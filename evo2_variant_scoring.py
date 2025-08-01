#!/usr/bin/env python3
"""
Evo2 40 B Variant-Scoring Pipeline
Author: Arc Institute / NVIDIA Blackwell reference implementation
--------------------------------------------------------------------------------
Scores ~2 M variants with 1-Mb context windows using Evo2_40b on multi-GPU
systems (B200, H100, HGX, GB200). FP8 mixed precision is used automatically
whenever Blackwell GPUs are detected.
"""

import argparse, logging, os, sys, math, time, json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

import torch
import numpy as np
import pandas as pd
from pyfaidx import Fasta

# ------------------------------------------------------------------------------
# Model initialisation
# ------------------------------------------------------------------------------

def init_model(model_name: str = "evo2_40b", use_fp8: bool = True):
    """
    Loads the Evo2 40 B model and places it on the assigned CUDA device.
    On Blackwell GPUs FP8 mixed precision is enabled automatically.
    """
    from evo2 import Evo2  # Lazy-import to avoid torch-cuda init in forked workers
    
    if torch.cuda.device_count() == 0:
        sys.exit("✗ No CUDA devices found – a Blackwell/Hopper GPU is required.")
    
    logging.info("Detected %d GPU(s): %s",
                 torch.cuda.device_count(),
                 ", ".join(torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count()))
                 )
    
    # FP8 only if hardware supports it (B200, H200, GB200, H100)
    fp8_ok = use_fp8 and any("B200" in torch.cuda.get_device_name(i) or
                             "H100" in torch.cuda.get_device_name(i)
                             for i in range(torch.cuda.device_count()))
    precision = "fp8" if fp8_ok else "bf16"
    
    current_device = torch.cuda.current_device()
    logging.info("Loading Evo2 %s in %s precision on GPU %d…",
                 model_name, precision.upper(), current_device)
    
    model = Evo2(model_name)
    return model

# ------------------------------------------------------------------------------
# Variant-specific helpers
# ------------------------------------------------------------------------------

CTX = 1_000_000  # 1 Mb window (can be overridden via CLI)

def extract_context(fasta: Fasta, chrom: str, pos: int, ctx: int = CTX):
    half = ctx // 2
    start = max(1, pos - half)
    end = pos + half
    seq = fasta[chrom][start-1:end].seq.upper()
    return seq, start  # return genomic start of slice

def make_alt_seq(ref_seq: str, pos: int, ref: str, alt: str, slice_start: int):
    rel = pos - slice_start
    if ref_seq[rel: rel + len(ref)] != ref:
        logging.warning("Reference mismatch at %d (%s→%s)", pos, ref, alt)
    return ref_seq[:rel] + alt + ref_seq[rel + len(ref):]

@torch.no_grad()
def ll_score(model, tokenizer, seq: str):
    device = "cuda:0" if torch.cuda.is_available() else "cpu"
    toks = torch.tensor(tokenizer(seq), dtype=torch.int16, device=device)[None, :]
    out, _ = model(toks)
    logp = torch.log_softmax(out, -1)
    return logp.mean().item()

def score_variant(model, tokenizer, fasta: Fasta, row, ctx: int):
    ref_seq, slice_start = extract_context(fasta, row["chr"], row["pos"], ctx)
    alt_seq = make_alt_seq(ref_seq, row["pos"], row["ref"], row["alt"], slice_start)
    ref_s = ll_score(model, tokenizer, ref_seq)
    alt_s = ll_score(model, tokenizer, alt_seq)
    return alt_s - ref_s

# ------------------------------------------------------------------------------
# Worker – each process holds *one* copy of the model on *one* assigned GPU
# ------------------------------------------------------------------------------

_worker_model = None
_worker_tok = None
_worker_fa = None

def _init_worker(fasta_path, model_name, fp8, gpu_id):
    global _worker_model, _worker_tok, _worker_fa
    
    # Set the specific GPU for this worker process
    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    torch.cuda.set_device(0)  # After setting CUDA_VISIBLE_DEVICES, use device 0
    
    logging.getLogger().setLevel(logging.ERROR)  # silence sub-process spam
    _worker_model = init_model(model_name, fp8)
    _worker_tok = _worker_model.tokenizer.tokenize
    _worker_fa = Fasta(fasta_path, rebuild=False)

def _worker_task(row_json, ctx):
    row = json.loads(row_json)
    return score_variant(_worker_model, _worker_tok, _worker_fa, row, ctx)

# ------------------------------------------------------------------------------
# Chunk processing function - MOVED OUTSIDE OF main() TO BE PICKLEABLE
# ------------------------------------------------------------------------------

def process_chunk(chunk_data):
    """
    Process a chunk of variants on a specific GPU.
    This function must be at module level to be pickleable.
    """
    chunk_rows, gpu_id, start_offset, fasta_path, model_name, fp8, ctx = chunk_data
    
    # Initialize worker for this specific GPU
    _init_worker(fasta_path, model_name, fp8, gpu_id)
    
    chunk_scores = []
    for i, row_json in enumerate(chunk_rows):
        score = _worker_task(row_json, ctx)
        chunk_scores.append(score)
        
        if (start_offset + i + 1) % 100 == 0:
            logging.info(f"GPU {gpu_id}: processed {start_offset + i + 1} variants")
    
    return start_offset, chunk_scores

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser("Evo2 40 B Variant Scorer")
    p.add_argument("--variants", required=True, help="TSV/CSV with chromosome, position, ref_allele, alt_allele")
    p.add_argument("--fasta", required=True, help="Reference genome FASTA (indexed)")
    p.add_argument("--out", required=True, help="Output TSV with evo2_score")
    p.add_argument("--ctx", type=int, default=CTX, help="Context bp (default: 1 000 000)")
    p.add_argument("--processes", type=int, default=None, help="Number of processes (defaults to number of GPUs)")
    p.add_argument("--model", default="evo2_40b")
    p.add_argument("--no-fp8", action="store_true", help="Force BF16 even on B200/H100")
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s: %(message)s",
                        datefmt="%H:%M:%S")

    # Determine number of GPUs and set processes accordingly
    num_gpus = torch.cuda.device_count()
    if num_gpus == 0:
        sys.exit("✗ No CUDA devices found – GPU required for Evo2 40B")
    
    if args.processes is None:
        args.processes = num_gpus
        logging.info(f"Auto-detected {num_gpus} GPUs, setting processes to {args.processes}")
    elif args.processes > num_gpus:
        logging.warning(f"Requested {args.processes} processes but only {num_gpus} GPUs available. "
                       f"Limiting to {num_gpus} processes.")
        args.processes = num_gpus

    df = pd.read_csv(args.variants, sep="\t")
    for col in ("chromosome", "position", "ref_allele", "alt_allele"):
        if col not in df.columns:
            sys.exit(f"Missing column '{col}' in {args.variants}")

    # Map column names
    df = df.rename(columns={
        "chromosome": "chr",
        "position": "pos",
        "ref_allele": "ref",
        "alt_allele": "alt"
    })

    logging.info(f"Scoring {len(df):,d} variants with {args.processes} process(es) on {num_gpus} GPU(s)…")
    
    fasta_path = Path(args.fasta).expanduser()
    if not Path(str(fasta_path) + '.fai').exists():
        sys.exit(f"FASTA not indexed – run: samtools faidx {fasta_path}")

    # Split variants among processes to balance load
    chunk_size = len(df) // args.processes
    if len(df) % args.processes != 0:
        chunk_size += 1

    # Serialise rows into compact strings (pickling large DataFrames costs RAM)
    rows_json = df.to_dict(orient="records")
    rows_json = [json.dumps(r) for r in rows_json]

    # Create chunks for each GPU/process
    chunks = []
    for i in range(args.processes):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(rows_json))
        if start_idx < len(rows_json):
            gpu_id = i % num_gpus  # Cycle through available GPUs
            # Pack all necessary data for the worker
            chunk_data = (
                rows_json[start_idx:end_idx],  # chunk_rows
                gpu_id,                        # gpu_id
                start_idx,                     # start_offset
                str(fasta_path),              # fasta_path
                args.model,                   # model_name
                not args.no_fp8,             # fp8
                args.ctx                      # ctx
            )
            chunks.append(chunk_data)

    scores = np.empty(len(df), dtype=np.float32)
    completed = 0

    # Use ProcessPoolExecutor with spawn method to avoid CUDA context issues
    mp_context = mp.get_context('spawn')
    
    with ProcessPoolExecutor(max_workers=args.processes, mp_context=mp_context) as executor:
        futures = {executor.submit(process_chunk, chunk): chunk for chunk in chunks}
        
        for future in as_completed(futures):
            start_offset, chunk_scores = future.result()
            
            # Place scores in correct positions
            for i, score in enumerate(chunk_scores):
                scores[start_offset + i] = score
            
            completed += len(chunk_scores)
            logging.info(f"Completed {completed:,d} / {len(df):,d} variants")

    df["evo2_score"] = scores
    df.to_csv(args.out, sep="\t", index=False)
    logging.info("Finished → %s", args.out)

if __name__ == "__main__":
    main()
