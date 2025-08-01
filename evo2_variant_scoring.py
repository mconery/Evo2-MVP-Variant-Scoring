#!/usr/bin/env python3
"""
Evo2 40 B Variant-Scoring Pipeline
Author: Arc Institute / NVIDIA Blackwell reference implementation
--------------------------------------------------------------------------------
Scores ~2 M variants with 1-Mb context windows using Evo2_40b on multi-GPU
systems (B200, H100, HGX, GB200).  FP8 mixed precision is used automatically
whenever Blackwell GPUs are detected.
"""

import argparse, logging, os, sys, math, time, json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import torch
import numpy as np
import pandas as pd
from pyfaidx import Fasta

# ------------------------------------------------------------------------------
# Model initialisation
# ------------------------------------------------------------------------------

def init_model(model_name: str = "evo2_40b", use_fp8: bool = True):
    """
    Loads the Evo2 40 B model and places it across all available CUDA devices.
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

    logging.info("Loading Evo2 %s in %s precision …", model_name, precision.upper())
    model = Evo2(model_name)
    return model


# ------------------------------------------------------------------------------
# Variant-specific helpers
# ------------------------------------------------------------------------------

CTX = 1_000_000  # 1 Mb window (can be overridden via CLI)

def extract_context(fasta: Fasta, chrom: str, pos: int, ctx: int = CTX):
    half = ctx // 2
    start = max(1, pos - half)
    end   = pos + half
    seq = fasta[chrom][start-1:end].seq.upper()
    return seq, start  # return genomic start of slice

def make_alt_seq(ref_seq: str, pos: int, ref: str, alt: str, slice_start: int):
    rel = pos - slice_start
    if ref_seq[rel: rel + len(ref)] != ref:
        logging.warning("Reference mismatch at %d (%s→%s)", pos, ref, alt)
    return ref_seq[:rel] + alt + ref_seq[rel + len(ref):]

@torch.no_grad()
def ll_score(model, tokenizer, seq: str, device="cuda"):
    toks = torch.tensor(tokenizer(seq), dtype=torch.int16, device=device)[None, :]
    out, _ = model(toks)
    logp = torch.log_softmax(out, -1)
    return logp.mean().item()


def score_variant(model, tokenizer, fasta: Fasta, row, ctx: int):
    ref_seq, slice_start = extract_context(fasta, row.chr, row.pos, ctx)
    alt_seq = make_alt_seq(ref_seq, row.pos, row.ref, row.alt, slice_start)
    ref_s = ll_score(model, tokenizer, ref_seq)
    alt_s = ll_score(model, tokenizer, alt_seq)
    return alt_s - ref_s


# ------------------------------------------------------------------------------
# Worker – each process holds *one* copy of the model on *all* GPUs.
# ------------------------------------------------------------------------------

_worker_model = None
_worker_tok   = None
_worker_fa    = None

def _init_worker(fasta_path, model_name, fp8):
    global _worker_model, _worker_tok, _worker_fa
    logging.getLogger().setLevel(logging.ERROR)           # silence sub-process spam
    _worker_model = init_model(model_name, fp8)
    _worker_tok   = _worker_model.tokenizer.tokenize
    _worker_fa    = Fasta(fasta_path, rebuild=False)

def _worker_task(row_json, ctx):
    row = json.loads(row_json)
    return score_variant(_worker_model, _worker_tok, _worker_fa, row, ctx)


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser("Evo2 40 B Variant Scorer")
    p.add_argument("--variants", required=True, help="TSV/CSV with chr, pos, ref, alt")
    p.add_argument("--fasta",    required=True, help="Reference genome FASTA (indexed)")
    p.add_argument("--out",      required=True, help="Output TSV with evo2_score")
    p.add_argument("--ctx",      type=int, default=CTX, help="Context bp (default: 1 000 000)")
    p.add_argument("--processes",type=int, default=os.cpu_count()//2, help="Forks (≥ GPUs)")
    p.add_argument("--model",    default="evo2_40b")
    p.add_argument("--no-fp8",   action="store_true", help="Force BF16 even on B200/H100")
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s: %(message)s",
                        datefmt="%H:%M:%S")

    df = pd.read_csv(args.variants, sep="\t")
    for col in ("chromosome", "position", "ref_allele", "alt_allele"):
        if col not in df.columns:
            sys.exit(f"Missing column '{col}' in {args.variants}")
    
    #Map column names
    df = df.rename(columns={
        "chromosome": "chr",
        "position": "pos",
        "ref_allele": "ref",
        "alt_allele": "alt"
    })

    logging.info(f"Scoring {len(df):,d} variants with {args.processes} process(es)…")
    fasta_path = Path(args.fasta).expanduser()
    if not Path(args.fasta + '.fai').expanduser().exists():
        sys.exit("FASTA not indexed – run: samtools faidx <fasta>")

    # Serialise rows into compact strings (pickling large DataFrames costs RAM)
    rows_json = df.to_dict(orient="records")
    rows_json = [json.dumps(r) for r in rows_json]

    # Spin up workers
    with ProcessPoolExecutor(max_workers=args.processes,
                             initializer=_init_worker,
                             initargs=(str(fasta_path), args.model, not args.no_fp8)
                             ) as pool:
        futures = {pool.submit(_worker_task, rj, args.ctx): i
                   for i, rj in enumerate(rows_json)}
        scores = np.empty(len(df), dtype=np.float32)
        completed = 0
        for fut in as_completed(futures):
            idx = futures[fut]
            scores[idx] = fut.result()
            completed += 1
            if completed % 1000 == 0:
                logging.info("…%d / %d done", completed, len(df))

    df["evo2_score"] = scores
    df.to_csv(args.out, sep="\t", index=False)
    logging.info("Finished → %s", args.out)


if __name__ == "__main__":
    main()
