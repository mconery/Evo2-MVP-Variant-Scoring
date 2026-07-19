"""
Step 3a: Evo2 variant scorer for T2D fine-mapping loci.

Adapted from MVP_Low-High_Variant_Test/mvp_variants_test.py with these changes:
  - Reads per-locus TSV produced by extract_locus_variants.py
    (columns: SNP_ID, CHR, POS, REF, ALT, ...)
  - Uses GRCh37/hg19 FASTA (chromosomes named chr1, chr2, ...)
  - Context window fixed at 16,384 bp
  - No PIP-class logic; scores all variants and writes delta scores

Two modes:
  prepare  -- extract sequences, write temp FASTAs, print shell variables
  process  -- load .pt predictions, compute deltas, append to output CSV
"""

import glob
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import argparse

from bionemo.core.data.load import load
from pyfaidx import Fasta


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def extract_context(fasta: Fasta, chrom: str, pos: int, ctx: int = 16_384):
    half = ctx // 2
    start = max(1, pos - half)
    end = pos + half
    seq = fasta[chrom][start - 1:end].seq.upper()
    return seq, start


def make_alt_seq(ref_seq: str, pos: int, ref: str, alt: str, slice_start: int):
    rel = pos - slice_start
    if ref_seq[rel: rel + len(ref)] != ref:
        print(f"Reference mismatch at {pos} ({ref}→{alt})", file=sys.stderr)
    return ref_seq[:rel] + alt + ref_seq[rel + len(ref):]


def extract_variant_sequences(row, fasta, ctx):
    try:
        chrom_key = "chr" + str(row["CHR"])
        ref_seq, slice_start = extract_context(fasta, chrom_key, int(row["POS"]), ctx)
        alt_seq = make_alt_seq(ref_seq, int(row["POS"]), row["REF"], row["ALT"], slice_start)
        return (row.name, (ref_seq, alt_seq))
    except Exception as e:
        print(f"Error at {row['CHR']}:{row['POS']}: {e}", file=sys.stderr)
        return (row.name, (None, None))


def check_fp8_support():
    if not torch.cuda.is_available():
        return False, "CUDA not available"
    p = torch.cuda.get_device_properties(0)
    supported = (p.major > 8) or (p.major == 8 and p.minor >= 9)
    return supported, f"Device: {p.name}, CC: {p.major}.{p.minor}"


def load_variant_df(variant_file):
    df = pd.read_csv(variant_file, sep="\t", dtype={"CHR": str})
    df.index = df["SNP_ID"]
    return df


def get_checkpoint_path(model_size, checkpoint_dir):
    if model_size == "7b_arc_longcontext":
        return load("evo2/7b-1m:1.0")
    elif model_size == "40b_arc_longcontext":
        return load("evo2/40b-1m-fp8-bf16:1.0")
    elif model_size == "1b":
        return load("evo2/1b-8k-bf16:1.0")
    else:
        ckpt = Path(checkpoint_dir) / f"nemo2_evo2_{model_size}_8k"
        if not ckpt.exists() or not any(ckpt.iterdir()):
            os.system(
                f"evo2_convert_to_nemo2 --model-path hf://arcinstitute/savanna_evo2_{model_size}_base "
                f"--model-size {model_size} --output-dir {ckpt}"
            )
        return ckpt


# ------------------------------------------------------------------
# prepare mode
# ------------------------------------------------------------------

def mode_prepare(args):
    variant_file = args.variants
    fasta_path = args.fasta
    out_file = args.out
    chunk_size = args.chunk_size
    ctx = args.ctx
    model_size = args.model
    checkpoint_dir = args.checkpoint_dir
    tp_size = args.tensor_parallel_size
    cp_size = args.context_parallel_size
    chunk_start_arg = args.chunk_start

    if not Path(str(fasta_path) + ".fai").exists():
        sys.exit(f"FASTA not indexed: {fasta_path}")

    variant_df = load_variant_df(variant_file)
    num_variants = len(variant_df)
    print(f"Processing {num_variants} variants ...", file=sys.stderr)

    checkpoint_path = get_checkpoint_path(model_size, checkpoint_dir)
    fp8_supported, gpu_info = check_fp8_support()
    fp8_option = "--fp8" if fp8_supported else ""
    print(f"FP8 support: {fp8_supported}  ({gpu_info})", file=sys.stderr)

    # Check if already complete
    if os.path.exists(out_file):
        try:
            existing = pd.read_csv(out_file)
            if len(existing) >= num_variants:
                print("Output already complete. Exiting.", file=sys.stderr)
                sys.exit(1)
        except Exception:
            pass

    # Determine resume point
    resume_from = 0
    if chunk_start_arg == 0 and os.path.exists(out_file):
        try:
            existing = pd.read_csv(out_file)
            n_done = len(existing)
            if n_done > 0:
                resume_from = (n_done // chunk_size) * chunk_size
                if resume_from > 0 and resume_from < n_done:
                    existing.iloc[:resume_from].to_csv(out_file, index=False)
                if resume_from > 0:
                    print(f"Resuming from variant {resume_from}", file=sys.stderr)
        except Exception:
            pass

    actual_chunk_start = resume_from if chunk_start_arg == 0 else chunk_start_arg
    if actual_chunk_start >= num_variants:
        sys.exit(1)

    chunk_end = min(actual_chunk_start + chunk_size, num_variants)

    # Probe call: just return variables without writing FASTAs
    if chunk_start_arg == 0 and resume_from > 0:
        print(f"FP8_FLAG={fp8_option}")
        print(f"CHUNK_END={chunk_end}")
        print(f"RESUME_FROM={actual_chunk_start}")
        print(f"CHECKPOINT_PATH={checkpoint_path}")
        sys.exit(0)

    # Extract sequences for this chunk
    chunk_rows = variant_df.iloc[actual_chunk_start:chunk_end]
    fasta = Fasta(str(fasta_path))

    with ThreadPoolExecutor(max_workers=os.cpu_count()) as ex:
        futures = {
            ex.submit(extract_variant_sequences, row, fasta, ctx): row.name
            for _, row in chunk_rows.iterrows()
        }
        results = {}
        for future in as_completed(futures):
            idx, (ref_seq, alt_seq) = future.result()
            results[idx] = (ref_seq, alt_seq)

    # Truncate for parallelism alignment
    divisor = 1
    if cp_size > 1:
        divisor *= 2 * cp_size
    if tp_size > 1:
        divisor *= tp_size
    if fp8_option:
        divisor = max(divisor, 8 * cp_size * tp_size)

    if divisor > 1:
        for v in results:
            r, a = results[v]
            if r is not None:
                tlen = (len(r) // divisor) * divisor
                if tlen > 0:
                    results[v] = (r[:tlen], a[:tlen])

    # Write FASTA files
    output_dir = Path(os.path.dirname(out_file))
    ref_fa = output_dir / f"temp_ref.{model_size}.{ctx}bp.fa"
    var_fa = output_dir / f"temp_var.{model_size}.{ctx}bp.fa"

    with open(ref_fa, "w") as f:
        for v in results:
            f.write(f">{v}_ref\n{results[v][0]}\n")
    with open(var_fa, "w") as f:
        for v in results:
            f.write(f">{v}_alt\n{results[v][1]}\n")

    (output_dir / f"reference_predictions.{model_size}.{ctx}bp").mkdir(parents=True, exist_ok=True)
    (output_dir / f"variant_predictions.{model_size}.{ctx}bp").mkdir(parents=True, exist_ok=True)

    print(f"FP8_FLAG={fp8_option}")
    print(f"CHUNK_END={chunk_end}")
    print(f"RESUME_FROM={actual_chunk_start}")
    print(f"CHECKPOINT_PATH={checkpoint_path}")
    sys.exit(0)


# ------------------------------------------------------------------
# process mode
# ------------------------------------------------------------------

def mode_process(args):
    variant_file = args.variants
    out_file = args.out
    chunk_size = args.chunk_size
    ctx = args.ctx
    model_size = args.model
    chunk_start = args.chunk_start
    chunk_end = args.chunk_end

    if chunk_end is None:
        sys.exit("--chunk-end is required for process mode")

    variant_df = load_variant_df(variant_file)
    num_variants = len(variant_df)

    chunk_df = variant_df.iloc[chunk_start:chunk_end].copy()
    chunk_df["ref_fasta_name"] = [f"{v}_ref" for v in chunk_df.index]
    chunk_df["var_fasta_name"] = [f"{v}_alt" for v in chunk_df.index]

    output_dir = Path(os.path.dirname(out_file))
    predict_ref_dir = output_dir / f"reference_predictions.{model_size}.{ctx}bp"
    predict_var_dir = output_dir / f"variant_predictions.{model_size}.{ctx}bp"
    ref_fa = output_dir / f"temp_ref.{model_size}.{ctx}bp.fa"
    var_fa = output_dir / f"temp_var.{model_size}.{ctx}bp.fa"

    ref_files = sorted(glob.glob(str(predict_ref_dir / "predictions__rank_*.pt")))
    var_files = sorted(glob.glob(str(predict_var_dir / "predictions__rank_*.pt")))
    if not ref_files:
        sys.exit(f"No prediction files in {predict_ref_dir}")
    if not var_files:
        sys.exit(f"No prediction files in {predict_var_dir}")

    with open(predict_ref_dir / "seq_idx_map.json") as f:
        ref_idx_map = json.load(f)
    with open(predict_var_dir / "seq_idx_map.json") as f:
        var_idx_map = json.load(f)

    ref_preds = torch.load(ref_files[0])
    var_preds = torch.load(var_files[0])

    ref_lps, var_lps = [], []
    for _, row in chunk_df.iterrows():
        ref_lps.append(ref_preds["log_probs_seqs"][ref_idx_map[row["ref_fasta_name"]]].item())
        var_lps.append(var_preds["log_probs_seqs"][var_idx_map[row["var_fasta_name"]]].item())

    chunk_df["ref_log_probs"] = ref_lps
    chunk_df["var_log_probs"] = var_lps
    chunk_df["evo2_delta_score"] = chunk_df["var_log_probs"] - chunk_df["ref_log_probs"]

    out_cols = ["SNP_ID", "CHR", "POS", "REF", "ALT", "ref_log_probs", "var_log_probs", "evo2_delta_score"]
    write_header = True
    if os.path.exists(out_file):
        try:
            if len(pd.read_csv(out_file)) > 0:
                write_header = False
        except Exception:
            pass
    mode = "w" if write_header else "a"
    chunk_df[out_cols].to_csv(out_file, mode=mode, header=write_header, index=False)
    print(f"Saved variants {chunk_start + 1}–{chunk_end} of {num_variants}", file=sys.stderr)

    # Cleanup prediction files
    for f in glob.glob(str(predict_ref_dir / "predictions__rank_*.pt")):
        os.remove(f)
    for f in glob.glob(str(predict_var_dir / "predictions__rank_*.pt")):
        os.remove(f)
    for p in [predict_ref_dir / "seq_idx_map.json", predict_var_dir / "seq_idx_map.json"]:
        if p.exists():
            p.unlink()
    if chunk_end >= num_variants:
        for fa in [ref_fa, var_fa]:
            if fa.exists():
                fa.unlink()

    sys.exit(0)


# ------------------------------------------------------------------
# main
# ------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser("Evo2 T2D locus variant scorer")
    parser.add_argument("--mode", required=True, choices=["prepare", "process"])
    parser.add_argument("--variants", required=True, help="Per-locus variant TSV")
    parser.add_argument("--fasta", required=True, help="GRCh37 reference FASTA (must have .fai index)")
    parser.add_argument("--out", required=True, help="Output scores CSV")
    parser.add_argument("--ctx", type=int, default=16_384, help="Context window in bp (default: 16384)")
    parser.add_argument("--chunk-size", type=int, default=10,
                        help="Variants per scoring batch (default: 10; increase if GPU memory allows)")
    parser.add_argument("--model", default="7b_arc_longcontext")
    parser.add_argument("--checkpoint-dir",
                        default="/lus/grand/projects/GeomicVar/mconery/tools/bionemo_cache")
    parser.add_argument("--tensor-parallel-size", type=int, default=4)
    parser.add_argument("--context-parallel-size", type=int, default=1)
    parser.add_argument("--chunk-start", type=int, default=0)
    parser.add_argument("--chunk-end", type=int, default=None)
    args = parser.parse_args()

    if args.mode == "prepare":
        mode_prepare(args)
    else:
        mode_process(args)


if __name__ == "__main__":
    main()
