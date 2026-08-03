"""
Evo2 variant scorer for T2D fine-mapping loci — runs on PARCC Betty.

Adapted from MVP_Low-High_Variant_Test/mvp_variants_test.py with these changes:
  - Reads the merged variant TSV produced by merge_variant_lists.py
    (columns: uid, locus_id, SNP_ID, CHR, POS, REF, ALT, ...), indexed on the
    collision-safe `uid` column so scoring/fasta-naming never gets two loci
    that happen to share an rsID
  - Uses GRCh37/hg19 FASTA (chromosomes named chr1, chr2, ...)
  - No PIP-class logic; scores all variants and writes a restricted set of
    delta-score columns, carrying `locus_id` through to the output so
    results can be split back into per-locus files after transfer back to
    Polaris

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

def extract_context(fasta: Fasta, chrom: str, pos: int, ctx: int = 8_192):
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
    """
    Helper worker function to extract ref and alt sequences for a single variant.
    Returns: (index, (ref_seq, alt_seq)) or (index, (None, None)) on error
    """
    try:
        ref_seq, slice_start = extract_context(fasta, "chr" + str(row["CHR"]), row["POS"], ctx)
        alt_seq = make_alt_seq(ref_seq, row["POS"], row["REF"], row["ALT"], slice_start)
        return (row.name, (ref_seq, alt_seq))
    except Exception as e:
        print(f"Error extracting sequences for variant {row['CHR']}:{row['POS']}: {e}", file=sys.stderr)
        return (row.name, (None, None))


def check_fp8_support():
    """Check if FP8 is supported on the current GPU.

    FP8 requires compute capability 8.9+ (Ada Lovelace/Hopper architecture or newer).
    """
    if not torch.cuda.is_available():
        return False, "CUDA not available"

    device_props = torch.cuda.get_device_properties(0)
    compute_capability = f"{device_props.major}.{device_props.minor}"
    device_name = device_props.name

    # FP8 is supported on compute capability 8.9+ (Ada Lovelace/Hopper architecture)
    is_supported = (device_props.major > 8) or (device_props.major == 8 and device_props.minor >= 9)

    return is_supported, f"Device: {device_name}, Compute Capability: {compute_capability}"


def load_variant_df(variant_file):
    df = pd.read_csv(variant_file, sep="\t", dtype={"CHR": str})
    df.index = df["uid"]
    return df


def get_checkpoint_path(model_size, checkpoint_dir):
    """Resolve or download checkpoint; return Path."""
    if model_size == "1b":
        return load("evo2/1b-8k-bf16:1.0")
    elif model_size == "40b_arc_longcontext":
        return load("evo2/40b-1m-fp8-bf16:1.0")
    elif model_size == "7b_arc_longcontext":
        return load("evo2/7b-1m:1.0")
    else:
        ckpt = Path(checkpoint_dir) / f"nemo2_evo2_{model_size}_8k"
        if not ckpt.exists() or not any(ckpt.iterdir()):
            os.system(
                f"evo2_convert_to_nemo2 --model-path hf://arcinstitute/savanna_evo2_{model_size}_base "
                f"--model-size {model_size} --output-dir {ckpt}"
            )
        else:
            print("Checkpoint directory is not empty. Skipping download.", file=sys.stderr)
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

    # Check FASTA index
    if not Path(str(fasta_path) + ".fai").exists():
        sys.exit(f"FASTA not indexed: {fasta_path}")

    # Load variant data
    variant_df = load_variant_df(variant_file)
    num_variants = len(variant_df)
    print(f"Processing {num_variants} variants...", file=sys.stderr)

    # Resolve/download checkpoint
    checkpoint_path = get_checkpoint_path(model_size, checkpoint_dir)

    # Check FP8 support
    fp8_supported, gpu_info = check_fp8_support()
    fp8_option = "--fp8" if fp8_supported else ""
    print(f"FP8 Support: {fp8_supported}", file=sys.stderr)
    print(gpu_info, file=sys.stderr)

    # Check if output is already complete
    if os.path.exists(out_file):
        try:
            existing_df = pd.read_csv(out_file)
            if len(existing_df) >= num_variants:
                print(f"Output file already complete ({len(existing_df)} variants). Exiting.", file=sys.stderr)
                sys.exit(1)
        except Exception:
            pass

    # Determine effective chunk_start (handle resume when called with chunk_start=0)
    resume_from = 0
    if chunk_start_arg == 0 and os.path.exists(out_file):
        try:
            existing_df = pd.read_csv(out_file)
            n_done = len(existing_df)
            if n_done > 0:
                resume_from = (n_done // chunk_size) * chunk_size
                if resume_from > 0:
                    # Truncate to full-chunk boundary in case a partial chunk was written
                    if resume_from < n_done:
                        existing_df.iloc[:resume_from].to_csv(out_file, index=False)
                    print(f"Resuming from variant {resume_from} of {num_variants}", file=sys.stderr)
        except Exception:
            pass  # Corrupt/empty file: start fresh

    actual_chunk_start = resume_from if chunk_start_arg == 0 else chunk_start_arg

    # Check if all chunks are done
    if actual_chunk_start >= num_variants:
        sys.exit(1)

    chunk_end = min(actual_chunk_start + chunk_size, num_variants)

    # If this is a probe call (bash passed chunk_start=0 but resume detected),
    # output variables without writing FASTAs. Bash will re-call with actual start.
    if chunk_start_arg == 0 and resume_from > 0:
        print(f"FP8_FLAG={fp8_option}")
        print(f"CHUNK_END={chunk_end}")
        print(f"RESUME_FROM={actual_chunk_start}")
        print(f"CHECKPOINT_PATH={checkpoint_path}")
        sys.exit(0)

    # --- Extract sequences for this chunk ---
    chunk_rows = variant_df.iloc[actual_chunk_start:chunk_end]
    fasta = Fasta(str(fasta_path))

    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = {
            executor.submit(extract_variant_sequences, row, fasta, ctx): row.name
            for _, row in chunk_rows.iterrows()
        }
        results_by_index = {}
        for future in as_completed(futures):
            idx, (ref_seq, alt_seq) = future.result()
            results_by_index[idx] = (ref_seq, alt_seq)

    # Truncate sequences for context parallelism compatibility.
    # extract_context produces 2*(ctx//2)+1 bases (always odd).
    # Megatron's get_batch_on_this_cp_rank requires seq length divisible by 2*cp_size*tp_size.
    # Additionally, FP8 GEMM (TransformerEngine) requires each GPU's local sequence shard
    # to be divisible by 8. The Hyena dense_projection splits the sequence across
    # cp_size * tp_size GPUs, so the full sequence must be divisible by 8 * cp_size * tp_size.
    divisor = 1
    if cp_size > 1:
        divisor *= 2 * cp_size     # Megatron CP requirement
    if tp_size > 1:
        divisor *= tp_size         # ensure local dim is divisible by TP
    if fp8_option:                 # FP8 GEMM: each rank's shard must be divisible by 8
        divisor = max(divisor, 8 * cp_size * tp_size)

    if divisor > 1:
        for variant in results_by_index:
            ref_seq, alt_seq = results_by_index[variant]
            if ref_seq is not None:
                trunc_len = (len(ref_seq) // divisor) * divisor
                if trunc_len == 0:
                    continue
                results_by_index[variant] = (ref_seq[:trunc_len], alt_seq[:trunc_len])

    # Write FASTA files
    output_dir = Path(os.path.dirname(out_file))
    ref_fasta_path = output_dir / f"temp_ref.{model_size}.{ctx}bp.fa"
    var_fasta_path = output_dir / f"temp_var.{model_size}.{ctx}bp.fa"

    ref_entries = []
    var_entries = []
    for variant in results_by_index.keys():
        ref_entries.append(f">{variant}_ref\n{results_by_index[variant][0]}\n")
        var_entries.append(f">{variant}_alt\n{results_by_index[variant][1]}\n")
    with open(ref_fasta_path, "w") as f:
        f.writelines(ref_entries)
    with open(var_fasta_path, "w") as f:
        f.writelines(var_entries)

    # Create prediction output dirs
    predict_ref_dir = output_dir / f"reference_predictions.{model_size}.{ctx}bp"
    predict_var_dir = output_dir / f"variant_predictions.{model_size}.{ctx}bp"
    predict_ref_dir.mkdir(parents=True, exist_ok=True)
    predict_var_dir.mkdir(parents=True, exist_ok=True)

    # Print metadata for bash to consume (stdout only — captured by bash $(...))
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

    # Load variant data
    variant_df = load_variant_df(variant_file)
    num_variants = len(variant_df)

    # Get the chunk and build fasta name columns
    chunk_df = variant_df.iloc[chunk_start:chunk_end].copy()
    chunk_df["ref_fasta_name"] = [f"{v}_ref" for v in chunk_df.index]
    chunk_df["var_fasta_name"] = [f"{v}_alt" for v in chunk_df.index]

    # Reconstruct prediction paths (deterministic naming, no state files needed)
    output_dir = Path(os.path.dirname(out_file))
    predict_ref_dir = output_dir / f"reference_predictions.{model_size}.{ctx}bp"
    predict_var_dir = output_dir / f"variant_predictions.{model_size}.{ctx}bp"
    ref_fasta_path = output_dir / f"temp_ref.{model_size}.{ctx}bp.fa"
    var_fasta_path = output_dir / f"temp_var.{model_size}.{ctx}bp.fa"

    # Find and load prediction files (rank 0 holds the full results)
    ref_pred_files = sorted(glob.glob(str(predict_ref_dir / "predictions__rank_*.pt")))
    var_pred_files = sorted(glob.glob(str(predict_var_dir / "predictions__rank_*.pt")))

    if not ref_pred_files:
        sys.exit(f"No prediction files found in {predict_ref_dir}")
    if not var_pred_files:
        sys.exit(f"No prediction files found in {predict_var_dir}")

    # Load sequence ID maps
    with open(predict_ref_dir / "seq_idx_map.json") as f:
        ref_seq_idx_map = json.load(f)
    with open(predict_var_dir / "seq_idx_map.json") as f:
        var_seq_idx_map = json.load(f)

    # Load predictions (rank 0 file has all sequence log-probs)
    ref_preds = torch.load(ref_pred_files[0])
    var_preds = torch.load(var_pred_files[0])

    # Calculate change in log-likelihoods
    ref_log_probs = []
    var_log_probs = []
    for _, row in chunk_df.iterrows():
        ref_name = row["ref_fasta_name"]
        var_name = row["var_fasta_name"]
        ref_log_probs.append(ref_preds["log_probs_seqs"][ref_seq_idx_map[ref_name]].item())
        var_log_probs.append(var_preds["log_probs_seqs"][var_seq_idx_map[var_name]].item())

    chunk_df["ref_log_probs"] = ref_log_probs
    chunk_df["var_log_probs"] = var_log_probs
    # ideally probability of a broken variant is lower than a good one. So a bad var - good ref is negative.
    chunk_df["evo2_delta_score"] = chunk_df["var_log_probs"] - chunk_df["ref_log_probs"]

    out_cols = ["locus_id", "SNP_ID", "CHR", "POS", "REF", "ALT", "ref_log_probs", "var_log_probs", "evo2_delta_score"]

    # Append to output CSV (write_header only if file has no data yet)
    write_header = True
    if os.path.exists(out_file):
        try:
            existing_df = pd.read_csv(out_file)
            if len(existing_df) > 0:
                write_header = False
        except Exception:
            pass
    mode = "w" if write_header else "a"
    chunk_df[out_cols].to_csv(out_file, mode=mode, header=write_header, index=False)
    print(f"Saved variants {chunk_start + 1}-{chunk_end} of {num_variants}", file=sys.stderr)

    # Clean up prediction files so stale data never affects subsequent chunks
    for f in glob.glob(str(predict_ref_dir / "predictions__rank_*.pt")):
        os.remove(f)
    seq_map_ref = predict_ref_dir / "seq_idx_map.json"
    if seq_map_ref.exists():
        seq_map_ref.unlink()
    for f in glob.glob(str(predict_var_dir / "predictions__rank_*.pt")):
        os.remove(f)
    seq_map_var = predict_var_dir / "seq_idx_map.json"
    if seq_map_var.exists():
        seq_map_var.unlink()

    # Remove temp FASTA files after the final chunk
    if chunk_end >= num_variants:
        if ref_fasta_path.exists():
            ref_fasta_path.unlink()
        if var_fasta_path.exists():
            var_fasta_path.unlink()

    sys.exit(0)


# ------------------------------------------------------------------
# main
# ------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser("Evo2 T2D locus variant scorer")
    parser.add_argument("--mode", required=True, choices=["prepare", "process"],
                        help="prepare: extract sequences and write FASTAs; process: load predictions and write CSV")
    parser.add_argument("--chunk-start", type=int, default=0, help="Starting variant index for this chunk")
    parser.add_argument("--chunk-end", type=int, default=None, help="Ending variant index (process mode only)")
    parser.add_argument("--variants", required=True, help="Merged variant TSV (see merge_variant_lists.py)")
    parser.add_argument("--fasta", required=True, help="GRCh37 reference FASTA (must have .fai index)")
    parser.add_argument("--out", required=True, help="Output scores CSV")
    parser.add_argument("--ctx", type=int, default=8_192, help="Context window in bp (default: 8192)")
    parser.add_argument("--chunk-size", type=int, default=100,
                        help="Number of variants per scoring batch (default: 100; increase if GPU memory allows)")
    parser.add_argument("--model", default="7b_arc_longcontext", help="Model name")
    parser.add_argument("--checkpoint-dir", default="/vast/projects/anuragv/cohort/mconery/bionemo",
                        help="Checkpoint Directory")
    parser.add_argument("--tensor-parallel-size", default=8, type=int, help="Tensor Parallel Size")
    parser.add_argument("--context-parallel-size", default=1, type=int, help="Context Parallel Size")
    args = parser.parse_args()

    if args.mode == "prepare":
        mode_prepare(args)
    elif args.mode == "process":
        mode_process(args)


if __name__ == "__main__":
    main()
