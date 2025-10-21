# Required imports
from Bio import SeqIO
import gzip
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns
from sklearn.metrics import roc_auc_score
from concurrent.futures import ThreadPoolExecutor, as_completed
from pyfaidx import Fasta
import argparse
from pathlib import Path
import torch
import sys
from bionemo.core.utils.subprocess_utils import run_subprocess_safely
import glob
import json

def extract_context(fasta: Fasta, chrom: str, pos: int, ctx: int = 1_000_000):
    half = ctx // 2
    start = max(1, pos - half)
    end = pos + half
    seq = fasta[chrom][start-1:end].seq.upper()
    return seq, start

def make_alt_seq(ref_seq: str, pos: int, ref: str, alt: str, slice_start: int):
    rel = pos - slice_start
    if ref_seq[rel: rel + len(ref)] != ref:
        print(f"Reference mismatch at {pos} ({ref}→{alt})")
    return ref_seq[:rel] + alt + ref_seq[rel + len(ref):]

def extract_variant_sequences(row, fasta, ctx):
    """
    Helper worker function to extract ref and alt sequences for a single variant.
    Returns: (index, (ref_seq, alt_seq)) or (index, (None, None)) on error
    """
    try:
        ref_seq, slice_start = extract_context(fasta, "chr" + str(row["CHR"]), row["BP38"], ctx)
        alt_seq = make_alt_seq(ref_seq, row["BP38"], row["REF"], row["ALT"], slice_start)
        return (row.name, (ref_seq, alt_seq))
    except Exception as e:
        print(f"Error extracting sequences for variant {row['CHR']}:{row['BP38']}: {e}")
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

def main():
    #Deploy argument parsing
    parser = argparse.ArgumentParser("Evo2 MVP Variant Scorer")
    parser.add_argument("--variants", required=False, help="Variant CSV file", default="/vast/projects/anuragv/cohort/mconery/mvp_variant_test/MVP_matched_variants.csv")
    parser.add_argument("--fasta", required=False, help="Reference FASTA file", default="/vast/projects/anuragv/cohort/mconery/genome_assemblies/GRCh38.primary_assembly.genome.fa")
    parser.add_argument("--out", required=True, help="Output file")
    parser.add_argument("--ctx", type=int, default=1_000_000, help="Context size (default: 1M)")
    parser.add_argument("--chunk-size", type=int, default=16, help="Number of variants per scoring batch (default: 16)")
    parser.add_argument("--model", default="7b", help="Model name")
    parser.add_argument("--checkpoint-dir", default="/vast/projects/anuragv/cohort/mconery/bionemo", help="Checkpoint Directory")
    parser.add_argument("--tensor-parallel-size", default=1, type=int, help="Tensor Parallel Size")
    parser.add_argument("--context-parallel-size", default=1, type=int, help="Context Parallel Size")
    args = parser.parse_args()
    #Copy args to variables to simplify debugging
    variant_file = args.variants
    fasta_path = args.fasta
    out_file = args.out
    chunk_size = args.chunk_size
    window_size = args.ctx
    MODEL_SIZE = args.model
    checkpoint_dir = args.checkpoint_dir
    tp_size = args.tensor_parallel_size
    cp_size = args.context_parallel_size
        
    # Define checkpoint path and download model weights if needed
    if MODEL_SIZE == "1b":
        from bionemo.core.data.load import load
        #  This line will download the checkpoint from NGC to your $HOME/.cache/bionemo directory and return the path.
        #  To do the same from the command line, use `CHECKPOINT_PATH=$(download_bionemo_data evo2/1b-8k-bf16:1.0)`
        checkpoint_path = load("evo2/1b-8k-bf16:1.0")
    else:
        checkpoint_path = Path(f"{checkpoint_dir}/nemo2_evo2_{MODEL_SIZE}_8k")
        # Check if the directory does not exist or is empty
        if not checkpoint_path.exists() or not any(checkpoint_path.iterdir()):
            os.system(f'evo2_convert_to_nemo2 --model-path hf://arcinstitute/savanna_evo2_{MODEL_SIZE}_base --model-size {MODEL_SIZE} --output-dir {checkpoint_dir}/nemo2_evo2_{MODEL_SIZE}_8k')
        else:
            print("Checkpoint directory is not empty. Skipping command.")
    
    #Check fp8 support
    fp8_supported, gpu_info = check_fp8_support()
    fp8_option = "--fp8" if fp8_supported else ""
    print(f"FP8 Support: {fp8_supported}")
    print(gpu_info)
    
    #Make temp fasta files for variant scoring
    ref_fasta_path = os.path.join(os.path.dirname(out_file), "temp_ref.fa")
    var_fasta_path = os.path.join(os.path.dirname(out_file), "temp_var.fa")
    
    #Read in variants file
    variant_df = pd.read_csv(variant_file)
    #Split off REF and ALT Alleles into new columns
    variant_df['REF'] = variant_df['MVP ID'].apply(lambda x: x.split(':')[2])
    variant_df['ALT'] = variant_df['MVP ID'].apply(lambda x: x.split(':')[3])
    #Set row indices
    variant_df.index = variant_df['MVP ID']
            
    # Check FASTA
    if not Path(str(fasta_path) + '.fai').exists():
        sys.exit(f"FASTA not indexed: {fasta_path}")
    
    # Initialize model and FASTA
    fasta = Fasta(str(fasta_path))
    
    print(f"Processing {len(variant_df)} variants...")
    num_variants = len(variant_df)
    scored_chunk_dfs = []  # preallocate list of score dataframed
    
    for chunk_start in range(0, num_variants, chunk_size):
        chunk_end = min(chunk_start + chunk_size, num_variants)
        chunk_rows = variant_df.iloc[chunk_start:chunk_end]
        
        # Extract sequences for all variants in chunk
        with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
            futures = {executor.submit(extract_variant_sequences, row, fasta, window_size): row.name for _, row in chunk_rows.iterrows()}
            # Use a dictionary to store results indexed by variant index to preserve order
            results_by_index = {}
            for future in as_completed(futures):
                idx, (ref_seq, alt_seq) = future.result()
                results_by_index[idx] = (ref_seq, alt_seq)
        
        #Write the sequences to fasta files
        ref_entries = []
        var_entries = []
        ref_names = []
        var_names = []
        for variant in results_by_index.keys():
            ref_names.append(f"{variant}_ref")
            var_names.append(f"{variant}_alt")
            ref_entries.append(f">{variant}_ref\n{results_by_index[variant][0]}\n")
            var_entries.append(f">{variant}_alt\n{results_by_index[variant][1]}\n")
        # Write unique sequences to FASTA files
        with open(ref_fasta_path, "w") as f:
            f.writelines(ref_entries)
        with open(var_fasta_path, "w") as f:
            f.writelines(var_entries)
        
        #Append fasta names to df
        chunk_df = variant_df.iloc[chunk_start:chunk_end].copy()
        chunk_df["ref_fasta_name"] = ref_names
        chunk_df["var_fasta_name"] = var_names
              
        #Set output directories
        output_dir = Path(os.path.dirname(out_file))
        predict_ref_dir = output_dir / "reference_predictions"
        predict_var_dir = output_dir / "variant_predictions"
        predict_ref_dir.mkdir(parents=True, exist_ok=True)
        predict_var_dir.mkdir(parents=True, exist_ok=True)
        # Create prediction commands
        predict_ref_command = (
            f"predict_evo2 --fasta {ref_fasta_path} --ckpt-dir {checkpoint_path} "
            f"--output-dir {predict_ref_dir} --model-size {MODEL_SIZE} --tensor-parallel-size {tp_size} "
            f"--pipeline-model-parallel-size 1 --context-parallel-size {cp_size} --output-log-prob-seqs {fp8_option}"
        )
        predict_var_command = (
            f"predict_evo2 --fasta {var_fasta_path} --ckpt-dir {checkpoint_path} "
            f"--output-dir {predict_var_dir} --model-size {MODEL_SIZE} --tensor-parallel-size {tp_size} "
            f"--pipeline-model-parallel-size 1 --context-parallel-size {cp_size} --output-log-prob-seqs {fp8_option}"
        )
        #Score reference sequences
        print(f"Running command: {predict_ref_command}")
        result = run_subprocess_safely(predict_ref_command)
        assert result["returncode"] == 0, result
        #Score alternate sequences
        print(f"Running command: {predict_var_command}")
        result = run_subprocess_safely(predict_var_command)
        assert result["returncode"] == 0, result

        # Find and load prediction files
        ref_pred_files = glob.glob(os.path.join(predict_ref_dir, "predictions__rank_*.pt"))
        var_pred_files = glob.glob(os.path.join(predict_var_dir, "predictions__rank_*.pt"))
        # Load sequence ID maps (maps sequence ID -> prediction index)
        with open(os.path.join(predict_ref_dir, "seq_idx_map.json"), "r") as f:
            ref_seq_idx_map = json.load(f)
        with open(os.path.join(predict_var_dir, "seq_idx_map.json"), "r") as f:
            var_seq_idx_map = json.load(f)
        # Load predictions
        ref_preds = torch.load(ref_pred_files[0])
        var_preds = torch.load(var_pred_files[0])
        # next, calculate change in likelihoods
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
        #Append chunk df to scored chunk dfs
        scored_chunk_dfs.append(chunk_df)
    
    #Concatenate all dataframes together
    scored_df = pd.concat(scored_chunk_dfs)
    scored_df['class'] = np.where(scored_df['Overall PIP'] > 0.95, 'High PIP', 'Low PIP')
    scored_df.to_csv(out_file, index=False)
    
    #Clean up temporary fasta files
    os.remove(ref_fasta_path)
    os.remove(var_fasta_path)

if __name__ == "__main__":
    main()
