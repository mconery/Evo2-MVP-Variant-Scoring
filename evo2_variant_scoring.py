#!/usr/bin/env python3
"""
Evo2 40 B Variant-Scoring Pipeline (Single Process Version)
Author: Arc Institute / NVIDIA Blackwell reference implementation
--------------------------------------------------------------------------------
Single-threaded version for debugging CUDA memory issues
"""

import argparse, logging, os, sys, time
from pathlib import Path
import torch
import numpy as np
import pandas as pd
from pyfaidx import Fasta

# Set environment variables for better CUDA error reporting
os.environ['CUDA_LAUNCH_BLOCKING'] = '1'
os.environ['TORCH_USE_CUDA_DSA'] = '1'

def init_model(model_name: str = "evo2_40b"):
    """Load Evo2 model on single GPU"""
    from evo2 import Evo2
    
    if torch.cuda.device_count() == 0:
        sys.exit("✗ No CUDA devices found")
    
    logging.info(f"Loading {model_name} on GPU 0...")
    torch.cuda.empty_cache()
    
    model = Evo2(model_name)
    
    # Log memory usage
    if torch.cuda.is_available():
        allocated = torch.cuda.memory_allocated() / 1e9
        reserved = torch.cuda.memory_reserved() / 1e9
        logging.info(f"GPU memory: {allocated:.1f}GB allocated, {reserved:.1f}GB reserved")
    
    return model

CTX = 1_000_000

def extract_context(fasta: Fasta, chrom: str, pos: int, ctx: int = CTX):
    half = ctx // 2
    start = max(1, pos - half)
    end = pos + half
    seq = fasta[chrom][start-1:end].seq.upper()
    return seq, start

def make_alt_seq(ref_seq: str, pos: int, ref: str, alt: str, slice_start: int):
    rel = pos - slice_start
    if ref_seq[rel: rel + len(ref)] != ref:
        logging.warning(f"Reference mismatch at {pos} ({ref}→{alt})")
    return ref_seq[:rel] + alt + ref_seq[rel + len(ref):]

@torch.no_grad()
def ll_score(model, tokenizer, seq: str):
    """Score a sequence with extensive error handling"""
    device = "cuda:0"
    
    # Limit sequence length
    MAX_LEN = 500_000  # Start with shorter sequences
    if len(seq) > MAX_LEN:
        logging.warning(f"Truncating sequence from {len(seq)} to {MAX_LEN}")
        seq = seq[:MAX_LEN]
    
    try:
        # Tokenize
        tokens = tokenizer(seq)
        if len(tokens) > MAX_LEN:
            tokens = tokens[:MAX_LEN]
            
        toks = torch.tensor(tokens, dtype=torch.long, device=device)[None, :]
        
        logging.debug(f"Processing sequence: {len(seq)} bp, {toks.size(1)} tokens")
        
        # Forward pass
        out, _ = model(toks)
        logp = torch.log_softmax(out, -1)
        score = logp.mean().item()
        
        logging.debug(f"Score computed: {score}")
        return score
        
    except Exception as e:
        logging.error(f"Error in ll_score: {e}")
        torch.cuda.empty_cache()
        return float('nan')

def score_variant(model, tokenizer, fasta: Fasta, row, ctx: int):
    """Score a single variant"""
    try:
        chrom = row["chr"]
        pos = row["pos"]
        ref = row["ref"]
        alt = row["alt"]
        
        logging.debug(f"Scoring {chrom}:{pos} {ref}>{alt}")
        
        # Extract sequences
        ref_seq, slice_start = extract_context(fasta, chrom, pos, ctx)
        alt_seq = make_alt_seq(ref_seq, pos, ref, alt, slice_start)
        
        # Score both sequences
        ref_score = ll_score(model, tokenizer, ref_seq)
        if np.isnan(ref_score):
            return float('nan')
            
        alt_score = ll_score(model, tokenizer, alt_seq)
        if np.isnan(alt_score):
            return float('nan')
        
        delta = alt_score - ref_score
        logging.debug(f"Variant {chrom}:{pos}: ref={ref_score:.6f}, alt={alt_score:.6f}, delta={delta:.6f}")
        
        return delta
        
    except Exception as e:
        logging.error(f"Error scoring variant {row}: {e}")
        return float('nan')

def main():
    parser = argparse.ArgumentParser("Evo2 40B Single-Process Variant Scorer")
    parser.add_argument("--variants", required=True, help="Variant TSV file")
    parser.add_argument("--fasta", required=True, help="Reference FASTA file")
    parser.add_argument("--out", required=True, help="Output file")
    parser.add_argument("--ctx", type=int, default=500_000, help="Context size (default: 500K)")
    parser.add_argument("--model", default="evo2_40b", help="Model name")
    parser.add_argument("--test", type=int, default=5, help="Test with N variants")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S"
    )
    
    # Load data
    df = pd.read_csv(args.variants, sep="\t")
    required_cols = ["chromosome", "position", "ref_allele", "alt_allele"]
    for col in required_cols:
        if col not in df.columns:
            sys.exit(f"Missing column: {col}")
    
    # Rename columns
    df = df.rename(columns={
        "chromosome": "chr",
        "position": "pos", 
        "ref_allele": "ref",
        "alt_allele": "alt"
    })
    
    # Test mode
    if args.test:
        df = df.head(args.test)
        logging.info(f"Test mode: processing {len(df)} variants")
    
    # Check FASTA
    fasta_path = Path(args.fasta)
    if not Path(str(fasta_path) + '.fai').exists():
        sys.exit(f"FASTA not indexed: {fasta_path}")
    
    # Initialize model and FASTA
    logging.info("Initializing model...")
    model = init_model(args.model)
    tokenizer = model.tokenizer.tokenize
    fasta = Fasta(str(fasta_path))
    
    logging.info(f"Processing {len(df)} variants...")
    
    # Process variants one by one
    scores = []
    for i, (_, row) in enumerate(df.iterrows()):
        logging.info(f"Processing variant {i+1}/{len(df)}: {row['chr']}:{row['pos']}")
        
        score = score_variant(model, tokenizer, fasta, row, args.ctx)
        scores.append(score)
        
        if (i + 1) % 10 == 0:
            valid = sum(1 for s in scores if not np.isnan(s))
            logging.info(f"Completed {i+1}/{len(df)}, valid scores: {valid}")
    
    # Save results
    df["evo2_score"] = scores
    df.to_csv(args.out, sep="\t", index=False)
    
    valid_scores = sum(1 for s in scores if not np.isnan(s))
    logging.info(f"Finished: {valid_scores}/{len(scores)} valid scores saved to {args.out}")

if __name__ == "__main__":
    main()
