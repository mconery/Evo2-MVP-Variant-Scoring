#!/usr/bin/env python3
"""
Parallel VCF Chunk Processor

Reads VCF files one at a time, parallelizes parsing over chunks of lines in each file,
groups variants by chromosome, sorts them, and writes out chunks of up to `m` variants each.

Outputs gzipped files named:
chr{chromosome}.{start_pos}.{end_pos}.topmed_freeze8_variants.txt.gz

Usage:
    python parallel_vcf_chunk_processor.py /path/to/vcf_dir /output/chunks_dir [-m 2000000] [-c 16]

Requires Python 3.6+.
"""

import os
import glob
import gzip
import argparse
import multiprocessing
from collections import defaultdict
from itertools import islice

def parse_variant_lines(lines):
    """
    Parse chunk of VCF lines.
    Returns list of tuples: (chrom, pos:int, ref, alt)
    Emits one tuple per alt allele if multiple.
    """
    variants = []
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        cols = line.strip().split('\t')
        if len(cols) < 5:
            continue
        chrom = cols[0]
        try:
            pos = int(cols[1])
        except ValueError:
            continue  # skip lines with bad pos
        ref = cols[3]
        alt_alleles = cols[4].split(',')
        for alt in alt_alleles:
            variants.append( (chrom, pos, ref, alt) )
    return variants

def chunked_iterable(iterable, size):
    """Yield fixed-size chunks from an iterable."""
    iterator = iter(iterable)
    while True:
        chunk = list(islice(iterator, size))
        if not chunk:
            break
        yield chunk

def write_chunks(chrom, variant_list, max_variants_per_file, output_dir):
    """
    Given sorted variants for a chromosome,
    write out chunks (gzipped) of up to max_variants_per_file variants.
    File names: chr{chrom}.{start_pos}.{end_pos}.topmed_freeze8_variants.txt.gz
    """
    total = len(variant_list)
    for start_idx in range(0, total, max_variants_per_file):
        chunk_vars = variant_list[start_idx:start_idx+max_variants_per_file]
        start_pos = chunk_vars[0][1]
        end_pos = chunk_vars[-1][1]
        filename = f"{chrom}.{start_pos}.{end_pos}.topmed_freeze8_variants.txt.gz"
        filepath = os.path.join(output_dir, filename)
        with gzip.open(filepath, 'wt') as out:
            out.write("chromosome\tposition\tref_allele\talt_allele\n")
            for v in chunk_vars:
                out.write(f"{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\n")

def process_single_vcf(vcf_file, max_variants_per_file, output_dir, n_cores):
    """Process one VCF file, parsing in parallel chunks and writing output chunk files."""

    print(f"Processing file: {vcf_file}")

    # Read all lines (to enable chunking)
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        lines = f.readlines()
    
    # Filter out header lines first (they start with '#')
    body_lines = [l for l in lines if not l.startswith('#') and l.strip()]
    print(f"Total variant lines (excluding headers): {len(body_lines):,}")

    if not body_lines:
        print("No variant lines found, skipping.")
        return

    # Chunk lines into pieces for parallel parsing (adjust size for chunking speed)
    chunk_size = max(100000, len(body_lines) // n_cores)  # at least 100k lines per chunk or roughly equal split
    chunks = list(chunked_iterable(body_lines, chunk_size))

    print(f"Parsing VCF body lines in {len(chunks)} chunks with {n_cores} cores...")

    # Parallel parse variant lines
    with multiprocessing.Pool(processes=n_cores) as pool:
        results = pool.map(parse_variant_lines, chunks)

    # Flatten all variants (list of lists -> list)
    all_variants = [v for sublist in results for v in sublist]
    print(f"Total variants parsed: {len(all_variants):,}")

    # Group by chromosome
    chrom_variants = defaultdict(list)
    for variant in all_variants:
        chrom_variants[variant[0]].append(variant)

    print(f"Variants grouped into {len(chrom_variants)} chromosomes.")

    # Sort variants in each chromosome by position
    for chrom in chrom_variants:
        chrom_variants[chrom].sort(key=lambda x: x[1])

    # Write chunked gzipped files chromosome by chromosome
    os.makedirs(output_dir, exist_ok=True)
    print(f"Writing chunked output files to '{output_dir}' (max {max_variants_per_file:,} variants per file)...")
    for chrom in sorted(chrom_variants.keys()):
        write_chunks(chrom, chrom_variants[chrom], max_variants_per_file, output_dir)

    print(f"Finished processing {vcf_file}.\n")

def main():
    parser = argparse.ArgumentParser(description="Parallel VCF Chunk Processor")
    parser.add_argument("vcf_dir", help="Directory with input VCF (.vcf or .vcf.gz) files")
    parser.add_argument("output_dir", help="Directory to save output chunked variant files")
    parser.add_argument("-m", "--max_variants", type=int, default=2_000_000,
                        help="Maximum number of variants per output chunk (default: 2,000,000)")
    parser.add_argument("-c", "--cores", type=int, default=None,
                        help="Number of CPU cores to use (default: use all available cores)")

    args = parser.parse_args()

    vcf_files = glob.glob(os.path.join(args.vcf_dir, "*.vcf")) + \
                glob.glob(os.path.join(args.vcf_dir, "*.vcf.gz"))

    if not vcf_files:
        print(f"No VCF files found in directory: {args.vcf_dir}")
        return

    n_cores = args.cores if args.cores else multiprocessing.cpu_count()
    print(f"Using {n_cores} CPU cores for parallel processing.\n")

    for vcf_file in sorted(vcf_files):
        process_single_vcf(vcf_file, args.max_variants, args.output_dir, n_cores)

    print("All files processed.")

if __name__ == "__main__":
    main()
