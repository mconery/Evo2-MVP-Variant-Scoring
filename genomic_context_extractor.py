#!/usr/bin/env python3
"""
Genomic Context Extractor for Variants

Given a variant chunk file (tsv or tsv.gz) with columns:
chromosome    position    ref_allele    alt_allele

and a GRCh38 reference FASTA file (indexed with pyfaidx),

this script extracts two fasta sequences per variant:
 - Reference allele context sequence
 - Alternate allele context sequence

Each sequence is 1 million nucleotides (default; configurable) centered on the variant.
For indels, the window is adjusted to keep the variant centered.
Boundaries near chromosome ends are handled gracefully.

Outputs fasta files to an output directory (default: "variant_contexts"),
named as:
 chr{chromosome}_{position}_{ref}_{alt}_{ref|alt}.fa

Usage:
    python genomic_context_extractor.py variants_chunk.tsv.gz /path/to/GRCh38.fa \
        -w 500000 \
        -o variant_contexts
"""

import os
import sys
import argparse
import gzip
from pyfaidx import Fasta
import re

def parse_variants(file_path):
    """
    Parses the variant file (tsv or tsv.gz).
    Yields tuples: (chrom, int(position), ref_allele, alt_allele)
    """
    opener = gzip.open if file_path.endswith(".gz") else open
    with opener(file_path, 'rt') as f:
        header = f.readline()
        header_cols = header.strip().split('\t')
        expected_cols = ["chromosome", "position", "ref_allele", "alt_allele"]
        if header_cols[:4] != expected_cols:
            print(f"Warning: Unexpected header columns in {file_path}: {header_cols}", file=sys.stderr)

        for line in f:
            if not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 4:
                continue
            chrom, pos_str, ref, alt = cols[:4]
            try:
                pos = int(pos_str)
            except ValueError:
                print(f"Skipping line with invalid position: {line.strip()}", file=sys.stderr)
                continue
            yield chrom, pos, ref, alt

def sanitize_chrom_name(chrom):
    """
    Normalize chromosome name to match reference fasta.
    Accepts chrom names with or without 'chr', returns version matching fasta.
    """
    # We will try several variants if needed during lookup.
    return chrom

def adjust_window(pos, ref, alt, window_radius, chrom_len):
    """
    Calculate start and end coordinates of context window.

    - pos is 1-based variant start position.
    - ref and alt are allele strings.
    - window_radius is half-window size (i.e., 500,000 for 1M total).
    - For indels, window is extended to keep variant centered by adjusting start or end.
    - Coordinates are constrained to [1, chrom_len].

    Returns (start, end) inclusive.
    """

    # Variant length difference (alt - ref)
    indel_len = len(alt) - len(ref)

    # Window boundaries initially (1-based coordinates)
    start = pos - window_radius
    end = pos + window_radius - 1  # inclusive window; length will be end - start + 1

    # Adjust window so the variant is centered, considering indels
    # If alt is longer => extend window end, shift window start left
    # If alt shorter => extend window start, shift window end right
    # Goal: keep variant evenly centered in window, adjusting for allele length difference

    ## Calculate the padding shift needed:
    # For indels, shift half indel_len bases from one side to the other
    # We'll shift start and end accordingly to keep center near variant middle

    if indel_len > 0:
        # insertion or longer alt allele: shift end by indel_len
        end += indel_len // 2 + indel_len % 2
        start -= indel_len // 2
    elif indel_len < 0:
        # deletion or shorter alt allele: shift start lower by half length 
        shift = (-indel_len) // 2 + (-indel_len) % 2
        start -= shift
        end += shift

    # Clamp to chromosome boundaries
    if start < 1:
        shift_up = 1 - start
        start = 1
        end = min(chrom_len, end + shift_up)
    if end > chrom_len:
        shift_down = end - chrom_len
        end = chrom_len
        start = max(1, start - shift_down)

    # Final length check (optional)
    # length = end - start + 1
    # If length != 2 * window_radius, this is due to chromosome boundary.
    return start, end

def build_alt_sequence(ref_seq, pos, ref, alt, window_start):
    """
    Build the alternate allele sequence by applying the variant substitution
    into the reference sequence string.

    Parameters:
    - ref_seq: str reference sequence from window_start to window_end (1-based)
    - pos: 1-based variant position (start)
    - ref: reference allele string
    - alt: alternate allele string
    - window_start: genomic position of first base in ref_seq (1-based)
    
    Returns:
    - alt_seq: str with alt allele sequence incorporated
    """

    # Variant offset in ref_seq (0-based)
    offset = pos - window_start

    # Confirm that ref sequence matches expected ref allele
    ref_in_seq = ref_seq[offset:offset + len(ref)].upper()
    if ref_in_seq != ref.upper():
        print(f"Warning: Reference allele mismatch at pos {pos}. Expected {ref}, found {ref_in_seq} in fasta.", file=sys.stderr)
        # We proceed anyway, substituting to alt

    # Construct alt sequence
    alt_seq = ref_seq[:offset] + alt + ref_seq[offset + len(ref):]

    return alt_seq

def sanitize_filename(s):
    """Make string safe for filenames by removing/replacing problematic characters."""
    return re.sub(r'[^A-Za-z0-9_.-]', '_', s)

def write_fasta(file_path, header, sequence):
    """
    Write a FASTA file with the given header and sequence.
    Wrap lines at 60 characters.
    """
    with open(file_path, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + '\n')

def main():
    parser = argparse.ArgumentParser(description="Extract genomic sequence windows around variants with ref and alt alleles")
    parser.add_argument("variant_file", help="Input variant chunk file (.tsv or .tsv.gz)")
    parser.add_argument("reference_fasta", help="FASTA file of GRCh38 reference genome (indexed with pyfaidx)")
    parser.add_argument("-w", "--half_window", type=int, default=500000,
                        help="Half window size in bp (default 500,000 for 1M total context)")
    parser.add_argument("-o", "--output_dir", default="variant_contexts",
                        help="Output directory to save fasta files")
    parser.add_argument("--skip_ambiguous", action='store_true',
                        help="Skip variants with ambiguous REF or ALT (e.g. containing 'N' or other non-ACGT bases)")
    args = parser.parse_args()

    fasta = Fasta(args.reference_fasta, sequence_always_upper=True)
    os.makedirs(args.output_dir, exist_ok=True)

    variant_count = 0
    skipped_count = 0

    for chrom, pos, ref, alt in parse_variants(args.variant_file):
        # Normalize chromosome naming:
        chrom_name = chrom
        if chrom_name not in fasta:
            # try removing or adding 'chr' prefix
            if chrom_name.startswith("chr") and chrom_name[3:] in fasta:
                chrom_name = chrom_name[3:]
            elif ("chr" + chrom_name) in fasta:
                chrom_name = "chr" + chrom_name
            else:
                print(f"ERROR: Chromosome '{chrom}' not found in reference fasta. Skipping variant at {chrom}:{pos}", file=sys.stderr)
                skipped_count += 1
                continue

        chrom_len = len(fasta[chrom_name])

        # Skip variants that contain ambiguous bases if requested
        if args.skip_ambiguous:
            if (not re.match(r'^[ACGTNacgtn]+$', ref)) or (not re.match(r'^[ACGTNacgtn]+$', alt)):
                print(f"Skipping variant with ambiguous allele: {chrom}:{pos} {ref}>{alt}", file=sys.stderr)
                skipped_count += 1
                continue

        # Calculate window coordinates for this variant
        window_start, window_end = adjust_window(pos, ref, alt, args.half_window, chrom_len)

        # Extract reference sequence window (pyfaidx is 0-based internally, slice is inclusive:exclusive)
        # FASTA is 1-based inclusive indexing, pyfaidx supports 1-based slice.
        ref_seq = fasta[chrom_name][window_start - 1:window_end].seq

        # Build alt allele sequence
        alt_seq = build_alt_sequence(ref_seq, pos, ref, alt, window_start)

        # Prepare filenames and headers
        safe_chrom = sanitize_filename(chrom_name)
        safe_ref = sanitize_filename(ref)
        safe_alt = sanitize_filename(alt)
        basename = f"{safe_chrom}_{pos}_{safe_ref}_{safe_alt}"

        ref_fasta_path = os.path.join(args.output_dir, f"{basename}_ref.fa")
        alt_fasta_path = os.path.join(args.output_dir, f"{basename}_alt.fa")

        ref_header = f"{chrom_name}:{window_start}-{window_end}_ref_{ref}"
        alt_header = f"{chrom_name}:{window_start}-{window_end}_alt_{alt}"

        # Write fasta files
        write_fasta(ref_fasta_path, ref_header, ref_seq)
        write_fasta(alt_fasta_path, alt_header, alt_seq)

        variant_count += 1
        if variant_count % 1000 == 0:
            print(f"Processed {variant_count} variants...", file=sys.stderr)

    print(f"Completed processing {variant_count} variants.")
    if skipped_count > 0:
        print(f"Skipped {skipped_count} variants due to errors or ambiguous bases.", file=sys.stderr)


if __name__ == '__main__':
    main()
