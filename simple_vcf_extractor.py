#!/usr/bin/env python3
"""
Simple VCF Variant Extractor
Extracts chromosome, position, and allele information from VCF files.
"""

import os
import glob
import gzip

def extract_variants_from_directory(vcf_directory: str, output_file: str = "variants.tsv.gz"):
    vcf_files = glob.glob(os.path.join(vcf_directory, "*.vcf"))
    vcf_gz_files = glob.glob(os.path.join(vcf_directory, "*.vcf.gz"))
    all_files = vcf_files + vcf_gz_files

    if not all_files:
        print(f"No VCF files found in {vcf_directory}")
        return

    print(f"Found {len(all_files)} VCF files")
    variant_count = 0

    with gzip.open(output_file, 'wt') as out:
        out.write("chromosome\tposition\tref_allele\talt_allele\n")
        for vcf_file in all_files:
            print(f"Processing: {os.path.basename(vcf_file)}")
            if vcf_file.endswith('.gz'):
                file_handle = gzip.open(vcf_file, 'rt')
            else:
                file_handle = open(vcf_file, 'r')
            with file_handle as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    columns = line.strip().split('\t')
                    if len(columns) >= 5:
                        chromosome = columns[0]
                        position = columns[1]
                        ref_allele = columns[3]
                        alt_alleles = columns[4].split(',')
                        for alt_allele in alt_alleles:
                            out.write(f"{chromosome}\t{position}\t{ref_allele}\t{alt_allele}\n")
                            variant_count += 1

    print(f"\nExtraction complete!")
    print(f"Total variants extracted: {variant_count:,}")
    print(f"Output saved to: {output_file}")

if __name__ == "__main__":
    import sys
    
    # Get directory path from command line or prompt user
    if len(sys.argv) > 1:
        directory = sys.argv[1]
    else:
        directory = input("Enter path to VCF directory: ").strip()
    
    # Get output file name (optional)
    output = sys.argv[2] if len(sys.argv) > 2 else "variants.tsv"
    
    # Run extraction
    extract_variants_from_directory(directory, output)
