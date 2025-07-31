#!/usr/bin/env python3
"""
Simple VCF Variant Extractor
Extracts chromosome, position, and allele information from VCF files.
"""
import os
import glob
import gzip
import multiprocessing
from tempfile import TemporaryDirectory

def extract_variants_from_vcf(input_vcf, out_path):
    """Extract variants from a single VCF file and write results to out_path."""
    opener = gzip.open if input_vcf.endswith('.gz') else open
    with opener(input_vcf, 'rt') as f, open(out_path, 'w') as out:
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

def main(vcf_directory, output_file="variants.tsv.gz"):
    vcf_files = glob.glob(os.path.join(vcf_directory, "*.vcf")) + \
                glob.glob(os.path.join(vcf_directory, "*.vcf.gz"))

    if not vcf_files:
        print(f"No VCF files found in {vcf_directory}")
        return

    n_cores = multiprocessing.cpu_count()
    print(f"Detected {n_cores} CPU cores for parallel processing.")

    with TemporaryDirectory() as tmpdir:
        # Prepare args: (input_vcf, out_path for variant results)
        job_args = []
        for i, vcf_file in enumerate(vcf_files):
            tmp_out = os.path.join(tmpdir, f"chunk_{i}.tsv")
            job_args.append((vcf_file, tmp_out))

        # Define wrapper for starmap
        def worker(args):
            extract_variants_from_vcf(*args)

        # Process in parallel
        with multiprocessing.Pool(n_cores) as pool:
            pool.map(worker, job_args)

        # Write header and concatenate all chunks to gzipped output
        with gzip.open(output_file, 'wt') as gzout:
            gzout.write("chromosome\tposition\tref_allele\talt_allele\n")
            for _, chunk_path in job_args:
                with open(chunk_path, 'rt') as chunk_file:
                    for line in chunk_file:
                        gzout.write(line)

    print(f"Parallel extraction complete. Output saved to: {output_file}")

if __name__ == "__main__":
    import sys
    vcf_dir = sys.argv[1] if len(sys.argv) > 1 else input("Enter path to VCF directory: ").strip()
    output = sys.argv[2] if len(sys.argv) > 2 else "variants.tsv.gz"
    main(vcf_dir, output)
