"""
Step 4a: Generate plink2 --r-unphased square commands to compute EUR LD matrices
for all 105 T2D loci from the 1000 Genomes Phase 3 plink2 files (hg19/GRCh37).

Uses --r-unphased square to produce a dense square correlation matrix and a
companion .vcor.vars file listing the variant IDs in matrix order.
Only variants present in the per-locus GWAS summary statistics are included
(via --extract), keeping matrix sizes manageable for CARMA.

One command per locus is written to ld_commands.txt, which is then consumed
by run_ld_matrices.sh.

Output:
  {BASE}/ld_matrices/extract/  {locus_id}.extract.txt  -- per-locus RSID lists
  {BASE}/ld_matrices/scripts/  ld_commands.txt          -- one shell command per line
"""

import sys
from pathlib import Path

import pandas as pd

BASE = "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping"
LOCI_FILE = Path(BASE) / "loci_definition" / "t2d_eur_loci.tsv"
VARIANT_DIR = Path(BASE) / "variant_lists"
OUT_DIR = Path(BASE) / "ld_matrices"
SCRIPTS_DIR = OUT_DIR / "scripts"
RESULTS_DIR = OUT_DIR / "results"
EXTRACT_DIR = OUT_DIR / "extract"
COMMANDS_FILE = SCRIPTS_DIR / "ld_commands.txt"

PFILE = "/grand/GeomicVar/mconery/resources/1000G/all_phase3"
KEEP_FILE = "/grand/GeomicVar/mconery/resources/1000G/1000G.EUR.txt"
PLINK2 = "/grand/GeomicVar/mconery/tools/bin/plink2"

THREADS = 8

for d in [SCRIPTS_DIR, RESULTS_DIR, EXTRACT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

loci = pd.read_csv(LOCI_FILE, sep="\t", dtype={"chr": str})
print(f"Generating commands for {len(loci)} loci ...", file=sys.stderr)

commands = []
skipped_no_variants = 0

for _, row in loci.iterrows():
    locus_id = row["locus_id"]
    out_prefix = RESULTS_DIR / locus_id

    # Load per-locus variant list to build extract file
    variant_file = VARIANT_DIR / f"{locus_id}.variants.tsv"
    if not variant_file.exists():
        print(f"WARNING: no variant file for {locus_id}; skipping", file=sys.stderr)
        skipped_no_variants += 1
        continue

    variants = pd.read_csv(variant_file, sep="\t", usecols=["SNP_ID"])
    extract_file = EXTRACT_DIR / f"{locus_id}.extract.txt"
    variants["SNP_ID"].to_csv(extract_file, index=False, header=False)

    # Skip if already complete (dense square format produces .vcor.vars)
    vcor_vars = Path(str(out_prefix) + ".vcor.vars")
    skip_check = (
        f"[ -f '{vcor_vars}' ] && "
        f"echo 'SKIP {locus_id}' || "
    )

    cmd = (
        f"{skip_check}"
        f"{PLINK2} "
        f"--pfile {PFILE} "
        f"--keep {KEEP_FILE} "
        f"--extract '{extract_file}' "
        f"--r-unphased square "
        f"--out '{out_prefix}' "
        f"--threads {THREADS} "
        f"--silent"
    )
    commands.append(cmd)

with open(COMMANDS_FILE, "w") as f:
    f.write("\n".join(commands) + "\n")

print(
    f"Written {len(commands)} commands to {COMMANDS_FILE}  "
    f"({skipped_no_variants} skipped: no variant file)",
    file=sys.stderr,
)
