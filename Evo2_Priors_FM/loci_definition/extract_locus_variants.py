"""
Step 2: For each T2D EUR locus, extract all variants from the MVP EUR GWAS
summary statistics (hg19 / GRCh37) and write a per-locus TSV.

The hg19 sumstats already contain BETA and SE; no OR/CI conversion needed.
RSID is used as SNP_ID to match the 1000G pvar variant identifiers.
Variants without an rsID (RSID == ".") are excluded.

Output per-locus TSV columns (tab-separated):
  SNP_ID  - rsID (matches 1000G pvar IDs)
  CHR     - chromosome (numeric, no chr prefix)
  POS     - position (hg19)
  REF     - reference allele
  ALT     - alternative allele
  BETA    - effect size (log-OR scale)
  SE      - standard error of BETA
  N       - total sample size
  PVAL    - association p-value
"""

import sys
from pathlib import Path

import pandas as pd

SUMSTATS = (
    "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping/summary_stats/"
    "MVP_R4.1000G_AGR.Phe_250_2.EUR.GIA.dbGaP.sumstats.hg19.txt.gz"
)
BASE = "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping"
LOCI_FILE = Path(BASE) / "loci_definition" / "t2d_eur_loci.tsv"
OUT_DIR = Path(BASE) / "variant_lists"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------
# Load loci
# ------------------------------------------------------------------
loci = pd.read_csv(LOCI_FILE, sep="\t", dtype={"chr": str})
loci["start"] = loci["start"].astype(int)
loci["end"] = loci["end"].astype(int)
print(f"Loaded {len(loci)} loci", file=sys.stderr)

# ------------------------------------------------------------------
# Load sumstats in one shot
# The file has a comment-style header (#ID RSID CHR POS ...)
# pandas handles the leading '#' in the column name automatically.
# ------------------------------------------------------------------
print(f"Reading {SUMSTATS} ...", file=sys.stderr)
ss = pd.read_csv(
    SUMSTATS, sep="\t", compression="gzip",
    dtype={"CHR": str},
    usecols=["RSID", "CHR", "POS", "REF", "ALT", "BETA", "SE", "N", "P"],
    comment=None,   # '#' in header is fine; pandas reads it as part of col name
)

# The first column may be read as "#ID" depending on pandas version;
# rename defensively.
ss.columns = [c.lstrip("#") for c in ss.columns]

ss["POS"] = ss["POS"].astype(int)
ss["BETA"] = pd.to_numeric(ss["BETA"], errors="coerce")
ss["SE"]   = pd.to_numeric(ss["SE"],   errors="coerce")
ss["N"]    = pd.to_numeric(ss["N"],    errors="coerce")
ss["P"]    = pd.to_numeric(ss["P"],    errors="coerce")

# Drop variants without a valid rsID or without BETA/SE
ss = ss[ss["RSID"].notna() & (ss["RSID"] != ".")]
ss = ss.dropna(subset=["BETA", "SE"]).copy()
print(f"Loaded {len(ss):,} variants with valid rsID and BETA/SE", file=sys.stderr)

# ------------------------------------------------------------------
# Assign variants to loci using vectorised range checks
# ------------------------------------------------------------------
written = 0
empty = 0

for _, lrow in loci.iterrows():
    locus_id = lrow["locus_id"]
    chrom    = str(lrow["chr"])   # numeric, no "chr" prefix — matches ss["CHR"]
    start    = lrow["start"]
    end      = lrow["end"]

    mask = (ss["CHR"] == chrom) & (ss["POS"] >= start) & (ss["POS"] <= end)
    sub = ss.loc[mask, ["RSID", "CHR", "POS", "REF", "ALT", "BETA", "SE", "N", "P"]].copy()
    sub = sub.rename(columns={"RSID": "SNP_ID", "P": "PVAL"})

    if sub.empty:
        print(f"WARNING: no variants for locus {locus_id}", file=sys.stderr)
        empty += 1
        continue

    out_path = OUT_DIR / f"{locus_id}.variants.tsv"
    sub.to_csv(out_path, sep="\t", index=False)
    written += 1

print(
    f"Wrote {written} locus variant files to {OUT_DIR}  ({empty} empty loci)",
    file=sys.stderr,
)
