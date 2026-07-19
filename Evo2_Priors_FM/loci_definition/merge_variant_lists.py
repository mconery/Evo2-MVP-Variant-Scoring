"""
Step 3: Merge all per-locus variant TSVs into a single file for transfer to
the Evo2 scoring cluster (PARCC Betty).

Reads every ${BASE}/variant_lists/{locus_id}.variants.tsv written by
extract_locus_variants.py and concatenates them into one file, adding:
  locus_id  - parsed from the filename
  uid       - f"{locus_id}__{SNP_ID}", a collision-safe key used by the Evo2
              scoring script (two loci could in principle share an rsID)

Output (tab-separated): ${BASE}/variant_lists/merged_variants.tsv
  uid  locus_id  SNP_ID  CHR  POS  REF  ALT  BETA  SE  N  PVAL
"""

import sys
from pathlib import Path

import pandas as pd

BASE = "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping"
VARIANT_DIR = Path(BASE) / "variant_lists"
OUT_FILE = VARIANT_DIR / "merged_variants.tsv"

variant_files = sorted(VARIANT_DIR.glob("*.variants.tsv"))
if not variant_files:
    sys.exit(f"No per-locus variant files found in {VARIANT_DIR}")

frames = []
for f in variant_files:
    locus_id = f.name[: -len(".variants.tsv")]
    df = pd.read_csv(f, sep="\t", dtype={"CHR": str})
    df.insert(0, "locus_id", locus_id)
    df.insert(0, "uid", locus_id + "__" + df["SNP_ID"].astype(str))
    frames.append(df)

merged = pd.concat(frames, ignore_index=True)
merged.to_csv(OUT_FILE, sep="\t", index=False)

print(
    f"Merged {len(variant_files)} locus files into {OUT_FILE} "
    f"({len(merged):,} variants total)",
    file=sys.stderr,
)
