"""
Step 1: Extract the 105 T2D EUR fine-mapped loci from the MVP phenome-wide
supplementary table (S11).

Criteria:
  - phenotype == "Phe_250_2"
  - EUR.best_variant is non-empty and non-NA (at least one EUR fine-mapped signal)

De-duplicates on (chr, start, end) to yield one row per genomic locus.

Output columns:
  locus_id  - e.g.  chr10.114250001.115250000  (matches S11 locus column)
  chr       - numeric, no "chr" prefix  (e.g. 10)   -- used by plink2 --chr
  chr_full  - with "chr" prefix         (e.g. chr10) -- used in file names / FASTA lookups
  start     - integer locus start (hg19)
  end       - integer locus end   (hg19)
"""

import pandas as pd
import sys
from pathlib import Path

S11 = "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping/mapping_loci/Supplementary_Table-S11.txt"
BASE = "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping"
OUT_DIR = Path(BASE) / "loci_definition"
OUT_FILE = OUT_DIR / "t2d_eur_loci.tsv"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# The file has two header rows; the first row is the real header.
# Read with header=0 then drop any row where phenotype == "phenotype" (duplicate header).
df = pd.read_csv(S11, sep="\t", header=0, dtype=str, low_memory=False)
df = df[df["phenotype"] != "phenotype"].copy()

# Filter: T2D AND EUR signal present
t2d = df[df["phenotype"] == "Phe_250_2"].copy()
eur_mapped = t2d[
    t2d["EUR.best_variant"].notna() & (t2d["EUR.best_variant"] != "")
].copy()

# De-duplicate on genomic locus
eur_mapped["start"] = eur_mapped["start"].astype(int)
eur_mapped["end"] = eur_mapped["end"].astype(int)

unique_loci = (
    eur_mapped[["locus", "chr", "start", "end"]]
    .drop_duplicates(subset=["locus"])
    .copy()
)
unique_loci = unique_loci.rename(columns={"locus": "locus_id", "chr": "chr_full"})

# chr_full is e.g. "chr10"; strip prefix for plink2
unique_loci["chr"] = unique_loci["chr_full"].str.replace("chr", "", regex=False)
unique_loci = unique_loci.sort_values(["chr", "start"]).reset_index(drop=True)

print(f"Found {len(unique_loci)} unique T2D EUR loci", file=sys.stderr)
assert len(unique_loci) > 0, "No loci found — check column indices in S11"

unique_loci[["locus_id", "chr", "chr_full", "start", "end"]].to_csv(
    OUT_FILE, sep="\t", index=False
)
print(f"Written: {OUT_FILE}", file=sys.stderr)
