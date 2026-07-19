"""
Step 5: Convert per-locus Evo2 delta scores into per-variant Bayesian prior
weights for CARMA.

Prior weighting model
---------------------
Evo2 delta score convention: delta = log_prob(variant) - log_prob(reference).
A more negative delta means the variant more strongly disrupts the DNA
"language model" — i.e., the variant is more likely to be functionally
consequential and therefore has higher prior probability of being causal.

For each variant i in a locus:

    w_i = exp(alpha * (-delta_i))

where alpha controls prior strength.  With alpha = 0.1 the maximum weight
ratio across a typical delta range of ~10 units is exp(0.1 * 10) ≈ 2.7 —
a mild, interpretable signal that will not overwhelm the likelihood.

Weights are then normalised to sum to the number of variants in the locus
so that CARMA's internal expectation of causal variants is preserved on
average.

Output per-locus TSV (tab-separated):
  SNP_ID        - variant identifier (matches CARMA sumstats VARIANT_ID)
  prior_weight  - normalized weight (mean ≈ 1.0)
  delta_score   - raw Evo2 delta score
"""

import math
import sys
from pathlib import Path

import pandas as pd

BASE = "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping"
LOCI_FILE = Path(BASE) / "loci_definition" / "t2d_eur_loci.tsv"
SCORES_DIR = Path(BASE) / "evo2_scoring" / "results"
OUT_DIR = Path(BASE) / "priors"

ALPHA = 0.1   # prior strength; increase to sharpen, decrease to flatten

OUT_DIR.mkdir(parents=True, exist_ok=True)

loci = pd.read_csv(LOCI_FILE, sep="\t", dtype={"chr": str})
written = 0
missing = 0

for _, row in loci.iterrows():
    locus_id = row["locus_id"]
    scores_file = SCORES_DIR / f"{locus_id}.evo2_scores.csv"

    if not scores_file.exists():
        print(f"WARNING: scores not found for {locus_id} — skipping", file=sys.stderr)
        missing += 1
        continue

    df = pd.read_csv(scores_file)
    if df.empty or "evo2_delta_score" not in df.columns:
        print(f"WARNING: empty or malformed scores for {locus_id} — skipping", file=sys.stderr)
        missing += 1
        continue

    df = df.dropna(subset=["evo2_delta_score"]).copy()
    if df.empty:
        print(f"WARNING: all delta scores are NA for {locus_id}", file=sys.stderr)
        missing += 1
        continue

    # w_i = exp(alpha * (-delta_i)); more negative delta -> higher weight
    df["raw_weight"] = df["evo2_delta_score"].apply(lambda d: math.exp(ALPHA * (-d)))

    # Normalize so weights sum to n_variants (preserves CARMA's causal count expectation)
    n = len(df)
    total = df["raw_weight"].sum()
    df["prior_weight"] = df["raw_weight"] * (n / total)

    out_path = OUT_DIR / f"{locus_id}.prior_weights.tsv"
    df[["SNP_ID", "prior_weight", "evo2_delta_score"]].to_csv(
        out_path, sep="\t", index=False
    )
    written += 1

print(
    f"Wrote prior weights for {written} loci  ({missing} missing scores)",
    file=sys.stderr,
)
