"""
Post-processing step (Polaris side): split the merged Evo2 scores file
produced on Betty (run_evo2_worker.sh -> merged_evo2_scores.csv) back into
per-locus files, so generate_priors.py and everything downstream keeps
reading exactly the same per-locus format it always has.

Run this after transferring merged_evo2_scores.csv back from Betty.

Output per-locus CSV: ${BASE}/evo2_scoring/results/{locus_id}.evo2_scores.csv
  SNP_ID, CHR, POS, REF, ALT, ref_log_probs, var_log_probs, evo2_delta_score
"""

import sys
from pathlib import Path

import pandas as pd

BASE = "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping"
RESULTS_DIR = Path(BASE) / "evo2_scoring" / "results"
MERGED_FILE = RESULTS_DIR / "merged_evo2_scores.csv"

OUT_COLS = ["SNP_ID", "CHR", "POS", "REF", "ALT", "ref_log_probs", "var_log_probs", "evo2_delta_score"]

if not MERGED_FILE.exists():
    sys.exit(f"Merged scores file not found: {MERGED_FILE}")

merged = pd.read_csv(MERGED_FILE, dtype={"CHR": str})

written = 0
for locus_id, group in merged.groupby("locus_id"):
    out_path = RESULTS_DIR / f"{locus_id}.evo2_scores.csv"
    group[OUT_COLS].to_csv(out_path, index=False)
    written += 1

print(f"Split {len(merged):,} scored variants into {written} per-locus files in {RESULTS_DIR}", file=sys.stderr)
