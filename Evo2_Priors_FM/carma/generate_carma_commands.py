"""
Step 6b: Generate CARMA command pairs (with and without Evo2 priors)
for all T2D EUR loci that have both an LD matrix and a variant sumstats file.

Output:
  {BASE}/carma_results/scripts/carma_commands.txt  -- one command per line
"""

import sys
from pathlib import Path

import pandas as pd

BASE = "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping"
LOCI_FILE = Path(BASE) / "loci_definition" / "t2d_eur_loci.tsv"
LD_DIR = Path(BASE) / "ld_matrices" / "results"
SUMSTATS_DIR = Path(BASE) / "variant_lists"
PRIORS_DIR = Path(BASE) / "priors"
OUT_NO_PRIOR = Path(BASE) / "carma_results" / "without_priors"
OUT_WITH_PRIOR = Path(BASE) / "carma_results" / "with_priors"
SCRIPTS_DIR = Path(BASE) / "carma_results" / "scripts"
COMMANDS_FILE = SCRIPTS_DIR / "carma_commands.txt"

REPO = "/lus/grand/projects/GeomicVar/mconery/Evo2-MVP-Variant-Scoring"
R_SCRIPT = f"{REPO}/Evo2_Priors_FM/carma/finemap_CARMA_t2d.R"

for d in [OUT_NO_PRIOR, OUT_WITH_PRIOR, SCRIPTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

loci = pd.read_csv(LOCI_FILE, sep="\t", dtype={"chr": str})

commands = []
skipped_no_ld = 0
skipped_complete = 0

for _, row in loci.iterrows():
    locus_id = row["locus_id"]

    vcor1 = LD_DIR / f"{locus_id}.vcor"
    sumstats = SUMSTATS_DIR / f"{locus_id}.variants.tsv"
    prior_file = PRIORS_DIR / f"{locus_id}.prior_weights.tsv"

    if not vcor1.exists() or not Path(str(vcor1) + ".vars").exists():
        print(f"SKIP (no LD matrix): {locus_id}", file=sys.stderr)
        skipped_no_ld += 1
        continue

    if not sumstats.exists():
        print(f"SKIP (no sumstats): {locus_id}", file=sys.stderr)
        skipped_no_ld += 1
        continue

    # Without priors
    out_no = OUT_NO_PRIOR / f"{locus_id}.carma.tsv"
    if not out_no.exists():
        commands.append(
            f"Rscript {R_SCRIPT} "
            f"'{sumstats}' '{vcor1}' '{out_no}'"
        )
    else:
        skipped_complete += 1

    # With priors (only if prior file exists)
    out_with = OUT_WITH_PRIOR / f"{locus_id}.carma.tsv"
    if prior_file.exists():
        if not out_with.exists():
            commands.append(
                f"Rscript {R_SCRIPT} "
                f"'{sumstats}' '{vcor1}' '{out_with}' "
                f"--use-priors --prior-file '{prior_file}' --input-alpha 0.1"
            )
        else:
            skipped_complete += 1
    else:
        print(f"NOTE: no prior file for {locus_id}; only uniform run queued", file=sys.stderr)

with open(COMMANDS_FILE, "w") as f:
    f.write("\n".join(commands) + "\n")

print(
    f"Written {len(commands)} CARMA commands to {COMMANDS_FILE}  "
    f"(skipped: {skipped_no_ld} missing LD/sumstats, {skipped_complete} already complete)",
    file=sys.stderr,
)
