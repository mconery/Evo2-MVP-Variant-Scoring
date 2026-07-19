#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=01:00:00
#PBS -l filesystems=grand:home
#PBS -q preemptable
#PBS -A GeomicVar
#PBS -N t2d_collation
#PBS -o /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/collation/logs/collation.log
#PBS -e /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/collation/logs/collation.err

###############################################################################
# Step 7: Collate CARMA fine-mapping results and compare with S11 mapping.
#
# Uses RSAIGE_GPU_V2 conda env, which has dplyr, readr, ggplot2, purrr,
# stringr (all required by compare_finemapping.R).
#
# Prerequisites:
#   - CARMA runs complete (carma_results/without_priors/*.carma.tsv present)
#   - Evo2-prior CARMA runs complete (carma_results/with_priors/*.carma.tsv)
###############################################################################

BASE=/grand/GeomicVar/mconery/evo2_variant_scoring_mapping
REPO=/lus/grand/projects/GeomicVar/mconery/Evo2-MVP-Variant-Scoring
CONDA_BASE=/grand/GeomicVar/mconery/tools/miniconda3

mkdir -p ${BASE}/collation/logs ${BASE}/collation/plots

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate RSAIGE_GPU_V2

echo "=== Step 7: Collation and comparison ===" >&2
echo "Without-prior results: $(ls ${BASE}/carma_results/without_priors/*.carma.tsv 2>/dev/null | wc -l)" >&2
echo "With-prior results:    $(ls ${BASE}/carma_results/with_priors/*.carma.tsv 2>/dev/null | wc -l)" >&2

Rscript ${REPO}/Evo2_Priors_FM/collation/compare_finemapping.R

echo "=== Collation complete ===" >&2
echo "Outputs in ${BASE}/collation/" >&2
