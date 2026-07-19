#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=02:00:00
#PBS -l filesystems=grand:home
#PBS -q preemptable
#PBS -A GeomicVar
#PBS -N t2d_setup
#PBS -o /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/logs/setup.log
#PBS -e /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/logs/setup.err

###############################################################################
# One-shot setup job: runs Steps 1-3 and 4a in sequence.
#   Step 1: define T2D EUR loci
#   Step 2: extract per-locus variant lists from GWAS sumstats
#   Step 3: merge per-locus variant lists into one file for transfer to Betty
#   Step 4a: generate LD matrix plink2 commands
#
# Run this first, before submitting Evo2 scoring (on Betty) or LD matrix
# (on Polaris) jobs.
###############################################################################

BASE=/grand/GeomicVar/mconery/evo2_variant_scoring_mapping
REPO=/lus/grand/projects/GeomicVar/mconery/Evo2-MVP-Variant-Scoring
CONDA_BASE=/grand/GeomicVar/mconery/tools/miniconda3

mkdir -p ${BASE}/logs

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate finemap

echo "=== Step 1: Define T2D EUR loci ===" >&2
python ${REPO}/Evo2_Priors_FM/loci_definition/define_t2d_loci.py
echo "Done: $(wc -l < ${BASE}/loci_definition/t2d_eur_loci.tsv) rows in loci file" >&2

echo "=== Step 2: Extract per-locus variant lists ===" >&2
python ${REPO}/Evo2_Priors_FM/loci_definition/extract_locus_variants.py
echo "Done: $(ls ${BASE}/variant_lists/*.tsv 2>/dev/null | wc -l) variant files written" >&2

echo "=== Step 3: Merge per-locus variant lists for Betty transfer ===" >&2
python ${REPO}/Evo2_Priors_FM/loci_definition/merge_variant_lists.py
echo "Done: $(wc -l < ${BASE}/variant_lists/merged_variants.tsv) rows in merged_variants.tsv" >&2

echo "=== Step 4a: Generate LD matrix commands ===" >&2
python ${REPO}/Evo2_Priors_FM/ld_matrices/generate_ld_commands.py
echo "Done: $(wc -l < ${BASE}/ld_matrices/scripts/ld_commands.txt) commands written" >&2

echo "=== Setup complete ===" >&2
echo "Next steps:" >&2
echo "  qsub ${REPO}/Evo2_Priors_FM/ld_matrices/run_ld_matrices.sh   (on Polaris)" >&2
echo "  Transfer to Betty: ${BASE}/variant_lists/merged_variants.tsv" >&2
echo "                     ${BASE}/reference/GRCh37.p13.genome.fa(.fai)" >&2
echo "  On Betty:          bash Evo2_Priors_FM/evo2_scoring/launch_evo2_jobs.sh" >&2
echo "  After scoring, transfer merged_evo2_scores.csv back and run:" >&2
echo "                     python ${REPO}/Evo2_Priors_FM/evo2_scoring/split_scores_by_locus.py" >&2
