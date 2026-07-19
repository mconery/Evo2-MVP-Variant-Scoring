#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=08:00:00
#PBS -l filesystems=grand:home
#PBS -q preemptable
#PBS -A GeomicVar
#PBS -N carma_t2d
#PBS -o /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/carma_results/logs/carma.log
#PBS -e /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/carma_results/logs/carma.err

###############################################################################
# Step 6c: Run all CARMA fine-mapping commands (with and without priors) in
# parallel on a single Polaris CPU node.
#
# Prerequisites:
#   python carma/generate_carma_commands.py
###############################################################################

BASE=/grand/GeomicVar/mconery/evo2_variant_scoring_mapping
COMMANDS_FILE=${BASE}/carma_results/scripts/carma_commands.txt
CONDA_BASE=/grand/GeomicVar/mconery/tools/miniconda3

mkdir -p ${BASE}/carma_results/logs

if [ ! -f "${COMMANDS_FILE}" ]; then
    echo "ERROR: Commands file not found: ${COMMANDS_FILE}" >&2
    echo "Run: python carma/generate_carma_commands.py" >&2
    exit 1
fi

# Activate CARMA conda environment
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate carma_env

N_COMMANDS=$(wc -l < "${COMMANDS_FILE}")
echo "Running ${N_COMMANDS} CARMA commands ..." >&2

# Each CARMA run is single-threaded; Polaris CPU node has 32 cores.
# With ~210 commands (105 loci × 2), 32 parallel workers finishes in ~30 min.
PARALLEL_JOBS=32

if command -v parallel &>/dev/null; then
    parallel --jobs ${PARALLEL_JOBS} --eta < "${COMMANDS_FILE}"
else
    cat "${COMMANDS_FILE}" | xargs -P ${PARALLEL_JOBS} -I{} bash -c '{}'
fi

echo "CARMA runs complete." >&2

# Report completion
no_prior_dir=${BASE}/carma_results/without_priors
with_prior_dir=${BASE}/carma_results/with_priors
n_done_no=$(ls "${no_prior_dir}"/*.carma.tsv 2>/dev/null | wc -l)
n_done_with=$(ls "${with_prior_dir}"/*.carma.tsv 2>/dev/null | wc -l)
echo "Completed: ${n_done_no} without-prior runs, ${n_done_with} with-prior runs" >&2
