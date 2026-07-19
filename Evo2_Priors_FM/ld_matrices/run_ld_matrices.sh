#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=04:00:00
#PBS -l filesystems=grand:home
#PBS -q preemptable
#PBS -A GeomicVar
#PBS -N ld_matrices_t2d
#PBS -o /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/ld_matrices/logs/ld_matrices.log
#PBS -e /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/ld_matrices/logs/ld_matrices.err

###############################################################################
# Step 4b: Run all plink2 LD matrix commands in parallel.
#
# Prerequisites:
#   python loci_definition/define_t2d_loci.py
#   python ld_matrices/generate_ld_commands.py
###############################################################################

BASE=/grand/GeomicVar/mconery/evo2_variant_scoring_mapping
COMMANDS_FILE=${BASE}/ld_matrices/scripts/ld_commands.txt

mkdir -p ${BASE}/ld_matrices/logs

if [ ! -f "${COMMANDS_FILE}" ]; then
    echo "ERROR: Commands file not found: ${COMMANDS_FILE}" >&2
    echo "Run: python ld_matrices/generate_ld_commands.py" >&2
    exit 1
fi

N_COMMANDS=$(wc -l < "${COMMANDS_FILE}")
echo "Running ${N_COMMANDS} plink2 commands with GNU parallel ..." >&2

# Polaris CPU nodes have 32 physical cores; plink2 uses 8 threads each,
# so run 4 jobs in parallel to fill the node.
PARALLEL_JOBS=4
if command -v parallel &>/dev/null; then
    parallel --jobs ${PARALLEL_JOBS} --eta < "${COMMANDS_FILE}"
else
    # Fallback: xargs
    cat "${COMMANDS_FILE}" | xargs -P ${PARALLEL_JOBS} -I{} bash -c '{}'
fi

echo "LD matrix generation complete." >&2

# Report any missing outputs
RESULTS_DIR=${BASE}/ld_matrices/results
LOCI_FILE=${BASE}/loci_definition/t2d_eur_loci.tsv
missing=0
while IFS=$'\t' read -r locus_id _rest; do
    [ "${locus_id}" = "locus_id" ] && continue
    if [ ! -f "${RESULTS_DIR}/${locus_id}.vcor1" ]; then
        echo "MISSING: ${locus_id}" >&2
        missing=$((missing + 1))
    fi
done < "${LOCI_FILE}"

echo "Missing matrices: ${missing}" >&2
