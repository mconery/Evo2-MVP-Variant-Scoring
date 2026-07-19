#!/bin/bash
# Launch Evo2 scoring PBS jobs in batches to stay within the 20-job queue limit.
#
# Usage:
#   bash launch_evo2_jobs.sh [--batch-size N] [--chunk-size N] [--dry-run]
#
# Options:
#   --batch-size N   Loci per PBS job (default: 7 → ~15 jobs for 105 loci)
#   --chunk-size N   Variants per Evo2 inference chunk (default: 10)
#   --dry-run        Print qsub commands without submitting

BASE=/grand/GeomicVar/mconery/evo2_variant_scoring_mapping
REPO=/lus/grand/projects/GeomicVar/mconery/Evo2-MVP-Variant-Scoring
LOCI_FILE=${BASE}/loci_definition/t2d_eur_loci.tsv
WORKER=${REPO}/Evo2_Priors_FM/evo2_scoring/run_evo2_worker.sh
RESULTS_DIR=${BASE}/evo2_scoring/results
LOG_DIR=${BASE}/evo2_scoring/logs

BATCH_SIZE=7
CHUNK_SIZE=10
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --batch-size) BATCH_SIZE="$2"; shift 2 ;;
        --chunk-size) CHUNK_SIZE="$2"; shift 2 ;;
        --dry-run)    DRY_RUN=1; shift ;;
        *) shift ;;
    esac
done

mkdir -p "${RESULTS_DIR}" "${LOG_DIR}"

if [ ! -f "${LOCI_FILE}" ]; then
    echo "ERROR: Loci file not found: ${LOCI_FILE}"
    exit 1
fi

# Collect loci that still need scoring
pending=()
skipped=0

while IFS=$'\t' read -r locus_id chr chr_full start end; do
    [ "${locus_id}" = "locus_id" ] && continue

    out_file="${RESULTS_DIR}/${locus_id}.evo2_scores.csv"
    variant_file="${BASE}/variant_lists/${locus_id}.variants.tsv"

    if [ -f "${out_file}" ] && [ -f "${variant_file}" ]; then
        n_vars=$(( $(wc -l < "${variant_file}") - 1 ))
        n_done=$(( $(wc -l < "${out_file}") - 1 ))
        if [ "${n_done}" -ge "${n_vars}" ]; then
            skipped=$((skipped + 1))
            continue
        fi
    fi

    pending+=("${locus_id}")
done < "${LOCI_FILE}"

echo "Pending loci: ${#pending[@]}  Already complete: ${skipped}"

if [ ${#pending[@]} -eq 0 ]; then
    echo "Nothing to submit."
    exit 0
fi

# Split pending loci into batches and submit one job per batch
batch_num=0
submitted=0
i=0

while [ $i -lt ${#pending[@]} ]; do
    batch_num=$((batch_num + 1))
    batch=("${pending[@]:$i:$BATCH_SIZE}")
    i=$((i + BATCH_SIZE))

    # Join locus IDs with colons for passing via -v
    locus_ids=$(IFS=:; echo "${batch[*]}")
    first_locus="${batch[0]}"
    last_locus="${batch[-1]}"
    n_in_batch="${#batch[@]}"

    job_name="evo2_b${batch_num}"

    cmd="qsub \
        -v LOCUS_IDS=${locus_ids},chunk_size=${CHUNK_SIZE},tp_size=4,cp_size=1,window_size=16384,MODEL_SIZE=7b_arc_longcontext \
        -N ${job_name} \
        ${WORKER}"

    if [ "${DRY_RUN}" = "1" ]; then
        echo "[DRY RUN] Batch ${batch_num} (${n_in_batch} loci: ${first_locus} ... ${last_locus})"
        echo "  ${cmd}"
    else
        job_id=$(eval "${cmd}" 2>&1)
        echo "Batch ${batch_num} (${n_in_batch} loci, ${first_locus}..${last_locus}): ${job_id}"
        submitted=$((submitted + 1))
    fi
done

echo ""
echo "Submitted: ${submitted} jobs covering ${#pending[@]} loci  |  Skipped: ${skipped} already complete"
