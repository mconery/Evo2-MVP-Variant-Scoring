#!/bin/bash
# Submit the chunk-size / GPU-memory saturation sweep for the Evo2 T2D
# fine-mapping priors pipeline on PARCC Betty.
#
# Scores real production variants (merged_variants.tsv) into
# results/merged_evo2_scores.csv while checking whether chunk_size=100 has
# saturated GPU memory -- walks chunk_size up through a candidate list
# (default 100/200/400/800/1600/3200), running --chunks-per-phase real
# chunks at each size before advancing. Fixed at ctx=8192bp, TP=8, CP=1,
# model=7b_arc_longcontext (matches production).
#
# On OOM (or any other failure) the job stops immediately -- no
# auto-fallback -- and leaves the remaining variants unscored. Review
# results/chunk_size_sweep_log.csv, then resume full production scoring
# with run_evo2_worker.sh at the last chunk_size that succeeded.
#
# Usage:
#   bash launch_evo2_chunk_sweep.sh [--chunk-sizes "100 200 400 800 1600 3200"] \
#                                    [--chunks-per-phase N] [--time HH:MM:SS] \
#                                    [--exclude NODE1,NODE2] [--dry-run]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CHUNK_SIZES="100 200 400 800 1600 3200"
CHUNKS_PER_PHASE=2
TIME_LIMIT="8:00:00"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --chunk-sizes)      CHUNK_SIZES="$2"; shift 2 ;;
        --chunks-per-phase) CHUNKS_PER_PHASE="$2"; shift 2 ;;
        --time)             TIME_LIMIT="$2"; shift 2 ;;
        --exclude)          EXCLUDE_NODES="$2"; shift 2 ;;
        --dry-run)          DRY_RUN=1; shift ;;
        *) shift ;;
    esac
done

JOB_NAME="evo2_chunk_sweep_7b_arc_longcontext_tp8_cp1_8192bp"

cmd="sbatch --parsable --export=ALL,CHUNK_SIZES='${CHUNK_SIZES}',CHUNKS_PER_PHASE=${CHUNKS_PER_PHASE} \
    --nodes=1 \
    --gpus=8 \
    --time=${TIME_LIMIT} \
    --job-name=${JOB_NAME} \
    ${SCRIPT_DIR}/run_evo2_chunk_sweep.sh"

echo "chunk_sizes=[${CHUNK_SIZES}] chunks_per_phase=${CHUNKS_PER_PHASE} time=${TIME_LIMIT} exclude=${EXCLUDE_NODES}"

if [ "${DRY_RUN}" = "1" ]; then
    echo "[DRY RUN] ${cmd}"
else
    job_id=$(eval "${cmd}" 2>&1)
    echo "Submitted: ${job_id}"
fi
