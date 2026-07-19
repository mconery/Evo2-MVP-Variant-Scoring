#!/bin/bash
# Submit the Evo2 scoring job for T2D fine-mapping priors on PARCC Betty.
#
# Scores the whole merged variant list (merge_variant_lists.py, run on
# Polaris and transferred here) in a single SLURM job.
#
# Usage:
#   bash launch_evo2_jobs.sh [--chunk-size N] [--tp-size N] [--cp-size N] \
#                             [--window-size N] [--model NAME] [--dry-run]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CHUNK_SIZE=100
TP_SIZE=8
CP_SIZE=1
WINDOW_SIZE=8192
MODEL_SIZE=7b_arc_longcontext
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --chunk-size)  CHUNK_SIZE="$2"; shift 2 ;;
        --tp-size)     TP_SIZE="$2"; shift 2 ;;
        --cp-size)     CP_SIZE="$2"; shift 2 ;;
        --window-size) WINDOW_SIZE="$2"; shift 2 ;;
        --model)       MODEL_SIZE="$2"; shift 2 ;;
        --dry-run)     DRY_RUN=1; shift ;;
        *) shift ;;
    esac
done

NUM_GPUS=$((TP_SIZE * CP_SIZE))
JOB_NAME="evo2_priors_${MODEL_SIZE}_${WINDOW_SIZE}bp"

cmd="sbatch --export=ALL,MODEL_SIZE=${MODEL_SIZE},window_size=${WINDOW_SIZE},tp_size=${TP_SIZE},cp_size=${CP_SIZE},chunk_size=${CHUNK_SIZE} \
    --nodes=1 \
    --gpus=${NUM_GPUS} \
    --job-name=${JOB_NAME} \
    ${SCRIPT_DIR}/run_evo2_worker.sh"

echo "MODEL_SIZE=${MODEL_SIZE} window_size=${WINDOW_SIZE}bp TP=${TP_SIZE} CP=${CP_SIZE} GPUs=${NUM_GPUS} chunk_size=${CHUNK_SIZE}"

if [ "${DRY_RUN}" = "1" ]; then
    echo "[DRY RUN] ${cmd}"
else
    job_id=$(eval "${cmd}" 2>&1)
    echo "Submitted: ${job_id}"
fi
