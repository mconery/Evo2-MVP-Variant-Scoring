#!/bin/bash

################################################################################################################
########################################## Step 1: Sample variants #############################################
################################################################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SIF_PATH=/vast/projects/anuragv/cohort/mconery/bionemo/bionemo-nightly.sif
BIND_PATH=/vast/projects/anuragv/cohort/mconery:/vast/projects/anuragv/cohort/mconery

SAMPLE_JOB=$(sbatch --parsable \
    --job-name=TIMING_sample_variants \
    --partition=genoa-std-mem \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=4 \
    --mem=16G \
    --time=0:30:00 \
    --output="${SCRIPT_DIR}/sample_variants_%j.log" \
    --wrap="bash -c 'source /etc/profile.d/modules.sh && module load apptainer/1.4.1 && apptainer exec --bind ${BIND_PATH} ${SIF_PATH} python ${SCRIPT_DIR}/sample_variants.py'")

echo "Submitted sample_variants job: ${SAMPLE_JOB}"

################################################################################################################
########## Step 2: 32 combinations (8192bp context) — 2 models x 4 chunk sizes x 4 TP/CP pairs ################
################################################################################################################

MODELS=("7b_arc_longcontext" "40b_arc_longcontext")
CHUNK_SIZES=(10 20 50 100)
TP_CP_PAIRS=("8 1" "4 2" "2 4" "1 8")
WINDOW=8192
TIME_LIMIT="6:00:00"

for MODEL in "${MODELS[@]}"; do
    for CHUNK in "${CHUNK_SIZES[@]}"; do
        for TP_CP in "${TP_CP_PAIRS[@]}"; do
            TP=$(echo $TP_CP | awk '{print $1}')
            CP=$(echo $TP_CP | awk '{print $2}')
            JOB_NAME="TIMING_${MODEL}_tp${TP}_cp${CP}_chunk${CHUNK}_${WINDOW}bp"
            echo "Submitting: ${JOB_NAME}"
            sbatch --export=ALL,MODEL_SIZE=${MODEL},window_size=${WINDOW},tp_size=${TP},cp_size=${CP},chunk_size=${CHUNK} \
                --nodes=1 \
                --gpus=8 \
                --time=${TIME_LIMIT} \
                --job-name=${JOB_NAME} \
                --dependency=afterok:${SAMPLE_JOB} \
                "${SCRIPT_DIR}/run_timing_evo2_worker.sh"
        done
    done
done

################################################################################################################
############## Step 3: 4 context-length jobs — 7b_arc_longcontext, TP=4, CP=2, chunk=100 ######################
################################################################################################################

MODEL=7b_arc_longcontext
TP=4
CP=2
CHUNK_SIZE=100
CONTEXT_SIZES=(16384 65536 131072 524288)

for CTX in "${CONTEXT_SIZES[@]}"; do
    JOB_NAME="TIMING_${MODEL}_tp${TP}_cp${CP}_chunk${CHUNK_SIZE}_${CTX}bp"
    echo "Submitting: ${JOB_NAME}"
    sbatch --export=ALL,MODEL_SIZE=${MODEL},window_size=${CTX},tp_size=${TP},cp_size=${CP},chunk_size=${CHUNK_SIZE} \
        --nodes=1 \
        --gpus=8 \
        --time=${TIME_LIMIT} \
        --job-name=${JOB_NAME} \
        --dependency=afterok:${SAMPLE_JOB} \
        "${SCRIPT_DIR}/run_timing_evo2_worker.sh"
done

echo "All timing jobs submitted. Sample job ID: ${SAMPLE_JOB}"
