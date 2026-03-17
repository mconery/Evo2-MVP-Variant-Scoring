#!/bin/bash

################################################################################################################
########################################## Define parameter combinations #######################################
################################################################################################################

# Define arrays of model sizes and context window sizes to test
MODEL_SIZES=("7b" "40b" "40b_arc_longcontext" "7b_arc_longcontext")
WINDOW_SIZES=(8192 16384 65536 131072 524288)

################################################################################################################
############################ Define per-model parallelism configurations #######################################
################################################################################################################

# get_model_config MODEL_SIZE -> sets TP, CP, NUM_GPUS
# All models use the same parallelism config
get_model_config() {
    TP=8; CP=1; NUM_GPUS=8
}

# get_chunk_size WINDOW_SIZE -> sets CHUNK_SIZE
get_chunk_size() {
    CHUNK_SIZE=40
}

# get_time_limit WINDOW_SIZE -> sets TIME_LIMIT
get_time_limit() {
    local ws="$1"
    case "$ws" in
        8192)    TIME_LIMIT="8:00:00"  ;;
        16384)   TIME_LIMIT="16:00:00" ;;
        65536)   TIME_LIMIT="24:00:00" ;;
        131072)  TIME_LIMIT="24:00:00" ;;
        524288)  TIME_LIMIT="36:00:00" ;;
        1000000) TIME_LIMIT="48:00:00" ;;
        *)       TIME_LIMIT="48:00:00" ;;
    esac
}

################################################################################################################
##################################### Submit jobs for each combination #########################################
################################################################################################################

for MODEL_SIZE in "${MODEL_SIZES[@]}"; do
    get_model_config "$MODEL_SIZE"
    for window_size in "${WINDOW_SIZES[@]}"; do
        get_time_limit "$window_size"
        get_chunk_size "$window_size"
        echo "Submitting job for MODEL_SIZE=$MODEL_SIZE, window_size=$window_size, TP=$TP, CP=$CP, GPUs=$NUM_GPUS, TIME=$TIME_LIMIT, CHUNK_SIZE=$CHUNK_SIZE"
        # Create unique job name
        JOB_NAME="MVP_${MODEL_SIZE}_${window_size}bp"
        # Submit the job
        sbatch --export=ALL,MODEL_SIZE=$MODEL_SIZE,window_size=$window_size,tp_size=$TP,cp_size=$CP,chunk_size=$CHUNK_SIZE \
            --nodes=1 \
	    --gpus=$NUM_GPUS \
            --time=$TIME_LIMIT \
            --job-name=$JOB_NAME \
            run_MVP_evo2_worker.sh
    done
done

echo "All jobs submitted successfully"
