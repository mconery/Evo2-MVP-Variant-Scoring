#!/bin/bash

#SBATCH -J launch_MVP_jobs
#SBATCH -o launch_MVP_jobs.log
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH -t 0:10:00

################################################################################################################
########################################## Define parameter combinations #######################################
################################################################################################################

# Define arrays of model sizes and context window sizes to test
MODEL_SIZES=("1b" "7b" "40b" "7b_arc_longcontext" "40b_arc_longcontext")
WINDOW_SIZES=(8192 16284 50000 100000 500000 1000000)

################################################################################################################
############################ Define per-model parallelism configurations #######################################
################################################################################################################

# get_model_config MODEL_SIZE -> sets TP, CP, NUM_GPUS
# All models use the same parallelism config
get_model_config() {
    TP=4; CP=2; NUM_GPUS=8
}

# get_time_limit WINDOW_SIZE -> sets TIME_LIMIT
get_time_limit() {
    local ws="$1"
    case "$ws" in
        8192)    TIME_LIMIT="8:00:00"  ;;
        16284)   TIME_LIMIT="16:00:00" ;;
        50000)   TIME_LIMIT="24:00:00" ;;
        100000)  TIME_LIMIT="24:00:00" ;;
        500000)  TIME_LIMIT="48:00:00" ;;
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
        echo "Submitting job for MODEL_SIZE=$MODEL_SIZE, window_size=$window_size, TP=$TP, CP=$CP, GPUs=$NUM_GPUS, TIME=$TIME_LIMIT"
        # Create unique job name
        JOB_NAME="MVP_${MODEL_SIZE}_${window_size}bp"
        # Submit the job
        sbatch --export=MODEL_SIZE=$MODEL_SIZE,window_size=$window_size,tp_size=$TP,cp_size=$CP \
            --gpus=$NUM_GPUS \
            --time=$TIME_LIMIT \
            --job-name=$JOB_NAME \
            run_MVP_evo2_worker.sh
    done
done

echo "All jobs submitted successfully"
