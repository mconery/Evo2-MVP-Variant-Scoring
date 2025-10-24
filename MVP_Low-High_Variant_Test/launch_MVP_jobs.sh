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
MODEL_SIZES=(1b 7b)
WINDOW_SIZES=(8192 16284 50000 100000)

################################################################################################################
##################################### Submit jobs for each combination #########################################
################################################################################################################

for MODEL_SIZE in "${MODEL_SIZES[@]}"; do
    for window_size in "${WINDOW_SIZES[@]}"; do
	#Check if file exists already
	out_file=/vast/projects/anuragv/cohort/mconery/mvp_variant_test/MVP_variant_scores."$MODEL_SIZE"_model."$window_size"bp_context.csv
	if [[ ! -f "$out_file" ]]; then	
        	echo "Submitting job for MODEL_SIZE=$MODEL_SIZE, window_size=$window_size"
        	# Create unique job name
       		JOB_NAME="MVP_${MODEL_SIZE}_${window_size}bp"
		# Submit the job
        	sbatch --export=MODEL_SIZE=$MODEL_SIZE,window_size=$window_size \
               		--job-name=$JOB_NAME \
               		run_MVP_evo2_worker.sh
    done
done

echo "All jobs submitted successfully"
