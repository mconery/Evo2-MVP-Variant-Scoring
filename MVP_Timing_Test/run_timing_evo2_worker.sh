#!/bin/bash

#SBATCH -o run_timing_evo2_%x.log
#SBATCH -p dgx-b200
#SBATCH --gpus=8
#SBATCH -t 6:00:00
#SBATCH -N 1

################################################################################################################
########################################## Define key model parameters #########################################
################################################################################################################

# MODEL_SIZE, window_size, tp_size, cp_size, chunk_size are passed via --export from launcher script
: "${chunk_size:=100}"
: "${tp_size:=8}"
: "${cp_size:=1}"

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

output_dir=/vast/projects/anuragv/cohort/mconery/mvp_timing_test

out_file=${output_dir}/timing_scores.${MODEL_SIZE}_tp${tp_size}_cp${cp_size}_chunk${chunk_size}_${window_size}bp.csv

evotwo_script="/vast/home/m/mconery/Evo2-TopMed-Variant-Scoring/MVP_Timing_Test/timing_variants_test.py"

SIF_PATH=/vast/projects/anuragv/cohort/mconery/bionemo/bionemo-nightly.sif
BIND_PATH=/vast/projects/anuragv/cohort/mconery:/vast/projects/anuragv/cohort/mconery
APPTAINER_CMD="apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH}"

TOTAL_GPUS=$((tp_size * cp_size))
GPUS_PER_NODE=$((TOTAL_GPUS / ${SLURM_NNODES:-1}))

ref_fasta_path=${output_dir}/temp_ref.${MODEL_SIZE}.tp${tp_size}_cp${cp_size}_chunk${chunk_size}.${window_size}bp.fa
var_fasta_path=${output_dir}/temp_var.${MODEL_SIZE}.tp${tp_size}_cp${cp_size}_chunk${chunk_size}.${window_size}bp.fa
predict_ref_dir=${output_dir}/reference_predictions.${MODEL_SIZE}.tp${tp_size}_cp${cp_size}_chunk${chunk_size}.${window_size}bp
predict_var_dir=${output_dir}/variant_predictions.${MODEL_SIZE}.tp${tp_size}_cp${cp_size}_chunk${chunk_size}.${window_size}bp

################################################################################################################
############################################# Call Evo2 Scoring Script #########################################
################################################################################################################

#Load Modules and set caches
module load cuda/12.8.1 cudnn/8.9.7.29-12 apptainer/1.4.1

export HF_HOME=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NEMO_CACHE_DIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export BIONEMO_CACHE_DIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NEMO_MODELS_CACHE=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NGC_CLI_API_KEY="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"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

chunk_start=0
FP8_FLAG=""
FIRST_CHUNK=1
JOB_STATUS=SUCCESS
FAIL_REASON=""

START_TIME=$(date +%s)

while true; do
    # Prepare FASTAs and detect FP8 / resume point for this chunk
    prep_output=$(${APPTAINER_CMD} python $evotwo_script \
        --mode prepare --chunk-start ${chunk_start} \
        --out $out_file --ctx $window_size --model $MODEL_SIZE \
        --tensor-parallel-size $tp_size --context-parallel-size $cp_size \
        --chunk-size $chunk_size)

    exit_code=$?
    [ $exit_code -ne 0 ] && break   # exit code 1 = all chunks done

    eval "$prep_output"             # sets FP8_FLAG, CHUNK_END, RESUME_FROM, CHECKPOINT_PATH

    # On first iteration, update chunk_start to the actual resume point if needed
    if [ "${FIRST_CHUNK}" = "1" ]; then
        FIRST_CHUNK=0
        chunk_start=${RESUME_FROM:-0}
        # Re-prepare with the actual starting chunk if we resumed mid-way
        if [ "${chunk_start}" != "0" ]; then
            prep_output=$(${APPTAINER_CMD} python $evotwo_script \
                --mode prepare --chunk-start ${chunk_start} \
                --out $out_file --ctx $window_size --model $MODEL_SIZE \
                --tensor-parallel-size $tp_size --context-parallel-size $cp_size \
                --chunk-size $chunk_size)
            exit_code=$?
            [ $exit_code -ne 0 ] && break
            eval "$prep_output"
        fi
    fi

    # Run distributed inference for REF sequences
    srun --ntasks=${TOTAL_GPUS} --ntasks-per-node=${GPUS_PER_NODE} \
        apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
        predict_evo2 --fasta ${ref_fasta_path} \
        --ckpt-dir ${CHECKPOINT_PATH} \
        --output-dir ${predict_ref_dir} \
        --model-size ${MODEL_SIZE} \
        --tensor-parallel-size ${tp_size} \
        --pipeline-model-parallel-size 1 \
        --context-parallel-size ${cp_size} \
        --output-log-prob-seqs ${FP8_FLAG}
    srun_exit=$?
    if [ $srun_exit -ne 0 ]; then
        JOB_STATUS=FAILED
        FAIL_REASON="srun predict_evo2 (REF) failed with exit code ${srun_exit} at chunk_start=${chunk_start}"
        break
    fi

    # Run distributed inference for ALT sequences
    srun --ntasks=${TOTAL_GPUS} --ntasks-per-node=${GPUS_PER_NODE} \
        apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
        predict_evo2 --fasta ${var_fasta_path} \
        --ckpt-dir ${CHECKPOINT_PATH} \
        --output-dir ${predict_var_dir} \
        --model-size ${MODEL_SIZE} \
        --tensor-parallel-size ${tp_size} \
        --pipeline-model-parallel-size 1 \
        --context-parallel-size ${cp_size} \
        --output-log-prob-seqs ${FP8_FLAG}
    srun_exit=$?
    if [ $srun_exit -ne 0 ]; then
        JOB_STATUS=FAILED
        FAIL_REASON="srun predict_evo2 (ALT) failed with exit code ${srun_exit} at chunk_start=${chunk_start}"
        break
    fi

    # Process results and write to CSV
    ${APPTAINER_CMD} python $evotwo_script \
        --mode process --chunk-start ${chunk_start} --chunk-end ${CHUNK_END} \
        --out $out_file --ctx $window_size --model $MODEL_SIZE \
        --tensor-parallel-size $tp_size --context-parallel-size $cp_size \
        --chunk-size $chunk_size

    chunk_start=${CHUNK_END}
done

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

################################################################################################################
############################################### Record timing results ##########################################
################################################################################################################

TIMING_CSV="${output_dir}/timing_results.csv"
HEADER="model,context_bp,tp_size,cp_size,chunk_size,start_time,end_time,duration_seconds,status,fail_reason"
ROW="${MODEL_SIZE},${window_size},${tp_size},${cp_size},${chunk_size},${START_TIME},${END_TIME},${DURATION},${JOB_STATUS},${FAIL_REASON}"

# Atomic append with flock; write header only if file does not yet exist
(
    flock -x 200
    if [ ! -f "${TIMING_CSV}" ]; then
        echo "${HEADER}" >> "${TIMING_CSV}"
    fi
    echo "${ROW}" >> "${TIMING_CSV}"
) 200>"${TIMING_CSV}.lock"
