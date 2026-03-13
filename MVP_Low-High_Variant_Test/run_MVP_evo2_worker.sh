#!/bin/bash

#SBATCH -o run_MVP_evo2_%x.log
#SBATCH -p dgx-b200
#SBATCH --gpus=16
#SBATCH -t 48:00:00
#SBATCH -N 1

################################################################################################################
########################################## Define key model parameters #########################################
################################################################################################################

# MODEL_SIZE, window_size, tp_size, cp_size, chunk_size are passed via --export from launcher script
: "${chunk_size:=100}"
: "${tp_size:=1}"
: "${cp_size:=1}"

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

out_file=/vast/projects/anuragv/cohort/mconery/mvp_variant_test/MVP_variant_scores."$MODEL_SIZE"_model."$window_size"bp_context.csv

evotwo_script="/vast/home/m/mconery/Evo2-TopMed-Variant-Scoring/MVP_Low-High_Variant_Test/mvp_variants_test.py"

SIF_PATH=/vast/projects/anuragv/cohort/mconery/bionemo/bionemo-nightly.sif
BIND_PATH=/vast/projects/anuragv/cohort/mconery:/vast/projects/anuragv/cohort/mconery
APPTAINER_CMD="apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH}"

TOTAL_GPUS=$((tp_size * cp_size))
GPUS_PER_NODE=$((TOTAL_GPUS / ${SLURM_NNODES:-1}))

output_dir=/vast/projects/anuragv/cohort/mconery/mvp_variant_test
ref_fasta_path=${output_dir}/temp_ref.${MODEL_SIZE}.${window_size}bp.fa
var_fasta_path=${output_dir}/temp_var.${MODEL_SIZE}.${window_size}bp.fa
predict_ref_dir=${output_dir}/reference_predictions.${MODEL_SIZE}.${window_size}bp
predict_var_dir=${output_dir}/variant_predictions.${MODEL_SIZE}.${window_size}bp

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

    # Run distributed inference across both nodes for REF sequences
    srun --ntasks=${TOTAL_GPUS} --ntasks-per-node=${GPUS_PER_NODE} \
        apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
        predict_evo2 --fasta ${ref_fasta_path} \
        --ckpt-dir ${CHECKPOINT_PATH} \
        --output-dir ${predict_ref_dir} \
        --model-size ${MODEL_SIZE} \
        --tensor-parallel-size ${tp_size} \
        --pipeline-model-parallel-size 1 \
        --context-parallel-size ${cp_size} \
        --output-log-prob-seqs ${FP8_FLAG} || break

    # Run distributed inference across both nodes for ALT sequences
    srun --ntasks=${TOTAL_GPUS} --ntasks-per-node=${GPUS_PER_NODE} \
        apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
        predict_evo2 --fasta ${var_fasta_path} \
        --ckpt-dir ${CHECKPOINT_PATH} \
        --output-dir ${predict_var_dir} \
        --model-size ${MODEL_SIZE} \
        --tensor-parallel-size ${tp_size} \
        --pipeline-model-parallel-size 1 \
        --context-parallel-size ${cp_size} \
        --output-log-prob-seqs ${FP8_FLAG} || break

    # Process results and write to CSV
    ${APPTAINER_CMD} python $evotwo_script \
        --mode process --chunk-start ${chunk_start} --chunk-end ${CHUNK_END} \
        --out $out_file --ctx $window_size --model $MODEL_SIZE \
        --tensor-parallel-size $tp_size --context-parallel-size $cp_size \
        --chunk-size $chunk_size

    chunk_start=${CHUNK_END}
done
