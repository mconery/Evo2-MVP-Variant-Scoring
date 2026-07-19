#!/bin/bash

#SBATCH -o run_evo2_priors_%x.log
#SBATCH -p dgx-b200
#SBATCH --gpus=8
#SBATCH -t 48:00:00
#SBATCH -N 1

###############################################################################
# Evo2 scoring worker for T2D fine-mapping priors — runs on PARCC Betty.
#
# Scores the entire merged variant list (see merge_variant_lists.py) in one
# continuous chunked loop. chunk_size/tp_size/cp_size/window_size/MODEL_SIZE
# are normally passed via `sbatch --export=...` from launch_evo2_jobs.sh.
###############################################################################
: "${chunk_size:=100}"
: "${tp_size:=8}"
: "${cp_size:=1}"
: "${window_size:=8192}"
: "${MODEL_SIZE:=7b_arc_longcontext}"

###############################################################################
# Paths
###############################################################################
BASE=/vast/projects/anuragv/cohort/mconery/evo2_variant_scoring_mapping
REPO=/vast/home/m/mconery/Evo2-TopMed-Variant-Scoring

FASTA=${BASE}/reference/GRCh37.p13.genome.fa
SCORING_SCRIPT=${REPO}/Evo2_Priors_FM/evo2_scoring/score_locus_variants.py
SIF_PATH=/vast/projects/anuragv/cohort/mconery/bionemo/bionemo-nightly.sif
BIND_PATH=/vast/projects/anuragv/cohort/mconery:/vast/projects/anuragv/cohort/mconery

APPTAINER_CMD="apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH}"
TOTAL_GPUS=$((tp_size * cp_size))

VARIANT_FILE=${BASE}/variant_lists/merged_variants.tsv
OUT_FILE=${BASE}/evo2_scoring/results/merged_evo2_scores.csv

mkdir -p "${BASE}/evo2_scoring/results"

###############################################################################
# Environment setup
###############################################################################
module load cuda/12.8.1 cudnn/8.9.7.29-12 apptainer/1.4.1

export HF_HOME=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NEMO_CACHE_DIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export BIONEMO_CACHE_DIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NEMO_MODELS_CACHE=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NGC_CLI_API_KEY="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"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

###############################################################################
# Validate
###############################################################################
if [ ! -f "${VARIANT_FILE}" ]; then
    echo "ERROR: Merged variant file not found: ${VARIANT_FILE}" >&2
    echo "       Run merge_variant_lists.py on Polaris and transfer it here first." >&2
    exit 1
fi
if [ ! -f "${FASTA}" ]; then
    echo "ERROR: FASTA not found: ${FASTA}" >&2
    echo "       Transfer the decompressed+indexed GRCh37.p13 FASTA from Polaris first." >&2
    exit 1
fi

echo "=== Scoring merged variant list: ${VARIANT_FILE} ===" >&2
echo "chunk_size=${chunk_size}  tp=${tp_size}  cp=${cp_size}  ctx=${window_size}bp  GPUs=${TOTAL_GPUS}" >&2

###############################################################################
# Chunked scoring loop
###############################################################################
chunk_start=0
FP8_FLAG=""
FIRST_CHUNK=1

while true; do
    prep_output=$(${APPTAINER_CMD} python ${SCORING_SCRIPT} \
        --mode prepare \
        --chunk-start ${chunk_start} \
        --variants "${VARIANT_FILE}" \
        --fasta "${FASTA}" \
        --out "${OUT_FILE}" \
        --ctx ${window_size} \
        --model ${MODEL_SIZE} \
        --tensor-parallel-size ${tp_size} \
        --context-parallel-size ${cp_size} \
        --chunk-size ${chunk_size})

    exit_code=$?
    [ $exit_code -ne 0 ] && break   # exit code 1 = all chunks done

    eval "$prep_output"    # sets FP8_FLAG, CHUNK_END, RESUME_FROM, CHECKPOINT_PATH

    # On first iteration, honour any resume offset
    if [ "${FIRST_CHUNK}" = "1" ]; then
        FIRST_CHUNK=0
        chunk_start=${RESUME_FROM:-0}
        if [ "${chunk_start}" != "0" ]; then
            prep_output=$(${APPTAINER_CMD} python ${SCORING_SCRIPT} \
                --mode prepare \
                --chunk-start ${chunk_start} \
                --variants "${VARIANT_FILE}" \
                --fasta "${FASTA}" \
                --out "${OUT_FILE}" \
                --ctx ${window_size} \
                --model ${MODEL_SIZE} \
                --tensor-parallel-size ${tp_size} \
                --context-parallel-size ${cp_size} \
                --chunk-size ${chunk_size})
            exit_code=$?
            [ $exit_code -ne 0 ] && break
            eval "$prep_output"
        fi
    fi

    output_dir=$(dirname "${OUT_FILE}")
    ref_fasta_path=${output_dir}/temp_ref.${MODEL_SIZE}.${window_size}bp.fa
    var_fasta_path=${output_dir}/temp_var.${MODEL_SIZE}.${window_size}bp.fa
    predict_ref_dir=${output_dir}/reference_predictions.${MODEL_SIZE}.${window_size}bp
    predict_var_dir=${output_dir}/variant_predictions.${MODEL_SIZE}.${window_size}bp

    # REF inference
    srun --ntasks=${TOTAL_GPUS} --ntasks-per-node=${TOTAL_GPUS} \
        apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
        predict_evo2 \
            --fasta "${ref_fasta_path}" \
            --ckpt-dir "${CHECKPOINT_PATH}" \
            --output-dir "${predict_ref_dir}" \
            --model-size ${MODEL_SIZE} \
            --tensor-parallel-size ${tp_size} \
            --pipeline-model-parallel-size 1 \
            --context-parallel-size ${cp_size} \
            --output-log-prob-seqs ${FP8_FLAG} || break

    # ALT inference
    srun --ntasks=${TOTAL_GPUS} --ntasks-per-node=${TOTAL_GPUS} \
        apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
        predict_evo2 \
            --fasta "${var_fasta_path}" \
            --ckpt-dir "${CHECKPOINT_PATH}" \
            --output-dir "${predict_var_dir}" \
            --model-size ${MODEL_SIZE} \
            --tensor-parallel-size ${tp_size} \
            --pipeline-model-parallel-size 1 \
            --context-parallel-size ${cp_size} \
            --output-log-prob-seqs ${FP8_FLAG} || break

    # Process results
    ${APPTAINER_CMD} python ${SCORING_SCRIPT} \
        --mode process \
        --chunk-start ${chunk_start} \
        --chunk-end ${CHUNK_END} \
        --variants "${VARIANT_FILE}" \
        --fasta "${FASTA}" \
        --out "${OUT_FILE}" \
        --ctx ${window_size} \
        --model ${MODEL_SIZE} \
        --tensor-parallel-size ${tp_size} \
        --context-parallel-size ${cp_size} \
        --chunk-size ${chunk_size} || break

    chunk_start=${CHUNK_END}
done

N_SCORED=0
[ -f "${OUT_FILE}" ] && N_SCORED=$(( $(wc -l < "${OUT_FILE}") - 1 ))
echo "=== Done: ${N_SCORED} variants scored, written to ${OUT_FILE} ===" >&2
