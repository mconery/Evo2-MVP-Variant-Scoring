#!/bin/bash
#PBS -l select=1:ngpus=4:system=polaris
#PBS -l walltime=24:00:00
#PBS -l filesystems=grand:home
#PBS -q preemptable
#PBS -A GeomicVar
#PBS -o /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/evo2_scoring/logs/worker_${PBS_JOBID}.log
#PBS -e /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/evo2_scoring/logs/worker_${PBS_JOBID}.err

###############################################################################
# Parameters — passed via -v when submitting:
#   LOCUS_IDS   colon-separated list of locus IDs to process in sequence
#   chunk_size  variants per Evo2 inference chunk (default: 10)
#   tp_size     tensor parallel degree (default: 4)
#   cp_size     context parallel degree (default: 1)
#   window_size context window in bp (default: 16384)
#   MODEL_SIZE  model checkpoint name (default: 7b_arc_longcontext)
###############################################################################
: "${chunk_size:=10}"
: "${tp_size:=4}"
: "${cp_size:=1}"
: "${window_size:=16384}"
: "${MODEL_SIZE:=7b_arc_longcontext}"

###############################################################################
# Paths
###############################################################################
BASE=/grand/GeomicVar/mconery/evo2_variant_scoring_mapping
REPO=/lus/grand/projects/GeomicVar/mconery/Evo2-MVP-Variant-Scoring

FASTA=${BASE}/reference/GRCh37.p13.genome.fa
SCORING_SCRIPT=${REPO}/Evo2_Priors_FM/evo2_scoring/score_locus_variants.py
SIF_PATH=/lus/grand/projects/GeomicVar/mconery/tools/bionemo-nightly.sif
BIND_PATH=/grand/GeomicVar:/grand/GeomicVar,/lus/grand/projects/GeomicVar:/lus/grand/projects/GeomicVar

APPTAINER_CMD="apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH}"
TOTAL_GPUS=$((tp_size * cp_size))

TIMING_LOG=${BASE}/evo2_scoring/results/timing_log.csv
LOCK_FILE=${BASE}/evo2_scoring/results/timing_log.lock

###############################################################################
# Environment setup
###############################################################################
module load cuda cudnn apptainer 2>/dev/null || true

export HF_HOME=/lus/grand/projects/GeomicVar/mconery/tools/bionemo_cache
export NEMO_CACHE_DIR=/lus/grand/projects/GeomicVar/mconery/tools/bionemo_cache
export BIONEMO_CACHE_DIR=/lus/grand/projects/GeomicVar/mconery/tools/bionemo_cache
export NEMO_MODELS_CACHE=/lus/grand/projects/GeomicVar/mconery/tools/bionemo_cache
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# NGC credentials for model weight download (first run only)
NGC_CREDS=/grand/GeomicVar/mconery/tools/ngc_credentials.sh
if [ -f "${NGC_CREDS}" ]; then
    source "${NGC_CREDS}"
else
    echo "WARNING: NGC credentials file not found at ${NGC_CREDS}" >&2
    echo "         Model weight download will fail if weights are not already cached." >&2
fi

mkdir -p ${BASE}/evo2_scoring/results ${BASE}/evo2_scoring/logs

###############################################################################
# Validate
###############################################################################
if [ -z "${LOCUS_IDS}" ]; then
    echo "ERROR: LOCUS_IDS not set. Submit with: qsub -v LOCUS_IDS=locus1:locus2:..." >&2
    exit 1
fi
if [ ! -f "${FASTA}" ]; then
    echo "ERROR: FASTA not found: ${FASTA}. Run setup/decompress_fasta.sh first." >&2
    exit 1
fi

GPU_MODEL=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || echo "unknown")

###############################################################################
# Loop over loci
###############################################################################
IFS=':' read -ra LOCI_ARRAY <<< "${LOCUS_IDS}"
echo "=== Batch: ${#LOCI_ARRAY[@]} loci  GPU: ${GPU_MODEL} ===" >&2

for LOCUS_ID in "${LOCI_ARRAY[@]}"; do

    VARIANT_FILE=${BASE}/variant_lists/${LOCUS_ID}.variants.tsv
    OUT_FILE=${BASE}/evo2_scoring/results/${LOCUS_ID}.evo2_scores.csv

    # Skip if already fully scored
    if [ -f "${OUT_FILE}" ] && [ -f "${VARIANT_FILE}" ]; then
        n_vars=$(( $(wc -l < "${VARIANT_FILE}") - 1 ))
        n_done=$(( $(wc -l < "${OUT_FILE}") - 1 ))
        if [ "${n_done}" -ge "${n_vars}" ]; then
            echo "SKIP (complete): ${LOCUS_ID} (${n_done}/${n_vars} variants)" >&2
            continue
        fi
    fi

    if [ ! -f "${VARIANT_FILE}" ]; then
        echo "SKIP (no variant file): ${LOCUS_ID}" >&2
        continue
    fi

    echo "=== Scoring: ${LOCUS_ID} ===" >&2
    echo "chunk_size=${chunk_size}  tp=${tp_size}  cp=${cp_size}  ctx=${window_size}bp" >&2

    chunk_start=0
    FP8_FLAG=""
    FIRST_CHUNK=1
    LOCUS_STATUS=SUCCESS
    REF_SEC=NA
    VAR_SEC=NA

    LOCUS_START_TIME=$(date +%s)

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
        if [ $exit_code -ne 0 ]; then
            LOCUS_STATUS=FAILED_PREPARE
            break
        fi

        eval "$prep_output"    # sets FP8_FLAG, CHUNK_END, RESUME_FROM, CHECKPOINT_PATH

        # On first chunk, honour any resume offset
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
                if [ $exit_code -ne 0 ]; then LOCUS_STATUS=FAILED_PREPARE; break; fi
                eval "$prep_output"
            fi
        fi

        output_dir=$(dirname "${OUT_FILE}")
        ref_fasta_path=${output_dir}/temp_ref.${MODEL_SIZE}.${window_size}bp.fa
        var_fasta_path=${output_dir}/temp_var.${MODEL_SIZE}.${window_size}bp.fa
        predict_ref_dir=${output_dir}/reference_predictions.${MODEL_SIZE}.${window_size}bp
        predict_var_dir=${output_dir}/variant_predictions.${MODEL_SIZE}.${window_size}bp

        # REF inference
        REF_START=$(date +%s)
        mpiexec -n ${TOTAL_GPUS} --ppn ${TOTAL_GPUS} \
            apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
            predict_evo2 \
                --fasta "${ref_fasta_path}" \
                --ckpt-dir "${CHECKPOINT_PATH}" \
                --output-dir "${predict_ref_dir}" \
                --model-size ${MODEL_SIZE} \
                --tensor-parallel-size ${tp_size} \
                --pipeline-model-parallel-size 1 \
                --context-parallel-size ${cp_size} \
                --output-log-prob-seqs ${FP8_FLAG}
        if [ $? -ne 0 ]; then LOCUS_STATUS=FAILED_REF; break; fi
        REF_END=$(date +%s)
        REF_SEC=$((REF_END - REF_START))

        # VAR inference
        VAR_START=$(date +%s)
        mpiexec -n ${TOTAL_GPUS} --ppn ${TOTAL_GPUS} \
            apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
            predict_evo2 \
                --fasta "${var_fasta_path}" \
                --ckpt-dir "${CHECKPOINT_PATH}" \
                --output-dir "${predict_var_dir}" \
                --model-size ${MODEL_SIZE} \
                --tensor-parallel-size ${tp_size} \
                --pipeline-model-parallel-size 1 \
                --context-parallel-size ${cp_size} \
                --output-log-prob-seqs ${FP8_FLAG}
        if [ $? -ne 0 ]; then LOCUS_STATUS=FAILED_VAR; break; fi
        VAR_END=$(date +%s)
        VAR_SEC=$((VAR_END - VAR_START))

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
            --chunk-size ${chunk_size}
        if [ $? -ne 0 ]; then LOCUS_STATUS=FAILED_PROCESS; break; fi

        chunk_start=${CHUNK_END}
    done

    LOCUS_END_TIME=$(date +%s)
    TOTAL_SEC=$((LOCUS_END_TIME - LOCUS_START_TIME))

    N_SCORED=0
    [ -f "${OUT_FILE}" ] && N_SCORED=$(( $(wc -l < "${OUT_FILE}") - 1 ))

    # Write timing row atomically
    (
      flock -x 200
      if [ ! -f "${TIMING_LOG}" ]; then
          echo "locus_id,n_variants,chunk_size,tp,context_bp,ref_inference_sec,var_inference_sec,total_sec,gpu_model,status" > "${TIMING_LOG}"
      fi
      echo "${LOCUS_ID},${N_SCORED},${chunk_size},${tp_size},${window_size},${REF_SEC},${VAR_SEC},${TOTAL_SEC},${GPU_MODEL},${LOCUS_STATUS}" >> "${TIMING_LOG}"
    ) 200>"${LOCK_FILE}"

    echo "Done: ${LOCUS_ID} — ${N_SCORED} variants scored in ${TOTAL_SEC}s (${LOCUS_STATUS})" >&2

done

echo "=== Batch complete ===" >&2
