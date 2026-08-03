#!/bin/bash

#SBATCH -o run_evo2_chunk_sweep_%x_%j.log
#SBATCH -p dgx-b200
#SBATCH --gpus=8
#SBATCH -t 8:00:00
#SBATCH -N 1

###############################################################################
# Chunk-size / GPU-memory saturation sweep for the Evo2 T2D fine-mapping
# priors pipeline -- runs on PARCC Betty.
#
# Production is fixed at ctx=8192bp, TP=8, CP=1, model=7b_arc_longcontext,
# chunk_size=100. This job checks whether chunk_size=100 has saturated GPU
# memory by walking chunk_size up through a candidate list (default: 100,
# 200, 400, 800, 1600, 3200; override via CHUNK_SIZES), running
# CHUNKS_PER_PHASE (default 2) real chunks of the LIVE production variant
# list at each size before advancing.
#
# Unlike a throwaway timing test, every chunk scored here is written into
# the real production output file (merged_evo2_scores.csv) via the
# unmodified score_locus_variants.py -- so calibration work is never wasted
# and never needs to be rerun. Sweep-only metrics (duration, peak GPU
# memory, status) are appended separately to chunk_size_sweep_log.csv.
#
# If any chunk_size OOMs (or otherwise fails), the sweep stops immediately
# -- no auto-fallback -- and the OOM is flagged in the log. The remaining
# variants are left unscored so you can review results and resume full
# production scoring manually with run_evo2_worker.sh at the last
# chunk_size that succeeded.
#
# CHUNK_SIZES and CHUNKS_PER_PHASE are normally passed via
# `sbatch --export=...` from launch_evo2_chunk_sweep.sh.
###############################################################################
: "${CHUNK_SIZES:=100 200 400 800 1600 3200}"
: "${CHUNKS_PER_PHASE:=2}"
tp_size=8
cp_size=1
window_size=8192
MODEL_SIZE=7b_arc_longcontext

###############################################################################
# Paths (match run_evo2_worker.sh / launch_evo2_jobs.sh conventions)
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
SWEEP_LOG=${BASE}/evo2_scoring/results/chunk_size_sweep_log.csv

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
    exit 1
fi
if [ ! -f "${FASTA}" ]; then
    echo "ERROR: FASTA not found: ${FASTA}" >&2
    exit 1
fi

###############################################################################
# Helpers
###############################################################################

# Number of variants already scored in the production output (0 if absent).
# Always fed back into --chunk-start explicitly so score_locus_variants.py's
# own resume_from = (n_done // chunk_size) * chunk_size recompute -- which
# would truncate the file if chunk_size differs from whatever wrote the last
# rows -- never runs. Explicit --chunk-start bypasses that path entirely.
count_done() {
    if [ -f "${OUT_FILE}" ]; then
        n=$(($(wc -l < "${OUT_FILE}") - 1))
        [ "$n" -lt 0 ] && n=0
        echo "$n"
    else
        echo 0
    fi
}

# Background GPU memory sampler: appends one memory.used (MiB) reading per
# GPU every 2s to $1 until killed.
#
# Must background directly in the caller's shell (not via `$(start_mem_monitor
# ...)`) -- wrapping this in a command substitution runs it in a subshell that
# exits as soon as it echoes the PID back, orphaning the backgrounded loop.
# Since bash never re-execs a background job, `ps` then shows the orphaned
# loop's cmd as this script's own path, indistinguishable from a second copy
# of the job script -- this is what the 2026-07-19 hang report misread as a
# duplicate/orphaned batch-script launch. Setting a global PID instead avoids
# the extra subshell entirely.
start_mem_monitor() {
    local outfile=$1
    : > "${outfile}"
    ( while true; do
          nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits >> "${outfile}" 2>/dev/null
          sleep 2
      done ) &
    MEM_MONITOR_PID=$!
}

stop_mem_monitor() {
    kill "${MEM_MONITOR_PID}" 2>/dev/null
    wait "${MEM_MONITOR_PID}" 2>/dev/null
}

peak_mem_mib() {
    local infile=$1
    if [ -s "${infile}" ]; then
        sort -n "${infile}" | tail -1
    fi
}

# Grep a captured srun log for common CUDA OOM signatures.
looks_like_oom() {
    grep -qiE "out of memory|OutOfMemoryError|CUBLAS_STATUS_ALLOC_FAILED|CUDA error" "$1" 2>/dev/null
}

append_sweep_row() {
    local chunk=$1 start_idx=$2 end_idx=$3 n_chunks=$4 t0=$5 t1=$6 dur=$7 peak=$8 status=$9 reason=${10}
    local header="model,context_bp,tp_size,cp_size,chunk_size,phase_start_idx,phase_end_idx,n_chunks_run,start_time,end_time,duration_seconds,peak_mem_mib,status,fail_reason"
    local row="${MODEL_SIZE},${window_size},${tp_size},${cp_size},${chunk},${start_idx},${end_idx},${n_chunks},${t0},${t1},${dur},${peak},${status},${reason}"
    (
        flock -x 200
        [ ! -f "${SWEEP_LOG}" ] && echo "${header}" >> "${SWEEP_LOG}"
        echo "${row}" >> "${SWEEP_LOG}"
    ) 200>"${SWEEP_LOG}.lock"
}

###############################################################################
# Sweep
###############################################################################
echo "=== Chunk-size sweep starting. Candidates: [${CHUNK_SIZES}] (${CHUNKS_PER_PHASE} chunks/phase) ===" >&2
echo "=== Scoring against live production list: ${VARIANT_FILE} -> ${OUT_FILE} ===" >&2

for CHUNK in ${CHUNK_SIZES}; do
    phase_start=$(count_done)
    phase_status=SUCCESS
    phase_reason=""
    n_chunks_run=0
    chunk_start=${phase_start}

    mem_log=${BASE}/evo2_scoring/results/gpu_mem.${MODEL_SIZE}.tp${tp_size}_cp${cp_size}.chunk${CHUNK}.${window_size}bp.log
    phase_t0=$(date +%s)

    for i in $(seq 1 ${CHUNKS_PER_PHASE}); do
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
            --chunk-size ${CHUNK})
        exit_code=$?
        if [ ${exit_code} -ne 0 ]; then
            echo "=== No variants remain to score at chunk_size=${CHUNK}; stopping sweep (production list exhausted). ===" >&2
            phase_status=DONE
            break 2
        fi
        eval "$prep_output"   # sets FP8_FLAG, CHUNK_END, RESUME_FROM, CHECKPOINT_PATH

        output_dir=$(dirname "${OUT_FILE}")
        ref_fasta_path=${output_dir}/temp_ref.${MODEL_SIZE}.${window_size}bp.fa
        var_fasta_path=${output_dir}/temp_var.${MODEL_SIZE}.${window_size}bp.fa
        predict_ref_dir=${output_dir}/reference_predictions.${MODEL_SIZE}.${window_size}bp
        predict_var_dir=${output_dir}/variant_predictions.${MODEL_SIZE}.${window_size}bp

        start_mem_monitor "${mem_log}"

        ref_log=$(mktemp)
        srun --ntasks=${TOTAL_GPUS} --ntasks-per-node=${TOTAL_GPUS} \
            apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
            predict_evo2 --fasta "${ref_fasta_path}" \
                --ckpt-dir "${CHECKPOINT_PATH}" \
                --output-dir "${predict_ref_dir}" \
                --model-size ${MODEL_SIZE} \
                --tensor-parallel-size ${tp_size} \
                --pipeline-model-parallel-size 1 \
                --context-parallel-size ${cp_size} \
                --output-log-prob-seqs ${FP8_FLAG} > "${ref_log}" 2>&1
        ref_exit=$?
        cat "${ref_log}" >&2

        if [ ${ref_exit} -ne 0 ]; then
            stop_mem_monitor
            if looks_like_oom "${ref_log}"; then
                phase_status=OOM
                phase_reason="OOM during REF inference at chunk_size=${CHUNK}, chunk_start=${chunk_start}"
            else
                phase_status=FAILED
                phase_reason="REF srun failed (exit ${ref_exit}) at chunk_size=${CHUNK}, chunk_start=${chunk_start}"
            fi
            rm -f "${ref_log}"
            break
        fi
        rm -f "${ref_log}"

        var_log=$(mktemp)
        srun --ntasks=${TOTAL_GPUS} --ntasks-per-node=${TOTAL_GPUS} \
            apptainer exec --nv --bind ${BIND_PATH} ${SIF_PATH} \
            predict_evo2 --fasta "${var_fasta_path}" \
                --ckpt-dir "${CHECKPOINT_PATH}" \
                --output-dir "${predict_var_dir}" \
                --model-size ${MODEL_SIZE} \
                --tensor-parallel-size ${tp_size} \
                --pipeline-model-parallel-size 1 \
                --context-parallel-size ${cp_size} \
                --output-log-prob-seqs ${FP8_FLAG} > "${var_log}" 2>&1
        var_exit=$?
        cat "${var_log}" >&2
        stop_mem_monitor

        if [ ${var_exit} -ne 0 ]; then
            if looks_like_oom "${var_log}"; then
                phase_status=OOM
                phase_reason="OOM during ALT inference at chunk_size=${CHUNK}, chunk_start=${chunk_start}"
            else
                phase_status=FAILED
                phase_reason="ALT srun failed (exit ${var_exit}) at chunk_size=${CHUNK}, chunk_start=${chunk_start}"
            fi
            rm -f "${var_log}"
            break
        fi
        rm -f "${var_log}"

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
            --chunk-size ${CHUNK}
        process_exit=$?
        if [ ${process_exit} -ne 0 ]; then
            phase_status=FAILED
            phase_reason="process step failed (exit ${process_exit}) at chunk_size=${CHUNK}, chunk_start=${chunk_start}"
            break
        fi

        n_chunks_run=$((n_chunks_run + 1))
        chunk_start=${CHUNK_END}
    done

    phase_t1=$(date +%s)
    phase_end=$(count_done)
    phase_dur=$((phase_t1 - phase_t0))
    peak=$(peak_mem_mib "${mem_log}")

    append_sweep_row "${CHUNK}" "${phase_start}" "${phase_end}" "${n_chunks_run}" \
        "${phase_t0}" "${phase_t1}" "${phase_dur}" "${peak}" "${phase_status}" "${phase_reason}"

    echo "=== chunk_size=${CHUNK}: ${n_chunks_run} chunk(s) run, variants ${phase_start}->${phase_end}, ${phase_dur}s, peak_mem=${peak} MiB, status=${phase_status} ===" >&2

    if [ "${phase_status}" != "SUCCESS" ]; then
        if [ "${phase_status}" = "OOM" ]; then
            echo "=== OOM at chunk_size=${CHUNK}. Stopping sweep -- do NOT run production at this chunk_size or higher. Review ${SWEEP_LOG}, then resume production with run_evo2_worker.sh at the last chunk_size that succeeded. ===" >&2
        elif [ "${phase_status}" = "FAILED" ]; then
            echo "=== Non-OOM failure at chunk_size=${CHUNK} (${phase_reason}). Stopping sweep for review. ===" >&2
        fi
        exit 1
    fi
done

echo "=== Sweep complete. Per-chunk-size duration/memory: ${SWEEP_LOG} ===" >&2
echo "=== ${OUT_FILE} now has $(count_done) real variants scored -- resume production with run_evo2_worker.sh (set chunk_size to your chosen value) to score the remainder. ===" >&2
