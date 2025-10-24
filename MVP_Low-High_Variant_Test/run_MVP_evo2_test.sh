#!/bin/bash
#SBATCH -J run_MVP_test
#SBATCH -o run_MVP_evo2_test.log
#SBATCH --gpus 8
#SBATCH -t 6:00:00

################################################################################################################
########################################## Define key model parameters #########################################
################################################################################################################
chunk_size=100
window_size=50000
MODEL_SIZE="7b"
tp_size=1
cp_size=1

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
out_file=/vast/projects/anuragv/cohort/mconery/mvp_variant_test/MVP_variant_scores."$MODEL_SIZE"_model."$window_size"bp_context.csv
evotwo_script="/vast/home/m/mconery/Evo2-TopMed-Variant-Scoring/MVP_Low-High_Variant_Test/mvp_variants_test.py"

################################################################################################################
############################################# Call Evo2 Scoring Script #########################################
################################################################################################################
#Load Modules and set caches
module load cuda/12.8.1 cudnn/8.9.7.29-12 apptainer/1.4.1
export HF_HOME=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NEMO_CACHE_DIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export BIONEMO_CACHE_DIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NEMO_MODELS_CACHE=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
apptainer exec --nv --bind /vast/projects/anuragv/cohort/mconery:/vast/projects/anuragv/cohort/mconery /vast/projects/anuragv/cohort/mconery/bionemo/bionemo-nightly.sif python $evotwo_script --out $out_file --chunk-size $chunk_size --ctx $window_size --model $MODEL_SIZE --tensor-parallel-size $tp_size --context-parallel-size $cp_size 
