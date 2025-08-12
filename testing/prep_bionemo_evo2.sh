########################## Download Bionemo Container #########################
#Load modules and pull container
module load cuda/12.8.1 cudnn/8.9.7.29-12 apptainer/1.4.1
NVIDIA_API_KEY="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"
export SINGULARITY_DOCKER_USERNAME='$oauthtoken'
export SINGULARITY_DOCKER_PASSWORD=$NVIDIA_API_KEY
cd /vast/projects/anuragv/cohort/mconery/bionemo
apptainer pull bionemo-framework_2.6.3.sif docker://nvcr.io/nvidia/clara/bionemo-framework:2.6.3

######################### Install Evo2 for Bionemo ############################
#Activate container shell and install evo2 inside it
apptainer shell --nv --bind /vast/projects/anuragv/cohort/mconery/bionemo:/vast/projects/anuragv/cohort/mconery/bionemo /vast/projects/anuragv/cohort/mconery/bionemo/bionemo-framework_2.6.3.sif
cd /vast/projects/anuragv/cohort/mconery/bionemo/bionemo-framework/sub-packages/bionemo-evo2
pip install -e .

######################### Download and Convert Savanna Checkpoint ###################
cd /vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export APPTAINER_CACHEDIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export MPLCONFIGDIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export HF_HOME=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NEMO_CACHE_DIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export BIONEMO_CACHE_DIR=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
export NEMO_MODELS_CACHE=/vast/projects/anuragv/cohort/mconery/bionemo/hf_cache
evo2_convert_to_nemo2 --model-path hf://arcinstitute/savanna_evo2_40b --model-size 40b_arc_longcontext --output-dir nemo2_evo2_40b_1m