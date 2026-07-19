#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=01:00:00
#PBS -l filesystems=grand:home
#PBS -q preemptable
#PBS -A GeomicVar
#PBS -N pull_bionemo
#PBS -o /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/logs/pull_bionemo.log
#PBS -e /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/logs/pull_bionemo.err

ml use /soft/modulefiles
ml spack-pe-base
ml apptainer

SIF_PATH=/lus/grand/projects/GeomicVar/mconery/tools/bionemo-nightly.sif

mkdir -p /lus/grand/projects/GeomicVar/mconery/tools/apptainer_tmp
mkdir -p /lus/grand/projects/GeomicVar/mconery/tools/apptainer_cache
export APPTAINER_TMPDIR=/lus/grand/projects/GeomicVar/mconery/tools/apptainer_tmp
export APPTAINER_CACHEDIR=/lus/grand/projects/GeomicVar/mconery/tools/apptainer_cache

# Source NGC credentials
NGC_CREDS=/grand/GeomicVar/mconery/tools/ngc_credentials.sh
if [ -f "${NGC_CREDS}" ]; then
    source "${NGC_CREDS}"
    export APPTAINER_DOCKER_USERNAME='$oauthtoken'
    export APPTAINER_DOCKER_PASSWORD="${NGC_CLI_API_KEY}"
fi

echo "Pulling BioNeMo container to ${SIF_PATH} ..."
apptainer pull "${SIF_PATH}" docker://nvcr.io/nvidia/clara/bionemo-framework:nightly

echo "Done. SIF size: $(du -sh ${SIF_PATH} 2>/dev/null || echo 'not found')"
