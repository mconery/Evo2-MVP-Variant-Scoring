#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=00:30:00
#PBS -l filesystems=grand:home
#PBS -q preemptable
#PBS -A GeomicVar
#PBS -N decompress_fasta
#PBS -o /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/logs/decompress_fasta.log
#PBS -e /grand/GeomicVar/mconery/evo2_variant_scoring_mapping/logs/decompress_fasta.err

# Decompress the GRCh37 FASTA and index it for pyfaidx.
# Run once before Evo2 scoring. Takes ~5 minutes.

FASTA_GZ=/grand/GeomicVar/mconery/resources/assemblies/GRCh37.p13.genome.fa.gz
BASE=/grand/GeomicVar/mconery/evo2_variant_scoring_mapping
OUT_DIR=${BASE}/reference
OUT_FASTA=${OUT_DIR}/GRCh37.p13.genome.fa
CONDA_BASE=/grand/GeomicVar/mconery/tools/miniconda3

mkdir -p "${OUT_DIR}"

if [ -f "${OUT_FASTA}.fai" ]; then
    echo "FASTA already decompressed and indexed at ${OUT_FASTA}"
    exit 0
fi

echo "Decompressing ${FASTA_GZ} ..."
gunzip -c "${FASTA_GZ}" > "${OUT_FASTA}"

echo "Indexing with pyfaidx ..."
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate finemap
python -c "import pyfaidx; pyfaidx.Faidx('${OUT_FASTA}')"

echo "Done: ${OUT_FASTA}"
