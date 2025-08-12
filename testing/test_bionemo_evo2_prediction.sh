bionemo_container="/vast/projects/anuragv/cohort/mconery/bionemo/bionemo-framework_2.6.3.sif"
module load cuda/12.8.1 cudnn/8.9.7.29-12 apptainer/1.4.1
apptainer exec --nv --bind /vast/projects/anuragv/cohort/mconery:/vast/projects/anuragv/cohort/mconery $bionemo_container predict_evo2 \
  --fasta /vast/projects/anuragv/cohort/mconery/sequence_1kb.fasta \
  --ckpt-dir /vast/projects/anuragv/cohort/mconery/bionemo/nemo2_evo2_40b_1m \
  --output-dir /vast/projects/anuragv/cohort/mconery/1kb_test.txt \
  --model-size 40b_arc_longcontext \
  --tensor-parallel-size 1 \
  --pipeline-model-parallel-size 8 \
  --context-parallel-size 1 \
  --output-log-prob-seqs