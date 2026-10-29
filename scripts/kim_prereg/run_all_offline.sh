#!/bin/bash
# After a producer case finished: submit my record-machinery arms and Kim's own-code arms.
# Usage: bash run_all_offline.sh <case tag> <ranks> <walltime>
set -eu
TAG=$1; NR=$2; WALL=$3
S=/oscar/scratch/eaguerov/kim_pipeline/prereg_src
for v in off_m1 off_m1f0 off_m3 off_m3f0 off_m4; do
  env -u SLURM_CPUS_PER_TASK -u SLURM_TRES_PER_TASK sbatch --parsable -n $NR -t $WALL -J ${TAG}_$v \
    -o $S/${TAG}_$v-%j.out --export=ALL,TAG=$TAG,VAR=$v $S/run_offline.sbatch
done
bash $S/kim_his_arms.sh $TAG $NR 20 $WALL
