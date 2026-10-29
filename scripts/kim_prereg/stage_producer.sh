#!/bin/bash
# Re-patch Kim's scratch copy, build one producer case with his pipeline, build my offline variants in its
# staged offline dir, and submit the producer. Prints the producer job id.
# Usage: bash stage_producer.sh <case tag> <N> <ranks> <walltime>
set -eu
TAG=$1; NG=$2; NR=$3; WALL=$4
S=/oscar/scratch/eaguerov/kim_pipeline/prereg_src
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
uv run python $S/patch_kim.py 2>&1 | grep -v warning
cd /oscar/scratch/eaguerov/kim_pipeline/prereg/oscar
module load openmpi >/dev/null 2>&1 || true
rm -rf builds/$TAG ../../prereg_out/$TAG
bash ./submit_case.sh $S/cases/$TAG.sh --build-only > $S/build_$TAG.log 2>&1
bash $S/build_offline.sh $TAG $NG 2>&1 | grep -v "warning: stencils"
env -u SLURM_CPUS_PER_TASK -u SLURM_TRES_PER_TASK sbatch --parsable --exclude=node2356 -J $TAG \
  -o $S/$TAG-%j.out -e $S/$TAG-%j.err --export=ALL,RUN=$PWD/builds/$TAG -N 1 -n $NR --mem 16G -t $WALL \
  -p batch A2_producer.slurm
