#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --mem 32G
#SBATCH --exclude=node2356
# One offline test variant (diary 2026-10-29) on the record of one producer run.
# Usage (ranks = the producer's rank count):
#   env -u SLURM_CPUS_PER_TASK -u SLURM_TRES_PER_TASK sbatch -n <ranks> -t <walltime> -J <name> \
#     -o <log> --export=ALL,TAG=<run tag>,VAR=<off_m1|...> run_offline.sbatch
set -eu
set +u; module load openmpi >/dev/null 2>&1; set -u
P=/oscar/scratch/eaguerov/kim_pipeline/prereg_out/$TAG/producer
B=/oscar/scratch/eaguerov/kim_pipeline/prereg/oscar/builds/$TAG/offline/$VAR
D=/oscar/scratch/eaguerov/kim_pipeline/prereg_out/$TAG/$VAR
. /oscar/scratch/eaguerov/kim_pipeline/prereg/oscar/builds/$TAG/run_env.sh
rm -rf "$D"; mkdir -p "$D"/Data_all; cd "$D"
for f in "$P"/prereg_rec.*.bin; do ln -s "$f" .; done
echo "host $(hostname) job ${SLURM_JOB_ID:-} tag $TAG var $VAR args $ARGS"
/usr/bin/time -f "wall %e s" srun --mpi=pmix "$B" $ARGS
