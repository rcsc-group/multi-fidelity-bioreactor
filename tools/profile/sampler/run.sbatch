#!/bin/bash
#SBATCH --job-name=mfb-sampler
#SBATCH --partition=batch
#SBATCH --mem-per-cpu=4G
#SBATCH --exclude=node2356,node2410,node2411
#SBATCH --output=mfb-sampler-%j.out
# One sampler-check arm: build SRC_DIR (production mpi-lean) and run the L6 snapshot case into OUT_DIR.
#   sbatch --ntasks=N --time=.. tools/profile/sampler/run.sbatch SRC_DIR OUT_DIR
set -euo pipefail
SRC=$1; RD=$2
cd "$SLURM_SUBMIT_DIR"
module load openmpi 2>/dev/null || true
export OMP_NUM_THREADS=1
BAS=/oscar/data/dharri15/eaguerov/basilisk/src
rm -rf "$RD"; mkdir -p "$RD/build"
cp -r "$SRC"/. "$RD/build/"
( cd "$RD/build" && CC99='mpicc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1' \
  "$BAS/qcc" -O2 -w -Wall -autolink -D_MPI=1 -DEXTRA_TRACERS=0 BioReactor.c -o BioReactor -L"$BAS/gl" -lglutils -lfb_tiny -lm )
uv run --project . python tools/profile/sampler/params.py "$RD/params.json"
cd "$RD"
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU SLURM_TRES_PER_TASK SLURM_CPUS_PER_TASK
/usr/bin/time -v srun --mpi=pmix -n "$SLURM_NTASKS" build/BioReactor params.json > stdout.txt 2> stderr.txt
echo "done $(date)"
