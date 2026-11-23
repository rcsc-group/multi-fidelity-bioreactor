#!/bin/bash
#SBATCH --job-name=mfb-cmprobe
#SBATCH --partition=batch
#SBATCH --mem-per-cpu=4G
#SBATCH --exclude=node2356,node2410,node2411
#SBATCH --output=mfb-cmprobe-%j.out
# One cm-probe arm (diary 2026-11-23, pre-registration "cm restriction"). From the repository root:
#   sbatch --ntasks=N --time=.. tools/profile/cm_probe/run.sbatch ARM LEVEL PERIODS WORK_DIR
set -euo pipefail
ARM=$1; LEVEL=$2; PERIODS=$3; WORK=$4
cd "$SLURM_SUBMIT_DIR"
module load openmpi 2>/dev/null || true
export OMP_NUM_THREADS=1
BAS=/oscar/data/dharri15/eaguerov/basilisk/src
RD="$WORK/cmprobe_${ARM}_l${LEVEL}_n${SLURM_NTASKS}"
mkdir -p "$RD"
uv run --project . python tools/profile/cm_probe/prepare.py "$ARM" "$RD/build"
( cd "$RD/build" && CC99='mpicc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1' \
  "$BAS/qcc" -O2 -w -Wall -autolink -D_MPI=1 -DEXTRA_TRACERS=0 BioReactor.c -o BioReactor -L"$BAS/gl" -lglutils -lfb_tiny -lm )
uv run --project . python tools/profile/params.py "$LEVEL" "$PERIODS" "$RD/params.json"
cd "$RD"
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU SLURM_TRES_PER_TASK SLURM_CPUS_PER_TASK
srun --mpi=pmix -n "$SLURM_NTASKS" build/BioReactor params.json > stdout.txt 2> stderr.txt
echo "done $(date)"
