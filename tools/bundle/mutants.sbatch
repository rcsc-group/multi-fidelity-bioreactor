#!/bin/bash
#SBATCH --job-name=bundle-mutant
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --exclude=node2356,node2410
#SBATCH --output=bundle-mutant-%x-%j.out
# One planted defect on a compute node. Submit from the repository root:
#   sbatch -J <NAME> tools/bundle/mutants.sbatch NAME WORK_DIR
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
uv run --project . python tools/bundle/mutants.py run "$@"
