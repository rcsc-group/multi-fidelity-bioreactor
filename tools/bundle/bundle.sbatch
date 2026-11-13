#!/bin/bash
#SBATCH --job-name=feature-bundle
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --exclude=node2356,node2410
#SBATCH --output=feature-bundle-%j.out
# Feature bundle on a compute node. Submit from the repository root:
#   sbatch tools/bundle/bundle.sbatch SRC_DIR WORK_DIR [--only CASE ...]
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
uv run --project . python tools/bundle/bundle.py run "$@"
