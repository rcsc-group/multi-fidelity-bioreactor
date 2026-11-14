#!/bin/bash
#SBATCH --job-name=feature-bundle-mpi
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --exclude=node2356,node2410
#SBATCH --output=feature-bundle-mpi-%j.out
# The MPI case of the feature bundle (4 ranks, one node). Submit from the repository root:
#   sbatch tools/bundle/bundle_mpi.sbatch SRC_DIR WORK_DIR
# Writes WORK_DIR/mpi4.json next to the serial job's bundle.json. Not on batch2 (no MPI there).
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load openmpi
cd "$SLURM_SUBMIT_DIR"
uv run --project . python tools/bundle/bundle.py mpi "$@"
