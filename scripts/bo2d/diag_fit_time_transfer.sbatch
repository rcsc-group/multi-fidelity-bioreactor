#!/bin/bash
#SBATCH -J diag-fitt
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 00:40:00
#SBATCH --exclude=node2356
#SBATCH -o logs/diag_fit_time_%j.out
# fit timing and transfer-function fairness.
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
uv run python scripts/bo2d/diag_fit_time_transfer.py
