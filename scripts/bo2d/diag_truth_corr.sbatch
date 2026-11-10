#!/bin/bash
#SBATCH -J diag-tcorr
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --exclude=node2356
#SBATCH -o logs/diag_truth_corr_%j.out
# Truth fairness check.
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
uv run python scripts/bo2d/diag_truth_corr.py
