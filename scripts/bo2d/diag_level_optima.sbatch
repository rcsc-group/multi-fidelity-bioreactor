#!/bin/bash
#SBATCH -J diag-optima
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --exclude=node2356
#SBATCH -o logs/diag_level_optima_%j.out
# Go/no-go for "search cheap, confirm at L9" (OOM-killed at 16G in interactive 8746011, 2026-11-09 23:36).
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
uv run python scripts/bo2d/diag_level_optima.py
