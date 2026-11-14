#!/bin/bash
#SBATCH -J diag-margin
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 00:30:00
#SBATCH --exclude=node2356
#SBATCH -o logs/diag_infeas_margin_%j.out
# infeasible margin check.
cd ${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
uv run python scripts/bo2d/diag_infeas_margin.py
