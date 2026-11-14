#!/bin/bash
#SBATCH -J diag-ctau
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --exclude=node2356
#SBATCH -o logs/diag_chain_tau_%j.out
# chain tau bias check.
cd ${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
uv run python scripts/bo2d/diag_chain_tau.py
