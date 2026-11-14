#!/bin/bash
#SBATCH -J bo-theta
#SBATCH -p batch
#SBATCH -n 1 -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --exclude=node2356
#SBATCH -o /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/logs/bo_theta_%j.out
cd ${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}
uv run python scripts/bo2d/diag_theta_basis.py
