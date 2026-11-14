#!/bin/bash
#SBATCH --job-name=c2d-diag-unc
#SBATCH --partition=batch
#SBATCH --account=default
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --exclude=node2356
#SBATCH --output=logs/c2d_diag_unc_%j.out
#SBATCH --error=logs/c2d_diag_unc_%j.err
set -euo pipefail
cd ${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}
uv run python scripts/campaign2d/diag_uncertainty.py
