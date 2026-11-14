#!/bin/bash
#SBATCH --job-name=c2d-fit-next
#SBATCH --partition=batch
#SBATCH --account=default
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --exclude=node2356
#SBATCH --output=logs/c2d_fit_next_%j.out
#SBATCH --error=logs/c2d_fit_next_%j.err
set -euo pipefail
cd ${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}
uv run python scripts/campaign2d/campaign.py fit
uv run python scripts/campaign2d/campaign.py next --k 4
