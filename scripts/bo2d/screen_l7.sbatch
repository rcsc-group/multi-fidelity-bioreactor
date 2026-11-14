#!/bin/bash
#SBATCH -J screen-l7
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 12:00:00
#SBATCH --exclude=node2356
#SBATCH -o logs/screen_l7_%j.out
# Throttled submitter for the dense L7 screen (screen_l7.py); <= 30 L7 jobs in flight.
cd ${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}
uv run python scripts/bo2d/screen_l7.py run --go
