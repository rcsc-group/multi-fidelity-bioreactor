#!/bin/bash
#SBATCH -J bo-screen
#SBATCH -p batch
#SBATCH -n 1 -c 4
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --exclude=node2356
#SBATCH -o /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/logs/bo_screen_%j.out
# screen-then-L9 (synth_screen.py, pre-registered). $1 = kla|mix, $2 = seed; env SCREEN_N.
cd ${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}
uv run python scripts/bo2d/synth_screen.py "$1" "$2"
