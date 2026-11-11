#!/bin/bash
#SBATCH -J bo-taumf
#SBATCH -p batch
#SBATCH -n 1 -c 4
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --exclude=node2356
#SBATCH -o /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/logs/bo_taumf_%j.out
# MF on the tau constraint only (synth_taumf.py, pre-registered in its docstring). $1 = kla|mix, $2 = seed.
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
uv run python scripts/bo2d/synth_taumf.py "$1" "$2"
