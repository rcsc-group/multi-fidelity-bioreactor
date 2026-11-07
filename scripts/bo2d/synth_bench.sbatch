#!/bin/bash
#SBATCH -J bo-synth
#SBATCH -p batch
#SBATCH -n 1 -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --exclude=node2356
#SBATCH -o /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/logs/bo_synth_%j.out
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
uv run python scripts/bo2d/synth_bench.py "${1:-0.1}" "${2:-dev}"
