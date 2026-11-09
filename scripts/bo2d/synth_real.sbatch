#!/bin/bash
#SBATCH -J bo-real
#SBATCH -p batch
#SBATCH -n 1 -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --exclude=node2356
#SBATCH -o /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/logs/bo_real_%j.out
# Synthetic test on a truth built from the REAL data (2026-11-09). $1 = kla|mix, $2 = seed, $3 = arm.
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
BO_OBJ=$1 BO_METHOD=mf BO_THETA=1 BO_FAIRSEED=1 uv run python scripts/bo2d/synth_bench.py 0.014 real "$2" "$3"
