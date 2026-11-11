#!/bin/bash -l
#SBATCH -J nest2d-submit
#SBATCH -p batch
#SBATCH -n 1 -c 1
#SBATCH --mem=4G
#SBATCH -t 00:20:00
#SBATCH --exclude=node2356
#SBATCH --export=NONE
#SBATCH -o /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/logs/nest2d_submit_%j.out
# Submits nest2d runs from a CLEAN login environment (2026-11-11: sbatch from the interactive shell passed its
# module list to the jobs; Lmod failed at 'module load openmpi' and 35 L8 jobs died in < 10 s). Same pattern as
# autopilot.sbatch.   sbatch --export=NONE scripts/bo2d/nest2d_submit.sbatch l8|l9
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
NEST_N="${2:-999}" uv run python scripts/bo2d/nest2d.py submit "$1" --go
