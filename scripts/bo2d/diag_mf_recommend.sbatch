#!/bin/bash
#SBATCH -J diag-mfrec
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH --exclude=node2356
#SBATCH --array=0-15
#SBATCH -o logs/diag_mfrec_%A_%a.out
# 2 objectives x 8 replicates; one replicate per task.
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
OBJ=$([ $SLURM_ARRAY_TASK_ID -lt 8 ] && echo kla || echo mix); REP=$((SLURM_ARRAY_TASK_ID % 8))
uv run python scripts/bo2d/diag_mf_recommend.py $OBJ $REP
