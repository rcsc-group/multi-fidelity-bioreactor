#!/bin/bash
#SBATCH -J mf-ablation
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --exclude=node2356
#SBATCH --array=0-8
#SBATCH -o logs/mf_ablation_%A_%a.out
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
Q=(kLa_1T_10 dtmix_0.25 tau_95); N=(4 8 16)
uv run python scripts/bo2d/mf_ablation.py ${Q[$((SLURM_ARRAY_TASK_ID / 3))]} ${N[$((SLURM_ARRAY_TASK_ID % 3))]}
