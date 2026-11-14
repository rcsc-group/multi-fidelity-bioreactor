#!/bin/bash
#SBATCH -J mf-abl-fair
#SBATCH -p batch
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH --exclude=node2356
#SBATCH --array=0-23
#SBATCH -o logs/mf_ablation_fair_%A_%a.out
# Fair ablation grid: qoi (3) x n9 {4,8} x on7 {0,1} x l8line {keep,drop} = 24 tasks
cd ${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}
Q=(kLa_1T_10 dtmix_0.25 tau_95); N=(4 8); O=(0 1); L=(keep drop); i=$SLURM_ARRAY_TASK_ID
uv run python scripts/bo2d/mf_ablation_fair.py ${Q[$((i / 8))]} ${N[$((i / 4 % 2))]} ${O[$((i / 2 % 2))]} ${L[$((i % 2))]}
