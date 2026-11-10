#!/bin/bash
#SBATCH -J bo-autopilot
#SBATCH -p batch
#SBATCH -n 1 -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH --exclude=node2356
#SBATCH -o /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/logs/bo_autopilot_%j.out
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
scripts/bo2d/autopilot_wed.sh
