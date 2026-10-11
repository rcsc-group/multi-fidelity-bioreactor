#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -t 02:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -J lagr
#SBATCH -o logs/lagrangian/%x_%j.out
# usage: sbatch scripts/lagrangian_gpu.sbatch <script.py> [args...]
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
uv run python "$@"
