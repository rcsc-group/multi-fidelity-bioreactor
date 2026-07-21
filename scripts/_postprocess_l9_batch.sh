#!/bin/bash
#SBATCH --job-name=pp-l9-batch
#SBATCH --account=mbessa-gcondo
#SBATCH --partition=mbessa-gcondo
#SBATCH --qos=mbessa-gcondo
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/oscar/data/dharri15/eaguerov/Github/BioReactor3D/dev/rocking-bioreactor-2d/logs/pp_l9_batch_%j.out
#SBATCH --error=/oscar/data/dharri15/eaguerov/Github/BioReactor3D/dev/rocking-bioreactor-2d/logs/pp_l9_batch_%j.err

cd /oscar/data/dharri15/eaguerov/Github/BioReactor3D/dev/rocking-bioreactor-2d
uv run python -m scripts._postprocess_l9_batch
