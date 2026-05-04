#!/bin/bash
#SBATCH --job-name=bioreactor_{run_id}
#SBATCH --output=runs/{run_id}/slurm_%j.out
#SBATCH --error=runs/{run_id}/slurm_%j.err
#SBATCH --time={walltime}
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elvis_alexander_aguero_vera@brown.edu

module load basilisk

cd {project_root}

OMP_NUM_THREADS=4 ./build/BioReactor runs/{run_id}/params.json
