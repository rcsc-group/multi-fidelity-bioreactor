#!/bin/bash
# Golden-run gate on a compute node: build + run all cases from one source tree.
#   sbatch tools/golden/golden.sbatch SRC_DIR WORK_DIR [run|compile]
#SBATCH -J golden
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 03:00:00
#SBATCH --exclude=node2356,node2410
#SBATCH -o %x-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
uv run --project . python tools/golden/golden.py "${3:-run}" "$1" "$2"
