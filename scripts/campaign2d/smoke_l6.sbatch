#!/bin/bash
# Campaign 2D smoke: one L6 three-cycle run, 4 MPI ranks, one node.
# Partition and account are given on the sbatch command line (batch/default or batch2/batch2).
# Usage: sbatch -p <part> -A <acct> --exclude=node2356 scripts/campaign2d/smoke_l6.sbatch <OUTDIR> <BINARY> <PARAMS>
#SBATCH --job-name=c2d-smoke-l6
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/c2d_smoke_%j.out
#SBATCH --error=logs/c2d_smoke_%j.err
set -euo pipefail
OUT="$1"; BIN="$2"; PARAMS="$3"
module load openmpi
mkdir -p "$OUT"
cp "$PARAMS" "$OUT/params.json"
cd "$OUT"
echo "partition=${SLURM_JOB_PARTITION:-?} account=${SLURM_JOB_ACCOUNT:-?} binary=$BIN host=$(hostname)"
srun --mpi=pmix -n 4 "$BIN" "$OUT/params.json"
