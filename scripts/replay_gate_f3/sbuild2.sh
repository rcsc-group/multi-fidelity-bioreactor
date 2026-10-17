#!/bin/bash
# Usage: sbuild2.sh DIR OUT [qcc flags]  -- compile DIR/BioReactor.c to OUT on a compute node (default account).
# BUILD_MPI=1 in the environment gives the MPI build of the Makefile's build-mpi-lean kind (without EXTRA_TRACERS change).
#SBATCH --job-name=replay-build
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/oscar/scratch/eaguerov/replay_gate/build_%j.log
set -euo pipefail
D="$1"; OUT="$2"; shift 2
BASILISK=/oscar/data/dharri15/eaguerov/basilisk/src
export BASILISK
cd "$D"
if [ "${BUILD_MPI:-0}" = "1" ]; then
  module load openmpi
  CC99='mpicc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1' "$BASILISK/qcc" -O2 -w -fopenmp -Wall -autolink -D_MPI=1 "$@" BioReactor.c -o "$OUT" -L"$BASILISK/gl" -lglutils -lfb_tiny -lm
else
  "$BASILISK/qcc" -O2 -w -fopenmp -Wall -autolink "$@" BioReactor.c -o "$OUT" -L"$BASILISK/gl" -lglutils -lfb_tiny -lm
fi
echo built "$OUT"
