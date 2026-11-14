#!/bin/bash
#SBATCH --job-name=henry-restart-repro
#SBATCH --output=logs/henry_restart_%j.out
#SBATCH --error=logs/henry_restart_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=elvis_alexander_aguero_vera@brown.edu

# Minimal reproducer: canonical henry.h + EMBED + MPI + restart.
# Tests whether NaN appears in tracer diffusion on checkpoint restart.
# See mpi_henry_restart.c for experimental design.
#
# Submit from the project root:
#   mkdir -p logs
#   sbatch sandbox/mpi_henry_restart.sh

set -euo pipefail

SPACK_BASILISK="/oscar/rt/9.6/25/spack/x86_64_v3/basilisk-2023.11.11-sgeoea5j7kth4wyjndwfitpekmopej4k"
SCRATCH_QCC="$HOME/scratch/basilisk/src/qcc"

# ── Rebuild qcc with correct BASILISK path if not present ──────────────────
# CPP99="" forces qcc to use CC99 (mpicc) as preprocessor, which:
#   1. Forwards -D_MPI=1 so _MPI-guarded code (balance.h, tree-mpi.h) is processed
#   2. Adds openmpi include paths so mpi.h can be found
# The spack include.o has a dead LIBDIR baked in — must recompile from source.
# The interpreter/ subdir provides ast_check_dimensions / ast_run symbols.
if [ ! -f "$SCRATCH_QCC" ]; then
  echo "Rebuilding qcc with correct BASILISK/LIBDIR paths..."
  SCRATCH_SRC="$HOME/scratch/basilisk/src"
  mkdir -p "$SCRATCH_SRC/ast"
  BPATH="$SPACK_BASILISK"
  # include.o: must recompile — spack version has dead LIBDIR
  cc -O2 -DLIBDIR=\"$BPATH\" -I"$BPATH" \
     -c "$BPATH/include.c" -o "$SCRATCH_SRC/include.o"
  # AST object files
  for f in "$BPATH/ast/"*.c; do
    base=$(basename "$f" .c)
    cc -O2 -DBASILISK=\"$BPATH\" -I"$BPATH" -I"$BPATH/ast" \
       -c "$f" -o "$SCRATCH_SRC/ast/${base}.o"
  done
  # interpreter/ is a separate subdirectory providing ast_check_dimensions
  cc -O2 -DBASILISK=\"$BPATH\" -I"$BPATH" -I"$BPATH/ast" \
     -c "$BPATH/ast/interpreter/interpreter.c" -o "$SCRATCH_SRC/ast/interpreter.o"
  ar rcs "$SCRATCH_SRC/ast/libast_new.a" "$SCRATCH_SRC/ast/"*.o
  # Link qcc — CPP99="" so mpicc is used as preprocessor (required for _MPI builds)
  cc -O2 \
     -DLIBDIR=\"$BPATH\" \
     -DCC99="\"mpicc -std=c99 -D_GNU_SOURCE=1 -D_XOPEN_SOURCE=700\"" \
     -DCPP99="\"\"" \
     -DCADNACC="\"\"" \
     -DBASILISK="\"$BPATH\"" \
     "$BPATH/qcc.c" "$SCRATCH_SRC/include.o" "$BPATH/postproc.o" \
     -o "$SCRATCH_QCC" -L"$SCRATCH_SRC/ast" -last_new -lm
  echo "qcc rebuilt at $SCRATCH_QCC"
fi

module load openmpi

cd /tmp
WORK="/tmp/henry_restart_$$"
mkdir -p "$WORK"
cp "$SLURM_SUBMIT_DIR/sandbox/mpi_henry_restart.c" "$WORK/"
cd "$WORK"

echo "=== Build ==="
CC99="mpicc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1" \
  "$SCRATCH_QCC" -D_MPI=1 -O2 -w \
  mpi_henry_restart.c -o mpi_henry_restart -lm
echo "Binary built: $(ls -la mpi_henry_restart)"

echo ""
echo "=== Run 1: Fresh simulation (t=0 → $(grep T_FRESH mpi_henry_restart.c | head -1)) ==="
rm -f checkpoint.dump
srun --mpi=pmix --mem=0 -n 4 ./mpi_henry_restart 2>&1
echo "Fresh run exit: $?"

echo ""
echo "=== Run 2: Restart simulation from checkpoint ==="
if [ ! -f checkpoint.dump ]; then
  echo "ERROR: checkpoint.dump not found after fresh run — aborting"
  exit 1
fi
srun --mpi=pmix --mem=0 -n 4 ./mpi_henry_restart restart 2>&1
EXIT_CODE=$?
echo "Restart run exit: $EXIT_CODE"

echo ""
if [ $EXIT_CODE -eq 0 ]; then
  echo "========================================================="
  echo "RESULT: PASS — canonical henry.h + EMBED + MPI + restart"
  echo "  is safe. Bug was specific to henry_oxy2.h modifications."
  echo "========================================================="
else
  echo "========================================================="
  echo "RESULT: FAIL — canonical henry.h + EMBED + MPI + restart"
  echo "  shows NaN. This is a Basilisk upstream bug."
  echo "  Report at: https://basilisk.fr/sandbox/bugs/README"
  echo "========================================================="
fi

cp "$WORK/mpi_henry_restart.c" "$SLURM_SUBMIT_DIR/sandbox/"
echo "Work dir: $WORK (kept for inspection)"
