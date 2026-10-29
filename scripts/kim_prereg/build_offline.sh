#!/bin/bash
# Build the offline test variants (diary 2026-10-29) in the staged offline dir of pre_A_N256.
# Usage: bash /oscar/scratch/eaguerov/kim_pipeline/prereg_src/build_offline.sh <build tag> <N>
set -eu
module load openmpi >/dev/null 2>&1 || true
export BASILISK=/oscar/data/dharri15/eaguerov/basilisk/src
export PATH="$BASILISK:$PATH"
BUILD=${1:?build tag, e.g. pre_A_N256}; NG=${2:?grid N}
cd /oscar/scratch/eaguerov/kim_pipeline/prereg/oscar/builds/$BUILD/offline
# qcc does not find the absolute include from this dir: use a local copy of the header
cp -f /oscar/scratch/eaguerov/kim_pipeline/prereg_src/prereg_MK.h .
sed -i 's|^#include "/oscar/scratch/eaguerov/kim_pipeline/prereg_src/prereg_MK.h"|#include "prereg_MK.h"|' TracerTest_test.c
rm -rf .qcc*
COMMON="-DNN_GRID=$NG -DT_CYCLE=0.6074 -DT_END=2.0 -DMOVIES=0 -DFIGURES=0 -DOUT_FILES=0 -DOUT_INTERFACE=0 -DPREREG_DIAG=1"
build () {   # name, extra flags
  CC99="mpicc -std=gnu99 -D_GNU_SOURCE" qcc -D_MPI=1 -O2 -w $COMMON $2 \
      TracerTest_test.c -o "$1" -L"$BASILISK"/gl -lglutils -lfb_tiny -lm
  echo "built $1 ($2)"
}
build off_m1   "-DPREREG_MODE=1 -DPREREG_FRESET=1"
build off_m1f0 "-DPREREG_MODE=1 -DPREREG_FRESET=0"
build off_m2   "-DPREREG_MODE=2 -DPREREG_FRESET=1"
build off_m3   "-DPREREG_MODE=3 -DPREREG_FRESET=1"
build off_m3f0 "-DPREREG_MODE=3 -DPREREG_FRESET=0"
build off_m4   "-DPREREG_MODE=4 -DPREREG_FRESET=1"
