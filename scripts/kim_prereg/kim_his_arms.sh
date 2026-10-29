#!/bin/bash
# Kim's own offline code (diag switch only) on snapshot subsets of one dense producer dump (diary amendment 2).
#   his_S    : every snapshot (NDUMP per cycle, 3 cycles), no wrap
#   his_H    : every STRIDE-th snapshot (100 per cycle), 3 cycles, no wrap
#   his_ship : every STRIDE-th snapshot of the first cycle, looped 3 times (his workflow)
# Usage: bash kim_his_arms.sh <producer tag> <ranks> <stride> <walltime>
set -eu
TAG=$1; NR=$2; STRIDE=$3; WALL=$4
OUT=/oscar/scratch/eaguerov/kim_pipeline/prereg_out/$TAG
BIN=$OUT/producer/Data_dump/Data_bin
SUB=$OUT/subsets
rm -rf "$SUB"
mkdir -p "$SUB"/all/Data_bin "$SUB"/strided/Data_bin "$SUB"/cycle1/Data_bin
mapfile -t SNAPS < <(find "$BIN" -maxdepth 1 -name 'Data_dump_*_tot.bin' | awk -F_ '{t=$(NF-1); print t, $0}' | sort -g | cut -d' ' -f2)
N=${#SNAPS[@]}
PER=$(( N / 3 ))                       # snapshots per cycle in the dense dump
for d in all strided cycle1; do ln -s "$BIN/flowgrid.bin" "$SUB/$d/Data_bin/flowgrid.bin"; done
for ((k = 0; k < N; k++)); do
  ln -s "${SNAPS[$k]}" "$SUB/all/Data_bin/"
  if (( k % STRIDE == 0 )); then
    ln -s "${SNAPS[$k]}" "$SUB/strided/Data_bin/"
    (( k < PER )) && ln -s "${SNAPS[$k]}" "$SUB/cycle1/Data_bin/"
  fi
done
echo "snapshots: all $N, strided $(ls "$SUB"/strided/Data_bin | grep -c tot), cycle1 $(ls "$SUB"/cycle1/Data_bin | grep -c tot)"

cd /oscar/scratch/eaguerov/kim_pipeline/prereg/oscar
module load openmpi >/dev/null 2>&1 || true
COMMONENV="OUTROOT=$OUT BASILISK=/oscar/data/dharri15/eaguerov/basilisk/src TPP=/oscar/scratch/eaguerov/kim_pipeline/TPP_headers"
TT="TT_IMPORT_METHOD=2 TT_PREREG_DIAG=1 TT_MOVIES=0 TT_FIGURES=0 TT_OUT_FILES=0 TT_OUT_INTERFACE=0"
for spec in "his_S:all:" "his_H:strided:" "his_ship:cycle1:3"; do
  arm=${spec%%:*}; rest=${spec#*:}; sub=${rest%%:*}; ncyc=${rest#*:}
  RUNB=/oscar/scratch/eaguerov/kim_pipeline/prereg/oscar/builds_offline/${TAG}_$arm
  rm -rf "$RUNB" "$OUT/$arm"
  if [ -n "$ncyc" ]; then
    env $COMMONENV $TT PRODUCER="$SUB/$sub" OFFTAG=$arm RUN=$RUNB NCYCLE=$ncyc bash B1_build_offline.sh > "$OUT/build_$arm.log" 2>&1
  else
    env $COMMONENV $TT PRODUCER="$SUB/$sub" OFFTAG=$arm RUN=$RUNB bash B1_build_offline.sh > "$OUT/build_$arm.log" 2>&1
  fi
  B=$RUNB
  echo "$arm built in $B: $(grep -E 'snapshots|t_ini|replays' "$OUT/build_$arm.log" | tr -s ' ' | tr '\n' ';')"
  env -u SLURM_CPUS_PER_TASK -u SLURM_TRES_PER_TASK sbatch --parsable -p batch -N 1 -n "$NR" --mem 32G -t "$WALL" \
    --exclude=node2356 -J "$arm" -o "$OUT/$arm-%j.out" -D "$OUT/$arm" \
    --wrap "module load openmpi; srun --mpi=pmix $B/offline_binary 0.25 7 20"
done
