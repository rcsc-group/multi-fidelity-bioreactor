#!/bin/bash
# Kim's own offline code (diag switch only): EVERY snapshot of the first cycle of the dense dump, looped 3 times
# (flux rebuilt from u, no hold, looped). Same build/run path as kim_his_arms.sh.
# Usage: bash kim_his_loop.sh <producer tag> <ranks> <walltime>
set -eu
TAG=$1; NR=$2; WALL=$3
OUT=/oscar/scratch/eaguerov/kim_pipeline/prereg_out/$TAG
BIN=$OUT/producer/Data_dump/Data_bin
SUB=$OUT/subsets/cycle1all
rm -rf "$SUB"
mkdir -p "$SUB"/Data_bin
mapfile -t SNAPS < <(find "$BIN" -maxdepth 1 -name 'Data_dump_*_tot.bin' | awk -F_ '{t=$(NF-1); print t, $0}' | sort -g | cut -d' ' -f2)
N=${#SNAPS[@]}
PER=$(( N / 3 ))
ln -s "$BIN/flowgrid.bin" "$SUB/Data_bin/flowgrid.bin"
for ((k = 0; k < PER; k++)); do ln -s "${SNAPS[$k]}" "$SUB/Data_bin/"; done
echo "snapshots: cycle1all $(ls "$SUB"/Data_bin | grep -c tot)"

cd /oscar/scratch/eaguerov/kim_pipeline/prereg/oscar
module load openmpi >/dev/null 2>&1 || true
COMMONENV="OUTROOT=$OUT BASILISK=/oscar/data/dharri15/eaguerov/basilisk/src TPP=/oscar/scratch/eaguerov/kim_pipeline/TPP_headers"
TT="TT_IMPORT_METHOD=2 TT_PREREG_DIAG=1 TT_MOVIES=0 TT_FIGURES=0 TT_OUT_FILES=0 TT_OUT_INTERFACE=0"
arm=his_loop
RUNB=/oscar/scratch/eaguerov/kim_pipeline/prereg/oscar/builds_offline/${TAG}_$arm
rm -rf "$RUNB" "$OUT/$arm"
env $COMMONENV $TT PRODUCER="$SUB" OFFTAG=$arm RUN=$RUNB NCYCLE=3 bash B1_build_offline.sh > "$OUT/build_$arm.log" 2>&1
echo "$arm built in $RUNB: $(grep -E 'snapshots|t_ini|replays' "$OUT/build_$arm.log" | tr -s ' ' | tr '\n' ';')"
env -u SLURM_CPUS_PER_TASK -u SLURM_TRES_PER_TASK sbatch --parsable -p batch -N 1 -n "$NR" --mem 32G -t "$WALL" \
  --exclude=node2356 -J "$arm" -o "$OUT/$arm-%j.out" -D "$OUT/$arm" \
  --wrap "module load openmpi; srun --mpi=pmix $RUNB/offline_binary 0.25 7 20"
