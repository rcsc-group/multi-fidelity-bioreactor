#!/bin/bash
# build_mpi.sh TAG flags...  : copy worktree src, MPI build on a compute node (default account)
TAG="$1"; shift
W=${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}/.claude/worktrees/agent-a0a362873e179658b
D=${MFB_SCRATCH:-/oscar/scratch/eaguerov}/replay_gate/build_$TAG
rm -rf "$D"; mkdir -p "$D"
cp "$W"/src/*.c "$W"/src/*.h "$D"/
BUILD_MPI=1 sbatch --wait --export=ALL ${MFB_SCRATCH:-/oscar/scratch/eaguerov}/replay_gate/sbuild2.sh "$D" ${MFB_DATA:-/oscar/data/dharri15/eaguerov}/bin/replay-gate/BioReactor-$TAG "$@" | tail -1
