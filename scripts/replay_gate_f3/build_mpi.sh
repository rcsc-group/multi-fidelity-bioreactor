#!/bin/bash
# build_mpi.sh TAG flags...  : copy worktree src, MPI build on a compute node (default account)
TAG="$1"; shift
W=/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/.claude/worktrees/agent-a0a362873e179658b
D=/oscar/scratch/eaguerov/replay_gate/build_$TAG
rm -rf "$D"; mkdir -p "$D"
cp "$W"/src/*.c "$W"/src/*.h "$D"/
BUILD_MPI=1 sbatch --wait --export=ALL /oscar/scratch/eaguerov/replay_gate/sbuild2.sh "$D" /oscar/data/dharri15/eaguerov/bin/replay-gate/BioReactor-$TAG "$@" | tail -1
