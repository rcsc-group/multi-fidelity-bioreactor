#!/bin/bash
#SBATCH --job-name=mfb-profile
#SBATCH --partition=batch
#SBATCH --mem-per-cpu=4G
#SBATCH --exclude=node2356,node2410,node2411
#SBATCH --output=mfb-profile-%j.out
# Profile the production build (mpi-lean) with Basilisk's built-in profiler (-DTRACE=2) and the multigrid log
# (-DDIAGNOSTICS=1). Submit from the repository root with --ntasks=N --time=...:
#   [FLAGS="-D..."] [TAG=name] sbatch --ntasks=N --time=HH:MM:SS tools/profile/profile.sbatch LEVEL PERIODS WORK_DIR
set -euo pipefail
LEVEL=$1; PERIODS=$2; WORK=$3
cd "$SLURM_SUBMIT_DIR"
OMPI_FALLBACK=/oscar/rt/9.6/25/spack/x86_64_v3/openmpi-5.0.8-q6ygsavgbbzcuxp2sog73dhpke6ttgfd
if ! module load openmpi 2>/dev/null; then
    export PATH="$OMPI_FALLBACK/bin:$PATH"
    export LD_LIBRARY_PATH="$OMPI_FALLBACK/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
export OMP_NUM_THREADS=1
BAS=/oscar/data/dharri15/eaguerov/basilisk/src
RD="$WORK/l${LEVEL}_n${SLURM_NTASKS}${TAG:+_$TAG}"
mkdir -p "$RD/build"
cp -r src/. "$RD/build/"
( cd "$RD/build" && CC99='mpicc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1' \
  "$BAS/qcc" -O2 -w -Wall -autolink -D_MPI=1 -DEXTRA_TRACERS=0 ${FLAGS:--DTRACE=2 -DDIAGNOSTICS=1} \
  BioReactor.c -o BioReactor -L"$BAS/gl" -lglutils -lfb_tiny -lm )
uv run --project . python tools/profile/params.py "$LEVEL" "$PERIODS" "$RD/params.json" ${NMIX:-1}
cd "$RD"
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU SLURM_TRES_PER_TASK SLURM_CPUS_PER_TASK   # inherited from the submitting allocation; srun refuses the mix
/usr/bin/time -v srun --mpi=pmix -n "$SLURM_NTASKS" build/BioReactor params.json > stdout.txt 2> stderr.txt
echo "done $(date)"
