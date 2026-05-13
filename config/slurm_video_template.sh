#!/bin/bash
# SLURM template for video runs (VIDEOS=1 binary).
# Produces vorticity3.mp4, volume_fraction3.mp4, oxygen3.mp4, tracer*.mp4
# in the run directory alongside the standard output files.
# Submit via scripts/submit_video_run.py.
#
#SBATCH --job-name=bioreactor-video
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elvis_alexander_aguero_vera@brown.edu

set -euo pipefail

# ffmpeg is required by Basilisk's save("*.mp4") to encode video frames.
# imagemagick is required for save("*.gif").
# Do NOT load the 'basilisk' module — the cluster's qcc is broken; we use the
# scratch-built binary at build/BioReactor-video which embeds TinyGL (fb_tiny)
# and has no runtime dependency on libGL or libGLX.
module load ffmpeg

if [ -z "${PARAMS:-}" ]; then
    echo "ERROR: PARAMS env var not set. Submit with: sbatch --export=PARAMS=<path> $0" >&2
    exit 1
fi

RUN_DIR="$(dirname "$PARAMS")"
PROJECT_ROOT="$(dirname "$(dirname "$RUN_DIR")")"

echo "Project root : $PROJECT_ROOT"
echo "Run dir      : $RUN_DIR"
echo "params.json  : $PARAMS"

mkdir -p "$RUN_DIR" "$PROJECT_ROOT/logs"

cd "$RUN_DIR"
# DISPLAY must be unset — fb_tiny renders to a memory buffer (no X server needed).
unset DISPLAY
OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" \
    "$PROJECT_ROOT/build/BioReactor-video" params.json

"$PROJECT_ROOT/.venv/bin/python" "$PROJECT_ROOT/scripts/postprocess.py" "$RUN_DIR"

echo "Done. Videos and results written to $RUN_DIR"
