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

# Basilisk's save("*.mp4") pipes PPM frames through ppm2mp4 (a helper script
# in the scratch Basilisk tree) which in turn calls ffmpeg.  Both must be
# findable via which().  We add the scratch Basilisk bin dir to PATH for the
# helper scripts, and load the ffmpeg module for the encoder itself.
# Do NOT load the 'basilisk' module — the cluster's qcc is broken.
module load ffmpeg
export PATH="$HOME/scratch/basilisk/src:$PATH"

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
