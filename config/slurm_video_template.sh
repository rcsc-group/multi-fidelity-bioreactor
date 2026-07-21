#!/bin/bash
# SLURM template for video runs (VIDEOS=1 binary).
# BioReactor-video only dumps raw binary frames to run_dir/frames/;
# render_videos.py (called below) is what actually renders them, producing
# volume_fraction.mp4 (body frame) and volume_fraction_lab.mp4 (lab frame)
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

# render_videos.py invokes ffmpeg directly to encode the rendered frames.
# Do NOT load the 'basilisk' module — the cluster's qcc is broken; use the
# persistent build instead (see docs_site/setup.md).
module load ffmpeg
export PATH="/oscar/data/dharri15/eaguerov/basilisk/src:$PATH"

if [ -z "${PARAMS:-}" ]; then
    echo "ERROR: PARAMS env var not set. Submit with: sbatch --export=PARAMS=<path> $0" >&2
    exit 1
fi

RUN_DIR="$(realpath "$(dirname "$PARAMS")")"
PROJECT_ROOT="$(realpath "$(dirname "$(dirname "$RUN_DIR")")")"

echo "Project root : $PROJECT_ROOT"
echo "Run dir      : $RUN_DIR"
echo "params.json  : $PARAMS"

mkdir -p "$RUN_DIR" "$PROJECT_ROOT/logs"

cd "$RUN_DIR"
# DISPLAY must be unset — fb_tiny renders to a memory buffer (no X server needed).
unset DISPLAY
if [ -n "${DUMP:-}" ]; then
    OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" \
        "$PROJECT_ROOT/build/BioReactor-video" params.json "$DUMP"
else
    OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" \
        "$PROJECT_ROOT/build/BioReactor-video" params.json
fi

"$PROJECT_ROOT/.venv/bin/python" "$PROJECT_ROOT/scripts/postprocess.py" "$RUN_DIR"
"$PROJECT_ROOT/.venv/bin/python" "$PROJECT_ROOT/scripts/render_videos.py" "$RUN_DIR"

echo "Done. Videos and results written to $RUN_DIR"
