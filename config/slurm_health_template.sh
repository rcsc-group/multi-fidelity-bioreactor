#!/bin/bash
# SLURM template for health-check runs (DIAGNOSTICS=1 binary).
# Writes pressure_diag.dat in addition to standard outputs.
# Submit via scripts/submit_health_runs.py.
#
#SBATCH --job-name=bioreactor-health
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elvis_alexander_aguero_vera@brown.edu

set -euo pipefail

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
OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" \
    "$PROJECT_ROOT/build/BioReactor-health" params.json

"$PROJECT_ROOT/.venv/bin/python" "$PROJECT_ROOT/scripts/postprocess.py" "$RUN_DIR"

echo "Done. Results written to $RUN_DIR/results.json"
