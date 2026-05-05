#!/bin/bash
# SLURM job template for a single BioReactor simulation.
# Submitted by simulate.py with --export=PARAMS=<absolute path to params.json>
#
# Key SLURM directives — override on the sbatch command line with --time, --mem, etc.
#SBATCH --job-name=bioreactor
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elvis_alexander_aguero_vera@brown.edu

set -euo pipefail

# PARAMS must be an absolute path to params.json (set by --export on sbatch)
if [ -z "${PARAMS:-}" ]; then
    echo "ERROR: PARAMS env var not set. Submit with: sbatch --export=PARAMS=<path> $0" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
RUN_DIR="$(dirname "$PARAMS")"

echo "Project root : $PROJECT_ROOT"
echo "Run dir      : $RUN_DIR"
echo "params.json  : $PARAMS"

mkdir -p "$RUN_DIR" "$PROJECT_ROOT/logs"

# Run simulation (binary must already be compiled via 'make build')
cd "$RUN_DIR"
OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" \
    "$PROJECT_ROOT/build/BioReactor" params.json

# Extract kLa and write results.json
"$PROJECT_ROOT/.venv/bin/python" "$PROJECT_ROOT/scripts/postprocess.py" "$RUN_DIR"

echo "Done. Results written to $RUN_DIR/results.json"
