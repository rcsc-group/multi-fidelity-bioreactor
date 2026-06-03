#!/bin/bash
# SLURM template for MPI-parallel BioReactor runs.
# Uses srun --mpi=pmix to launch across ntasks MPI ranks.
# Build the binary first: make build-mpi
#
#SBATCH --job-name=bioreactor-mpi
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elvis_alexander_aguero_vera@brown.edu

set -euo pipefail

module load openmpi
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
echo "MPI ranks    : ${SLURM_NTASKS:-16}"

mkdir -p "$RUN_DIR" "$PROJECT_ROOT/logs"

unset DISPLAY
cd "$RUN_DIR"
if [ -n "${DUMP:-}" ]; then
    srun --mpi=pmix \
        "$PROJECT_ROOT/build/BioReactor-mpi" params.json "$DUMP"
else
    srun --mpi=pmix \
        "$PROJECT_ROOT/build/BioReactor-mpi" params.json
fi

uv run python "$PROJECT_ROOT/scripts/postprocess.py" "$RUN_DIR"

echo "Done. Results written to $RUN_DIR/results.json"

# Self-submitting chain: submit next segment on completion
NEXT_RUN=$(python3 -c "
import json, sys
try:
    p = json.load(open(sys.argv[1]))
    print(p.get('next_run_id', ''))
except:
    print('')
" "$PARAMS" 2>/dev/null)

if [ -n "$NEXT_RUN" ]; then
    NEXT_PARAMS="$PROJECT_ROOT/runs/$NEXT_RUN/params.json"
    NEXT_DUMP="$RUN_DIR/checkpoint.dump"
    WALLTIME=$(python3 -c "
import json, sys
try:
    p = json.load(open(sys.argv[1]))
    print(p.get('_walltime', '04:00:00'))
except:
    print('04:00:00')
" "$NEXT_PARAMS" 2>/dev/null)
    MEM=$(python3 -c "
import json, sys
try:
    p = json.load(open(sys.argv[1]))
    print(p.get('_mem', '4G'))
except:
    print('4G')
" "$NEXT_PARAMS" 2>/dev/null)
    NTASKS=$(python3 -c "
import json, sys
try:
    p = json.load(open(sys.argv[1]))
    print(p.get('_ntasks', 16))
except:
    print(16)
" "$NEXT_PARAMS" 2>/dev/null)
    sbatch --no-requeue \
        --time="$WALLTIME" \
        --mem-per-cpu="$MEM" \
        --ntasks="$NTASKS" \
        --cpus-per-task=1 \
        --export="NONE,PARAMS=$NEXT_PARAMS,DUMP=$NEXT_DUMP" \
        "$PROJECT_ROOT/config/slurm_mpi_template.sh"
    echo "Submitted next segment: $NEXT_RUN"
fi
