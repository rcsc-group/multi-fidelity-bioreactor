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

# Derive project root from PARAMS (always an absolute path to runs/<id>/params.json).
# Do NOT use BASH_SOURCE — SLURM stages the script under /var/spool/slurmd, so
# dirname(BASH_SOURCE) resolves to the SLURM daemon's directory, not the repo.
RUN_DIR="$(dirname "$PARAMS")"
PROJECT_ROOT="$(dirname "$(dirname "$RUN_DIR")")"

echo "Project root : $PROJECT_ROOT"
echo "Run dir      : $RUN_DIR"
echo "params.json  : $PARAMS"

mkdir -p "$RUN_DIR" "$PROJECT_ROOT/logs"

# Run simulation (binary must already be compiled via 'make build').
# DUMP, if set, is the absolute path to a checkpoint.dump from a previous
# segment; the binary will restore it and start from that flow field.
cd "$RUN_DIR"
if [ -n "${DUMP:-}" ]; then
    OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" \
        "$PROJECT_ROOT/build/BioReactor" params.json "$DUMP"
else
    OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" \
        "$PROJECT_ROOT/build/BioReactor" params.json
fi

# Extract kLa and write results.json
"$PROJECT_ROOT/.venv/bin/python" "$PROJECT_ROOT/scripts/postprocess.py" "$RUN_DIR"

echo "Done. Results written to $RUN_DIR/results.json"

# Self-submitting chain: if this segment has a successor, submit it now.
# This avoids upfront dependency graphs (afterok:) which OSCAR's SLURM
# cancels instead of queuing when the CPU cap is hit simultaneously.
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
    print(p.get('_mem', '12G'))
except:
    print('12G')
" "$NEXT_PARAMS" 2>/dev/null)
    sbatch --no-requeue \
        --time="$WALLTIME" \
        --mem="$MEM" \
        --cpus-per-task="${SLURM_CPUS_PER_TASK:-4}" \
        --export="NONE,PARAMS=$NEXT_PARAMS,DUMP=$NEXT_DUMP" \
        "$PROJECT_ROOT/config/slurm_template.sh"
    echo "Submitted next segment: $NEXT_RUN"
fi
