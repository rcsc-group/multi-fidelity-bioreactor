#!/bin/bash
# SLURM template for MPI-parallel BioReactor runs.
# Uses srun --mpi=pmix to launch across ntasks MPI ranks.
# Build the binary first: make build-mpi
#
# IMPORTANT: submit via simulate.py with template=slurm_mpi_template.sh.
# simulate.py stages params.json to /oscar/scratch before submitting so
# PARAMS points to scratch (accessible from MPI compute nodes).
#
#SBATCH --job-name=bioreactor-mpi
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=04:00:00

#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elvis_alexander_aguero_vera@brown.edu

set -euo pipefail

module load openmpi
module load ffmpeg
export PATH="/oscar/data/dharri15/eaguerov/basilisk/src:$PATH"

# Set HOME so srun workers can initialise OpenMPI's opal layer.
if [ -z "${HOME:-}" ]; then
    export HOME=$(getent passwd "$(id -u)" | cut -d: -f6)
fi

if [ -z "${PARAMS:-}" ]; then
    echo "ERROR: PARAMS env var not set (must point to /oscar/scratch)." >&2
    exit 1
fi

SCRATCH_RUN="$(dirname "$PARAMS")"
echo "Scratch run  : $SCRATCH_RUN"
echo "params.json  : $PARAMS"
echo "MPI ranks    : ${SLURM_NTASKS:-16}"
echo "HOME         : $HOME"

# Canonical Lustre run dir (for results collection after the simulation)
CANON_RUN=$(python3 -c "
import json, sys
try:
    p = json.load(open(sys.argv[1]))
    print(p.get('_canonical_run_dir', ''))
except:
    print('')
" "$PARAMS" 2>/dev/null)

# Derive CANON_RUN from _experiment_dir when not set (chain-submitted jobs lack _canonical_run_dir)
if [ -z "$CANON_RUN" ]; then
    _EXP_DIR=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_experiment_dir',''))" "$PARAMS" 2>/dev/null)
    if [ -n "$_EXP_DIR" ]; then
        _RUN_ID=$(basename "$SCRATCH_RUN")
        CANON_RUN="$(dirname "$(dirname "$_EXP_DIR")")/runs/$_RUN_ID"
    fi
fi

# Binary must be in /oscar/scratch (accessible from compute nodes)
BINARY="/oscar/scratch/eaguerov/BioReactor-mpi-video"
_BINARY_OVERRIDE=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1])).get('_binary',''))" "$PARAMS" 2>/dev/null)
[ -n "$_BINARY_OVERRIDE" ] && BINARY="$_BINARY_OVERRIDE"
if [ ! -f "$BINARY" ]; then
    echo "ERROR: $BINARY not found. Run: cp build/BioReactor-mpi-video $BINARY" >&2
    exit 1
fi

unset DISPLAY
if [ -n "${DUMP:-}" ]; then
    srun --mpi=pmix --mem=0 --chdir="$SCRATCH_RUN" --export=HOME \
        "$BINARY" "$PARAMS" "$DUMP"
else
    srun --mpi=pmix --mem=0 --chdir="$SCRATCH_RUN" --export=HOME \
        "$BINARY" "$PARAMS"
fi

echo "Simulation complete. Syncing results..."

# Copy output files to canonical Lustre path for postprocessing
if [ -n "$CANON_RUN" ] && [ "$CANON_RUN" != "$SCRATCH_RUN" ]; then
    mkdir -p "$CANON_RUN"
    rsync -a --exclude="params.json" "$SCRATCH_RUN/" "$CANON_RUN/"
    PROJECT_ROOT="$(dirname "$(dirname "$CANON_RUN")")"
    uv run python "$PROJECT_ROOT/scripts/postprocess.py" "$CANON_RUN"
    echo "Done. Results written to $CANON_RUN/results.json"
else
    PROJECT_ROOT="$(cd "$SCRATCH_RUN/../../.." && pwd)"
    uv run python "$PROJECT_ROOT/scripts/postprocess.py" "$SCRATCH_RUN"
    echo "Done. Results written to $SCRATCH_RUN/results.json"
fi

# Self-submitting chain (next segment must also be staged to scratch)
NEXT_RUN=$(python3 -c "
import json, sys
try:
    p = json.load(open(sys.argv[1]))
    print(p.get('next_run_id', ''))
except:
    print('')
" "$PARAMS" 2>/dev/null)

if [ -n "$NEXT_RUN" ]; then
    # Derive canon runs/ dir from _experiment_dir in params (fallback: dirname of CANON_RUN)
    RUNS_ROOT=$(python3 -c "
import json, sys, os
p = json.load(open(sys.argv[1]))
exp = p.get('_experiment_dir', '')
if exp:
    print(os.path.join(os.path.dirname(os.path.dirname(exp)), 'runs'))
" "$PARAMS" 2>/dev/null)
    [ -z "$RUNS_ROOT" ] && [ -n "$CANON_RUN" ] && RUNS_ROOT="$(dirname "$CANON_RUN")"
    NEXT_CANON="$RUNS_ROOT/$NEXT_RUN"
    NEXT_PARAMS_CANON="$NEXT_CANON/params.json"
    if [ -f "$NEXT_PARAMS_CANON" ]; then
        NEXT_SCRATCH="/oscar/scratch/eaguerov/mpi_runs/$NEXT_RUN"
        mkdir -p "$NEXT_SCRATCH"
        # checkpoint.dump is in the current seg's output dir; fall back to scratch
        CURR_CKPT="${CANON_RUN:+$CANON_RUN/checkpoint.dump}"
        [ -z "$CURR_CKPT" ] && CURR_CKPT="$SCRATCH_RUN/checkpoint.dump"
        cp "$CURR_CKPT" "$NEXT_SCRATCH/checkpoint.dump" 2>/dev/null || true
        cp "$NEXT_PARAMS_CANON" "$NEXT_SCRATCH/params.json"
        WALLTIME=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_walltime','04:00:00'))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        MEM=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_mem','4G'))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        NTASKS=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_ntasks',16))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        MAIL_USER=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_mail_user',''))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        MAIL_TYPE=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_mail_type','FAIL'))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        NEXT_DUMP_ARG=""
        if [ -f "$NEXT_SCRATCH/checkpoint.dump" ]; then
            NEXT_DUMP_ARG="DUMP=$NEXT_SCRATCH/checkpoint.dump"
        fi
        # Resolve project root from RUNS_ROOT (runs/ is one level below project root)
        TEMPLATE="$RUNS_ROOT/../config/slurm_mpi_template.sh"
        [ ! -f "$TEMPLATE" ] && [ -n "$CANON_RUN" ] && TEMPLATE="$(dirname "$CANON_RUN")/../config/slurm_mpi_template.sh"
        MAIL_ARGS=()
        if [ -n "$MAIL_USER" ]; then
            MAIL_ARGS=(--mail-type="$MAIL_TYPE" --mail-user="$MAIL_USER")
        fi
        NEXT_JID=$(sbatch --no-requeue \
            --time="$WALLTIME" \
            --mem-per-cpu="$MEM" \
            --ntasks="$NTASKS" \
            --cpus-per-task=1 \
            "${MAIL_ARGS[@]}" \
            --export="NONE,PARAMS=$NEXT_SCRATCH/params.json${NEXT_DUMP_ARG:+,$NEXT_DUMP_ARG}" \
            "$TEMPLATE" | awk '{print $NF}')
        echo "$NEXT_JID" > "$NEXT_SCRATCH/.slurm_jid"
        echo "Submitted next segment: $NEXT_RUN (job $NEXT_JID)"
    fi
fi
