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
# [PROJECT ADDED, 2026-09-15] Pin every MPI job to ONE node. Multi-node MPI
# has a real, unfixed output-corruption bug in this project: the 240-task /
# 5-node L10 scaling probe (job 6190727) completed its 10 cycles per
# Basilisk's own end-of-run print, but all four output files were
# null-byte-padded from just after the header to a 4096-byte boundary -- a
# parallel-filesystem write race on the first block, never seen in any
# single-node job (diary.md 2026-09-10 (2)). Re-verified 2026-09-15: that
# run's shear_stress.dat still raises ValueError in np.loadtxt, i.e. the
# data is permanently unusable. Without this line SLURM is free to spread a
# large --ntasks across nodes and silently corrupt the run, which matters
# now that 48/64-rank requests are being used. Remove only when that write
# race is actually fixed and re-tested.
#SBATCH --nodes=1
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

# [PROJECT ADDED, 2026-09-21] Tell the solver how much wall clock it has, so
# it can checkpoint and stop instead of being killed mid-run.
#
# event dump_checkpoint fires exactly ONCE, at t_end. Before this, a job that
# ran out of clock first lost everything it had computed: on 2026-09-21
# fig9_l9_rpm20 was 49% through 24 h of L9 with no checkpoint on disk, and
# SLURM refuses to raise TimeLimit on a RUNNING job, so there was no recovery.
#
# Read from SLURM's own EndTime rather than from --time, because a job that
# was requeued or had its limit changed while pending has a --time that no
# longer matches when it will actually be killed. Computed immediately before
# srun so the figure is the time the SOLVER has, not the time the script had.
JOB_END=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | grep -oP 'EndTime=\K\S+' | head -1)
if [ -n "${JOB_END:-}" ] && [ "$JOB_END" != "Unknown" ]; then
    _END_EPOCH=$(date -d "$JOB_END" +%s 2>/dev/null || echo "")
    if [ -n "$_END_EPOCH" ]; then
        export BIOREACTOR_MAXRUNTIME_S=$(( _END_EPOCH - $(date +%s) ))
        echo "Walltime left: ${BIOREACTOR_MAXRUNTIME_S}s (solver reserves 300s to checkpoint)"
    fi
fi
# --export is a whitelist: a variable not named here does not reach the ranks,
# which would leave the deadline unset and silently restore the old behaviour.
SRUN_EXPORT="HOME"
[ -n "${BIOREACTOR_MAXRUNTIME_S:-}" ] && SRUN_EXPORT="HOME,BIOREACTOR_MAXRUNTIME_S"

if [ -n "${DUMP:-}" ]; then
    srun --mpi=pmix --mem=0 --chdir="$SCRATCH_RUN" --export="$SRUN_EXPORT" \
        "$BINARY" "$PARAMS" "$DUMP"
else
    srun --mpi=pmix --mem=0 --chdir="$SCRATCH_RUN" --export="$SRUN_EXPORT" \
        "$BINARY" "$PARAMS"
fi

echo "Simulation complete. Syncing results..."

# PROJECT_ROOT is hardcoded rather than derived from CANON_RUN/SCRATCH_RUN path
# arithmetic (2026-08-05, diary.md): that derivation silently resolved to
# /oscar/scratch whenever _canonical_run_dir was absent from the staged
# params.json (e.g. a raw sbatch submission that bypassed
# scripts.simulate._prepare_run_dir) -- hit twice this session (jobs 4645673,
# 4631100), both times the simulation completed fine but this postprocessing
# step failed with "can't open file '/oscar/scratch/scripts/postprocess.py'".
# This repo lives at a single fixed OSCAR path (see CLAUDE.md); no derivation
# needed.
PROJECT_ROOT="/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"

# Copy output files to canonical Lustre path for postprocessing
if [ -n "$CANON_RUN" ] && [ "$CANON_RUN" != "$SCRATCH_RUN" ]; then
    mkdir -p "$CANON_RUN"
    rsync -a --exclude="params.json" "$SCRATCH_RUN/" "$CANON_RUN/"
    uv run python "$PROJECT_ROOT/scripts/postprocess.py" "$CANON_RUN"
    echo "Done. Results written to $CANON_RUN/results.json"
else
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
    # RUNS_ROOT used to be derived from _experiment_dir/CANON_RUN, which is
    # empty for any self-submitted segment (only Python's submit_slurm()
    # stamps _canonical_run_dir, and only for segment 0 -- every later
    # segment's params.json is a plain `cp` of the pre-written canonical
    # copy, never annotated). That silently no-op'd this whole block from
    # segment 1 onward: NEXT_CANON became a bogus root-level path like
    # "/$NEXT_RUN", the `-f` check below failed, and the chain died with no
    # error (2026-09-03, diary.md -- same fragile-derivation bug the
    # 2026-08-05 PROJECT_ROOT fix addressed at the OTHER call site in this
    # file, never propagated here). PROJECT_ROOT is a single fixed path
    # (CLAUDE.md) -- just use it directly, no derivation needed.
    RUNS_ROOT="$PROJECT_ROOT/runs"
    NEXT_CANON="$RUNS_ROOT/$NEXT_RUN"
    NEXT_PARAMS_CANON="$NEXT_CANON/params.json"
    if [ -f "$NEXT_PARAMS_CANON" ]; then
        NEXT_SCRATCH="/oscar/scratch/eaguerov/mpi_runs/$NEXT_RUN"
        mkdir -p "$NEXT_SCRATCH"
        # checkpoint.dump is in the current seg's output dir; fall back to scratch
        CURR_CKPT="${CANON_RUN:+$CANON_RUN/checkpoint.dump}"
        [ -z "$CURR_CKPT" ] && CURR_CKPT="$SCRATCH_RUN/checkpoint.dump"
        cp "$CURR_CKPT" "$NEXT_SCRATCH/checkpoint.dump" 2>/dev/null || true
        # checkpoint.phase must travel WITH the dump: without it the next
        # segment falls back to inferring the phase by rounding to the
        # nearest period boundary, which is wrong by up to half a period for
        # any dump not taken at a zero-crossing (i.e. every emergency one).
        cp "$(dirname "$CURR_CKPT")/checkpoint.phase" \
           "$NEXT_SCRATCH/checkpoint.phase" 2>/dev/null || true
        # Stamp _canonical_run_dir into the copy so THIS segment's own
        # results-copy-back (top of this script, next time it runs) and its
        # own self-submission of the segment after it both work too --
        # otherwise the bug just recurs one hop later.
        python3 -c "
import json, sys
p = json.load(open(sys.argv[1]))
p['_canonical_run_dir'] = sys.argv[3]
json.dump(p, open(sys.argv[2], 'w'), indent=2)
" "$NEXT_PARAMS_CANON" "$NEXT_SCRATCH/params.json" "$NEXT_CANON"
        WALLTIME=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_walltime','04:00:00'))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        MEM=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_mem','4G'))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        NTASKS=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_ntasks',16))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        MAIL_USER=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_mail_user',''))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        MAIL_TYPE=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_mail_type','FAIL'))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        # [PROJECT FIXED, 2026-09-08] --exclude was applied to segment 0's
        # own submit_slurm() call (chain.py) but never propagated to THIS
        # self-submission, so segments 1+ of a long chain could still land
        # on a known-flaky node -- silently stalling the whole chain (a
        # killed segment never reaches "Simulation complete", so it never
        # self-submits the next one). Same bug class as the ntasks/mem
        # misses already fixed here; see chain.py's "_exclude" stamping.
        EXCLUDE=$(python3 -c "import json,sys; p=json.load(open(sys.argv[1])); print(p.get('_exclude',''))" "$NEXT_PARAMS_CANON" 2>/dev/null)
        NEXT_DUMP_ARG=""
        if [ -f "$NEXT_SCRATCH/checkpoint.dump" ]; then
            NEXT_DUMP_ARG="DUMP=$NEXT_SCRATCH/checkpoint.dump"
        fi
        TEMPLATE="$PROJECT_ROOT/config/slurm_mpi_template.sh"
        MAIL_ARGS=()
        if [ -n "$MAIL_USER" ]; then
            MAIL_ARGS=(--mail-type="$MAIL_TYPE" --mail-user="$MAIL_USER")
        fi
        EXCLUDE_ARGS=()
        if [ -n "$EXCLUDE" ]; then
            EXCLUDE_ARGS=(--exclude="$EXCLUDE")
        fi
        NEXT_JID=$(sbatch --no-requeue \
            --time="$WALLTIME" \
            --mem-per-cpu="$MEM" \
            --ntasks="$NTASKS" \
            --cpus-per-task=1 \
            "${MAIL_ARGS[@]}" \
            "${EXCLUDE_ARGS[@]}" \
            --export="NONE,PARAMS=$NEXT_SCRATCH/params.json${NEXT_DUMP_ARG:+,$NEXT_DUMP_ARG}" \
            "$TEMPLATE" | awk '{print $NF}')
        echo "$NEXT_JID" > "$NEXT_SCRATCH/.slurm_jid"
        echo "Submitted next segment: $NEXT_RUN (job $NEXT_JID)"
    fi
fi
