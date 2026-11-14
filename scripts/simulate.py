"""Run a BioReactor simulation locally or via SLURM.

Public API
----------
run_local(params, project_root, runs_root)
    Write params.json, execute binary directly, return run_dir Path.

submit_slurm(params, project_root, runs_root, walltime, **sbatch_kwargs)
    Write params.json, submit via sbatch, return SLURM job_id string.

wait_for_result(run_dir, timeout, poll)
    Block until results.json appears in run_dir; return its contents as dict.
    Raises TimeoutError if the file does not appear within `timeout` seconds.

Typical HPC workflow
--------------------
    run_dir  = run_local(params, project_root)      # sets up directory
    job_id   = submit_slurm(params, project_root)   # submit to SLURM
    results  = wait_for_result(run_dir)              # block until done
    kLa      = results["kLa_25"]
"""
from __future__ import annotations

import json
import math
import re
import subprocess
import time
import warnings
from pathlib import Path
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root
from scripts import paths  # noqa: E402  (machine-specific folders)
from scripts import run_key as _rk  # noqa: E402


class DuplicateRun(Exception):
    """The index has a run with the same key (same physics, binary and parent): use its results.
    ``run_dir`` is that run's folder; ``finished`` tells if its results.json exists."""

    def __init__(self, run_id: str, run_dir: Path, finished: bool):
        self.run_id, self.run_dir, self.finished = run_id, run_dir, finished
        state = "finished" if finished else "submitted, no results yet (if it failed, resubmit it under its own run_id)"
        super().__init__(f"same run key as {run_id} ({state}): {run_dir}. "
                         f"Use its results, or pass allow_duplicate=True for a deliberate repeat.")

_DEFAULT_TEMPLATE = Path(__file__).parents[1] / "config" / "slurm_mpi_template.sh"
_DEFAULT_MPI_SCRATCH_ROOT = Path(f"{paths.SCRATCH}/mpi_runs")
_DEFAULT_BINARY   = Path(__file__).parents[1] / "build" / "BioReactor"


def _t_mix_nd(params: dict) -> float:
    """Compute the non-dim time at which oxygen/tracer transfer starts.

    Mirrors BioReactor.c: t_mix = T_per_st * n_mix_cycles.
    """
    omega_b      = params.get("omega_b", 3.93)
    L            = params.get("geometry", {}).get("a", 0.25)
    H            = 2 * params.get("geometry", {}).get("b", 0.03575)  # full height = 2x half-height semi-axis (BioReactor.c:295, fixed 2026-08-03)
    th           = math.radians(params.get("theta_max", [7.0])[0])
    n_mix_cycles = params.get("n_mix_cycles", 80)

    T_per   = 2 * math.pi / omega_b
    V       = L / 4 * (H + 0.5 * L * math.tan(th))
    U       = V / (H * 0.5) / T_per
    T_bio   = L / U
    T_per_st = T_per / T_bio
    return T_per_st * n_mix_cycles


RESTART_PREV_KEYS = ("omega_b_prev", "theta_max_prev")
# Binaries built from src/ since 2026-11-13 take a restart's source forcing from the dump's label
# (checkpoint.phase, src/forcing_state.h) and refuse *_prev keys. Older binaries need *_prev. The new solver
# prints "restart refused", so the string marks a binary that reads the label.
LABEL_MARKER = b"restart refused"


def binary_reads_label(binary) -> bool:
    """True if `binary` takes the source forcing from the checkpoint label. None or a missing file means the
    binary built from this src/, which does."""
    if binary is None or not Path(binary).is_file():
        return True
    return LABEL_MARKER in Path(binary).read_bytes()


def source_forcing_keys(src: dict, binary) -> dict:
    """The *_prev keys a restart from a dump of run `src` needs with `binary`: none for a binary that reads the
    label, src's forcing for an older binary."""
    if binary_reads_label(binary):
        return {}
    return {"omega_b_prev": src["omega_b"], "theta_max_prev": list(src["theta_max"]),
            "phi_angular_prev": list(src["phi_angular"]), "amplitude_h_prev": list(src["amplitude_h"]),
            "phi_horizontal_prev": list(src["phi_horizontal"]), "omega_h_prev": src["omega_h"]}


def _check_restart_forcing(params: dict, checkpoint=None) -> None:
    """Refuse a restart that does not say, in exactly one place, which forcing it continues.

    Every restart ramps the forcing from its source to params.json over 3 periods. A binary that reads the
    label takes the source from checkpoint.phase next to the dump; *_prev keys next to such a label are a second
    source and are refused; an old label (phase only) still needs *_prev. An older binary needs *_prev: with
    none it ramps from zero rocking while the liquid still sloshes (fig9_l9_rpm37.5_ext1: velocity +49 % after
    the seam). A replay run (rec_dir set) is not driven by the forcing.
    Pitfall "restart-ramp" in docs_site/developer/solver-pitfalls.md.
    """
    if params.get("t_checkpoint", 0) <= 0 or "rec_dir" in params:
        return
    rid = params.get("run_id")
    prev = sorted(k for k in params if k.endswith("_prev"))
    if not binary_reads_label(params.get("_binary")):
        missing = [k for k in RESTART_PREV_KEYS if k not in params]
        if missing:
            raise ValueError(f"{rid}: restart without {missing}: this binary would ramp the forcing again from "
                             f"zero (pitfall restart-ramp, docs_site/developer/solver-pitfalls.md). Build the "
                             f"params with scripts/submit_continue.build_params.")
        return
    if not checkpoint or not Path(checkpoint).exists():
        return   # a dump that a chained job has yet to write: the solver checks its label
    label = Path(checkpoint).parent / "checkpoint.phase"
    if not label.exists():
        raise ValueError(f"{rid}: no checkpoint.phase (the dump's label) next to {checkpoint}: the solver would "
                         f"refuse the restart (pitfall restart-ramp). Stage it with submit_continue.stage_phase.")
    has_forcing = "omega_b" in json.loads(label.read_text())
    if has_forcing and prev:
        raise ValueError(f"{rid}: {prev} in params AND a label with the forcing next to the dump: two sources "
                         f"(pitfall restart-ramp). Remove the *_prev keys.")
    if not has_forcing and not prev:
        raise ValueError(f"{rid}: old label (phase only) next to the dump and no *_prev keys: the source "
                         f"forcing is unknown (pitfall restart-ramp).")


def _check_t_end_vs_t_mix(params: dict) -> None:
    """Warn if t_end is set and is less than t_mix — run will produce no kLa data."""
    t_end = params.get("t_end")
    if t_end is None:
        return  # solver default (250) is always > t_mix
    t_mix = _t_mix_nd(params)
    if t_end <= t_mix:
        warnings.warn(
            f"t_end={t_end:.1f} <= t_mix≈{t_mix:.1f} (n_mix_cycles="
            f"{params.get('n_mix_cycles', 80)} rocking cycles). "
            "The run will complete but tr_oxy.dat will be all zeros and kLa will be NaN.",
            UserWarning,
            stacklevel=3,
        )


def _prepare_run_dir(params: dict, runs_root: Path) -> Path:
    """Create run directory and write params.json; return run_dir."""
    run_id  = params.get("run_id", "unnamed")
    run_dir = runs_root / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "params.json").write_text(json.dumps(params, indent=2))
    return run_dir


def run_local(
    params: dict,
    project_root: Path | str | None = None,
    runs_root: Path | str | None = None,
    timeout: int = 86400,
) -> Path:
    """Execute BioReactor binary directly (blocks until completion or timeout).

    Parameters
    ----------
    params       : BioReactor parameter dict (will be written as params.json)
    project_root : repo root; defaults to two directories above this file
    runs_root    : parent directory for run dirs; defaults to project_root/runs
    timeout      : subprocess timeout in seconds (default 24 h)

    Returns
    -------
    Path of the run directory.
    """
    project_root = Path(project_root) if project_root else Path(__file__).parents[1]
    runs_root    = Path(runs_root) if runs_root else project_root / "runs"
    binary       = project_root / "build" / "BioReactor"

    _check_restart_forcing(params)
    _check_t_end_vs_t_mix(params)
    if not binary.exists():
        raise FileNotFoundError(f"BioReactor binary not found at {binary}; run 'make build'")

    run_dir = _prepare_run_dir(params, runs_root)
    try:
        subprocess.run(
            [str(binary.resolve()), "params.json"],
            cwd=run_dir, check=False, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        pass  # long-running sim; caller can check output files

    from scripts.postprocess import main as _postprocess
    _postprocess(str(run_dir))
    return run_dir


def submit_slurm(
    params: dict,
    project_root: Path | str | None = None,
    runs_root: Path | str | None = None,
    walltime: str = "04:00:00",
    template: Path | str | None = None,
    mem: str = "12G",
    cpus: int = 4,
    ntasks: int | None = None,
    checkpoint: str | None = None,
    dependency: str | None = None,
    begin: str | None = None,
    exclude: str | None = None,
    mpi_scratch_root: Path | str | None = None,
    account: str | None = None,
    qos: str | None = None,
    partition: str | None = None,
    allow_duplicate: bool = False,
    allow_unknown_binary: bool = False,
) -> str:
    """Write params.json and submit a SLURM job via sbatch.

    Before it submits, it computes the run key (scripts/run_key.py: physics, binary, parent run) and refuses
    (DuplicateRun) a run whose key is already in runs_root/_index.jsonl under another run_id, and (ValueError) a
    binary that is not on disk, and a run_id whose finished results came from other physics.

    Parameters
    ----------
    params       : BioReactor parameter dict
    project_root : repo root; defaults to two directories above this file
    runs_root    : parent directory for run dirs; defaults to project_root/runs
    walltime     : SLURM wall-clock limit (HH:MM:SS)
    template     : path to SLURM script; defaults to config/slurm_template.sh
    mem          : memory request (e.g. "12G")
    cpus         : CPUs per task
    ntasks       : MPI rank count; overrides #SBATCH --ntasks in template
    checkpoint   : absolute path to checkpoint.dump for restart runs; if set,
                   DUMP env var is exported so the binary receives it as argv[2]
    dependency   : SLURM dependency string, e.g. "afterok:12345"; passed as
                   --dependency to sbatch (enables chained job submission)
    exclude      : comma-separated node names to pass as --exclude, e.g.
                   "node2336" -- for routing around a node that reproducibly
                   kills jobs instantly (RunTime~1s, empty output, no
                   documented reason) rather than gambling on the scheduler
                   picking elsewhere (diary.md 2026-09-07)
    mpi_scratch_root : where to stage params/checkpoint for MPI jobs (real MPI
                   compute nodes can't see project_root); defaults to the
                   real OSCAR scratch path. Overridable so tests exercising
                   the MPI staging branch don't touch the real filesystem.
    allow_duplicate : submit even if the key is already indexed (a deliberate repeat, e.g. a determinism check)
    allow_unknown_binary : submit with a binary that is not on disk here; the run gets no usable key

    Returns
    -------
    SLURM job ID as a string (e.g. "123456").
    """
    project_root = Path(project_root) if project_root else Path(__file__).parents[1]
    runs_root    = Path(runs_root) if runs_root else project_root / "runs"
    template     = Path(template) if template else _DEFAULT_TEMPLATE
    _check_restart_forcing(params, checkpoint)
    _check_t_end_vs_t_mix(params)
    params, key_entry = _keyed(params, runs_root, template, project_root, checkpoint,
                               allow_duplicate, allow_unknown_binary)

    run_dir      = _prepare_run_dir(params, runs_root)
    params_path  = (run_dir / "params.json").resolve()
    logs_dir     = project_root / "logs"
    logs_dir.mkdir(exist_ok=True)

    # Use NONE instead of ALL to avoid SLURM "user env retrieval failed" on nodes
    # with LDAP issues.  The SLURM script only needs PARAMS and DUMP; all
    # executables use absolute paths so PATH is not required.  SLURM_* vars
    # (SLURM_CPUS_PER_TASK etc.) are always injected regardless of --export.
    #
    # MPI exception: /oscar/data/dharri15/ is not mounted on MPI compute nodes.
    # When the MPI template is used, stage params.json and checkpoint to
    # /oscar/scratch at submission time (login node has full filesystem access),
    # then pass the scratch paths as PARAMS/DUMP so srun workers can read them.
    _mpi_template = (project_root / "config" / "slurm_mpi_template.sh").resolve()
    _using_mpi    = Path(template).resolve() == _mpi_template
    if _using_mpi:
        import shutil as _shutil
        _scratch_root = Path(mpi_scratch_root) if mpi_scratch_root else _DEFAULT_MPI_SCRATCH_ROOT
        scratch_base = _scratch_root / params["run_id"]
        scratch_base.mkdir(parents=True, exist_ok=True)
        # Store canonical path so the MPI job can write results back to Lustre
        canon_params = json.loads(params_path.read_text())
        canon_params["_canonical_run_dir"] = str(run_dir.resolve())
        scratch_params = scratch_base / "params.json"
        scratch_params.write_text(json.dumps(canon_params, indent=2))
        effective_params = scratch_params
        if checkpoint:
            scratch_ck = scratch_base / "checkpoint.dump"
            _shutil.copy2(checkpoint, scratch_ck)
            label = Path(checkpoint).parent / "checkpoint.phase"   # the dump's label travels with it
            if label.exists() and label.resolve() != (scratch_base / "checkpoint.phase").resolve():
                _shutil.copy2(label, scratch_base / "checkpoint.phase")
            checkpoint = str(scratch_ck)
    else:
        effective_params = params_path
        if checkpoint:   # the solver runs in run_dir and reads the dump's label there
            import shutil as _shutil
            label = Path(checkpoint).parent / "checkpoint.phase"
            if label.exists() and label.resolve() != (run_dir / "checkpoint.phase").resolve():
                _shutil.copy2(label, run_dir / "checkpoint.phase")

    # The job sees only these variables (NONE): the submitting checkout and this machine's folders (scripts/paths.py)
    export_str = (f"NONE,PARAMS={effective_params},MFB_REPO={Path(project_root).resolve()},"
                  f"MFB_SCRATCH={paths.SCRATCH},BASILISK={paths.BASILISK}")
    if checkpoint:
        export_str += f",DUMP={checkpoint}"

    # MPI jobs use --mem-per-cpu (per rank) instead of --mem (total per node)
    mem_flag = "--mem-per-cpu" if _using_mpi else "--mem"
    cmd = [
        "sbatch",
        "--no-requeue",           # prevent SLURM from re-running on node failure
        f"--time={walltime}",
        f"{mem_flag}={mem}",
        f"--cpus-per-task={cpus}",
        f"--export={export_str}",
        str(template.resolve()),
    ]
    if ntasks is not None:
        cmd.insert(1, f"--ntasks={ntasks}")
    if dependency:
        cmd.insert(1, f"--dependency={dependency}")
    if begin:
        cmd.insert(1, f"--begin={begin}")
    if exclude:
        cmd.insert(1, f"--exclude={exclude}")
    # [PROJECT ADDED, 2026-09-10] account/qos override -- ONLY for an
    # explicitly user-authorized, narrowly-scoped exception (e.g.
    # mbessa-condo for one named experiment; see feedback_mbessa_condo_
    # scope.md). Never set these as a default; every call site in this
    # project omits them unless the user has named that specific job.
    if account:
        cmd.insert(1, f"--account={account}")
    if qos:
        cmd.insert(1, f"--qos={qos}")
    # [PROJECT ADDED, 2026-11-04] partition for the Campaign 2D split (batch 312 cpu cap, batch2 512 cpu cap).
    if partition:
        cmd.insert(1, f"--partition={partition}")

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # sbatch stdout: "Submitted batch job 123456"
    match = re.search(r"(\d+)", result.stdout)
    if not match:
        raise RuntimeError(f"Could not parse job ID from sbatch output: {result.stdout!r}")
    _rk.append_index(runs_root, key_entry | {"job": match.group(1)})
    return match.group(1)


def _keyed(params, runs_root, template, project_root, checkpoint, allow_duplicate, allow_unknown_binary):
    """The run key checks of submit_slurm. Returns (params with _run_key, index entry); raises on a refusal."""
    rid = params.get("run_id", "unnamed")
    binary = params.get("_binary") or _rk.default_binary(template, project_root)
    bid = _rk.binary_identity(binary)
    if bid is None and not allow_unknown_binary:
        raise ValueError(f"{rid}: binary {binary} is not on disk, so the run would have no known code. Build it, "
                         f"set _binary, or pass allow_unknown_binary=True.")
    parent_run = Path(checkpoint).parent.name if checkpoint else None
    known = {e["run_id"]: e for e in _rk.read_index(runs_root)}
    e = _rk.entry(rid, params, bid, parent_run, known.get(parent_run))
    finished = Path(runs_root) / rid / "results.json"
    if finished.exists():
        old = json.loads((finished.parent / "params.json").read_text()).get("_run_key")
        if old != e["key"]:
            raise ValueError(f"{rid}: {finished} exists and came from other params (key {old}); a new run "
                             f"would mix with it. Use a new run_id.")
        if not allow_duplicate:
            raise DuplicateRun(rid, finished.parent, True)
    dup = None if allow_duplicate else _rk.find_duplicate(runs_root, e)
    if dup:
        d = Path(runs_root) / dup["run_id"]
        raise DuplicateRun(dup["run_id"], d, (d / "results.json").exists())
    return {**params, "_run_key": e["key"]}, e


def get_job_elapsed_seconds(job_id: str) -> float | None:
    """Query sacct for a completed job's wall-clock Elapsed time, in seconds.

    SLURM's Elapsed format is [D-]HH:MM:SS (the day component and its dash
    are only present once a job has run past 24h).

    Returns None if sacct has no record for this job yet (e.g. queried too
    soon after completion) -- wall time is optional metadata, so callers
    should not fail a run just because it's momentarily unavailable.
    """
    result = subprocess.run(
        ["sacct", "-j", job_id, "--format=Elapsed", "--noheader", "-X"],
        capture_output=True, text=True,
    )
    elapsed_str = result.stdout.strip()
    if not elapsed_str:
        return None

    day_part, _, clock_part = elapsed_str.rpartition("-")
    days = int(day_part) if day_part else 0
    hours, minutes, seconds = (int(x) for x in clock_part.split(":"))
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def wait_for_result(
    run_dir: Path | str,
    timeout: float = 7200,
    poll: float = 30,
) -> dict:
    """Block until results.json appears in run_dir; return its contents.

    Parameters
    ----------
    run_dir : completed run directory (contains params.json)
    timeout : maximum seconds to wait before raising TimeoutError
    poll    : polling interval in seconds

    Returns
    -------
    dict with kLa_10, kLa_25, kLa_50.

    Raises
    ------
    TimeoutError if results.json does not appear within `timeout` seconds.
    """
    run_dir     = Path(run_dir)
    results_file = run_dir / "results.json"
    deadline    = time.monotonic() + timeout

    while time.monotonic() < deadline:
        if results_file.exists():
            return json.loads(results_file.read_text())
        time.sleep(poll)

    raise TimeoutError(
        f"results.json not found in {run_dir} after {timeout:.0f}s"
    )


def run_trial(
    params: dict,
    backend: str = "slurm",
    kla_key: str = "kLa_25",
    walltime: str = "04:00:00",
    timeout: float = 7200,
    project_root: Path | str | None = None,
    runs_root: Path | str | None = None,
) -> float:
    """Black-box trial: run one simulation and return a single kLa scalar.

    Parameters
    ----------
    params       : BioReactor parameter dict
    backend      : "slurm" (submit + wait) or "local" (blocking subprocess)
    kla_key      : which saturation level to return — "kLa_10", "kLa_25", or "kLa_50"
    walltime     : SLURM wall-clock limit (HH:MM:SS); ignored for local backend
    timeout      : seconds to wait for results.json (SLURM) or subprocess (local)
    project_root : repo root; defaults to two directories above this file
    runs_root    : parent for run dirs; defaults to project_root/runs

    Returns
    -------
    float — kLa value at the requested saturation level.
    Returns float('nan') if the key is missing or NaN (failed/incomplete run).
    Raises TimeoutError if the result does not appear within `timeout` seconds.
    """
    project_root = Path(project_root) if project_root else Path(__file__).parents[1]
    runs_root    = Path(runs_root) if runs_root else project_root / "runs"

    if backend == "local":
        run_dir = run_local(params, project_root=project_root,
                            runs_root=runs_root, timeout=int(timeout))
        results = json.loads((run_dir / "results.json").read_text())
    elif backend == "slurm":
        submit_slurm(params, project_root=project_root, runs_root=runs_root,
                     walltime=walltime)
        run_dir = runs_root / params.get("run_id", "unnamed")
        results = wait_for_result(run_dir, timeout=timeout)
    else:
        raise ValueError(f"backend must be 'slurm' or 'local', got {backend!r}")

    value = results.get(kla_key, math.nan)
    if value is None or value != value:   # None or NaN
        return math.nan
    return float(value)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run BioReactor simulation")
    parser.add_argument("params", help="Path to params.json")
    parser.add_argument("--slurm", action="store_true", help="Submit via SLURM instead of running locally")
    parser.add_argument("--walltime", default="04:00:00")
    parser.add_argument("--wait", action="store_true", help="Wait for results.json after SLURM submission")
    parser.add_argument("--template", default=None, help="SLURM script (default: config/slurm_mpi_template.sh)")
    parser.add_argument("--runs-root", default=None, help="Folder of run directories (default: runs/)")
    parser.add_argument("--allow-duplicate", action="store_true",
                        help="Submit even if a run with the same key exists (a deliberate repeat)")
    args = parser.parse_args()

    p = json.loads(Path(args.params).read_text())
    root = Path(__file__).parents[1]
    runs_root = Path(args.runs_root) if args.runs_root else root / "runs"

    if args.slurm:
        try:
            job_id = submit_slurm(p, project_root=root, runs_root=runs_root, walltime=args.walltime,
                                  template=args.template, allow_duplicate=args.allow_duplicate)
        except DuplicateRun as e:
            print(e)
            sys.exit(3)
        print(f"Submitted job {job_id}")
        if args.wait:
            run_dir = runs_root / p.get("run_id", "unnamed")
            results = wait_for_result(run_dir)
            print(json.dumps(results, indent=2))
    else:
        run_dir = run_local(p, project_root=root)
        print(f"Completed: {run_dir}")
