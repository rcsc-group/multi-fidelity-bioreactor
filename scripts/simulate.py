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

_DEFAULT_TEMPLATE = Path(__file__).parents[1] / "config" / "slurm_mpi_template.sh"
_DEFAULT_BINARY   = Path(__file__).parents[1] / "build" / "BioReactor"


def _t_mix_nd(params: dict) -> float:
    """Compute the non-dim time at which oxygen/tracer transfer starts.

    Mirrors BioReactor.c: t_mix = T_per_st * n_mix_cycles.
    """
    omega_b      = params.get("omega_b", 3.93)
    L            = params.get("geometry", {}).get("a", 0.25)
    H            = params.get("geometry", {}).get("b", 0.071)
    th           = math.radians(params.get("theta_max", [7.0])[0])
    n_mix_cycles = params.get("n_mix_cycles", 80)

    T_per   = 2 * math.pi / omega_b
    V       = L / 4 * (H + 0.5 * L * math.tan(th))
    U       = V / (H * 0.5) / T_per
    T_bio   = L / U
    T_per_st = T_per / T_bio
    return T_per_st * n_mix_cycles


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
) -> str:
    """Write params.json and submit a SLURM job via sbatch.

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

    Returns
    -------
    SLURM job ID as a string (e.g. "123456").
    """
    project_root = Path(project_root) if project_root else Path(__file__).parents[1]
    runs_root    = Path(runs_root) if runs_root else project_root / "runs"
    template     = Path(template) if template else _DEFAULT_TEMPLATE
    _check_t_end_vs_t_mix(params)

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
        scratch_base = Path("/oscar/scratch/eaguerov/mpi_runs") / params["run_id"]
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
            checkpoint = str(scratch_ck)
    else:
        effective_params = params_path

    export_str = f"NONE,PARAMS={effective_params}"
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

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    # sbatch stdout: "Submitted batch job 123456"
    match = re.search(r"(\d+)", result.stdout)
    if not match:
        raise RuntimeError(f"Could not parse job ID from sbatch output: {result.stdout!r}")
    return match.group(1)


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
    args = parser.parse_args()

    p = json.loads(Path(args.params).read_text())
    root = Path(__file__).parents[1]

    if args.slurm:
        job_id = submit_slurm(p, project_root=root, walltime=args.walltime)
        print(f"Submitted job {job_id}")
        if args.wait:
            run_dir = root / "runs" / p.get("run_id", "unnamed")
            results = wait_for_result(run_dir)
            print(json.dumps(results, indent=2))
    else:
        run_dir = run_local(p, project_root=root)
        print(f"Completed: {run_dir}")
