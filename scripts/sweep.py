"""JSON-driven multi-parameter sweep runner.

Any param in a params JSON can be "swept" by giving it a list of values.
The system expands to individual simulation configs using zip or cartesian
product, then either submits each run independently (default) or groups them
into checkpoint-restart chains (opt-in).

Sweep detection
---------------
- Scalar or dict value → fixed param, same in every simulation.
- List of scalars for a non-vector param → sweep list.
- List of scalars for a known 3-vector param (theta_max, phi_angular,
  amplitude_h, phi_horizontal) → fixed vector, NOT a sweep.
- List of lists or list of dicts → sweep list (nested form for vector params).

Expansion
---------
- All swept lists the same length N → zip: N simulations, element-wise.
- Swept lists with different lengths → cartesian product.

Checkpoint grouping (opt-in via "chain": true in _sweep)
-------------------
Checkpoint restart requires identical Basilisk grid structure (same fidelity
and geometry).  Simulations are clustered by (fidelity, a, b, n); each cluster
is submitted as an independent chain.  NOTE: MPI checkpoint restart is broken
in the current Basilisk version (stale coarse MPI ghost cells → SIGFPE on first
multigrid solve after injection).  Serial checkpoint restart works.

Sweep control options (in "_sweep" key of the JSON):
  chain:               false (default) → each combination is an independent
                       fresh MPI run; true → checkpoint-restart chains (serial)
  n_mix_cycles:        mixing cycles for each run / seg-0 of each chain (default 80)
  n_transition_cycles: mixing cycles for restart segments when chain=true (default 10)
  t_buffer:            non-dim kLa measurement window (default 150.0)
  walltime:            SLURM walltime string (default "04:00:00")
  submit:              true → sbatch; false → write params only (dry run)

Usage
-----
    python scripts/sweep.py config/sweep_example.json
"""
from __future__ import annotations

import itertools
import json
import math
import re
import sys
from pathlib import Path
from uuid import uuid4

_PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(_PROJECT_ROOT))


# ── walltime estimation ───────────────────────────────────────────────────────

_LOGSTATS_PAT = re.compile(
    r"t:\s*([\d.]+).*Wall clock time \(s\):\s*([\d.]+)"
)

def _collect_timing_data(runs_root: Path) -> list[tuple[int, int, float]]:
    """Return [(fidelity, ntasks, wall_s_per_t_nd), ...] from all logstats.dat.

    Includes both completed and timed-out runs — the last recorded row gives
    the instantaneous rate, which is valid regardless of whether the run
    finished.  Deduplicates by resolved logstats path so scratch-mirrored
    files are not double-counted.
    """
    scratch = Path("/oscar/scratch/eaguerov/mpi_runs")
    search_dirs = [runs_root, scratch] if scratch.exists() else [runs_root]

    seen: set[str] = set()
    data: list[tuple[int, int, float]] = []

    for root in search_dirs:
        for d in root.iterdir():
            if not d.is_dir():
                continue
            lf = d / "logstats.dat"
            pf = d / "params.json"
            key = str(lf.resolve())
            if key in seen or not lf.exists() or not pf.exists():
                continue
            seen.add(key)
            try:
                p = json.loads(pf.read_text())
                fid    = p.get("fidelity")
                ntasks = p.get("_ntasks", 1)
                if not fid:
                    continue
                rows = _LOGSTATS_PAT.findall(lf.read_text())
                if len(rows) < 10:
                    continue
                t_sim  = float(rows[-1][0])
                wall_s = float(rows[-1][1])
                if t_sim <= 0 or wall_s <= 0:
                    continue
                data.append((int(fid), int(ntasks), wall_s / t_sim))
            except Exception:
                pass
    return data


def _fit_timing_model(
    data: list[tuple[int, int, float]],
) -> tuple[float, float, float] | None:
    """OLS fit of log(rate) = alpha*fidelity + beta*log(ntasks) + gamma.

    Returns (alpha, beta, gamma) or None if fewer than 3 distinct points.
    Uses a simple numpy-free normal-equations solve to avoid adding a hard dep.
    """
    # deduplicate to one median per (fid, ntasks) bucket so large sweeps
    # don't dominate the fit
    import statistics
    from collections import defaultdict
    buckets: dict[tuple, list] = defaultdict(list)
    for fid, ntasks, rate in data:
        buckets[(fid, ntasks)].append(rate)
    if len(buckets) < 3:
        return None

    rows = []
    for (fid, ntasks), rates in buckets.items():
        # use 75th-percentile (conservative — faster jobs are the tail risk)
        s = sorted(rates)
        p75 = s[int(0.75 * len(s))]
        rows.append((fid, math.log(max(ntasks, 1)), math.log(p75)))

    # design matrix X (n×3): [fid, log(ntasks), 1]
    n = len(rows)
    Xt = [[r[0] for r in rows], [r[1] for r in rows], [1.0] * n]
    y  = [r[2] for r in rows]

    # normal equations: (X^T X) b = X^T y  — 3×3 system, solve by Cramer
    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    XtX = [[dot(Xt[i], Xt[j]) for j in range(3)] for i in range(3)]
    Xty = [dot(Xt[i], y) for i in range(3)]

    # Gaussian elimination with partial pivoting on 3×3
    A = [row[:] + [Xty[i]] for i, row in enumerate(XtX)]
    for col in range(3):
        pivot = max(range(col, 3), key=lambda r: abs(A[r][col]))
        A[col], A[pivot] = A[pivot], A[col]
        if abs(A[col][col]) < 1e-12:
            return None
        for row in range(col + 1, 3):
            f = A[row][col] / A[col][col]
            A[row] = [A[row][k] - f * A[col][k] for k in range(4)]
    coeffs = [0.0] * 3
    for i in range(2, -1, -1):
        coeffs[i] = (A[i][3] - sum(A[i][j] * coeffs[j] for j in range(i + 1, 3))) / A[i][i]
    return float(coeffs[0]), float(coeffs[1]), float(coeffs[2])


def estimate_walltime(
    fidelity: int,
    t_end: float,
    ntasks: int = 1,
    runs_root: Path | None = None,
    headroom: float = 1.20,
) -> str:
    """Predict SLURM walltime for a run, fitted from empirical logstats data.

    Strategy:
      1. If the requested (fidelity, ntasks) bucket has direct observations,
         use the 75th-percentile of those rates — avoids poor model extrapolation.
      2. Otherwise fall back to an OLS model:
           log(rate) = alpha*fidelity + beta*log(ntasks) + gamma
         fit on 75th-percentile rates per bucket.
      3. If fewer than 3 buckets exist, return "24:00:00" (conservative default).

    Output: HH:MM:SS rounded UP to the nearest 30 minutes after applying headroom.
    """
    if runs_root is None:
        runs_root = _PROJECT_ROOT / "runs"

    data = _collect_timing_data(runs_root)
    if not data:
        return "24:00:00"

    # Build per-bucket 75th-percentile rates
    from collections import defaultdict
    buckets: dict[tuple, list] = defaultdict(list)
    for fid, n, rate in data:
        buckets[(fid, n)].append(rate)

    def _p75(rates: list) -> float:
        s = sorted(rates)
        return s[int(0.75 * len(s))]

    # ── Direct observation ───────────────────────────────────────────────────
    key = (fidelity, ntasks)
    if key in buckets:
        rate_s_per_tnd = _p75(buckets[key])
    else:
        # ── Model fallback ───────────────────────────────────────────────────
        model_data = [(fid, n, _p75(r)) for (fid, n), r in buckets.items()]
        coeffs = _fit_timing_model(model_data)
        if coeffs is None:
            return "24:00:00"
        alpha, beta, gamma = coeffs
        rate_s_per_tnd = math.exp(
            alpha * fidelity + beta * math.log(max(ntasks, 1)) + gamma
        )

    predicted_s = rate_s_per_tnd * t_end * headroom

    # Round UP to nearest 30 minutes
    half_hours = math.ceil(predicted_s / 1800)
    total_s    = half_hours * 1800
    h, rem     = divmod(int(total_s), 3600)
    m          = rem // 60
    return f"{h:02d}:{m:02d}:00"

import scripts.simulate as simulate
from scripts.chain import _t_period_nd
from scripts.postprocess import validate_params

_DEFAULT_TEMPLATE = _PROJECT_ROOT / "config" / "slurm_template.sh"

# Known fixed 3-element vector params — a plain list value for these is the
# vector itself, not a sweep.  A list-of-lists value IS a sweep.
_KNOWN_3VECTORS = {"theta_max", "phi_angular", "amplitude_h", "phi_horizontal"}


# ── public API ────────────────────────────────────────────────────────────────

def detect_sweep_params(raw: dict) -> dict[str, list]:
    """Return {param_name: [val, ...]} for every sweep param in raw.

    Parameters with list values are classified as sweeps unless they are one
    of the known 3-vector params given as a plain list of scalars.
    The special '_sweep' key is always excluded.
    """
    sweep: dict[str, list] = {}
    for key, val in raw.items():
        if key.startswith("_"):   # skip _sweep, _comment, and any other metadata keys
            continue
        if not isinstance(val, list):
            continue
        if key in _KNOWN_3VECTORS:
            # Sweep only if elements are themselves lists or dicts
            if val and isinstance(val[0], (list, dict)):
                sweep[key] = val
        else:
            sweep[key] = val
    return sweep


def expand_combinations(sweep_params: dict[str, list]) -> list[dict]:
    """Expand sweep params into one dict per simulation.

    If all swept lists have the same length → zip (element-wise pairing).
    Otherwise → cartesian product.
    """
    if not sweep_params:
        return [{}]

    lengths = [len(v) for v in sweep_params.values()]
    keys = list(sweep_params.keys())

    if len(set(lengths)) == 1:
        # All same length → zip
        return [dict(zip(keys, vals)) for vals in zip(*sweep_params.values())]
    else:
        # Different lengths → cartesian
        return [
            dict(zip(keys, combo))
            for combo in itertools.product(*sweep_params.values())
        ]


def group_by_checkpoint_key(params_list: list[dict]) -> dict[tuple, list[dict]]:
    """Cluster params dicts by (fidelity, geometry_a, geometry_b, geometry_n).

    Only runs sharing the same Basilisk grid structure can be checkpoint-chained.
    Order within each group is preserved.
    """
    groups: dict[tuple, list[dict]] = {}
    for p in params_list:
        key = _checkpoint_key(p)
        groups.setdefault(key, []).append(p)
    return groups


def build_segment_list(
    group: list[dict],
    options: dict,
    swept_keys: frozenset[str] = frozenset(),
) -> list[dict]:
    """Add chain-restart fields to each params dict in the group.

    Segment 0 is a fresh run; segments 1+ restore from the previous segment.
    n_mix_cycles is set by chain logic (n_mix_cycles for seg 0,
    n_transition_cycles for segs 1+) UNLESS "n_mix_cycles" is in swept_keys,
    meaning the user explicitly swept over different n_mix_cycles values and
    those per-combo values should be preserved as-is.

    Fields added / set:
      run_id, n_mix_cycles (unless swept), t_end, and for restarts:
      t_checkpoint, omega_b_prev, theta_max_prev, phi_angular_prev,
      amplitude_h_prev, phi_horizontal_prev, omega_h_prev.
    """
    n_mix_default   = int(options.get("n_mix_cycles", 80))
    n_trans_default = int(options.get("n_transition_cycles", 10))
    t_buffer        = float(options.get("t_buffer", 150.0))

    segments: list[dict] = []
    t_checkpoint = 0.0
    prev: dict | None = None

    for k, raw_p in enumerate(group):
        p = {key: (list(val) if isinstance(val, list) else val)
             for key, val in raw_p.items()}

        p["run_id"] = uuid4().hex[:8]

        # Apply chain default for n_mix_cycles unless it was explicitly swept
        if "n_mix_cycles" not in swept_keys:
            p["n_mix_cycles"] = n_trans_default if k > 0 else n_mix_default

        T_per = _t_period_nd(p)
        p["t_end"] = p["n_mix_cycles"] * T_per + t_buffer

        if k > 0 and prev is not None:
            p["t_checkpoint"]        = t_checkpoint
            p["omega_b_prev"]        = float(prev.get("omega_b", 3.93))
            p["theta_max_prev"]      = list(prev.get("theta_max", [7.0, 0.0, 0.0]))
            p["phi_angular_prev"]    = list(prev.get("phi_angular", [0.0, 0.0, 0.0]))
            p["amplitude_h_prev"]    = list(prev.get("amplitude_h", [0.0, 0.0, 0.0]))
            p["phi_horizontal_prev"] = list(prev.get("phi_horizontal", [0.0, 0.0, 0.0]))
            p["omega_h_prev"]        = float(prev.get("omega_h", 0.0))

        segments.append(p)

        # Advance checkpoint time to next full period boundary after t_end
        t_end_abs = t_checkpoint + p["t_end"]
        n_per = int(t_end_abs / T_per) + 1
        t_checkpoint = n_per * T_per
        prev = p

    return segments


def parse_sweep_config(path: str | Path) -> tuple[list[dict], dict]:
    """Parse a JSON sweep config.

    Returns
    -------
    params_list : one fully-merged params dict per simulation
    options     : sweep control settings from the '_sweep' key
    """
    raw = json.loads(Path(path).read_text())
    options = dict(raw.pop("_sweep", {}))
    # Remove remaining metadata keys (_comment, etc.) — they are not sim params
    for k in [k for k in list(raw) if k.startswith("_")]:
        raw.pop(k)

    sweep_params = detect_sweep_params(raw)
    combos = expand_combinations(sweep_params)

    fixed = {k: v for k, v in raw.items() if k not in sweep_params}

    params_list: list[dict] = []
    for combo in combos:
        p = {k: (list(v) if isinstance(v, list) else v) for k, v in fixed.items()}
        p.update(combo)
        params_list.append(p)

    return params_list, options


def _init_experiment_store(path: Path, options: dict,
                            params_list: list[dict]) -> str:
    """Create the experiment store for this sweep and return its path.

    The store lives at experiments/<config_stem>/ and is the canonical diary
    for all runs produced by this sweep config.  Each completed segment will
    append its inputs and all 10 KPIs via postprocess._register_to_experiment_store.

    Always creates the directory and writes _sweep_metadata.json regardless of
    whether f3dasm is available — run_ids are appended to that file after
    submission so plot_convergence.py can load exactly the right runs.
    """
    import time as _time

    exp_name = path.stem                  # e.g. "sweep_fb_theta_l7"
    exp_dir  = _PROJECT_ROOT / "experiments" / exp_name
    exp_dir.mkdir(parents=True, exist_ok=True)

    meta = {
        "config_file":    str(path.resolve()),
        "created_at":     _time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime()),
        "n_runs_planned": len(params_list),
        "run_ids":        [],          # populated by submit_sweep after submission
        "sweep_options":  {k: v for k, v in options.items()
                           if not k.startswith("_")},
    }
    (exp_dir / "_sweep_metadata.json").write_text(json.dumps(meta, indent=2))

    # f3dasm ExperimentData store (optional — only if f3dasm is installed)
    try:
        import pandas as pd   # noqa: F401
        from f3dasm import ExperimentData  # noqa: F401
    except ImportError:
        pass

    return str(exp_dir)


def submit_sweep(path: str | Path) -> list[str]:
    """End-to-end: parse → expand → submit.

    Default mode (chain=false in _sweep): every combination is submitted as a
    fully independent fresh MPI run — no checkpoint dependency, all jobs can
    run in parallel.  Returns one job ID per combination.

    Chain mode (chain=true in _sweep): combinations sharing the same grid
    structure (fidelity + geometry) are grouped into checkpoint-restart chains.
    Only seg-0 of each chain is submitted; the SLURM script self-submits
    subsequent segments on completion.  Returns one job ID per chain.
    NOTE: MPI checkpoint restart is currently broken (Basilisk coarse-ghost
    bug); chain mode should only be used with the serial template.

    An ExperimentData store is created at experiments/<config_stem>/ in both
    modes and is the canonical provenance record for all runs in this sweep.
    """
    params_list, options = parse_sweep_config(path)

    runs_root    = _PROJECT_ROOT / "runs"
    walltime_cfg = options.get("walltime", "04:00:00")   # may be "auto"
    cpus         = int(options.get("cpus", 4))
    mem          = str(options.get("mem", "12G"))
    submit       = bool(options.get("submit", True))
    use_chain    = bool(options.get("chain", False))   # default: independent
    use_mpi      = bool(options.get("mpi", False))
    ntasks       = int(options.get("ntasks", 16)) if use_mpi else None
    _mpi_template = _PROJECT_ROOT / "config" / "slurm_mpi_template.sh"
    _ser_template = _DEFAULT_TEMPLATE
    active_template = _mpi_template if use_mpi else _ser_template

    experiment_dir = _init_experiment_store(Path(path), options, params_list)
    job_ids: list[str] = []

    if not use_chain:
        # ── Independent mode (default) ──────────────────────────────────────
        # Each combination is a self-contained fresh run with n_mix_cycles from
        # the config (or default 80).  All jobs are submitted upfront and are
        # fully independent — no next_run_id, no checkpoint dependency.
        n_mix = int(options.get("n_mix_cycles", 80))
        t_buf = float(options.get("t_buffer", 150.0))

        mpi_label = f"MPI ntasks={ntasks}" if use_mpi else "serial"
        print(f"  Experiment store: {experiment_dir or '(f3dasm unavailable)'}")
        print(f"  Mode: independent ({mpi_label}) — {len(params_list)} runs\n")

        auto_walltime = walltime_cfg == "auto"
        if auto_walltime:
            print("  Walltime: auto (fitted from empirical logstats + 20% headroom)\n")

        _submitted_run_ids: list[str] = []
        for i, raw_p in enumerate(params_list):
            p = {k: (list(v) if isinstance(v, list) else v) for k, v in raw_p.items()}
            p["run_id"]       = uuid4().hex[:8]
            p["n_mix_cycles"] = p.get("n_mix_cycles", n_mix)
            T_per             = _t_period_nd(p)
            p["t_end"]        = p["n_mix_cycles"] * T_per + t_buf

            walltime = (
                estimate_walltime(
                    fidelity=int(p["fidelity"]),
                    t_end=float(p["t_end"]),
                    ntasks=ntasks or 1,
                    runs_root=runs_root,
                )
                if auto_walltime
                else walltime_cfg
            )

            p["_walltime"]    = walltime
            p["_mem"]         = mem
            if ntasks is not None:
                p["_ntasks"] = ntasks   # used by chain self-submission & provenance
            if experiment_dir:
                p["_experiment_dir"] = experiment_dir

            run_dir = runs_root / p["run_id"]
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "params.json").write_text(json.dumps(p, indent=2))

            print(
                f"  [{i:3d}] run={p['run_id']}"
                f"  omega_b={p.get('omega_b', '?'):.4g}"
                f"  n_mix={p['n_mix_cycles']}"
                f"  t_end≈{p['t_end']:.1f}"
                f"  walltime={walltime}",
                end="",
            )

            _submitted_run_ids.append(p["run_id"])
            if submit:
                validate_params(p)
                job_id = simulate.submit_slurm(
                    p,
                    project_root=_PROJECT_ROOT,
                    runs_root=runs_root,
                    walltime=walltime,
                    template=active_template,
                    cpus=cpus,
                    mem=mem,
                    ntasks=ntasks,
                )
                job_ids.append(job_id)
                print(f"  → job {job_id}")
            else:
                job_ids.append(f"dry-{i}")
                print(f"  → dry run")

        # Record the submitted run IDs in the experiment metadata so that
        # plot_convergence.py --experiment can load exactly these runs.
        meta_path = Path(experiment_dir) / "_sweep_metadata.json"
        meta = json.loads(meta_path.read_text())
        meta["run_ids"] = _submitted_run_ids
        meta_path.write_text(json.dumps(meta, indent=2))

        print(f"\nDone. {len(job_ids)} independent run(s) submitted: {job_ids}")
        return job_ids

    # ── Chain mode (opt-in, chain=true) ─────────────────────────────────────
    # Preserved for serial+checkpoint sweeps.  Groups by grid structure and
    # submits seg-0 only; the SLURM script self-submits subsequent segments.
    groups     = group_by_checkpoint_key(params_list)
    swept_keys = frozenset(detect_sweep_params(
        json.loads(Path(path).read_text())
    ).keys())

    print(f"  Experiment store: {experiment_dir or '(f3dasm unavailable)'}")
    print(f"  Mode: chained (chain=true) — {len(groups)} chain(s)\n")

    for g_idx, (key, group) in enumerate(groups.items()):
        segments = build_segment_list(group, options, swept_keys=swept_keys)

        print(f"\nGroup {g_idx} (fidelity={key[0]}, a={key[1]}, b={key[2]}, n={key[3]})"
              f" — {len(segments)} segment(s)")

        prev_run_id: str | None = None
        for k, params in enumerate(segments):
            next_run_id = segments[k + 1]["run_id"] if k + 1 < len(segments) else None
            checkpoint  = (
                str((runs_root / prev_run_id / "checkpoint.dump").resolve())
                if prev_run_id is not None else None
            )
            params["next_run_id"] = next_run_id
            params["_walltime"]   = walltime
            params["_mem"]        = mem
            if experiment_dir:
                params["_experiment_dir"] = experiment_dir

            run_dir = runs_root / params["run_id"]
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "params.json").write_text(json.dumps(params, indent=2))
            if checkpoint:
                (run_dir / "checkpoint_path.txt").write_text(checkpoint)

            print(
                f"  [seg {k}] run={params['run_id']}"
                f"  omega_b={params.get('omega_b', '?'):.4g}"
                f"  n_mix={params['n_mix_cycles']}"
                f"  t_end≈{params['t_end']:.1f}"
                f"  → {'next:' + next_run_id[:8] if next_run_id else 'last'}",
                end="",
            )
            prev_run_id = params["run_id"]
            print()

        seg0_params = segments[0]
        if submit:
            validate_params(seg0_params)
            job_id = simulate.submit_slurm(
                seg0_params,
                project_root=_PROJECT_ROOT,
                runs_root=runs_root,
                walltime=walltime,
                cpus=cpus,
                mem=mem,
                checkpoint=None,
            )
            job_ids.append(job_id)
            print(f"  → submitted seg-0 as job {job_id} (chain self-submits from here)")
        else:
            job_ids.append(f"dry-{g_idx}")
            print(f"  → dry run (seg-0 params written, chain would self-submit)")

    print(f"\nDone. {len(job_ids)} chain(s) started: {job_ids}")
    return job_ids


def submit_sweep_videos(run_dirs: list[Path | str]) -> list[str]:
    """Submit independent video jobs for a list of completed run directories.

    Each video job depends only on itself having a valid params.json and
    checkpoint.dump (for restart segments).  Jobs are submitted with NO
    inter-job dependency so a single failure cannot cascade.

    Call this after all simulation jobs from submit_sweep have completed.

    Returns list of submitted SLURM job IDs.
    """
    video_template = _PROJECT_ROOT / "config" / "slurm_video_template.sh"
    if not video_template.exists():
        raise FileNotFoundError(f"Video template not found: {video_template}")

    runs_root = _PROJECT_ROOT / "runs"
    job_ids: list[str] = []

    for run_dir in run_dirs:
        run_dir = Path(run_dir)
        params_path = run_dir / "params.json"
        if not params_path.exists():
            print(f"  SKIP {run_dir.name}: no params.json")
            continue

        params = json.loads(params_path.read_text())

        # Pass checkpoint for restart segments so BioReactor-video can restore state
        checkpoint: str | None = None
        if params.get("t_checkpoint", 0.0) > 0.0:
            ck_path = run_dir / "checkpoint.dump"
            if ck_path.exists():
                checkpoint = str(ck_path.resolve())
            else:
                print(f"  SKIP {run_dir.name}: restart segment but checkpoint.dump missing")
                continue

        job_id = simulate.submit_slurm(
            params,
            project_root=_PROJECT_ROOT,
            runs_root=runs_root,
            walltime="00:30:00",
            template=video_template,
            checkpoint=checkpoint,
            dependency=None,   # fully independent — no cascade risk
        )
        job_ids.append(job_id)
        print(f"  {run_dir.name}  → video job {job_id}")

    return job_ids


# ── helpers ───────────────────────────────────────────────────────────────────

def _checkpoint_key(params: dict) -> tuple:
    g = params.get("geometry", {})
    if isinstance(g, dict):
        a, b, n = g.get("a", 0.25), g.get("b", 0.071), g.get("n", 8.0)
    else:
        a, b, n = 0.25, 0.071, 8.0
    # fill_level changes the interface height → incompatible checkpoint fields
    fill = round(float(params.get("fill_level", 0.5)), 4)
    # theta_max[0] groups chains by angle for parallel execution efficiency
    th = round(float((params.get("theta_max") or [7.0])[0]), 4)
    return (params.get("fidelity", 7), a, b, n, fill, th)


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/sweep.py <sweep_config.json>", file=sys.stderr)
        sys.exit(1)
    ids = submit_sweep(sys.argv[1])
    print(f"\nDone. {len(ids)} segment(s) submitted: {ids}")
