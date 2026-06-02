"""JSON-driven multi-parameter sweep runner.

Any param in a params JSON can be "swept" by giving it a list of values.
The system expands to individual simulation configs using zip or cartesian
product, groups them by checkpoint compatibility (fidelity + geometry), and
submits each group as a checkpoint-restart chain.

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

Checkpoint grouping
-------------------
Checkpoint restart requires identical Basilisk grid structure (same fidelity
and geometry).  Simulations are clustered by (fidelity, a, b, n); each cluster
is submitted as an independent chain.

Sweep control options (in "_sweep" key of the JSON):
  n_mix_cycles:        mixing cycles for segment 0 of each group (default 80)
  n_transition_cycles: mixing cycles for restart segments (default 10)
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
import sys
from pathlib import Path
from uuid import uuid4

_PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(_PROJECT_ROOT))

import scripts.simulate as simulate
from scripts.chain import _t_period_nd

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


def submit_sweep(path: str | Path) -> list[str]:
    """End-to-end: parse → expand → group → chain → submit.

    Returns a flat list of all simulation SLURM job IDs (or dry-run tags).

    Videos are NOT submitted inline.  Submit them after all simulations
    complete via submit_sweep_videos() to avoid co-scheduling video and sim
    jobs on the same node — a video job failure would otherwise CANCEL the
    next sim segment via SLURM's dependency propagation.
    """
    params_list, options = parse_sweep_config(path)
    groups = group_by_checkpoint_key(params_list)

    runs_root    = _PROJECT_ROOT / "runs"
    walltime     = options.get("walltime", "04:00:00")
    submit       = bool(options.get("submit", True))
    # Stagger seg-0 start times so chains don't pile up on the SLURM scheduler.
    # When N chains all release their seg-1 dependency at the same second, SLURM
    # cancels the excess with Reason=Dependency even if the parent completed 0:0.
    # 300s between chain starts (5 min) is enough for the previous chain's seg-0
    # to have started before the next chain's seg-0 is released.
    chain_stagger_s = int(options.get("chain_stagger_s", 10))

    all_job_ids: list[str] = []

    swept_keys = frozenset(detect_sweep_params(
        json.loads(Path(path).read_text())
    ).keys())

    for g_idx, (key, group) in enumerate(groups.items()):
        segments = build_segment_list(group, options, swept_keys=swept_keys)
        prev_run_id: str | None = None

        print(f"\nGroup {g_idx} (fidelity={key[0]}, a={key[1]}, b={key[2]}, n={key[3]})"
              f" — {len(segments)} segment(s)")

        for k, params in enumerate(segments):
            checkpoint = (
                str((runs_root / prev_run_id / "checkpoint.dump").resolve())
                if prev_run_id is not None else None
            )
            dependency = f"afterok:{all_job_ids[-1]}" if k > 0 else None
            # Stagger seg-0 start times to prevent simultaneous dependency release.
            # Chain g_idx starts g_idx * chain_stagger_s seconds from now.
            begin = f"now+{g_idx * chain_stagger_s}" if (k == 0 and g_idx > 0) else None

            print(
                f"  [seg {k}] run={params['run_id']}"
                f"  omega_b={params.get('omega_b', '?'):.4g}"
                f"  n_mix={params['n_mix_cycles']}"
                f"  t_end≈{params['t_end']:.1f}",
                end="",
            )

            if submit:
                job_id = simulate.submit_slurm(
                    params,
                    project_root=_PROJECT_ROOT,
                    runs_root=runs_root,
                    walltime=walltime,
                    checkpoint=checkpoint,
                    dependency=dependency,
                    begin=begin,
                )
                all_job_ids.append(job_id)
                print(f"  → job {job_id}")
            else:
                run_dir = runs_root / params["run_id"]
                run_dir.mkdir(parents=True, exist_ok=True)
                (run_dir / "params.json").write_text(json.dumps(params, indent=2))
                if checkpoint:
                    (run_dir / "checkpoint_path.txt").write_text(checkpoint)
                print("  → params written (dry run)")
                all_job_ids.append(f"dry-{k}")

            prev_run_id = params["run_id"]

    return all_job_ids


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
