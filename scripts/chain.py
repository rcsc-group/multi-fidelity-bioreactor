"""Chained simulation sweep via checkpoint restart.

Each segment in the sweep is a separate SLURM job.  Segment k+1 restores the
checkpoint.dump written by segment k, so transients from the parameter change
decay in ~n_transition_cycles rocking cycles instead of the full n_mix_cycles.

Driven by a YAML config file — see config/chain_config.yaml.

Usage
-----
    python scripts/chain.py config/chain_config.yaml

Supported sweep parameters (scalar motion params):
    omega_b, omega_h, theta_max_0/1/2, amplitude_h_0/1/2,
    phi_angular_1/2, phi_horizontal_0/1/2

The sweep parameter name ``theta_max_0`` maps to ``theta_max[0]`` in the
params dict; analogously for other indexed vector params.

Results land in runs/<run_id>/results.json for each segment independently.
"""
from __future__ import annotations

import sys
from pathlib import Path
from uuid import uuid4

import yaml

_PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(_PROJECT_ROOT))

import scripts.simulate as simulate

_DEFAULT_TEMPLATE       = _PROJECT_ROOT / "config" / "slurm_template.sh"
_DEFAULT_VIDEO_TEMPLATE = _PROJECT_ROOT / "config" / "slurm_video_template.sh"

# Map sweep parameter names that encode a vector index (e.g. "theta_max_0")
# to the vector field name and integer index in the params dict.
_VECTOR_PARAMS: dict[str, tuple[str, int]] = {
    f"{field}_{i}": (field, i)
    for field in ("theta_max", "amplitude_h", "phi_angular", "phi_horizontal")
    for i in range(3)
}


def _t_period_nd(params: dict) -> float:
    """Non-dimensional rocking period T_per / T_bio."""
    return simulate._t_mix_nd({**params, "n_mix_cycles": 1})


def _apply_sweep_param(params: dict, sweep_param: str, value: float) -> dict:
    """Return a copy of params with the named sweep parameter set to value."""
    p = {k: (list(v) if isinstance(v, list) else v) for k, v in params.items()}

    if sweep_param in _VECTOR_PARAMS:
        field, idx = _VECTOR_PARAMS[sweep_param]
        vec = list(p.get(field, [0.0, 0.0, 0.0]))
        vec[idx] = float(value)
        p[field] = vec
    else:
        p[sweep_param] = float(value)
    return p


def build_chain(cfg: dict) -> list[dict]:
    """Return one params dict per sweep value.

    Normal mode (no initial_checkpoint):
      Segment 0 is a fresh run (n_mix_cycles).
      Segments 1+ restart from the previous segment's checkpoint (n_transition_cycles).

    Restart-from-existing mode (initial_checkpoint present):
      ALL segments restart — segment 0 from the provided checkpoint,
      subsequent segments from their predecessor.  n_transition_cycles is used
      for every segment's n_mix_cycles.
    """
    motion          = cfg["motion"]
    sweep_param     = cfg["sweep"]["parameter"]
    sweep_values    = cfg["sweep"]["values"]
    fidelity        = int(cfg["fidelity"])
    geometry        = dict(cfg["geometry"])
    fill_level      = float(cfg["fill_level"])
    n_mix_cycles    = int(cfg["n_mix_cycles"])
    n_transition    = int(cfg["n_transition_cycles"])
    t_buffer        = float(cfg["t_buffer"])

    initial_ck      = cfg.get("initial_checkpoint")   # optional

    # T_per_nd is constant for fixed geometry/theta_max (independent of omega_b).
    # Use seg-0 params to compute it once for all checkpoint time calculations.
    base0 = {
        "omega_b":    float(motion.get("omega_b", 1.0)),
        "theta_max":  list(motion.get("theta_max", [5.0, 0.0, 0.0])),
    }
    base0 = _apply_sweep_param(base0, sweep_param, sweep_values[0])
    T_per_nd = _t_period_nd({**base0, "geometry": geometry})

    # Seed checkpoint state from initial_checkpoint or zero (fresh start).
    if initial_ck:
        t_checkpoint = float(initial_ck["t_dump"])
        prev_omega_b: float | None = float(initial_ck["omega_b"])
        prev_motion: dict | None   = {
            "theta_max":      list(initial_ck.get("theta_max",      [5.0, 0.0, 0.0])),
            "phi_angular":    list(initial_ck.get("phi_angular",    [0.0, 0.0, 0.0])),
            "omega_h":        float(initial_ck.get("omega_h",       0.0)),
            "amplitude_h":    list(initial_ck.get("amplitude_h",    [0.0, 0.0, 0.0])),
            "phi_horizontal": list(initial_ck.get("phi_horizontal", [0.0, 0.0, 0.0])),
        }
    else:
        t_checkpoint = 0.0
        prev_omega_b = None
        prev_motion  = None

    chain: list[dict] = []

    for k, val in enumerate(sweep_values):
        # Base motion params with sweep value applied
        base = {
            "omega_b":        float(motion.get("omega_b", 1.0)),
            "n_harmonics":    int(motion.get("n_harmonics", 1)),
            "theta_max":      list(motion.get("theta_max", [5.0, 0.0, 0.0])),
            "phi_angular":    list(motion.get("phi_angular", [0.0, 0.0, 0.0])),
            "omega_h":        float(motion.get("omega_h", 0.0)),
            "amplitude_h":    list(motion.get("amplitude_h", [0.0, 0.0, 0.0])),
            "phi_horizontal": list(motion.get("phi_horizontal", [0.0, 0.0, 0.0])),
        }
        base = _apply_sweep_param(base, sweep_param, val)

        is_restart = (k > 0) or (initial_ck is not None)
        n_mix = n_transition if is_restart else n_mix_cycles
        t_end = n_mix * T_per_nd + t_buffer   # relative to this segment's start

        params = {
            "run_id":       uuid4().hex[:8],
            "fidelity":     fidelity,
            "geometry":     geometry,
            "fill_level":   fill_level,
            "n_mix_cycles": n_mix,
            "t_end":        t_end,
            **base,
        }
        if is_restart:
            # Restart segment: checkpoint time, omega_b rescaling, and full
            # prev-motion params for smooth-step parameter interpolation in C.
            params["t_checkpoint"]        = t_checkpoint
            params["omega_b_prev"]        = prev_omega_b
            params["theta_max_prev"]      = list(prev_motion["theta_max"])
            params["phi_angular_prev"]    = list(prev_motion["phi_angular"])
            params["amplitude_h_prev"]    = list(prev_motion["amplitude_h"])
            params["phi_horizontal_prev"] = list(prev_motion["phi_horizontal"])
            params["omega_h_prev"]        = float(prev_motion["omega_h"])

        chain.append(params)

        # Compute this segment's checkpoint time (for the next segment's t_checkpoint).
        t_end_abs = t_checkpoint + t_end
        n_per = int(t_end_abs / T_per_nd) + 1
        t_checkpoint = n_per * T_per_nd
        prev_omega_b = float(base["omega_b"])
        prev_motion  = base

    return chain


def submit_chain(cfg: dict) -> list[str]:
    """Build and submit the chain.  Returns list of SLURM job IDs.

    Set ``videos: true`` in the config to use the video binary
    (BioReactor-video + render_videos.py) for every segment.
    Omit or set ``videos: false`` for the standard kLa-only binary.
    """
    chain       = build_chain(cfg)
    runs_root   = _PROJECT_ROOT / "runs"
    walltime    = cfg.get("walltime", "02:00:00")
    submit      = bool(cfg.get("submit", True))
    videos      = bool(cfg.get("videos", False))
    template    = _DEFAULT_VIDEO_TEMPLATE if videos else _DEFAULT_TEMPLATE

    job_ids: list[str] = []
    prev_run_id: str | None = None
    initial_ck  = cfg.get("initial_checkpoint")

    for k, params in enumerate(chain):
        if k == 0 and initial_ck:
            # First segment restarts from an externally supplied checkpoint.
            checkpoint = str(Path(initial_ck["checkpoint_path"]).resolve())
        elif prev_run_id is not None:
            checkpoint = str((runs_root / prev_run_id / "checkpoint.dump").resolve())
        else:
            checkpoint = None
        dependency = f"afterok:{job_ids[k-1]}" if k > 0 else None

        print(
            f"  [seg {k}] run={params['run_id']}  "
            f"{cfg['sweep']['parameter']}={params.get(cfg['sweep']['parameter'], '?'):.3g}  "
            f"n_mix={params['n_mix_cycles']}  t_end≈{params['t_end']:.1f}"
            f"{'  [video]' if videos else ''}",
            end="",
        )

        if submit:
            job_id = simulate.submit_slurm(
                params,
                project_root=_PROJECT_ROOT,
                runs_root=runs_root,
                walltime=walltime,
                template=template,
                checkpoint=checkpoint,
                dependency=dependency,
            )
            job_ids.append(job_id)
            print(f"  → job {job_id}")
        else:
            import json
            run_dir = runs_root / params["run_id"]
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "params.json").write_text(json.dumps(params, indent=2))
            if checkpoint:
                (run_dir / "checkpoint_path.txt").write_text(checkpoint)
            print("  → params written (not submitted)")
            job_ids.append(f"dry-run-{k}")

        prev_run_id = params["run_id"]

    return job_ids


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/chain.py <config.yaml>", file=sys.stderr)
        sys.exit(1)
    cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
    ids = submit_chain(cfg)
    print(f"\nChain submitted: {ids}")
