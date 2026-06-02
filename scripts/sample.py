"""Generate and (optionally) submit a batch of BioReactor simulation runs.

Driven by a YAML config file — see config/sample_config.yaml.

Usage
-----
    python scripts/sample.py config/sample_config.yaml

Supported sampling strategies (via f3dasm):
    latin   — Latin Hypercube Sampling (default, best space-filling)
    random  — uniform random
    grid    — full Cartesian grid (n_samples rounded to nearest perfect power)
    sobol   — Sobol quasi-random sequence

t_end is computed per run as t_mix(params) + t_buffer, so it adapts to
each parameter combination automatically.  No fixed t_end in the config.

Results are collected in an f3dasm ExperimentData store at experiment_dir.
Job IDs are printed as they are submitted.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import yaml
from f3dasm import ExperimentData, create_sampler

_PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(_PROJECT_ROOT))

import scripts.simulate as simulate
from scripts.loop import _build_domain, _compute_t_end

_STRATEGY_ALIASES = {
    "latin":  "latin_sampler",
    "random": "random_sampler",
    "grid":   "grid",
    "sobol":  "sobol_sampler",
}


def _row_to_params(row: pd.Series, fidelity: int, t_buffer: float,
                   n_max: int) -> dict:
    """Convert one ExperimentData input row to a nested params dict."""
    def _vec(base: str) -> list:
        return [float(row.get(f"{base}_{i}", 0.0)) for i in range(n_max)]

    params = {
        "run_id":       uuid4().hex[:8],
        "fidelity":     fidelity,
        "omega_b":      float(row["omega_b"]),
        "n_harmonics":  max(1, int(round(float(row["n_harmonics"])))),
        "theta_max":    _vec("theta_max"),
        "phi_angular":  [0.0] + [float(row.get(f"phi_angular_{i}", 0.0))
                                  for i in range(1, n_max)],
        "phi_horizontal": _vec("phi_horizontal"),
        "omega_h":      float(row["omega_h"]),
        "amplitude_h":  _vec("amplitude_h"),
        "geometry": {
            "a": float(row["geometry_a"]),
            "b": float(row["geometry_b"]),
            "n": float(row["geometry_n"]),
        },
        "fill_level":   float(row["fill_level"]),
    }
    params["t_end"] = _compute_t_end(params, t_buffer)
    return params


def run_sampling(cfg: dict) -> None:
    exp_dir  = _PROJECT_ROOT / cfg["experiment_dir"]
    exp_dir.mkdir(parents=True, exist_ok=True)

    param_space_path = _PROJECT_ROOT / "config" / "param_space.yaml"
    spec     = yaml.safe_load(param_space_path.read_text())
    n_max    = spec["N_max"]

    fidelity  = int(cfg["fidelity"])
    t_buffer  = float(cfg.get("t_buffer", 150.0))
    n_samples = int(cfg["sampling"]["n_samples"])
    seed      = int(cfg["sampling"].get("seed", 0))
    strategy  = cfg["sampling"]["strategy"].lower()
    submit    = bool(cfg.get("submit", True))
    walltime  = cfg.get("walltime", "01:00:00")

    f3dasm_name = _STRATEGY_ALIASES.get(strategy)
    if f3dasm_name is None:
        raise ValueError(
            f"Unknown sampling strategy {strategy!r}. "
            f"Choose from: {list(_STRATEGY_ALIASES)}")

    # Generate design points
    domain = _build_domain(spec)
    ed     = ExperimentData(domain=domain)
    ed     = create_sampler(f3dasm_name, seed=seed).call(ed, n_samples=n_samples)
    inp_df, _ = ed.to_pandas()

    print(f"Sampling: {strategy} | n={n_samples} | fidelity={fidelity} | "
          f"t_buffer={t_buffer} | submit={submit}")
    print(f"Experiment dir: {exp_dir}")

    job_ids: list[str] = []
    params_list: list[dict] = []

    for i, (_, row) in enumerate(inp_df.iterrows()):
        params = _row_to_params(row, fidelity, t_buffer, n_max)
        params_list.append(params)

        run_dir = _PROJECT_ROOT / "runs" / params["run_id"]
        print(f"  [{i+1:3d}/{n_samples}] run={params['run_id']}  "
              f"omega_b={params['omega_b']:.2f}  "
              f"fill={params['fill_level']:.2f}  "
              f"t_end={params['t_end']:.1f}", end="")

        if submit:
            job_id = simulate.submit_slurm(
                params,
                project_root=_PROJECT_ROOT,
                runs_root=_PROJECT_ROOT / "runs",
                walltime=walltime,
            )
            job_ids.append(job_id)
            print(f"  → job {job_id}")
        else:
            # Write params.json without submitting
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "params.json").write_text(json.dumps(params, indent=2))
            print(f"  → params written (not submitted)")

    # Persist the sampling plan to ExperimentData (inputs only; outputs added later)
    ed_store_path = exp_dir / "experiment_data"
    if ed_store_path.exists():
        existing = ExperimentData.from_file(str(exp_dir))
        combined = existing + ed
        combined.store(str(exp_dir))
    else:
        ed.store(str(exp_dir))

    if submit:
        print(f"\nSubmitted {len(job_ids)} jobs: {job_ids}")
    print(f"Sampling plan saved to {exp_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/sample.py <config.yaml>", file=sys.stderr)
        sys.exit(1)
    cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
    run_sampling(cfg)
