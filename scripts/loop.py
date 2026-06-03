"""Multi-fidelity Bayesian Optimisation loop for rocking bioreactor conditions.

Driven by a YAML config file (see config/bo_config.yaml for all options).

Usage
-----
    python scripts/loop.py config/bo_config.yaml

The loop can be interrupted and resumed: it checks how many runs are already
finished in the ExperimentData store and skips the initial DoE if complete.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
import yaml
from f3dasm import ExperimentData, create_sampler
from f3dasm.design import Domain

_PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(_PROJECT_ROOT))

import scripts.simulate as simulate
import scripts.suggest  as suggest_mod
import scripts.train_surrogate as ts_mod


# ── helpers ───────────────────────────────────────────────────────────────────

def _build_domain(spec: dict) -> Domain:
    """Build f3dasm Domain from param_space.yaml spec."""
    pspace = spec["parameters"]
    n_max  = spec["N_max"]
    d      = Domain()

    def _add(name: str, key: str | None = None) -> None:
        lo, hi = pspace[key or name]["bounds"]
        d.add_float(name, lo, hi)

    _add("omega_b")
    d.add_float("n_harmonics",
                pspace["n_harmonics"]["bounds"][0],
                pspace["n_harmonics"]["bounds"][1])

    for i in range(n_max):
        d.add_float(f"theta_max_{i}",    *pspace["theta_max"]["bounds"])
    for i in range(n_max):
        lo, hi = (0.0, 0.0) if i == 0 else pspace["phi_angular"]["bounds"]
        d.add_float(f"phi_angular_{i}", lo, hi)
    _add("omega_h")
    for i in range(n_max):
        d.add_float(f"amplitude_h_{i}",    *pspace["amplitude_h"]["bounds"])
    for i in range(n_max):
        d.add_float(f"phi_horizontal_{i}", *pspace["phi_horizontal"]["bounds"])
    for sub in ("a", "b", "n"):
        _add(f"geometry_{sub}", key=f"geometry.{sub}")
    _add("fill_level")
    d.add_float("fidelity", 1, 10)
    return d


def _compute_t_end(params: dict, t_buffer: float) -> float:
    """t_end = t_mix_nd(params) + t_buffer — minimum run length for kLa extraction."""
    return simulate._t_mix_nd(params) + t_buffer


def _ed_to_params(row: pd.Series, fidelity: int,
                  t_end: float | None, n_max: int) -> dict:
    """Convert a flat ExperimentData input row to a nested params dict.

    t_end=None means the caller will set it after computing t_mix.
    """
    def _vec(base: str) -> list:
        return [float(row.get(f"{base}_{i}", 0.0)) for i in range(n_max)]

    run_id = uuid4().hex[:8]
    d = {
        "run_id":       run_id,
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
    if t_end is not None:
        d["t_end"] = t_end
    return d


def _submit_and_wait(params: dict, cfg: dict) -> dict | None:
    """Submit one SLURM job and block until results.json appears."""
    runs_root  = _PROJECT_ROOT / "runs"
    job_id     = simulate.submit_slurm(
        params,
        project_root=_PROJECT_ROOT,
        runs_root=runs_root,
        walltime=cfg["walltime"],
    )
    # run_dir is created by simulate._prepare_run_dir using params["run_id"]
    run_dir    = runs_root / params["run_id"]
    print(f"    submitted job {job_id} → {run_dir.name}")
    try:
        results = simulate.wait_for_result(
            run_dir, timeout=cfg["job_timeout"], poll=30)
        return results
    except TimeoutError:
        print(f"    WARNING: job {job_id} timed out after {cfg['job_timeout']} s")
        return None


def _append_to_ed(ed: ExperimentData, params: dict, results: dict,
                  fidelity: int, phase: str = "unknown",
                  bo_iteration: int | None = None) -> ExperimentData:
    """Append one completed run to ExperimentData and return the updated object.

    Stores all 10 postprocessing KPIs plus provenance metadata (run_id,
    timestamp, phase, bo_iteration) so the ED is the complete record of
    every run rather than a subset.
    """
    import time
    row_in = {
        "run_id":          params.get("run_id", ""),
        "omega_b":         params["omega_b"],
        "n_harmonics":     params["n_harmonics"],
        "theta_max_0":     params["theta_max"][0],
        "theta_max_1":     params["theta_max"][1],
        "theta_max_2":     params["theta_max"][2],
        "phi_angular_0":   params["phi_angular"][0],
        "phi_angular_1":   params["phi_angular"][1],
        "phi_angular_2":   params["phi_angular"][2],
        "omega_h":         params["omega_h"],
        "amplitude_h_0":   params["amplitude_h"][0],
        "amplitude_h_1":   params["amplitude_h"][1],
        "amplitude_h_2":   params["amplitude_h"][2],
        "phi_horizontal_0": params["phi_horizontal"][0],
        "phi_horizontal_1": params["phi_horizontal"][1],
        "phi_horizontal_2": params["phi_horizontal"][2],
        "geometry_a":      params["geometry"]["a"],
        "geometry_b":      params["geometry"]["b"],
        "geometry_n":      params["geometry"]["n"],
        "fill_level":      params["fill_level"],
        "fidelity":        fidelity,
        "phase":           phase,
        "bo_iteration":    bo_iteration if bo_iteration is not None else -1,
        "completed_at":    time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    new_ed = ExperimentData(
        input_data=pd.DataFrame([row_in]),
        output_data=pd.DataFrame([{
            "kLa_10":       results.get("kLa_10",       float("nan")),
            "kLa_25":       results.get("kLa_25",       float("nan")),
            "kLa_50":       results.get("kLa_50",       float("nan")),
            "kLa_inst_10":  results.get("kLa_inst_10",  float("nan")),
            "kLa_inst_25":  results.get("kLa_inst_25",  float("nan")),
            "kLa_inst_50":  results.get("kLa_inst_50",  float("nan")),
            "dtmix_0.50":   results.get("dtmix_0.50",   float("nan")),
            "dtmix_0.75":   results.get("dtmix_0.75",   float("nan")),
            "dtmix_0.95":   results.get("dtmix_0.95",   float("nan")),
            "vor_mean":     results.get("vor_mean",     float("nan")),
        }]),
    )
    return ed + new_ed


# ── main loop ─────────────────────────────────────────────────────────────────

def run_loop(cfg: dict) -> None:
    exp_dir = _PROJECT_ROOT / cfg["experiment_dir"]
    exp_dir.mkdir(parents=True, exist_ok=True)

    param_space_path = _PROJECT_ROOT / "config" / "param_space.yaml"
    spec  = yaml.safe_load(param_space_path.read_text())
    n_max = spec["N_max"]

    lf_fidelity  = int(cfg["lf_fidelity"])
    hf_fidelity  = int(cfg["hf_fidelity"])
    n_lf_init    = int(cfg["n_lf_init"])
    n_hf_init    = int(cfg["n_hf_init"])
    n_iter       = int(cfg["n_iter"])
    kla_key      = cfg["kla_key"]
    n_candidates = int(cfg.get("n_candidates", 2000))
    t_buffer     = float(cfg.get("t_buffer", 150.0))

    # Load or create ExperimentData
    ed_store = exp_dir / "experiment_data"
    if ed_store.exists():
        ed = ExperimentData.from_file(str(exp_dir))
        print(f"Resuming from {exp_dir}: {len(ed.to_numpy()[0])} existing rows")
    else:
        domain = _build_domain(spec)
        ed     = ExperimentData(domain=domain, project_dir=str(exp_dir))
        print(f"New experiment at {exp_dir}")

    inp_df, out_df = ed.to_numpy()
    n_done = inp_df.shape[0]
    n_init_needed = n_lf_init + n_hf_init

    # ── Phase 1: initial DoE ─────────────────────────────────────────────────
    if n_done < n_init_needed:
        print(f"\nPhase 1: initial DoE ({n_lf_init} LF + {n_hf_init} HF runs)")
        domain = _build_domain(spec)
        # Generate LHS candidates
        ed_lhs = ExperimentData(domain=domain)
        ed_lhs = create_sampler("latin_sampler", seed=0).call(
            ed_lhs, n_samples=n_lf_init + n_hf_init)
        X_lhs, _ = ed_lhs.to_numpy()
        inp_lhs_df, _ = ed_lhs.to_pandas()

        for run_i in range(n_done, n_init_needed):
            fid = lf_fidelity if run_i < n_lf_init else hf_fidelity
            row = inp_lhs_df.iloc[run_i]
            params = _ed_to_params(row, fid, None, n_max)
            params["t_end"] = _compute_t_end(params, t_buffer)
            print(f"  DoE run {run_i + 1}/{n_init_needed}: "
                  f"fidelity={fid}, omega_b={params['omega_b']:.2f}, "
                  f"t_end={params['t_end']:.1f}")
            results = _submit_and_wait(params, cfg)
            if results is None:
                continue
            phase_tag = "doe_lf" if fid == lf_fidelity else "doe_hf"
            ed = _append_to_ed(ed, params, results, fid, phase=phase_tag)
            ed.store(str(exp_dir))
            print(f"    {kla_key}={results.get(kla_key, 'nan'):.4f}")

    # ── Phase 2: BO iterations ───────────────────────────────────────────────
    print(f"\nPhase 2: {n_iter} BO iterations (maximising {kla_key})")
    _, out_arr = ed.to_numpy()
    best_so_far = float(out_arr[:, list(
        ed.to_pandas()[1].columns).index(kla_key)].max())

    for it in range(n_iter):
        print(f"\n  Iteration {it + 1}/{n_iter}  (best so far: {best_so_far:.5f})")

        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            tmp_model = f.name

        ts_mod.main(str(exp_dir), tmp_model,
                    lf_fidelity=lf_fidelity, hf_fidelity=hf_fidelity,
                    kla_key=kla_key)

        candidate = suggest_mod.main(
            str(exp_dir), str(param_space_path),
            model_path=tmp_model,
            kla_key=kla_key,
            lf_fidelity=lf_fidelity,
            hf_fidelity=hf_fidelity,
            n_candidates=n_candidates,
            seed=it,
        )
        candidate["run_id"]   = uuid4().hex[:8]
        candidate["fidelity"] = hf_fidelity
        candidate["t_end"]    = _compute_t_end(candidate, t_buffer)
        print(f"  Suggested: omega_b={candidate['omega_b']:.2f}, "
              f"fill_level={candidate['fill_level']:.2f}, "
              f"theta_max[0]={candidate['theta_max'][0]:.1f}")

        results = _submit_and_wait(candidate, cfg)
        if results is None:
            continue
        kla_val = results.get(kla_key, float("nan"))
        print(f"  Result: {kla_key}={kla_val:.5f}")

        ed = _append_to_ed(ed, candidate, results, hf_fidelity,
                           phase="bo", bo_iteration=it)
        ed.store(str(exp_dir))

        if kla_val > best_so_far:
            best_so_far = kla_val
            print(f"  *** New best: {kla_key}={best_so_far:.5f} ***")

    # ── Summary ──────────────────────────────────────────────────────────────
    inp_final, out_final = ed.to_pandas()
    out_cols = list(out_final.columns)
    k_idx    = out_cols.index(kla_key)
    hf_mask  = inp_final["fidelity"] == hf_fidelity
    best_row = out_final.loc[hf_mask, kla_key].idxmax()

    print("\n" + "=" * 60)
    print(f"Optimisation complete.  Best {kla_key} = {best_so_far:.5f}")
    print("Best parameters:")
    print(json.dumps({
        k: float(inp_final.loc[best_row, k]) for k in inp_final.columns
    }, indent=2))
    print("=" * 60)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/loop.py <config.yaml>", file=sys.stderr)
        sys.exit(1)
    cfg_path = Path(sys.argv[1])
    cfg      = yaml.safe_load(cfg_path.read_text())
    run_loop(cfg)
