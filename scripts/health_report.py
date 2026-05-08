"""Numerical health report for a completed BioReactor run.

Usage
-----
    python scripts/health_report.py runs/health_l6/ runs/health_l7/

Prints three KPIs per run:
  1. VOF mass drift      — (max-min)/mean of f_liq_sum [%]; threshold < 0.1%
  2. Interface stability — mean f_liq_interf in 2nd half / 1st half;   threshold < 3
  3. Velocity quasi-steady — mean post-ramp vel_rms in 2nd half / 1st; threshold < 3
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


# ── loaders (mirrors conftest.py) ────────────────────────────────────────────

def _load(path: Path, skip_header: bool = True) -> np.ndarray:
    lines = [
        l for l in path.read_text().splitlines()
        if l.strip() and not (skip_header and l.strip().startswith("i"))
    ]
    return np.array([[float(x) for x in l.split()] for l in lines])


def load_vol_frac(run_dir: Path) -> np.ndarray:
    """(N,6): i t f_liq_sum f_liq_interf posY_max posY_min"""
    return _load(run_dir / "vol_frac_interf.dat")


def load_normf(run_dir: Path) -> np.ndarray:
    """(N,14): i t Omega_avg Omega_rms Omega_vol Omega_max ux_avg ux_rms ux_vol ux_max uy_avg uy_rms uy_vol uy_max"""
    return _load(run_dir / "normf.dat")


# ── non-dim ramp end time (mirrors test_quasi_steady_flow.py) ────────────────

def _t_ramp_nd(params: dict) -> float:
    omega_b = params["omega_b"]
    L = params["geometry"]["a"]
    H = params["geometry"]["b"]
    th = math.radians(params["theta_max"][0])
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(th))
    U = V / (H * 0.5) / T_per
    T_bio = L / U
    T_per_st = T_per / T_bio
    return 3 * T_per_st  # N_RAMP_CYCLES = 3


# ── KPI computations ─────────────────────────────────────────────────────────

def kpi_mass_drift(vf: np.ndarray) -> tuple[float, str]:
    """VOF mass drift: (max-min)/mean of f_liq_sum [%]."""
    f = vf[:, 2]
    drift_pct = (f.max() - f.min()) / f.mean() * 100
    status = "OK" if drift_pct < 0.1 else "FAIL"
    return drift_pct, status


def kpi_interface_stability(vf: np.ndarray) -> tuple[float, str]:
    """Interface area ratio: mean 2nd half / mean 1st half of f_liq_interf."""
    iface = vf[:, 3]
    # skip first row (t=0, interface not yet initialised)
    iface = iface[iface > 0]
    if len(iface) < 4:
        return float("nan"), "SKIP"
    mid = len(iface) // 2
    ratio = iface[mid:].mean() / (iface[:mid].mean() + 1e-30)
    status = "OK" if ratio < 3.0 else "FAIL"
    return ratio, status


def kpi_velocity_steady(normf: np.ndarray, t_ramp: float) -> tuple[float, str]:
    """Velocity RMS ratio: 2nd half / 1st half of post-ramp run."""
    t = normf[:, 1]
    vel_rms = np.sqrt(normf[:, 7] ** 2 + normf[:, 11] ** 2)
    post = vel_rms[t > t_ramp]
    if len(post) < 4:
        return float("nan"), "SKIP"
    mid = len(post) // 2
    ratio = post[mid:].mean() / (post[:mid].mean() + 1e-30)
    status = "OK" if ratio < 3.0 else "FAIL"
    return ratio, status


# ── report ───────────────────────────────────────────────────────────────────

def report(run_dir: Path) -> None:
    run_dir = Path(run_dir)
    params_path = run_dir / "params.json"

    if not params_path.exists():
        print(f"[{run_dir.name}] ERROR: params.json not found")
        return
    params = json.loads(params_path.read_text())

    vf_path    = run_dir / "vol_frac_interf.dat"
    norm_path  = run_dir / "normf.dat"

    if not vf_path.exists() or not norm_path.exists():
        print(f"[{run_dir.name}] ERROR: output files not found — simulation may not have completed")
        return

    vf   = load_vol_frac(run_dir)
    nf   = load_normf(run_dir)
    t_ramp = _t_ramp_nd(params)

    rows_vf   = len(vf)
    t_final   = vf[-1, 1] if rows_vf > 0 else 0.0

    drift,     drift_status = kpi_mass_drift(vf)
    iface_rat, iface_status = kpi_interface_stability(vf)
    vel_rat,   vel_status   = kpi_velocity_steady(nf, t_ramp)

    print(f"\n{'='*60}")
    print(f"  Run: {run_dir.name}   (fidelity={params.get('fidelity','?')}, "
          f"t_end={params.get('t_end','?')}, rows={rows_vf}, t_final={t_final:.1f})")
    print(f"{'='*60}")
    print(f"  KPI 1 — VOF mass drift        : {drift:.4f}%   [{drift_status}]  (threshold < 0.1%)")
    print(f"  KPI 2 — Interface area ratio  : {iface_rat:.3f}      [{iface_status}]  (threshold < 3.0)")
    print(f"  KPI 3 — Velocity RMS ratio    : {vel_rat:.3f}      [{vel_status}]  (threshold < 3.0)")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(description="Numerical health report for BioReactor runs")
    parser.add_argument("run_dirs", nargs="+", help="One or more run directories")
    args = parser.parse_args()

    for d in args.run_dirs:
        report(Path(d))


if __name__ == "__main__":
    main()
