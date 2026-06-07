"""Numerical health report for a completed BioReactor run.

Usage
-----
    python scripts/health_report.py runs/health_l6/ runs/health_l7/

KPIs per run:
  1. VOF mass drift        — (max-min)/mean of f_liq_sum [%];       threshold < 0.5%
  2. Interface stability   — f_liq_interf 2nd/1st half mean;        threshold < 1.5
  3. Velocity quasi-steady — vel_rms 2nd/1st half mean (post-ramp); threshold [0.7, 1.5]
  4. CFL estimate          — dt * U_max / dx;                        threshold < 0.6
  5. KE quasi-steady       — KE 2nd/1st half (post-ramp);           threshold [0.7, 1.5]
  6. Poisson convergence   — max mgp.i from pressure_diag.dat
                             (only if DIAGNOSTICS=1 binary was used); threshold < NITERMAX (1000)
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

_PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(_PROJECT_ROOT))
import scripts.simulate as simulate


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

_N_RAMP_CYCLES = 3  # must match BioReactor.c: #define N_RAMP_CYCLES 3

def _t_ramp_nd(params: dict) -> float:
    return simulate._t_mix_nd({**params, "n_mix_cycles": _N_RAMP_CYCLES})


# ── KPI computations ─────────────────────────────────────────────────────────

def kpi_mass_drift(vf: np.ndarray) -> tuple[float, str]:
    """VOF mass drift: (max-min)/mean of f_liq_sum [%]."""
    f = vf[:, 2]
    drift_pct = (f.max() - f.min()) / f.mean() * 100
    status = "OK" if drift_pct < 0.5 else "FAIL"
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
    status = "OK" if ratio < 1.5 else "FAIL"
    return ratio, status


def kpi_velocity_steady(normf: np.ndarray, t_ramp: float) -> tuple[float, str]:
    """Velocity RMS ratio: 2nd half / 1st half of post-ramp run. Threshold [0.7, 1.5]."""
    t = normf[:, 1]
    vel_rms = np.sqrt(normf[:, 7] ** 2 + normf[:, 11] ** 2)
    post = vel_rms[t > t_ramp]
    if len(post) < 4:
        return float("nan"), "SKIP"
    mid = len(post) // 2
    ratio = post[mid:].mean() / (post[:mid].mean() + 1e-30)
    status = "OK" if 0.7 <= ratio <= 1.5 else "FAIL"
    return ratio, status


def _parse_logstats(run_dir: Path) -> list[tuple[int, float, float]]:
    """Parse logstats.dat → list of (i, t, dt) triples.

    Each line has format: 'i: X t: Y dt: Z #Cells: ...'
    """
    pattern = re.compile(r"i:\s*(\d+)\s+t:\s*(\S+)\s+dt:\s*(\S+)")
    result = []
    for line in (run_dir / "logstats.dat").read_text().splitlines():
        m = pattern.match(line.strip())
        if m:
            result.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    return result


def kpi_cfl(
    logstats: list[tuple[int, float, float]],
    normf: np.ndarray,
    fidelity: int,
) -> tuple[float, str]:
    """Estimated CFL = dt * U_max / dx. Threshold < 0.6.

    dx = 1/2^fidelity (non-dim uniform grid).
    U_max = max(|ux_liq_max|, |uy_liq_max|) from normf columns 9 and 13.
    dt matched from logstats by step index i (both files output at t+=0.1).
    """
    dx = 1.0 / (2**fidelity)
    dt_by_step = {row[0]: row[2] for row in logstats}

    cfl_vals = []
    for row in normf:
        step_i = int(row[0])
        if step_i not in dt_by_step:
            continue
        dt = dt_by_step[step_i]
        u_max = max(abs(row[9]), abs(row[13]))
        if u_max > 0:
            cfl_vals.append(dt * u_max / dx)

    if not cfl_vals:
        return float("nan"), "SKIP"

    cfl_max = float(np.max(cfl_vals))
    status = "OK" if cfl_max < 0.6 else "FAIL"
    return cfl_max, status


def kpi_kinetic_energy(normf: np.ndarray, t_ramp: float) -> tuple[float, str]:
    """KE quasi-steady ratio: 2nd/1st half mean post-ramp. Threshold [0.5, 2.0].

    KE ∝ ux_rms² + uy_rms² (V_liq constant by mass conservation).
    Lower bound catches stagnant flow; upper bound catches blowup.
    """
    t = normf[:, 1]
    ke = 0.5 * (normf[:, 7] ** 2 + normf[:, 11] ** 2)
    post = ke[t > t_ramp]
    if len(post) < 4:
        return float("nan"), "SKIP"
    mid = len(post) // 2
    ratio = float(post[mid:].mean() / (post[:mid].mean() + 1e-30))
    status = "OK" if 0.7 <= ratio <= 1.5 else "FAIL"
    return ratio, status


_NITERMAX = 1000  # must match BioReactor.c: NITERMAX = 1000


def kpi_pressure_residual(run_dir: Path) -> tuple[float, str]:
    """Max Poisson solver iteration count from pressure_diag.dat.

    project() passes tolerance=TOLERANCE/sq(dt), so mgp.resa scales with dt
    and is not directly comparable to a fixed threshold.  The correct health
    signal is max(mgp.i) < NITERMAX: if the solver exhausts its iteration
    budget, Basilisk prints a warning and the divergence-free constraint is
    not satisfied.  Only written by the DIAGNOSTICS=1 build (BioReactor-health).
    Returns (nan, 'SKIP') if the file is absent.
    """
    pdiag = Path(run_dir) / "pressure_diag.dat"
    if not pdiag.exists():
        return float("nan"), "SKIP"

    lines = [l for l in pdiag.read_text().splitlines()
             if l.strip() and not l.strip().startswith("i")]
    if not lines:
        return float("nan"), "SKIP"

    rows = np.array([[float(x) for x in l.split()] for l in lines])
    max_i = float(rows[:, 4].max())   # col 4: mgp_i
    status = "OK" if max_i < _NITERMAX else "FAIL"
    return max_i, status


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
    ke_rat,    ke_status    = kpi_kinetic_energy(nf, t_ramp)

    logstats_path = run_dir / "logstats.dat"
    if logstats_path.exists():
        logstats = _parse_logstats(run_dir)
        cfl_max, cfl_status = kpi_cfl(logstats, nf, params.get("fidelity", 8))
    else:
        cfl_max, cfl_status = float("nan"), "SKIP"

    max_iter, iter_status = kpi_pressure_residual(run_dir)

    print(f"\n{'='*60}")
    print(f"  Run: {run_dir.name}   (fidelity={params.get('fidelity','?')}, "
          f"t_end={params.get('t_end','?')}, rows={rows_vf}, t_final={t_final:.1f})")
    print(f"{'='*60}")
    print(f"  KPI 1 — VOF mass drift        : {drift:.4f}%       [{drift_status}]  (< 0.5%)")
    print(f"  KPI 2 — Interface area ratio  : {iface_rat:.3f}         [{iface_status}]  (< 3.0)")
    print(f"  KPI 3 — Velocity RMS ratio    : {vel_rat:.3f}         [{vel_status}]  ([0.7, 1.5])")
    print(f"  KPI 4 — CFL estimate          : {cfl_max:.3f}         [{cfl_status}]  (< 0.6)")
    print(f"  KPI 5 — KE quasi-steady ratio : {ke_rat:.3f}         [{ke_status}]  ([0.7, 1.5])")
    iter_str = f"{int(max_iter):4d}     " if not math.isnan(max_iter) else "  n/a  "
    print(f"  KPI 6 — Poisson max iters     : {iter_str}     [{iter_status}]  (< {_NITERMAX}; needs -health build)")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(description="Numerical health report for BioReactor runs")
    parser.add_argument("run_dirs", nargs="+", help="One or more run directories")
    args = parser.parse_args()

    for d in args.run_dirs:
        report(Path(d))


if __name__ == "__main__":
    main()
