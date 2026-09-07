"""Delta-theta_max monotonicity study, decreasing-direction mirror v2,
target=theta=4deg (diary.md 2026-09-07). Same noise-floor-aware
methodology as analyze_dtheta_monotonicity_decreasing.py (calibrate
tolerance from the trivial Delta_theta=0 identity-restart control's own
full-length deviation, not a short reference tail).

Usage:
    uv run python scripts/analyze_dtheta_th4_target.py
"""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
TOL = 0.05
MANIFEST = json.loads((ROOT / "experiments" / "dtheta_monotonicity_th4target_manifest.json").read_text())
RPM = MANIFEST["rpm"]
TARGET_THETA = MANIFEST["target_theta"]


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    th = math.radians(theta)
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th))
    U_bio = V_bio / (H_bio * 0.5) / T_per
    T_bio = L_bio / U_bio
    return T_per / T_bio


T_ND = t_per_nd(RPM, TARGET_THETA)


def per_cycle_peak(t, tau_mean, t0, n_cycles):
    cycles = (t - t0) / T_ND
    out = np.full(n_cycles, np.nan)
    for c in range(n_cycles):
        mask = (cycles >= c) & (cycles < c + 1)
        if mask.sum() > 0:
            out[c] = tau_mean[mask].max()
    return out


def find_settle(peaks, target, tol):
    reldiff = np.abs(peaks - target) / abs(target)
    for start in range(len(reldiff)):
        if np.all(reldiff[start:] < tol):
            return start
    return None


baseline = np.loadtxt(RUNS / "settling_ref_L6_rpm32.5_th4" / "shear_stress.dat", skiprows=1)
n_cycles_base = int(baseline[-1, 1] / T_ND)
base_peaks = per_cycle_peak(baseline[:, 1], baseline[:, 5], 0.0, n_cycles_base)
target_steady = float(np.nanmean(base_peaks[-5:]))

_trivial_run_id = MANIFEST["edges"]["4"]["run_id"]
_trivial_params = json.loads((RUNS / _trivial_run_id / "params.json").read_text())
_t_ckpt0 = _trivial_params["t_checkpoint"]
_trivial = np.loadtxt(RUNS / _trivial_run_id / "shear_stress.dat", skiprows=1)
_n0 = int((_trivial[-1, 1] - _t_ckpt0) / T_ND)
_trivial_peaks = per_cycle_peak(_trivial[:, 1], _trivial[:, 5], _t_ckpt0, _n0)
_trivial_mean = float(np.nanmean(_trivial_peaks))
_max_reldev = float(np.nanmax(np.abs(_trivial_peaks - _trivial_mean))) / _trivial_mean
tol_eff = max(TOL, 1.2 * _max_reldev)
offset = abs(_trivial_mean - target_steady) / target_steady
print(f"trivial-run mean={_trivial_mean:.4e}  reference tail mean={target_steady:.4e}  "
      f"offset={100*offset:.1f}%")
print(f"trivial-run noise calibration: max single-cycle deviation={100*_max_reldev:.1f}%  "
      f"-> effective tolerance={100*tol_eff:.1f}% (flat TOL={100*TOL:.0f}%)")

cold_settle = find_settle(base_peaks, target_steady, tol_eff)

rows = []
for key, e in MANIFEST["edges"].items():
    theta_source = e["theta_source"]
    run_id = e["run_id"]
    d_theta = theta_source - TARGET_THETA

    params = json.loads((RUNS / run_id / "params.json").read_text())
    t_ckpt = params["t_checkpoint"]
    tr = np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)
    n_cycles_tr = int((tr[-1, 1] - t_ckpt) / T_ND)
    tr_peaks = per_cycle_peak(tr[:, 1], tr[:, 5], t_ckpt, n_cycles_tr)

    settle = find_settle(tr_peaks, target_steady, tol_eff)
    rows.append({"theta_source": theta_source, "d_theta": d_theta, "t_settle": settle,
                 "n_cycles_available": n_cycles_tr})
    settle_str = str(settle) if settle is not None else "NEVER"
    print(f"theta_source={theta_source:>4g}  d_theta={d_theta:>4.1f}  t_settle={settle_str:>6}  "
          f"(of {n_cycles_tr} recorded)")

rows.sort(key=lambda r: r["d_theta"])
print(f"\ncold start's own settling time (self-referential): {cold_settle}")

resolved = [r for r in rows if r["t_settle"] is not None]
d_thetas = [r["d_theta"] for r in resolved]
settles = [r["t_settle"] for r in resolved]
is_monotonic = all(settles[i] <= settles[i + 1] for i in range(len(settles) - 1))
zero_at_zero = next((r["t_settle"] for r in rows if r["d_theta"] == 0.0), "MISSING")

print(f"\n(1) monotonic non-decreasing in d_theta: {is_monotonic}")
print(f"    settling times (sorted by d_theta): {settles}")
print(f"(2) settling time at d_theta=0: {zero_at_zero}")

fig, ax = plt.subplots(figsize=(5.5, 4))
ax.scatter(d_thetas, settles, s=50, color="#2166ac", zorder=3, label="checkpoint restart")
never = [r["d_theta"] for r in rows if r["t_settle"] is None]
if never:
    ax.scatter(never, [max(settles, default=0) * 1.1] * len(never), marker="^",
               color="#b2182b", zorder=3, label="not settled in window")
if cold_settle is not None:
    ax.axhline(cold_settle, color="#555555", linestyle="--", linewidth=1.2, label="cold start")

ax.set_xlabel(r"$\Delta\theta_\mathrm{max}$ (deg)")
ax.set_ylabel("Settling time (cycles)")
ax.legend(frameon=False)
fig.tight_layout()

out_dir = ROOT / "experiments" / "dtheta_monotonicity"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "settling_vs_dtheta_th4target.png"
fig.savefig(out_path, dpi=150)
print(f"\nSaved {out_path}")
