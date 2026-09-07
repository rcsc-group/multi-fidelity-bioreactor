"""L6 vs L7 resolution-independence check for the Delta_theta_max
settling-time pattern (diary.md 2026-09-07): does the near-flat/steep-onset
shape found at L6 replicate at L7 (real phenomenon), or shift/soften
(discretization artifact)?

Usage:
    uv run python scripts/plot_dtheta_L6_vs_L7.py
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


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    th = math.radians(theta)
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th))
    U_bio = V_bio / (H_bio * 0.5) / T_per
    T_bio = L_bio / U_bio
    return T_per / T_bio


def per_cycle_peak(t, tau_mean, t0, T_nd, n_cycles):
    cycles = (t - t0) / T_nd
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


def load_series(manifest_path, baseline_run_id, rpm, target_theta, tol=TOL):
    manifest = json.loads(manifest_path.read_text())
    T_nd = t_per_nd(rpm, target_theta)
    baseline = np.loadtxt(RUNS / baseline_run_id / "shear_stress.dat", skiprows=1)
    n_base = int(baseline[-1, 1] / T_nd)
    base_peaks = per_cycle_peak(baseline[:, 1], baseline[:, 5], 0.0, T_nd, n_base)
    target_steady = float(np.nanmean(base_peaks[-5:]))
    cold_settle = find_settle(base_peaks, target_steady, tol)

    rows = []
    for key, e in manifest["edges"].items():
        theta_source = e["theta_source"]
        run_id = e["run_id"]
        d_theta = abs(target_theta - theta_source)
        params = json.loads((RUNS / run_id / "params.json").read_text())
        t_ckpt = params["t_checkpoint"]
        tr = np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)
        n_cycles_tr = int((tr[-1, 1] - t_ckpt) / T_nd)
        tr_peaks = per_cycle_peak(tr[:, 1], tr[:, 5], t_ckpt, T_nd, n_cycles_tr)
        settle = find_settle(tr_peaks, target_steady, tol)
        rows.append({"d_theta": d_theta, "t_settle": settle})
    return rows, cold_settle


rows_l6, cold_l6 = load_series(
    ROOT / "experiments" / "dtheta_monotonicity_manifest.json",
    "settling_baseline_L6", 32.5, 7.0,
)
rows_l7, cold_l7 = load_series(
    ROOT / "experiments" / "dtheta_L7_check_manifest.json",
    "settling_baseline_L7", 32.5, 7.0,
)

print("L6:", [(r["d_theta"], r["t_settle"]) for r in rows_l6])
print("L7:", [(r["d_theta"], r["t_settle"]) for r in rows_l7])
print("cold start: L6=", cold_l6, " L7=", cold_l7)

fig, ax = plt.subplots(figsize=(6.5, 4.5))

d_l6 = [r["d_theta"] for r in rows_l6 if r["t_settle"] is not None]
s_l6 = [r["t_settle"] for r in rows_l6 if r["t_settle"] is not None]
d_l7 = [r["d_theta"] for r in rows_l7 if r["t_settle"] is not None]
s_l7 = [r["t_settle"] for r in rows_l7 if r["t_settle"] is not None]
never_l7 = [r["d_theta"] for r in rows_l7 if r["t_settle"] is None]

y_never = 65

ax.scatter(d_l6, s_l6, s=70, marker="o", color="#2166ac", zorder=3, label="L6, settled")
ax.scatter(d_l7, s_l7, s=110, marker="D", facecolors="none", edgecolors="#b2182b",
           linewidths=2.2, zorder=4, label="L7, settled")
if never_l7:
    ax.scatter(never_l7, [y_never] * len(never_l7), s=110, marker="^",
               facecolors="none", edgecolors="#b2182b", linewidths=2.2, zorder=4,
               label="L7, not settled in window")
ax.axhline(cold_l6, color="#555555", linestyle="--", linewidth=1.2, label="cold start (L6)")

ax.set_xlabel(r"$\Delta\theta_\mathrm{max}$ (deg)")
ax.set_ylabel("Settling time (cycles)")
ax.set_ylim(-3, y_never + 5)
ax.legend(frameon=False, loc="center right")
fig.tight_layout()

out_path = ROOT / "experiments" / "dtheta_monotonicity" / "settling_vs_dtheta_L6_vs_L7.png"
fig.savefig(out_path, dpi=150)
print(f"\nSaved {out_path}")
