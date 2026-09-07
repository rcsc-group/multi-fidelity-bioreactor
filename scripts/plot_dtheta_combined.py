"""Combined settling-time-vs-Delta_theta_max figure (diary.md 2026-09-07):
both directions on one axis -- target=theta=7deg (sources approaching from
below) and target=theta=4deg (sources approaching from above) -- plus the
finer-resolution Test A point (Delta_theta=0.05, target=4.05) that first
showed the dependence is steep-but-continuous, not a step function.

Usage:
    uv run python scripts/plot_dtheta_combined.py
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


rows_inc, cold_inc = load_series(
    ROOT / "experiments" / "dtheta_monotonicity_manifest.json",
    "settling_baseline_L6", 32.5, 7.0,
)
rows_dec, cold_dec = load_series(
    ROOT / "experiments" / "dtheta_monotonicity_th4target_manifest.json",
    "settling_ref_L6_rpm32.5_th4", 32.5, 4.0,
)

fig, ax = plt.subplots(figsize=(6.5, 4.5))

d_inc = [r["d_theta"] for r in rows_inc if r["t_settle"] is not None]
s_inc = [r["t_settle"] for r in rows_inc if r["t_settle"] is not None]
never_inc = [r["d_theta"] for r in rows_inc if r["t_settle"] is None]

d_dec = [r["d_theta"] for r in rows_dec if r["t_settle"] is not None]
s_dec = [r["t_settle"] for r in rows_dec if r["t_settle"] is not None]
never_dec = [r["d_theta"] for r in rows_dec if r["t_settle"] is None]

y_never = 65

ax.scatter(d_inc, s_inc, s=60, marker="o", color="#2166ac", zorder=3,
           label=r"target $\theta_\mathrm{max}=7°$")
ax.scatter(d_dec, s_dec, s=60, marker="s", color="#b2182b", zorder=3,
           label=r"target $\theta_\mathrm{max}=4°$")
if never_inc:
    ax.scatter(never_inc, [y_never] * len(never_inc), s=70, marker="^",
               facecolors="none", edgecolors="#2166ac", zorder=3)
if never_dec:
    ax.scatter(never_dec, [y_never] * len(never_dec), s=70, marker="^",
               facecolors="none", edgecolors="#b2182b", zorder=3)
ax.axhline(cold_inc, color="#555555", linestyle="--", linewidth=1.2, label="cold start")

ax.set_xlabel(r"$\Delta\theta_\mathrm{max}$ (deg)")
ax.set_ylabel("Settling time (cycles)")
ax.set_ylim(-3, y_never + 5)
ax.legend(frameon=False, loc="upper right")
fig.tight_layout()

out_dir = ROOT / "experiments" / "dtheta_monotonicity"
out_path = out_dir / "settling_vs_dtheta_combined.png"
fig.savefig(out_path, dpi=150)
print(f"Saved {out_path}")
print(f"cold start settling (both targets): inc={cold_inc}, dec={cold_dec}")
