"""Delta-theta_max monotonicity study, analysis + plot (diary.md 2026-09-07).

Tests, for a single fixed target condition (32.5rpm, theta_max=7deg, L6),
whether checkpoint-restart settling time behaves sanely as a function of
the SIZE of the warm-start perturbation:
  (1) settling time should decrease monotonically as Delta_theta_max -> 0
  (2) settling time should be exactly 0 at Delta_theta_max = 0 (source ==
      target, "already converged")
Cold start's own settling time (self-referential: its own early cycles vs
its own converged last-5-cycle value) is plotted as a reference line for
context -- it answers a different question (distance from REST, not from a
nearby converged state) so isn't on the same x-axis, just a fixed y-value.

Usage:
    uv run python scripts/analyze_dtheta_monotonicity.py
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
MANIFEST = json.loads((ROOT / "experiments" / "dtheta_monotonicity_manifest.json").read_text())
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


def find_settle(peaks, target, tol=TOL):
    reldiff = np.abs(peaks - target) / abs(target)
    for start in range(len(reldiff)):
        if np.all(reldiff[start:] < tol):
            return start
    return None


# Target's own converged steady value, from its cold start.
baseline = np.loadtxt(RUNS / "settling_baseline_L6" / "shear_stress.dat", skiprows=1)
n_cycles_base = int(baseline[-1, 1] / T_ND)
base_peaks = per_cycle_peak(baseline[:, 1], baseline[:, 5], 0.0, n_cycles_base)
target_steady = float(np.nanmean(base_peaks[-5:]))

# Cold start's own settling time (self-referential), for the reference line.
cold_settle = find_settle(base_peaks, target_steady)

rows = []
for key, e in MANIFEST["edges"].items():
    theta_source = e["theta_source"]
    run_id = e["run_id"]
    d_theta = TARGET_THETA - theta_source

    params = json.loads((RUNS / run_id / "params.json").read_text())
    t_ckpt = params["t_checkpoint"]
    tr = np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)
    n_cycles_tr = int((tr[-1, 1] - t_ckpt) / T_ND)
    tr_peaks = per_cycle_peak(tr[:, 1], tr[:, 5], t_ckpt, n_cycles_tr)

    settle = find_settle(tr_peaks, target_steady)
    rows.append({"theta_source": theta_source, "d_theta": d_theta, "t_settle": settle,
                 "n_cycles_available": n_cycles_tr})
    settle_str = str(settle) if settle is not None else "NEVER"
    print(f"theta_source={theta_source:>4g}  d_theta={d_theta:>4.1f}  t_settle={settle_str:>6}  "
          f"(of {n_cycles_tr} recorded)")

rows.sort(key=lambda r: r["d_theta"])
print(f"\ncold start's own settling time (self-referential): {cold_settle}")

# ── checks ───────────────────────────────────────────────────────────────
resolved = [r for r in rows if r["t_settle"] is not None]
d_thetas = [r["d_theta"] for r in resolved]
settles = [r["t_settle"] for r in resolved]

is_monotonic = all(settles[i] <= settles[i + 1] for i in range(len(settles) - 1))
zero_at_zero = next((r["t_settle"] for r in rows if r["d_theta"] == 0.0), "MISSING")

print(f"\n(1) monotonic non-decreasing in d_theta: {is_monotonic}")
print(f"    settling times (sorted by d_theta): {settles}")
print(f"(2) settling time at d_theta=0: {zero_at_zero}")

# ── plot ─────────────────────────────────────────────────────────────────
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
out_path = out_dir / "settling_vs_dtheta.png"
fig.savefig(out_path, dpi=150)
print(f"\nSaved {out_path}")
