"""Overlay Figure 13(a) from Kim et al. (2024): shear stress vs RPM.

Kim's figure shows:
  - Solid blue circles  : tau_liq_max  (absolute spatial max over one period)
  - Hollow blue circles : tau_liq_mean (max of spatially-averaged tau over one period)

We overlay our L8 theta sweep (theta ∈ {2,3,4,5,6,7}°, RPM ∈ {15..37.5}):
  - Upward triangles, filled
  - Marker size proportional to theta_max
  - tau_100_max on the y-axis  (direct analog of Kim's tau_liq_max)

...and our L9 tau/RPM sweep at fixed theta=7° (hydrodynamics-only runs, no
oxygen injection -- t_end is far shorter than t_mix by design, so only the
shear-stress KPIs are meaningful here):
  - Black squares (filled: tau_100_max, hollow: tau_mean_max), same dual
    solid/hollow convention as Kim's own two series.

...and our L10 (fidelity 10) tau/RPM sweep, same theta=7° conditions --
PARTIAL DATA: only seg0 of each 3-4 segment chain has completed as of this
writing, so these are transient (not yet quasi-steady-periodic) values and
should be read as a preliminary sanity check, not a validated comparison:
  - Gray diamonds (filled: tau_100_max, hollow: tau_mean_max)

Both L9 and L10 tau_100_max/tau_mean_max are computed in postprocess.py's
_compute_tau98_kpis(), restricted to the QSS window (t_ramp < t < t_inject)
to match Kim's "over one period" methodology (Main.tex Sec. 4.5): tau_100_max
is the absolute max over space+time (solid-circle analog), tau_mean_max is
the max over time of the spatially-averaged tau (hollow-circle analog).

Usage:
    uv run python scripts/plot_kim_overlay_tau.py
"""

import json
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────
HERE       = Path(__file__).parent.parent          # dev/rocking-bioreactor-2d/
KIM_CSV    = HERE / "docs/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv"
EXP_DIR    = HERE / "experiments/sweep_fb_theta_l8_mpi_ckpt/experiment_data"
L9_META    = HERE / "experiments/sweep_tau_theta7_l9/_sweep_metadata.json"
L10_MANIFEST = HERE / "experiments/sweep_tau_theta7_l10/_chain_manifest.json"
SCRATCH    = Path("/oscar/scratch/eaguerov/mpi_runs")
OUT_DIR    = HERE / "experiments/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── load Kim data ──────────────────────────────────────────────────────────────
# CSV has two header rows; drop the second (units) row
kim_raw = pd.read_csv(KIM_CSV, header=0)
kim = kim_raw[kim_raw["RPM"].apply(lambda x: str(x).replace('.','').isdigit()
              or (isinstance(x, float) and not np.isnan(x)))].copy()
kim["RPM"]          = pd.to_numeric(kim["RPM"])
kim["tau_liq_max"]  = pd.to_numeric(kim["tau_liq_max"])
kim["tau_liq_mean"] = pd.to_numeric(kim["tau_liq_mean"])
kim = kim.sort_values("RPM")

# ── load our data ──────────────────────────────────────────────────────────────
inp = pd.read_csv(EXP_DIR / "input.csv").reset_index(drop=True)
out = pd.read_csv(EXP_DIR / "output.csv", index_col=0).reset_index(drop=True)

ours = pd.DataFrame({
    "rpm":       inp["omega_b"] * 60.0 / (2 * math.pi),
    "theta":     inp["theta_max_0"],
    "tau_max":   out["tau_100_max"],
}).dropna(subset=["tau_max"])

THETAS = sorted(ours["theta"].unique())   # [2, 3, 4, 5, 6, 7]

# ── load L9 tau/RPM sweep (theta=7 fixed, scratch-resident results) ────────────
def _load_tau_sweep(rpms, run_ids):
    rows = []
    for rpm, run_id in zip(rpms, run_ids):
        results_path = SCRATCH / run_id / "results.json"
        if not results_path.exists():
            continue
        r = json.loads(results_path.read_text())
        tau_max, tau_mean_max = r.get("tau_100_max"), r.get("tau_mean_max")
        if tau_max is None or not math.isfinite(tau_max):
            continue
        rows.append({
            "rpm": rpm,
            "tau_max": tau_max,
            "tau_mean_max": tau_mean_max if tau_mean_max and math.isfinite(tau_mean_max) else None,
        })
    return pd.DataFrame(rows)

l9_meta = json.loads(L9_META.read_text())
l9 = _load_tau_sweep(l9_meta["rpms"], l9_meta["run_ids"])

# ── load L10 tau/RPM sweep -- latest completed segment in each chain; as of
# this writing that's still only seg0 for most conditions, so these are
# PARTIAL/transient values, not a converged quasi-steady comparison.
def _latest_completed_segment(chain):
    for run_id in reversed(chain):
        if (SCRATCH / run_id / "results.json").exists():
            return run_id
    return chain[0]

l10_manifest = json.loads(L10_MANIFEST.read_text())
l10 = _load_tau_sweep(
    [float(rpm) for rpm in l10_manifest.keys()],
    [_latest_completed_segment(chain) for chain in l10_manifest.values()],
)

# Red colormap: darker shade for larger theta
_cmap = plt.colormaps["Reds"]
_n    = len(THETAS)
# sample 0.40-0.90 so the lightest shade is still visible
c_map = {th: _cmap(0.40 + 0.50 * i / (_n - 1)) for i, th in enumerate(THETAS)}

# Marker sizes: linear in theta area
S_MIN, S_MAX = 30, 200
s_map = {th: S_MIN + (S_MAX - S_MIN) * (th - THETAS[0]) / (THETAS[-1] - THETAS[0])
         for th in THETAS}

# ── figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5.5, 4.0))

# Kim — solid circles (tau_liq_max)
ax.plot(kim["RPM"], kim["tau_liq_max"],
        color="royalblue", marker="o", ms=7, lw=1.2,
        label=r"Kim $\tau_\mathrm{max}$ (2024)", zorder=3)

# Kim — hollow circles (tau_liq_mean = max of <tau>)
ax.plot(kim["RPM"], kim["tau_liq_mean"],
        color="royalblue", marker="o", ms=7, lw=1.2, ls="--",
        markerfacecolor="white", markeredgecolor="royalblue",
        label=r"Kim $\langle\tau\rangle_\mathrm{max}$ (2024)", zorder=3)

# Ours — triangles: color = red scale by theta, size proportional to theta
for _, row in ours.iterrows():
    ax.scatter(row["rpm"], row["tau_max"],
               marker="^", s=s_map[row["theta"]],
               color=c_map[row["theta"]], alpha=0.70,
               edgecolors="none", zorder=4)

# Ours — L9 tau/RPM sweep at theta=7°: black squares (filled=tau_100_max, hollow=tau_mean_max)
ax.plot(l9["rpm"], l9["tau_max"],
        color="black", marker="s", ms=6, lw=1.0, ls="-",
        markerfacecolor="black", alpha=0.85, zorder=4)
l9_mean = l9.dropna(subset=["tau_mean_max"])
ax.plot(l9_mean["rpm"], l9_mean["tau_mean_max"],
        color="black", marker="s", ms=6, lw=1.0, ls="--",
        markerfacecolor="white", markeredgecolor="black", alpha=0.85, zorder=4)

# Ours — L10 tau/RPM sweep at theta=7°, PARTIAL (mostly seg0 only, transient):
# gray diamonds (filled=tau_100_max, hollow=tau_mean_max)
ax.plot(l10["rpm"], l10["tau_max"],
        color="dimgray", marker="D", ms=6, lw=1.0, ls="-",
        markerfacecolor="dimgray", alpha=0.85, zorder=5)
l10_mean = l10.dropna(subset=["tau_mean_max"])
ax.plot(l10_mean["rpm"], l10_mean["tau_mean_max"],
        color="dimgray", marker="D", ms=6, lw=1.0, ls="--",
        markerfacecolor="white", markeredgecolor="dimgray", alpha=0.85, zorder=5)

# ── legend: Kim series + size guide for theta ─────────────────────────────────
kim_max_handle  = mlines.Line2D([], [], color="royalblue", marker="o",
                                ms=7, lw=1.2,
                                label=r"Kim $\tau_\mathrm{max}$")
kim_mean_handle = mlines.Line2D([], [], color="royalblue", marker="o",
                                ms=7, lw=1.2, ls="--",
                                markerfacecolor="white",
                                markeredgecolor="royalblue",
                                label=r"Kim $\langle\tau\rangle_\mathrm{max}$")

theta_handles = [
    mlines.Line2D([], [], color=c_map[th], marker="^", lw=0,
                  ms=math.sqrt(s_map[th]),
                  label=rf"$\theta={int(th)}°$")
    for th in THETAS
]

l9_max_handle  = mlines.Line2D([], [], color="black", marker="s", ms=6, lw=1.0,
                                label=r"L9 $\tau_\mathrm{max}$")
l9_mean_handle = mlines.Line2D([], [], color="black", marker="s", ms=6, lw=1.0, ls="--",
                                markerfacecolor="white", markeredgecolor="black",
                                label=r"L9 $\langle\tau\rangle_\mathrm{max}$")
l10_max_handle  = mlines.Line2D([], [], color="dimgray", marker="D", ms=6, lw=1.0,
                                 label=r"L10 $\tau_\mathrm{max}$ (partial)")
l10_mean_handle = mlines.Line2D([], [], color="dimgray", marker="D", ms=6, lw=1.0, ls="--",
                                 markerfacecolor="white", markeredgecolor="dimgray",
                                 label=r"L10 $\langle\tau\rangle_\mathrm{max}$ (partial)")

ax.legend(handles=[kim_max_handle, kim_mean_handle] + theta_handles
                   + [l9_max_handle, l9_mean_handle, l10_max_handle, l10_mean_handle],
          fontsize=6.8, framealpha=0.9, loc="lower left",
          handlelength=1.4, handletextpad=0.5)

# ── axes ──────────────────────────────────────────────────────────────────────
ax.set_xlabel(r"Rocking frequency $f_b$ (rpm)", fontsize=11)
ax.set_ylabel(r"Shear stress (Pa)", fontsize=11)
ax.set_xlim(13, 39)
ax.set_yscale("log")
ax.set_ylim(5e-4, 3.0)
ax.tick_params(which="both", direction="in", top=True, right=True)
ax.grid(True, which="major", ls=":", alpha=0.4)

fig.tight_layout()

out_path = OUT_DIR / "overlay_tau_rpm.pdf"
fig.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"Saved: {out_path}")

out_path_png = OUT_DIR / "overlay_tau_rpm.png"
fig.savefig(out_path_png, dpi=150, bbox_inches="tight")
print(f"Saved: {out_path_png}")
