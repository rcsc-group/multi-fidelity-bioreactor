"""Regenerate docs/kimetal2024/figure_replicas/replicated_Fig13.png from
currently-valid data only (diary.md 2026-09-02).

The version pushed 2026-08-05 (commit 7f87e46) baked in:
  - L6/L8 sweeps run before the H_bio nondim factor-of-2 fix (052e9e4,
    2026-08-20) and the tau-histogram OpenMP data-race fix (a648ca2,
    2026-08-21) -- both change tau/EDR postprocessing outputs.
  - A single L10 point recovered from `l10_kim_seg2`, one chain segment
    explicitly documented (2026-08-05 diary entry, plot_kim_overlay_tau.py
    docstring) as PARTIAL/transient, not a converged quasi-steady value --
    exactly the restart-transient-biased class of run this project later
    characterized as unreliable.

This script plots only what is currently valid:
  - Kim et al.'s published curve (unchanged).
  - Our L8 fig13a_rampmatch sweep (2026-09-01): current bug-fixed driver,
    upstream's own ramp matched, 9 independent cold starts. Supersedes the
    old L8 series outright, zero new compute.

L6 and L10 are omitted rather than shown stale -- neither has a rerun on
the current driver yet. See diary.md 2026-09-02.
"""
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).parent.parent
KIM_CSV = HERE / "docs/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv"
RUNS_DIR = HERE / "runs"
OUT_PATH = HERE / "docs/kimetal2024/figure_replicas/replicated_Fig13.png"

RPMS = [17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5]

kim = pd.read_csv(KIM_CSV, skiprows=[1])  # row 1 is a units/label row, not data
kim["RPM"] = pd.to_numeric(kim["RPM"])
kim = kim.sort_values("RPM")

rows = []
for rpm in RPMS:
    d = json.loads((RUNS_DIR / f"fig13a_rampmatch_rpm{rpm:g}" / "results.json").read_text())
    rows.append({"rpm": rpm, "tau_max": d["tau_100_max"], "tau_mean_max": d["tau_mean_max"]})
ours = pd.DataFrame(rows).sort_values("rpm")

fig, ax = plt.subplots(figsize=(6, 4.3))

ax.plot(kim["RPM"], kim["tau_liq_max"], color="royalblue", marker="o", ms=7, lw=1.3,
        label=r"Kim $\tau_\mathrm{max}$")
ax.plot(kim["RPM"], kim["tau_liq_mean"], color="royalblue", marker="o", ms=7, lw=1.3, ls="--",
        markerfacecolor="white", label=r"Kim $\langle\tau\rangle_\mathrm{max}$")

ax.plot(ours["rpm"], ours["tau_max"], color="darkred", marker="^", ms=8, lw=1.3,
        label=r"Our fork $\tau_\mathrm{max}$ (L8)")
ax.plot(ours["rpm"], ours["tau_mean_max"], color="darkred", marker="^", ms=8, lw=1.3, ls="--",
        markerfacecolor="white", label=r"Our fork $\langle\tau\rangle_\mathrm{max}$ (L8)")

ax.set_xlabel(r"Rocking frequency $f_b$ (rpm)", fontsize=11)
ax.set_ylabel("Shear stress (Pa)", fontsize=11)
ax.set_yscale("log")
ax.set_title(r"$\theta=7°$", fontsize=10)
ax.tick_params(which="both", direction="in", top=True, right=True)
ax.grid(True, which="major", ls=":", alpha=0.4)
ax.legend(fontsize=8, framealpha=0.9, loc="lower left")
fig.tight_layout()
fig.savefig(OUT_PATH, dpi=150)
print(f"saved {OUT_PATH}")
