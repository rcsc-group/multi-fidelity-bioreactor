"""Regenerate experiments/kimetal2024/figure_replicas/replicated_Fig13.png from
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
  - Our L6 fig13a_l6 sweep (2026-09-02): same driver/ramp/methodology,
    fidelity=6 instead of 8.
  - [ADDED 2026-09-10] Our L9 CHAINED sweep: 17.5rpm independent cold
    start, then a checkpoint chain 20->22.5->...->37.5rpm, each hop
    warm-started from the previous rpm's own converged state -- the
    production chained-sweep protocol this week's phase-offset/su/g fixes
    targeted, not independent cold starts like L6/L8. Every restart
    segment passed scripts/validate_run.py (bit-exact restore against its
    parent, correct parent linkage, settled tail) before being trusted
    here; the one exception is the 17.5rpm source's own restore fidelity,
    unverifiable rather than confirmed-bad (its checkpoint predates the
    diagnostic instrumentation -- see diary.md 2026-09-08/09). Also the
    run that caught and required fixing a real bug in postprocess.py's
    QSS-window logic (_ramp_end_nd ignored t_checkpoint, so a restart
    segment's own post-restart transient was silently included in what
    was supposed to be its converged window) -- fixed before these numbers
    were generated, not after.

L10 is omitted rather than shown stale -- no rerun on the current driver
yet. See diary.md 2026-09-02.
"""
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).parent.parent
KIM_CSV = HERE / "experiments/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv"
RUNS_DIR = HERE / "runs"
OUT_PATH = HERE / "experiments/kimetal2024/figure_replicas/replicated_Fig13.png"

RPMS = [17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5]

kim = pd.read_csv(KIM_CSV, skiprows=[1])  # row 1 is a units/label row, not data
kim["RPM"] = pd.to_numeric(kim["RPM"])
kim = kim.sort_values("RPM")


def _ediss_peak(run_id: str) -> float:
    """Max over the settled window of the per-timestep domain-mean EDR.

    Mirrors how tau_mean_max is defined (max over time of a domain mean), so
    it is the analogue of Kim's Ediss_liq_mean. Kim's Ediss_liq_MAX is a
    SPATIAL max over the liquid, which this project does not log as a KPI --
    only ediss_field in the video frames carries it -- so that series is not
    plotted rather than substituted with a different quantity.
    """
    import numpy as np
    a = np.loadtxt(RUNS_DIR / run_id / "shear_stress.dat", skiprows=1)
    t = a[:, 1]
    tail = t >= t[0] + 0.5 * (t[-1] - t[0])   # settled half
    return float(np.max(a[tail, 9]))


def _load(run_id_fmt):
    rows = []
    for rpm in RPMS:
        d = json.loads((RUNS_DIR / run_id_fmt.format(rpm=rpm) / "results.json").read_text())
        rid = run_id_fmt.format(rpm=rpm)
        rows.append({"rpm": rpm, "tau_max": d["tau_100_max"],
                     "tau_mean_max": d["tau_mean_max"],
                     "ediss": _ediss_peak(rid)})
    return pd.DataFrame(rows).sort_values("rpm")


l8 = _load("fig13a_rampmatch_rpm{rpm:g}")
l6 = _load("fig13a_l6_rpm{rpm:g}")

L9_RUN_IDS = {
    17.5: "l9_sweep_rpm17.5", 20.0: "20a14369", 22.5: "255e0b87",
    25.0: "bce29aa5", 27.5: "9ec56180", 30.0: "60d09a80",
    32.5: "a34fc4d4", 35.0: "a281a16f", 37.5: "0f0ad3ea",
}
l9 = pd.DataFrame(
    [{"rpm": rpm,
      "tau_max": json.loads((RUNS_DIR / rid / "results.json").read_text())["tau_100_max"],
      "tau_mean_max": json.loads((RUNS_DIR / rid / "results.json").read_text())["tau_mean_max"],
      "ediss": _ediss_peak(rid)}
     for rpm, rid in L9_RUN_IDS.items()]
).sort_values("rpm")

# [ADDED 2026-09-15] L10, the two points that exist so far (32-rank runs).
# These were submitted as 5-cycle TIMING calibrations, and with a 3-cycle
# ramp they carry only ~2 post-ramp cycles -- per-cycle tau was still
# climbing at the final cycle (17.5rpm +12%, 37.5rpm +62%), i.e. neither has
# reached quasi-steady. They are plotted at the user's explicit request;
# the caveat lives here and in diary.md 2026-09-15, not in the figure (see
# the project's figure rules: no status/caveat text baked into the image).
# Replace with the converged L9->L10 warm-started sweep when it lands.
L10_RUN_IDS = {
    17.5: "l10_fig13a_rpm17.5_32rank_calib",
    37.5: "l10_fig13a_rpm37.5_32rank_calib",
}
l10 = pd.DataFrame(
    [{"rpm": rpm,
      "tau_max": json.loads((RUNS_DIR / rid / "results.json").read_text())["tau_100_max"],
      "tau_mean_max": json.loads((RUNS_DIR / rid / "results.json").read_text())["tau_mean_max"],
      "ediss": _ediss_peak(rid)}
     for rpm, rid in L10_RUN_IDS.items()]
).sort_values("rpm")

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 4.6))

ax.plot(kim["RPM"], kim["tau_liq_max"], color="royalblue", marker="o", ms=7, lw=1.3,
        label=r"Kim $\tau$")
ax.plot(kim["RPM"], kim["tau_liq_mean"], color="royalblue", marker="o", ms=7, lw=1.3, ls="--",
        markerfacecolor="white", label=r"Kim $\langle\tau\rangle$")

ax.plot(l8["rpm"], l8["tau_max"], color="darkred", marker="^", ms=8, lw=1.3,
        label=r"Ours $\tau$ (L8)")
ax.plot(l8["rpm"], l8["tau_mean_max"], color="darkred", marker="^", ms=8, lw=1.3, ls="--",
        markerfacecolor="white", label=r"Ours $\langle\tau\rangle$ (L8)")

ax.plot(l6["rpm"], l6["tau_max"], color="darkorange", marker="s", ms=7, lw=1.3,
        label=r"Ours $\tau$ (L6)")
ax.plot(l6["rpm"], l6["tau_mean_max"], color="darkorange", marker="s", ms=7, lw=1.3, ls="--",
        markerfacecolor="white", label=r"Ours $\langle\tau\rangle$ (L6)")

ax.plot(l9["rpm"], l9["tau_max"], color="seagreen", marker="D", ms=6, lw=1.3,
        label=r"Ours $\tau$ (L9)")
ax.plot(l9["rpm"], l9["tau_mean_max"], color="seagreen", marker="D", ms=6, lw=1.3, ls="--",
        markerfacecolor="white", label=r"Ours $\langle\tau\rangle$ (L9)")

ax.plot(l10["rpm"], l10["tau_max"], color="purple", marker="v", ms=7, lw=1.3,
        label=r"Ours $\tau$ (L10)")
ax.plot(l10["rpm"], l10["tau_mean_max"], color="purple", marker="v", ms=7, lw=1.3, ls="--",
        markerfacecolor="white", label=r"Ours $\langle\tau\rangle$ (L10)")

ax.set_xlabel(r"Rocking frequency $f_b$ (rpm)", fontsize=11)
ax.set_ylabel("Peak shear stress (Pa)", fontsize=11)
ax.set_yscale("log")
ax.tick_params(which="both", direction="in", top=True, right=True)
ax.grid(True, which="major", ls=":", alpha=0.4)
ax.legend(fontsize=8, framealpha=0.9, loc="upper left", bbox_to_anchor=(1.02, 1.0))

# ── (b) energy dissipation rate ──────────────────────────────────────────
ax2.plot(kim["RPM"], kim["Ediss_liq_mean"], color="royalblue", marker="o", ms=7, lw=1.3,
         label=r"Kim $\langle\epsilon\rangle$")
for df, col, mk, ms_, lbl in [
        (l6, "darkorange", "s", 7, "L6"), (l8, "darkred", "^", 8, "L8"),
        (l9, "seagreen", "D", 6, "L9"), (l10, "purple", "v", 7, "L10")]:
    ax2.plot(df["rpm"], df["ediss"], color=col, marker=mk, ms=ms_, lw=1.3,
             label=rf"Ours $\langle\epsilon\rangle$ ({lbl})")

ax2.set_xlabel(r"Rocking frequency $f_b$ (rpm)", fontsize=11)
ax2.set_ylabel(r"Energy dissipation rate (W/m$^3$)", fontsize=11)
ax2.set_yscale("log")
ax2.tick_params(which="both", direction="in", top=True, right=True)
ax2.grid(True, which="major", ls=":", alpha=0.4)
ax2.legend(fontsize=8, framealpha=0.9, loc="upper left", bbox_to_anchor=(1.02, 1.0))

fig.suptitle(r"$\theta=7°$", fontsize=11)
fig.tight_layout()
fig.savefig(OUT_PATH, dpi=150)
print(f"saved {OUT_PATH}")
