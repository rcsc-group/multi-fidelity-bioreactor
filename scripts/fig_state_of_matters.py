"""Overview figures of the scalar-QoI evidence (kLa, mixing time), plain axes: relative change in percent, shared y within a figure.

fig_state_levels   kLa_1T_10 and dtmix_0.50 vs rpm, theta 7, per level (L9 split by protocol).
fig_state_scatter  L8 0.625-rpm grids at theta 4, 7, 9: percent deviation of each run from the quadratic trend; outliers ringed.
fig_state_dose     change of the result (percent) for a tiny change of rpm, base 22.5 and 30 rpm.
fig_state_r12      change of the result (percent) when the tracer is released after 12 instead of 80 cycles.
fig_state_e9       difference (percent) between two runs 0.2 rpm apart, L8 and L9.
Usage: uv run python scripts/fig_state_of_matters.py
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts.outliers_l8 import zscores  # noqa: E402
from scripts.test_e8_theta_scatter import GRID, run_id  # noqa: E402

MF = ROOT / "experiments/multifidelity"
QL = {"kLa_1T_10": (fs.MK["cstar_10"], r"$k_La$, $C^*\!=\!0.10$"),
      "kLa_1T_25": (fs.MK["cstar_25"], r"$k_La$, $C^*\!=\!0.25$"),
      "dtmix_0.50": (fs.MK["chi_0.50"], r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$"),
      "dtmix_0.95": (fs.MK["chi_0.95"], r"$t_{\mathrm{mix}}$, $\chi\!=\!0.95$")}
YLAB = {"kLa_1T_10": r"$k_La$, $C^*\!=\!0.10$ (1/h)", "dtmix_0.50": r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$ (s)"}


def val(run, q):
    f = ROOT / "runs" / run / "results.json"
    if not f.exists():
        return math.nan
    v = json.loads(f.read_text()).get(q, math.nan)
    return float(v) if v is not None and float(v) > 0 else math.nan


def rpm_of(run):
    return json.loads((ROOT / "runs" / run / "params.json").read_text())["omega_b"] * 60 / (2 * math.pi)


def save(fig, name, extra=()):
    for ext in ("png", "pdf"):
        fig.savefig(MF / f"{name}.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=list(extra))
    print("wrote", MF / f"{name}.png")
    plt.close(fig)


def mkw(L, q):
    kw = fs.series_kw(fs.level_colour(L), QL[q][0], stat="mean", ours=True)
    kw.pop("ls")
    return kw


def levels():
    rows = [r for r in csv.DictReader((MF / "kla_mix_dataset.csv").open()) if float(r["theta"]) == 7.0]
    l8 = [(float(r["rpm"]), r["run"]) for r in rows if r["level"] == "8"]
    l9_old = [(float(r["rpm"]), r["run"]) for r in rows if r["level"] == "9"]
    e9 = [(rpm_of(r.name), r.name) for r in (ROOT / "runs").glob("e9_l9_rpm*_th7")]
    l10 = [(rpm_of(n), n) for n in ("l10c_rpm32.5_seg1", "v_l10_rpm32.5_th7")]
    l10_old = [(rpm_of(n), n) for n in ("fig11_l10_rpm35", "fig11_l10_rpm37.5")]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.9), sharex=True)
    for ax, q in zip(axes, ("kLa_1T_10", "dtmix_0.50")):
        for L, items, a in ((8, l8, 1.0), (9, l9_old, 0.4), (9, e9, 1.0), (10, l10, 1.0), (10, l10_old, 0.35)):
            xy = np.array([(x, val(r, q)) for x, r in items if np.isfinite(val(r, q))])
            if len(xy):
                ax.plot(xy[:, 0], xy[:, 1], ls="none", alpha=a, **mkw(L, q))
        ax.set_yscale("log")
        ax.set_ylabel(YLAB[q])
        ax.set_xlabel("rocking speed (rpm)")
        ax.grid(**fs.GRID_KW)
    h = [Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(8), ms=5, label="L8"),
         Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(9), ms=5, label="L9, same protocol as L8"),
         Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(9), alpha=0.4, ms=5, label="L9, cold start, other binary"),
         Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(10), ms=5, label="L10, 32.5 rpm"),
         Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(10), alpha=0.35, ms=5, label="L10, earlier runs")]
    lg = axes[1].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    save(fig, "fig_state_levels", [lg])


def scatter():
    ids = json.loads((MF / "test_t5_roughness.json").read_text())["run_ids"]
    fig, axes = plt.subplots(4, 3, figsize=(10.0, 9.0), sharex=True, sharey=True)
    for j, th in enumerate((4.0, 7.0, 9.0)):
        runs = [ids[f"{r:g}"] for r in GRID] if th == 7.0 else [run_id(r, th) for r in GRID]
        for i, q in enumerate(QL):
            ax = axes[i, j]
            y = np.array([val(rn, q) for rn in runs])
            x = np.array(GRID)
            ok = np.isfinite(y)
            c = np.polyfit(x[ok], np.log(y[ok]), 2)
            dev = 100 * (y / np.exp(np.polyval(c, x)) - 1)
            ax.plot(x[ok], dev[ok], ls="none", **mkw(8, q))
            z, _ = zscores(y)
            out = ok & (np.abs(z) > 3)
            ax.plot(x[out], dev[out], ls="none", marker="o", ms=12, mfc="none", mec="0.35", mew=1.0)
            ax.axhline(0, color=fs.GUIDE, lw=1.0)
            ax.grid(**fs.GRID_KW)
            if j == 0:
                ax.set_ylabel(QL[q][1] + "\n" + "deviation from trend (%)", fontsize=9)
            if i == 0:
                ax.set_title(rf"$\theta$ = {th:g}$^\circ$", fontsize=10)
            if i == 3:
                ax.set_xlabel("rocking speed (rpm)")
    ax.set_ylim(-60, 130)
    fig.tight_layout()
    save(fig, "fig_state_scatter")


def dose():
    ids = json.loads((MF / "test_t5_roughness.json").read_text())["run_ids"]
    eps = (1e-4, 1e-3, 3e-3, 1e-2, 3e-2)
    pairs = {22.5: ("t5_l8_rpm22.4_th7", "t5_l8_rpm22.6_th7"), 30.0: ("t5_l8_rpm29.9_th7", "t5_l8_rpm30.1_th7")}
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.0), sharex=True, sharey=True)
    for i, base in enumerate((22.5, 30.0)):
        for j, q in enumerate(("dtmix_0.50", "dtmix_0.95", "kLa_1T_10")):
            ax = axes[i, j]
            q0 = val(ids[f"{base:g}"], q)
            xs = [-0.1] + list(eps) + [0.1]
            ys = [100 * (val(pairs[base][0], q) / q0 - 1)] + [100 * (val(f"dose_l8_rpm{base + e:g}_th7", q) / q0 - 1) for e in eps] \
                + [100 * (val(pairs[base][1], q) / q0 - 1)]
            kw = mkw(8, q)
            ax.plot(xs, ys, ls=":", lw=0.8, color=fs.level_colour(8), marker=kw["marker"], mfc="w", mec=fs.level_colour(8), ms=6, mew=1.1)
            ax.axhline(0, color=fs.GUIDE, lw=1.0)
            ax.set_xscale("symlog", linthresh=1e-4, linscale=0.4)
            ax.grid(**fs.GRID_KW)
            if i == 0:
                ax.set_title(QL[q][1], fontsize=10)
            if i == 1:
                ax.set_xlabel("change of rocking speed (rpm)")
            if j == 0:
                ax.set_ylabel(f"base {base:g} rpm\nchange of result (%)", fontsize=9)
    axes[0, 0].set_ylim(-50, 60)
    axes[0, 0].set_xticks([-0.1, 1e-4, 1e-3, 1e-2, 0.1])
    axes[0, 0].set_xticklabels(["-0.1", "1e-4", "1e-3", "1e-2", "0.1"])
    fig.tight_layout()
    save(fig, "fig_state_dose")


def r12():
    from scripts.submit_r12 import RPMS
    ids = json.loads((MF / "test_t5_roughness.json").read_text())["run_ids"]
    names = ("kLa_1T_10", "dtmix_0.50", "kLa_1T_25", "dtmix_0.95")
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for k, q in enumerate(names):
        y = np.array([100 * (val(f"r12_l8_rpm{r:g}_th7", q) / val(ids[f"{r:g}"], q) - 1) for r in RPMS])
        jit = np.linspace(-0.18, 0.18, len(y))
        ax.plot(k + jit, y, ls="none", **mkw(8, q))
        ax.plot([k - 0.28, k + 0.28], [y.mean()] * 2, color=fs.level_colour(8), lw=1.6)
    ax.axhline(0, color="0.2", lw=1.0)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([QL[q][1] for q in names], fontsize=9)
    ax.set_ylabel("change of result, release after 12 vs 80 cycles (%)")
    ax.grid(**fs.GRID_KW, axis="y")
    fig.tight_layout()
    save(fig, "fig_state_r12")


def e9():
    pairs = {"22.4 vs 22.6 rpm": ("22.4", "22.6"), "29.9 vs 30.1 rpm": ("29.9", "30.1")}
    qs = ("kLa_1T_10", "dtmix_0.50", "dtmix_0.95")
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), sharey=True)
    for ax, (title, (a, b)) in zip(axes, pairs.items()):
        for k, (L, pre) in enumerate(((8, "t5_l8_rpm{}_th7"), (9, "e9_l9_rpm{}_th7"))):
            for i, q in enumerate(qs):
                va, vb = val(pre.format(a), q), val(pre.format(b), q)
                x = i + (k - 0.5) * 0.25
                if np.isfinite(va) and np.isfinite(vb):
                    ax.plot(x, 100 * abs(va - vb) / ((va + vb) / 2), ls="none", **mkw(L, q))
                elif np.isfinite(va):    # 22.6 rpm at L9 never reached chi 0.95: the true gap is larger than the one shown
                    vb = 374.0           # time reached by the run (diary 2026-10-18)
                    ax.plot(x, 100 * abs(va - vb) / ((va + vb) / 2), ls="none", **mkw(L, q))
                    ax.annotate("", xy=(x, 100 * abs(va - vb) / ((va + vb) / 2) + 28), xytext=(x, 100 * abs(va - vb) / ((va + vb) / 2) + 4),
                                arrowprops=dict(arrowstyle="->", color=fs.level_colour(L)))
        ax.set_xticks(range(len(qs)))
        ax.set_xticklabels([QL[q][1] for q in qs], fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.grid(**fs.GRID_KW, axis="y")
    axes[0].set_ylabel("difference between the two runs (%)")
    axes[0].set_ylim(-2, 100)
    h = [Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(L), ms=6, label=f"L{L}") for L in (8, 9)]
    lg = axes[1].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    save(fig, "fig_state_e9", [lg])


if __name__ == "__main__":
    plt.rcParams.update(fs.rcparams())
    for fn in (levels, scatter, dose, r12, e9):
        fn()
