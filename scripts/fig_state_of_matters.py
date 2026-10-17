"""Overview figures of the scalar-QoI evidence (kLa, mixing time) for the DFD plan, from runs/*/results.json and the test jsons.

fig_state_levels   kLa_1T_10 and dtmix_0.50 vs rpm, theta 7, per level (L9 split by protocol).
fig_state_scatter  L8 0.625-rpm grids at theta 4, 7, 9, four QoI, quadratic trend (tests E0, E8).
fig_state_dose     test D: |d log q| vs the rpm perturbation, base 22.5 and 30 rpm.
fig_state_r12      test R12: log ratio of release at cycle 12 vs 80 per rpm.
fig_state_e9       pairs 0.2 rpm apart: |d log q| at L8 and L9 (same binary and protocol).
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
from scripts.test_e8_theta_scatter import GRID, run_id  # noqa: E402

MF = ROOT / "experiments/multifidelity"
QL = {"kLa_1T_10": (fs.MK["cstar_10"], r"$k_La$, $C^*\!=\!0.10$ (1/h)"),
      "kLa_1T_25": (fs.MK["cstar_25"], r"$k_La$, $C^*\!=\!0.25$ (1/h)"),
      "dtmix_0.50": (fs.MK["chi_0.50"], r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$ (s)"),
      "dtmix_0.95": (fs.MK["chi_0.95"], r"$t_{\mathrm{mix}}$, $\chi\!=\!0.95$ (s)")}


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


def levels():
    rows = [r for r in csv.DictReader((MF / "kla_mix_dataset.csv").open()) if float(r["theta"]) == 7.0]
    l8 = [(float(r["rpm"]), r) for r in rows if r["level"] == "8"]
    l9_old = [(float(r["rpm"]), r) for r in rows if r["level"] == "9"]
    e9 = [r for r in (ROOT / "runs").glob("e9_l9_rpm*_th7")]
    l10 = [("l10c_rpm32.5_seg1", 1.0), ("v_l10_rpm32.5_th7", 1.0), ("fig11_l10_rpm35", 0.35), ("fig11_l10_rpm37.5", 0.35)]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.9), sharex=True)
    for ax, q in zip(axes, ("kLa_1T_10", "dtmix_0.50")):
        mk = QL[q][0]
        def pts(L, items, alpha=1.0, label=None):
            xy = np.array([(x, v) for x, v in items if np.isfinite(v)])
            if len(xy):
                kw = fs.series_kw(fs.level_colour(L), mk, stat="mean", ours=True)
                kw.pop("ls")
                ax.plot(xy[:, 0], xy[:, 1], ls="none", alpha=alpha, **kw)
        pts(8, [(x, val(r["run"], q)) for x, r in l8])
        pts(9, [(x, val(r["run"], q)) for x, r in l9_old], alpha=0.4)
        pts(9, [(rpm_of(r.name), val(r.name, q)) for r in e9])
        pts(10, [(rpm_of(n), val(n, q)) for n, a in l10 if a == 1.0])
        pts(10, [(rpm_of(n), val(n, q)) for n, a in l10 if a < 1.0], alpha=0.35)
        ax.set_yscale("log")
        ax.set_ylabel(QL[q][1])
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
    fig, axes = plt.subplots(4, 3, figsize=(10.0, 10.0), sharex=True)
    for j, th in enumerate((4.0, 7.0, 9.0)):
        runs = [ids[f"{r:g}"] for r in GRID] if th == 7.0 else [run_id(r, th) for r in GRID]
        for i, q in enumerate(QL):
            ax = axes[i, j]
            y = np.array([val(rn, q) for rn in runs])
            ok = np.isfinite(y)
            x = np.array(GRID)
            kw = fs.series_kw(fs.level_colour(8), QL[q][0], stat="mean", ours=True)
            kw.pop("ls")
            ax.plot(x[ok], y[ok], ls="none", **kw)
            c = np.polyfit(x[ok], np.log(y[ok]), 2)
            xx = np.linspace(20, 32.5, 100)
            ax.plot(xx, np.exp(np.polyval(c, xx)), color=fs.level_colour(8), ls=":", lw=1.0)
            ax.set_yscale("log")
            ax.grid(**fs.GRID_KW)
            if j == 0:
                ax.set_ylabel(QL[q][1])
            if i == 0:
                ax.set_title(rf"$\theta$ = {th:g}$^\circ$", fontsize=10)
            if i == 3:
                ax.set_xlabel("rocking speed (rpm)")
    fig.tight_layout()
    save(fig, "fig_state_scatter")


def dose():
    ids = json.loads((MF / "test_t5_roughness.json").read_text())["run_ids"]
    eps = (1e-4, 1e-3, 3e-3, 1e-2, 3e-2)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
    pairs = {22.5: ("t5_l8_rpm22.4_th7", "t5_l8_rpm22.6_th7"), 30.0: ("t5_l8_rpm29.9_th7", "t5_l8_rpm30.1_th7")}
    ls = {22.5: "-", 30.0: "--"}
    for ax, q in zip(axes, ("dtmix_0.50", "dtmix_0.95", "kLa_1T_10")):
        for base in (22.5, 30.0):
            q0 = val(ids[f"{base:g}"], q)
            d = [abs(np.log(val(f"dose_l8_rpm{base + e:g}_th7", q)) - np.log(q0)) for e in eps]
            d01 = [abs(np.log(val(p, q)) - np.log(q0)) for p in pairs[base]]
            kw = fs.series_kw(fs.level_colour(8), QL[q][0], stat="mean", ours=True)
            kw["ls"] = ls[base]
            ax.plot(eps, np.maximum(d, 3e-4), **kw)
            ax.plot([0.1, 0.1], d01, ls="none", marker=QL[q][0], mfc="w", mec=fs.level_colour(8), ms=5, mew=1.1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(3e-4, 1.0)
        ax.set_xlabel("rpm perturbation")
        ax.set_title(QL[q][1], fontsize=10)
        ax.grid(**fs.GRID_KW)
    axes[0].set_ylabel(r"$|\Delta \ln q|$")
    h = [Line2D([], [], color=fs.level_colour(8), ls="-", label="base 22.5 rpm"),
         Line2D([], [], color=fs.level_colour(8), ls="--", label="base 30 rpm")]
    lg = axes[2].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    save(fig, "fig_state_dose", [lg])


def r12():
    d = json.loads((MF / "test_r12_release.json").read_text())
    from scripts.submit_r12 import RPMS
    names = ("kLa_1T_10", "dtmix_0.50", "kLa_1T_25", "dtmix_0.95", "tau_95_t", "edr_mean_t")
    labs = {**{k: v[1] for k, v in QL.items()}, "tau_95_t": r"$\tau_{95}$", "edr_mean_t": "EDR"}
    mks = {**{k: v[0] for k, v in QL.items()}, "tau_95_t": fs.MK["tau"], "edr_mean_t": fs.MK["ediss"]}
    fig, axes = plt.subplots(2, 3, figsize=(10.0, 5.6), sharex=True)
    for ax, q in zip(axes.ravel(), names):
        y = np.array(d[q]["d"])
        kw = fs.series_kw(fs.level_colour(8), mks[q], stat="mean", ours=True)
        kw.pop("ls")
        ax.plot(RPMS, y, ls="none", **kw)
        ax.axhline(0, color=fs.GUIDE, lw=1.0)
        ax.axhline(d[q]["mean"], color=fs.level_colour(8), ls=":", lw=1.0)
        ax.set_title(labs[q], fontsize=10)
        ax.grid(**fs.GRID_KW)
    for ax in axes[1]:
        ax.set_xlabel("rocking speed (rpm)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\ln(q_{12}/q_{80})$")
    fig.tight_layout()
    save(fig, "fig_state_r12")


def e9():
    pairs = {"22.4 / 22.6 rpm": ("22.4", "22.6"), "29.9 / 30.1 rpm": ("29.9", "30.1")}
    qs = ("kLa_1T_10", "dtmix_0.50", "dtmix_0.95")
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.7), sharey=True)
    for ax, (title, (a, b)) in zip(axes, pairs.items()):
        for k, (L, pre) in enumerate(((8, "t5_l8_rpm{}_th7"), (9, "e9_l9_rpm{}_th7"))):
            for i, q in enumerate(qs):
                va, vb = val(pre.format(a), q), val(pre.format(b), q)
                d = abs(np.log(va) - np.log(vb))
                if np.isfinite(d):
                    ax.plot(i + (k - 0.5) * 0.25, d, ls="none", marker=QL[q][0], mfc="w", mec=fs.level_colour(L), ms=8, mew=1.3)
        ax.set_xticks(range(len(qs)))
        ax.set_xticklabels([r"$k_La$ 0.10", r"$t_{\mathrm{mix}}$ 0.50", r"$t_{\mathrm{mix}}$ 0.95"], fontsize=9)
        ax.set_title(title, fontsize=10)
        ax.grid(**fs.GRID_KW)
    axes[0].set_ylabel(r"$|\Delta \ln q|$ between the two runs")
    h = [Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(L), ms=6, label=f"L{L}") for L in (8, 9)]
    lg = axes[1].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    save(fig, "fig_state_e9", [lg])


if __name__ == "__main__":
    plt.rcParams.update(fs.rcparams())
    for fn in (levels, scatter, dose, r12, e9):
        fn()
