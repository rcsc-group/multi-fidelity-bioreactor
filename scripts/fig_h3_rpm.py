"""Figure: L10 QoIs vs rocking speed, simulation and multi-fidelity prediction, per rocking angle.

Rows: tau_95, mean tau, EDR. Columns: theta = 4, 7, 9 deg. x: rocking speed (rpm).
  filled/hollow markers: simulations at L8 (13 x... all L8 runs), L9 and L10; L10 bars = spread of the last three cycles.
  dashed lines: L8 fit (KRR) and L9 fit (GP).
  grey line + band: three-level prediction trained on all 13 L10 runs, with 95% interval incl. the L9 uncertainty.
  grey bars: the same model trained on the 12 other L10 runs, predicting that run (each L10 run is left out once).
Method: scripts/test_h3_mfbml.py (log QoI, FullBayesKRRLRGP), band propagation: scripts/h3_band_check.py (variant B).
Usage: uv run python scripts/fig_h3_rpm.py
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import figstyle as fs  # noqa: E402
from scripts.h3_band_check import DELTA, Offset, sd_log  # noqa: E402
from scripts.hydro_dataset_v2 import COL_EDISS, RHO, qois  # noqa: E402
from scripts.postprocess import _load_dat, _t_scales  # noqa: E402
from scripts.test_h1_mfbml import (BOUNDS, LF_THETAS, RPMS, hf_points, hf_run, lf_run,  # noqa: E402
                                   load_points, train_lf2d)
from scripts.test_h2_mfbml import l9_points, l9_run  # noqa: E402

MF = ROOT / "experiments/multifidelity"
THETAS = (4.0, 7.0, 9.0)
ROWS = (("tau_95_t", fs.MK["tau"], "max", r"$\tau_{95}$ (Pa)"),
        ("tau_mean_t", fs.MK["tau"], "mean", r"$\overline{\tau}$ (Pa)"),
        ("edr_mean_t", fs.MK["ediss"], "mean", r"EDR (W/m$^3$)"))
GRID = np.linspace(17.5, 32.5, 31)


def cycle_sd(run, q):
    """Relative sd of the QoI over the three whole cycles of the QoI window (same window as hydro_dataset_v2.qois)."""
    h = qois(run)
    d = ROOT / "runs" / run
    import math
    p = json.load(open(d / "params.json"))
    arr = _load_dat(d / "shear_stress.dat")
    T_bio, T = _t_scales(p)
    U = p.get("geometry", {}).get("a", 0.25) / T_bio
    L = p.get("geometry", {}).get("a", 0.25)
    cyc = arr[:, 1] / T
    c0, c1 = h["window"]
    vals = []
    for k in range(int(round(c1 - c0))):
        m = (cyc >= c0 + k) & (cyc < c0 + k + 1)
        a = arr[m]
        vals.append({"tau_mean_t": a[:, 5].mean(), "tau_95_t": np.median(a[:, 2]), "edr_mean_t": a[:, COL_EDISS].mean(),
                     "tau_amp_t": 0.5 * (a[:, 10].max() - a[:, 10].min()), "edr_max_t": a[:, COL_EDISS].max()}[q])
    vals = np.array(vals)
    return float(vals.std(ddof=1) / vals.mean())


def full_model(q, lf_pts, l9_pts, hf_pts):
    from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP
    k8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
    lf = train_lf2d(np.array(k8), np.log([lf_pts[k][q] for k in k8]))
    k9 = sorted(l9_pts, key=lambda k: (k[1], k[0]))
    m1 = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(np.array(k9), np.log([l9_pts[k][q] for k in k9]))
    kh = sorted(hf_pts, key=lambda k: (k[1], k[0]))
    ad = Offset(m1)
    m2 = FullBayesKRRLRGP(BOUNDS, lf_model=ad, basis="linear").train(np.array(kh), np.log([hf_pts[k][q] for k in kh]))

    def pred(theta):
        X = np.column_stack([GRID, np.full_like(GRID, theta)])
        mean = np.ravel(m2.predict_mean(X))
        qq = m2.quantile(X, [0.025, 0.975])
        ad.off = DELTA
        rho = (np.ravel(m2.predict_mean(X)) - mean) / DELTA
        ad.off = 0.0
        hw = np.sqrt(((qq[:, 1] - qq[:, 0]) / 2) ** 2 + (1.96 * rho * sd_log(m1, X)) ** 2)
        return {8: np.exp(np.ravel(lf.predict(X))),
                9: np.exp(np.ravel(m1.predict_mean(X))),
                10: (np.exp(mean), np.exp(mean - hw), np.exp(mean + hw))}
    return pred


def main():
    warnings.filterwarnings("ignore")
    lf_pts, _ = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
    l9_pts, _ = load_points(l9_points(), l9_run)
    hf_pts, _ = load_points(hf_points(), hf_run)
    loo = json.loads((MF / "h3_band_check.json").read_text())
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(3, 3, figsize=(9.2, 7.6), sharex=True, sharey=False)
    c10, cp = fs.level_colour(10), "0.25"
    for i, (q, mk, stat, lab) in enumerate(ROWS):
        pred = full_model(q, lf_pts, l9_pts, hf_pts)
        for j, th in enumerate(THETAS):
            ax = axes[i, j]
            curves = pred(th)
            for L, src in ((8, lf_pts), (9, l9_pts)):
                cL = fs.level_colour(L)
                ax.plot(GRID, curves[L], color=cL, lw=1.0, ls="--")
                ks = sorted(k for k in src if k[1] == th and GRID[0] <= k[0] <= GRID[-1])
                ax.plot([k[0] for k in ks], [src[k][q] for k in ks], ls="none",
                        **{k: v for k, v in fs.series_kw(cL, mk, stat=stat, ours=True).items() if k != "ls"})
            m, lo, hi = curves[10]
            ax.fill_between(GRID, lo, hi, color=cp, alpha=0.15, lw=0)
            ax.plot(GRID, m, color=cp, lw=1.0)
            pts = [p for p in loo[q]["points_B"] if p["theta"] == th]
            r = np.array([p["rpm"] for p in pts])
            ax.errorbar(r + 0.35, [np.exp(p["mean"]) for p in pts],
                        yerr=[[np.exp(p["mean"]) - np.exp(p["lo"]) for p in pts], [np.exp(p["hi"]) - np.exp(p["mean"]) for p in pts]],
                        ls="none", color=cp, marker="_", ms=6, elinewidth=1.0, capsize=2)
            obs = np.array([np.exp(p["obs"]) for p in pts])
            sd = np.array([cycle_sd(hf_run(p["rpm"], th), q) for p in pts])
            ax.errorbar(r, obs, yerr=obs * sd, ls="none", color=c10, elinewidth=1.0, capsize=0, zorder=3)
            ax.plot(r, obs, ls="none", **{k: v for k, v in fs.series_kw(c10, mk, stat=stat, ours=True).items() if k != "ls"})
            ax.grid(**fs.GRID_KW)
            if i == 0:
                ax.set_title(rf"$\theta$ = {th:g}$^\circ$", fontsize=10)
            if j == 0:
                ax.set_ylabel(lab)
            if i == 2:
                ax.set_xlabel("rocking speed (rpm)")
    handles = [Line2D([], [], color=fs.level_colour(L), marker="o", mfc=fs.level_colour(L), ls="none", ms=5, label=f"L{L}")
               for L in (8, 9, 10)] + \
              [Patch(facecolor=cp, alpha=0.3, label="prediction"),
               Line2D([], [], color=cp, marker="_", ls="none", ms=8, mew=1.0, label="prediction without that run")]
    lg = axes[0, 2].legend(handles=handles, fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.04, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(MF / f"fig_h3_rpm.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])
        print("wrote", MF / f"fig_h3_rpm.{ext}")


if __name__ == "__main__":
    main()
