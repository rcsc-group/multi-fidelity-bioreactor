"""Figure: L10 QoIs vs rocking speed, simulation and H3 prediction with 95% uncertainty, per rocking angle.

Rows: tau_95, mean tau, EDR. Columns: theta = 4, 7, 9 deg. x: rocking speed (rpm).
  markers   = L10 simulation (13 runs).
  bars      = H3 prediction at that point with the point LEFT OUT of the training (leave-one-out, as test H3).
  line+band = H3 prediction trained on all 13 L10 points, posterior mean and 95% interval.
H3 = three-level chain on log(QoI): KRR on 40 L8 runs -> GP correction to 16 L9 runs -> GP correction to the
L10 runs (FullBayesKRRLRGP, mfbml_local). Data/models: scripts/test_h1_mfbml.py, test_h2_mfbml.py, test_h3_mfbml.py.
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
from scripts.test_h1_mfbml import (BOUNDS, LF_THETAS, RPMS, hf_points, hf_run, lf_run,  # noqa: E402
                                   load_points, train_lf2d)
from scripts.test_h2_mfbml import MeanAdapter, l9_points, l9_run  # noqa: E402

MF = ROOT / "experiments/multifidelity"
THETAS = (4.0, 7.0, 9.0)
ROWS = (("tau_95_t", fs.MK["tau"], "max", r"$\tau_{95}$ (Pa)"),
        ("tau_mean_t", fs.MK["tau"], "mean", r"$\overline{\tau}$ (Pa)"),
        ("edr_mean_t", fs.MK["ediss"], "mean", r"EDR (W/m$^3$)"))
GRID = np.linspace(17.5, 32.5, 31)


def full_model(q, lf_pts, l9_pts, hf_pts):
    """H3 trained on all points; returns a function rpm-grid, theta -> (mean, lo, hi) on the raw scale."""
    from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP
    k8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
    lf = train_lf2d(np.array(k8), np.log([lf_pts[k][q] for k in k8]))
    k9 = sorted(l9_pts, key=lambda k: (k[1], k[0]))
    m1 = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(
        np.array(k9), np.log([l9_pts[k][q] for k in k9]))
    kh = sorted(hf_pts, key=lambda k: (k[1], k[0]))
    m2 = FullBayesKRRLRGP(BOUNDS, lf_model=MeanAdapter(m1), basis="linear").train(
        np.array(kh), np.log([hf_pts[k][q] for k in kh]))

    def pred(theta):
        X = np.column_stack([GRID, np.full_like(GRID, theta)])
        qq = m2.quantile(X, [0.025, 0.975])
        return np.exp(np.ravel(m2.predict_mean(X))), np.exp(qq[:, 0]), np.exp(qq[:, 1])
    return pred


def main():
    warnings.filterwarnings("ignore")
    lf_pts, _ = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
    l9_pts, _ = load_points(l9_points(), l9_run)
    hf_pts, _ = load_points(hf_points(), hf_run)
    loo = json.loads((MF / "test_h3_mfbml.json").read_text())["results"]
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(3, 3, figsize=(10.0, 7.6), sharex=True, sharey="row")
    c10, cp = fs.level_colour(10), "0.25"
    for i, (q, mk, stat, lab) in enumerate(ROWS):
        pred = full_model(q, lf_pts, l9_pts, hf_pts)
        for j, th in enumerate(THETAS):
            ax = axes[i, j]
            m, lo, hi = pred(th)
            ax.fill_between(GRID, lo, hi, color=cp, alpha=0.15, lw=0)
            ax.plot(GRID, m, color=cp, lw=1.0)
            pts = [p for p in loo[q]["points"] if p["theta"] == th]
            r = np.array([p["rpm"] for p in pts])
            ax.errorbar(r + 0.35, [p["mf_pred"] for p in pts],
                        yerr=[[p["mf_pred"] - p["mf_lo"] for p in pts], [p["mf_hi"] - p["mf_pred"] for p in pts]],
                        ls="none", color=cp, marker="_", ms=6, elinewidth=1.0, capsize=2)
            ax.plot(r, [p["obs"] for p in pts], ls="none",
                    **{k: v for k, v in fs.series_kw(c10, mk, stat=stat, ours=True).items() if k != "ls"})
            ax.grid(**fs.GRID_KW)
            if i == 0:
                ax.set_title(rf"$\theta$ = {th:g}$^\circ$", fontsize=10)
            if j == 0:
                ax.set_ylabel(lab)
            if i == 2:
                ax.set_xlabel("rocking speed (rpm)")
    handles = [Line2D([], [], color=c10, marker="o", mfc=c10, ls="none", ms=5, label="L10 simulation"),
               Line2D([], [], color=cp, marker="_", ls="none", ms=8, mew=1.0,
                      label="prediction, point left out\n(95% interval)"),
               Patch(facecolor=cp, alpha=0.15, label="prediction, all points\n(95% interval)"),
               Line2D([], [], color=cp, lw=1.0, label="prediction, all points\n(mean)")]
    lg = axes[0, 2].legend(handles=handles, fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.04, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = MF / f"fig_h3_rpm.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])
        print("wrote", out)


if __name__ == "__main__":
    main()
