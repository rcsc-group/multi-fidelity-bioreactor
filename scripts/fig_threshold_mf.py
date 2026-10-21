"""PI idea 1 as a multi-fidelity problem (pre-registered, diary 2026-10-21): the fidelity is the C* threshold at which a run may be stopped.
LF = kLa at the lower threshold (known at all L8 runs), HF = kLa at the higher threshold (known at n runs). Learning curves of three predictors on the other runs:
MF (FullBayesKRRLRGP), GPR with a good prior (quadratic mean + ARD GP on the residual), ratio (ln HF = ln LF + one offset).
Panel 0 shows the problem: C*(t) of one run with the four thresholds. Usage: uv run python scripts/fig_threshold_mf.py
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from scripts import figstyle as fs  # noqa: E402
from scripts import postprocess as pp  # noqa: E402
from scripts.test_h1_mfbml import BOUNDS, train_lf2d  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

PAIRS = ((10, 50), (25, 50), (25, 75), (50, 75))
COL = {10: 2, 25: 3, 50: 4, 75: 5}
SIZES = (4, 6, 8, 12, 16, 24, 32)
DRAWS = 15


def load():
    rows = [r for r in csv.reader((ROOT / "experiments/multifidelity/kla_thresholds.csv").open())][1:]
    A = np.array([[float(x) for x in r[2:]] for r in rows if r[0] == "8"], float)
    return A[np.all(np.isfinite(A) & (A > 0), axis=1)]


def quad(X):
    x, t = (X[:, 0] - 25) / 7.5, (X[:, 1] - 6.5) / 2.5
    return np.column_stack([np.ones_like(x), x, t, x * x, x * t, t * t])


def lin(X):
    return quad(X)[:, :3]


def gpr_good_prior(Xtr, ytr, Xte, basis=lin):
    """GP on the residual of a low-order mean (linear in rpm and theta: 3 coefficients, fair at n >= 4)."""
    c, *_ = np.linalg.lstsq(basis(Xtr), ytr, rcond=None)
    res = ytr - basis(Xtr) @ c
    Z = lambda X: np.column_stack([(X[:, 0] - 25) / 7.5, (X[:, 1] - 6.5) / 2.5])   # noqa: E731
    k = ConstantKernel(0.02, (1e-4, 1.0)) * RBF(np.ones(2), (0.2, 20.0)) + WhiteKernel(0.02, (1e-4, 0.5))
    gp = GaussianProcessRegressor(k, normalize_y=False, n_restarts_optimizer=0, random_state=0).fit(Z(Xtr), res)
    return basis(Xte) @ c + gp.predict(Z(Xte))


def gpr_lf_prior(Xtr, ytr, ltr, Xte, lte):
    """GP whose prior mean is the lower-threshold value plus one offset; the GP models the residual."""
    off = np.mean(ytr - ltr)
    res = ytr - ltr - off
    Z = lambda X: np.column_stack([(X[:, 0] - 25) / 7.5, (X[:, 1] - 6.5) / 2.5])   # noqa: E731
    k = ConstantKernel(0.02, (1e-4, 1.0)) * RBF(np.ones(2), (0.2, 20.0)) + WhiteKernel(0.02, (1e-4, 0.5))
    gp = GaussianProcessRegressor(k, normalize_y=False, n_restarts_optimizer=0, random_state=0).fit(Z(Xtr), res)
    return lte + off + gp.predict(Z(Xte))


def main():
    warnings.filterwarnings("ignore")
    A = load()
    X = A[:, :2]
    N = len(A)
    print(f"{N} L8 runs with all four thresholds")
    rng = np.random.default_rng(7)
    out = {}
    for lo, hi in PAIRS:
        ylf, yhf = np.log(A[:, COL[lo]]), np.log(A[:, COL[hi]])
        lf = train_lf2d(X, ylf)
        flf = np.ravel(lf.predict(X))
        curve = {}
        for n in SIZES:
            acc = {"MF": [], "GPR": [], "GPR_LF": [], "ratio": []}
            for _ in range(DRAWS):
                tr = np.zeros(N, bool)
                tr[rng.choice(N, n, replace=False)] = True
                # drop runs within 0.7 rpm at the same theta from the test set of each training run (near-replicates)
                near = np.zeros(N, bool)
                for i in np.where(tr)[0]:
                    near |= (np.abs(X[:, 0] - X[i, 0]) <= 0.7) & (X[:, 1] == X[i, 1])
                te = ~near
                if te.sum() < 5:
                    continue
                m = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(X[tr], yhf[tr])
                acc["MF"].append(np.sqrt(np.mean((np.ravel(m.predict_mean(X[te])) - yhf[te]) ** 2)))
                acc["GPR"].append(np.sqrt(np.mean((gpr_good_prior(X[tr], yhf[tr], X[te]) - yhf[te]) ** 2)))
                acc["GPR_LF"].append(np.sqrt(np.mean((gpr_lf_prior(X[tr], yhf[tr], ylf[tr], X[te], ylf[te]) - yhf[te]) ** 2)))
                acc["ratio"].append(np.sqrt(np.mean((ylf[te] + np.mean(yhf[tr] - ylf[tr]) - yhf[te]) ** 2)))
            curve[n] = {k: float(np.mean(v)) for k, v in acc.items()}
        out[f"{lo}->{hi}"] = curve
        print(f"{lo / 100:.2f} -> {hi / 100:.2f}: " + " | ".join(f"n={n}: MF {c['MF']:.2f} GPR(lin) {c['GPR']:.2f} GPR(LF prior) {c['GPR_LF']:.2f} ratio {c['ratio']:.2f}" for n, c in curve.items() if n in (4, 8, 16, 32)))
    adv = sum(1 for k, c in out.items() if all(c[n]["MF"] <= 0.8 * min(c[n]["GPR"], c[n]["GPR_LF"], c[n]["ratio"]) for n in (4, 6, 8)))
    verdict = "REAL ADVANTAGE" if adv >= 3 else "NO ADVANTAGE OVER A GPR WITH A GOOD PRIOR"
    print("verdict:", verdict, f"({adv} of 4 pairs meet the 20 % rule at n <= 8)")
    out["verdict"] = verdict
    (ROOT / "experiments/multifidelity/threshold_mf.json").write_text(json.dumps(out, indent=1))

    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 5, figsize=(17.0, 3.5))
    ax = axes[0]
    run = ROOT / "runs" / "kmix_l8_rpm25"
    p = json.loads((run / "params.json").read_text())
    t, c = pp._compute_c_star(run)
    T_bio, T_nd = pp._t_scales(p)
    tt = (t - t[np.argmax(c > 1e-6)]) * T_bio
    ax.plot(tt, c, color=fs.level_colour(8), lw=1.0)
    for x in (0.10, 0.25, 0.50, 0.75):
        ax.axhline(x, color=fs.GUIDE, lw=1.0, ls=":")
    ax.set_xlim(0, 700)
    ax.set_xlabel("time after release (s)")
    ax.set_ylabel(r"$C^*$")
    ax.set_title("the problem: L8, 25 rpm, 7$^\\circ$", fontsize=10)
    ax.grid(**fs.GRID_KW)
    sty = {"MF": (fs.level_colour(8), "o", "-"), "GPR": ("0.25", "s", "-"), "GPR_LF": ("0.55", "D", "--"), "ratio": (fs.level_colour(9), "^", ":")}
    for ax, (lo, hi) in zip(axes[1:], PAIRS):
        cv = out[f"{lo}->{hi}"]
        for k, (col, mk, ls) in sty.items():
            ax.plot(SIZES, [cv[n][k] for n in SIZES], color=col, marker=mk, mfc="w", mec=col, ls=ls, lw=1.0, ms=5)
        ax.set_xscale("log")
        ax.set_ylim(0, None)
        ax.set_xticks([4, 8, 16, 32])
        ax.set_xticklabels(["4", "8", "16", "32"])
        ax.set_xlabel("runs with the higher threshold")
        ax.set_title(rf"$C^*$ {lo / 100:.2f} $\to$ {hi / 100:.2f}", fontsize=10)
        ax.grid(**fs.GRID_KW)
    axes[1].set_ylabel("error on unseen runs (RMS of ln ratio)")
    h = [Line2D([], [], color=sty["MF"][0], marker="o", mfc="w", label="mfbml"), Line2D([], [], color=sty["GPR"][0], marker="s", mfc="w", label="GPR, linear mean prior"),
         Line2D([], [], color=sty["GPR_LF"][0], ls="--", marker="D", mfc="w", label="GPR, lower threshold as prior"),
         Line2D([], [], color=sty["ratio"][0], ls=":", marker="^", mfc="w", label="lower threshold times a constant")]
    lg = axes[-1].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / f"experiments/multifidelity/fig_threshold_mf.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])


if __name__ == "__main__":
    main()
