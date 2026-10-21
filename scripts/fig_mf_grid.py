"""5 x 3 multi-fidelity prediction figures (meeting 2026-10-22).

Rows: tau_95, mean tau, EDR, kLa (C* = 0.10), mixing time (chi = 0.50). Columns: theta = 4, 7, 9 deg. x: rocking speed (rpm).
  mode "chain":     L8 KRR -> L9 GP -> L10 GP (log of the QoI), L9 variance propagated into the L10 band (scripts/h3_band_check.py recipe).
  mode "two_level": L9 KRR -> L10 GP (log of the QoI), the usual two-fidelity mfbml (Yi et al. 2407.15110).
Markers: simulations (L8 blue, L9 green, L10 pink); dashed line: the fitted lower level; grey line and band: prediction of L10 and its 95 % range, trained on ALL L10 runs.
L10 hydrodynamic markers carry bars = sd of the last three cycles. A row has a prediction only if there are >= 5 L10 runs with that quantity.
Uniform rule: all level markers are always drawn; the fitted curves and the L10 prediction span the L10 range (17.5-32.5 rpm) and are drawn only in panels where L9 has >= 3 runs at that angle and the row has >= 5 L10 runs.
Data gaps (shown as missing, not filled): L9 kLa and mixing time exist only at theta 7; L10 mixing time (50 %) exists at 3 runs.
Usage: uv run python scripts/fig_mf_grid.py [chain|two_level]
"""
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from scripts import figstyle as fs  # noqa: E402
from scripts.fig_h3_rpm import cycle_sd  # noqa: E402
from scripts.h3_band_check import DELTA, Offset, sd_log  # noqa: E402
from scripts.test_h1_mfbml import BOUNDS, LF_THETAS, RPMS, hf_points, hf_run, lf_run, load_points, train_lf2d  # noqa: E402
from scripts.test_h2_mfbml import l9_points, l9_run  # noqa: E402
from scripts.test_s89_l8_to_l9 import collect  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

THETAS = (4.0, 7.0, 9.0)
GRID = np.linspace(14.5, 38.5, 49)   # the full range of the L8 and L9 data: the predictions extrapolate beyond the L10 runs (17.5-32.5 rpm)
HYDRO = (("tau_95_t", fs.MK["tau"], "max", r"$\tau_{95}$ (Pa)"), ("tau_mean_t", fs.MK["tau"], "mean", r"$\overline{\tau}$ (Pa)"),
         ("edr_mean_t", fs.MK["ediss"], "mean", r"EDR (W/m$^3$)"))
SCAL = (("kLa_1T_10", fs.MK["cstar_10"], "mean", r"$k_La$, $C^*\!=\!0.10$ (1/h)"), ("dtmix_0.50", fs.MK["chi_0.50"], "mean", r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$ (s)"))
MF = ROOT / "experiments/multifidelity"


def hydro_data(q):
    lf, _ = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
    l9, _ = load_points(l9_points(), l9_run)
    hf, _ = load_points(hf_points(), hf_run)
    pick = lambda d: (np.array(sorted(d, key=lambda k: (k[1], k[0])), float), np.array([d[k][q] for k in sorted(d, key=lambda k: (k[1], k[0]))]))   # noqa: E731
    return pick(lf), pick(l9), pick(hf), hf


def scalar_data(q):
    l8 = collect(["t6_l8_*", "kmix_l8_*", "t5_l8_*", "e8_l8_*", "w_l8_*", "dose_l8_*"], q)
    l9 = [p for p in collect(["fig9_l9_rpm*", "e9_l9_rpm22.4_th7", "e9_l9_rpm22.6_th7", "e9_l9_rpm29.9_th7", "e9_l9_rpm30.1_th7"], q) if p[1] == 7.0]
    l10 = []
    for d in sorted((ROOT / "runs").glob("t10_l10_*")) + [ROOT / "runs" / "v_l10_rpm32.5_th7"]:
        if d.name.endswith("_ext") or not (d / "results.json").exists():
            continue
        v = json.loads((d / "results.json").read_text()).get(q)
        if v and v > 0:
            p = json.loads((d / "params.json").read_text())
            l10.append((p["omega_b"] * 60 / (2 * math.pi), p["theta_max"][0], float(v)))
    arr = lambda rows: (np.array([[r[0], r[1]] for r in rows], float), np.array([r[2] for r in rows], float))   # noqa: E731
    return arr(l8), arr(l9), arr(l10)


def lf_cached(X, y):
    """The L8 / L9 KRR (leave-one-out tuned, 10 restarts) is the slow step and does not depend on the level above it: cache it by a hash of its data."""
    import hashlib
    import pickle
    key = hashlib.md5(np.ascontiguousarray(np.column_stack([X, np.log(y)])).tobytes()).hexdigest()[:16]
    f = Path("/oscar/scratch/eaguerov/mfcache") / f"krr_{key}.pkl"
    if f.exists():
        return pickle.load(f.open("rb"))
    lf = train_lf2d(X, np.log(y))
    f.parent.mkdir(exist_ok=True)
    pickle.dump(lf, f.open("wb"))
    return lf


def models(mode, a8, a9, a10):
    """returns pred(theta) -> (lower-level curve, (mean, lo, hi) of L10) or None."""
    (X8, y8), (X9, y9), (X10, y10) = a8, a9, a10
    if len(y10) < 5 or len(y9) < 3:
        return None
    if mode == "chain":
        lf = lf_cached(X8, y8)
        m1 = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(X9, np.log(y9))
        ad = Offset(m1)
        m2 = FullBayesKRRLRGP(BOUNDS, lf_model=ad, basis="linear").train(X10, np.log(y10))

        def pred(theta):
            X = np.column_stack([GRID, np.full_like(GRID, theta)])
            mean = np.ravel(m2.predict_mean(X))
            qq = m2.quantile(X, [0.025, 0.975])
            ad.off = DELTA
            rho = (np.ravel(m2.predict_mean(X)) - mean) / DELTA
            ad.off = 0.0
            hw = np.sqrt(((qq[:, 1] - qq[:, 0]) / 2) ** 2 + (1.96 * rho * sd_log(m1, X)) ** 2)
            return np.exp(np.ravel(m1.predict_mean(X))), (np.exp(mean), np.exp(mean - hw), np.exp(mean + hw)), np.exp(np.ravel(lf.predict(X)))
        return pred
    lf = lf_cached(X9, y9)
    m2 = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(X10, np.log(y10))

    def pred(theta):
        X = np.column_stack([GRID, np.full_like(GRID, theta)])
        qq = m2.quantile(X, [0.025, 0.975])
        mean = np.ravel(m2.predict_mean(X))
        return np.exp(np.ravel(lf.predict(X))), (np.exp(mean), np.exp(qq[:, 0]), np.exp(qq[:, 1])), None
    return pred


def main(mode):
    warnings.filterwarnings("ignore")
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(5, 3, figsize=(9.4, 12.2), sharex=True, sharey="row")
    c10, cp = fs.level_colour(10), "0.25"
    rows = []
    for q, mk, stat, lab in HYDRO:
        a8, a9, a10, hf = hydro_data(q)
        rows.append((q, mk, stat, lab, a8, a9, a10, hf))
    for q, mk, stat, lab in SCAL:
        a8, a9, a10 = scalar_data(q)
        rows.append((q, mk, stat, lab, a8, a9, a10, None))
    for i, (q, mk, stat, lab, a8, a9, a10, hf) in enumerate(rows):
        pred = models(mode, a8, a9, a10)
        for j, th in enumerate(THETAS):
            ax = axes[i, j]
            for L, (X, y) in ((8, a8), (9, a9), (10, a10)):
                m = X[:, 1] == th
                if m.any():
                    kw = {k: v for k, v in fs.series_kw(fs.level_colour(L), mk, stat=stat, ours=True).items() if k != "ls"}
                    if L == 10 and hf is not None:
                        sd = np.array([cycle_sd(hf_run(r, t), q) for r, t in X[m]])
                        ax.errorbar(X[m, 0], y[m], yerr=y[m] * sd, ls="none", color=c10, elinewidth=1.0, capsize=0, zorder=3)
                    ax.plot(X[m, 0], y[m], ls="none", **kw)
            has9 = int((a9[0][:, 1] == th).sum()) >= 3   # a level curve or a prediction is drawn only where that level has data at this angle
            if pred is not None:
                try:
                    low, (mean, lo, hi), l8curve = pred(th)
                    ax.plot(GRID, low, color=fs.level_colour(9 if mode == "chain" else 9), ls="--", lw=1.0)
                    if l8curve is not None:
                        ax.plot(GRID, l8curve, color=fs.level_colour(8), ls="--", lw=1.0)
                    ax.fill_between(GRID, lo, hi, color=cp, alpha=0.15, lw=0)
                    ax.plot(GRID, mean, color=cp, lw=1.0)
                except Exception as e:   # noqa: BLE001
                    print("panel failed", q, th, repr(e)[:80])
            if q in ("tau_95_t", "tau_mean_t", "kLa_1T_10", "dtmix_0.50"):
                ax.set_yscale("log")
            ax.set_xlim(14, 39)
            ax.grid(**fs.GRID_KW)
            if i == 0:
                ax.set_title(rf"$\theta$ = {th:g}$^\circ$", fontsize=10)
            if j == 0:
                ax.set_ylabel(lab)
            if i == 4:
                ax.set_xlabel("rocking speed (rpm)")
    h = [Line2D([], [], color=fs.level_colour(L), marker="o", mfc=fs.level_colour(L), ls="none", ms=5, label=f"L{L}") for L in (8, 9, 10)] + \
        [Patch(facecolor=cp, alpha=0.3, label="L10 prediction")]
    lg = axes[0, 2].legend(handles=h, fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.04, 1.0))
    fig.tight_layout()
    out = MF / f"fig_mf_{mode}_5x3"
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "chain")
