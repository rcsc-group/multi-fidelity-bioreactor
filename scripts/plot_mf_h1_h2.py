"""Figures for the two cheaper-replication hypotheses (rpm sweep, 7 deg).

H1 (time horizon): predict the LATE threshold from the EARLY one, same level (L9).
   panels: kLa 10% -> 50%, dtmix chi 0.50 -> 0.95.
H2 (grid): predict L9 from L7, same threshold, on Kim's six numbers.
   panels: kLa at C* = 10/25/50%, dtmix at chi = 0.50/0.75/0.95.
Each panel shows ONE training set (the lowest and highest available rpm + 25). The held-out points are
drawn (x markers), because judging the test needs them. Errors over all
training sets are in diag_mf_early_late.py / diag_mf_l7_l9_kim.py.
Model as in the Fig 13 MF panels (LOO-KRR on LF + universal kriging with
HF sigma = 5% of y, REML); baselines: HF-only GP and LF x mean ratio.
Kim's values are drawn where the panel's quantity is one of his.

Usage: uv run python scripts/plot_mf_h1_h2.py h1|h2
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs                              # noqa: E402
from scripts.plot_mf_fig13 import LooKRR                        # noqa: E402
from scripts.plot_mf_fig13_unc import HeteroUK                  # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
REL_NOISE = 0.05
KIM = pd.read_csv(ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_frequency.csv").set_index("RPM")
KIM_COL = {"kLa_1T_10": "kLa_exp5pts_10", "kLa_1T_25": "kLa_exp5pts_25",
           "kLa_1T_50": "kLa_exp5pts_50", "dtmix_0.50": "dtmix_strict_0.5",
           "dtmix_0.75": "dtmix_strict_0.75", "dtmix_0.95": "dtmix_strict_0.95"}
LABEL = {"kLa_1T_10": r"$k_La$ at $C^*=10\%$", "kLa_1T_25": r"$k_La$ at $C^*=25\%$",
         "kLa_1T_50": r"$k_La$ at $C^*=50\%$", "dtmix_0.50": r"$\Delta t_{\mathrm{mix}}$ at $\chi=0.50$",
         "dtmix_0.75": r"$\Delta t_{\mathrm{mix}}$ at $\chi=0.75$",
         "dtmix_0.95": r"$\Delta t_{\mathrm{mix}}$ at $\chi=0.95$"}
UNIT = {"k": r"(h$^{-1}$)", "d": "(s)"}


def val(run, key):
    f = ROOT / "runs" / run / "results.json"
    return json.load(open(f)).get(key, math.nan) if f.exists() else math.nan


def fit(x, lf, hf, train):
    lo, hi = x.min(), x.max()
    krr = LooKRR(design_space=np.array([[lo, hi]]), params_optimize=True,
                 noise_data=True, optimizer_restart=10, seed=42)
    krr.train(x.reshape(-1, 1), lf)
    f_lf = lambda z: np.ravel(krr.predict(np.asarray(z, float).reshape(-1, 1)))
    sig = REL_NOISE * hf[train]
    mf = HeteroUK(lo, hi, lambda z: np.column_stack([np.ones_like(z), f_lf(z)])).fit(x[train], hf[train], sig)
    sf = HeteroUK(lo, hi, lambda z: np.ones((len(z), 1))).fit(x[train], hf[train], sig)
    ratio = float(np.mean(hf[train] / lf[train]))
    return f_lf, mf, sf, ratio


def panel(ax, x, lf, hf, lf_lab, hf_lab, c_lf, c_hf, kim_key):
    # endpoints of the AVAILABLE data plus 25 rpm (a missing endpoint would
    # otherwise leave two training points and force extrapolation)
    train = np.isin(x, [x.min(), 25.0, x.max()])
    f_lf, mf, sf, ratio = fit(x, lf, hf, train)
    xg = np.linspace(x.min(), x.max(), 200)
    mu, sd = mf.predict(xg)
    mu_s, _ = sf.predict(xg)
    ax.fill_between(xg, np.maximum(mu - 1.96 * sd, 1e-6), mu + 1.96 * sd,
                    color=c_hf, alpha=0.15, lw=0, label="MF 95%")
    ax.plot(xg, mu, color=c_hf, lw=1.4, ls="--", label="MF")
    ax.plot(x, lf, ls="none", marker="o", color=c_lf, ms=4, label=lf_lab)
    ax.errorbar(x[train], hf[train], yerr=1.96 * REL_NOISE * hf[train], fmt="o",
                color=c_hf, ms=5, capsize=2, label=f"{hf_lab}, training")
    ax.plot(x[~train], hf[~train], ls="none", marker="x", color=c_hf, ms=6,
            mew=1.4, label=f"{hf_lab}, held out")
    vals = np.concatenate([lf, hf])
    if kim_key in KIM_COL:
        k = KIM[KIM_COL[kim_key]].reindex(RPMS)
        ax.plot(k.index, k.values, color=fs.KIM, lw=fs.LW_KIM, label="Kim et al.")
        vals = np.concatenate([vals, k.dropna().values])
    ax.set_yscale("log")
    # limits from the data, so a wide band cannot drag the axis decades down
    ax.set_ylim(0.5 * vals.min(), 2.0 * vals.max())
    ax.grid(**fs.GRID_KW)


def data(pairs):
    out = []
    for key_lf, run_lf, key_hf, run_hf in pairs:
        rows = [(r, val(run_lf.format(r), key_lf), val(run_hf.format(r), key_hf)) for r in RPMS]
        a = np.array([q for q in rows if all(np.isfinite(q))], float)
        out.append(a)
    return out


def main() -> None:
    which = sys.argv[1]
    plt.rcParams.update(fs.rcparams())
    if which == "h1":
        pairs = [("kLa_1T_10", "fig9_l9_rpm{:g}", "kLa_1T_50", "fig9_l9_rpm{:g}"),
                 ("dtmix_0.50", "fig9_l9_rpm{:g}", "dtmix_0.95", "fig9_l9_rpm{:g}")]
        titles = [("early", "late", "kLa_1T_50", r"$C^*$: 10% $\to$ 50%"),
                  ("early", "late", "dtmix_0.95", r"$\chi$: 0.50 $\to$ 0.95")]
        fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.9))
        for ax, a, (l, h, key, ttl) in zip(axes, data(pairs), titles):
            panel(ax, a[:, 0], a[:, 1], a[:, 2], l, h, "0.6", fs.level_colour(9), key)
            ax.set_title(ttl, loc="left", fontsize=10)
            ax.set_ylabel(LABEL[key].split(" at")[0] + " " + UNIT[key[0]])
            ax.set_xlabel("rocking frequency [rpm]")
        fig.suptitle("L9", x=0.02, ha="left")
        out = "mf_h1_early_to_late.png"
    else:
        keys = ["kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50", "dtmix_0.75", "dtmix_0.95"]
        lfl = int(sys.argv[2]) if len(sys.argv) > 2 else 7
        pairs = [(k, f"kmix_l{lfl}_rpm{{:g}}", k, "fig9_l9_rpm{:g}") for k in keys]
        fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.0))
        for ax, a, key in zip(axes.ravel(), data(pairs), keys):
            if len(a) < 5:
                ax.set_title(f"{LABEL[key]}: not enough data", loc="left", fontsize=9)
                continue
            panel(ax, a[:, 0], a[:, 1], a[:, 2], f"L{lfl}", "L9", fs.level_colour(lfl),
                  fs.level_colour(9), key)
            ax.set_title(LABEL[key], loc="left", fontsize=10)
            ax.set_ylabel(UNIT[key[0]])
        for ax in axes[-1]:
            ax.set_xlabel("rocking frequency [rpm]")
        out = f"mf_h2_l{lfl}_to_l9.png"
    h, lab = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper left", bbox_to_anchor=(1.0, 0.95), frameon=False, fontsize=8)
    fig.tight_layout()
    p = ROOT / "experiments/multifidelity" / out
    fig.savefig(p, dpi=170, bbox_inches="tight")
    print("saved", p)


if __name__ == "__main__":
    main()
