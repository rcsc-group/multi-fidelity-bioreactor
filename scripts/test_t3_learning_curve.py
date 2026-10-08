"""T3 (pre-registered, diary 2026-10-08): learning curve of MF vs baselines on held-out L10.

Pre-registration (diary.md, quoted):
  "T3 (zero compute): learning curve, n_HF = 2..6 (all designs with both endpoints when feasible,
   else 20 random designs), tau_mean/tau_95/EDR, MF vs HF-only GP vs LF x ratio. Descriptive; no pass
   criterion. Decides L9 vs L10 as HF for T6."
Task spec (supersedes "20"): designs = all subsets of the 9 L10 rpm containing 17.5 and 37.5 when their
count is <= 35, else 35 random such subsets (seed 0); n_HF = 2 -> the endpoints only. QoIs tau_mean_t,
tau_95_t, edr_mean_t; LF L8 and L9. Methods: MF (reference), MF-Bayes (T2, a0=2, b0=0.1), HF-only GP,
LF x ratio. Metric: median and IQR (25-75%) of held-out relative RMSE over designs; pooled 95% coverage
for MF, MF-Bayes, HF-only GP.

Counts of designs: n=2:1, 3:7, 4:21, 5:35, 6:35 -> never > 35, so the random branch (seed 0) is
implemented but not exercised.

Choices where the spec was silent:
  * Speed: the LF KRR (LooKRR, seed 42) is deterministic given the LF data, so its trained state is
    cached per (QoI, LF) and re-used (CachedLooKRR); main() asserts that cached and uncached fits give the
    same predictions (max abs diff printed). MF-Bayes is a view of the trained MF model
    (as_bayes on a shallow copy) -- same hyperparameters, as T2 verified with an independent fit.
  * Failed fits are counted and skipped for that design in all methods.
  * Coverage for MF-Bayes uses t_{0.975, nu} * scale (noise inside the amplitude, as T2).
  * Figure: one panel per QoI x LF; median as line+marker, IQR as error bar; y log scale.
Usage: uv run python scripts/test_t3_learning_curve.py
"""
from __future__ import annotations

import copy
import itertools
import json
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import figstyle as fs                                         # noqa: E402
from scripts.plot_mf_fig13 import LooKRR                                   # noqa: E402
from scripts.test_mf_l10_holdout import load, col, rel_rmse, ENDS, DS, LFS, OUT  # noqa: E402
from scripts.test_t2_bayes_band import as_bayes                            # noqa: E402
from mfbml_local.krr_lr_gpr import KernelRidgeLinearGaussianProcess        # noqa: E402

QOIS = ("tau_mean_t", "tau_95_t", "edr_mean_t")
NS = (2, 3, 4, 5, 6)
METHODS = ("MF", "MF-Bayes", "HF-only GP", "LF x ratio")
_LF_CACHE: dict = {}


class CachedLooKRR(LooKRR):
    """LooKRR whose training (deterministic given data and seed) is memoised."""

    def train(self, X, Y, portion_test=0.1):
        key = (np.asarray(X).tobytes(), np.asarray(Y).tobytes(), portion_test, self.seed)
        if key not in _LF_CACHE:
            super().train(X, Y, portion_test=portion_test)
            _LF_CACHE[key] = dict(self.__dict__)
        else:
            self.__dict__.update(_LF_CACHE[key])


def fit_cached(poly, xh, yh, xl, yl):
    lf = CachedLooKRR(design_space=DS, params_optimize=True, noise_data=True,
                      optimizer_restart=10, seed=42)
    m = KernelRidgeLinearGaussianProcess(design_space=DS, optimizer_restart=10, lf_poly_order=poly,
                                         seed=42, lf_model=lf, noise_prior=None)
    m.train(X=[col(xh), col(xl)], Y=[np.asarray(yh, float), np.asarray(yl, float)])
    return m


def bayes_view(m):
    return as_bayes(copy.copy(m))


def designs_for(n, hrpm):
    interior = [r for r in hrpm if r not in ENDS]
    combos = list(itertools.combinations(interior, n - 2))
    if len(combos) > 35:
        rng = np.random.default_rng(0)
        combos = [combos[i] for i in sorted(rng.choice(len(combos), 35, replace=False))]
    return [list(ENDS) + list(c) for c in combos]


def evaluate(xl, yl, hrpm, yh_all, tr):
    msk = np.isin(hrpm, tr)
    xh, yh, xt, yt = hrpm[msk], yh_all[msk], hrpm[~msk], yh_all[~msk]
    mf, sf = fit_cached("linear", xh, yh, xl, yl), fit_cached("ordinary", xh, yh, xl, yl)
    mum, sdm = (a.ravel() for a in mf.predict(col(xt), return_std=True))
    mus, sds = (a.ravel() for a in sf.predict(col(xt), return_std=True))
    mub, scale, nu = bayes_view(mf).predict_dist(col(xt))
    lf_at = lambda x: np.interp(x, xl, yl)
    ratio = float(np.mean(yh / lf_at(xh)))
    return dict(
        rel={"MF": rel_rmse(mum, yt), "MF-Bayes": rel_rmse(mub, yt), "HF-only GP": rel_rmse(mus, yt),
             "LF x ratio": rel_rmse(ratio * lf_at(xt), yt)},
        hits={"MF": int(np.sum(np.abs(mum - yt) <= 1.96 * sdm)),
              "MF-Bayes": int(np.sum(np.abs(mub - yt) <= scale * stats.t.ppf(0.975, nu))),
              "HF-only GP": int(np.sum(np.abs(mus - yt) <= 1.96 * sds))},
        n_test=int(len(yt)))


def main():
    t0 = time.time()
    warnings.filterwarnings("ignore")
    D = load()
    hrpm = np.array(sorted(D[10]))
    # equivalence of cached and uncached training
    from scripts.test_mf_l10_holdout import fit as fit_ref
    q0, tr0 = "tau_95_t", [17.5, 27.5, 37.5]
    xl0 = np.array(sorted(D[8])); yl0 = np.array([D[8][r][q0] for r in xl0])
    yh0 = np.array([D[10][r][q0] for r in hrpm]); k0 = np.isin(hrpm, tr0)
    a = fit_ref("linear", hrpm[k0], yh0[k0], xl0, yl0).predict(col(hrpm))
    b = fit_cached("linear", hrpm[k0], yh0[k0], xl0, yl0).predict(col(hrpm))
    b = fit_cached("linear", hrpm[k0], yh0[k0], xl0, yl0).predict(col(hrpm))   # second call = cache hit
    eq = float(np.max(np.abs(a - b)))
    print("cached vs uncached max abs diff:", eq)
    assert eq < 1e-12

    res = {}
    for q in QOIS:
        yh_all = np.array([D[10][r][q] for r in hrpm])
        for lfl in LFS:
            xl = np.array(sorted(D[lfl]))
            yl = np.array([D[lfl][r][q] for r in xl])
            for n in NS:
                rows, fails = [], 0
                for tr in designs_for(n, hrpm):
                    try:
                        rows.append(evaluate(xl, yl, hrpm, yh_all, tr))
                    except Exception as e:  # noqa: BLE001
                        fails += 1
                        print("  FAIL", q, lfl, tr, e)
                summ = {}
                for m in METHODS:
                    v = np.array([r["rel"][m] for r in rows])
                    summ[m] = dict(median=float(np.median(v)), q25=float(np.percentile(v, 25)),
                                   q75=float(np.percentile(v, 75)))
                for m in ("MF", "MF-Bayes", "HF-only GP"):
                    summ[m]["coverage95"] = sum(r["hits"][m] for r in rows) / sum(r["n_test"] for r in rows)
                res[f"{q}|L{lfl}|n{n}"] = dict(summary=summ, n_designs=len(rows), n_fail=fails,
                                               all_rel=[r["rel"] for r in rows])
    lines = [f"{'case':18s}{'n':>2s}{'#d':>4s}  " + "  ".join(f"{m:>17s}" for m in METHODS) + "   cov MF/MFB/HF"]
    for k, r in res.items():
        q, lf, n = k.split("|")
        s = r["summary"]
        lines.append(f"{q + ' ' + lf:18s}{n[1:]:>2s}{r['n_designs']:4d}  "
                     + "  ".join(f"{s[m]['median']:6.1%}[{s[m]['q25']:5.1%},{s[m]['q75']:5.1%}]" for m in METHODS)
                     + f"   {s['MF']['coverage95']:.0%}/{s['MF-Bayes']['coverage95']:.0%}/{s['HF-only GP']['coverage95']:.0%}")
    table = "\n".join(lines)
    print(table)

    plt.rcParams.update(fs.rcparams())
    f, axes = plt.subplots(3, 2, figsize=(8.0, 9.0), sharex=True)
    ch = fs.level_colour(10)
    sty = {"MF": dict(color=ch, ls="--"), "MF-Bayes": dict(color=ch, ls=(0, (6, 1.5, 1, 1.5))),
           "HF-only GP": dict(color=fs.NEUTRAL, ls="-."), "LF x ratio": dict(color=fs.NEUTRAL, ls=(0, (1, 3)))}
    ylab = {"tau_mean_t": "mean wall shear stress", "tau_95_t": "95th-percentile wall shear stress",
            "edr_mean_t": "mean dissipation rate"}
    for i, q in enumerate(QOIS):
        mk = fs.MK["ediss"] if q.startswith("edr") else fs.MK["tau"]
        stat = "max" if q == "tau_95_t" else "mean"
        for j, lfl in enumerate(LFS):
            ax = axes[i, j]
            for mi, m in enumerate(METHODS):
                x = np.array(NS) + (mi - 1.5) * 0.08
                med = np.array([res[f"{q}|L{lfl}|n{n}"]["summary"][m]["median"] * 100 for n in NS])
                lo = np.array([res[f"{q}|L{lfl}|n{n}"]["summary"][m]["q25"] * 100 for n in NS])
                hi = np.array([res[f"{q}|L{lfl}|n{n}"]["summary"][m]["q75"] * 100 for n in NS])
                kw = fs.series_kw(sty[m]["color"], mk, stat=stat, ours=True, ls=sty[m]["ls"], ms=4.5)
                ax.errorbar(x, med, yerr=[med - lo, hi - med], capsize=2, elinewidth=0.7, label=m, **kw)
            ax.set_yscale("log")
            ax.set_title(f"{ylab[q]}, LF L{lfl}", fontsize=9)
            ax.grid(**fs.GRID_KW)
            ax.set_xticks(NS)
            if j == 0:
                ax.set_ylabel("relative RMSE (%)")
            if i == 2:
                ax.set_xlabel("number of L10 runs")
            if i == 0 and j == 1:
                ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    f.tight_layout()
    f.savefig(OUT / "test_t3_learning_curve.png", dpi=200, bbox_inches="tight")
    rt = time.time() - t0
    json.dump(dict(results=res, table=table, cache_equivalence_max_abs_diff=eq, runtime_s=rt),
              open(OUT / "test_t3_learning_curve.json", "w"), indent=1)
    print(f"runtime {rt:.1f} s")


if __name__ == "__main__":
    main()
