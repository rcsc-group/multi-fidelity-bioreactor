"""TEST A (pre-registered, diary 2026-10-08): does MF beat HF-only GP and LF x ratio on held-out L10?

Data: experiments/multifidelity/hydro_dataset_v2.json. QoIs tau_mean_t, tau_95_t.
LF = L8 (variant L9), all 10 rpm. HF = L10 (9 rpm). Designs: n_HF=3 -> {17.5, 37.5}+1 interior
(7 designs); n_HF=4 -> {17.5, 37.5}+2 interior (21 designs). Score on the held-out L10 rpm.

Methods
  MF       Yi et al. KRR-LR-GPR (mfbml_local), linear basis, noise_prior=None (the method
           estimates its own noise). LF model = LooKRR exactly as plot_mf_fig13.py (LOO-CV
           hyperparameters; optimizer_restart=10, seed=42). Input = raw rpm, design space
           [15, 37.5] (the model scales it internally).
  HF-only  Same class, lf_poly_order="ordinary" (constant basis, HF points only), noise_prior=None
           (noise estimated too, same treatment as MF). Same LF model object is trained but unused
           by the basis. (Choice: the reference class's "ordinary" option instead of HeteroUK, whose
           noise must be supplied; this keeps the noise treatment identical.)
  LF x r   LF(x) * mean(HF/LF over training rpm).
  LF       LF(x) as is.
Metric: relative RMSE sqrt(mean(((pred-obs)/obs)^2)) on held-out points; median (min/max) over designs.
Secondary: pooled 95% coverage |pred-obs| <= 1.96*sd_total (sd_total includes estimated noise).
PASS: MF median < both HF-only and LF x r in >= 3 of 4 (QoI x LF) cases at n_HF=3.

Choices made where the spec was silent: representative figure design = the n_HF=3 design whose MF
rel RMSE is the median (7 designs, so the 4th ranked), chosen per panel (it can differ per case; the
rpm is stated in the output JSON). Failed fits (Cholesky errors) are recorded and skipped, counted.

Usage: uv run python scripts/test_mf_l10_holdout.py
"""
from __future__ import annotations

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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import figstyle as fs                                      # noqa: E402
from scripts.plot_mf_fig13 import LooKRR                                # noqa: E402
from mfbml_local.krr_lr_gpr import KernelRidgeLinearGaussianProcess     # noqa: E402

OUT = ROOT / "experiments/multifidelity"
DS = np.array([[15.0, 37.5]])
QOIS = ("tau_mean_t", "tau_95_t")
LFS = (8, 9)
METHODS = ("MF", "HF-only GP", "LF x ratio", "LF alone")
ENDS = (17.5, 37.5)


def load():
    rows = json.load(open(OUT / "hydro_dataset_v2.json"))
    d = {}
    for r in rows:
        d.setdefault(r["level"], {})[r["rpm"]] = r
    return d


def col(a):
    return np.asarray(a, float).reshape(-1, 1)


def fit(poly, xh, yh, xl, yl):
    lf = LooKRR(design_space=DS, params_optimize=True, noise_data=True,
                optimizer_restart=10, seed=42)
    m = KernelRidgeLinearGaussianProcess(design_space=DS, optimizer_restart=10,
                                         lf_poly_order=poly, seed=42, lf_model=lf,
                                         noise_prior=None)
    m.train(X=[col(xh), col(xl)], Y=[np.asarray(yh, float), np.asarray(yl, float)])
    return m


def pred(m, x):
    mu, sd = m.predict(col(x), return_std=True)
    return mu.ravel(), sd.ravel()


def rel_rmse(p, o):
    return float(np.sqrt(np.mean(((np.asarray(p) - np.asarray(o)) / np.asarray(o)) ** 2)))


def run_design(xl, yl, xh_all, yh_all, train_rpm):
    tr = np.isin(xh_all, train_rpm)
    xh, yh, xt, yt = xh_all[tr], yh_all[tr], xh_all[~tr], yh_all[~tr]
    mf, sf = fit("linear", xh, yh, xl, yl), fit("ordinary", xh, yh, xl, yl)
    out = {"train": list(map(float, train_rpm)), "test": list(map(float, xt))}
    lf_at = lambda x: np.interp(x, xl, yl)   # LF grid contains every HF rpm
    ratio = float(np.mean(yh / lf_at(xh)))
    mum, sdm = pred(mf, xt)
    mus, sds = pred(sf, xt)
    out["rel"] = {"MF": rel_rmse(mum, yt), "HF-only GP": rel_rmse(mus, yt),
                  "LF x ratio": rel_rmse(ratio * lf_at(xt), yt),
                  "LF alone": rel_rmse(lf_at(xt), yt)}
    out["cov_hits"] = {"MF": int(np.sum(np.abs(mum - yt) <= 1.96 * sdm)),
                       "HF-only GP": int(np.sum(np.abs(mus - yt) <= 1.96 * sds))}
    out["n_test"] = int(len(xt))
    out["noise"] = {"MF": float(mf.noise), "HF-only GP": float(sf.noise)}
    return out, (mf, sf, ratio)


def main():
    t0 = time.time()
    warnings.filterwarnings("ignore")
    D = load()
    hrpm = np.array(sorted(D[10]))
    interior = [r for r in hrpm if r not in ENDS]
    res = {}
    models = {}
    for q in QOIS:
        yh_all = np.array([D[10][r][q] for r in hrpm])
        for lfl in LFS:
            xl = np.array(sorted(D[lfl]))
            yl = np.array([D[lfl][r][q] for r in xl])
            for n in (3, 4):
                designs, fails = [], 0
                for mid in itertools.combinations(interior, n - 2):
                    tr = list(ENDS) + list(mid)
                    try:
                        o, mods = run_design(xl, yl, hrpm, yh_all, tr)
                    except Exception as e:     # noqa: BLE001
                        fails += 1
                        print(f"  FAIL {q} L{lfl} {tr}: {e}")
                        continue
                    designs.append(o)
                    models[(q, lfl, n, tuple(tr))] = mods
                key = f"{q}|L{lfl}|n{n}"
                summ = {}
                for m in METHODS:
                    v = np.array([d["rel"][m] for d in designs])
                    summ[m] = {"median": float(np.median(v)), "min": float(v.min()),
                               "max": float(v.max())}
                for m in ("MF", "HF-only GP"):
                    h = sum(d["cov_hits"][m] for d in designs)
                    t = sum(d["n_test"] for d in designs)
                    summ[m]["coverage95"] = h / t
                res[key] = {"summary": summ, "designs": designs, "n_fail": fails}

    # PASS
    wins = {}
    for q in QOIS:
        for lfl in LFS:
            s = res[f"{q}|L{lfl}|n3"]["summary"]
            wins[f"{q}|L{lfl}"] = bool(s["MF"]["median"] < s["HF-only GP"]["median"]
                                       and s["MF"]["median"] < s["LF x ratio"]["median"])
    passed = sum(wins.values()) >= 3

    # table
    lines = [f"{'case':22s}{'nHF':>4s}  " + "  ".join(f"{m:>22s}" for m in METHODS) + "  cov MF/HF"]
    for key, r in res.items():
        q, lfl, n = key.split("|")
        s = r["summary"]
        cells = "  ".join(f"{s[m]['median']:7.1%} [{s[m]['min']:5.1%},{s[m]['max']:6.1%}]" for m in METHODS)
        lines.append(f"{q + ' ' + lfl:22s}{n[1:]:>4s}  {cells}  "
                     f"{s['MF']['coverage95']:.0%}/{s['HF-only GP']['coverage95']:.0%}"
                     + (f"  fails={r['n_fail']}" if r["n_fail"] else ""))
    table = "\n".join(lines)
    print(table)
    print("MF wins (n_HF=3):", wins)
    print("PASS" if passed else "FAIL", f"({sum(wins.values())}/4 cases)")

    # representative designs and figure
    rep = {}
    plt.rcParams.update(fs.rcparams())
    f, axes = plt.subplots(2, 2, figsize=(8.0, 6.4), sharex=True)
    xg = np.linspace(15, 37.5, 200)
    for i, q in enumerate(QOIS):
        for j, lfl in enumerate(LFS):
            r = res[f"{q}|L{lfl}|n3"]
            order = np.argsort([d["rel"]["MF"] for d in r["designs"]])
            d = r["designs"][order[len(order) // 2]]
            tr = tuple(d["train"])
            rep[f"{q}|L{lfl}"] = {"train": d["train"], "mf_rel_rmse": d["rel"]["MF"]}
            mf, sf, ratio = models[(q, lfl, 3, tr)]
            xl = np.array(sorted(D[lfl]))
            yl = np.array([D[lfl][x][q] for x in xl])
            yh = np.array([D[10][x][q] for x in hrpm])
            trm = np.isin(hrpm, tr)
            mu, sd = pred(mf, xg)
            mus, _ = pred(sf, xg)
            lfg = np.interp(xg, xl, yl)
            ax = axes[i, j]
            cl, ch, mk = fs.level_colour(lfl), fs.level_colour(10), fs.MK["tau"]
            stat = "mean" if q == "tau_mean_t" else "max"
            ax.fill_between(xg, mu - 1.96 * sd, mu + 1.96 * sd, color=ch, alpha=0.18, lw=0,
                            label="MF 95%")
            ax.plot(xg, mu, color=ch, lw=1.4, ls="--", label="MF")
            ax.plot(xg, mus, color=fs.NEUTRAL, lw=1.2, ls="-.", label="HF-only GP")
            ax.plot(xg, ratio * lfg, color=fs.NEUTRAL, lw=1.2, ls=(0, (1, 3)), label="LF x ratio")
            ax.plot(xl, yl, **fs.series_kw(cl, mk, stat=stat, ours=True, ls="none"),
                    label=f"L{lfl}")
            ax.plot(hrpm[trm], yh[trm], **fs.series_kw(ch, mk, stat=stat, ours=True, ls="none"),
                    label="L10 training")
            ax.plot(hrpm[~trm], yh[~trm],
                    **fs.series_kw(ch, mk, stat=stat, ours=True, ls="none", ms=3.0, alpha=0.6),
                    label="L10 held out")
            ax.set_title(f"LF L{lfl}", fontsize=10)
            ax.grid(**fs.GRID_KW)
            if j == 0:
                ax.set_ylabel(("mean wall shear stress (Pa)" if i == 0 else
                               r"95th-percentile wall shear stress (Pa)"))
            if i == 1:
                ax.set_xlabel("rocking speed (rpm)")
            if i == 0 and j == 1:
                ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    f.tight_layout()
    f.savefig(OUT / "test_mf_l10_holdout.png", dpi=200, bbox_inches="tight")
    plt.close(f)

    f, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), sharey=True)
    cases = [(q, l) for q in QOIS for l in LFS]
    cl_ = {"MF": fs.level_colour(10), "HF-only GP": fs.NEUTRAL, "LF x ratio": "0.2", "LF alone": "0.7"}
    mks = {"MF": "o", "HF-only GP": "s", "LF x ratio": "^", "LF alone": "v"}
    for k, n in enumerate((3, 4)):
        ax = axes[k]
        for mi, m in enumerate(METHODS):
            xs = np.arange(len(cases)) + (mi - 1.5) * 0.15
            med = [res[f"{q}|L{l}|n{n}"]["summary"][m]["median"] * 100 for q, l in cases]
            lo = [res[f"{q}|L{l}|n{n}"]["summary"][m]["min"] * 100 for q, l in cases]
            hi = [res[f"{q}|L{l}|n{n}"]["summary"][m]["max"] * 100 for q, l in cases]
            ax.errorbar(xs, med, yerr=[np.array(med) - lo, np.array(hi) - med], fmt=mks[m],
                        color=cl_[m], ms=5, lw=0.8, capsize=2, label=m)
        ax.set_xticks(range(len(cases)))
        ax.set_xticklabels([("mean" if q == "tau_mean_t" else "95th") + f"\nLF L{l}" for q, l in cases],
                           fontsize=8)
        ax.set_yscale("log")
        ax.set_title(rf"$n_{{HF}}={n}$", fontsize=10)
        ax.grid(**fs.GRID_KW)
        if k == 0:
            ax.set_ylabel("relative RMSE, median over designs (%)")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    f.tight_layout()
    f.savefig(OUT / "test_mf_l10_holdout_summary.png", dpi=200, bbox_inches="tight")

    rt = time.time() - t0
    json.dump({"pass": passed, "mf_wins_n3": wins, "representative_designs_n3": rep,
               "results": res, "table": table, "runtime_s": rt},
              open(OUT / "test_mf_l10_holdout.json", "w"), indent=1)
    print("representative n3 designs:", rep)
    print(f"runtime {rt:.1f} s")


if __name__ == "__main__":
    main()
