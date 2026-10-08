"""T2b (pre-registered, diary 2026-10-08): fully Bayesian KRR-LR-GPR, re-run of test A.

Pre-registration (diary.md, quoted):
  "T2b PRE-REGISTERED (own follow-up of the T2 root cause; T2 approved by the user): fully Bayesian
   KRR-LR-GPR residual GP.
   - Keep: LF KRR, basis [1, f_l], rho flat, amplitude IG(2, 0.1) integrated analytically (as T2).
   - New: integrate the kernel length scale and the noise ratio on a 2-D grid (60 x 40) with priors:
     length scale in scaled input units log l ~ N(log 0.3, 0.75^2); noise ratio (noise sd / HF std)
     log-uniform on [1e-3, 0.3]. Marginal likelihood per grid node in closed form:
     |K|^-1/2 |F'K^-1F|^-1/2 (b0 + S/2)^-(a0 + (n-p)/2). Predictive = mixture of Student-t over the grid.
   - The SAME priors and machinery for the HF-only GP (constant basis), so the baseline is no longer
     a collapsed MLE fit.
   - PASS = on test A (tau_mean, tau_95; LF L8, L9; n_HF = 3): pooled 95% coverage >= 85% AND MF
     median rel RMSE <= 1.2x the test-A value in every case. Reported, not criteria: the fair
     HF-only GP, T1 (EDR), and T3 at n_HF = 2..6."

Model, derivation and the l <-> RBF theta conversion: scripts/mfbml_local/krr_lr_gpr_fullbayes.py.
  theta = log10(1/(2 l^2)) because RBF is exp(-10^theta d^2) on inputs scaled to [0,1].

Protocol = test A (scripts/test_mf_l10_holdout.py) for tau_mean_t, tau_95_t, LF L8/L9, n_HF = 3 (the
two end rpm + one interior = 7 designs), held-out L10 rpm scored. Metric: relative RMSE of the
predictive mean (weighted mean of node locations); coverage: hit if the observation lies in the
central 95% interval of the mixture (quantiles 0.025/0.975 by brentq CDF inversion).
PASS reference value = MF "median" in experiments/multifidelity/test_mf_l10_holdout.json
(key "<qoi>|L<lf>|n3").  Pooled coverage = hits / held-out points over all designs of the 4 cases.

Reported, not criteria: HF-only full-Bayes GP and LF x ratio (same designs); T1 (edr_mean_t, LF L8/L9,
n_HF = 3); T3 learning curve n_HF = 2..6 for tau_mean_t, tau_95_t, edr_mean_t with designs_for() of
scripts/test_t3_learning_curve.py (all subsets containing both endpoints; <= 35 per n).

Choices where the spec was silent (closest to the spec):
  * LF KRR: LooKRR (seed 42, optimizer_restart 10, portion_test 0.2 = reference lf_portion) trained
    once per (QoI, LF level) on all 10 LF rpm and reused for every design (the LF data do not depend
    on the design; T3 verified the training is deterministic).
  * HF-only GP = same class with basis "ordinary" ([1] only), the same grid and priors.
  * Input scaling: raw rpm, design space [15, 37.5], as test A.
  * Grid: l log-spaced on exp(log 0.3 +- 3.29 * 0.75) (central 99.9%), 60 nodes, grid weight = prior
    density at the node; eta log-spaced on [1e-3, 0.3], 40 nodes, equal prior weight.
  * Predictive includes the noise term eta^2 inside the amplitude-scaled variance (as T2).
  * Coverage band = central 95% mixture quantiles (not a symmetric +- band).
  * Hyperparameter posterior reported for the median-error test-A design of each case (marginal
    median and 90% central interval of l in scaled units and of eta).
  * Figure (a): y axis limited to the data range (L10 and LF) extended by 40% of the span so a wide
    HF-only band does not squash the data; bands are clipped by the axes. Representative design =
    the 4th-ranked of 7 by full-Bayes MF rel RMSE. HF-only band is a NEUTRAL fill.
  * Figure (b): rows = QoI, columns = (rel RMSE, coverage) for each LF level; RMSE = median over
    designs with IQR bars (as T3), coverage = pooled over designs, 0.95 shown as a NEUTRAL guide.

Usage: uv run python scripts/test_t2b_fullbayes.py
"""
from __future__ import annotations

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
from scripts import figstyle as fs                                            # noqa: E402
from scripts.test_mf_l10_holdout import load, col, rel_rmse, ENDS, DS, LFS, OUT  # noqa: E402
from scripts.test_t3_learning_curve import CachedLooKRR, designs_for          # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP                 # noqa: E402

QOIS_A = ("tau_mean_t", "tau_95_t")
QOIS_T3 = ("tau_mean_t", "tau_95_t", "edr_mean_t")
NS = (2, 3, 4, 5, 6)
LF_PORTION = 0.2


def train_lf(xl, yl):
    lf = CachedLooKRR(design_space=DS, params_optimize=True, noise_data=True,
                      optimizer_restart=10, seed=42)
    lf.train(col(xl), np.asarray(yl, float), portion_test=LF_PORTION)
    return lf


def evaluate(lf, xl, yl, hrpm, yh_all, tr, keep=False):
    msk = np.isin(hrpm, tr)
    xh, yh, xt, yt = hrpm[msk], yh_all[msk], hrpm[~msk], yh_all[~msk]
    out = dict(train=list(map(float, tr)), n_test=int(len(yt)), rel={}, hits={}, post={})
    models = {}
    for name, basis in (("MF", "linear"), ("HF-only GP", "ordinary")):
        m = FullBayesKRRLRGP(DS, lf_model=lf, basis=basis).train(col(xh), yh)
        mu = m.predict_mean(col(xt))
        qq = m.quantile(col(xt), [0.025, 0.975])
        out["rel"][name] = rel_rmse(mu, yt)
        out["hits"][name] = int(np.sum((yt >= qq[:, 0]) & (yt <= qq[:, 1])))
        out["post"][name] = m.posterior_summary(0.9)
        models[name] = m
    ratio = float(np.mean(yh / np.interp(xh, xl, yl)))
    out["rel"]["LF x ratio"] = rel_rmse(ratio * np.interp(xt, xl, yl), yt)
    models["ratio"] = ratio
    return out, models


def summarise(rows, methods=("MF", "HF-only GP", "LF x ratio")):
    s = {}
    for m in methods:
        v = np.array([r["rel"][m] for r in rows])
        s[m] = dict(median=float(np.median(v)), q25=float(np.percentile(v, 25)),
                    q75=float(np.percentile(v, 75)))
    for m in ("MF", "HF-only GP"):
        s[m]["coverage95"] = sum(r["hits"][m] for r in rows) / sum(r["n_test"] for r in rows)
    return s


def run_case(D, hrpm, q, lfl, n, lf_cache):
    key = (q, lfl)
    xl = np.array(sorted(D[lfl]))
    yl = np.array([D[lfl][r][q] for r in xl])
    if key not in lf_cache:
        lf_cache[key] = train_lf(xl, yl)
    yh_all = np.array([D[10][r][q] for r in hrpm])
    rows, mods, fails = [], {}, 0
    for tr in designs_for(n, hrpm):
        try:
            o, m = evaluate(lf_cache[key], xl, yl, hrpm, yh_all, tr)
        except Exception as e:  # noqa: BLE001
            fails += 1
            print("  FAIL", q, lfl, tr, e)
            continue
        rows.append(o)
        mods[tuple(tr)] = m
    return dict(summary=summarise(rows), rows=rows, n_fail=fails), mods


def main():
    t0 = time.time()
    warnings.filterwarnings("ignore")
    D = load()
    hrpm = np.array(sorted(D[10]))
    old = json.load(open(OUT / "test_mf_l10_holdout.json"))["results"]
    lf_cache = {}

    # ---------------- test A (the pre-registered criterion)
    resA, modsA, hits_t, n_t = {}, {}, 0, 0
    for q in QOIS_A:
        for lfl in LFS:
            r, mods = run_case(D, hrpm, q, lfl, 3, lf_cache)
            ref = old[f"{q}|L{lfl}|n3"]["summary"]["MF"]["median"]
            s = r["summary"]
            r["testA_MF_median"] = ref
            r["ratio_to_testA"] = s["MF"]["median"] / ref
            hits_t += sum(d["hits"]["MF"] for d in r["rows"])
            n_t += sum(d["n_test"] for d in r["rows"])
            order = np.argsort([d["rel"]["MF"] for d in r["rows"]])
            rep = r["rows"][order[len(order) // 2]]
            r["representative"] = dict(train=rep["train"], post=rep["post"], rel=rep["rel"])
            resA[f"{q}|L{lfl}"] = r
            modsA[(q, lfl)] = mods
    pooled = hits_t / n_t
    cov_ok = pooled >= 0.85
    rmse_ok = all(r["ratio_to_testA"] <= 1.2 for r in resA.values())
    passed = bool(cov_ok and rmse_ok)

    lines = [f"{'case':18s} {'MF':>7s} {'testA MF':>9s} {'ratio':>6s} {'HF-only':>8s} {'LFxr':>7s} "
             f"{'cov MF':>7s} {'cov HF':>7s}"]
    for k, r in resA.items():
        s = r["summary"]
        lines.append(f"{k:18s} {s['MF']['median']:7.1%} {r['testA_MF_median']:9.1%} {r['ratio_to_testA']:6.2f} "
                     f"{s['HF-only GP']['median']:8.1%} {s['LF x ratio']['median']:7.1%} "
                     f"{s['MF']['coverage95']:7.0%} {s['HF-only GP']['coverage95']:7.0%}")
    lines.append(f"pooled MF coverage {pooled:.1%} ({hits_t}/{n_t})")
    table = "\n".join(lines)
    print(table)
    print("PASS" if passed else "FAIL", f"(coverage>=85%: {cov_ok}; rmse ratio<=1.2 in all: {rmse_ok})")
    for k, r in resA.items():
        for m in ("MF", "HF-only GP"):
            p = r["representative"]["post"][m]
            print(f"  post {k} {m} design {r['representative']['train']}: l {p['l_med']:.3f} "
                  f"[{p['l_lo']:.3f},{p['l_hi']:.3f}], eta {p['eta_med']:.3g} [{p['eta_lo']:.3g},{p['eta_hi']:.3g}]")

    # ---------------- T1 (EDR, n_HF = 3), reported
    resT1 = {}
    for lfl in LFS:
        r, _ = run_case(D, hrpm, "edr_mean_t", lfl, 3, lf_cache)
        resT1[f"edr_mean_t|L{lfl}"] = r
        s = r["summary"]
        print(f"T1 edr L{lfl}: MF {s['MF']['median']:.1%}  HF-only {s['HF-only GP']['median']:.1%}  "
              f"ratio {s['LF x ratio']['median']:.1%}  cov MF {s['MF']['coverage95']:.0%} "
              f"HF {s['HF-only GP']['coverage95']:.0%}")

    # ---------------- T3 learning curve, reported
    resLC = {}
    for q in QOIS_T3:
        for lfl in LFS:
            for n in NS:
                if n == 3 and q in resA_keys(resA, q, lfl):
                    resLC[f"{q}|L{lfl}|n3"] = dict(summary=resA[f"{q}|L{lfl}"]["summary"],
                                                   n_designs=len(resA[f"{q}|L{lfl}"]["rows"]))
                    continue
                r, _ = run_case(D, hrpm, q, lfl, n, lf_cache)
                resLC[f"{q}|L{lfl}|n{n}"] = dict(summary=r["summary"], n_designs=len(r["rows"]), n_fail=r["n_fail"])
    lc = [f"{'case':18s}{'n':>2s}  {'MF':>20s} {'HF-only':>20s} {'LF x ratio':>20s}   cov MF/HF"]
    for k, r in resLC.items():
        q, lf, n = k.split("|")
        s = r["summary"]
        lc.append(f"{q + ' ' + lf:18s}{n[1:]:>2s}  " + " ".join(
            f"{s[m]['median']:7.1%}[{s[m]['q25']:5.1%},{s[m]['q75']:5.1%}]" for m in ("MF", "HF-only GP", "LF x ratio"))
            + f"   {s['MF']['coverage95']:.0%}/{s['HF-only GP']['coverage95']:.0%}")
    lc_table = "\n".join(lc)
    print(lc_table)

    # ---------------- figure (a)
    plt.rcParams.update(fs.rcparams())
    f, axes = plt.subplots(2, 2, figsize=(8.0, 6.4), sharex=True)
    xg = np.linspace(15, 37.5, 120)
    for i, q in enumerate(QOIS_A):
        for j, lfl in enumerate(LFS):
            r = resA[f"{q}|L{lfl}"]
            tr = tuple(r["representative"]["train"])
            m = modsA[(q, lfl)][tr]
            mf, hf, ratio = m["MF"], m["HF-only GP"], m["ratio"]
            xl = np.array(sorted(D[lfl]))
            yl = np.array([D[lfl][x][q] for x in xl])
            yh = np.array([D[10][x][q] for x in hrpm])
            trm = np.isin(hrpm, tr)
            mu, qm = mf.predict_mean(col(xg)), mf.quantile(col(xg), [0.025, 0.975])
            muh, qh = hf.predict_mean(col(xg)), hf.quantile(col(xg), [0.025, 0.975])
            ax = axes[i, j]
            cl, ch, mk = fs.level_colour(lfl), fs.level_colour(10), fs.MK["tau"]
            stat = "mean" if q == "tau_mean_t" else "max"
            ax.fill_between(xg, qh[:, 0], qh[:, 1], color=fs.NEUTRAL, alpha=0.12, lw=0, label="HF-only GP 95%")
            ax.fill_between(xg, qm[:, 0], qm[:, 1], color=ch, alpha=0.18, lw=0, label="MF 95%")
            ax.plot(xg, mu, color=ch, lw=1.4, ls="--", label="MF")
            ax.plot(xg, muh, color=fs.NEUTRAL, lw=1.2, ls="-.", label="HF-only GP")
            ax.plot(xg, ratio * np.interp(xg, xl, yl), color="0.2", lw=1.2, ls=(0, (1, 3)), label="LF x ratio")
            ax.plot(xl, yl, **fs.series_kw(cl, mk, stat=stat, ours=True, ls="none"), label=f"L{lfl}")
            ax.plot(hrpm[trm], yh[trm], **fs.series_kw(ch, mk, stat=stat, ours=True, ls="none"),
                    label="L10 training")
            ax.plot(hrpm[~trm], yh[~trm],
                    **fs.series_kw(ch, mk, stat=stat, ours=True, ls="none", ms=3.0, alpha=0.6),
                    label="L10 held out")
            allv = np.concatenate([yl, yh])
            span = allv.max() - allv.min()
            ax.set_ylim(allv.min() - 0.4 * span, allv.max() + 0.4 * span)
            ax.set_title(f"LF L{lfl}", fontsize=10)
            ax.grid(**fs.GRID_KW)
            if j == 0:
                ax.set_ylabel("mean wall shear stress (Pa)" if i == 0 else "95th-percentile wall shear stress (Pa)")
            if i == 1:
                ax.set_xlabel("rocking speed (rpm)")
            if i == 0 and j == 1:  # one legend for both LF columns, so L8 and L9 are both labelled
                hl = dict(zip(*axes[0, 0].get_legend_handles_labels()[::-1]))
                hl.update(zip(*ax.get_legend_handles_labels()[::-1]))
                ax.legend(list(hl.values()), list(hl), loc="upper left", bbox_to_anchor=(1.02, 1.0),
                          frameon=False, fontsize=8)
    f.tight_layout()
    f.savefig(OUT / "test_t2b_fullbayes.png", dpi=200, bbox_inches="tight")
    plt.close(f)

    # ---------------- figure (b)
    f, axes = plt.subplots(3, 4, figsize=(12.5, 8.6), sharex=True)
    ch = fs.level_colour(10)
    sty = {"MF": dict(color=ch, ls="--"), "HF-only GP": dict(color=fs.NEUTRAL, ls="-."),
           "LF x ratio": dict(color="0.2", ls=(0, (1, 3)))}
    ylab = {"tau_mean_t": "mean wall shear stress", "tau_95_t": "95th-percentile wall shear stress",
            "edr_mean_t": "mean dissipation rate"}
    for i, q in enumerate(QOIS_T3):
        mk = fs.MK["ediss"] if q.startswith("edr") else fs.MK["tau"]
        stat = "max" if q == "tau_95_t" else "mean"
        for j, lfl in enumerate(LFS):
            axr, axc = axes[i, 2 * j], axes[i, 2 * j + 1]
            for mi, m in enumerate(sty):
                x = np.array(NS) + (mi - 1) * 0.08
                S = [resLC[f"{q}|L{lfl}|n{n}"]["summary"][m] for n in NS]
                med = np.array([s["median"] * 100 for s in S])
                lo = np.array([s["q25"] * 100 for s in S])
                hi = np.array([s["q75"] * 100 for s in S])
                kw = fs.series_kw(sty[m]["color"], mk, stat=stat, ours=True, ls=sty[m]["ls"], ms=4.5)
                axr.errorbar(x, med, yerr=[med - lo, hi - med], capsize=2, elinewidth=0.7, label=m, **kw)
                if m != "LF x ratio":
                    axc.plot(x, [s["coverage95"] for s in S], **kw)
            axc.axhline(0.95, color=fs.NEUTRAL, lw=0.8, ls="-")
            axc.set_ylim(0, 1.02)
            axr.set_yscale("log")
            axr.set_title(f"{ylab[q]}, LF L{lfl}", fontsize=9)
            for ax in (axr, axc):
                ax.grid(**fs.GRID_KW)
                ax.set_xticks(NS)
            if j == 0:
                axr.set_ylabel("relative RMSE (%)")
            axc.set_ylabel("95% coverage (fraction)")
            if i == 2:
                axr.set_xlabel("number of L10 runs")
                axc.set_xlabel("number of L10 runs")
            if i == 0 and j == 1:
                axc.legend(*axr.get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(1.02, 1.0),
                           frameon=False, fontsize=8)
    f.tight_layout()
    f.savefig(OUT / "test_t2b_learning_curve.png", dpi=200, bbox_inches="tight")
    plt.close(f)

    rt = time.time() - t0
    json.dump(dict(pass_=passed, coverage_ok=cov_ok, rmse_ok=rmse_ok, pooled_cov_mf=pooled,
                   cases=resA, t1_edr=resT1, table=table, runtime_s=rt),
              open(OUT / "test_t2b_fullbayes.json", "w"), indent=1)
    json.dump(dict(results=resLC, table=lc_table, runtime_s=rt),
              open(OUT / "test_t2b_learning_curve.json", "w"), indent=1)
    print(f"runtime {rt:.1f} s")


def resA_keys(resA, q, lfl):
    """Names of QoIs of test A available for re-use at n = 3 (EDR is not in test A)."""
    return (q,) if f"{q}|L{lfl}" in resA else ()


if __name__ == "__main__":
    main()
