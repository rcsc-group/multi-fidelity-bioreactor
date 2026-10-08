"""T7 (pre-registered, diary 2026-10-08): can mfbml carry dtmix and kLa from L8 to L9?

Pre-registration (diary.md, quoted verbatim):
  "- T7 PRE-REGISTERED (zero compute): test A protocol with L8 -> L9 as a proxy for L8 -> L10, on dtmix and kLa.
  - QoIs: dtmix_0.50, dtmix_0.75, dtmix_0.95, kLa_1T_10, kLa_1T_25, kLa_1T_50 (whole-period kLa, as Figs 11/12).
  - LF = kmix_l8_rpm* (10 rpm). HF = fig9_l9_rpm* (10 rpm); dtmix_0.95 at 37.5 from fig9_l9_rpm37.5_ext1; kLa_1T_50 at 30 rpm is nan -> that rpm is dropped for that QoI only.
  - Fit on log y (positive QoIs that span up to 60x); score on exp of the predictive median. Same metric otherwise (rel RMSE over held-out L9 rpm; 95% coverage from mixture quantiles).
  - Model: T2b full-Bayes MF (FullBayesKRRLRGP, basis linear) and the same machinery with basis ordinary (HF-only GP); LF x ratio (ratio = mean of yh/yl over training rpm) as the honest baseline. n_HF = 3, designs = both endpoints (15, 37.5) + one interior.
  - PASS per QoI: MF median rel RMSE < LF x ratio AND < HF-only GP. Overall pass = >= 4 of 6 QoIs pass AND pooled MF coverage >= 85%.
  - Reported, not criteria: n_HF = 2..6 learning curve; the raw L8 vs L9 values (does L8 even rank the rpm like L9?), Spearman rank correlation L8 vs L9 per QoI.
  - Caveat fixed now: a pass says mfbml can carry dtmix/kLa from L8 to L9, not to L10. L8/L9/L10 dtmix_0.95 at 32.5 rpm = 68/128/191 s: not converged, so L9 is not a stand-in for L10 in value."

Implementation notes from the task brief (not silent): data paths runs/kmix_l8_rpm<r>/results.json (LF) and
runs/fig9_l9_rpm<r>/results.json (HF), r = 15, 17.5, ..., 37.5 formatted with :g; designs use endpoints
(15, 37.5), all interior combinations when <= 35 else 35 random (seed 0), restricted to the rpm
available for that QoI; predictive median / 95% quantiles are exp of the log-space mixture quantiles
(exact for a monotone map); LF x ratio uses ratio = mean(yh / yl_interp) on RAW values over training rpm.

Choices where the spec was silent (closest to scripts/test_t2b_fullbayes.py):
  * LF KRR: CachedLooKRR (seed 42, optimizer_restart 10, portion_test 0.2), trained once per QoI on the
    log of all valid LF rpm of that QoI (an LF nan would be dropped; none is expected).
  * Input scaling: raw rpm, design space [15, 37.5] (DS of test A); grid, priors, a0, b0 = T2b defaults.
  * Rel RMSE is computed on raw (back-transformed) values: sqrt(mean(((exp(median)-y)/y)^2)).
  * Coverage hit: observation inside [exp q0.025, exp q0.975]; pooled = hits / held-out points over all
    designs and all 6 QoIs at n_HF = 3.
  * HF-only GP = same class, basis "ordinary", same grid and priors (as T2b).
  * "median over designs" and IQR (25-75%) over designs, as T2b/T3.
  * Failed fits are counted, printed and skipped for that design (as T2b).
  * Representative design for figure (a) = the design ranked len//2 by MF rel RMSE (4th of 7 as T2b, 0-based).
  * When a QoI loses an rpm (nan), that rpm is removed from LF, from the LF-x-ratio interpolation
    (LF lookup is np.interp over the valid LF rpm), from HF and from the Spearman correlation.
  * HF dtmix_0.95 at 37.5 from fig9_l9_rpm37.5_ext1 for that QoI only; every other QoI at 37.5 uses
    fig9_l9_rpm37.5. Both values are printed in the table.
  * Spearman: scipy.stats.spearmanr(L8, L9) on the rpm available for the QoI.
  * Markers: each panel uses the figstyle threshold marker of its own QoI (MK["chi_0.50/0.75/0.95"],
    MK["cstar_10/25/50"]).
    Fill: hollow (mean) for every QoI, because neither mixing time nor whole-period kLa is an
    extremum.
  * Axis units: dtmix in s (collect_results.py: "dimensional mixing time (s)"); kLa in 1/h
    (plot_fig11_fig12.py y label k_La (h^-1)).
  * Figures: legends outside axes (right of the top-right panel in (a); right of the top-right panel
    in (b)); titles = the quantity only; MF line in the L9 colour dashed with its 95% band; LF x ratio
    dotted dark grey; HF-only GP and the baselines in (b) drawn as in T2b (MF = L9 colour dashed,
    HF-only GP neutral dash-dot, ratio dark-grey dotted). Learning curve: x offset +-0.08 per method.

Usage: uv run python scripts/test_t7_mix_kla_l8_l9.py
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
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import figstyle as fs                                  # noqa: E402
from scripts.test_mf_l10_holdout import col, rel_rmse, DS, OUT      # noqa: E402
from scripts.test_t3_learning_curve import CachedLooKRR             # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP       # noqa: E402

QOIS = ("dtmix_0.50", "dtmix_0.75", "dtmix_0.95", "kLa_1T_10", "kLa_1T_25", "kLa_1T_50")
RPMS = [15 + 2.5 * i for i in range(10)]
ENDS7 = (15.0, 37.5)
NS = (2, 3, 4, 5, 6)
LF_PORTION = 0.2
Q3 = [0.025, 0.5, 0.975]
METHODS = ("MF", "HF-only GP", "LF x ratio")
TITLE = {"dtmix_0.50": r"$\Delta t_{mix}$, $\chi=0.50$", "dtmix_0.75": r"$\Delta t_{mix}$, $\chi=0.75$",
         "dtmix_0.95": r"$\Delta t_{mix}$, $\chi=0.95$", "kLa_1T_10": r"$k_La$, $C^*=10\%$",
         "kLa_1T_25": r"$k_La$, $C^*=25\%$", "kLa_1T_50": r"$k_La$, $C^*=50\%$"}
YLAB = {"dtmix": "mixing time (s)", "kLa": r"$k_La$ (h$^{-1}$)"}
MARK = {"dtmix_0.50": fs.MK["chi_0.50"], "dtmix_0.75": fs.MK["chi_0.75"], "dtmix_0.95": fs.MK["chi_0.95"],
        "kLa_1T_10": fs.MK["cstar_10"], "kLa_1T_25": fs.MK["cstar_25"], "kLa_1T_50": fs.MK["cstar_50"]}


# ----------------------------------------------------------------- data
def load_t7():
    lf, hf = {}, {}
    for r in RPMS:
        lf[r] = json.load(open(ROOT / f"runs/kmix_l8_rpm{r:g}/results.json"))
        hf[r] = json.load(open(ROOT / f"runs/fig9_l9_rpm{r:g}/results.json"))
    ext = json.load(open(ROOT / "runs/fig9_l9_rpm37.5_ext1/results.json"))
    return lf, hf, ext


def series(q, lf, hf, ext):
    """Valid (rpm, yl, yh) arrays for QoI q (nan rows dropped)."""
    rows = []
    for r in RPMS:
        yh = ext[q] if (q == "dtmix_0.95" and r == 37.5) else hf[r][q]
        rows.append((r, float(lf[r][q]), float(yh)))
    a = np.array(rows)
    ok = np.isfinite(a[:, 1]) & np.isfinite(a[:, 2])
    return a[ok, 0], a[ok, 1], a[ok, 2], a[~ok, 0]


def designs7(n, rpm):
    interior = [r for r in rpm if r not in ENDS7]
    combos = list(itertools.combinations(interior, n - 2))
    if len(combos) > 35:
        rng = np.random.default_rng(0)
        combos = [combos[i] for i in sorted(rng.choice(len(combos), 35, replace=False))]
    return [list(ENDS7) + list(c) for c in combos]


# ----------------------------------------------------------------- model
def train_lf(x, ly):
    lf = CachedLooKRR(design_space=DS, params_optimize=True, noise_data=True,
                      optimizer_restart=10, seed=42)
    lf.train(col(x), np.asarray(ly, float), portion_test=LF_PORTION)
    return lf


def evaluate(lf, x, yl, yh, tr):
    msk = np.isin(x, tr)
    xh, yhh, xt, yt = x[msk], yh[msk], x[~msk], yh[~msk]
    out = dict(train=list(map(float, tr)), n_test=int(len(yt)), rel={}, hits={})
    models = {}
    for name, basis in (("MF", "linear"), ("HF-only GP", "ordinary")):
        m = FullBayesKRRLRGP(DS, lf_model=lf, basis=basis).train(col(xh), np.log(yhh))
        qq = np.exp(m.quantile(col(xt), Q3))
        out["rel"][name] = rel_rmse(qq[:, 1], yt)
        out["hits"][name] = int(np.sum((yt >= qq[:, 0]) & (yt <= qq[:, 2])))
        models[name] = m
    ratio = float(np.mean(yhh / np.interp(xh, x, yl)))
    out["rel"]["LF x ratio"] = rel_rmse(ratio * np.interp(xt, x, yl), yt)
    models["ratio"] = ratio
    return out, models


def summarise(rows):
    s = {}
    for m in METHODS:
        v = np.array([r["rel"][m] for r in rows])
        s[m] = dict(median=float(np.median(v)), q25=float(np.percentile(v, 25)),
                    q75=float(np.percentile(v, 75)))
    for m in ("MF", "HF-only GP"):
        s[m]["hits"] = int(sum(r["hits"][m] for r in rows))
        s[m]["coverage95"] = s[m]["hits"] / sum(r["n_test"] for r in rows)
    s["n_test_total"] = int(sum(r["n_test"] for r in rows))
    return s


def run_case(lf, x, yl, yh, n, keep=False):
    rows, mods, fails = [], {}, 0
    for tr in designs7(n, list(x)):
        try:
            o, m = evaluate(lf, x, yl, yh, tr)
        except Exception as e:  # noqa: BLE001
            fails += 1
            print("  FAIL", n, tr, e)
            continue
        rows.append(o)
        if keep:
            mods[tuple(tr)] = m
    return dict(summary=summarise(rows), rows=rows, n_fail=fails, n_designs=len(rows)), mods


# ----------------------------------------------------------------- main
def main():
    t0 = time.time()
    warnings.filterwarnings("ignore")
    lf_d, hf_d, ext = load_t7()
    data = {q: series(q, lf_d, hf_d, ext) for q in QOIS}

    # full raw table
    print("RAW TABLE (LF = L8 kmix, HF = L9 fig9; nan rpm dropped per QoI)")
    for q in QOIS:
        x, yl, yh, dropped = data[q]
        print(f"\n{q}  (dropped rpm: {list(dropped) if len(dropped) else 'none'})")
        print(f"  {'rpm':>6s} {'L8':>10s} {'L9':>10s} {'L9/L8':>7s}")
        for a, b, c in zip(x, yl, yh):
            print(f"  {a:6g} {b:10.4g} {c:10.4g} {c / b:7.2f}")
    print("\ndtmix_0.95 at 37.5: fig9_l9_rpm37.5 =", hf_d[37.5]["dtmix_0.95"], " ext1 =", ext["dtmix_0.95"])

    spear = {}
    for q in QOIS:
        x, yl, yh, _ = data[q]
        rho, p = spearmanr(yl, yh)
        spear[q] = dict(rho=float(rho), p=float(p), n=int(len(x)))

    res, mods_rep, rep = {}, {}, {}
    hits_t = n_t = 0
    for q in QOIS:
        x, yl, yh, _ = data[q]
        lf = train_lf(x, np.log(yl))
        for n in NS:
            r, mods = run_case(lf, x, yl, yh, n, keep=(n == 3))
            res[f"{q}|n{n}"] = r
            if n == 3:
                order = np.argsort([d["rel"]["MF"] for d in r["rows"]])
                d = r["rows"][order[len(order) // 2]]
                rep[q] = dict(train=d["train"], rel=d["rel"])
                mods_rep[q] = mods[tuple(d["train"])]
                hits_t += r["summary"]["MF"]["hits"]
                n_t += r["summary"]["n_test_total"]
    pooled = hits_t / n_t

    wins = {}
    lines = [f"{'QoI':12s} {'MF':>7s} {'HF-only':>8s} {'LFxratio':>9s} {'covMF':>6s} {'covHF':>6s} "
             f"{'#d':>3s} {'win':>4s} {'Spearman':>9s}"]
    for q in QOIS:
        s = res[f"{q}|n3"]["summary"]
        wins[q] = bool(s["MF"]["median"] < s["LF x ratio"]["median"] and s["MF"]["median"] < s["HF-only GP"]["median"])
        lines.append(f"{q:12s} {s['MF']['median']:7.1%} {s['HF-only GP']['median']:8.1%} "
                     f"{s['LF x ratio']['median']:9.1%} {s['MF']['coverage95']:6.0%} {s['HF-only GP']['coverage95']:6.0%} "
                     f"{res[f'{q}|n3']['n_designs']:3d} {'PASS' if wins[q] else 'FAIL':>4s} {spear[q]['rho']:9.3f}")
    n_pass = sum(wins.values())
    cov_ok = pooled >= 0.85
    passed = bool(n_pass >= 4 and cov_ok)
    lines.append(f"QoIs passing: {n_pass}/6 (need >= 4); pooled MF coverage {pooled:.1%} ({hits_t}/{n_t}), need >= 85%")
    lines.append("VERDICT: " + ("PASS" if passed else "FAIL"))
    table = "\n".join(lines)
    print("\n" + table)

    lc = [f"{'QoI':12s}{'n':>2s}{'#d':>4s}  " + " ".join(f"{m:>22s}" for m in METHODS) + "  cov MF/HF"]
    for q in QOIS:
        for n in NS:
            r = res[f"{q}|n{n}"]
            s = r["summary"]
            lc.append(f"{q:12s}{n:>2d}{r['n_designs']:4d}  " + " ".join(
                f"{s[m]['median']:7.1%}[{s[m]['q25']:5.1%},{s[m]['q75']:6.1%}]" for m in METHODS)
                + f"  {s['MF']['coverage95']:.0%}/{s['HF-only GP']['coverage95']:.0%}")
    lc_table = "\n".join(lc)
    print("\nLEARNING CURVE (reported only)\n" + lc_table)

    # ------------------------------------------------------- figure (a)
    plt.rcParams.update(fs.rcparams())
    f, axes = plt.subplots(2, 3, figsize=(11.0, 6.4), sharex=True)
    xg = np.linspace(15, 37.5, 120)
    cl, ch = fs.level_colour(8), fs.level_colour(9)
    for k, q in enumerate(QOIS):
        ax = axes[k // 3, k % 3]
        x, yl, yh, _ = data[q]
        tr = tuple(rep[q]["train"])
        mf, ratio = mods_rep[q]["MF"], mods_rep[q]["ratio"]
        qm = np.exp(mf.quantile(col(xg), Q3))
        trm = np.isin(x, tr)
        mk = MARK[q]
        ax.fill_between(xg, qm[:, 0], qm[:, 2], color=ch, alpha=0.18, lw=0, label="MF 95%")
        ax.plot(xg, qm[:, 1], color=ch, lw=1.4, ls="--", label="MF median")
        ax.plot(xg, ratio * np.interp(xg, x, yl), color="0.2", lw=1.2, ls=(0, (1, 3)), label="LF x ratio")
        ax.plot(x, yl, **fs.series_kw(cl, mk, stat="mean", ours=True, ls="none"), label="L8")
        ax.plot(x[trm], yh[trm], **fs.series_kw(ch, mk, stat="mean", ours=True, ls="none"), label="L9 training")
        ax.plot(x[~trm], yh[~trm], **fs.series_kw(ch, mk, stat="mean", ours=True, ls="none", ms=3.0, alpha=0.6),
                label="L9 held out")
        ax.set_yscale("log")
        ax.set_title(TITLE[q], fontsize=10)
        ax.grid(**fs.GRID_KW)
        if k % 3 == 0:
            ax.set_ylabel(YLAB["dtmix"] if k < 3 else YLAB["kLa"])
        if k >= 3:
            ax.set_xlabel("rocking speed (rpm)")
        if k == 2:
            ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    f.tight_layout()
    f.savefig(OUT / "test_t7_mix_kla_l8_l9.png", dpi=200, bbox_inches="tight")
    plt.close(f)

    # ------------------------------------------------------- figure (b)
    f, axes = plt.subplots(2, 3, figsize=(11.0, 6.2), sharex=True)
    sty = {"MF": dict(color=ch, ls="--"), "HF-only GP": dict(color=fs.NEUTRAL, ls="-."),
           "LF x ratio": dict(color="0.2", ls=(0, (1, 3)))}
    for k, q in enumerate(QOIS):
        ax = axes[k // 3, k % 3]
        for mi, m in enumerate(METHODS):
            xs = np.array(NS) + (mi - 1) * 0.08
            S = [res[f"{q}|n{n}"]["summary"][m] for n in NS]
            med = np.array([s["median"] * 100 for s in S])
            lo = np.array([s["q25"] * 100 for s in S])
            hi = np.array([s["q75"] * 100 for s in S])
            kw = fs.series_kw(sty[m]["color"], MARK[q], stat="mean", ours=True, ls=sty[m]["ls"], ms=4.5)
            ax.errorbar(xs, med, yerr=[med - lo, hi - med], capsize=2, elinewidth=0.7, label=m, **kw)
        ax.set_yscale("log")
        ax.set_title(TITLE[q], fontsize=10)
        ax.grid(**fs.GRID_KW)
        ax.set_xticks(NS)
        if k % 3 == 0:
            ax.set_ylabel("relative RMSE (%)")
        if k >= 3:
            ax.set_xlabel("number of L9 runs")
        if k == 2:
            ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    f.tight_layout()
    f.savefig(OUT / "test_t7_learning_curve.png", dpi=200, bbox_inches="tight")
    plt.close(f)

    rt = time.time() - t0
    slim = {k: dict(summary=v["summary"], n_designs=v["n_designs"], n_fail=v["n_fail"],
                    rows=v["rows"] if k.endswith("|n3") else None) for k, v in res.items()}
    json.dump(dict(pass_=passed, n_qoi_pass=n_pass, wins=wins, pooled_cov_mf=pooled, hits=hits_t, n_test=n_t,
                   spearman=spear, representative=rep,
                   raw={q: dict(rpm=data[q][0].tolist(), L8=data[q][1].tolist(), L9=data[q][2].tolist(),
                                dropped=data[q][3].tolist()) for q in QOIS},
                   results=slim, table=table, learning_curve_table=lc_table, runtime_s=rt),
              open(OUT / "test_t7_mix_kla_l8_l9.json", "w"), indent=1)
    print(f"runtime {rt:.1f} s")


if __name__ == "__main__":
    main()
