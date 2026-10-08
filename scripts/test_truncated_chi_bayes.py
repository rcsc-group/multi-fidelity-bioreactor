"""Test B (pre-registered, diary 2026-10-08): predict dtmix_0.95 from the chi(t) curve up to chi = 0.50, Bayesian.

The user's "whole curve" idea: a run is observed only until chi = 0.50 (about 100+ output samples); predict the
time to chi = 0.95. Target: the L9 rpm sweep (fig9_l9_rpm*), whose full curves are held out.

Model. After chi = 0.50 the variance decays as 1 - chi(t) = 0.5 exp(-r lam_e (t - t_50)), so
    t_95 = t_50 + ln(10) / (r lam_e).
  lam_e: early decay rate, least-squares slope of ln(1 - chi) on chi in [0.25, 0.50]. It is treated as known
         (assumption: its error is small next to the spread of r; the residuals are correlated within a cycle).
  r:     late/early rate ratio, the only unknown. For a run with a full curve, r is defined exactly by
         r = ln(10) / ((t_95 - t_50) lam_e), so the model reproduces t_95 when r is known.
Priors on log r:
  P1 physics      N(0, 0.6^2): the late rate is within a factor ~3 of the early rate.
  P2 population   N(m, s^2 (1 + 1/n)): m, s from the full L6-L8 curves, all rpm (cheap runs).
  P3 grid MF      N(log r_L8(rpm), s3^2), s3 = sd over rpm of log r_L8 - log r_L7.
Reference: the deterministic single exponential (r = 1), scripts/pilot_truncated_chi.py.
PASS per prior: 95% coverage >= 8/10 AND median |rel err| of the median < 20% AND median (hi - lo)/median <= 0.6.

Writes experiments/multifidelity/test_truncated_chi_bayes.{json,png} and ..._summary.png.
"""
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts.pilot_truncated_chi import chi_curve  # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
RUN = {6: "kmix_l6_rpm{:g}", 7: "kmix_l7_rpm{:g}", 8: "kmix_l8_rpm{:g}", 9: "fig9_l9_rpm{:g}"}
LN10 = math.log(10.0)
NDRAW = 20000


def features(run):
    """t_50, lam_e from the curve up to chi = 0.50; t_95 (nan if not reached); the observed part of the curve."""
    t, chi = chi_curve(run)
    k50 = int(np.argmax(chi >= 0.50))
    if chi[k50] < 0.50:
        return None
    m = (chi >= 0.25) & (np.arange(len(t)) <= k50)
    b, _ = np.polyfit(t[m], np.log(1 - chi[m]), 1)
    k95 = int(np.argmax(chi >= 0.95))
    t95 = t[k95] if chi[k95] >= 0.95 else math.nan
    return dict(t50=float(t[k50]), lam=float(-b), t95=float(t95), t=t, chi=chi, k50=k50, n_obs=int(k50 + 1))


def log_r(f):
    return math.log(LN10 / ((f["t95"] - f["t50"]) * f["lam"]))


def main():
    rng = np.random.default_rng(0)
    F = {lev: {r: features(RUN[lev].format(r)) for r in RPMS} for lev in RUN}
    lr = {lev: {r: log_r(f) for r, f in F[lev].items() if f and f["t95"] == f["t95"]} for lev in RUN}
    pop = np.array([v for lev in (6, 7, 8) for v in lr[lev].values()])
    m2, s2 = pop.mean(), pop.std(ddof=1) * math.sqrt(1 + 1 / len(pop))
    d87 = np.array([lr[8][r] - lr[7][r] for r in RPMS if r in lr[8] and r in lr[7]])
    s3 = float(np.std(d87, ddof=1))
    print(f"population log r (L6-L8, n={len(pop)}): mean {m2:+.2f}, sd {pop.std(ddof=1):.2f}; "
          f"L9 log r: {np.round([lr[9].get(r, np.nan) for r in RPMS], 2)}")
    print(f"P3 sd (log r_L8 - log r_L7 over rpm, n={len(d87)}): {s3:.2f}")

    priors = {"P1": lambda r: (0.0, 0.6), "P2": lambda r: (m2, s2),
              "P3": lambda r: (lr[8][r], s3) if r in lr[8] else (m2, s2)}
    res = {k: [] for k in list(priors) + ["exp1"]}
    for r in RPMS:
        f = F[9][r]
        if f is None or f["t95"] != f["t95"]:
            print(f"{r:5g} rpm: no measured t_95 at L9, skipped")
            continue
        for name, pr in priors.items():
            mu, sd = pr(r)
            rr = np.exp(mu + sd * rng.standard_normal(NDRAW))
            t95 = f["t50"] + LN10 / (rr * f["lam"])
            lo, med, hi = np.quantile(t95, [0.025, 0.5, 0.975])
            res[name].append(dict(rpm=r, obs=f["t95"], lo=lo, med=med, hi=hi, mu=mu, sd=sd))
        p = f["t50"] + LN10 / f["lam"]
        res["exp1"].append(dict(rpm=r, obs=f["t95"], lo=p, med=p, hi=p))
    summ = {}
    for name, rows in res.items():
        obs = np.array([x["obs"] for x in rows])
        med = np.array([x["med"] for x in rows])
        lo, hi = np.array([x["lo"] for x in rows]), np.array([x["hi"] for x in rows])
        cov = int(np.sum((obs >= lo) & (obs <= hi)))
        err = float(np.median(np.abs(med / obs - 1)))
        wid = float(np.median((hi - lo) / med))
        ok = cov >= 0.8 * len(rows) and err < 0.20 and wid <= 0.6
        summ[name] = dict(n=len(rows), coverage=cov, median_abs_rel_err=err, median_rel_width=wid,
                          passed=bool(ok) if name != "exp1" else None)
        print(f"{name:5s} coverage {cov}/{len(rows)}  median|err| {err:5.1%}  median width {wid:5.2f}  "
              f"{'' if name == 'exp1' else ('PASS' if ok else 'FAIL')}")
    out = dict(summary=summ, rows=res, log_r={str(k): v for k, v in lr.items()},
               n_obs_points={str(r): F[9][r]["n_obs"] for r in RPMS if F[9][r]})
    json.dump(out, open(ROOT / "experiments/multifidelity/test_truncated_chi_bayes.json", "w"), indent=1,
              default=float)
    plot_curves(F[9], priors, rng)
    plot_summary(res)


def plot_curves(F9, priors, rng):
    """Per rpm: 1 - chi on a log axis (an exponential tail is a straight line). Observed part bold, hidden part
    faint; fans = 95% band of the prediction under P1 (neutral) and P3 (L9 hue)."""
    plt.rcParams.update(fs.rcparams())
    c9 = fs.level_colour(9)
    fig, axes = plt.subplots(2, 5, figsize=(15, 5.6), sharey=True)
    for ax, r in zip(axes.ravel(), RPMS):
        f = F9[r]
        if f is None:
            ax.set_visible(False)
            continue
        t, chi, k = f["t"], f["chi"], f["k50"]
        ax.plot(t[k:], 1 - chi[k:], color=c9, ls=":", lw=0.8, alpha=0.45, label="L9, hidden")
        ax.plot(t[:k + 1], 1 - chi[:k + 1], color=c9, ls=":", lw=1.8, label=r"L9, observed ($\chi \leq 0.5$)")
        tmax = (f["t95"] if f["t95"] == f["t95"] else t[-1]) * 1.6
        tg = np.linspace(f["t50"], tmax, 200)
        for name, colour, alpha in (("P1", fs.NEUTRAL, 0.18), ("P3", c9, 0.28)):
            mu, sd = priors[name](r)
            rr = np.exp(mu + sd * rng.standard_normal(4000))
            curves = 0.5 * np.exp(-np.outer(rr * f["lam"], tg - f["t50"]))
            lo, hi = np.quantile(curves, [0.025, 0.975], axis=0)
            ax.fill_between(tg, lo, hi, color=colour, alpha=alpha, lw=0, label=f"prediction, prior {name}")
        if f["t95"] == f["t95"]:
            ax.plot([f["t95"]], [0.05], **fs.series_kw(c9, fs.MK["chi_0.95"], stat="max", ours=True),
                    label=r"measured $\Delta t_{0.95}$")
        ax.axhline(0.05, color=fs.GUIDE, lw=0.8)
        ax.set_yscale("log")
        ax.set_ylim(0.01, 1.05)
        ax.set_xlim(0, tmax)
        ax.set_title(f"{r:g} rpm", fontsize=10)
        ax.grid(**fs.GRID_KW)
    for ax in axes[1]:
        ax.set_xlabel("time after release (s)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$1 - \chi$")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper left", bbox_to_anchor=(1.0, 0.9), frameon=False)
    fig.tight_layout()
    fig.savefig(ROOT / "experiments/multifidelity/test_truncated_chi_bayes.png", dpi=150, bbox_inches="tight")


def plot_summary(res):
    """Predicted against measured dtmix_0.95, 95% bars, one panel per prior."""
    plt.rcParams.update(fs.rcparams())
    c9 = fs.level_colour(9)
    names = ["exp1", "P1", "P2", "P3"]
    titles = {"exp1": "single exponential", "P1": "prior P1 (physics)", "P2": "prior P2 (L6-L8 population)",
              "P3": "prior P3 (L8 at same rpm)"}
    fig, axes = plt.subplots(1, 4, figsize=(15, 4), sharex=True, sharey=True)
    for ax, n in zip(axes, names):
        rows = res[n]
        obs = np.array([x["obs"] for x in rows])
        med = np.array([x["med"] for x in rows])
        err = np.array([med - [x["lo"] for x in rows], [x["hi"] for x in rows] - med])
        ax.errorbar(obs, med, yerr=err, ls="none", ecolor=c9, elinewidth=1,
                    **{k: v for k, v in fs.series_kw(c9, fs.MK["chi_0.95"], stat="max", ours=True).items()
                       if k != "ls"})
        g = np.array([50, 2000])
        ax.plot(g, g, color=fs.GUIDE, lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*g)
        ax.set_ylim(*g)
        ax.set_title(titles[n], fontsize=10)
        ax.set_xlabel(r"measured $\Delta t_{0.95}$ (s)")
        ax.grid(**fs.GRID_KW)
    axes[0].set_ylabel(r"predicted $\Delta t_{0.95}$ (s)")
    fig.tight_layout()
    fig.savefig(ROOT / "experiments/multifidelity/test_truncated_chi_bayes_summary.png", dpi=150,
                bbox_inches="tight")


if __name__ == "__main__":
    main()
