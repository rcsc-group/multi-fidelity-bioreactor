"""T4 (pre-registered, diary 2026-10-08): retrospective pool-based constrained BO on the 9 L10 rpm.

Pre-registration (diary.md, quoted):
  "T4 (zero compute): retrospective pool-based BO on the 10-rpm grid. Objective: maximise EDR at L10
   subject to tau_95(L10) <= the median L10 tau_95 (a stated proxy for 'mixing vs shear'; kLa at L10 has
   3 points only). MF-BO (L8 known at all rpm, L10 acquired one at a time from 2 initial endpoints) vs
   single-fidelity BO (L10 only, same start). Metric: L10 evaluations to find the true best feasible
   rpm, over all start pairs/seeds. Descriptive."
Task spec: pool = the 9 L10 rpm; truth = observed L10 values (hydro_dataset_v2.json, edr_mean_t,
tau_95_t); limit = median of the 9 L10 tau_95_t (5th value; feasible = tau_95 <= limit, so 5 rpm are
feasible); true best feasible = argmax EDR among feasible. Start: all 36 pairs of L10 rpm. Each step
acquires the pool rpm maximising constrained EI = EI_EDR * P(tau_95 <= limit) (incumbent = best observed
feasible EDR; if none observed feasible, P(feasible) alone), both from the surrogate predictive
distribution (Student-t for MF-Bayes, normal otherwise; total predictive incl. noise).
Arms: (i) MF-Bayes LF L8, (ii) HF-only GP (L10 only; reference class with constant basis, normal),
(iii) MF-Bayes LF L9, (iv) random acquisition (mean over 50 seeds per start pair).
Metric: number of L10 evaluations (starts included) until the true best feasible rpm has been evaluated.

Choices where the spec was silent:
  * One surrogate per QoI (EDR, tau_95), each refit on all observed L10 points at each step; raw units.
  * The MF-Bayes models are the reference-trained model viewed as the Bayes subclass (same trained
    state, T2) with LF training cached (T3 helpers).
  * Ties / all-zero acquisition (EI underflow): lowest rpm among the maximisers (np.argmax). The number of
    such steps is reported per arm. A failed fit (Cholesky) falls back to a uniformly random pool pick
    (seeded by start pair); count reported.
  * Student-t EI closed form (derived here): for Y ~ t_nu(m, s^2), z = (m-b)/s,
    E[(Y-b)+] = s [ z T_nu(z) + (nu + z^2)/(nu - 1) t_nu(z) ].
  * If a start pair already contains the true best feasible rpm the count is 2.
  * Random arm: per start pair the mean count over 50 seeds (np.random.default_rng(seed)).
Figure: histogram-free summary -- per arm, the empirical distribution (jittered dots) of the count over the
36 start pairs with the median marked.
Usage: uv run python scripts/test_t4_retro_bo.py
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
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import figstyle as fs                                          # noqa: E402
from scripts.test_mf_l10_holdout import load, col, OUT                      # noqa: E402
from scripts.test_t3_learning_curve import fit_cached, bayes_view           # noqa: E402

ARMS = ("MF-Bayes L8", "HF-only GP", "MF-Bayes L9", "random")


def ei_normal(m, s, b):
    s = np.maximum(s, 1e-300)
    z = (m - b) / s
    return s * (z * stats.norm.cdf(z) + stats.norm.pdf(z))


def ei_t(m, s, nu, b):
    s = np.maximum(s, 1e-300)
    z = (m - b) / s
    return s * (z * stats.t.cdf(z, nu) + (nu + z ** 2) / (nu - 1) * stats.t.pdf(z, nu))


def predictive(arm, D, rpm, obs_idx, qty):
    """Return (m, s, nu or None) at all pool rpm for QoI qty, trained on L10 points obs_idx."""
    xh = rpm[obs_idx]
    yh = np.array([D[10][r][qty] for r in xh])
    if arm == "HF-only GP":
        xl = np.array(sorted(D[8])); yl = np.array([D[8][r][qty] for r in xl])   # unused by constant basis
        m = fit_cached("ordinary", xh, yh, xl, yl)
        mu, sd = m.predict(col(rpm), return_std=True)
        return mu.ravel(), sd.ravel(), None
    lfl = 8 if arm.endswith("L8") else 9
    xl = np.array(sorted(D[lfl])); yl = np.array([D[lfl][r][qty] for r in xl])
    m = fit_cached("linear", xh, yh, xl, yl)
    mu, scale, nu = bayes_view(m).predict_dist(col(rpm))
    return mu, scale, nu


def run_arm(arm, D, rpm, edr, tau, limit, best_i, start, rng_seed):
    obs = list(start)
    stats_ = dict(fallback=0, zero_acq=0)
    while best_i not in obs:
        pool = [i for i in range(len(rpm)) if i not in obs]
        feas_obs = [i for i in obs if tau[i] <= limit]
        try:
            me, se, ne = predictive(arm, D, rpm, obs, "edr_mean_t")
            mt, st_, nt = predictive(arm, D, rpm, obs, "tau_95_t")
            pf = (stats.t.cdf((limit - mt) / np.maximum(st_, 1e-300), nt) if nt is not None
                  else stats.norm.cdf((limit - mt) / np.maximum(st_, 1e-300)))
            if feas_obs:
                b = max(edr[i] for i in feas_obs)
                ei = ei_t(me, se, ne, b) if ne is not None else ei_normal(me, se, b)
                acq = ei * pf
            else:
                acq = pf
            acq = np.nan_to_num(acq[pool], nan=0.0)
            if acq.max() <= 1e-300:
                stats_["zero_acq"] += 1
            nxt = pool[int(np.argmax(acq))]
        except Exception:  # noqa: BLE001
            stats_["fallback"] += 1
            nxt = pool[int(np.random.default_rng(rng_seed + len(obs)).integers(len(pool)))]
        obs.append(nxt)
    return len(obs), stats_


def run_random(rpm_n, best_i, start, n_seeds=50):
    out = []
    for s in range(n_seeds):
        rng = np.random.default_rng(s)
        obs = list(start)
        while best_i not in obs:
            pool = [i for i in range(rpm_n) if i not in obs]
            obs.append(pool[int(rng.integers(len(pool)))])
        out.append(len(obs))
    return float(np.mean(out))


def main():
    t0 = time.time()
    warnings.filterwarnings("ignore")
    D = load()
    rpm = np.array(sorted(D[10]))
    edr = np.array([D[10][r]["edr_mean_t"] for r in rpm])
    tau = np.array([D[10][r]["tau_95_t"] for r in rpm])
    limit = float(np.median(tau))
    feas = tau <= limit
    best_i = int(np.argmax(np.where(feas, edr, -np.inf)))
    print(f"limit tau_95 = {limit:.5f}; feasible rpm {rpm[feas].tolist()}; true best feasible rpm {rpm[best_i]}")
    pairs = list(itertools.combinations(range(len(rpm)), 2))
    counts = {a: [] for a in ARMS}
    diag = {a: dict(fallback=0, zero_acq=0) for a in ARMS}
    for pi, pr in enumerate(pairs):
        for a in ARMS[:3]:
            c, d = run_arm(a, D, rpm, edr, tau, limit, best_i, pr, rng_seed=pi * 1000)
            counts[a].append(c)
            for k in d:
                diag[a][k] += d[k]
        counts["random"].append(run_random(len(rpm), best_i, pr))
    summ = {}
    for a in ARMS:
        v = np.array(counts[a], float)
        summ[a] = dict(median=float(np.median(v)), mean=float(v.mean()), q25=float(np.percentile(v, 25)),
                       q75=float(np.percentile(v, 75)), min=float(v.min()), max=float(v.max()),
                       hist={str(k): int(np.sum(v == k)) for k in range(2, 10)} if a != "random" else None,
                       **diag[a])
    lines = [f"{'arm':14s} median  mean  IQR        min max  zero_acq fallback"]
    for a in ARMS:
        s = summ[a]
        lines.append(f"{a:14s} {s['median']:6.2f} {s['mean']:5.2f} [{s['q25']:.2f},{s['q75']:.2f}]  {s['min']:.0f}  {s['max']:.0f}"
                     f"   {s['zero_acq']:5d} {s['fallback']:5d}")
    table = "\n".join(lines)
    print(table)
    for a in ARMS[:3]:
        print(a, "count histogram (2..9):", summ[a]["hist"])

    plt.rcParams.update(fs.rcparams())
    f, ax = plt.subplots(figsize=(5.6, 3.6))
    cols = {"MF-Bayes L8": fs.level_colour(10), "MF-Bayes L9": fs.level_colour(10),
            "HF-only GP": fs.NEUTRAL, "random": "0.7"}
    mks = {"MF-Bayes L8": "s", "MF-Bayes L9": "s", "HF-only GP": "s", "random": "s"}
    rng = np.random.default_rng(1)
    for xi, a in enumerate(ARMS):
        v = np.array(counts[a], float)
        mk = fs.series_kw(cols[a], fs.MK["ediss"], stat="mean", ours=True, ls="none", ms=3.5, alpha=0.5)
        ax.plot(xi + rng.uniform(-0.18, 0.18, len(v)), v + rng.uniform(-0.12, 0.12, len(v)), **mk)
        ax.hlines(np.median(v), xi - 0.3, xi + 0.3, color="k", lw=1.6)
    ax.set_xticks(range(len(ARMS)))
    ax.set_xticklabels(["MF-Bayes\nLF L8", "HF-only\nGP", "MF-Bayes\nLF L9", "random"], fontsize=8)
    ax.set_ylabel("L10 evaluations to best feasible rpm")
    ax.set_yticks(range(2, 10))
    ax.grid(axis="y", **fs.GRID_KW)
    f.tight_layout()
    f.savefig(OUT / "test_t4_retro_bo.png", dpi=200, bbox_inches="tight")
    rt = time.time() - t0
    json.dump(dict(limit_tau95=limit, feasible_rpm=rpm[feas].tolist(), best_feasible_rpm=float(rpm[best_i]),
                   start_pairs=[[float(rpm[i]), float(rpm[j])] for i, j in pairs], counts=counts,
                   summary=summ, table=table, runtime_s=rt), open(OUT / "test_t4_retro_bo.json", "w"), indent=1)
    print(f"runtime {rt:.1f} s")


if __name__ == "__main__":
    main()
