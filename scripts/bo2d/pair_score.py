"""Score of the paired L8 -> L9 test, exactly as pre-registered in pair_l8.py (2026-11-10).

Sites: the 12 pair sites (L8 = pr8_l8_*, L9 = bo1_l9_* at the same x); secondary set: bo4 chain sites with L8 and L9.
Outer leave-one-site-out. L8 value known at the held-out site (nested design).
  (a) L9-only GP (basis ordinary)   (b) mfbml T = [1, f8]   (c) T = [1, f8, f8^2]   (d) T = [1, f8, theta]
all FullBayesKRRLRGP with the GP discrepancy Delta(x), ln units, mixing as the rate 1/t.
'best MF' = form chosen by INNER LOO on the 11 training sites (no outer peeking).
Rule: L8 carries L9 information off-line if RMSE(best MF) <= 0.8 x RMSE(a), for kLa_1T_10 (primary).
Also: tau_95 and dtmix_0.25; Spearman rho(f8, y9) descriptive only.
  uv run python scripts/bo2d/pair_score.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as LP  # noqa: E402
from scripts.bo2d import bo  # noqa: E402
from scripts.mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

C = LP.C
import os  # noqa: E402
if os.environ.get("PAIR_NOGATE") == "1":   # SENSITIVITY ONLY (2026-11-10 22:20): the 1 % wave-height gate drops 8/12
    C.warm_run_gate = lambda d: ("PASS", {})   # pair L8 runs whose tau/EDR drift is < 2 %; not the pre-registered score


class Lookup:
    def __init__(self, X, f):
        self.d = {tuple(np.round(x, 2)): v for x, v in zip(X, f)}

    def predict(self, X):
        return np.array([self.d[tuple(np.round(x, 2))] for x in np.atleast_2d(X)])


class Form(FullBayesKRRLRGP):
    extra = None   # "quad" | "theta" | None

    def _basis(self, X):
        B = super()._basis(X)
        if self.extra == "quad":
            return np.hstack([B, B[:, 1:2] ** 2])
        if self.extra == "theta":
            X = np.atleast_2d(X)
            return np.hstack([B, ((X[:, 1] - 4.5) / 2.5)[:, None]])
        return B


def fit_pred(form, Xtr, ytr, ftr_model, Xte, q):
    grid = C.hf_grid(C.noise_floor(q), float(np.std(ytr)))
    if form == "a":
        m = FullBayesKRRLRGP(C.DOMAIN, basis="ordinary", grid=grid)
    else:
        m = Form(C.DOMAIN, lf_model=ftr_model, basis="linear", grid=grid)
        m.extra = {"b": None, "c": "quad", "d": "theta"}[form]
    m.train(Xtr, ytr)
    return float(np.ravel(m.predict_mean(Xte))[0])


def loo(forms, X, y, lk, q):
    err = {k: [] for k in forms}
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        for k in forms:
            try:
                err[k].append(fit_pred(k, X[m], y[m], lk, X[i:i + 1], q) - y[i])
            except np.linalg.LinAlgError:
                err[k].append(np.nan)
    return {k: float(np.sqrt(np.nanmean(np.square(v)))) for k, v in err.items()}, err


def pairs(q, pre8, pre9, sites=None):
    d8 = LP.as_objective(LP.rows(q, [pre8], {8}), q)
    d9 = LP.as_objective(LP.rows(q, [pre9], {9}), q)
    out = []
    for _, r in d9.iterrows():
        y9 = float(r.y)
        if not (np.isfinite(y9) and y9 > 0):
            continue
        if sites is not None and not any(np.isclose([r.rpm, r.theta], s, atol=0.02).all() for s in sites):
            continue
        m = d8[np.isclose(d8.rpm, r.rpm, atol=0.02) & np.isclose(d8.theta, r.theta, atol=0.02)]
        m = m[np.isfinite(m.y.astype(float)) & (m.y.astype(float) > 0)]
        if len(m):
            out.append((r.rpm, r.theta, np.log(float(m.y.iloc[0])), np.log(y9)))
    return np.array(out, float).reshape(-1, 4)


def score(name, P, q):
    X, f, y = P[:, :2], P[:, 2], P[:, 3]
    lk = Lookup(X, f)
    outer, _ = loo("abcd", X, y, lk, q)
    # inner LOO form choice, per outer fold
    e_best = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        inner, _ = loo("bcd", X[m], y[m], lk, q)
        k = min(inner, key=inner.get)
        e_best.append(fit_pred(k, X[m], y[m], lk, X[i:i + 1], q) - y[i])
    best = float(np.sqrt(np.mean(np.square(e_best))))
    rho = spearmanr(f, y).correlation
    verdict = "CARRIES" if best <= 0.8 * outer["a"] else "does NOT carry"
    print(f"{name:10s} {q:11s} n={len(y):2d} L9-only {outer['a']:.3f} | lin {outer['b']:.3f} quad {outer['c']:.3f} "
          f"+theta {outer['d']:.3f} | best(inner) {best:.3f} ratio {best / outer['a']:.2f} -> {verdict}  "
          f"Spearman {rho:+.2f}  sd(y9) {np.std(y):.3f}", flush=True)


if __name__ == "__main__":
    from scripts.bo2d.pair_l8 import SITES
    for q in ("kLa_1T_10", "tau_95", "dtmix_0.25"):
        P = pairs(q, "pr8_", "bo1_", SITES)
        score("pair12", P, q)
        P2 = pairs(q, "bo4_", "bo4_")
        if len(P2) >= 5:
            score("bo4chain", P2, q)
