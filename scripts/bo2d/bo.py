"""Constrained BO for the DFD comparison (Elvis, 2026-11-07).

Problem: max kLa_1T_10 subject to tau_95 <= CAP, on (rpm, theta) in the campaign domain. Both in ln units.
Arms:
  l9only : one fully Bayesian GP per QoI on L9 rows (FullBayesKRRLRGP, constant basis).
  chain  : L7 KRR (LF, as the campaign) -> L8 FullBayes with basis [1, f7] -> L9 FullBayes with basis [1, mean f8].
           The L8 model variance is not propagated into the L9 basis (BACKLOG H3): the chain's band is too narrow
           between design points. Stated, not fixed.
Acquisition: cEI(x) = EI of ln kLa above the best feasible L9 observation, times P(ln tau_95 <= ln CAP).
  Predictive of each model = Gaussian with the mixture mean and variance (as campaign.crit).
  Batch of q: greedy on a dense grid, each pick removes a disc of radius R_EXCL (scaled units) around it.
Chain cascade: a proposed x runs the lowest level not yet run at x (L7, then L8, then L9). A point reaches L9 only
  if it is still a top proposal after its L7 and L8 rows are in. That screening is the saving we test.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts")); sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))
import campaign as C  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

CAP = None              # tau_95 cap [nondim tau_95 units of the database]; set by choose_cap() and frozen in state
R_EXCL = 0.08
SAME_X = 0.02           # scaled distance under which two rows are "the same point"
MAP_N = 101
Z_REC = 0.84          # recommendation needs P(feasible) >= 0.8


class _Mean:
    """Wrap a FullBayes model so it can be the LF of the next level (its predictive mean, ln units)."""
    def __init__(self, m):
        self.m = m

    def predict(self, X):
        return self.m.predict_mean(np.atleast_2d(X))


def gauss(m, X, chunk: int = 300):
    """Mixture mean and sd (ln units), chunked over X: predict_nodes holds (G, n, m) arrays."""
    X = np.atleast_2d(X); mo = np.empty(len(X)); so = np.empty(len(X))
    for i in range(0, len(X), chunk):
        mu, sc, nu, w = m.predict_dist(X[i:i + chunk])
        wm = w @ mu
        v = (w[:, None] * sc ** 2 * (nu / (nu - 2.0))).sum(0) + (w[:, None] * (mu - wm) ** 2).sum(0)
        mo[i:i + chunk], so[i:i + chunk] = wm, np.sqrt(v)
    return mo, so


def recommend(mk, mt, d9, lncap):
    """Recommended point = the evaluated L9 site with the highest model-mean ln kLa among sites whose model-mean
    ln tau_95 is <= ln cap. Not the best noisy observation: that one can be infeasible in truth."""
    X = d9[["rpm", "theta"]].drop_duplicates().to_numpy(float)
    k, _ = gauss(mk, X); t, st = gauss(mt, X)
    ok = t + Z_REC * st <= lncap          # feasible with probability >= 80 % under the model
    return (X[ok][np.argmax(k[ok])], float(k[ok].max())) if ok.any() else (None, None)


def _xy(d):
    return d[["rpm", "theta"]].to_numpy(float), np.log(d.y.to_numpy(float))


def fit_single(d9, qoi):
    X, y = _xy(d9)
    m = FullBayesKRRLRGP(C.DOMAIN, basis="ordinary", grid=C.hf_grid(C.noise_floor(qoi), float(np.std(y))))
    return m.train(X, y)


def fit_chain(d7, d8, d9, qoi):
    floor = C.noise_floor(qoi)
    X7, y7 = _xy(d7)
    lf7 = C.fit_lf_selected(X7, y7.reshape(-1, 1), floor)[0]
    X8, y8 = _xy(d8)
    m8 = FullBayesKRRLRGP(C.DOMAIN, lf_model=lf7, grid=C.hf_grid(floor, float(np.std(y8)))).train(X8, y8)
    X9, y9 = _xy(d9)
    m9 = FullBayesKRRLRGP(C.DOMAIN, lf_model=_Mean(m8), grid=C.hf_grid(floor, float(np.std(y9)))).train(X9, y9)
    return m9, m8


def grid():
    r = np.linspace(*C.RPM, MAP_N); t = np.linspace(*C.THETA, MAP_N)
    RR, TT = np.meshgrid(r, t, indexing="ij")
    return np.column_stack([RR.ravel(), TT.ravel()])


def scaled(X):
    return (np.atleast_2d(X) - C.DOMAIN[:, 0]) / (C.DOMAIN[:, 1] - C.DOMAIN[:, 0])


def cei(mk, sk, mt, st, best, lncap):
    pof = stats.norm.cdf((lncap - mt) / np.maximum(st, 1e-12))
    if best is None:
        return pof
    z = (mk - best) / np.maximum(sk, 1e-12)
    ei = (mk - best) * stats.norm.cdf(z) + sk * stats.norm.pdf(z)
    return ei * pof


def best_feasible(dk9, dt9, cap):
    """Best observed ln kLa among L9 rows whose observed tau_95 <= cap (same run)."""
    m = dk9.merge(dt9[["run", "y"]], on="run", suffixes=("_k", "_t"))
    f = m[m.y_t <= cap]
    return (float(np.log(f.y_k.max())), f.loc[f.y_k.idxmax(), "run"]) if len(f) else (None, None)


def pick_batch(score, G, q):
    Gs = scaled(G); s = score.copy(); out = []
    for _ in range(q):
        i = int(np.argmax(s))
        if not np.isfinite(s[i]) or s[i] <= 0:
            break
        out.append(G[i]); s[np.linalg.norm(Gs - Gs[i], axis=1) < R_EXCL] = -np.inf
    return np.array(out)


def next_level(x, rows_by_level):
    """Lowest level (7, 8, 9) with no row within SAME_X of x."""
    xs = scaled(x)[0]
    for lev in (7, 8, 9):
        d = rows_by_level.get(lev)
        if d is None or len(d) == 0 or np.min(np.linalg.norm(scaled(d[["rpm", "theta"]].to_numpy(float)) - xs, axis=1)) > SAME_X:
            return lev
    return None
