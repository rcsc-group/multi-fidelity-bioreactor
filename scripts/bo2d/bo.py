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


# ---------------------------------------------------------------- augmented EI (2026-11-07, Elvis chose it)
# Huang, Allen, Notz, Miller (2006) MFSKO, read in full in its restatement Schenk et al., WSC 2005, eqs. (7)-(9):
#   EI(x, l) = E[max(f_m(x) - f_m(x*), 0)] * a1 * a2 * a3        (written here for maximisation)
#   a1 = corr[f_l^p(x), f_m^p(x)]   posterior correlation of level l with the top level at x
#   a2 = 1 - s_eps,l / sqrt(s_l^2(x) + s_eps,l^2)                  noise discount at level l
#   a3 = C_m / C_l                                                cost ratio (C_l includes the warm-start prerequisites)
#   x* = argmax over evaluated points (any level) of u(x) = f_m(x) - c s_m(x), c = 1.
# Their model is additive (rho = 1). Ours is Kennedy-O'Hagan with a fitted rho per level (basis [1, f_{l-1}]); the
# variance is propagated up the chain assuming the level posteriors independent:
#   V8 = rho8^2 V7 + Vd8,  V9 = rho9^2 V8 + Vd9,  cov(f8, f9) = rho9 V8,  cov(f7, f9) = rho9 rho8 V7.
# Constraint: the score is multiplied by P(tau_95 <= cap) from the tau chain at L9.

def latent(m, X):
    """(mean, latent variance, noise variance) in ln units for one FullBayes model; noise = eta^2 term."""
    X = np.atleast_2d(X)
    mu, sd = gauss(m, X)
    noise = float(m.w @ (m.b_n / m.a_n * m.eta ** 2)) * m.yh_std ** 2
    return mu, np.maximum(sd ** 2 - noise, 1e-12), noise


def rho_of(m):
    return float(m.w @ m._rho[:, 1]) if m._rho.shape[1] > 1 else 0.0


class _ThetaBasis(FullBayesKRRLRGP):
    """Basis [1, f_lower, theta_n]: lets the L8->L9 discrepancy carry a linear theta trend (flat prior, integrated).
    Why (2026-11-07): the L8 kLa error grows as theta falls (real data ~4x at low theta); with basis [1, f8] the
    chain inherited L8's preference for theta ~ 2 on the dev synthetic set."""
    def _basis(self, X):
        F = super()._basis(X)
        return np.hstack([F, self.normalize_input(X)[:, 1:2]])


CHAIN_THETA = False   # set by callers; True = _ThetaBasis for the L8 and L9 layers


def fit_chain3(d7, d8, d9, qoi):
    """L7 FullBayes (constant basis) -> L8 (basis [1, mean f7]) -> L9 (basis [1, mean f8]). Returns (m7, m8, m9)."""
    floor = C.noise_floor(qoi)
    m7 = fit_single(d7, qoi)
    X8, y8 = _xy(d8)
    K = _ThetaBasis if CHAIN_THETA else FullBayesKRRLRGP
    m8 = K(C.DOMAIN, lf_model=_Mean(m7), grid=C.hf_grid(floor, float(np.std(y8)))).train(X8, y8)
    X9, y9 = _xy(d9)
    m9 = K(C.DOMAIN, lf_model=_Mean(m8), grid=C.hf_grid(floor, float(np.std(y9)))).train(X9, y9)
    return m7, m8, m9


def chain_moments(ms, X):
    """Per level l: mean, propagated latent variance V_l, noise var, and corr(f_l, f_9). ms = (m7, m8, m9)."""
    m7, m8, m9 = ms
    mu7, V7, n7 = latent(m7, X)
    mu8, Vd8, n8 = latent(m8, X)
    mu9, Vd9, n9 = latent(m9, X)
    r8, r9 = rho_of(m8), rho_of(m9)
    V8 = r8 ** 2 * V7 + Vd8
    V9 = r9 ** 2 * V8 + Vd9
    corr = {7: r9 * r8 * V7 / np.sqrt(V7 * V9), 8: r9 * V8 / np.sqrt(V8 * V9), 9: np.ones_like(V9)}
    return {7: (mu7, V7, n7), 8: (mu8, V8, n8), 9: (mu9, V9, n9)}, corr


def single_moments(m9, X):
    mu, V, n = latent(m9, X)
    return {9: (mu, V, n)}, {9: np.ones_like(V)}


def aei(mom_k, corr_k, mom_t, best, lncap, cost):
    """Score per level l: EI of ln kLa at L9 (latent) * a1 * a2 * a3 * P(feasible). cost[l] = core-h incl. prerequisites."""
    mu9, V9, _ = mom_k[9]
    s9 = np.sqrt(V9)
    z = (mu9 - best) / s9
    ei = (mu9 - best) * stats.norm.cdf(z) + s9 * stats.norm.pdf(z)
    mt, Vt, nt = mom_t[9]
    pof = stats.norm.cdf((lncap - mt) / np.sqrt(Vt + nt))
    out = {}
    for l, (mu, V, n) in mom_k.items():
        a2 = 1.0 - np.sqrt(n) / np.sqrt(V + n)
        out[l] = ei * np.abs(corr_k[l]) * a2 * (cost[9] / cost[l]) * pof
    return out


def effective_best(mom_k_at, mom_t_at, lncap):
    """x* rule: max over evaluated points of mu9 - s9 among points with P(feasible) >= 0.8. Returns ln value or None."""
    mu9, V9, _ = mom_k_at[9]
    mt, Vt, nt = mom_t_at[9]
    ok = mt + Z_REC * np.sqrt(Vt + nt) <= lncap
    u = mu9 - np.sqrt(V9)
    return float(u[ok].max()) if ok.any() else None


def pick_budget(scores, G, budget, cost, max_picks: int = 4):
    """Greedy (x, level) picks by score: at most max_picks, and within the core-h budget; each pick removes a disc
    around x. max_picks: aEI is a one-at-a-time rule; filling a budget with it floods the cheap level (77 L7 picks in
    one round of the synthetic test, 2026-11-07)."""
    Gs = scaled(G); S = {l: s.copy() for l, s in scores.items()}; out = []; spent = 0.0
    while True:
        l, i = max(((l, int(np.argmax(s))) for l, s in S.items()), key=lambda li: S[li[0]][li[1]])
        if len(out) >= max_picks or not np.isfinite(S[l][i]) or S[l][i] <= 0 or spent + cost[l] > budget:
            break
        out.append((G[i], l)); spent += cost[l]
        mask = np.linalg.norm(Gs - Gs[i], axis=1) < R_EXCL
        for s in S.values():
            s[mask] = -np.inf
    return out


def recommend_m(Mk, Mt, d9, lncap):
    """recommend() for moment functions: the evaluated L9 site with the highest mu9 among sites with P(feasible)>=0.8."""
    X = d9[["rpm", "theta"]].drop_duplicates().to_numpy(float)
    k, _ = Mk(X); t, _ = Mt(X)
    mu, (mt, Vt, nt) = k[9][0], t[9]
    ok = mt + Z_REC * np.sqrt(Vt + nt) <= lncap
    return (X[ok][np.argmax(mu[ok])], float(mu[ok].max())) if ok.any() else (None, None)


# --- MF-GP-UCB level rule (Kandasamy, Dasarathy, Oliva, Schneider, Poczos, NeurIPS 2016), with constrained EI ----
# Step 1: x from the top-level acquisition only (constrained EI of ln kLa at L9).
# Step 2: at x, run the lowest level l < 9 whose latent sd exceeds its threshold gamma_l; else run L9.
# gamma_l = noise sd of level l (the level is "learned" at x once its latent sd falls under its own noise).
# Why (2026-11-07): aEI's cost factor C9/C7 ~ 20 made L7 win nearly every round and pulled x to theta ~ 2, where
# L8 over-predicts kLa ~4x; the chain stalled at 9.5 % regret on the dev synthetic set.
def mf_pick(mom_k, mom_t, best, lncap, G, q: int = 4):
    mu9, V9, _ = mom_k[9]
    s9 = np.sqrt(V9)
    z = (mu9 - best) / s9
    ei = (mu9 - best) * stats.norm.cdf(z) + s9 * stats.norm.pdf(z)
    mt, Vt, nt = mom_t[9]
    score = ei * stats.norm.cdf((lncap - mt) / np.sqrt(Vt + nt))
    B = pick_batch(score, G, q)
    out = []
    for x in B:
        i = int(np.argmin(np.linalg.norm(G - x, axis=1)))
        lev = 9
        for l in sorted(k for k in mom_k if k < 9):
            _, V, n = mom_k[l]
            if np.sqrt(V[i]) > np.sqrt(np.broadcast_to(n, V.shape)[i]):
                lev = l
                break
        out.append((x, lev))
    return out
