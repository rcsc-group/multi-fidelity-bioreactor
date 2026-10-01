"""Numerical uncertainty from a grid-refinement study, Eca & Hoekstra (2014).

L. Eca, M. Hoekstra, "A procedure for the estimation of the numerical
uncertainty of CFD calculations based on grid refinement studies",
J. Comput. Phys. 262 (2014) 104-130. Implemented from Appendix A (procedure)
and Appendix B (least-squares solutions), read from the publisher PDF on
marin.nl on 2026-10-01. Equation numbers below are the paper's.

Error estimators (phi_i - phi_o as a function of cell size h_i):
    RE  alpha h^p                  Eq. (1), p free
    1   alpha h                    Eq. (5)
    2   alpha h^2                  Eq. (6)
    12  alpha1 h + alpha2 h^2      Eq. (7)
Each is fitted by least squares, unweighted (w_i = 1) and weighted
(w_i = (1/h_i)/sum_j(1/h_j), Eq. 16). Standard deviation of a fit (App. B):
    sigma = sqrt( sum nw_i (phi_i - fit_i)^2 / (n_g - n_par) ),
    nw_i = 1 (unweighted) or n_g w_i (weighted).

Procedure (App. A):
 1. Fit RE with and without weights; discard fits with p <= 0.
    - if a fit has 0.5 <= p <= 2: use RE (the smallest sigma among those);
    - elif the observed p > 2: best of {1, 2} x {unweighted, weighted};
    - else (p < 0.5, or no p > 0: "anomalous"): best of {1, 2, 12} x {u, w}.
 2. Data range Delta_phi = (max phi - min phi) / (n_g - 1)      Eq. (19)
 3. F_s = 1.25 if 0.5 <= p < 2.1 and sigma < Delta_phi, else 3.
 4. eps = |fit(h_i) - phi_o|, the error estimate at grid i.
    sigma <  Delta_phi:  U = F_s eps + sigma + |phi_i - fit(h_i)|          Eq. (20)
    sigma >= Delta_phi:  U = 3 (sigma/Delta_phi)(eps + sigma + |phi_i - fit(h_i)|)  Eq. (21)
U is the 95% half-width: phi_i - U <= phi_exact <= phi_i + U  (Eq. 18).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar

P_SEARCH = (-8.0, 8.0)   # bracket for the observed order; p = 0 is excluded


def weights(h: np.ndarray) -> np.ndarray:
    """Eq. (16): w_i = (1/h_i) / sum_j (1/h_j)."""
    inv = 1.0 / np.asarray(h, float)
    return inv / inv.sum()


def _linear_fit(h, phi, w, exps):
    """Weighted LS for phi = phi_o + sum_k a_k h^exps[k]; returns coefficients."""
    A = np.column_stack([np.ones_like(h)] + [h ** e for e in exps])
    sw = np.sqrt(w)
    c, *_ = np.linalg.lstsq(A * sw[:, None], phi * sw, rcond=None)
    return c


def _make(c, exps):
    c = np.array(c, float)
    exps = tuple(float(e) for e in exps)
    return lambda x: c[0] + sum(ck * np.asarray(x, float) ** e for ck, e in zip(c[1:], exps))


def _sigma(h, phi, nw, f, n_par):
    dof = len(h) - n_par
    return float(np.sqrt(np.sum(nw * (phi - f(h)) ** 2) / dof)) if dof > 0 else 0.0


def _fit_re(h, phi, w, nw):
    """Eq. (1) with p free: for fixed p the problem is linear; minimise over p."""
    def s_of(p):
        c = _linear_fit(h, phi, w, (p,))
        return float(np.sum(w * (phi - _make(c, (p,))(h)) ** 2))
    best = None
    for lo, hi in ((P_SEARCH[0], -1e-3), (1e-3, P_SEARCH[1])):
        grid = np.linspace(lo, hi, 400)
        p0 = grid[int(np.argmin([s_of(p) for p in grid]))]
        step = (hi - lo) / 399
        r = minimize_scalar(s_of, bounds=(max(lo, p0 - step), min(hi, p0 + step)),
                            method="bounded")
        if best is None or r.fun < best.fun:
            best = r
    p = float(best.x)
    c = _linear_fit(h, phi, w, (p,))
    f = _make(c, (p,))
    return dict(estimator="RE", p=p, phi0=float(c[0]), fit=f,
                sigma=_sigma(h, phi, nw, f, 3))


def _fit_fixed(h, phi, w, nw, name):
    exps = {"1": (1.0,), "2": (2.0,), "12": (1.0, 2.0)}[name]
    c = _linear_fit(h, phi, w, exps)
    f = _make(c, exps)
    return dict(estimator=name, p=None, phi0=float(c[0]), fit=f,
                sigma=_sigma(h, phi, nw, f, 1 + len(exps)))


def numerical_uncertainty(h, phi, index: int = 0) -> dict:
    """U for phi[index] (default: the first entry, which should be the finest grid)."""
    h, phi = np.asarray(h, float), np.asarray(phi, float)
    ng = len(h)
    if ng < 4:
        raise ValueError(f"Eca & Hoekstra need at least 4 grids, got {ng}")
    variants = [("u", np.ones(ng), np.ones(ng)),
                ("w", weights(h), ng * weights(h))]

    re = []
    for tag, w, nw in variants:
        r = _fit_re(h, phi, w, nw)
        r["weighted"] = tag == "w"
        if r["p"] > 0:
            re.append(r)
    in_range = [r for r in re if 0.5 <= r["p"] <= 2.0]
    p_obs = min(re, key=lambda r: r["sigma"])["p"] if re else None

    if in_range:
        cands = in_range
    else:
        names = ("1", "2") if (p_obs is not None and p_obs > 2.0) else ("1", "2", "12")
        cands = []
        for name in names:
            for tag, w, nw in variants:
                r = _fit_fixed(h, phi, w, nw, name)
                r["weighted"] = tag == "w"
                cands.append(r)
    best = min(cands, key=lambda r: r["sigma"])
    p = best["p"] if best["estimator"] == "RE" else p_obs

    data_range = (phi.max() - phi.min()) / (ng - 1)
    sigma, f = best["sigma"], best["fit"]
    fs = 1.25 if (p is not None and 0.5 <= p < 2.1 and sigma < data_range) else 3.0
    eps = abs(float(f(h[index])) - best["phi0"])
    resid = abs(phi[index] - float(f(h[index])))
    if sigma < data_range:
        U = fs * eps + sigma + resid
    else:
        U = 3.0 * sigma / data_range * (eps + sigma + resid)
    return dict(U=float(U), Fs=fs, p=p, phi0=best["phi0"], sigma=sigma, eps=eps,
                resid=float(resid), data_range=float(data_range), fit=f,
                estimator=best["estimator"], weighted=best["weighted"],
                candidates=[(c["estimator"], c["weighted"], c["sigma"]) for c in cands])
