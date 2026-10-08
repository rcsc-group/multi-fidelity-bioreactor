"""Fully Bayesian KRR-LR-GPR residual GP (T2b, diary 2026-10-08).

Pre-registration quoted from diary.md (2026-10-08, "T2b PRE-REGISTERED"):
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

Kernel parameterisation (read from scripts/mfbml_local/kernels.py)
  RBF stores param = 10**theta and evaluates k(x, y) = exp(-param * |x - y|^2) on inputs scaled to
  [0, 1] by the design space (normalize_input: (x - lo) / (hi - lo)).
  The conventional length scale l is defined by k = exp(-d^2 / (2 l^2)), hence
        param = 1 / (2 l^2),   theta = log10(param) = -log10(2 l^2).
  The prior "log l ~ N(log 0.3, 0.75^2) in scaled input units" is therefore placed on
  l = sqrt(1 / (2 param)); the grid is uniform in log l over log(0.3) +- z * 0.75 with
  z = Phi^-1(0.9995) = 3.29 (central 99.9%), i.e. l in [0.0255, 3.53]. Equal spacing in log l means
  the grid weight of a node is the prior density N(log l; log 0.3, 0.75^2) (the spacing is constant).
  The noise ratio eta = noise sd / HF std (HF std = np.std of the training outputs, as in the
  reference class) is log-spaced on [1e-3, 0.3]; log-uniform prior = equal prior weight per node.

Model (standardised units; y = (HF - mean) / std)
  y = F rho + Z,  Z ~ GP(0, s2 K),  K = k_l(X, X) + (eta^2 + 1e-10) I     (1e-10 = reference jitter)
  F = [1, f_l(X)] (MF, p = 2) or [1] (HF-only, p = 1); f_l = (LF KRR prediction - mean) / std.
  rho flat, s2 ~ IG(a0, b0), (l, eta) ~ prior above, on a grid.

Derivation of the marginal likelihood p(y | l, eta) (rho and s2 integrated out), per node
  A = F'K^-1F, rho_hat = A^-1 F'K^-1 y, r = y - F rho_hat, S = r'K^-1 r.
  (y - F rho)'K^-1(y - F rho) = S + (rho - rho_hat)' A (rho - rho_hat). Integrating the flat prior in
  rho (p dims) gives (2 pi s2)^(p/2) |A|^(-1/2), so
     p(y | s2, l, eta) = (2 pi)^(-(n-p)/2) |K|^(-1/2) |A|^(-1/2) s2^(-(n-p)/2) exp(-S / (2 s2)).
  With p(s2) = b0^a0 / Gamma(a0) s2^(-a0-1) exp(-b0 / s2):
     p(y | l, eta) = c * |K|^(-1/2) |A|^(-1/2) * Gamma(a_n) / (b0 + S/2)^(a_n),   a_n = a0 + (n-p)/2,
  where c = (2 pi)^(-(n-p)/2) b0^a0 / Gamma(a0) and Gamma(a_n) do not depend on (l, eta).  Hence
     log p(y | l, eta) = -1/2 log|K| - 1/2 log|A| - a_n log(b0 + S/2) + const.
  Posterior grid weight  w_g = prior_g * p(y | l_g, eta_g) / sum_g'(...), normalised by logsumexp.

Predictive at x* (per node, exactly the T2 algebra, eqs. (2)-(4) of krr_lr_gpr_bayes.py)
  y* | y, l, eta ~ t_nu(m*, (b_n / a_n) v*),  nu = n - p + 2 a0,
  m* = f*' rho_hat + k*' K^-1 r,  v* = 1 - k*' K^-1 k* + R'A^-1R + eta^2,  R = f* - F'K^-1 k*.
  Full predictive = sum_g w_g t_nu(m*_g, scale_g^2); the mean is sum_g w_g m*_g (nu > 1 always since
  nu >= 2 a0 = 4). Quantiles: root of the mixture CDF by scipy.optimize.brentq (in standardised units,
  then rescaled), or sampling (sample()).

The LF KRR is passed in already trained (train once, reuse), the same object the reference class
trains (LooKRR, seed 42, portion_test = 0.2 = reference lf_portion).
"""
from __future__ import annotations

import numpy as np
from scipy import special, stats
from scipy.optimize import brentq

JITTER = 1e-10


def theta_from_l(l):
    """log10 theta of mfbml_local.kernels.RBF for length scale l (k = exp(-d^2 / (2 l^2)))."""
    return -np.log10(2.0 * np.asarray(l, float) ** 2)


def l_from_theta(theta):
    return np.sqrt(1.0 / (2.0 * 10.0 ** np.asarray(theta, float)))


def default_grid(n_l=60, n_eta=40, l_med=0.3, l_sd=0.75, eta_lo=1e-3, eta_hi=0.3, mass=0.999):
    """Grid nodes (l, eta, log prior weight), flattened with l the slow index."""
    z = stats.norm.ppf(0.5 + mass / 2.0)
    logl = np.linspace(np.log(l_med) - z * l_sd, np.log(l_med) + z * l_sd, n_l)
    eta = np.exp(np.linspace(np.log(eta_lo), np.log(eta_hi), n_eta))
    LL, EE = np.meshgrid(logl, eta, indexing="ij")
    logprior = -0.5 * ((LL - np.log(l_med)) / l_sd) ** 2
    return np.exp(LL).ravel(), EE.ravel(), logprior.ravel()


class FullBayesKRRLRGP:
    def __init__(self, design_space, lf_model=None, basis="linear", a0=2.0, b0=0.1, grid=None,
                 **grid_kw):
        """lf_model: trained LF model with .predict(X) (raw inputs); None only for basis='ordinary'.
        grid: optional (l, eta, logprior) arrays to override the default 60 x 40 grid."""
        self.bounds = np.asarray(design_space, float)
        self.lf_model, self.basis = lf_model, basis
        self.a0, self.b0 = float(a0), float(b0)
        self.l, self.eta, self.logprior = (np.asarray(a, float) for a in (grid or default_grid(**grid_kw)))
        self.n_l = len(np.unique(self.l)) if grid is None else None
        self.n_eta = len(self.l) // self.n_l if self.n_l else None

    # -- helpers
    def normalize_input(self, X):
        return (np.atleast_2d(X) - self.bounds[:, 0]) / (self.bounds[:, 1] - self.bounds[:, 0])

    def _basis(self, X):
        X = np.atleast_2d(X)
        one = np.ones((X.shape[0], 1))
        if self.basis == "ordinary":
            return one
        f = np.asarray(self.lf_model.predict(X), float).reshape(-1, 1)
        return np.hstack([one, (f - self.yh_mean) / self.yh_std])

    @staticmethod
    def _d2(X, Y):
        return np.sum((X[:, None, :] - Y[None, :, :]) ** 2, axis=2)

    def _kern(self, d2):
        param = 1.0 / (2.0 * self.l ** 2)                       # (G,)
        return np.exp(-param[:, None, None] * d2[None])

    # -- training
    def train(self, xh, yh):
        xh = np.atleast_2d(np.asarray(xh, float))
        yh = np.asarray(yh, float).ravel()
        self.yh_mean, self.yh_std = float(np.mean(yh)), float(np.std(yh))
        self.xs = self.normalize_input(xh)
        self.y = (yh - self.yh_mean) / self.yh_std
        self.F = self._basis(xh)
        n, p = self.F.shape
        self.n, self.p = n, p
        G = len(self.l)
        K = self._kern(self._d2(self.xs, self.xs)) + (self.eta ** 2 + JITTER)[:, None, None] * np.eye(n)
        L = np.linalg.cholesky(K)
        Linv = np.linalg.solve(L, np.broadcast_to(np.eye(n), (G, n, n)))
        Kinv = np.swapaxes(Linv, 1, 2) @ Linv
        KF = Kinv @ self.F
        A = np.swapaxes(self.F[None], 1, 2) @ KF + JITTER * np.eye(p)
        Ainv = np.linalg.inv(A)
        rho = np.einsum("gpq,gq->gp", Ainv, np.einsum("np,gnm,m->gp", self.F, Kinv, self.y))
        r = self.y[None] - np.einsum("np,gp->gn", self.F, rho)
        gamma = np.einsum("gnm,gm->gn", Kinv, r)
        S = np.einsum("gn,gn->g", r, gamma)
        a_n, b_n = self.a0 + (n - p) / 2.0, self.b0 + S / 2.0
        logdetK = 2.0 * np.sum(np.log(np.diagonal(L, axis1=1, axis2=2)), axis=1)
        logdetA = np.linalg.slogdet(A)[1]
        self.logml = -0.5 * logdetK - 0.5 * logdetA - a_n * np.log(b_n)
        lw = self.logprior + self.logml
        self.logw = lw - special.logsumexp(lw)
        self.w = np.exp(self.logw)
        self._Kinv, self._Ainv, self._rho, self._gamma = Kinv, Ainv, rho, gamma
        self.a_n, self.b_n, self.S = a_n, b_n, S
        self.nu = 2.0 * a_n
        return self

    # -- prediction
    def predict_nodes(self, X):
        """Per-node (location, scale) in standardised units, each (G, m); and nu."""
        X = np.atleast_2d(np.asarray(X, float))
        Xn = self.normalize_input(X)
        ks = self._kern(self._d2(self.xs, Xn))                   # (G,n,m)
        fs = self._basis(X)                                      # (m,p)
        Kk = self._Kinv @ ks
        R = fs.T[None] - np.einsum("np,gnm->gpm", self.F, Kk)    # (G,p,m)
        v = (1.0 - np.einsum("gnm,gnm->gm", ks, Kk)
             + np.einsum("gpm,gpq,gqm->gm", R, self._Ainv, R) + (self.eta ** 2)[:, None])
        mu = np.einsum("mp,gp->gm", fs, self._rho) + np.einsum("gnm,gn->gm", ks, self._gamma)
        scale = np.sqrt(np.maximum((self.b_n / self.a_n)[:, None] * v, 0.0))
        return mu, scale, self.nu

    def predict_dist(self, X):
        """(mu_nodes, scale_nodes, nu, w) in original units."""
        mu, sc, nu = self.predict_nodes(X)
        return mu * self.yh_std + self.yh_mean, sc * self.yh_std, nu, self.w

    def predict_mean(self, X):
        mu, _, _, w = self.predict_dist(X)
        return w @ mu

    def cdf(self, x, X1):
        """Mixture CDF at original-units values x (m,) for inputs X1 (m, d)."""
        mu, sc, nu, w = self.predict_dist(X1)
        return np.einsum("g,gm->m", w, special.stdtr(nu, (np.atleast_1d(x)[None] - mu) / sc))

    def quantile(self, X, qs):
        """Mixture quantiles by brentq inversion of the CDF; returns (m, len(qs)) original units."""
        mu, sc, nu, w = self.predict_dist(X)
        qs = np.atleast_1d(qs)
        tc = stats.t.ppf(1.0 - 1e-12, nu)
        out = np.empty((mu.shape[1], len(qs)))
        for i in range(mu.shape[1]):
            m, s = mu[:, i], sc[:, i]
            lo, hi = np.min(m - s * tc), np.max(m + s * tc)
            xt = 1e-13 * max(abs(lo), abs(hi), s.max())
            for j, q in enumerate(qs):
                out[i, j] = brentq(lambda x: float(w @ special.stdtr(nu, (x - m) / s)) - q, lo, hi,
                                   xtol=xt, rtol=1e-14, maxiter=500)
        return out

    def sample(self, X, N, rng):
        """N draws from the mixture predictive for each row of X: (N, m)."""
        mu, sc, nu, w = self.predict_dist(X)
        g = rng.choice(len(w), size=N, p=w)
        t = rng.standard_t(nu, size=(N, mu.shape[1]))
        return mu[g] + sc[g] * t

    # -- hyperparameter posterior
    def posterior_summary(self, level=0.9):
        """Marginal posterior of l (scaled input units) and eta: median and central `level` interval."""
        def wq(v, w, qs):
            o = np.argsort(v)
            c = np.cumsum(w[o])
            return [float(v[o][min(np.searchsorted(c, q), len(v) - 1)]) for q in qs]
        qs = [0.5 - level / 2, 0.5, 0.5 + level / 2]
        lq = wq(self.l, self.w, qs)
        eq = wq(self.eta, self.w, qs)
        return dict(l_lo=lq[0], l_med=lq[1], l_hi=lq[2], eta_lo=eq[0], eta_med=eq[1], eta_hi=eq[2],
                    l_mean_log=float(np.exp(self.w @ np.log(self.l))))
