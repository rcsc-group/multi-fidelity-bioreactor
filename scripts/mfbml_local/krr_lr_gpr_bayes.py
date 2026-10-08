"""KRR-LR-GPR with a Bayesian (inverse-gamma) prior on the GP amplitude (T2, diary 2026-10-08).

Pre-registration quoted from diary.md (2026-10-08):
  "T2 (zero compute): prior on the residual amplitude in KRR-LR-GPR (fully Bayesian over the GP
   amplitude, replacing the plug-in sigma2 that collapses with 1 residual dof). PASS = pooled 95%
   coverage >= 85% on test A (tau, all 4 cases, n_HF = 3) AND median rel RMSE <= 1.2x the test-A
   value in every case."

Problem: the reference predict() multiplies the kriging variance by the plug-in
sigma2_hat = r'K^-1 r / n. With n = 3 and p = 2 basis terms one residual dof remains and
sigma2_hat -> ~0, so the band collapses.

Model (standardised units, as in the reference class; y is the standardised HF response)
  y = F rho + Z,  Z ~ GP(0, s2 * K),  K = k(X,X) + nr^2 I,   F = [1, f_l(X)]  (n x p, p = 2)
  rho ~ flat (improper),  s2 ~ IG(a0, b0)  (a0 = 2, b0 = 0.1 by default)
  Kernel hyperparameters and the noise ratio nr = noise / yh_std are those fitted by the reference
  eq-23 optimiser; train() is inherited unchanged. The reference estimate of beta is the GLS estimate.

Derivation (derived here)
  Let A = F' K^-1 F, rho_hat = A^-1 F' K^-1 y (GLS), r = y - F rho_hat, S = r' K^-1 r.
  Likelihood: p(y | rho, s2) = (2 pi s2)^(-n/2) |K|^(-1/2) exp(-(y-F rho)' K^-1 (y-F rho) / (2 s2)).
  Completing the square in rho:
     (y-F rho)' K^-1 (y-F rho) = S + (rho - rho_hat)' A (rho - rho_hat).
  Integrating the flat prior over rho (p dims) gives (2 pi s2)^(p/2) |A|^(-1/2), so
     p(y | s2) ∝ s2^(-(n-p)/2) exp(-S / (2 s2)).
  With the IG(a0, b0) prior  p(s2) ∝ s2^(-a0-1) exp(-b0/s2):
     s2 | y ~ IG(a_n, b_n),   a_n = a0 + (n-p)/2,   b_n = b0 + S/2.             (1)
  Conditional on s2, rho | y, s2 ~ N(rho_hat, s2 A^-1). For a new input x with kernel vector
  k* (n), basis f* (p), the new (noisy) response y*, with rho integrated out (universal kriging):
     y* | y, s2 ~ N(m*, s2 v*),
     m*  = f*' rho_hat + k*' K^-1 r,                                              (2)
     v*  = 1 - k*' K^-1 k* + R' A^-1 R + nr^2,   R = f* - F' K^-1 k*.            (3)
  (1 = k(x,x) prior variance; the R term is the basis-uncertainty term; nr^2 is the observation
  noise, which lives inside K and so is also multiplied by s2.)
  Mixing a Normal scale family over s2 ~ IG(a_n, b_n) gives a Student-t:
     y* | y ~ t_nu( m*, scale^2 = (b_n / a_n) v* ),   nu = 2 a_n = n - p + 2 a0.  (4)
  (Proof: integral of N(m, s2 v) IG(s2; a, b) ds2 is the t density with 2a dof and scale^2 = b v / a.)
  Location m* equals the reference mean exactly (inherited predict()); only the band changes.
  Limit: a0 -> inf with b0 = a0 * s2_hat gives b_n/a_n -> s2_hat = S/n, nu -> inf, so (4) -> the
  reference band (up to the noise convention: the reference adds noise^2 unscaled by s2_hat, here
  nr^2 is inside v*; the `include_noise=False` option removes it for comparison).
Scale back to original units: location * yh_std + yh_mean, scale * yh_std.

Usage: model = KernelRidgeLinearGaussianProcessBayes(design_space=..., lf_model=..., a0=2, b0=0.1);
model.train(...); mu, scale, nu = model.predict_dist(X); quantile = mu + scale * scipy.stats.t.ppf(q, nu).
"""
from __future__ import annotations

import numpy as np
from numpy.linalg import solve
from scipy import stats

from .krr_lr_gpr import KernelRidgeLinearGaussianProcess


class KernelRidgeLinearGaussianProcessBayes(KernelRidgeLinearGaussianProcess):
    def __init__(self, *args, a0: float = 2.0, b0: float = 0.1, **kw) -> None:
        super().__init__(*args, **kw)
        self.a0, self.b0 = float(a0), float(b0)

    @property
    def posterior_amplitude(self):
        """(a_n, b_n, S) of eq. (1), standardised units."""
        n, p = self.f.shape
        r = self.sample_yh_scaled - np.dot(self.f, self.beta)
        S = float(np.dot(r.ravel(), np.asarray(self.gamma).ravel()))
        return self.a0 + (n - p) / 2.0, self.b0 + S / 2.0, S

    def predict_dist(self, X: np.ndarray, include_noise: bool = True):
        """Return (location, scale, nu) in original units; y* ~ t_nu(location, scale^2)."""
        mu = super().predict(X).ravel()                      # reference mean, eq. (2)
        Xn = np.atleast_2d(self.normalize_input(X))
        knew = self.kernel.get_kernel_matrix(self.sample_xh_scaled, Xn)
        f = self._basis_function(X, poly_order=self.lf_poly_order)
        delta = solve(self.L.T, solve(self.L, knew))
        R = f.T - np.dot(self.f.T, delta)
        v = (1 - np.diag(np.dot(knew.T, delta))
             + np.diag(R.T.dot(solve(self.ld.T, solve(self.ld, R)))))
        if include_noise:
            v = v + (self.noise / self.yh_std) ** 2
        a_n, b_n, _ = self.posterior_amplitude
        scale = np.sqrt(np.maximum((b_n / a_n) * v, 0.0)) * self.yh_std
        return mu, scale, 2.0 * a_n

    def predict_interval_halfwidth(self, X, level: float = 0.95):
        mu, scale, nu = self.predict_dist(X)
        return mu, scale * stats.t.ppf(0.5 + level / 2.0, nu)
