"""LooKRR with the closed-form leave-one-out error (one matrix inverse per evaluation instead of n Cholesky solves).
For a kernel ridge / GP mean with K = k(X, X) + noise I, the leave-one-out prediction is mu_{-i} = y_i - [K^-1 y]_i / [K^-1]_ii (exact).
Same objective as scripts/plot_mf_fig13.LooKRR; test: scripts/fast_loo_check.py (values agree to 1e-8, fits agree). Usage: from scripts.fast_loo_krr import FastLooKRR
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from scripts.plot_mf_fig13 import LooKRR


class FastLooKRR(LooKRR):
    def _loo_objective(self, params, X, Y, n):
        self._set_kernel_params(params=params)
        K = self._training_kernel_matrix(scaled_x=X, scaled_noise_std=self.noise_std / self.y_std)
        Kinv = np.linalg.inv(K)
        resid = (Kinv @ Y) / np.diag(Kinv)[:, None] if np.ndim(Y) > 1 else (Kinv @ Y) / np.diag(Kinv)
        return float(np.sum(np.square(resid)) / n)

    def _optimize_kernel_params(self) -> None:
        X, Y, n = self.sample_x_scaled, self.sample_y_scaled, len(self.sample_x_scaled)
        bounds = self._bound_definition_for_optimization()
        best, rng = None, np.random.default_rng(self.seed)
        for _ in range(self.optimizer_restart + 1):
            r = minimize(lambda p: self._loo_objective(p, X, Y, n), rng.uniform(bounds[:, 0], bounds[:, 1]), method="l-bfgs-b", bounds=bounds.tolist())
            if best is None or r.fun < best.fun:
                best = r
        self._set_kernel_params(params=best.x)
