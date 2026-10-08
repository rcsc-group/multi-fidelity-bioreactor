"""T2 unit tests for the Bayesian-amplitude KRR-LR-GPR (scripts/mfbml_local/krr_lr_gpr_bayes.py)."""
import pathlib
import sys

import numpy as np
import pytest
from scipy import stats

ROOT = pathlib.Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from mfbml_local.krr_lr_gpr import KernelRidgeLinearGaussianProcess          # noqa: E402
from mfbml_local.krr_lr_gpr_bayes import KernelRidgeLinearGaussianProcessBayes  # noqa: E402
from scripts.plot_mf_fig13 import LooKRR                                       # noqa: E402

DS = np.array([[0.0, 1.0]])


def problem(n_hf=5, n_lf=12):
    xl = np.linspace(0, 1, n_lf).reshape(-1, 1)
    xh = np.linspace(0.05, 0.95, n_hf).reshape(-1, 1)
    lf = lambda x: np.sin(3 * x) + 0.3 * x
    hf = lambda x: 1.7 * lf(x) + 0.4 + 0.05 * np.cos(9 * x)
    return xh, hf(xh).ravel(), xl, lf(xl).ravel()


def make(cls, **kw):
    lf = LooKRR(design_space=DS, params_optimize=True, noise_data=True, optimizer_restart=5, seed=1)
    return cls(design_space=DS, optimizer_restart=5, lf_poly_order="linear", seed=1,
               lf_model=lf, noise_prior=None, **kw)


def trained(cls, **kw):
    xh, yh, xl, yl = problem()
    m = make(cls, **kw)
    m.train(X=[xh, xl], Y=[yh, yl])
    return m


XQ = np.linspace(0.0, 1.0, 11).reshape(-1, 1)


def test_location_equals_reference_mean():
    ref, bay = trained(KernelRidgeLinearGaussianProcess), trained(KernelRidgeLinearGaussianProcessBayes)
    np.testing.assert_allclose(bay.predict_dist(XQ)[0], ref.predict(XQ).ravel(), rtol=0, atol=1e-12)


def test_band_converges_to_reference_when_prior_is_tight():
    ref = trained(KernelRidgeLinearGaussianProcess)
    ref.predict(XQ, return_std=True)
    s2hat = float(np.squeeze(ref.sigma2))
    a0 = 1e7
    bay = trained(KernelRidgeLinearGaussianProcessBayes, a0=a0, b0=a0 * s2hat)
    mu, scale, nu = bay.predict_dist(XQ, include_noise=False)
    assert nu > 1e6
    # reference epistemic std (no noise), original units
    np.testing.assert_allclose(scale, ref.epistemic.ravel(), rtol=1e-3)
    # 95% half-widths: t -> normal
    np.testing.assert_allclose(scale * stats.t.ppf(0.975, nu), 1.96 * ref.epistemic.ravel(), rtol=2e-3)


def test_monte_carlo_matches_closed_form_student_t():
    """Sample s2 ~ IG(a_n, b_n); rho | s2 ~ N(rho_hat, s2 A^-1); y* | rho, y, s2 by direct
    Gaussian conditioning (own matrix algebra, independent of eq. (3)); compare quantiles."""
    bay = trained(KernelRidgeLinearGaussianProcessBayes, a0=2.0, b0=0.1)
    x = np.array([[0.33], [0.5], [0.98]])
    mu, scale, nu = bay.predict_dist(x)
    a_n, b_n, _ = bay.posterior_amplitude
    rng = np.random.default_rng(0)
    N = 400_000
    K, F = bay.K, bay.f
    Ki = np.linalg.inv(K)
    A = F.T @ Ki @ F
    rho_hat = np.linalg.solve(A, F.T @ Ki @ bay.sample_yh_scaled)
    Lc = np.linalg.cholesky(np.linalg.inv(A))
    nr2 = (bay.noise / bay.yh_std) ** 2
    s2 = b_n / rng.gamma(a_n, size=N)                       # IG(a_n, b_n)
    qs = [0.025, 0.25, 0.5, 0.75, 0.975]
    for i in range(len(x)):
        xi = x[i:i + 1]
        kn = bay.kernel.get_kernel_matrix(bay.sample_xh_scaled, bay.normalize_input(xi)).ravel()
        fs = bay._basis_function(xi, poly_order="linear").ravel()
        rho = rho_hat[None, :] + np.sqrt(s2)[:, None] * (rng.standard_normal((N, 2)) @ Lc.T)
        m = rho @ fs + (bay.sample_yh_scaled[None, :] - rho @ F.T) @ Ki @ kn
        var = 1.0 + nr2 - kn @ Ki @ kn
        ys = (m + np.sqrt(s2 * var) * rng.standard_normal(N)) * bay.yh_std + bay.yh_mean
        emp = np.quantile(ys, qs)
        cf = mu[i] + scale[i] * stats.t.ppf(qs, nu)
        tol = 0.02 * scale[i] * (stats.t.ppf(0.975, nu))     # heavy tails: MC error well below 2% of 95% half-width
        np.testing.assert_allclose(emp, cf, atol=tol)
