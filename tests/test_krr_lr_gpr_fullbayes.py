"""T2b unit tests for the fully Bayesian KRR-LR-GPR (scripts/mfbml_local/krr_lr_gpr_fullbayes.py)."""
import pathlib
import sys

import numpy as np
from scipy import stats

ROOT = pathlib.Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from mfbml_local.krr_lr_gpr_bayes import KernelRidgeLinearGaussianProcessBayes     # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import (FullBayesKRRLRGP, l_from_theta,       # noqa: E402
                                              theta_from_l, default_grid)
from scripts.plot_mf_fig13 import LooKRR                                            # noqa: E402

DS = np.array([[0.0, 1.0]])
XQ = np.linspace(0.0, 1.0, 11).reshape(-1, 1)


def problem(n_hf=5, n_lf=12):
    xl = np.linspace(0, 1, n_lf).reshape(-1, 1)
    xh = np.linspace(0.05, 0.95, n_hf).reshape(-1, 1)
    lf = lambda x: np.sin(3 * x) + 0.3 * x
    hf = lambda x: 1.7 * lf(x) + 0.4 + 0.05 * np.cos(9 * x)
    return xh, hf(xh).ravel(), xl, lf(xl).ravel()


def test_theta_conversion_roundtrip_and_kernel():
    l = np.array([0.05, 0.3, 2.0])
    np.testing.assert_allclose(l_from_theta(theta_from_l(l)), l, rtol=1e-12)
    # RBF value at distance d equals exp(-d^2 / (2 l^2))
    from mfbml_local.kernels import RBF
    k = RBF(theta=np.array([theta_from_l(0.3)]))
    d = 0.4
    kv = k.get_kernel_matrix(np.array([[0.0]]), np.array([[d]]))[0, 0]
    np.testing.assert_allclose(kv, np.exp(-d ** 2 / (2 * 0.3 ** 2)), rtol=1e-8)  # reference adds 1e-10 jitter


def test_default_grid_shape_and_range():
    l, eta, lp = default_grid()
    assert l.size == eta.size == lp.size == 2400
    np.testing.assert_allclose([l.min(), l.max()], 0.3 * np.exp(np.array([-1, 1]) * 3.2905 * 0.75), rtol=1e-3)
    np.testing.assert_allclose([eta.min(), eta.max()], [1e-3, 0.3])


def test_single_node_reproduces_T2():
    """One grid node at the T2 model's own (theta, noise ratio) -> identical Student-t."""
    xh, yh, xl, yl = problem()
    lf = LooKRR(design_space=DS, params_optimize=True, noise_data=True, optimizer_restart=5, seed=1)
    t2 = KernelRidgeLinearGaussianProcessBayes(
        design_space=DS, optimizer_restart=5, lf_poly_order="linear", seed=1, lf_model=lf,
        noise_prior=0.02 * np.std(yh), a0=2.0, b0=0.1)
    t2.train(X=[xh, xl], Y=[yh, yl])
    mu2, sc2, nu2 = t2.predict_dist(XQ)
    l = float(l_from_theta(np.log10(np.ravel(t2.kernel.param)[0])))
    eta = t2.noise / t2.yh_std
    fb = FullBayesKRRLRGP(DS, lf_model=t2.lf_model, basis="linear", a0=2.0, b0=0.1,
                          grid=(np.array([l]), np.array([eta]), np.array([0.0])))
    fb.train(xh, yh)
    mu, sc, nu, w = fb.predict_dist(XQ)
    assert w.shape == (1,) and abs(w[0] - 1.0) < 1e-15
    np.testing.assert_allclose(mu[0], mu2, rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose(sc[0], sc2, rtol=1e-6)
    assert nu == nu2
    np.testing.assert_allclose(fb.predict_mean(XQ), mu2, rtol=1e-7, atol=1e-9)
    qs = [0.025, 0.5, 0.975]
    np.testing.assert_allclose(fb.quantile(XQ, qs), mu2[:, None] + sc2[:, None] * stats.t.ppf(qs, nu2),
                               rtol=1e-6, atol=1e-8)


def _gp_data(l_true, n, seed, eta=0.01):
    rng = np.random.default_rng(seed)
    x = np.sort(rng.uniform(0, 1, n))
    K = np.exp(-(x[:, None] - x[None]) ** 2 / (2 * l_true ** 2)) + 1e-8 * np.eye(n)
    y = np.linalg.cholesky(K) @ rng.standard_normal(n)
    return x.reshape(-1, 1), y + eta * np.std(y) * rng.standard_normal(n)


def test_posterior_over_l_concentrates_near_truth():
    """n = 15 GP draws with l_true = 0.15 and 0.5: posterior median of l within a factor 2, and the
    posterior is narrower than the prior (90% interval spans < 0.75 * prior 90% log-range)."""
    prior_logrange = 2 * 1.645 * 0.75
    for l_true, seed in ((0.15, 3), (0.5, 4)):
        x, y = _gp_data(l_true, 15, seed)
        fb = FullBayesKRRLRGP(DS, basis="ordinary").train(x, y)
        s = fb.posterior_summary(0.9)
        assert l_true / 2 <= s["l_med"] <= l_true * 2, (l_true, s)
        assert np.log(s["l_hi"] / s["l_lo"]) < 0.75 * prior_logrange, (l_true, s)


def test_mixture_quantiles_match_sampling():
    x, y = _gp_data(0.2, 8, 7)
    fb = FullBayesKRRLRGP(DS, basis="ordinary").train(x, y)
    xq = np.array([[0.1], [0.5], [0.93]])
    qs = [0.025, 0.5, 0.975]
    inv = fb.quantile(xq, qs)
    N = 400_000
    draws = fb.sample(xq, N, np.random.default_rng(0))
    # CDF residual of the sampling quantile: |F(emp) - q| must be within 5 sigma of the binomial MC error.
    for i in range(len(xq)):
        emp = np.quantile(draws[:, i], qs)
        Fe = np.array([fb.cdf(e, xq[i:i + 1])[0] for e in emp])
        se = np.sqrt(np.array(qs) * (1 - np.array(qs)) / N)
        assert np.all(np.abs(Fe - np.array(qs)) < 5 * se), (Fe, qs)
        # and the inverted quantile really inverts the CDF
        Fi = np.array([fb.cdf(v, xq[i:i + 1])[0] for v in inv[i]])
        np.testing.assert_allclose(Fi, qs, atol=1e-9)
