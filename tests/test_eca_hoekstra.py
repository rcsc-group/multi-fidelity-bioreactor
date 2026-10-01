"""scripts/eca_hoekstra.py must implement Eca & Hoekstra (2014) Appendix A.

Each test cites the passage of J. Comput. Phys. 262:104-130 it checks. The
paper was read 2026-10-01 from the publisher PDF on marin.nl.
The previous implementation (diag_num_uncertainty.fits) failed several of
these: no weighted fits, no p-dependent choice of estimators, F_s = 1.25 for
fixed-order fits, epsilon taken from the data instead of the fit, and a
late-bound exponent that evaluated a selected p=1 fit as h^2.
"""
import numpy as np
import pytest

from scripts.eca_hoekstra import numerical_uncertainty

H = np.array([1.0, 2.0, 4.0, 8.0])          # finest first


def test_requires_at_least_four_grids():
    # App. A: "n_g is supposed to be at least 4"
    with pytest.raises(ValueError):
        numerical_uncertainty(H[:3], 1 + 0.1 * H[:3])


def test_exact_power_law_in_range_is_gci_with_fs_1_25():
    # Sec. 3: monotonic smooth data with 0.5 <= p < 2.1 -> "reduces to the GCI",
    # F_s = 1.25, U = F_s * |alpha h_1^p| (sigma = 0, no residual).
    r = numerical_uncertainty(H, 1 + 0.2 * H ** 1.5)
    assert r["estimator"] == "RE"
    assert r["p"] == pytest.approx(1.5, rel=1e-3)
    assert r["phi0"] == pytest.approx(1.0, abs=1e-6)
    assert r["Fs"] == 1.25
    assert r["U"] == pytest.approx(1.25 * 0.2, rel=1e-3)


def test_p_above_two_uses_the_four_fixed_order_fits_and_fs_3():
    # App. A step 1: p > 2 -> delta_1 and delta_2, with and without weights (4 fits).
    # Step 3: F_s = 1.25 only if 0.5 <= p < 2.1.
    r = numerical_uncertainty(H, 1 + 0.01 * H ** 3)
    assert r["p"] > 2.1
    assert {c[0] for c in r["candidates"]} == {"1", "2"}
    assert len(r["candidates"]) == 4
    assert r["Fs"] == 3.0


def test_p_below_half_uses_the_six_fits_including_two_term():
    # App. A step 1: p < 0.5 -> delta_1, delta_2, delta_12, with and without weights.
    r = numerical_uncertainty(H, 1 + 0.3 * H ** 0.3)
    assert {c[0] for c in r["candidates"]} == {"1", "2", "12"}
    assert len(r["candidates"]) == 6
    assert r["Fs"] == 3.0


def test_divergent_data_is_anomalous_and_uses_six_fits():
    # Sec. 2.4: "If both give negative p, the data behaviour is classified as anomalous"
    r = numerical_uncertainty(H, 1 + 0.5 / H)
    assert r["p"] is None
    assert len(r["candidates"]) == 6
    assert r["Fs"] == 3.0


def test_error_estimate_comes_from_the_fit_not_the_data():
    # Eq. (1): eps(phi_i) = phi_i - phi_o := alpha h_i^p, i.e. the FIT at h_i;
    # Eq. (20) adds |phi_i - phi_fit| separately.
    rng = np.random.default_rng(0)
    q = 1 + 0.2 * H ** 1.5 + rng.normal(0, 0.003, H.size)
    r = numerical_uncertainty(H, q)
    fit1 = r["fit"](H[0])
    assert r["eps"] == pytest.approx(abs(fit1 - r["phi0"]))
    assert r["resid"] == pytest.approx(abs(q[0] - fit1))


def test_selected_fixed_order_fit_keeps_its_own_exponent():
    # The bug that motivated this module: a p=1 fit evaluated as h^2 after return.
    q = 1 + 0.1 * H ** 0.2          # p < 0.5 -> fixed-order candidates
    r = numerical_uncertainty(H, q)
    f = r["fit"]
    # a fit of the form c0 + c1 h (+ c2 h^2) is exactly reproduced by its own coefficients
    assert f(0.0) == pytest.approx(r["phi0"])
    for h in H:
        assert np.isfinite(f(h))
    if r["estimator"] == "1":
        assert f(2.0) - f(1.0) == pytest.approx(f(1.0) - f(0.0))


def test_bad_fit_branch_inflates_uncertainty():
    # Eq. (21): sigma >= Delta_phi -> U = 3 (sigma/Delta_phi)(eps + sigma + |phi_i - phi_fit|)
    q = np.array([1.00, 1.10, 0.95, 1.08, 0.97])
    h = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    r = numerical_uncertainty(h, q)
    if r["sigma"] >= r["data_range"]:
        expect = 3 * r["sigma"] / r["data_range"] * (r["eps"] + r["sigma"] + r["resid"])
        assert r["U"] == pytest.approx(expect)


def test_uncertainty_at_a_coarser_grid_uses_that_grid():
    # Sec. 3: "the determination of U_phi, usually for the finest grid ... but
    # essentially for any phi_i of a data set"
    q = 1 + 0.2 * H ** 1.5
    r1 = numerical_uncertainty(H, q, index=1)
    assert r1["U"] == pytest.approx(1.25 * 0.2 * 2.0 ** 1.5, rel=1e-3)


def test_weights_are_inverse_cell_size_normalised():
    # Eq. (16)-(17): w_i = (1/h_i) / sum(1/h_j), sum w = 1
    from scripts.eca_hoekstra import weights
    w = weights(H)
    assert w.sum() == pytest.approx(1.0)
    assert w[0] / w[1] == pytest.approx(2.0)
