"""Unit tests for scripts/test_h2_mfbml.py (synthetic data, no runs, no real models)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.test_h1_mfbml import BOUNDS, h1_verdict  # noqa: E402
from scripts.test_h2_mfbml import (  # noqa: E402
    MeanAdapter, h2_verdict, l9_points, l9_run, loo_h2,
)


class StubLF:
    """Stands in for the L8 KRR: smooth function of (rpm, theta)."""

    def predict(self, X):
        X = np.atleast_2d(X)
        return 0.05 * X[:, 0] + 0.3 * X[:, 1]


class StubModel:
    def predict_mean(self, X):
        return np.atleast_2d(X)[:, 0] * 2.0


def test_adapter_shape_and_values():
    X = np.array([[15.0, 7.0], [20.0, 4.0], [30.0, 9.0]])
    out = MeanAdapter(StubModel()).predict(X)
    assert out.shape == (3,)
    assert np.allclose(out, [30.0, 40.0, 60.0])
    assert MeanAdapter(StubModel()).predict(X[:1]).shape == (1,)


def test_adapter_works_inside_fullbayes():
    from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP
    X = np.array([[15, 4], [20, 7], [25, 5.5], [30, 9], [35, 7], [37.5, 4.0]], float)
    y = 1.0 + 0.1 * X[:, 0]
    m = FullBayesKRRLRGP(BOUNDS, lf_model=MeanAdapter(StubModel()), basis="linear").train(X, y)
    assert m.predict_mean(X[:2]).shape == (2,)


def test_l9_run_ids():
    assert l9_run(15, 7) == "fig9_l9_rpm15"
    assert l9_run(17.5, 7.0) == "fig9_l9_rpm17.5"
    assert l9_run(37.5, 7) == "fig9_l9_rpm37.5"
    assert l9_run(25, 4.0) == "x8_l9_rpm25_th4"
    assert l9_run(32.5, 9) == "x8_l9_rpm32.5_th9"


def test_l9_points():
    pts = l9_points()
    assert len(pts) == 10 + 6
    assert (15.0, 7.0) in pts and (37.5, 7.0) in pts
    assert (17.5, 4.0) in pts and (32.5, 9.0) in pts
    assert len(set(pts)) == len(pts)


def test_verdict_is_h1_rule():
    q = dict(mf_rmse=0.05, hf_rmse=0.20, ratio_rmse=0.25, mf_cover=12, n=13)
    s = dict(tau_95_t=dict(q), edr_mean_t=dict(q))
    assert h2_verdict(s) == h1_verdict(s)
    s["edr_mean_t"]["mf_rmse"] = 0.2
    assert h2_verdict(s) == h1_verdict(s) and h2_verdict(s)["verdict"] == "FAIL"
    q["n"] = 12
    assert h2_verdict(dict(tau_95_t=dict(q), edr_mean_t=dict(q)))["verdict"].startswith("PENDING")


def test_loo_exact_scaling_gives_small_error():
    rng = np.random.default_rng(0)
    hf_keys = [(17.5, 7.0), (22.5, 7.0), (27.5, 7.0), (32.5, 7.0), (25.0, 4.0), (25.0, 9.0)]
    extra = [(15.0, 7.0), (20.0, 7.0), (30.0, 7.0), (37.5, 7.0)]
    q = "tau_95_t"
    f = lambda k: 1.0 + 0.04 * k[0] + 0.1 * k[1]  # noqa: E731
    l9 = {k: {q: f(k) * (1 + 0.0 * rng.normal())} for k in hf_keys + extra}
    lf = {k: {q: 0.8 * f(k)} for k in hf_keys + extra}
    hf = {k: {q: 1.1 * l9[k][q]} for k in hf_keys}
    s, rows = loo_h2(q, StubLF(), l9, lf, hf)
    assert s["n"] == 6 and len(rows) == 6
    assert s["mf_rmse"] < 0.05
    assert s["l9_ratio_rmse"] < 1e-9          # L10 = 1.1 * L9 exactly -> ratio baseline is exact
    assert s["l9_rmse"] == pytest.approx(abs(1 / 1.1 - 1.0), rel=1e-6)
    assert s["n_l9_model1"] == len(l9)
    for key in ("mf_cover", "hf_cover", "hf_rmse", "ratio_rmse"):
        assert key in s
