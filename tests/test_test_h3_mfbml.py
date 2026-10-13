"""Unit tests for scripts/test_h3_mfbml.py (synthetic data, no runs)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.test_h1_mfbml import h1_verdict  # noqa: E402
from scripts.test_h2_mfbml import loo_h2  # noqa: E402
from scripts.test_h3_mfbml import UNSEEN, h3_verdict, unseen_stats  # noqa: E402

BASE = dict(mf_rmse=0.05, hf_rmse=0.20, ratio_rmse=0.25, mf_cover=12, n=13)
UNS = dict(mf_rmse=0.05, n=5)


def _stats(**over):
    return dict(tau_95_t=dict(BASE), edr_mean_t=dict(BASE), tau_mean_t=dict(BASE)), \
        dict(tau_95_t=dict(UNS), edr_mean_t=dict(UNS), tau_mean_t=dict(UNS))


class StubLF:
    def predict(self, X):
        X = np.atleast_2d(X)
        return np.log(0.8) + 0.04 * X[:, 0] + 0.1 * X[:, 1]


def _synthetic():
    hf_keys = [(17.5, 7.0), (22.5, 7.0), (27.5, 7.0), (32.5, 7.0), (25.0, 4.0), (25.0, 9.0)]
    extra = [(15.0, 7.0), (20.0, 7.0), (30.0, 7.0), (37.5, 7.0)]
    q = "tau_95_t"
    f = lambda k: np.exp(0.04 * k[0] + 0.1 * k[1])  # noqa: E731
    l9 = {k: {q: f(k)} for k in hf_keys + extra}
    lf = {k: {q: 0.8 * f(k)} for k in hf_keys + extra}
    hf = {k: {q: 1.1 * f(k)} for k in hf_keys}
    return q, l9, lf, hf


def test_unseen_constant():
    assert set(UNSEEN) == {(17.5, 4.0), (25.0, 4.0), (32.5, 4.0), (17.5, 9.0), (25.0, 9.0)}
    assert len(UNSEEN) == 5


def test_verdict_pass():
    s, u = _stats()
    assert h3_verdict(s, u) == dict(verdict="PASS", failed=[])


def test_verdict_h2_bar_failure_propagates():
    s, u = _stats()
    s["edr_mean_t"]["mf_rmse"] = 0.2
    v = h3_verdict(s, u)
    assert v["verdict"] == "FAIL" and v["failed"] == h1_verdict(s)["failed"]


def test_verdict_unseen_failure():
    s, u = _stats()
    u["tau_95_t"]["mf_rmse"] = 0.11
    v = h3_verdict(s, u)
    assert v["verdict"] == "FAIL" and ("tau_95_t", "unseen rmse<=10%") in v["failed"]
    s, u = _stats()
    u["edr_mean_t"]["mf_rmse"] = 0.101
    assert ("edr_mean_t", "unseen rmse<=10%") in h3_verdict(s, u)["failed"]


def test_verdict_unseen_boundary_passes():
    s, u = _stats()
    u["tau_95_t"]["mf_rmse"] = 0.10
    assert h3_verdict(s, u)["verdict"] == "PASS"


def test_verdict_unseen_wrong_n_fails():
    s, u = _stats()
    u["tau_95_t"]["n"] = 4
    v = h3_verdict(s, u)
    assert v["verdict"] == "FAIL" and ("tau_95_t", "unseen n==5") in v["failed"]


def test_verdict_pending():
    s, u = _stats()
    for q in s:
        s[q]["n"] = 12
    u["tau_95_t"]["mf_rmse"] = 0.5   # ignored while PENDING
    assert h3_verdict(s, u)["verdict"].startswith("PENDING")
    assert h3_verdict(s, None)["verdict"].startswith("PENDING")


def test_h2_default_path_unchanged():
    """loo_h2 default (no log) must equal loo_h2(log=False) and keep its original numbers."""
    q, l9, lf, hf = _synthetic()
    a, ra = loo_h2(q, StubLF(), l9, lf, hf)
    b, rb = loo_h2(q, StubLF(), l9, lf, hf, log=False)
    assert a == b and ra == rb


def test_log_round_trip_recovers_constant_ratio_truth():
    q, l9, lf, hf = _synthetic()
    s, rows = loo_h2(q, StubLF(), l9, lf, hf, log=True)
    assert s["n"] == 6
    assert s["mf_rmse"] < 0.05
    for r in rows:
        assert r["mf_lo"] <= r["mf_pred"] <= r["mf_hi"]
        assert r["mf_lo"] > 0 and r["obs"] == pytest.approx(1.1 * np.exp(0.04 * r["rpm"] + 0.1 * r["theta"]))
        assert r["ratio_pred"] > 0 and r["l9"] == pytest.approx(np.exp(0.04 * r["rpm"] + 0.1 * r["theta"]))
    # non-decisive columns stay on the raw scale
    assert s["l9_rmse"] == pytest.approx(abs(1 / 1.1 - 1.0), rel=1e-6)


def test_unseen_stats_subset():
    rows = [dict(rpm=r, theta=t, obs=1.0, mf_pred=1.0 + (0.1 if (r, t) in UNSEEN else 0.0))
            for r, t in [(17.5, 7.0), (25.0, 7.0)] + sorted(UNSEEN)]
    u = unseen_stats(rows)
    assert u["n"] == 5 and u["mf_rmse"] == pytest.approx(0.1)


def test_h2_default_matches_direct_fullbayes_calls():
    """Reference built straight from FullBayesKRRLRGP (raw scale), independent of _fit_predict."""
    from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP
    from scripts.test_h1_mfbml import BOUNDS
    from scripts.test_h2_mfbml import MeanAdapter
    q, l9, lf, hf = _synthetic()
    k9 = sorted(l9, key=lambda k: (k[1], k[0]))
    m1 = FullBayesKRRLRGP(BOUNDS, lf_model=StubLF(), basis="linear").train(
        np.array(k9), np.array([l9[k][q] for k in k9]))
    hk = sorted(hf, key=lambda k: (k[1], k[0]))
    X = np.array(hk)
    yh = np.array([hf[k][q] for k in hk])
    tr = np.arange(len(hk)) != 0
    m2 = FullBayesKRRLRGP(BOUNDS, lf_model=MeanAdapter(m1), basis="linear").train(X[tr], yh[tr])
    qq = m2.quantile(X[:1], [0.025, 0.975])[0]
    _, rows = loo_h2(q, StubLF(), l9, lf, hf)
    assert rows[0]["mf_pred"] == float(m2.predict_mean(X[:1])[0])
    assert (rows[0]["mf_lo"], rows[0]["mf_hi"]) == (float(qq[0]), float(qq[1]))


def test_log_mode_trains_every_stage_on_log_targets(monkeypatch):
    """Rule: model 1 on log L9, model 2 and HF-only on log L10. Spy on every FullBayesKRRLRGP.train call."""
    from mfbml_local import krr_lr_gpr_fullbayes as fb
    seen = []
    orig = fb.FullBayesKRRLRGP.train

    def spy(self, xh, yh):
        seen.append(np.asarray(yh, float).copy())
        return orig(self, xh, yh)

    monkeypatch.setattr(fb.FullBayesKRRLRGP, "train", spy)
    q, l9, lf, hf = _synthetic()
    loo_h2(q, StubLF(), l9, lf, hf, log=True)
    k9 = sorted(l9, key=lambda k: (k[1], k[0]))
    np.testing.assert_allclose(seen[0], np.log([l9[k][q] for k in k9]))      # model 1
    log_hf = np.log([hf[k][q] for k in sorted(hf, key=lambda k: (k[1], k[0]))])
    for y in seen[1:]:                                                        # model 2 / HF-only folds
        assert np.all(np.isin(np.round(y, 12), np.round(log_hf, 12)))
