"""Unit tests for the pre-registered test K rule (pure function, synthetic data)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from test_k_diffusivity import evaluate_rule  # noqa: E402

TH = ("dtmix_0.50", "dtmix_0.75", "dtmix_0.95")


def run(v):
    return {t: v for t in TH}


def full(l7s10=1.0, l8s10=1.0, l9s5=1.0, l9s10=1.0):
    """dtmix = 10 everywhere unless overridden by a multiplier on the s=10 / L9 s=1e5 run."""
    d = {7: {0: run(10.0), 1: run(10.0 * l7s10)},
         8: {0: run(10.0), 1: run(10.0 * l8s10), 5: run(4.0)},
         9: {0: run(12.0), 1: run(12.0 * l9s10), 3: run(11.0), 5: run(4.0 * l9s5)}}
    return d


def test_pass():
    r = evaluate_rule(full())
    assert r["verdict"] == "NOT WELL-POSED AT PHYSICAL D"
    assert r["cond_i"][7]["holds"] and r["cond_i"][8]["holds"]
    assert r["cond_i"][9]["holds"]  # amended 2026-10-10: k_l9_s1e1 added
    assert r["cond_ii"]["holds"]


def test_reject():
    r = evaluate_rule(full(l8s10=1.4))
    assert r["verdict"] == "REJECTED"


def test_pending_missing_runs():
    d = full()
    del d[9]
    r = evaluate_rule(d)
    assert r["verdict"] == "PENDING"
    assert any("L9" in m for m in r["missing"])


def test_pending_null_value():
    d = full()
    d[8][5]["dtmix_0.95"] = None
    assert evaluate_rule(d)["verdict"] == "PENDING"


def test_inconclusive_names_failing_condition():
    r = evaluate_rule(full(l9s5=1.5))  # levels do not agree at s=1e5
    assert r["verdict"] == "INCONCLUSIVE"
    assert "(ii)" in " ".join(r["failing"])
    r = evaluate_rule(full(l7s10=1.15))  # 15% change at L7: (i) fails, but < 30% so not rejected
    assert r["verdict"] == "INCONCLUSIVE"
    assert any("(i) L7" in f for f in r["failing"])


def test_boundaries():
    # (i) strict <0.10, (ii) <=0.10
    assert evaluate_rule(full(l7s10=1.10))["cond_i"][7]["holds"] is False
    assert evaluate_rule(full(l9s5=1.10))["cond_ii"]["holds"] is True


def test_l9_s10_fails_inconclusive():
    r = evaluate_rule(full(l9s10=1.2))
    assert r["verdict"] == "INCONCLUSIVE"
    assert "(i) L9" in r["failing"]


def test_l9_s10_missing_pending():
    d = full()
    del d[9][1]
    assert evaluate_rule(d)["verdict"] == "PENDING"
