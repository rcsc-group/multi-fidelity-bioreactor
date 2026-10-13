"""Unit tests for the pure logic of scripts/test_h1_mfbml.py (no data, no models)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.test_h1_mfbml import h1_verdict, ratio_baseline, rel_rmse_ratio  # noqa: E402


def good(n=13):
    q = dict(mf_rmse=0.05, hf_rmse=0.20, ratio_rmse=0.25, mf_cover=12, n=n)
    return dict(tau_95_t=dict(q), edr_mean_t=dict(q))


def test_pass():
    v = h1_verdict(good())
    assert v["verdict"] == "PASS" and v["failed"] == []


@pytest.mark.parametrize("key,val,cond", [
    ("mf_rmse", 0.11, "rmse<=10%"),
    ("hf_rmse", 0.04, "beats HF-only"),
    ("ratio_rmse", 0.04, "beats ratio"),
    ("mf_cover", 10, "coverage>=11/13"),
])
def test_each_failure(key, val, cond):
    s = good()
    s["edr_mean_t"][key] = val
    v = h1_verdict(s)
    assert v["verdict"] == "FAIL"
    assert v["failed"] == [("edr_mean_t", cond)]


def test_boundaries_pass():
    s = good()
    s["tau_95_t"].update(mf_rmse=0.10, mf_cover=11)
    s["tau_95_t"]["hf_rmse"] = s["tau_95_t"]["ratio_rmse"] = 0.1000001
    assert h1_verdict(s)["verdict"] == "PASS"


def test_pending():
    v = h1_verdict(good(12))
    assert v["verdict"] == "PENDING (12/13 HF points)"


def test_ratio_baseline_toy():
    # train HF/LF ratios: 2/1 and 6/3 -> mean 2 ; held-out LF 5 -> 10
    assert ratio_baseline(np.array([2.0, 6.0]), np.array([1.0, 3.0]), 5.0) == pytest.approx(10.0)
    # unequal ratios 2 and 3 -> mean 2.5
    assert ratio_baseline(np.array([2.0, 6.0]), np.array([1.0, 2.0]), 4.0) == pytest.approx(10.0)


def test_rel_rmse():
    assert rel_rmse_ratio(np.array([1.1, 0.9]), np.array([1.0, 1.0])) == pytest.approx(0.1)
