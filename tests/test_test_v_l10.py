import numpy as np

from scripts.test_v_l10 import binned, verdict


def test_verdict_thresholds():
    assert verdict(10.0) == "NUMERICAL-DOMINATED"
    assert verdict(2.0) == "PHYSICAL"
    assert verdict(5.0) == "MIXED"


def test_binned_uses_window_and_skips_nan():
    chi = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.2])
    r = np.array([999.0, 10.0, 20.0, 30.0, 999.0, np.nan])
    b = binned(chi, r)
    assert b["window"] == 20.0 and b["n_window"] == 3
    assert b["bins"] == [10.0, 20.0, 30.0]
