"""Unit tests for test I (islands): a50 and the verdict rule (pure functions, synthetic data)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from test_i_islands import a50, island_rule, jaccard  # noqa: E402


def test_a50_uniform_magnitude():
    c = np.zeros((20, 20))
    c[::2, :] = 1.0                      # c' = +-0.5 everywhere
    assert a50(c, np.ones_like(c)) == pytest.approx(0.5, abs=1e-9)


def test_a50_compact_blob():
    c = np.zeros((100, 100))
    c[:20, :20] = 1.0                    # 4% of the area
    assert a50(c, np.ones_like(c)) == pytest.approx(0.02, abs=0.003)


def test_a50_ignores_zero_weight_cells():
    c = np.zeros((10, 10))
    c[:5] = 1.0
    w = np.zeros_like(c)
    w[:, :] = 1.0
    w[5:] = 0.0                          # only the c = 1 half is liquid: c' = 0 there
    assert a50(c, w) == 0.0 or np.isnan(a50(c, w))


def S(a90, a95):
    return {"0.90": a90, "0.95": a95}


def test_rule_island():
    r = island_rule(S(0.1, 0.1), [S(0.3, 0.3), S(0.25, 0.4)], 0.8)
    assert r == "ISLAND"


def test_rule_island_needs_persistence():
    assert island_rule(S(0.1, 0.1), [S(0.3, 0.3), S(0.3, 0.3)], 0.2) == "INCONCLUSIVE"


def test_rule_island_needs_both_chi_levels():
    assert island_rule(S(0.1, 0.25), [S(0.3, 0.3), S(0.3, 0.3)], 0.8) == "INCONCLUSIVE"


def test_rule_not_island():
    assert island_rule(S(0.3, 0.31), [S(0.3, 0.3), S(0.28, 0.33)], 0.5) == "NOT ISLAND"


def test_rule_not_island_boundary_fails_one_fast():
    assert island_rule(S(0.3, 0.3), [S(0.3, 0.3), S(0.15, 0.3)], 0.5) == "INCONCLUSIVE"


def test_rule_pending():
    assert island_rule(None, [S(0.3, 0.3), S(0.3, 0.3)], 0.05) == "PENDING"
    assert island_rule(S(0.1, 0.1), [S(0.3, 0.3), None], 0.05) == "PENDING"


def test_rule_nan_is_inconclusive():
    assert island_rule(S(float("nan"), 0.1), [S(0.3, 0.3), S(0.3, 0.3)], 0.05) == "INCONCLUSIVE"


def test_jaccard():
    a = np.zeros((10, 10), bool)
    b = np.zeros((10, 10), bool)
    a[:4] = True
    b[2:6] = True
    assert jaccard(a, a) == 1.0
    assert jaccard(a, b) == pytest.approx(20 / 60)
    assert jaccard(a, np.zeros_like(a)) == 0.0


def test_rule_island_symmetric_blobs_need_overlap_not_centroid():
    # persistence is now the min Jaccard overlap (>= 0.5 needed); a low overlap blocks ISLAND
    assert island_rule(S(0.1, 0.1), [S(0.3, 0.3), S(0.3, 0.3)], 0.3) == "INCONCLUSIVE"
