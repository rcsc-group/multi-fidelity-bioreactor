"""Chained segments must be stitched automatically, not by hand.

A continuation segment restores an already-evolved state, so on its own it is
not a complete experiment: a later segment's oxygen curve starts near
saturation and crosses every kLa threshold at its first row, and a later
segment's tracer is already partly mixed so chi never starts at 0. Computing
either KPI from one segment in isolation yields a number that looks ordinary
and is meaningless.

postprocess therefore walks `_parent_run` backwards and concatenates the raw
series before computing anything. These tests pin that walk, because the
failure it prevents is silent -- the wrong answer has the right shape.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.postprocess import _compute_c_star, _segment_chain  # noqa: E402

PARAMS = {"omega_b": 2 * math.pi * 32.5 / 60.0,
          "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
          "theta_max": [7.0, 0.0, 0.0]}
V_LIQ, ALPHA = 0.28, 1.0 / 30.0


def _seg(root: Path, name: str, t, c_star, parent: str | None = None) -> Path:
    """One segment directory whose oxygen series encodes a prescribed C*."""
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    p = dict(PARAMS, run_id=name)
    if parent:
        p["_parent_run"] = parent
    (d / "params.json").write_text(json.dumps(p))
    oxy = np.asarray(c_star) * ALPHA * V_LIQ
    rows = [f"{k} {tk:.17g} {o:.17g} 0 0 0 0 0 0 0 0 0"
            for k, (tk, o) in enumerate(zip(t, oxy))]
    (d / "tr_oxy.dat").write_text(
        "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 c1_liq_sum "
        "c1_liq_sum2 c2_liq_sum c2_liq_sum2 c3_liq_sum c3_liq_sum2\n"
        + "\n".join(rows) + "\n")
    (d / "vol_frac_interf.dat").write_text(
        "i t f_liq_sum f_liq_interf posY_max posY_min\n"
        + "\n".join(f"{k} {tk:.17g} {V_LIQ:.17g} 0 0 0" for k, tk in enumerate(t))
        + "\n")
    return d


def test_chain_walks_back_to_the_root(tmp_path):
    t = np.linspace(0, 1, 5)
    _seg(tmp_path, "s1", t, 0.1 * np.ones(5))
    _seg(tmp_path, "s2", t + 1, 0.2 * np.ones(5), parent="s1")
    s3 = _seg(tmp_path, "s3", t + 2, 0.3 * np.ones(5), parent="s2")
    assert [d.name for d in _segment_chain(s3)] == ["s1", "s2", "s3"]


def test_unchained_run_is_its_own_chain(tmp_path):
    solo = _seg(tmp_path, "solo", np.linspace(0, 1, 5), 0.1 * np.ones(5))
    assert [d.name for d in _segment_chain(solo)] == ["solo"]


def test_c_star_is_joined_across_segments(tmp_path):
    """A later segment alone starts mid-curve; joined, it starts at 0."""
    t1 = np.linspace(0.0, 5.0, 60)
    t2 = np.linspace(5.0, 10.0, 60)
    c1 = 1.0 - np.exp(-t1 / 2.0)          # rises from 0
    c2 = 1.0 - np.exp(-t2 / 2.0)          # already well up the curve
    _seg(tmp_path, "a", t1, c1)
    seg_b = _seg(tmp_path, "b", t2, c2, parent="a")

    t_joined, c_joined = _compute_c_star(seg_b)
    assert c_joined[0] == pytest.approx(c1[0], abs=1e-6), (
        "stitched C* must begin at the ROOT segment's first sample; starting "
        "mid-curve is exactly the failure that makes a chained kLa meaningless")
    assert t_joined[-1] == pytest.approx(t2[-1], rel=1e-9)
    assert np.all(np.diff(t_joined) > 0), "joined series must be time-ordered"
    # the duplicated boundary sample appears once
    assert len(t_joined) == len(t1) + len(t2) - 1
