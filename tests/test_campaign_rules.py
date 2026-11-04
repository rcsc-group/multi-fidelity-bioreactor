"""Campaign 2D (2026-11-04): kLa_1T_50 reported-only rule, replay saving under that rule, and the `next` selection.

Pure functions only (no data, no fit). Each expected value is stated in the test.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))

import campaign as c  # noqa: E402

CRIT = ["tau_95", "edr_95", "kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50"]


def test_kla50_is_reported_not_coupled_below_35_rpm():
    assert c.scalar_source(8, 20.0, "kLa_1T_50") == "reported"
    assert c.scalar_source(9, 34.0, "kLa_1T_50") == "reported"


def test_kla50_does_not_force_a_coupled_tail():
    pure_scalars = ["kLa_1T_50", "kLa_1T_10", "kLa_1T_25", "dtmix_0.50"]
    assert not c.point_needs_coupled(8, 20.0, pure_scalars)
    assert not c.point_needs_coupled(9, 30.0, pure_scalars)


def test_tau_and_edr_still_need_the_coupled_tail():
    assert c.point_needs_coupled(8, 20.0, ["tau_95"])
    assert c.point_needs_coupled(8, 20.0, ["edr_95"])


def test_replay_scalars_below_35_rpm():
    assert c.scalar_source(8, 20.0, "kLa_1T_10") == "replay"
    assert c.scalar_source(8, 35.0, "kLa_1T_10") == "coupled"   # rpm >= 35 keeps the coupled tail


def test_saving_is_zero_under_the_decided_rule():
    pts = [(20.0, 3.0), (25.0, 5.0), (30.0, 6.0)]
    assert c.replay_saving(pts, level=8) == 0.0   # tau_95 and edr_95 keep the coupled tail


def test_saving_if_tau_and_edr_were_replay_sourced():
    # Alternative (not decided): only the scalars plus kLa_1T_50 reported. Points below 35 rpm then save
    # (L8 cold minus L8 replay), points at or above 35 rpm do not.
    scal = [q for q in CRIT if q not in c.COUPLED_ONLY_QOIS]
    pts = [(20.0, 3.0), (25.0, 5.0), (36.0, 6.0)]
    per = c.COST["L8_cold"][0] - c.COST["L8_replay"][0]
    assert c.replay_saving(pts, level=8, qois=scal) == 2 * per


def test_select_next_greedy_with_separation():
    sites = [(20.0, 3.0), (20.5, 3.05), (30.0, 6.0), (16.0, 2.2)]
    scores = {sites[0]: 0.9, sites[1]: 0.85, sites[2]: 0.5, sites[3]: 0.1}
    out = c.select_next(sites, scores, k=2)
    # the second site is a near-duplicate of the first, so the third-best distinct site is taken
    assert [p[:2] for p in out] == [(20.0, 3.0), (30.0, 6.0)]


def test_select_next_returns_at_most_k_and_sorted():
    sites = [(15.0, 2.0), (25.0, 4.0), (35.0, 6.0)]
    scores = {sites[0]: 0.2, sites[1]: 0.7, sites[2]: 0.4}
    out = c.select_next(sites, scores, k=2)
    assert len(out) == 2
    assert out[0][2] >= out[1][2] == 0.4
