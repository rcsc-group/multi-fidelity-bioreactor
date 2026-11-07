"""Warm-start settled-state gate for the Campaign 2D driver (2026-11-04).

Synthetic series only. One settled run must PASS; one drifting run must FAIL-GATE.
Also covers PENDING (too few cycles) and the fs-sanity .err line.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))

from campaign import settled_gate  # noqa: E402

PERIOD = 1.0
WARM_END = 3.0            # warm-up = first 3 cycles (synthetic; plan uses 12)
FS_OK = "fs sanity [post-solid]: OK"


def _series(n_cycles, drift=0.0, noise=0.002, seed=0, per_cycle=40):
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, n_cycles * PERIOD, n_cycles * per_cycle, endpoint=False)
    k = np.floor(t / PERIOD)
    base_tau = 1.0 * (1.0 + drift) ** k
    base_edr = 0.2 * (1.0 + drift) ** k
    tau = base_tau * (1.0 + noise * rng.standard_normal(t.size))
    edr = base_edr * (1.0 + noise * rng.standard_normal(t.size))
    return t, tau, edr


def test_settled_run_passes():
    t, tau, edr = _series(n_cycles=13, drift=0.0)
    status, det = settled_gate(t, tau, edr, PERIOD, WARM_END, release_t=WARM_END,
                               height_range=0.301, src_height_range=0.300, err_text=FS_OK)
    assert status == "PASS", det


def test_drifting_run_fails_gate():
    t, tau, edr = _series(n_cycles=13, drift=0.03)   # 3 % per cycle, above the 2 % limit
    status, det = settled_gate(t, tau, edr, PERIOD, WARM_END, release_t=WARM_END,
                               height_range=0.300, src_height_range=0.300, err_text=FS_OK)
    assert status == "FAILED-GATE", det
    assert det["drift_ok"] is False


def test_release_height_mismatch_fails_gate():
    t, tau, edr = _series(n_cycles=13, drift=0.0)
    status, det = settled_gate(t, tau, edr, PERIOD, WARM_END, release_t=WARM_END,
                               height_range=0.310, src_height_range=0.300, err_text=FS_OK)   # 3.3 %
    assert status == "FAILED-GATE", det
    assert det["height_ok"] is False


def test_missing_fs_sanity_line_fails_gate():
    t, tau, edr = _series(n_cycles=13, drift=0.0)
    status, det = settled_gate(t, tau, edr, PERIOD, WARM_END, release_t=WARM_END,
                               height_range=0.300, src_height_range=0.300, err_text="no sanity line here")
    assert status == "FAILED-GATE", det
    assert det["fs_ok"] is False


def test_too_few_cycles_is_pending():
    t, tau, edr = _series(n_cycles=4, drift=0.0)     # 1 cycle after warm-up: need 3 blocks of 3
    status, _ = settled_gate(t, tau, edr, PERIOD, WARM_END, release_t=WARM_END,
                             height_range=0.300, src_height_range=0.300, err_text=FS_OK)
    assert status == "PENDING"


def test_period3_limit_cycle_passes_drift():
    """2026-11-07: 19.14/2.89, 32.19/4.51, 34.81/4.26 at L9 repeat every 3 cycles (cycle means ~1, 1, 0.96).
    That is a steady 3T limit cycle, not drift. Drift compares 3-cycle block means, so it must PASS."""
    n = 30
    t = np.linspace(0.0, n * PERIOD, n * 40, endpoint=False)
    k = np.floor(t / PERIOD).astype(int)
    amp = np.where(k % 3 == 2, 0.96, 1.0)
    st, det = settled_gate(t, amp, 0.2 * amp, PERIOD, WARM_END, None, None, None, FS_OK)
    assert det["drift_ok"] is True, det


def test_slow_drift_over_blocks_fails():
    t, tau, edr = _series(30, drift=0.01)          # 1 %/cycle -> ~3 % per 3-cycle block
    st, det = settled_gate(t, tau, edr, PERIOD, WARM_END, None, None, None, FS_OK)
    assert det["drift_ok"] is False, det
