"""The solver stops on the period boundary scripts/periods.py predicts.

A request for exactly k periods, rounded to 6 decimals as every submit script
does, must stop at k periods, not k+1. The old rule floor(t/T)+1 gave 11
periods for 10 and 20 for 20 on the same binary, which put the two sides of
test_restart_continue 2 periods apart (diary 2026-09-30). A request between
boundaries must round UP to the next one.
"""
from __future__ import annotations

import math
import sys
import pathlib

import pytest

PROJECT_ROOT = pathlib.Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from tests.conftest import run_bioreactor                       # noqa: E402
from tests.verification.test_restart_continue import _params, T_PER_ND  # noqa: E402
from scripts.dump_fields import fields as dump_fields           # noqa: E402
from scripts.periods import next_period_boundary                # noqa: E402


@pytest.mark.medium
@pytest.mark.parametrize("cycles", [10, 10.4])
def test_fresh_run_stops_on_the_predicted_boundary(tmp_path, cycles):
    d = run_bioreactor(_params(f"pb_{cycles:g}", cycles), tmp_path, timeout=900,
                       require_complete=True)
    t_dump = dump_fields(d / "checkpoint.dump")[0]["t"]
    n_expect, t_expect = next_period_boundary(round(cycles * T_PER_ND, 6), T_PER_ND)
    assert n_expect == math.ceil(cycles)
    assert t_dump == pytest.approx(t_expect, abs=1e-6 * T_PER_ND), (
        f"requested {cycles} periods; solver dumped at {t_dump / T_PER_ND:.6f} "
        f"periods, scripts/periods.py predicts {n_expect}")
