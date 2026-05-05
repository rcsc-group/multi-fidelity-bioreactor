"""Verification: dissolved O2 in liquid (oxy_liq_sum) must trend upward over the run.

Physical basis: Henry's law drives net O2 flux from air (high) to liquid (low).
OXYGEN_AIR=1 initialises oxy in the air phase only; no degassing path exists.

Why not strict step-by-step monotonicity: oxy_liq_sum = Σ oxy·f·dV is a
VOF-weighted integral. As the free surface sloshes, interface cells cycle between
f≈1 (liquid) and f≈0 (air), causing the sum to oscillate at the rocking
frequency even while the *net* dissolved O2 increases. Step-by-step decreases
are therefore expected and not a numerical pathology.

Correct invariant: the linear trend of oxy_liq_sum over the full simulation
must be strictly positive.  A zero or negative trend means Henry's law is not
transferring oxygen into the liquid — the most likely cause being a sign error
in the diffusion coefficient or an O2 leak through the embed boundary.
"""
import numpy as np
import pytest
from tests.conftest import CANONICAL_PARAMS, run_bioreactor, load_tr_oxy


@pytest.mark.medium
def test_oxy_liq_sum_trends_upward(tmp_path):
    params = {**CANONICAL_PARAMS, "run_id": "oxy_trend"}
    run_dir = run_bioreactor(params, tmp_path)

    data = load_tr_oxy(run_dir)
    assert len(data) >= 10, (
        f"tr_oxy.dat has only {len(data)} rows — need at least 10 for trend test"
    )

    t   = data[:, 1]    # time
    oxy = data[:, 2]    # oxy_liq_sum

    slope, _ = np.polyfit(t, oxy, 1)
    assert slope > 0, (
        f"oxy_liq_sum linear trend is {slope:.3e} (≤ 0): "
        "Henry's law is not transferring O2 into liquid over the run"
    )
