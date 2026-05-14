"""Verification: interface perimeter (f_liq_interf) must not grow monotonically.

Physical basis: for a rocking bioreactor without atomisation, the interface
oscillates quasi-periodically at the forcing frequency. Monotonic growth of
f_liq_interf signals spurious droplet/filament generation from numerical
fragmentation — a sign of under-resolution or a broken PLIC reconstruction.

Test: the time-mean of f_liq_interf in the second half of the run must not
exceed 1.5× the mean of the first half.

Threshold calibrated from L6 reference run (health_l6_video): ratio = 0.998.
1.5 gives a large margin while firmly catching fragmentation (old threshold 3.0
was explicitly labelled 'deliberately generous' in the prior version).
"""
import pytest
from tests.conftest import CANONICAL_PARAMS, run_bioreactor, load_vol_frac


@pytest.mark.medium
def test_interface_area_does_not_blow_up(tmp_path):
    params = {**CANONICAL_PARAMS, "run_id": "iface_bound"}
    run_dir = run_bioreactor(params, tmp_path)

    data = load_vol_frac(run_dir)
    assert len(data) >= 10, (
        f"vol_frac_interf.dat has only {len(data)} rows — sim may not have produced output"
    )

    iface = data[:, 3]          # column index 3: f_liq_interf
    mid   = len(iface) // 2
    mean_first  = iface[:mid].mean()
    mean_second = iface[mid:].mean()

    assert mean_first > 0, "First-half interface area is zero — unexpected geometry"
    ratio = mean_second / mean_first

    assert ratio < 1.5, (
        f"Interface area grew by factor {ratio:.2f} between first and second half "
        f"(first_mean={mean_first:.4g}, second_mean={mean_second:.4g})"
    )
