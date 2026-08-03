"""Verification: VOF liquid volume (f_liq_sum) must not drift over a run.

Physical basis: navier-stokes/conserving.h implements a conservative VOF scheme.
f_liq_sum = statsf2(f).sum = integral of liquid volume fraction (true liquid volume).
At fidelity=3 (8×8 cells) with an embedded boundary, O(0.2%) variation is
expected from VOF reconstruction at cut cells. A drift > 0.5% signals a broken
VOF reconstruction or an embed mask leak.

Threshold calibrated from L6 reference run (health_l6_video): 0.04% measured
(recalibrated 2026-08-03 after the geometry.b fix, see below).

fidelity=4, not CANONICAL_PARAMS' default fidelity=3 (2026-08-03, diary.md):
geometry.b's fix (0.071->0.03575, see diary.md) halved the bag's fraction of
the L0=1 domain box, which at any FIXED fidelity halves the number of cells
resolving the bag's height. At fidelity=3 that dropped from ~4.5 to ~2.3
cells across the whole bag -- too coarse to keep VOF drift under 0.5%
(measured 2.6% post-fix at fidelity=3 in CI). fidelity=4 restores ~4.6 cells,
matching the effective resolution this test was originally calibrated
against; measured drift there is back down to ~0.09%.
"""
import pytest
from tests.conftest import CANONICAL_PARAMS, run_bioreactor, load_vol_frac


@pytest.mark.medium
def test_liquid_volume_conserved(tmp_path):
    params = {**CANONICAL_PARAMS, "run_id": "mass_cons", "fidelity": 4}
    run_dir = run_bioreactor(params, tmp_path)

    data = load_vol_frac(run_dir)
    assert len(data) >= 5, (
        f"vol_frac_interf.dat has only {len(data)} rows — sim may not have produced output"
    )

    f_liq_sum = data[:, 2]     # column index 2: f_liq_sum
    mean_f    = f_liq_sum.mean()
    drift     = (f_liq_sum.max() - f_liq_sum.min()) / mean_f

    assert drift < 5e-3, (
        f"VOF mass drift {drift:.2e} exceeds 0.5% "
        f"(min={f_liq_sum.min():.6g}, max={f_liq_sum.max():.6g}, mean={mean_f:.6g})"
    )
