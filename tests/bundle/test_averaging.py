"""Guard S11 on a real run: averages are taken over the liquid INSIDE the bag.

The 2026-09-15 defect averaged over every cell with f > 0.5, which is half the domain (0.50) instead of the
liquid in the bag (0.1425 at L10; tests/verification/test_shear_averaging_region.py). The only area the
solver prints is f_liq_sum (vol_frac_interf.dat), so this checks that one; the tau/EDR mask itself keeps its
static guard in test_shear_averaging_region.py.
"""
import numpy as np

from tests.bundle import _series as s
from tests.bundle.conftest import case

BAG_LIQUID = 0.1425   # liquid area / domain area at fill 0.5 (test_shear_averaging_region.py docstring)


def test_s11_liquid_area_is_the_bag_liquid_not_half_the_domain(bundle):
    _, area = s.column(case(bundle, "cold")["dir"], "vol_frac_interf.dat", "f_liq_sum")
    assert np.all(np.abs(area - BAG_LIQUID) <= 0.05 * BAG_LIQUID), (
        f"f_liq_sum in [{area.min():.4f}, {area.max():.4f}], expected {BAG_LIQUID} +- 5 % "
        f"(0.50 means the outside of the bag is counted)")
