"""Regression guard: spatial averages of shear stress and EDR must be taken
over the LIQUID INSIDE THE BAG, not over every cell with f > 0.5.

Why this test exists
--------------------
2026-09-15 (diary.md 2026-09-15 (2)). `f` is initialised with
`fraction(f, y_init - y)`, which fills the lower half of the WHOLE DOMAIN --
including everything outside the embedded bag, where there is no fluid.
The tau/EDR diagnostic masked on `f[] > 0.5` alone, with no `cs[]` check, so
those cells entered the spatial averages: they contribute (near-)nothing to
the numerator but their full area to the denominator.

Measured at L10 (runs/57f68830, 1024^2): cells with f>0.5 covered 0.5000 of
the domain while liquid actually inside the bag covered 0.1425 -- a 3.51x
inflation. Every reported tau_mean / ediss_mean was consequently ~3x too
small, which was the entire standing "~3x below Kim et al." discrepancy:

    quantity (32.5rpm, th=7)   as reported   bag-masked      Kim
    tau_liq_mean   [Pa]          5.39e-4      1.685e-3     1.803e-3
    Ediss_liq_mean [W/m3]        8.69e-2      0.2773       0.2866

0.30x -> 0.93x and 0.30x -> 0.97x of Kim. The physics was never wrong; the
diagnostic was.

Kim's subscript `w` means WATER, not wall -- Main.tex, sec. "Shear stress and
energy dissipation rate": "The shear stress in the water is defined as
tau_w' = mu_w(du_x/dy + du_y/dx)", "where mu_w is the viscosity of water",
and Fig. 8's caption says the histograms are "in water". So the correct
averaging region genuinely is the water; the only error was including
domain cells outside the bag.

Why it hid: spatial MAXIMA are unaffected (u=0 outside the bag => tau=0
there, which cannot change a max), so peak tau always looked sane. And being
an area/mask error it is resolution-independent, so it survived every
grid-convergence check -- the signed mean was flat to <5% from L6 to L10.

This is a static guard rather than a CFD run: the defect is a one-token
omission in a mask, it is cheap and deterministic to assert directly, and it
therefore runs in the default fast suite rather than behind `-m medium`.
"""
from __future__ import annotations

import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[2]
SRC = PROJECT_ROOT / "src" / "BioReactor.c"

# Real `if (...)` masks that select liquid cells. Comments mentioning the
# mask are deliberately excluded -- only executable conditions count.
# Currently: the kLa/interface loop (which already required cs[]==1), the
# tau/EDR averaging pass, the histogram pass, and the video field pass.
_EXPECTED_MASKED_LOOPS = 4


def _mask_lines() -> list[str]:
    text = SRC.read_text()
    out = []
    for ln in text.splitlines():
        st = ln.strip()
        if st.startswith("//"):
            continue
        if "f[] > 0.5" in st and st.startswith("if ("):
            out.append(st)
    return out


def test_every_liquid_mask_also_checks_cs():
    """No tau/EDR mask may select on `f[] > 0.5` without also requiring cs[] > 0."""
    masks = _mask_lines()
    assert masks, f"no `f[] > 0.5` masks found in {SRC} -- has the diagnostic moved?"
    assert len(masks) == _EXPECTED_MASKED_LOOPS, (
        f"expected {_EXPECTED_MASKED_LOOPS} liquid-masked loops, found {len(masks)}:\n"
        + "\n".join(masks)
        + "\nIf a loop was added or removed, update this guard deliberately."
    )
    bad = [m for m in masks if "cs[]" not in m]
    assert not bad, (
        "A tau/EDR loop masks on `f[] > 0.5` without requiring `cs[] > 0.`:\n"
        + "\n".join(bad)
        + "\n\n`f` is filled over the whole lower half-DOMAIN, so this sweeps in "
        "cells outside the embedded bag and dilutes every spatial mean "
        "(measured 3.51x at L10). Use `cs[] > 0. && f[] > 0.5`. "
        "See diary.md 2026-09-15 (2)."
    )


def test_area_accumulators_are_weighted_by_cs():
    """Areas must be cs-weighted: a cut cell is only partly fluid."""
    text = SRC.read_text()
    # the averaging pass defines dA once and uses it for every accumulation
    assert re.search(r"double\s+dA\s*=\s*cs\[\]\s*\*\s*\(\s*Delta\s*\*\s*Delta\s*\)", text), (
        "expected `double dA = cs[]*(Delta*Delta);` in the tau/EDR averaging pass -- "
        "cut cells are only partially fluid, so their area must be weighted by cs[]. "
        "See diary.md 2026-09-15 (2)."
    )
    offenders = [
        ln.strip()
        for ln in text.splitlines()
        if re.search(r"(tau_sum|tau_vol|ediss_sum|tau_sum_signed|tau_sum_strict|tau_vol_strict)\s*\+=", ln)
        and "dA" not in ln
    ]
    assert not offenders, (
        "tau/EDR area accumulation not using the cs-weighted area `dA`:\n"
        + "\n".join(offenders)
        + "\nSee diary.md 2026-09-15 (2)."
    )
