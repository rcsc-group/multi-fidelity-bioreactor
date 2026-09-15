"""Regression guard: the video-frame cadence must NOT divide the rocking period.

Why this exists (diary.md 2026-09-15 (3))
-----------------------------------------
`dt_video` used to be `T_per_st / N`, an exact divisor of the rocking period.
Every frame then lands on the SAME N phases of the cycle forever, so recording
more frames adds no phase coverage at all and the true peak of any oscillating
quantity is never sampled -- only a max-over-N-samples lower bound.

Measured on real runs before the fix: a34fc4d4 (L9, 130 frames) and 57f68830
(L10, 50 frames) each sampled only FIVE distinct phases (0.0, 0.2, 0.4, 0.6,
0.8). Every <tau'_w> and <eps'_w> peak this project reported from frames was
therefore understated. l8_coldstart_vid escaped only by accident, because its
timestep drifted relative to its period -- which is why L8 was the one level
that showed a real waveform while the better-resolved L9/L10 showed 5 dots.

Kim et al. designed around exactly this (Main.tex): "a constant time gap of
0.05 simulation time ... approximately 13 simulation points per cycle but is
intentionally misaligned with the period to ensure convergence."

The fix offsets the divisor by the golden-ratio conjugate (0.6180339887...),
the irrational hardest to approximate by rationals, giving low-discrepancy
(near-uniform) phase coverage as frames accumulate.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

SRC = Path(__file__).parents[2] / "src" / "BioReactor.c"
PARAMS_H = Path(__file__).parents[2] / "src" / "params_read.h"


def test_dt_video_is_not_an_exact_divisor_of_the_period():
    text = SRC.read_text()
    m = re.search(r"dt_video\s*=\s*T_per_st\s*/\s*\(([^;]+)\);", text)
    assert m, "could not find the dt_video assignment in BioReactor.c"
    expr = m.group(1)
    assert "0.618" in expr, (
        "dt_video divides the rocking period exactly:\n"
        f"    dt_video = T_per_st / ({expr});\n\n"
        "That locks every frame onto the same N phases forever, so the peak of "
        "any oscillating quantity is never sampled. Offset the divisor by the "
        "golden-ratio conjugate (0.6180339887498949). See diary.md 2026-09-15 (3)."
    )


def test_offset_cadence_actually_scans_the_phase():
    """The property that matters, checked numerically rather than textually."""
    Tp = 0.607329
    n_periods, N = 60, 13
    for off, min_distinct in ((0.0, None), (0.6180339887498949, 500)):
        dt = Tp / (N + off)
        ph = (np.arange(0, n_periods * Tp, dt) / Tp) % 1.0
        distinct = len(np.unique(np.round(ph, 4)))
        if off == 0.0:
            assert distinct <= N + 1, (
                f"sanity check on the test: an exact divisor should give <= {N} "
                f"phases, got {distinct}"
            )
        else:
            assert distinct >= min_distinct, (
                f"offset cadence only reached {distinct} distinct phases over "
                f"{n_periods} periods -- it is not scanning the cycle."
            )


def test_default_cadence_matches_kim():
    assert re.search(r"p\.frames_per_period\s*=\s*13\b", PARAMS_H.read_text()), (
        "default frames_per_period should be 13 (~Kim's stated cadence of "
        "approximately 13 points per cycle). See diary.md 2026-09-15 (3)."
    )
