"""Where a run stops: the period-boundary rule shared by the solver and Python.

A run always ends (and dumps its checkpoint) on a rocking-period boundary, at
theta = 0. The rule is the SMALLEST boundary at or after t, with t/T snapped to
an integer when it is within PERIOD_TOL of one:

    n = max(1, ceil(t/T - PERIOD_TOL)),   t_stop = n*T

The old rule, floor(t/T) + 1, gave k or k+1 periods for a request of exactly k,
depending on the last bit of t/T (diary 2026-09-30). The C code in
src/BioReactor.c (fresh and restarted runs) implements the same expression.
Keep the two in step: chain.py and sweep.py use this to predict the
checkpoint time the solver will write.

PERIOD_TOL = 1e-5 period is far above the round-off of t/T, and above the
~1e-6 period error of round(t_end, 6) in the submit scripts. It is far below
any duration anyone would request.
"""
from __future__ import annotations

import math

PERIOD_TOL = 1e-5


def next_period_boundary(t: float, T: float) -> tuple[int, float]:
    """(n, n*T) for the first period boundary at or after t."""
    n = max(1, math.ceil(t / T - PERIOD_TOL))
    return n, n * T
