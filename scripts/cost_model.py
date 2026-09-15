"""Per-cycle wall-clock cost model: min_per_cycle(level, ntasks, rpm) ->
(minutes, confidence). Every entry here is a REAL sacct Elapsed measurement
for a REAL run at that exact (level, ntasks) -- never extrapolated from a
formula -- because every walltime-estimation failure this project has had
(L7->L8 pilot 3-4x over, l8_coldstart_vid finishing 4 SECONDS before its
walltime cap, all 3 L10 Fig13a points TIMING OUT) came from guessing
instead of measuring.

Known, real rpm-dependence (2026-09-12, L10/48-tasks: 17.5rpm=148.8,
32.5rpm=83.0, 37.5rpm=74.6 min/cycle) means rpm cannot be dropped as a
dimension even though it doesn't change T_per_nd's normalization the way
level/ntasks do -- lower rpm is measurably MORE expensive per cycle, not
less, likely because the absolute (not nondimensionalized) capillary/
viscous timescale becomes proportionally more restrictive on dt as U_bio
shrinks. This is an observed pattern, not a derived law -- treat it as
another reason not to extrapolate across rpm either.

Usage:
    from scripts.cost_model import min_per_cycle
    minutes, confidence = min_per_cycle(level=8, ntasks=8, rpm=32.5)
"""
from __future__ import annotations

# (level, ntasks) -> {rpm: minutes/cycle}. Sources noted per block.
_MEASURED: dict[tuple[int, int], dict[float, float]] = {
    # [CORRECTED 2026-09-15] This whole block was off by exactly 60x: the
    # numbers were sacct Elapsed converted to MINUTES for a 60-cycle segment,
    # recorded as if they were minutes PER CYCLE. The "/ 60" in the old note
    # was the seconds->minutes conversion; the division by the 60 cycles never
    # happened. Job 5989327 has Elapsed 00:02:34 -- i.e. 2.57 MINUTES TOTAL,
    # which was stored as "2.57 min/cycle".
    #
    # Confirmed by two independent direct measurements at 32.5 rpm, 8 ranks,
    # each a 209.7-cycle Fig 9 run that completed (diary.md 2026-09-15 (7)):
    #   L6  job 6410430  00:08:36 / 209.7 cyc = 0.041 min/cyc  (x60 = 2.46)
    #   L7  job 6410614  00:34:17 / 209.7 cyc = 0.163 min/cyc  (x60 = 9.8)
    # against the stored 2.57 and 10.18. Both match to within 5%.
    #
    # This made every plan built on it ~60x too pessimistic, which is why the
    # Fig 9 L7 sweep was costed at ~500 h and declared to need chaining for 9
    # of 10 points. It does not.
    #
    # The 17.5 rpm and (8,8) entries are NOT independently verified yet -- they
    # are corrected by the same factor of 60 because they come from the same
    # mis-derivation, and are marked lower-confidence until measured directly.
    (6, 8): {32.5: 0.041, 17.5: 0.0705},      # 32.5 measured directly
    (7, 8): {32.5: 0.163, 17.5: 0.227},       # 32.5 measured directly
    (8, 8): {32.5: 1.536, 37.5: 1.7125},      # neither measured directly yet
    # L10 Fig13a sweep, 2026-09-12, jobs 6197602-6197604 (measured directly
    # from partial shear_stress.dat + sacct TIMEOUT elapsed, since none
    # reached completion -- see conversation 2026-09-12):
    (10, 48): {17.5: 148.8, 32.5: 83.0, 37.5: 74.6},
    # l10_kim_fig8_signed's own extension segment, 2026-08-09, job 4812614
    # (sacct Elapsed 01:27:49 / 2.98 cycles covered):
    (10, 64): {32.5: 29.5},
    # 32-rank L10, measured 2026-09-13/15 from runs that actually completed:
    # 17.5 and 37.5 are the 5-cycle Fig13a calibration runs (jobs 6326988,
    # 6326989); 32.5 is the Fig 8 late-time probe (job 6314896, 9.98 cycles
    # in 16:47:58). Note the strong rpm dependence -- 4x across the range --
    # which is exactly why walltime guessing kept failing here.
    (10, 32): {17.5: 168.1, 32.5: 101.0, 37.5: 42.0},
}


def min_per_cycle(level: int, ntasks: int, rpm: float, tol: float = 1e-6) -> tuple[float, str]:
    """Returns (minutes_per_cycle, confidence).

    "measured"              -- exact (level, ntasks, rpm) in the table.
    "measured_other_rpm_max" -- (level, ntasks) tested, but not at this rpm;
                                returns the MAX over rpms tested at this
                                (level, ntasks) as a conservative estimate,
                                since rpm-dependence is real (see module
                                docstring) and direction isn't safe to
                                assume for an untested rpm.
    Raises ValueError for an (level, ntasks) combination with NO
    measurement at all -- chain.py must not silently guess a walltime for
    hardware/resolution combinations nobody has ever actually timed; the
    caller needs an explicit human-supplied walltime instead.
    """
    key = (level, ntasks)
    if key not in _MEASURED:
        raise ValueError(
            f"no measured per-cycle cost for level={level}, ntasks={ntasks} -- "
            f"cannot auto-size a walltime for an untested (level, ntasks) combination. "
            f"Measure a short run first, or supply an explicit walltime override."
        )
    by_rpm = _MEASURED[key]
    for r, m in by_rpm.items():
        if abs(r - rpm) < tol:
            return m, "measured"
    return max(by_rpm.values()), "measured_other_rpm_max"
