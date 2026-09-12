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
    # settling-study grid extension runs, 2026-09-12, jobs 5989327-5989350
    # and 6287058-6287064 (60-cycle segments, sacct Elapsed / 60):
    (6, 8): {32.5: 2.57, 17.5: 4.23},
    (7, 8): {32.5: 10.18, 17.5: 13.63},
    (8, 8): {32.5: 92.15, 37.5: 102.75},  # 37.5 averaged over 2 runs (91.62, 113.9)
    # L10 Fig13a sweep, 2026-09-12, jobs 6197602-6197604 (measured directly
    # from partial shear_stress.dat + sacct TIMEOUT elapsed, since none
    # reached completion -- see conversation 2026-09-12):
    (10, 48): {17.5: 148.8, 32.5: 83.0, 37.5: 74.6},
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
