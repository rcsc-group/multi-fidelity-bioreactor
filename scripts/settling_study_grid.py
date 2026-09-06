"""Shared grid definition for the settling-time study (diary.md 2026-09-06).

3x3 factorial in (omega_b, theta_max), baseline at (32.5rpm, 7deg), at
three fidelity levels. Used by submit_settling_study_baselines.py,
submit_settling_study_references.py, and submit_settling_study_transitions.py
so all three scripts share one source of truth for the grid and run_id
naming convention.
"""
import math

LEVELS = {6: 8, 7: 8, 8: 8}  # fidelity -> ntasks (uniform 8 ranks per user's request)

BASELINE = {"rpm": 32.5, "theta": 7.0}

OMEGA_RPMS = [17.5, 32.5, 37.5]
THETAS = [2.0, 4.0, 7.0]

GRID = [
    {"rpm": rpm, "theta": theta}
    for rpm in OMEGA_RPMS
    for theta in THETAS
]
OUT_POINTS = [p for p in GRID if not (p["rpm"] == BASELINE["rpm"] and p["theta"] == BASELINE["theta"])]
assert len(OUT_POINTS) == 8, f"expected 8 non-baseline points, got {len(OUT_POINTS)}"


def omega_b_of(rpm: float) -> float:
    return rpm * 2 * math.pi / 60.0


def point_tag(p: dict) -> str:
    return f"rpm{p['rpm']:g}_th{p['theta']:g}"


def baseline_run_id(fidelity: int) -> str:
    return f"settling_baseline_L{fidelity}"


def reference_run_id(fidelity: int, p: dict) -> str:
    return f"settling_ref_L{fidelity}_{point_tag(p)}"


T_PER_ND = 0.607329  # nondim period per cycle -- independent of omega_b (fixed geometry/theta)
CYCLES_PER_RUN = 25  # uniform for baseline/reference/transition -- mainline's ramp is a
                     # FIXED 3 cycles (N_RAMP_CYCLES) regardless of RPM/theta, unlike the
                     # rampmatched binary's RPM-dependent upstream ramp, so 25 gives ~22
                     # genuine post-ramp cycles either way -- comfortable margin given the
                     # observed settling range (0-9 cycles) for both baseline/reference cold
                     # starts and post-checkpoint transitions.
T_END = CYCLES_PER_RUN * T_PER_ND

PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
TEMPLATE = f"{PROJECT_ROOT}/config/slurm_mpi_template.sh"
MAINLINE_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi"  # not rampmatched -- see 2026-09-05 (2)
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
FILL_LEVEL = 0.5
