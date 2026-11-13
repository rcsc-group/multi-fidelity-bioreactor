"""The cases of the feature bundle: short L6 runs that the guard tests in tests/bundle/ read.

Each case may restart from a file of its parent (a checkpoint, or the record's start dump), so one run's
output serves several tests and nothing is computed twice. Design: docs/specs/2026-11-13-feature-guardrails-design.md.

Time: `t_end` of a fresh run is absolute; on a restart it is the time run AFTER the checkpoint (as in the
solver and tools/golden/golden.py). T is one rocking period in the solver's non-dimensional time.
"""
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.postprocess import _t_scales  # noqa: E402  (the solver's own period, as golden.py uses it)

RPM_A, RPM_B = 32.5, 25.0
ZERO3 = [0.0, 0.0, 0.0]

# golden's five builds, plus the cross-level restart and the MPI lean build (made by a separate MPI job)
BUILDS = {"default": [], "diag": ["-DDIAGNOSTICS=1"], "lean": ["-DEXTRA_TRACERS=0"],
          "record": ["-DRECORD_FLUX=1"], "replay": ["-DREPLAY_FLUX=1"],
          "crosslevel": ["-DCROSS_LEVEL_WARMSTART=1"], "mpi": ["-D_MPI=1", "-DEXTRA_TRACERS=0"]}


@dataclass(frozen=True)
class Case:
    name: str
    build: str
    parent: str | None          # case whose output this one needs
    restart: str | None         # file of the parent to restart from ("checkpoint.dump", "replay_start.dump")
    level: int
    ranks: int
    env: dict = field(default_factory=dict)   # extra environment (e.g. a wall-clock limit)


CASES = {c.name: c for c in [
    Case("cold", "record", None, None, 6, 1),            # reference; also records the flux for 1 period
    Case("cold2", "record", None, None, 6, 1),           # same again: the noise floor
    Case("half", "record", None, None, 6, 1),            # stops at 2T and dumps; history must equal cold's
    Case("cont", "default", "half", "checkpoint.dump", 6, 1),
    # wall-clock limit set so that the early dump fires near 1.5 T; value from cold's wall time (Task 2)
    Case("kill", "default", None, None, 6, 1, {"BIOREACTOR_MAXRUNTIME_S": "0"}),
    Case("kill_resume", "default", "kill", "checkpoint.dump", 6, 1),
    Case("replay", "replay", "cold", "replay_start.dump", 6, 1),
    Case("loop", "replay", "cold", "replay_start.dump", 6, 1),
    Case("lean", "lean", None, None, 6, 1),
    Case("l5", "default", None, None, 5, 1),
    Case("xl6", "crosslevel", "l5", "checkpoint.dump", 6, 1),
    Case("caseB_cold", "default", None, None, 6, 1),
    Case("caseB_warm", "default", "half", "checkpoint.dump", 6, 1),
    Case("mpi4", "mpi", None, None, 6, 4),
]}


def order(cases: dict) -> list:
    """Case names with every parent before its children. ValueError on an unknown parent or a cycle."""
    for c in cases.values():
        if c.parent is not None and c.parent not in cases:
            raise ValueError(f"case {c.name}: unknown parent {c.parent}")
    done, seq = set(), []
    while len(seq) < len(cases):
        ready = [n for n, c in cases.items() if n not in done and (c.parent is None or c.parent in done)]
        if not ready:
            raise ValueError(f"cycle among cases {sorted(set(cases) - done)}")
        for n in ready:
            done.add(n)
            seq.append(n)
    return seq


def omega(rpm: float) -> float:
    return rpm * 2 * math.pi / 60


def base(run_id: str, rpm: float = RPM_A, level: int = 6) -> dict:
    """tools/golden/golden.py base_params, with the rate and level as arguments."""
    return {"run_id": run_id, "fidelity": level, "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
            "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": ZERO3, "omega_b": omega(rpm),
            "omega_h": 0.0, "amplitude_h": ZERO3, "phi_horizontal": ZERO3, "n_mix_cycles": 1}


def period(rpm: float = RPM_A) -> float:
    return _t_scales(base("x", rpm))[1]


def _prev(rpm: float) -> dict:
    """Forcing of the run a restart comes from (chain.py keys)."""
    return {"omega_b_prev": omega(rpm), "theta_max_prev": [7.0, 0.0, 0.0], "phi_angular_prev": ZERO3,
            "amplitude_h_prev": ZERO3, "phi_horizontal_prev": ZERO3, "omega_h_prev": 0.0}


def params(name: str, runs: Path) -> dict:
    """params.json of a case; `runs` is the folder that holds one sub-folder per case."""
    T = period(RPM_A)
    rec = {"n_rec_periods": 1, "rec_dir": "."}
    if name in ("cold", "cold2"):
        return base(name) | {"t_end": 3 * T} | rec
    if name == "half":
        return base(name) | {"t_end": 2 * T} | rec
    if name == "cont":
        return base(name) | {"t_checkpoint": 2 * T, "restart_continue": 1, "n_mix_cycles": 0, "t_end": T}
    if name == "kill":
        return base(name) | {"t_end": 3 * T}
    if name == "kill_resume":   # estimate; the runner sets t_checkpoint = the dump's time, t_end = 3T - that
        return base(name) | {"t_checkpoint": 1.5 * T, "restart_continue": 1, "n_mix_cycles": 0, "t_end": 1.5 * T}
    if name == "replay":        # golden's replay keys: the record covers [T, 2T]
        return base(name) | {"t_checkpoint": T, "restart_continue": 1, "n_mix_cycles": 0, "t_end": T,
                             "rec_dir": str(runs / "cold")}
    if name == "loop":          # the one recorded period, looped for 2 periods: [T, 3T]
        return base(name) | {"t_checkpoint": T, "restart_continue": 1, "n_mix_cycles": 0, "t_end": 2 * T,
                             "n_rec_periods": 1, "replay_loop": 1, "rec_dir": str(runs / "cold")}
    if name == "lean":
        return base(name) | {"t_end": 3 * T}
    if name == "l5":
        return base(name, level=5) | {"t_end": 2 * T}
    if name == "xl6":           # as tests/verification/test_cross_level_warmstart.py
        return base(name) | {"t_checkpoint": 2 * T, "t_end": T} | _prev(RPM_A)
    if name == "caseB_cold":
        return base(name, RPM_B) | {"t_end": 3 * period(RPM_B)}
    if name == "caseB_warm":    # condition A's state at 2T -> condition B, 1 transition cycle, then release
        return base(name, RPM_B) | {"t_checkpoint": 2 * T, "t_end": 3 * period(RPM_B),
                                    "n_transition_cycles": 1} | _prev(RPM_A)
    if name == "mpi4":
        return base(name) | {"t_end": 3 * T}
    raise KeyError(name)
