"""Two runs that together decide whether the escape is restart-specific and
whether it terminates.

Context (diary.md 2026-09-07). The 2x2 showed the escape tracks the restored
STATE, not the ramp code path: ramp live + bit-exact state -> 0.999, ramp off
+ state 0.3% off target -> 1.246. Meanwhile a 1e-4 UNIFORM velocity rescale
on a pure cold start decayed. A uniform rescale is a poor probe (div-free,
symmetry-preserving, nearly tangent to the periodic-orbit family), so that
null result does not establish stability against a generic perturbation.

  (1) kicknoise -- cold start, no checkpoint anywhere, with a one-shot
      per-cell random u *= (1 + 3e-3*noise()) at cycle 40, matching the
      ~0.3% magnitude of the state mismatch that does escape. 200 cycles.
        escapes  -> the reference state is linearly stable but has a finite
                    basin, 0.3% is outside it, and NOTHING here is
                    restart-specific: the restart is just a convenient way
                    to deliver a finite perturbation.
        decays   -> a 0.3% generic perturbation is inside the basin, so the
                    restart path is still injecting something extra, and the
                    remaining defect is on that path.

  (2) longrun -- case C (ramp off, state off target, the cleanest escaping
      configuration) rerun to 400 cycles. The escaping runs were all still
      climbing at cycle 120 (tau_100 0.0075 -> 0.0211 over the last four
      20-cycle blocks), so "a second attractor" is not established.
        saturates -> a genuine second periodic state; report its separation.
        grows without bound -> a numerical instability, and the L7 grid's
                    ability to hold the reference state is the real issue.

Usage:
    uv run python scripts/submit_kick_noise_and_longrun.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.simulate import submit_slurm
from scripts.settling_study_grid import (GEOMETRY, FILL_LEVEL, PROJECT_ROOT,
                                         TEMPLATE, omega_b_of)

chain.validate_params = lambda params: None

T_PER_ND = 0.607329
RPM, THETA, FIDELITY = 32.5, 7.0, 7
ROOT = Path(PROJECT_ROOT)
runs_dir = ROOT / "runs"

# --- (1) generic finite-amplitude kick on a pure cold start ----------------
N = 200
kick = {
    "run_id": "kicknoise_L7_th7", "fidelity": FIDELITY,
    "geometry": GEOMETRY, "fill_level": FILL_LEVEL, "n_harmonics": 1,
    "theta_max": [THETA, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega_b_of(RPM), "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": N * T_PER_ND, "n_mix_cycles": N,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-kicknoise",
}
job = submit_slurm(kick, project_root=ROOT, walltime="04:00:00",
                   template=Path(TEMPLATE), cpus=1, ntasks=8, mem="4G",
                   exclude="node2336")
print(f"kicknoise (3e-3 random kick at cycle 40, {N} cycles): job={job}")

# --- (2) does the escape saturate? ----------------------------------------
src_dir = runs_dir / "dtheta_src_L7_th6.9"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
cfg = {
    "motion": {"omega_b": omega_b_of(RPM), "theta_max": [THETA, 0.0, 0.0]},
    "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": 400, "n_transition_cycles": 400, "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM)]},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(RPM),
        "theta_max": [THETA, 0.0, 0.0],          # ramp off, state off target
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "08:00:00",
    "binary": "/oscar/scratch/eaguerov/BioReactor-mpi-testc",
    "submit": True, "exclude": "node2336",
}
print(f"longrun (case C, 400 cycles): {chain.submit_chain(cfg)}")
