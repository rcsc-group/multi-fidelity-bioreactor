"""Re-verification of the phase fix after the offset formula changed.

v1 derived t_phase_offset from the U_bio ratio between segments. That fixed
single restarts (Delta_theta 0.1/3/5 -> 1.0001/1.0009/0.9999) but failed on
chains: a second segment that changes nothing computes a ratio of 1 and so
an offset of 0, while its restored t is already displaced by the first
segment's offset. Measured on theta 2 -> 7 -> 7: segment 0 at 1.0003,
segment 1 at 1.2564.

v2 derives the offset from the restored t alone -- nearest period boundary
-- which is provenance-free and self-heals over any chain depth. Since the
formula changed, the single-segment cases have to be re-run too, not just
the chain.

Usage:
    uv run python scripts/verify_phase_fix_v2.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix2"
RPM, TARGET, FIDELITY = 32.5, 7.0, 7
runs_dir = Path(PROJECT_ROOT) / "runs"

SOURCES = {7.0: "settling_baseline_L7", 6.9: "dtheta_src_L7_th6.9",
           4.0: "settling_ref_L7_rpm32.5_th4", 2.0: "settling_ref_L7_rpm32.5_th2"}


def build(source_run, theta_source, cycles, n_segments):
    src_dir = runs_dir / source_run
    t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
    return {
        "motion": {"omega_b": omega_b_of(RPM), "theta_max": [TARGET, 0.0, 0.0]},
        "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
        "n_mix_cycles": cycles, "n_transition_cycles": cycles, "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM)] * n_segments},
        "initial_checkpoint": {
            "t_dump": t_dump, "omega_b": omega_b_of(RPM),
            "theta_max": [theta_source, 0.0, 0.0],
            "checkpoint_path": str(src_dir / "checkpoint.dump"),
        },
        "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
        "binary": BINARY, "submit": True, "exclude": "node2336",
    }


for theta_source, source_run in SOURCES.items():
    ids = chain.submit_chain(build(source_run, theta_source, 120, 1))
    print(f"single dtheta={TARGET - theta_source:g}: {ids}")

print(f"chain theta 2 -> 7 -> 7: "
      f"{chain.submit_chain(build('settling_ref_L7_rpm32.5_th2', 2.0, 90, 2))}")
