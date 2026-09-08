"""A/B the theta_max contribution to the restart velocity rescale.

Commit 646b7b3 changed su from omega_b_prev/omega_b to the full
U_bio_prev/U_bio ratio and called it a confirmed bug fix. It was not
confirmed: its one verification run moved the metric from 5.8% to 5.4%,
which is nothing, and the change was argued from dimensional analysis
alone. The restart velocity rescale itself has only ever been touched twice
(fbc500d introduced it, 646b7b3 changed it), but this project's
nondimensionalization has been corrected wrongly or partially many times
(H_bio off by 2x in 0482dd7 -> d035047 -> c8c21e1, again in postprocessing
in 052e9e4, unit conversions in 456bc4a and 02519a7, T_bio duplication in
699df66), so "it is dimensionally correct" is not sufficient grounds.

It also now sits stacked with the restart phase-offset fix, so the
verification sweep cannot attribute a result to either one.

This runs two theta-changing restarts with the phase fix present and the su
term REVERTED (-DLEGACY_SU=1), against the same two cases in the
verification sweep which have it applied:

  theta 6.9 -> 7 : su term is 0.26% on velocity. Predicted to make no
      difference at all -- the kicknoise run applied a deliberate 0.3%
      per-cell random velocity perturbation to a converged flow and it
      decayed completely (everything back to 1.000 within 0.1% over 160
      cycles). If a larger, less structured perturbation vanishes, this one
      cannot matter.
  theta 2 -> 7   : su term is 14% on velocity. This is the case where it
      could plausibly matter, and the only one that can supply actual
      evidence for the change.

Usage:
    uv run python scripts/ab_test_su_rescale.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix-legacysu"
RPM, TARGET, FIDELITY, CYCLES = 32.5, 7.0, 7, 120
runs_dir = Path(PROJECT_ROOT) / "runs"

CASES = [(6.9, "dtheta_src_L7_th6.9"), (2.0, "settling_ref_L7_rpm32.5_th2")]

for theta_source, source_run in CASES:
    src_dir = runs_dir / source_run
    t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
    cfg = {
        "motion": {"omega_b": omega_b_of(RPM), "theta_max": [TARGET, 0.0, 0.0]},
        "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
        "n_mix_cycles": CYCLES, "n_transition_cycles": CYCLES, "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM)]},
        "initial_checkpoint": {
            "t_dump": t_dump, "omega_b": omega_b_of(RPM),
            "theta_max": [theta_source, 0.0, 0.0],
            "checkpoint_path": str(src_dir / "checkpoint.dump"),
        },
        "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
        "binary": BINARY, "submit": True, "exclude": "node2336",
    }
    print(f"legacy-su, source theta={theta_source:g}: {chain.submit_chain(cfg)}")
