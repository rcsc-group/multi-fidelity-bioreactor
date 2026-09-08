"""Is the omega escape caused by the mis-scaled g field, or is it a genuine
basin crossing?

g combines the pressure gradient with the full acceleration field, whose
four terms scale as su^2 (gravity), su (Coriolis) and 1 (centrifugal,
Euler). The restart scales all of g by su^2, so the su-invariant terms are
wrong by |su^2-1|: 24% at theta 2->7 (converges), 33% and 71% at the two
omega cases (both escape). Consistent -- but so is a plain magnitude
threshold, since the clean/escaping split also sits between su 0.874 and
1.154. These separate the two.

Runs the worst escaping case, omega 17.5 -> 32.5 (su=0.538, 71% g error,
currently 1.2308), under:
  mode 1  g left exactly as restored (no rescale)
  mode 2  g zeroed -- which is precisely what a COLD start does, since
          centered.h's event init never sets g and the first predictor
          half-step therefore runs with g==0

If either converges, g is the culprit and the fix is to stop scaling it. If
both still escape, g is not what pushes it over and the escape is a genuine
basin crossing driven by the size of the omega change.

Usage:
    uv run python scripts/test_g_restart_modes.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

src_dir = Path(PROJECT_ROOT) / "runs" / "settling_ref_L7_rpm17.5_th7"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])

for mode in (1, 2):
    cfg = {
        "motion": {"omega_b": omega_b_of(32.5), "theta_max": [7.0, 0.0, 0.0]},
        "fidelity": 7, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
        "n_mix_cycles": 120, "n_transition_cycles": 120, "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(32.5)]},
        "initial_checkpoint": {
            "t_dump": t_dump, "omega_b": omega_b_of(17.5),
            "theta_max": [7.0, 0.0, 0.0],
            "checkpoint_path": str(src_dir / "checkpoint.dump"),
        },
        "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
        "binary": f"/oscar/scratch/eaguerov/BioReactor-mpi-gmode{mode}",
        "submit": True, "exclude": "node2336",
    }
    print(f"g mode {mode}: {chain.submit_chain(cfg)}")
