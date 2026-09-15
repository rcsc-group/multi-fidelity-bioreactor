"""Rerun the L6 and L8 Fig13a sweeps on the mask-fixed binary (diary.md
2026-09-15 (2)).

Every tau_mean / ediss_mean in the existing L6/L8/L9 sweeps was computed
with the broken averaging mask (`f[] > 0.5` alone, which reaches outside the
embedded bag and inflates the averaging area by ~3.5x). The correction is
NOT a post-hoc multiplier -- the excluded cells carried spurious velocity
gradients, so different quantities shift by different factors (measured at
L7: 2.06x on mean|tau|, 1.66x on EDR, 0.56x on the signed mean) -- so the
sweeps have to be rerun rather than rescaled.

New run_ids (`*_mf_*`) so the original runs are preserved as the record of
what the buggy diagnostic produced.

NOTE ON QUEUE: submits to the normal batch queue. scripts/submit_fig13a_l6.py
hardcodes `--account=mbessa-condo --qos=mbessa-condo` from a one-off
permission that was revoked 2026-09-04; do not reuse that script.

Usage:
    uv run python scripts/submit_fig13a_maskfix_reruns.py [--dry-run]
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-maskfix"
DRY = "--dry-run" in sys.argv
RPMS = [17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5]

# (fidelity, run_id prefix, t_end, n_mix_cycles, ntasks, walltime, mem)
SWEEPS = [
    (6, "fig13a_l6_mf_rpm",  20.0, 80,  4, "02:00:00", "2G"),
    (8, "fig13a_rm_mf_rpm",  20.0, 80, 16, "06:00:00", "4G"),
]

jobs = []
for fid, prefix, t_end, n_mix, ntasks, walltime, mem in SWEEPS:
    for rpm in RPMS:
        omega = rpm * 2 * math.pi / 60.0
        run_id = f"{prefix}{rpm:g}"
        params = {
            "run_id": run_id, "fidelity": fid,
            "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
            "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0],
            "phi_angular": [0.0, 0.0, 0.0],
            "omega_b": omega, "omega_h": 0.0,
            "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
            "t_end": t_end, "n_mix_cycles": n_mix,
            "_binary": BINARY,
        }
        if DRY:
            job = "DRY"
        else:
            job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                               walltime=walltime,
                               template=ROOT / "config" / "slurm_mpi_template.sh",
                               cpus=1, ntasks=ntasks, mem=mem)
        jobs.append((run_id, job))
        print(f"  L{fid} rpm={rpm:<5g} {run_id:<26} job={job}")

print(f"\n{len(jobs)} jobs submitted on the normal batch queue")
