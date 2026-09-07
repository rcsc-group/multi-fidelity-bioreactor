"""Resubmit the full-length p-fix verification (run_id=dea0e98d),
excluding node2336.

Usage:
    uv run python scripts/resubmit_p_fix_verification.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
params = json.loads((root / "runs" / "dea0e98d" / "params.json").read_text())
src_dir = root / "runs" / "dtheta_src_L6_th4.1"

job_id = submit_slurm(
    params, project_root=root, runs_root=root / "runs", walltime="06:00:00",
    template=root / "config" / "slurm_mpi_template.sh",
    checkpoint=str((src_dir / "checkpoint.dump").resolve()),
    cpus=1, ntasks=8, mem="4G", exclude="node2336",
)
print(f"Resubmitted: job {job_id}  run_id=dea0e98d  (node2336 excluded)")
