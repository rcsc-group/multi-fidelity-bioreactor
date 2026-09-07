"""Resubmit theta_source=6.9 -> target=4 (run_id=35e2e274), excluding
node2336 again (diary.md 2026-09-07).

Usage:
    uv run python scripts/resubmit_dtheta_th4_th69_exclnode.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
params = json.loads((root / "runs" / "35e2e274" / "params.json").read_text())
src_dir = root / "runs" / "dtheta_src_L6_th6.9"

job_id = submit_slurm(
    params, project_root=root, runs_root=root / "runs", walltime="06:00:00",
    template=root / "config" / "slurm_mpi_template.sh",
    checkpoint=str((src_dir / "checkpoint.dump").resolve()),
    cpus=1, ntasks=8, mem="4G", exclude="node2336",
)
print(f"Submitted theta_source=6.9 -> target=4: job {job_id}  run_id=35e2e274  (node2336 excluded)")
