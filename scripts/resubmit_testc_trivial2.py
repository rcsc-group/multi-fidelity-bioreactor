"""Resubmit Test C's trivial restart, verification run (run_id=54c658b0),
excluding node2336.

Usage:
    uv run python scripts/resubmit_testc_trivial2.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
params = json.loads((root / "runs" / "54c658b0" / "params.json").read_text())
src_dir = root / "runs" / "testc_source_th4.1"

job_id = submit_slurm(
    params, project_root=root, runs_root=root / "runs", walltime="01:00:00",
    template=root / "config" / "slurm_mpi_template.sh",
    checkpoint=str((src_dir / "checkpoint.dump").resolve()),
    cpus=1, ntasks=8, mem="4G", exclude="node2336",
)
print(f"Resubmitted: job {job_id}  run_id=54c658b0  (node2336 excluded)")
