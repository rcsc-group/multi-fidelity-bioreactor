"""Set up a BioReactor run directory and generate its SLURM script.

Does NOT submit — call simulate.submit_slurm() for that.

Public API
----------
main(params_file, runs_root) -> dict
    Returns {"run_id": str, "run_dir": str, "slurm_script": str}
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from uuid import uuid4

_TEMPLATE = Path(__file__).parents[1] / "config" / "slurm_template.sh"


def main(
    params_file: str,
    runs_root: str | None = None,
    template: str | None = None,
) -> dict:
    """Create run directory, copy params.json, write a SLURM script.

    Parameters
    ----------
    params_file : path to the source params.json
    runs_root   : parent directory for run dirs (default: <project>/runs)
    template    : SLURM template path (default: config/slurm_template.sh)

    Returns
    -------
    dict with keys: run_id, run_dir, slurm_script
    """
    params_path = Path(params_file)
    params      = json.loads(params_path.read_text())
    fidelity    = params.get("fidelity", 7)

    project_root = Path(__file__).parents[1]
    runs_root    = Path(runs_root) if runs_root else project_root / "runs"
    template     = Path(template) if template else _TEMPLATE

    run_id  = uuid4().hex[:8]
    run_dir = runs_root / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    # Copy params.json into run directory
    dest_params = run_dir / "params.json"
    shutil.copy2(params_path, dest_params)

    # Write a SLURM script stamped with the fidelity level
    template_text = template.read_text()
    # Insert LEVEL=<fidelity> after the shebang line so it's visible in the script
    lines = template_text.splitlines(keepends=True)
    insert_at = 1  # after #!/bin/bash
    lines.insert(insert_at, f"LEVEL={fidelity}  # fidelity level for this run\n")
    slurm_script = run_dir / "run.sh"
    slurm_script.write_text("".join(lines))

    return {
        "run_id":       run_id,
        "run_dir":      str(run_dir),
        "slurm_script": str(slurm_script),
    }


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: launch.py <params_file> [runs_root]", file=sys.stderr)
        sys.exit(1)
    runs_root = sys.argv[2] if len(sys.argv) > 2 else None
    result = main(sys.argv[1], runs_root=runs_root)
    print(json.dumps(result, indent=2))
