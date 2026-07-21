"""Submit a single run with VIDEOS=1 binary to produce mp4s alongside outputs.

Usage:
    python scripts/submit_video_run.py [params.json]

If no params.json is given, submits the canonical L6 health parameters as
'health_l6_video' so the health run folder is not overwritten.

Build the video binary first: make build-video
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.simulate import submit_slurm

root     = Path(__file__).parents[1]
template = root / "config" / "slurm_video_template.sh"

if len(sys.argv) > 1:
    params = json.loads(Path(sys.argv[1]).read_text())
else:
    params = {
        "run_id":   "health_l6_video",
        "fidelity": 6,
        "omega_b": 3.93, "n_harmonics": 1,
        "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_h": 0.0, "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
        "fill_level": 0.5,
        "t_end": 100.0,
    }

job_id = submit_slurm(params, project_root=root, walltime="02:00:00", template=template)
print(f"Submitted job {job_id}  run_id={params['run_id']}  run_dir=runs/{params['run_id']}/")
print("Expected outputs: volume_fraction.mp4, volume_fraction_lab.mp4")
