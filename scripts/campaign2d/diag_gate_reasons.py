"""Print the settled-gate details (drift, height change vs the L8 source) for every campaign c2d_l9 run."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import campaign as C  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root
from scripts import paths  # noqa: E402  (machine-specific folders)
for d in sorted(Path(f"{paths.SCRATCH}/mpi_runs").glob("c2d_l9_*")):
    v, info = C.warm_run_gate(d)
    print(d.name, v, json.dumps({k: (round(x, 4) if isinstance(x, float) else x) for k, x in info.items()}))
