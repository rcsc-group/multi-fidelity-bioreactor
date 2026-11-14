"""Machine-specific folders, in one place.

Each folder is an environment variable. The default is the value on Brown's OSCAR cluster, where this
project runs; on another machine, set the variables (for example in your shell profile):

    export MFB_SCRATCH=/big/fast/disk/me     # run folders, staged MPI runs, binaries built for jobs
    export MFB_DATA=/kept/disk/me            # kept data: flux records, shared binaries
    export BASILISK=/path/to/basilisk/src    # Basilisk's src/ (qcc lives here)

Python code uses the attributes below (`from scripts import paths`; `f"{paths.SCRATCH}/mpi_runs"`).
Shell and SLURM scripts use the same variables with the same defaults (`${MFB_SCRATCH:-/oscar/scratch/eaguerov}`).
Config files (YAML/JSON) may write `$MFB_SCRATCH/...`; the code that reads them calls `expand()`.
"""
import os
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]                                     # this checkout
SCRATCH = Path(os.environ.get("MFB_SCRATCH", "/oscar/scratch/eaguerov"))       # large, purged after ~30 days
DATA = Path(os.environ.get("MFB_DATA", "/oscar/data/dharri15/eaguerov"))       # kept
BASILISK = Path(os.environ.get("BASILISK", "/oscar/data/dharri15/eaguerov/basilisk/src"))
# The first 2D run database (2,915 runs), in the older BioReactor3D repository; figures of the Kim et al.
# comparison read it.
LEGACY_2D = Path(os.environ.get("MFB_LEGACY_2D",
                                "/oscar/data/dharri15/eaguerov/Github/BioReactor3D/dev/rocking-bioreactor-2d"))

VARS = {"MFB_REPO": REPO, "MFB_SCRATCH": SCRATCH, "MFB_DATA": DATA, "BASILISK": BASILISK,
        "MFB_LEGACY_2D": LEGACY_2D}


def expand(s: str) -> str:
    """`$MFB_SCRATCH/x` or `${MFB_SCRATCH}/x` in a config value -> the folder of this machine."""
    return re.sub(r"\$\{?(" + "|".join(VARS) + r")\}?", lambda m: str(VARS[m.group(1)]), s)
