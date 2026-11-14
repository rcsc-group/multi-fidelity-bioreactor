"""Wave height (posY_max - posY_min), cycle means every 3 cycles from the checkpoint, for given campaign runs.
Question: is the L9 wave still changing after the 12 warm cycles (unsettled), or steady (grid offset only)?"""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.postprocess import _load_dat, _t_scales  # noqa: E402
from scripts import paths  # noqa: E402  (machine-specific folders)
for r in sys.argv[1:]:
    d = Path(f"{paths.SCRATCH}/mpi_runs") / r
    p = json.loads((d / "params.json").read_text()); _, P = _t_scales(p)
    v = _load_dat(d / "vol_frac_interf.dat"); h = v[:, 4] - v[:, 5]; t = v[:, 1]
    t0 = float(p["t_checkpoint"]); n = int((t.max() - t0) / P)
    print(f"{r:24s}", " ".join(f"{h[(t >= t0 + k * P) & (t < t0 + (k + 1) * P)].mean():.4f}" for k in range(0, n, 3)))
