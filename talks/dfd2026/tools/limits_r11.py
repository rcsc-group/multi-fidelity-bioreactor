"""Round 11: colour limits = 10th and 85th percentiles of pooled liquid |u| (common 256 grid, each level equal weight).
Saturation = fraction of each level's full-resolution liquid cells above vmax / below vmin.
Output: data/limits_r11.json. Frames: phase 0.25 as in fig_s6_row1_r10.py."""
import json, struct, sys
from pathlib import Path
import numpy as np
REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from scripts.fields import load_snapshot, snapshots, liquid_mask
from scripts.postprocess import _t_scales
TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
RUNS = [("L8", "flow_l8_rpm32.5_th7"), ("L9", "scal_l9_rpm32.5_th7"), ("L10", "v_l10_rpm32.5_th7")]
spd, msk, picked = {}, {}, {}
for lab, run in RUNS:
    d = REPO / "runs" / run
    _, T = _t_scales(json.load(open(d / "params.json")))
    files = snapshots(d); tt = []
    for f in files:
        with open(f, "rb") as fh:
            fh.read(8); (t,) = struct.unpack("<d", fh.read(8)); tt.append(t)
    ph = (np.array(tt) / T) % 1.0
    dph = np.minimum(np.abs(ph - 0.25), 1 - np.abs(ph - 0.25))
    picked[lab] = files[int(np.argmin(dph))].name
    s = load_snapshot(files[int(np.argmin(dph))])
    msk[lab] = liquid_mask(s)
    spd[lab] = np.where(msk[lab], np.hypot(s["ux"], s["uy"]), np.nan)
def block(a, m, k):
    n = a.shape[0] // k
    A = (a * m).reshape(n, k, n, k).sum(axis=(1, 3)); M = m.reshape(n, k, n, k).sum(axis=(1, 3))
    return np.where(M >= 0.5 * k * k, A / np.maximum(M, 1), np.nan)
pool = []
for lab in spd:
    k = spd[lab].shape[0] // 256
    pool.append(block(np.nan_to_num(spd[lab]), msk[lab], k)[np.isfinite(block(np.nan_to_num(spd[lab]), msk[lab], k))])
pool = np.concatenate(pool)
PHI = float(sys.argv[1]) if len(sys.argv) > 1 else 85.0
vmin, vmax = float(np.percentile(pool, 10)), float(np.percentile(pool, PHI))
sat = {}
for lab in spd:
    v = spd[lab][msk[lab]]
    sat[lab] = dict(n_liquid=int(v.size), above_vmax=float(np.mean(v > vmax)), below_vmin=float(np.mean(v < vmin)))
    print(lab, sat[lab])
print(f"vmin {vmin:.4f} vmax {vmax:.4f}")
json.dump(dict(vmin=vmin, vmax=vmax, percentiles="10th and 85th of pooled liquid |u|, common 256 grid, each level equal weight",
               saturation=sat, picked=picked), open(TALK / "data" / "limits_r11.json", "w"), indent=1)
