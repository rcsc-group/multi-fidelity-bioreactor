"""Audit of the time window behind tau_mean_max / tau_95_qss for each run of hydro_dataset.json.

postprocess._compute_tau98_kpis uses the window t_ramp (3 periods) < t < t_inject. This prints, per run, the
window in rocking cycles (start, end, length) and the per-cycle mean of the spatially averaged tau over the
last cycles, so that drift (spin-up or restart transients) and window-length effects are visible.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _load_dat, _ramp_end_nd, _t_scales  # noqa: E402


def main():
    rows = json.load(open(ROOT / "experiments/multifidelity/hydro_dataset.json"))
    for r in rows:
        d = ROOT / "runs" / r["run"]
        p = json.load(open(d / "params.json"))
        arr = _load_dat(d / "shear_stress.dat") if (d / "shear_stress.dat").exists() else None
        if arr is None:
            print(f"L{r['level']} {r['rpm']:5} {r['run']}: no shear_stress.dat")
            continue
        _, T = _t_scales(p)
        t = arr[:, 1]
        cyc = t / T
        nmix = p.get("n_mix_cycles", math.nan)
        c0 = p.get("t_checkpoint", 0.0) / T
        ramp = _ramp_end_nd(p) / T
        # per-cycle means of spatial-mean tau over the last 6 whole cycles before release (or end)
        c_end = min(cyc.max(), c0 + nmix) if np.isfinite(nmix) else cyc.max()
        per = []
        for k in range(6, 0, -1):
            m = (cyc >= c_end - k) & (cyc < c_end - k + 1)
            per.append(arr[m, 5].mean() if m.any() else np.nan)
        per = np.array(per) / per[-1]
        print(f"L{r['level']} {r['rpm']:5} cycles {cyc.min():6.1f}-{cyc.max():6.1f} ramp_end {ramp:5.1f} "
              f"release {c0 + nmix:6.1f}  last-6-cycle means / last = {np.round(per, 3)}  {r['run']}")


if __name__ == "__main__":
    main()
