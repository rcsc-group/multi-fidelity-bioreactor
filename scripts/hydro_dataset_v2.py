"""Protocol-consistent hydrodynamic QoIs from the raw shear_stress.dat series (replaces hydro_dataset.py).

Audit (scripts/hydro_windows.py, diary 2026-10-05): the stored tau_mean_max / tau_95_qss use the window
"3 periods < t < tracer release", i.e. 77 cycles including spin-up for L6-L9, but only 3 cycles after an rpm
change for the chained L10 runs; tau_mean_max is also a MAXIMUM over that window, so it grows with its length.
Here every QoI is computed over the SAME window: the last N_CYC whole cycles before the tracer release:
  tau_mean_t  = time mean over the window of the spatially averaged shear stress
  tau_95_t    = time median over the window of the per-step 95th percentile
Runs (post 2026-09-15 mask fix): L6-L8 lean (continuous from rest), L9 fig9 (prod). L10: the chained
cross-level runs are kept but flagged (window right after an rpm change).
Writes experiments/multifidelity/hydro_dataset_v2.json.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _load_dat, _t_scales  # noqa: E402

N_CYC = 3
RHO = 1.0e3


def qois(run):
    d = ROOT / "runs" / run
    p = json.load(open(d / "params.json"))
    arr = _load_dat(d / "shear_stress.dat")
    T_bio, T = _t_scales(p)
    U = p.get("geometry", {}).get("a", 0.25) / T_bio
    scale = RHO * U**2
    cyc = arr[:, 1] / T
    c_rel = p.get("t_checkpoint", 0.0) / T + p["n_mix_cycles"]
    c_end = min(math.floor(c_rel + 1e-6), math.floor(cyc.max() + 1e-6))
    m = (cyc >= c_end - N_CYC) & (cyc < c_end)
    return dict(tau_mean_t=float(arr[m, 5].mean() * scale), tau_95_t=float(np.median(arr[m, 2]) * scale),
                window=[float(c_end - N_CYC), float(c_end)], n_steps=int(m.sum()))


def main():
    old = json.load(open(ROOT / "experiments/multifidelity/hydro_dataset.json"))
    out = []
    for r in old:
        q = qois(r["run"])
        flag = "chained_rpm_change" if r["level"] == 10 else "clean"
        out.append(dict(level=r["level"], rpm=r["rpm"], run=r["run"], flag=flag, **q))
        print(f"L{r['level']} {r['rpm']:5} tau_mean_t={q['tau_mean_t']:.5f} (was max {r['tau_mean_max']:.5f}) "
              f"tau95_t={q['tau_95_t']:.5f} (was {r['tau_95_qss']:.5f}) window {q['window']} {flag}")
    json.dump(out, open(ROOT / "experiments/multifidelity/hydro_dataset_v2.json", "w"), indent=1)


if __name__ == "__main__":
    main()
