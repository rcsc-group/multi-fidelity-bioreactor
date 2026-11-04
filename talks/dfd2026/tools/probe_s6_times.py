"""Probe: header times of snap_*.bin for the L8/L9/L10 triple at 32.5 rpm, theta 7.

Reads only the 16-byte header of each frame (n, nplanes, t). Used to choose
phase-matched frames for fig_s6_chain.py. Data (mtime in stat output):
  runs/flow_l8_rpm32.5_th7/fields/snap_*.bin
  runs/scal_l9_rpm32.5_th7/fields/snap_*.bin
  runs/v_l10_rpm32.5_th7/fields/snap_*.bin
"""
import struct
from pathlib import Path
import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs")
for run in ["flow_l8_rpm32.5_th7", "scal_l9_rpm32.5_th7", "v_l10_rpm32.5_th7"]:
    files = sorted((ROOT / run / "fields").glob("snap_*.bin"))
    ts, ns = [], None
    for p in files:
        with open(p, "rb") as fh:
            n, npl = struct.unpack("<ii", fh.read(8))
            (t,) = struct.unpack("<d", fh.read(8))
        ts.append(t); ns = (n, npl)
    ts = np.array(ts)
    dt = np.diff(ts)
    print(f"{run}: frames={len(ts)} n={ns[0]} nplanes={ns[1]} t0={ts[0]:.5f} t1={ts[-1]:.5f} "
          f"dt_median={np.median(dt):.6f} dt_min={dt.min():.6f} dt_max={dt.max():.6f}")
    print("  first 4 t:", np.round(ts[:4], 5))
