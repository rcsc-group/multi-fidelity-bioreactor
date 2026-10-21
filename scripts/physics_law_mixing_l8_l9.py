"""Mixing-time law across levels (zero compute): ln dtmix_0.50 = c_L + s ln EDR with the SLOPE s fitted on L8 (88 runs) and only the offset c_9 fitted on L9 runs.
Leave-the-rpm-group-out over the 14 L9 theta-7 runs (runs within 0.25 rpm removed from the offset fit). Compared with the S89 chain (0.222), L9-only GP (0.376) and L8 x ratio (0.333).
Usage: uv run python scripts/physics_law_mixing_l8_l9.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402


def rows(patterns):
    out = []
    for pat in patterns:
        for d in sorted((ROOT / "runs").glob(pat)):
            f = d / "results.json"
            if "_ext" in d.name or not f.exists():
                continue
            v = json.loads(f.read_text()).get("dtmix_0.50")
            try:
                h = qois(d.name)
            except Exception:
                continue
            if v and v > 0 and h:
                p = json.loads((d / "params.json").read_text())
                out.append((p["omega_b"] * 60 / (2 * math.pi), p["theta_max"][0], h["edr_mean_t"], v))
    return np.array(out)


def main():
    l8 = rows(("t6_l8_*", "t5_l8_*", "e8_l8_*"))
    l9 = rows(("fig9_l9_rpm*", "e9_l9_rpm22.4_th7", "e9_l9_rpm22.6_th7", "e9_l9_rpm29.9_th7", "e9_l9_rpm30.1_th7"))
    s8, c8 = np.polyfit(np.log(l8[:, 2]), np.log(l8[:, 3]), 1)
    s9, c9 = np.polyfit(np.log(l9[:, 2]), np.log(l9[:, 3]), 1)
    print(f"L8: slope {s8:.2f} (n {len(l8)}); L9: slope {s9:.2f} (n {len(l9)}); offset L9 - L8 at the same EDR: x{np.exp(np.mean(np.log(l9[:, 3]) - s8 * np.log(l9[:, 2])) - c8):.2f}")
    e_law, e_free, e_const = [], [], []
    for i in range(len(l9)):
        tr = np.abs(l9[:, 0] - l9[i, 0]) > 0.25
        off = np.mean(np.log(l9[tr, 3]) - s8 * np.log(l9[tr, 2]))
        e_law.append(off + s8 * np.log(l9[i, 2]) - np.log(l9[i, 3]))
        b = np.polyfit(np.log(l9[tr, 2]), np.log(l9[tr, 3]), 1)
        e_free.append(np.polyval(b, np.log(l9[i, 2])) - np.log(l9[i, 3]))
        e_const.append(np.mean(np.log(l9[tr, 3])) - np.log(l9[i, 3]))
    rms = lambda e: float(np.sqrt(np.mean(np.square(e))))   # noqa: E731
    print(f"L9 leave-group-out RMS ln error: law with L8 slope + L9 offset {rms(e_law):.3f}; law with L9 slope and offset {rms(e_free):.3f}; L9 constant {rms(e_const):.3f}  (chain from L8: 0.222, L9-only GP 0.376, L8 x ratio 0.333)")


if __name__ == "__main__":
    main()
