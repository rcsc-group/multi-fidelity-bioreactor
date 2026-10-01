"""Eca & Hoekstra (2014) grid uncertainty of the L9 kLa at 32.5 rpm, 7 deg.

The only condition with kLa at 4 levels: L6 fig9_validate_l6, L7 fig9_ladder_l7,
L8 fig9_ladder_l8, L9 fig9_l9 (all release at cycle 80). The L10 anchor
(fig9_l10_seg2) is excluded: its kLa is contaminated by a restart transient at
the seg1->seg2 seam (diary 2026-10-01). h_N = 2^(10-N).
Question: is U(L9)/kLa(L9) large enough that the L9 kLa curve cannot serve as
"truth" for the MF tests, and that the H1 failure could be grid noise?
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.eca_hoekstra import numerical_uncertainty  # noqa: E402

RUNS = {9: "fig9_l9_rpm32.5", 8: "fig9_ladder_l8_rpm32.5",
        7: "fig9_ladder_l7_rpm32.5", 6: "fig9_validate_l6_rpm32.5"}


def main() -> None:
    for key in ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "kLa_25", "kLa_50",
                "dtmix_0.50", "dtmix_0.95"):
        lv = [N for N in (9, 8, 7, 6)
              if np.isfinite(json.load(open(ROOT / "runs" / RUNS[N] / "results.json")).get(key, np.nan))]
        q = np.array([json.load(open(ROOT / "runs" / RUNS[N] / "results.json"))[key] for N in lv])
        line = f"{key:11s} " + " ".join(f"L{N}={v:8.3g}" for N, v in zip(lv, q))
        if len(lv) < 4:
            print(line + "   (fewer than 4 levels)")
            continue
        h = np.array([2.0 ** (10 - N) for N in lv])
        r = numerical_uncertainty(h, q, index=0)
        p = "-" if r["p"] is None else f"{r['p']:.2f}"
        print(line + f"   est {r['estimator']}{'w' if r['weighted'] else ''} p {p} Fs {r['Fs']} "
              f"phi0 {r['phi0']:.3g}  U(L9)/L9 = {r['U'] / q[0]:.0%}")


if __name__ == "__main__":
    main()
