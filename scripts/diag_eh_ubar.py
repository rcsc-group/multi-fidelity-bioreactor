"""U(L9)/q(L9) by Eca & Hoekstra (2014) on the 7-deg rpm sweep, per rpm and pooled.

Same statistics as the Fig 13 MF panels (plot_mf_fig13_unc.QTY): <tau> is the
half peak-to-peak of tau_mean_signed, <eps> the peak of ediss_mean, cycle-averaged.
Level sets: L7-L10 (4 grids, the minimum E&H allow) and L6-L10 (5 grids).
h_N = 2^(10-N), finest first.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.eca_hoekstra import numerical_uncertainty          # noqa: E402
from scripts.plot_mf_fig13_unc import stat, CONV_RPMS          # noqa: E402

RUNS = {10: ("l10_fig13a_xlevel_rpm{:g}", 0.0), 9: ("fig9_l9_rpm{:g}", 0.6),
        8: ("fig13a_rm_mf_rpm{:g}", 0.6), 7: ("mf_l7_rpm{:g}", 0.6),
        6: ("fig13a_l6_mf_rpm{:g}", 0.6)}
SETS = {"L7-L10": (10, 9, 8, 7), "L6-L10": (10, 9, 8, 7, 6)}


def main() -> None:
    for kind in ("tau", "ediss"):
        print(f"\n== {kind}")
        for label, levels in SETS.items():
            rel = []
            print(f"  {label}")
            for rpm in CONV_RPMS:
                q = np.array([stat(RUNS[N][0].format(rpm), kind, RUNS[N][1])[0] for N in levels])
                h = np.array([2.0 ** (10 - N) for N in levels])
                r = numerical_uncertainty(h, q, index=levels.index(9))
                u = r["U"] / q[levels.index(9)]
                rel.append(u)
                p = "-" if r["p"] is None else f"{r['p']:.2f}"
                print(f"    {rpm:5g} rpm  q = {np.array2string(q, precision=4)}  "
                      f"est {r['estimator']:>2}{'w' if r['weighted'] else ' '} p {p:>5} "
                      f"Fs {r['Fs']:.2f}  sig/Dq {r['sigma'] / r['data_range']:.2f}  U9/q9 {u:.3f}")
            print(f"    pooled RMS u_bar = {math.sqrt(np.mean(np.square(rel))):.3f}")


if __name__ == "__main__":
    main()
