"""kL measured vs the most a resolved diffusive layer can carry (diagnostic for kLa; diary 2026-10-12).

kL = kLa_1T_25 / a, a = interface length / liquid area (vol_frac_interf.dat: f_liq_interf = interface_area(f),
f_liq_sum = liquid area; nondim, L0 = 1 = L_bio), time-mean over the run after release. Bound on a resolved
diffusive flux per unit interface: the whole saturation jump across one cell, k_res = 2 D_oxy / dx (the factor 2 is
generous: half-cell distance). Physical reference: penetration theory kL = 2 sqrt(D / (pi T)) with the rocking
period T as the renewal time (an order-of-magnitude estimate, not a model of this flow).
Runs: k_l<L>_s1e0 (L6-L9, physical D, 32.5 rpm) and l10c_rpm32.5_seg1 (L10).
Usage: uv run python scripts/diag_kl_bound.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _load_dat, _t_scales  # noqa: E402

D_OXY, L_BIO = 1.9e-9, 0.25
RUNS = {6: "k_l6_s1e0", 7: "k_l7_s1e0", 8: "k_l8_s1e0", 9: "k_l9_s1e0", 10: "l10c_rpm32.5_seg1"}


def main():
    rows = []
    for L, run in RUNS.items():
        d = ROOT / "runs" / run
        p = json.loads((d / "params.json").read_text())
        res = json.loads((d / "results.json").read_text())
        v = _load_dat(d / "vol_frac_interf.dat")
        T_bio, T_nd = _t_scales(p)
        m = v[:, 2] > 0
        a_nd = np.mean(v[m, 3] / v[m, 2])           # 1/L_bio
        a = a_nd / L_BIO                             # 1/m
        kla = res["kLa_1T_25"] / 3600.0              # 1/s
        kl = kla / a                                 # m/s
        dx = L_BIO / 2 ** L
        k_res = 2 * D_OXY / dx
        k_pen = 2 * np.sqrt(D_OXY / (np.pi * T_nd * T_bio))
        rows.append(dict(level=L, run=run, a_per_m=a, kLa_1T_25_per_h=res["kLa_1T_25"], kL=kl, dx=dx,
                         k_res=k_res, ratio=kl / k_res, k_pen=k_pen))
        print(f"L{L:2d} {run:20s} a {a:6.1f} 1/m  kLa {res['kLa_1T_25']:7.2f} 1/h  kL {kl:.2e} m/s  "
              f"2D/dx {k_res:.2e}  kL/(2D/dx) {kl / k_res:6.1f}  penetration {k_pen:.2e}")
    (ROOT / "experiments/multifidelity/diag_kl_bound.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
