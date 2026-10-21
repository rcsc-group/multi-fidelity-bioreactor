"""Split kLa = kL x a and test the small-eddy scaling on each factor (zero compute).
a = time-mean interface length / liquid area (vol_frac_interf.dat col 3 / col 2, 1/m with L_bio = 0.25 m); kL = kLa(C* = 0.10) / a.
Regress ln a, ln kL and ln kLa on ln EDR per level. Small-eddy theory (Lamont-Scott): kL ~ (eps nu)^(1/4) Sc^(-1/2), slope 0.25 on ln EDR for kL.
Usage: uv run python scripts/physics_law_split.py
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
from scripts.postprocess import _load_dat  # noqa: E402
from scripts.test_h1_mfbml import hf_run  # noqa: E402

L_BIO = 0.25


def interface_a(run):
    v = _load_dat(ROOT / "runs" / run / "vol_frac_interf.dat")
    m = v[:, 2] > 0
    n = int(m.sum())
    sel = np.where(m)[0][n // 5:]
    return float(np.mean(v[sel, 3] / v[sel, 2]) / L_BIO)


def kla(run):
    f = ROOT / "runs" / run / "results.json"
    v = json.loads(f.read_text()).get("kLa_1T_10") if f.exists() else None
    return float(v) / 3600.0 if v and v > 0 else np.nan


def rows_for(level):
    out = []
    if level == 8:
        runs = [(d.name, d.name) for pat in ("t6_l8_*", "t5_l8_*", "e8_l8_*") for d in sorted((ROOT / "runs").glob(pat))]
    elif level == 9:
        runs = [(d.name, d.name) for d in sorted((ROOT / "runs").glob("fig9_l9_rpm*")) if "_ext" not in d.name]
    else:
        runs = []
        for d in sorted((ROOT / "runs").glob("t10_l10_*")) + [ROOT / "runs" / "v_l10_rpm32.5_th7"]:
            p = json.loads((d / "params.json").read_text())
            r, t = p["omega_b"] * 60 / (2 * math.pi), p["theta_max"][0]
            runs.append((d.name, f"x8_l10_rpm{r:g}_th5.5" if t == 5.5 else hf_run(r, t)))
    for name, src in runs:
        k = kla(name)
        try:
            h = qois(src)
            a = interface_a(name)
        except Exception:
            continue
        if np.isfinite(k) and h:
            out.append((a, k / a, k, h["edr_mean_t"]))
    return np.array(out)


def slope(x, y):
    A = np.column_stack([np.ones(len(x)), x])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return c[1], np.corrcoef(x, y)[0, 1]


def main():
    for L in (8, 9, 10):
        T = rows_for(L)
        e = np.log(T[:, 3])
        sa, ra = slope(e, np.log(T[:, 0]))
        sk, rk = slope(e, np.log(T[:, 1]))
        sl, rl = slope(e, np.log(T[:, 2]))
        print(f"L{L}: n {len(T)}  a mean {T[:, 0].mean():.1f} 1/m (cv {T[:, 0].std() / T[:, 0].mean():.2f});  slope on ln EDR: a {sa:+.2f} (r {ra:+.2f}), kL {sk:+.2f} (r {rk:+.2f}), kLa {sl:+.2f} (r {rl:+.2f});  kL mean {1000 * T[:, 1].mean():.3f} mm/s")


if __name__ == "__main__":
    main()
