"""Metrics M1-M3 of the offline-vs-DNS test (diary 2026-10-29), Kim's solver: every run's prereg_diag.dat against the
producer's (DNS A) at the producer's step times. Columns: i t dt S0 S1 S2 O (S1 = sum dv c2 f, S2 = sum dv c2^2 f,
S0 = sum dv f, O = sum dv oxy f). Release row = first row with S1 > 0. Time in cycles of T_PER from the release.
Usage: uv run python compare_kim.py <run tag> [floor tag]   (from any dir; writes <prereg_out>/<tag>/metrics.txt and series.csv)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

OUTROOT = Path("/oscar/scratch/eaguerov/kim_pipeline/prereg_out")
T_PER = 0.607329   # the producer's T_per_st at 20 rpm, theta 7 (printed by pre_open as the record span / 3); checked below
VARS = ("off_m1", "off_m1f0", "off_m3", "off_m3f0", "off_m4", "his_S", "his_H", "his_ship", "his_loop")


def load(p: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    a = np.loadtxt(p, skiprows=1)
    k = int(np.argmax(a[:, 4] > 0))
    a = a[k:]
    s0, s1, s2, o = a[:, 3], a[:, 4], a[:, 5], a[:, 6]
    var = s2 / s0 - (s1 / s0) ** 2
    return a[:, 1], s1 / s1[0], 1 - var / var[0], o


def main() -> None:
    tag = sys.argv[1]
    base = OUTROOT / tag
    ta, ma, xa, oa = load(base / "producer/prereg_diag.dat")
    t0 = ta[0]
    rows, cols = [], {"t_cyc": (ta - t0) / T_PER, "A_m": ma, "A_chi": xa, "A_oxy": oa}
    runs = [(v, base / v / "prereg_diag.dat") for v in VARS]
    if len(sys.argv) > 2:
        runs.insert(0, ("floor_B", OUTROOT / sys.argv[2] / "producer/prereg_diag.dat"))
    for name, p in runs:
        if not p.exists():
            rows.append(f"{name:9s} missing")
            continue
        t, m, x, o = load(p)
        w = (ta > ta[0]) & (ta <= t[-1] + 1e-12) & (ta - t0 <= 3 * T_PER + 1e-12)
        mi, xi, oi = (np.interp(ta[w], t, v) for v in (m, x, o))
        m1 = np.max(np.abs(mi - ma[w]))
        m2 = np.max(np.abs(xi - xa[w]))
        w3 = (ta[w] - t0) / T_PER >= 0.25        # M3 window [0.25, 3] cycles (diary amendment 2)
        m3 = np.nanmax(np.abs(oi[w3] / oa[w][w3] - 1)) if w3.any() else np.nan   # run ended before 0.25 cycle
        dt0 = (t[0] - t0) / T_PER
        rows.append(f"{name:9s} M1 {m1:.3e}  M2 {m2:.3e}  M3 {m3:.3e}  to {(t[-1] - t0) / T_PER:.3f} cyc  "
                    f"start offset {dt0:+.2e} cyc  rows {len(t)}")
        for k, v in (("m", m), ("chi", x), ("oxy", o)):
            cols[f"{name}_{k}"] = np.interp(ta, t, v, left=np.nan, right=np.nan)
    (base / "metrics.txt").write_text("\n".join(rows) + "\n")
    keys = list(cols)
    np.savetxt(base / "series.csv", np.column_stack([cols[k] for k in keys]), delimiter=",",
               header=",".join(keys), comments="")
    print(f"producer rows {len(ta)}, span {(ta[-1] - t0) / T_PER:.4f} cycles")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
