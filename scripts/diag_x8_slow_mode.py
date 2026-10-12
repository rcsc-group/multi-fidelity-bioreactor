"""What is slow after an L8 -> L10 warm start at 17.5 rpm? (diagnostic, not pre-registered; diary 2026-10-12)

Hypothesis H_stream: the slow drift of tau_mean after the restart is the period-averaged mean circulation (streaming)
adjusting, not the oscillatory part. Proxy: the signed mean wall shear tau_mean_signed (shear_stress.dat col 11, the
wall average of signed tau), whose period average is the mean-circulation wall shear (the oscillatory part averages
out). Prediction if H_stream: per-cycle tau_mean_signed drifts with the same time scale as tau_mean, while the
x8 run approaches the xlevel reference value. If tau_mean_signed is flat while tau_mean drifts, H_stream is rejected.

Per cycle (cycles counted from the L8 -> L10 restart; x8 then x8e): tau_mean, tau_mean_signed, EDR (nondimensional,
as logged). Reference: the xlevel run's last 3 cycles before release. Usage: uv run python scripts/diag_x8_slow_mode.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.test_warmstart_x8 import load  # noqa: E402
from scripts.postprocess import _load_dat  # noqa: E402

COLS = {"tau_mean": 5, "tau_mean_signed": 10, "ediss": 9}


def per_cycle(r, c0, n):
    out = []
    for k in range(n):
        m = (r["cyc"] >= c0 + k) & (r["cyc"] < c0 + k + 1)
        if m.sum() < 25:
            break
        a = r["arr"][m]
        out.append({q: float(a[:, j].mean()) for q, j in COLS.items()})
    return out


NCOLS = {"Omega_liq_avg": 2, "ux_liq_avg": 6}  # normf.dat: period mean of Omega_liq_avg = net mean circulation


def normf_cycles(run, r, c0, n):
    """Per-cycle means of normf.dat columns (period average of the liquid-averaged vorticity = mean circulation
    of the liquid / liquid area, by Stokes' theorem; the oscillatory part averages out over one period)."""
    f = ROOT / "runs" / run / "normf.dat"
    if not f.exists():  # running job: read the live file in scratch
        f = Path("/oscar/scratch/eaguerov/mpi_runs") / run / "normf.dat"
    a = _load_dat(f)
    if a is None or len(a) == 0:
        return []
    cyc = a[:, 1] / r["T"]
    out = []
    for k in range(n):
        m = (cyc >= c0 + k) & (cyc < c0 + k + 1)
        if m.sum() < 25:
            break
        out.append({q: float(a[m, j].mean()) for q, j in NCOLS.items()})
    return out


def main():
    for rpm in ("17.5", "25"):
        x8, xe, xl = load(f"x8_l10_rpm{rpm}_th7"), load(f"x8e_l10_rpm{rpm}_th7"), load(f"l10_fig13a_xlevel_rpm{rpm}")
        if x8 is not None and xl is not None:
            nf = normf_cycles(f"x8_l10_rpm{rpm}_th7", x8, x8["p"]["t_checkpoint"] / x8["T"], 6)
            pe = ROOT / "runs" / f"x8e_l10_rpm{rpm}_th7" / "params.json"
            if pe.exists():  # continuation, finished or still running (same T as x8)
                te = json.loads(pe.read_text())["t_checkpoint"] / x8["T"]
                nf += normf_cycles(f"x8e_l10_rpm{rpm}_th7", x8, te, 7)
            c_rel = xl["p"]["t_checkpoint"] / xl["T"] + xl["p"]["n_mix_cycles"]
            nref = normf_cycles(f"l10_fig13a_xlevel_rpm{rpm}", xl, np.floor(c_rel) - 3, 3)
            nxl = normf_cycles(f"l10_fig13a_xlevel_rpm{rpm}", xl, xl["p"]["t_checkpoint"] / xl["T"], 8)
            src = f"kmix_l8_rpm{rpm}"
            ls = load(src)
            nsrc = normf_cycles(src, ls, np.floor(x8["p"]["t_checkpoint"] / ls["T"]) - 3, 3) if ls else []
            fmt = lambda s: ", ".join(f"{c['Omega_liq_avg']:.4g}" for c in s)  # noqa: E731
            print(f"=== {rpm} rpm, period-mean Omega_liq_avg (mean circulation / liquid area)")
            print(f"   L8 source, last 3 cycles before dump: {fmt(nsrc)}")
            print(f"   x8 (+x8e) cycles 1..: {fmt(nf)}")
            print(f"   xlevel after its L9 -> L10 restart: {fmt(nxl)}")
            print(f"   xlevel reference (last 3): {fmt(nref)}")
        if x8 is None or xl is None:
            print(rpm, "missing")
            continue
        seq = per_cycle(x8, x8["p"]["t_checkpoint"] / x8["T"], 6)
        if xe is not None:
            seq += per_cycle(xe, xe["p"]["t_checkpoint"] / xe["T"], 7)
        c_rel = xl["p"]["t_checkpoint"] / xl["T"] + xl["p"]["n_mix_cycles"]
        ref = per_cycle(xl, np.floor(c_rel) - 3, 3)
        refm = {q: np.mean([c[q] for c in ref]) for q in COLS}
        print(f"--- {rpm} rpm; xlevel ref (last 3 cycles): " + ", ".join(f"{q} {v:.4g}" for q, v in refm.items()))
        print(f"{'cyc':>4s} " + " ".join(f"{q:>16s} {'/ref':>6s}" for q in COLS))
        for i, c in enumerate(seq, 1):
            print(f"{i:4d} " + " ".join(f"{c[q]:16.5g} {c[q] / refm[q]:6.3f}" for q in COLS))
        # xlevel's own transient after ITS restart (L9 -> L10), for comparison
        xs = per_cycle(xl, xl["p"]["t_checkpoint"] / xl["T"], 8)
        print("   xlevel after its own L9 -> L10 restart, tau_mean/ref: "
              + ", ".join(f"{c['tau_mean'] / refm['tau_mean']:.3f}" for c in xs)
              + " | tau_mean_signed/ref: " + ", ".join(f"{c['tau_mean_signed'] / refm['tau_mean_signed']:.3f}" for c in xs))


if __name__ == "__main__":
    main()
