"""When is steady streaming actually established?

Kim releases the tracer at t/T_p = 80. Every Fig 9/10/11/12 run therefore
pays 80 cycles of spin-up before the measurement even starts -- 38% of a
32.5 rpm run and 75% of a 37.5 rpm one.

We already own eight converged L9 checkpoints (the Fig 13a sweep, 20-37.5
rpm) and one at L10, but they stop at t_end = 15.18 non-dim = 25 cycles.
Warm-starting Fig 9 from them would skip most of that spin-up -- but only if
the flow is genuinely settled by then FOR MIXING PURPOSES.

That caveat is the point. This project's settling studies measured tau, a
wall quantity that settles in ~3-7 cycles. Mixing is driven by STEADY
STREAMING, a slow secondary flow that can take far longer to establish --
which is presumably why Kim waits 80 cycles rather than 10.

This script answers the question from data already on disk, with no new
compute: it reads the liquid-phase mean vorticity from a completed run's
normf.dat and reports its cycle-averaged value over time. If that plateaus
well before cycle 80, an early tracer release is defensible and the whole
warm-start saving unlocks. If it is still drifting at cycle 25, Kim's 80
cycles must be honoured and only the cross-level part of the saving is real.

Usage:  uv run python scripts/probe_streaming_settle.py <run_id> [run_id ...]
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
_COL_T = 1
_COL_OMEGA_AVG = 2


def t_per_nd(params: dict) -> float:
    omega_b = params["omega_b"]
    L = params["geometry"]["a"]
    H = 2 * params["geometry"]["b"]
    th = math.radians(params["theta_max"][0])
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(th))
    U = V / (H * 0.5) / T_per
    return T_per / (L / U)


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    for run_id in sys.argv[1:]:
        rd = ROOT / "runs" / run_id
        nf = rd / "normf.dat"
        if not nf.exists():
            print(f"{run_id}: no normf.dat")
            continue
        params = json.loads((rd / "params.json").read_text())
        Tp = t_per_nd(params)
        rows = [l for l in nf.read_text().splitlines()
                if l.strip() and not l.strip().startswith("i")]
        a = np.array([[float(x) for x in l.split()] for l in rows])
        t, om = a[:, _COL_T], a[:, _COL_OMEGA_AVG]
        cyc = t / Tp

        # cycle-averaged |omega|: the oscillation averages out, what is left
        # is the slow secondary (streaming) component
        last = int(np.floor(cyc.max()))
        per_cycle = []
        for k in range(last):
            m = (cyc >= k) & (cyc < k + 1)
            if m.sum() > 3:
                per_cycle.append((k, float(np.abs(om[m]).mean())))
        if not per_cycle:
            print(f"{run_id}: too short")
            continue
        ks = np.array([p[0] for p in per_cycle])
        vs = np.array([p[1] for p in per_cycle])
        ref = vs[-1]

        print(f"\n{run_id}  (fidelity {params.get('fidelity')}, "
              f"{params['omega_b']*60/(2*math.pi):.1f} rpm, {last} cycles)")
        print(f"  cycle-mean |omega| at final cycle = {ref:.6g}")
        for probe in (5, 10, 25, 40, 60, 80):
            if probe < len(vs):
                dev = 100 * (vs[probe] - ref) / ref
                print(f"  cycle {probe:>3}: {vs[probe]:.6g}  ({dev:+.2f}% vs final)")
        # first cycle after which it stays within 1% of final
        within = np.abs(vs - ref) / abs(ref) < 0.01
        idx = len(vs)
        for i in range(len(vs)):
            if within[i:].all():
                idx = ks[i]
                break
        print(f"  --> settles within 1% of final from cycle {idx}")


if __name__ == "__main__":
    main()
