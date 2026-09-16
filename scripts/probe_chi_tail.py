"""Is dtmix_0.95 measuring mixing, or a numerical floor under sigma^2?

dt95/Kim jumped 0.154 -> 0.371 from L7 to L8 while dt50/Kim moved only
0.133 -> 0.192. Three thresholds of the same curve diverging that way points
at dt95 being ill-conditioned rather than at three different physics.

Candidate: postprocess normalises the tracer integrals by `f_mean`, the
TIME-AVERAGE of the liquid volume, not its instantaneous value. Liquid
volume oscillates as the bag rocks, so c_mean = c_sum/f_mean carries a
spurious oscillation and sigma^2 = <c^2> - <c>^2 gets a spurious FLOOR.
Once true sigma^2 falls near that floor, chi stops climbing, and the time at
which it crosses 0.95 becomes hypersensitive -- or never happens. dt50 sits
far above the floor and is barely affected.

This prints, per run, the chi curve computed BOTH ways -- with the
time-averaged f (what postprocess does) and with the instantaneous f -- plus
the noise amplitude in the tail, so the floor is visible if it exists.

Usage:  uv run python scripts/probe_chi_tail.py <run_id> [run_id ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
_COL_T, _COL_C2_SUM, _COL_C2_SUM2 = 1, 8, 9
_COL_F_LIQ = 2


def load(p: Path) -> np.ndarray:
    rows = [l for l in p.read_text().splitlines()
            if l.strip() and not l.strip().startswith("i")]
    return np.array([[float(x) for x in l.split()] for l in rows])


def chi_curve(c_sum, c_sum2, f):
    c1 = c_sum / f
    c2 = c_sum2 / f
    return c2 - c1**2


def report(run_id: str) -> None:
    rd = ROOT / "runs" / run_id
    if not (rd / "tr_oxy.dat").exists():
        print(f"{run_id}: no tr_oxy.dat")
        return
    tr = load(rd / "tr_oxy.dat")
    vf = load(rd / "vol_frac_interf.dat")
    t = tr[:, _COL_T]
    c_sum, c_sum2 = tr[:, _COL_C2_SUM], tr[:, _COL_C2_SUM2]

    # align vol_frac to tr_oxy by time (same event, same cadence, but guard)
    f_inst = np.interp(t, vf[:, _COL_T], vf[:, _COL_F_LIQ])
    f_mean = float(vf[:, _COL_F_LIQ].mean())

    nz = np.where(c_sum > 1e-10 * f_mean)[0]
    if len(nz) == 0:
        print(f"{run_id}: tracer never injected")
        return
    t0 = int(nz[0])

    print(f"\n{run_id}  (fidelity {json.loads((rd/'params.json').read_text()).get('fidelity')})")
    print(f"  liquid volume: mean={f_mean:.6g} "
          f"min={f_inst.min():.6g} max={f_inst.max():.6g} "
          f"peak-to-peak={100*(f_inst.max()-f_inst.min())/f_mean:.2f}% of mean")

    for label, f in (("time-averaged f (postprocess)", f_mean),
                     ("instantaneous f", f_inst)):
        s2 = chi_curve(c_sum, c_sum2, f)
        s2max = float(s2[t0])
        chi = np.clip(1.0 - s2 / s2max, 0.0, 1.0)
        tail = chi[t0:]
        # last 20% of the record
        k = int(0.8 * len(tail))
        print(f"  -- {label}")
        print(f"     sigma2_max={s2max:.6f}  chi_max={tail.max():.4f}  "
              f"chi_final={tail[-1]:.4f}")
        print(f"     tail (last 20%): mean={tail[k:].mean():.4f} "
              f"std={tail[k:].std():.4f}  min={tail[k:].min():.4f}")
        for thr in (0.50, 0.75, 0.95):
            idx = int(np.argmax(tail >= thr))
            if tail[idx] < thr:
                print(f"     chi={thr:.2f}: NEVER REACHED")
            else:
                # how many times does it cross back below?
                after = tail[idx:]
                recross = int(np.sum(np.diff((after >= thr).astype(int)) < 0))
                print(f"     chi={thr:.2f}: first at t={t[t0+idx]-t[t0]:.3f} "
                      f"nd, drops back below {recross}x afterwards")


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    for r in sys.argv[1:]:
        report(r)


if __name__ == "__main__":
    main()
