"""Falsifiable probe: does a restart segment RE-INJECT the passive tracer?

Hypothesis (from src/BioReactor.c:476, `t_mix = params.t_checkpoint +
T_per_st * params.n_mix_cycles`): on restart t_mix is recomputed relative to
the checkpoint, so `event tracer(t = t_mix)` fires again inside the restart
segment, overwriting the partially-mixed tracer field with a fresh
segregated top-half state -- and resetting chi to 0 mid-experiment.

This is correct behaviour for chain.py's intended use (each segment is a
DIFFERENT rocking condition, warm-started, and wants its own tracer
experiment). It is fatal for Kim Fig. 9, where 9 of the 10 L7 points exceed
the 48 h walltime cap and must be continued across segments as ONE physical
experiment with the tracer injected once.

Prediction if the hypothesis holds: in segment 2, sigma^2 of the tracer
jumps back up to ~0.25 (fully segregated) at t = t_checkpoint + 2 periods,
after having decayed below that in segment 1.
Prediction if it is refuted: sigma^2 in segment 2 continues its monotone
decay with no jump.

Runs at L5 so the whole probe is a couple of minutes.

Usage:  uv run python scripts/probe_tracer_reinjection.py
"""
from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))

BINARY = ROOT / "build" / "BioReactor-mpi"
# mpirun is not on PATH until `module load openmpi`; resolve it once rather
# than relying on the caller's shell having done so.
MPIRUN = (shutil.which("mpirun")
          or "/oscar/rt/9.6/25/spack/x86_64_v3/openmpi-5.0.8-"
             "q6ygsavgbbzcuxp2sog73dhpke6ttgfd/bin/mpirun")
WORK = Path("/oscar/scratch/eaguerov/tracer_reinject_probe")
RPM = 32.5
LEVEL = 5
N_MIX = 2          # inject after 2 cycles
SEG_CYCLES = 4     # each segment spans 4 cycles

GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
THETA = [7.0, 0.0, 0.0]


def t_per_nd() -> float:
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA[0])))
    U = V / (H * 0.5) / T_per
    return T_per / (L / U)


def base_params(**over) -> dict:
    p = {
        "run_id": "probe", "fidelity": LEVEL,
        "geometry": GEOMETRY, "fill_level": 0.5,
        "n_harmonics": 1, "theta_max": THETA,
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": RPM * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "n_mix_cycles": N_MIX,
    }
    p.update(over)
    return p


def run(seg_dir: Path, params: dict, restart: str | None = None) -> None:
    """restart_file is argv[2] (src/BioReactor.c:398); t_checkpoint is read
    back out of the dump itself (line 463), not from params.json."""
    seg_dir.mkdir(parents=True, exist_ok=True)
    (seg_dir / "params.json").write_text(json.dumps(params, indent=2))
    cmd = [MPIRUN, "-n", "4", str(BINARY), "params.json"]
    if restart:
        cmd.append(restart)
    r = subprocess.run(cmd, cwd=seg_dir, capture_output=True, text=True, timeout=3600)
    (seg_dir / "run.out").write_text(r.stdout)
    (seg_dir / "run.err").write_text(r.stderr)
    if not (seg_dir / "tr_oxy.dat").exists():
        raise SystemExit(f"{seg_dir.name}: no tr_oxy.dat\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")


def sigma2(seg_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    """(t, sigma^2) of the VERTICAL_MIXUP tracer, cols 8/9 of tr_oxy.dat."""
    rows = [l for l in (seg_dir / "tr_oxy.dat").read_text().splitlines()
            if l.strip() and not l.strip().startswith("i")]
    a = np.array([[float(x) for x in l.split()] for l in rows])
    vf = [l for l in (seg_dir / "vol_frac_interf.dat").read_text().splitlines()
          if l.strip() and not l.strip().startswith("i")]
    V = np.array([[float(x) for x in l.split()] for l in vf])[:, 2].mean()
    return a[:, 1], a[:, 9] / V - (a[:, 8] / V) ** 2


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)
    Tp = t_per_nd()

    print(f"L{LEVEL}, {RPM} rpm, T_per_nd={Tp:.4f}, inject at {N_MIX} cycles")
    print("\n--- segment 1 (fresh) ---")
    s1 = WORK / "seg1"
    run(s1, base_params(t_end=round(SEG_CYCLES * Tp, 6)))
    t1, v1 = sigma2(s1)
    inj = np.where(v1 > 1e-9)[0]
    print(f"  injection at t={t1[inj[0]]:.4f} (cycle {t1[inj[0]]/Tp:.2f}), "
          f"sigma^2={v1[inj[0]]:.5f}   [analytic 0.25]")
    print(f"  sigma^2 at end of segment 1: {v1[-1]:.5f}")

    ck = s1 / "checkpoint.dump"
    if not ck.exists():
        raise SystemExit(f"segment 1 wrote no checkpoint.dump; found "
                         f"{[p.name for p in s1.iterdir()]}")
    t_ck = float(t1[-1])
    ok = True

    for mode, label in [(0, "restart_continue=0 (NEW EXPERIMENT, chain.py sweeps)"),
                        (1, "restart_continue=1 (CONTINUATION, Fig 9/10/11/12)")]:
        print(f"\n--- segment 2, {label} ---")
        seg = WORK / f"seg2_mode{mode}"
        seg.mkdir(parents=True, exist_ok=True)
        shutil.copy(ck, seg / "restart.dump")
        run(seg, base_params(t_end=round(SEG_CYCLES * Tp, 6), t_checkpoint=t_ck,
                             restart_continue=mode), restart="restart.dump")
        t2, v2 = sigma2(seg)
        d = seg / "restart_diagnostic_post_init.txt"
        c2_post_init = "n/a"
        if d.exists():
            for ln in d.read_text().splitlines():
                if ln.startswith("c2_sum "):
                    c2_post_init = ln.split()[1]
        print(f"  c2_sum after init: {c2_post_init}")
        print(f"  sigma^2 at start of segment 2: {v2[0]:.5f}  "
              f"(segment 1 ended at {v1[-1]:.5f})")
        print(f"  sigma^2 max over segment 2:    {v2.max():.5f} "
              f"at t={t2[np.argmax(v2)]:.4f}")

        reinjected = v2.max() > 1.5 * v1[-1] and v2.max() > 0.15
        preserved = abs(v2[0] - v1[-1]) < 0.25 * max(v1[-1], 1e-12)
        if mode == 0:
            good = reinjected and v2[0] < 1e-9
            print(f"  EXPECT wipe + re-injection (sweep semantics): "
                  f"{'OK' if good else 'FAIL'}")
        else:
            good = preserved and not reinjected
            print(f"  EXPECT tracer preserved, no re-injection: "
                  f"{'OK' if good else 'FAIL'}")
        ok = ok and good

    print("\n--- verdict ---")
    print("  ALL MODES OK" if ok else "  FAILURE -- see above")
    print(f"\nartifacts: {WORK}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
