"""Restart-history test at L8, 32.5 rpm: is the L10 point (fig9_l10b_seg2) comparable?

The L10 point released its tracer at absolute cycle 80 (Kim's protocol), but
through checkpoint restarts:
  57f68830      restart at cycle 37 with *_prev set, ran to cycle 47
  fig9_l10_seg1 restart at cycle 47 WITHOUT *_prev, release at 47+33 = 80, ran to ~85
  fig9_l10b_seg2 restart at ~85 with *_prev, restart_continue=1, ran to the end
A missing *_prev injects a seam transient (diary 2026-10-01: ~8 cycles), and
the spin-up test (plot_spinup_test.py) shows the release state matters at the
10-25% level. This test repeats that history at L8, where the continuous run
kmix_l8_rpm32.5 and its replicates (release 80/82/85) agree within 1%:
  A  cold start, cycles 0-37, no release
  B  restart at 37 with *_prev, to 47                     (as 57f68830)
  C  restart at 47 WITHOUT *_prev, release at 80, to 85   (as fig9_l10_seg1)
  Cp same as C but WITH *_prev                            (clean warm restart)
  D  restart_continue from C at 85 with *_prev, to 230    (as fig9_l10b_seg2)
  Dp same from Cp
Falsifier of "restart history is harmless": dtmix or kLa of D (or Dp) differ
from the continuous run by more than ~3x the replicate spread (~2%).

Usage: uv run python scripts/submit_restart_test_l8.py --stage A|B|C|Cp|D|Dp [--dry-run]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm            # noqa: E402
from scripts.dump_fields import fields as _dump_fields  # noqa: E402

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-f1c11e0"   # kmix_l8 binary
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
THETA, RPM, LEVEL, NTASKS = 7.0, 32.5, 8, 16
ZEROS = [0.0, 0.0, 0.0]
# stage -> (parent stage, nominal start cycle, ABSOLUTE end cycle, n_mix_cycles, restart_continue,
#           set *_prev). n_mix_cycles "rel80" = release at ABSOLUTE cycle 80 from the actual dump
#           cycle (runs stop at period boundaries, so a dump can land up to one cycle late).
STAGES = {"A": (None, 0, 37, 80, None, None), "B": ("A", 37, 47, 10, 0, True),
          "C": ("B", 47, 85, "rel80", 0, False), "Cp": ("B", 47, 85, "rel80", 0, True),
          "D": ("C", 85, 230, 1, 1, True), "Dp": ("Cp", 85, 230, 1, 1, True)}


def t_nd_per_cycle():
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA)))
    return T_per / (L / (V / (H * 0.5) / T_per))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=sorted(STAGES))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    parent, c0, c1, n_mix, cont, prev = STAGES[a.stage]
    T = t_nd_per_cycle()
    omega = RPM * 2 * math.pi / 60.0
    p = {"run_id": f"rt_l8_{a.stage}", "fidelity": LEVEL, "geometry": GEOMETRY,
         "fill_level": 0.5, "n_harmonics": 1, "theta_max": [THETA, 0.0, 0.0],
         "phi_angular": ZEROS, "omega_b": omega, "omega_h": 0.0, "amplitude_h": ZEROS,
         "phi_horizontal": ZEROS, "n_mix_cycles": n_mix,
         "t_end": round((c1 - c0) * T, 4), "_binary": BINARY}
    ckpt = None
    if parent:
        ckpt = ROOT / "runs" / f"rt_l8_{parent}" / "checkpoint.dump"
        if not ckpt.exists():
            raise SystemExit(f"checkpoint not found: {ckpt}")
        t_ck = _dump_fields(ckpt)[0]["t"]
        c_dump = t_ck / T
        if abs(c_dump - c0) > 1.5:
            raise SystemExit(f"parent dump at cycle {c_dump:.2f}, expected about {c0}")
        if n_mix == "rel80":
            n_mix = int(round(80 - c_dump))
            p["n_mix_cycles"] = n_mix
        p["t_end"] = round((c1 - c_dump) * T, 4)          # run to the ABSOLUTE end cycle
        print(f"  parent dump at absolute cycle {c_dump:.2f}")
        p.update(t_checkpoint=t_ck, restart_continue=cont)
        if prev:
            p.update(omega_b_prev=omega, theta_max_prev=[THETA, 0.0, 0.0],
                     phi_angular_prev=ZEROS, amplitude_h_prev=ZEROS,
                     phi_horizontal_prev=ZEROS, omega_h_prev=0.0)
        if cont == 1:
            p["_parent_run"] = f"rt_l8_{parent}"
    print(f"  {p['run_id']}: cycles {c0}-{c1}, n_mix={n_mix}, continue={cont}, *_prev={prev}")
    if a.dry_run:
        return
    job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime="06:00:00",
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=NTASKS, mem="4G",
                       checkpoint=str(ckpt) if ckpt else None)
    print(f"    job={job}")


if __name__ == "__main__":
    main()
