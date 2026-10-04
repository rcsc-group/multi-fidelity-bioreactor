"""A clean L10 point at 32.5 rpm: cold start, Kim's protocol, release at cycle 80.

Why. The only L10 mixing point (fig9_l10b_seg2) is contaminated: its chain had
restarts without *_prev (at cycles ~11, ~22 and 47). The L8 restart-history test
(scripts/submit_restart_test_l8.py, diary 2026-10-03) showed that a seam without
*_prev at cycle 47 shifts dtmix by +28..57% and kLa by up to -39% at release 80,
while the same history WITH *_prev reproduces the continuous run within 2%.
So: a continuous cold start; any later segment is a same-condition
continuation with *_prev set and restart_continue=1 (the validated case).

Segment 1 runs as far as the 96 h cap allows (graceful checkpoint before it).
Later segments continue until chi reaches 0.95 (check results.json first).

Usage:
    uv run python scripts/submit_l10_cold.py --segment 1
    uv run python scripts/submit_l10_cold.py --segment 2 --from-run l10c_rpm32.5_seg1
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm               # noqa: E402
from scripts.dump_fields import fields as _dump_fields  # noqa: E402
from scripts.cost_model import min_per_cycle            # noqa: E402

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-f1c11e0"
RPM, THETA, NTASKS, WALL_H = 32.5, 7.0, 64, 96
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
ZEROS = [0.0, 0.0, 0.0]
SEG_CYCLES = 260          # nominal; the graceful checkpoint stops earlier at the wall-time cap


def t_nd_per_cycle():
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA)))
    return T_per / (L / (V / (H * 0.5) / T_per))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--segment", type=int, required=True)
    ap.add_argument("--from-run")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if not Path(BINARY).exists():
        sys.exit(f"binary missing: {BINARY}")
    T = t_nd_per_cycle()
    omega = RPM * 2 * math.pi / 60.0
    mc, conf = min_per_cycle(10, NTASKS, RPM)
    p = {"run_id": f"l10c_rpm{RPM:g}_seg{a.segment}", "fidelity": 10, "geometry": GEOMETRY,
         "fill_level": 0.5, "n_harmonics": 1, "theta_max": [THETA, 0.0, 0.0],
         "phi_angular": ZEROS, "omega_b": omega, "omega_h": 0.0, "amplitude_h": ZEROS,
         "phi_horizontal": ZEROS, "t_end": round(SEG_CYCLES * T, 4), "_binary": BINARY}
    ckpt = None
    if a.segment == 1:
        p["n_mix_cycles"] = 80                      # release at absolute cycle 80
    else:
        if not a.from_run:
            sys.exit("--from-run required for segments > 1")
        ckpt = ROOT / "runs" / a.from_run / "checkpoint.dump"
        if not ckpt.exists():
            sys.exit(f"checkpoint not found: {ckpt}")
        t_ck = _dump_fields(ckpt)[0]["t"]
        p.update(t_checkpoint=t_ck, restart_continue=1, n_mix_cycles=1, _parent_run=a.from_run,
                 omega_b_prev=omega, theta_max_prev=[THETA, 0.0, 0.0], phi_angular_prev=ZEROS,
                 amplitude_h_prev=ZEROS, phi_horizontal_prev=ZEROS, omega_h_prev=0.0)
        print(f"  continuation from {a.from_run} at absolute cycle {t_ck / T:.2f}")
    print(f"  {p['run_id']}: {mc:.1f} min/cycle ({conf}) -> about {WALL_H * 60 / mc:.0f} cycles "
          f"in {WALL_H} h at {NTASKS} ranks, {mc / 60 * NTASKS:.1f} core-h per cycle")
    if a.dry_run:
        return
    job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=f"{WALL_H}:00:00",
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=NTASKS, mem="4G", checkpoint=str(ckpt) if ckpt else None)
    print(f"    job={job}")


if __name__ == "__main__":
    main()
