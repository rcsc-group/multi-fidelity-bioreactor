"""L7 rpm sweep at 7 deg: the LF side of the Fig 13(a) multi-fidelity demo.

Kim's Fig 13 means are the SIGNED tau / EDR on his fully-liquid, uncut
mask (tau_kim_mean, ediss_kim_mean), which only the lean binary at f1c11e0
logs. The L8/L10 rpm runs predate it; the HF side here, fig9_l9_rpm*, has it.
Same binary, so LF and HF differ in level only, as in the angle sweep
(submit_mf_l7_angle.py).

Cold start, 80 cycles (Kim's spin-up), no tracer or oxygen.

Usage:
    uv run python scripts/submit_mf_l7_rpm.py --dry-run
    uv run python scripts/submit_mf_l7_rpm.py
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm    # noqa: E402

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-f1c11e0"
THETA, LEVEL, NTASKS, CYCLES = 7.0, 7, 8, 80
RPMS = [17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}


def t_scales(rpm: float, theta: float = THETA) -> tuple[float, float]:
    omega_b = rpm * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T_per, L / (V / (H * 0.5) / T_per)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for rpm in RPMS:
        T_per, T_bio = t_scales(rpm)
        params = {
            "run_id": f"mf_l7_rpm{rpm:g}", "fidelity": LEVEL,
            "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
            "theta_max": [THETA, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
            "omega_b": rpm * 2 * math.pi / 60.0, "omega_h": 0.0,
            "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
            "t_end": round(CYCLES * T_per / T_bio, 4),
            "n_mix_cycles": CYCLES + 10,
            "_binary": BINARY,
        }
        print(f"  {params['run_id']}: {CYCLES} cycles, t_end={params['t_end']}")
        if a.dry_run:
            continue
        job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                           walltime="01:00:00",
                           template=ROOT / "config" / "slurm_mpi_template.sh",
                           cpus=1, ntasks=NTASKS, mem="4G")
        print(f"    submitted job={job}")


if __name__ == "__main__":
    main()
