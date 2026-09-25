"""L7 angle sweep at 32.5 rpm: the LF side of the Fig 13(b) multi-fidelity demo.

Every existing L7 angle run predates the bag-mask fix (2026-09-15), so its
spatial means are diluted ~3.5x, and none logs Kim's exact ediss/tau columns
(ediss_kim_max, ...). The HF side, fig10_l9_th*, was run on the lean binary
at f1c11e0 and has both. Same binary here, so LF and HF differ in level only.

Cold start, 80 cycles (Kim's spin-up), no tracer or oxygen (n_mix_cycles past
t_end). At (7, 8) 0.163 min/cycle measured at 32.5 rpm: ~13 min per angle.

Usage:
    uv run python scripts/submit_mf_l7_angle.py --dry-run
    uv run python scripts/submit_mf_l7_angle.py
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
RPM, LEVEL, NTASKS, CYCLES = 32.5, 7, 8, 80
ANGLES = [7.0, 6.0, 5.0, 4.0, 3.0, 2.0]
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}


def t_scales(theta: float) -> tuple[float, float]:
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T_per, L / (V / (H * 0.5) / T_per)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for th in ANGLES:
        T_per, T_bio = t_scales(th)
        params = {
            "run_id": f"mf_l7_th{th:g}", "fidelity": LEVEL,
            "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
            "theta_max": [th, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
            "omega_b": RPM * 2 * math.pi / 60.0, "omega_h": 0.0,
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
