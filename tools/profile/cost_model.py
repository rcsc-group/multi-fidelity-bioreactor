"""Cost model per level and rank count from profile runs: wall time and core-h per rocking period.

    cost_model.py RUN_DIR [RUN_DIR ...]   (dirs written by tools/profile/profile.sbatch, named l{L}_n{N}[_tag])

s/step = slope of wall-clock time over steps between logstats rows (from the 2nd row when there are >= 3 rows,
so start-up is excluded; else rows 1-2). Steps per period = T / mean dt over the same rows (T = rocking period
of the profile base case), so a run stopped by its time limit counts too. Diary 2026-11-23 (profiling).
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.bundle import cases as bc  # noqa: E402

ROW = re.compile(r"i: (\d+) t: (\S+) dt: \S+ #Cells: (\d+) Wall clock time \(s\): (\S+)")


def rows(run: Path) -> list:
    return [(int(m[1]), float(m[2]), int(m[3]), float(m[4]))
            for m in map(ROW.search, (run / "logstats.dat").read_text().splitlines()) if m]


def measure(run: Path) -> dict:
    lev, n = map(int, re.match(r"l(\d+)_n(\d+)", run.name).groups())
    r = rows(run)
    a, b = (r[1], r[-1]) if len(r) >= 3 else (r[0], r[-1])
    s_step = (b[3] - a[3]) / (b[0] - a[0])
    steps_period = bc.period(bc.RPM_A) * (b[0] - a[0]) / (b[1] - a[1])
    return {"run": run.name, "level": lev, "ranks": n, "cells_rank": r[-1][2], "s_step": s_step,
            "us_cell_step": s_step / r[-1][2] * 1e6, "steps_period": steps_period,
            "wall_h_period": steps_period * s_step / 3600, "core_h_period": steps_period * s_step / 3600 * n}


def table(runs: list) -> list:
    return sorted((measure(p) for p in runs), key=lambda m: (m["level"], m["ranks"]))


if __name__ == "__main__":
    print(f"{'run':22s} {'L':>2s} {'ranks':>5s} {'cells/rk':>8s} {'us/c-s':>7s} {'s/step':>7s} {'steps/T':>7s} "
          f"{'wall h/T':>8s} {'core-h/T':>8s}")
    for m in table([Path(p) for p in sys.argv[1:]]):
        print(f"{m['run']:22s} {m['level']:2d} {m['ranks']:5d} {m['cells_rank']:8d} {m['us_cell_step']:7.2f} "
              f"{m['s_step']:7.4f} {m['steps_period']:7.0f} {m['wall_h_period']:8.3f} {m['core_h_period']:8.2f}")
