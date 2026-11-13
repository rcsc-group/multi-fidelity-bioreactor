"""Pre-registered confound check (diary 2026-10-31): the recording DNS (pr4-rec) must reproduce the database flow
(older binary) before the release: mean free-surface height max and mean tau_95 over cycles 70-80, within 1 %.
Usage: uv run python scripts/analyse/replay_chain_flowcheck.py 8
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts.replay_chain import DB, rid  # noqa: E402
from scripts.test_h1_mfbml import RPMS  # noqa: E402


def window_mean(run: Path, fname: str, col: int, T: float) -> float:
    a = np.loadtxt(run / fname, skiprows=1)
    m = (a[:, 1] >= 70 * T) & (a[:, 1] <= 80 * T)
    return float(a[m, col].mean())


def main(level: int) -> None:
    worst = 0.0
    for rpm in RPMS:
        db = ROOT / "runs" / DB[level].format(rpm)
        rec = rg.RUNS / rid(level, rpm, "rec")
        _, T = pp._t_scales(json.loads((db / "params.json").read_text()))
        row = []
        for name, f, c in (("posY_max", "vol_frac_interf.dat", 4), ("tau_95", "shear_stress.dat", 2)):
            a, b = window_mean(db, f, c, T), window_mean(rec, f, c, T)
            e = 100 * (b / a - 1)
            worst = max(worst, abs(e))
            row.append(f"{name} {e:+.4f} %")
        print(f"{rpm:5g} rpm  " + "  ".join(row))
    print(f"worst {worst:.4f} %: {'PASS' if worst <= 1 else 'FAIL'} (limit 1 %)")


if __name__ == "__main__":
    main(int(sys.argv[1]))
