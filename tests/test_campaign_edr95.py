"""Campaign 2D (2026-11-05): edr_95 for the fit = median over the last 3 whole cycles of the run's edr_pct.dat,
scaled to W/m^3 (same scale as the hydro QoIs). Synthetic run folder only.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))
sys.path.insert(0, str(ROOT))

import campaign as c  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402


def _folder(root: Path, name: str, rpm: float = 20.0):
    d = root / name
    d.mkdir(parents=True)
    p = {"fidelity": 8, "omega_b": rpm * 2 * math.pi / 60.0, "theta_max": [4.0, 0.0, 0.0],
         "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "t_checkpoint": 0.0, "n_mix_cycles": 3}
    (d / "params.json").write_text(json.dumps(p))
    return d, p


def test_edr95_is_median_of_last_three_cycles_scaled(tmp_path):
    d, p = _folder(tmp_path, "run_a")
    _, T = _t_scales(p)
    a = p["geometry"]["a"]
    T_bio, _ = _t_scales(p)
    escale = 1000.0 * (a / T_bio) ** 3 / a
    t = np.arange(0.0, 10.0 * T, 0.02 * T)
    v = np.where(t >= 7.0 * T, 2.0, 5.0)              # last 3 cycles (t >= 7T of 10T) hold 2.0, earlier 5.0
    with open(d / "edr_pct.dat", "w") as f:
        f.write("i t edr_95 edr_max\n")
        for i, (tt, vv) in enumerate(zip(t, v)):
            f.write(f"{i} {tt:.10g} {vv:.10g} {vv * 2:.10g}\n")
    got = c.edr95_late(d)
    assert got == pytest_approx(2.0 * escale)


def test_edr95_missing_file_gives_none(tmp_path):
    d, _ = _folder(tmp_path, "run_b")
    assert c.edr95_late(d) is None


def pytest_approx(x):
    import pytest
    return pytest.approx(x, rel=1e-9)
