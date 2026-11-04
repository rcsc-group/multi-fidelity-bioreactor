"""Unit test: postprocess._edr95_kpis reads edr_pct.dat over the QSS window (2026-11-04).

Synthetic file, exact expected values. Fails before the reader exists, passes after.
"""
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from postprocess import _edr95_kpis  # noqa: E402


def _write(run: Path, rows) -> None:
    lines = ["i t edr_95 edr_max "] + [f"{i} {t:g} {p:g} {m:g} " for i, t, p, m in rows]
    (run / "edr_pct.dat").write_text("\n".join(lines) + "\n")


def test_edr95_is_median_over_qss_window(tmp_path):
    # t = 0.0 .. 0.9 in 0.1 steps; QSS window is 0.25 < t < 0.85 -> t = 0.3 .. 0.8 (6 samples)
    rows = [(i, round(0.1 * i, 6), 0.01 * (i + 1), 0.1 * (i + 1)) for i in range(10)]
    _write(tmp_path, rows)
    out = _edr95_kpis(tmp_path, t_ramp=0.25, t_inject=0.85, ediss_scale=2.0)
    in_win = [r for r in rows if 0.25 < r[1] < 0.85]
    expect_edr95 = float(np.median([r[2] for r in in_win]) * 2.0)
    expect_max = float(max(r[3] for r in in_win) * 2.0)
    assert math.isclose(out["edr_95"], expect_edr95, rel_tol=1e-12)
    assert math.isclose(out["edr_max_qss"], expect_max, rel_tol=1e-12)


def test_edr95_nan_when_file_absent(tmp_path):
    out = _edr95_kpis(tmp_path, t_ramp=0.0, t_inject=None, ediss_scale=1.0)
    assert math.isnan(out["edr_95"]) and math.isnan(out["edr_max_qss"])


def test_edr95_nan_when_window_too_short(tmp_path):
    _write(tmp_path, [(0, 0.0, 1.0, 2.0), (1, 0.1, 1.0, 2.0)])
    out = _edr95_kpis(tmp_path, t_ramp=0.5, t_inject=None, ediss_scale=1.0)
    assert math.isnan(out["edr_95"])
