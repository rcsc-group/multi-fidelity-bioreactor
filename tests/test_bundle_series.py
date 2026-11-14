"""tests/bundle/_series.compare with the x8 bar max(tol, 2 sd): synthetic run folders, no solver."""
import struct

import numpy as np

from tests.bundle import _series as s
from tools.bundle import cases as bc

T = bc.period(bc.RPM_B)


def _run(path, per_period, n_rows=20):
    """A folder with edr_pct.dat whose edr_95 is constant inside each period (values per_period, oldest first)
    and a checkpoint.dump at the end of the last period."""
    path.mkdir()
    rows = ["i t edr_95"]
    for k, v in enumerate(per_period):
        for j in range(n_rows):
            t = (k + (j + 0.5) / n_rows) * T
            rows.append(f"0 {t!r} {v!r}")
    (path / "edr_pct.dat").write_text("\n".join(rows) + "\n")
    (path / "checkpoint.dump").write_bytes(struct.pack("<d", len(per_period) * T) + b"\0" * 8)
    return path


BAR = {"window": "last", "tol": {"edr_mean": 0.05}}


def test_without_sd_periods_the_bar_is_the_fixed_tolerance(tmp_path):
    ref = _run(tmp_path / "ref", [1.0, 1.1, 0.9, 1.0])
    case = _run(tmp_path / "case", [1.0, 1.0, 1.0, 1.06])
    (a, b, d, tol), = s.compare(BAR, ref, case, rpm=bc.RPM_B).values()
    assert tol == 0.05 and abs(d - 0.06 / 1.06) < 1e-9
    assert s.failures(s.compare(BAR, ref, case, rpm=bc.RPM_B))


def test_sd_periods_widens_the_bar_to_two_sd_of_the_reference(tmp_path):
    ref = _run(tmp_path / "ref", [1.0, 1.1, 0.9, 1.0])     # sd (ddof 1) = 0.0816, mean 1.0 -> 2 sd = 16.3 %
    case = _run(tmp_path / "case", [1.0, 1.0, 1.0, 1.06])
    r = s.compare(BAR | {"sd_periods": 4}, ref, case, rpm=bc.RPM_B)
    (_, _, d, tol), = r.values()
    assert abs(tol - 2 * np.std([1.0, 1.1, 0.9, 1.0], ddof=1)) < 1e-9
    assert not s.failures(r)


def test_sd_periods_never_narrows_the_fixed_tolerance(tmp_path):
    ref = _run(tmp_path / "ref", [1.0, 1.0, 1.0, 1.0])      # no scatter: the 5 % floor stays
    case = _run(tmp_path / "case", [1.0, 1.0, 1.0, 1.06])
    r = s.compare(BAR | {"sd_periods": 4}, ref, case, rpm=bc.RPM_B)
    (_, _, _, tol), = r.values()
    assert tol == 0.05
    assert s.failures(r)
