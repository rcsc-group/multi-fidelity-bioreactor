"""Unit tests for scripts/kim_convert.py (snapshots -> Kim's *_tot.txt dumps)."""
import pathlib
import struct
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))
from scripts import kim_convert as kc  # noqa: E402
from scripts import fields  # noqa: E402

NAMES = fields.SNAP_PLANES


def _write(path, n, t, base):
    """Snapshot in the solver's layout; plane k holds base*(k+1) + cell index."""
    with open(path, "wb") as fh:
        fh.write(struct.pack("<iid", n, len(NAMES), t))
        for k, _ in enumerate(NAMES):
            plane = (k + 1) * 1000.0 + np.arange(n * n, dtype=float).reshape(n, n)
            fh.write((plane * base).astype("<f4").tobytes())


def test_columns_order_and_cell_centres(tmp_path):
    n, L0, x0, y0 = 7, 0.3, -0.15, -0.15
    snap_p = tmp_path / "snap_000000.bin"
    _write(snap_p, n, 1.0, 1.0)
    snap = fields.load_snapshot(snap_p)
    rows = kc.snapshot_rows(snap, x0, y0, L0)
    assert rows.shape == (n * n, 12)
    # row index r = j*n + i  (j over y, i over x, as the solver writes)
    dx = L0 / n
    j, i = 3, 5
    r = rows[j * n + i]
    assert r[0] == np.float64(x0 + (i + 0.5) * dx)
    assert r[1] == np.float64(y0 + (j + 0.5) * dx)
    order = ["ux", "uy", "f", "c", "cs", "oxy", "omega", "c1", "c2", "c3"]
    for col, name in enumerate(order, start=2):
        assert r[col] == snap[name][j, i], (col, name)


def test_grid_and_filename():
    g = kc.kim_grid(48.5863329, 1.0633, 48.5863329 + 1.0633 * 3 + 1e-9)
    assert len(g) == 4
    assert abs(g[2] - (48.5863329 + 2 * 1.0633)) < 1e-8
    assert kc.kim_filename(256, 48.5863329) == "Data_all_256_48.5863329_tot.txt"
    assert kc.kim_filename(256, 1.0) == "Data_all_256_1_tot.txt"


def test_convert_nearest_and_text_format(tmp_path):
    n, L0, x0, y0 = 7, 0.3, -0.15, -0.15
    fd = tmp_path / "fields"
    fd.mkdir()
    for k, t in enumerate([10.0, 10.4, 10.9, 11.6]):
        _write(fd / f"snap_{k:06d}.bin", n, t, 1.0 + k)
    out = tmp_path / "out"
    log = kc.convert(fd, out, N=7, ini=10.0, dt=0.8, fin=11.6,
                     x0=x0, y0=y0, L0=L0)
    # targets 10.0, 10.8, 11.6 -> snapshots at 10.0, 10.9 (|0.1| < |0.4|), 11.6
    assert [round(e["snap_t"], 6) for e in log] == [10.0, 10.9, 11.6]
    names = sorted(p.name for p in out.iterdir())
    assert names == sorted(kc.kim_filename(7, t) for t in (10.0, 10.0 + 0.8, 10.0 + 1.6))
    text = (out / kc.kim_filename(7, 10.0 + 0.8)).read_text().splitlines()
    assert len(text) == n * n
    toks = text[0].split()
    assert len(toks) == 12                      # no header, 12 columns
    # x of first cell is -0.15+0.5*0.3/7 = -0.128571428...; >= 9 significant digits
    digits = len(toks[0].lstrip("-").replace(".", "").lstrip("0"))
    assert digits >= 9
    snap = fields.load_snapshot(fd / "snap_000002.bin")
    arr = np.loadtxt(out / kc.kim_filename(7, 10.0 + 0.8))
    assert np.allclose(arr[:, 10], snap["c2"].ravel(), rtol=1e-8)


def test_grid_matches_matlab_materialised_colon():
    # MATLAB's [a:d:b] computes the upper half backwards from the end point, so
    # element 326 differs from a + k*d in the 9th digit (found by running
    # bio_mix_elvis.m: it asked for ..63.0980034.. where a + k*d prints ..035).
    g = kc.kim_grid(48.5863329, 0.044651294, 74.528734714)
    assert len(g) == 582
    assert f"{g[325]:.9g}" == "63.0980034"
    assert f"{g[375]:.9g}" == "65.3305681"
    assert f"{g[0]:.9g}" == "48.5863329" and f"{g[-1]:.9g}" == "74.528734714"[:10]
