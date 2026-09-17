"""The field files are a binary contract between the solver and the analysis.

Nothing in either half declares the layout, so a plane inserted on the C side
without updating SNAP_PLANES silently shifts every field after it: vorticity
reads as the interface, the tracer reads as the geometry, and every figure
built on them is wrong in a way that still looks like a field. These tests pin
the header layout, the plane order, and the fact that the plane COUNT comes
from the file rather than from the reader's own list.
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.fields import (  # noqa: E402
    SNAP_PLANES, cell_centres, liquid_mask, load_snapshot,
    load_streaming_field, load_uv, nearest_snapshot,
)


def _write_snap(path: Path, n: int, t: float, planes: list[np.ndarray]):
    with open(path, "wb") as fh:
        fh.write(struct.pack("<ii", n, len(planes)))
        fh.write(struct.pack("<d", t))
        for p in planes:
            fh.write(p.astype("<f4").tobytes())


def test_snapshot_planes_come_back_in_the_written_order(tmp_path):
    n = 4
    planes = [np.full((n, n), float(k)) for k in range(len(SNAP_PLANES))]
    _write_snap(tmp_path / "snap_000000.bin", n, 1.25, planes)
    s = load_snapshot(tmp_path / "snap_000000.bin")
    assert s["t"] == 1.25 and s["n"] == n
    for k, name in enumerate(SNAP_PLANES):
        assert s[name][0, 0] == float(k)


def test_plane_count_is_read_from_the_file_not_assumed(tmp_path):
    """An older recording has fewer planes. It must load as far as it goes,
    not read past the end of the file into whatever follows."""
    n = 4
    planes = [np.full((n, n), float(k)) for k in range(5)]
    _write_snap(tmp_path / "snap_000000.bin", n, 0.5, planes)
    s = load_snapshot(tmp_path / "snap_000000.bin")
    assert SNAP_PLANES[4] in s
    assert SNAP_PLANES[5] not in s


def test_streaming_field_carries_its_window_length(tmp_path):
    """f_acc is an integral over the averaging window, so a reader cannot tell
    a cell that held water throughout from one the interface swept unless the
    window length travels with it."""
    n = 4
    p = tmp_path / "streaming_field.bin"
    with open(p, "wb") as fh:
        fh.write(struct.pack("<ii", n, 4))
        fh.write(struct.pack("<dd", 9.0, 0.61))
        for k in range(4):
            fh.write(np.full((n, n), float(k), dtype="<f4").tobytes())
    d = load_streaming_field(p)
    assert d["window_dt"] == pytest.approx(0.61)
    assert d["omega_bar"][0, 0] == 2.0 and d["f_acc"][0, 0] == 3.0


def test_uv_field_has_its_own_shorter_header(tmp_path):
    n = 4
    p = tmp_path / "uv_field.bin"
    with open(p, "wb") as fh:
        fh.write(struct.pack("<i", n))
        fh.write(struct.pack("<d", 2.5))
        for k in range(4):
            fh.write(np.full((n, n), float(k), dtype="<f4").tobytes())
    d = load_uv(p)
    assert d["t"] == 2.5 and d["cs"][0, 0] == 3.0


def test_nearest_snapshot_picks_by_time_not_by_filename(tmp_path):
    """Frames are named by index and the instant we want comes from a threshold
    crossing in the time series, so the pick has to read each frame's own t."""
    (tmp_path / "fields").mkdir()
    for idx, t in enumerate((10.0, 12.0, 14.0)):
        _write_snap(tmp_path / "fields" / f"snap_{idx:06d}.bin", 2, t,
                    [np.zeros((2, 2))])
    assert nearest_snapshot(tmp_path, 12.4).name == "snap_000001.bin"
    assert nearest_snapshot(tmp_path, 13.9).name == "snap_000002.bin"


def test_liquid_mask_excludes_water_outside_the_bag():
    """`f` spans the whole lower half-domain; only cs distinguishes the bag."""
    snap = {"f": np.array([[1.0, 1.0]]), "cs": np.array([[1.0, 0.0]])}
    assert liquid_mask(snap).tolist() == [[True, False]]


def test_cell_centres_match_the_writers_own_sampling():
    X, Y = cell_centres(4)
    assert X[0, 0] == pytest.approx(-0.375)
    assert Y[0, 0] == pytest.approx(-0.375)
    assert X[0, 1] == pytest.approx(-0.125)
