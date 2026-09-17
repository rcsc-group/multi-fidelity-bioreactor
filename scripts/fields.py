"""Readers for the solver's uniform-grid field files.

Three formats, all little-endian, all written by BioReactor.c:

``fields/snap_%06d.bin``  (write_snapshot)
    ``int n; int nplanes; double t;`` then ``nplanes`` planes of ``n*n``
    float32, in the order ux, uy, omega, f, cs, c, c1, c2, c3, oxy.

``streaming_field.bin``  (event streaming_vorticity)
    ``int n; int nplanes; double t; double window_dt;`` then ubar_x, ubar_y,
    omega_bar, f_acc.

``uv_field.bin``  (write_uv_field, at the checkpoint instant)
    ``int n; double t;`` then ux, uy, f, cs.

The plane ORDER is part of the format and the plane COUNT is in the header;
read the count from the file rather than assuming it, because the snapshot
format gained planes once already and older recordings will not have them.

Every reader returns planes as ``float64`` arrays indexed ``[j, i]`` with j
running over y from the bottom of the domain, matching the write loop.
"""
from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

# Plane order of fields/snap_*.bin, as written.
SNAP_PLANES = ["ux", "uy", "omega", "f", "cs", "c", "c1", "c2", "c3", "oxy"]
STREAM_PLANES = ["ubar_x", "ubar_y", "omega_bar", "f_acc"]
# Tracer key -> Kim's name for that initial configuration (Figs. 5 and 6).
TRACER_CONFIG = {"c2": "top half", "c1": "left half",
                 "c3": "circle", "c": "line"}


def _planes(fh, n: int, count: int) -> list[np.ndarray]:
    return [np.frombuffer(fh.read(n * n * 4), dtype="<f4")
            .reshape(n, n).astype(float) for _ in range(count)]


def load_snapshot(path: Path) -> dict:
    """{'t', 'n', and one entry per plane} from fields/snap_*.bin."""
    with open(path, "rb") as fh:
        n, nplanes = struct.unpack("<ii", fh.read(8))
        (t,) = struct.unpack("<d", fh.read(8))
        data = _planes(fh, n, nplanes)
    out = {"t": t, "n": n}
    out.update(dict(zip(SNAP_PLANES[:nplanes], data)))
    return out


def load_streaming_field(path: Path) -> dict:
    """{'t', 'n', 'window_dt', 'ubar_x', 'ubar_y', 'omega_bar', 'f_acc'}."""
    with open(path, "rb") as fh:
        n, nplanes = struct.unpack("<ii", fh.read(8))
        t, window_dt = struct.unpack("<dd", fh.read(16))
        data = _planes(fh, n, nplanes)
    out = {"t": t, "n": n, "window_dt": window_dt}
    out.update(dict(zip(STREAM_PLANES[:nplanes], data)))
    return out


def load_uv(path: Path) -> dict:
    """{'t', 'n', 'ux', 'uy', 'f', 'cs'} from uv_field.bin."""
    with open(path, "rb") as fh:
        (n,) = struct.unpack("<i", fh.read(4))
        (t,) = struct.unpack("<d", fh.read(8))
        ux, uy, f, cs = _planes(fh, n, 4)
    return {"t": t, "n": n, "ux": ux, "uy": uy, "f": f, "cs": cs}


def snapshots(run_dir: Path) -> list[Path]:
    return sorted((Path(run_dir) / "fields").glob("snap_*.bin"))


def nearest_snapshot(run_dir: Path, t_target: float) -> Path | None:
    """The recorded frame closest in time to `t_target` (non-dimensional).

    Kim's field figures are at instants defined by a THRESHOLD -- chi = 0.5,
    C* = 0.10/0.25/0.50 -- which the solver does not know in advance. The
    thresholds are crossed in the time series, and the frame is picked here.
    """
    paths = snapshots(run_dir)
    if not paths:
        return None
    times = np.array([load_snapshot(p)["t"] for p in paths])
    return paths[int(np.argmin(np.abs(times - t_target)))]


def cell_centres(n: int, x0: float = -0.5, length: float = 1.0):
    """(X, Y) meshgrid of cell centres, matching the writer's own loop."""
    a = x0 + (np.arange(n) + 0.5) * (length / n)
    return np.meshgrid(a, a)


def liquid_mask(snap: dict, f_key: str = "f") -> np.ndarray:
    """Water inside the embedded bag -- `f` alone spans the whole lower half."""
    return (snap[f_key] > 0.5) & (snap["cs"] > 0.5)
