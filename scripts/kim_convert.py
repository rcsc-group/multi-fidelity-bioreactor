"""Convert solver snapshots (fields/snap_*.bin) to Kim's full-field text dumps.

bio_mix_elvis.m (dev/postprocessing) reads `Data_all_<N>_<t>_tot.txt`, 12
whitespace-separated columns, no header, one row per cell:

    x y ux uy f c cs oxy omega c1 c2 c3

(the layout is inferred from the script's own indexing, textscan at line 98 and
the "tracer: 6,10,11,12" comment; tracer = column 11 = our c2).  Row order is
irrelevant, the script sorts.  Time t in the file name is non-dimensional,
printed with %.9g like MATLAB's num2str(t, 9), on the grid Ini_t:dt_File:Fin_t.
Each target time is served by the snapshot nearest in header time.

Usage: uv run python scripts/kim_convert.py FIELDS_DIR OUT_DIR N INI DT FIN
           [--max-files K]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import fields  # noqa: E402

# snapshot plane for each output column 3..12 (columns 1, 2 are x, y)
COLUMN_PLANES = ["ux", "uy", "f", "c", "cs", "oxy", "omega", "c1", "c2", "c3"]


def kim_grid(ini: float, dt: float, fin: float) -> np.ndarray:
    """MATLAB's materialised [Ini:dt:Fin].

    n = floor((Fin-Ini)/dt) steps; the first half is Ini + k*dt, the second half
    is computed BACKWARDS from the end point (Fin if Ini + n*dt reaches it to
    rounding, else Ini + n*dt).  The two halves differ in the last bits, which
    shows in the 9th digit of the file name (checked against MATLAB for the KP
    grids, tests/test_kim_convert.py).  Note `for i = a:d:b` without brackets
    does NOT do this; bio_mix_elvis.m writes File_t = [Ini_t:dt_File:Fin_t].
    """
    n = int(math.floor((fin - ini) / dt * (1 + 1e-12) + 1e-12))
    k = np.arange(n + 1)
    last = ini + n * dt
    end = fin if abs(last - fin) < 1e-9 * max(1.0, abs(fin)) else last
    return np.where(k <= n / 2, ini + k * dt, end - (n - k) * dt)


def kim_filename(N: int, t: float) -> str:
    return f"Data_all_{N}_{t:.9g}_tot.txt"


def snapshot_rows(snap: dict, x0: float, y0: float, L0: float) -> np.ndarray:
    """(n*n, 12) array; row = j*n + i with j over y (from the bottom)."""
    n = snap["n"]
    dx = L0 / n
    xs = x0 + (np.arange(n) + 0.5) * dx
    ys = y0 + (np.arange(n) + 0.5) * dx
    X, Y = np.meshgrid(xs, ys)          # [j, i]
    cols = [X.ravel(), Y.ravel()] + [snap[p].ravel() for p in COLUMN_PLANES]
    return np.column_stack(cols)


def convert(fields_dir, out_dir, N, ini, dt, fin, x0=-0.5, y0=-0.5, L0=1.0,
            max_files=None):
    """Write one tot.txt per target time; returns a log of dicts."""
    files = sorted(Path(fields_dir).glob("snap_*.bin"))
    times = np.array([_header_t(p)
                      for p in files])
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    grid = kim_grid(ini, dt, fin)
    if max_files is not None:
        grid = grid[:max_files]
    log = []
    for tk in grid:
        k = int(np.argmin(np.abs(times - tk)))
        rows = snapshot_rows(fields.load_snapshot(files[k]), x0, y0, L0)
        name = kim_filename(N, float(tk))
        np.savetxt(out / name, rows, fmt="%.9g")
        log.append({"name": name, "target_t": float(tk),
                    "snap": files[k].name, "snap_t": float(times[k])})
    return log


def _header_t(path) -> float:
    import struct
    with open(path, "rb") as fh:
        _, _, t = struct.unpack("<iid", fh.read(16))
    return t


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("fields_dir")
    ap.add_argument("out_dir")
    ap.add_argument("N", type=int)
    ap.add_argument("ini", type=float)
    ap.add_argument("dt", type=float)
    ap.add_argument("fin", type=float)
    ap.add_argument("--max-files", type=int, default=None)
    a = ap.parse_args(argv)
    log = convert(a.fields_dir, a.out_dir, a.N, a.ini, a.dt, a.fin,
                  max_files=a.max_files)
    errs = [abs(e["snap_t"] - e["target_t"]) for e in log]
    print(f"{len(log)} files; max |snap_t - target_t| = {max(errs):.4g}")


if __name__ == "__main__":
    main()
