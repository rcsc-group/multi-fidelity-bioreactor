"""List the field names stored in a Basilisk checkpoint.dump.

Format (basilisk/src/output.h:1018-1074): struct DumpHeader {double t; long
len; int i, depth, npe, version; coord n;} followed by `len` entries of
{unsigned name_len; char name[name_len]} and then {X0,Y0,Z0,L0} as doubles.

Useful for answering "was this field ever written?" without guessing -- e.g.
whether the passive tracers and oxygen survive a checkpoint, which decides
whether a mixing or kLa run can be continued across a walltime boundary.

Usage:  uv run python scripts/dump_fields.py <path/to/checkpoint.dump>
"""
from __future__ import annotations

import struct
import sys
from pathlib import Path


def fields(path: Path) -> tuple[dict, list[str]]:
    with open(path, "rb") as fh:
        t, n_fields, i, depth, npe, version = struct.unpack("=dqiiii", fh.read(8 + 8 + 16))  # q: C long is 8 bytes on LP64
        dims = struct.unpack("=3d", fh.read(24))
        names = []
        for _ in range(n_fields):
            (ln,) = struct.unpack("=I", fh.read(4))
            names.append(fh.read(ln).decode("ascii", "replace"))
    meta = {"t": t, "n_fields": n_fields, "iter": i, "depth": depth,
            "npe": npe, "version": version, "n": dims}
    return meta, names


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    path = Path(sys.argv[1])
    meta, names = fields(path)
    print(f"{path}")
    print(f"  t={meta['t']:.6f}  iter={meta['iter']}  depth={meta['depth']}  "
          f"npe={meta['npe']}  version={meta['version']}  n_fields={meta['n_fields']}")
    print(f"  fields: {', '.join(names)}")
    for want in ("c", "c1", "c2", "c3", "oxy"):
        print(f"  {want:>4} present: {want in names}")


if __name__ == "__main__":
    main()
