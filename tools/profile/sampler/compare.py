"""Bitwise comparison of every .bin output (and the .dat time series) between two sampler-check runs.

    compare.py RUN_A RUN_B      exit 0 only if every file is byte-identical and both runs wrote snapshots
"""
import filecmp
import re
import sys
from pathlib import Path


def outputs(d: Path) -> list:
    return sorted(str(p.relative_to(d)) for p in d.rglob("*") if p.suffix in (".bin", ".dat") and "build" not in p.parts)


TIMES = re.compile(r"Wall clock time \(s\): \S+ CPU time \(s\): \S+")


def same(fa: Path, fb: Path) -> bool:
    """Byte-identical; logstats.dat is compared without its wall-clock and CPU-time fields."""
    if fa.name == "logstats.dat":
        return TIMES.sub("", fa.read_text()) == TIMES.sub("", fb.read_text())
    return filecmp.cmp(fa, fb, shallow=False)


def main(a: Path, b: Path) -> int:
    fa, fb = outputs(a), outputs(b)
    bad = 0
    if fa != fb:
        print("file sets differ:", sorted(set(fa) ^ set(fb)))
        bad += 1
    nsnap = sum(f.startswith("fields/") for f in fa)
    for f in sorted(set(fa) & set(fb)):
        ok = same(a / f, b / f)
        bad += not ok
        print(("SAME " if ok else "DIFF ") + f)
    print(f"{len(fa)} files, {nsnap} snapshots, {bad} differences")
    return 0 if bad == 0 and nsnap > 0 and "uv_field.bin" in fa and "streaming_field.bin" in fa else 1


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2])))
