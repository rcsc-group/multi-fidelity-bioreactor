"""One-off check of the frames_tau binary format before writing the
comparison-video renderer against it.

Usage:
    uv run python scripts/check_frame_format.py <frame.bin>
"""
import struct
import sys

import numpy as np

p = sys.argv[1] if len(sys.argv) > 1 else "runs/smoke_video_test/frames_tau/frame_000010.bin"
with open(p, "rb") as f:
    n, = struct.unpack("i", f.read(4))
    t, = struct.unpack("d", f.read(8))
    Th, = struct.unpack("d", f.read(8))
    xh, = struct.unpack("d", f.read(8))
    fbuf = np.frombuffer(f.read(n * n * 4), dtype="float32")
    tbuf = np.frombuffer(f.read(n * n * 4), dtype="float32")
    ebuf = np.frombuffer(f.read(n * n * 4), dtype="float32")
print("n=", n, "t=", t, "Th=", Th, "xh=", xh)
print("f range", fbuf.min(), fbuf.max())
print("tau range", tbuf.min(), tbuf.max())
mask = fbuf > 0.5
print("mean|tau| over f>0.5:", np.abs(tbuf[mask]).mean() if mask.any() else "n/a")
