"""Animated trajectory video: per-cycle peak tau_mean, cold start vs
cross-level warm-start, revealed cycle by cycle. This is the actual
evidence for the two-basin hypothesis shown directly -- two curves that
track together early, then diverge into two distinct, persistent
plateaus -- rather than a spatial field video where the divergence isn't
visually obvious (checked directly: it isn't).

y-axis nondimensionalized as a ratio to the cold start's own converged
value, per this session's standing figure convention. No baked-in title/
caption -- context goes in prose when the file is sent.

Usage:
    uv run python scripts/render_two_basin_trajectory.py <warm_run> <cold_run> <out.mp4>
"""
import argparse
import math
import subprocess
import tempfile
from pathlib import Path

import imageio_ffmpeg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def t_per_nd(rpm=32.5, theta=7.0):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)


def cycle_peak(run):
    d = np.loadtxt(Path("runs") / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, 5]
    T_ND = t_per_nd()
    c = (t - t[0]) / T_ND
    n = int(c[-1])
    return np.array([y[(c >= k) & (c < k + 1)].max() for k in range(n)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("warm_run")
    ap.add_argument("cold_run")
    ap.add_argument("out_mp4")
    ap.add_argument("--fps", type=float, default=6.0)
    args = ap.parse_args()

    warm = cycle_peak(args.warm_run)
    cold = cycle_peak(args.cold_run)
    ref = float(np.mean(cold[-15:]))
    warm, cold = warm / ref, cold / ref
    n = min(len(warm), len(cold))
    warm, cold = warm[:n], cold[:n]

    ymax = max(1.15, 1.05 * max(warm.max(), cold.max()))
    ymin = min(-0.05, 1.05 * min(0, warm.min(), cold.min()))

    plt.rcParams.update({"font.size": 13, "axes.spines.top": False,
                         "axes.spines.right": False})

    with tempfile.TemporaryDirectory() as tmp:
        for k in range(2, n + 1):
            fig, ax = plt.subplots(figsize=(7.2, 4.2))
            x = np.arange(k)
            ax.plot(x, cold[:k], color="#3d3d3d", lw=1.8, label="cold start")
            ax.plot(x, warm[:k], color="#c44536", lw=1.8, label="cross-level warm start")
            ax.axhline(1.0, color="#bbbbbb", lw=0.8, zorder=0)
            ax.set_xlim(0, n)
            ax.set_ylim(ymin, ymax)
            ax.set_xlabel("rocking cycle")
            ax.set_ylabel(r"per-cycle peak $\tau_{\rm mean}$ / cold-start converged value")
            ax.legend(frameon=False, loc="upper right")
            fig.tight_layout()
            fig.savefig(f"{tmp}/f{k:05d}.png", dpi=130)
            plt.close(fig)
        # hold the last frame briefly
        for extra in range(int(args.fps * 2)):
            import shutil
            shutil.copy(f"{tmp}/f{n:05d}.png", f"{tmp}/f{n + 1 + extra:05d}.png")

        ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
        subprocess.run([ffmpeg, "-y", "-r", str(args.fps), "-i", f"{tmp}/f%05d.png",
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", str(args.out_mp4)],
                       check=True, capture_output=True)
    print(f"wrote {args.out_mp4}")


if __name__ == "__main__":
    main()
