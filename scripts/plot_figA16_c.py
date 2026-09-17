"""Kim et al. (2024) Fig. A.16(c) -- L2 norm of the velocity error vs resolution.

Kim's definition (Main.tex, appendix `app:grid`):

    e_{i,L2} = sqrt( (1/N_e) * sum_e (u_i' - u_i,ref')^2 ) / <u_i,ref'>

with i = x, y and the reference taken at n_L = 2^11.

Two departures from his panel, both forced and both stated rather than
absorbed:

  * THE REFERENCE IS OUR FINEST, not 2^11. One level above Kim's production
    mesh is ~11x his per-cycle cost (cost_model.PER_LEVEL_WORK_FACTOR); the
    reference here is n_L = 2^10, so the curve reads as convergence toward
    our finest grid rather than toward his. The shape of the convergence is
    the statement; its absolute floor is set by the reference and is not
    comparable to his.

  * THE INSTANT IS A PERIOD BOUNDARY, not his t/T_p = 29.77 (maximum rocking
    angle). uv_field.bin is written at the checkpoint, which the solver pins
    to a zero-crossing by construction. A convergence measure needs the SAME
    phase across levels, which that guarantees; an arbitrary mid-cycle
    instant would need a solver change and buys nothing here.

The reference field is block-averaged down to each coarse grid before
differencing -- both grids are uniform and the refinement ratio is a power of
two, so the coarse cell's value is exactly the mean of the fine cells inside
it. Interpolating the coarse field UP instead would compare the interpolant's
smoothness rather than the solution's error.

Only liquid inside the bag enters the norm (f > 0.5 and cs > 0.5 on both
grids), for the same reason every other KPI in this project is bag-masked:
the region outside carries no fluid and dilutes any mean taken over it.

Usage:  uv run python scripts/plot_figA16_c.py
"""
from __future__ import annotations

import math
import struct
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
RUNS = ROOT / "runs"
OUT = ROOT / "experiments/kimetal2024/figure_replicas/replicated_FigA16_c.png"

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

# level -> run_id, written by scripts/submit_figA16_c.py
SERIES = {6: "figA16c_l6", 7: "figA16c_l7", 8: "figA16c_l8",
          9: "figA16c_l9", 10: "figA16c_l10"}
REF_LEVEL = 10


def load_uv(path: Path):
    """(t, ux, uy, f, cs) from uv_field.bin -- four n x n float32 planes."""
    with open(path, "rb") as fh:
        (n,) = struct.unpack("i", fh.read(4))
        (t,) = struct.unpack("d", fh.read(8))
        planes = [np.frombuffer(fh.read(n * n * 4), dtype=np.float32)
                  .reshape(n, n).astype(float) for _ in range(4)]
    return (t, *planes)


def block_mean(a: np.ndarray, factor: int) -> np.ndarray:
    """Exact coarse-cell average of a uniform field, factor a power of two."""
    n = a.shape[0] // factor
    return a.reshape(n, factor, n, factor).mean(axis=(1, 3))


def l2_error(coarse, ref, factor: int) -> tuple[float, float]:
    """(e_x, e_y) of `coarse` against `ref` block-averaged onto its grid."""
    t_c, ux_c, uy_c, f_c, cs_c = coarse
    t_r, ux_r, uy_r, f_r, cs_r = ref
    ux_r, uy_r = block_mean(ux_r, factor), block_mean(uy_r, factor)
    f_r, cs_r = block_mean(f_r, factor), block_mean(cs_r, factor)
    m = (f_c > 0.5) & (cs_c > 0.5) & (f_r > 0.5) & (cs_r > 0.5)
    if not m.any():
        return math.nan, math.nan
    out = []
    for a, b in ((ux_c, ux_r), (uy_c, uy_r)):
        denom = abs(b[m].mean())
        rms = math.sqrt(float(((a[m] - b[m]) ** 2).mean()))
        out.append(rms / denom if denom > 0 else math.nan)
    return out[0], out[1]


def main() -> None:
    fields = {}
    for level, run in SERIES.items():
        p = RUNS / run / "uv_field.bin"
        if p.exists():
            fields[level] = load_uv(p)
    if REF_LEVEL not in fields:
        raise SystemExit(
            f"no reference field: runs/{SERIES[REF_LEVEL]}/uv_field.bin is "
            f"missing. Submit the set with scripts/submit_figA16_c.py first; "
            f"levels present: {sorted(fields)}")

    ref = fields[REF_LEVEL]
    times = {lv: d[0] for lv, d in fields.items()}
    if max(times.values()) - min(times.values()) > 1e-6:
        raise SystemExit(
            f"fields are not at a common instant -- an L2 norm across levels "
            f"at different phases measures the phase, not the error: {times}")

    rows = []
    for level in sorted(fields):
        if level == REF_LEVEL:
            continue
        ex, ey = l2_error(fields[level], ref, 2 ** (REF_LEVEL - level))
        rows.append((level, ex, ey))
        print(f"  n_L=2^{level}   e_x={ex:.4f}   e_y={ey:.4f}")
    if not rows:
        raise SystemExit("only the reference level is present; nothing to plot")

    n_l = [2 ** r[0] for r in rows]
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    ax.plot(n_l, [r[1] for r in rows], color="crimson", marker="o", ms=6,
            lw=1.3, label=r"$e_{x,L_2}$")
    ax.plot(n_l, [r[2] for r in rows], color="mediumorchid", marker="s", ms=6,
            lw=1.3, label=r"$e_{y,L_2}$")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel(r"$n_L$", fontsize=13)
    ax.set_ylabel(r"$e_{i,L_2}$", fontsize=13)
    ax.tick_params(which="both", direction="in")
    ax.grid(True, which="major", ls=":", alpha=0.4)
    ax.legend(fontsize=10, frameon=False, loc="upper left",
              bbox_to_anchor=(1.02, 1.0))
    ax.text(-0.18, 1.02, r"$(c)$", transform=ax.transAxes, fontsize=15,
            style="italic")
    fig.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
