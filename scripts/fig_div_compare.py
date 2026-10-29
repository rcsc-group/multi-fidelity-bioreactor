"""How divergence-free is the face velocity that advects the scalars? Same metric in both solvers, per step over one rocking cycle:
    max over fully open cells (cs = 1, four open faces) of |(uf.x[1]-uf.x[]) + (uf.y[0,1]-uf.y[])|  /  max |uf|
 - ours, true flux      : the face flux recorded from the Navier-Stokes projection (replay_diag, src/BioReactor.c), L8, 32.5 rpm, theta 7
 - ours, Kim-style      : fm*face_value(u.x,0) of the same coupled run in the same record step (REC_KIM_FLUX), same replay_diag
 - Kim's own solver     : codes/TracerTest_test.c line 815, his producer at N = 128 (divdiag.dat of the agent run in /oscar/scratch/eaguerov/kim_pipeline)
Colour = whose data (figstyle: Kim black, our L8 blue). Usage: uv run python scripts/fig_div_compare.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts import postprocess as pp  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts import replay_lf as rl  # noqa: E402

OUT = ROOT / "docs/meetings/2026-10-29"
KIM_DIV = Path("/oscar/scratch/eaguerov/kim_pipeline/out/off_real_c1/divdiag.dat")
KIM_T = 0.6074      # his non-dimensional rocking period (T_CYCLE in TracerTest_test.c)
RPM, THETA = 32.5, 7.0


def ours(tag: str) -> tuple[np.ndarray, np.ndarray]:
    sp = json.loads((rl.SRC_RUNS / rl.source_id(RPM, THETA) / "params.json").read_text())
    _, T = pp._t_scales(sp)
    a = np.loadtxt(rg.scratch_dir(f"divdiag_{tag}_l8_rpm{RPM:g}_th{THETA:g}") / "replay_diag.txt")
    return (a[:, 1] - a[0, 1]) / T, a[:, 2] / a[:, 4]


plt.rcParams.update(fs.rcparams())
fig, ax = plt.subplots(figsize=(6.4, 3.8))
tt, yt = ours("true")
tk, yk = ours("kim")
b = np.loadtxt(KIM_DIV, skiprows=1)
tb, yb = b[:, 1] / KIM_T, b[:, 2] / b[:, 3]
c8 = fs.level_colour(8)
ax.semilogy(tt, yt, color=c8, lw=2.0)
ax.semilogy(tk, yk, color=c8, lw=2.0, ls=(0, (4, 2)))
ax.semilogy(tb, yb, color=fs.KIM, lw=2.0)
for x, y, s, col in ((tt, yt, "ours", c8), (tk, yk, "rebuilt", c8), (tb, yb, "Kim", fs.KIM)):
    ax.annotate(s, (x[-1], np.median(y)), xytext=(6, 0), textcoords="offset points", color=col, va="center", fontsize=10)
ax.set_xlim(0, 1)
ax.set_xlabel("rocking cycles")
ax.set_ylabel(r"$\max\,|\nabla\!\cdot u_f|\,\Delta\,/\,\max|u_f|$")
ax.spines[["top", "right"]].set_visible(False)
fig.subplots_adjust(right=0.84)
OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"fig_divergence_compare.{ext}", dpi=200)
for n, y in (("ours true", yt), ("ours rebuilt", yk), ("Kim", yb)):
    print(f"{n:13s} median {np.median(y):.2e}  max {y.max():.2e}  steps {len(y)}")
