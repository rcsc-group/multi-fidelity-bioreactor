"""Oxygen transfer (the curve kLa is fitted from) and mixing against time, reference flux vs rebuilt face field, in each solver.
 ours (blue): L8, 32.5 rpm, theta 7: faithful replay lf_l8_* (projected flux, solid) and the replay with the Kim-style rebuild kim_l8_* (dashed)
 Kim (black): his producer (solid, coupled run) and his offline solver (dashed), N = 128, 20 rpm, theta 7, his default case (agent run, scratch kim_pipeline)
Left: liquid oxygen total, divided by the reference value at 3 cycles of the same solver. Right: chi = 1 - var/var0 of the upper-half tracer c2
(var of the tracer over the liquid, var0 at release). Time = rocking cycles from the release. The two solvers run different flows (rpm, grid).
Usage: uv run python scripts/fig_kla_mix_compare.py
"""
from __future__ import annotations

import csv
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
KIM_CSV = Path("/oscar/scratch/eaguerov/kim_pipeline/out/kim_off_vs_prod_3cycles.csv")
RPM, THETA, NCYC = 32.5, 7.0, 3.0


def ours(rid: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    run = rg.RUNS / rid
    sp = json.loads((run / "params.json").read_text())
    _, T = pp._t_scales(sp)
    t, oxy, c, c2 = pp._joined_cols(run, "tr_oxy.dat", [2, 8, 9])
    _, fl = pp._joined_cols(run, "vol_frac_interf.dat", [pp._COL_F_LIQ])
    fm = fl.mean()
    mean = c / fm
    var = c2 / fm - mean ** 2
    chi = 1 - var / var[0]
    cyc = (t - t[0]) / T
    m = cyc <= NCYC
    return cyc[m], oxy[m] / fm, chi[m]


def kim(src: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = [r for r in csv.DictReader(open(KIM_CSV)) if r["src"] == src]
    cyc = np.array([float(r["cycle"]) for r in rows])
    oxy = np.array([float(r["oxy_liq_prod"]) for r in rows])   # the producer's liquid-oxygen definition, applied to both fields
    var = np.array([float(r["c2_var_liq"]) for r in rows])
    return cyc, oxy, 1 - var / var[0]


plt.rcParams.update(fs.rcparams())
fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.4), sharex=True)
c8 = fs.level_colour(8)
sets = ((c8, ours(f"lf_l8_rpm{RPM:g}_th{THETA:g}"), ours(f"kim_l8_rpm{RPM:g}_th{THETA:g}"), "ours"),
        (fs.KIM, kim("producer"), kim("offline"), "Kim"))
for col, ref, reb, name in sets:
    norm = np.interp(NCYC if ref[0][-1] >= NCYC - 1e-6 else ref[0][-1], ref[0], ref[1])
    for (t, o, chi), ls in ((ref, "-"), (reb, (0, (4, 2)))):
        axs[0].plot(t, o / norm, color=col, lw=2.0, ls=ls)
        axs[1].plot(t, chi, color=col, lw=2.0, ls=ls)
    axs[1].annotate(name, (ref[0][-1], ref[2][-1]), xytext=(5, 0), textcoords="offset points", color=col, va="center", fontsize=10)
axs[0].set_ylabel("liquid oxygen")
axs[1].set_ylabel(r"mixing $\chi$")
for ax in axs:
    ax.set_xlabel("rocking cycles")
    ax.set_xlim(0, NCYC)
    ax.spines[["top", "right"]].set_visible(False)
axs[0].plot([], [], color=fs.NEUTRAL, lw=2.0, label="projected")
axs[0].plot([], [], color=fs.NEUTRAL, lw=2.0, ls=(0, (4, 2)), label="rebuilt")
axs[0].legend(frameon=False, loc="upper left", fontsize=10)
fig.tight_layout()
OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"fig_oxygen_mixing_compare.{ext}", dpi=200)
for col, ref, reb, name in sets:
    i, j = -1, -1
    print(f"{name:5s} at {ref[0][-1]:.2f} cycles: oxygen rebuilt/reference {reb[1][j] / np.interp(reb[0][j], ref[0], ref[1]):.3f}, "
          f"chi reference {ref[2][i]:.3f} rebuilt {np.interp(ref[0][i], reb[0], reb[2]):.3f}")
