"""L9 warm starts: change of each QoI when the restart ramp starts from the source forcing instead of from rest.

Pairs: c2d_l9_<point> (old rule, ramp from rest) vs c2d_l9_<point>_fix (forcing continued; ab9prev_l9_rpm34.81_th4.26
for 34.81 rpm). Same binary, source run and checkpoint (diary 2026-11-13). dtmix from results.json, else
experiments/multifidelity/dtmix_cache.json.

    uv run python scripts/fig_l9_restart_fix.py OUT.png
"""
import json
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import figstyle as fs  # noqa: E402

QOIS = [("dtmix_0.25", "mixing time\n" r"$\chi=0.25$ (%)"), ("edr_95", "EDR, 95th pct (%)"),
        ("tau_95_qss", r"$\tau$, 95th pct (%)")]
CACHE = json.loads((ROOT / "experiments/multifidelity/dtmix_cache.json").read_text())


def qoi(name, q):
    f = ROOT / "runs" / name / "results.json"
    v = json.loads(f.read_text()).get(q) if f.exists() else None
    if v is None:
        v = (CACHE.get(name) or {}).get(q)
    return float(v) if v is not None and math.isfinite(float(v)) else math.nan


def pairs():
    out = []
    for od in sorted((ROOT / "runs").glob("c2d_l9_rpm*_th*")):
        if od.name.endswith("_fix"):
            continue
        new = "ab9prev_l9_rpm34.81_th4.26" if od.name == "c2d_l9_rpm34.81_th4.26" else od.name + "_fix"
        if not (ROOT / "runs" / new / "results.json").exists():
            continue
        rpm, th = map(float, re.search(r"rpm([\d.]+)_th([\d.]+)", od.name).groups())
        out.append((rpm, th, od.name, new))
    return out


def main(png):
    plt.rcParams.update(fs.rcparams())
    P = pairs()
    fig, axs = plt.subplots(len(QOIS), 1, figsize=(4.2, 6.4), sharex=True)
    c = fs.level_colour(9)
    for ax, (q, lab) in zip(axs, QOIS):
        th = [p[1] for p in P]
        d = [100 * (qoi(p[3], q) / qoi(p[2], q) - 1) for p in P]
        ax.axhline(0, color=fs.GUIDE, lw=1)
        ax.plot(th, d, ls="none", **{k: v for k, v in fs.series_kw(c, fs.LEVEL_MARKER[9], ours=True).items() if k != "ls"})
        ax.set_ylabel(lab)
        fs.qoi_axes(ax, q)
        ax.grid(**fs.GRID_KW)
    axs[-1].set_xlabel(r"rocking amplitude $\theta_{max}$ (deg)")
    fig.tight_layout()
    fig.savefig(png, dpi=200)
    for p in P:
        print(f"{p[2]:28s} " + "  ".join(f"{q} {100*(qoi(p[3], q)/qoi(p[2], q)-1):+6.1f}%" for q, _ in QOIS))


if __name__ == "__main__":
    main(sys.argv[1])
