"""Pre-registered chain test (diary 2026-10-31): the chained mfbml (L8 KRR -> L9 GP -> L10 GP, fig_mf_grid.py chain)
trained with the theta-7 L8 DNS points (kmix_l8_rpm*) replaced by their looped replays (rc8_rpm*_rep). PASS if the
replay-trained L10 mean stays inside the DNS-trained 95 % band at every grid rpm, for each quantity.
Also prints replay/DNS per rpm for the chain quantities. Figure: docs/meetings/prereg/replay_chain_L8.png
Usage: uv run python scripts/analyse/replay_chain_mf.py
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts import fig_mf_grid as g  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts.test_s89_l8_to_l9 import collect  # noqa: E402

OUT = ROOT / "docs/meetings/prereg"
L8_PATS = ["t6_l8_*", "kmix_l8_*", "t5_l8_*", "e8_l8_*", "w_l8_*", "dose_l8_*"]


def l8_points(q: str, replay: bool):
    pts = collect(L8_PATS, q)
    if replay:
        out = []
        for rpm, th, v, name in pts:
            if name.startswith("kmix_l8_rpm"):
                r = json.loads((rg.RUNS / f"rc8_rpm{rpm:g}_rep" / "results.json").read_text())[q]
                print(f"  {q} {rpm:5g} rpm  DNS {v:8.3f}  replay {r:8.3f}  ({100 * (r / v - 1):+5.1f} %)")
                v = r
            out.append((rpm, th, v, name))
        pts = out
    return np.array([[p[0], p[1]] for p in pts], float), np.array([p[2] for p in pts], float)


def main() -> None:
    warnings.filterwarnings("ignore")
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for a, (q, lab) in zip(ax, (("kLa_1T_10", r"$k_La$, $C^*\!=\!0.10$ (1/h)"),
                                ("dtmix_0.50", r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$ (s)"))):
        _, a9, a10, _ = g.scalar_data(q)
        res = {}
        for kind in ("DNS", "replay"):
            a8 = l8_points(q, kind == "replay")
            res[kind] = (a8, g.models("chain", a8, a9, a10)(7.0))
        (_, (m_d, lo, hi), _), (_, (m_r, _, _), _) = res["DNS"][1], res["replay"][1]
        out = (m_r < lo) | (m_r > hi)
        print(f"{q}: replay-trained L10 mean outside the DNS 95 % band at {out.sum()} of {len(g.GRID)} grid rpm "
              f"-> {'PASS' if not out.any() else 'FAIL'}; max |mean shift| {100 * np.max(np.abs(m_r / m_d - 1)):.1f} %")
        a.fill_between(g.GRID, lo, hi, color="0.85", lw=0)
        a.plot(g.GRID, m_d, color="0.3", lw=3.0, label="L10 prediction, L8 from DNS")
        a.plot(g.GRID, m_r, color=fs.KIM if hasattr(fs, "KIM") else "k", lw=1.3, ls="--",
               label="L10 prediction, L8 from replay")
        for kind, mk in (("DNS", "o"), ("replay", "x")):
            X8, y8 = res[kind][0]
            k = X8[:, 1] == 7.0
            a.plot(X8[k, 0], y8[k], mk, color="tab:blue", ms=6, mfc="none" if mk == "o" else None,
                   label=f"L8 {kind}")
        X10, y10 = a10
        k = X10[:, 1] == 7.0
        a.plot(X10[k, 0], y10[k], "s", color="tab:pink", ms=6, label="L10")
        a.set_xlabel("rocking speed (rpm)"); a.set_ylabel(lab); a.grid(**fs.GRID_KW)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"replay_chain_L8.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "replay_chain_L8.png")


if __name__ == "__main__":
    main()
