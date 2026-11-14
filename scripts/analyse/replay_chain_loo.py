"""The established mfbml protocols, run twice: L8 from the database DNS, and L8 (kmix_l8_rpm*, theta 7) from the
looped replays rc8_rpm*_rep (diary 2026-10-31). Nothing else changes (L9, L10, models, leave-out rules, baselines).
  step L8 -> L9   scripts/experiments/test_s89_l8_to_l9.py: leave-one-rpm-group-out over the theta-7 L9 runs; kLa_1T_10, dtmix_0.50
  chain L8 -> L9 -> L10   scripts/experiments/test_chain_l10_theta7.py: leave-one-out over the 7 L10 theta-7 runs; kLa_1T_10
     (variant A: L9 at the held-out rpm; C: L9 within 0.25 rpm removed). L10 mixing time exists at 3 runs only: no chain test.
Figure: internal/docs/meetings/prereg/replay_chain_loo.png. JSON: experiments/multifidelity/replay_chain_loo.json
Usage: uv run python scripts/analyse/replay_chain_loo.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from scripts import figstyle as fs  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts.experiments import test_chain_l10_theta7 as tc  # noqa: E402
from scripts.experiments import test_s89_l8_to_l9 as s89  # noqa: E402
from scripts.hypotheses.test_h1_mfbml import RPMS  # noqa: E402

EXP = ROOT / "experiments/multifidelity"
OUT = ROOT / "internal/docs/meetings/prereg"


def swap(q: str) -> dict:
    return {f"kmix_l8_rpm{r:g}": json.loads((rg.RUNS / f"rc8_rpm{r:g}_rep" / "results.json").read_text())[q] for r in RPMS}


def main() -> None:
    res = {}
    for kind in ("dns", "replay"):
        print(f"===== L8 from {kind}")
        sw = {q: swap(q) for q in s89.QS} if kind == "replay" else None
        step = s89.main(l8_swap=sw, out_path=EXP / f"replay_chain_loo_s89_{kind}.json")
        chain = tc.main(q="kLa_1T_10", l8_swap=sw["kLa_1T_10"] if sw else None,
                        out_path=EXP / f"replay_chain_loo_chain_{kind}.json")
        res[kind] = dict(step=step, chain=chain)
    summary = {k: dict(step={q: v["step"][q]["rms"] | dict(cover=v["step"][q]["cover"], verdict=v["step"][q]["verdict"])
                             for q in s89.QS},
                       chain={a: dict(rms=v["chain"][a]["rms"], cover=v["chain"][a]["cover"]) for a in "AC"}
                       | dict(baselines=v["chain"]["baselines"], verdict=v["chain"]["verdict"]))
               for k, v in res.items()}
    (EXP / "replay_chain_loo.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    figure(res)


def figure(res: dict) -> None:
    """Held-out predictions against the observed value, DNS-trained (thick, hollow) vs replay-trained (thin, filled)."""
    plt.rcParams.update(fs.rcparams())
    panels = [("step", "kLa_1T_10", r"$k_La$, L8$\to$L9 (1/h)", fs.COL[9] if hasattr(fs, "COL") else "tab:green"),
              ("step", "dtmix_0.50", r"$t_{\mathrm{mix}}$, L8$\to$L9 (s)", fs.COL[9] if hasattr(fs, "COL") else "tab:green"),
              ("chain", "A", r"$k_La$, L8$\to$L9$\to$L10 (1/h)", fs.COL[10] if hasattr(fs, "COL") else "tab:pink")]
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.9))
    for a, (kind, key, lab, col) in zip(ax, panels):
        for src, mfc, lw, ms, off in (("dns", "none", 2.4, 8, -0.15), ("replay", col, 1.0, 5, 0.15)):
            rows = res[src][kind][key]["rows"]
            x = np.array([r["rpm"] for r in rows]) + off
            o = np.exp([r["obs"] for r in rows])
            m = np.exp([r["mf" if kind == "step" else "mean"] for r in rows])
            lo, hi = np.exp([r["lo"] for r in rows]), np.exp([r["hi"] for r in rows])
            a.errorbar(x, m, yerr=[m - lo, hi - m], fmt="o", color=col, mfc=mfc, ms=ms, lw=lw / 2, capsize=0,
                       label=f"prediction, L8 from {'DNS' if src == 'dns' else 'replay'}")
        a.plot(x - 0.15, o, "_", color="k", ms=14, mew=2.2, label="held-out run")
        a.set_xlabel("rocking speed (rpm)"); a.set_ylabel(lab); a.grid(**fs.GRID_KW)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"replay_chain_loo.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "replay_chain_loo.png")


if __name__ == "__main__":
    main()
