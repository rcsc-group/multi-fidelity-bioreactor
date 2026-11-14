"""Feature ledger: error each mechanism adds in production, from its clean pairs (tools/ledger, feature_ledger.json).

One row per mechanism with clean pairs (confounded rows have no number to show), one panel per quantity.
y = with / without - 1 in %, one point per clean pair; colour and marker = mesh level of the pair (figstyle);
the grey line is the median, the number the ledger reports. Status and saved core-h stay in the JSON.
Usage: uv run python scripts/figures/fig_feature_ledger.py [LEDGER.json] [OUT.png]
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

LABEL = {"dtmix_0.25": r"$t_{mix}$, $\chi=0.25$", "dtmix_0.50": r"$t_{mix}$, $\chi=0.50$",
         "edr_95": "EDR, 95th pct", "tau_95_qss": r"$\tau$, 95th pct", "kLa_1T_10": r"$k_La$"}


def level(run: str) -> int:
    m = re.search(r"_l(\d+)", run)
    return int(m.group(1)) if m else 9


def clean_rows(ledger: dict) -> dict:
    return {k: r for k, r in ledger.items() if r["n_clean"] > 0}


def main(src: Path, out: Path) -> None:
    rows = clean_rows(json.loads(src.read_text()))
    ncol = max(len(r["error"]) for r in rows.values())
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(len(rows), ncol, figsize=(2.6 * ncol, 2.6 * len(rows)), squeeze=False)
    for i, (mech, r) in enumerate(rows.items()):
        clean = [p for p in r["runs"] if p["pair"] == "same-source"]
        for j in range(ncol):
            ax = axes[i, j]
            if j >= len(r["error"]):
                ax.set_axis_off()
                continue
            q = list(r["error"])[j]
            pts = [(k, 100 * p["error"][q], level(p["with"])) for k, p in enumerate(clean)
                   if p["error"][q] is not None and math.isfinite(p["error"][q])]
            for k, y, lv in pts:
                ax.plot(k, y, ls="none", marker=fs.LEVEL_MARKER.get(lv, "h"), color=fs.level_colour(lv), ms=6)
            if r["error"][q] is not None:
                ax.axhline(100 * r["error"][q], color=fs.NEUTRAL, lw=1.2)
            ax.axhline(0, color=fs.GUIDE, lw=0.8)
            ax.set_xticks([])
            ax.set_ylabel(f"{LABEL.get(q, q)} (%)")
            fs.qoi_axes(ax, q)
            ax.grid(axis="y", **fs.GRID_KW)
        axes[i, 0].set_title(mech.replace("_", " "), loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(Path(a[0]) if a else ROOT / "experiments" / "feature_ledger.json",
         Path(a[1]) if len(a) > 1 else ROOT / "experiments" / "feature_ledger.png")
