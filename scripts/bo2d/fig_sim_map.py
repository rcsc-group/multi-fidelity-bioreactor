"""Where are the simulations? Every real run in the frozen snapshots (union over kLa_1T_10, dtmix_0.25, tau_95), one
panel per level, in (rpm, theta). Colour = level (figstyle). Marker = how the point was chosen:
  o  space-filling design (campaign2d L8 maximin LHS; L7 seeding at the L8 sites: bo_l7)
  ^  max-uncertainty acquisition (campaign2d `next`: max 2 sigma/mean of the MF model, L8 sites only)
  s  1D tests (mostly theta = 7) (fig9, T5, T6, kmix tests)
  D  BO, L9-only arm (constrained EI; its L8 spin-up runs too)
  v  BO, chain arm (cascade constrained EI, L7 -> L8 -> L9)
No annotations.  uv run python scripts/bo2d/fig_sim_map.py OUT.png"""
import re
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

S = ROOT / "scripts/bo2d/state"
d = pd.concat([pd.read_csv(S / f"truth_snapshot_{q}.csv") for q in ("kLa_1T_10", "dtmix_0.25", "tau_95")])
d = d[d.level.isin([7, 8, 9])].drop_duplicates("run")
ORIGIN = [(r"^c2d_l8|^bo_l7", "o", "space-filling design"), (r"^c2d_l9", "^", "max-uncertainty acquisition"),
          (r"^(fig9|t5|t6|kmix)_", "s", "1D tests (mostly theta = 7)"), (r"^bo1_", "D", "BO, L9-only arm"),
          (r"^bo4_", "v", "BO, chain arm")]


def origin(run):
    for pat, mk, lab in ORIGIN:
        if re.search(pat, run):
            return mk, lab
    return "x", "other"


d[["mk", "lab"]] = [origin(r) for r in d.run]
plt.rcParams.update(fs.rcparams())
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.3), sharey=True)
for ax, lev in zip(axs, (7, 8, 9)):
    s = d[d.level == lev]
    for _, mk, lab in ORIGIN + [(None, "x", "other")]:
        t = s[s.lab == lab]
        if len(t):
            ax.scatter(t.rpm, t.theta, marker=mk, s=26, facecolors="none" if mk in "o^" else fs.LEVEL_COLOUR[lev],
                       edgecolors=fs.LEVEL_COLOUR[lev], lw=1.0, label=lab)
    ax.set_title(f"L{lev}  (n = {len(s)})"); ax.set_xlabel("rpm")
axs[0].set_ylabel("theta (deg)")
from matplotlib.lines import Line2D
h = [Line2D([], [], ls="", marker=mk, mec="0.3", mfc="none" if mk in "o^" else "0.3", label=lab) for _, mk, lab in ORIGIN]
axs[-1].legend(handles=h, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=200, bbox_inches="tight")
print(d.groupby(["level", "lab"]).size().to_string())
print("theta of 1D-line runs:", sorted(d[d.lab == "1D tests (mostly theta = 7)"].theta.round(2).unique()))
