"""Fig 9 with the grid-extrapolated mixing time (E&H, L6-L9 per rpm) vs Kim.

Kim's encoding (colour/marker = chi). Kim: markers only. L9: small markers,
dotted. Extrapolated h -> 0: hollow-edged markers with the 95% E&H bound.
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs                       # noqa: E402
from scripts.diag_dtmix_extrapolation import table, KIM   # noqa: E402

ENC = {0.95: ("#1f5fa8", "o"), 0.75: ("#d62728", "^"), 0.50: ("#2ca02c", "v")}
t = table()
plt.rcParams.update(fs.rcparams())
fig, ax = plt.subplots(figsize=(5.8, 4.4))
for thr, (c, m) in ENC.items():
    s = t[t.thr == thr]
    ax.plot(KIM.index, KIM[f"dtmix_strict_{thr:g}"], ls="none", marker=m, color=c, ms=7)
    ax.plot(s.rpm, s.L9, ls=":", lw=1.0, marker=m, color=c, ms=3.5, alpha=0.6)
    ok = np.isfinite(s.phi0)
    lo = np.minimum(s.U[ok], 0.95 * s.phi0[ok])
    ax.errorbar(s.rpm[ok] * 1.01, s.phi0[ok], yerr=[lo, s.U[ok]], fmt=m, mfc="w", mec=c,
                ecolor=c, elinewidth=0.9, capsize=2, ms=6, mew=1.3)
ax.set_xscale("log"); ax.set_yscale("log"); ax.tick_params(which="both", direction="in")
from matplotlib.ticker import NullFormatter
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlim(10, 100); ax.set_ylim(1, 3000)
ax.set_xlabel(r"$f_b$ (rpm)"); ax.set_ylabel(r"$\Delta t_{\mathrm{mix}}$ (s)")
from matplotlib.lines import Line2D
q = [Line2D([], [], ls="none", marker=ENC[k][1], color=ENC[k][0], ms=6.5, label=rf"$\chi={k:.2f}$") for k in ENC]
d = [Line2D([], [], ls="none", marker="o", color="0.3", ms=6.5, label="Kim et al."),
     Line2D([], [], ls=":", marker="o", color="0.3", ms=3.5, alpha=0.6, label="L9"),
     Line2D([], [], ls="none", marker="o", mfc="w", mec="0.3", mew=1.3, ms=6, label=r"$h\to0$ (L6-L9), 95%")]
l1 = fig.legend(handles=q, loc="upper left", bbox_to_anchor=(1.0, 0.95), frameon=False, fontsize=9)
fig.add_artist(l1)
fig.legend(handles=d, loc="upper left", bbox_to_anchor=(1.0, 0.62), frameon=False, fontsize=9)
fig.tight_layout()
out = ROOT / "experiments/kimetal2024/figure_replicas/replicated_Fig09_extrapolated.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved", out)
