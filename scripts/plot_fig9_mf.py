"""Fig 9 (mixing time vs rpm, 7 deg) with a threshold-fidelity MF prediction.

Kim's encoding as in plot_fig9.py: chi = 0.95 blue circles, 0.75 red
up-triangles, 0.50 green down-triangles, log-log axes over his ranges.
  Kim et al.: markers only.
  L9: markers joined by a dotted line.  L8: smaller markers, no line.
  MF: the "fidelity" is the threshold. LF = L9 dtmix at chi=0.50 (cheap, every
      rpm). HF = L9 dtmix at chi=0.75 and 0.95 at 3 rpm only (lowest, 25, highest
      available). Dashed curve + 95% band per HF threshold, same model as the
      Fig 13 MF panels (LOO-KRR + universal kriging, HF sigma 5%).
"""
import json, math, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs                    # noqa: E402
from scripts.plot_mf_h1_h2 import fit, REL_NOISE      # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
ENC = {0.95: ("#1f5fa8", "o"), 0.75: ("#d62728", "^"), 0.50: ("#2ca02c", "v")}
KIM = pd.read_csv(ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_frequency.csv").sort_values("RPM")


def val(run, key):
    f = ROOT / "runs" / run / "results.json"
    return json.load(open(f)).get(key, math.nan) if f.exists() else math.nan


plt.rcParams.update(fs.rcparams())
fig, ax = plt.subplots(figsize=(5.8, 4.4))
for thr, (c, m) in ENC.items():
    ax.plot(KIM["RPM"], KIM[f"dtmix_strict_{thr:g}"], ls="none", marker=m, color=c, ms=7)
    y9 = np.array([val(f"fig9_l9_rpm{r:g}", f"dtmix_{thr:.2f}") for r in RPMS])
    y8 = np.array([val(f"kmix_l8_rpm{r:g}", f"dtmix_{thr:.2f}") for r in RPMS])
    ax.plot(RPMS, y9, ls=":", lw=1.1, marker=m, color=c, ms=5)
    ax.plot(RPMS, y8, ls="none", marker=m, color=c, ms=3.5, alpha=0.55)
    if thr == 0.50:
        continue
    lf_all = np.array([val(f"fig9_l9_rpm{r:g}", "dtmix_0.50") for r in RPMS])
    ok = np.isfinite(lf_all) & np.isfinite(y9)
    x, lf, hf = np.array(RPMS, float)[ok], lf_all[ok], y9[ok]
    train = np.isin(x, [x.min(), 25.0, x.max()])
    _, mf, _, _ = fit(x, lf, hf, train)
    xg = np.linspace(x.min(), x.max(), 200)
    mu, sd = mf.predict(xg)
    ax.fill_between(xg, np.maximum(mu - 1.96 * sd, 1.0), mu + 1.96 * sd, color=c, alpha=0.13, lw=0)
    ax.plot(xg, mu, color=c, lw=1.3, ls="--")
    ax.plot(x[train], hf[train], ls="none", marker=m, mfc="none", mec="k", mew=1.0, ms=9)

for a_ in (ax,):
    a_.set_xscale("log"); a_.set_yscale("log"); a_.tick_params(which="both", direction="in")
from matplotlib.ticker import NullFormatter
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlim(10, 100); ax.set_ylim(1, 1000)
ax.set_xlabel(r"$f_b$ (rpm)"); ax.set_ylabel(r"$\Delta t_{\mathrm{mix}}$ (s)")

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
q = [Line2D([], [], ls="none", marker=ENC[t][1], color=ENC[t][0], ms=6.5, label=rf"$\chi={t:.2f}$") for t in ENC]
d = [Line2D([], [], ls="none", marker="o", color="0.3", ms=6.5, label="Kim et al."),
     Line2D([], [], ls=":", marker="o", color="0.3", ms=5, label="L9"),
     Line2D([], [], ls="none", marker="o", color="0.3", ms=3.5, alpha=0.55, label="L8"),
     Line2D([], [], ls="--", color="0.3", label=r"MF from $\chi=0.50$"),
     Line2D([], [], ls="none", marker="o", mfc="none", mec="k", ms=9, label="MF training")]
l1 = fig.legend(handles=q, loc="upper left", bbox_to_anchor=(1.0, 0.95), frameon=False, fontsize=9)
fig.add_artist(l1)
fig.legend(handles=d, loc="upper left", bbox_to_anchor=(1.0, 0.62), frameon=False, fontsize=9)
fig.tight_layout()
out = ROOT / "experiments/kimetal2024/figure_replicas/replicated_Fig09_mf.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved", out)
