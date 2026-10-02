"""chi = 0.75 mixing rate 1/dtmix: joint-model fit (shared p) and its h -> 0 limit.

Left: 1/dtmix vs cell size h/h10 at each rpm (L6-L9 markers, fitted
r_inf + C h^p curves). Right: extrapolated physical rate r_inf vs rpm with 95%
bounds, against Kim's 1/dtmix and our L9. Also prints, per rpm, the share of
the L9 rate that the fit attributes to the grid term, C h9^p / (r_inf + C h9^p).
That share is a MODEL INFERENCE, not a measured numerical diffusivity.
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.argv.append("--rate")
sys.path.insert(0, str(ROOT))
import scripts.joint_grid_model as jm   # noqa: E402
from scripts import figstyle as fs     # noqa: E402

thr, lv = 0.75, (6, 7, 8, 9)
d = jm.data(thr, lv)
out, p, pU, s, dof, _ = jm.fit(thr, lv)
# refit to get C (fit() returns q_inf only): least squares again with p fixed
js = sorted(set(d[:, 0].astype(int)))
r_inf, C = [], []
for j in js:
    m = d[:, 0] == j
    A = np.column_stack([np.ones(m.sum()), d[m, 1] ** p])
    c, *_ = np.linalg.lstsq(A, d[m, 2], rcond=None)
    r_inf.append(c[0]); C.append(c[1])
r_inf, C = np.array(r_inf), np.array(C)
share = C * 2.0 ** p / (r_inf + C * 2.0 ** p)
for j, sh in zip(js, share):
    print(f"rpm {jm.RPMS[j]:5g}: grid term = {sh:.0%} of the L9 rate")

plt.rcParams.update(fs.rcparams())
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
cmap = plt.get_cmap("viridis")
hg = np.linspace(0, 16, 100)
for k, j in enumerate(js):
    m = d[:, 0] == j
    col = cmap(k / (len(js) - 1))
    a1.plot(d[m, 1], d[m, 2], "o", color=col, ms=4)
    a1.plot(hg, r_inf[k] + C[k] * hg ** p, color=col, lw=0.9, label=f"{jm.RPMS[j]:g} rpm")
a1.set_xlabel(r"cell size $h/h_{10}$"); a1.set_ylabel(r"$1/\Delta t_{\mathrm{mix}}$ (1/s)")
a1.set_title(rf"$\chi=0.75$, shared $p={p:.2f}\pm{pU:.2f}$", loc="left", fontsize=10)
a1.grid(**fs.GRID_KW)
rpm = out.rpm.values
lo = r_inf - 1.0 / (out.q_inf + out.U95).values  # not used; bounds below from rate sd
sd_rate = (out.U95 / out.q_inf ** 2).values * 0 + 0
# rate-space 95% bound: delta method inverted back (U95_time = U_rate / r^2)
U_rate = out.U95.values / out.q_inf.values ** 2 * 0 + (out.U95.values * r_inf ** 2)
a2.errorbar(rpm, r_inf, yerr=U_rate, fmt="o", mfc="w", mec="#d62728", ecolor="#d62728",
            capsize=2, label=r"$r_\infty$ (h$\to$0), 95%")
a2.plot(rpm, 1.0 / out.kim.values, "^", color="#d62728", ms=7, label="Kim et al.")
l9 = d[d[:, 1] == 2.0]
a2.plot([jm.RPMS[int(j)] for j in l9[:, 0]], l9[:, 2], ":", marker="^", color="#d62728", ms=4,
        alpha=0.6, label="L9")
a2.axhline(0, color="0.5", lw=0.8)
a2.set_xlabel("rocking frequency [rpm]"); a2.set_ylabel(r"$1/\Delta t_{\mathrm{mix}}$ (1/s)")
a2.grid(**fs.GRID_KW)
a1.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=7)
a2.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=8)
fig.tight_layout()
out_p = ROOT / "experiments/multifidelity/joint_rate_chi075.png"
fig.savefig(out_p, dpi=160, bbox_inches="tight")
print("saved", out_p)
