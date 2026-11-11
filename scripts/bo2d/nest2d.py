"""Nested space-filling 2D design (user go 2026-11-11 ~10:30: "20 L8 + 8 L9; reuse simulations we have if
pseudo-random; nobody should raise an eyebrow at the scatter; L9-only may need more L9 than the chain").

Design (fixed before any result):
- L7: 40 of the 200 scr_l7 screen runs (maximin LHS, pseudo-random; no new compute). L8: 20 WITHIN the 40. L9: 8
  of the 20. Each subset: min coverage radius over 4000 random subsets of the parent (v1 greedy maximin was
  edge-biased: no L8 below 20 rpm, half the L9 in corners).
- New runs: L8 warm from the scr_l7 run at the same site (30 settle cycles, release, 40), tag ns; L9 warm from that
  ns L8. Training sets: chain = 40 L7 / 20 L8 / 8 L9; L9-only GP = the 8 L9 (more space-filling L9 later if needed).
- Test set: the 33 old L9 runs (bo1 / bo4 / fig9), never trained on (they sit near the tau cap and theta = 7).
Accounts: default (priority QoS, 312 CPU per user) first; overflow to mbessa-condo up to 160 CPU (user, 2026-11-11).
  uv run python scripts/bo2d/nest2d.py design        -> state/nest2d_design.json, figures/bo2d/fig_nest2d_design.png
  uv run python scripts/bo2d/nest2d.py submit l8|l9 [--go]
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.bo2d import submit as S  # noqa: E402

C = S.C
DESIGN = S.STATE / "nest2d_design.json"
TAG = "ns"
CONDO_CAP = 160


def unit(P):
    return (np.asarray(P) - [C.RPM[0], C.THETA[0]]) / [C.RPM[1] - C.RPM[0], C.THETA[1] - C.THETA[0]]


GRID = np.stack(np.meshgrid(np.linspace(0, 1, 41), np.linspace(0, 1, 41)), -1).reshape(-1, 2)


def cover(P, n, tries=4000, seed=0):
    """Subset of n rows of P minimising the coverage radius (max over a 41x41 grid of the distance to the nearest
    chosen point): spread out WITHOUT the edge bias of greedy maximin (v1 of this design had no L8 below 20 rpm)."""
    U = unit(P); rng = np.random.default_rng(seed); best, bi = np.inf, None
    for _ in range(tries):
        i = rng.choice(len(U), n, replace=False)
        r = np.sqrt(np.min(np.sum((GRID[:, None] - U[i][None]) ** 2, axis=2), axis=1)).max()
        if r < best:
            best, bi = r, i
    return [int(k) for k in bi], float(best)


def screen_sites():
    out = []
    for d in sorted((ROOT / "runs").glob("scr_l7_*")):
        if (d / "checkpoint.dump").exists() and (d / "results.json").exists():
            p = json.loads((d / "params.json").read_text())
            out.append((d.name, round(float(p["omega_b"]) * 60 / (2 * np.pi), 4), float(p["theta_max"][0])))
    return out


def design():
    s = screen_sites(); P = np.array([[r, t] for _, r, t in s])
    i7, r7 = cover(P, 40); P7 = P[i7]
    i8, r8 = cover(P7, 20); i9, r9 = cover(P7[i8], 8)
    print(f"coverage radius (unit box): L7 {r7:.3f}  L8 {r8:.3f}  L9 {r9:.3f}")
    l7 = [s[i] for i in i7]; l8 = [l7[i] for i in i8]; l9 = [l8[i] for i in i9]
    DESIGN.write_text(json.dumps(dict(n_screen=len(s), l7=l7, l8=l8, l9=l9), indent=1))
    print(f"screen {len(s)} -> L7 {len(l7)} L8 {len(l8)} L9 {len(l9)}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scripts import figstyle as fs
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 4.2), sharey=True)
    import pandas as pd
    t = pd.read_csv(S.STATE / "truth_snapshot_v2_tau_95.csv"); t = t[t.level == 9]   # the 33 test L9 (truth v2)
    for L, pts, ms in ((7, l7, 4), (8, l8, 7), (9, l9, 11)):
        c = fs.level_colour(L)
        ax[0].plot([p[1] for p in pts], [p[2] for p in pts], "o", ms=ms, color=c, mfc="none", mew=1.4,
                   ls="none", label=f"L{L} ({len(pts)})")
    ax[0].set_title("training (nested)")
    c9 = fs.level_colour(9)
    ax[1].plot(t.rpm, t.theta, "o", ms=7, color=c9, mfc="w", mew=1.4, ls="none", label=f"L9 test ({len(t)})")
    ax[1].set_title("test (never trained on)")
    for a in ax:
        a.set_xlim(14.5, 38); a.set_ylim(1.8, 7.2); a.set_xticks([15, 20, 25, 30, 35]); a.grid(**fs.GRID_KW); a.set_xlabel("rocking speed (rpm)")
        a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False)
    ax[0].set_ylabel(r"rocking angle $\theta$ (deg)")
    fig.tight_layout(); out = ROOT / "figures/bo2d/fig_nest2d_design.png"; fig.savefig(out, dpi=170); print("wrote", out)


def my_cpus(account):
    o = subprocess.run(["squeue", "-u", "eaguerov", "-h", "-A", account, "-o", "%C"], capture_output=True, text=True).stdout
    return sum(int(x) for x in o.split())


def submit(level, go):
    d = json.loads(DESIGN.read_text()); ranks = S.RANKS[level]
    used = {"default": my_cpus("default"), "mbessa-condo": my_cpus("mbessa-condo")}
    for name, rpm, th in d[f"l{level}"][:int(os.environ.get("NEST_N", "999"))]:   # NEST_N: submit only the first n (tests)
        src = name if level == 8 else C.run_id(8, rpm, th, tag=TAG)
        if level == 9 and not (ROOT / "runs" / src / "shear_stress.dat").exists():
            print("L8 not finished:", src); continue
        acct = "default" if used["default"] + ranks <= 312 else ("mbessa-condo" if used["mbessa-condo"] + ranks <= CONDO_CAP else None)
        if acct is None:
            print("no room now for", rpm, th); continue
        jid = S.warm(level, src, rpm, th, arm="nest", go=go, tag=TAG, account=acct)
        if jid:
            used[acct] += ranks; print(f"submitted L{level} {rpm} {th} on {acct} job {jid}")


if __name__ == "__main__":
    {"design": lambda: design(), "submit": lambda: submit(int(sys.argv[2][1]), "--go" in sys.argv)}[sys.argv[1]]()
