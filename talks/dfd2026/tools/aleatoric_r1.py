"""Aleatoric spread of the slide-7 top statistics (tau_amp_t, edr_max_t) at L10, theta 7.

Definition used: spread of a QoI between estimates at the same (rpm, theta, level) that differ only
in what we do not control (which cycles are in the window, start phase, start-up history).

Two estimates:
 (a) window-to-window: slide the N_CYC=3 cycle window (qois() protocol) over every complete
     post-transient cycle of each run; sd of ln(statistic) across non-overlapping windows.
 (b) run-to-run: every pair of runs at the same rpm, theta 7, L10 (different start histories);
     |ln ratio| of the qois() values.
Output: data/aleatoric_r1.json and a printed table.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import scripts.hydro_dataset_v2 as H  # noqa: E402

RUNS = REPO / "runs"
SKIP_CYC = 2  # cycles dropped as start-up transient after the restart


def per_window(run):
    d = RUNS / run
    p = json.load(open(d / "params.json"))
    arr = H._load_dat(d / "shear_stress.dat")
    T_bio, T = H._t_scales(p)
    a = p.get("geometry", {}).get("a", 0.25)
    U = a / T_bio
    scale, escale = H.RHO * U**2, H.RHO * U**3 / a
    cyc = arr[:, 1] / T
    c0 = math.ceil(p.get("t_checkpoint", 0.0) / T - 1e-6) + SKIP_CYC
    c1 = math.floor(cyc.max() + 1e-6)
    out = []
    for s in range(c0, c1 - H.N_CYC + 1, H.N_CYC):  # non-overlapping windows
        m = (cyc >= s) & (cyc < s + H.N_CYC)
        if m.sum() < 10:
            continue
        out.append((0.5 * (arr[m, 10].max() - arr[m, 10].min()) * scale, arr[m, H.COL_EDISS].max() * escale))
    return np.array(out)


def main():
    cands = sorted(r.name for r in RUNS.iterdir()
                   if (r / "shear_stress.dat").exists() and (r / "params.json").exists()
                   and (r.name.startswith("l10_fig13a_xlevel_rpm") and not r.name.endswith("_e95")
                        or (r.name.startswith(("t10_l10_rpm", "x8_l10_rpm")) and "_th7" in r.name and not r.name.endswith("_ext"))))
    res = {}
    for run in cands:
        p = json.load(open(RUNS / run / "params.json"))
        rpm = round(p["omega_b"] * 60 / (2 * math.pi), 3)
        if abs(p["theta_max"][0] - 7.0) > 1e-6:
            continue
        try:
            q = H.qois(run)
            w = per_window(run)
        except Exception as e:  # noqa: BLE001
            print(f"skip {run}: {e}")
            continue
        lw = np.log(w) if len(w) else np.empty((0, 2))
        res[run] = dict(rpm=rpm, tau_amp_t=q["tau_amp_t"], edr_max_t=q["edr_max_t"], window=q["window"],
                        n_win=int(len(w)), sd_ln_tau=float(lw[:, 0].std(ddof=1)) if len(w) > 1 else None,
                        sd_ln_edr=float(lw[:, 1].std(ddof=1)) if len(w) > 1 else None)
        r = res[run]
        print(f"{run:34s} rpm {rpm:5.2f} win {r['n_win']:2d}  sd ln tau {r['sd_ln_tau']}  sd ln edr {r['sd_ln_edr']}")
    pairs = []
    names = sorted(res)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if abs(res[a]["rpm"] - res[b]["rpm"]) < 1e-6:
                pairs.append(dict(a=a, b=b, rpm=res[a]["rpm"],
                                  dln_tau=float(math.log(res[a]["tau_amp_t"] / res[b]["tau_amp_t"])),
                                  dln_edr=float(math.log(res[a]["edr_max_t"] / res[b]["edr_max_t"]))))
                print(f"pair {a} / {b}: dln tau {pairs[-1]['dln_tau']:+.3f}  dln edr {pairs[-1]['dln_edr']:+.3f}")
    json.dump(dict(runs=res, pairs=pairs, skip_cyc=SKIP_CYC, n_cyc=H.N_CYC), open(TALK / "data/aleatoric_r1.json", "w"), indent=1)


if __name__ == "__main__":
    main()
