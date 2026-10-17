import json
import sys
import numpy as np
sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/.claude/worktrees/agent-a0a362873e179658b")
from scripts import replay_gate as rg

def mass_drift(path):
    d = np.loadtxt(path)
    # columns: i t cf cu cc2 co cog cogn nz ntiny nneg cf2 mc2 mo mf
    m = d[:, 12]
    return float(np.max(np.abs(m - m[0])) / abs(m[0])), float(m[0]), float(m[-1]), float(np.max(np.abs(d[:, 14] - d[0, 14])) / d[0, 14])

rows = {}
res = rg.load_results(rg.run_dir("f3p_coupled"))
md = mass_drift(rg.scratch_dir("f3p_coupled") / "rec_checksum.txt")
rows["coupled"] = {"dtmix_0.50": res["dtmix_0.50"], "kLa_1T_25": res["kLa_1T_25"], "c2_total_mass_drift": md[0],
                   "c2_mass_first_last": md[1:3], "f_total_drift": md[3]}
for ls, pr in ((6, 0), (7, 0), (7, 1), (8, 0), (8, 1)):
    rid = f"f3p_ls{ls}_p{pr}"
    d = np.loadtxt(rg.scratch_dir(rid) / "replay_diag.txt")
    far, allc, um = d[:, 2], d[:, 3], d[:, 4]
    vf = rg.load_tr_oxy(rg.scratch_dir(rid) / "vol_frac_interf.dat")
    r = rg.load_results(rg.run_dir(rid))
    m = mass_drift(rg.scratch_dir(rid) / "rec_checksum.txt")
    rows[rid] = {"div_far": float(far.max() / um.max()), "div_all": float(allc.max() / um.max()),
                 "liq_vol_drift_f_liq_sum": float(np.max(np.abs(vf[:, 2] - vf[0, 2])) / vf[0, 2]),
                 "f_total_drift": m[3], "c2_total_mass_drift": m[0], "dtmix_0.50": r["dtmix_0.50"], "kLa_1T_25": r["kLa_1T_25"]}
    try:
        pj = np.loadtxt(rg.scratch_dir(rid) / "replay_project.txt")
        rows[rid]["project_iters_mean"], rows[rid]["project_resa_max"] = float(pj[:, 1].mean()), float(pj[:, 3].max())
    except OSError:
        pass
print(json.dumps(rows, indent=1))
for ls in (7, 8):
    a, b = rows[f"f3p_ls{ls}_p0"], rows[f"f3p_ls{ls}_p1"]
    print(ls, "proj/noproj - 1:", {k: b[k] / a[k] - 1 for k in ("dtmix_0.50", "kLa_1T_25")})
print("Ls=6 vs coupled - 1:", {k: rows["f3p_ls6_p0"][k] / rows["coupled"][k] - 1 for k in ("dtmix_0.50", "kLa_1T_25")})
