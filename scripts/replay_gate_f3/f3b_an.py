"""F3b analysis: Ls=8 vs coupled (1e-6), QoIs per Ls and the diary rule. Usage: f3b_an.py [run ids...]"""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import replay_gate as rg

if len(sys.argv) > 1 and sys.argv[1] == "check8":
    a = rg.load_tr_oxy(rg.run_dir("f3b_coupled") / "tr_oxy.dat")
    b = rg.load_tr_oxy(rg.run_dir("f3b_ls8") / "tr_oxy.dat")
    d = rg.rel_diff_by_time(a, b, [2, 8, 9])
    print("Ls=8 vs coupled: common", d["n_common"], "of", d["n_b"], "max rel", d["max"], "t_last", b[-1, 1])
    sys.exit()

q = ("kLa_1T_10", "kLa_1T_25", "dtmix_0.50")
res = {}
for ls in (8, 9, 10):
    res[ls] = rg.load_results(rg.run_dir(f"f3b_ls{ls}"))
    print(ls, {k: res[ls][k] for k in q})
r98 = {k: res[9][k] / res[8][k] for k in q}
r109 = {k: res[10][k] / res[9][k] for k in q}
print("q(9)/q(8):", r98)
print("q(10)/q(9):", r109)
dev = {k: abs(v - 1) for k, v in r109.items()}
n_conv = sum(1 for v in dev.values() if v <= 0.10)
n_not = sum(1 for v in dev.values() if v >= 0.25)
print("|q(10)/q(9)-1|:", dev)
print("VERDICT:", "CONVERGING" if n_conv >= 2 else "NOT CONVERGING" if n_not >= 2 else "MIXED")
