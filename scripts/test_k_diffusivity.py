"""Test K (diffusivity scan): pre-registered rule, table, descriptive s*(L), summary JSON, figure.

Pre-registration (diary.md, "PRE-REGISTERED test K"): s = 10^k multiplies tracer and oxygen
diffusivity; L6-L8 k = 0..5, L9 k in {0, 1, 3, 5} (k = 1 at L9 added 2026-10-10 before any L8/L9 result:
the rule needs (i) at L9 and the original design had no L9 s = 10 run).
  (i)  grid-controlled at L: |dtmix(s=10)/dtmix(s=1) - 1| < 0.10 for all three thresholds (L7, L8, L9).
       L9 s = 1e3 vs s = 1 is reported only.
  (ii) physics-controlled: |dtmix_L9/dtmix_L8 - 1| <= 0.10 at s = 1e5, each threshold.
  Rejection: dtmix at L8 changes > 30% between s = 1 and 10 (any threshold).
Verdict precedence implemented: REJECTED (decidable from L8 s=1, s=10 alone) > PENDING (a needed
value missing or null) > NOT WELL-POSED AT PHYSICAL D > INCONCLUSIVE.

Usage: uv run python scripts/test_k_diffusivity.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "experiments/multifidelity"

THRESH = ("dtmix_0.50", "dtmix_0.75", "dtmix_0.95")
QOIS = THRESH + ("kLa_1T_25", "kLa_25")
GRID = {6: range(6), 7: range(6), 8: range(6), 9: (0, 1, 3, 5)}
TOL = 1e-9  # float guard so a ratio of exactly 1.10 sits on the stated boundary


def _val(data, L, k, q):
    """value, or None if run missing / null."""
    v = data.get(L, {}).get(k, {}).get(q)
    return None if v is None or (isinstance(v, float) and not math.isfinite(v)) else float(v)


def _rel(data, L, k, q, k0=0):
    a, b = _val(data, L, k, q), _val(data, L, k0, q)
    return None if a is None or b is None else a / b - 1.0


def evaluate_rule(data):
    """Pure rule. data[L][k][qoi] = float|None; absent run = absent key."""
    missing, failing = [], []

    def need(L, k, q=None):
        for t in ([q] if q else THRESH):
            if _val(data, L, k, t) is None:
                m = f"L{L} s=1e{k} {t}"
                if m not in missing:
                    missing.append(m)

    cond_i = {}
    for L in (7, 8, 9):
        need(L, 0); need(L, 1)
        ch = {t: _rel(data, L, 1, t) for t in THRESH}
        ok = None if any(v is None for v in ch.values()) else all(round(abs(v), 12) < 0.10 - TOL for v in ch.values())
        cond_i[L] = dict(changes=ch, holds=ok)
        if ok is False:
            failing.append(f"(i) L{L}")
    need(9, 3)
    cond_i[9]["reported_s1e3_vs_s1"] = {t: _rel(data, 9, 3, t) for t in THRESH}

    need(8, 5); need(9, 5)
    rat = {}
    for t in THRESH:
        a, b = _val(data, 9, 5, t), _val(data, 8, 5, t)
        rat[t] = None if a is None or b is None else a / b - 1.0
    ok2 = None if any(v is None for v in rat.values()) else all(round(abs(v), 12) <= 0.10 + TOL for v in rat.values())
    cond_ii = dict(changes=rat, holds=ok2)
    if ok2 is False:
        failing.append("(ii) L9 vs L8 at s=1e5")

    ch8 = cond_i[8]["changes"]
    known = [v for v in ch8.values() if v is not None]
    rejected = any(round(abs(v), 12) > 0.30 + TOL for v in known)

    if rejected:
        verdict = "REJECTED"
    elif missing:
        verdict = "PENDING"
    elif cond_i[7]["holds"] and cond_i[8]["holds"] and cond_i[9]["holds"] and cond_ii["holds"]:
        verdict = "NOT WELL-POSED AT PHYSICAL D"
    else:
        verdict = "INCONCLUSIVE"
    return dict(verdict=verdict, cond_i=cond_i, cond_ii=cond_ii, rejected=rejected,
                missing=missing, failing=failing)


def s_star(data, L, q, thr=0.10):
    """Smallest scanned s with |log(q(s)/q(1))| > thr; None if never (or q(1) missing); 'n/a' if no data."""
    base = _val(data, L, 0, q)
    if base is None or base <= 0:
        return "n/a"
    for k in sorted(data.get(L, {})):
        v = _val(data, L, k, q)
        if v is not None and v > 0 and abs(math.log(v / base)) > thr:
            return 10.0 ** k
    return None


def load():
    data, pending = {}, []
    for L, ks in GRID.items():
        for k in ks:
            p = ROOT / f"runs/k_l{L}_s1e{k}/results.json"
            if p.exists():
                d = json.load(open(p))
                data.setdefault(L, {})[k] = {q: d.get(q) for q in QOIS}
            else:
                pending.append(f"L{L} s=1e{k}")
    return data, pending


def fmt(v, nd=3):
    return "null" if v is None else f"{v:.{nd}f}"


def fpct(v):
    return "null" if v is None else f"{100 * v:+.1f}%"


def print_table(data, pending):
    print("TABLE: level x s x QoI (value; change from s=1 of same level)")
    print(f"{'L':>2} {'s':>6} " + " ".join(f"{q:>22}" for q in QOIS))
    for L, ks in GRID.items():
        for k in ks:
            if k not in data.get(L, {}):
                print(f"{L:>2} {'1e%d' % k:>6} pending")
                continue
            cells = [f"{fmt(_val(data, L, k, q), 2):>10} ({fpct(_rel(data, L, k, q)):>9})" for q in QOIS]
            print(f"{L:>2} {'1e%d' % k:>6} " + " ".join(f"{c:>22}" for c in cells))
    print(f"\nPending runs ({len(pending)}): " + (", ".join(pending) if pending else "none"))


def print_rule(r):
    print("\nPRE-REGISTERED RULE")
    for L in (7, 8, 9):
        c = r["cond_i"][L]
        nums = ", ".join(f"{t}: {fpct(v)}" for t, v in c["changes"].items())
        print(f"  (i) L{L} grid-controlled, |dtmix(10)/dtmix(1)-1| < 10%: {c['holds']}  [{nums}]")
    nums = ", ".join(f"{t}: {fpct(v)}" for t, v in r["cond_i"][9]["reported_s1e3_vs_s1"].items())
    print(f"  L9 s=1e3 vs s=1 (reported only): {nums}")
    nums = ", ".join(f"{t}: {fpct(v)}" for t, v in r["cond_ii"]["changes"].items())
    print(f"  (ii) |dtmix_L9/dtmix_L8 - 1| <= 10% at s=1e5: {r['cond_ii']['holds']}  [{nums}]")
    ch8 = ", ".join(f"{t}: {fpct(v)}" for t, v in r["cond_i"][8]["changes"].items())
    print(f"  Rejection (L8 change s=1 -> 10 > 30%, any threshold): {r['rejected']}  [{ch8}]")
    if r["missing"]:
        print("  Missing/null needed values: " + "; ".join(r["missing"]))
    if r["failing"]:
        print("  Failing conditions: " + "; ".join(r["failing"]))
    print(f"  VERDICT: {r['verdict']}")


def sstar_all(data):
    out = {}
    print("\ns*(L) (descriptive, not pre-registered): smallest scanned s with |log(q(s)/q(1))| > 0.10")
    for L in GRID:
        out[L] = {q: s_star(data, L, q) for q in ("dtmix_0.95", "kLa_1T_25")}
        txt = ", ".join(f"{q}: " + ("n/a (no s=1 value)" if v == "n/a" else "not reached in scan" if v is None else f"{v:g}")
                        for q, v in out[L].items())
        print(f"  L{L}: {txt}")
    return out


def figure(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT))
    from scripts import figstyle as fs
    plt.rcParams.update(fs.rcparams())
    f, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for ax, q, mk, yl in ((axes[0], "dtmix_0.95", fs.MK["chi_0.95"], r"$\Delta t_{mix,0.95}$ (s)"),
                          (axes[1], "kLa_1T_25", fs.MK["cstar_25"], r"$k_La$ (h$^{-1}$)")):
        for L in GRID:
            pts = [(10.0 ** k, _val(data, L, k, q)) for k in sorted(data.get(L, {}))]
            pts = [(x, y) for x, y in pts if y is not None and y > 0]
            if pts:
                ax.plot(*zip(*pts), label=f"L{L}",
                        **fs.series_kw(fs.level_colour(L), mk, stat="mean", ours=True))
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$D / D_{phys}$"); ax.set_ylabel(yl)
        ax.grid(**fs.GRID_KW)
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    f.tight_layout()
    path = OUT / "test_k_diffusivity.png"
    f.savefig(path, dpi=200, bbox_inches="tight")
    return path


def main():
    data, pending = load()
    print_table(data, pending)
    r = evaluate_rule(data)
    print_rule(r)
    ss = sstar_all(data)
    OUT.mkdir(parents=True, exist_ok=True)
    summary = dict(
        runs={f"L{L}": {f"1e{k}": v for k, v in ks.items()} for L, ks in data.items()},
        pending=pending, rule={**r, "cond_i": {f"L{L}": v for L, v in r["cond_i"].items()}},
        s_star_descriptive_not_preregistered={f"L{L}": v for L, v in ss.items()})
    p = OUT / "test_k_summary.json"
    json.dump(summary, open(p, "w"), indent=2)
    print(f"\nWrote {p}\nWrote {figure(data)}")


if __name__ == "__main__":
    main()
