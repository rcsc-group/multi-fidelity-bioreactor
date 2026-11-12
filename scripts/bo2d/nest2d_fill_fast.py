"""Genoa duplicates of fill L9 runs that landed on slow Intel Cascade Lake nodes (2026-11-11 ~21:10).
Evidence: ns_l9 15.42/4.66 ran 939 s per time unit on node1818 (cascade), its Genoa copy 470 incl. start-up (test
passed, slow run cancelled). Fill runs at 21:05: cascade 794-1037 s/unit, Genoa 202-315 s/unit. Same rule as
nest2d_dup_fast.py (user: "send in parallel on a fast node, if it is truly faster, then cancel the slow one"):
identical warm start from the same ns_l8 parent, tag nfg, pinned to Features=genoa; keep the faster, cancel the other.
  uv run python scripts/bo2d/nest2d_fill_fast.py <rpm,theta> ... [--go]
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.bo2d import nest2d as N  # noqa: E402

S, C = N.S, N.C
go = "--go" in sys.argv
for arg in [a for a in sys.argv[1:] if a != "--go"]:
    rpm, th = map(float, arg.split(","))
    src = C.run_id(8, rpm, th, tag=N.TAG)
    acct = "default" if N.my_cpus("default") + 32 <= 312 else (
        "mbessa-condo" if N.my_cpus("mbessa-condo") + 32 <= N.CONDO_CAP else None)
    if acct is None:
        print("no room for", rpm, th); continue
    jid = S.warm(9, src, rpm, th, arm="fill-fast", go=go, tag="nfg", account=acct)
    if jid:
        subprocess.run(["scontrol", "update", f"JobId={jid}", "Features=genoa"], check=True)
        print(f"submitted nfg {rpm} {th} on {acct} job {jid} (genoa)")
