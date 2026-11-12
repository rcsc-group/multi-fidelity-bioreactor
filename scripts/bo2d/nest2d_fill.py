"""L9 fill for the L9-only model (2026-11-11 ~19:30, user: "we will need more L9 cases for the L9-only GPR, because it
is looking quite sparse ... Can we do that overnight? ... send a more filling-space L9 set (try to reuse whatever we
already had just like the L7 ones you just salvaged)").
Reuse: the 20 nested L8 runs (ns_l8_) are already a space-filling set (min coverage radius inside the 40 L7, nest2d.py),
and 8 of them already carry an L9. The other 12 get an L9 warm start from their finished L8 (same recipe as the nested
L9: 30 settle cycles, release, 40), tag nf. Then L9 only = 20 L9 on the L8 sites; the chain keeps its 4 L9.
Test set unchanged (the 33 old cold L9). Accounts as nest2d.submit: default up to 312 CPU, mbessa-condo up to 160.
  uv run python scripts/bo2d/nest2d_fill.py [--go]
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.bo2d import nest2d as N  # noqa: E402

S, C = N.S, N.C
TAG = "nf"


def main(go):
    d = json.loads(N.DESIGN.read_text()); ranks = S.RANKS[9]
    have = {(round(r, 2), round(t, 2)) for _, r, t in d["l9"]}
    used = {"default": N.my_cpus("default"), "mbessa-condo": N.my_cpus("mbessa-condo")}
    print("CPU in use:", used)
    todo = [(r, t) for _, r, t in d["l8"] if (round(r, 2), round(t, 2)) not in have]
    print(len(todo), "L8 sites without an L9")
    for rpm, th in todo:
        src = C.run_id(8, rpm, th, tag=N.TAG)
        if not (ROOT / "runs" / src / "shear_stress.dat").exists():
            print("L8 not finished:", src); continue
        acct = ("default" if used["default"] + ranks <= 312
                else ("mbessa-condo" if used["mbessa-condo"] + ranks <= N.CONDO_CAP else None))
        if acct is None:
            print("no room now for", rpm, th); continue
        jid = S.warm(9, src, rpm, th, arm="fill", go=go, tag=TAG, account=acct)
        if jid:
            used[acct] += ranks; print(f"submitted L9 {rpm} {th} on {acct} job {jid}")


if __name__ == "__main__":
    main("--go" in sys.argv)
