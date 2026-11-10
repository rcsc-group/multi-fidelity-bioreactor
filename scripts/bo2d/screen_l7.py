"""Dense L7 screen (user-approved 2026-11-10: "screen L7 very well, then some L8, then L9").
200 cold L7 runs on a maximin LHS of the campaign domain (seed 20261110), skipping points within 0.02 (scaled) of an
existing L7 run. Same cold protocol and binary as the bo4_l7 / bo_l7 runs (submit.cold: 80-cycle spin-up, release,
150 cycles; 8 ranks; measured 5.1 core-h, 0.62 h wall, n = 38). Tag "scr" -> run ids scr_l7_rpmR_thT.
Throttle: at most MAX_INFLIGHT of these jobs pending or running at once (shared cluster: 30 x 8 ranks = 240 cores).
  uv run python scripts/bo2d/screen_l7.py plan        # print the points, submit nothing
  uv run python scripts/bo2d/screen_l7.py run --go    # throttled submit loop (run inside a SLURM job)"""
import subprocess
import sys
import time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import submit as S  # noqa: E402
C = S.C
N, SEED, MAX_INFLIGHT, MIN_D = 200, 20261110, 30, 0.02


def points():
    P = C.lhs_maximin(N, seed=SEED, tries=60)
    old = [(float(r), float(t)) for r, t in C.load_db("tau_95").query("level == 7")[["rpm", "theta"]].to_numpy()]
    sc = lambda x: (np.asarray(x) - C.DOMAIN[:, 0]) / (C.DOMAIN[:, 1] - C.DOMAIN[:, 0])
    if old:
        dmin = np.min(np.linalg.norm(sc(P)[:, None] - sc(old)[None], axis=2), axis=1)
        P = P[dmin > MIN_D]
    return [(round(float(r), 3), round(float(t), 3)) for r, t in P]


def inflight(jobs):
    if not jobs:
        return 0
    out = subprocess.run(["squeue", "-h", "-j", ",".join(jobs), "-o", "%i"], capture_output=True, text=True).stdout
    return len(out.split())


if __name__ == "__main__":
    pts = points()
    if sys.argv[1] == "plan":
        print(len(pts), "points"); [print(p) for p in pts[:10]]; sys.exit()
    go = "--go" in sys.argv; jobs = []
    for r, t in pts:
        while go and inflight(jobs) >= MAX_INFLIGHT:
            time.sleep(60)
        jid = S.cold(7, r, t, arm="screen", go=go, tag="scr")
        if jid:
            jobs.append(jid); print(time.strftime("%H:%M"), "submitted", jid, r, t, flush=True)
    while go and inflight(jobs):
        time.sleep(120)
    print("SCREEN DONE", len(jobs), "jobs", flush=True)
