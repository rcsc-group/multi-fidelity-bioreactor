"""Advance prediction of bio_mix_elvis.m dtmix_0.95 from snapshots (no MATLAB).
Faithful to the script: cs==1 cells, liquid = f>0 (plain) / f>1-1e-10 (strict),
equal weights, sigma0 = first sample, dtmix = sample NEAREST chi=0.95, seconds
via Kim's T_char (H = L/3.5).  Also the chi_scale.py 'kim' block=1 estimator."""
import json, math, sys, glob
from pathlib import Path
import numpy as np
ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts import kim_convert as kc, fields, chi_scale

TH = 7 * math.pi / 180; L = 0.25; RPM = 32.5
Tper = 60 / RPM; H = L / 3.5
V = L / 4 * (H + 0.5 * L * math.tan(TH)); U = V / (H * 0.5) / Tper; TCH = L / U

def chi_series(files, times, grid, strict):
    sig = []; sigck = []
    for tk in grid:
        k = int(np.argmin(np.abs(times - tk)))
        s = chi_scale.read_snapshot(files[k], ["f", "cs", "c2"])
        m = (s["cs"] == 1) & ((s["f"] > 1 - 1e-10) if strict else (s["f"] > 0))
        x = s["c2"][m]; sig.append((x * x).mean() - x.mean() ** 2)
        mk = s["f"] > 0 if not strict else s["f"] > 1 - 1e-10
        sigck.append(chi_scale.sigma2_coarse(s["f"], s["cs"], s["c2"], 1, "kim") if not strict else np.nan)
    sig = np.array(sig); sigck = np.array(sigck)
    return 1 - sig / sig[0], 1 - sigck / sigck[0]

def dtmix(chi, grid):
    if chi.max() < 0.95: return float("nan")
    i = int(np.argmin(np.abs(chi - 0.95)))
    return (grid[i] - grid[0]) * TCH

def run(rid, ini, dt, fin, nmax=None):
    fd = ROOT / "runs" / rid / "fields"
    files = sorted(fd.glob("snap_*.bin"))
    times = np.array([kc._header_t(p) for p in files])
    grid = kc.kim_grid(ini, dt, fin)
    if nmax: grid = grid[:nmax]
    out = {}
    for strict in (False, True):
        chi, chick = chi_series(files, times, grid, strict)
        out["strict" if strict else "plain"] = dtmix(chi, grid)
        if not strict:
            out["chi_scale_kim_block1"] = dtmix(chick, grid)
            first = int(np.argmax(chi >= 0.95)) if chi.max() >= 0.95 else None
            out["first_idx_plain"] = first
    out["nfiles"] = len(grid)
    out["ours_results_json"] = json.load(open(ROOT / "runs" / rid / "results.json"))["dtmix_0.95"]
    return out, times

if __name__ == "__main__":
    mode = sys.argv[1]
    ini = 48.5863329
    if mode == "kp1":
        for rid in ("scal_l8_rpm32.5_th7", "scal_l9_rpm32.5_th7"):
            fd = ROOT / "runs" / rid / "fields"
            tl = kc._header_t(sorted(fd.glob("snap_*.bin"))[-1])
            dt = 0.1519 * 7
            n = int(math.floor((tl - ini) / dt))
            fin = ini + n * dt
            print(rid, "fin", fin, run(rid, ini, dt, fin)[0])
    else:
        rid = "scal_l8_rpm32.5_th7"
        fd = ROOT / "runs" / rid / "fields"
        times = np.array([kc._header_t(p) for p in sorted(fd.glob("snap_*.bin"))])
        dt = round((times[-1] - times[0]) / (len(times) - 1), 9)
        print("mean snapshot spacing dt_File =", dt)
        full, _ = run(rid, ini, dt, times[-1])
        print("full-series", full)
        n = full["first_idx_plain"] + 31
        fin = ini + (n - 1) * dt
        print("truncate at nfiles", n, "fin", repr(fin))
        print("TRUNCATED prediction", run(rid, ini, dt, fin)[0])
