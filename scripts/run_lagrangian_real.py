"""Run the Lagrangian estimator on a recorded run (fields/snap_*.bin).

Usage: uv run python scripts/run_lagrangian_real.py RUN_DIR [--npc 16 64] [--periods 10]
Release = first snapshot with t >= n_mix_cycles * T_nd.  The recorded flow from
release is used directly (no reuse) for `periods` periods.  Writes
OUT/summary.json, OUT/chi_<npc>.npz and OUT/chi_l32.png.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import jax

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts import lagrangian_mix as lm  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--npc", type=int, nargs="+", default=[16, 64])
    ap.add_argument("--periods", type=int, default=10)
    ap.add_argument("--substeps", type=int, default=2)
    ap.add_argument("--nb", type=int, default=64, help="fine boxes per side (l = L0/nb)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    run = Path(args.run_dir)
    params = json.load(open(run / "params.json"))
    cu = lm.code_units(params)
    T = cu["T_nd"]
    out = Path(args.out or f"/oscar/scratch/eaguerov/lagrangian/{run.name}")
    out.mkdir(parents=True, exist_ok=True)
    t_mix = params["n_mix_cycles"] * T
    print("code units", cu, flush=True)

    w = lm.load_window(run / "fields", planes=("ux", "uy", "f", "cs", "c2"),
                       t_min=t_mix - 1e-9, t_max=t_mix + args.periods * T + 2 * T / 100)
    times, n = w["times"], w["n"]
    S = len(times)
    print(f"{S} snapshots, n = {n}, t in [{times[0]:.4f}, {times[-1]:.4f}], "
          f"release t_mix = {t_mix:.4f} (= cycle {times[0] / T:.4f})", flush=True)
    if times[-1] < t_mix + args.periods * T - 1e-3 * T:
        raise SystemExit(f"snapshots cover only {(times[-1] - times[0]) / T:.2f} periods")
    dt_snap = float(np.median(np.diff(times)))
    domain = (-0.5, -0.5, 1.0, 1.0)
    flow = lm.SnapshotFlow(times, w["ux"], w["uy"], w["f"], w["cs"], domain,
                           window_start=times[0], window_len=0.0)
    f0, cs0, c20 = w["f"][0], w["cs"][0], w["c2"][0]
    print("release c2: min %.3g max %.3g mean(liquid) %.3g" %
          (c20.min(), c20.max(), c20[(f0 >= .5) & (cs0 >= .5)].mean()), flush=True)

    summary = dict(run=run.name, code_units=cu, n=n, snapshots=S, dt_snap=dt_snap,
                   t_release=float(times[0]), periods=args.periods, substeps=args.substeps)

    # divergence diagnostic on snapshots spread over the first period
    dd = []
    for k in np.linspace(0, min(S - 1, int(100)), 6).astype(int):
        d = lm.divergence_diagnostic(w["ux"][k], w["uy"][k], w["f"][k], w["cs"][k], 1.0 / n)
        dd.append(d)
    summary["divergence"] = dict(per_snapshot=dd,
                                 ratio_mean=float(np.mean([d["ratio"] for d in dd])))
    print("divergence ratio rms(div)/rms(|grad u|):", [round(d["ratio"], 4) for d in dd], flush=True)

    # snapshot c2 coarse-grained, chi vs time
    wt = (w["f"] * w["cs"]).astype(np.float64)
    snap_chi = {nb: lm.blockmean_chi_series(w["c2"].astype(np.float64), wt, nb)
                for nb in (16, 32, 64)}
    rel_t = times - times[0]

    steps_per_snap = args.substeps
    dt = dt_snap / args.substeps
    n_steps = int(round(args.periods * T / dt))
    n_steps -= n_steps % steps_per_snap
    K = {16: args.nb // 16, 32: args.nb // 32, 64: args.nb // 64}
    summary["runs"] = {}
    for npc in args.npc:
        key = jax.random.PRNGKey(npc)
        x, y = lm.init_liquid_particles(key, f0, cs0, npc, domain)
        lab_tr = lm.labels_from_field(c20, x, y, domain)
        lab_rand = (np.random.default_rng(npc).random(x.size) < 0.5).astype(np.float32)
        kw = dict(dt=dt, n_steps=n_steps, D=cu["D_nd"], record_every=steps_per_snap,
                  nbx=args.nb, nby=args.nb, key=jax.random.PRNGKey(7))
        t0 = time.time()
        o = lm.run_particles(flow, x, y, lab_tr, **kw)
        o["x"].block_until_ready()
        wall1 = time.time() - t0
        t0 = time.time()
        o = lm.run_particles(flow, x, y, lab_tr, **kw)
        o["x"].block_until_ready()
        wall = time.time() - t0
        orand = lm.run_particles(flow, x, y, lab_rand, **kw)
        t_rec = o["t_rec"]
        r = dict(particles=int(x.size), wall_first_s=wall1, wall_s=wall,
                 wall_s_per_period=wall / args.periods)
        # leakage per period
        per = np.minimum((t_rec // T).astype(int), args.periods - 1)
        l1 = [int(o["leak1"][(per == p) & (t_rec > 0)].sum()) for p in range(args.periods)]
        l2 = [int(o["leak2"][(per == p) & (t_rec > 0)].sum()) for p in range(args.periods)]
        r["leak1_per_period"], r["leak2_per_period"] = l1, l2
        r["leak1_per_particle_per_period_mean"] = float(np.mean(l1) / x.size)
        # uniform checks at l = L0/32 and L0/16: fully liquid boxes at the snapshot
        # nearest in time; one-sided (counts start at zero variance, see below)
        uni = {}
        for nb in (16, 32):
            k = args.nb // nb
            zs, dis = [], []
            for p in range(1, args.periods + 1):
                ridx = int(np.argmin(np.abs(t_rec - p * T)))
                sidx = int(np.argmin(np.abs(rel_t - t_rec[ridx])))
                full = ((w["f"][sidx] > .99) & (w["cs"][sidx] > .99)).reshape(nb, n // nb, nb, n // nb)
                full = full.all(axis=(1, 3))
                Nb = lm._coarsen(o["N"][ridx:ridx + 1], k)[0]
                # expected count of a box: n/cell x cells; Poisson reference mean
                Nsel = Nb.T[full]  # N is (ix, iy); full is (j=y, i=x)
                zs.append(lm.dispersion_z(Nsel))
                dis.append(float(Nsel.var() / Nsel.mean()))
            uni[f"l=L0/{nb}"] = dict(z_end_of_period=zs, dispersion_index=dis)
        r["uniform_density_check"] = uni
        # random-label excess variance relative to cbar(1-cbar) (should be ~0)
        cb = orand["S"][0].sum() / orand["N"][0].sum()
        for nb in (16, 32):
            v = lm.excess_variance_series(orand["N"], orand["S"], args.nb // nb) / (cb * (1 - cb))
            r[f"random_label_excess_var_over_pq_l{nb}"] = dict(
                mean=float(v.mean()), max_abs=float(np.abs(v).max()))
        # chi vs snapshot c2 coarse-grained
        cmp = {}
        for nb in (16, 32, 64):
            chi_p = lm.chi_series(o["N"], o["S"], coarsen=K[nb])
            chi_s = np.interp(t_rec, rel_t, snap_chi[nb])
            cmp[f"l=L0/{nb}"] = dict(
                particles_end_of_period=[float(chi_p[int(np.argmin(np.abs(t_rec - p * T)))])
                                         for p in range(1, args.periods + 1)],
                snapshot_end_of_period=[float(chi_s[int(np.argmin(np.abs(t_rec - p * T)))])
                                        for p in range(1, args.periods + 1)],
                rms_diff=float(np.sqrt(np.mean((chi_p - chi_s) ** 2))),
                max_abs_diff=float(np.abs(chi_p - chi_s).max()))
            if nb == 32:
                np.savez(out / f"chi_l32_npc{npc}.npz", t=t_rec / T, chi_particles=chi_p,
                         chi_snapshot=chi_s)
        r["chi_vs_snapshot_c2"] = cmp
        summary["runs"][str(npc)] = r
        print(json.dumps({str(npc): r}, indent=1), flush=True)
        json.dump(summary, open(out / "summary.json", "w"), indent=1)

    # figure: chi at l = L0/32, particles vs snapshot c2
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 3.2))
    for npc in args.npc:
        d = np.load(out / f"chi_l32_npc{npc}.npz")
        ax.plot(d["t"], d["chi_particles"], lw=1, label=f"particles, {npc} per cell")
    ax.plot(d["t"], d["chi_snapshot"], "k--", lw=1, label="c2 coarse-grained")
    ax.set_xlabel("time since release (rocking periods)")
    ax.set_ylabel(r"$\chi_l$, $l = L_0/32$")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.savefig(out / "chi_l32.png", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
