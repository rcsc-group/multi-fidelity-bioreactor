"""GPU smoke/timing test of lagrangian_mix on a synthetic L8-sized flow
(256^2 cells, 1000 snapshots, liquid band |y| < 0.14, divergence-free
double-gyre-like flow in the band).  Not a benchmark of correctness."""
import json
import sys
import time
from pathlib import Path

import numpy as np
import jax

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts import lagrangian_mix as lm  # noqa: E402

n, S, periods, T = 256, 1000, 10, 0.6073
dt_snap = T / 100.618
xc = -0.5 + (np.arange(n) + 0.5) / n
X, Y = np.meshgrid(xc, xc)
band = (np.abs(Y) < 0.14).astype(np.float32)
times = np.arange(S) * dt_snap
u = np.empty((S, n, n), np.float32)
v = np.empty_like(u)
for k, t in enumerate(times):
    ph = 2 * np.pi * t / T
    psi_x = 0.05 * np.sin(2 * np.pi * (X + 0.3 * np.sin(ph))) * np.cos(np.pi * Y / 0.14)
    u[k] = -0.3 * np.sin(np.pi * Y / 0.28 * 2) * np.cos(2 * np.pi * (X + 0.3 * np.sin(ph))) * band
    v[k] = 0.3 * np.cos(np.pi * Y / 0.28 * 2) * 0.0 + 0.0
f = np.broadcast_to(band, (S, n, n))
flow = lm.SnapshotFlow(times, u, v, f, np.ones_like(f), (-0.5, -0.5, 1.0, 1.0))
print("devices", jax.devices(), flush=True)
res = {}
for npc in (16, 64):
    key = jax.random.PRNGKey(0)
    x, y = lm.init_liquid_particles(key, band, np.ones_like(band), npc, flow.domain)
    lab = np.asarray(x > 0).astype(np.float32)
    kw = dict(dt=dt_snap / 2, n_steps=2 * 100 * periods, D=3.1e-8, record_every=2,
              nbx=64, nby=64, key=key)
    for rep in range(2):
        t0 = time.time()
        out = lm.run_particles(flow, x, y, lab, **kw)
        out["x"].block_until_ready()
        dtw = time.time() - t0
        print(f"N/cell {npc} particles {x.size} rep {rep} wall {dtw:.1f} s", flush=True)
    res[npc] = dict(particles=int(x.size), wall_s_10_periods=dtw, wall_s_per_period=dtw / periods,
                    leak1_total=int(out["leak1"].sum()), leak2_total=int(out["leak2"].sum()))
print(json.dumps(res, indent=1))
