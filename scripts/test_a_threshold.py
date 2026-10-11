"""Detection threshold of benchmark (a): double gyre + weak compressive part
(kappa); reports div/|grad u| (rms over the box, t = 0 and averaged over a
period) and the dispersion statistics after 10 periods."""
import os, sys
os.environ.setdefault("JAX_PLATFORMS", "cpu")
from pathlib import Path
import numpy as np
from scipy import stats
import jax
sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts import lagrangian_mix as lm

T, A, eps = 10.0, 0.1, 0.25
om = 2 * np.pi / T
dom = (0.0, 0.0, 2.0, 1.0)


def ratio(kappa):
    n = 400
    x = (np.arange(2 * n) + .5) / n; y = (np.arange(n) + .5) / n
    X, Y = np.meshgrid(x, y)
    r = []
    for t in np.linspace(0, T, 11)[:-1]:
        u, v = [np.asarray(a) for a in lm._double_gyre_compressive((A, eps, om, kappa), X, Y, t)]
        h = 1.0 / n
        ux = np.gradient(u, h, axis=1); uy = np.gradient(u, h, axis=0)
        vx = np.gradient(v, h, axis=1); vy = np.gradient(v, h, axis=0)
        r.append((np.mean((ux + vy) ** 2), np.mean(ux**2 + uy**2 + vx**2 + vy**2)))
    r = np.mean(r, axis=0)
    return float(np.sqrt(r[0] / r[1]))


def run(kappa, seed, n=100000, dt=0.05):
    k0, k1 = jax.random.split(jax.random.PRNGKey(seed))
    x, y = lm.init_uniform_box(k0, n, dom)
    fl = lm.AnalyticFlow("double_gyre_compressive", (A, eps, om, kappa), dom)
    per = int(round(T / dt))
    o = lm.run_particles(fl, x, y, np.zeros(n, np.float32), dt=dt, n_steps=10 * per, D=0.0,
                         record_every=10 * per, nbx=20, nby=10, key=k1)
    N = o["N"].astype(float)
    out = []
    for r in (0, -1):
        c = N[r].ravel(); m = c.mean()
        z = (c - m) / np.sqrt(m)
        out.append((lm.dispersion_z(c), stats.kstest(z, "norm").pvalue))
    return out


if __name__ == "__main__":
    print("kappa  div/|grad u|  z(t=0)  z(10T)  KSp(0) KSp(10T)  [3 seeds]")
    for kappa in [0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]:
        rs = ratio(kappa)
        res = [run(kappa, s) for s in range(3)]
        print(f"{kappa:7.0e} {rs:.2e} ", " | ".join(
            f"{a[0][0]:+.1f} {a[1][0]:+.1f} {a[0][1]:.2f} {a[1][1]:.3f}" for a in res), flush=True)
