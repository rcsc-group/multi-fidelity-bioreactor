"""Correctness benchmarks for scripts/lagrangian_mix.py (CPU, small N).

(a) uniform density is invariant under a divergence-free flow (double gyre);
(b) a synthetic NON-divergence-free cell-centred field clusters particles;
(c) pure diffusion in a reflecting box: variance decays at the analytic rate;
(d) reader round-trip on a synthetic snapshot file.
"""
import os
import struct
import sys
from pathlib import Path

os.environ.setdefault("JAX_PLATFORMS", "cpu")

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
import jax  # noqa: E402
from scripts import lagrangian_mix as lm  # noqa: E402


def _dispersion_z(counts):
    """z-score of the index of dispersion var/mean of box counts vs Poisson (=1)."""
    c = np.asarray(counts, float).ravel()
    B = c.size
    d = ((c - c.mean()) ** 2).sum() / c.mean()  # ~ chi2(B-1) for Poisson counts
    return (d - (B - 1)) / np.sqrt(2.0 * (B - 1))


# --------------------------------------------------------------------------- (a)
def test_a_uniform_density_invariant_double_gyre():
    T = 10.0
    flow = lm.AnalyticFlow("double_gyre", (0.1, 0.25, 2 * np.pi / T), (0.0, 0.0, 2.0, 1.0))
    key = jax.random.PRNGKey(1)
    k0, k1 = jax.random.split(key)
    n = 100000
    x, y = lm.init_uniform_box(k0, n, flow.domain)
    label = np.zeros(n, np.float32)
    dt = 0.05
    per = int(round(T / dt))
    out = lm.run_particles(flow, x, y, label, dt=dt, n_steps=10 * per, D=0.0,
                           record_every=per, nbx=20, nby=10, key=k1)
    N = np.asarray(out["N"])  # (11, 20, 10)
    assert N.shape == (11, 20, 10)
    z = np.array([_dispersion_z(N[r]) for r in range(N.shape[0])])
    print("double gyre dispersion z per period:", np.round(z, 2))
    assert np.all(np.abs(z) < 3.0)


# --------------------------------------------------------------------------- (b)
def test_b_non_divergence_free_cell_centred_field_clusters():
    n = 64
    xc = (np.arange(n) + 0.5) / n
    X, Y = np.meshgrid(xc, xc)  # (j=y, i=x)
    u = 0.3 * np.sin(2 * np.pi * X)  # converges on x = 0.5, div u != 0
    v = np.zeros_like(u)
    ones = np.ones((2, n, n), np.float32)
    flow = lm.SnapshotFlow(times=np.array([0.0, 1000.0], np.float32),
                           u=np.stack([u, u]).astype(np.float32),
                           v=np.stack([v, v]).astype(np.float32),
                           f=ones, cs=ones, domain=(0.0, 0.0, 1.0, 1.0))
    key = jax.random.PRNGKey(2)
    k0, k1 = jax.random.split(key)
    N0 = 100000
    x, y = lm.init_uniform_box(k0, N0, flow.domain)
    dt = 0.05
    out = lm.run_particles(flow, x, y, np.zeros(N0, np.float32), dt=dt, n_steps=200, D=0.0,
                           record_every=200, nbx=20, nby=10, key=k1)
    N = np.asarray(out["N"])
    z0, z1 = _dispersion_z(N[0]), _dispersion_z(N[-1])
    print("clustering z initial/final:", z0, z1)
    assert abs(z0) < 3.0
    assert z1 > 50.0


# --------------------------------------------------------------------------- (c)
def test_c_pure_diffusion_reflecting_box_first_mode_rate():
    L, D = 1.0, 1.0e-3
    flow = lm.AnalyticFlow("zero", (), (0.0, 0.0, L, L))
    key = jax.random.PRNGKey(3)
    k0, k1 = jax.random.split(key)
    n = 300000
    x, y = lm.init_uniform_box(k0, n, flow.domain)
    label = (np.asarray(x) < 0.5 * L).astype(np.float32)
    dt, rec = 0.5, 10
    n_steps = 300
    out = lm.run_particles(flow, x, y, label, dt=dt, n_steps=n_steps, D=D,
                           record_every=rec, nbx=8, nby=8, key=k1)
    var = lm.excess_variance_series(out["N"], out["S"])  # Poisson-corrected
    t = np.asarray(out["t_rec"])
    rate_theory = np.pi ** 2 * D / L ** 2  # amplitude (std) rate; variance rate = 2x
    # fit window: higher modes (rate 9x) are gone, noise still small
    m = (t >= 40.0) & (t <= 150.0)
    slope = np.polyfit(t[m], 0.5 * np.log(var[m]), 1)[0]
    rate = -slope
    print("fitted std rate", rate, "theory", rate_theory, "ratio", rate / rate_theory)
    assert abs(rate / rate_theory - 1.0) < 0.05


# --------------------------------------------------------------------------- (d)
def _write_snapshot(path, n, t, planes):
    with open(path, "wb") as fp:
        fp.write(struct.pack("<iid", n, len(planes), t))
        fp.write(np.asarray(planes, np.float32).tobytes())


def test_d_reader_round_trip(tmp_path):
    rng = np.random.default_rng(0)
    n, npl = 16, 10
    snaps = {}
    d = tmp_path / "fields"
    d.mkdir()
    for k, t in enumerate([0.1, 0.35, 0.9]):
        pl = rng.standard_normal((npl, n, n)).astype(np.float32)
        snaps[k] = (t, pl)
        _write_snapshot(d / f"snap_{k:06d}.bin", n, t, pl)
    t, pl = lm.read_snapshot(d / "snap_000001.bin")
    assert t == pytest.approx(0.35)
    assert set(pl) == set(lm.PLANES)
    for i, name in enumerate(lm.PLANES):
        np.testing.assert_array_equal(pl[name], snaps[1][1][i])
    w = lm.load_window(d, planes=("ux", "uy", "f", "cs"))
    np.testing.assert_allclose(w["times"], [0.1, 0.35, 0.9])
    assert w["ux"].shape == (3, n, n)
    np.testing.assert_array_equal(w["f"][2], snaps[2][1][3])
    # origin row (j = 0) is y = Y0 bottom: layout is row-major (j, i)
    assert w["n"] == n


# ------------------------------------------------------------ extras (not benchmarks)
def test_time_interpolation_and_periodic_window():
    n = 8
    times = np.array([0.0, 0.3, 1.0, 1.6, 2.0])  # irregular
    u = np.stack([np.full((n, n), 2.0 * t) for t in times]).astype(np.float32)  # u = 2t
    z = np.zeros_like(u)
    o = np.ones_like(u)
    flow = lm.SnapshotFlow(times, u, z, o, o, (0.0, 0.0, 1.0, 1.0), window_start=0.0, window_len=2.0)
    for s, expect in [(0.15, 0.3), (0.7, 1.4), (1.9, 3.8), (2.5, 1.0), (4.25, 0.5)]:
        ux, _ = flow.frame(np.float32(s)).vel(jnp_arr(0.5), jnp_arr(0.5))
        assert float(ux[0]) == pytest.approx(expect, abs=1e-5)


def jnp_arr(v):
    import jax.numpy as jnp
    return jnp.full((1,), v, jnp.float32)


def test_code_units_matches_formula():
    p = {"geometry": {"a": 0.25, "b": 0.03575}, "theta_max": [7.0], "omega_b": 2.345722514680379}
    cu = lm.code_units(p)
    assert cu["T_nd"] == pytest.approx(cu["T_per"] * cu["U_bio"] / 0.25)
    assert cu["D_nd"] == pytest.approx(0.44e-9 / (cu["U_bio"] * 0.25))
