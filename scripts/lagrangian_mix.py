"""Lagrangian (particle) estimator of mixing from recorded flow snapshots.

Library only; the driver for the real run is scripts/run_lagrangian_real.py.
Implements docs/dfd_plan.md section "Estimator": passive particles advected
(RK2 midpoint) in the recorded, time-interpolated velocity field plus a
Brownian step at prescribed diffusivity D, confined to the liquid; labels give
a tracer concentration; chi_l(t) is computed from particle counts in boxes of
side l with explicit Poisson shot-noise subtraction.

--------------------------------------------------------------------------
Facts read from src/BioReactor.c (not recalled)
--------------------------------------------------------------------------
* Snapshot file (write_snapshot, ~line 1891): int32 n, int32 np, float64 t_nd,
  then np planes of n*n float32, row-major (j = y, i = x), in the order
  ux, uy, omega, f, cs, c, c1, c2, c3, oxy.  Header is 16 bytes, no padding.
* Cell centres: x_i = X0 + (i + 0.5) dx, y_j = Y0 + (j + 0.5) dx, dx = L0/n,
  with L0 = 1 and origin(-L0/2, -L0/2) (main(), ~line 389): X0 = Y0 = -0.5.
  Lengths are in units of L_bio = geometry.a (the bag half-length, m).
* Frame: the velocity is in the rocking (tank) frame.  The tank is fixed in
  the grid (cs is static) and the rocking enters only as the fictitious
  acceleration (event acceleration(), line ~1670, with Coriolis/centrifugal/
  Euler terms).  So particles are advected with u as recorded, no frame
  transformation.
* Non-dimensionalisation (main(), lines ~399-426), with a = geometry.a,
  b = geometry.b, theta = theta_max[0] (rad), T_per = 2 pi / omega_b:
      H_bio = 2 b;  V_bio = (a/4) (H_bio + 0.5 a tan theta)
      U_bio = V_bio / (H_bio / 2) / T_per;  T_bio = a / U_bio
  Lengths / a, times / T_bio, velocities / U_bio.  Hence
      period (code time units)   T_nd = T_per / T_bio = T_per U_bio / a
      tracer diffusivity (code)  D_nd = 1 / Pe_tracer = D / (U_bio a),
  with D = D_tracer_1 = 0.44e-9 m^2/s (water; line 201; c.D1 = 1/Pe_tracer_1).
  See code_units().
* Snapshots are written every dt_video = T_nd / (N + 0.618...) (incommensurate
  with the period); this module interpolates linearly in time between the two
  bracketing snapshots, so no phase resampling is done.

--------------------------------------------------------------------------
Concentration estimator
--------------------------------------------------------------------------
Boxes b of side l tile the domain.  N_b = particles in b, S_b = tracer-labelled
particles in b (label fixed per particle), c_b = S_b / N_b over boxes with
N_b >= 2, weights w_b = N_b (particles are uniform in liquid volume, so this is
the liquid-volume weight), cbar = sum S_b / sum N_b.
    var_obs = sum_b w_b (c_b - cbar)^2 / sum_b w_b
    shot    = sum_b w_b c_b (1 - c_b) / (N_b - 1) / sum_b w_b
    var_l   = var_obs - shot
(c_b(1-c_b)/(N_b-1) is the unbiased estimate of the binomial sampling variance
of a box fraction from N_b particles; it is the Poisson shot noise for 0/1
labels.)  chi_l(t) = 1 - var_l(t) / var_l(t0), t0 = release (box mode,
default), or 1 - var_l(t) / (cbar (1 - cbar)) (binary mode, the cell-scale
initial variance of a 0/1 field).

Uniform-density check: index of dispersion var(N_b)/mean(N_b) of the box counts
is 1 for Poisson counts; (sum (N_b-m)^2/m - (B-1)) / sqrt(2 (B-1)) is the
z-score (chi^2 with B-1 degrees of freedom).

--------------------------------------------------------------------------
Confinement
--------------------------------------------------------------------------
A step x -> x_adv + sqrt(2 D dt) xi is valid if the landing cell has
(time-interpolated) f >= 0.5 and cs >= 0.5 (and lies inside the domain).
Event 1: the first try is invalid (COUNTED: leak1).  Retry once with a fresh
Brownian draw; if still invalid (COUNTED: leak2) the first landing point is
projected to the nearest liquid cell centre within a 7x7 cell window; if there
is none the particle stays where it was.  Leakage per period = leak1 summed
over the steps of a period (also reported per particle).

Choices where the spec was silent
--------------------------------------------------------------------------
1. "Reject the Brownian part and retry once": retry = a fresh Brownian draw
   added to the same advected point; the projection fallback acts on the
   first landing point.  leak1 counts first-try failures, leak2 counts
   projections (both returned).
2. Validity lookup uses the cell containing the landing point, with f and cs
   linearly interpolated in time (not in space).
3. Velocity is bilinear in space between cell centres with constant
   extrapolation outside the outermost centres; velocity of air/solid cells is
   used as recorded (no masking), so particles near the interface see a mix.
4. Time stepping: dt = dt_snapshot / substeps; chi recorded at fixed step
   multiples; comparisons with snapshots interpolate the particle chi to the
   snapshot times.  Beyond the recorded time range the flow is held at the
   first/last snapshot.
5. Periodic reuse is a window of m whole periods, s -> t_start + (s mod m T);
   m = 1 repeats one period, m > 1 cycles through m different recorded
   periods (window_periods()).  The flow is not exactly periodic, so there is
   a discontinuity at each wrap (measured by comparing different m).
6. Particle initialisation: n per liquid cell (f >= 0.5 and cs >= 0.5 at
   release), uniform within the cell (cells with 0.5 <= f < 1 get the same n).
   Label = (c2 > 0.5) in the cell containing the particle.
7. chi_l reference variance is the variance of the same estimator at release
   (so chi_l(0) = 0 at every l); 'binary' mode is available.  Boxes with
   N_b < 2 are dropped.
8. Pure-diffusion test: the variance of the box fractions decays at twice the
   amplitude rate; the test fits the standard deviation (sqrt var), whose
   analytic rate is pi^2 D / L^2 (variance rate 2 pi^2 D / L^2).
9. Snapshot comparison: block mean of c2 weighted by w = f * cs (cs also
   removes solid cells), blocks of n/nbox cells; the same box-variance chi
   with its own release as reference (the shot-noise term is zero).
10. Divergence diagnostic: central differences, liquid mask f > 0.99 and
    cs > 0.99 eroded by `erode` cells (default 3) in every direction;
    |grad u|^2 = ux_x^2 + ux_y^2 + uy_x^2 + uy_y^2 on the same mask.
11. float32 everywhere (JAX default), no x64.
"""
from __future__ import annotations

import json
import math
import struct
from dataclasses import dataclass
from functools import partial
from pathlib import Path

import numpy as np
import jax
import jax.numpy as jnp

PLANES = ("ux", "uy", "omega", "f", "cs", "c", "c1", "c2", "c3", "oxy")
_HEADER = 16  # int32 n, int32 np, float64 t


# ----------------------------------------------------------------------------
# Reader
# ----------------------------------------------------------------------------
def read_snapshot(path, planes=None):
    """Return (t_nd, {name: array (n, n)}) for one snapshot file."""
    with open(path, "rb") as fp:
        n, npl = struct.unpack("<ii", fp.read(8))
        (t,) = struct.unpack("<d", fp.read(8))
    names = PLANES[:npl]
    mm = np.memmap(path, dtype="<f4", mode="r", offset=_HEADER, shape=(npl, n, n))
    want = names if planes is None else planes
    return t, {k: np.array(mm[names.index(k)]) for k in want}


def snapshot_files(fields_dir):
    return sorted(Path(fields_dir).glob("snap_*.bin"))


def load_window(fields_dir, planes=("ux", "uy", "f", "cs"), t_min=None, t_max=None,
                files=None):
    """Load all snapshots with t_min <= t <= t_max.  Returns dict with 'times' (S,),
    one (S, n, n) float32 array per plane, and 'n'."""
    files = snapshot_files(fields_dir) if files is None else files
    times, data, n = [], {k: [] for k in planes}, None
    for p in files:
        with open(p, "rb") as fp:
            n_, _ = struct.unpack("<ii", fp.read(8))
            (t,) = struct.unpack("<d", fp.read(8))
        if (t_min is not None and t < t_min) or (t_max is not None and t > t_max):
            continue
        t, pl = read_snapshot(p, planes)
        n = n_
        times.append(t)
        for k in planes:
            data[k].append(pl[k])
    out = {k: np.stack(v) for k, v in data.items()}
    out["times"] = np.array(times, dtype=np.float64)
    out["n"] = n
    return out


# ----------------------------------------------------------------------------
# Units
# ----------------------------------------------------------------------------
def code_units(params, D_phys=0.44e-9):
    """Code-unit scales from a params.json dict (formulas of src/BioReactor.c)."""
    a = params["geometry"]["a"]
    b = params["geometry"]["b"]
    th = math.radians(params["theta_max"][0])
    T_per = 2 * math.pi / params["omega_b"]
    H = 2 * b
    V = a / 4 * (H + 0.5 * a * math.tan(th))
    U = V / (H / 2) / T_per
    T_bio = a / U
    return dict(L_bio=a, U_bio=U, T_bio=T_bio, T_per=T_per, T_nd=T_per / T_bio,
                D_phys=D_phys, D_nd=D_phys / (U * a), Pe=U * a / D_phys)


def window_periods(times, T_nd, t_start, m):
    """Window (t_start, length m T_nd) of m whole periods; checks the snapshots
    cover it (one snapshot beyond each end)."""
    t_end = t_start + m * T_nd
    if times[0] > t_start or times[-1] < t_end - 1e-6 * T_nd:
        raise ValueError("snapshots do not cover the requested window of %d periods" % m)
    return float(t_start), float(m * T_nd)


# ----------------------------------------------------------------------------
# Frames (a flow evaluated at one instant) and flows (pytrees passed to jit)
# ----------------------------------------------------------------------------
class _GridFrame:
    def __init__(self, u, v, liq, domain):
        self.u, self.v, self.liq = u, v, liq
        self.x0, self.y0, self.Lx, self.Ly = domain
        self.ny, self.nx = u.shape
        self.dx, self.dy = self.Lx / self.nx, self.Ly / self.ny

    def vel(self, x, y):
        gx = (x - self.x0) / self.dx - 0.5
        gy = (y - self.y0) / self.dy - 0.5
        i0 = jnp.clip(jnp.floor(gx), 0, self.nx - 2)
        j0 = jnp.clip(jnp.floor(gy), 0, self.ny - 2)
        wx = jnp.clip(gx - i0, 0.0, 1.0)
        wy = jnp.clip(gy - j0, 0.0, 1.0)
        i0 = i0.astype(jnp.int32)
        j0 = j0.astype(jnp.int32)

        def bil(a):
            return ((1 - wy) * ((1 - wx) * a[j0, i0] + wx * a[j0, i0 + 1])
                    + wy * ((1 - wx) * a[j0 + 1, i0] + wx * a[j0 + 1, i0 + 1]))

        return bil(self.u), bil(self.v)

    def _cell(self, x, y):
        i = jnp.clip(jnp.floor((x - self.x0) / self.dx), 0, self.nx - 1).astype(jnp.int32)
        j = jnp.clip(jnp.floor((y - self.y0) / self.dy), 0, self.ny - 1).astype(jnp.int32)
        return i, j

    def valid(self, x, y):
        inside = ((x >= self.x0) & (x < self.x0 + self.Lx)
                  & (y >= self.y0) & (y < self.y0 + self.Ly))
        i, j = self._cell(x, y)
        return inside & self.liq[j, i]

    def project(self, xl, yl, xold, yold):
        """Nearest liquid cell centre in a 7x7 window around the landing cell;
        particles with none stay at (xold, yold)."""
        i, j = self._cell(xl, yl)
        best = jnp.full(xl.shape, jnp.inf)
        bx, by = xold, yold
        for dj in range(-3, 4):
            for di in range(-3, 4):
                ii = jnp.clip(i + di, 0, self.nx - 1)
                jj = jnp.clip(j + dj, 0, self.ny - 1)
                cx = self.x0 + (ii + 0.5) * self.dx
                cy = self.y0 + (jj + 0.5) * self.dy
                d = (cx - xl) ** 2 + (cy - yl) ** 2
                ok = self.liq[jj, ii] & (d < best)
                best = jnp.where(ok, d, best)
                bx = jnp.where(ok, cx, bx)
                by = jnp.where(ok, cy, by)
        return bx, by


def _double_gyre(params, x, y, t):
    A, eps, om = params
    a = eps * jnp.sin(om * t)
    b = 1.0 - 2.0 * a
    fx = a * x * x + b * x
    dfx = 2.0 * a * x + b
    u = -jnp.pi * A * jnp.sin(jnp.pi * fx) * jnp.cos(jnp.pi * y)
    v = jnp.pi * A * jnp.cos(jnp.pi * fx) * jnp.sin(jnp.pi * y) * dfx
    return u, v


def _double_gyre_compressive(params, x, y, t):
    """Double gyre plus u_x += kappa sin(pi x) (zero normal flux at x = 0, 2);
    div = kappa pi cos(pi x)."""
    A, eps, om, kappa = params
    u, v = _double_gyre((A, eps, om), x, y, t)
    return u + kappa * jnp.sin(jnp.pi * x), v


class _AnalyticFrame:
    def __init__(self, kind, params, domain, t):
        self.kind, self.params, self.t = kind, params, t
        self.x0, self.y0, self.Lx, self.Ly = domain

    def vel(self, x, y):
        if self.kind == "double_gyre":
            return _double_gyre(self.params, x, y, self.t)
        if self.kind == "double_gyre_compressive":
            return _double_gyre_compressive(self.params, x, y, self.t)
        return jnp.zeros_like(x), jnp.zeros_like(y)

    def valid(self, x, y):
        return ((x >= self.x0) & (x <= self.x0 + self.Lx)
                & (y >= self.y0) & (y <= self.y0 + self.Ly))

    def project(self, xl, yl, xold, yold):
        return (jnp.clip(xl, self.x0, self.x0 + self.Lx),
                jnp.clip(yl, self.y0, self.y0 + self.Ly))


@jax.tree_util.register_pytree_node_class
@dataclass(frozen=True)
class AnalyticFlow:
    """Analytic flow in a rectangular box (reflecting walls).
    kind 'double_gyre': params (A, eps, omega), domain [0,2]x[0,1];
    kind 'zero': no flow."""
    kind: str
    params: tuple
    domain: tuple

    def frame(self, s):
        return _AnalyticFrame(self.kind, self.params, self.domain, s)

    def tree_flatten(self):
        return (), (self.kind, self.params, self.domain)

    @classmethod
    def tree_unflatten(cls, aux, children):
        return cls(*aux)


@jax.tree_util.register_pytree_node_class
class SnapshotFlow:
    """Cell-centred snapshots (S, ny, nx) of u, v, f, cs at `times` (S,),
    linearly interpolated in time.  Simulation time s starts at window_start
    (default times[0]); if window_len > 0 the flow is reused periodically:
    t = window_start + (s mod window_len).  Out-of-range t is clamped."""

    def __init__(self, times, u, v, f, cs, domain, window_start=None, window_len=0.0):
        times = np.asarray(times, np.float64)
        self.window_start = float(times[0] if window_start is None else window_start)
        self.window_len = float(window_len)
        self.domain = tuple(float(d) for d in domain)
        self.times = jnp.asarray(times - self.window_start, jnp.float32)
        self.u, self.v = jnp.asarray(u, jnp.float32), jnp.asarray(v, jnp.float32)
        self.f, self.cs = jnp.asarray(f, jnp.float32), jnp.asarray(cs, jnp.float32)

    def frame(self, s):
        t = jnp.mod(s, self.window_len) if self.window_len > 0 else s
        S = self.times.shape[0]
        k = jnp.clip(jnp.searchsorted(self.times, t, side="right") - 1, 0, S - 2)
        t0, t1 = self.times[k], self.times[k + 1]
        a = jnp.clip((t - t0) / (t1 - t0), 0.0, 1.0)

        def lerp(A):
            return (1 - a) * A[k] + a * A[k + 1]

        liq = (lerp(self.f) >= 0.5) & (lerp(self.cs) >= 0.5)
        return _GridFrame(lerp(self.u), lerp(self.v), liq, self.domain)

    def tree_flatten(self):
        return ((self.times, self.u, self.v, self.f, self.cs),
                (self.domain, self.window_start, self.window_len))

    @classmethod
    def tree_unflatten(cls, aux, ch):
        o = object.__new__(cls)
        o.times, o.u, o.v, o.f, o.cs = ch
        o.domain, o.window_start, o.window_len = aux
        return o


# ----------------------------------------------------------------------------
# Particles
# ----------------------------------------------------------------------------
def init_uniform_box(key, n, domain):
    x0, y0, Lx, Ly = domain
    kx, ky = jax.random.split(key)
    return (x0 + Lx * jax.random.uniform(kx, (n,)), y0 + Ly * jax.random.uniform(ky, (n,)))


def init_liquid_particles(key, f, cs, n_per_cell, domain):
    """n_per_cell particles uniform in every liquid cell (f >= 0.5, cs >= 0.5)."""
    x0, y0, Lx, Ly = domain
    ny, nx = f.shape
    dx, dy = Lx / nx, Ly / ny
    jj, ii = np.nonzero((f >= 0.5) & (cs >= 0.5))
    jj, ii = np.repeat(jj, n_per_cell), np.repeat(ii, n_per_cell)
    kx, ky = jax.random.split(key)
    rx = np.asarray(jax.random.uniform(kx, (jj.size,)))
    ry = np.asarray(jax.random.uniform(ky, (jj.size,)))
    return (jnp.asarray(x0 + (ii + rx) * dx, jnp.float32),
            jnp.asarray(y0 + (jj + ry) * dy, jnp.float32))


def labels_from_field(c, x, y, domain, thresh=0.5):
    """Label 1.0 where the cell containing the particle has c > thresh."""
    x0, y0, Lx, Ly = domain
    ny, nx = c.shape
    i = np.clip(np.floor((np.asarray(x) - x0) / (Lx / nx)), 0, nx - 1).astype(int)
    j = np.clip(np.floor((np.asarray(y) - y0) / (Ly / ny)), 0, ny - 1).astype(int)
    return (np.asarray(c)[j, i] > thresh).astype(np.float32)


def _counts(x, y, label, domain, nbx, nby):
    x0, y0, Lx, Ly = domain
    ix = jnp.clip(jnp.floor((x - x0) / Lx * nbx), 0, nbx - 1).astype(jnp.int32)
    iy = jnp.clip(jnp.floor((y - y0) / Ly * nby), 0, nby - 1).astype(jnp.int32)
    idx = ix * nby + iy
    N = jnp.bincount(idx, length=nbx * nby).reshape(nbx, nby)
    S = jnp.bincount(idx, weights=label, length=nbx * nby).reshape(nbx, nby)
    return N, S


@partial(jax.jit, static_argnames=("dt", "n_steps", "D", "record_every", "nbx", "nby"))
def _run(flow, x, y, label, key, t0, *, dt, n_steps, D, record_every, nbx, nby):
    domain = flow.domain
    sig = math.sqrt(2.0 * D * dt)
    nrec = n_steps // record_every

    def step(i, carry):
        x, y, e1, e2, istep = carry
        s = t0 + istep * dt
        u1, v1 = flow.frame(s).vel(x, y)
        xm, ym = x + 0.5 * dt * u1, y + 0.5 * dt * v1
        u2, v2 = flow.frame(s + 0.5 * dt).vel(xm, ym)
        xa, ya = x + dt * u2, y + dt * v2
        fr = flow.frame(s + dt)
        if D > 0.0:
            z = jax.random.normal(jax.random.fold_in(key, istep), (4,) + x.shape)
            x1, y1 = xa + sig * z[0], ya + sig * z[1]
            x2, y2 = xa + sig * z[2], ya + sig * z[3]
        else:
            x1, y1, x2, y2 = xa, ya, xa, ya
        ok1 = fr.valid(x1, y1)
        ok2 = fr.valid(x2, y2)
        bad = (~ok1) & (~ok2)
        xp, yp = jax.lax.cond(jnp.any(bad), lambda: fr.project(x1, y1, x, y),
                              lambda: (x1, y1))
        xn = jnp.where(ok1, x1, jnp.where(ok2, x2, xp))
        yn = jnp.where(ok1, y1, jnp.where(ok2, y2, yp))
        return xn, yn, e1 + jnp.sum(~ok1), e2 + jnp.sum(bad), istep + 1

    def block(carry, _):
        x, y, istep = carry
        x, y, e1, e2, istep = jax.lax.fori_loop(0, record_every, step,
                                                (x, y, 0, 0, istep))
        N, S = _counts(x, y, label, domain, nbx, nby)
        return (x, y, istep), (N, S, e1, e2)

    N0, S0 = _counts(x, y, label, domain, nbx, nby)
    (xf, yf, _), (N, S, e1, e2) = jax.lax.scan(block, (x, y, jnp.int32(0)), None, length=nrec)
    zero = jnp.zeros((1,), jnp.int32)
    return dict(x=xf, y=yf,
                N=jnp.concatenate([N0[None], N]), S=jnp.concatenate([S0[None], S]),
                leak1=jnp.concatenate([zero, e1]), leak2=jnp.concatenate([zero, e2]))


def run_particles(flow, x, y, label, *, dt, n_steps, D, record_every, nbx, nby, key, t0=0.0):
    """Advance particles n_steps of size dt (simulation time starts at t0,
    measured from the flow's window start).  Returns dict with final x, y;
    N, S box counts (R+1, nbx, nby) at the initial time and every record_every
    steps; leak1/leak2 (R+1,) event counts in each record interval; t_rec."""
    if n_steps % record_every:
        raise ValueError("n_steps must be a multiple of record_every")
    out = _run(flow, jnp.asarray(x, jnp.float32), jnp.asarray(y, jnp.float32),
               jnp.asarray(label, jnp.float32), key, jnp.float32(t0), dt=float(dt),
               n_steps=int(n_steps), D=float(D), record_every=int(record_every),
               nbx=int(nbx), nby=int(nby))
    out = {k: np.asarray(v) if k in ("N", "S", "leak1", "leak2") else v
           for k, v in out.items()}
    out["t_rec"] = t0 + np.arange(out["N"].shape[0]) * record_every * dt
    return out


# ----------------------------------------------------------------------------
# Concentration statistics
# ----------------------------------------------------------------------------
def _coarsen(A, k):
    if k == 1:
        return A
    R, bx, by = A.shape
    return A.reshape(R, bx // k, k, by // k, k).sum(axis=(2, 4))


def excess_variance_series(N, S, coarsen=1):
    """Shot-noise-subtracted box variance var_l per record (see module docstring)."""
    N = _coarsen(np.asarray(N, float), coarsen).reshape(len(N), -1)
    S = _coarsen(np.asarray(S, float), coarsen).reshape(len(S), -1)
    out = np.empty(len(N))
    for r in range(len(N)):
        m = N[r] >= 2
        n, s = N[r][m], S[r][m]
        c = s / n
        cbar = s.sum() / n.sum()
        obs = (n * (c - cbar) ** 2).sum() / n.sum()
        shot = (n * c * (1 - c) / (n - 1)).sum() / n.sum()
        out[r] = obs - shot
    return out


def chi_series(N, S, coarsen=1, sigma0_mode="box"):
    """chi_l(t) = 1 - var_l(t)/var_ref; var_ref = var_l(record 0) ('box') or
    cbar (1 - cbar) at record 0 ('binary')."""
    var = excess_variance_series(N, S, coarsen)
    if sigma0_mode == "box":
        ref = var[0]
    else:
        cb = np.asarray(S[0], float).sum() / np.asarray(N[0], float).sum()
        ref = cb * (1 - cb)
    return 1.0 - var / ref


def dispersion_z(counts):
    """z-score of the index of dispersion of box counts against Poisson."""
    c = np.asarray(counts, float).ravel()
    B = c.size
    d = ((c - c.mean()) ** 2).sum() / c.mean()
    return float((d - (B - 1)) / math.sqrt(2.0 * (B - 1)))


def blockmean_chi_series(c, w, nbox):
    """chi of the snapshot field c (S, n, n) coarse-grained to nbox x nbox
    blocks, block mean weighted by w (e.g. f*cs); weighted variance of the
    block means relative to snapshot 0."""
    S, n, _ = c.shape
    k = n // nbox
    W = w.reshape(S, nbox, k, nbox, k).sum(axis=(2, 4))
    C = (w * c).reshape(S, nbox, k, nbox, k).sum(axis=(2, 4))
    ok = W > 0
    cb = np.where(ok, C / np.where(ok, W, 1), 0.0)
    cbar = C.sum(axis=(1, 2)) / W.sum(axis=(1, 2))
    var = (W * (cb - cbar[:, None, None]) ** 2).sum(axis=(1, 2)) / W.sum(axis=(1, 2))
    return 1.0 - var / var[0]


# ----------------------------------------------------------------------------
# Divergence diagnostic
# ----------------------------------------------------------------------------
def divergence_diagnostic(ux, uy, f, cs, dx, erode=3):
    """rms of the discrete divergence of the cell-centred u over liquid cells
    away from the interface and solid, relative to the rms of |grad u|.
    Central differences; mask = (f > 0.99 & cs > 0.99) eroded `erode` cells."""
    ux, uy = np.asarray(ux, float), np.asarray(uy, float)
    mask = (np.asarray(f) > 0.99) & (np.asarray(cs) > 0.99)
    for _ in range(erode):
        m = mask.copy()
        m[1:, :] &= mask[:-1, :]
        m[:-1, :] &= mask[1:, :]
        m[:, 1:] &= mask[:, :-1]
        m[:, :-1] &= mask[:, 1:]
        mask = m
    mask[0, :] = mask[-1, :] = False
    mask[:, 0] = mask[:, -1] = False

    def ddx(a):
        g = np.zeros_like(a)
        g[:, 1:-1] = (a[:, 2:] - a[:, :-2]) / (2 * dx)
        return g

    def ddy(a):
        g = np.zeros_like(a)
        g[1:-1, :] = (a[2:, :] - a[:-2, :]) / (2 * dx)
        return g

    uxx, uxy, uyx, uyy = ddx(ux), ddy(ux), ddx(uy), ddy(uy)
    div = (uxx + uyy)[mask]
    g2 = (uxx ** 2 + uxy ** 2 + uyx ** 2 + uyy ** 2)[mask]
    rms_div = float(np.sqrt(np.mean(div ** 2))) if div.size else float("nan")
    rms_grad = float(np.sqrt(np.mean(g2))) if div.size else float("nan")
    return dict(n_cells=int(mask.sum()), rms_div=rms_div, rms_grad_u=rms_grad,
                ratio=rms_div / rms_grad if div.size else float("nan"))
