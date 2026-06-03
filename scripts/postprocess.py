"""Extract KPIs from a completed BioReactor simulation run.

What this script does
---------------------
After a simulation finishes, this script reads the raw time-series output files
and computes ten numbers that characterise how well the bioreactor is mixing and
oxygenating the liquid.  These numbers are saved to ``results.json`` in the run
directory and are what appear in the heatmap figures.

Quick-start for interns
-----------------------
To run on a finished simulation::

    python scripts/postprocess.py runs/<run_id>/

The script is also called automatically at the end of every SLURM job
(``config/slurm_template.sh``), so you normally never need to run it manually.

The ten KPIs (Key Performance Indicators)
------------------------------------------
All KPIs live in ``results.json`` next to the simulation outputs.

**Oxygen transfer rate — kLa  (units: 1/non-dim-time)**

kLa (pronounced "kay-el-ay") measures how fast dissolved oxygen accumulates in
the liquid.  A higher kLa means the bioreactor is better at supplying oxygen to
cells.  It is defined by the first-order model::

    dC*(t)/dt = kLa × (1 − C*(t))

where C* is the dimensionless dissolved oxygen level (0 = none, 1 = fully
saturated).  We compute kLa at three saturation crossing points:

- ``kLa_10``  — value of kLa when C* first reaches 10 %
- ``kLa_25``  — value of kLa when C* first reaches 25 %   ← most commonly cited
- ``kLa_50``  — value of kLa when C* first reaches 50 %

Each is estimated two ways:

- **5-point log-linear fit** (``kLa_10``, ``kLa_25``, ``kLa_50``): fits a
  straight line to ln(1−C*) vs time over the 5 data points closest to the
  crossing.  Smoother, less sensitive to noise.

- **Instantaneous** (``kLa_inst_10``, ``kLa_inst_25``, ``kLa_inst_50``):
  computes dC*/dt at the crossing point via finite differences, then divides by
  (1−C*).  Faster but noisier.

Both methods should agree within ~20 %.  Large disagreements suggest a noisy
oxygen curve (e.g., t_buffer too short).

**Mixing time — dtmix  (units: seconds)**

The passive tracer is injected into the top half of the liquid at t_mix.  The
degree of mixing χ(t) rises from 0 (completely segregated) to 1 (fully
homogeneous) as the tracer spreads::

    χ(t) = 1 − σ²(t) / σ²_max

where σ²(t) is the spatial variance of tracer concentration.  We record how
long it takes to reach three levels of homogeneity:

- ``dtmix_0.50``  — time (s) to reach 50 % mixing
- ``dtmix_0.75``  — time (s) to reach 75 % mixing
- ``dtmix_0.95``  — time (s) to reach 95 % mixing

A faster mixer reaches 95 % sooner.  Values of NaN mean the simulation did not
run long enough (increase t_buffer in the sweep config).

**Steady streaming vorticity — vor_mean  (units: 1/s)**

When a bag rocks back and forth, the time-averaged flow forms a slow "streaming"
circulation (two counter-rotating vortices) that is responsible for most of the
long-term mixing and oxygen transport.  ``vor_mean`` is the spatial average of
the absolute vorticity |ξ| of this time-averaged flow, in dimensional units of
1/s.

Higher vor_mean → stronger steady streaming → faster mixing and higher kLa.
It is the hydrodynamic root cause connecting operating conditions to
bioreactor performance.

Input files read
-----------------
All files are in the run directory (``runs/<run_id>/``):

``tr_oxy.dat``
    Time series of spatially integrated quantities.  Columns (0-indexed):

    0 i        — timestep index
    1 t        — non-dimensional simulation time
    2 oxy_liq_sum   — ∫ C_oxy dV over liquid phase (oxygen integral)
    3 oxy_liq_sum2  — ∫ C_oxy² dV (not used here)
    4–5 c_liq_sum/2  — horizontal-left tracer (not used; only set if HORIZONTAL_MIXL=1)
    6–7 c1_liq_*     — horizontal-right tracer (not used)
    8 c2_liq_sum    — ∫ C_tracer dV over liquid, VERTICAL_MIXUP configuration ← **used**
    9 c2_liq_sum2   — ∫ C_tracer² dV                                           ← **used**
    10–11 c3_liq_*   — vertical-down tracer (not used)

``vol_frac_interf.dat``
    Time series of the liquid volume fraction and interface geometry.  The mean
    of column 2 (``f_liq_sum``) gives the total liquid volume, used to
    normalise the oxygen and tracer integrals into mean concentrations.

``normf.dat``
    Time series of spatially averaged velocity and vorticity statistics.
    Column 2 (``Omega_liq_avg``) is the liquid-phase mean absolute vorticity
    in non-dimensional units; we convert to 1/s using T_bio.

``params.json``
    Simulation parameters — used only for dimensional time conversion
    (T_bio = L_bio / U_bio) so that mixing times come out in seconds.

Non-dimensionalisation
----------------------
The simulation runs in dimensionless units.  The characteristic scales are::

    T_per  = 2π / omega_b          [rocking period, seconds]
    V_char = L/4 × (H + ½ L tan θ) [characteristic swept volume, m²]
    U_bio  = V_char / (H/2) / T_per [characteristic velocity, m/s]
    T_bio  = L / U_bio              [characteristic time, seconds]

To convert a non-dimensional time τ to seconds: t_seconds = τ × T_bio.
To convert a non-dimensional vorticity ω to 1/s: ω_dim = ω_nd / T_bio.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


MIN_WINDOW = 5   # minimum number of time points for a valid kLa estimate


# ── time-scale helpers ────────────────────────────────────────────────────────

def _t_scales(params: dict) -> tuple[float, float]:
    """Return (T_bio [s], T_per_nd [-]) from simulation params.

    T_bio is the characteristic dimensional time scale used to convert
    non-dimensional simulation time to real seconds.  T_per_nd is one rocking
    period expressed in those same non-dimensional units.

    Parameters
    ----------
    params : dict
        Contents of params.json.  Reads omega_b, geometry.a, geometry.b,
        and theta_max[0].

    Returns
    -------
    T_bio : float
        Characteristic time [seconds].
    T_per_nd : float
        One rocking period in non-dimensional time units [-].
    """
    omega_b  = params.get("omega_b", 3.93)
    L        = params.get("geometry", {}).get("a", 0.25)
    H        = params.get("geometry", {}).get("b", 0.071)
    th       = math.radians(params.get("theta_max", [7.0])[0])
    T_per    = 2 * math.pi / omega_b
    V        = L / 4 * (H + 0.5 * L * math.tan(th))
    U        = V / (H * 0.5) / T_per
    T_bio    = L / U
    T_per_nd = T_per / T_bio
    return T_bio, T_per_nd


def _load_params(run_dir: Path) -> dict:
    """Load params.json from run_dir, return empty dict if absent."""
    p = run_dir / "params.json"
    return json.loads(p.read_text()) if p.exists() else {}


# ── raw data loaders ──────────────────────────────────────────────────────────

def _load_dat(path: Path) -> np.ndarray:
    """Load a space-delimited Basilisk .dat file into a 2-D numpy array.

    Skips the header line (which starts with 'i') and any blank lines.
    Each subsequent line becomes one row; columns are floats.
    """
    lines = [
        l for l in path.read_text().splitlines()
        if l.strip() and not l.strip().startswith("i")
    ]
    return np.array([[float(x) for x in l.split()] for l in lines])


def _load_col(path: Path, col: int) -> np.ndarray:
    """Load a single column (0-indexed) from a .dat file."""
    return _load_dat(path)[:, col]


# ── kLa helpers ──────────────────────────────────────────────────────────────

def _compute_c_star(run_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    """Compute the dimensionless dissolved-oxygen saturation C*(t).

    C* ranges from 0 (no dissolved oxygen) to 1 (fully saturated).  It is
    calculated as::

        C*(t) = oxy_liq_sum(t) / mean(f_liq_sum)

    where oxy_liq_sum is the total dissolved oxygen in the liquid (col 2 of
    tr_oxy.dat) and f_liq_sum is the liquid volume (col 2 of
    vol_frac_interf.dat, averaged over the whole run to remove interface noise).

    Returns
    -------
    t : np.ndarray
        Non-dimensional simulation time at each output step.
    c_star : np.ndarray
        C*(t) at each output step.
    """
    tr_path = run_dir / "tr_oxy.dat"
    vf_path = run_dir / "vol_frac_interf.dat"
    t           = _load_col(tr_path, 1)
    oxy_liq_sum = _load_col(tr_path, 2)
    f_liq_sum   = _load_col(vf_path, 2)
    f_mean = f_liq_sum.mean()
    if f_mean <= 0:
        raise ValueError("f_liq_sum mean is zero — VOF field may be empty")
    return t, oxy_liq_sum / f_mean


def _kla_5pt_at_threshold(t: np.ndarray, c_star: np.ndarray,
                           threshold: float) -> float:
    """Estimate kLa using a 5-point log-linear fit centred on a C* crossing.

    From the first-order model dC*/dt = kLa(1−C*), taking the log gives::

        ln(1 − C*(t)) = −kLa × t + const

    We fit this line over the 5 output steps closest to the first time C*
    crosses ``threshold``.  The slope gives −kLa.

    Returns NaN if C* never reaches the threshold or fewer than 5 points exist.
    """
    if len(t) < MIN_WINDOW:
        return math.nan
    idx = np.argmax(c_star >= threshold)
    if c_star[idx] < threshold:
        return math.nan
    half = MIN_WINDOW // 2
    lo = max(0, idx - half)
    hi = min(len(t), lo + MIN_WINDOW)
    lo = max(0, hi - MIN_WINDOW)
    t_win  = t[lo:hi]
    cs_win = np.clip(c_star[lo:hi], 0.0, 1.0 - 1e-10)
    y = np.log(1.0 - cs_win)
    slope, _ = np.polyfit(t_win, y, 1)
    return float(-slope)


def _kla_inst_at_threshold(t: np.ndarray, c_star: np.ndarray,
                            threshold: float) -> float:
    """Estimate kLa instantaneously at the moment C* crosses a threshold.

    Computes kLa = (dC*/dt) / (1 − C*) at the crossing point using a central
    finite difference for dC*/dt.  Falls back to a one-sided difference at
    array boundaries.

    Returns NaN if C* never reaches the threshold or fewer than 3 points exist.
    """
    if len(t) < 3:
        return math.nan
    idx = int(np.argmax(c_star >= threshold))
    if c_star[idx] < threshold:
        return math.nan
    if idx == 0:
        dc_dt = (float(c_star[1]) - float(c_star[0])) / (float(t[1]) - float(t[0]))
    elif idx == len(t) - 1:
        dc_dt = (float(c_star[-1]) - float(c_star[-2])) / (float(t[-1]) - float(t[-2]))
    else:
        dc_dt = (float(c_star[idx + 1]) - float(c_star[idx - 1])) / \
                (float(t[idx + 1])      - float(t[idx - 1]))
    denom = 1.0 - float(c_star[idx])
    if denom <= 1e-10:
        return math.nan
    return float(dc_dt / denom)


# ── mixing time ───────────────────────────────────────────────────────────────

def _compute_mixing_metrics(run_dir: Path, params: dict) -> dict:
    """Compute dimensional mixing times in seconds at χ = 0.50, 0.75, and 0.95.

    The passive tracer (rhodamine-like dye) is injected into the **top half**
    of the liquid at t_mix (after n_mix_cycles rocking cycles).  Before
    injection all tracer values are zero; at injection the top half is set to 1
    and the bottom half stays 0.

    The degree of mixing χ(t) is defined as::

        χ(t) = 1 − σ²(t) / σ²_max

    where σ²(t) is the spatial variance of the tracer concentration in the
    liquid::

        σ²(t) = ⟨c²⟩_liquid − ⟨c⟩²_liquid

    and σ²_max = σ²(t_inject) is the maximum variance (at the moment of
    injection, before any mixing occurs).  χ = 0 means completely segregated;
    χ = 1 means perfectly uniform.

    The tracer data come from column 8 (c2_liq_sum) and column 9
    (c2_liq_sum2) of tr_oxy.dat, which correspond to the VERTICAL_MIXUP
    configuration compiled into BioReactor.c.

    Parameters
    ----------
    run_dir : Path
        Run directory containing tr_oxy.dat and vol_frac_interf.dat.
    params : dict
        Simulation parameters (used to convert non-dimensional time to seconds
        via T_bio).

    Returns
    -------
    dict with keys ``dtmix_0.50``, ``dtmix_0.75``, ``dtmix_0.95`` in seconds.
    NaN is returned for a threshold when it is never reached within the
    simulation duration.
    """
    nan_result = {"dtmix_0.50": math.nan, "dtmix_0.75": math.nan,
                  "dtmix_0.95": math.nan}
    tr_path = run_dir / "tr_oxy.dat"
    vf_path = run_dir / "vol_frac_interf.dat"
    if not tr_path.exists() or not vf_path.exists():
        return nan_result

    arr_tr = _load_dat(tr_path)
    arr_vf = _load_dat(vf_path)
    if arr_tr.shape[0] < 5 or arr_tr.shape[1] < 10:
        return nan_result

    t      = arr_tr[:, 1]
    c_sum  = arr_tr[:, 8]   # c2_liq_sum  — VERTICAL_MIXUP tracer (top-half init)
    c_sum2 = arr_tr[:, 9]   # c2_liq_sum2
    f_mean = float(arr_vf[:, 2].mean())
    if f_mean <= 0:
        return nan_result

    c_mean  = c_sum  / f_mean
    c_mean2 = c_sum2 / f_mean
    sigma2  = c_mean2 - c_mean ** 2

    # Detect injection: first timestep where tracer becomes non-zero
    nonzero = np.where(c_sum > 1e-10 * f_mean)[0]
    if len(nonzero) == 0:
        return nan_result
    t0_idx     = int(nonzero[0])
    sigma2_max = float(sigma2[t0_idx])
    if sigma2_max <= 0:
        return nan_result

    chi = np.clip(1.0 - sigma2 / sigma2_max, 0.0, 1.0)
    chi[:t0_idx] = 0.0   # before injection χ is undefined → set to 0

    T_bio, _    = _t_scales(params)
    t_inject_nd = float(t[t0_idx])

    result = {}
    for threshold in (0.50, 0.75, 0.95):
        sub = chi[t0_idx:]
        idx = int(np.argmax(sub >= threshold))
        if sub[idx] < threshold:
            result[f"dtmix_{threshold:.2f}"] = math.nan
        else:
            t_cross_nd = float(t[t0_idx + idx]) - t_inject_nd
            result[f"dtmix_{threshold:.2f}"] = t_cross_nd * T_bio
    return result


# ── vorticity (steady streaming) ─────────────────────────────────────────────

def _compute_vor_mean(run_dir: Path, params: dict) -> float:
    """Compute the period-averaged mean absolute vorticity in 1/s.

    When the bag rocks sinusoidally, the *instantaneous* flow field oscillates
    at the rocking frequency.  However, the *time-averaged* flow (averaged over
    one full rocking period) reveals a slow secondary circulation called
    "steady streaming": two large counter-rotating vortices that drive net
    transport of tracers and oxygen over many cycles.

    This function computes::

        vor_mean = ⟨|Ω_liquid|⟩_time  [1/s]

    where Ω_liquid is the spatially averaged absolute vorticity over the liquid
    phase (column 2 of normf.dat, non-dimensional), and the time average is
    taken over the quasi-steady post-ramp portion of the run (t > 3 rocking
    periods after the ramp).

    The conversion from non-dimensional to dimensional uses::

        vor_mean [1/s] = Omega_nd [-] / T_bio [s]

    A larger vor_mean indicates stronger steady streaming, which correlates
    directly with faster mixing and higher kLa in the laminar regime.

    Parameters
    ----------
    run_dir : Path
        Run directory containing normf.dat.
    params : dict
        Simulation parameters for dimensional conversion.

    Returns
    -------
    float : vor_mean in 1/s, or NaN if normf.dat is absent or too short.
    """
    normf_path = run_dir / "normf.dat"
    if not normf_path.exists():
        return math.nan
    arr = _load_dat(normf_path)
    if arr.shape[0] < 10 or arr.shape[1] < 3:
        return math.nan

    t         = arr[:, 1]
    omega_avg = np.abs(arr[:, 2])   # Omega_liq_avg (non-dimensional)

    T_bio, T_per_nd = _t_scales(params)
    t_ramp = 3.0 * T_per_nd         # skip the soft-start ramp (3 rocking cycles)
    mask   = t > t_ramp
    if mask.sum() < 5:
        return math.nan

    vor_nd = float(np.mean(omega_avg[mask]))
    return vor_nd / T_bio


# ── main entry point ──────────────────────────────────────────────────────────

def main(run_dir: str, params: dict | None = None) -> dict:
    """Compute all ten KPIs and write them to results.json.

    This is the top-level function called by the SLURM script at the end of
    each simulation job.  It reads the raw output files, computes all metrics,
    and writes ``results.json`` to the run directory.

    Parameters
    ----------
    run_dir : str or Path
        Path to a completed simulation run directory (e.g. ``runs/abc12345/``).
        Must contain ``tr_oxy.dat``, ``vol_frac_interf.dat``, and ``normf.dat``.
    params : dict, optional
        Simulation parameters dict.  If None, loaded automatically from
        ``run_dir/params.json``.  Needed only for dimensional conversions;
        if absent, dtmix and vor_mean will be NaN.

    Returns
    -------
    dict
        All ten KPIs.  Any metric that cannot be computed (file missing,
        simulation too short, etc.) is set to NaN rather than raising.

        Keys and units:

        ============== =============================== ============
        Key            Description                     Unit
        ============== =============================== ============
        kLa_10         O2 transfer rate at C*=10%      1/t_nd
        kLa_25         O2 transfer rate at C*=25%      1/t_nd
        kLa_50         O2 transfer rate at C*=50%      1/t_nd
        kLa_inst_10    Instantaneous kLa at C*=10%     1/t_nd
        kLa_inst_25    Instantaneous kLa at C*=25%     1/t_nd
        kLa_inst_50    Instantaneous kLa at C*=50%     1/t_nd
        dtmix_0.50     Time to 50% mixing              seconds
        dtmix_0.75     Time to 75% mixing              seconds
        dtmix_0.95     Time to 95% mixing              seconds
        vor_mean       Mean absolute vorticity         1/s
        ============== =============================== ============

    Side effects
    ------------
    Writes ``results.json`` in run_dir (overwrites if exists).
    """
    path = Path(run_dir)
    if params is None:
        params = _load_params(path)

    nan_base = {
        "kLa_10": math.nan, "kLa_25": math.nan, "kLa_50": math.nan,
        "kLa_inst_10": math.nan, "kLa_inst_25": math.nan, "kLa_inst_50": math.nan,
        "dtmix_0.50": math.nan, "dtmix_0.75": math.nan, "dtmix_0.95": math.nan,
        "vor_mean": math.nan,
    }

    try:
        t, c_star = _compute_c_star(path)
    except (FileNotFoundError, ValueError):
        (path / "results.json").write_text(json.dumps(nan_base, indent=2))
        return nan_base

    results = {
        "kLa_10":      _kla_5pt_at_threshold(t, c_star, 0.10),
        "kLa_25":      _kla_5pt_at_threshold(t, c_star, 0.25),
        "kLa_50":      _kla_5pt_at_threshold(t, c_star, 0.50),
        "kLa_inst_10": _kla_inst_at_threshold(t, c_star, 0.10),
        "kLa_inst_25": _kla_inst_at_threshold(t, c_star, 0.25),
        "kLa_inst_50": _kla_inst_at_threshold(t, c_star, 0.50),
    }
    results.update(_compute_mixing_metrics(path, params))
    results["vor_mean"] = _compute_vor_mean(path, params)

    (path / "results.json").write_text(json.dumps(results, indent=2))
    return results


def validate_params(params: dict) -> None:
    """Check that a params dict is within the declared parameter-space bounds.

    Raises ValueError with a descriptive message if any parameter is out of
    range or structurally malformed (e.g. wrong vector length, phi_angular[0]
    ≠ 0).  Bounds are read from ``config/param_space.yaml``.

    Parameters
    ----------
    params : dict
        Contents of a params.json file.

    Raises
    ------
    ValueError
        If any parameter fails its bound or structural check.
    """
    import yaml
    _PARAM_SPACE = Path(__file__).parents[1] / "config" / "param_space.yaml"
    spec   = yaml.safe_load(_PARAM_SPACE.read_text())
    n_max  = spec["N_max"]
    pspace = spec["parameters"]

    def _check_scalar(name: str, val: float) -> None:
        lo, hi = pspace[name]["bounds"]
        if not (lo <= val <= hi):
            raise ValueError(f"{name}={val!r} out of bounds [{lo}, {hi}]")

    def _check_vector(name: str, vec: list, key: str | None = None,
                      n_active: int | None = None) -> None:
        if len(vec) != n_max:
            raise ValueError(f"{name} must have length {n_max}, got {len(vec)}")
        bounds = pspace[key or name]["bounds"]
        active = n_active if n_active is not None else n_max
        for i, v in enumerate(vec):
            if i < active:
                lo, hi = bounds
                if not (lo <= v <= hi):
                    raise ValueError(f"{name}[{i}]={v!r} out of bounds [{lo}, {hi}]")
            else:
                if v != 0.0:
                    raise ValueError(f"{name}[{i}] must be 0.0 (padding), got {v!r}")

    _check_scalar("omega_b", params["omega_b"])
    n_harm = params.get("n_harmonics", 1)
    lo, hi = pspace["n_harmonics"]["bounds"]
    if not (lo <= n_harm <= hi):
        raise ValueError(f"n_harmonics={n_harm!r} out of bounds [{lo}, {hi}]")
    _check_vector("theta_max", params["theta_max"], n_active=n_harm)

    phi_ang = params["phi_angular"]
    if len(phi_ang) != n_max:
        raise ValueError(f"phi_angular must have length {n_max}, got {len(phi_ang)}")
    if phi_ang[0] != 0.0:
        raise ValueError(f"phi_angular[0] must be 0.0, got {phi_ang[0]!r}")
    for i, v in enumerate(phi_ang[1:n_harm], start=1):
        lo, hi = pspace["phi_angular"]["bounds"]
        if not (lo <= v <= hi):
            raise ValueError(f"phi_angular[{i}]={v!r} out of bounds [{lo}, {hi}]")
    for i, v in enumerate(phi_ang[n_harm:], start=n_harm):
        if v != 0.0:
            raise ValueError(f"phi_angular[{i}] must be 0.0, got {v!r}")

    _check_scalar("omega_h", params["omega_h"])
    _check_vector("amplitude_h", params["amplitude_h"], n_active=n_harm)
    _check_vector("phi_horizontal", params["phi_horizontal"], n_active=n_harm)

    geom = params["geometry"]
    for sub in ("a", "b", "n"):
        key = f"geometry.{sub}"
        lo, hi = pspace[key]["bounds"]
        v = geom[sub]
        if not (lo <= v <= hi):
            raise ValueError(f"geometry.{sub}={v!r} out of bounds [{lo}, {hi}]")

    _check_scalar("fill_level", params["fill_level"])


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: postprocess.py <run_dir>", file=sys.stderr)
        sys.exit(1)
    res = main(sys.argv[1])
    print(json.dumps(res, indent=2))
