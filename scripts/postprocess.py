"""Extract kLa, mixing time, and vorticity metrics from BioReactor output files.

Input files (in run_dir):
  tr_oxy.dat:          i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 ...  (12 cols)
  vol_frac_interf.dat: i t f_liq_sum f_liq_interf posY_max posY_min
  normf.dat:           i t Omega_liq_avg Omega_liq_rms ... (14 cols)
  params.json:         simulation parameters (used for dimensional conversion)

Metrics computed
----------------
kLa_10, kLa_25, kLa_50          5-point log-linear fit at C*=0.10/0.25/0.50
kLa_inst_10, kLa_inst_25, kLa_inst_50   instantaneous dC*/dt / (1-C*) at crossing
dtmix_0.50, dtmix_0.75, dtmix_0.95      dimensional mixing time (s) at χ=0.50/0.75/0.95
vor_mean                                 dimensional period-averaged |ξ| (1/s)

Notes
-----
- C* = oxy_liq_sum / f_liq_mean  (dimensionless dissolved O2, 0→1 at saturation)
- χ = 1 - σ²(t)/σ²_max  (degree of mixing, 0=segregated, 1=mixed)
- σ²(t) = c_liq_sum2/f_mean - (c_liq_sum/f_mean)²
- σ²_max is taken at the first non-zero tracer time step (injection event)
- Shear stress and EDR are NOT computed here — they require BioReactor.c diagnostics
  not present in the current output files.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


MIN_WINDOW = 5   # minimum number of time points for a valid kLa estimate


# ── time-scale helpers ────────────────────────────────────────────────────────

def _t_scales(params: dict) -> tuple[float, float]:
    """Return (T_bio [s], T_per_nd [-]) from params.

    Mirrors BioReactor.c non-dimensionalisation:
      T_per   = 2π / omega_b
      V_char  = L/4 * (H + 0.5*L*tan(theta_max[0]))
      U_bio   = V_char / (H/2) / T_per
      T_bio   = L / U_bio
      T_per_nd = T_per / T_bio
    """
    omega_b = params.get("omega_b", 3.93)
    L       = params.get("geometry", {}).get("a", 0.25)
    H       = params.get("geometry", {}).get("b", 0.071)
    th      = math.radians(params.get("theta_max", [7.0])[0])
    T_per   = 2 * math.pi / omega_b
    V       = L / 4 * (H + 0.5 * L * math.tan(th))
    U       = V / (H * 0.5) / T_per
    T_bio   = L / U
    T_per_nd = T_per / T_bio
    return T_bio, T_per_nd


def _load_params(run_dir: Path) -> dict:
    p = run_dir / "params.json"
    return json.loads(p.read_text()) if p.exists() else {}


# ── raw data loaders ──────────────────────────────────────────────────────────

def _load_dat(path: Path) -> np.ndarray:
    """Load space-delimited .dat file, skipping header line starting with 'i'."""
    lines = [
        l for l in path.read_text().splitlines()
        if l.strip() and not l.strip().startswith("i")
    ]
    return np.array([[float(x) for x in l.split()] for l in lines])


def _load_col(path: Path, col: int) -> np.ndarray:
    return _load_dat(path)[:, col]


# ── kLa helpers ──────────────────────────────────────────────────────────────

def _compute_c_star(run_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (t, C*) from tr_oxy.dat and vol_frac_interf.dat."""
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
    """kLa from 5-point log-linear fit of ln(1-C*) vs t centred on crossing."""
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
    """Instantaneous kLa = (dC*/dt) / (1-C*) at first crossing via central diff."""
    if len(t) < 3:
        return math.nan
    idx = int(np.argmax(c_star >= threshold))
    if c_star[idx] < threshold:
        return math.nan
    # Use central diff where possible; fall back to one-sided at boundaries
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
    """Compute dimensional mixing times (s) at χ = 0.50, 0.75, 0.95.

    χ = 1 - σ²(t)/σ²_max
    σ²(t) = c_liq_sum2/f_mean - (c_liq_sum/f_mean)²

    Returns dict with keys dtmix_0.50, dtmix_0.75, dtmix_0.95.
    NaN when the threshold is never reached or data are missing.
    """
    nan_result = {"dtmix_0.50": math.nan, "dtmix_0.75": math.nan,
                  "dtmix_0.95": math.nan}
    tr_path = run_dir / "tr_oxy.dat"
    vf_path = run_dir / "vol_frac_interf.dat"
    if not tr_path.exists() or not vf_path.exists():
        return nan_result

    arr_tr = _load_dat(tr_path)
    arr_vf = _load_dat(vf_path)
    if arr_tr.shape[0] < 5 or arr_tr.shape[1] < 6:
        return nan_result

    t      = arr_tr[:, 1]
    c_sum  = arr_tr[:, 8]   # c2_liq_sum  (VERTICAL_MIXUP tracer; cols 4-5 are HORIZONTAL_MIXL)
    c_sum2 = arr_tr[:, 9]   # c2_liq_sum2
    f_mean = float(arr_vf[:, 2].mean())
    if f_mean <= 0:
        return nan_result

    c_mean  = c_sum  / f_mean
    c_mean2 = c_sum2 / f_mean
    sigma2  = c_mean2 - c_mean ** 2

    # find first time step where tracer is injected (c_sum goes non-zero)
    nonzero = np.where(c_sum > 1e-10 * f_mean)[0]
    if len(nonzero) == 0:
        return nan_result
    t0_idx      = int(nonzero[0])
    sigma2_max  = float(sigma2[t0_idx])
    if sigma2_max <= 0:
        return nan_result

    chi = np.clip(1.0 - sigma2 / sigma2_max, 0.0, 1.0)
    chi[:t0_idx] = 0.0   # before injection chi is undefined, set to 0

    T_bio, _ = _t_scales(params)
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
    """Period-averaged mean absolute vorticity in dimensional units (1/s).

    Averages Omega_liq_avg (col 2 of normf.dat) over post-ramp quasi-steady
    data (t > 3 * T_per_nd), then converts to dimensional 1/s via T_bio.
    """
    normf_path = run_dir / "normf.dat"
    if not normf_path.exists():
        return math.nan
    arr = _load_dat(normf_path)
    if arr.shape[0] < 10 or arr.shape[1] < 3:
        return math.nan

    t         = arr[:, 1]
    omega_avg = np.abs(arr[:, 2])   # Omega_liq_avg

    T_bio, T_per_nd = _t_scales(params)
    t_ramp = 3.0 * T_per_nd
    mask = t > t_ramp
    if mask.sum() < 5:
        return math.nan

    vor_nd = float(np.mean(omega_avg[mask]))
    return vor_nd / T_bio           # dimensional 1/s


# ── main entry point ──────────────────────────────────────────────────────────

def main(run_dir: str, params: dict | None = None) -> dict:
    """Compute all postprocessing metrics and write results.json.

    Parameters
    ----------
    run_dir : str | Path
        Completed BioReactor run directory.
    params : dict, optional
        Simulation parameters.  If None, loaded from run_dir/params.json.
        Required for dimensional conversions (mixing time, vorticity).

    Returns
    -------
    dict with keys:
      kLa_10, kLa_25, kLa_50          — 5-point log-linear fit
      kLa_inst_10, kLa_inst_25, kLa_inst_50  — instantaneous fit
      dtmix_0.50, dtmix_0.75, dtmix_0.95     — dimensional mixing time (s)
      vor_mean                                — mean absolute vorticity (1/s)
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
        "kLa_10":       _kla_5pt_at_threshold(t, c_star, 0.10),
        "kLa_25":       _kla_5pt_at_threshold(t, c_star, 0.25),
        "kLa_50":       _kla_5pt_at_threshold(t, c_star, 0.50),
        "kLa_inst_10":  _kla_inst_at_threshold(t, c_star, 0.10),
        "kLa_inst_25":  _kla_inst_at_threshold(t, c_star, 0.25),
        "kLa_inst_50":  _kla_inst_at_threshold(t, c_star, 0.50),
    }
    results.update(_compute_mixing_metrics(path, params))
    results["vor_mean"] = _compute_vor_mean(path, params)

    (path / "results.json").write_text(json.dumps(results, indent=2))
    return results


def validate_params(params: dict) -> None:
    """Raise ValueError if params are out of bounds or structurally malformed."""
    import yaml
    _PARAM_SPACE = Path(__file__).parents[1] / "config" / "param_space.yaml"
    spec = yaml.safe_load(_PARAM_SPACE.read_text())
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
