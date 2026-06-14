"""Convergence time-series plot for a parameter sweep.

Produces a 3-panel figure (shared x-axis = non-dimensional time):

  Panel 1 — vel_rms(t) = sqrt(u'x_rms² + u'y_rms²) [normf.dat]
             Global flow convergence to quasi-steady state
             (feeds into the vel_rms_qss KPI, Appendix A of Kim et al.)

  Panel 2 — C*(t) dissolved-oxygen saturation [tr_oxy.dat]
             Local kLa convergence; slope at any crossing = kLa
             (feeds into kLa_25, kla_fit_rmse_25)

  Panel 3 — χ(t) mixing homogeneity [tr_oxy.dat]
             1 − σ²(t)/σ²_max; rises from 0 to 1 after tracer injection

Line encoding (auto-detected from the actual sweep):
  Color     — the swept parameter with the most unique values
  Linewidth — the swept parameter with the second-most unique values
              (if only one parameter varies, linewidth is fixed)

Sweep detection: all numeric params in params.json are extracted; any param
whose value differs across runs is considered swept.  The top-2 by number of
distinct values get color and linewidth respectively.

Usage
-----
    python scripts/plot_convergence.py [--fidelity 5] [--out path.pdf]
    python scripts/plot_convergence.py --fidelity 5 --color rpm --lw fill_level

If multiple runs share the same parameter combination at the requested fidelity,
the most-recently modified run is used (later segment of a chain wins).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

_PROJECT_ROOT = Path(__file__).parents[1]
_RUNS_ROOT    = _PROJECT_ROOT / "runs"
_FIG_DIR      = _PROJECT_ROOT / "experiments" / "figures"

# normf.dat: i t Omega_avg Omega_rms Omega_vol Omega_max ux_avg ux_rms ux_vol ux_max uy_avg uy_rms uy_vol uy_max
_N_T      = 1
_N_UX_RMS = 7
_N_UY_RMS = 11
# tr_oxy.dat: i t oxy_liq_sum oxy_liq_sum2 ... c2_liq_sum c2_liq_sum2 ...
_T_T       = 1
_T_OXY_SUM = 2
_T_C2_SUM  = 8
_T_C2_SUM2 = 9
# vol_frac_interf.dat: i t f_liq_sum f_liq_interf posY_max posY_min
_V_F_LIQ = 2

_LW_MIN, _LW_MAX = 0.6, 2.8


def _load_dat(path: Path) -> np.ndarray | None:
    try:
        lines = [l for l in path.read_text().splitlines()
                 if l.strip() and not l.strip().startswith("i")]
        if not lines:
            return None
        return np.array([[float(x) for x in l.split()] for l in lines])
    except Exception:
        return None


# ── Parameter extraction ──────────────────────────────────────────────────────

def _flat_params(params: dict) -> dict[str, float]:
    """Flatten params.json into named scalars for sweep detection."""
    def _get(d, key, default=float("nan")):
        return float(d.get(key, default)) if d.get(key) is not None else default

    geom  = params.get("geometry") or {}
    th    = params.get("theta_max")    or []
    ph    = params.get("phi_angular")  or []
    amp   = params.get("amplitude_h")  or []
    phh   = params.get("phi_horizontal") or []

    omega_b = _get(params, "omega_b")
    flat = {
        "rpm":          round(omega_b * 60.0 / (2 * math.pi), 2),
        "omega_b":      omega_b,
        "fill_level":   _get(params, "fill_level"),
        "n_harmonics":  _get(params, "n_harmonics", 1),
        "omega_h":      _get(params, "omega_h", 0.0),
        "geometry_a":   float(geom.get("a", float("nan"))),
        "geometry_b":   float(geom.get("b", float("nan"))),
        "geometry_n":   float(geom.get("n", float("nan"))),
        "theta_max_0":  float(th[0]) if len(th) > 0 else float("nan"),
        "theta_max_1":  float(th[1]) if len(th) > 1 else 0.0,
        "theta_max_2":  float(th[2]) if len(th) > 2 else 0.0,
        "phi_angular_0": float(ph[0]) if len(ph) > 0 else 0.0,
        "phi_angular_1": float(ph[1]) if len(ph) > 1 else 0.0,
        "phi_angular_2": float(ph[2]) if len(ph) > 2 else 0.0,
        "amplitude_h_0": float(amp[0]) if len(amp) > 0 else 0.0,
        "amplitude_h_1": float(amp[1]) if len(amp) > 1 else 0.0,
        "phi_horizontal_0": float(phh[0]) if len(phh) > 0 else 0.0,
        "phi_horizontal_1": float(phh[1]) if len(phh) > 1 else 0.0,
    }
    return flat


def _detect_sweep_dims(records: list[dict]) -> list[str]:
    """Return param names that vary across runs, sorted by n_unique desc.

    'rpm' and 'omega_b' encode the same thing — keep only 'rpm' if both vary.
    """
    all_keys = [k for k in records[0]["flat"].keys()]
    n_unique = {}
    for k in all_keys:
        vals = {r["flat"][k] for r in records if not math.isnan(r["flat"][k])}
        if len(vals) > 1:
            n_unique[k] = len(vals)

    # drop omega_b if rpm is already there (same physical quantity)
    if "rpm" in n_unique and "omega_b" in n_unique:
        del n_unique["omega_b"]

    return sorted(n_unique, key=lambda k: -n_unique[k])


def _param_label(name: str) -> str:
    labels = {
        "rpm":            "RPM",
        "fill_level":     "Fill level",
        "theta_max_0":    r"$\theta_{max,0}$ (deg)",
        "theta_max_1":    r"$\theta_{max,1}$ (deg)",
        "omega_h":        r"$\omega_h$ (nd)",
        "amplitude_h_0":  r"$A_{h,0}$",
        "geometry_a":     "Bag length $a$ (m)",
        "geometry_b":     "Bag height $b$ (m)",
        "geometry_n":     "Bag taper $n$",
        "n_harmonics":    "N harmonics",
        "phi_angular_1":  r"$\phi_{\angle,1}$ (rad)",
    }
    return labels.get(name, name)


# ── Run loading ───────────────────────────────────────────────────────────────

def _load_runs(runs_root: Path, fidelity: int) -> list[dict]:
    """Load all fidelity-N runs; deduplicate by full param fingerprint."""
    best: dict[tuple, tuple[float, dict]] = {}

    for results_path in runs_root.glob("*/results.json"):
        run_dir     = results_path.parent
        params_path = run_dir / "params.json"
        if not params_path.exists():
            continue
        try:
            params  = json.loads(params_path.read_text())
            results = json.loads(results_path.read_text())
        except Exception:
            continue
        if params.get("fidelity") != fidelity:
            continue

        flat  = _flat_params(params)
        # fingerprint = all numeric params rounded to 4 sig figs
        key   = tuple(round(v, 4) for v in flat.values())
        mtime = results_path.stat().st_mtime
        rec   = {"run_dir": run_dir, "params": params, "results": results, "flat": flat}

        if key not in best or mtime > best[key][0]:
            best[key] = (mtime, rec)

    return [v for _, v in best.values()]


# ── Time-series loaders ───────────────────────────────────────────────────────

def _vel_rms_series(run_dir: Path) -> tuple[np.ndarray, np.ndarray] | None:
    arr = _load_dat(run_dir / "normf.dat")
    if arr is None or arr.shape[1] < 12:
        return None
    return arr[:, _N_T], np.sqrt(arr[:, _N_UX_RMS] ** 2 + arr[:, _N_UY_RMS] ** 2)


def _c_star_series(run_dir: Path) -> tuple[np.ndarray, np.ndarray] | None:
    tr = _load_dat(run_dir / "tr_oxy.dat")
    vf = _load_dat(run_dir / "vol_frac_interf.dat")
    if tr is None or vf is None or tr.shape[1] < 3 or vf.shape[1] < 3:
        return None
    f_mean = float(vf[:, _V_F_LIQ].mean())
    if f_mean <= 0:
        return None
    return tr[:, _T_T], tr[:, _T_OXY_SUM] / f_mean


def _chi_series(run_dir: Path) -> tuple[np.ndarray, np.ndarray] | None:
    tr = _load_dat(run_dir / "tr_oxy.dat")
    vf = _load_dat(run_dir / "vol_frac_interf.dat")
    if tr is None or vf is None or tr.shape[1] < 10 or vf.shape[1] < 3:
        return None
    f_mean = float(vf[:, _V_F_LIQ].mean())
    if f_mean <= 0:
        return None
    t, c_sum, c_sum2 = tr[:, _T_T], tr[:, _T_C2_SUM], tr[:, _T_C2_SUM2]
    c_mean  = c_sum  / f_mean
    sigma2  = c_sum2 / f_mean - c_mean ** 2
    nonzero = np.where(c_sum > 1e-10 * f_mean)[0]
    if not len(nonzero):
        return None
    t0 = int(nonzero[0])
    s2max = float(sigma2[t0])
    if s2max <= 0:
        return None
    chi = np.clip(1.0 - sigma2 / s2max, 0.0, 1.0)
    chi[:t0] = 0.0
    return t, chi


# ── Visual encoding helpers ───────────────────────────────────────────────────

def _linear_map(val: float, vals: list[float], lo: float, hi: float) -> float:
    mn, mx = min(vals), max(vals)
    if mx == mn:
        return (lo + hi) / 2
    return lo + (val - mn) / (mx - mn) * (hi - lo)


# ── Main plot ─────────────────────────────────────────────────────────────────

def plot(fidelity: int = 5,
         color_param: str | None = None,
         lw_param: str | None = None,
         out_path: Path | None = None) -> Path:

    records = _load_runs(_RUNS_ROOT, fidelity)
    if not records:
        raise RuntimeError(f"No fidelity-{fidelity} runs found in {_RUNS_ROOT}")

    sweep_dims = _detect_sweep_dims(records)

    # Resolve color and lw params
    if color_param is None:
        color_param = sweep_dims[0] if sweep_dims else "rpm"
    if lw_param is None:
        lw_param = sweep_dims[1] if len(sweep_dims) >= 2 else None

    color_vals = sorted({r["flat"][color_param] for r in records
                         if not math.isnan(r["flat"].get(color_param, float("nan")))})
    lw_vals    = sorted({r["flat"][lw_param] for r in records
                         if lw_param and not math.isnan(r["flat"].get(lw_param, float("nan")))}) \
                 if lw_param else []

    cmap      = matplotlib.colormaps["plasma"]
    color_norm = mcolors.Normalize(vmin=min(color_vals), vmax=max(color_vals))

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True, constrained_layout=True)
    ax_vel, ax_cstar, ax_chi = axes

    plotted = 0
    for rec in sorted(records, key=lambda r: r["flat"].get(color_param, 0)):
        run_dir = rec["run_dir"]
        flat    = rec["flat"]

        cv  = flat.get(color_param, float("nan"))
        lwv = flat.get(lw_param, float("nan")) if lw_param else float("nan")

        if math.isnan(cv):
            continue

        color = cmap(color_norm(cv))
        lw    = _linear_map(lwv, lw_vals, _LW_MIN, _LW_MAX) if lw_vals and not math.isnan(lwv) \
                else (_LW_MIN + _LW_MAX) / 2

        kw = dict(color=color, linewidth=lw, alpha=0.75)

        vel_data   = _vel_rms_series(run_dir)
        cstar_data = _c_star_series(run_dir)
        chi_data   = _chi_series(run_dir)

        if vel_data   is not None: ax_vel.plot(vel_data[0],   vel_data[1],   **kw)
        if cstar_data is not None: ax_cstar.plot(cstar_data[0], cstar_data[1], **kw)
        if chi_data   is not None: ax_chi.plot(chi_data[0],   chi_data[1],   **kw)
        plotted += 1

    color_lbl = _param_label(color_param)
    lw_lbl    = _param_label(lw_param) if lw_param else "fixed"

    ax_vel.set_ylabel(r"$u'_{rms}(t)$ [nd]", fontsize=10)
    ax_vel.set_title(
        f"Convergence diagnostics — fidelity {fidelity}  "
        f"({plotted} runs)   color={color_lbl}   lw={lw_lbl}",
        fontsize=10,
    )
    ax_vel.set_yscale("log")
    ax_vel.grid(True, which="both", alpha=0.3)

    ax_cstar.set_ylabel(r"$C^*(t)$", fontsize=10)
    ax_cstar.set_ylim(-0.02, 1.05)
    ax_cstar.axhline(0.25, color="gray", lw=0.8, ls="--", alpha=0.5, label="C*=0.25")
    ax_cstar.grid(True, alpha=0.3)
    ax_cstar.legend(fontsize=8, loc="upper left")

    ax_chi.set_ylabel(r"$\chi(t)$", fontsize=10)
    ax_chi.set_xlabel("Non-dimensional time $t$ [-]", fontsize=10)
    ax_chi.set_ylim(-0.02, 1.05)
    ax_chi.axhline(0.95, color="gray", lw=0.8, ls="--", alpha=0.5, label=r"$\chi=0.95$")
    ax_chi.grid(True, alpha=0.3)
    ax_chi.legend(fontsize=8, loc="upper left")

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=color_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02)
    cbar.set_label(color_lbl, fontsize=10)

    # Linewidth legend
    if lw_vals:
        from matplotlib.lines import Line2D
        lw_handles = [
            Line2D([0], [0], color="gray",
                   linewidth=_linear_map(v, lw_vals, _LW_MIN, _LW_MAX),
                   label=f"{lw_lbl}={v:.3g}")
            for v in lw_vals
        ]
        ax_vel.legend(handles=lw_handles, fontsize=7, loc="upper right",
                      title="Line width", title_fontsize=7)

    if out_path is None:
        _FIG_DIR.mkdir(parents=True, exist_ok=True)
        out_path = _FIG_DIR / f"convergence_sweep_f{fidelity}.pdf"

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fidelity", type=int, default=5)
    parser.add_argument("--color",  dest="color_param", default=None,
                        help="Param name to encode as color (default: auto-detected)")
    parser.add_argument("--lw",     dest="lw_param",    default=None,
                        help="Param name to encode as linewidth (default: auto-detected)")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    out = plot(fidelity=args.fidelity,
               color_param=args.color_param,
               lw_param=args.lw_param,
               out_path=args.out)
    print(f"Saved: {out}")
