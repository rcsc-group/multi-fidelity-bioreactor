"""Convergence diagnostics plot for a parameter sweep.

Three panels, shared x-axis = time relative to oxygen injection (t - t_inject):

  Panel 1 — Cycle-averaged u_rms(t)   [normf.dat]
             Rolling mean over one rocking period T = 2π/omega_b.
             Removes per-cycle oscillations; plateau = QSS reached (global).

  Panel 2 — Rolling kLa(t)            [tr_oxy.dat + vol_frac_interf.dat]
             Instantaneous kLa from slope of ln(1 - C*(t)) over a sliding
             window of W rocking periods.  Plateau = KPI has converged (local).

  Panel 3 — Normalised C*(t) = C_liq(t) / C_sat   [tr_oxy.dat]
             C_sat ≈ 1.0 (confirmed from longest runs).  Reference lines at
             C*=0.10 and C*=0.25 (the kLa measurement windows).

Time axis: t - t_inject so all curves align at the injection event.

Line encoding (auto-detected from the actual sweep, or overridden with flags):
  Color     — swept param with the most unique values
  Linewidth — swept param with the second-most unique values

Usage
-----
    # Load by experiment (preferred — exact run set from one submission):
    python scripts/plot_convergence.py --experiment experiments/sweep_fb_theta_l7

    # Load by fidelity (all matching runs in runs/):
    python scripts/plot_convergence.py --fidelity 7

    # Override encoding:
    python scripts/plot_convergence.py --experiment ... --color rpm --lw theta_max_0

Outputs: experiments/figures/convergence_<stem>.pdf  (experiment mode)
      or experiments/figures/convergence_sweep_f<N>.pdf (fidelity mode)
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

# normf.dat columns
_N_T      = 1
_N_UX_RMS = 7
_N_UY_RMS = 11

# tr_oxy.dat columns
_T_T       = 1
_T_OXY_SUM = 2

# vol_frac_interf.dat columns
_V_T     = 1
_V_F_LIQ = 2

_C_SAT   = 1.0      # confirmed from longest runs (oxy tracer in [0,1] in gas)
_LW_MIN, _LW_MAX = 0.6, 2.8
_KLA_WINDOW_PERIODS = 2   # rolling window for kLa estimate, in rocking cycles


# ── Raw file loading ──────────────────────────────────────────────────────────

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
    def _get(d, key, default=float("nan")):
        return float(d.get(key, default)) if d.get(key) is not None else default

    geom = params.get("geometry") or {}
    th   = params.get("theta_max")   or []
    ph   = params.get("phi_angular") or []
    amp  = params.get("amplitude_h") or []
    phh  = params.get("phi_horizontal") or []

    omega_b = _get(params, "omega_b")
    return {
        "rpm":             round(omega_b * 60.0 / (2 * math.pi), 2),
        "omega_b":         omega_b,
        "fill_level":      _get(params, "fill_level"),
        "n_harmonics":     _get(params, "n_harmonics", 1),
        "omega_h":         _get(params, "omega_h", 0.0),
        "geometry_a":      float(geom.get("a", float("nan"))),
        "geometry_b":      float(geom.get("b", float("nan"))),
        "geometry_n":      float(geom.get("n", float("nan"))),
        "theta_max_0":     float(th[0])  if len(th) > 0 else float("nan"),
        "theta_max_1":     float(th[1])  if len(th) > 1 else 0.0,
        "theta_max_2":     float(th[2])  if len(th) > 2 else 0.0,
        "phi_angular_0":   float(ph[0])  if len(ph) > 0 else 0.0,
        "phi_angular_1":   float(ph[1])  if len(ph) > 1 else 0.0,
        "amplitude_h_0":   float(amp[0]) if len(amp) > 0 else 0.0,
        "phi_horizontal_0":float(phh[0]) if len(phh) > 0 else 0.0,
    }


def _detect_sweep_dims(records: list[dict]) -> list[str]:
    all_keys = list(records[0]["flat"].keys())
    n_unique = {}
    for k in all_keys:
        vals = {r["flat"][k] for r in records if not math.isnan(r["flat"][k])}
        if len(vals) > 1:
            n_unique[k] = len(vals)
    if "rpm" in n_unique and "omega_b" in n_unique:
        del n_unique["omega_b"]
    return sorted(n_unique, key=lambda k: -n_unique[k])


def _param_label(name: str) -> str:
    labels = {
        "rpm":           "RPM",
        "fill_level":    "Fill level",
        "theta_max_0":   r"$\theta_{max}$ (deg)",
        "omega_h":       r"$\omega_h$ (nd)",
        "amplitude_h_0": r"$A_{h,0}$",
        "geometry_a":    "Bag length $a$ (m)",
        "geometry_b":    "Bag height $b$ (m)",
    }
    return labels.get(name, name)


# ── Run loading ───────────────────────────────────────────────────────────────

def _load_runs(
    runs_root: Path,
    fidelity: int | None = None,
    run_ids: set[str] | None = None,
) -> list[dict]:
    """Load runs, deduplicating by param fingerprint (keep most recently modified).

    If run_ids is given, only consider those run directories (experiment mode).
    If fidelity is given as well, additionally filter by fidelity level.
    """
    candidates = (
        [runs_root / rid for rid in run_ids if (runs_root / rid).is_dir()]
        if run_ids is not None
        else [d for d in runs_root.iterdir() if d.is_dir()]
    )

    best: dict[tuple, tuple[float, dict]] = {}
    for run_dir in candidates:
        params_path  = run_dir / "params.json"
        results_path = run_dir / "results.json"
        if not params_path.exists() or not results_path.exists():
            continue
        try:
            params  = json.loads(params_path.read_text())
            results = json.loads(results_path.read_text())
        except Exception:
            continue
        if fidelity is not None and params.get("fidelity") != fidelity:
            continue

        flat  = _flat_params(params)
        key   = tuple(round(v, 4) for v in flat.values())
        mtime = results_path.stat().st_mtime
        rec   = {"run_dir": run_dir, "params": params, "results": results, "flat": flat}
        if key not in best or mtime > best[key][0]:
            best[key] = (mtime, rec)

    return [v for _, v in best.values()]


# ── Time-series computation ───────────────────────────────────────────────────

def _t_rock(params: dict) -> float:
    """Non-dimensional rocking period T = 2π / omega_b."""
    return 2 * math.pi / params.get("omega_b", 3.14)


def _t_inject(run_dir: Path) -> float | None:
    """Return absolute injection time from tr_oxy.dat, or None if not found."""
    tr = _load_dat(run_dir / "tr_oxy.dat")
    vf = _load_dat(run_dir / "vol_frac_interf.dat")
    if tr is None or vf is None or tr.shape[1] < 3 or vf.shape[1] < 3:
        return None
    f_mean = float(vf[:, _V_F_LIQ].mean())
    if f_mean <= 0:
        return None
    c = tr[:, _T_OXY_SUM] / f_mean
    inj = np.where(c > 1e-6)[0]
    return float(tr[inj[0], _T_T]) if len(inj) else None


def _cycle_avg_urms(
    run_dir: Path, T_rock: float, t_inj: float
) -> tuple[np.ndarray, np.ndarray] | None:
    """Rolling mean of sqrt(ux_rms² + uy_rms²) over one rocking period.

    Returns (t_rel, urms_smooth) where t_rel = t - t_inj, clipped to t_rel >= 0.
    Aligns with the C* and kLa panels so all three share the same time origin.
    """
    arr = _load_dat(run_dir / "normf.dat")
    if arr is None or arr.shape[1] < 12:
        return None
    t    = arr[:, _N_T]
    urms = np.sqrt(arr[:, _N_UX_RMS] ** 2 + arr[:, _N_UY_RMS] ** 2)

    dt_uniform = float(np.median(np.diff(t)))
    if dt_uniform <= 0:
        return None
    half = max(1, int(round(0.5 * T_rock / dt_uniform)))
    kernel = np.ones(2 * half + 1) / (2 * half + 1)
    smooth = np.convolve(urms, kernel, mode="same")

    # Trim the half-kernel edge region where zero-padding corrupts the average
    smooth = smooth[half:-half] if half < len(smooth) // 2 else smooth
    t_trim = t[half:-half] if half < len(t) // 2 else t

    t_rel = t_trim - t_inj
    mask  = t_rel >= 0
    return t_rel[mask], smooth[mask]


_C_STAR_LO = 0.05   # below this C* the kLa signal is too noisy (early)
_C_STAR_HI = 0.80   # above this C* → 1 makes ln(1-C*) blow up (late)


def _rolling_kla(
    run_dir: Path, T_rock: float
) -> tuple[np.ndarray, np.ndarray] | None:
    """Instantaneous kLa from slope of ln(1 - C*(t)) over a rolling window.

    Window = _KLA_WINDOW_PERIODS * T_rock.
    Returns (t_rel, kla) where t_rel is time relative to injection.
    Only covers t > t_inject; returns None if fewer than 2 window-widths of
    post-injection data are available.
    """
    tr = _load_dat(run_dir / "tr_oxy.dat")
    vf = _load_dat(run_dir / "vol_frac_interf.dat")
    if tr is None or vf is None or tr.shape[1] < 3 or vf.shape[1] < 3:
        return None

    f_mean = float(vf[:, _V_F_LIQ].mean())
    if f_mean <= 0:
        return None

    t   = tr[:, _T_T]
    c   = tr[:, _T_OXY_SUM] / f_mean   # C_liq (≡ C* since C_sat=1)

    inj = np.where(c > 1e-6)[0]
    if len(inj) == 0:
        return None
    i0 = inj[0]
    t_inj = t[i0]

    t_post = t[i0:] - t_inj
    c_post = np.clip(c[i0:] / _C_SAT, 0.0, 1.0 - 1e-6)
    lnterm = -np.log(1.0 - c_post)     # = kLa * (t - t_inj)

    dt_uniform = float(np.median(np.diff(t_post))) if len(t_post) > 1 else 1.0
    if dt_uniform <= 0:
        return None
    W = max(2, int(round(_KLA_WINDOW_PERIODS * T_rock / dt_uniform)))
    if len(t_post) < 2 * W:
        return None

    # Rolling OLS slope of lnterm vs t_post
    kla = np.full(len(t_post), np.nan)
    for k in range(W, len(t_post) - W):
        sl = slice(k - W, k + W + 1)
        x  = t_post[sl] - t_post[k]
        y  = lnterm[sl] - lnterm[k]
        sx2 = float(np.dot(x, x))
        if sx2 > 0:
            kla[k] = float(np.dot(x, y)) / sx2

    # Mask unreliable C* range
    in_range = (c_post >= _C_STAR_LO) & (c_post <= _C_STAR_HI)
    valid = ~np.isnan(kla) & in_range
    if valid.sum() < 5:
        return None
    return t_post[valid], kla[valid]


def _c_star_series(
    run_dir: Path,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Normalised C*(t) = C_liq / C_sat, time relative to injection."""
    tr = _load_dat(run_dir / "tr_oxy.dat")
    vf = _load_dat(run_dir / "vol_frac_interf.dat")
    if tr is None or vf is None or tr.shape[1] < 3 or vf.shape[1] < 3:
        return None

    f_mean = float(vf[:, _V_F_LIQ].mean())
    if f_mean <= 0:
        return None

    t = tr[:, _T_T]
    c = tr[:, _T_OXY_SUM] / f_mean / _C_SAT

    inj = np.where(c > 1e-6)[0]
    t_inj = t[inj[0]] if len(inj) else 0.0
    return t - t_inj, np.clip(c, 0.0, 1.0)


# ── Visual encoding ───────────────────────────────────────────────────────────

def _linear_map(val: float, vals: list[float], lo: float, hi: float) -> float:
    mn, mx = min(vals), max(vals)
    if mx == mn:
        return (lo + hi) / 2
    return lo + (val - mn) / (mx - mn) * (hi - lo)


# ── Main plot ─────────────────────────────────────────────────────────────────

def plot(
    fidelity: int | None = None,
    experiment: Path | None = None,
    color_param: str | None = None,
    lw_param: str | None = None,
    out_path: Path | None = None,
) -> Path:
    if experiment is None and fidelity is None:
        raise ValueError("Provide --experiment or --fidelity")

    run_ids: set[str] | None = None
    if experiment is not None:
        meta_path = Path(experiment) / "_sweep_metadata.json"
        if not meta_path.exists():
            raise FileNotFoundError(f"No _sweep_metadata.json in {experiment}")
        meta    = json.loads(meta_path.read_text())
        run_ids = set(meta.get("run_ids", []))
        if not run_ids:
            raise RuntimeError(f"_sweep_metadata.json has no run_ids in {experiment}")

    records = _load_runs(_RUNS_ROOT, fidelity=fidelity, run_ids=run_ids)
    if not records:
        raise RuntimeError("No matching runs found")

    sweep_dims = _detect_sweep_dims(records)
    if color_param is None:
        color_param = sweep_dims[0] if sweep_dims else "rpm"
    if lw_param is None:
        lw_param = sweep_dims[1] if len(sweep_dims) >= 2 else None

    color_vals = sorted({r["flat"][color_param] for r in records
                         if not math.isnan(r["flat"].get(color_param, float("nan")))})
    lw_vals    = sorted({r["flat"][lw_param] for r in records
                         if lw_param and not math.isnan(r["flat"].get(lw_param, float("nan")))}) \
                 if lw_param else []

    cmap       = matplotlib.colormaps["plasma"]
    color_norm = mcolors.Normalize(vmin=min(color_vals), vmax=max(color_vals))

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True, constrained_layout=True)
    ax_urms, ax_kla, ax_cstar = axes

    fid_label = str(fidelity) if fidelity is not None else Path(experiment).name
    plotted = 0
    for rec in sorted(records, key=lambda r: r["flat"].get(color_param, 0)):
        run_dir = rec["run_dir"]
        flat    = rec["flat"]
        params  = rec["params"]

        cv  = flat.get(color_param, float("nan"))
        lwv = flat.get(lw_param,    float("nan")) if lw_param else float("nan")
        if math.isnan(cv):
            continue

        color = cmap(color_norm(cv))
        lw    = _linear_map(lwv, lw_vals, _LW_MIN, _LW_MAX) if lw_vals and not math.isnan(lwv) \
                else (_LW_MIN + _LW_MAX) / 2
        kw    = dict(color=color, linewidth=lw, alpha=0.7)

        T_rock  = _t_rock(params)
        t_inj   = _t_inject(run_dir)
        if t_inj is None:
            continue

        urms_data  = _cycle_avg_urms(run_dir, T_rock, t_inj)
        kla_data   = _rolling_kla(run_dir, T_rock)
        cstar_data = _c_star_series(run_dir)

        if urms_data  is not None: ax_urms.plot( urms_data[0],  urms_data[1],  **kw)
        if kla_data   is not None: ax_kla.plot(  kla_data[0],   kla_data[1],   **kw)
        if cstar_data is not None: ax_cstar.plot(cstar_data[0], cstar_data[1], **kw)
        plotted += 1

    color_lbl = _param_label(color_param)
    lw_lbl    = _param_label(lw_param) if lw_param else "fixed"

    # Set left margin to -5% of whatever xmax matplotlib chose from the data
    _, x_right = ax_urms.get_xlim()
    for ax in axes:
        ax.set_xlim(left=-0.05 * x_right)

    ax_urms.set_ylabel(r"$\langle u_{rms} \rangle_T$ [nd]", fontsize=10)
    ax_urms.set_title(
        f"Convergence diagnostics — {fid_label}  ({plotted} runs)"
        f"   color={color_lbl}   lw={lw_lbl}",
        fontsize=10,
    )
    ax_urms.set_yscale("log")
    ax_urms.grid(True, which="both", alpha=0.3)

    ax_kla.set_ylabel(r"Rolling $kLa$ (nd)", fontsize=10)
    ax_kla.set_ylim(bottom=0)
    ax_kla.axhline(0, color="gray", lw=0.6, ls="--", alpha=0.4)
    ax_kla.grid(True, alpha=0.3)

    ax_cstar.set_ylabel(r"$C^*(t) = C_{liq}/C_{sat}$", fontsize=10)
    ax_cstar.set_xlabel(r"$t - t_{inject}$ [nd]", fontsize=10)
    ax_cstar.set_ylim(-0.02, 1.05)
    ax_cstar.axhline(0.10, color="gray", lw=0.8, ls="--", alpha=0.5, label=r"$C^*=0.10$")
    ax_cstar.axhline(0.25, color="gray", lw=0.8, ls=":",  alpha=0.5, label=r"$C^*=0.25$")
    ax_cstar.grid(True, alpha=0.3)
    ax_cstar.legend(fontsize=8, loc="upper left")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=color_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02)
    cbar.set_label(color_lbl, fontsize=10)

    if lw_vals:
        from matplotlib.lines import Line2D
        lw_handles = [
            Line2D([0], [0], color="gray",
                   linewidth=_linear_map(v, lw_vals, _LW_MIN, _LW_MAX),
                   label=f"{lw_lbl}={v:.3g}")
            for v in lw_vals
        ]
        ax_urms.legend(handles=lw_handles, fontsize=7, loc="lower right",
                       title="Line width", title_fontsize=7)

    if out_path is None:
        _FIG_DIR.mkdir(parents=True, exist_ok=True)
        if experiment is not None:
            stem = Path(experiment).name
        else:
            stem = f"sweep_f{fidelity}"
        out_path = _FIG_DIR / f"convergence_{stem}.pdf"

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--experiment", type=Path, default=None,
                     help="Path to experiment dir (e.g. experiments/sweep_fb_theta_l7)")
    grp.add_argument("--fidelity", type=int, default=None,
                     help="Load all runs at this fidelity from runs/")
    parser.add_argument("--color",  dest="color_param", default=None)
    parser.add_argument("--lw",     dest="lw_param",    default=None)
    parser.add_argument("--out",    type=Path,           default=None)
    args = parser.parse_args()
    out = plot(
        fidelity=args.fidelity,
        experiment=args.experiment,
        color_param=args.color_param,
        lw_param=args.lw_param,
        out_path=args.out,
    )
    print(f"Saved: {out}")
