"""Readers and metrics of the guard tests. Each metric is defined in the header of tools/bundle/bars.yaml."""
import math
import struct
from pathlib import Path

import numpy as np

from scripts import periods
from scripts.replay_gate import rel_diff_by_time
from tools.bundle import cases as bc

EPS = 1e-9


def table(rd, fname: str) -> tuple[list, np.ndarray]:
    """Header names and the numbers of one solver .dat file (first row is the header)."""
    lines = (Path(rd) / fname).read_text().splitlines()
    rows = [[float(x) for x in ln.split()] for ln in lines[1:] if ln.strip()]
    return lines[0].split(), np.array(rows)


def column(rd, fname: str, name: str) -> tuple[np.ndarray, np.ndarray]:
    head, d = table(rd, fname)
    return d[:, head.index("t")], d[:, head.index(name)]


def dump_time(path) -> float:
    return struct.unpack("<d", Path(path).read_bytes()[:8])[0]


def window(spec, rpm: float = bc.RPM_A, rd=None) -> tuple[float, float]:
    """bars.yaml window as absolute times: [a, b] in periods, or "last" / "lastN" = the 1 / N periods before
    the run's own final dump (a warm start begins at another absolute time than its cold reference)."""
    T = bc.period(rpm)
    if isinstance(spec, str) and spec.startswith("last"):
        n = int(spec[4:] or 1)
        t1 = dump_time(Path(rd) / "checkpoint.dump")
        return t1 - n * T, t1
    return spec[0] * T, spec[1] * T


def _in(t, w):
    return (t >= w[0] - EPS) & (t <= w[1] + EPS)


def _mean(rd, fname, name, w):
    t, v = column(rd, fname, name)
    m = _in(t, w)
    if m.sum() < 3:
        raise AssertionError(f"{Path(rd).name}/{fname}: {m.sum()} rows in window {w}")
    return t[m], v[m]


def metric(name: str, rd, w) -> float:
    """One run's value of a bars.yaml metric over the window w = (t0, t1)."""
    if name == "vel_rms":
        _, ux = _mean(rd, "normf.dat", "ux_liq_rms", w)
        _, uy = _mean(rd, "normf.dat", "uy_liq_rms", w)
        return float(np.mean(np.sqrt(ux ** 2 + uy ** 2)))
    if name in ("tau_rms", "ediss_rms"):
        _, v = _mean(rd, "shear_stress.dat", "tau_mean" if name == "tau_rms" else "ediss_mean", w)
        return float(np.sqrt(np.mean(v ** 2)))
    if name == "tau_95":
        return float(np.mean(_mean(rd, "shear_stress.dat", "tau_95", w)[1]))
    if name == "edr_mean":
        return float(np.mean(_mean(rd, "edr_pct.dat", "edr_95", w)[1]))
    if name in ("oxy_last", "c2_last"):
        return float(_mean(rd, "tr_oxy.dat", "oxy_liq_sum" if name == "oxy_last" else "c2_liq_sum", w)[1][-1])
    raise KeyError(name)


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def series_max(col: str, rd_ref, rd, w) -> tuple[float, int]:
    """Max relative difference of a tr_oxy.dat column at the times both runs logged in w; and that count."""
    head, a = table(rd_ref, "tr_oxy.dat")
    _, b = table(rd, "tr_oxy.dat")
    a, b = a[_in(a[:, 1], w)], b[_in(b[:, 1], w)]
    d = rel_diff_by_time(a, b, [head.index(col)])
    return d["max"][head.index(col)], d["n_common"]


def compare(bar: dict, rd_ref, rd, rpm: float = bc.RPM_A) -> dict:
    """{metric: (ref value, value, relative difference, tolerance)} for a `kind: metric` bar."""
    w_ref, w = window(bar["window"], rpm, rd_ref), window(bar["window"], rpm, rd)
    out = {}
    for m, tol in bar["tol"].items():
        if m.endswith("_max"):
            assert w_ref == w, "a series comparison needs one window"
            d, n = series_max(m[:-4], rd_ref, rd, w)
            assert n >= 3, f"{m}: only {n} common times in window {w}"
            out[m] = (math.nan, math.nan, d, tol)
        else:
            a, b = metric(m, rd_ref, w_ref), metric(m, rd, w)
            if bar.get("sd_periods"):
                tol = max(tol, 2 * period_sd(m, rd_ref, w_ref, bar["sd_periods"], rpm))
            out[m] = (a, b, rel(a, b), tol)
    return out


def period_sd(m: str, rd, w, n: int, rpm: float) -> float:
    """Relative period-to-period scatter of metric m in run rd: sample sd / |mean| of its value over the n
    periods that end at w[1] (the x8 bar's "2 sd": a single period of an unsteady flow scatters by this much)."""
    T = bc.period(rpm)
    v = np.array([metric(m, rd, (w[1] - (k + 1) * T, w[1] - k * T)) for k in range(n)])
    return float(np.std(v, ddof=1) / abs(np.mean(v)))


def failures(result: dict) -> list:
    """Metrics over their bar, as readable lines (a NaN difference is a failure, never a pass)."""
    return [f"{m}: ref {a:.6g}, case {b:.6g}, rel diff {d:.3g} > {tol}"
            for m, (a, b, d, tol) in result.items() if not d <= tol]


def stop_time(t_request: float, rpm: float = bc.RPM_A) -> float:
    """Where the solver must stop for a requested absolute end time (scripts/periods.py rule)."""
    return periods.next_period_boundary(t_request, bc.period(rpm))[1]
