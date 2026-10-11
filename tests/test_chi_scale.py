"""Unit tests for scripts/chi_scale.py on synthetic snapshots."""
import math
import struct
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))
from scripts import chi_scale as cs_mod  # noqa: E402

NAMES = cs_mod.PLANE_NAMES


def _write(path, n, t, planes):
    """Write a snapshot in the solver's binary layout (BioReactor.c write_snapshot)."""
    np_ = len(NAMES)
    with open(path, "wb") as fh:
        fh.write(struct.pack("<iid", n, np_, t))
        for name in NAMES:
            fh.write(np.asarray(planes.get(name, np.zeros((n, n))),
                                dtype="<f4").tobytes())


def _fields(n, seed=0, fractional=False):
    rng = np.random.default_rng(seed)
    f = np.zeros((n, n))
    f[: n // 2 + 3, :] = 1.0                      # liquid in the lower part, free surface above
    if fractional:
        f[n // 2 + 3, :] = rng.uniform(0.2, 0.8, n)
    cs = np.ones((n, n))
    cs[:2, :] = 0.0                               # solid wall row
    cs[2, :] = 0.5                                # cut cell
    c = rng.uniform(0, 1, (n, n))
    return f, cs, c


def test_reader_roundtrip(tmp_path):
    n = 16
    f, cs, c = _fields(n)
    p = tmp_path / "snap_000000.bin"
    _write(p, n, 3.25, {"f": f, "cs": cs, "c2": c})
    s = cs_mod.read_snapshot(p)
    assert s["n"] == n and s["np"] == 10 and s["t"] == pytest.approx(3.25)
    np.testing.assert_allclose(s["c2"], c.astype("f4"))
    np.testing.assert_allclose(s["f"], f)
    assert s["c2"].shape == (n, n)
    # j (y) is the slow index: row 0 is the bottom
    c_marker = np.zeros((n, n)); c_marker[1, 5] = 1.0
    _write(p, n, 0.0, {"c2": c_marker})
    assert cs_mod.read_snapshot(p, planes=["c2"])["c2"][1, 5] == 1.0
    assert "f" not in cs_mod.read_snapshot(p, planes=["c2"])


def test_iter_snapshots_time_order(tmp_path):
    n = 8
    for k, t in enumerate([0.0, 1.0, 2.0]):
        _write(tmp_path / f"snap_{k:06d}.bin", n, t, {})
    ts = [s["t"] for s in cs_mod.iter_snapshots(tmp_path, planes=["f"])]
    assert ts == [0.0, 1.0, 2.0]


def test_uniform_field_is_fully_mixed():
    n = 32
    f, cs, _ = _fields(n)
    c_seg = np.where(np.arange(n)[:, None] < n // 4, 1.0, 0.0) * np.ones((n, n))
    c_uni = np.full((n, n), 0.5)
    for block in (1, 4, 8):
        for mode in ("vol", "kim"):
            s0 = cs_mod.sigma2_coarse(f, cs, c_seg, block, mode=mode)
            s1 = cs_mod.sigma2_coarse(f, cs, c_uni, block, mode=mode)
            assert s0 > 0
            assert 1.0 - s1 / s0 == pytest.approx(1.0, abs=1e-12)


def test_coarse_graining_conserves_weighted_mean():
    n = 32
    f, cs, c = _fields(n, fractional=True)
    for block in (1, 2, 4, 8):
        Vb, Cb = cs_mod.coarse_grain(f, cs, c, block)
        total = (cs * f * c).sum()
        assert (Vb * Cb).sum() == pytest.approx(total, rel=1e-12)
        assert Vb.sum() == pytest.approx((cs * f).sum(), rel=1e-12)
    # boxes with zero liquid are excluded (zero weight), not NaN
    Vb, Cb = cs_mod.coarse_grain(f, cs, c, 8)
    assert np.all(np.isfinite(Cb))
    assert (Vb == 0).any()


def test_block1_reproduces_native_sigma2():
    n = 32
    f, cs, c = _fields(n, fractional=False)       # sharp interface: f in {0,1}
    nat = cs_mod.sigma2_native(f, cs, c)
    assert cs_mod.sigma2_coarse(f, cs, c, 1, mode="vol") == pytest.approx(nat, rel=1e-12)
    # native uses the solver's (c f) quantity; with f in {0,1} it equals the c-based value
    f2, cs2, c2 = _fields(n, fractional=True)
    nat2 = cs_mod.sigma2_native(f2, cs2, c2)
    vol2 = cs_mod.sigma2_coarse(f2, cs2, c2, 1, mode="vol")
    assert nat2 != pytest.approx(vol2, rel=1e-6)  # documented interface-cell difference
    assert abs(nat2 - vol2) / nat2 < 0.1


def test_crossing_interpolates_linearly():
    t = np.array([0.0, 1.0, 2.0, 3.0])
    chi = np.array([0.0, 0.2, 0.6, 0.9])
    assert cs_mod.crossing_time(t, chi, 0.4) == pytest.approx(1.5)
    assert np.isnan(cs_mod.crossing_time(t, chi, 0.95))


def test_chi_series_end_to_end_matches_solver_pipeline(tmp_path, monkeypatch):
    """Synthetic exponential mixing: snapshots + the .dat files the solver would
    write; native dtmix from snapshots must reproduce postprocess's dtmix."""
    import json
    n, tau = 32, 4.0
    f, cs, _ = _fields(n)
    seg = np.where(np.arange(n)[:, None] >= 11, 1.0, 0.0) * np.ones((n, n))
    run = tmp_path / "runs" / "fake"
    (run / "fields").mkdir(parents=True)
    params = {"omega_b": 3.4, "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
              "theta_max": [7.0, 0, 0], "n_mix_cycles": 1}
    (run / "params.json").write_text(json.dumps(params))
    from scripts import postprocess as pp
    T_bio, _ = pp._t_scales(params)
    t0 = 1.0
    dv = cs / n ** 2
    f_liq = float((dv * f).sum())

    def field(t):
        return 0.5 + (seg - 0.5) * np.exp(-max(t - t0, 0) / tau)

    tr_rows, vf_rows = [], []
    for k, t in enumerate(np.arange(0.0, t0 + 14 * tau, 0.05)):
        c = field(t) if t >= t0 else np.zeros((n, n))
        a = c * f
        tr_rows.append(f"{k} {t} 0 0 0 0 0 0 {(dv*a).sum()} {(dv*a*a).sum()} 0 0")
        vf_rows.append(f"{k} {t} {f_liq} 0 0 0 0 0")
    hdr = "i t a b c d e f g h\n"
    (run / "tr_oxy.dat").write_text(hdr + "\n".join(tr_rows) + "\n")
    (run / "vol_frac_interf.dat").write_text(hdr + "\n".join(vf_rows) + "\n")
    for k, t in enumerate(np.arange(t0, t0 + 12 * tau, 0.4)):
        _write(run / "fields" / f"snap_{k:06d}.bin", n, t, {"f": f, "cs": cs, "c2": field(t)})
    monkeypatch.setattr(cs_mod, "ROOT", tmp_path)
    monkeypatch.setattr(cs_mod, "SCRATCH", tmp_path / "none")
    s = cs_mod.chi_series("fake", divisors=(4,))
    assert s["status"] == "ok"
    ref = pp._compute_mixing_metrics(run, params)
    for thr in ("0.50", "0.75", "0.95"):
        exact = -tau / 2 * math.log(1 - float(thr)) * T_bio
        assert s["dtmix"]["native"][thr] == pytest.approx(exact, rel=1.5e-2)  # linear-in-time interpolation error at 0.4 cadence
        assert s["dtmix"]["native"][thr] == pytest.approx(ref[f"dtmix_{thr}"], rel=0.02)
        # uniform decay of a pattern: coarse chi follows the same law up to coarse-graining
        assert np.isfinite(s["dtmix"]["vol_4"][thr])
    assert s["release_offset_s"] == pytest.approx(0.0, abs=1e-9)


def _fake(dt_native, dt_l):
    keys = ("0.50", "0.75", "0.95")
    return {"status": "ok", "solver": {"dtmix": {}},
            "dtmix": {"native": dict(zip(keys, dt_native)), "vol_32": dict(zip(keys, dt_l))}}


def test_test_l_amended_rule():
    # 32.5 rpm: L9/L8 native ratio e^0.3 at all chi; l-scale ratio ~1 => artefact
    s = {}
    for rpm in (32.5, 37.5):
        s[(8, rpm)] = _fake((10, 20, 60), (10, 20, 60))
        s[(9, rpm)] = _fake((10 * math.e**0.3, 20 * math.e**0.3, 60 * math.e**0.3),
                            (10.1, 20.1, 60.1))
    v = cs_mod.test_l_verdict(s)
    assert v["n_counted"] == 6 and v["verdict"] == "SCALE ARTEFACT"
    s[(9, 37.5)] = _fake((10.1, 20.1, 60.1), (10.1, 20.1, 60.1))   # d_grid small: not counted
    assert cs_mod.test_l_verdict(s)["n_counted"] == 3
    s[(9, 32.5)] = _fake((10.1, 20.1, 60.1), (10.1, 20.1, 60.1))
    assert cs_mod.test_l_verdict(s)["verdict"].startswith("inconclusive")
    s[(9, 32.5)] = {"status": "missing"}
    assert "missing" in cs_mod.test_l_verdict(s)["verdict"]


def test_test_p_reading():
    assert cs_mod.test_p_reading({0.1: .25, 1e-2: .2, 1e-3: .15, 1e-6: .2})["verdict"] == "NOISE"
    assert cs_mod.test_p_reading({0.1: .25, 1e-2: .02, 1e-3: .002, 1e-6: 0.0})["verdict"] == "SMOOTH"
    assert cs_mod.test_p_reading({0.1: .25, 1e-2: .2, 1e-3: .05, 1e-6: .01})["verdict"] == "mixed"


def test_gradient_scale_sinusoid():
    n, periods = 256, 8
    k = 2 * np.pi * periods / n               # per cell
    x = np.arange(n)
    c = np.tile(0.5 + 0.5 * np.sin(k * x), (n, 1))
    f = np.ones((n, n)); cs = np.ones((n, n))
    g = cs_mod.gradient_scale(f, cs, c)
    assert g["lambda_dx"] == pytest.approx(1 / k, rel=0.02)
    assert 0.5 < g["frac_mid"] < 1.0
    # gas / solid cells are excluded: halving the liquid must not change lambda
    f2 = f.copy(); f2[n // 2:, :] = 0.0
    assert cs_mod.gradient_scale(f2, cs, c)["lambda_dx"] == pytest.approx(1 / k, rel=0.02)


def test_test_l_lists_unreached_cases():
    s = {}
    for rpm in (32.5, 37.5):
        s[(8, rpm)] = _fake((10, 20, 60), (10, 20, 60))
        s[(9, rpm)] = _fake((10 * math.e**0.3, 20 * math.e**0.3, 60 * math.e**0.3), (10, 20, 60))
    s[(9, 37.5)] = _fake((10 * math.e**0.3, 20 * math.e**0.3, math.nan), (10, 20, math.nan))
    s[(9, 37.5)]["chi"] = {"native": np.array([0.0, 0.8, 0.93]), "vol_32": np.array([0.0, 0.7, 0.9])}
    v = cs_mod.test_l_verdict(s)
    un = [c for c in v["cases"] if c["status"] == "unreached"]
    assert len(un) == 1 and un[0]["rpm"] == 37.5 and un[0]["chi"] == "0.95"
    assert {(u["level"], u["estimator"]): u["max_chi_reached"] for u in un[0]["unreached"]} == {
        (9, "native"): pytest.approx(0.93), (9, "vol_32"): pytest.approx(0.9)}
    assert v["n_counted"] == 5            # unreached is listed but not counted
