"""Replay flow gate (plan 3.1, amended 2026-11-04): tau_95 and posY_max of a replay run must agree with the
recording DNS within REPLAY_FLOW_TOL (2 %) before its scalars are released. Synthetic run folders only.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))
sys.path.insert(0, str(ROOT / "scripts"))

import campaign as c  # noqa: E402
import replay_gate as rg  # noqa: E402

T = 1.0


def _write_run(folder: Path, tau_scale: float, posY_scale: float, n_cycles: int = 12, per_cycle: int = 20):
    folder.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0.0, n_cycles * T, n_cycles * per_cycle, endpoint=False)
    tau = 1.0 + 0.05 * np.sin(2 * np.pi * t / T)
    pos = 0.3 + 0.02 * np.cos(2 * np.pi * t / T)
    idx = np.arange(t.size)
    with open(folder / "shear_stress.dat", "w") as f:
        f.write("i t tau_95 tau_98\n")
        for i in range(t.size):
            f.write(f"{idx[i]} {t[i]:.8g} {tau[i] * tau_scale:.8g} 0\n")
    with open(folder / "vol_frac_interf.dat", "w") as f:
        f.write("i t f_liq_sum f_liq_interf posY_max posY_min posY_left posY_right\n")
        for i in range(t.size):
            f.write(f"{idx[i]} {t[i]:.8g} 0 0 {pos[i] * posY_scale:.8g} 0 0 0\n")
    return folder


def test_matching_replay_passes(tmp_path):
    rec = _write_run(tmp_path / "rec", 1.0, 1.0)
    rep = _write_run(tmp_path / "rep", 1.0, 1.0)
    ok, dev = rg.flow_gate(rec, rep, T, t_lo=11 * T)
    assert ok, dev
    assert abs(dev["tau_95"]) < 1e-9 and abs(dev["posY_max"]) < 1e-9


def test_tau_offset_within_2_percent_passes(tmp_path):
    rec = _write_run(tmp_path / "rec", 1.0, 1.0)
    rep = _write_run(tmp_path / "rep", 1.015, 1.0)
    ok, dev = rg.flow_gate(rec, rep, T, t_lo=11 * T)
    assert ok, dev
    assert abs(dev["tau_95"] - 1.5) < 0.05 and abs(dev["posY_max"]) < 1e-6


def test_posY_offset_above_2_percent_fails(tmp_path):
    rec = _write_run(tmp_path / "rec", 1.0, 1.0)
    rep = _write_run(tmp_path / "rep", 1.0, 0.97)
    ok, dev = rg.flow_gate(rec, rep, T, t_lo=11 * T)
    assert not ok
    assert abs(dev["posY_max"] + 3.0) < 0.05


def test_flow_check_ok_is_the_limit_used():
    assert c.flow_check_ok(0.019, -0.019)
    assert not c.flow_check_ok(0.021, 0.0)
    assert not c.flow_check_ok(0.0, -0.025)
