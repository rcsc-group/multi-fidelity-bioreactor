"""Campaign 2D (2026-11-04): finished campaign runs (runs/c2d_*) reach the fit through campaign_row.

The database builders read fixed run lists, so the fit adds the campaign folders itself. Synthetic run folder only.
A run without results.json is skipped (not yet finished). A kLa value comes from results.json; a hydro value is
the same qois() as hydro_dataset_v2.py, so it is not re-derived here.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))

import campaign as c  # noqa: E402


def _make_run(folder: Path, finished: bool) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    rpm = 20.0
    omega = rpm * 2 * np.pi / 60.0
    p = {"fidelity": 8, "omega_b": omega, "theta_max": [4.0, 0.0, 0.0], "n_mix_cycles": 3,
         "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "t_checkpoint": 0.0}
    (folder / "params.json").write_text(json.dumps(p))
    if finished:
        (folder / "results.json").write_text(json.dumps({"kLa_1T_10": 12.5}))
    return folder


def test_unfinished_run_is_skipped(tmp_path):
    d = _make_run(tmp_path / "c2d_l8_rpm20.00_th4.00", finished=False)
    assert c.campaign_row(d, "kLa_1T_10", "kla") is None


def test_kla_value_read_from_results(tmp_path):
    d = _make_run(tmp_path / "c2d_l8_rpm20.00_th4.00", finished=True)
    row = c.campaign_row(d, "kLa_1T_10", "kla")
    assert row["y"] == 12.5 and row["level"] == 8
    assert abs(row["rpm"] - 20.0) < 1e-9 and abs(row["theta"] - 4.0) < 1e-9
    assert row["run"] == d.name


def test_duplicate_run_rows_are_dropped_once():
    import pandas as pd
    d = pd.DataFrame({"level": [8, 8, 9], "rpm": [20.0, 20.0, 25.0], "theta": [7.0, 7.0, 4.0],
                      "run": ["kmix_l8_rpm20", "kmix_l8_rpm20", "x"], "y": [1.0, 1.0, 2.0]})
    out = c.dedup_runs(d)
    assert list(out.run) == ["kmix_l8_rpm20", "x"]


def test_database_fit_rows_are_one_per_run():
    # kla CSV lists six runs twice (sets 'design' and 't5_grid'); the fit must see each run once
    d = c.load_db("kLa_1T_10")
    assert d.run.is_unique
