"""Suggest the next experiment via Expected Improvement on the MF surrogate.

Public API
----------
main(experiment_dir, param_space, model_path, kla_key, lf_fidelity, hf_fidelity,
     n_candidates) -> dict
    Returns a nested params dict (same structure as params.json) for the
    highest-EI candidate sampled from param_space.
"""
from __future__ import annotations

import pickle
import tempfile
from pathlib import Path

import numpy as np
import yaml
from scipy.stats import norm

import scripts.train_surrogate as _ts

_N_MAX = 3
_FEATURE_COLS = _ts._FEATURE_COLS   # 18-element list


def _load_param_space(path: str) -> dict:
    return yaml.safe_load(Path(path).read_text())


def _sample_candidates(spec: dict, n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw n random candidates in the 18-D feature space (uniform per bound)."""
    pspace = spec["parameters"]
    bounds = []
    for col in _FEATURE_COLS:
        if col.startswith("theta_max") or col.startswith("phi_angular") \
                or col.startswith("amplitude_h") or col.startswith("phi_horizontal"):
            # vector param: extract base name
            base = col.rsplit("_", 1)[0]
            idx  = int(col.rsplit("_", 1)[1])
            if idx == 0 and base == "phi_angular":
                bounds.append((0.0, 0.0))   # always fixed to 0
                continue
            bounds.append(tuple(pspace[base]["bounds"]))
        elif col.startswith("geometry_"):
            sub  = col[len("geometry_"):]   # a, b, or n
            key  = f"geometry.{sub}"
            bounds.append(tuple(pspace[key]["bounds"]))
        else:
            bounds.append(tuple(pspace[col]["bounds"]))

    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    X  = rng.uniform(lo, hi, size=(n, len(bounds)))
    if "phi_angular_0" in _FEATURE_COLS:
        X[:, _FEATURE_COLS.index("phi_angular_0")] = 0.0
    return X


def _expected_improvement(mean: np.ndarray, std: np.ndarray,
                           y_best: float, xi: float = 0.01) -> np.ndarray:
    """Standard EI acquisition (maximisation)."""
    improvement = mean - y_best - xi
    with np.errstate(divide="ignore"):
        Z   = np.where(std > 0, improvement / std, 0.0)
    ei  = np.where(std > 0,
                   improvement * norm.cdf(Z) + std * norm.pdf(Z),
                   0.0)
    return ei


def _flat_to_nested(flat: np.ndarray, spec: dict, hf_fidelity: int) -> dict:
    """Convert 18-element feature vector back to a nested params dict."""
    pspace = spec["parameters"]
    n_max  = spec["N_max"]
    cols   = _FEATURE_COLS
    v      = dict(zip(cols, flat))

    def _vec(base: str) -> list:
        vals = [float(v.get(f"{base}_{i}", 0.0)) for i in range(n_max)]
        return vals

    return {
        "fidelity":        hf_fidelity,
        "omega_b":         float(v["omega_b"]),
        "n_harmonics":     max(1, int(round(float(v["n_harmonics"])))),
        "theta_max":       _vec("theta_max"),
        "phi_angular":     [0.0] + [float(v.get(f"phi_angular_{i}", 0.0))
                                    for i in range(1, n_max)],
        "phi_horizontal":  _vec("phi_horizontal"),
        "omega_h":         float(v["omega_h"]),
        "amplitude_h":     _vec("amplitude_h"),
        "geometry": {
            "a": float(v["geometry_a"]),
            "b": float(v["geometry_b"]),
            "n": float(v["geometry_n"]),
        },
        "fill_level":      float(v["fill_level"]),
    }


def main(
    experiment_dir: str,
    param_space: str,
    model_path: str | None = None,
    kla_key: str = "kLa_25",
    lf_fidelity: int = 5,
    hf_fidelity: int = 7,
    n_candidates: int = 2000,
    seed: int = 0,
) -> dict:
    """Train surrogate (or load model_path) and return highest-EI candidate.

    Parameters
    ----------
    experiment_dir : path to f3dasm ExperimentData directory
    param_space    : path to config/param_space.yaml
    model_path     : if given, load existing model instead of training
    kla_key        : output column to maximise
    lf_fidelity    : fidelity tag for LF rows in ExperimentData
    hf_fidelity    : fidelity tag for the suggested next run
    n_candidates   : random candidates evaluated by acquisition function
    seed           : RNG seed for reproducibility

    Returns
    -------
    Nested params dict compatible with params.json (same structure).
    """
    spec = _load_param_space(param_space)

    if model_path is None:
        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            tmp_model = f.name
        _ts.main(experiment_dir, tmp_model,
                 lf_fidelity=lf_fidelity, hf_fidelity=hf_fidelity,
                 kla_key=kla_key)
        model_path = tmp_model

    with open(model_path, "rb") as f:
        model = pickle.load(f)

    # Best observed HF value
    from f3dasm import ExperimentData
    ed = ExperimentData.from_file(str(experiment_dir))
    inp_df, out_df = ed.to_pandas()
    hf_mask = inp_df["fidelity"] == hf_fidelity
    if hf_mask.sum() > 0:
        y_best = float(out_df.loc[hf_mask, kla_key].max())
    else:
        y_best = float(out_df[kla_key].max())

    rng        = np.random.default_rng(seed)
    X_cand     = _sample_candidates(spec, n_candidates, rng)
    mean, _var = model.predict(X_cand)
    std        = np.sqrt(np.maximum(_var, 0.0))
    ei         = _expected_improvement(mean.ravel(), std.ravel(), y_best)
    best_idx   = int(np.argmax(ei))

    return _flat_to_nested(X_cand[best_idx], spec, hf_fidelity)


if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment_dir")
    ap.add_argument("param_space")
    ap.add_argument("--model-path")
    ap.add_argument("--kla-key", default="kLa_25")
    ap.add_argument("--lf-fidelity", type=int, default=5)
    ap.add_argument("--hf-fidelity", type=int, default=7)
    ap.add_argument("--n-candidates", type=int, default=2000)
    args = ap.parse_args()
    result = main(args.experiment_dir, args.param_space,
                  model_path=args.model_path, kla_key=args.kla_key,
                  lf_fidelity=args.lf_fidelity, hf_fidelity=args.hf_fidelity,
                  n_candidates=args.n_candidates)
    print(json.dumps(result, indent=2))
