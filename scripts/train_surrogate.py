"""Train a multi-fidelity KRR-LR-GPR surrogate from an f3dasm ExperimentData store.

Public API
----------
main(experiment_dir, model_path, lf_fidelity, hf_fidelity, kla_key) -> None
    Load ExperimentData, split by fidelity, train KRR-LR-GPR, pickle a
    wrapper whose .predict(X) -> (mean, var) matches the test interface.
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
from f3dasm import ExperimentData

# Use vendored copy (numpy-only, no torch) with numpy-compat fixes applied.
sys.path.insert(0, str(Path(__file__).parent))
from mfbml_local.krr_lr_gpr import KernelRidgeLinearGaussianProcess

# Input columns present in ExperimentData; 'fidelity' and 'phi_angular_0' are
# excluded from the surrogate feature vector (tracking column and fixed zero).
_ALL_INPUT_COLS = [
    "omega_b", "n_harmonics",
    "theta_max_0", "theta_max_1", "theta_max_2",
    "phi_angular_0", "phi_angular_1", "phi_angular_2",
    "omega_h",
    "amplitude_h_0", "amplitude_h_1", "amplitude_h_2",
    "phi_horizontal_0", "phi_horizontal_1", "phi_horizontal_2",
    "geometry_a", "geometry_b", "geometry_n",
    "fill_level", "fidelity",
]
_FEATURE_COLS = [c for c in _ALL_INPUT_COLS
                 if c not in ("phi_angular_0", "fidelity")]  # 18 features


class _WrappedModel:
    """Thin wrapper so pickle gives back a stable .predict(X) -> (mean, var)."""

    def __init__(self, model: KernelRidgeLinearGaussianProcess,
                 design_space: np.ndarray) -> None:
        self._model        = model
        self._design_space = design_space

    def predict(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return (mean, variance) arrays of shape (n_samples,)."""
        mean, std = self._model.predict(X, return_std=True)
        return mean.ravel(), (std.ravel() ** 2)


def main(
    experiment_dir: str,
    model_path: str,
    lf_fidelity: int = 5,
    hf_fidelity: int = 7,
    kla_key: str = "kLa_25",
) -> None:
    """Load ExperimentData, train KRR-LR-GPR, save pickled _WrappedModel.

    Parameters
    ----------
    experiment_dir : path to the f3dasm ExperimentData directory
    model_path     : where to write model.pkl
    lf_fidelity    : fidelity value identifying low-fidelity rows
    hf_fidelity    : fidelity value identifying high-fidelity rows
    kla_key        : output column to use as the surrogate target
    """
    ed = ExperimentData.from_file(str(experiment_dir))
    inp_df, out_df = ed.to_pandas()

    # Build full feature matrix from the 18 surrogate columns
    missing = [c for c in _FEATURE_COLS if c not in inp_df.columns]
    if missing:
        raise ValueError(f"ExperimentData missing input columns: {missing}")
    X_all = inp_df[_FEATURE_COLS].to_numpy(dtype=float)

    fidelity_col = inp_df["fidelity"].to_numpy(dtype=float)
    y_all        = out_df[kla_key].to_numpy(dtype=float)

    lf_mask = fidelity_col == lf_fidelity
    hf_mask = fidelity_col == hf_fidelity

    if lf_mask.sum() == 0:
        raise ValueError(f"No rows with lf_fidelity={lf_fidelity} in experiment data")
    if hf_mask.sum() == 0:
        raise ValueError(f"No rows with hf_fidelity={hf_fidelity} in experiment data")

    X_hf, y_hf = X_all[hf_mask], y_all[hf_mask]
    X_lf, y_lf = X_all[lf_mask], y_all[lf_mask]

    # design_space: (n_features, 2) with [min, max] per column
    lo = X_all.min(axis=0)
    hi = X_all.max(axis=0)
    # Guard against degenerate columns (all same value)
    degenerate = lo == hi
    hi[degenerate] = lo[degenerate] + 1.0
    design_space = np.column_stack([lo, hi])

    model = KernelRidgeLinearGaussianProcess(
        design_space=design_space,
        optimizer_restart=5,
        seed=42,
    )
    # API: X[0]=HF, X[1]=LF; Y[0]=HF responses, Y[1]=LF responses
    model.train(X=[X_hf, X_lf], Y=[y_hf, y_lf])

    wrapped = _WrappedModel(model, design_space)
    with open(model_path, "wb") as f:
        pickle.dump(wrapped, f)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("experiment_dir")
    ap.add_argument("model_path")
    ap.add_argument("--lf-fidelity", type=int, default=5)
    ap.add_argument("--hf-fidelity", type=int, default=7)
    ap.add_argument("--kla-key", default="kLa_25")
    args = ap.parse_args()
    main(args.experiment_dir, args.model_path,
         lf_fidelity=args.lf_fidelity, hf_fidelity=args.hf_fidelity,
         kla_key=args.kla_key)
