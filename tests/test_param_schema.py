"""Tests for parameter schema validation."""
import json
from pathlib import Path

import pytest

from scripts.postprocess import validate_params  # shared validator

VALID_PARAMS = {
    "fidelity": 7,
    "omega_b": 3.93,
    "n_harmonics": 2,
    "theta_max": [7.0, 2.0, 0.0],
    "phi_angular": [0.0, 1.57, 0.0],
    "omega_h": 2.1,
    "amplitude_h": [0.01, 0.005, 0.0],
    "phi_horizontal": [0.0, 0.78, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 2.0},
    "fill_level": 0.5,
}


def test_valid_params_pass_validation():
    validate_params(VALID_PARAMS)  # should not raise


def test_out_of_bounds_param_rejected():
    bad = {**VALID_PARAMS, "omega_b": 99.0}
    with pytest.raises(ValueError, match="omega_b"):
        validate_params(bad)


def test_phi_angular_first_element_must_be_zero():
    bad = {**VALID_PARAMS, "phi_angular": [1.0, 1.57, 0.0]}
    with pytest.raises(ValueError, match="phi_angular"):
        validate_params(bad)


def test_vectors_must_have_length_nmax():
    bad = {**VALID_PARAMS, "theta_max": [7.0, 2.0]}  # length 2, not N_max=3
    with pytest.raises(ValueError, match="theta_max"):
        validate_params(bad)
