"""Smallest concentration-streak width that the flow makes before diffusion erases it, vs the L10 cell size.

A streak of width w in a flow with strain rate gamma thins as w ~ exp(-gamma t); diffusion widens it as
sqrt(D t). The two balance at the Batchelor width lambda = sqrt(D / gamma). The strain rate comes from the
measured dissipation: EDR = mu <2 S:S> = mu gamma^2, so gamma = sqrt(EDR / mu).
Constants: src/BioReactor.c (D_oxy_1 = 1.90e-9 m^2/s, mu_w = 1.0e-3 Pa s, L_bio = 0.25 m, NN = 2^level on
L0 = 1 L_bio, uniform grid). EDR: the L10 points of test_h3_mfbml.json (time-mean domain-mean, W/m^3).
Usage: uv run python scripts/diffusion_length.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
D, MU, L_BIO = 1.90e-9, 1.0e-3, 0.25


def main():
    pts = json.loads((ROOT / "experiments/multifidelity/test_h3_mfbml.json").read_text())["results"]["edr_mean_t"]["points"]
    dx = L_BIO / 2 ** 10
    print(f"L10 cell: dx = {L_BIO} m / 2^10 = {dx * 1e3:.3f} mm")
    print(f"{'rpm':>5} {'theta':>5} {'EDR W/m3':>9} {'gamma 1/s':>9} {'lambda mm':>9} {'dx/lambda':>9}")
    lam = []
    for p in pts:
        g = np.sqrt(p["obs"] / MU)
        lb = np.sqrt(D / g)
        lam.append(lb)
        print(f"{p['rpm']:5g} {p['theta']:5g} {p['obs']:9.4f} {g:9.2f} {lb * 1e3:9.4f} {dx / lb:9.1f}")
    lam = np.array(lam)
    print(f"lambda range {lam.min() * 1e3:.3f}-{lam.max() * 1e3:.3f} mm; dx/lambda {dx / lam.max():.1f}-{dx / lam.min():.1f}")


if __name__ == "__main__":
    main()
