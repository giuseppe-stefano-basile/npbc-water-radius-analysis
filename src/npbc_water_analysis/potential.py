"""Frozen Gaussian mean-field reconstruction."""

from __future__ import annotations

import numpy as np


def gaussian_mean_field(
    r: np.ndarray,
    amplitudes: np.ndarray,
    centers: np.ndarray,
    widths: np.ndarray,
    scale_eV: float = 5.0e-4,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``U(r)`` and the analytic outward force ``-dU/dr`` in eV units."""
    r = np.asarray(r, dtype=float)
    a = np.asarray(amplitudes, dtype=float)
    c = np.asarray(centers, dtype=float)
    s = np.asarray(widths, dtype=float)
    if not (a.shape == c.shape == s.shape):
        raise ValueError("amplitudes, centers, and widths must have identical shapes")
    if a.ndim != 1 or not len(a):
        raise ValueError("at least one Gaussian layer is required")
    if np.any(s <= 0):
        raise ValueError("all Gaussian widths must be positive")
    delta = r[..., None] - c
    weights = np.exp(-0.5 * (delta / s) ** 2)
    shifted = a - a[0]
    potential = scale_eV * np.sum(shifted * weights, axis=-1)
    outward_force = scale_eV * np.sum(shifted * weights * delta / s**2, axis=-1)
    return potential, outward_force
