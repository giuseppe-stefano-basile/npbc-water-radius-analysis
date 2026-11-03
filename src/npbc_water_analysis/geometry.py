"""Independent geometry kernels used to explain and test the estimators."""

from __future__ import annotations

import math

import numpy as np


def spherical_shell_volume(r_inner: float, r_outer: float) -> float:
    if r_inner < 0 or r_outer < r_inner:
        raise ValueError("require 0 <= r_inner <= r_outer")
    return 4.0 * math.pi * (r_outer**3 - r_inner**3) / 3.0


def sphere_intersection_volume(R: float, s: float, d: float) -> float:
    """Volume shared by spheres of radii ``R`` and ``s`` at separation ``d``."""
    if min(R, s, d) < 0:
        raise ValueError("sphere radii and center separation must be non-negative")
    if d >= R + s:
        return 0.0
    if d <= abs(R - s):
        return 4.0 * math.pi * min(R, s) ** 3 / 3.0
    numerator = math.pi * (R + s - d) ** 2
    bracket = d**2 + 2.0 * d * (s + R) - 3.0 * (s - R) ** 2
    return numerator * bracket / (12.0 * d)


def accessible_pair_shell_volume(
    cavity_radius: float,
    reference_radius: float,
    r_inner: float,
    r_outer: float,
) -> float:
    """Exact pair-shell volume inside a spherical cavity."""
    outer = sphere_intersection_volume(cavity_radius, r_outer, reference_radius)
    inner = sphere_intersection_volume(cavity_radius, r_inner, reference_radius)
    value = outer - inner
    if value < -1.0e-10:
        raise ArithmeticError("accessible shell volume became negative")
    return max(0.0, value)


def radial_rdf_from_counts(
    counts: np.ndarray,
    accessible_volumes: np.ndarray,
    number_of_oxygens: int,
    cavity_radius: float,
) -> np.ndarray:
    """Normalize directed O--O counts with the exact distinct-particle density."""
    counts = np.asarray(counts, dtype=float)
    accessible_volumes = np.asarray(accessible_volumes, dtype=float)
    if counts.shape != accessible_volumes.shape:
        raise ValueError("counts and accessible_volumes must have identical shapes")
    if number_of_oxygens < 2:
        raise ValueError("at least two oxygens are required")
    cavity_volume = 4.0 * math.pi * cavity_radius**3 / 3.0
    rho_distinct = (number_of_oxygens - 1) / cavity_volume
    denominator = rho_distinct * accessible_volumes
    return np.divide(counts, denominator, out=np.full_like(counts, np.nan), where=denominator > 0)
