"""Small, transparent statistical kernels used by the teaching material."""

from __future__ import annotations

import numpy as np
from scipy.spatial.distance import jensenshannon
from scipy.stats import wasserstein_distance


def normalize_histogram(counts: np.ndarray) -> np.ndarray:
    values = np.asarray(counts, dtype=float)
    if np.any(values < 0):
        raise ValueError("histogram counts must be non-negative")
    total = values.sum()
    if total <= 0:
        raise ValueError("histogram must contain positive mass")
    return values / total


def volume_weighted_distribution(shell_probabilities: np.ndarray) -> np.ndarray:
    """Average normalized per-shell distributions with equal shell weight."""
    values = np.asarray(shell_probabilities, dtype=float)
    if values.ndim != 2:
        raise ValueError("expected shape (shell, bin)")
    rows = np.array([normalize_histogram(row) for row in values])
    return normalize_histogram(rows.mean(axis=0))


def molecule_weighted_distribution(shell_counts: np.ndarray) -> np.ndarray:
    """Pool raw counts, making each molecule contribute equally."""
    values = np.asarray(shell_counts, dtype=float)
    if values.ndim != 2:
        raise ValueError("expected shape (shell, bin)")
    return normalize_histogram(values.sum(axis=0))


def first_temporal_sd(values: np.ndarray, axis: int = 0) -> np.ndarray:
    return np.nanstd(np.asarray(values, dtype=float), axis=axis, ddof=1)


def distribution_distances(p: np.ndarray, q: np.ndarray, bin_midpoints: np.ndarray) -> tuple[float, float]:
    p = normalize_histogram(p)
    q = normalize_histogram(q)
    x = np.asarray(bin_midpoints, dtype=float)
    if x.shape != p.shape:
        raise ValueError("bin_midpoints must match the histogram shape")
    js = float(jensenshannon(p, q, base=2.0))
    w1 = float(wasserstein_distance(x, x, u_weights=p, v_weights=q))
    return js, w1


def standardized_mean_difference(mu_a: float, mu_b: float, sd_a: float, sd_b: float) -> float:
    pooled = np.sqrt((sd_a**2 + sd_b**2) / 2.0)
    if pooled == 0:
        return 0.0 if mu_a == mu_b else float(np.sign(mu_a - mu_b) * np.inf)
    return float((mu_a - mu_b) / pooled)
