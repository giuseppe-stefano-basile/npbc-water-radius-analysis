import math

import numpy as np

from npbc_water_analysis.statistics import (
    distribution_distances,
    first_temporal_sd,
    molecule_weighted_distribution,
    standardized_mean_difference,
    volume_weighted_distribution,
)


def test_weighting_conventions_are_intentionally_different():
    counts = np.array([[90.0, 10.0], [1.0, 9.0]])
    molecule = molecule_weighted_distribution(counts)
    volume = volume_weighted_distribution(counts)
    np.testing.assert_allclose(molecule, [91.0 / 110.0, 19.0 / 110.0])
    np.testing.assert_allclose(volume, [0.5, 0.5])


def test_identical_distribution_distances_are_zero():
    p = np.array([0.1, 0.3, 0.6])
    js, w1 = distribution_distances(p, p, np.array([0.0, 1.0, 2.0]))
    assert js == 0.0
    assert w1 == 0.0


def test_temporal_sd_uses_sample_definition():
    assert math.isclose(first_temporal_sd(np.array([1.0, 2.0, 3.0])), 1.0)


def test_standardized_mean_difference():
    assert math.isclose(standardized_mean_difference(1.2, 1.0, 1.0, 1.0), 0.2)
