import math

import numpy as np

from npbc_water_analysis.geometry import (
    accessible_pair_shell_volume,
    radial_rdf_from_counts,
    sphere_intersection_volume,
    spherical_shell_volume,
)


def test_intersection_limiting_cases():
    assert sphere_intersection_volume(3.0, 1.0, 5.0) == 0.0
    assert math.isclose(
        sphere_intersection_volume(3.0, 1.0, 1.0),
        4.0 * math.pi / 3.0,
    )
    assert math.isclose(
        sphere_intersection_volume(2.0, 2.0, 0.0),
        4.0 * math.pi * 8.0 / 3.0,
    )


def test_intersection_is_symmetric():
    assert math.isclose(
        sphere_intersection_volume(4.0, 2.5, 3.0),
        sphere_intersection_volume(2.5, 4.0, 3.0),
        rel_tol=1.0e-14,
    )


def test_centered_accessible_shell_is_full_shell():
    observed = accessible_pair_shell_volume(15.0, 0.0, 2.0, 2.5)
    expected = spherical_shell_volume(2.0, 2.5)
    assert math.isclose(observed, expected, rel_tol=1.0e-14)


def test_accessible_shell_never_exceeds_full_shell():
    full = spherical_shell_volume(4.0, 4.5)
    for reference_radius in np.linspace(0.0, 15.5, 64):
        value = accessible_pair_shell_volume(15.0, reference_radius, 4.0, 4.5)
        assert 0.0 <= value <= full + 1.0e-12


def test_distinct_particle_normalization_returns_ideal_gas_unity():
    accessible = np.array([10.0, 20.0, 30.0])
    n_oxygen = 526
    radius = 15.0
    rho_distinct = (n_oxygen - 1) / (4.0 * math.pi * radius**3 / 3.0)
    counts = rho_distinct * accessible
    observed = radial_rdf_from_counts(counts, accessible, n_oxygen, radius)
    np.testing.assert_allclose(observed, 1.0)
