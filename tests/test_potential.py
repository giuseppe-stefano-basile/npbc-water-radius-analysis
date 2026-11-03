import numpy as np

from npbc_water_analysis.potential import gaussian_mean_field


def test_analytic_force_matches_negative_numerical_derivative():
    r = np.linspace(0.1, 9.9, 5001)
    amplitudes = np.array([0.0, 2.0, -1.0])
    centers = np.array([2.0, 5.0, 8.0])
    widths = np.array([0.7, 1.2, 0.9])
    potential, force = gaussian_mean_field(r, amplitudes, centers, widths)
    numerical = -np.gradient(potential, r)
    np.testing.assert_allclose(force[5:-5], numerical[5:-5], rtol=2.0e-4, atol=1.0e-8)
