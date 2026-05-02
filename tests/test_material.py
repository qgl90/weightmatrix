"""
Tests for the material layer module (material.py).
"""

import numpy as np
import pytest

from weightmatrix.material import MaterialLayer
from weightmatrix.matrix import N_PARAMS


SILICON_X0 = 93.7     # mm
SILICON_THICKNESS = 0.3  # mm


@pytest.fixture
def silicon_layer():
    return MaterialLayer(thickness=SILICON_THICKNESS, x0=SILICON_X0, name="Si_layer")


class TestMaterialLayerProperties:
    def test_x_over_x0(self, silicon_layer):
        expected = SILICON_THICKNESS / SILICON_X0
        assert silicon_layer.x_over_x0 == pytest.approx(expected)

    def test_repr(self, silicon_layer):
        r = repr(silicon_layer)
        assert "Si_layer" in r
        assert "MaterialLayer" in r


class TestHighlandFormula:
    def test_theta0_positive(self, silicon_layer):
        theta0 = silicon_layer.highland_theta0(p_gev=1.0)
        assert theta0 > 0

    def test_theta0_decreases_with_momentum(self, silicon_layer):
        """Higher momentum → less scattering."""
        theta_low = silicon_layer.highland_theta0(p_gev=0.5)
        theta_high = silicon_layer.highland_theta0(p_gev=5.0)
        assert theta_low > theta_high

    def test_zero_thickness(self):
        layer = MaterialLayer(thickness=0.0, x0=SILICON_X0)
        assert layer.highland_theta0(p_gev=1.0) == pytest.approx(0.0)

    def test_zero_momentum(self, silicon_layer):
        assert silicon_layer.highland_theta0(p_gev=0.0) == pytest.approx(0.0)

    def test_theta0_scales_with_charge(self, silicon_layer):
        """Double charge → double scattering angle."""
        t1 = silicon_layer.highland_theta0(p_gev=1.0, q=1.0)
        t2 = silicon_layer.highland_theta0(p_gev=1.0, q=2.0)
        assert t2 == pytest.approx(2 * t1)

    def test_theta0_value_1gev_silicon(self, silicon_layer):
        """Cross-check against known value for 1 GeV/c through 300 µm of Si."""
        t = SILICON_THICKNESS / SILICON_X0
        expected = (13.6e-3 / 1.0) * np.sqrt(t) * (1.0 + 0.038 * np.log(t))
        assert silicon_layer.highland_theta0(p_gev=1.0) == pytest.approx(expected)


class TestScatteringWeightCorrection:
    def test_returns_5x5(self, silicon_layer):
        Q = silicon_layer.scattering_weight_correction(p_gev=1.0)
        assert Q.shape == (N_PARAMS, N_PARAMS)

    def test_only_angular_elements_nonzero(self, silicon_layer):
        Q = silicon_layer.scattering_weight_correction(p_gev=1.0)
        # Only (1,1) and (2,2) should be non-zero
        for i in range(N_PARAMS):
            for j in range(N_PARAMS):
                if (i, j) in [(1, 1), (2, 2)]:
                    assert Q[i, j] > 0
                else:
                    assert Q[i, j] == pytest.approx(0.0)

    def test_symmetry(self, silicon_layer):
        Q = silicon_layer.scattering_weight_correction(p_gev=1.0)
        np.testing.assert_array_equal(Q, Q.T)

    def test_values_decrease_with_momentum(self, silicon_layer):
        Q1 = silicon_layer.scattering_weight_correction(p_gev=0.5)
        Q10 = silicon_layer.scattering_weight_correction(p_gev=10.0)
        assert Q1[1, 1] > Q10[1, 1]
