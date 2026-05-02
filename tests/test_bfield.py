"""
Tests for the B-field module (bfield.py).
"""

import numpy as np
import pytest

from weightmatrix.bfield import BField
from weightmatrix.matrix import N_PARAMS


class TestBFieldMagnitude:
    def test_default_is_2T(self):
        b = BField()
        assert b.magnitude == pytest.approx(2.0)

    def test_zero_field(self):
        b = BField(bx=0, by=0, bz=0)
        assert b.magnitude == pytest.approx(0.0)

    def test_combined_field(self):
        b = BField(bx=1.0, by=0.0, bz=0.0)
        assert b.magnitude == pytest.approx(1.0)

    def test_vector_shape(self):
        b = BField(bx=1.0, by=2.0, bz=3.0)
        assert b.vector.shape == (3,)
        np.testing.assert_array_equal(b.vector, [1.0, 2.0, 3.0])


class TestBFieldCurvature:
    def test_curvature_sign_positive_charge(self):
        b = BField()
        kappa = b.curvature(qop=1.0)  # q/p > 0 → positive curvature
        assert kappa > 0

    def test_curvature_sign_negative_charge(self):
        b = BField()
        kappa = b.curvature(qop=-1.0)
        assert kappa < 0

    def test_zero_field_zero_curvature(self):
        b = BField(bx=0, by=0, bz=0)
        assert b.curvature(qop=1.0) == pytest.approx(0.0)

    def test_curvature_value(self):
        # kappa = 0.3e-3 * qop * |B| = 0.3e-3 * 1.0 * 2.0 = 6e-4 /mm
        b = BField(bz=2.0)
        assert b.curvature(qop=1.0) == pytest.approx(6e-4)


class TestTransportJacobian:
    def test_jacobian_shape(self):
        b = BField()
        F = b.transport_jacobian_simple(qop=1.0, step_length=100.0)
        assert F.shape == (N_PARAMS, N_PARAMS)

    def test_zero_field_is_identity_like(self):
        """With B=0 the curvature-dependent off-diagonal terms should be zero."""
        b = BField(bx=0, by=0, bz=0)
        F = b.transport_jacobian_simple(qop=1.0, step_length=100.0)
        # Curvature-driven terms: d(phi)/d(qop) and d(y)/d(qop) must be zero
        assert F[2, 0] == pytest.approx(0.0)
        assert F[3, 0] == pytest.approx(0.0)

    def test_identity_diagonal(self):
        b = BField()
        F = b.transport_jacobian_simple(qop=1.0, step_length=100.0)
        np.testing.assert_array_equal(np.diag(F), np.ones(N_PARAMS))

    def test_jacobian_nonzero_offdiag_with_field(self):
        b = BField(bz=2.0)
        F = b.transport_jacobian_simple(qop=1.0, step_length=100.0)
        # With non-zero field there should be off-diagonal terms
        off_diag_sum = np.sum(np.abs(F - np.eye(N_PARAMS)))
        assert off_diag_sum > 0
