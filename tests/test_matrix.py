"""
Tests for the core weight matrix formalism (matrix.py).
"""

import numpy as np
import pytest

from weightmatrix.matrix import (
    WeightMatrix,
    TrackState,
    propagate_weight_matrix,
    N_PARAMS,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def identity_weight():
    return WeightMatrix()


@pytest.fixture
def diagonal_variances():
    return np.array([1e-6, 1e-4, 1e-4, 0.01, 0.01])


@pytest.fixture
def diagonal_weight(diagonal_variances):
    return WeightMatrix.diagonal(diagonal_variances)


@pytest.fixture
def simple_track_state():
    params = np.array([0.5, 0.1, 0.2, 0.0, 0.0])   # q/p=0.5, lambda=0.1, phi=0.2
    variances = np.array([1e-6, 1e-4, 1e-4, 0.04, 0.04])
    return TrackState(params=params, covariance=np.diag(variances))


# ---------------------------------------------------------------------------
# WeightMatrix tests
# ---------------------------------------------------------------------------

class TestWeightMatrixIdentity:
    def test_default_is_identity(self, identity_weight):
        np.testing.assert_array_equal(identity_weight.matrix, np.eye(N_PARAMS))

    def test_covariance_of_identity_is_identity(self, identity_weight):
        np.testing.assert_array_almost_equal(
            identity_weight.covariance, np.eye(N_PARAMS)
        )


class TestWeightMatrixFromCovariance:
    def test_roundtrip(self, diagonal_variances):
        cov = np.diag(diagonal_variances)
        wm = WeightMatrix.from_covariance(cov)
        np.testing.assert_array_almost_equal(wm.covariance, cov)

    def test_inverse_relation(self, diagonal_variances):
        cov = np.diag(diagonal_variances)
        wm = WeightMatrix.from_covariance(cov)
        product = wm.matrix @ cov
        np.testing.assert_array_almost_equal(product, np.eye(N_PARAMS), decimal=10)


class TestWeightMatrixDiagonal:
    def test_diagonal_entries(self, diagonal_variances):
        wm = WeightMatrix.diagonal(diagonal_variances)
        expected_diag = 1.0 / diagonal_variances
        np.testing.assert_array_almost_equal(np.diag(wm.matrix), expected_diag)

    def test_off_diagonal_zero(self, diagonal_variances):
        wm = WeightMatrix.diagonal(diagonal_variances)
        m = wm.matrix
        off_diag = m - np.diag(np.diag(m))
        np.testing.assert_array_equal(off_diag, np.zeros((N_PARAMS, N_PARAMS)))

    def test_wrong_size_raises(self):
        with pytest.raises(ValueError):
            WeightMatrix.diagonal(np.array([1.0, 2.0, 3.0]))


class TestWeightMatrixAdd:
    def test_identity_plus_identity(self, identity_weight):
        result = identity_weight + identity_weight
        np.testing.assert_array_equal(result.matrix, 2 * np.eye(N_PARAMS))

    def test_addition_is_commutative(self, identity_weight, diagonal_weight):
        ab = identity_weight + diagonal_weight
        ba = diagonal_weight + identity_weight
        np.testing.assert_array_almost_equal(ab.matrix, ba.matrix)


class TestWeightMatrixConstructorValidation:
    def test_wrong_shape_raises(self):
        with pytest.raises(ValueError):
            WeightMatrix(np.eye(3))

    def test_non_square_raises(self):
        with pytest.raises(ValueError):
            WeightMatrix(np.ones((5, 3)))


# ---------------------------------------------------------------------------
# TrackState tests
# ---------------------------------------------------------------------------

class TestTrackState:
    def test_default_weight_is_identity(self):
        params = np.zeros(N_PARAMS)
        ts = TrackState(params=params)
        np.testing.assert_array_equal(ts.weight, np.eye(N_PARAMS))

    def test_covariance_roundtrip(self, simple_track_state):
        cov = np.diag([1e-6, 1e-4, 1e-4, 0.04, 0.04])
        np.testing.assert_array_almost_equal(simple_track_state.covariance, cov)

    def test_sigmas_shape(self, simple_track_state):
        assert simple_track_state.sigmas.shape == (N_PARAMS,)

    def test_sigmas_values(self, simple_track_state):
        expected = np.sqrt([1e-6, 1e-4, 1e-4, 0.04, 0.04])
        np.testing.assert_array_almost_equal(simple_track_state.sigmas, expected)

    def test_weight_from_weight(self):
        w = 3 * np.eye(N_PARAMS)
        ts = TrackState(params=np.zeros(N_PARAMS), weight=w)
        np.testing.assert_array_almost_equal(ts.weight, w)

    def test_both_cov_and_weight_raises(self):
        with pytest.raises(ValueError, match="either"):
            TrackState(
                params=np.zeros(N_PARAMS),
                covariance=np.eye(N_PARAMS),
                weight=np.eye(N_PARAMS),
            )

    def test_wrong_param_shape_raises(self):
        with pytest.raises(ValueError):
            TrackState(params=np.zeros(3))

    def test_repr(self, simple_track_state):
        r = repr(simple_track_state)
        assert "TrackState" in r


# ---------------------------------------------------------------------------
# propagate_weight_matrix tests
# ---------------------------------------------------------------------------

class TestPropagateWeightMatrix:
    def test_identity_jacobian_unchanged(self):
        W = 4 * np.eye(N_PARAMS)
        F = np.eye(N_PARAMS)
        W_new = propagate_weight_matrix(W, F)
        np.testing.assert_array_almost_equal(W_new, W)

    def test_scale_jacobian(self):
        """Scaling all coordinates by factor s should scale weight by 1/s^2."""
        s = 2.0
        W = np.eye(N_PARAMS)
        F = s * np.eye(N_PARAMS)
        W_new = propagate_weight_matrix(W, F)
        expected = (1.0 / s**2) * np.eye(N_PARAMS)
        np.testing.assert_array_almost_equal(W_new, expected)

    def test_wrong_shapes_raise(self):
        with pytest.raises(ValueError):
            propagate_weight_matrix(np.eye(3), np.eye(5))

    def test_roundtrip_through_inverse_jacobian(self):
        """Propagating forward then backward should return to original W."""
        rng = np.random.default_rng(42)
        # Random symmetric positive definite W
        A = rng.standard_normal((N_PARAMS, N_PARAMS))
        W = A @ A.T + np.eye(N_PARAMS)

        F = np.eye(N_PARAMS) + 0.01 * rng.standard_normal((N_PARAMS, N_PARAMS))
        F_inv = np.linalg.inv(F)

        W_fwd = propagate_weight_matrix(W, F)
        W_back = propagate_weight_matrix(W_fwd, F_inv)
        np.testing.assert_array_almost_equal(W_back, W, decimal=10)
