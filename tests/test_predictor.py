"""
Integration tests for the ResolutionPredictor.

Tests that the full propagation pipeline (B-field + material layers) produces
physically sensible track-resolution estimates.
"""

import numpy as np
import pytest

from weightmatrix.matrix import TrackState, N_PARAMS
from weightmatrix.bfield import BField
from weightmatrix.material import MaterialLayer
from weightmatrix.predictor import ResolutionPredictor


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_detector():
    """A simple 5-layer silicon barrel tracker in a 2 T solenoid."""
    bfield = BField(bz=2.0)
    layers = [
        (MaterialLayer(thickness=0.3, x0=93.7, name=f"Si_{i}"), 50.0 * (i + 1))
        for i in range(5)
    ]
    return bfield, layers


@pytest.fixture
def initial_state():
    """1 GeV/c track at the interaction point with generous initial uncertainties."""
    params = np.array([1.0, 0.0, 0.0, 0.0, 0.0])   # q/p=1, others=0
    variances = np.array([1e-4, 1e-3, 1e-3, 1.0, 1.0])
    return TrackState(params=params, covariance=np.diag(variances))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestResolutionPredictorBasic:
    def test_returns_one_state_per_layer(self, simple_detector, initial_state):
        bfield, layers = simple_detector
        predictor = ResolutionPredictor(bfield, layers)
        states = predictor.predict(initial_state, p_gev=1.0)
        assert len(states) == len(layers)

    def test_all_states_are_track_states(self, simple_detector, initial_state):
        bfield, layers = simple_detector
        predictor = ResolutionPredictor(bfield, layers)
        states = predictor.predict(initial_state, p_gev=1.0)
        for state in states:
            assert isinstance(state, TrackState)

    def test_covariance_remains_positive_definite(self, simple_detector, initial_state):
        """All eigenvalues of the covariance must be positive."""
        bfield, layers = simple_detector
        predictor = ResolutionPredictor(bfield, layers)
        states = predictor.predict(initial_state, p_gev=1.0)
        for state in states:
            eigvals = np.linalg.eigvalsh(state.covariance)
            assert np.all(eigvals > 0), f"Non-positive eigenvalues: {eigvals}"

    def test_covariance_is_symmetric(self, simple_detector, initial_state):
        bfield, layers = simple_detector
        predictor = ResolutionPredictor(bfield, layers)
        states = predictor.predict(initial_state, p_gev=1.0)
        for state in states:
            np.testing.assert_array_almost_equal(
                state.covariance, state.covariance.T, decimal=12
            )

    def test_sigmas_shape(self, simple_detector, initial_state):
        bfield, layers = simple_detector
        predictor = ResolutionPredictor(bfield, layers)
        states = predictor.predict(initial_state, p_gev=1.0)
        for state in states:
            assert state.sigmas.shape == (N_PARAMS,)


class TestResolutionPredictorPhysics:
    def test_higher_momentum_better_resolution(self):
        """Higher momentum → smaller scattering → better angular resolution."""
        bfield = BField(bz=2.0)
        layers = [
            (MaterialLayer(thickness=0.3, x0=93.7), 100.0)
        ]
        params = np.zeros(N_PARAMS)
        params[0] = 1.0  # q/p
        cov = np.diag([1e-6, 1e-4, 1e-4, 0.01, 0.01])

        predictor = ResolutionPredictor(bfield, layers)

        state_low_p = TrackState(params.copy(), covariance=cov.copy())
        state_high_p = TrackState(params.copy(), covariance=cov.copy())

        # q/p=2 → p=0.5 GeV, q/p=0.2 → p=5 GeV
        states_low = predictor.predict(state_low_p, p_gev=0.5)
        states_high = predictor.predict(state_high_p, p_gev=5.0)

        # sigma_lambda (idx 1) should be larger for lower momentum
        sigma_lam_low = states_low[-1].sigmas[1]
        sigma_lam_high = states_high[-1].sigmas[1]
        assert sigma_lam_low > sigma_lam_high

    def test_empty_layers(self, initial_state):
        """Predictor with no layers returns empty list."""
        bfield = BField()
        predictor = ResolutionPredictor(bfield, [])
        states = predictor.predict(initial_state, p_gev=1.0)
        assert states == []

    def test_zero_material_no_scattering_growth(self):
        """Vacuum layer: resolution only changes via B-field propagation."""
        bfield = BField(bz=0.0)  # no field either → pure drift
        vacuum = MaterialLayer(thickness=0.0, x0=93.7, name="vacuum")
        layers = [(vacuum, 100.0)]

        params = np.zeros(N_PARAMS)
        params[0] = 1.0
        initial_cov = np.diag([1e-6, 1e-4, 1e-4, 0.01, 0.01])

        predictor = ResolutionPredictor(bfield, layers)
        state0 = TrackState(params.copy(), covariance=initial_cov.copy())
        states = predictor.predict(state0, p_gev=1.0)

        # Angular resolutions (lambda, phi) should be unchanged (no scattering, no B)
        np.testing.assert_allclose(
            states[0].sigmas[1], np.sqrt(initial_cov[1, 1]), rtol=1e-6
        )
        np.testing.assert_allclose(
            states[0].sigmas[2], np.sqrt(initial_cov[2, 2]), rtol=1e-6
        )
