"""
High-level detector resolution predictor.

Combines weight matrix propagation, B-field transport, and material
scattering to parametrically predict track-state resolutions through
a sequence of detector layers.
"""

from __future__ import annotations

import numpy as np
from typing import List, Tuple

from .matrix import TrackState, WeightMatrix, propagate_weight_matrix, N_PARAMS
from .bfield import BField
from .material import MaterialLayer


class ResolutionPredictor:
    """
    Predict track-parameter resolutions propagating through detector layers.

    Parameters
    ----------
    bfield : BField
        Magnetic field configuration.
    layers : list of (MaterialLayer, float)
        Sequence of (layer, step_length_mm) tuples describing the detector.
        Each tuple gives the material layer and the propagation distance
        (arc-length in mm) to reach that layer.
    """

    def __init__(self, bfield: BField, layers: List[Tuple[MaterialLayer, float]]):
        self.bfield = bfield
        self.layers = layers

    def predict(
        self,
        initial_state: TrackState,
        p_gev: float,
        q: float = 1.0,
        beta: float = 1.0,
    ) -> List[TrackState]:
        """
        Propagate *initial_state* through all detector layers and return the
        predicted track state (with updated weight/covariance) at each layer.

        Parameters
        ----------
        initial_state : TrackState
            Starting track state (e.g. at the interaction point).
        p_gev : float
            Particle momentum [GeV].
        q : float
            Charge number.
        beta : float
            Relativistic beta.

        Returns
        -------
        list of TrackState
            One TrackState per layer, in order.
        """
        states: List[TrackState] = []
        current_weight = initial_state.weight.copy()
        current_params = initial_state.params.copy()
        qop = current_params[0]
        lambda_ = current_params[1]

        for layer, step_mm in self.layers:
            # --- Transport Jacobian (B-field deflection) ---
            F = self.bfield.transport_jacobian_simple(
                qop=qop, step_length=step_mm, lambda_=lambda_
            )
            current_weight = propagate_weight_matrix(current_weight, F)

            # --- Multiple scattering: add noise Q to covariance, re-invert ---
            Q = layer.scattering_weight_correction(p_gev, q=q, beta=beta)
            cov = np.linalg.inv(current_weight)
            cov_new = cov + Q
            current_weight = np.linalg.inv(cov_new)

            # Update params (simple straight-line update, angles shift slightly)
            current_params = F @ current_params

            states.append(
                TrackState(params=current_params, weight=current_weight)
            )

        return states
