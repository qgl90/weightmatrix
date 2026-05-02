"""
Core weight matrix formalism for track-state resolution prediction.

The weight matrix W is the inverse of the covariance matrix C (W = C^{-1}).
Track states are represented in a 5-parameter local coordinate system:
  (q/p, lambda, phi, y, z)
where:
  q/p    - charge over momentum [1/GeV]
  lambda - dip angle [rad]
  phi    - azimuthal angle [rad]
  y      - local transverse position [mm]
  z      - local longitudinal position [mm]
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


# Indices for the 5 track parameters
IDX_QOP = 0   # q/p  (charge over momentum)
IDX_LAM = 1   # lambda (dip angle)
IDX_PHI = 2   # phi (azimuthal angle)
IDX_Y = 3     # local y
IDX_Z = 4     # local z

N_PARAMS = 5


@dataclass
class TrackState:
    """
    A local 5-parameter track state with an associated weight (inverse covariance) matrix.

    Parameters
    ----------
    params : array-like of shape (5,)
        Track parameters [q/p, lambda, phi, y, z].
    covariance : array-like of shape (5, 5), optional
        Covariance matrix C. Mutually exclusive with ``weight``.
    weight : array-like of shape (5, 5), optional
        Weight matrix W = C^{-1}.  Mutually exclusive with ``covariance``.

    If neither is provided the weight matrix defaults to the 5×5 identity matrix.
    """

    params: np.ndarray
    _weight: np.ndarray = field(init=False, repr=False)

    def __init__(
        self,
        params,
        covariance=None,
        weight=None,
    ):
        self.params = np.asarray(params, dtype=float)
        if self.params.shape != (N_PARAMS,):
            raise ValueError(
                f"params must have shape ({N_PARAMS},), got {self.params.shape}"
            )

        if covariance is not None and weight is not None:
            raise ValueError("Provide either covariance or weight, not both.")

        if covariance is not None:
            cov = np.asarray(covariance, dtype=float)
            if cov.shape != (N_PARAMS, N_PARAMS):
                raise ValueError(
                    f"covariance must have shape ({N_PARAMS}, {N_PARAMS})."
                )
            self._weight = np.linalg.inv(cov)
        elif weight is not None:
            w = np.asarray(weight, dtype=float)
            if w.shape != (N_PARAMS, N_PARAMS):
                raise ValueError(
                    f"weight must have shape ({N_PARAMS}, {N_PARAMS})."
                )
            self._weight = w.copy()
        else:
            self._weight = np.eye(N_PARAMS)

    @property
    def weight(self) -> np.ndarray:
        """Weight matrix W = C^{-1}."""
        return self._weight

    @property
    def covariance(self) -> np.ndarray:
        """Covariance matrix C = W^{-1}."""
        return np.linalg.inv(self._weight)

    @property
    def sigmas(self) -> np.ndarray:
        """Standard deviations (sqrt of diagonal covariance elements)."""
        return np.sqrt(np.diag(self.covariance))

    def __repr__(self) -> str:
        return (
            f"TrackState(params={self.params}, "
            f"sigma_qop={self.sigmas[IDX_QOP]:.4g}, "
            f"sigma_y={self.sigmas[IDX_Y]:.4g} mm, "
            f"sigma_z={self.sigmas[IDX_Z]:.4g} mm)"
        )


class WeightMatrix:
    """
    Utility class for constructing and manipulating weight matrices.

    The weight matrix encapsulates the precision of a track measurement.
    It is the inverse of the covariance matrix.
    """

    def __init__(self, matrix: Optional[np.ndarray] = None):
        if matrix is None:
            self._w = np.eye(N_PARAMS)
        else:
            w = np.asarray(matrix, dtype=float)
            if w.shape != (N_PARAMS, N_PARAMS):
                raise ValueError(
                    f"matrix must have shape ({N_PARAMS}, {N_PARAMS}), "
                    f"got {w.shape}"
                )
            self._w = w.copy()

    @classmethod
    def from_covariance(cls, covariance) -> "WeightMatrix":
        """Construct a WeightMatrix from a covariance matrix C."""
        cov = np.asarray(covariance, dtype=float)
        return cls(np.linalg.inv(cov))

    @classmethod
    def diagonal(cls, variances) -> "WeightMatrix":
        """Construct a WeightMatrix from a vector of variances (diagonal C)."""
        variances = np.asarray(variances, dtype=float)
        if variances.shape != (N_PARAMS,):
            raise ValueError(
                f"variances must have shape ({N_PARAMS},), got {variances.shape}"
            )
        return cls(np.diag(1.0 / variances))

    @property
    def matrix(self) -> np.ndarray:
        """The 5×5 weight matrix array."""
        return self._w

    @property
    def covariance(self) -> np.ndarray:
        """The 5×5 covariance matrix C = W^{-1}."""
        return np.linalg.inv(self._w)

    def add(self, other: "WeightMatrix") -> "WeightMatrix":
        """
        Return the combined weight matrix after adding information from *other*.

        In the weight formalism the combined precision is simply W_total = W_1 + W_2
        (equivalent to combining two independent measurements).
        """
        return WeightMatrix(self._w + other._w)

    def __add__(self, other: "WeightMatrix") -> "WeightMatrix":
        return self.add(other)

    def __repr__(self) -> str:
        return f"WeightMatrix(\n{self._w}\n)"


def propagate_weight_matrix(
    weight: np.ndarray,
    jacobian: np.ndarray,
) -> np.ndarray:
    """
    Propagate a weight matrix through a linear transformation.

    Given the Jacobian F of the propagation (from state k to state k+1),
    the new weight matrix is:

        W_{k+1} = (F · C_k · F^T)^{-1}
                = (F^{-T}) · W_k · (F^{-1})

    where C_k = W_k^{-1} is the covariance at step k.

    Parameters
    ----------
    weight : ndarray of shape (5, 5)
        Weight matrix at the current surface.
    jacobian : ndarray of shape (5, 5)
        Transport Jacobian F from current to next surface.

    Returns
    -------
    ndarray of shape (5, 5)
        Weight matrix at the next surface.
    """
    weight = np.asarray(weight, dtype=float)
    F = np.asarray(jacobian, dtype=float)
    if weight.shape != (N_PARAMS, N_PARAMS) or F.shape != (N_PARAMS, N_PARAMS):
        raise ValueError("Both weight and jacobian must be 5×5 matrices.")

    cov = np.linalg.inv(weight)
    new_cov = F @ cov @ F.T
    return np.linalg.inv(new_cov)
