"""
Magnetic field (B-field) representation and track curvature helper.

The B-field is represented as a 3-component vector (Bx, By, Bz) in Tesla.
Helper utilities compute the track curvature in the field and the resulting
contribution to the transport Jacobian.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass


@dataclass
class BField:
    """
    A uniform magnetic field vector B = (Bx, By, Bz) in Tesla.

    Parameters
    ----------
    bx, by, bz : float
        Field components in Tesla.
    """

    bx: float = 0.0
    by: float = 0.0
    bz: float = 2.0  # default: 2 T along z (typical LHC-like solenoid)

    @property
    def magnitude(self) -> float:
        """Field magnitude |B| in Tesla."""
        return float(np.sqrt(self.bx**2 + self.by**2 + self.bz**2))

    @property
    def vector(self) -> np.ndarray:
        """Field as a 3-element numpy array [Bx, By, Bz]."""
        return np.array([self.bx, self.by, self.bz])

    def curvature(self, qop: float) -> float:
        """
        Compute the track curvature kappa = q·|B|/p in 1/mm.

        Parameters
        ----------
        qop : float
            Charge over momentum q/p in 1/GeV.

        Returns
        -------
        float
            Signed curvature kappa [1/mm].  Uses the relation:
            kappa = (0.3 · q · |B|) / p
            where 0.3 = c / (10^9) converts SI units so that
            p is in GeV and B in Tesla gives kappa in 1/m;
            then / 1000 → 1/mm.
        """
        # 0.3 T·m/GeV factor (qop already carries the sign)
        return 0.3e-3 * qop * self.magnitude  # 1/mm

    def transport_jacobian_simple(
        self,
        qop: float,
        step_length: float,
        lambda_: float = 0.0,
    ) -> np.ndarray:
        """
        Return a simplified 5×5 transport Jacobian for a step in a uniform B-field.

        This uses a first-order approximation valid for small deflections.
        The full Jacobian would require a proper Runge-Kutta stepper; this
        linear approximation is sufficient for testing and parametric studies.

        Parameters
        ----------
        qop : float
            Charge over momentum [1/GeV].
        step_length : float
            Arc length of the step [mm].
        lambda_ : float
            Dip angle [rad] (default 0 = transverse track).

        Returns
        -------
        ndarray of shape (5, 5)
        """
        F = np.eye(5)
        kappa = self.curvature(qop)  # 1/mm

        cos_lam = np.cos(lambda_)
        # Deflection in the transverse plane over step_length
        delta_phi = kappa * step_length * cos_lam

        # d(phi) / d(q/p):   dphi/dqop = step * |B| * 0.3e-3 * cos_lam
        F[2, 0] = 0.3e-3 * self.magnitude * step_length * cos_lam

        # d(y) / d(phi) and d(y) / d(q/p) from circular motion
        F[3, 2] = step_length * cos_lam       # dy/dphi
        F[3, 0] = 0.5 * step_length**2 * cos_lam * 0.3e-3 * self.magnitude

        # d(z) / d(lambda)
        F[4, 1] = step_length                 # dz/dlambda  (approx)

        # Phi itself shifts
        F[2, 2] = 1.0  # identity already set; no first-order correction needed
        _ = delta_phi  # suppress unused warning

        return F
