"""
Material layer representation and multiple-scattering / energy-loss contributions.

Uses the Highland formula for multiple Coulomb scattering and the Bethe-Bloch
formula for mean energy loss to compute the contributions to track resolution
from a material layer (e.g. a GDML-described detector volume).
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass


# Physical constants
ELECTRON_MASS_GEV = 0.000511   # GeV
PROTON_MASS_GEV = 0.938272     # GeV (default particle mass)
C_LIGHT = 299.792458            # mm/ns


@dataclass
class MaterialLayer:
    """
    A material layer described by its radiation length X0 and thickness.

    Parameters
    ----------
    thickness : float
        Physical thickness of the layer [mm].
    x0 : float
        Radiation length of the material [mm].
        Typical values: Si ~ 93.7 mm, Be ~ 353 mm, Air ~ 304_000 mm.
    name : str
        Human-readable identifier (e.g. GDML volume name).
    """

    thickness: float
    x0: float
    name: str = "layer"

    @property
    def x_over_x0(self) -> float:
        """Material budget t/X0 (dimensionless)."""
        return self.thickness / self.x0

    def highland_theta0(self, p_gev: float, q: float = 1.0, beta: float = 1.0) -> float:
        """
        RMS scattering angle from the Highland formula [rad].

        theta_0 = (13.6 MeV / (beta * c * p)) * |z| * sqrt(t/X0)
                  * [1 + 0.038 * ln(t/X0)]

        Parameters
        ----------
        p_gev : float
            Particle momentum [GeV].
        q : float
            Particle charge number (default 1).
        beta : float
            Relativistic beta = v/c (default 1 = ultra-relativistic).

        Returns
        -------
        float
            RMS projected scattering angle theta_0 [rad].
        """
        t = self.x_over_x0
        if t <= 0 or p_gev <= 0:
            return 0.0
        theta0 = (
            (13.6e-3 / (beta * p_gev))
            * abs(q)
            * np.sqrt(t)
            * (1.0 + 0.038 * np.log(t))
        )
        return float(theta0)

    def scattering_weight_correction(
        self, p_gev: float, q: float = 1.0, beta: float = 1.0
    ) -> np.ndarray:
        """
        Return the 5×5 additive weight-matrix correction from multiple scattering.

        In the weight-matrix formalism a thin scatterer adds a noise term Q
        to the covariance matrix.  The weight matrix contribution is subtracted:

            W_new = (C + Q)^{-1}

        This function returns the noise covariance Q for the scattering angles
        (lambda and phi directions), assuming a thin-scatterer approximation.

        Parameters
        ----------
        p_gev : float
            Momentum [GeV].
        q : float
            Charge number.
        beta : float
            Relativistic beta.

        Returns
        -------
        ndarray of shape (5, 5)
            Noise covariance Q to be *added* to the track covariance matrix.
        """
        theta0 = self.highland_theta0(p_gev, q=q, beta=beta)
        sigma2 = theta0**2

        Q = np.zeros((5, 5))
        # Scattering affects the angular parameters lambda (idx 1) and phi (idx 2)
        Q[1, 1] = sigma2   # lambda scatter
        Q[2, 2] = sigma2   # phi scatter
        return Q

    def __repr__(self) -> str:
        return (
            f"MaterialLayer(name={self.name!r}, "
            f"thickness={self.thickness:.3g} mm, "
            f"X0={self.x0:.3g} mm, "
            f"t/X0={self.x_over_x0:.4g})"
        )
