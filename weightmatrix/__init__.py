"""
weightmatrix: Weight matrix formalism for particle detector track resolution prediction.

Uses weight matrix (inverse covariance) formalism combined with B-field (xyz components)
and GDML material input to parametrically predict track-state resolution in a detector.
"""

from .matrix import WeightMatrix, TrackState, propagate_weight_matrix
from .bfield import BField
from .material import MaterialLayer

__all__ = [
    "WeightMatrix",
    "TrackState",
    "propagate_weight_matrix",
    "BField",
    "MaterialLayer",
]
