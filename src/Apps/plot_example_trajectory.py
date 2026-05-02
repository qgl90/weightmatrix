from __future__ import annotations

from pathlib import Path

import matplotlib

# Headless-safe backend (CI / batch jobs)
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from ROOT import TLorentzVector

from Core.BField import BField
from Core.ParabolicExtrapolator import ParabolicExtrapolator
from Core.TrackState import TrackState


def plot_example_trajectory(
    *,
    settings: dict,
    fourvector: TLorentzVector,
    charge: int,
    origin_position: tuple[float, float, float],
    out_png: str,
    z_values_mm: list[float] | None = None,
) -> None:
    """
    Plot an example propagated trajectory (x(z), y(z), and x-y projection).

    This uses the same propagator as the resolution calculation:
    - field map: settings["fieldMap"] (+ scale factors)
    - propagator: `Core.ParabolicExtrapolator.ParabolicExtrapolator`

    Parameters
    ----------
    settings:
        Resolved settings dict.
    fourvector:
        ROOT TLorentzVector for the particle.
    charge:
        Particle charge (+1/-1).
    origin_position:
        (x,y,z) in mm.
    out_png:
        Output PNG path.
    z_values_mm:
        Optional list of z locations to sample. If omitted, uses a default coarse list.
    """
    field = BField(settings["fieldMap"], settings["fieldShrinkFactor"], settings["fieldMagFactor"])
    prop = ParabolicExtrapolator(field)

    unit = fourvector.Vect().Unit()
    tx = float(unit.X() / unit.Z())
    ty = float(unit.Y() / unit.Z())
    qop = float(charge) / float(fourvector.P())

    state0 = TrackState(
        x=float(origin_position[0]),
        y=float(origin_position[1]),
        z=float(origin_position[2]),
        tx=tx,
        ty=ty,
        qop=qop,
    )

    if z_values_mm is None:
        z0 = state0.Z()
        z_values_mm = [z for z in [z0 + 100, 900, 1500, 2500, 4000, 6000, 8000, 10000, 12000, 15000, 18000] if z > z0]

    states, curling = prop.CreateStatesAlongTrajectory_ZVALUES(state0, z_values_mm)

    z = np.asarray([s.Z() for s in states], dtype=float)
    x = np.asarray([s.X() for s in states], dtype=float)
    y = np.asarray([s.Y() for s in states], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    axes[0].plot(z, x, marker="o", ms=3)
    axes[0].set_xlabel("z [mm]")
    axes[0].set_ylabel("x [mm]")
    axes[0].set_title("Example trajectory: x(z)")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(z, y, marker="o", ms=3, color="tab:orange")
    axes[1].set_xlabel("z [mm]")
    axes[1].set_ylabel("y [mm]")
    axes[1].set_title("Example trajectory: y(z)")
    axes[1].grid(True, alpha=0.25)

    axes[2].plot(x, y, marker="o", ms=3, color="tab:green")
    axes[2].set_xlabel("x [mm]")
    axes[2].set_ylabel("y [mm]")
    axes[2].set_title("Example trajectory: x-y")
    axes[2].grid(True, alpha=0.25)
    axes[2].axis("equal")

    title = f"p={fourvector.P()/1000.0:.1f} GeV, eta={fourvector.Eta():.2f}, q={charge}"
    if curling:
        title += " (curling cut)"
    fig.suptitle(title)

    out_path = Path(out_png)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)

