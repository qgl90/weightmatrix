from __future__ import annotations
from unicodedata import name
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap

def preliminary() : 
    import sys
    from pathlib import Path
    # def add_src_to_path() -> None:
    """
    Allow running scripts directly (e.g. `python3 scripts/run_it.py`) without
    requiring an editable install by prepending `<repo>/src` to `sys.path`.
    """
    repo_root = Path(__file__).resolve().parents[2]    
    src_dir = repo_root / "src"
    if src_dir.is_dir():
        print(f"Append to sys.path {src_dir}")
        sys.path.insert(0, str(src_dir))
preliminary()


from alive_progress import alive_bar
import ROOT as r
from ROOT import TLorentzVector
r.DisableImplicitMT() 
from Utils.Logger import Logger
from Inputs.Loader import LoadInfoTracks
from Settings.Settings import load_settings
from Core.MomentumResolutionCalculator import MomentumResolutionCalculator
from Core.Plane import Plane, PicklePlane
from Utils.Logger import Logger
import mplhep as hep 
hep.style.use("LHCb2")
if __name__ == "__main__":
    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(10*4, 8), sharex=True, sharey=True)    
    for ax, plane_setup in zip( axs.flatten(), [
        "Blake_Pixel_x",
        "Blake_Fibre_x",
        "Blake_Fibre_u",
        "Blake_Fibre_v"
        ]):         
        plane = Plane(zval=0, name=plane_setup, geomLayer = plane_setup, thickness = 0.01, dzMat = 1.0, sigmaX = 0.01, sigmaY = 0.01, angle=0)
        n_points = 500
        x_values = np.linspace(-4000, 4000, n_points)
        y_values = np.linspace(-4000, 4000, n_points)

        # acceptance_map[y_index, x_index] = 1 for in-acceptance, 0 for out-of-acceptance
        acceptance_map = np.zeros((n_points, n_points), dtype=np.uint8)

        with alive_bar(n_points, title="Evaluating plane acceptance") as bar:
            for ix, x in enumerate(x_values):
                for iy, y in enumerate(y_values):
                    if plane.PointInAcceptance(float(x), float(y)):
                        acceptance_map[iy, ix] = 1
                bar()

        # red = in acceptance, blue = out of acceptance
        cmap = ListedColormap(["royalblue", "red"])

        ax.imshow(
            acceptance_map,
            origin="lower",
            extent=[x_values[0], x_values[-1], y_values[0], y_values[-1]],
            cmap=cmap,
            alpha = 0.5,
            interpolation="nearest",
            aspect="equal",
        )
        ax.set_xlim( -4000,4000)
        ax.set_ylim( -4000,4000)
        ax.set_title(f"Acceptance map: {plane.name}")
        ax.set_xlabel("x [mm]", loc='center')
        ax.set_ylabel("y [mm]", loc='center')

        legend_handles = [
            Patch(facecolor="red", edgecolor="black", label="In acceptance"),
            Patch(facecolor="royalblue", edgecolor="black", label="Out of acceptance"),
        ]
        ax.legend(handles=legend_handles, loc="upper right")
    Logger.info(f"Saving acceptance map figure to tests/AcceptanceMap_Blake_Layout.png")
    fig.savefig(f"slides/plots/AcceptanceMap_Blake_Layout.png")
    