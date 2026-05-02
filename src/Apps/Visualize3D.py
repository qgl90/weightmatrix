from __future__ import annotations
from unicodedata import name
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from mpl_toolkits.mplot3d import art3d

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
import pandas as pd
from Core.Paths import expand_env_var

def plot_design(design):
    """
    Plot the design of a given tracker downstream layout used
    """    
    Logger.info(f"Processing design: {design}")
    settings  = load_settings(f"layouts/{design}.yaml")

    df                = pd.concat((pd.read_csv(expand_env_var(f), comment ="#", skip_blank_lines=True) for f in settings["detectionlayers"]), ignore_index=True)
    dfUse             = df[ df['Use'] ==1]
    planes            = []
    # MyField           = BField( settings["fieldMap"],     settings["fieldShrinkFactor"] , settings["fieldMagFactor"])        
    # TrajectoryCreator = ParabolicExtrapolator(MyField)
    """
    Add custom Planes Here for the tracking, with their material budget and resolution if they are measurement planes.
    The planes are sorted in z and then used for the tracking and resolution estimation.        
    """
    
    for  z, thick,  sigmaX ,sigmaY, angle, geomLayer, hiteff in zip( dfUse['Z'],
                                                                        dfUse['thickness'],
                                                                        dfUse['sigmaX'],
                                                                        dfUse['sigmaY'],
                                                                        dfUse['angle'] ,
                                                                        dfUse['xyshape'].replace(" ",""),
                                                                        dfUse['hiteff']) :
        planes.append( Plane(
            zval     = float(z),
            thickness= float(thick), dzMat= 0.0,
            sigmaX   = float(sigmaX) if not( pd.isnull( sigmaX)) else None ,
            sigmaY   = float(sigmaY) if not( pd.isnull( sigmaY)) else None ,
            angle    = float(angle)  if not( pd.isnull( sigmaX)) else None ,
            geomLayer= geomLayer,
            hitefficiency = hiteff
    ))
    planes = sorted(planes)

    fig = plt.figure(figsize=(15,15))
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.1)
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    ax_table = fig.add_subplot(gs[1, 0])
    ax_table.axis("off")

    color_measurement   = "steelblue"
    color_passive       = "tomato"

    def _ring_vertices_codes(ring_coords):
        ring_arr = np.asarray(ring_coords, dtype=float)
        if ring_arr.shape[0] < 3:
            return [], []
        if not np.allclose(ring_arr[0], ring_arr[-1]):
            ring_arr = np.vstack([ring_arr, ring_arr[0]])
        verts = ring_arr.tolist()
        codes = [Path.MOVETO] + [Path.LINETO] * (len(verts) - 2) + [Path.CLOSEPOLY]
        return verts, codes

    def _polygon_path(poly):
        verts_all = []
        codes_all = []
        v_ext, c_ext = _ring_vertices_codes(poly.exterior.coords)
        verts_all.extend(v_ext)
        codes_all.extend(c_ext)
        for interior in poly.interiors:
            v_int, c_int = _ring_vertices_codes(interior.coords)
            verts_all.extend(v_int)
            codes_all.extend(c_int)
        return Path(np.asarray(verts_all, dtype=float), codes_all)

    def _add_geometry_patch_3d(ax3d, geom, z_value, color_value, alpha_value):
        if geom is None:
            return False
        geom_type = getattr(geom, "geom_type", "")

        if geom_type == "Polygon":
            patch = PathPatch(_polygon_path(geom), facecolor=color_value,
                            edgecolor=color_value, lw=0.4, alpha=alpha_value)
            ax3d.add_patch(patch)
            # Axis remap for 3D view: old Z -> new X, old X -> new Y, old Y -> new Z.
            art3d.pathpatch_2d_to_3d(patch, z=z_value, zdir="x")
            return True

        if geom_type == "MultiPolygon":
            drew = False
            for poly in geom.geoms:
                patch = PathPatch(_polygon_path(poly), facecolor=color_value,
                                edgecolor=color_value, lw=0.4, alpha=alpha_value)
                ax3d.add_patch(patch)
                art3d.pathpatch_2d_to_3d(patch, z=z_value, zdir="x")
                drew = True
            return drew

        return False

    def _geom_points_for_bounds(geom):
        if geom is None:
            return []
        geom_type = getattr(geom, "geom_type", "")
        if geom_type == "Polygon":
            return [np.asarray(geom.exterior.coords, dtype=float)]
        if geom_type == "MultiPolygon":
            return [np.asarray(poly.exterior.coords, dtype=float) for poly in geom.geoms]
        return []

    table_rows = []
    projection_rows = []

    for p in planes:
        if p.geomLayer == "any":
            continue

        zLocation = p.Z
        thickness = p.Thickness
        xy_geom   = p.GetPolygon()
        alpha_face = 0.3 if p.IsMeasurement else 0.05
        if p.IsMeasurement and p.IsPixelMeasurement():
            color = "purple"
            detail = f"errX={p.sigmaX*1000:.1f} um, errY={p.sigmaY*1000:.1f} um"
            
        elif p.IsMeasurement:
            color = "steelblue"
            detail = f"errX={p.sigmaX*1000:.1f} um, stereo={p.angle} deg"
        else: 
            color = 'red'
            detail = "passive"
        
        if xy_geom is None:
            continue

        # Visual slab depth in mm; thickness is in X0 so this is a display scaling.
        dz = max(thickness * 10.0, 0.1)
        drew_front = _add_geometry_patch_3d(ax, xy_geom, zLocation, color, alpha_face)
        # _add_geometry_patch_3d(ax, xy_geom, zLocation + dz, color, alpha_face * 0.55)
        if not drew_front:
            continue

        bound_parts = _geom_points_for_bounds(xy_geom)
        if not bound_parts:
            continue
        bound_xy = np.vstack(bound_parts)
        x_min, x_max = float(np.min(bound_xy[:, 0])), float(np.max(bound_xy[:, 0]))
        y_min, y_max = float(np.min(bound_xy[:, 1])), float(np.max(bound_xy[:, 1]))
        projection_rows.append((x_min, x_max, y_min, y_max, zLocation, zLocation + dz, color, thickness * 100.0))

        # ---- thickness annotation at centroid ----
        if hasattr(xy_geom, "centroid"):
            cx, cy = float(xy_geom.centroid.x), float(xy_geom.centroid.y)
        else:
            cx, cy = float(np.mean(bound_xy[:, 0])), float(np.mean(bound_xy[:, 1]))
        # ax.text(zLocation, cx, cy, f"x/x0={thickness*100:.2f} %", fontsize=5,
        #         color=color, ha="center", va="bottom", alpha=0.9)

        table_rows.append([
            p.geomLayer,
            "pixel" if p.IsMeasurement and p.IsPixelMeasurement() else ("strip" if p.IsMeasurement else "passive"),
            f"{zLocation:.1f}",
            f"{thickness*100:.2f}",
            detail,
        ])

        Logger.info(p)

    # Remove pane fills and borders so only the drawn polygons are visible
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor("none")
    ax.grid(False)
    z_view_min, z_view_max = 6500.0, 10000.0
    ax.set_xlim(z_view_min, z_view_max)
    ax.margins(x=0)
    ax.set_xlabel("Z [mm]", fontsize=9)
    ax.set_ylabel("X [mm]", fontsize=9)
    ax.set_zlabel("Y [mm]", fontsize=9)
    ax.set_title("3D MightyTracker Layout", fontsize=11)
    ax.tick_params(axis="x", labelsize=7)
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(axis="z", labelsize=7)
    ax.set_box_aspect((4.0, 1.0, 1.0))
    ax.view_init(elev=30, azim=-120)

    # Fix axes limits from the planes' bounding boxes
    all_xy_parts = []
    all_z = []
    for p in planes:
        if p.geomLayer == "any":
            continue
        geom = p.GetPolygon()
        parts = _geom_points_for_bounds(geom)
        if not parts:
            continue
        all_xy_parts.extend(parts)
        all_z.append(p.Z)

    if all_xy_parts and all_z:
        all_xy = np.vstack(all_xy_parts)
        ax.set_ylim(all_xy[:, 0].min(), all_xy[:, 0].max())
        ax.set_zlim(all_xy[:, 1].min(), all_xy[:, 1].max())

    # Build an informative table on the right side of the viewer.
    if table_rows:
        max_rows = 30
        shown_rows = table_rows[:max_rows]
        if len(table_rows) > max_rows:
            shown_rows.append(["...", "...", "...", "...", f"+{len(table_rows) - max_rows} more rows"])

        table = ax_table.table(
            cellText=shown_rows,
            colLabels=["Layer", "Type", "Z [mm]", "x/x0 [%]", "Details"],
            loc="center",
            cellLoc="left",
            colLoc="left",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.0, 1.00)

        # Style header and row backgrounds for readability.
        n_cols = 5
        for c in range(n_cols):
            header_cell = table[(0, c)]
            header_cell.set_facecolor("#2F2F2F")
            header_cell.get_text().set_color("white")

        for r_idx, row in enumerate(shown_rows, start=1):
            row_type = row[1]
            if row_type == "pixel":
                bg = "#E9D8FD"
            elif row_type == "strip":
                bg = "#DBEAFE"
            elif row_type == "passive":
                bg = "#FEE2E2"
            else:
                bg = "#F3F4F6"
            for c in range(n_cols):
                table[(r_idx, c)].set_facecolor(bg)
    fig.subplots_adjust(top=0.98, bottom=0.04, left=0.03, right=0.99, hspace=0.02)
    # Force a wide top panel so the remapped Z axis spans almost all figure width.
    ax.set_position([0.03, 0.26, 0.96, 0.70])
    ax_table.set_position([0.03, 0.03, 0.96, 0.19])
    Logger.info(f"Saving file to slides/plots/geometry_3d_{design}.pdf")
    fig.savefig(f"slides/plots/geometry_3d_{design}.pdf", bbox_inches='tight')
    # XZ and YZ projection view in a dedicated figure.
    fig_proj, (ax_xz, ax_yz) = plt.subplots(1, 2, figsize=(30, 12))
    xz_labels = []
    yz_labels = []
    for x_min, x_max, y_min, y_max, z0, z1, color, x_over_x0_pct in projection_rows:
        z_span = max(z1 - z0, 0.1)
        x_span = max(x_max - x_min, 0.1)
        y_span = max(y_max - y_min, 0.1)

        ax_xz.add_patch(Rectangle((z0, x_min), z_span, x_span,
                                facecolor=color, edgecolor=color, alpha=1.0, lw=0.4))
        ax_yz.add_patch(Rectangle((z0, y_min), z_span, y_span,
                                facecolor=color, edgecolor=color, alpha=1.0, lw=0.4))

        text_label = f"x/x0={x_over_x0_pct:.1f}%"
        xz_labels.append((z0 + 0.5 * z_span, x_min + 0.5 * x_span, text_label))
        yz_labels.append((z0 + 0.5 * z_span, y_min + 0.5 * y_span, text_label))

    ax_xz.set_title("Projection on XZ") #, fontsize=11)
    ax_yz.set_title("Projection on YZ") #, fontsize=11)
    ax_xz.set_xlabel("Z [mm]") #, fontsize=10)
    ax_yz.set_xlabel("Z [mm]") #, fontsize=10)
    ax_xz.set_ylabel("X [mm]") #, fontsize=10)
    ax_yz.set_ylabel("Y [mm]") #, fontsize=10)

    for ax2d in (ax_xz, ax_yz):
        ax2d.grid(True, alpha=0.25)
        # ax2d.tick_params(labelsize=8)
        ax2d.set_xlim(z_view_min, z_view_max)

    ax_xz.autoscale_view(scalex=False, scaley=True)
    ax_yz.autoscale_view(scalex=False, scaley=True)

    
    # Rotate labels and stagger them from low to high vertical coordinates to reduce overlaps.
    def _add_staggered_labels(axis, labels):
        if not labels:
            return
        y_lo, y_hi = axis.get_ylim()
        y_span = max(abs(y_hi - y_lo), 1.0)
        step = 0.018 * y_span
        sorted_labels = sorted(labels, key=lambda t: t[1])
        n = len(sorted_labels)
        for idx, (x_c, y_c, txt) in enumerate(sorted_labels):
            y_shift = (idx - 0.5 * (n - 1)) * step
            # axis.text(x_c, y_c + y_shift, txt,
            #           fontsize=6, color="black", ha="center", va="center",
            #           rotation=90, clip_on=True)

    _add_staggered_labels(ax_xz, xz_labels)
    _add_staggered_labels(ax_yz, yz_labels)
    fig_proj.tight_layout()
    Logger.info(f"Saving file to slides/plots/geometry_projections_{design}.pdf")
    fig_proj.savefig(f"slides/plots/geometry_projections_{design}.pdf", bbox_inches='tight')
    
if __name__ == "__main__":
    Logger.info("Running script to visualize 3D trajectories")        
    for design in "u2_baseline", "u2_uli":
        plot_design(design)
    