from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from Settings.Settings import load_settings
from Core.BField import BField


def plot_bfield_z_scan(settings_file: str, out_png: str, *, x_mm: float, y_mm: float, n: int) -> None:
    settings = load_settings(settings_file)
    field = BField(settings["fieldMap"], settings["fieldShrinkFactor"], settings["fieldMagFactor"])

    zmin, zmax = field.z_min_max
    z = np.linspace(zmin, zmax, n)
    bx = np.empty_like(z)
    by = np.empty_like(z)
    bz = np.empty_like(z)
    for i, zi in enumerate(z):
        bx[i], by[i], bz[i] = field.GetFieldVector(x_mm, y_mm, float(zi))

    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.plot(z, bx, label="Bx [T]")
    ax.plot(z, by, label="By [T]")
    ax.plot(z, bz, label="Bz [T]")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("B [T]")
    ax.set_title(f"B-field z-scan at x={x_mm:.1f} mm, y={y_mm:.1f} mm")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")

    out_path = Path(out_png)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)


def main() -> None:
    ap = ArgumentParser(description="Quick B-field sanity check: scan Bx/By/Bz vs z at fixed (x,y).")
    ap.add_argument("--settings_file", default="layouts/u2_debug_material.yaml")
    ap.add_argument("--out_png", default="outputs/bfield_z_scan.png")
    ap.add_argument("--x_mm", type=float, default=0.0)
    ap.add_argument("--y_mm", type=float, default=0.0)
    ap.add_argument("--n", type=int, default=400)
    args = ap.parse_args()
    plot_bfield_z_scan(args.settings_file, args.out_png, x_mm=args.x_mm, y_mm=args.y_mm, n=args.n)


if __name__ == "__main__":
    main()

