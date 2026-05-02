from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def add_src_to_path() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    src_dir = repo_root / "src"
    if src_dir.is_dir():
        sys.path.insert(0, str(src_dir))

add_src_to_path()

from Core.BField import BField


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sample B-field at x=0, y=0 while stepping z and produce plots."
    )
    parser.add_argument(
        "fieldmap",
        type=str,
        help="Path to the CSV field map file.",
    )
    parser.add_argument(
        "--zmin",
        type=float,
        default=-200.0,
        help="Start z value in mm (default: -200).",
    )
    parser.add_argument(
        "--zmax",
        type=float,
        default=15000.0,
        help="End z value in mm (default: 15000).",
    )
    parser.add_argument(
        "--zstep",
        type=float,
        default=10.0,
        help="z step in mm (default: 10).",
    )
    parser.add_argument(
        "--x",
        type=float,
        default=0.0,
        help="x position in mm (default: 0).",
    )
    parser.add_argument(
        "--y",
        type=float,
        default=0.0,
        help="y position in mm (default: 0).",
    )
    parser.add_argument(
        "--grid-factor",
        type=float,
        default=1.0,
        help="GridFactor passed to BField (default: 1.0).",
    )
    parser.add_argument(
        "--vector-factor",
        type=float,
        default=1.0,
        help="VectorFactor passed to BField (default: 1.0).",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="plots/bfield",
        help="Output directory for plot and CSV (default: plots/bfield).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.zstep <= 0:
        raise ValueError("--zstep must be > 0")

    z_values = np.arange(args.zmin, args.zmax + 0.5 * args.zstep, args.zstep)

    bfield = BField(
        csvfilemap=args.fieldmap,
        GridFactor=args.grid_factor,
        VectorFactor=args.vector_factor,
    )

    bx_values: list[float] = []
    by_values: list[float] = []
    bz_values: list[float] = []

    for z in z_values:
        bx, by, bz = bfield.GetFieldVector(args.x, args.y, float(z))
        bx_values.append(float(bx))
        by_values.append(float(by))
        bz_values.append(float(bz))

    bx_arr = np.asarray(bx_values)
    by_arr = np.asarray(by_values)
    bz_arr = np.asarray(bz_values)
    bmag_arr = np.sqrt(bx_arr * bx_arr + by_arr * by_arr + bz_arr * bz_arr)

    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(
        {
            "x": np.full_like(z_values, args.x, dtype=float),
            "y": np.full_like(z_values, args.y, dtype=float),
            "z": z_values,
            "Bx": bx_arr,
            "By": by_arr,
            "Bz": bz_arr,
            "|B|": bmag_arr,
        }
    )
    csv_path = output_dir / "bfield_scan_x0_y0.csv"
    df.to_csv(csv_path, index=False)

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

    axes[0].plot(z_values, bx_arr, label="Bx", linewidth=1.5)
    axes[0].plot(z_values, by_arr, label="By", linewidth=1.5)
    axes[0].plot(z_values, bz_arr, label="Bz", linewidth=1.5)
    axes[0].set_ylabel("Field component")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].plot(z_values, bmag_arr, color="black", label="|B|", linewidth=1.5)
    axes[1].set_xlabel("z [mm]")
    axes[1].set_ylabel("|B|")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    fig.suptitle(f"B-field scan at x={args.x} mm, y={args.y} mm")
    fig.tight_layout()

    plot_path = output_dir / "bfield_scan_x0_y0.png"
    fig.savefig(plot_path, dpi=160)
    plt.close(fig)

    print(f"Saved samples to: {csv_path}")
    print(f"Saved plot to:    {plot_path}")


if __name__ == "__main__":
    main()