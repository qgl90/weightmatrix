from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from Settings.Settings import load_settings
from Core.Paths import expand_env_var


def _load_detectionlayers_csvs(files: list[str]) -> pd.DataFrame:
    frames = []
    for file in files:
        frames.append(pd.read_csv(expand_env_var(file), comment="#", skip_blank_lines=True))
    if not frames:
        raise ValueError("No detectionlayers CSV files provided.")
    return pd.concat(frames, ignore_index=True)


def plot_detector_layout(settings_file: str, out_png: str) -> None:
    settings = load_settings(settings_file)
    df = _load_detectionlayers_csvs(settings["detectionlayers"])

    if "Use" in df.columns:
        df_use = df[df["Use"] == 1].copy()
    else:
        df_use = df.copy()

    fig, ax = plt.subplots(figsize=(11, 4.5))

    # Passive/material planes (all)
    ax.scatter(df["Z"], df.get("thickness", 0.0) * 100.0, s=12, alpha=0.35, label="All planes (x/X0 %)")

    # Measurement planes (subset where sigmaX/sigmaY are present)
    is_meas = (~df_use.get("sigmaX").isna()) | (~df_use.get("sigmaY").isna())
    meas = df_use[is_meas]
    if not meas.empty:
        ax.scatter(
            meas["Z"],
            meas.get("thickness", 0.0) * 100.0,
            s=26,
            marker="s",
            label="Measurement planes (from CSV)",
        )

    ax.set_xlabel("z [mm]")
    ax.set_ylabel("x/X0 [%] (per plane)")
    ax.set_title(f"Detector layout from {Path(settings_file).name}")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")

    out_path = Path(out_png)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)


def main() -> None:
    ap = ArgumentParser(description="Plot detector layout from detectionlayers CSV(s) referenced by a YAML settings file.")
    ap.add_argument("--settings_file", default="layouts/u2_debug_material.yaml")
    ap.add_argument("--out_png", default="outputs/detector_layout.png")
    args = ap.parse_args()
    plot_detector_layout(args.settings_file, args.out_png)


if __name__ == "__main__":
    main()

