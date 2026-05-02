from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ROOT


def _read_branch(tree, name: str) -> np.ndarray:
    vals = []
    for entry in tree:
        vals.append(getattr(entry, name))
    return np.asarray(vals, dtype=float)


def plot_resolution(root_file: str, out_png: str, *, pres_branch: str = "fwd_pres_LastMeasurement") -> None:
    f = ROOT.TFile.Open(root_file, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {root_file}")

    tree = f.Get("ntuple")
    if not tree:
        raise RuntimeError("Cannot find TNtuple named 'ntuple' in file.")

    p = _read_branch(tree, "p")
    pres = _read_branch(tree, pres_branch)

    # pres is typically a resolution-like observable; keep only sane values
    mask = np.isfinite(p) & np.isfinite(pres) & (p > 0) & (pres > 0)
    p = p[mask]
    pres = pres[mask]

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.hist2d(p / 1000.0, pres, bins=(80, 80), cmap="viridis")
    ax.set_xlabel("p [GeV]")
    ax.set_ylabel(pres_branch)
    ax.set_title("Momentum resolution observable vs p")
    fig.colorbar(ax.collections[0], ax=ax, label="counts")

    out_path = Path(out_png)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)


def main() -> None:
    ap = ArgumentParser(description="Plot a quick 2D histogram of momentum resolution observable vs momentum from the output TNtuple.")
    ap.add_argument("--root_file", required=True, help="ROOT file produced by run_sampled_tracks.py")
    ap.add_argument("--out_png", default="outputs/momentum_resolution_vs_p.png")
    ap.add_argument("--pres_branch", default="fwd_pres_LastMeasurement")
    args = ap.parse_args()
    plot_resolution(args.root_file, args.out_png, pres_branch=args.pres_branch)


if __name__ == "__main__":
    main()

