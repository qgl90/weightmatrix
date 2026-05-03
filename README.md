# Weightmatrix estimator of resolutions along tracks

The code uses the weight formalism to evaluate what the kalman filter for tracks would do in the transport of tracks in a generic 3D magnetic field, in presence of measurements with specific position resolutions and material budget. It also navigates existing materials from full geometries to identify for a projectile tracks the amount of crossed material and introduce noise terms to resolutions. The tracks states are expressed in (x,y,tx,ty,q/p) as in LHCb experiment for a forward tracker layout. Planes and states are defined at specific z locations. Custom acceptance masks and layouts can be parametrically injected with csv tables inputs.

## What it does

- Loads a **B-field map** (CSV grid) and interpolates `Bx/By/Bz`.
- Loads a **detector layout** from one or more CSV files (planes with material + optional measurement resolutions).
- Optionally navigates a **GDML / ROOT TGeo** geometry to add passive material (and some measurement proxies) along track steps.
- Runs a parabolic propagation and multiple-scattering model to estimate **momentum resolution** and writes results to a ROOT `TNtuple`.

## Quick start (local)

1) Create an environment that has **PyROOT/ROOT** available.
2) Install Python deps:

```bash
python3 -m pip install -r requirements.txt
```

3) Run:

```bash
PYTHONPATH=$PWD/src:$PYTHONPATH python3 src/Apps/run_sampled_tracks.py \
  --kinematic_file samples/bsmumu_rdf.parquet \
  --settings_file layouts/u2_debug_material.yaml \
  --output_file outputs/debug_baseline_bsmumu.root \
  --start_entry 0 --stop_entry 100
```

## Plot helpers

All write PNGs under `outputs/` by default:

```bash
PYTHONPATH=$PWD/src:$PYTHONPATH python3 src/Apps/plot_detector_layout.py --settings_file layouts/u2_debug_material.yaml
PYTHONPATH=$PWD/src:$PYTHONPATH python3 src/Apps/plot_bfield_scan.py --settings_file layouts/u2_debug_material.yaml
PYTHONPATH=$PWD/src:$PYTHONPATH python3 src/Apps/plot_resolution_from_root.py --root_file outputs/debug_baseline_bsmumu_0_100.root
```

## Notebook

- `notebooks/momentum_resolution_quicklook.ipynb`: runs 100 events and makes a quick resolution-vs-p plot.

## Notes

- Outputs are written under `outputs/` (gitignored).
- A concise pipeline description is in `docs/momentum_resolution_pipeline.md`.
