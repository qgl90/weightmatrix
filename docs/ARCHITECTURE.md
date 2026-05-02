# Architecture overview

This project estimates track momentum resolution by propagating a 5-parameter track state through a sequence of planes:

- **Measurement planes** (with an associated measurement weight matrix `WMeas`)
- **Passive material planes** (contribute via multiple scattering “process noise”)

Planes come from two sources:

1. **CSV detector parameterization** (`settings["detectionlayers"]`)
2. **GDML/ROOT TGeo navigation** (`settings["gdmlInput"]`) for additional passive material (and optionally measurement proxies)

The runner then writes a per-track summary to a ROOT `TNtuple`.

## Packages / modules

- `src/Apps/`
  - `run_sampled_tracks.py`: main runner used to produce ROOT output for later comparisons
  - plotting helpers: `plot_detector_layout.py`, `plot_bfield_scan.py`, `plot_resolution_from_root.py`
- `src/Settings/`
  - `Settings.py`: YAML settings loader (`load_settings`)
- `src/Inputs/`
  - `Loader.py`: ROOT/Parquet loader utilities (not used by current bsmumu runner path)
- `src/Core/`
  - `MomentumResolutionCalculator.py`: pipeline orchestrator (plane construction, propagation, scattering, output writing)
  - `Plane.py`: plane definition (geometry acceptance, measurement weight matrix creation)
  - `TrackState.py`: 5D state `(x, y, tx, ty, q/p)` at a `z`
  - `ParabolicExtrapolator.py`: field-dependent propagation + transport matrix
  - `BField.py`: loads field CSV → ROOT `TH3D` → interpolated `(Bx,By,Bz)`
  - `MaterialLocator.py`: steps ROOT TGeo through volumes and returns layer segments
  - `Propagator.py`, `Matrix.py`, `Precision.py`: algebra + multiple scattering update

## Data flow (high level)

1. **Runner** (`Apps/run_sampled_tracks.py`)
   - loads settings YAML
   - reads parquet events
   - per muon: calls `MomentumResolutionCalculator.ComputeResolutions(...)`
2. **MomentumResolutionCalculator**
   - builds baseline planes from CSV detector description
   - constructs a `BField` and a `ParabolicExtrapolator`
   - constructs a `MaterialLocator` and inserts passive planes from geometry along each propagation step
   - transports covariance, applies measurement updates and multiple-scattering noise
   - fills `TNtuple` branches and writes output ROOT file

## Conventions / units (as used in code)

- Positions are in **mm**
- Track parameters are `tx=dx/dz`, `ty=dy/dz`
- Momentum is treated in **MeV/c** in several places (see `ParabolicExtrapolator.c_light` comment)
- Field map values are interpreted as **Tesla**

## Known hotspots / technical debt (for later refactor)

- Very large `MomentumResolutionCalculator.ComputeResolutions` (multiple responsibilities: I/O, geometry stepping, physics model, output)
- Implicit global configuration (e.g. `_GEOMETRY_FOLDER_` in `Plane.py`)
- Mixed “script style” and “library style”; many modules import heavy dependencies at import time (ROOT, shapely, matplotlib)
- Settings are untyped and path resolution is currently “best effort” (YAML values are used as-is)

