# Refactor proposal (structure + new classes)

This document proposes a restructuring that keeps physics logic intact but makes
the codebase easier to understand, test, and evolve.

The current code works as a *single pipeline script* centered on
`Core/MomentumResolutionCalculator.py`. The main goal of the refactor is to
separate concerns:

- configuration / settings
- inputs (kinematics)
- detector model (CSV-defined planes)
- geometry/material sampling (TGeo)
- field + propagation (transport matrices)
- filtering / scattering model
- output writing

## Proposed package layout

Create a real Python package (keeping `src/` layout), e.g.:

```
src/weightmatrix/
  __init__.py
  config.py
  units.py
  field/
    __init__.py
    csv_fieldmap.py
  geometry/
    __init__.py
    tgeo_navigator.py
  detector/
    __init__.py
    layers.py
    acceptance.py
  tracking/
    __init__.py
    state.py
    plane.py
    propagators.py
    filter.py
    resolution.py
  io/
    __init__.py
    parquet_kinematics.py
    root_output.py
  apps/
    run_sampled_tracks.py
    plot_*.py
```

This can be introduced incrementally by re-exporting the old module names to
avoid breaking existing scripts.

## Core domain objects (new classes)

### 1) `Settings` (typed config)

Current: settings is a plain dict loaded from YAML.

Proposed:

- `@dataclass Settings` with typed fields
- `Settings.from_yaml(path)` that:
  - resolves relative paths relative to YAML location
  - validates required keys
  - normalizes booleans and lists

### 2) `FieldMap` interface

Current: `Core/BField.py` is tightly coupled to ROOT histograms and CSV.

Proposed:

- `class FieldMap(Protocol): Get(x,y,z)->np.ndarray`
- `CsvRootFieldMap(FieldMap)`: current implementation
- (optional) `ConstantFieldMap(FieldMap)` for testing

Benefits:
- unit tests for propagation without requiring a full field CSV

### 3) `GeometryNavigator` interface

Current: `MaterialLocator` returns low-level arrays and `MomentumResolutionCalculator`
converts them into planes inline.

Proposed:

- `GeometryNavigator.sample_segments(start_state, end_state) -> list[MaterialSegment]`
- `MaterialSegment` dataclass:
  - `entry_z`, `exit_z`, `radlen_mm`, `node_id`
  - computed `dz_mm`, `x_over_x0`

Optional extension:
- `GeometryNavigator.sample_measurements(...)` for active sensors in GDML

### 4) `DetectorModel` (CSV-defined plane set)

Current: `MomentumResolutionCalculator` reads CSVs and instantiates `Plane` objects directly.

Proposed:

- `DetectorModel.from_csvs(list[path]) -> DetectorModel`
- Stores:
  - list of `PlaneDefinition` (z, x/X0, measurement model, acceptance tag, efficiency)
  - utilities for plane subsets (e.g. “measurement-only”)

### 5) `Trajectory` / `Propagator`

Current: `ParabolicExtrapolator` mixes propagation and state list creation.

Proposed:

- `Propagator.propagate(state, z_new) -> (state_new, transport_matrix)`
- `TrajectorySampler.sample(state0, z_values) -> list[state]`

### 6) `ResolutionEngine` (the filter)

Current: `ComputeResolutions` is a very large method with multiple nested loops.

Proposed main abstraction:

- `ResolutionEngine.compute(track: TrackInput) -> ResolutionResult`

Where:

- `TrackInput` dataclass:
  - four-vector
  - charge
  - origin position, PV position
- `ResolutionResult` dataclass:
  - per-plane states (requested storing planes)
  - summary counters
  - final resolution metrics (the current `pres` etc.)

### 7) `OutputWriter`

Current: ROOT `TNtuple` definition/branch building is embedded in `MomentumResolutionCalculator`.

Proposed:

- `RootTupleWriter(schema).fill(result).close()`
- schema is built once (plane names known up-front)

Benefits:
- algorithm code becomes testable without ROOT (by swapping writer)

## Phased migration plan (safe + incremental)

1. Introduce typed settings + path resolution, keep old dict API for compatibility.
2. Extract small pure helpers from `ComputeResolutions`:
   - plane list construction
   - conversion of geometry segments → passive planes
   - branch filling helpers
3. Wrap output writing into `RootTupleWriter`.
4. Create `ResolutionEngine` and keep `MomentumResolutionCalculator` as a thin wrapper.
5. Move plotting utilities under `apps/` and keep stable CLIs.

## Immediate “clean up” targets (low risk)

- Remove unused imports and duplicated imports (many modules import the same things multiple times).
- Make all imports package-relative (`from Core...` not `from src.Core...`).
- Move heavyweight imports (ROOT, shapely, matplotlib) into function scope where possible so importing the package is cheap.
- Add type hints for public functions and dataclasses for structured records.
- Add a small set of unit tests:
  - transport matrix shape/invertibility in trivial cases
  - field map interpolation boundary behavior
  - geometry segment sampling returns monotonic z segments

