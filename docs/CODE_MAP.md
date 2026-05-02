# Code map (what to read first)

## Main pipeline (momentum resolution)

1. `src/Apps/run_sampled_tracks.py`
2. `src/Core/MomentumResolutionCalculator.py`
3. `src/Core/ParabolicExtrapolator.py`
4. `src/Core/Propagator.py`
5. `src/Core/Plane.py`, `src/Core/TrackState.py`

## Geometry/material navigation

- `src/Core/MaterialLocator.py` (ROOT TGeo stepping)
- `src/Core/MaterialAndMeasurementLocator.py` (alternative locator with active-element patterns; not currently wired into the runner path)

## Inputs / settings

- `src/Settings/Settings.py` (YAML)
- `src/Inputs/Loader.py` (ROOT/Parquet conversion helper; runner currently reads parquet directly)

## Output

- `src/Core/MomentumResolutionCalculator.py` (ROOT `TNtuple` writing)
- `src/Core/RNTupleMaker.py` (unused by default; experimental)

