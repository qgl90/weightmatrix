# Roadmap / TODOs

This is a suggested worklist to evolve the current code into the refactor plan.

## Documentation

- [ ] Add per-module “public API” sections (what is safe to call from Apps)
- [ ] Add a glossary for units (mm, MeV/c, Tesla) and branch definitions

## Safety / correctness

- [ ] Audit units in `ParabolicExtrapolator` (the numerical factors assume a specific momentum unit)
- [ ] Consolidate and validate settings path resolution (YAML-relative paths)
- [ ] Add deterministic seeding where random choices are used

## Code structure

- [ ] Extract a `Settings` dataclass with validation
- [ ] Extract `RootTupleWriter` from `MomentumResolutionCalculator`
- [ ] Extract geometry-to-plane conversion into a helper module
- [ ] Split `ComputeResolutions` into smaller private methods (≤100 LOC each)

## Tests

- [ ] Add a “no ROOT” unit test suite for pure helpers
- [ ] Keep the current ROOT-based smoke test for end-to-end regressions (`tests/test_smoke_100_events.py`)

