# Output writer refactor: steps + extension points

This project writes per-track outputs into a ROOT `TNtuple` called `ntuple`.

Historically the output writing was mixed into `Core/MomentumResolutionCalculator.py`:

- schema creation (`branches=":".join(branch_list)`)
- ROOT file/tuple creation
- per-track filling (`TNtuple.Fill(array('f', ...))`)
- final write/close

## New component: `RootTupleWriter`

Implemented in:

- `src/weightmatrix/io/root_output.py`

### Step 1: Build a schema

Schema is pure data and defines a stable branch order:

```python
from weightmatrix.io.root_output import TupleSchema
schema = TupleSchema(branches=(...,))
```

### Step 2: Open writer

```python
from weightmatrix.io.root_output import RootTupleWriter
writer = RootTupleWriter.open("outputs/out.root", schema)
```

### Step 3: Fill rows

The writer intentionally fills *rows* (ordered float sequences) to keep the
core algorithm free from ROOT details:

```python
writer.fill_row([1.0, 2.0, ...])
```

### Step 4: Close

```python
writer.close()
```

## Why this helps

- Output branch schema is now explicit and reusable.
- Adding new properties becomes “schema + fill one more column”, without touching
  ROOT file/tuple creation code.
- Future upgrades (switch to `RNTuple`, write multiple trees, add metadata
  objects) can be added behind this interface.

## Planned next improvements (optional)

- Add a `RowBuilder` helper that fills rows from dictionaries (with defaults).
- Add a `PlaneBlock` helper to generate per-plane branch blocks and keep naming
  conventions in one place.
- Add a `RootMetadataWriter` to store settings YAML, git hash, and schema into
  the output ROOT file.

