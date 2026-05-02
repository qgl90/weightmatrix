from __future__ import annotations

"""
ROOT output writer (modular TNtuple writing).

This is a small abstraction layer over ROOT `TFile` + `TNtuple` to keep output
schema definition and filling logic out of the physics pipeline.

Design goals
------------
- Define schema once (branch order is explicit and stable)
- Fill rows by passing a list of floats in schema order
- Allow incremental extension:
  - add extra scalar branches
  - add per-plane blocks
  - add new blocks without rewriting the calculator
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import array


@dataclass(frozen=True)
class TupleSchema:
    """
    Defines the column order for the output TNtuple.

    Keep this as pure data so schema building is testable without ROOT.
    """

    branches: tuple[str, ...]

    def as_root_branch_string(self) -> str:
        return ":".join(self.branches)


@dataclass
class RootTupleWriter:
    """
    Writes a ROOT `TNtuple` named `ntuple` to a ROOT file.

    Typical usage:

        schema = TupleSchema(branches=(...,))
        writer = RootTupleWriter.open("out.root", schema)
        writer.fill_row([1.0, 2.0, ...])
        writer.close()
    """

    # ROOT types are intentionally kept as `object` to allow importing this
    # module in environments where PyROOT is not installed.
    file: object
    tuple: object
    schema: TupleSchema
    _closed: bool = field(default=False, init=False)

    @classmethod
    def open(cls, out_file: str | Path, schema: TupleSchema, *, tuple_name: str = "ntuple", title: str = "") -> "RootTupleWriter":
        import ROOT  # type: ignore
        from ROOT import TFile, TNtuple  # type: ignore

        out_path = Path(out_file)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        tf = TFile(str(out_path), "RECREATE")
        tn = TNtuple(tuple_name, title or tuple_name, schema.as_root_branch_string())
        return cls(file=tf, tuple=tn, schema=schema)

    def fill_row(self, values: Sequence[float]) -> None:
        if self._closed:
            raise RuntimeError("Cannot fill: writer is closed.")
        if len(values) != len(self.schema.branches):
            raise ValueError(f"Row length {len(values)} does not match schema length {len(self.schema.branches)}")
        arr = array.array("f", [float(v) for v in values])
        self.tuple.Fill(arr)

    def fill_rows(self, rows: Iterable[Sequence[float]]) -> None:
        for row in rows:
            self.fill_row(row)

    def close(self) -> None:
        if self._closed:
            return
        # Store tuple and close file.
        self.file.WriteObjectAny(self.tuple, self.tuple.ClassName(), self.tuple.GetName())
        self.file.Close()
        self._closed = True
