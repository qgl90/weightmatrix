from __future__ import annotations

"""
Typed settings with path resolution.

The historical code uses a plain `dict` from YAML. This module introduces a
dataclass wrapper that:

- validates required keys exist
- resolves file paths relative to the YAML file location
- expands env vars in paths (`$VAR` / `${VAR}`)

This is designed to be adopted incrementally:
existing code can continue to accept `dict`, while Apps can start using
`Settings.from_yaml(...).as_dict()`.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from Core.Paths import expand_env_var


_PATH_KEYS = {
    "fieldMap",
    "gdmlInput",
    "storing_planes",
    "path_geom_pickle",
}

_PATH_LIST_KEYS = {
    "detectionlayers",
}


def _resolve_path(value: str, *, base_dir: Path) -> str:
    expanded = Path(expand_env_var(value)).expanduser()
    if expanded.is_absolute():
        return str(expanded)

    # Prefer paths relative to the YAML file location, but fall back to a
    # best-effort "project root" resolution to support historical configs where
    # settings live under `layouts/` but paths are repo-root-relative.
    cand1 = (base_dir / expanded).resolve()
    if cand1.exists():
        return str(cand1)

    repo_root = _find_project_root(base_dir)
    if repo_root is not None:
        cand2 = (repo_root / expanded).resolve()
        if cand2.exists():
            return str(cand2)

    # Last resort: keep it relative to base_dir.
    return str(cand1)


def _find_project_root(start: Path) -> Path | None:
    for candidate in [start, *start.parents]:
        if (candidate / "pyproject.toml").exists() and (candidate / "src").is_dir():
            return candidate
    return None


@dataclass(frozen=True)
class Settings:
    raw: dict[str, Any]
    base_dir: Path

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Settings":
        settings_path = Path(path).expanduser().resolve()
        with settings_path.open("r", encoding="utf-8") as f:
            raw = yaml.safe_load(f) or {}
        if not isinstance(raw, dict):
            raise ValueError("Settings YAML must parse to a mapping/dict.")
        return cls(raw=raw, base_dir=settings_path.parent)

    def resolved(self) -> dict[str, Any]:
        """
        Return a copy of the settings dict with resolved/expanded paths.
        """
        out: dict[str, Any] = dict(self.raw)

        for key in _PATH_KEYS:
            if key in out and isinstance(out[key], str):
                out[key] = _resolve_path(out[key], base_dir=self.base_dir)

        for key in _PATH_LIST_KEYS:
            if key in out and isinstance(out[key], list):
                out[key] = [_resolve_path(str(v), base_dir=self.base_dir) for v in out[key]]

        return out

    def as_dict(self) -> dict[str, Any]:
        return self.resolved()
