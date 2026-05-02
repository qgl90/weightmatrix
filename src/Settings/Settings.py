from __future__ import annotations

import json
from pathlib import Path
from typing import Any
# _PATH_KEYS = {
#     "fieldMap",
#     "gdmlInput",
#     "storing_planes",
#     "path_geom_pickle",
# }

# _PATH_LIST_KEYS = {
#     "detectionlayers",
# }


# def _find_repo_root(start: Path) -> Path | None:
#     for candidate in [start, *start.parents]:
#         if (candidate / ".git").is_dir():
#             return candidate
#     return None
# def _reanchor_absolute_path(path: Path, *, repo_root: Path) -> Path | None:
#     """
#     Best-effort portability helper: if a settings file contains an absolute path
#     pointing to a different checkout, try to map it onto `<repo>/inputs/...`.
#     """
#     anchors = ("planes", 
#                "magfield", 
#                "gdmls", 
#                "xyshapes", 
#                "jsons", 
#                "kinematics")
#     parts = path.parts
#     for anchor in anchors:
#         if anchor in parts:
#             idx = parts.index(anchor)
#             subpath = Path(*parts[idx:])
#             candidate = repo_root / "inputs" / subpath
#             if candidate.exists():
#                 return candidate
#     return None


# def _candidate_paths(value: Path, *, base_dir: Path, repo_root: Path | None) -> list[Path]:
#     candidates = [base_dir / value]
#     if repo_root is not None:
#         candidates.append(repo_root / value)
#         candidates.append(repo_root / "inputs" / value)
#     return candidates


# def _resolve_path(value: str, *, base_dir: Path, repo_root: Path | None) -> str:
#     p = Path(value).expanduser()
#     if p.is_absolute():
#         if p.exists():
#             return str(p)
#         if repo_root is not None:
#             remapped = _reanchor_absolute_path(p, repo_root=repo_root)
#             if remapped is not None:
#                 return str(remapped.resolve())
#         return str(p)
#     candidates = _candidate_paths(p, base_dir=base_dir, repo_root=repo_root)
#     for candidate in candidates:
#         if candidate.exists():
#             return str(candidate.resolve())
#     # Fall back to the canonical "relative to json file" resolution
#     return str(candidates[0].resolve())


# def resolve_settings_paths(settings: dict[str, Any], *, base_dir: str | Path) -> dict[str, Any]:
#     print( f"resolve_settings_paths -- base_dir = {base_dir}")
#     """
#     Return a copy of the settings dict where file paths are resolved relative
#     to the JSON file location (instead of the current working directory).
#     """
#     from .Paths import expand_env_var    
#     resolved: dict[str, Any] = dict(settings)
#     for key in _PATH_KEYS:
#         if key in resolved and isinstance(resolved[key], str):            
#             resolved[key] = expand_env_var( resolved[key])            
#     for key in _PATH_LIST_KEYS:
#         if key in resolved and isinstance(resolved[key], list):
#             resolved[key] = [                
#                 expand_env_var(item) for item in resolved[key]              
#             ]
#     import os 
#     for _ in resolved[key]:
#         if not(os.path.exists( _)):
#             print(f"key = {resolved[key]} \n {_} not existing... something went bad")
#             exit()            

#     return resolved
import yaml

def load_settings(path: str | Path) -> dict[str, Any]:
    """
    Historical settings loader (returns a raw dict).

    For new code prefer `weightmatrix.config.Settings.from_yaml(...).as_dict()`,
    which also resolves paths relative to the YAML location.
    """
    settings_path = Path(path).expanduser().resolve()
    with settings_path.open("r", encoding="utf-8") as f:
        settings = yaml.safe_load(f) or {}
    if not isinstance(settings, dict):
        raise ValueError("Settings YAML must parse to a mapping/dict.")
    return settings
