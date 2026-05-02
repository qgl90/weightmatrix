from __future__ import annotations

import sys
from pathlib import Path

def add_src_to_path() -> None:
    """
    Allow running scripts directly (e.g. `python3 scripts/run_it.py`) without
    requiring an editable install by prepending `<repo>/src` to `sys.path`.
    """
    repo_root = Path(__file__).resolve().parents[2]
    src_dir = repo_root / "src"
    if src_dir.is_dir():        
        sys.path.insert(0, str(src_dir))