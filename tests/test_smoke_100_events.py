from __future__ import annotations

import os
import subprocess
from pathlib import Path


def test_smoke_100_events(tmp_path: Path) -> None:
    """
    CI smoke test:
    - runs 100 events through the runner
    - verifies a ROOT output is produced and contains an 'ntuple'
    """
    try:
        import ROOT  # type: ignore
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "PyROOT is required for this test. In CI use a ROOT-enabled environment/container."
        ) from e

    repo = Path(__file__).resolve().parents[1]
    out_root = tmp_path / "smoke.root"

    env = dict(os.environ)
    env["PYTHONPATH"] = str(repo / "src") + os.pathsep + env.get("PYTHONPATH", "")

    cmd = [
        "python3",
        str(repo / "src" / "Apps" / "run_sampled_tracks.py"),
        "--kinematic_file",
        str(repo / "samples" / "bsmumu_rdf.parquet"),
        "--settings_file",
        str(repo / "layouts" / "u2_debug_material.yaml"),
        "--output_file",
        str(out_root),
        "--start_entry",
        "0",
        "--stop_entry",
        "100",
    ]
    subprocess.run(cmd, check=True, cwd=repo, env=env)

    produced = Path(str(out_root).replace(".root", "_0_100.root"))
    assert produced.exists(), f"Expected output ROOT file not found: {produced}"

    f = ROOT.TFile.Open(str(produced), "READ")
    assert f and not f.IsZombie()
    t = f.Get("ntuple")
    assert t is not None, "Missing TNtuple 'ntuple' in output"
    assert int(t.GetEntries()) > 0

