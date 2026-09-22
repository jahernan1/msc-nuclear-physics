"""Shared fixtures."""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_ROOT = REPO_ROOT / "results"


@pytest.fixture(scope="session")
def results_root() -> Path:
    if not RESULTS_ROOT.is_dir():
        pytest.skip(f"results tree not present at {RESULTS_ROOT}")
    return RESULTS_ROOT


@pytest.fixture(scope="session")
def ca40_potential(results_root: Path) -> Path:
    path = results_root / "potentials" / "Ca40_fsugold_potential.dat"
    if not path.is_file():
        pytest.skip(f"missing reference potential: {path}")
    return path
