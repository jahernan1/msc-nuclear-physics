"""Shared fixtures.

RESULTS_ROOT points at the legacy `nuc-th/` tree until Task 9 migrates the
data files into `results/`. Task 9 updates the constant below; no test body
needs to change.
"""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_ROOT = REPO_ROOT / "nuc-th"


@pytest.fixture(scope="session")
def results_root() -> Path:
    if not RESULTS_ROOT.is_dir():
        pytest.skip(f"results tree not present at {RESULTS_ROOT}")
    return RESULTS_ROOT


@pytest.fixture(scope="session")
def ca40_potential(results_root: Path) -> Path:
    path = results_root / "Potentials" / "Ca40_fsugold_potential.dat"
    if not path.is_file():
        pytest.skip(f"missing reference potential: {path}")
    return path
