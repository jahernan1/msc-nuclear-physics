"""Data-root resolution.

Every path in this package resolves through here. No other module contains a
filesystem path literal.
"""

from __future__ import annotations

import os
from pathlib import Path

_PACKAGE_DATA = Path(__file__).resolve().parent / "data"
_REPO_DATA = Path(__file__).resolve().parents[2] / "data"


def data_root() -> Path:
    """Directory holding GPN.txt, Masses2016.txt and parameters/.

    Resolution order:

    1. ``NUCTH_DATA`` environment variable, if set.
    2. ``data/`` bundled inside an installed wheel.
    3. ``data/`` at the repository root, for an editable/source checkout.

    Raises
    ------
    FileNotFoundError
        If no candidate directory exists.
    """
    env = os.environ.get("NUCTH_DATA")
    if env:
        candidate = Path(env).expanduser().resolve()
        if not candidate.is_dir():
            raise FileNotFoundError(f"NUCTH_DATA points at a missing directory: {candidate}")
        return candidate

    for candidate in (_PACKAGE_DATA, _REPO_DATA):
        if candidate.is_dir():
            return candidate

    raise FileNotFoundError(
        "Could not locate the nucth data directory. Set NUCTH_DATA to the "
        "directory containing GPN.txt and Masses2016.txt."
    )


DATA_ROOT = data_root()