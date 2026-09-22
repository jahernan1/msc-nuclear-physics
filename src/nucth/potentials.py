"""Convert RMF meson fields into the nucleon scalar and vector potentials.

The RMF solver (see ``external/rmf/``) writes ``potential0.txt`` with columns

    r [fm], Phi [MeV], W [MeV], B [MeV], A [MeV]

where Phi is the sigma (scalar) field, W the omega (isoscalar-vector) field,
B the rho (isovector-vector) field, and A the Coulomb field. A nucleon of a
given isospin sees

    S(p) = S(n) = -Phi
    V(n) = W - B/2
    V(p) = W + B/2 + A

Ported from `nuc-th/convertPotential.py`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

#: Column header written to converted potential files.
POTENTIAL_HEADER = "raxis[fm], SN, VN, SP, VP [MeV]"

#: Column indices of a converted potential array.
COL_R, COL_SN, COL_VN, COL_SP, COL_VP = range(5)


def convert(raw: np.ndarray) -> np.ndarray:
    """Map raw meson fields to nucleon scalar/vector potentials.

    Parameters
    ----------
    raw : array, shape (n, 5)
        Columns ``r [fm], Phi, W, B, A`` in MeV, as written by the RMF solver.

    Returns
    -------
    array, shape (n, 5)
        Columns ``r [fm], S_n, V_n, S_p, V_p`` in MeV.
    """
    raw = np.asarray(raw, dtype=float)
    if raw.ndim != 2 or raw.shape[1] != 5:
        raise ValueError(f"expected a raw potential of shape (n, 5), got {raw.shape}")

    radius = raw[:, 0]
    phi, w, b, a = raw[:, 1], raw[:, 2], raw[:, 3], raw[:, 4]

    scalar = -phi                 # identical for both isospins
    vector_n = w - 0.5 * b
    vector_p = w + 0.5 * b + a    # Coulomb acts on protons only

    return np.stack((radius, scalar, vector_n, scalar, vector_p), axis=-1)


def load_raw(path: str | Path) -> np.ndarray:
    """Load a five-column raw RMF potential file."""
    raw = np.loadtxt(Path(path))
    if raw.ndim != 2 or raw.shape[1] != 5:
        raise ValueError(f"{path}: expected 5 columns, got shape {raw.shape}")
    return raw


def save_potential(path: str | Path, potential: np.ndarray) -> None:
    """Write a converted potential array with the standard header."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.asarray(potential, dtype=float), header=POTENTIAL_HEADER)
