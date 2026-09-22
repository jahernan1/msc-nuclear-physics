"""Spherical (sine) Fourier transform pair for radially symmetric densities.

For a spherically symmetric rho(r) the three-dimensional transform collapses to

    F(q) = (4*pi/q) * int_0^inf r sin(q r) rho(r) dr

and its inverse

    rho(r) = (1/(2 pi^2 r)) * int_0^inf q sin(q r) F(q) dq

Both integrals are evaluated with Simpson's rule on the supplied grid. The
q = 0 and r = 0 endpoints are regulated with a small offset rather than
special-cased, matching the original implementation's numerics.

Ported from `nuc-th/fourier.py`. Behavioral change: the inputs are copied, so
neither function mutates the caller's arrays.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import simpson

#: Offset substituted for a zero endpoint to avoid division by zero.
ZERO_REGULATOR = 1.0e-6


def _regulated(axis: np.ndarray) -> np.ndarray:
    """Return a copy of ``axis`` with a zero first element nudged off zero."""
    out = np.array(axis, dtype=float, copy=True)
    if out[0] == 0.0:
        out[0] = ZERO_REGULATOR
    return out


def fourier(q: np.ndarray, r: np.ndarray, pr: np.ndarray) -> np.ndarray:
    """Forward transform: r-space density to q-space form factor.

    Parameters
    ----------
    q : array
        Momentum-transfer grid [fm^-1].
    r : array
        Radial grid [fm].
    pr : array
        Density sampled on ``r`` [fm^-3].

    Returns
    -------
    array
        Form factor sampled on ``q``. ``F(0)`` equals the volume integral
        ``4*pi*int r^2 pr dr`` (i.e. the total charge or nucleon number).
    """
    q = _regulated(q)
    r = np.asarray(r, dtype=float)
    pr = np.asarray(pr, dtype=float)

    pq = np.array([simpson(r * np.sin(qval * r) * pr, x=r) for qval in q])
    return (4.0 * np.pi / q) * pq


def fourier_back(q: np.ndarray, r: np.ndarray, pq: np.ndarray) -> np.ndarray:
    """Inverse transform: q-space form factor to r-space density.

    Parameters
    ----------
    q : array
        Momentum-transfer grid [fm^-1].
    r : array
        Radial grid [fm] on which to evaluate the result.
    pq : array
        Form factor sampled on ``q``.

    Returns
    -------
    array
        Density sampled on ``r`` [fm^-3].
    """
    q = np.asarray(q, dtype=float)
    r = _regulated(r)
    pq = np.asarray(pq, dtype=float)

    pr = np.array([simpson(q * np.sin(q * rval) * pq, x=q) for rval in r])
    return (1.0 / (2.0 * np.pi**2 * r)) * pr
