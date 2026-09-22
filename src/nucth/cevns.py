"""Coherent elastic neutrino-nucleus scattering (CEvNS) cross sections.

For a neutrino of energy E scattering off a nucleus of mass M with nuclear
recoil energy T, the differential cross section is

    dsigma/dT = (Gf^2 / 4 pi) M Qw^2 [1 - T/E - M T / (2 E^2)] |F_w(Q)|^2

with Q = sqrt(2 M T) the momentum transfer and

    Qw = N - (1 - 4 sin^2 theta_W) Z

the weak charge of the nucleus. Because 1 - 4 sin^2 theta_W is approximately
0.076, the proton contribution is strongly suppressed and Qw is essentially the
neutron number: CEvNS is a neutron-distribution probe, and the weak form factor
F_w carries the neutron-skin dependence this thesis measures.

Ported from `nuc-th/CEvNS.py`.

Prefactor note
--------------
The original wrote ``(Gf^2 / 8 pi) M (2 - 2T/E - M T/E^2)``, which factors
exactly into the ``(Gf^2 / 4 pi) M (1 - T/E - M T / 2E^2)`` used here. An older
copy of the same file in the `JesseDocs` backup paired ``Gf^2 / 4 pi`` with
the *unfactored* bracket, double-counting by a factor of two. The form below is
the correct one.

Sign note
---------
The original also wrote ``Qw = -N + (1 - 4 sin^2) Z``, the negative of the
conventional definition. Qw enters only squared, so no published number
changes, but the sign is corrected here so that the weak charge and the weak
density share a consistent convention and no plot needs a compensating minus.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import splev, splrep

from nucth.paths import data_root

#: Fermi constant [MeV^-2].
G_FERMI = 1.1663787e-11
#: Weak mixing angle, sin^2(theta_W), at low momentum transfer.
SIN2_THETA_W = 0.231
#: Conversion constant [fm MeV].
HBARC = 197.326
#: (hbar c)^2 in cm^2 MeV^2, for converting MeV^-2 to cm^2.
MEV_INV2_TO_CM2 = HBARC**2 * 1e-26

#: Proton and neutron masses [MeV].
PROTON_MASS = 938.28
NEUTRON_MASS = 939.57


def weak_charge(N: int, Z: int) -> float:
    """Weak charge ``Qw = N - (1 - 4 sin^2 theta_W) Z``."""
    return N - (1.0 - 4.0 * SIN2_THETA_W) * Z


def binding_energy_per_nucleon(N: int, Z: int, root: str | Path | None = None) -> float | None:
    """Binding energy per nucleon [MeV] from ``Masses2016.txt``, or None."""
    path = (Path(root) if root is not None else data_root()) / "Masses2016.txt"
    table = np.loadtxt(path)
    match = table[(table[:, 0] == N) & (table[:, 1] == Z)]
    return float(match[0, 2]) if len(match) else None


def nuclear_mass(N: int, Z: int, root: str | Path | None = None) -> float:
    """Nuclear mass [MeV] from the constituent masses less the binding energy."""
    binding = binding_energy_per_nucleon(N, Z, root)
    if binding is None:
        raise KeyError(f"no binding energy tabulated for N={N}, Z={Z}")
    return N * NEUTRON_MASS + Z * PROTON_MASS - (N + Z) * binding


class WeakFormFactor:
    """Cubic-spline interpolant of ``F_w(q)`` that accepts ``Q`` in MeV.

    Parameters
    ----------
    q_axis : array
        Momentum-transfer grid in fm^-1, as written by
        :func:`nucth.formfactors.compute`.
    values : array
        ``F_w`` on that grid, normalized to 1 at q = 0.
    """

    def __init__(self, q_axis: np.ndarray, values: np.ndarray) -> None:
        self._tck = splrep(np.asarray(q_axis, float), np.asarray(values, float))

    def __call__(self, q_mev: float | np.ndarray) -> float | np.ndarray:
        """Evaluate at momentum transfer ``q_mev`` [MeV]."""
        return splev(np.asarray(q_mev, float) / HBARC, self._tck)

    @classmethod
    def from_file(cls, path: str | Path) -> WeakFormFactor:
        """Load from a ``*_WeakFF.dat`` file: columns ``q [fm^-1], F_w, F_ch``."""
        table = np.loadtxt(Path(path))
        return cls(table[:, 0], table[:, 1])


def max_recoil(M: float, E_nu: float) -> float:
    """Maximum nuclear recoil energy ``T_max = 2 E^2 / (2E + M)`` [MeV]."""
    return 2.0 * E_nu**2 / (2.0 * E_nu + M)


def differential_cross_section(
    M: float,
    E_nu: float,
    T: float,
    Qw: float,
    form_factor: Callable[[float], float],
) -> float:
    """``dsigma/dT`` in cm^2/MeV.

    Parameters
    ----------
    M : float
        Nuclear mass [MeV].
    E_nu : float
        Incoming neutrino energy [MeV].
    T : float
        Nuclear recoil energy [MeV].
    Qw : float
        Weak charge of the nucleus.
    form_factor : callable
        ``F_w(Q)`` with Q in MeV — e.g. a :class:`WeakFormFactor`.
    """
    q = np.sqrt(2.0 * M * T)
    kinematic = 1.0 - T / E_nu - (M * T) / (2.0 * E_nu**2)
    natural = (G_FERMI**2 / (4.0 * np.pi)) * M * Qw**2 * kinematic * form_factor(q) ** 2
    return float(natural * MEV_INV2_TO_CM2)


def total_cross_section(
    M: float,
    E_nu: float,
    Qw: float,
    form_factor: Callable[[float], float],
    *,
    n_T: int = 500,
    T_min: float = 1.0e-6,
) -> float:
    """``sigma(E_nu)`` in cm^2, integrated over recoil energy."""
    t_axis = np.linspace(T_min, max_recoil(M, E_nu), n_T)
    dsig = np.array(
        [differential_cross_section(M, E_nu, t, Qw, form_factor) for t in t_axis]
    )
    return float(simpson(dsig, x=t_axis))
