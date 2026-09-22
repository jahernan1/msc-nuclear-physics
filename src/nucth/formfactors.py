"""Charge and weak form factors and densities of a nucleus.

The point proton and neutron densities from the Dirac solver are Fourier
transformed and folded with the single-nucleon electric form factors Gep, Gen
(``data/GPN.txt``) to give the nuclear charge form factor

    F_ch(q) = Gep(q) F_p(q) + Gen(q) F_n(q)

and, with the weak vector charges of the nucleons, the weak form factor

    Gwp = gvp Gep + gvn Gen        gvp = +0.0712
    Gwn = gvn Gep + gvp Gen        gvn = -0.9877
    F_w(q) = Gwp(q) F_p(q) + Gwn(q) F_n(q)

The neutron carries essentially all of the weak charge, which is why coherent
neutrino scattering is sensitive to the neutron distribution and hence to the
neutron skin r_n - r_p.

Ported from `nuc-th/Density.py`. Results are returned in a dataclass rather
than smuggled out through module-level globals.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

from nucth.paths import data_root
from nucth.transforms import fourier, fourier_back

#: Weak vector charge of the proton.
WEAK_CHARGE_PROTON = 0.0712
#: Weak vector charge of the neutron.
WEAK_CHARGE_NEUTRON = -0.9877

#: Number of q points retained from GPN.txt. The original used 1501.
DEFAULT_N_Q = 1501


def load_nucleon_form_factors(root: str | Path | None = None):
    """Load ``(q, Gep, Gen)`` from ``GPN.txt``.

    Columns: momentum transfer [fm^-1], proton electric form factor, neutron
    electric form factor.
    """
    path = Path(root) if root is not None else data_root()
    table = np.loadtxt(path / "GPN.txt")
    return table[:, 0], table[:, 1], table[:, 2]


def weak_nucleon_form_factors(gep: np.ndarray, gen: np.ndarray):
    """Fold the electric form factors with the weak vector charges.

    Returns
    -------
    (gwn, gwp)
        Neutron and proton weak form factors, in that order.
    """
    gwp = WEAK_CHARGE_PROTON * gep + WEAK_CHARGE_NEUTRON * gen
    gwn = WEAK_CHARGE_NEUTRON * gep + WEAK_CHARGE_PROTON * gen
    return gwn, gwp


def _rms_radius(r: np.ndarray, density: np.ndarray, normalization: float) -> float:
    """sqrt( 4 pi / Q * int r^4 rho dr )."""
    return float(np.sqrt((4.0 * np.pi / normalization) * simpson(r**4 * density, x=r)))


def _total_charge(r: np.ndarray, density: np.ndarray) -> float:
    """4 pi * int r^2 rho dr."""
    return float(4.0 * np.pi * simpson(r**2 * density, x=r))


@dataclass(frozen=True)
class NuclearFormFactors:
    """Charge and weak structure of one nucleus in one mean field."""

    r: np.ndarray            # [fm]
    q: np.ndarray            # [fm^-1]
    charge_density: np.ndarray   # [fm^-3]
    weak_density: np.ndarray     # [fm^-3]
    charge_ff: np.ndarray        # normalized to 1 at q = 0
    weak_ff: np.ndarray          # normalized to 1 at q = 0
    Q_charge: float              # = Z
    Q_weak: float
    r_p: float
    r_n: float
    r_ch: float
    r_w: float

    @property
    def neutron_skin(self) -> float:
        """r_n - r_p [fm]."""
        return self.r_n - self.r_p

    @property
    def weak_skin(self) -> float:
        """r_w - r_ch [fm]."""
        return self.r_w - self.r_ch


def compute(
    nucleon_density_path: str | Path,
    *,
    n_q: int = DEFAULT_N_Q,
    root: str | Path | None = None,
) -> NuclearFormFactors:
    """Charge and weak form factors from a nucleon density file.

    Parameters
    ----------
    nucleon_density_path : path
        Five-column file ``r [fm], rho_s(n), rho_v(n), rho_s(p), rho_v(p)``
        as written by :func:`nucth.cli.run_structure`.
    n_q : int
        Number of q points to retain from ``GPN.txt``.
    root : path, optional
        Data directory override; defaults to :func:`nucth.paths.data_root`.
    """
    table = np.loadtxt(Path(nucleon_density_path))
    r = table[:, 0]
    rho_n = table[:, 2]   # neutron vector density
    rho_p = table[:, 4]   # proton vector density

    q_full, gep_full, gen_full = load_nucleon_form_factors(root)
    gwn_full, gwp_full = weak_nucleon_form_factors(gep_full, gen_full)

    q = q_full[:n_q]
    gep, gen = gep_full[:n_q], gen_full[:n_q]
    gwn, gwp = gwn_full[:n_q], gwp_full[:n_q]

    f_p = fourier(q, r, rho_p)
    f_n = fourier(q, r, rho_n)
    Z, N = f_p[0], f_n[0]   # F(0) is the particle number

    f_ch = gep * f_p + gen * f_n
    f_w = gwp * f_p + gwn * f_n

    rho_ch = fourier_back(q, r, f_ch)
    rho_w = fourier_back(q, r, f_w)

    Q_charge = _total_charge(r, rho_ch)
    Q_weak = _total_charge(r, rho_w)

    return NuclearFormFactors(
        r=r,
        q=q,
        charge_density=rho_ch,
        weak_density=rho_w,
        charge_ff=f_ch / Q_charge,
        weak_ff=f_w / Q_weak,
        Q_charge=Q_charge,
        Q_weak=Q_weak,
        r_p=_rms_radius(r, rho_p, Z),
        r_n=_rms_radius(r, rho_n, N),
        r_ch=_rms_radius(r, rho_ch, Q_charge),
        r_w=_rms_radius(r, rho_w, Q_weak),
    )
