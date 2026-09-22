"""Coupled radial Dirac equations in a relativistic mean field.

For a nucleon of mass M in scalar potential S(r) and vector potential V(r), the
upper and lower radial components g, f obey

    dg/dx = (E - V + S + 1) f - (kappa / x) g
    df/dx = (kappa / x) f - (E - V - S - 1) g

in dimensionless form: x = r M / hbarc, energies in units of M, potentials in
units of M. kappa is the Dirac quantum number (j = |kappa| - 1/2; l = kappa for
kappa > 0, l = -1 - kappa for kappa < 0).

Bound states are found by shooting and matching: integrate outward from the
origin to the matching radius with RK4, inward from the box boundary back to the
same radius, and root-find the determinant g_out f_in - g_in f_out with Brent's
method. Occupied levels are filled to N and Z by ascending energy, and the
scalar and vector densities follow as

    rho_v(r) = sum_occ (2j+1) / (4 pi r^2) * (g^2 + f^2)
    rho_s(r) = sum_occ (2j+1) / (4 pi r^2) * (g^2 - f^2)

Ported from `nuc-th/NuclearClass.py`.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import splev, splrep
from scipy.optimize import brentq

NEUTRON = "neutron"
PROTON = "proton"
_SPECIES = (NEUTRON, PROTON)

#: Selects the single-species branch of ``_func`` (the neutron/proton flag the
#: original class called ``n`` when integrating one species' wavefunction).
_FUNC_ARG = {NEUTRON: 0, PROTON: 1}

#: Upper bound of the bound-state search, in units of the nucleon mass.
E_MAX = 0.994
#: Initial energy of the bound-state scan, in units of the nucleon mass.
E_START = 0.9
#: Coarse step of the bracketing scan, in units of the nucleon mass.
E_STEP = 0.012


def binding_energies(levels: Sequence[Sequence[float]], mass: float = 939.0) -> np.ndarray:
    """Convert ``[kappa, E/M]`` level pairs to binding energies ``E - M`` in MeV."""
    return np.array([level[1] * mass - mass for level in levels], dtype=float)


class Nuclear:
    """Dirac solver for one nucleus in one mean field.

    Parameters
    ----------
    potential_path : path
        Five-column file ``r [fm], S_n, V_n, S_p, V_p [MeV]`` as written by
        :func:`nucth.potentials.save_potential`.
    N, Z : int
        Neutron and proton numbers.
    r_match : float, default 7.0
        Matching radius [fm]. Increase for heavy nuclei.
    r_box : float, default 15.0
        Outer box boundary [fm]. Must exceed ``r_match``.
    """

    MASS = 939.0     # [MeV]
    HBARC = 197.326  # [fm MeV]

    def __init__(
        self,
        potential_path: str | Path,
        N: int,
        Z: int,
        *,
        r_match: float = 7.0,
        r_box: float = 15.0,
    ) -> None:
        if r_box <= r_match:
            raise ValueError(f"r_box ({r_box}) must exceed r_match ({r_match})")

        self.potential_path = Path(potential_path)
        self.N = int(N)
        self.Z = int(Z)
        self.kappa = -1  # per-instance, not shared across instances

        table = np.loadtxt(self.potential_path)
        self.raxis = table[:, 0]
        self.xaxis = table[:, 0] * (self.MASS / self.HBARC)
        self.neutronS = table[:, 1] / self.MASS
        self.neutronV = table[:, 2] / self.MASS
        self.protonS = table[:, 3] / self.MASS
        self.protonV = table[:, 4] / self.MASS

        self.tck_NS = splrep(self.xaxis, self.neutronS)
        self.tck_NV = splrep(self.xaxis, self.neutronV)
        self.tck_PS = splrep(self.xaxis, self.protonS)
        self.tck_PV = splrep(self.xaxis, self.protonV)

        self.h = self.xaxis[1] - self.xaxis[0]
        self.a = 0.0
        self.b = r_match * (self.MASS / self.HBARC)
        self.c = r_box * (self.MASS / self.HBARC)
        self.axis = self.xaxis[: int(self.c / self.h)]
        self.midindex = int(self.b / self.h) - 1

    @property
    def radial_axis(self) -> np.ndarray:
        """Solver axis in fm."""
        return self.axis * (self.HBARC / self.MASS)

    # ------------------------------------------------------------------
    # numerical core (bodies unchanged from NuclearClass.py except naming)
    # ------------------------------------------------------------------

    def _func(self, r, x, *args):
        if x == 0:
            x = 1.0e-6
        vn = splev(x, self.tck_NV)
        sn = splev(x, self.tck_NS)
        vp = splev(x, self.tck_PV)
        sp = splev(x, self.tck_PS)

        if args[0] == 0:
            enn = args[1]
            gn = r[0]
            fn = r[1]
            dgn = (enn - vn + sn + 1) * fn - (self.kappa / x) * gn
            dfn = (self.kappa / x) * fn - (enn - vn - sn - 1) * gn
            return np.array([dgn, dfn], float)

        elif args[0] == 1:
            enp = args[1]
            gp = r[0]
            fp = r[1]
            dgp = (enp - vp + sp + 1) * fp - (self.kappa / x) * gp
            dfp = (self.kappa / x) * fp - (enp - vp - sp - 1) * gp
            return np.array([dgp, dfp], float)

        else:
            enn = args[0]
            enp = args[1]

            gn = r[0]
            fn = r[1]
            gp = r[2]
            fp = r[3]

            dgn = (enn - vn + sn + 1) * fn - (self.kappa / x) * gn
            dfn = (self.kappa / x) * fn - (enn - vn - sn - 1) * gn
            dgp = (enp - vp + sp + 1) * fp - (self.kappa / x) * gp
            dfp = (self.kappa / x) * fp - (enp - vp - sp - 1) * gp

            return np.array([dgn, dfn, dgp, dfp], float)

    def _matchpnt(self, first, last, step, enn, enp):
        if self.kappa < 0:
            r = np.array([0.0001, 0.0, 0.0001, 0.0], float)
        else:
            r = np.array([0.0, 0.0001, 0.0, 0.0001], float)

        for x in np.arange(first, last, step):
            k1 = step * self._func(r, x, enn, enp)
            k2 = step * self._func(r + (0.5 * k1), x + 0.5 * step, enn, enp)
            k3 = step * self._func(r + 0.5 * k2, x + 0.5 * step, enn, enp)
            k4 = step * self._func(r + k3, x + step, enn, enp)
            r += (k1 + 2 * (k2 + k3) + k4) / 6

        return r

    def _det(self, enn, enp):
        rint = self._matchpnt(self.a, self.b, self.h, enn, enp)
        rext = self._matchpnt(self.c, self.b, -self.h, enn, enp)

        gn1, gp1 = rint[0], rint[2]
        gn2, gp2 = rext[0], rext[2]
        fn1, fp1 = rint[1], rint[3]
        fn2, fp2 = rext[1], rext[3]

        ndet, pdet = (gn2 * fn1) - (gn1 * fn2), (gp2 * fp1) - (gp1 * fp2)

        return ndet, pdet

    @staticmethod
    def _sign_change(det1: float, det2: float) -> bool:
        """True when the determinant brackets a root between two energies."""
        return det1 * det2 < 0.0

    def _neutron_det(self, en, partner):
        return self._det(en, partner)[0]

    def _proton_det(self, en, partner):
        return self._det(partner, en)[1]

    def _species_number(self, species: str) -> int:
        if species == NEUTRON:
            return self.N
        if species == PROTON:
            return self.Z
        raise ValueError(f"species must be one of {_SPECIES}, got {species!r}")

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------

    def _scan_species(self, species: str, kappa: int) -> list:
        """Bracket-and-refine bound states of one species at one kappa.

        The neutron and proton determinants in ``_det`` never actually depend
        on each other's trial energy (each branch of ``_func`` only reads its
        own species' potential), so the energy passed as the "other" species'
        argument is a fixed dummy rather than live loop state from a coupled
        scan (edit 8: the original let the neutron loop's leftover energy
        leak into the proton scan's determinant calls).
        """
        self.kappa = kappa
        det_fn = self._neutron_det if species == NEUTRON else self._proton_det
        dummy = E_MAX

        levels = []
        en = E_START
        while en < E_MAX:
            en2 = en + E_STEP
            det1, det2 = det_fn(en, dummy), det_fn(en2, dummy)

            while not self._sign_change(det1, det2):
                if en > E_MAX:
                    break
                en += E_STEP
                en2 = en + E_STEP
                det1 = det2
                det2 = det_fn(en2, dummy)

            if en > E_MAX:
                break

            en_root = brentq(det_fn, en, en2, args=(dummy,))
            levels.append([kappa, en_root])
            en = en2

        return levels

    def spectrum(self, kappa_values: Sequence[int]) -> tuple[list, list]:
        """Bound-state spectrum over a range of kappa values.

        Returns
        -------
        (neutron_levels, proton_levels)
            Each a list-of-lists: one list of ``[kappa, E/M]`` per kappa.
        """
        neutron_levels = []
        proton_levels = []
        for kappa in kappa_values:
            neutron_levels.append(self._scan_species(NEUTRON, kappa))
            proton_levels.append(self._scan_species(PROTON, kappa))
        return neutron_levels, proton_levels

    @staticmethod
    def _sorted_levels(levels) -> list:
        flat = [level for kappa_levels in levels for level in kappa_levels]
        flat.sort(key=lambda level: level[1])
        return flat

    def occupied(self, species: str, levels) -> list:
        """Fill levels by ascending energy up to N (neutron) or Z (proton)."""
        target = self._species_number(species)
        occupied_levels = []
        filled = 0.0

        for kappa, en in self._sorted_levels(levels):
            if filled >= target:
                break
            j = abs(kappa) - 0.5
            filled += 2 * j + 1
            occupied_levels.append([kappa, en])

        return occupied_levels

    def wavefunction(self, species: str, level) -> tuple[list, np.ndarray, np.ndarray]:
        """Radial wavefunction ``(g, f)`` for one occupied level."""
        kappa, en = level
        self.kappa = kappa
        arg0 = _FUNC_ARG[species]

        g_wf = np.zeros(len(self.axis) + 1, float)
        f_wf = np.zeros(len(self.axis) + 1, float)

        if self.kappa < 0:
            r = np.array([0.0001, 0.0], float)
        else:
            r = np.array([0.0, 0.0001], float)

        i = 0
        for xint in np.arange(self.a, self.b + self.h, self.h):
            g_wf[i] = r[0]
            f_wf[i] = r[1]
            i += 1

            k1 = self.h * self._func(r, xint, arg0, en)
            k2 = self.h * self._func(r + (0.5 * k1), xint + 0.5 * self.h, arg0, en)
            k3 = self.h * self._func(r + 0.5 * k2, xint + 0.5 * self.h, arg0, en)
            k4 = self.h * self._func(r + k3, xint + self.h, arg0, en)
            r += (k1 + 2 * (k2 + k3) + k4) / 6

        if self.kappa < 0:
            r = np.array([0.0001, 0.0], float)
        else:
            r = np.array([0.0, 0.0001], float)

        j = len(self.axis)

        for xext in np.arange(self.c, self.b - self.h, -self.h):
            g_wf[j] = r[0]
            f_wf[j] = r[1]
            j -= 1

            k1 = -self.h * self._func(r, xext, arg0, en)
            k2 = -self.h * self._func(r + (0.5 * k1), xext + 0.5 * -self.h, arg0, en)
            k3 = -self.h * self._func(r + 0.5 * k2, xext + 0.5 * -self.h, arg0, en)
            k4 = -self.h * self._func(r + k3, xext + -self.h, arg0, en)
            r += (k1 + 2 * (k2 + k3) + k4) / 6

        matchcoef = g_wf[self.midindex] / g_wf[self.midindex + 1]

        g_wf[self.midindex :] *= matchcoef
        f_wf[self.midindex :] *= matchcoef

        g_wf = np.delete(g_wf, self.midindex)
        f_wf = np.delete(f_wf, self.midindex)

        norm_const = np.sqrt(1 / simpson(g_wf**2 + f_wf**2, x=self.axis))

        g_wf = norm_const * g_wf
        f_wf = norm_const * f_wf

        return [self.kappa, en], g_wf, f_wf

    def wavefunctions(self, species: str, levels) -> list:
        return [self.wavefunction(species, level) for level in levels]

    def density(self, species: str, wavefunctions) -> tuple[np.ndarray, np.ndarray]:
        """Scalar and vector densities in fm^-3, filled to N or Z."""
        target = self._species_number(species)

        # Regulate r = 0 on a local copy; never mutate self.axis (review item 7).
        axis = np.array(self.axis, dtype=float, copy=True)
        if axis[0] == 0.0:
            axis[0] = 1.0e-6

        scalar = np.zeros_like(axis)
        vector = np.zeros_like(axis)
        filled = 0.0

        for (kappa, _en), g_wf, f_wf in wavefunctions:
            j = abs(kappa) - 0.5
            n_occ = 2.0 * j + 1.0
            filled += n_occ
            if filled > target:              # partially filled outermost shell
                n_occ = target - (filled - n_occ)

            weight = n_occ / (4.0 * np.pi * axis**2)
            vector += weight * (g_wf**2 + f_wf**2)
            scalar += weight * (g_wf**2 - f_wf**2)

        unit = (self.MASS / self.HBARC) ** 3
        return scalar * unit, vector * unit
