"""Characterization tests for the coupled Dirac solver.

`nuc-th/Energy/Ca40_fsugold_energy.txt` was produced by the code being ported.
If the port changes any single-particle energy, these tests fail.
"""

import numpy as np
import pytest

from nucth.dirac import Nuclear, binding_energies

# From nuc-th/Energy/Ca40_fsugold_energy.txt, neutron column, MeV.
CA40_FSUGOLD_NEUTRON_MEV = [
    -54.63107378641769,
    -38.90327656808108,
    -35.167732356029774,
    -23.247723394421882,
    -16.938560924251874,
    -16.911425584187327,
]
CA40_FSUGOLD_NEUTRON_KAPPA = [-1, -2, 1, -3, 2, -1]

# Proton column of the same file, MeV.
CA40_FSUGOLD_PROTON_MEV = [
    -46.20828495117985,
    -30.874232640446166,
    -27.107982395624163,
    -15.603630774916496,
    -9.361886319340215,
    -9.339538039115041,
]
CA40_FSUGOLD_PROTON_KAPPA = [-1, -2, 1, -3, -1, 2]


@pytest.fixture(scope="module")
def ca40(ca40_potential):
    return Nuclear(ca40_potential, N=20, Z=20)


def test_radial_axis_starts_at_the_origin_and_is_uniform(ca40):
    axis = ca40.radial_axis
    assert axis[0] == pytest.approx(0.0, abs=1e-9)
    spacing = np.diff(axis)
    np.testing.assert_allclose(spacing, spacing[0], rtol=1e-9)


@pytest.mark.slow
def test_ca40_fsugold_neutron_spectrum_matches_the_committed_reference(ca40):
    neutron_levels, _ = ca40.spectrum(list(range(-4, 0)) + list(range(1, 4)))
    occupied = ca40.occupied("neutron", neutron_levels)

    kappas = [level[0] for level in occupied]
    energies = binding_energies(occupied)

    assert kappas == CA40_FSUGOLD_NEUTRON_KAPPA
    np.testing.assert_allclose(energies, CA40_FSUGOLD_NEUTRON_MEV, rtol=1e-6)


@pytest.mark.slow
def test_ca40_fsugold_proton_spectrum_matches_the_committed_reference(ca40):
    _, proton_levels = ca40.spectrum(list(range(-4, 0)) + list(range(1, 4)))
    occupied = ca40.occupied("proton", proton_levels)

    kappas = [level[0] for level in occupied]
    energies = binding_energies(occupied)

    assert kappas == CA40_FSUGOLD_PROTON_KAPPA
    np.testing.assert_allclose(energies, CA40_FSUGOLD_PROTON_MEV, rtol=1e-6)


@pytest.mark.slow
def test_neutron_vector_density_integrates_to_the_neutron_number(ca40):
    neutron_levels, _ = ca40.spectrum(list(range(-4, 0)) + list(range(1, 4)))
    occupied = ca40.occupied("neutron", neutron_levels)
    wfs = ca40.wavefunctions("neutron", occupied)
    _, vector = ca40.density("neutron", wfs)

    from scipy.integrate import simpson

    r = ca40.radial_axis
    particle_number = 4.0 * np.pi * simpson(r**2 * vector, x=r)
    assert particle_number == pytest.approx(20.0, rel=2e-2)


@pytest.mark.slow
def test_scalar_density_is_smaller_than_vector_density_everywhere(ca40):
    """rho_s = sum (g^2 - f^2) and rho_v = sum (g^2 + f^2), so rho_s <= rho_v."""
    neutron_levels, _ = ca40.spectrum(list(range(-4, 0)) + list(range(1, 4)))
    occupied = ca40.occupied("neutron", neutron_levels)
    wfs = ca40.wavefunctions("neutron", occupied)
    scalar, vector = ca40.density("neutron", wfs)

    assert np.all(scalar <= vector + 1e-12)


def test_two_instances_do_not_share_kappa_state(ca40_potential):
    """kappa was a class attribute in the original; it must now be per-instance."""
    a = Nuclear(ca40_potential, N=20, Z=20)
    b = Nuclear(ca40_potential, N=20, Z=20)

    a.kappa = -3
    assert b.kappa == -1
