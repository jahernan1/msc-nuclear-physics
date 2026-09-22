import numpy as np
import pytest

from nucth.cevns import (
    G_FERMI,
    SIN2_THETA_W,
    WeakFormFactor,
    differential_cross_section,
    max_recoil,
    nuclear_mass,
    total_cross_section,
    weak_charge,
)

UNITY_FF = WeakFormFactor(np.linspace(0.0, 5.0, 100), np.ones(100))


def test_weak_charge_is_positive_and_neutron_dominated():
    """Qw = N - (1 - 4 sin^2 theta_W) Z. The Z term is strongly suppressed."""
    qw = weak_charge(N=22, Z=18)
    assert qw > 0
    assert qw == pytest.approx(22 - (1 - 4 * SIN2_THETA_W) * 18, rel=1e-12)
    # Z contributes only ~7% of what an equal neutron count would.
    assert abs(1 - 4 * SIN2_THETA_W) < 0.08


def test_argon40_mass_is_near_37_gev():
    """Ar40: 18 protons, 22 neutrons, ~8.6 MeV/nucleon binding."""
    mass = nuclear_mass(N=22, Z=18)
    assert 37_000 < mass < 37_400


def test_max_recoil_grows_quadratically_at_low_neutrino_energy():
    """T_max = 2 E^2 / (2E + M) -> 2E^2/M when E << M."""
    mass = 37_200.0
    assert max_recoil(mass, 30.0) == pytest.approx(2 * 30.0**2 / (2 * 30.0 + mass), rel=1e-12)
    assert max_recoil(mass, 60.0) / max_recoil(mass, 30.0) == pytest.approx(4.0, rel=1e-2)


def test_cross_section_vanishes_at_maximum_recoil():
    """The kinematic bracket 1 - T/E - MT/(2E^2) is zero at T = T_max."""
    mass, energy = 37_200.0, 30.0
    t_max = max_recoil(mass, energy)
    qw = weak_charge(22, 18)

    at_max = differential_cross_section(mass, energy, t_max, qw, UNITY_FF)
    assert at_max == pytest.approx(0.0, abs=1e-45)


def test_cross_section_is_positive_below_maximum_recoil():
    mass, energy = 37_200.0, 30.0
    qw = weak_charge(22, 18)
    t_max = max_recoil(mass, energy)

    for frac in (0.01, 0.25, 0.5, 0.9):
        assert differential_cross_section(mass, energy, frac * t_max, qw, UNITY_FF) > 0


def test_total_cross_section_scales_as_energy_squared():
    """sigma ~ Gf^2 Qw^2 E^2 / (4 pi) in the coherent, low-energy limit."""
    mass, qw = 37_200.0, weak_charge(22, 18)

    sigma_20 = total_cross_section(mass, 20.0, qw, UNITY_FF)
    sigma_40 = total_cross_section(mass, 40.0, qw, UNITY_FF)

    assert sigma_40 / sigma_20 == pytest.approx(4.0, rel=5e-2)


def test_total_cross_section_has_the_expected_order_of_magnitude():
    """Ar40 at 30 MeV: order 1e-39 cm^2, the canonical CEvNS scale."""
    sigma = total_cross_section(37_200.0, 30.0, weak_charge(22, 18), UNITY_FF)
    assert 1e-40 < sigma < 1e-38


def test_form_factor_suppresses_the_cross_section_at_large_momentum_transfer():
    """A realistic falling F_w gives a smaller sigma than a flat F_w = 1."""
    q = np.linspace(0.0, 5.0, 200)
    falling = WeakFormFactor(q, np.exp(-(q**2)))

    mass, qw = 37_200.0, weak_charge(22, 18)
    assert total_cross_section(mass, 40.0, qw, falling) < total_cross_section(
        mass, 40.0, qw, UNITY_FF
    )


def test_fermi_constant_is_the_pdg_value():
    assert G_FERMI == pytest.approx(1.1663787e-11, rel=1e-12)
