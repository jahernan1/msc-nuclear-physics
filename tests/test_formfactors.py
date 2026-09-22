import numpy as np
import pytest

from nucth.formfactors import (
    WEAK_CHARGE_NEUTRON,
    WEAK_CHARGE_PROTON,
    compute,
    load_nucleon_form_factors,
    weak_nucleon_form_factors,
)

# nuc-th/Skins/Ca40_neutron_skins.txt, "### fsugold" row.
CA40_FSUGOLD = {
    "r_p": 3.33771,
    "r_n": 3.28643,
    "neutron_skin": -0.05128,
    "r_ch": 3.42651,
    "r_w": 3.37267,
    "weak_skin": -0.05384,
}


@pytest.fixture(scope="module")
def ca40_density(results_root):
    path = results_root / "densities" / "Ca40_fsugold_nucleon_density.dat"
    if not path.is_file():
        pytest.skip(f"missing reference density: {path}")
    return path


def test_nucleon_form_factors_load_with_the_expected_shape():
    q, gep, gen = load_nucleon_form_factors()
    assert q.shape == gep.shape == gen.shape == (8001,)
    assert gep[0] == pytest.approx(1.0, abs=1e-12)   # Gep(0) = 1 (proton charge)
    assert gen[0] == pytest.approx(0.0, abs=1e-12)   # Gen(0) = 0 (neutron neutral)


def test_weak_nucleon_form_factors_apply_the_weak_vector_charges():
    q, gep, gen = load_nucleon_form_factors()
    gwn, gwp = weak_nucleon_form_factors(gep, gen)

    # At q = 0 the form factors reduce to the bare weak vector charges.
    assert gwp[0] == pytest.approx(WEAK_CHARGE_PROTON, abs=1e-9)
    assert gwn[0] == pytest.approx(WEAK_CHARGE_NEUTRON, abs=1e-9)


def test_weak_charge_of_the_neutron_dominates_the_proton():
    """The neutron carries nearly all the weak charge: this is why CEvNS probes rn."""
    assert abs(WEAK_CHARGE_NEUTRON) > 10.0 * abs(WEAK_CHARGE_PROTON)


@pytest.mark.slow
def test_ca40_fsugold_radii_match_the_committed_reference(ca40_density):
    ff = compute(ca40_density)

    assert ff.r_p == pytest.approx(CA40_FSUGOLD["r_p"], abs=5e-5)
    assert ff.r_n == pytest.approx(CA40_FSUGOLD["r_n"], abs=5e-5)
    assert ff.r_ch == pytest.approx(CA40_FSUGOLD["r_ch"], abs=5e-5)
    assert ff.r_w == pytest.approx(CA40_FSUGOLD["r_w"], abs=5e-5)


@pytest.mark.slow
def test_ca40_fsugold_skins_match_the_committed_reference(ca40_density):
    ff = compute(ca40_density)

    assert ff.neutron_skin == pytest.approx(CA40_FSUGOLD["neutron_skin"], abs=5e-5)
    assert ff.weak_skin == pytest.approx(CA40_FSUGOLD["weak_skin"], abs=5e-5)


@pytest.mark.slow
def test_ca40_charge_is_the_proton_number(ca40_density):
    """Ca40 has Z = 20, so the electric charge must come out at 20."""
    ff = compute(ca40_density)
    assert ff.Q_charge == pytest.approx(20.0, rel=1e-3)


@pytest.mark.slow
def test_form_factors_are_normalized_to_unity_at_zero_momentum(ca40_density):
    ff = compute(ca40_density)
    assert ff.charge_ff[0] == pytest.approx(1.0, rel=1e-6)
    assert ff.weak_ff[0] == pytest.approx(1.0, rel=1e-6)
