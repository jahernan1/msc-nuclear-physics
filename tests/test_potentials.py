import numpy as np

from nucth.potentials import POTENTIAL_HEADER, convert, load_raw, save_potential


def _raw_row(r, phi, w, b, a):
    return [r, phi, w, b, a]


def test_convert_applies_the_meson_field_to_nucleon_potential_mapping():
    # raw columns: r, Phi, W, B, A
    raw = np.array([_raw_row(0.0, 400.0, 350.0, -4.0, 10.0)])

    out = convert(raw)

    r, sn, vn, sp, vp = out[0]
    assert r == 0.0
    assert sn == -400.0                    # S(n) = -Phi
    assert sp == -400.0                    # S(p) = -Phi
    assert vn == 350.0 - 0.5 * -4.0        # V(n) = W - B/2  = 352.0
    assert vp == 350.0 + 0.5 * -4.0 + 10.0 # V(p) = W + B/2 + A = 358.0


def test_convert_preserves_shape_and_radial_axis():
    rng = np.random.default_rng(0)
    raw = rng.normal(size=(50, 5))
    raw[:, 0] = np.linspace(0.0, 12.0, 50)

    out = convert(raw)

    assert out.shape == (50, 5)
    np.testing.assert_array_equal(out[:, 0], raw[:, 0])


def test_isoscalar_limit_makes_proton_and_neutron_vector_potentials_differ_only_by_coulomb():
    """With B = 0 (no rho field), V(p) - V(n) is exactly the Coulomb term A."""
    raw = np.array([_raw_row(1.0, 300.0, 340.0, 0.0, 7.5)])

    _, _, vn, _, vp = convert(raw)[0]

    assert vp - vn == 7.5


def test_round_trip_through_disk(tmp_path):
    raw = np.array([
        _raw_row(0.000, 424.906, 346.759, -3.82645, 10.5182),
        _raw_row(0.025, 424.906, 346.759, -3.82645, 10.5182),
        _raw_row(0.050, 424.592, 346.265, -3.82123, 10.5176),
    ])
    out = convert(raw)

    path = tmp_path / "Test_model_potential.dat"
    save_potential(path, out)

    assert path.is_file()
    assert path.read_text().splitlines()[0].lstrip("# ").strip() == POTENTIAL_HEADER
    np.testing.assert_allclose(np.loadtxt(path), out, rtol=1e-12)


def test_load_raw_reads_a_five_column_file(tmp_path):
    path = tmp_path / "raw_Test_model.dat"
    path.write_text("0 424.906 346.759 -3.82645 10.5182\n0.025 424.9 346.7 -3.8 10.5\n")

    raw = load_raw(path)

    assert raw.shape == (2, 5)
    assert raw[0, 1] == 424.906
