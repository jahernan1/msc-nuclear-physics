import numpy as np
import pytest
from scipy.integrate import simpson

from nucth.transforms import fourier, fourier_back

EPS = 1.0e-6


def test_forward_does_not_mutate_its_inputs():
    q = np.linspace(0.0, 5.0, 64)
    r = np.linspace(0.0, 10.0, 256)
    pr = np.exp(-(r**2))
    q_before, r_before = q.copy(), r.copy()

    fourier(q, r, pr)

    np.testing.assert_array_equal(q, q_before)
    np.testing.assert_array_equal(r, r_before)


def test_backward_does_not_mutate_its_inputs():
    q = np.linspace(0.0, 5.0, 64)
    r = np.linspace(0.0, 10.0, 256)
    pq = np.exp(-(q**2))
    q_before, r_before = q.copy(), r.copy()

    fourier_back(q, r, pq)

    np.testing.assert_array_equal(q, q_before)
    np.testing.assert_array_equal(r, r_before)


def test_forward_at_zero_momentum_gives_the_normalization():
    """F(q -> 0) equals the volume integral 4*pi*int r^2 rho(r) dr."""
    r = np.linspace(0.0, 12.0, 2000)
    rho = np.exp(-(r**2) / 2.0)
    q = np.linspace(0.0, 1.0, 8)

    norm = 4.0 * np.pi * simpson(r**2 * rho, x=r)
    assert fourier(q, r, rho)[0] == pytest.approx(norm, rel=1e-4)


def test_round_trip_recovers_a_gaussian_density():
    r = np.linspace(0.0, 12.0, 1200)
    q = np.linspace(0.0, 25.0, 1200)
    rho = np.exp(-(r**2) / 2.0)

    recovered = fourier_back(q, r, fourier(q, r, rho))

    # Compare on the interior; both endpoints carry the 1/r and 1/q regulators.
    interior = slice(50, 600)
    np.testing.assert_allclose(recovered[interior], rho[interior], rtol=2e-2, atol=2e-3)


def test_zero_endpoint_is_regulated_not_infinite():
    r = np.linspace(0.0, 10.0, 500)
    q = np.linspace(0.0, 10.0, 500)
    rho = np.exp(-r)

    assert np.all(np.isfinite(fourier(q, r, rho)))
    assert np.all(np.isfinite(fourier_back(q, r, rho)))
