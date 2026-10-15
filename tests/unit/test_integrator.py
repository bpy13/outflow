"""Cash-Karp Runge-Kutta step and the paper code's error measure and step control."""

import math

import pytest

from outflow import integrator


def test_cash_karp_exponential():
    """Adaptive integration of y' = y from 0 to 1 reaches e to better than 1e-10."""
    tolerance = 1e-12
    r, y, h = 0.0, (1.0,), 1e-3
    while r < 1.0:
        h = min(h, 1.0 - r)
        slope = (y[0],)
        new_y, error = integrator.cash_karp_step(lambda _r, state: (state[0],), h, r, y, slope)
        delta = integrator.relative_error(new_y, error)
        if delta <= tolerance:
            r, y = r + h, new_y
            h = integrator.next_step_size(h, delta, tolerance, accepted=True)
        else:
            h = integrator.next_step_size(h, delta, tolerance, accepted=False)
    assert abs(y[0] - math.e) < 1e-10


def test_cash_karp_single_step_is_fourth_order():
    """Halving h reduces the one-step error by about 2^5 (local error of a 4th-order method)."""

    def derivative(_r, state):
        return (math.cos(_r),)

    errors = []
    for h in (0.1, 0.05):
        new_y, _ = integrator.cash_karp_step(derivative, h, 0.0, (0.0,), (1.0,))
        errors.append(abs(new_y[0] - math.sin(h)))
    assert errors[0] / errors[1] == pytest.approx(32, rel=0.2)


@pytest.mark.parametrize(
    "values, errors",
    [
        ((0.0, 2.0, 4.0), (1.0, 2e-12, 4e-12)),
        ((2.0, 0.0, 4.0), (2e-12, 1.0, 4e-12)),
        ((2.0, 4.0, 0.0), (2e-12, 4e-12, 1.0)),
    ],
)
def test_error_norm_skips_zero_components_anywhere(values, errors):
    assert integrator.relative_error(values, errors) == pytest.approx(1e-12)


def test_error_norm_is_largest_relative_error():
    assert integrator.relative_error((1.0, -10.0), (1e-3, 1.0)) == pytest.approx(0.1)


def test_step_size_control_matches_paper_code():
    tolerance = 1e-12
    assert integrator.next_step_size(1.0, 1e-13, tolerance, accepted=True) == pytest.approx(
        0.9 * 10**0.2
    )
    assert integrator.next_step_size(1.0, 1e-11, tolerance, accepted=False) == pytest.approx(
        0.9 * 0.1**0.25
    )
