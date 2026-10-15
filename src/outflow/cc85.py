"""Analytic thermally driven wind of Chevalier & Clegg (1985), paper §2.1.1 (eqs. 7-9).

With x = r / r_sb and u = v / v_inf, the solution is given implicitly by paper eq. 8:

    x = A u (1 + (5 gamma + 1)/(gamma - 1) u^2)^(-(3 gamma + 1)/(5 gamma + 1))   for x <= 1,
    x = B u^(-1/2) (1 - u^2)^(-1/(2 (gamma - 1)))                                   for x >= 1,

with A and B fixed by the sonic condition at r_sb, eq. 9: u(1)^2 = (gamma - 1)/(gamma + 1).
The inverse u(x) has no closed form, so it is found numerically (bracketed root finding).
"""

import math

import numpy as np
from scipy.optimize import brentq

from outflow.constants import GAMMA_GAS


def _inner_x(u: float, gamma: float) -> float:
    """x(u) / A for the subsonic branch of eq. 8."""
    return u * (1 + (5 * gamma + 1) / (gamma - 1) * u**2) ** (-(3 * gamma + 1) / (5 * gamma + 1))


def _outer_x(u: float, gamma: float) -> float:
    """x(u) / B for the supersonic branch of eq. 8."""
    return u ** (-0.5) * (1 - u**2) ** (-1 / (2 * (gamma - 1)))


def sonic_velocity_ratio(gamma: float = GAMMA_GAS) -> float:
    """u at r_sb from eq. 9: sqrt((gamma - 1)/(gamma + 1)) (1/2 for gamma = 5/3)."""
    return math.sqrt((gamma - 1) / (gamma + 1))


def _dimensionless_velocity(x: float, gamma: float) -> float:
    u_sonic = sonic_velocity_ratio(gamma)
    if x == 0:
        return 0.0
    if x == 1:
        return u_sonic
    if x < 1:
        a = 1 / _inner_x(u_sonic, gamma)
        return brentq(lambda u: a * _inner_x(u, gamma) - x, 0.0, u_sonic, xtol=1e-300, rtol=1e-15)
    b = 1 / _outer_x(u_sonic, gamma)
    # Solve in log form: x grows without bound as u -> 1.
    upper = 1.0 - 1e-15
    return brentq(
        lambda u: math.log(b * _outer_x(u, gamma)) - math.log(x),
        u_sonic,
        upper,
        xtol=1e-300,
        rtol=1e-15,
    )


def cc85_dimensionless_velocity(x, gamma: float = GAMMA_GAS):
    """u = v / v_inf of the CC85 wind at x = r / r_sb (scalar or array), paper eq. 8."""
    if np.ndim(x) == 0:
        return _dimensionless_velocity(float(x), gamma)
    return np.array([_dimensionless_velocity(float(xi), gamma) for xi in np.ravel(x)]).reshape(
        np.shape(x)
    )


def cc85_velocity_profile(radius, starburst_radius: float, terminal_velocity: float):
    """CC85 wind velocity (cm/s) at ``radius`` (cm, scalar or array), paper eqs. 7-9.

    ``terminal_velocity`` is v_inf = sqrt(2 E-dot / M-dot) (see
    :attr:`WindParameters.terminal_velocity`).
    """
    return terminal_velocity * cc85_dimensionless_velocity(
        np.asarray(radius, dtype=float) / starburst_radius
    )
