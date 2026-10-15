"""Thermal-only limit: the solver against the analytic CC85 solution (paper eqs. 7-9)."""

import numpy as np
import pytest

from outflow import WindParameters, solve_wind
from outflow.cc85 import cc85_dimensionless_velocity, cc85_velocity_profile

THERMAL_ONLY = WindParameters.paper_default().replace(
    thermal_fraction=1.0, cr_fraction=0.0, opacity=0.0, magnetic_field=0.0, halo_mass=0.0
)


def test_sonic_point_at_starburst_radius():
    # Eq. 9: v^2(r_sb) = (gamma - 1)/(gamma + 1) v_inf^2, i.e. u = 1/2 for gamma = 5/3.
    assert cc85_dimensionless_velocity(1.0) == pytest.approx(0.5, rel=1e-14)
    assert cc85_dimensionless_velocity(0.0) == 0.0


def test_analytic_solution_satisfies_eq8_and_limits():
    x = np.array([1e-3, 0.3, 0.9, 1.1, 3.0, 100.0])
    u = cc85_dimensionless_velocity(x)
    assert np.all(np.diff(u) > 0) and np.all(u < 1)
    # Deep inside, u grows linearly: x = A u with A = 2 * 4.5^(9/14).
    assert u[0] == pytest.approx(1e-3 / (2 * 4.5 ** (9 / 14)), rel=1e-5)
    # Far out, u -> 1.
    assert u[-1] > 0.99


def test_solver_matches_cc85():
    r_sb = THERMAL_ONLY.starburst_radius
    profile = solve_wind(THERMAL_ONLY, outer_radius=10 * r_sb)
    mask = profile.r >= 0.2 * r_sb
    analytic = cc85_velocity_profile(profile.r[mask], r_sb, THERMAL_ONLY.terminal_velocity)
    difference = np.max(np.abs(profile.v[mask] / analytic - 1))
    print(f"CC85 max relative difference in v: {difference:.3e}")
    assert difference <= 1e-4
