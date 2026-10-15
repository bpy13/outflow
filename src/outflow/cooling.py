"""Experimental radiative cooling (paper §2.1.2, eqs. 10 and 14). Off by default.

The paper's published runs (Figs. 6-9) never used cooling, and this option is not validated
against the paper. It uses the Cloudy cooling function Lambda(T) (solar metallicity), shipped
as ``outflow/data/cloudy_cooling.csv`` and loaded on first use, interpolated with the same
cubic spline as the paper code (``plotcool.py``: ``scipy.interpolate.splrep``).

With cooling, the central temperature T_0 sets the central density through paper eq. 14,

    rho_0 = mu sqrt((Q_th - q c_s0^2 / (gamma_g - 1)) / Lambda(T_0)),
    with c_s0^2 = gamma_g k_B T_0 / mu,

and T_0 is searched for as in ``legacy/outflow.py`` (``initial_condition``): find the T_0 that
maximises the sonic radius; if that radius exceeds r_sb, solve for the T_0 that puts the sonic
point at r_sb. Otherwise no cooling solution exists and the wind is solved without cooling.
(The paper code's own search, ``Tsolve``, has a branch that returns nothing.)
"""

import functools
import math
from importlib import resources

import numpy as np
from scipy import interpolate, optimize

from outflow import boundary, physics
from outflow.constants import BOLTZMANN, GAMMA_GAS, MEAN_PARTICLE_MASS
from outflow.equations import WindEquations
from outflow.parameters import WindParameters


@functools.lru_cache(maxsize=1)
def load_cooling_table() -> tuple[np.ndarray, np.ndarray]:
    """The Cloudy table (T in K, Lambda in erg cm^3 s^-1), read from the package data once."""
    resource = resources.files("outflow") / "data" / "cloudy_cooling.csv"
    with resources.as_file(resource) as path:
        table = np.loadtxt(path, delimiter=",", comments="#", skiprows=2)
    return table[:, 0], table[:, 1]


@functools.lru_cache(maxsize=1)
def _spline():
    temperature, rate = load_cooling_table()
    return interpolate.splrep(temperature, rate)


def cooling_function(temperature):
    """Lambda(T) in erg cm^3 s^-1 (cubic spline through the Cloudy table), scalar or array."""
    value = interpolate.splev(temperature, _spline())
    return float(value) if np.ndim(temperature) == 0 else value


def _scalar_cooling_function(temperature: float) -> float:
    return float(interpolate.splev(temperature, _spline()))


def maximum_central_temperature(params: WindParameters) -> float:
    """T_0 at which the central density of eq. 14 drops to zero (no energy left to radiate).

    T_max = (gamma_g - 1)/gamma_g mu/k_B Q_th/q, the central temperature without cooling.
    """
    return (
        (GAMMA_GAS - 1.0)
        / GAMMA_GAS
        * MEAN_PARTICLE_MASS
        / BOLTZMANN
        * params.thermal_energy_injection_rate
        / params.mass_injection_rate
    )


def cooling_central_density(params: WindParameters, temperature: float) -> float:
    """rho_0 from paper eq. 14 for a central temperature T_0, or 0 if there is none."""
    sound_speed_squared = GAMMA_GAS * BOLTZMANN * temperature / MEAN_PARTICLE_MASS
    numerator = (
        params.thermal_energy_injection_rate
        - params.mass_injection_rate * sound_speed_squared / (GAMMA_GAS - 1.0)
    )
    rate = _scalar_cooling_function(temperature)
    if numerator <= 0 or rate <= 0:
        return 0.0
    return MEAN_PARTICLE_MASS * math.sqrt(numerator / rate)


def central_state_for_temperature(params: WindParameters, temperature: float):
    """Central values for T_0 with cooling (eq. 14 for rho_0; eq. 44 for P_C,0), or None."""
    density = cooling_central_density(params, temperature)
    if not density > 0 or not math.isfinite(density):
        return None
    velocity_gradient, alfven_gradient = physics.central_velocity_gradients(
        density, params.mass_injection_rate, params.magnetic_field, params.starburst_radius
    )
    pressure = density * BOLTZMANN * temperature / MEAN_PARTICLE_MASS
    # Same P_C,0 as without cooling (eq. 44): it depends only on rho_0 and the gradients.
    cr_pressure = boundary.central_state(params, density).cr_pressure
    return boundary.CentralState(density, pressure, cr_pressure, velocity_gradient, alfven_gradient)


def cooling_central_state(params: WindParameters):
    """Central state with cooling and the cooling function, or (None, None) if there is none.

    Follows ``legacy/outflow.py``'s ``initial_condition``: a bounded maximisation of the
    sonic radius over T_0 in (0, T_max), then a bracketed root find for r_s = r_sb.
    """
    equations = WindEquations(params, _scalar_cooling_function)
    starburst_radius = params.starburst_radius
    t_max = maximum_central_temperature(params)

    def sonic_radius(temperature: float) -> float:
        state = central_state_for_temperature(params, temperature)
        if state is None:
            return 0.0
        return boundary.sonic_trial(equations, state).radius

    peak = optimize.minimize_scalar(
        lambda t: -sonic_radius(t), bounds=(0.0, t_max), method="bounded"
    )
    if not -peak.fun > starburst_radius:
        return None, None
    try:
        temperature = optimize.brentq(lambda t: sonic_radius(t) - starburst_radius, peak.x, t_max)
    except ValueError:
        return None, None
    state = central_state_for_temperature(params, temperature)
    if state is None:
        return None, None
    return state, _scalar_cooling_function
