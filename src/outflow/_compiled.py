"""Compiled copy (Numba) of the trial integration used by the central-density search.

This module is optional: it is used only when Numba is installed, and only for winds without
cooling. It repeats, line by line, the arithmetic of :mod:`outflow.physics`,
:mod:`outflow.equations`, :mod:`outflow.integrator` and ``boundary._sonic_trial``; those
modules remain the readable reference, and the tests check that this copy gives bit-for-bit
the same floats.

Two details keep the results identical to the Python code:

- **Powers.** Python evaluates ``x**2`` by calling the C library's ``pow``, which is not
  always the exactly rounded ``x*x``. A compiler that sees a constant exponent replaces
  ``x**2`` (even ``x**2.0``) by ``x*x``. So every exponent here is read from the parameter
  array at run time (``p[TWO]``, ``p[THREE]``, ...): the compiler cannot see its value and calls
  ``pow``, as Python does.
- **Sums.** Since Python 3.12, ``sum()`` adds floats with Neumaier's compensated summation;
  ``_compensated_sum`` repeats that algorithm.
"""

import math

import numpy as np
from numba import njit

from outflow.constants import FOUR_PI, GAMMA_CR, GAMMA_GAS, GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT
from outflow.integrator import _A, _B, _C4, _ERROR_WEIGHTS
from outflow.physics import NFW_SERIES_LIMIT, NFW_SERIES_TERMS

# Cash-Karp tableau as arrays (row `stage` of _STAGE_WEIGHTS holds _B[stage], zero-padded).
_STAGE_NODES = np.array(_A)
_STAGE_WEIGHTS = np.zeros((6, 5))
for _stage, _row in enumerate(_B):
    _STAGE_WEIGHTS[_stage, : len(_row)] = _row
_FOURTH_ORDER_WEIGHTS = np.array(_C4)
_ERROR_WEIGHTS_ARRAY = np.array(_ERROR_WEIGHTS)

# Positions of the inputs in the packed parameter array (see pack_parameters).
(
    STARBURST_RADIUS,
    MASS_RATE,
    TERMINAL_VELOCITY,
    MASS_INJECTION_RATE,
    THERMAL_INJECTION_RATE,
    CR_INJECTION_RATE,
    OPACITY,
    MAGNETIC_FIELD,
    NFW_DENSITY,
    NFW_SCALE_RADIUS,
    # Exponents, passed at run time so that powers call pow() as Python does (see above).
    ONE,
    TWO,
    THREE,
    GROWTH_EXPONENT,
    SHRINK_EXPONENT,
) = range(15)


def pack_parameters(params) -> np.ndarray:
    """The WindParameters values the equations use, and the exponents, as a float array."""
    return np.array(
        [
            params.starburst_radius,
            params.mass_rate,
            params.terminal_velocity,
            params.mass_injection_rate,
            params.thermal_energy_injection_rate,
            params.cr_energy_injection_rate,
            params.opacity,
            params.magnetic_field,
            params.nfw_characteristic_density,
            params.nfw_scale_radius,
            1.0,
            2.0,
            3.0,
            0.2,  # step growth exponent of integrator.next_step_size
            0.25,  # step shrink exponent of integrator.next_step_size
        ]
    )


@njit(cache=True)
def _compensated_sum(values, count):
    """sum(values[:count]) exactly as Python 3.12's built-in sum() adds floats."""
    total = 0.0
    compensation = 0.0
    for i in range(count):
        x = values[i]
        t = total + x
        if abs(total) >= abs(x):
            compensation += (total - t) + x
        else:
            compensation += (x - t) + total
        total = t
    if compensation != 0.0 and math.isfinite(compensation):
        total += compensation
    return total


# --- physics.py and equations.py --------------------------------------------------------------


@njit(cache=True)
def _density(r, v, p, searching):
    starburst_radius = p[STARBURST_RADIUS]
    if searching or r < starburst_radius:
        return p[MASS_RATE] * r / (FOUR_PI * starburst_radius ** p[THREE] * v)
    return p[MASS_RATE] / (FOUR_PI * r ** p[TWO] * v)


@njit(cache=True)
def _alfven_speed(r, density, p, searching):
    starburst_radius = p[STARBURST_RADIUS]
    if searching or r <= starburst_radius:
        field = p[MAGNETIC_FIELD] * r / starburst_radius
    else:
        field = p[MAGNETIC_FIELD] * starburst_radius ** p[TWO] / r ** p[TWO]
    return field / math.sqrt(FOUR_PI * density)


@njit(cache=True)
def _sound_speed_squared_from(P, Pc, rho, v, v_alfven):
    cr_factor = (
        (2.0 * v + v_alfven) * (v - (GAMMA_GAS - 1.0) * v_alfven) / (2.0 * v * (v + v_alfven))
    )
    return GAMMA_GAS * P / rho + GAMMA_CR * Pc * cr_factor / rho


@njit(cache=True)
def _sound_speed_squared(r, v, P, Pc, p, searching):
    rho = _density(r, v, p, searching)
    return _sound_speed_squared_from(P, Pc, rho, v, _alfven_speed(r, rho, p, searching))


@njit(cache=True)
def _luminosity(v, P, Pc, rho, cr_velocity_ratio, p):
    return p[MASS_RATE] * (
        p[TERMINAL_VELOCITY] ** p[TWO] / 2.0
        - v ** p[TWO] / 2.0
        - GAMMA_GAS * P / ((GAMMA_GAS - 1.0) * rho)
        - GAMMA_CR * cr_velocity_ratio * Pc / (rho * (GAMMA_CR - 1.0))
    )


@njit(cache=True)
def _radiative_force(r, luminosity, p, inner):
    opacity, starburst_radius = p[OPACITY], p[STARBURST_RADIUS]
    if inner:
        return opacity * luminosity * r / (FOUR_PI * SPEED_OF_LIGHT * starburst_radius ** p[THREE])
    return opacity * luminosity / (FOUR_PI * SPEED_OF_LIGHT * r ** p[TWO])


@njit(cache=True)
def _radiative_power_factor(r, luminosity, p):
    opacity, starburst_radius = p[OPACITY], p[STARBURST_RADIUS]
    if r == 0:
        return 3.0 * opacity * luminosity / (FOUR_PI * starburst_radius ** p[TWO] * SPEED_OF_LIGHT)
    if r == starburst_radius:
        return (
            3.0
            * opacity
            * luminosity
            / (8.0 * math.pi * starburst_radius ** p[TWO] * SPEED_OF_LIGHT)
        )
    shape = 2.0 * starburst_radius * r - (r ** p[TWO] - starburst_radius ** p[TWO]) * math.log(
        (r + starburst_radius) / abs(r - starburst_radius)
    )
    return (
        3.0
        * opacity
        * luminosity
        * shape
        / (16.0 * math.pi * starburst_radius ** p[THREE] * SPEED_OF_LIGHT * r)
    )


@njit(cache=True)
def _nfw_gravity(r, p):
    density, scale_radius = p[NFW_DENSITY], p[NFW_SCALE_RADIUS]
    if r > scale_radius * NFW_SERIES_LIMIT:
        shape = math.log((scale_radius + r) / scale_radius) - r / (scale_radius + r)
        return (
            FOUR_PI
            * GRAVITATIONAL_CONSTANT
            * density
            * scale_radius ** p[THREE]
            * shape
            / r ** p[TWO]
        )
    x = -r / scale_radius
    terms = np.empty(NFW_SERIES_TERMS)
    for i in range(NFW_SERIES_TERMS):
        terms[i] = (i + 1) * x ** (i * p[ONE]) / (i + 2)  # x**i, exponent opaque to compiler
    series = _compensated_sum(terms, NFW_SERIES_TERMS)
    return FOUR_PI * GRAVITATIONAL_CONSTANT * density * scale_radius * series


@njit(cache=True)
def _derivatives(r, v, P, Pc, p, searching, out):
    """WindEquations.derivatives without cooling (r > 0), written into ``out``."""
    gamma, gamma_cr = GAMMA_GAS, GAMMA_CR
    rho = _density(r, v, p, searching)
    v_alfven = _alfven_speed(r, rho, p, searching)
    cooling = 0.0
    cs2 = _sound_speed_squared_from(P, Pc, rho, v, v_alfven)
    luminosity = _luminosity(v, P, Pc, rho, (v + v_alfven) / v, p)
    power_factor = _radiative_power_factor(r, luminosity, p)
    gravity = _nfw_gravity(r, p)
    subsonic = v ** p[TWO] <= cs2
    force = _radiative_force(r, luminosity, p, subsonic)
    external = (
        (gamma - 1.0) * (cooling / rho - power_factor * v) + gamma * force * v - gravity * v
    ) / (v ** p[TWO] - cs2)
    cr_factor = gamma_cr * Pc * (2.0 * v + v_alfven) / (2.0 * (v + v_alfven))
    gas_factor = gamma * P - gamma_cr * Pc * (gamma - 1.0) * v_alfven * (2.0 * v + v_alfven) / (
        2.0 * v * (v + v_alfven)
    )
    if subsonic:
        source_cr = -3.0 * gamma_cr * Pc / r + (gamma_cr - 1.0) * p[CR_INJECTION_RATE] / (
            v + v_alfven
        )
        source_gas = (
            -(gamma - 1.0) * v_alfven * source_cr / v
            + (gamma - 1.0) * p[THERMAL_INJECTION_RATE] / v
            - 3.0 * gamma * P / r
            + 3.0 * (gamma - 1.0) * rho * v ** p[TWO] / (2.0 * r)
        )
        dv_cr = cs2 * v / (r * (cs2 - v ** p[TWO])) + v * (
            source_gas + source_cr + p[MASS_INJECTION_RATE] * v
        ) / (rho * (cs2 - v ** p[TWO]))
        dv_dr = dv_cr + external
        log_density_gradient = 1.0 / r - dv_dr / v
        dpc_dr = log_density_gradient * cr_factor + source_cr
        dp_dr = log_density_gradient * gas_factor + source_gas
    else:
        dv_dr = 2.0 * cs2 * v / (r * (v ** p[TWO] - cs2)) + external
        log_density_gradient = -2.0 / r - dv_dr / v
        dpc_dr = log_density_gradient * cr_factor
        dp_dr = log_density_gradient * gas_factor
    dp_dr = dp_dr + (gamma - 1.0) * (rho * power_factor - rho * force - cooling / v)
    out[0], out[1], out[2] = dv_dr, dp_dr, dpc_dr


# --- integrator.py ----------------------------------------------------------------------------


@njit(cache=True)
def _cash_karp_step(h, r, y, first_slope, p, searching, new_y, error, stages, work, slope):
    """integrator.cash_karp_step for the 3-component wind state; results in new_y, error."""
    for i in range(3):
        stages[i, 0] = h * first_slope[i]
    state = np.empty(3)
    for stage in range(1, 6):
        for i in range(3):
            for j in range(stage):
                work[j] = _STAGE_WEIGHTS[stage, j] * stages[i, j]
            state[i] = y[i] + _compensated_sum(work, stage)
        _derivatives(r + _STAGE_NODES[stage] * h, state[0], state[1], state[2], p, searching, slope)
        for i in range(3):
            stages[i, stage] = h * slope[i]
    for i in range(3):
        for j in range(6):
            work[j] = _FOURTH_ORDER_WEIGHTS[j] * stages[i, j]
        new_y[i] = y[i] + _compensated_sum(work, 6)
        for j in range(6):
            work[j] = _ERROR_WEIGHTS_ARRAY[j] * stages[i, j]
        error[i] = _compensated_sum(work, 6)


@njit(cache=True)
def _relative_error(values, errors):
    largest = 0.0
    for i in range(3):
        if values[i] != 0.0:
            ratio = abs(errors[i] / values[i])
            if ratio > largest:
                largest = ratio
    return largest


@njit(cache=True)
def _next_step_size(h, delta, tolerance, accepted, growth_when_error_is_zero, p):
    if delta == 0.0:
        return h * growth_when_error_is_zero
    exponent = p[GROWTH_EXPONENT] if accepted else p[SHRINK_EXPONENT]
    return 0.9 * h * (tolerance / delta) ** exponent


# --- boundary._sonic_trial --------------------------------------------------------------------


@njit(cache=True)
def sonic_trial(
    p,
    pressure,
    cr_pressure,
    velocity_gradient,
    central_cs2,
    initial_step,
    step_tolerance,
    first_step_tolerance,
    smallest_step_fraction,
    growth_when_error_is_zero,
):
    """boundary._sonic_trial: returns (stopping radius, dv/dr there, or inf if it failed)."""
    starburst_radius = p[STARBURST_RADIUS]
    r = 0.0
    y = np.array([0.0, pressure, cr_pressure])
    slope = np.array([velocity_gradient, 0.0, 0.0])
    h = initial_step
    cs2 = central_cs2
    gradient = velocity_gradient
    new_y, error, scratch_slope = np.empty(3), np.empty(3), np.empty(3)
    stages, work = np.empty((3, 6)), np.empty(6)
    while cs2 >= y[0] ** p[TWO] and gradient > 0:
        _cash_karp_step(h, r, y, slope, p, True, new_y, error, stages, work, scratch_slope)
        delta = _relative_error(new_y, error)
        if h < starburst_radius * smallest_step_fraction:
            return r, math.inf
        if delta <= (first_step_tolerance if r == 0 else step_tolerance):
            r = r + h
            y[:] = new_y
            cs2 = _sound_speed_squared(r, y[0], y[1], y[2], p, True)
            if r + h == r and y[0] ** p[TWO] < 0.9 * cs2:
                return r, math.inf
            _derivatives(r, y[0], y[1], y[2], p, True, slope)
            gradient = slope[0]
            h = _next_step_size(h, delta, step_tolerance, True, growth_when_error_is_zero, p)
        else:
            h = _next_step_size(h, delta, step_tolerance, False, growth_when_error_is_zero, p)
    return r, gradient
