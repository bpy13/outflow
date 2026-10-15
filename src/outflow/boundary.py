"""Central boundary conditions and the search for the central density rho_0.

Ported from the paper code (``CRline.py``: ``rsonic()``, ``gss()`` and the start of
``CRline()``), in cgs. Paper §2.3.2 and §2.4 (eqs. 44-45): the central pressures follow from
rho_0, and rho_0 is varied until the sonic radius r_s (where v = c*) equals r_sb.

The search, as in the paper code:

1. Start at rho_0 = 0.1 Msun/pc^3 and multiply by 10 until a trial solution either turns over
   (dv/dr <= 0 while still subsonic: a "breeze", rho_0 too high) or reaches r_sb.
2. Bisect between 0 and that rho_0 for the critical value that separates winds that turn
   supersonic from breezes (the "peak"); this gives the largest possible sonic radius.
3. If that sonic radius lies beyond r_sb, find the rho_0 between 0 and the peak value whose
   sonic radius is r_sb. The paper code bisects; by default this port uses Brent's method,
   which reaches the same precision with fewer trials (``method="bisection"`` restores the
   paper code's search).

Trial solutions use the inner forms of density and Alfven speed at every radius
(``searching=True``, the paper code's tag ``'find0'``).
"""

import functools
import math
from dataclasses import dataclass

from outflow import physics
from outflow.constants import GAMMA_CR, GAMMA_GAS, PARSEC, SOLAR_MASS
from outflow.equations import WindEquations
from outflow.errors import NoSolutionError
from outflow.integrator import (
    _GROWTH_WHEN_ERROR_IS_ZERO,
    cash_karp_step,
    next_step_size,
    relative_error,
)
from outflow.parameters import WindParameters

# Numerical settings of the paper code, converted from its natural units (pc, Msun).
INITIAL_STEP = 1e-5 * PARSEC  # first step size, cm
STEP_TOLERANCE = 1e-12  # relative error per step
SEARCH_TOLERANCE = 1e-12  # bisection stops when the bracket is this fraction of its start
FIRST_DENSITY_TRIAL = 0.1 * SOLAR_MASS / PARSEC**3  # first rho_0 trial, g/cm^3
# A trial gives up when the step falls below this fraction of r_sb (paper code: 1e-50).
SMALLEST_STEP_FRACTION = 1e-50
# Safety limit (not in the paper code): stop multiplying rho_0 by 10 after this many tries.
MAX_DENSITY_INCREASES = 40
# Brent search (not in the paper code): halve rho_0 at most this many times looking for a sonic
# radius inside r_sb before falling back to the bisection. 40 halvings reach the bisection's
# own precision (2^-40 ~ 1e-12).
MAX_BRACKET_HALVINGS = 40
# Fallback for the first step from r = 0 (not in the paper code). Near r = 0 the equations are
# unstable within a step: round-off in the Runge-Kutta stages grows by a factor of ~40 per
# stage, whatever the step size, so the error estimate of the first step has a floor of about
# 1e-12 to 5e-12. Usually the floor is below STEP_TOLERANCE; for some strongly magnetised
# CR-driven winds (Fig. 7, B_0 = 100 uG) it is not, and the step size collapses. Only then is
# the integration restarted, accepting the first step at this looser tolerance.
FIRST_STEP_FALLBACK_TOLERANCE = 100 * STEP_TOLERANCE

# Use the compiled copy of the trial integration (outflow._compiled) when Numba is installed.
# It gives the same results, about 30x faster. Set to False to run the plain Python version.
USE_COMPILED_TRIALS = True

# Marker for a trial whose step size collapsed before it became sonic or turned over (the
# paper code returns the string 'RKfail'). It is treated as larger than any gradient, which is
# how Python 2 compares a string with a number; see docs/open-items.md.
TRIAL_FAILED = math.inf


@dataclass(frozen=True)
class CentralState:
    """Values at r = 0 for a given central density (paper eqs. 44-45)."""

    density: float  # rho_0, g/cm^3
    pressure: float  # P_0, erg/cm^3
    cr_pressure: float  # P_C,0, erg/cm^3
    velocity_gradient: float  # dv/dr at r = 0, 1/s
    alfven_gradient: float  # dv_A/dr at r = 0, 1/s


@dataclass(frozen=True)
class SonicTrial:
    """Result of one trial integration from r = 0 (the paper code's ``rsonic()``).

    ``radius`` is where the trial stopped: the sonic radius if it turned supersonic, or where
    it turned over. ``velocity_gradient`` is dv/dr there: positive for a wind that went
    supersonic, zero or negative for a breeze, ``TRIAL_FAILED`` if the integration broke down.
    """

    radius: float
    velocity_gradient: float


def central_state(params: WindParameters, central_density: float) -> CentralState:
    """Central pressures for a given rho_0 without cooling, paper eqs. 44-45.

    P_0 = (gamma_g - 1)/gamma_g rho_0 Q_th/q and
    P_C,0 = (gamma_C - 1)/gamma_C rho_0 Q_CR/q * v/(v + v_A), with v/(v + v_A) at r = 0 from the
    L'Hopital gradients. Q_th/q = f_th v_inf^2 / 2 and Q_CR/q = f_CR v_inf^2 / 2.
    """
    v_inf = params.terminal_velocity
    velocity_gradient, alfven_gradient = physics.central_velocity_gradients(
        central_density, params.mass_injection_rate, params.magnetic_field, params.starburst_radius
    )
    pressure = (
        params.thermal_fraction * central_density * (GAMMA_GAS - 1.0) * v_inf**2 / (2.0 * GAMMA_GAS)
    )
    cr_pressure = (
        params.cr_fraction
        * central_density
        * (GAMMA_CR - 1.0)
        * v_inf**2
        * velocity_gradient
        / (2.0 * GAMMA_CR * (velocity_gradient + alfven_gradient))
    )
    return CentralState(central_density, pressure, cr_pressure, velocity_gradient, alfven_gradient)


def central_sound_speed_squared(state: CentralState) -> float:
    """c*^2 at r = 0 (eq. 35 with the L'Hopital ratio v_A / v)."""
    return physics.effective_sound_speed_squared(
        state.pressure,
        state.cr_pressure,
        state.density,
        state.velocity_gradient,
        state.alfven_gradient,
    )


def sonic_trial(equations: WindEquations, state: CentralState) -> SonicTrial:
    """Integrate from r = 0 until the flow turns sonic or turns over (paper code ``rsonic()``).

    If the integration breaks down on its very first step (see
    ``FIRST_STEP_FALLBACK_TOLERANCE``), it is repeated with the looser first-step tolerance.
    """
    trial = _sonic_trial(equations, state, STEP_TOLERANCE)
    if trial.velocity_gradient == TRIAL_FAILED and trial.radius == 0:
        trial = _sonic_trial(equations, state, FIRST_STEP_FALLBACK_TOLERANCE)
    return trial


def _sonic_trial(
    equations: WindEquations, state: CentralState, first_step_tolerance: float
) -> SonicTrial:
    """One trial: the compiled copy if Numba is installed and there is no cooling, else Python.

    Both give bit-for-bit the same result (tests/unit/test_compiled.py).
    """
    if USE_COMPILED_TRIALS and equations.cooling_function is None:
        kernel = _compiled_kernel()
        if kernel is not None:
            radius, gradient = kernel.sonic_trial(
                kernel.pack_parameters(equations.params),
                state.pressure,
                state.cr_pressure,
                state.velocity_gradient,
                central_sound_speed_squared(state),
                INITIAL_STEP,
                STEP_TOLERANCE,
                first_step_tolerance,
                SMALLEST_STEP_FRACTION,
                _GROWTH_WHEN_ERROR_IS_ZERO,
            )
            return SonicTrial(radius, gradient)
    return _sonic_trial_python(equations, state, first_step_tolerance)


def compiled_trials_available() -> bool:
    """True if the central-density search runs compiled (Numba installed and switched on)."""
    return USE_COMPILED_TRIALS and _compiled_kernel() is not None


@functools.cache
def _compiled_kernel():
    """The Numba-compiled trial integration (outflow._compiled), or None without Numba.

    Imported on first use, not when the package is imported (Numba takes a moment to load).
    The first use after installing compiles it (a few seconds); the compiled code is then
    cached on disk.
    """
    try:
        from outflow import _compiled
    except ImportError:
        return None
    return _compiled


def _sonic_trial_python(
    equations: WindEquations, state: CentralState, first_step_tolerance: float
) -> SonicTrial:
    """The trial integration in plain Python: the readable reference for outflow._compiled."""
    starburst_radius = equations.starburst_radius

    def derivative(r, y):
        return equations.derivatives(r, y[0], y[1], y[2], True)

    r, y = 0.0, (0.0, state.pressure, state.cr_pressure)
    h = INITIAL_STEP
    slope = (state.velocity_gradient, 0.0, 0.0)  # dv/dr at r = 0, as the paper code's du
    cs2 = central_sound_speed_squared(state)
    velocity_gradient = state.velocity_gradient
    while cs2 >= y[0] ** 2 and velocity_gradient > 0:
        new_y, error = cash_karp_step(derivative, h, r, y, slope)
        delta = relative_error(new_y, error)
        if h < starburst_radius * SMALLEST_STEP_FRACTION:
            return SonicTrial(r, TRIAL_FAILED)  # unstable iteration
        if delta <= (first_step_tolerance if r == 0 else STEP_TOLERANCE):
            r, y = r + h, new_y
            cs2 = equations.sound_speed_squared(r, y[0], y[1], y[2], True)
            if r + h == r and y[0] ** 2 < 0.9 * cs2:
                return SonicTrial(r, TRIAL_FAILED)  # step below round-off while subsonic
            slope = derivative(r, y)  # also the first stage of the next step
            velocity_gradient = slope[0]
            h = next_step_size(h, delta, STEP_TOLERANCE, accepted=True)
        else:
            h = next_step_size(h, delta, STEP_TOLERANCE, accepted=False)
    return SonicTrial(r, velocity_gradient)


def _bisect(equations: WindEquations, params, low: float, high: float, find_peak: bool):
    """The paper code's ``gss()`` without cooling: bisection on rho_0.

    ``find_peak``: bracket the critical rho_0 between winds (dv/dr > 0 at the sonic point,
    rho_0 too low) and breezes (dv/dr <= 0, too high). Otherwise: bracket the rho_0 whose sonic
    radius is r_sb (a sonic radius beyond r_sb means rho_0 too high).
    Returns (stopping radius of the last trial, lower end of the final bracket).
    """
    width = high - low
    trial = None
    while high - low > SEARCH_TOLERANCE * width:
        middle = (low + high) / 2.0
        trial = sonic_trial(equations, central_state(params, middle))
        if find_peak:
            indicator = trial.velocity_gradient
        else:
            indicator = params.starburst_radius - trial.radius
        if indicator < 0:
            high = middle
        elif indicator > 0:
            low = middle
        else:
            break  # exactly on target (the paper code would loop forever here)
    if trial is None:
        raise NoSolutionError("The central-density search could not start.")
    return trial.radius, low


def _brent(equations: WindEquations, params, peak_density: float) -> float:
    """rho_0 in (0, peak_density) whose sonic radius is r_sb, by Brent's method.

    Replaces the paper code's second bisection (decision of 2026-10-15, docs/decisions.md):
    the sonic radius changes smoothly with rho_0 below the peak, so Brent's method (bisection
    combined with secant and inverse quadratic steps) needs about a quarter of the trials. It
    stops at the same absolute precision as the bisection, SEARCH_TOLERANCE * peak_density.
    Falls back to the bisection if the bracket cannot be set up.
    """
    from scipy.optimize import brentq  # imported here: only needed for this search

    def sonic_radius_mismatch(density: float) -> float:
        trial = sonic_trial(equations, central_state(params, density))
        return trial.radius - params.starburst_radius  # > 0: rho_0 too high

    high = peak_density
    if not sonic_radius_mismatch(high) > 0:
        return _bisect(equations, params, 0.0, peak_density, find_peak=False)[1]
    # Lower end of the bracket: halve rho_0 until the sonic radius lies inside r_sb (these are
    # the same first trials the bisection makes).
    low = peak_density / 2.0
    for _ in range(MAX_BRACKET_HALVINGS):
        mismatch = sonic_radius_mismatch(low)
        if mismatch == 0:
            return low
        if mismatch < 0:
            return brentq(sonic_radius_mismatch, low, high, xtol=SEARCH_TOLERANCE * peak_density)
        high, low = low, low / 2.0
    return _bisect(equations, params, 0.0, peak_density, find_peak=False)[1]


def find_central_density(
    params: WindParameters, equations: WindEquations | None = None, method: str = "brent"
):
    """rho_0 for which the wind's sonic radius is r_sb, as in the paper code's ``CRline()``.

    Args:
        params: model inputs.
        equations: the wind equations for ``params`` (built if not given).
        method: how the last stage finds the rho_0 with sonic radius r_sb: ``"brent"``
            (default; the whole solve is about 1.6x faster) or ``"bisection"`` (the paper
            code's method).
            Both stop at the same precision; the first two stages are the paper code's.

    Raises:
        NoSolutionError: if no trial density turns the flow over within the safety limit.
    """
    if method not in ("brent", "bisection"):
        raise ValueError(f"method must be 'brent' or 'bisection', got {method!r}")
    if equations is None:
        equations = WindEquations(params)
    density = FIRST_DENSITY_TRIAL
    trial = sonic_trial(equations, central_state(params, density))
    increases = 0
    while trial.velocity_gradient > 0 and trial.radius < params.starburst_radius:
        increases += 1
        if increases > MAX_DENSITY_INCREASES:
            raise NoSolutionError(
                "No stationary wind solution was found for these inputs: the central density "
                "search did not converge."
            )
        density *= 10.0
        trial = sonic_trial(equations, central_state(params, density))
    sonic_radius, peak_density = _bisect(equations, params, 0.0, density, find_peak=True)
    if sonic_radius > params.starburst_radius:
        if method == "brent":
            density = _brent(equations, params, peak_density)
        else:
            _, density = _bisect(equations, params, 0.0, peak_density, find_peak=False)
    # Otherwise the paper code keeps the last trial of step 1 (see docs/open-items.md).
    return density
