"""Solve for a steady wind: :func:`solve_wind`.

Ported from the paper code's ``CRline()`` (Yu et al. 2020, §2.4): find the central density
rho_0 (see :mod:`outflow.boundary`), then integrate outward from r = 0 with the Cash-Karp
method until the outer radius, or until the wind stalls (falls back to subsonic beyond
1.1 r_sb, e.g. in a heavy halo).
"""

import numpy as np

from outflow import boundary
from outflow.constants import BOLTZMANN, MEAN_PARTICLE_MASS
from outflow.equations import WindEquations
from outflow.errors import NoSolutionError
from outflow.integrator import cash_karp_step, next_step_size, relative_error
from outflow.parameters import WindParameters
from outflow.profile import WindProfile

DEFAULT_OUTER_RADIUS_IN_SB = 50.0  # default outer radius, in units of r_sb
# Beyond this many r_sb, a wind that falls back to subsonic has stalled (paper code: 1.1).
STALL_CHECK_START_IN_SB = 1.1
# Safety limit (not in the paper code) on the number of integration steps.
MAX_STEPS = 2_000_000


def solve_wind(
    params: WindParameters, outer_radius: float | None = None, cooling: bool = False
) -> WindProfile:
    """Solve the steady wind equations for ``params`` out to ``outer_radius`` (cm).

    Args:
        params: model inputs (see :class:`~outflow.parameters.WindParameters`).
        outer_radius: where to stop (cm). Default 50 r_sb. The last point can lie slightly
            beyond it (the last step is not shortened, as in the paper code).
        cooling: experimental radiative cooling (not used in the paper's runs). If no cooling
            solution exists, the wind is solved without cooling and ``cooling_used`` is False.

    Returns:
        The profile from r = 0. If the wind stalls, it ends at the stall radius and
        ``stalled`` is True.

    Raises:
        NoSolutionError: if no stationary wind solution is found.
    """
    if outer_radius is None:
        outer_radius = DEFAULT_OUTER_RADIUS_IN_SB * params.starburst_radius
    if not outer_radius > 0:
        raise ValueError("outer_radius must be greater than zero")

    if cooling:
        # Imported here so that the cooling table is only loaded when cooling is asked for.
        from outflow import cooling as cooling_module

        state, cooling_function = cooling_module.cooling_central_state(params)
        if state is not None:
            try:
                return integrate_outward(
                    WindEquations(params, cooling_function), state, outer_radius
                )
            except NoSolutionError:
                pass  # the cooled wind breaks down; fall back to no cooling
    equations = WindEquations(params)
    central_density = boundary.find_central_density(params, equations)
    state = boundary.central_state(params, central_density)
    return integrate_outward(equations, state, outer_radius)


class _FirstStepFailed(Exception):
    """The step size collapsed on the first step from r = 0 (round-off floor)."""


def integrate_outward(
    equations: WindEquations, state: boundary.CentralState, outer_radius: float
) -> WindProfile:
    """Integrate from r = 0 with the given central state (the main loop of ``CRline()``).

    If the first step from r = 0 cannot meet the tolerance (see
    ``boundary.FIRST_STEP_FALLBACK_TOLERANCE``), the integration is repeated with the looser
    first-step tolerance.
    """
    try:
        return _integrate_outward(equations, state, outer_radius, boundary.STEP_TOLERANCE)
    except _FirstStepFailed:
        return _integrate_outward(
            equations, state, outer_radius, boundary.FIRST_STEP_FALLBACK_TOLERANCE
        )


def _integrate_outward(
    equations: WindEquations,
    state: boundary.CentralState,
    outer_radius: float,
    first_step_tolerance: float,
) -> WindProfile:
    params = equations.params
    starburst_radius = params.starburst_radius
    stall_check_start = STALL_CHECK_START_IN_SB * starburst_radius

    def derivative(r, y):
        return equations.derivatives(r, y[0], y[1], y[2], False)

    r, y = 0.0, (0.0, state.pressure, state.cr_pressure)
    slope = (state.velocity_gradient, 0.0, 0.0)  # dv/dr at r = 0, as the paper code's du
    h = boundary.INITIAL_STEP
    radii, states = [r], [y]
    densities, alfven_speeds = [state.density], [0.0]
    steps = 0
    stalled = False
    while r < outer_radius:
        if r > stall_check_start:
            cs2 = equations.sound_speed_squared(r, y[0], y[1], y[2], False)
            if y[0] ** 2 <= cs2:
                stalled = True  # fell back to subsonic: the wind stalls
                break
        steps += 1
        if steps > MAX_STEPS or h < starburst_radius * boundary.SMALLEST_STEP_FRACTION:
            if r == 0 and first_step_tolerance == boundary.STEP_TOLERANCE:
                raise _FirstStepFailed
            if r <= stall_check_start:
                raise NoSolutionError(
                    "No stationary wind solution was found for these inputs: the integration "
                    "broke down inside the starburst region.",
                    details=f"stopped at r = {r:.6e} cm after {steps} steps, step {h:.3e} cm",
                )
            stalled = True  # the step collapsed at a singular point beyond r_sb
            break
        new_y, error = cash_karp_step(derivative, h, r, y, slope)
        delta = relative_error(new_y, error)
        if delta <= (first_step_tolerance if r == 0 else boundary.STEP_TOLERANCE):
            r, y = r + h, new_y
            rho = equations.density(r, y[0], False)
            radii.append(r)
            states.append(y)
            densities.append(rho)
            alfven_speeds.append(equations.alfven_speed(r, rho, False))
            slope = derivative(r, y)
            h = next_step_size(h, delta, boundary.STEP_TOLERANCE, accepted=True)
        else:
            h = next_step_size(h, delta, boundary.STEP_TOLERANCE, accepted=False)

    values = np.array(states, dtype=float)
    rho = np.array(densities, dtype=float)
    pressure = values[:, 1]
    return WindProfile(
        r=np.array(radii, dtype=float),
        v=values[:, 0],
        rho=rho,
        P=pressure,
        Pc=values[:, 2],
        T=MEAN_PARTICLE_MASS * pressure / (BOLTZMANN * rho),
        vA=np.array(alfven_speeds, dtype=float),
        params=params,
        stalled=stalled,
        cooling_used=equations.cooling_function is not None,
        outer_radius=outer_radius,
    )
