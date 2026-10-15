"""Local physics of the combined outflow model (Yu et al. 2020, §2.4), in cgs.

Each function is a single formula of the paper, written for scalar floats (they sit in the
inner loop of the integrator, so they avoid numpy overhead). Equation numbers refer to the
published paper (MNRAS 492, 3179).

Where the paper code (``CRline.py``) makes a choice the paper does not spell out, the choice is
copied and noted here:

- Inside the starburst the density and Alfven speed take their "inner" forms (q r / 3v and
  B_0 r / r_sb). The density switches to the outer form at r = r_sb, the Alfven speed only
  beyond r_sb; both forms agree at r_sb, so the switch is continuous. While searching for the
  central density the inner forms are used at every radius.
- Gravity uses a four-term series below 1e-3 R_s, where the closed form loses precision.
"""

import math

from outflow.constants import (
    FOUR_PI,
    GAMMA_CR,
    GAMMA_GAS,
    GRAVITATIONAL_CONSTANT,
    SPEED_OF_LIGHT,
)

# Below this fraction of the NFW scale radius, gravity is evaluated from its series expansion.
NFW_SERIES_LIMIT = 1e-3
NFW_SERIES_TERMS = 4


# --- switching rules of the paper code -------------------------------------------------------


def uses_inner_density(radius: float, starburst_radius: float, searching: bool = False) -> bool:
    """True where the density takes its inner form q r / (3 v): r < r_sb, or while searching."""
    return searching or radius < starburst_radius


def uses_inner_alfven_speed(
    radius: float, starburst_radius: float, searching: bool = False
) -> bool:
    """True where v_A takes its inner form (B = B_0 r / r_sb): r <= r_sb, or while searching."""
    return searching or radius <= starburst_radius


# --- local quantities -------------------------------------------------------------------------


def gas_density(
    radius: float, velocity: float, starburst_radius: float, mass_rate: float, inner: bool
) -> float:
    """Gas density from mass conservation, paper eq. 13 (with eq. 1, 4).

    Inner form: rho = q r / (3 v) = M-dot r / (4 pi r_sb^3 v).
    Outer form: rho = M-dot / (4 pi r^2 v).
    """
    if inner:
        return mass_rate * radius / (FOUR_PI * starburst_radius**3 * velocity)
    return mass_rate / (FOUR_PI * radius**2 * velocity)


def alfven_speed(
    radius: float, density: float, starburst_radius: float, magnetic_field: float, inner: bool
) -> float:
    """Alfven speed v_A = B / sqrt(4 pi rho), paper eq. 34 and footnote of §2.3.2.

    Inner form: B = B_0 r / r_sb (a tractable choice for the boundary condition, paper §2.3.2).
    Outer form: B = B_0 r_sb^2 / r^2 (magnetic flux 4 pi r^2 B conserved).
    """
    if inner:
        field = magnetic_field * radius / starburst_radius
    else:
        field = magnetic_field * starburst_radius**2 / radius**2
    return field / math.sqrt(FOUR_PI * density)


def central_velocity_gradients(
    central_density: float,
    mass_injection_rate: float,
    magnetic_field: float,
    starburst_radius: float,
) -> tuple[float, float]:
    """dv/dr and dv_A/dr at r = 0 (L'Hopital limits), paper §2.2.2 footnote and §2.3.2.

    dv/dr = q / (3 rho_0) and dv_A/dr = B_0 / (r_sb sqrt(4 pi rho_0)). Both v and v_A vanish
    linearly at r = 0, so ratios such as v_A / v take the ratio of these gradients there.
    """
    velocity_gradient = mass_injection_rate / (3.0 * central_density)
    alfven_gradient = magnetic_field / (starburst_radius * math.sqrt(FOUR_PI * central_density))
    return velocity_gradient, alfven_gradient


def effective_sound_speed_squared(
    pressure: float, cr_pressure: float, density: float, velocity: float, alfven: float
) -> float:
    """Effective sound speed squared c*^2 of gas plus cosmic rays, paper eq. 35.

    c*^2 = gamma_g P / rho
           + gamma_C P_C / rho * (2 v + v_A)(v - (gamma_g - 1) v_A) / (2 v (v + v_A)).

    The CR term depends on v and v_A only through v_A / v, so at r = 0 pass the gradients from
    :func:`central_velocity_gradients` instead (the L'Hopital limit used by the paper code).
    """
    cr_factor = (
        (2.0 * velocity + alfven)
        * (velocity - (GAMMA_GAS - 1.0) * alfven)
        / (2.0 * velocity * (velocity + alfven))
    )
    return GAMMA_GAS * pressure / density + GAMMA_CR * cr_pressure * cr_factor / density


def luminosity(
    velocity: float,
    pressure: float,
    cr_pressure: float,
    density: float,
    cr_velocity_ratio: float,
    mass_rate: float,
    terminal_velocity: float,
) -> float:
    """Radiation luminosity left after the wind's energy flux is taken out, paper eq. 49.

    L = E-dot - M-dot (v^2/2 + gamma_g/(gamma_g-1) P/rho
                       + gamma_C/(gamma_C-1) (v + v_A)/v P_C/rho),
    with E-dot = M-dot v_inf^2 / 2. ``cr_velocity_ratio`` is (v + v_A) / v; at r = 0 pass the
    ratio of the central gradients (and velocity = 0).
    """
    return mass_rate * (
        terminal_velocity**2 / 2.0
        - velocity**2 / 2.0
        - GAMMA_GAS * pressure / ((GAMMA_GAS - 1.0) * density)
        - GAMMA_CR * cr_velocity_ratio * cr_pressure / (density * (GAMMA_CR - 1.0))
    )


def radiative_force(
    radius: float, luminosity_value: float, opacity: float, starburst_radius: float, inner: bool
) -> float:
    """Radiative force per unit mass f_rad, paper eq. 25.

    Inner form: kappa L r / (4 pi r_sb^3 c). Outer form: kappa L / (4 pi r^2 c). The paper code
    uses the inner form while the flow is subsonic and the outer form once it is supersonic.
    """
    if inner:
        return (
            opacity * luminosity_value * radius / (FOUR_PI * SPEED_OF_LIGHT * starburst_radius**3)
        )
    return opacity * luminosity_value / (FOUR_PI * SPEED_OF_LIGHT * radius**2)


def radiative_power_factor(
    radius: float, luminosity_value: float, opacity: float, starburst_radius: float
) -> float:
    """F_rad, the radiative heating term per unit mass and velocity (F_rad v is the power), eq. 26.

    F_rad = 3 kappa L (2 r_sb r - (r^2 - r_sb^2) ln((r + r_sb)/|r - r_sb|)) / (16 pi r_sb^3 c r),
    with its limits 3 kappa L / (4 pi r_sb^2 c) at r = 0 and 3 kappa L / (8 pi r_sb^2 c) at
    r = r_sb, where the formula is 0/0 or has a log singularity times zero.
    """
    if radius == 0:
        return 3.0 * opacity * luminosity_value / (FOUR_PI * starburst_radius**2 * SPEED_OF_LIGHT)
    if radius == starburst_radius:
        return (
            3.0
            * opacity
            * luminosity_value
            / (8.0 * math.pi * starburst_radius**2 * SPEED_OF_LIGHT)
        )
    shape = 2.0 * starburst_radius * radius - (radius**2 - starburst_radius**2) * math.log(
        (radius + starburst_radius) / abs(radius - starburst_radius)
    )
    return (
        3.0
        * opacity
        * luminosity_value
        * shape
        / (16.0 * math.pi * starburst_radius**3 * SPEED_OF_LIGHT * radius)
    )


def nfw_mass_shape(radius: float, scale_radius: float) -> float:
    """NFW mass shape xi(r) = ln((R_s + r)/R_s) - r/(R_s + r), paper eq. 29.

    The enclosed halo mass is M(r) = 4 pi rho_M R_s^3 xi(r).
    """
    return math.log((scale_radius + radius) / scale_radius) - radius / (scale_radius + radius)


def nfw_gravity(radius: float, characteristic_density: float, scale_radius: float) -> float:
    """Magnitude of the NFW gravitational acceleration G M(r) / r^2, paper eq. 28-29.

    Written with the characteristic density rho_M: G M(r)/r^2 = 4 pi G rho_M R_s^3 xi(r)/r^2.
    Below 1e-3 R_s the closed form loses precision (xi ~ r^2/2 is a small difference of two
    terms), so the paper code's four-term series is used:
    4 pi G rho_M R_s sum_{i<4} (i+1)(-r/R_s)^i / (i+2).
    """
    if radius > scale_radius * NFW_SERIES_LIMIT:
        return (
            FOUR_PI
            * GRAVITATIONAL_CONSTANT
            * characteristic_density
            * scale_radius**3
            * nfw_mass_shape(radius, scale_radius)
            / radius**2
        )
    x = -radius / scale_radius
    series = sum((i + 1) * x**i / (i + 2) for i in range(NFW_SERIES_TERMS))
    return FOUR_PI * GRAVITATIONAL_CONSTANT * characteristic_density * scale_radius * series
