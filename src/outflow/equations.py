"""The combined wind equations (Yu et al. 2020, §2.4), ported from the paper code's ``CC85()``.

The state is (v, P, P_C) as a function of radius r. Two sets of equations are used, as in the
paper code:

- **Inner** (source-term) equations, paper eqs. 39-43 with the radiation and cooling terms of
  eqs. 50-52, while the flow is subsonic (v^2 <= c*^2). They include the injection terms q,
  Q_th and Q_CR.
- **Outer** equations, paper eqs. 36-38 with eqs. 50-52, once the flow is supersonic.

The paper describes the switch at r = r_sb; the paper code switches on the Mach number instead
(the two agree for a wind whose sonic point is at r_sb). This port follows the code.
"""

from collections.abc import Callable

from outflow import physics
from outflow.constants import BOLTZMANN, GAMMA_CR, GAMMA_GAS, MEAN_PARTICLE_MASS
from outflow.parameters import WindParameters

CoolingFunction = Callable[[float], float]


class WindEquations:
    """Right-hand side d(v, P, P_C)/dr of the wind equations for one set of parameters.

    Args:
        params: model inputs.
        cooling_function: Lambda(T) in erg cm^3 s^-1, or None for no radiative cooling (the
            paper's default; the published runs never used cooling).
    """

    def __init__(self, params: WindParameters, cooling_function: CoolingFunction | None = None):
        self.params = params
        self.cooling_function = cooling_function
        # Cache the scalars used at every evaluation.
        self.starburst_radius = params.starburst_radius
        self.mass_rate = params.mass_rate
        self.terminal_velocity = params.terminal_velocity
        self.mass_injection_rate = params.mass_injection_rate
        self.thermal_injection_rate = params.thermal_energy_injection_rate
        self.cr_injection_rate = params.cr_energy_injection_rate
        self.opacity = params.opacity
        self.magnetic_field = params.magnetic_field
        self.nfw_characteristic_density = params.nfw_characteristic_density
        self.nfw_scale_radius = params.nfw_scale_radius

    def density(self, r: float, v: float, searching: bool) -> float:
        """Gas density at (r, v) with the paper code's inner/outer switch (eq. 13)."""
        inner = physics.uses_inner_density(r, self.starburst_radius, searching)
        return physics.gas_density(r, v, self.starburst_radius, self.mass_rate, inner)

    def alfven_speed(self, r: float, density: float, searching: bool) -> float:
        """Alfven speed at r with the paper code's inner/outer switch (eq. 34)."""
        inner = physics.uses_inner_alfven_speed(r, self.starburst_radius, searching)
        return physics.alfven_speed(r, density, self.starburst_radius, self.magnetic_field, inner)

    def sound_speed_squared(self, r: float, v: float, P: float, Pc: float, searching: bool):
        """Effective sound speed squared c*^2 at r > 0 (eq. 35)."""
        rho = self.density(r, v, searching)
        v_alfven = self.alfven_speed(r, rho, searching)
        return physics.effective_sound_speed_squared(P, Pc, rho, v, v_alfven)

    def cooling_rate(self, density: float, pressure: float) -> float:
        """Radiative cooling rate per volume C = rho^2 Lambda(T) / mu^2 (eq. 10), or 0."""
        if self.cooling_function is None:
            return 0.0
        temperature = MEAN_PARTICLE_MASS * pressure / (BOLTZMANN * density)
        return self.cooling_function(temperature) * (density / MEAN_PARTICLE_MASS) ** 2

    def derivatives(
        self,
        r: float,
        v: float,
        P: float,
        Pc: float,
        searching: bool,
        central_velocity_gradient: float = 0.0,
    ) -> tuple[float, float, float]:
        """(dv/dr, dP/dr, dP_C/dr) at radius r, paper eqs. 36-43 and 50-52.

        Args:
            r, v, P, Pc: radius (cm), velocity (cm/s), gas and CR pressure (erg/cm^3).
            searching: True while searching for the central density (inner density and Alfven
                speed forms at every radius, as the paper code's tag ``'find0'``).
            central_velocity_gradient: dv/dr at r = 0, q / (3 rho_0). At r = 0 the paper code
                returns (dv/dr, 0, 0) with this L'Hopital value; it is only used there.
        """
        if r == 0:
            # Centre: v = 0, and the gradients follow from the L'Hopital limit (paper §2.2.2).
            return (central_velocity_gradient, 0.0, 0.0)

        gamma, gamma_cr = GAMMA_GAS, GAMMA_CR
        rho = self.density(r, v, searching)
        v_alfven = self.alfven_speed(r, rho, searching)
        cooling = self.cooling_rate(rho, P)
        cs2 = physics.effective_sound_speed_squared(P, Pc, rho, v, v_alfven)
        luminosity = physics.luminosity(
            v, P, Pc, rho, (v + v_alfven) / v, self.mass_rate, self.terminal_velocity
        )
        power_factor = physics.radiative_power_factor(
            r, luminosity, self.opacity, self.starburst_radius
        )
        gravity = physics.nfw_gravity(r, self.nfw_characteristic_density, self.nfw_scale_radius)
        subsonic = v**2 <= cs2
        force = physics.radiative_force(
            r, luminosity, self.opacity, self.starburst_radius, inner=subsonic
        )
        # Radiation, cooling and gravity terms of eq. 52 (dv/dr), common to both regimes.
        external = (
            (gamma - 1.0) * (cooling / rho - power_factor * v) + gamma * force * v - gravity * v
        ) / (v**2 - cs2)
        # Coupling factors of the CR pressure (eqs. 36-37 and 39-40).
        cr_factor = gamma_cr * Pc * (2.0 * v + v_alfven) / (2.0 * (v + v_alfven))
        gas_factor = gamma * P - gamma_cr * Pc * (gamma - 1.0) * v_alfven * (2.0 * v + v_alfven) / (
            2.0 * v * (v + v_alfven)
        )

        if subsonic:
            # Inner equations with injection terms, eqs. 39-43.
            source_cr = -3.0 * gamma_cr * Pc / r + (gamma_cr - 1.0) * self.cr_injection_rate / (
                v + v_alfven
            )  # C_2
            source_gas = (
                -(gamma - 1.0) * v_alfven * source_cr / v
                + (gamma - 1.0) * self.thermal_injection_rate / v
                - 3.0 * gamma * P / r
                + 3.0 * (gamma - 1.0) * rho * v**2 / (2.0 * r)
            )  # C_1
            dv_cr = cs2 * v / (r * (cs2 - v**2)) + v * (
                source_gas + source_cr + self.mass_injection_rate * v
            ) / (rho * (cs2 - v**2))
            dv_dr = dv_cr + external
            log_density_gradient = 1.0 / r - dv_dr / v  # d ln(rho)/dr for rho = q r / (3 v)
            dpc_dr = log_density_gradient * cr_factor + source_cr
            dp_dr = log_density_gradient * gas_factor + source_gas
        else:
            # Outer equations, eqs. 36-38.
            dv_dr = 2.0 * cs2 * v / (r * (v**2 - cs2)) + external
            log_density_gradient = -2.0 / r - dv_dr / v  # d ln(rho)/dr for rho ~ 1 / (r^2 v)
            dpc_dr = log_density_gradient * cr_factor
            dp_dr = log_density_gradient * gas_factor
        # Radiative heating, radiative force and cooling in the gas pressure, eq. 51.
        dp_dr += (gamma - 1.0) * (rho * power_factor - rho * force - cooling / v)
        return (dv_dr, dp_dr, dpc_dr)
