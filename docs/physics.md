# Physics map: paper equations and code

Equation numbers are those of the published paper, Yu et al. 2020, MNRAS 492, 3179
(arXiv:2001.04384). Everything in the code is in cgs units. The code is a port of the paper's
own code (`CRline.py`), which produced Figs. 6–9; where that code makes a choice the paper does
not spell out, the port copies it and this page says so.

## Symbols

| Paper | Code | Meaning |
|---|---|---|
| r, r_sb | `radius`, `starburst_radius` | radius; radius of the injection region |
| v, ρ, P, P_C | `v`, `rho`, `P`, `Pc` | wind velocity, gas density, gas and cosmic-ray (CR) pressure |
| v_A, B₀ | `vA`, `magnetic_field` | Alfvén speed; field strength at r_sb |
| c\* | `sound_speed` | effective sound speed of gas plus CRs |
| Ṁ, Ė | `mass_rate`, `energy_rate` | mass and energy injection rates |
| q, Q_th, Q_CR | `mass_injection_rate`, `thermal_energy_injection_rate`, `cr_energy_injection_rate` | injection rates per volume inside r_sb |
| v∞ | `terminal_velocity` | √(2Ė/Ṁ) |
| κ, L | `opacity`, `luminosity` | opacity; radiation luminosity left in the wind |
| ρ_M, R_s, R_vir, M_tot | `nfw_characteristic_density`, `nfw_scale_radius`, `virial_radius`, `halo_mass` | NFW halo |
| γ_g, γ_C | `GAMMA_GAS`, `GAMMA_CR` | adiabatic indices 5/3 and 4/3 |
| μ | `MEAN_PARTICLE_MASS` | 1.4 m_H |

## Inputs and derived rates (`parameters.py`)

| Code | Paper |
|---|---|
| `WindParameters.mass_injection_rate` | q, eq. 4 |
| `WindParameters.energy_injection_rate`, `thermal_energy_injection_rate` | Q_th, eq. 5 (times the thermal fraction) |
| `WindParameters.cr_energy_injection_rate` | Q_CR, §2.3 / eq. 44 (times the CR fraction) |
| `WindParameters.terminal_velocity` | v∞² = 2Ė/Ṁ, §2.1.1 |
| `WindParameters.nfw_characteristic_density`, `nfw_total_mass` | ρ_M from M_tot, eqs. 28–29 |

## Local physics (`physics.py`)

| Function | Paper | Notes |
|---|---|---|
| `gas_density` | eq. 13 | inner form q r / (3v), outer form Ṁ / (4π r² v) |
| `alfven_speed` | eq. 34 | B = B₀ r / r_sb inside (footnote in §2.3.2), B₀ r_sb² / r² outside |
| `central_velocity_gradients` | limits of eqs. 13 and 34 at r = 0 | dv/dr = q/(3ρ₀), dv_A/dr = B₀/(r_sb √(4πρ₀)) (L'Hôpital) |
| `effective_sound_speed_squared` | eq. 35 | c\*²; at r = 0 the gradients replace v and v_A |
| `luminosity` | eq. 49 | L = Ė − Ṁ(v²/2 + enthalpy of gas and CRs) |
| `radiative_force` | eq. 25 | f_rad, inner form while subsonic (see below) |
| `radiative_power_factor` | eq. 26 | F_rad, with its limits at r = 0 and r = r_sb |
| `nfw_mass_shape` | eq. 29 | ξ(r) |
| `nfw_gravity` | eqs. 28–29 | G M(r)/r²; four-term series below 10⁻³ R_s |
| `uses_inner_density` | eq. 13 (code convention) | inner density form for r < r_sb |
| `uses_inner_alfven_speed` | eq. 34 (code convention) | inner Alfvén form for r ≤ r_sb |

## Wind equations (`equations.py`)

`WindEquations` holds one parameter set and returns d(v, P, P_C)/dr (paper §2.4, eqs. 36–52).

| Method | Paper | Notes |
|---|---|---|
| `WindEquations.derivatives` | eqs. 36–43 with 50–52 | outer eqs. 36–38 or inner eqs. 39–43, plus the radiation and cooling terms of eqs. 50–52 and gravity (eqs. 28–29) |
| `WindEquations.density` | eq. 13 | with the code's inner/outer switch |
| `WindEquations.alfven_speed` | eq. 34 | with the code's inner/outer switch |
| `WindEquations.sound_speed_squared` | eq. 35 | at r > 0 |
| `WindEquations.cooling_rate` | eq. 10 | C = ρ²Λ(T)/μ²; zero unless cooling is on |

**Inner/outer switch (code convention).** The paper uses the inner equations (39–43, with the
injection terms q, Q_th, Q_CR) for r ≤ r_sb and the outer ones (36–38) beyond. The paper code
switches on the Mach number instead: inner equations while v² ≤ c\*², outer once the flow is
supersonic. The two agree for a wind whose sonic point is at r_sb, which is what the solver
looks for. The radiative force uses the same switch. The density and Alfvén speed, by contrast,
switch at r_sb (density for r < r_sb, Alfvén speed for r ≤ r_sb). The port follows the code.

## Boundary conditions and the central-density search (`boundary.py`)

| Function | Paper | Notes |
|---|---|---|
| `central_state`, `CentralState` | eqs. 44–45 | P₀ and P_C,0 from ρ₀ |
| `central_sound_speed_squared` | eq. 35 at r = 0 | with the L'Hôpital ratio v_A/v |
| `sonic_trial`, `_sonic_trial`, `SonicTrial` | §2.3.2, §2.4 (code: `rsonic()`) | integrate from r = 0 until sonic or turning over |
| `_bisect` | §2.3.2 (code: `gss()`) | bisection on ρ₀ |
| `find_central_density` | §2.3.2, §2.4 (code: `CRline()`) | ρ₀ for which the sonic radius is r_sb |

The search: start at ρ₀ = 0.1 M☉/pc³, multiply by 10 until a trial turns over or reaches r_sb,
bisect for the critical ("peak") ρ₀, then bisect for the ρ₀ whose sonic radius is r_sb. Trials
use the inner density and Alfvén forms at every radius (the code's `'find0'` tag).

## Solver (`solver.py`, `integrator.py`)

- `solve_wind` finds ρ₀ and calls `integrate_outward`, which integrates from r = 0 to the outer
  radius, or stops when the wind falls back to subsonic beyond 1.1 r_sb (a stalled wind).
- `integrator.cash_karp_step`: Cash–Karp Runge–Kutta 4(5) (Press et al., *Numerical Recipes*);
  the fourth-order solution is kept and (5th − 4th) is the error estimate, as in the paper code.
- `integrator.relative_error`: the largest |error/value|, skipping components that are zero.
- `integrator.next_step_size`: h ← 0.9 h (tol/δ)^0.2 (accepted) or ^0.25 (rejected); tol = 1e-12.

## Other modules

| Code | Paper |
|---|---|
| `cc85.cc85_velocity_profile`, `cc85_dimensionless_velocity` | eqs. 7–9 (Chevalier & Clegg 1985 thermal wind) |
| `cooling.cooling_central_density` | eq. 14 (ρ₀ with cooling) |
| `cooling.cooling_function` | Λ(T), Cloudy, solar metallicity (§2.1.2) |
| `cooling.maximum_central_temperature` | eq. 45 with P₀ = ρ₀ k_B T₀/μ |
| `constants.PAPER_ENERGY_RATE` | Ė = 4.2e41 erg/s as in the paper code (Table 1 says 1e42; see known-issues) |
