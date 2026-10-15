"""Quickstart: solve the paper's default wind and a cosmic-ray-driven variant.

Run from the repository root:  python examples/quickstart.py
"""

from outflow import WindParameters, solve_wind
from outflow.constants import KM_S

# The paper's standard M82 parameters (Table 1, with the code's E-dot = 4.2e41 erg/s), in cgs.
thermal = WindParameters.paper_default()
# Same galaxy, but 80 % of the energy goes into cosmic rays (inputs in astronomer's units).
cosmic_ray = WindParameters.from_astro_units(thermal_fraction=0.1, cr_fraction=0.8)

for label, params in (("thermal", thermal), ("cosmic rays", cosmic_ray)):
    profile = solve_wind(params, outer_radius=5 * params.starburst_radius)
    sonic = profile.sound_speed
    print(
        f"{label:>12}: {len(profile)} radii; at 5 r_sb v = {profile.v[-1] / KM_S:.0f} km/s, "
        f"T = {profile.T[-1]:.3g} K, n = {profile.n[-1]:.3g} cm^-3, "
        f"Mach = {profile.v[-1] / sonic[-1]:.2f}"
    )
