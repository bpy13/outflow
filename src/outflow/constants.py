"""Physical constants and unit conversions, in cgs.

The values are the rounded constants of the paper code (``CRline.py``, Yu et al. 2020), not
CODATA values, so that results match the published figures exactly.
"""

import math

# Adiabatic indices: thermal gas (gamma_g) and relativistic cosmic rays (gamma_C), paper §2.1, §2.3.
GAMMA_GAS = 5.0 / 3.0
GAMMA_CR = 4.0 / 3.0

# Physical constants (cgs), rounded as in the paper code.
HYDROGEN_MASS = 1.67e-24  # g
MEAN_PARTICLE_MASS = 1.4 * HYDROGEN_MASS  # mu = 1.4 m_H, g (paper §2.1.2)
BOLTZMANN = 1.38e-16  # k_B, erg/K
GRAVITATIONAL_CONSTANT = 6.67e-8  # G, cm^3 g^-1 s^-2
SPEED_OF_LIGHT = 2.998e10  # c, cm/s

# Astronomical units, as in the paper code.
SOLAR_MASS = 2e33  # g
PARSEC = 3.086e18  # cm
KILOPARSEC = 1e3 * PARSEC  # cm
YEAR = 365.25 * 86400.0  # s
KM_S = 1e5  # cm/s
MICROGAUSS = 1e-6  # gauss

# The paper code's energy injection rate: dEdt = 1e12/1.543 in its natural units (Msun, pc,
# km/s; time unit pc/(km/s)). This is 4.2e41 erg/s, not the 1e42 erg/s of the paper's Table 1;
# the published figures used this value (docs/decisions.md).
PAPER_ENERGY_RATE = 1e12 / 1.543 * SOLAR_MASS * KM_S**2 / (PARSEC / KM_S)  # erg/s

# Characteristic NFW density of the M82 halo in the paper code: rho_M = 8.81e-3 Msun/pc^3.
PAPER_NFW_DENSITY = 8.81e-3 * SOLAR_MASS / PARSEC**3  # g/cm^3

FOUR_PI = 4.0 * math.pi
