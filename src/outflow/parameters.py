"""Model inputs: :class:`WindParameters`, in cgs, with the paper's defaults (Table 1).

All defaults follow the paper code (``CRline.py``), which produced the published Figs. 6-9.
Use :meth:`WindParameters.from_astro_units` to give inputs in pc, Msun/yr, erg/s, uG, Msun
and kpc.
"""

import dataclasses
import math
from dataclasses import dataclass

from outflow import constants as const
from outflow.physics import nfw_mass_shape


def nfw_total_mass(characteristic_density: float, scale_radius: float, virial_radius: float):
    """Total NFW halo mass M_tot = M(R_vir) = 4 pi rho_M R_s^3 xi(R_vir) (g), paper eq. 28-29."""
    return (
        const.FOUR_PI
        * characteristic_density
        * scale_radius**3
        * nfw_mass_shape(virial_radius, scale_radius)
    )


# Default halo mass: the value implied by the paper code's rho_M (about 5.557e11 Msun).
PAPER_HALO_MASS = nfw_total_mass(
    const.PAPER_NFW_DENSITY, 14.7 * const.KILOPARSEC, 164.0 * const.KILOPARSEC
)

# Tolerance for "thermal + CR fraction <= 1", so that e.g. 0.9 + 0.1 is accepted.
_FRACTION_TOLERANCE = 1e-12


@dataclass(frozen=True)
class WindParameters:
    """Inputs of the combined outflow model (paper §2.4), all in cgs.

    Attributes:
        starburst_radius: r_sb, radius of the uniform injection region (cm).
        energy_rate: total energy injection rate E-dot (erg/s), split into thermal, CR and
            radiation parts.
        mass_rate: mass injection rate M-dot (g/s).
        thermal_fraction: fraction of E-dot injected as thermal energy (Q_th).
        cr_fraction: fraction of E-dot injected as cosmic rays (Q_CR). The rest is radiation.
        opacity: mean opacity kappa of the wind to the starburst radiation (cm^2/g).
        magnetic_field: B_0, field strength at r_sb that sets the Alfven speed (gauss).
        halo_mass: total NFW dark-matter halo mass M_tot = M(R_vir) (g). Zero switches gravity
            off.
        nfw_scale_radius: NFW scale radius R_s (cm).
        virial_radius: virial radius R_vir (cm).
    """

    starburst_radius: float = 200.0 * const.PARSEC
    energy_rate: float = const.PAPER_ENERGY_RATE
    mass_rate: float = 2.6 * const.SOLAR_MASS / const.YEAR
    thermal_fraction: float = 0.8
    cr_fraction: float = 0.1
    opacity: float = 1e4
    magnetic_field: float = 50.0 * const.MICROGAUSS
    halo_mass: float = PAPER_HALO_MASS
    nfw_scale_radius: float = 14.7 * const.KILOPARSEC
    virial_radius: float = 164.0 * const.KILOPARSEC

    def __post_init__(self):
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{field.name} must be a finite number, got {value!r}")
            if value < 0:
                raise ValueError(f"{field.name} must not be negative, got {value!r}")
        for name in ("starburst_radius", "energy_rate", "mass_rate", "thermal_fraction"):
            if getattr(self, name) == 0:
                raise ValueError(f"{name} must be greater than zero")
        for name in ("nfw_scale_radius", "virial_radius"):
            if getattr(self, name) == 0:
                raise ValueError(f"{name} must be greater than zero")
        if self.thermal_fraction + self.cr_fraction > 1 + _FRACTION_TOLERANCE:
            raise ValueError(
                "thermal_fraction + cr_fraction must not exceed 1 (the rest is radiation); "
                f"got {self.thermal_fraction} + {self.cr_fraction}"
            )

    # --- constructors ------------------------------------------------------------------------

    @classmethod
    def paper_default(cls) -> "WindParameters":
        """The paper's standard M82 parameters (Table 1, with the code's E-dot = 4.2e41 erg/s)."""
        return cls()

    @classmethod
    def from_astro_units(
        cls,
        starburst_radius_pc: float = 200.0,
        energy_rate_erg_s: float = const.PAPER_ENERGY_RATE,
        mass_rate_msun_yr: float = 2.6,
        thermal_fraction: float = 0.8,
        cr_fraction: float = 0.1,
        opacity_cm2_g: float = 1e4,
        magnetic_field_microgauss: float = 50.0,
        halo_mass_msun: float = PAPER_HALO_MASS / const.SOLAR_MASS,
        nfw_scale_radius_kpc: float = 14.7,
        virial_radius_kpc: float = 164.0,
    ) -> "WindParameters":
        """Build parameters from astronomer's units (pc, erg/s, Msun/yr, cm^2/g, uG, Msun, kpc)."""
        return cls(
            starburst_radius=starburst_radius_pc * const.PARSEC,
            energy_rate=energy_rate_erg_s,
            mass_rate=mass_rate_msun_yr * const.SOLAR_MASS / const.YEAR,
            thermal_fraction=thermal_fraction,
            cr_fraction=cr_fraction,
            opacity=opacity_cm2_g,
            magnetic_field=magnetic_field_microgauss * const.MICROGAUSS,
            halo_mass=halo_mass_msun * const.SOLAR_MASS,
            nfw_scale_radius=nfw_scale_radius_kpc * const.KILOPARSEC,
            virial_radius=virial_radius_kpc * const.KILOPARSEC,
        )

    def to_astro_units(self) -> dict:
        """The inverse of :meth:`from_astro_units`: keyword arguments in astronomer's units."""
        return {
            "starburst_radius_pc": self.starburst_radius / const.PARSEC,
            "energy_rate_erg_s": self.energy_rate,
            "mass_rate_msun_yr": self.mass_rate * const.YEAR / const.SOLAR_MASS,
            "thermal_fraction": self.thermal_fraction,
            "cr_fraction": self.cr_fraction,
            "opacity_cm2_g": self.opacity,
            "magnetic_field_microgauss": self.magnetic_field / const.MICROGAUSS,
            "halo_mass_msun": self.halo_mass / const.SOLAR_MASS,
            "nfw_scale_radius_kpc": self.nfw_scale_radius / const.KILOPARSEC,
            "virial_radius_kpc": self.virial_radius / const.KILOPARSEC,
        }

    def replace(self, **changes) -> "WindParameters":
        """A copy with some fields changed (validated again)."""
        return dataclasses.replace(self, **changes)

    def as_dict(self) -> dict:
        """Field name -> value (cgs)."""
        return dataclasses.asdict(self)

    # --- derived quantities ------------------------------------------------------------------

    @property
    def radiation_fraction(self) -> float:
        """Fraction of E-dot injected as radiation: 1 - thermal - CR."""
        return max(0.0, 1.0 - self.thermal_fraction - self.cr_fraction)

    @property
    def terminal_velocity(self) -> float:
        """v_inf = sqrt(2 E-dot / M-dot) (cm/s), paper §2.1.1."""
        return math.sqrt(2.0 * self.energy_rate / self.mass_rate)

    @property
    def starburst_volume(self) -> float:
        """Volume of the injection region, 4/3 pi r_sb^3 (cm^3)."""
        return const.FOUR_PI * self.starburst_radius**3 / 3.0

    @property
    def mass_injection_rate(self) -> float:
        """q = 3 M-dot / (4 pi r_sb^3) inside r_sb (g cm^-3 s^-1), paper eq. 4."""
        return 3.0 * self.mass_rate / (const.FOUR_PI * self.starburst_radius**3)

    @property
    def energy_injection_rate(self) -> float:
        """Total 3 E-dot / (4 pi r_sb^3) inside r_sb (erg cm^-3 s^-1), paper eq. 5."""
        return 3.0 * self.energy_rate / (const.FOUR_PI * self.starburst_radius**3)

    @property
    def thermal_energy_injection_rate(self) -> float:
        """Q_th, thermal part of the energy injection rate density, paper eq. 5."""
        return self.energy_injection_rate * self.thermal_fraction

    @property
    def cr_energy_injection_rate(self) -> float:
        """Q_CR, cosmic-ray part of the energy injection rate density, paper eq. 44."""
        return self.energy_injection_rate * self.cr_fraction

    @property
    def nfw_characteristic_density(self) -> float:
        """rho_M of the NFW halo, from M_tot = 4 pi rho_M R_s^3 xi(R_vir) (g/cm^3), eq. 28-29."""
        return self.halo_mass / (
            const.FOUR_PI
            * self.nfw_scale_radius**3
            * nfw_mass_shape(self.virial_radius, self.nfw_scale_radius)
        )
