"""WindParameters: paper defaults, derived quantities, unit conversion and validation."""

import math

import pytest

from outflow import constants as const
from outflow.parameters import WindParameters


def test_paper_defaults_in_cgs():
    params = WindParameters.paper_default()
    assert params == WindParameters()
    assert params.starburst_radius == pytest.approx(200 * 3.086e18)
    assert params.energy_rate == pytest.approx(4.2e41, rel=1e-3)
    assert params.mass_rate == pytest.approx(2.6 * 2e33 / (365.25 * 86400))
    assert (params.thermal_fraction, params.cr_fraction) == (0.8, 0.1)
    assert params.radiation_fraction == pytest.approx(0.1)
    assert params.opacity == 1e4
    assert params.magnetic_field == pytest.approx(50e-6)
    assert params.nfw_scale_radius == pytest.approx(14.7e3 * 3.086e18)
    assert params.virial_radius == pytest.approx(164e3 * 3.086e18)
    # M_tot implied by the paper code's rho_M = 8.81e-3 Msun/pc^3 (decisions.md: 5.557e11 Msun).
    assert params.halo_mass / const.SOLAR_MASS == pytest.approx(5.557e11, rel=1e-3)
    assert params.nfw_characteristic_density == pytest.approx(
        8.81e-3 * 2e33 / 3.086e18**3, rel=1e-12
    )


def test_derived_injection_rates():
    params = WindParameters.paper_default()
    volume = 4 / 3 * math.pi * params.starburst_radius**3
    assert params.mass_injection_rate == pytest.approx(params.mass_rate / volume, rel=1e-14)
    assert params.thermal_energy_injection_rate == pytest.approx(
        0.8 * params.energy_rate / volume, rel=1e-14
    )
    assert params.cr_energy_injection_rate == pytest.approx(
        0.1 * params.energy_rate / volume, rel=1e-14
    )


def test_terminal_velocity_is_714_km_s():
    v_inf = WindParameters.paper_default().terminal_velocity
    assert v_inf / 1e5 == pytest.approx(714, abs=1)


def test_replace_returns_new_parameters():
    params = WindParameters.paper_default()
    changed = params.replace(opacity=0.0, magnetic_field=0.0)
    assert changed.opacity == 0 and changed.magnetic_field == 0
    assert params.opacity == 1e4  # the original is unchanged (frozen)
    with pytest.raises(AttributeError):
        params.opacity = 1.0  # type: ignore[misc]


def test_from_astro_units_defaults_match_paper():
    params = WindParameters.from_astro_units()
    default = WindParameters.paper_default()
    for name, value in default.as_dict().items():
        assert getattr(params, name) == pytest.approx(value, rel=1e-12), name


def test_astro_units_round_trip():
    params = WindParameters.paper_default().replace(opacity=3e3, halo_mass=1e45)
    again = WindParameters.from_astro_units(**params.to_astro_units())
    for name, value in params.as_dict().items():
        assert getattr(again, name) == pytest.approx(value, rel=1e-12), name


def test_from_astro_units_converts():
    params = WindParameters.from_astro_units(
        starburst_radius_pc=100, magnetic_field_microgauss=5, halo_mass_msun=1e12
    )
    assert params.starburst_radius == pytest.approx(100 * 3.086e18)
    assert params.magnetic_field == pytest.approx(5e-6)
    assert params.halo_mass == pytest.approx(1e12 * 2e33)


@pytest.mark.parametrize(
    "field",
    ["starburst_radius", "energy_rate", "mass_rate", "opacity", "magnetic_field", "halo_mass"],
)
def test_negative_values_rejected(field):
    with pytest.raises(ValueError, match=field):
        WindParameters.paper_default().replace(**{field: -1.0})


def test_fractions_must_leave_room_for_radiation():
    with pytest.raises(ValueError, match="fraction"):
        WindParameters(thermal_fraction=0.7, cr_fraction=0.4)
    with pytest.raises(ValueError, match="fraction"):
        WindParameters(cr_fraction=-0.1)
    # Exactly 1 is allowed (no radiation), and so are zero CR, opacity, field and halo.
    WindParameters(thermal_fraction=0.9, cr_fraction=0.1)
    WindParameters(cr_fraction=0, opacity=0, magnetic_field=0, halo_mass=0)
