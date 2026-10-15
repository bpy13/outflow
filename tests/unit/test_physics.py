"""Physics functions of the combined model (paper §2.4)."""

import math

import pytest

from outflow import constants as const
from outflow import physics
from outflow.parameters import WindParameters

PARAMS = WindParameters.paper_default()
R_SB = PARAMS.starburst_radius


def test_nfw_series_and_closed_form_agree_at_switch():
    rho_m, r_s = PARAMS.nfw_characteristic_density, PARAMS.nfw_scale_radius
    switch = 1e-3 * r_s
    series = physics.nfw_gravity(switch, rho_m, r_s)  # r <= 1e-3 R_s: series
    closed = physics.nfw_gravity(switch * (1 + 1e-12), rho_m, r_s)  # just above: closed form
    assert series == pytest.approx(closed, rel=1e-9)


def test_nfw_gravity_far_out_is_point_mass():
    r = 1e4 * PARAMS.virial_radius
    gravity = physics.nfw_gravity(r, PARAMS.nfw_characteristic_density, PARAMS.nfw_scale_radius)
    enclosed = PARAMS.halo_mass * (
        physics.nfw_mass_shape(r, PARAMS.nfw_scale_radius)
        / physics.nfw_mass_shape(PARAMS.virial_radius, PARAMS.nfw_scale_radius)
    )
    assert gravity == pytest.approx(const.GRAVITATIONAL_CONSTANT * enclosed / r**2, rel=1e-12)


def test_nfw_gravity_zero_without_halo():
    assert physics.nfw_gravity(R_SB, 0.0, PARAMS.nfw_scale_radius) == 0.0


@pytest.mark.parametrize("offset", [1e-8, 1e-6])
def test_radiative_power_continuous_across_starburst_radius(offset):
    luminosity, kappa = 1e41, 1e4
    at = physics.radiative_power_factor(R_SB, luminosity, kappa, R_SB)
    below = physics.radiative_power_factor(R_SB * (1 - offset), luminosity, kappa, R_SB)
    above = physics.radiative_power_factor(R_SB * (1 + offset), luminosity, kappa, R_SB)
    assert below == pytest.approx(at, rel=100 * offset * abs(math.log(offset)))
    assert above == pytest.approx(at, rel=100 * offset * abs(math.log(offset)))


def test_radiative_power_limits():
    luminosity, kappa = 1e41, 1e4
    # r -> 0: F_rad -> 3 kappa L / (4 pi r_sb^2 c); far out: F_rad -> kappa L / (4 pi r^2 c).
    centre = physics.radiative_power_factor(0.0, luminosity, kappa, R_SB)
    near = physics.radiative_power_factor(1e-4 * R_SB, luminosity, kappa, R_SB)
    assert centre == pytest.approx(near, rel=1e-7)
    far_r = 1e4 * R_SB
    far = physics.radiative_power_factor(far_r, luminosity, kappa, R_SB)
    force = physics.radiative_force(far_r, luminosity, kappa, R_SB, inner=False)
    assert far == pytest.approx(force, rel=1e-6)


def test_effective_sound_speed_finite_at_origin():
    rho0, pressure, cr_pressure = 3e-24, 2e-9, 1e-10
    dv, dva = physics.central_velocity_gradients(
        rho0, PARAMS.mass_injection_rate, PARAMS.magnetic_field, R_SB
    )
    at_origin = physics.effective_sound_speed_squared(pressure, cr_pressure, rho0, dv, dva)
    assert math.isfinite(at_origin) and at_origin > 0
    # L'Hopital: v = v'(0) r and v_A = v_A'(0) r near r = 0 give the same value.
    r = 1e-10 * R_SB
    near = physics.effective_sound_speed_squared(pressure, cr_pressure, rho0, dv * r, dva * r)
    assert near == pytest.approx(at_origin, rel=1e-12)


def test_effective_sound_speed_without_field():
    value = physics.effective_sound_speed_squared(2.0, 3.0, 5.0, 7.0, 0.0)
    assert value == pytest.approx((const.GAMMA_GAS * 2.0 + const.GAMMA_CR * 3.0) / 5.0)


def test_density_and_alfven_speed_switch_as_paper_code():
    # Density switches to the outer form at r = r_sb; the Alfven speed only beyond r_sb.
    assert physics.uses_inner_density(0.999 * R_SB, R_SB)
    assert not physics.uses_inner_density(R_SB, R_SB)
    assert physics.uses_inner_alfven_speed(R_SB, R_SB)
    assert not physics.uses_inner_alfven_speed(1.001 * R_SB, R_SB)
    # While searching for the central density, the inner forms are used at every radius.
    assert physics.uses_inner_density(10 * R_SB, R_SB, searching=True)
    assert physics.uses_inner_alfven_speed(10 * R_SB, R_SB, searching=True)


def test_density_and_alfven_speed_forms():
    v, mass_rate = 3e7, PARAMS.mass_rate
    r_in, r_out = 0.5 * R_SB, 2 * R_SB
    inner = physics.gas_density(r_in, v, R_SB, mass_rate, inner=True)
    assert inner == pytest.approx(PARAMS.mass_injection_rate * r_in / (3 * v), rel=1e-14)
    outer = physics.gas_density(r_out, v, R_SB, mass_rate, inner=False)
    assert outer == pytest.approx(mass_rate / (4 * math.pi * r_out**2 * v), rel=1e-14)
    # The two forms agree at r_sb, so the switch does not make a jump.
    assert physics.gas_density(R_SB, v, R_SB, mass_rate, inner=True) == pytest.approx(
        physics.gas_density(R_SB, v, R_SB, mass_rate, inner=False), rel=1e-14
    )
    b0, rho = PARAMS.magnetic_field, 1e-26
    assert physics.alfven_speed(R_SB, rho, R_SB, b0, inner=True) == pytest.approx(
        physics.alfven_speed(R_SB, rho, R_SB, b0, inner=False), rel=1e-14
    )
    assert physics.alfven_speed(0.0, rho, R_SB, b0, inner=True) == 0.0


# --- cross-check against the paper code -------------------------------------------------------

SAMPLE_STATES = [
    # (r / r_sb, v [cm/s], rho [g/cm^3], P [erg/cm^3], Pc [erg/cm^3])
    (1e-6, 1e3, 3e-24, 2e-9, 1e-10),
    (1e-4, 1e4, 3e-24, 2e-9, 1e-10),
    (0.5, 1.5e7, 1e-24, 1e-9, 5e-11),
    (1.0, 2.5e7, 4e-25, 6e-10, 3e-11),
    (3.0, 5e7, 2e-26, 1e-11, 1e-12),
    (60.0, 6e7, 1e-28, 1e-15, 1e-15),
    (200.0, 6.5e7, 1e-29, 1e-17, 1e-17),
]


@pytest.mark.legacy
@pytest.mark.parametrize("state", SAMPLE_STATES)
@pytest.mark.parametrize("tag", ["find0", "CRline"])
def test_matches_paper_code(crline, state, tag):
    units = crline.units
    params = PARAMS
    args = units.crline_arguments(params)
    r_sb_nat, energy_nat, v_inf_nat, _, _, _, b0_nat, r_s_nat, rho_m_nat = args
    x, v, _, pressure, cr_pressure = state
    r = x * R_SB
    searching = tag == "find0"

    # Gravity (eq. 28-29).
    ours = physics.nfw_gravity(r, params.nfw_characteristic_density, params.nfw_scale_radius)
    theirs = crline.fgrav(r / units.length, r_s_nat, rho_m_nat) * units.acceleration
    assert ours == pytest.approx(theirs, rel=1e-12)

    # Density and Alfven speed (eq. 13, 34), with the paper code's switching rules.
    rho = physics.gas_density(
        r, v, R_SB, params.mass_rate, inner=physics.uses_inner_density(r, R_SB, searching)
    )
    theirs_rho = crline.rhos(
        r / units.length, v / units.velocity, r_sb_nat, energy_nat, v_inf_nat, tag
    )
    assert rho == pytest.approx(theirs_rho * units.density, rel=1e-12)
    v_alfven = physics.alfven_speed(
        r,
        rho,
        R_SB,
        params.magnetic_field,
        inner=physics.uses_inner_alfven_speed(r, R_SB, searching),
    )
    theirs_va = crline.uAs(
        r / units.length, v / units.velocity, r_sb_nat, energy_nat, v_inf_nat, b0_nat, tag
    )
    assert v_alfven == pytest.approx(theirs_va * units.velocity, rel=1e-12)

    # Effective sound speed (eq. 35).
    ours_cs2 = physics.effective_sound_speed_squared(pressure, cr_pressure, rho, v, v_alfven)
    theirs_cs2 = crline.cs2s(
        r / units.length,
        v / units.velocity,
        pressure / units.pressure,
        cr_pressure / units.pressure,
        r_sb_nat,
        energy_nat,
        v_inf_nat,
        b0_nat,
        tag,
    )
    assert ours_cs2 == pytest.approx(theirs_cs2 * units.velocity**2, rel=1e-12)


@pytest.mark.legacy
def test_effective_sound_speed_at_origin_matches_paper_code(crline):
    units = crline.units
    args = units.crline_arguments(PARAMS)
    r_sb_nat, energy_nat, v_inf_nat, *_ = args
    b0_nat = args[6]
    rho0, pressure, cr_pressure = 3e-24, 2e-9, 1e-10
    dv, dva = physics.central_velocity_gradients(
        rho0, PARAMS.mass_injection_rate, PARAMS.magnetic_field, R_SB
    )
    ours = physics.effective_sound_speed_squared(pressure, cr_pressure, rho0, dv, dva)
    # The paper code takes du = dv/dr at r = 0 and recovers rho0 = q / (3 du) from it.
    theirs = crline.cs2s(
        0,
        0,
        pressure / units.pressure,
        cr_pressure / units.pressure,
        r_sb_nat,
        energy_nat,
        v_inf_nat,
        b0_nat,
        "CRline",
        du=dv / units.velocity * units.length,
    )
    assert ours == pytest.approx(theirs * units.velocity**2, rel=1e-12)
