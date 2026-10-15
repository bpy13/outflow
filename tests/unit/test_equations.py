"""The wind equations (paper eqs. 36-43 and 50-52) against the paper code's ``CC85()``."""

import math

import pytest

from outflow import physics
from outflow.equations import WindEquations
from outflow.parameters import WindParameters

PARAMS = WindParameters.paper_default()
R_SB = PARAMS.starburst_radius


def test_derivatives_finite_in_both_regimes():
    equations = WindEquations(PARAMS)
    subsonic = equations.derivatives(0.5 * R_SB, 1e7, 2e-9, 1e-10, searching=False)
    supersonic = equations.derivatives(2 * R_SB, 6e7, 1e-12, 1e-12, searching=False)
    for values in (subsonic, supersonic):
        assert len(values) == 3 and all(math.isfinite(x) for x in values)
    assert subsonic[0] > 0  # the subsonic wind accelerates inside the starburst


def test_origin_uses_central_gradient():
    equations = WindEquations(PARAMS)
    rho0 = 3e-24
    gradient, _ = physics.central_velocity_gradients(
        rho0, PARAMS.mass_injection_rate, PARAMS.magnetic_field, R_SB
    )
    values = equations.derivatives(0.0, 0.0, 2e-9, 1e-10, False, central_velocity_gradient=gradient)
    assert values == (gradient, 0.0, 0.0)


def test_no_cr_pressure_stays_zero():
    equations = WindEquations(PARAMS.replace(cr_fraction=0.0))
    for r, v, p in ((0.5 * R_SB, 1e7, 2e-9), (2 * R_SB, 6e7, 1e-12)):
        assert equations.derivatives(r, v, p, 0.0, searching=False)[2] == 0.0


# --- cross-check against the paper code -------------------------------------------------------

# (label, parameter changes, r / r_sb, v [cm/s], P [erg/cm^3], Pc [erg/cm^3], tag)
CASES = [
    ("subsonic", {}, 0.5, 1e7, 2e-9, 1e-10, "CRline"),
    ("subsonic search", {}, 1.3, 2e7, 5e-10, 4e-11, "find0"),
    ("supersonic", {}, 2.0, 6e7, 1e-12, 1e-12, "CRline"),
    ("supersonic inside", {}, 0.8, 4e7, 1e-11, 1e-11, "CRline"),
    ("origin", {}, 0.0, 0.0, 2e-9, 1e-10, "CRline"),
    ("no field", {"magnetic_field": 0.0}, 0.5, 1e7, 2e-9, 1e-10, "CRline"),
    ("no field supersonic", {"magnetic_field": 0.0}, 2.0, 6e7, 1e-12, 1e-12, "CRline"),
    ("no opacity", {"opacity": 0.0}, 0.5, 1e7, 2e-9, 1e-10, "CRline"),
    ("no opacity supersonic", {"opacity": 0.0}, 3.0, 6e7, 1e-13, 1e-13, "CRline"),
    ("no halo", {"halo_mass": 0.0}, 2.0, 6e7, 1e-12, 1e-12, "CRline"),
]


@pytest.mark.legacy
@pytest.mark.parametrize("case", CASES, ids=[c[0] for c in CASES])
def test_derivatives_match_paper_code(crline, case):
    _, changes, x, v, pressure, cr_pressure, tag = case
    params = PARAMS.replace(**changes)
    units = crline.units
    args = units.crline_arguments(params)
    r = x * R_SB
    equations = WindEquations(params)
    central_gradient = 0.0
    if r == 0:
        rho0 = 3e-24
        central_gradient, _ = physics.central_velocity_gradients(
            rho0, params.mass_injection_rate, params.magnetic_field, R_SB
        )
    ours = equations.derivatives(
        r, v, pressure, cr_pressure, tag == "find0", central_velocity_gradient=central_gradient
    )
    theirs = crline.CC85(
        r / units.length,
        v / units.velocity,
        pressure / units.pressure,
        cr_pressure / units.pressure,
        *args,
        tag,
        central_gradient * units.length / units.velocity,
        "off",
    )
    scales = (
        units.velocity / units.length,
        units.pressure / units.length,
        units.pressure / units.length,
    )
    for mine, paper, scale in zip(ours, theirs, scales, strict=True):
        assert mine == pytest.approx(paper * scale, rel=1e-10, abs=0)
