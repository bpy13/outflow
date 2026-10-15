"""Central boundary values (paper eqs. 44-45) and the first-step fallback at r = 0."""

import pytest

from outflow import boundary, paper_runs
from outflow.constants import GAMMA_CR, GAMMA_GAS
from outflow.equations import WindEquations
from outflow.parameters import WindParameters


def test_central_state_eqs_44_45():
    params = WindParameters.paper_default()
    rho0 = 3e-24
    state = boundary.central_state(params, rho0)
    q = params.mass_injection_rate
    # Eq. 45: gamma_g/(gamma_g - 1) P_0 / rho_0 = Q_th / q.
    assert GAMMA_GAS / (GAMMA_GAS - 1) * state.pressure / rho0 == pytest.approx(
        params.thermal_energy_injection_rate / q, rel=1e-14
    )
    # Eq. 44: gamma_C/(gamma_C - 1) (v + v_A)/v P_C,0 / rho_0 = Q_CR / q, with L'Hopital ratios.
    ratio = (state.velocity_gradient + state.alfven_gradient) / state.velocity_gradient
    assert GAMMA_CR / (GAMMA_CR - 1) * ratio * state.cr_pressure / rho0 == pytest.approx(
        params.cr_energy_injection_rate / q, rel=1e-14
    )
    assert state.velocity_gradient == pytest.approx(q / (3 * rho0), rel=1e-14)


def test_unknown_search_method_is_rejected():
    with pytest.raises(ValueError, match="brent"):
        boundary.find_central_density(WindParameters.paper_default(), method="secant")


@pytest.mark.slow
@pytest.mark.parametrize("run_name", ["fig6_thermal", "fig7_B0_100uG", "fig8_kappa_1e5"])
def test_brent_search_finds_the_bisection_density(run_name):
    """Brent's method (default) and the paper code's bisection agree on rho_0.

    Both stop when rho_0 is known to 1e-12 of the search range, so they agree to a few 1e-12.
    """
    run = next(r for figure in (6, 7, 8) for r in paper_runs(figure) if r.name == run_name)
    brent = boundary.find_central_density(run.params, method="brent")
    bisection = boundary.find_central_density(run.params, method="bisection")
    assert brent == pytest.approx(bisection, rel=1e-11)


def test_first_step_fallback_rescues_round_off_start():
    """A density where the first step from r = 0 sits on its round-off floor (Fig. 7, 100 uG).

    With the paper code's tolerance alone the step size collapsed at r = 0 (on Windows, Python
    3.12; the exact round-off can differ between platforms); the fallback lets the trial go on
    to its sonic point, just beyond r_sb.
    """
    run = next(r for r in paper_runs(7) if r.name == "fig7_B0_100uG")
    equations = WindEquations(run.params)
    state = boundary.central_state(run.params, 5.973568084651753e-24)
    trial = boundary.sonic_trial(equations, state)
    assert trial.radius == pytest.approx(run.params.starburst_radius, rel=1e-6)
    assert 0 < trial.velocity_gradient < boundary.TRIAL_FAILED
