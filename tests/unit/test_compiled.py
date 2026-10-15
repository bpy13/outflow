"""The Numba-compiled trial integration gives bit-for-bit the same results as the Python code."""

import random

import numpy as np
import pytest

pytest.importorskip("numba")

from outflow import _compiled, boundary, load_paper_profile, paper_runs, solve_wind  # noqa: E402
from outflow.equations import WindEquations  # noqa: E402
from outflow.parameters import WindParameters  # noqa: E402

ALL_RUNS = [run for figure in (6, 7, 8, 9) for run in paper_runs(figure)]


def test_compensated_sum_is_pythons_sum():
    rng = random.Random(3)
    for _ in range(2000):
        count = rng.randint(1, 6)
        values = [rng.uniform(-1.0, 1.0) * 10.0 ** rng.randint(-20, 20) for _ in range(count)]
        if rng.random() < 0.3:  # cancellation, where compensation matters most
            values.append(-sum(values[:-1]) * (1 + 1e-15))
        array = np.array(values)
        assert _compiled._compensated_sum(array, len(values)) == sum(values)


def trial_cases():
    """(name, params, central density): the published rho_0 of every paper run and a factor of
    3 either side (breezes and winds), plus inputs that switch physics off."""
    cases = []
    for run in ALL_RUNS:
        density = float(load_paper_profile(run).rho[0])
        for factor in (1 / 3, 1.0, 3.0):
            cases.append((f"{run.name} x{factor:.2f}", run.params, density * factor))
    default = WindParameters.paper_default()
    switched_off = {
        "no cosmic rays": dict(cr_fraction=0.0),
        "no radiation force": dict(opacity=0.0),
        "no magnetic field": dict(magnetic_field=0.0),
        "no halo": dict(halo_mass=0.0),
    }
    density = float(load_paper_profile(ALL_RUNS[0]).rho[0])
    for name, change in switched_off.items():
        cases.append((name, default.replace(**change), density))
    return cases


@pytest.mark.parametrize("name, params, density", trial_cases(), ids=lambda v: str(v)[:30])
@pytest.mark.parametrize(
    "first_step_tolerance",
    [boundary.STEP_TOLERANCE, boundary.FIRST_STEP_FALLBACK_TOLERANCE],
    ids=["paper tolerance", "first-step fallback"],
)
def test_compiled_trial_is_bit_identical(name, params, density, first_step_tolerance, monkeypatch):
    equations = WindEquations(params)
    state = boundary.central_state(params, density)
    python = boundary._sonic_trial_python(equations, state, first_step_tolerance)
    compiled = boundary._sonic_trial(equations, state, first_step_tolerance)
    assert boundary._compiled_kernel() is not None
    assert compiled == python, name


def test_compiled_trial_reproduces_a_failed_first_step():
    """Fig. 7, 100 uG: at this density the first step from r = 0 fails at the paper tolerance."""
    run = next(r for r in ALL_RUNS if r.name == "fig7_B0_100uG")
    equations = WindEquations(run.params)
    state = boundary.central_state(run.params, 5.973568084651753e-24)
    python = boundary._sonic_trial_python(equations, state, boundary.STEP_TOLERANCE)
    compiled = boundary._sonic_trial(equations, state, boundary.STEP_TOLERANCE)
    assert compiled == python


def test_switch_selects_python(monkeypatch):
    calls = []
    monkeypatch.setattr(boundary, "USE_COMPILED_TRIALS", False)
    monkeypatch.setattr(
        boundary,
        "_sonic_trial_python",
        lambda *args: calls.append(args) or boundary.SonicTrial(0.0, 1.0),
    )
    params = WindParameters.paper_default()
    boundary._sonic_trial(WindEquations(params), boundary.central_state(params, 1e-24), 1e-12)
    assert len(calls) == 1


@pytest.mark.slow
@pytest.mark.parametrize("run", ALL_RUNS, ids=lambda run: run.name)
def test_every_paper_run_is_bit_identical(run, monkeypatch):
    compiled = solve_wind(run.params, outer_radius=run.outer_radius)
    monkeypatch.setattr(boundary, "USE_COMPILED_TRIALS", False)
    python = solve_wind(run.params, outer_radius=run.outer_radius)
    for name in ("r", "v", "rho", "P", "Pc", "T", "vA"):
        np.testing.assert_array_equal(getattr(compiled, name), getattr(python, name), name)
