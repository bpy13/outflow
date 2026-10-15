"""Reproduce the published runs (Yu et al. 2020, Figs. 6-9) from the package data.

Comparison as in docs/decisions.md: the published profile is interpolated linearly in r onto
the new radii between 0.2 r_sb and the smaller outer radius; v, rho, P, P_C, T and v_A must
agree to 1e-3 relative.
"""

import time

import numpy as np
import pytest

from outflow import NoSolutionError, WindProfile, load_paper_profile, paper_runs, solve_wind
from outflow.constants import BOLTZMANN, MEAN_PARTICLE_MASS

TOLERANCE = 1e-3


def run_by_name(name):
    return next(run for run in paper_runs() if run.name == name)


def compare(name, outer_in_sb=None):
    run = run_by_name(name)
    reference = load_paper_profile(run)
    outer = run.outer_radius if outer_in_sb is None else outer_in_sb * run.params.starburst_radius
    start = time.perf_counter()
    profile = solve_wind(run.params, outer_radius=outer)
    seconds = time.perf_counter() - start
    differences = profile.max_relative_differences(reference)
    print(f"{name}: {seconds:.1f} s, max relative differences {differences}")
    return profile, reference, differences


def test_paper_run_parameters_match_saved_headers():
    for run in paper_runs():
        saved = load_paper_profile(run)
        for name, value in run.params.as_dict().items():
            assert getattr(saved.params, name) == pytest.approx(value, rel=1e-12), (run.name, name)
        assert saved.outer_radius == pytest.approx(run.outer_radius, rel=1e-12)


def test_paper_runs_by_figure():
    assert [len(paper_runs(f)) for f in (6, 7, 8, 9)] == [3, 4, 4, 3]
    assert len(paper_runs()) == 14
    with pytest.raises(ValueError):
        paper_runs(5)


def test_fig6_thermal_to_10_rsb():
    profile, _, differences = compare("fig6_thermal", outer_in_sb=10)
    assert profile.r[-1] >= 10 * profile.params.starburst_radius
    assert not profile.stalled
    for column, value in differences.items():
        assert value <= TOLERANCE, (column, value)


def test_temperature_from_each_radius_own_density():
    run = run_by_name("fig6_thermal")
    profile = solve_wind(run.params, outer_radius=2 * run.params.starburst_radius)
    expected = MEAN_PARTICLE_MASS * profile.P / (BOLTZMANN * profile.rho)
    np.testing.assert_allclose(profile.T, expected, rtol=1e-12, atol=0)
    assert np.all(np.diff(profile.T[1:]) != 0)  # not a constant (the outflow.py bug)


def test_csv_round_trip_is_exact(tmp_path):
    reference = load_paper_profile("fig7_B0_0uG")
    first = reference.to_csv(tmp_path / "first.csv")
    again = WindProfile.from_csv(first)
    for name in ("r", "v", "rho", "P", "Pc", "T", "vA"):
        np.testing.assert_array_equal(getattr(again, name), getattr(reference, name))
    assert again.params == reference.params
    assert again.outer_radius == reference.outer_radius
    second = again.to_csv(tmp_path / "second.csv")
    assert first.read_bytes() == second.read_bytes()


def test_computed_profile_csv_round_trip(tmp_path):
    run = run_by_name("fig6_thermal")
    profile = solve_wind(run.params, outer_radius=1.5 * run.params.starburst_radius)
    written = WindProfile.from_csv(profile.to_csv(tmp_path / "wind.csv"))
    assert written.params == profile.params
    assert (written.stalled, written.cooling_used) == (profile.stalled, profile.cooling_used)
    np.testing.assert_allclose(written.v, profile.v, rtol=1e-10)
    rewritten = WindProfile.from_csv(written.to_csv(tmp_path / "again.csv"))
    np.testing.assert_array_equal(rewritten.v, written.v)


def test_derived_columns():
    profile = load_paper_profile("fig6_thermal")
    np.testing.assert_allclose(profile.n, profile.rho / MEAN_PARTICLE_MASS)
    np.testing.assert_allclose(profile.ram_pressure, profile.rho * profile.v**2)
    assert np.all(np.isfinite(profile.sound_speed))
    # The thermal wind passes the sonic point near r_sb.
    r_sb = profile.params.starburst_radius
    inside = (profile.r > 0.5 * r_sb) & (profile.r < 0.9 * r_sb)
    outside = profile.r > 1.5 * r_sb
    assert np.all(profile.v[inside] < profile.sound_speed[inside])
    assert np.all(profile.v[outside] > profile.sound_speed[outside])


def test_no_solution_error_is_plain_language():
    error = NoSolutionError("No stationary wind solution was found.", details="x")
    assert str(error) == "No stationary wind solution was found."
    assert error.details == "x"


@pytest.mark.slow
def test_heavy_halo_stalls():
    profile, _, differences = compare("fig9_halo_10x")
    assert profile.stalled
    r_sb = profile.params.starburst_radius
    assert profile.stall_radius == pytest.approx(44.8 * r_sb, rel=0.02)
    assert profile.r[-1] < 50 * r_sb
