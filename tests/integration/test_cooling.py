"""Experimental radiative cooling (paper §2.1.2): smoke test and the cooling-function spline."""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from outflow import WindParameters, cooling, solve_wind

LEGACY_PAPER_DIR = Path(__file__).resolve().parents[2] / "legacy" / "Brian2019MNRAS"


def test_table_loads_lazily():
    code = "import outflow, outflow.cooling as c; print(c.load_cooling_table.cache_info().currsize)"
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "0"


def test_cooling_function_values():
    temperature, rate = cooling.load_cooling_table()
    assert temperature.shape == rate.shape == (10000,)
    # The spline passes through the table points.
    sample = slice(100, 10000, 997)
    np.testing.assert_allclose(
        cooling.cooling_function(temperature[sample]), rate[sample], rtol=1e-9
    )
    assert cooling.cooling_function(1e6) > 0


def test_solve_with_cooling_paper_default():
    params = WindParameters.paper_default()
    profile = solve_wind(params, cooling=True)
    assert isinstance(profile.cooling_used, bool)
    for name in ("r", "v", "rho", "P", "Pc", "T", "vA"):
        assert np.all(np.isfinite(getattr(profile, name))), name
    print(f"cooling used: {profile.cooling_used}, points {len(profile)}, T0 {profile.T[0]:.4e} K")


@pytest.mark.slow
def test_falls_back_without_cooling_solution():
    # A very weak thermal injection cannot balance cooling at the centre.
    params = WindParameters.paper_default().replace(
        thermal_fraction=0.01, cr_fraction=0.0, mass_rate=26 * 2e33 / (365.25 * 86400)
    )
    state, function = cooling.cooling_central_state(params)
    if state is None:
        assert function is None
        profile = solve_wind(params, outer_radius=2 * params.starburst_radius, cooling=True)
        assert profile.cooling_used is False
    else:
        pytest.skip("a cooling solution exists for these inputs")


@pytest.mark.legacy
def test_spline_matches_legacy_plotcool(monkeypatch):
    """Our spline equals the paper code's (plotcool.py on func_map.dta) to 1e-12."""
    monkeypatch.chdir(LEGACY_PAPER_DIR)  # plotcool reads func_map.dta from the working folder
    namespace: dict = {}
    source = (LEGACY_PAPER_DIR / "plotcool.py").read_text(encoding="utf-8")
    exec(compile(source, "plotcool.py", "exec"), namespace)
    from scipy import interpolate

    temperatures = np.concatenate([np.logspace(1, 9.9, 400), [1.234e4, 2e5, 8.0798e9]])
    theirs = interpolate.splev(temperatures, namespace["tck"])
    ours = cooling.cooling_function(temperatures)
    np.testing.assert_allclose(ours, theirs, rtol=1e-12, atol=0)
