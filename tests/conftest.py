"""Shared test setup.

- Skips ``legacy``-marked tests when the git-ignored legacy/ folder is absent.
- The ``crline`` fixture gives the paper code's functions (legacy/.../version1/CRline.py) for
  cross-checks, plus the factors that convert its natural units (Msun, pc, km/s) to cgs.
"""

import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
LEGACY = REPO / "legacy"
CRLINE = LEGACY / "Brian2019MNRAS" / "version1" / "CRline.py"


class NaturalUnits:
    """cgs value of one natural unit of the paper code (Msun, pc, km/s), its own rounding."""

    mass = 2e33  # Msun, g
    length = 3.086e18  # pc, cm
    velocity = 1e5  # km/s, cm/s
    time = length / velocity
    density = mass / length**3
    pressure = mass * velocity**2 / length**3
    acceleration = velocity**2 / length
    energy_rate = mass * velocity**2 / time
    opacity = length**2 / mass
    magnetic_field = pressure**0.5  # B^2 has units of pressure in Gaussian units

    @classmethod
    def crline_arguments(cls, params):
        """WindParameters -> CRline's (rSB, dEdt, v_inf, acr, ath, k, B0, rM, rhoM)."""
        energy_rate = params.energy_rate / cls.energy_rate
        mass_rate = params.mass_rate / (cls.mass / cls.time)
        return (
            params.starburst_radius / cls.length,
            energy_rate,
            (2 * energy_rate / mass_rate) ** 0.5,
            params.cr_fraction,
            params.thermal_fraction,
            params.opacity / cls.opacity,
            params.magnetic_field / cls.magnetic_field,
            params.nfw_scale_radius / cls.length,
            params.nfw_characteristic_density / cls.density,
        )


@pytest.fixture(scope="session")
def crline():
    """The paper code's function block, executed without its run loop (which takes minutes).

    The block is everything before ``# setting up parameter space``. Its module-level
    ``import plotcool`` would read ``func_map.dta`` from the working directory, so a stub module
    is put in its place; ``cool()`` is therefore not usable from this namespace.
    """
    source = CRLINE.read_text(encoding="utf-8")
    block = source.split("# setting up parameter space")[0]
    stub = types.ModuleType("plotcool")
    stub.tck = None
    saved = sys.modules.get("plotcool")
    sys.modules["plotcool"] = stub
    try:
        namespace: dict = {}
        exec(compile(block, str(CRLINE), "exec"), namespace)
    finally:
        if saved is None:
            sys.modules.pop("plotcool", None)
        else:
            sys.modules["plotcool"] = saved
    module = types.SimpleNamespace(**namespace)
    module.units = NaturalUnits
    return module


def pytest_collection_modifyitems(config, items):
    if LEGACY.is_dir():
        return
    skip = pytest.mark.skip(reason="needs the local legacy/ folder (paper-era code and data)")
    for item in items:
        if "legacy" in item.keywords:
            item.add_marker(skip)
