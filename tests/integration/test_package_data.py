"""Checks on the package data converted from the paper-era files (scripts/convert_paper_data.py)."""

from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "src" / "outflow" / "data"
LEGACY = REPO / "legacy"


def read_csv(path):
    """Return (header dict, column names, data array) of a package CSV file.

    Header lines start with ``#``; those of the form ``# key = value`` go into the dict. The
    first line without ``#`` holds the column names.
    """
    header, names, rows = {}, None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            key, sep, value = line[1:].partition("=")
            if sep:
                header[key.strip()] = value.strip()
        elif names is None:
            names = line.split(",")
        elif line:
            rows.append([float(x) for x in line.split(",")])
    return header, names, np.array(rows)


# --- cooling table ---------------------------------------------------------------------------


def test_cooling_table_layout():
    path = DATA / "cloudy_cooling.csv"
    first = path.read_text(encoding="utf-8").splitlines()[0]
    assert first.startswith("# Cloudy cooling function, solar metallicity")
    _, names, values = read_csv(path)
    assert names == ["T_K", "Lambda_erg_cm3_s"]
    assert values.shape == (10000, 2)


@pytest.mark.legacy
def test_cooling_table_equals_legacy_exactly():
    _, _, values = read_csv(DATA / "cloudy_cooling.csv")
    temperature, cooling = np.load(LEGACY / "Cloudy.npy")
    np.testing.assert_array_equal(values[:, 0], temperature)
    np.testing.assert_array_equal(values[:, 1], cooling)


# --- the 14 published runs -------------------------------------------------------------------

PAPER_RUNS = [
    "fig6_thermal",
    "fig6_cosmic_ray",
    "fig6_radiation",
    "fig7_B0_0uG",
    "fig7_B0_5uG",
    "fig7_B0_50uG",
    "fig7_B0_100uG",
    "fig8_kappa_0",
    "fig8_kappa_1e3",
    "fig8_kappa_1e4",
    "fig8_kappa_1e5",
    "fig9_halo_0x",
    "fig9_halo_1x",
    "fig9_halo_10x",
]
COLUMNS = ["r_cm", "v_cm_s", "rho_g_cm3", "P_erg_cm3", "Pc_erg_cm3", "T_K", "vA_cm_s"]
PARAMETERS = [
    "starburst_radius",
    "energy_rate",
    "mass_rate",
    "thermal_fraction",
    "cr_fraction",
    "opacity",
    "magnetic_field",
    "halo_mass",
    "nfw_scale_radius",
    "virial_radius",
    "outer_radius",
]
# Paper-code constants (cgs).
MEAN_PARTICLE_MASS = 1.4 * 1.67e-24
BOLTZMANN = 1.38e-16


def test_exactly_the_14_paper_runs_exist():
    found = sorted(p.stem for p in (DATA / "paper").glob("*.csv"))
    assert found == sorted(PAPER_RUNS)


@pytest.mark.parametrize("name", PAPER_RUNS)
def test_paper_run_layout_and_temperature(name):
    header, names, values = read_csv(DATA / "paper" / f"{name}.csv")
    assert names == COLUMNS
    for key in PARAMETERS:
        assert float(header[key]) >= 0, key
    assert values.ndim == 2 and values.shape[1] == 7 and len(values) > 1000
    assert np.all(np.isfinite(values))
    r, rho, pressure, temperature = values[:, 0], values[:, 2], values[:, 3], values[:, 5]
    # Non-decreasing, not strictly increasing: the paper data repeat a radius where the step
    # size fell below round-off (near r_sb, and where fig9_halo_10x stalls).
    assert np.all(np.diff(r) >= 0)
    # Each radius's temperature comes from its own density (T = mu P / (k_B rho)).
    expected = MEAN_PARTICLE_MASS * pressure / (BOLTZMANN * rho)
    np.testing.assert_allclose(temperature, expected, rtol=1e-9, atol=0)
