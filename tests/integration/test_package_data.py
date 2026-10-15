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
