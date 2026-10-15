"""Shared test setup: skip ``legacy``-marked tests when the git-ignored legacy/ folder is absent."""

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
LEGACY = REPO / "legacy"


def pytest_collection_modifyitems(config, items):
    if LEGACY.is_dir():
        return
    skip = pytest.mark.skip(reason="needs the local legacy/ folder (paper-era code and data)")
    for item in items:
        if "legacy" in item.keywords:
            item.add_marker(skip)
