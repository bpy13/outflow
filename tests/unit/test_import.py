"""The package imports cleanly and reports its version."""

import outflow


def test_version():
    assert outflow.__version__ == "0.1.0"
