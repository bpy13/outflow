"""Matplotlib helpers save figures without error (Agg backend, no display)."""

import matplotlib

matplotlib.use("Agg")

import pytest  # noqa: E402

from outflow import load_paper_profile, paper_runs  # noqa: E402
from outflow.plotting import plot_comparison, plot_profiles  # noqa: E402


@pytest.fixture(scope="module")
def fig7_profiles():
    return [load_paper_profile(run) for run in paper_runs(7)]


def test_plot_profiles_saves_png(tmp_path, fig7_profiles):
    figure = plot_profiles(fig7_profiles, labels=[r.label for r in paper_runs(7)], x_max=5)
    assert len(figure.axes) == 4
    path = tmp_path / "profiles.png"
    figure.savefig(path)
    assert path.stat().st_size > 10_000


def test_plot_comparison_saves_png(tmp_path, fig7_profiles):
    reference, *others = fig7_profiles
    figure = plot_comparison(reference, others, labels=["5 uG", "50 uG", "100 uG"], x_max=5)
    assert len(figure.axes) == 4  # ratio panels: v, n, T, pressures
    path = tmp_path / "comparison.png"
    figure.savefig(path)
    assert path.stat().st_size > 10_000


def test_too_many_profiles_rejected(fig7_profiles):
    with pytest.raises(ValueError):
        plot_profiles(fig7_profiles * 3)
