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


def test_curves_are_trimmed_to_the_shown_range():
    # fig6_thermal was computed to 100 r_sb; with x_max = 5 nothing far beyond the axis edge
    # may be drawn, or the y axes are scaled to values that are never shown.
    profile = load_paper_profile("fig6_thermal")
    for figure in (
        plot_profiles([profile], x_max=5),
        plot_comparison(profile, [load_paper_profile("fig6_cosmic_ray")], x_max=5),
    ):
        for axis in figure.axes:
            for line in axis.get_lines():
                x = line.get_xdata()
                if len(x) > 2:  # skip the reference line at 1 and legend proxies
                    assert max(x) < 6, max(x)
