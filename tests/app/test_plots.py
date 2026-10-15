"""Plotly figure builders of the app (app/plots.py)."""

import sys
from pathlib import Path

from outflow import load_paper_profile

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "app"))
import plots  # noqa: E402

RATIO_AXES = {"x2", "x4", "x6", "x8"}  # right-hand column of the 4 x 2 grid
COMPUTED_COLOR = "#123456"


def ratio_traces_with_color(figure, color):
    return [
        trace
        for trace in figure.data
        if trace.xaxis in RATIO_AXES and color in str(trace.to_plotly_json())
    ]


def test_ratios_never_fall_back_to_another_run_as_baseline():
    first = load_paper_profile("fig7_B0_0uG")
    second = load_paper_profile("fig7_B0_5uG")
    published = [
        plots.Series("B0 = 0", first, "#aa0000", published=True),
        plots.Series("B0 = 5", second, "#00aa00", published=True),
    ]
    # The recomputed baseline run ("B0 = 0") failed; only "B0 = 5" was recomputed.
    computed = [plots.Series("B0 = 5", second, COMPUTED_COLOR)]
    figure = plots.profile_and_ratio_figure(computed + published, baseline_label="B0 = 0")
    assert not ratio_traces_with_color(figure, COMPUTED_COLOR)
    assert ratio_traces_with_color(figure, "#00aa00")  # published ratios still drawn

    # With the baseline present, the recomputed run gets its ratio curves.
    computed_both = [plots.Series("B0 = 0", first, "#999999")] + computed
    figure = plots.profile_and_ratio_figure(computed_both + published, baseline_label="B0 = 0")
    assert ratio_traces_with_color(figure, COMPUTED_COLOR)
