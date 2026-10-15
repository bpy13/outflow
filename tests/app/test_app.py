"""The Streamlit app, driven headless with Streamlit's AppTest."""

import json
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from outflow import WindParameters
from outflow.solver import DEFAULT_OUTER_RADIUS_IN_SB

APP = Path(__file__).resolve().parents[2] / "app" / "streamlit_app.py"
WIDGET_TYPES = (
    "button",
    "checkbox",
    "color_picker",
    "date_input",
    "download_button",
    "multiselect",
    "number_input",
    "radio",
    "select_slider",
    "selectbox",
    "slider",
    "text_area",
    "text_input",
    "time_input",
    "toggle",
)


def start_app() -> AppTest:
    app = AppTest.from_file(str(APP), default_timeout=30)
    app.run()
    assert not app.exception, app.exception
    return app


def plotted_traces(app: AppTest, chart_index: int = 0) -> list[dict]:
    """The Plotly traces of one chart, read from the figure JSON the app sent."""
    chart = app.get("plotly_chart")[chart_index]
    return json.loads(chart.proto.spec)["data"]


def all_widgets(app: AppTest):
    for widget_type in WIDGET_TYPES:
        yield from app.get(widget_type)


def test_app_loads_with_six_tabs():
    app = start_app()
    labels = [tab.label for tab in app.tabs]
    assert len(labels) == 6
    assert labels[0] == "Explore" and labels[-1] == "About"


def test_explore_defaults_are_the_paper_defaults():
    app = start_app()
    expected = WindParameters.paper_default().to_astro_units()
    for name, value in expected.items():
        assert app.number_input(key=f"explore_{name}").value == value, name
    outer = app.number_input(key="explore_outer_radius_in_sb").value
    assert outer == DEFAULT_OUTER_RADIUS_IN_SB
    built = WindParameters.from_astro_units(
        **{name: app.number_input(key=f"explore_{name}").value for name in expected}
    )
    for name, value in WindParameters.paper_default().as_dict().items():
        assert getattr(built, name) == pytest.approx(value, rel=1e-12), name


def test_no_cooling_widget():
    app = start_app()
    labels = [widget.label for widget in all_widgets(app)]
    assert labels, "expected some widgets"
    assert not [label for label in labels if "cooling" in label.lower()]


def test_explore_run_shows_plot_and_download():
    app = start_app()
    app.number_input(key="explore_outer_radius_in_sb").set_value(5.0)
    app.button(key="explore_run").click()
    app.run(timeout=180)
    assert not app.exception, app.exception
    assert not app.error, [e.value for e in app.error]
    traces = plotted_traces(app)
    velocity = [t for t in traces if t["meta"]["panel"] == "v"]
    assert len(velocity) == 1
    assert velocity[0]["name"].startswith("Run 1: paper defaults")
    downloads = app.get("download_button")
    assert len(downloads) == 1
    # Still no cooling widget once results are shown.
    assert not [w.label for w in all_widgets(app) if "cooling" in w.label.lower()]


def test_no_solution_shows_plain_error_with_details(monkeypatch):
    import outflow
    from outflow import NoSolutionError

    def no_solution(params, outer_radius=None, cooling=False):
        raise NoSolutionError("No stationary wind solution was found.", details="step collapsed")

    # The app imports solve_wind from the package on every run, so patch it there.
    monkeypatch.setattr(outflow, "solve_wind", no_solution)
    app = start_app()
    app.button(key="explore_run").click()
    app.run()
    assert not app.exception, app.exception
    assert [e.value for e in app.error] == ["No stationary wind solution was found."]
    assert app.expander[0].label == "Technical details"
    assert not app.get("plotly_chart")


@pytest.mark.slow
def test_stalled_wind_shows_a_notice():
    app = start_app()
    heavy_halo = 10 * WindParameters.paper_default().to_astro_units()["halo_mass_msun"]
    app.number_input(key="explore_halo_mass_msun").set_value(heavy_halo)
    app.button(key="explore_run").click()
    app.run(timeout=300)
    assert not app.exception, app.exception
    assert any("stalls" in notice.value for notice in app.info)


def test_invalid_inputs_disable_run_with_a_plain_message():
    app = start_app()
    app.number_input(key="explore_cr_fraction").set_value(0.5)
    app.number_input(key="explore_thermal_fraction").set_value(0.8)
    app.run()
    assert not app.exception
    assert any("more than 1" in w.value for w in app.warning)
    assert app.button(key="explore_run").disabled
