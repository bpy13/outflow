"""Streamlit app for the galactic outflow model of Yu et al. 2020 (MNRAS 492, 3179), §2.4.

Run it with::

    streamlit run app/streamlit_app.py

Tabs: Explore (any inputs), Figs. 6-9 (the paper's run sets, recomputed and compared with the
published curves) and About. Nothing is cached: every run is computed fresh when its button is
pressed. Runs made in this session are kept in ``st.session_state`` so they can be overlaid.
Radiative cooling is deliberately not offered here (it is an experimental API option only).
"""

import math
import tempfile
import traceback
from dataclasses import dataclass
from pathlib import Path

import plots
import streamlit as st

from outflow import (
    NoSolutionError,
    PaperRun,
    WindParameters,
    WindProfile,
    load_paper_profile,
    paper_runs,
    solve_wind,
)
from outflow.solver import DEFAULT_OUTER_RADIUS_IN_SB

TAB_NAMES = (
    "Explore",
    "Fig. 6 · Driving mechanisms",
    "Fig. 7 · Magnetic field",
    "Fig. 8 · Opacity",
    "Fig. 9 · Gravity",
    "About",
)
# Explore keeps at most this many runs, one per colour of the package's palette.
MAX_KEPT_RUNS = len(plots.SERIES_COLORS)


# --- running the model -------------------------------------------------------------------------


@dataclass
class RunOutcome:
    """The result of one model run: a profile, or a plain-language error with details."""

    profile: WindProfile | None = None
    error: str = ""
    details: str = ""


def solve_safely(build_params, outer_radius_in_sb: float) -> RunOutcome:
    """Build the parameters and solve the wind; turn every failure into a plain message.

    ``build_params`` is a function returning :class:`WindParameters`, so that invalid inputs
    are reported the same way as solver failures.
    """
    try:
        params = build_params()
        outer_radius = outer_radius_in_sb * params.starburst_radius
        return RunOutcome(profile=solve_wind(params, outer_radius=outer_radius))
    except NoSolutionError as error:
        return RunOutcome(error=error.message, details=error.details)
    except ValueError as error:
        return RunOutcome(error="Please check the inputs.", details=str(error))
    except Exception:  # noqa: BLE001 - shown to the user as a message, never as a stack trace
        return RunOutcome(
            error="Something unexpected went wrong while solving the wind equations.",
            details=traceback.format_exc(),
        )


def show_error(error: str, details: str):
    """A plain error message, with the technical details folded away."""
    st.error(error)
    if details:
        with st.expander("Technical details"):
            st.code(details, language=None)


def stall_notice(label: str, profile: WindProfile) -> str:
    """Plain explanation of a stalled wind (a result, not an error)."""
    radius = profile.stall_radius / profile.params.starburst_radius
    return (
        f"{label}: the wind stalls at r ≈ {radius:.1f} r_sb. Beyond the starburst it slowed "
        "back below the sound speed, so there is no steady supersonic wind further out and the "
        "profile ends there."
    )


def profile_csv(profile: WindProfile) -> str:
    """The profile in the package's CSV format (parameter header plus one row per radius)."""
    with tempfile.TemporaryDirectory() as folder:
        path = profile.to_csv(Path(folder) / "profile.csv")
        return path.read_text(encoding="utf-8")


# --- Explore tab: inputs -----------------------------------------------------------------------


@dataclass(frozen=True)
class InputSpec:
    """One model input as shown in the Explore tab.

    ``name`` is the keyword of :meth:`WindParameters.from_astro_units`; ``symbol`` is used in
    run labels such as "B₀ = 5 µG".
    """

    name: str
    label: str
    symbol: str
    unit: str
    number_format: str
    step: float
    help: str


INPUT_GROUPS = {
    "Starburst": (
        InputSpec(
            "starburst_radius_pc",
            "Starburst radius r_sb [pc]",
            "r_sb",
            "pc",
            "%.1f",
            10.0,
            "Radius of the region where mass and energy are injected uniformly.",
        ),
        InputSpec(
            "energy_rate_erg_s",
            "Energy injection rate Ė [erg/s]",
            "Ė",
            "erg/s",
            "%.3e",
            1e40,
            "Total energy injected per second, split between gas heating, cosmic rays and "
            "radiation. The published figures used 4.2e41 erg/s (see About).",
        ),
        InputSpec(
            "mass_rate_msun_yr",
            "Mass injection rate Ṁ [M☉/yr]",
            "Ṁ",
            "M☉/yr",
            "%.2f",
            0.1,
            "Gas mass injected per year by supernovae and stellar winds.",
        ),
    ),
    "Energy split": (
        InputSpec(
            "thermal_fraction",
            "Thermal fraction",
            "thermal",
            "",
            "%.2f",
            0.05,
            "Share of Ė that heats the gas. Must be above zero.",
        ),
        InputSpec(
            "cr_fraction",
            "Cosmic-ray fraction",
            "CR",
            "",
            "%.2f",
            0.05,
            "Share of Ė given to cosmic rays. The rest of Ė is radiation.",
        ),
    ),
    "Radiation and magnetic field": (
        InputSpec(
            "opacity_cm2_g",
            "Opacity κ [cm²/g]",
            "κ",
            "cm²/g",
            "%.3g",
            1e3,
            "How strongly the wind absorbs the starburst's radiation. Zero switches radiation "
            "pressure off.",
        ),
        InputSpec(
            "magnetic_field_microgauss",
            "Magnetic field B₀ [µG]",
            "B₀",
            "µG",
            "%.1f",
            5.0,
            "Field strength at r_sb. It sets the Alfvén speed, which couples cosmic rays to the "
            "gas. Zero means cosmic rays move with the gas.",
        ),
    ),
    "Dark-matter halo (NFW)": (
        InputSpec(
            "halo_mass_msun",
            "Halo mass M_tot [M☉]",
            "M_tot",
            "M☉",
            "%.3e",
            1e10,
            "Total halo mass within the virial radius. Zero switches gravity off.",
        ),
        InputSpec(
            "nfw_scale_radius_kpc",
            "Scale radius R_s [kpc]",
            "R_s",
            "kpc",
            "%.1f",
            1.0,
            "Radius where the halo's density profile bends (NFW profile).",
        ),
        InputSpec(
            "virial_radius_kpc",
            "Virial radius R_vir [kpc]",
            "R_vir",
            "kpc",
            "%.1f",
            10.0,
            "Outer edge of the halo; with R_s it sets how concentrated the halo is.",
        ),
    ),
}
ALL_INPUTS = tuple(spec for group in INPUT_GROUPS.values() for spec in group)
OUTER_RADIUS_KEY = "explore_outer_radius_in_sb"


def widget_key(spec: InputSpec) -> str:
    return f"explore_{spec.name}"


def default_inputs() -> dict:
    """Paper defaults in astronomer's units, keyed by input name."""
    return WindParameters.paper_default().to_astro_units()


def reset_explore_inputs():
    """Put every Explore input back to the paper default (a button callback)."""
    for name, value in default_inputs().items():
        st.session_state[f"explore_{name}"] = value
    st.session_state[OUTER_RADIUS_KEY] = DEFAULT_OUTER_RADIUS_IN_SB


def init_explore_state():
    """First visit: fill the inputs with the paper defaults and start an empty run list."""
    if OUTER_RADIUS_KEY not in st.session_state:
        reset_explore_inputs()
    st.session_state.setdefault("explore_runs", [])
    st.session_state.setdefault("explore_run_count", 0)
    st.session_state.setdefault("explore_error", None)


def input_widgets() -> dict:
    """Draw the inputs, grouped in columns. Returns the values keyed by input name."""
    values = {}
    columns = st.columns(len(INPUT_GROUPS))
    for column, (group, specs) in zip(columns, INPUT_GROUPS.items(), strict=True):
        with column:
            st.markdown(f"**{group}**")
            for spec in specs:
                values[spec.name] = st.number_input(
                    spec.label,
                    min_value=0.0,
                    step=spec.step,
                    format=spec.number_format,
                    key=widget_key(spec),
                    help=spec.help,
                )
            if group == "Energy split":
                radiation = 1.0 - values["thermal_fraction"] - values["cr_fraction"]
                st.caption(f"Radiation fraction (the remainder): {radiation:.2f}")
    return values


def input_problems(values: dict, outer_radius_in_sb: float) -> list[str]:
    """Plain-language problems that would stop a run, checked before the Run button."""
    problems = []
    for spec in ALL_INPUTS:
        must_be_positive = spec.name in (
            "starburst_radius_pc",
            "energy_rate_erg_s",
            "mass_rate_msun_yr",
            "thermal_fraction",
            "nfw_scale_radius_kpc",
            "virial_radius_kpc",
        )
        if must_be_positive and not values[spec.name] > 0:
            problems.append(f"{spec.label} must be greater than zero.")
    if values["thermal_fraction"] + values["cr_fraction"] > 1 + 1e-12:
        problems.append("The thermal and cosmic-ray fractions add up to more than 1.")
    if not outer_radius_in_sb > 0:
        problems.append("The outer radius must be greater than zero.")
    return problems


def describe_changes(values: dict) -> str:
    """Short text listing the inputs that differ from the paper defaults."""
    defaults = default_inputs()
    changes = [
        f"{spec.symbol} = {spec.number_format % values[spec.name]} {spec.unit}".rstrip()
        for spec in ALL_INPUTS
        if not math.isclose(values[spec.name], defaults[spec.name], rel_tol=1e-9, abs_tol=0.0)
    ]
    return ", ".join(changes) if changes else "paper defaults"


# --- Explore tab: runs and results -------------------------------------------------------------


@dataclass(frozen=True)
class ExploreRun:
    """One kept Explore run: its number, label, colour, profile and CSV text."""

    number: int
    label: str
    color: str
    profile: WindProfile
    csv: str


def free_color(runs: list[ExploreRun]) -> str:
    """The first palette colour not used by a kept run (so a run keeps its colour)."""
    used = {run.color for run in runs}
    return next(color for color in plots.SERIES_COLORS if color not in used)


def keep_run(profile: WindProfile, values: dict):
    """Add a finished run to the session's list, dropping the oldest beyond the limit."""
    runs = st.session_state.explore_runs
    if len(runs) >= MAX_KEPT_RUNS:
        runs.pop(0)
    st.session_state.explore_run_count += 1
    number = st.session_state.explore_run_count
    label = f"Run {number}: {describe_changes(values)}"
    profile.metadata["app_run_label"] = label
    runs.append(ExploreRun(number, label, free_color(runs), profile, profile_csv(profile)))


def run_explore(values: dict, outer_radius_in_sb: float):
    """Solve with the current inputs (on the Run button) and record the result."""
    with st.spinner("Solving the wind equations (usually a few seconds)…"):
        outcome = solve_safely(
            lambda: WindParameters.from_astro_units(**values), outer_radius_in_sb
        )
    if outcome.profile is None:
        st.session_state.explore_error = outcome
    else:
        st.session_state.explore_error = None
        keep_run(outcome.profile, values)


def clear_runs():
    st.session_state.explore_runs = []
    st.session_state.explore_error = None


def show_explore_results():
    outcome = st.session_state.explore_error
    if outcome is not None:
        show_error(outcome.error, outcome.details)
    runs = st.session_state.explore_runs
    if not runs:
        st.caption("No runs yet. Set the inputs and press Run; each run is added to the plot.")
        return
    for run in runs:
        if run.profile.stalled:
            st.info(stall_notice(run.label, run.profile))
    log_r = st.toggle("Logarithmic r axis", key="explore_log_r")
    series = [plots.Series(run.label, run.profile, run.color) for run in runs]
    x_max = max(float(run.profile.r[-1] / run.profile.params.starburst_radius) for run in runs)
    st.plotly_chart(plots.profile_figure(series, x_max=x_max, log_r=log_r), key="explore_plot")
    st.markdown("**Download** (CSV: parameters in the header, then r, v, ρ, P, P_C, T, v_A in cgs)")
    for run in runs:
        st.download_button(
            f"{run.label}",
            data=run.csv,
            file_name=f"outflow_run_{run.number}.csv",
            mime="text/csv",
            key=f"explore_download_{run.number}",
            on_click="ignore",
            icon=":material/download:",
        )


def explore_tab():
    init_explore_state()
    st.markdown(
        "Solve the steady wind for any inputs. The defaults are the paper's M82 model. Press "
        "**Run** to compute; each run is kept and overlaid so you can compare them."
    )
    values = input_widgets()
    outer_radius_in_sb = st.number_input(
        "Outer radius of the computation [r_sb]",
        min_value=0.0,
        step=5.0,
        format="%.1f",
        key=OUTER_RADIUS_KEY,
        help="How far out to integrate, in units of the starburst radius.",
    )
    problems = input_problems(values, outer_radius_in_sb)
    for problem in problems:
        st.warning(problem)
    buttons = st.columns([1, 1, 1, 3])
    run_pressed = buttons[0].button(
        "Run", type="primary", disabled=bool(problems), key="explore_run"
    )
    buttons[1].button("Clear all runs", on_click=clear_runs, key="explore_clear")
    buttons[2].button("Reset inputs", on_click=reset_explore_inputs, key="explore_reset")
    if run_pressed:
        run_explore(values, outer_radius_in_sb)
    show_explore_results()


# --- figure tabs (Figs. 6-9) -------------------------------------------------------------------


@dataclass(frozen=True)
class FigureSpec:
    """One of the paper's figures: what it varies and how it is drawn."""

    number: int
    description: str
    x_max: float  # right edge of the r / r_sb axis, as in the published figure
    with_ratios: bool  # Figs. 7-9 add ratio panels normalised to the first run


FIGURES = (
    FigureSpec(
        6,
        "Three winds that put 80 % of the injected energy into one channel: gas heating "
        "(thermal), cosmic rays, or radiation. Computed to 100 r_sb, shown to 5 r_sb.",
        5.0,
        False,
    ),
    FigureSpec(
        7,
        "A cosmic-ray-driven wind with magnetic field B₀ = 0, 5, 50 and 100 µG. The field sets "
        "how fast cosmic rays stream through the gas. Ratios are relative to B₀ = 0.",
        5.0,
        True,
    ),
    FigureSpec(
        8,
        "A radiation-driven wind with opacity κ = 0, 1e3, 1e4 and 1e5 cm²/g, i.e. how strongly "
        "the gas absorbs the starburst's light. Ratios are relative to κ = 0.",
        5.0,
        True,
    ),
    FigureSpec(
        9,
        "A thermally driven wind in dark-matter haloes of 0, 1 and 10 times the M82 mass. "
        "Computed to 50 r_sb, shown to 30 r_sb; the heaviest halo stops the wind. Ratios are "
        "relative to the run without a halo.",
        30.0,
        True,
    ),
)


@dataclass(frozen=True)
class FigureRunResult:
    """One recomputed paper run, with its published profile for the overlay."""

    run: PaperRun
    outcome: RunOutcome
    published: WindProfile


def run_figure(spec: FigureSpec):
    """Recompute every run of one figure, showing "run k of n" while it works."""
    runs = paper_runs(spec.number)
    progress = st.progress(0.0, text=f"Run 1 of {len(runs)}: {runs[0].label}…")
    results = []
    for index, run in enumerate(runs):
        progress.progress(index / len(runs), text=f"Run {index + 1} of {len(runs)}: {run.label}…")
        outcome = solve_safely(lambda run=run: run.params, run.outer_radius_in_sb)
        results.append(FigureRunResult(run, outcome, load_paper_profile(run)))
    progress.empty()
    st.session_state[f"fig{spec.number}_results"] = results


def figure_series(results: list[FigureRunResult], show_published: bool) -> list[plots.Series]:
    """Recomputed curves, then (optionally) the published ones; one colour per paper run."""
    computed, published = [], []
    for index, result in enumerate(results):
        color = plots.series_color(index)
        if result.outcome.profile is not None:
            computed.append(plots.Series(result.run.label, result.outcome.profile, color))
        published.append(plots.Series(result.run.label, result.published, color, published=True))
    return computed + (published if show_published else [])


def show_figure_results(spec: FigureSpec, results: list[FigureRunResult]):
    for result in results:
        outcome = result.outcome
        if outcome.profile is None:
            show_error(f"{result.run.label}: {outcome.error}", outcome.details)
        elif outcome.profile.stalled:
            st.info(stall_notice(result.run.label, outcome.profile))
    show_published = st.toggle(
        "Show the published curves", value=True, key=f"fig{spec.number}_show_published"
    )
    series = figure_series(results, show_published)
    if not series:
        return
    if spec.with_ratios:
        figure = plots.profile_and_ratio_figure(series, x_max=spec.x_max)
    else:
        figure = plots.profile_figure(series, x_max=spec.x_max)
    st.plotly_chart(figure, key=f"fig{spec.number}_plot")
    st.caption(
        "Thin lines: recomputed now. Pale wide bands of the same colour: the published curves "
        "(the paper code's saved runs). Grey shading: inside the starburst region (r < r_sb)."
    )


def figure_tab(spec: FigureSpec):
    runs = paper_runs(spec.number)
    st.markdown(spec.description)
    pressed = st.button(
        f"Run the {len(runs)} cases of Fig. {spec.number}",
        type="primary",
        key=f"fig{spec.number}_run",
        help="Each case takes a few seconds; the whole figure about half a minute.",
    )
    if pressed:
        run_figure(spec)
    results = st.session_state.get(f"fig{spec.number}_results")
    if results is None:
        st.caption("Press the button to recompute the figure and compare it with the paper.")
        return
    show_figure_results(spec, results)


# --- page --------------------------------------------------------------------------------------


def main():
    st.set_page_config(page_title="Galactic outflow model", layout="wide")
    st.title("Galactic outflows driven by heat, cosmic rays and radiation")
    st.caption("The semi-analytic wind model of Yu et al. 2020 (MNRAS 492, 3179), section 2.4.")
    tabs = st.tabs(TAB_NAMES)
    with tabs[0]:
        explore_tab()
    for tab, spec in zip(tabs[1:5], FIGURES, strict=True):
        with tab:
            figure_tab(spec)
    with tabs[5]:
        st.caption("Coming soon.")


main()
