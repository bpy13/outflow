"""Plotly figures for the Streamlit app.

Two layouts, both against r / r_sb:

- :func:`profile_figure`: four panels like the paper's Fig. 6 (velocity, number density,
  temperature, pressures).
- :func:`profile_and_ratio_figure`: the same four panels on the left and, on the right, each
  profile divided by the first one, like the paper's Figs. 7-9.

Published curves (from the paper's saved runs) are drawn as wide, pale bands of the same colour
under the recomputed lines, so a match shows as a thin line running inside its band.

Every trace carries ``meta = {"panel", "kind", "label", "ratio"}`` so tests can count what is
plotted.
"""

from dataclasses import dataclass

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from outflow import WindProfile
from outflow.constants import KM_S, MEAN_PARTICLE_MASS
from outflow.plotting import SERIES_COLORS

COMPUTED = "computed"
PUBLISHED = "published"
STARBURST_SHADE = "rgba(150, 148, 140, 0.15)"  # the region inside r_sb, as in the paper
GREY = "#7a7975"


@dataclass(frozen=True)
class Series:
    """One curve set to draw: a profile, its legend label, colour and whether it is published."""

    label: str
    profile: WindProfile
    color: str
    published: bool = False

    @property
    def kind(self) -> str:
        return PUBLISHED if self.published else COMPUTED


@dataclass(frozen=True)
class Quantity:
    """One plotted quantity: which panel, how to get it from a profile, its line style."""

    panel: str
    name: str
    values: object  # callable WindProfile -> array
    dash: str = "solid"


# Panel key -> (y-axis title, log scale?). Order is the reading order of the paper's panels.
PANELS = {
    "v": ("v [km s⁻¹]", False),
    "n": ("n [cm⁻³]", True),
    "T": ("T [K]", True),
    "P": ("pressure [erg cm⁻³]", True),
}
# The same for the ratio panels. Pressure ratios can span decades, so that axis is logarithmic.
RATIO_PANELS = {
    "v": ("v / v of first run", False),
    "n": ("n / n of first run", False),
    "T": ("T / T of first run", False),
    "P": ("pressure / first run", True),
}
# Pressure line styles follow outflow.plotting and the paper: P solid, ρv² dashed, P_C dotted.
QUANTITIES = (
    Quantity("v", "v", lambda p: p.v / KM_S),
    Quantity("n", "n", lambda p: p.n),
    Quantity("T", "T", lambda p: p.T),
    Quantity("P", "P (gas)", lambda p: p.P),
    Quantity("P", "ρv² (ram)", lambda p: p.ram_pressure, "dash"),
    Quantity("P", "P_C (cosmic rays)", lambda p: p.Pc, "dot"),
)


def series_color(index: int) -> str:
    """The colour of the index-th series (the package's colourblind-checked order)."""
    return SERIES_COLORS[index % len(SERIES_COLORS)]


# --- public figure builders -----------------------------------------------------------------


def profile_figure(series: list[Series], x_max: float | None = None, log_r: bool = False):
    """Four panels (2 x 2): v, n, T and the pressures against r / r_sb.

    Args:
        series: what to draw; computed and published curves may be mixed.
        x_max: right edge of the r / r_sb axis (None: fit the data).
        log_r: logarithmic r axis (the point at r = 0 is then left out).
    """
    positions = {"v": (1, 1), "n": (1, 2), "T": (2, 1), "P": (2, 2)}
    figure = make_subplots(
        rows=2,
        cols=2,
        shared_xaxes=True,
        vertical_spacing=0.08,
        horizontal_spacing=0.09,
    )
    for item in series:
        _add_profile_traces(figure, item, positions, x_max, log_r)
    _style_axes(figure, positions, x_max, log_r, ratio=False)
    _fit_pressure_axis(figure, series, positions["P"], x_max)
    _add_pressure_key(figure, positions["P"])
    _finish(figure, series, height=640)
    return figure


def profile_and_ratio_figure(series: list[Series], x_max: float | None = None):
    """Profiles on the left and ratios to the first run on the right (4 rows x 2 columns).

    Recomputed curves are divided by the first recomputed curve; published curves by the first
    published curve. The reference is interpolated linearly in r onto each profile's radii.
    """
    positions = {"v": (1, 1), "n": (2, 1), "T": (3, 1), "P": (4, 1)}
    ratio_positions = {key: (row, 2) for key, (row, _) in positions.items()}
    figure = make_subplots(
        rows=4,
        cols=2,
        shared_xaxes=True,
        vertical_spacing=0.04,
        horizontal_spacing=0.09,
        column_titles=["Profiles", "Ratio to the first run"],
    )
    references = {kind: _first_of_kind(series, kind) for kind in (COMPUTED, PUBLISHED)}
    for item in series:
        _add_profile_traces(figure, item, positions, x_max, log_r=False)
        _add_ratio_traces(figure, item, references[item.kind], ratio_positions, x_max)
    _style_axes(figure, positions, x_max, log_r=False, ratio=False)
    _style_axes(figure, ratio_positions, x_max, log_r=False, ratio=True)
    _fit_pressure_axis(figure, series, positions["P"], x_max)
    _add_pressure_key(figure, positions["P"])
    _finish(figure, series, height=1100)
    return figure


# --- traces ---------------------------------------------------------------------------------


def _x_values(profile: WindProfile) -> np.ndarray:
    return profile.r / profile.params.starburst_radius


def _positive_or_nan(values) -> np.ndarray:
    """Log axes cannot show zero or negative values (e.g. P_C without cosmic rays, ρv² at r = 0)."""
    values = np.asarray(values, dtype=float)
    return np.where(values > 0, values, np.nan)


def _line_style(item: Series, dash: str) -> dict:
    if item.published:
        return {"color": item.color, "width": 7, "dash": dash}
    return {"color": item.color, "width": 1.8, "dash": dash}


def _trace(item: Series, quantity: Quantity, x, y, ratio: bool, show_legend: bool):
    name = f"{item.label} (published)" if item.published else item.label
    return go.Scatter(
        x=x,
        y=y,
        mode="lines",
        name=name,
        legendgroup=f"{item.kind}:{item.label}",
        showlegend=show_legend,
        line=_line_style(item, quantity.dash),
        opacity=0.3 if item.published else 1.0,
        meta={"panel": quantity.panel, "kind": item.kind, "label": item.label, "ratio": ratio},
        hovertemplate=f"{name}<br>{quantity.name}<br>r/r_sb = %{{x:.3g}}<br>%{{y:.4g}}"
        "<extra></extra>",
    )


def _shown(x: np.ndarray, x_max: float | None, log_r: bool = False) -> np.ndarray:
    """Mask of the points to draw: within the r axis (Plotly scales y to all data, even data
    beyond the x range), and without r = 0 on a log r axis."""
    keep = np.ones_like(x, dtype=bool)
    if x_max is not None:
        keep &= x <= x_max
    if log_r:
        keep &= x > 0
    return keep


def _add_profile_traces(figure, item: Series, positions: dict, x_max, log_r: bool):
    x = _x_values(item.profile)
    keep = _shown(x, x_max, log_r)
    for quantity in QUANTITIES:
        y = np.asarray(quantity.values(item.profile), dtype=float)
        if PANELS[quantity.panel][1]:
            y = _positive_or_nan(y)
        row, col = positions[quantity.panel]
        show_legend = quantity.panel == "v"
        figure.add_trace(
            _trace(item, quantity, x[keep], y[keep], ratio=False, show_legend=show_legend),
            row=row,
            col=col,
        )


def _ratio(profile: WindProfile, reference: WindProfile, quantity: Quantity):
    """profile / reference for one quantity, on the profile's radii that the reference covers.

    Returns the ratio and the mask of the profile's radii it belongs to. Where the reference
    value is zero (e.g. v at r = 0) the ratio is left undefined (NaN).
    """
    within = profile.r <= reference.r[-1]
    reference_at_radii = _InterpolatedProfile(reference.interpolate(profile.r[within]))
    mine = np.asarray(quantity.values(profile), dtype=float)[within]
    theirs = np.asarray(quantity.values(reference_at_radii), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(theirs != 0, mine / theirs, np.nan)
    return ratio, within


class _InterpolatedProfile:
    """Interpolated columns that answer like a WindProfile for the quantity functions."""

    def __init__(self, columns: dict):
        self.v = columns["v"]
        self.rho = columns["rho"]
        self.P = columns["P"]
        self.Pc = columns["Pc"]
        self.T = columns["T"]

    @property
    def n(self):
        return self.rho / MEAN_PARTICLE_MASS

    @property
    def ram_pressure(self):
        return self.rho * self.v**2


def _add_ratio_traces(figure, item: Series, reference: Series, positions: dict, x_max):
    x = _x_values(item.profile)
    for quantity in QUANTITIES:
        ratio, within = _ratio(item.profile, reference.profile, quantity)
        keep = _shown(x[within], x_max)
        row, col = positions[quantity.panel]
        figure.add_trace(
            _trace(item, quantity, x[within][keep], ratio[keep], ratio=True, show_legend=False),
            row=row,
            col=col,
        )


def _first_of_kind(series: list[Series], kind: str) -> Series | None:
    return next((item for item in series if item.kind == kind), None)


# --- layout ---------------------------------------------------------------------------------


def _fit_pressure_axis(figure, series: list[Series], position, x_max):
    """Limit the log pressure axis to values outside the centre (ρv² -> 0 at r = 0 would
    otherwise stretch the axis over twenty decades), as ``outflow.plotting`` does."""
    logs = []
    for item in series:
        x = _x_values(item.profile)
        outside = (x > 0.1) & _shown(x, x_max)
        for quantity in QUANTITIES:
            if quantity.panel != "P":
                continue
            values = np.asarray(quantity.values(item.profile), dtype=float)[outside]
            values = values[values > 0]
            if values.size:
                logs.append(np.log10(values))
    if logs:
        joined = np.concatenate(logs)
        row, col = position
        # Plotly takes the range of a log axis in powers of ten.
        figure.update_yaxes(range=[joined.min() - 0.3, joined.max() + 0.3], row=row, col=col)


def _style_axes(figure, positions: dict, x_max, log_r: bool, ratio: bool):
    bottom_row = max(row for row, _ in positions.values())
    for key, (row, col) in positions.items():
        title, log_y = RATIO_PANELS[key] if ratio else PANELS[key]
        figure.update_yaxes(title_text=title, type="log" if log_y else "linear", row=row, col=col)
        if log_r:
            figure.update_xaxes(type="log", row=row, col=col)
        elif x_max is not None:
            figure.update_xaxes(range=[0, x_max], row=row, col=col)
        if row == bottom_row:
            figure.update_xaxes(title_text="r / r_sb", row=row, col=col)
        if not log_r:
            figure.add_vrect(
                x0=0, x1=1, fillcolor=STARBURST_SHADE, line_width=0, layer="below", row=row, col=col
            )


def _add_pressure_key(figure, position):
    """Grey legend entries explaining the three pressure line styles."""
    row, col = position
    for quantity in QUANTITIES:
        if quantity.panel != "P":
            continue
        figure.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines",
                name=quantity.name,
                line={"color": GREY, "dash": quantity.dash, "width": 1.8},
                legendgroup="pressure key",
                legendgrouptitle_text="Pressure panel",
                meta={"panel": "key", "kind": "key", "label": quantity.name, "ratio": False},
            ),
            row=row,
            col=col,
        )


def _finish(figure, series: list[Series], height: int):
    if any(item.published for item in series):
        figure.add_annotation(
            text="Pale wide bands: published curves (Yu et al. 2020)",
            xref="paper",
            yref="paper",
            x=0,
            y=1.0,
            yanchor="bottom",
            showarrow=False,
            font={"size": 12, "color": GREY},
        )
    figure.update_layout(
        height=height,
        margin={"l": 60, "r": 20, "t": 50, "b": 50},
        legend={"orientation": "v", "x": 1.02, "y": 1},
        hovermode="closest",
    )
