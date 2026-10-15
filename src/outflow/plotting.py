"""Matplotlib helpers: Fig. 6-style profile panels and Figs. 7-9-style ratio panels.

These make paper-style figures, not pixel copies of the published ones. Each function returns
a :class:`matplotlib.figure.Figure` and never calls ``plt.show()``; save it with
``figure.savefig(...)``. Matplotlib is imported only when a function is called.
"""

import numpy as np

from outflow.profile import WindProfile

# Validated categorical order (colourblind-checked), assigned in order, never cycled.
SERIES_COLORS = (
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#e87ba4",  # magenta
    "#008300",  # green
    "#6250d6",  # violet
    "#e34948",  # red
)
# Pressure panel line styles, as in the paper: P solid, ram pressure dashed, P_C dotted.
PRESSURE_STYLES = (("P", "-", "P"), ("ram_pressure", "--", "ρv²"), ("Pc", ":", "P_C"))
STARBURST_SHADE = "#d9d8d4"  # the region inside r_sb is shaded, as in the paper's figures
TEXT_COLOR = "#52514e"


def _figure(rows: int, columns: int, size):
    from matplotlib.figure import Figure

    figure = Figure(figsize=size, layout="constrained")
    axes = figure.subplots(rows, columns, sharex=True).ravel()
    return figure, axes


def _check_count(profiles, labels, first_is_reference=True):
    if not profiles:
        raise ValueError("give at least one profile")
    if len(profiles) > len(SERIES_COLORS) - 1 + int(first_is_reference):
        raise ValueError("too many profiles for one figure (at most 8 including a reference)")
    if labels is None:
        labels = [f"run {i + 1}" for i in range(len(profiles))]
    if len(labels) != len(profiles):
        raise ValueError("give one label per profile")
    return list(labels)


def _style_axis(axis, x_max: float, ylabel: str):
    axis.axvspan(0, 1, color=STARBURST_SHADE, alpha=0.5, linewidth=0, zorder=0)
    axis.set_xlim(0, x_max)
    axis.set_ylabel(ylabel)
    axis.grid(True, color="#e8e7e3", linewidth=0.6)
    axis.tick_params(colors=TEXT_COLOR)
    for spine in ("top", "right"):
        axis.spines[spine].set_visible(False)


def _x(profile: WindProfile) -> np.ndarray:
    return profile.r / profile.params.starburst_radius


def plot_profiles(profiles, labels=None, x_max: float = 5.0):
    """Four panels like the paper's Fig. 6: v, log n, log T and log pressures against r/r_sb.

    Args:
        profiles: list of :class:`WindProfile` (at most 8).
        labels: legend label per profile.
        x_max: right edge of the r / r_sb axis.
    """
    profiles = list(profiles)
    labels = _check_count(profiles, labels)
    figure, axes = _figure(2, 2, (10, 6.5))
    ylabels = (
        "v [km s⁻¹]",
        "log₁₀ n [cm⁻³]",
        "log₁₀ T [K]",
        "log₁₀ P [erg cm⁻³]",
    )
    for axis, ylabel in zip(axes, ylabels, strict=True):
        _style_axis(axis, x_max, ylabel)
    with np.errstate(divide="ignore"):
        for profile, label, color in zip(profiles, labels, SERIES_COLORS, strict=False):
            x = _x(profile)
            axes[0].plot(x, profile.v / 1e5, color=color, linewidth=1.5, label=label)
            axes[1].plot(x, np.log10(profile.n), color=color, linewidth=1.5)
            axes[2].plot(x, np.log10(profile.T), color=color, linewidth=1.5)
            for name, style, _ in PRESSURE_STYLES:
                axes[3].plot(x, np.log10(getattr(profile, name)), style, color=color, linewidth=1.5)
    axes[0].legend(frameon=False, loc="lower right")
    _fit_pressure_axis(axes[3], profiles)
    _pressure_legend(axes[3])
    for axis in axes[2:]:
        axis.set_xlabel("r / r_sb")
    return figure


def _pressure_legend(axis):
    from matplotlib.lines import Line2D

    handles = [
        Line2D([], [], color=TEXT_COLOR, linestyle=style, label=label)
        for _, style, label in PRESSURE_STYLES
    ]
    axis.legend(handles=handles, frameon=False, loc="best")


def _fit_pressure_axis(axis, profiles):
    """Limit the log-pressure axis to the range outside the centre, where rho v^2 -> 0."""
    values = []
    for profile in profiles:
        outside = profile.r > 0.1 * profile.params.starburst_radius
        for name, _, _ in PRESSURE_STYLES:
            data = np.asarray(getattr(profile, name))[outside]
            data = data[data > 0]
            if data.size:
                values.append(np.log10(data))
    if values:
        joined = np.concatenate(values)
        axis.set_ylim(joined.min() - 0.3, joined.max() + 0.3)


def plot_comparison(reference: WindProfile, others, labels=None, x_max: float = 5.0):
    """Ratio panels like the paper's Figs. 7-9: each profile divided by ``reference``.

    Panels: v, n, T and the pressures (P, rho v^2, P_C), each as profile / reference, with the
    reference interpolated linearly in r onto each profile's radii.

    Args:
        reference: the profile everything is normalised to (the first run in the paper).
        others: profiles to compare (at most 7).
        labels: legend label per profile in ``others``.
        x_max: right edge of the r / r_sb axis.
    """
    others = list(others)
    labels = _check_count(others, labels, first_is_reference=False)
    figure, axes = _figure(2, 2, (10, 6.5))
    ylabels = ("v / v_ref", "n / n_ref", "T / T_ref", "pressure / reference")
    for axis, ylabel in zip(axes, ylabels, strict=True):
        _style_axis(axis, x_max, ylabel)
        axis.axhline(1.0, color=TEXT_COLOR, linewidth=0.8)
    with np.errstate(divide="ignore", invalid="ignore"):
        # Colours follow the entity: the reference takes the first colour, as in plot_profiles.
        for profile, label, color in zip(others, labels, SERIES_COLORS[1:], strict=False):
            x = _x(profile)
            ref = reference.interpolate(profile.r)
            axes[0].plot(x, profile.v / ref["v"], color=color, linewidth=1.5, label=label)
            axes[1].plot(x, profile.rho / ref["rho"], color=color, linewidth=1.5)
            axes[2].plot(x, profile.T / ref["T"], color=color, linewidth=1.5)
            reference_values = {
                "P": ref["P"],
                "Pc": ref["Pc"],
                "ram_pressure": ref["rho"] * ref["v"] ** 2,
            }
            for name, style, _ in PRESSURE_STYLES:
                ratio = getattr(profile, name) / reference_values[name]
                axes[3].plot(x, ratio, style, color=color, linewidth=1.5)
    axes[0].legend(frameon=False, loc="best")
    _pressure_legend(axes[3])
    for axis in axes[2:]:
        axis.set_xlabel("r / r_sb")
    return figure
