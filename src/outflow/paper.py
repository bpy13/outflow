"""The 14 published runs behind Figs. 6-9 of Yu et al. 2020, and their saved profiles.

The run settings were rebuilt from the paper code (``version1/CRline.py``); see "Reference run
map" in ``docs/decisions.md``. All runs share r_sb = 200 pc, E-dot = 4.2e41 erg/s,
M-dot = 2.6 Msun/yr, R_s = 14.7 kpc and R_vir = 164 kpc.
"""

from dataclasses import dataclass
from importlib import resources

from outflow.constants import MICROGAUSS
from outflow.parameters import PAPER_HALO_MASS, WindParameters
from outflow.profile import WindProfile


@dataclass(frozen=True)
class PaperRun:
    """One published run: its name (CSV file stem), figure, legend label and settings."""

    name: str
    figure: int
    label: str
    params: WindParameters
    outer_radius_in_sb: float

    @property
    def outer_radius(self) -> float:
        """Outer radius of the published computation (cm)."""
        return self.outer_radius_in_sb * self.params.starburst_radius


def _run(name, figure, label, thermal, cr, opacity, b0_microgauss, halo_factor, outer):
    params = WindParameters(
        thermal_fraction=thermal,
        cr_fraction=cr,
        opacity=opacity,
        magnetic_field=b0_microgauss * MICROGAUSS,
        halo_mass=halo_factor * PAPER_HALO_MASS,
    )
    return PaperRun(name, figure, label, params, float(outer))


# In the order of the paper code's saved runs (version1/CRline.npy, runs 0-13).
_PAPER_RUNS = (
    _run("fig6_thermal", 6, "Thermal", 0.8, 0.1, 1e4, 50.0, 1.0, 100),
    _run("fig6_cosmic_ray", 6, "Cosmic rays", 0.1, 0.8, 1e4, 50.0, 1.0, 100),
    _run("fig6_radiation", 6, "Radiation", 0.1, 0.1, 1e4, 50.0, 1.0, 100),
    _run("fig7_B0_0uG", 7, "B₀ = 0 µG", 0.1, 0.8, 1e4, 0.0, 1.0, 5),
    _run("fig7_B0_5uG", 7, "B₀ = 5 µG", 0.1, 0.8, 1e4, 5.0, 1.0, 5),
    _run("fig7_B0_50uG", 7, "B₀ = 50 µG", 0.1, 0.8, 1e4, 50.0, 1.0, 5),
    _run("fig7_B0_100uG", 7, "B₀ = 100 µG", 0.1, 0.8, 1e4, 100.0, 1.0, 5),
    _run("fig8_kappa_0", 8, "κ = 0 cm²/g", 0.1, 0.1, 0.0, 50.0, 1.0, 5),
    _run("fig8_kappa_1e3", 8, "κ = 1e3 cm²/g", 0.1, 0.1, 1e3, 50.0, 1.0, 5),
    _run("fig8_kappa_1e4", 8, "κ = 1e4 cm²/g", 0.1, 0.1, 1e4, 50.0, 1.0, 5),
    _run("fig8_kappa_1e5", 8, "κ = 1e5 cm²/g", 0.1, 0.1, 1e5, 50.0, 1.0, 5),
    _run("fig9_halo_0x", 9, "M_tot = 0", 0.8, 0.1, 1e4, 50.0, 0.0, 50),
    _run("fig9_halo_1x", 9, "M_tot = M_tot,M82", 0.8, 0.1, 1e4, 50.0, 1.0, 50),
    _run("fig9_halo_10x", 9, "M_tot = 10 × M_tot,M82", 0.8, 0.1, 1e4, 50.0, 10.0, 50),
)


def paper_runs(figure: int | None = None) -> list[PaperRun]:
    """The published runs of one figure (6, 7, 8 or 9), or all 14 if ``figure`` is None."""
    if figure is None:
        return list(_PAPER_RUNS)
    runs = [run for run in _PAPER_RUNS if run.figure == figure]
    if not runs:
        raise ValueError(f"figure must be 6, 7, 8 or 9, got {figure!r}")
    return runs


def load_paper_profile(run: str | PaperRun) -> WindProfile:
    """The published profile of a run (by name or :class:`PaperRun`), from the package data."""
    name = run.name if isinstance(run, PaperRun) else run
    names = {r.name for r in _PAPER_RUNS}
    if name not in names:
        raise ValueError(f"unknown paper run {name!r}; choose from {sorted(names)}")
    resource = resources.files("outflow") / "data" / "paper" / f"{name}.csv"
    with resources.as_file(resource) as path:
        return WindProfile.from_csv(path)
