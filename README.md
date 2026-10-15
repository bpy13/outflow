# outflow

A Python package and a small web app for the galactic wind model of
Yu, Owen, Wu & Ferreras (2020), *A hydrodynamical study of outflows in starburst galaxies with
different driving mechanisms*, MNRAS 492, 3179
([arXiv:2001.04384](https://arxiv.org/abs/2001.04384),
[MNRAS](https://academic.oup.com/mnras/article/492/3/3179/5697960),
DOI [10.1093/mnras/staa021](https://doi.org/10.1093/mnras/staa021)).

It reproduces the model of the paper's section 2.4 and its Figs. 6–9.

## What this models

Some galaxies form stars in a burst, in a small central region (here the galaxy M82, with a
"starburst" 200 parsecs across in radius). Exploding stars and stellar winds pour gas and energy
into that region, and the result is a hot wind that blows out of the galaxy at hundreds of
kilometres per second.

The model computes that wind as a function of distance from the centre: its speed, density,
temperature and pressures. It treats the wind as steady (not changing in time) and spherical
(depending on radius only). The injected energy can drive the wind in three ways, in any mix:

- **heat**: hot gas expands and pushes outwards;
- **cosmic rays**: very energetic particles that push on the gas through the magnetic field;
- **radiation**: starlight that pushes on the gas, more so the more opaque the gas is.

The wind also has to climb out of the gravity of the galaxy's dark-matter halo, and if gravity
wins it stalls. Mathematically this is a set of ordinary differential equations in radius,
integrated outward from the centre with a Runge–Kutta method; the one unknown, the central
density, is found by requiring the flow to pass smoothly through the speed of sound at the edge
of the starburst.

## Install

Python 3.12. With conda:

```powershell
conda create -n outflow python=3.12
conda activate outflow
pip install -e ".[dev]"
```

Any other Python 3.12 environment works the same way with `pip install -e ".[dev]"`.

## Quickstart

```python
from outflow import WindParameters, solve_wind

params = WindParameters.paper_default()  # the paper's M82 values, in cgs units
profile = solve_wind(params, outer_radius=5 * params.starburst_radius)

print(f"{len(profile)} radii, speed at 5 r_sb: {profile.v[-1] / 1e5:.0f} km/s")
profile.to_csv("wind.csv")  # r, v, rho, P, Pc, T, vA with a header of the inputs

# Inputs in astronomer's units (pc, erg/s, Msun/yr, cm^2/g, microgauss, Msun, kpc):
cr_driven = WindParameters.from_astro_units(thermal_fraction=0.1, cr_fraction=0.8)
```

The run takes a few seconds. Main entry points:

- `WindParameters` – the inputs (cgs), with `paper_default()`, `from_astro_units()` and `replace()`.
- `solve_wind(params, outer_radius=None)` – returns a `WindProfile`: arrays `r, v, rho, P, Pc,
  T, vA`, derived `n`, `ram_pressure`, `sound_speed`, a `stalled` flag, `to_csv()`/`from_csv()`.
  Default outer radius 50 r_sb. Raises `NoSolutionError` when no steady wind exists.
- `paper_runs(figure)` and `load_paper_profile(run)` – the 14 published runs of Figs. 6–9 and
  their saved profiles.
- `cc85_velocity_profile(...)` – the analytic thermal wind of Chevalier & Clegg (1985).
- `outflow.plotting.plot_profiles(...)` / `plot_comparison(...)` – paper-style matplotlib figures.

See `examples/` for a runnable quickstart and a script that redraws Fig. 6.

## The app

```powershell
streamlit run app/streamlit_app.py
```

Tabs: **Explore** (change any input, overlay runs, download CSV), **Fig. 6–9** (recompute the
paper's run sets and compare them with the published curves) and **About**. Every run is
computed when you press its button. A figure tab solves its runs at the same time, one per CPU
core, so it takes about as long as its slowest run (about 5 s on a multi-core computer).

## How the results compare with the paper

The physics and numerical method are a port of the code that made the published figures. All 14
published runs are recomputed and compared between 0.2 r_sb and the outer radius: the largest
relative differences are about 3e-5 in velocity and 2e-4 in pressure (the tests allow 1e-3). The
thermal-only case matches the analytic Chevalier & Clegg wind to 3e-6.

One step differs from the paper code, for speed: the last stage of the search for the central
density uses Brent's method instead of bisection, to the same precision (central densities agree
to 2e-12). `outflow.boundary.find_central_density(params, method="bisection")` gives the paper
code's search.

**A note on Ė.** The published Figs. 6–9 were made with an energy injection rate
Ė = 4.2e41 erg/s (the value in the paper's code), although the paper's Table 1 lists
1e42 erg/s. The defaults here use 4.2e41 erg/s so that they reproduce the figures.

Radiative cooling exists as an experimental option (`solve_wind(..., cooling=True)`); it is off
by default, not in the app, and not validated against the paper. Other limits are listed in
[`docs/known-issues.md`](docs/known-issues.md).

## Documentation

- [`docs/physics.md`](docs/physics.md) – each paper equation and the function that implements it
- [`docs/architecture.md`](docs/architecture.md) – modules and how data flows through them
- [`docs/known-issues.md`](docs/known-issues.md) – known limits and conventions
- [`docs/decisions.md`](docs/decisions.md) – what was decided and why

## Development

```powershell
pytest                                   # fast tests
pytest -m slow                           # all 14 paper runs and other long tests
ruff check . ; ruff format --check .     # lint
```

The paper-era scripts that produced the published figures are not included in this repository;
the published profiles ship as CSV package data in `src/outflow/data/paper/`. The tests that
compare against those published profiles run everywhere (including CI); a further set of
line-by-line cross-checks against the paper-era functions (marker `legacy`) runs only where
those scripts are present and is skipped otherwise.

## Citation

If you use this code, please cite the paper (see also [`CITATION.cff`](CITATION.cff)):

> Yu B. P. B., Owen E. R., Wu K., Ferreras I., 2020, MNRAS, 492, 3179,
> doi:10.1093/mnras/staa021

## Licence

MIT, see [`LICENSE`](LICENSE).
