# Architecture

Two parts: the package `outflow` (in `src/outflow/`, the physics and numerics, usable on its
own) and a Streamlit app (in `app/`) built only on the package's public API. There is no
command-line tool; scripts use the Python API.

## Data flow

```
WindParameters ──► solve_wind ──► boundary.find_central_density ──► boundary.sonic_trial (many)
   (inputs, cgs)       │                     │                              │
                       │                     └── equations.WindEquations ◄──┘
                       ▼                            │ uses physics.* formulas
               solver.integrate_outward ──► integrator.cash_karp_step
                       │
                       ▼
                  WindProfile ──► to_csv / plotting / app/plots.py
                       ▲
paper.load_paper_profile (published runs, package data CSV)
```

1. The user builds `WindParameters` (paper defaults, or `from_astro_units`).
2. `solve_wind` searches for the central density ρ₀ (`boundary`), each trial being an
   integration from r = 0 with `WindEquations` (`equations`) and the Cash–Karp step
   (`integrator`).
3. With ρ₀ known, `integrate_outward` integrates to the outer radius and returns a
   `WindProfile`.
4. The profile is saved as CSV, compared with a published run, or plotted.

## Package modules (`src/outflow/`)

| Module | Role |
|---|---|
| `__init__.py` | public API: `WindParameters`, `solve_wind`, `WindProfile`, `NoSolutionError`, `paper_runs`, `load_paper_profile`, `PaperRun`, `cc85_velocity_profile`. No side effects on import. |
| `constants.py` | physical constants and unit conversions, rounded as in the paper code; the paper code's Ė and ρ_M |
| `parameters.py` | `WindParameters` (frozen dataclass, cgs, validated) and derived rates q, Q_th, Q_CR, v∞, ρ_M |
| `physics.py` | single formulas of the paper (density, Alfvén speed, c\*², luminosity, radiative force and heating, NFW gravity) |
| `equations.py` | `WindEquations`: d(v, P, P_C)/dr, inner and outer equations |
| `integrator.py` | Cash–Karp step, error measure and step-size rule of the paper code |
| `boundary.py` | central values (eqs. 44–45) and the search for ρ₀ (last stage: Brent's method by default, or the paper code's bisection) |
| `solver.py` | `solve_wind` and the outward integration, stall detection |
| `profile.py` | `WindProfile`: arrays, derived quantities, interpolation, comparison, CSV read/write |
| `paper.py` | the 14 published runs (`paper_runs`) and their saved profiles (`load_paper_profile`) |
| `cc85.py` | analytic Chevalier & Clegg (1985) thermal wind, used to test the solver |
| `cooling.py` | experimental radiative cooling: Cloudy Λ(T) table, central temperature search |
| `errors.py` | `NoSolutionError` with a plain-language message and optional technical `details` |
| `plotting.py` | matplotlib figures: Fig. 6-style profiles and Figs. 7–9-style ratio panels |
| `data/` | package data: `cloudy_cooling.csv` and `paper/*.csv` (the 14 published runs) |

Dependencies point one way: `constants` ← `physics` ← `parameters` ← `equations` ←
`boundary` ← `solver`; `profile` sits beside them; `paper`, `cooling`, `plotting` and the app
sit on top.

## App (`app/`)

| Module | Role |
|---|---|
| `streamlit_app.py` | the page: tabs Explore, Fig. 6–9 and About; inputs, run buttons, progress, session runs, CSV download, plain error messages |
| `plots.py` | Plotly figure builders: four profile panels, and profile plus ratio panels for Figs. 7–9, with the published curves as pale bands |

Runs are computed fresh on each button press (no caching). Explore keeps up to 8 runs of the
session in `st.session_state` for overlay. Cooling is not offered in the app.

A figure tab solves its runs in parallel (`solve_runs_in_parallel`: a `ProcessPoolExecutor`
with one worker per run, up to the number of CPU cores), and falls back to solving them one
after another if worker processes cannot be started. Workers re-import the app script, so its
page code runs only under `if __name__ == "__main__"`.

## Other folders

| Path | Role |
|---|---|
| `tests/unit`, `tests/integration`, `tests/regression`, `tests/app` | tests; `slow` marks long runs (all 14 paper runs), `legacy` marks tests that need the git-ignored `legacy/` folder |
| `scripts/convert_paper_data.py` | one-off conversion of the paper code's saved data to the package CSVs (needs `legacy/`) |
| `examples/` | runnable examples: quickstart and a Fig. 6 redraw |
| `docs/` | decisions, plans, open items, physics map, this file, known issues |
| `legacy/` | paper-era code and the paper source; git-ignored, kept locally as reference only |
