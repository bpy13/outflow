# Plan: refactor the paper model into a package + Streamlit app

**Goal.** Turn the paper code for Yu et al. 2020 (MNRAS 492, 3179; arXiv:2001.04384), §2.4 /
Figs. 6–9, into a public, modular, tested Python package `outflow` with a Streamlit app, whose
results match the published figures.

**Approach.** Port the physics and numerics of `legacy/Brian2019MNRAS/version1/CRline.py`
(the code that made the published Figs. 6–9) into small modules in cgs units. Convert the
paper's saved runs (`version1/CRline.npy`) to CSV package data, and use them as the reference
for regression tests and app overlays. Build the app on top of the public API.

**Constraints** (from `docs/decisions.md`, verbatim where short):
- All defaults follow the paper code: r_sb = 200 pc, Ė = 4.2e41 erg/s, Ṁ = 2.6 M☉/yr,
  κ = 1e4 cm²/g, B₀ = 50 µG, R_s = 14.7 kpc, R_vir = 164 kpc, ρ_M = 8.81e-3 M☉/pc³,
  split (thermal, radiation, CR) = (0.8, 0.1, 0.1); constants μ = 1.4 × 1.67e-24 g,
  k_B = 1.38e-16, G = 6.67e-8, c = 2.998e10, M☉ = 2e33 g, pc = 3.086e18 cm, yr = 365.25 d.
- The numerical method is the paper code's: Cash–Karp RK4(5), tol 1e-12, initial step 1e-5 pc,
  CRline's ρ₀ search, inner equations while v² ≤ c*², stop on stall beyond 1.1 r_sb.
- Cooling is off by default, experimental in the API, and absent from the UI.
- No CLI. App: Streamlit tabs Explore / Fig. 6 / Fig. 7 / Fig. 8 / Fig. 9 / About. No caching;
  runs start on a button press; runs from the current session are kept for overlay.
- `legacy/` is git-ignored and never edited. Public repo `bpy13/outflow`, MIT licence.
- Full list and the builder-check calls (run map, CSV layout, API names, comparison method):
  `docs/decisions.md`.

**Acceptance criteria**
- A1. `pip install -e ".[dev]"` works in the conda env `outflow`; `import outflow` has no side
  effects (no file reads, no plots).
- A2. All 14 published runs are reproduced: relative difference ≤ 1e-3 in v, ρ, P, P_C, T, v_A
  from 0.2 r_sb to the smaller outer radius. T = μP/(k_B ρ) at every point.
- A3. Thermal-only limit matches analytic CC85 (eq. 8) to ≤ 1e-4.
- A4. The app runs locally. All six tabs work with paper defaults and have no cooling switch.
  Figure tabs overlay the published curves. Explore has session overlay, clear and CSV
  download. Errors show plain messages.
- A5. README, `docs/physics.md`, `docs/architecture.md`, `CITATION.cff`, `LICENSE` exist;
  `legacy/` is not in git.
- A6. `ruff check .` and `ruff format --check .` are clean; `pytest` passes locally (including
  `-m slow`); GitHub Actions passes on push.

Commands (env `outflow`): tests `pytest`, slow set `pytest -m slow`, lint
`ruff check . ; ruff format --check .`, app `streamlit run app/streamlit_app.py`.

---

## Phase 1: Setup and reference data (A1, A5 partly)

- [x] **Task 1: Environment and package skeleton**
  - Files: create `pyproject.toml`, `src/outflow/__init__.py` (`__version__ = "0.1.0"`),
    `tests/unit/test_import.py`
  - Done when: `conda run -n outflow python -c "import outflow; print(outflow.__version__)"`
    prints `0.1.0`, and `conda run -n outflow pytest` passes.
  - Notes: conda env `outflow`, Python 3.12 (python-env skill). Deps: numpy, scipy, matplotlib,
    streamlit, plotly; extra `dev`: pytest, ruff. pytest config: `testpaths = ["tests"]`,
    markers `slow` and `legacy` (needs `legacy/`, skipped when absent), `addopts = "-m 'not slow'"`.
    ruff: line length 100. Package data: `src/outflow/data/**/*.csv`.

- [x] **Task 2: Move paper-era files into `legacy/`**
  - Files: move `Brian2019MNRAS/`, `outflow.py`, `Cloudy.npy`, `Brian2019MNRAS.pdf`,
    `Brian2019MNRAS.tex`, `Brian2019MNRAS.bib` into `legacy/`
  - Done when: `Test-Path legacy\Brian2019MNRAS\version1\CRline.npy` is True and
    `git status --porcelain --ignored` lists `legacy/` only as ignored.
  - Notes: move, never delete. Do not edit anything inside `legacy/`.

- [x] **Task 3: Convert the cooling table**
  - Files: create `scripts/convert_paper_data.py` (cooling part), `src/outflow/data/cloudy_cooling.csv`,
    `tests/integration/test_package_data.py`
  - Done when: the script runs, and a `legacy`-marked test confirms the CSV values equal
    `legacy/Cloudy.npy` exactly (written with `%.17g`).
  - Notes: columns `T_K,Lambda_erg_cm3_s`; header line `# Cloudy cooling function, solar
    metallicity (Yu et al. 2020, §2.1.2). Source: func_map.dta`.

- [ ] **Task 4: Convert the 14 published runs**
  - Files: extend `scripts/convert_paper_data.py`; create 14 CSVs in `src/outflow/data/paper/`;
    extend `tests/integration/test_package_data.py`
  - Done when: 14 files exist. The tests (which need no legacy) check each has the 7 columns,
    a parameter header, and T = μP/(k_B ρ) to 1e-9.
  - Notes: file names `fig6_thermal, fig6_cosmic_ray, fig6_radiation, fig7_B0_0uG,
    fig7_B0_5uG, fig7_B0_50uG, fig7_B0_100uG, fig8_kappa_0, fig8_kappa_1e3, fig8_kappa_1e4,
    fig8_kappa_1e5, fig9_halo_0x, fig9_halo_1x, fig9_halo_10x` (= npy runs 0–13). Run settings:
    the table in decisions.md. npy layout per run: r/r_sb, v [km/s], P, P_C, ρ (CRline natural
    units), T [K], v_A [km/s]. Conversions: r × 200 × 3.086e18; v, v_A × 1e5;
    ρ × 2/(3.086³ × 1e21); P, P_C × 2/(3.086³ × 1e11). Load with `allow_pickle=True`.
    CSV layout as in decisions.md (`%.10e`, `# key = value` header). Header parameters in cgs.

## Phase 2: The package (A1, A2, A3)

Write each test before the code it checks. Every physics function's docstring cites the paper
equation. The cross-checks against legacy (marked `legacy`) import the legacy functions by
executing `version1/CRline.py`'s function block only: the code between its first `'''` and
`# setting up parameter space`. Never run its run loop, which takes minutes.

- [ ] **Task 5: Constants and parameters**
  - Files: `src/outflow/constants.py`, `src/outflow/parameters.py`, `tests/unit/test_parameters.py`
  - Done when: tests pass for paper defaults, the derived q, Q_th, Q_CR, v∞ (714 km/s ± 1),
    `replace()`, the `from_astro_units()` round trip, and rejection of negative values or
    fractions with thermal + CR > 1.
  - Notes: frozen dataclass `WindParameters` with fields `starburst_radius, energy_rate,
    mass_rate, thermal_fraction, cr_fraction, opacity, magnetic_field, halo_mass,
    nfw_scale_radius, virial_radius` (cgs). Property `radiation_fraction`. Default `halo_mass`
    is derived from ρ_M = 8.81e-3 M☉/pc³.

- [ ] **Task 6: Physics functions**
  - Files: `src/outflow/physics.py`, `tests/unit/test_physics.py`
  - Done when: tests pass for:
    - the NFW series and the closed form agreeing at r/R_s = 1e-3;
    - F_rad continuous across r_sb;
    - c*² finite at r = 0 via the L'Hôpital limit;
    - v_A and ρ switching at r_sb as the paper code does;
    - a `legacy` test matching CRline's `fgrav`, `cs2s` and `uAs` after unit conversion, to 1e-12.

- [ ] **Task 7: Equations and integrator**
  - Files: `src/outflow/equations.py`, `src/outflow/integrator.py`, `tests/unit/test_integrator.py`,
    `tests/unit/test_equations.py`
  - Done when: tests pass for:
    - Cash–Karp on y' = y over [0, 1] with error < 1e-10;
    - the error norm ignoring zero components regardless of their position;
    - a `legacy` test matching CRline `CC85()` derivatives at 5 sample states (subsonic,
      supersonic, r = 0 with du, B₀ = 0, κ = 0) to 1e-10.

- [ ] **Task 8: Boundary conditions, solver and profile**
  - Files: `src/outflow/boundary.py`, `src/outflow/solver.py`, `src/outflow/profile.py`,
    `src/outflow/paper.py` (`paper_runs`, `load_paper_profile`), `tests/regression/test_paper_runs.py`
  - Done when:
    - `pytest tests/regression -k fig6_thermal` passes, comparing to 10 r_sb with tolerance 1e-3;
    - a test checks T = μP/(k_B ρ);
    - a test that halts on stall returns `stalled=True` (fig9_halo_10x to 50 r_sb with
      a shortened check is fine as a slow test);
    - the `to_csv`/`from_csv` round trip is exact.
  - Notes: `solve_wind(params, outer_radius=None, cooling=False) -> WindProfile`. Default
    outer radius 50 r_sb. `WindProfile` fields `r, v, rho, P, Pc, T, vA` + `params`, `stalled`,
    `cooling_used`; properties `n`, `ram_pressure`, `sound_speed`. Error type
    `NoSolutionError(message)` with a plain-language message.

- [ ] **Task 9: All 14 paper runs and the analytic CC85**
  - Files: extend `tests/regression/test_paper_runs.py` (parametrised, `slow`),
    `src/outflow/cc85.py`, `tests/regression/test_cc85.py`
  - Done when: `pytest -m slow` passes for all 14 at ≤ 1e-3 (comparison as in decisions.md).
    The CC85 test (split 1,0,0; κ, B₀, M = 0) matches eq. 8 to ≤ 1e-4. Record the measured max
    differences in the phase report.

- [ ] **Task 10: Experimental cooling**
  - Files: `src/outflow/cooling.py`, extend `boundary.py`, `tests/integration/test_cooling.py`
  - Done when: `solve_wind(WindParameters.paper_default(), cooling=True)` returns finite
    profiles with `cooling_used` reported, and the spline matches legacy `plotcool` values to
    1e-12 (`legacy` test).
  - Notes: T₀ search follows `legacy/outflow.py`'s `initial_condition` (bounded maximisation of
    sonic radius, then brentq). If no cooling solution exists, fall back to no cooling and set
    `cooling_used=False`. The table loads lazily on first use.

- [ ] **Task 11: Matplotlib helpers**
  - Files: `src/outflow/plotting.py`, `tests/unit/test_plotting.py`
  - Done when: `plot_profiles([...])` (four panels) and `plot_comparison(reference, others)`
    (ratio panels) save PNGs in a temporary folder without error (Agg backend).

## Phase 3: Streamlit app (A4)

- [ ] **Task 12: App shell and Explore tab**
  - Files: `app/streamlit_app.py`, `app/plots.py` (Plotly figure builders), `tests/app/test_app.py`
  - Done when: a `streamlit.testing.v1.AppTest` test:
    - loads the app with no exception and shows 6 tabs;
    - sets the Explore outer radius to 5 r_sb, clicks Run, and sees a plot and a download button;
    - finds no widget whose label contains "cooling".

    Also, a manual `streamlit run` launch, with a screenshot or checked page, shows the
    Explore tab working.
  - Notes: inputs in astronomer's units: r_sb pc, Ṁ M☉/yr, Ė erg/s, thermal and CR fractions
    (radiation shown as the remainder), κ cm²/g, B₀ µG, M_tot M☉, R_s and R_vir kpc, outer
    radius in r_sb. Runs are kept in `st.session_state` with labels; there is a clear button.
    Panels: v [km/s], log n [cm⁻³], log T [K], log P / P_C / ρv² [erg cm⁻³] against r/r_sb
    (log-r toggle). Stall shows `st.info`; `NoSolutionError` shows `st.error` with the
    message and an expander holding the details.

- [ ] **Task 13: Figure tabs 6–9**
  - Files: extend `app/streamlit_app.py` and `app/plots.py`; extend `tests/app/test_app.py`
    (`slow`)
  - Done when: an AppTest clicks Run in the Fig. 7 tab and finds 4 runs plotted plus 4
    published overlays (slow test). A manual launch checks Fig. 6 against the paper by eye.
  - Notes: run sets from `paper_runs(figure)`. Progress shows "run k of n". The published
    curves are dashed, with a toggle. Figs. 7–9 add ratio panels (as in the paper) normalised to
    the first run. The x range is 0–5 r_sb (Figs. 6–8) or 0–30 r_sb (Fig. 9), with the
    computation outer radius as in the run map.

- [ ] **Task 14: About tab**
  - Files: extend `app/streamlit_app.py`
  - Done when: an AppTest finds the arXiv and MNRAS links and the Ė note in the About tab.
  - Notes: covers what the model is (paper §2.4, short), its assumptions (steady, spherical,
    inviscid, cooling off), the Ė note (published figures used 4.2e41 erg/s; Table 1 says 1e42),
    and the citation.

## Phase 4: Docs, examples, CI (A5, A6)

- [ ] **Task 15: README, licence, citation**
  - Files: `README.md`, `LICENSE` (MIT, "B. P. Brian Yu"), `CITATION.cff`
  - Done when: `cffconvert --validate` (or a YAML load plus the required keys) passes, and the
    README quickstart code runs as a test (`tests/integration/test_readme.py` extracts and runs
    the first python block with outer radius 5 r_sb).
  - Notes: the README explains what the model does in plain words, then covers install, the
    API quickstart, the app, the paper links (arXiv:2001.04384; MNRAS URL without the guest
    key), the Ė note, and the fact that the paper-era code is not included. DOI: take it from
    the arXiv/MNRAS page; don't guess. If it can't be confirmed, use the URL only.

- [ ] **Task 16: Physics map and architecture docs**
  - Files: `docs/physics.md`, `docs/architecture.md`, `docs/known-issues.md`
  - Done when: every function in `physics.py`, `equations.py` and `boundary.py` appears in
    `physics.md` with its paper equation number, and the architecture doc covers every module
    in `src/outflow/` and `app/`. Check with a short script that greps the function names.
  - Notes: equation numbers come from `legacy/Brian2019MNRAS.tex`. Known issues include: the
    inner/outer switch convention, the Table 1 vs code Ė, and that cooling is unvalidated.

- [ ] **Task 17: Examples**
  - Files: `examples/quickstart.py`, `examples/reproduce_fig6.py`
  - Done when: both run in the env (`python examples/...`) and `reproduce_fig6.py` writes
    `fig6.png`. Add `examples/*.png` to `.gitignore`.

- [ ] **Task 18: CI**
  - Files: `.github/workflows/ci.yml`
  - Done when: after push, `gh run watch` shows success. The workflow runs ubuntu with
    Python 3.12, `pip install -e ".[dev]"`, ruff check and format check, and `pytest` (slow and
    legacy tests excluded or skipped).
