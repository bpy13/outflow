# Decisions

What was decided, why, and when. Change these only with a reason, and note the change here
(mark reversed decisions as reversed rather than deleting them). "(Claude's call)" marks
minor decisions made without asking; the user can overturn them.

## Background (findings from kickoff, 2026-10-15)

- The paper's Figs. 6–9 were produced by `Brian2019MNRAS/version1/CRline.py` (data saved in
  `version1/CRline.npy`, 14 runs). The top-level `Brian2019MNRAS/CRline.py` has identical physics
  code but a later run set (heavy halo 5× instead of 10×, log axes in pc).
- `outflow.py` was a later attempt at a handover version. Compared with the paper code it:
  - stores a wrong temperature (`profile()` uses the central density ρ₀ for every T; the
    paper code is correct, verified: T = μP/(k_B ρ) to 6e-16 in all saved runs);
  - switches cooling on automatically when a cooling solution exists (the paper code always
    runs with cooling off), which made the thermal case 5–7 % slower than the paper;
  - uses a looser step tolerance (1e-8 vs 1e-12).
  With cooling forced off and tolerance 1e-12 it matches the paper data to ≤1.5e-4 relative,
  so its equations are a faithful copy of `CRline.py`.
- `Cloudy.npy` is exactly the cooling column of `Brian2019MNRAS/func_map.dta` (Cloudy output),
  parsed as `plotcool.py` does: Λ(T) in erg cm³ s⁻¹, T from 0 to 9.9e9 K, 10 000 rows.
- The paper code's energy injection rate is `dEdt = 10**12/1.543` in its natural units
  (M☉, pc, km/s), which is 4.2e41 erg/s, not the 1e42 erg/s in the paper's Table 1.

## Scope

- Refactor the combined model of paper §2.4 (Figs. 6–9) into a modular, documented, tested,
  public Python package plus a Streamlit app. (user, 2026-10-15)
- Out of scope for now: Figs. 1–5 scripts (other than analytic CC85) and the line-emission work
  of §3.5 (Figs. 10–11). They stay in `legacy/` for a later phase. (user, default accepted)
- Analytic CC85 (paper eqs. 7–9, Fig. 1) is included, mainly as a test of the solver's
  thermal-only limit. (user, default accepted)
- No command-line tool. The scriptable interface is the Python API; the UI is the only app.
  This overrides the "CLI/UI parity" app convention. (user)
- Speed improvements are out of scope. (Claude's call; never raised as a goal)

## Physics and numerics (follow the paper code)

- All defaults follow the paper / paper code. (user)
  - r_sb = 200 pc, Ṁ = 2.6 M☉/yr, κ = 1e4 cm²/g, B₀ = 50 µG, NFW R_s = 14.7 kpc,
    characteristic density ρ_M = 8.81e-3 M☉/pc³ (M_tot ≈ 5.54e11 M☉ within R_vir = 164 kpc),
    energy split (thermal, radiation, CR) = (0.8, 0.1, 0.1).
  - Ė = 4.2e41 erg/s, the value the paper code actually used (see Background). Confirmed by
    the user 2026-10-15. The app's About tab and the README note the difference from Table 1.
  - Physical constants exactly as the paper code: μ = 1.4 × 1.67e-24 g, k_B = 1.38e-16,
    G = 6.67e-8, c = 2.998e10 (cgs), M☉ = 2e33 g, pc = 3.086e18 cm, yr = 365.25 d.
- Heavy halo in the Fig. 9 comparison is 10 × M_tot,M82, as in the paper PDF and `version1`.
  (user: "stick with paper") Confirmed 2026-10-15: the paper's Fig. 7–9 captions (LaTeX source)
  say the r axis is r/r_sb, 0–5 r_sb for earlier plots and 0–30 for Fig. 9, and the legend
  says 10 ×. `version1/34.pdf` has an r/r_sb axis and a "10 ×" label; the top-level `34.pdf`
  has an r/pc axis and a "5 ×" label. So the published figures came from `version1`.
- Ė actually used by the published figures is 4.2e41 erg/s, confirmed from the results:
  v∞ = √(2Ė/Ṁ) is 714 km/s for 4.2e41 but 1100 km/s for 1e42. The paper's thermal wind
  reaches ~650 km/s (Fig. 6, Table 3: 633–651 km/s), which is possible only with 4.2e41.
- Numerical method is the paper code's (`CRline.py`), ported, not redesigned (user: "stick to
  the paper"):
  - Cash–Karp RK4(5) with relative error control, tolerance 1e-12, step update
    h ← 0.9 h (tol/δ)^0.2 on accept, ^0.25 on reject; initial step 1e-5 pc.
  - Central boundary values from paper eqs. 44–45 with the L'Hôpital limits at r = 0
    (dv/dr = q/3ρ₀, dv_A/dr = B₀/(r_sb √(4πρ₀))).
  - ρ₀ found as in `CRline.CRline()`: start at 0.1 (natural units), multiply by 10 until the
    trial solution turns over or reaches r_sb, bisect for the peak, then bisect for r_s = r_sb.
  - Source-term ("inner") equations are used while v² ≤ c*² (subsonic), as the paper code does;
    density and Alfvén speed switch at r_sb (ρ: r < r_sb; v_A: r ≤ r_sb). Documented in
    `docs/physics.md` as the paper code's convention.
  - Outward integration stops at the requested outer radius, or when the flow falls back to
    subsonic beyond 1.1 r_sb (wind stalls, e.g. heavy halo), as in `CRline.py`.
- The error measure must not depend on the order of array elements when a component is zero
  (e.g. P_C ≡ 0 when no CR injection): components with zero value are skipped explicitly.
  Results unchanged for the paper runs. (Claude's call)
- **First-step fallback at r = 0** (Claude's call, 2026-10-15, Phase 2 Task 9). Near r = 0
  the equations are unstable within one Runge–Kutta step (round-off grows ~40× per stage), so
  the error estimate of the first step has a floor of about 1e-12–5e-12, close to the 1e-12
  tolerance. For Fig. 7, B₀ = 100 µG, a few trial densities of the ρ₀ search hit that floor in
  the cgs port (the paper code, with its own round-off, did not), the step size collapsed, and
  the search ended 1.8e-9 away from the paper's ρ₀ and then could not start the outward
  integration. Fix: only when the integration breaks down on its very first step, it is
  repeated accepting that first step at 100 × the tolerance (1e-10). All other steps, and
  every trial that succeeds with the paper code's rule, are unchanged. With it, all 14 ρ₀ agree
  with the paper's to ≤ 1e-11 and all runs match to ≤ 2.2e-4.
- Cooling: off by default. Kept in the API as an experimental option
  (`cooling=True`) using the S04-style T₀ search; not exposed in the UI; not validated against
  the paper (the paper's §3 runs never used it). The paper code's own cooling search
  (`Tsolve`) has a branch that returns nothing, so the search follows `outflow.py`'s
  `initial_condition` logic instead. (user: off by default, no UI switch; Claude's call on
  the method)
- Units: cgs internally (standard for readers from both backgrounds), plus
  `WindParameters.from_astro_units(...)` taking pc, M☉/yr, erg/s, µG, cm²/g, M☉, kpc. No
  astropy dependency. (user accepted default)
- The temperature is computed from each radius's own density (the `outflow.py` bug is not
  carried over).

## Verification

- "Same results" = differences not visible by eye. Test threshold: every reference run matches
  to ≤ 1e-3 relative in v, ρ, P, P_C, T between 0.2 r_sb and the run's outer radius
  (measured agreement of the faithful port is ~1.5e-4). (user: "not catchable by eyeball";
  threshold is Claude's call)
- Reference data: the 14 runs of `version1/CRline.npy`, converted to cgs and stored as one
  CSV per run in `src/outflow/data/paper/` (package data, so the installed package, the app and
  the tests all read the same single copy). The same CSV format is what
  `WindProfile.to_csv()` writes and `WindProfile.from_csv()` reads. (user: not .npy; "formats
  the refactored module reuses"; CSV chosen by Claude as readable by anyone)
- Tests split into `tests/unit`, `tests/integration`, `tests/regression` (app convention). The
  full paper reproduction (14 runs × ~5–10 s) is marked `slow`; the default run includes one
  fast representative case.

## Builder-check calls (Claude's calls, 2026-10-15)

From the builder check on these documents. All minor; the user can overturn any.

- **Reference run map** (rebuilt from `version1/CRline.py`; the .npy stores no settings).
  Shared: r_sb 200 pc, Ė 4.2e41 erg/s, Ṁ 2.6 M☉/yr, R_s 14.7 kpc, R_vir 164 kpc.
  Energy split is (thermal, radiation, CR).
  | Runs | Figure | Varied | Fixed | Outer radius |
  |---|---|---|---|---|
  | 0–2 | Fig. 6 | split (0.8,0.1,0.1) / (0.1,0.1,0.8) / (0.1,0.8,0.1) | κ 1e4, B₀ 50 µG, M82 halo | 100 r_sb |
  | 3–6 | Fig. 7 | B₀ = 0, 5, 50, 100 µG | split (0.1,0.1,0.8), κ 1e4, M82 halo | 5 r_sb |
  | 7–10 | Fig. 8 | κ = 0, 1e3, 1e4, 1e5 cm²/g | split (0.1,0.8,0.1), B₀ 50 µG, M82 halo | 5 r_sb |
  | 11–13 | Fig. 9 | halo ρ_M × 0, 1, 10 | split (0.8,0.1,0.1), κ 1e4, B₀ 50 µG | 50 r_sb |
- **Unit conversion of the reference data** uses the paper code's own conversion factors
  (M☉ = 2e33 g, pc = 3.086e18 cm, km/s = 1e5 cm/s), as in `CRline.py`'s plotting lines.
  v_A is stored in the .npy (`uA_list`), converted and tested like the other columns; zero
  values (B₀ = 0) are skipped in relative comparisons.
- **CSV layout:** comma-separated, ASCII column names with units
  `r_cm, v_cm_s, rho_g_cm3, P_erg_cm3, Pc_erg_cm3, T_K, vA_cm_s`, values `%.10e`, header lines
  `# key = value` (parameters, units, provenance, package version).
- **Conversion script** is committed as `scripts/convert_paper_data.py`. It needs the
  git-ignored `legacy/` folder, says so at the top, and is not run in CI.
- **Halo parameters:** the API takes M_tot, R_s, R_vir (paper Table 1). The default M_tot is
  the value implied by the paper code's ρ_M = 8.81e-3 M☉/pc³ (5.557e11 M☉, shown as 5.56e11),
  so results match the code exactly. The Fig. 9 heavy halo is M_tot × 10 with R_s, R_vir fixed,
  identical to ρ_M × 10.
- **Internal units stay cgs.** Starting values from the paper code are converted: initial
  step 1e-5 pc, first ρ₀ trial 0.1 M☉/pc³. Relative error control does not depend on units;
  any change from floating-point round-off is far below the 1e-3 test threshold.
- **Inputs and defaults:** the energy fractions are thermal and CR; radiation = 1 − thermal − CR,
  which must be ≥ 0. κ = 0, B₀ = 0, M_tot = 0 and a CR fraction of 0 are allowed, which switches
  that physics off. Default outer radius is 50 r_sb (Explore tab); the paper presets use the
  table above.
- **Stalled wind** (falls back to subsonic beyond 1.1 r_sb) is a result, not an error: the
  shorter profile is returned with `stalled=True` and the stall radius; the app shows a
  notice. Errors are raised only when no solution is found.
- **Cooling is built now** as the experimental option. The table is converted from
  `Cloudy.npy` to package CSV. It gets a smoke test only (runs, finite output, reports whether
  cooling was used).
- **Public API:** `WindParameters` (with `.paper_default()`, `.replace()`,
  `.from_astro_units()`), `solve_wind(params, outer_radius=None, cooling=False) -> WindProfile`,
  `WindProfile` (arrays plus derived n, ram pressure, c*; `to_csv`/`from_csv`),
  `paper_runs(figure)` returning the parameter sets in the table above,
  `load_paper_profile(run)`, and `cc85_velocity_profile()` (analytic, eqs. 7–9). Module layout
  as in the kickoff proposal; `docs/architecture.md` is written in the build.
- **Comparison method:** reference values are interpolated linearly in r onto the new radii,
  over 0.2 r_sb to the smaller of the two outer radii.
- **Fast default test:** run 0 (Fig. 6, thermal) to 10 r_sb. CC85 test: thermal only
  (split 1,0,0; κ = 0; B₀ = 0; no halo) against eq. 8, to 1e-4 relative.
- **App, runs within a session:** every run is computed fresh, with no caching. Runs the user
  made in this session are kept in the page's memory so they can be overlaid, and can be
  cleared.
- **App, progress:** a progress indicator shows "run k of n" while a figure tab computes
  (about 20–40 s locally).
- **App, axes:** as in the published figures. x is r/r_sb on a linear axis (0–5 for
  Figs. 6–8, 0–30 for Fig. 9); v in km/s, n in cm⁻³ (log), T in K (log), P in erg cm⁻³ (log).
  The Explore tab has a log-r toggle.
- **If the user moves the default Ė to 1e42:** only the Explore defaults change. The paper
  tabs and the reference runs keep 4.2e41, so they still reproduce the figures.
- **The matplotlib helper** makes a four-panel Fig. 6-style plot and a residual comparison
  plot like Figs. 7–9, not pixel-identical copies of the published figures.
- **Docs written in this build:** `docs/physics.md`, `docs/architecture.md`, a known-issues
  file as issues are found, README and `CITATION.cff`. The plan goes in `docs/plans/` after the
  user says "go".

## Platform and hosting

- Python 3.12, new conda env `outflow`. Runtime deps: numpy, scipy, matplotlib, streamlit,
  plotly. Dev deps: pytest, ruff. (user: new env; deps Claude's call)
- Packaging: `pyproject.toml`, src layout (`src/outflow/`), installable with `pip install -e .`.
- App runs locally with `streamlit run`. Hosting on Streamlit Community Cloud is possible later
  but is publishing, so only on the user's say-so. (user accepted default)

## Data and secrets

- No secrets. No user data stored. The app offers CSV download of what is plotted.
- The Cloudy cooling table ships as package data (CSV, `src/outflow/data/cloudy_cooling.csv`)
  with a provenance note; read on first use.

## UI

- Streamlit, multiple tabs (user):
  - **Explore** – every model input with paper defaults, run button, Fig. 6-style panels
    (velocity, density, temperature, pressures), overlay of several runs, CSV download.
  - **Fig. 6 Driving mechanisms**, **Fig. 7 Magnetic field**, **Fig. 8 Opacity**,
    **Fig. 9 Gravity** – the paper's run sets, recomputed on demand, with the published curves
    (from `tests/reference/`) overlaid for comparison. (overlay is Claude's call)
  - **About** – what the model is, assumptions, citation, the Ė note.
- No cooling switch. (user)
- Always compute on the spot, no caching. Runs start on a button press, not on tab open, so
  opening a tab is instant. (user: no cache; button is Claude's call)
- Plain-language error messages (e.g. "no stationary wind solution for these inputs"); technical
  details in a collapsible section. A full reference-ID/log system is not needed for a
  stateless science app. (Claude's call, adapting the app convention)
- Interactive plots with Plotly in the app; paper-style matplotlib helpers in the package.
  (Claude's call)

## Repository

- Public GitHub repo `bpy13/outflow`, MIT licence; README asks users to cite the paper;
  `CITATION.cff` included. (user)
- `legacy/` holds `Brian2019MNRAS/` (whole folder), `outflow.py`, the paper PDF and the paper's
  LaTeX source (`Brian2019MNRAS.tex`, `.bib`), and is git-ignored. (user) The LaTeX source is
  the reference for equation numbers and wording in `docs/physics.md`; it matches the PDF
  (Table 1 Ė = 1e42, Fig. 9 heavy halo 10×, figure files 211–352).
- Paper links in public files: arXiv:2001.04384 and MNRAS 492(3), 3179 (2020),
  https://academic.oup.com/mnras/article/492/3/3179/5697960 (no guest-access key).
- Cooling-table provenance is a one-line pointer to the paper (Cloudy, solar metallicity); no
  further explanation. (user)
- CI: GitHub Actions running lint and the fast tests on every push. (Claude's call)
