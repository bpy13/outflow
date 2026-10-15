# Open items

As of 2026-10-15. Tick off, or move to `decisions.md` when settled. A bug that will be fixed goes
here; a limit or trap we live with goes in `known-issues.md`. Each item in one place only.

## Blockers

## Waiting on the user

- [x] Ė default = 4.2e41 erg/s, confirmed by the user (2026-10-15).
- [x] Fig. 9 heavy halo: 10× (user, 2026-10-15).
- [x] Paper links (user, 2026-10-15): arXiv:2001.04384; MNRAS 492(3), 3179 (2020),
      https://academic.oup.com/mnras/article/492/3/3179/5697960 (guest-access key in the link
      the user gave is left out of public files).
- [x] Cloudy settings: not needed; provenance note says "Cloudy cooling function, solar
      metallicity (paper §2.1.2, §3.7)". Only the experimental cooling option uses it.

## Pending work

- [ ] (Question for the user, non-blocking; found in Phase 2, Task 8.) Two branches of the
      paper code's ρ₀ search are ported as written but never happen in the 14 paper runs:
      1. A trial whose step size collapses returns the string `'RKfail'`, which the search then
         compares with numbers. Python 3 would crash there; the port treats it as larger than
         any number (what Python 2 did), i.e. "the wind turned supersonic".
      2. If even the critical ("peak") ρ₀ gives a sonic radius inside r_sb (strong radiation
         or weak thermal driving, the paper's r_s < r_sb case), the paper code keeps the last
         ρ₀ of its ×10 loop (0.1, 1, 10 … M☉/pc³), which is not a converged solution, and goes
         on. The port does the same. Default if you say "ok": keep both as is and list them in
         `docs/known-issues.md`; alternative: raise `NoSolutionError` in case 2 so the app
         says "no stationary wind solution" instead of showing an unconverged profile.
      Partly addressed after the pre-release review: a profile that never becomes supersonic
      now raises `NoSolutionError` instead of being reported as a "stall" (e.g. thermal 0.01,
      CR 0, κ = 1e7 cm²/g). Still open: whether to raise for every unconverged case 2, even
      when the profile does become supersonic somewhere. The current behaviour is described in
      `docs/known-issues.md`; update it there if the answer changes the code.

- [ ] (Phase 3, non-blocking.) The app has not been looked at in a real browser: the checks were
      AppTest (every tab, runs, errors, stall notice), a headless `streamlit run` answering on
      `/_stcore/health` and `/`, and the app's Plotly figures for Figs. 6, 7 and 9 redrawn with
      matplotlib. A quick look in a browser (layout, legend, hover text) is still worth doing.

- [x] For Task 8/9 (comparison with the paper runs): the reference radii are non-decreasing,
      not strictly increasing. Every paper run repeats a radius 35–107 times (where the step
      fell below round-off near r_sb, while v, P, ρ still change slightly), and
      `fig9_halo_10x` ends with many points at the same radius, 44.8 r_sb, where it stalls.
      Interpolation onto new radii must handle the repeats (e.g. keep the last point per radius).
      Found in Phase 1, Task 4. Done in Task 8: `WindProfile.interpolate` keeps the last point
      per radius.

- [x] (Phase 4, Task 17, non-blocking, cosmetic; fixed after the pre-release review: curves are
      trimmed to the shown range.) `outflow.plotting.plot_profiles` and
      `plot_comparison` scale the y axes over the whole profile, not only the shown r range.
      With profiles computed to 100 r_sb and `x_max=5` (as in `examples/reproduce_fig6.py`),
      the n, T and pressure panels leave empty space (e.g. the pressure axis reaches 1e-16.5
      although the curves shown stay above 1e-13). The app is not affected (it trims curves to
      the shown range). Fix: trim each curve to r ≤ x_max r_sb before plotting, as the app does.

- [ ] (Phase 4, Task 18, non-blocking.) CI shows a warning that `actions/checkout@v4` and
      `actions/setup-python@v5` target the deprecated Node.js 20 runtime (they still run).
      Bump them to versions built for Node.js 24 when convenient.

## Later

- Port Figs. 1–5 scripts and the §3.5 line-emission work (Figs. 10–11).
- Validate the experimental cooling option against `Silich.py` (Fig. 2).
- Possible erratum on Ė in Table 1 (user's decision).
- Hosting the app on Streamlit Community Cloud (needs user's go-ahead: publishing).
