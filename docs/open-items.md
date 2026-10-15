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

- [x] For Task 8/9 (comparison with the paper runs): the reference radii are non-decreasing,
      not strictly increasing. Every paper run repeats a radius 35–107 times (where the step
      fell below round-off near r_sb, while v, P, ρ still change slightly), and
      `fig9_halo_10x` ends with many points at the same radius, 44.8 r_sb, where it stalls.
      Interpolation onto new radii must handle the repeats (e.g. keep the last point per radius).
      Found in Phase 1, Task 4. Done in Task 8: `WindProfile.interpolate` keeps the last point
      per radius.

## Later

- Port Figs. 1–5 scripts and the §3.5 line-emission work (Figs. 10–11).
- Validate the experimental cooling option against `Silich.py` (Fig. 2).
- Possible erratum on Ė in Table 1 (user's decision).
- Hosting the app on Streamlit Community Cloud (needs user's go-ahead: publishing).
