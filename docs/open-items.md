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

- [ ] For Task 8/9 (comparison with the paper runs): the reference radii are non-decreasing,
      not strictly increasing. Every paper run repeats a radius 35–107 times (where the step
      fell below round-off near r_sb, while v, P, ρ still change slightly), and
      `fig9_halo_10x` ends with many points at the same radius, 44.8 r_sb, where it stalls.
      Interpolation onto new radii must handle the repeats (e.g. keep the last point per radius).
      Found in Phase 1, Task 4.

## Later

- Port Figs. 1–5 scripts and the §3.5 line-emission work (Figs. 10–11).
- Validate the experimental cooling option against `Silich.py` (Fig. 2).
- Possible erratum on Ė in Table 1 (user's decision).
- Hosting the app on Streamlit Community Cloud (needs user's go-ahead: publishing).
