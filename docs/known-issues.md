# Known issues and conventions

Limits and traps we live with. Bugs that will be fixed go in `open-items.md` instead.

## Ė in the paper's Table 1 differs from the code

The published Figs. 6–9 were made with Ė = 4.2e41 erg/s (the paper code's
`dEdt = 10**12/1.543` in its natural units), while Table 1 lists 1e42 erg/s. The results
confirm the code value: v∞ = √(2Ė/Ṁ) is 714 km/s for 4.2e41 erg/s and 1100 km/s for
1e42 erg/s, and the published thermal wind reaches about 650 km/s. All defaults use
4.2e41 erg/s so that the figures are reproduced. To model the Table 1 value, pass
`energy_rate=1e42`.

## Inner/outer equations switch on the Mach number, not at r_sb

The paper describes the inner (source-term) equations for r ≤ r_sb and the outer equations
beyond. The paper code, and this port, use the inner equations while the flow is subsonic
(v² ≤ c\*²) and the outer ones once it is supersonic; the radiative force switches the same
way. Density and Alfvén speed switch at r_sb instead (density for r < r_sb, Alfvén speed for
r ≤ r_sb). For the solutions the solver looks for (sonic point at r_sb) the two descriptions
agree. See `physics.md`.

## First-step fallback at r = 0 (not in the paper code)

Near r = 0 the equations amplify round-off inside one Runge–Kutta step, so the error estimate
of the very first step has a floor of about 1e-12 to 5e-12, close to the 1e-12 tolerance. For
Fig. 7, B₀ = 100 µG, some trials of the central-density search hit that floor and the step size
collapsed. Only when an integration breaks down on its very first step is it repeated with that
first step accepted at 1e-10 (100 × the tolerance); every other step, and every integration that
works with the paper code's rule, is unchanged. With it all 14 published runs are reproduced
(central densities to 1e-11). Other inputs could in principle meet the same floor later than the
first step; none of the published runs do.

## Cooling is experimental and unvalidated

`solve_wind(..., cooling=True)` adds radiative cooling with a Cloudy cooling function (solar
metallicity). The published runs never used cooling, and this option has not been compared with
the paper's cooling results (Fig. 2). Its central-temperature search follows the later
`outflow.py` handover script, not the paper code (whose own search has a branch that returns
nothing). If no cooling solution exists, or the cooled wind breaks down, the wind is solved
without cooling and `cooling_used` is False. Cooling is not offered in the app.

## Two branches of the central-density search are untested

Both are ported as the paper code has them and never occur in the 14 published runs. A decision
on them is open (see `open-items.md`); this is the current behaviour:

1. **A trial whose step size collapses** (the paper code returns the string `'RKfail'`) is
   treated as "the wind turned supersonic", which is how Python 2 compared that string with a
   number.
2. **If even the critical central density gives a sonic radius inside r_sb** (strong radiation
   or weak thermal driving), the paper code keeps the last density of its ×10 loop
   (0.1, 1, 10 … M☉/pc³), which is not a converged solution, and goes on. The port does the
   same. If the resulting profile never becomes supersonic, the solver raises
   `NoSolutionError` ("the flow never becomes supersonic") rather than reporting a stall. If it
   does become supersonic somewhere, the profile is returned, but its sonic point is not at r_sb.

## Stalled winds

A wind that falls back to subsonic beyond 1.1 r_sb (for example in a heavy halo) is returned as
a result with `stalled=True`, ending at the stall radius, not raised as an error. A step size
that collapses beyond 1.1 r_sb is also treated as a stall. Only a wind that was supersonic can
stall; a flow still subsonic beyond 1.1 r_sb that never was supersonic is `NoSolutionError`. The published `fig9_halo_10x` run
stalls at about 44.8 r_sb.

## Safety limits not in the paper code

The central-density search stops after 40 increases by ×10 (`NoSolutionError`). The outward
integration stops after 2 000 000 steps: inside 1.1 r_sb this raises `NoSolutionError`, beyond
it the profile is returned as stalled. Neither limit is reached by the published runs.

## Profiles are sensitive to ρ₀ at the 1e-4 level

Near the sonic point the solution is very sensitive to the central density: changing ρ₀ by
about 1e-12 (the search's own precision) moves the profiles by up to 2e-4 relative. This is why
Brent's method and the paper code's bisection, which agree on ρ₀ to 2e-12, give profiles that
differ by up to 2e-4, and why the match to the published runs cannot be much better than
that. Differences of this size are far below what is visible in a plot.

## Figure tabs on machines with few cores

The figure tabs solve their runs in parallel, one process per CPU core. On a machine or host
with one or two cores (e.g. free cloud hosting) the gain is correspondingly smaller. Where
worker processes cannot be started at all, the runs are solved one after another.

## Repeated radii in the published data

The published profiles repeat some radii (where the step size fell below round-off near r_sb,
and at the stall of `fig9_halo_10x`). `WindProfile.interpolate` keeps the last point at each
radius.
