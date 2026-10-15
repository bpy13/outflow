"""Recompute the paper's Fig. 6 (thermal, cosmic-ray and radiation-driven winds) and save it.

Run from the repository root:  python examples/reproduce_fig6.py
Writes examples/fig6.png and prints how far each run is from the published profile.
It takes about 20 s (three runs to 100 r_sb, as in the paper).
"""

from pathlib import Path

from outflow import load_paper_profile, paper_runs, solve_wind
from outflow.plotting import plot_profiles

OUTPUT = Path(__file__).with_name("fig6.png")

profiles, labels = [], []
for run in paper_runs(6):
    profile = solve_wind(run.params, outer_radius=run.outer_radius)
    differences = profile.max_relative_differences(load_paper_profile(run))
    worst = max(differences, key=differences.get)
    print(f"{run.label:>12}: largest difference from the paper {differences[worst]:.1e} ({worst})")
    profiles.append(profile)
    labels.append(run.label)

figure = plot_profiles(profiles, labels, x_max=5.0)
figure.savefig(OUTPUT, dpi=150)
print(f"saved {OUTPUT}")
