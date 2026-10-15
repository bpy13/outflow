# outflow – notes for Claude

A public, refactored Python package and Streamlit app for the semi-analytic galactic-outflow
model of Yu, Owen, Wu & Ferreras (2020, MNRAS 492, 3179), "A hydrodynamical study of outflows in starburst
galaxies with different driving mechanisms". The model solves the steady, spherically symmetric
wind equations with thermal, radiation and cosmic-ray driving plus an NFW halo (paper §2.4).
Audience: astrophysicists reading the paper and programmers with no astrophysics background.
The user is the paper's first author.

Read these before substantial changes:

- `docs/decisions.md` – what was decided and why (don't re-litigate without a reason)
- `docs/open-items.md` – pending work and questions waiting on the user
- `docs/plans/` – the latest plan
- `docs/physics.md` – paper equation ↔ code map
- `docs/architecture.md` – modules, data flow, where things live
- `docs/known-issues.md` – limits and traps found so far
- `README.md` – user-facing setup and commands

## Environment
- Conda env: `outflow` (Python 3.12)
- Install: `conda activate outflow; pip install -e ".[dev]"`
- Run the app: `streamlit run app/streamlit_app.py`
- Tests: `pytest` (fast); `pytest -m slow` for the full paper-reproduction set
- Lint: `ruff check . ; ruff format --check .`
- Docs coverage: `python scripts/check_docs.py` (every physics/equations/boundary function in
  `docs/physics.md`, every module in `docs/architecture.md`)
- Examples: `python examples/quickstart.py`, `python examples/reproduce_fig6.py` (writes
  `examples/fig6.png`, git-ignored)
- CI: `.github/workflows/ci.yml` runs lint and the fast tests on every push (`gh run watch`)

## Project rules
- `legacy/` (paper-era scripts `Brian2019MNRAS/`, the handover attempt `outflow.py`, the paper
  PDF) is git-ignored and is the physics reference. Never edit it; never commit it.
- "Same as the paper" means the paper code `legacy/Brian2019MNRAS/version1/CRline.py` and its
  output, now stored as package data in `src/outflow/data/paper/`, not `legacy/outflow.py` (which has a temperature
  bug and auto-enables cooling).
- Defaults follow the paper code exactly (values, rounded constants, tolerances, numerical
  method). Changing any of them is a decision for the user. One such decision so far: the last
  stage of the ρ₀ search uses Brent's method (`docs/decisions.md`, "Speed").
- A speed-up that changes results is a numerics change (the user's decision). To show a
  speed-up changes nothing, use `tests/unit/test_integrator.py` (bit-identical RK step) and
  compare all 14 paper runs before and after.
- The app script keeps its `if __name__ == "__main__": main()` guard: the figure tabs' worker
  processes re-import it.
- Cooling is off by default and is not exposed in the UI.

## Conventions
- Internal units are cgs; docstrings give the paper equation number for every physics formula.
- Variable names follow the paper's symbols spelled out (`v_inf`, `q_mass`, `sonic_radius`,
  `nfw_scale_radius`); never reuse one name for two quantities.
- Importing the package has no side effects (no file loads at import time beyond package data
  on first use, no plots).
