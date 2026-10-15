"""Semi-analytic galactic outflow model of Yu et al. 2020 (MNRAS 492, 3179), paper section 2.4.

Importing this package has no side effects: no files are read and nothing is plotted.

Main entry points::

    from outflow import WindParameters, solve_wind
    profile = solve_wind(WindParameters.paper_default())
"""

__version__ = "0.1.0"

from outflow.errors import NoSolutionError  # noqa: E402
from outflow.paper import PaperRun, load_paper_profile, paper_runs  # noqa: E402
from outflow.parameters import WindParameters  # noqa: E402
from outflow.profile import WindProfile  # noqa: E402
from outflow.solver import solve_wind  # noqa: E402

__all__ = [
    "NoSolutionError",
    "PaperRun",
    "WindParameters",
    "WindProfile",
    "__version__",
    "load_paper_profile",
    "paper_runs",
    "solve_wind",
]
