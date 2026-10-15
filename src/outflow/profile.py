"""Wind profiles: :class:`WindProfile`, its CSV format, and comparison of two profiles.

The CSV format is shared with the published runs in ``outflow/data/paper/``: header lines
``# key = value`` (parameters in cgs, provenance), one line of column names, then one row per
radius with values written as ``%.10e``.
"""

import dataclasses
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from outflow import __version__, physics
from outflow.constants import BOLTZMANN, MEAN_PARTICLE_MASS
from outflow.parameters import WindParameters

COLUMNS = ("r_cm", "v_cm_s", "rho_g_cm3", "P_erg_cm3", "Pc_erg_cm3", "T_K", "vA_cm_s")
ARRAY_FIELDS = ("r", "v", "rho", "P", "Pc", "T", "vA")
PARAMETER_NAMES = tuple(f.name for f in dataclasses.fields(WindParameters))
VALUE_FORMAT = "%.10e"


@dataclass
class WindProfile:
    """A steady wind solution on the radii chosen by the integrator (cgs).

    Attributes:
        r: radius (cm), starting at 0.
        v: velocity (cm/s).
        rho: gas density (g/cm^3).
        P: gas pressure (erg/cm^3).
        Pc: cosmic-ray pressure (erg/cm^3).
        T: temperature (K), T = mu P / (k_B rho) at each radius.
        vA: Alfven speed (cm/s).
        params: the inputs that produced it.
        stalled: True if the wind fell back to subsonic beyond 1.1 r_sb before the requested
            outer radius (e.g. a heavy halo); the profile then ends at the stall radius.
        cooling_used: True if radiative cooling was on for this solution.
        outer_radius: the requested outer radius (cm), if known.
        metadata: other header entries (for example the source of a published run).
    """

    r: np.ndarray
    v: np.ndarray
    rho: np.ndarray
    P: np.ndarray
    Pc: np.ndarray
    T: np.ndarray
    vA: np.ndarray
    params: WindParameters
    stalled: bool = False
    cooling_used: bool = False
    outer_radius: float | None = None
    metadata: dict = field(default_factory=dict)

    # --- derived quantities ------------------------------------------------------------------

    @property
    def n(self) -> np.ndarray:
        """Number density n = rho / mu (cm^-3), with mu = 1.4 m_H."""
        return self.rho / MEAN_PARTICLE_MASS

    @property
    def ram_pressure(self) -> np.ndarray:
        """Ram pressure rho v^2 (erg/cm^3)."""
        return self.rho * self.v**2

    @property
    def sound_speed(self) -> np.ndarray:
        """Effective sound speed c* of gas plus cosmic rays (cm/s), paper eq. 35.

        At r = 0, where v = v_A = 0, the ratio v_A / v takes its L'Hopital limit.
        """
        params = self.params
        velocity = self.v.astype(float).copy()
        alfven = self.vA.astype(float).copy()
        at_origin = self.r == 0
        if np.any(at_origin):
            dv, dva = physics.central_velocity_gradients(
                float(self.rho[at_origin][0]),
                params.mass_injection_rate,
                params.magnetic_field,
                params.starburst_radius,
            )
            velocity[at_origin], alfven[at_origin] = dv, dva
        squared = physics.effective_sound_speed_squared(self.P, self.Pc, self.rho, velocity, alfven)
        return np.sqrt(squared)

    @property
    def stall_radius(self) -> float | None:
        """Radius where the wind stalled (cm), or None if it reached the outer radius."""
        return float(self.r[-1]) if self.stalled else None

    def __len__(self) -> int:
        return len(self.r)

    # --- comparison --------------------------------------------------------------------------

    def interpolate(self, radii) -> dict:
        """Columns linearly interpolated in r onto ``radii``: {"v": ..., "rho": ..., ...}.

        Radii may repeat in a profile (where the step size fell below round-off, the state
        still changes slightly); the last point at each radius is used.
        """
        r = np.asarray(self.r, dtype=float)
        keep = np.append(np.diff(r) > 0, True)  # last point of each run of equal radii
        radii = np.asarray(radii, dtype=float)
        return {
            name: np.interp(radii, r[keep], np.asarray(getattr(self, name), dtype=float)[keep])
            for name in ARRAY_FIELDS[1:]
        }

    def max_relative_differences(self, reference: "WindProfile", inner_fraction: float = 0.2):
        """Largest |self / reference - 1| per column, as in docs/decisions.md.

        The reference is interpolated onto this profile's radii between ``inner_fraction``
        r_sb and the smaller of the two outer radii. Points where the reference value is zero
        (e.g. v_A for B_0 = 0) are skipped.
        """
        r_sb = self.params.starburst_radius
        r_max = min(self.r[-1], reference.r[-1])
        mask = (self.r >= inner_fraction * r_sb) & (self.r <= r_max)
        radii = self.r[mask]
        theirs = reference.interpolate(radii)
        result = {}
        for name in ARRAY_FIELDS[1:]:
            mine = np.asarray(getattr(self, name), dtype=float)[mask]
            nonzero = theirs[name] != 0
            if not np.any(nonzero):
                result[name] = 0.0
                continue
            result[name] = float(np.max(np.abs(mine[nonzero] / theirs[name][nonzero] - 1.0)))
        return result

    # --- CSV ---------------------------------------------------------------------------------

    def header_items(self) -> list[tuple[str, object]]:
        """Header entries (key, value) written at the top of the CSV file."""
        items: list[tuple[str, object]] = []
        for key in ("description", "source", "figure"):
            if key in self.metadata:
                items.append((key, self.metadata[key]))
        if "description" not in self.metadata:
            items.append(("description", "Wind profile from the outflow package (Yu et al. 2020)"))
        for name in PARAMETER_NAMES:
            items.append((name, getattr(self.params, name)))
        items.append(("radiation_fraction", self.params.radiation_fraction))
        if self.outer_radius is not None:
            items.append(("outer_radius", self.outer_radius))
        items.append(("stalled", self.stalled))
        items.append(("cooling_used", self.cooling_used))
        for key, value in self.metadata.items():
            if key not in ("description", "source", "figure", "package_version"):
                items.append((key, value))
        items.append(("units", "cgs (cm, g, s, erg, gauss, K)"))
        items.append(("package_version", self.metadata.get("package_version", __version__)))
        return items

    def to_csv(self, path) -> Path:
        """Write the profile as CSV (UTF-8, LF line endings). Returns the path."""
        path = Path(path)
        lines = []
        for key, value in self.header_items():
            # Numbers as repr() of a float so they read back exactly; flags as True/False.
            if isinstance(value, bool | str):
                text = str(value)
            else:
                text = repr(float(value))
            lines.append(f"# {key} = {text}")
        lines.append(",".join(COLUMNS))
        table = np.column_stack([getattr(self, name) for name in ARRAY_FIELDS])
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            np.savetxt(
                handle, table, fmt=VALUE_FORMAT, delimiter=",", header="\n".join(lines), comments=""
            )
        return path

    @classmethod
    def from_csv(cls, path) -> "WindProfile":
        """Read a profile written by :meth:`to_csv` (or a published run from the package data)."""
        path = Path(path)
        header: dict[str, str] = {}
        rows = []
        names = None
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    key, sep, value = line[1:].partition("=")
                    if sep:
                        header[key.strip()] = value.strip()
                elif names is None:
                    names = tuple(line.split(","))
                elif line:
                    rows.append([float(x) for x in line.split(",")])
        if names != COLUMNS:
            raise ValueError(f"{path} does not have the columns {','.join(COLUMNS)}")
        table = np.array(rows, dtype=float)
        params = WindParameters(**{name: float(header[name]) for name in PARAMETER_NAMES})
        known = set(PARAMETER_NAMES) | {
            "radiation_fraction",
            "outer_radius",
            "stalled",
            "cooling_used",
            "units",
        }
        metadata = {key: value for key, value in header.items() if key not in known}
        outer = float(header["outer_radius"]) if "outer_radius" in header else None
        return cls(
            *(table[:, i] for i in range(len(COLUMNS))),
            params=params,
            stalled=header.get("stalled", "False") == "True",
            cooling_used=header.get("cooling_used", "False") == "True",
            outer_radius=outer,
            metadata=metadata,
        )

    def temperature_consistent(self, rtol: float = 1e-9) -> bool:
        """True if T = mu P / (k_B rho) at every radius (a basic consistency check)."""
        expected = MEAN_PARTICLE_MASS * self.P / (BOLTZMANN * self.rho)
        return bool(np.allclose(self.T, expected, rtol=rtol, atol=0))
