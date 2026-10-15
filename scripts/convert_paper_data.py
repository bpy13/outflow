"""Convert the paper-era data files into the package's CSV data.

Needs the git-ignored ``legacy/`` folder (paper-era code and data, kept locally only), so it is
not run in CI. Run from anywhere:

    python scripts/convert_paper_data.py

It writes:

- ``src/outflow/data/cloudy_cooling.csv`` from ``legacy/Cloudy.npy``
  (the Cloudy cooling function used by the experimental cooling option);
- ``src/outflow/data/paper/<name>.csv``, one file per published run, from
  ``legacy/Brian2019MNRAS/version1/CRline.npy`` (the data behind the paper's Figs. 6-9).
"""

from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
LEGACY = REPO / "legacy"
DATA = REPO / "src" / "outflow" / "data"

# Paper-code constants (CRline.py), in cgs.
SOLAR_MASS = 2e33  # g
PARSEC = 3.086e18  # cm
KM_S = 1e5  # cm/s
YEAR = 86400 * 365.25  # s

COOLING_HEADER = (
    "# Cloudy cooling function, solar metallicity (Yu et al. 2020, §2.1.2). Source: func_map.dta"
)

# Shared settings of all 14 runs, as in version1/CRline.py ("setting up parameter space").
STARBURST_RADIUS_PC = 200.0
# The paper code's dEdt = 1e12/1.543 in natural units (M_sun, pc, km/s; time unit pc/(km/s)),
# which is 4.2e41 erg/s (not the 1e42 of the paper's Table 1; see docs/decisions.md).
ENERGY_RATE = 1e12 / 1.543 * SOLAR_MASS * KM_S**2 / (PARSEC / KM_S)
MASS_RATE = 2.6 * SOLAR_MASS / YEAR
NFW_SCALE_RADIUS = 14700.0 * PARSEC
VIRIAL_RADIUS = 164e3 * PARSEC
NFW_DENSITY_M82 = 8.81e-3 * SOLAR_MASS / PARSEC**3  # characteristic density rho_M, g/cm^3

# npy run index -> (file name, figure, thermal, CR, opacity cm^2/g, B0 uG, halo factor,
# outer radius in r_sb). Rebuilt from version1/CRline.py; see "Reference run map" in
# docs/decisions.md.
RUNS = [
    ("fig6_thermal", 6, 0.8, 0.1, 1e4, 50.0, 1.0, 100),
    ("fig6_cosmic_ray", 6, 0.1, 0.8, 1e4, 50.0, 1.0, 100),
    ("fig6_radiation", 6, 0.1, 0.1, 1e4, 50.0, 1.0, 100),
    ("fig7_B0_0uG", 7, 0.1, 0.8, 1e4, 0.0, 1.0, 5),
    ("fig7_B0_5uG", 7, 0.1, 0.8, 1e4, 5.0, 1.0, 5),
    ("fig7_B0_50uG", 7, 0.1, 0.8, 1e4, 50.0, 1.0, 5),
    ("fig7_B0_100uG", 7, 0.1, 0.8, 1e4, 100.0, 1.0, 5),
    ("fig8_kappa_0", 8, 0.1, 0.1, 0.0, 50.0, 1.0, 5),
    ("fig8_kappa_1e3", 8, 0.1, 0.1, 1e3, 50.0, 1.0, 5),
    ("fig8_kappa_1e4", 8, 0.1, 0.1, 1e4, 50.0, 1.0, 5),
    ("fig8_kappa_1e5", 8, 0.1, 0.1, 1e5, 50.0, 1.0, 5),
    ("fig9_halo_0x", 9, 0.8, 0.1, 1e4, 50.0, 0.0, 50),
    ("fig9_halo_1x", 9, 0.8, 0.1, 1e4, 50.0, 1.0, 50),
    ("fig9_halo_10x", 9, 0.8, 0.1, 1e4, 50.0, 10.0, 50),
]

COLUMNS = "r_cm,v_cm_s,rho_g_cm3,P_erg_cm3,Pc_erg_cm3,T_K,vA_cm_s"


def nfw_mass(characteristic_density, scale_radius, virial_radius):
    """Total NFW halo mass within the virial radius (g)."""
    concentration = virial_radius / scale_radius
    return (
        4
        * np.pi
        * characteristic_density
        * scale_radius**3
        * (np.log(1 + concentration) - concentration / (1 + concentration))
    )


def write_csv(path, table, fmt, header):
    """Write header lines then comma-separated rows, UTF-8 with LF line endings on every OS."""
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        np.savetxt(handle, table, fmt=fmt, delimiter=",", header=header, comments="")


def convert_cooling():
    """Write the Cloudy cooling table as CSV, with values that round-trip exactly."""
    temperature, cooling = np.load(LEGACY / "Cloudy.npy")
    out = DATA / "cloudy_cooling.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    write_csv(
        out,
        np.column_stack([temperature, cooling]),
        "%.17g",
        COOLING_HEADER + "\nT_K,Lambda_erg_cm3_s",
    )
    print(f"wrote {out.relative_to(REPO)} ({temperature.size} rows)")


def run_header(index, figure, thermal, cr, opacity, b0_microgauss, halo_factor, outer):
    """Header lines (``# key = value``) for one published run, parameters in cgs."""
    starburst_radius = STARBURST_RADIUS_PC * PARSEC
    halo_mass = halo_factor * nfw_mass(NFW_DENSITY_M82, NFW_SCALE_RADIUS, VIRIAL_RADIUS)
    items = [
        ("description", f"Published run for Fig. {figure} of Yu et al. 2020 (MNRAS 492, 3179)"),
        ("source", f"version1/CRline.npy run {index} (paper code), converted to cgs"),
        ("figure", str(figure)),
        ("starburst_radius", starburst_radius),
        ("energy_rate", ENERGY_RATE),
        ("mass_rate", MASS_RATE),
        ("thermal_fraction", thermal),
        ("cr_fraction", cr),
        ("radiation_fraction", round(1 - thermal - cr, 12)),
        ("opacity", opacity),
        ("magnetic_field", b0_microgauss / 1e6),  # gauss
        ("halo_mass", halo_mass),
        ("nfw_scale_radius", NFW_SCALE_RADIUS),
        ("virial_radius", VIRIAL_RADIUS),
        ("outer_radius", outer * starburst_radius),
        ("units", "cgs (cm, g, s, erg, gauss, K)"),
        ("package_version", "0.1.0"),
    ]
    # Numbers are written with repr() of a plain float, so they read back exactly.
    return "\n".join(
        f"# {key} = {value}" if isinstance(value, str) else f"# {key} = {float(value)!r}"
        for key, value in items
    )


def convert_paper_runs():
    """Write the 14 published runs of version1/CRline.npy as CSV files in cgs."""
    data = np.load(LEGACY / "Brian2019MNRAS" / "version1" / "CRline.npy", allow_pickle=True)
    r_over_rsb, v, pressure, cr_pressure, density, temperature, alfven = data
    out_dir = DATA / "paper"
    out_dir.mkdir(parents=True, exist_ok=True)
    # Conversion factors of the paper code's plotting lines (natural units -> cgs).
    density_factor = 2 / (3.086**3 * 1e21)
    pressure_factor = 2 / (3.086**3 * 1e11)
    for index, (name, figure, *settings) in enumerate(RUNS):
        table = np.column_stack(
            [
                np.asarray(r_over_rsb[index]) * STARBURST_RADIUS_PC * PARSEC,
                np.asarray(v[index]) * KM_S,
                np.asarray(density[index]) * density_factor,
                np.asarray(pressure[index]) * pressure_factor,
                np.asarray(cr_pressure[index]) * pressure_factor,
                np.asarray(temperature[index]),
                np.asarray(alfven[index]) * KM_S,
            ]
        )
        out = out_dir / f"{name}.csv"
        header = run_header(index, figure, *settings) + "\n" + COLUMNS
        write_csv(out, table, "%.10e", header)
        print(f"wrote {out.relative_to(REPO)} ({len(table)} rows)")


if __name__ == "__main__":
    convert_cooling()
    convert_paper_runs()
