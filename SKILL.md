---
name: salt-neutronics-tbr
description: >-
  Emulate a Shift Monte Carlo neutronics simulation of a FLiBe (LiF–BeF2)
  molten-salt tritium-breeding blanket WITHOUT running neutron transport on HPC.
  Given a salt composition (mole % BeF2) and a Li-6 enrichment, it interpolates a
  precomputed parameter study to return the tritium breeding ratio (TBR), and can
  also return neutron flux or nuclide number densities, sweep over compositions,
  and render plots. Use this skill whenever the user asks about tritium breeding
  ratio or TBR, FLiBe / LiF-BeF2 salt composition, mole percent of BeF2, Li-6
  (lithium-6) enrichment, the beryllium multiplier, neutron flux or
  isotopic/number densities in a breeding blanket, or wants to mock, emulate,
  predict, or interpolate what a Shift (or any Monte Carlo) neutronics run on HPC
  (e.g. OLCF Frontier) would produce for a given blanket composition — even if
  they do not name the tool. Also trigger it for sweeping TBR across composition
  or plotting TBR / flux / density surfaces.
---

# Salt Neutronics TBR (Shift FLiBe emulator)

## What this skill does

A full Shift Monte Carlo neutronics+depletion run for a fusion breeding blanket is
expensive — it runs on HPC and takes hours. This skill instead **interpolates a
precomputed parameter study** so you can answer "what tritium breeding ratio would
Shift have produced for *this* FLiBe composition?" in milliseconds, on a laptop or
a single HPC core.

The precomputed study (bundled at `data/neutronics_isotopics.h5`) scanned two
design parameters on a grid:

- **Li-6 enrichment** — atom fraction of ⁶Li in the lithium. Grid: `0.07, 0.20,
  0.50, 0.80, 1.00`.
- **Beryllium multiplier** — scale factor on the nominal beryllium atom density,
  a proxy for BeF₂ content. Grid: `0.9, 1.0, 1.1, 1.2, 1.4`.

The headline output is the **tritium breeding ratio (TBR)** — tritium atoms bred
per source neutron. A blanket needs TBR ≳ 1.0 (with margin for losses) to be
tritium self-sufficient.

## The key input mapping: BeF₂ mol % → beryllium multiplier

Users think in **salt composition**; the data is indexed by **beryllium
multiplier**. The skill converts between them:

```
beryllium_multiplier = (mol% BeF2) / 33.33
```

Nominal eutectic FLiBe is 2 LiF : 1 BeF₂ = **33.33 mol % BeF₂**, which maps to
multiplier **1.0**. The scanned multiplier grid (0.9–1.4) therefore covers roughly
**30–46.7 mol % BeF₂**. This assumes the beryllium atom density is proportional to
BeF₂ content; see `references/physics.md` for the full assumption and how to
override the nominal value (`--nominal-bef2`). Valid query ranges:

- BeF₂: **30–46.67 mol %**
- Li-6 enrichment: **0.07–1.0** (default 0.075 ≈ natural when unspecified)

Requests outside these ranges are rejected for single-point queries (you can pass
`--allow-extrapolation` to override) and returned as `NaN` within sweeps.

## Quickstart — the headline TBR query

Run from the skill root directory. No installation is required; the package is
importable in place.

```bash
# Human-readable
python -m salt_neutronics.cli tbr --bef2 33.3 --li6 0.075

# Machine-readable JSON (preferred when you will parse the result)
python -m salt_neutronics.cli tbr --bef2 40 --li6 0.5 --json
```

The JSON report contains everything needed to interpret and cite the result:

```json
{
  "quantity": "tritium_breeding_ratio",
  "input_composition": {"bef2_mol_percent": 40.0, "lif_mol_percent": 60.0, "li6_enrichment": 0.5},
  "derived": {"beryllium_multiplier": 1.2, "nominal_bef2_mol_percent": 33.33},
  "result": {"tbr": 1.3794, "interpretation": "TBR = 1.3794 >= 1.0: ... self-sufficient ..."},
  "provenance": {
    "method": "linear interpolation of a precomputed Shift Monte Carlo table",
    "is_interpolated": true, "is_extrapolated": false,
    "data_file": ".../neutronics_isotopics.h5",
    "grid": {"li6_enrichment": [...], "beryllium_multiplier": [...]},
    "source_strength_n_per_s": 7.53e17, "irradiation_time_s": 86400.0,
    "incident_neutrons": 6.5e22
  }
}
```

You can specify the beryllium multiplier directly instead of BeF₂ with
`--be-multiplier 1.0` (mutually exclusive with `--bef2`).

## Running on HPC

Use the portable launcher, which loads site modules if present, creates/activates
the conda env from `environment.yml` on first use, and forwards all arguments to
the CLI:

```bash
bash scripts/run_on_hpc.sh tbr --bef2 33.3 --li6 0.075 --json
```

This is a **CPU-only, single-node, sub-second** post-processing task — it does not
need GPUs or multiple nodes. To submit as a Slurm batch job, prepend `#SBATCH`
directives to the launcher (an example header is in the script comments). The
launcher works unchanged on a laptop (falls back to system Python).

## Other queries

```bash
# Sweep TBR over composition and write a CSV (+ optional surface plot)
python -m salt_neutronics.cli sweep --bef2-range 30 46 60 --li6 0.075 --plot
# 2-D sweep over both axes
python -m salt_neutronics.cli sweep --bef2-range 30 46 40 --li6-range 0.07 1.0 40 --output sweep.csv

# Neutron flux at a spatial position (cm) for one composition
python -m salt_neutronics.cli flux --position 15 --bef2 33.3 --li6 0.075 --json

# Number density of a nuclide by ZAID (1003 = tritium, 8016 = O-16, ...)
python -m salt_neutronics.cli density --zaid 1003 --position 15 --bef2 33.3 --li6 0.5 --json
```

`--bef2-range LOW HIGH [N]` builds `N` linearly spaced points (default 50). The
sweep CSV columns are `bef2_mol_percent, li6_enrichment, beryllium_multiplier, tbr`.

## Plotting

`sweep --plot` writes a TBR surface. For finer control, import the plotting
helpers (all save PNGs to an output dir and return the path; a headless Agg
backend is selected automatically so they work on compute nodes):

```python
from salt_neutronics import ShiftFlibeProcessor
from salt_neutronics.plotting import plot_tbr_surface, plot_flux_profiles, plot_nuclide_density_surface

proc = ShiftFlibeProcessor()
plot_tbr_surface(proc, mark=(40.0, 0.5), output_dir="figs")   # star marks a queried point
plot_flux_profiles(proc, output_dir="figs")
plot_nuclide_density_surface(proc, zaid=1003, position=15.0, output_dir="figs")
```

## Programmatic use

For multi-step analysis, use the library directly rather than shelling out:

```python
from salt_neutronics import ShiftFlibeProcessor, bef2_molpct_to_be_multiplier
proc = ShiftFlibeProcessor()                       # loads the bundled data file
mult = float(bef2_molpct_to_be_multiplier(40.0))   # -> 1.2
tbr = float(proc.tritium_ratio_interp(0.5, mult))  # -> 1.3794
```

`tritium_ratio_interp` accepts NumPy arrays for vectorized sweeps and returns
`NaN` for out-of-grid points (no exception). Use `proc.in_bounds(li6, mult)` to
test, or `proc.check_bounds(...)` to raise.

## How to interpret results for the user

Always tell the user **whether the value was interpolated or extrapolated**
(`provenance.is_extrapolated`) — an extrapolated number is not backed by any
simulation. Then read the physics:

- **TBR ≥ 1.0** is the self-sufficiency threshold; **≥ ~1.1** gives a comfortable
  margin once first-wall losses, structural absorption, and neutron streaming are
  accounted for.
- In this dataset TBR **rises with both** Li-6 enrichment and BeF₂ content, and
  Li-6 enrichment is the stronger lever (especially below ~50 %). State the trend,
  not just the single number, when it helps the user's decision.
- Remember this is an **emulation by interpolation**, not a fresh transport solve.
  Don't claim Monte Carlo statistical uncertainties; the limiting error is
  interpolation between grid points.

## Reference material

- `references/physics.md` — FLiBe chemistry, TBR definition, the composition→
  multiplier assumption and how to change it, modeling caveats.
- `references/data_schema.md` — exact layout of `neutronics_isotopics.h5`.
- `README.md` — human-oriented overview, install, and dev/test instructions.

## Tests

`python -m pytest tests/` validates physics consistency, interpolation exactness
at grid nodes, the composition mapping, and bounds handling. Run it after any
change to the processor or composition logic.
