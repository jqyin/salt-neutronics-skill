# salt-neutronics-tbr

Emulate a **Shift Monte Carlo neutronics simulation** of a FLiBe (LiF–BeF₂)
molten-salt tritium-breeding blanket — without running neutron transport on HPC.
Given a salt composition (mole % BeF₂) and a Li-6 enrichment, the tool interpolates
a precomputed parameter study and returns the **tritium breeding ratio (TBR)**
(plus neutron flux and nuclide number densities, sweeps, and plots).

This repository is also packaged as a **Claude skill** (`SKILL.md`) so an AI agent
can drive it directly.

## Why

A full Shift run is an HPC job that takes hours. For design exploration you often
just want "what TBR would Shift give for *this* composition?" The underlying physics
was already scanned on a grid; this tool blends those results in milliseconds, on a
laptop or a single HPC core.

## Install

No installation is strictly required — the package is importable from the repo root
and the CLI runs as a module. For a clean environment:

```bash
conda env create -f environment.yml   # creates 'salt-neutronics'
conda activate salt-neutronics
# or, editable install into any environment:
pip install -e .
```

Dependencies are light: `numpy`, `scipy`, `h5py`, `matplotlib`. No GPU, PyTorch,
OpenMM, or MACE (the original environment carried those for an unrelated MD workflow;
they have been removed).

## Quickstart

```bash
# Headline TBR query (human-readable)
python -m salt_neutronics.cli tbr --bef2 33.3 --li6 0.075

# JSON output for programmatic use
python -m salt_neutronics.cli tbr --bef2 40 --li6 0.5 --json

# Sweep over composition and plot the TBR surface
python -m salt_neutronics.cli sweep --bef2-range 30 46 60 --li6 0.075 --plot

# Run on HPC (loads site modules if present; CPU-only, single node)
bash scripts/run_on_hpc.sh tbr --bef2 33.3 --li6 0.075 --json
```

Inputs and valid ranges:

- **BeF₂ content**: 30–46.67 mol % (maps to beryllium multiplier 0.9–1.4; nominal
  eutectic FLiBe = 33.33 mol % → multiplier 1.0).
- **Li-6 enrichment**: 0.07–1.0 atom fraction (default 0.075 ≈ natural).

Outside these ranges, single-point queries are rejected and sweep points return
`NaN`; pass `--allow-extrapolation` to extrapolate (cautiously).

## Library API

```python
from salt_neutronics import ShiftFlibeProcessor, bef2_molpct_to_be_multiplier

proc = ShiftFlibeProcessor()                       # loads data/neutronics_isotopics.h5
mult = float(bef2_molpct_to_be_multiplier(40.0))   # 1.2
tbr  = float(proc.tritium_ratio_interp(0.5, mult)) # 1.3794
```

`tritium_ratio_interp`, `flux_interp`, and `zaid_cell_interp` accept scalars or
NumPy arrays and return `NaN` outside the grid. See `salt_neutronics/plotting.py` for
figure helpers.

## Layout

```
salt-neutronics-skill/
├── SKILL.md                  # Claude skill manifest (agent-facing)
├── README.md                 # this file
├── pyproject.toml            # packaging + console entry point
├── environment.yml           # slim conda environment
├── salt_neutronics/          # the package
│   ├── processor.py          # load + interpolate the precomputed study
│   ├── composition.py        # BeF2 mol% <-> beryllium multiplier
│   ├── plotting.py           # TBR / flux / density figures (headless-safe)
│   ├── cli.py                # command-line interface
│   └── __main__.py           # `python -m salt_neutronics`
├── scripts/run_on_hpc.sh     # portable HPC launcher
├── data/neutronics_isotopics.h5
├── references/               # physics + data-schema docs
└── tests/test_neutronics.py
```

## What changed from the original

The original `neutronics_agent/` scripts are preserved for provenance. The
modernization:

- Packaged the loose scripts into an importable `salt_neutronics` package with type
  hints, NumPy-style docstrings, and `pathlib` paths (no hardcoded absolute paths).
- **Fixed an off-by-one bug** in cell selection (`searchsorted` without `-1` selected
  the wrong mesh cell).
- Vectorized the TBR computation (single `tensordot` instead of a Python double loop)
  and cache the TBR interpolator.
- Added the **BeF₂ mol % → beryllium multiplier** mapping so composition is a
  first-class input.
- Added a structured CLI (`tbr`/`sweep`/`flux`/`density`) emitting JSON "mock
  simulation reports" with provenance and in/out-of-grid flags.
- Added bounds handling, headless plotting, a slim environment, a portable HPC
  runner, and a `pytest` suite.

## Tests

```bash
python -m pytest tests/
```

## Caveats

This is **interpolation, not simulation**: results emulate Shift and carry
interpolation error, not Monte Carlo statistics. The composition→multiplier mapping
is a first-order (constant total atom density) assumption. See
[`references/physics.md`](references/physics.md) for the full discussion.
