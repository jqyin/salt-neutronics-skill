# Physics reference

This document explains the physical quantities, the FLiBe composition mapping, and
the modeling assumptions/caveats behind the skill. Read it when you need to justify
or adjust a result, or when a user questions the composition mapping.

## Contents

- [The system: a FLiBe breeding blanket](#the-system-a-flibe-breeding-blanket)
- [Tritium breeding ratio (TBR)](#tritium-breeding-ratio-tbr)
- [The two scanned design parameters](#the-two-scanned-design-parameters)
- [Composition mapping: BeF2 mol % ↔ beryllium multiplier](#composition-mapping-bef2-mol--beryllium-multiplier)
- [How TBR is computed from the table](#how-tbr-is-computed-from-the-table)
- [Interpolation method](#interpolation-method)
- [Caveats and limits](#caveats-and-limits)

## The system: a FLiBe breeding blanket

FLiBe is a molten salt, a binary mixture of lithium fluoride (`LiF`) and beryllium
fluoride (`BeF₂`). In a fusion reactor it surrounds the plasma as a **breeding
blanket**: 14 MeV fusion neutrons enter the salt and breed tritium via

- `⁶Li(n, α)T` — the dominant, exothermic channel, strong at low neutron energy.
- `⁷Li(n, n'α)T` — a threshold reaction that also yields a neutron.

Beryllium acts as a **neutron multiplier** via `⁹Be(n, 2n)`, increasing the neutron
population available to breed tritium. More beryllium (more BeF₂) therefore tends to
raise the breeding ratio, as does enriching the lithium in ⁶Li.

## Tritium breeding ratio (TBR)

The TBR is the number of tritium atoms bred per source (fusion) neutron:

```
TBR = (total tritium atoms produced) / (total incident source neutrons)
```

A reactor must replace the tritium it burns plus losses, so a blanket design
targets **TBR > 1**, typically with margin (≈1.05–1.15) to cover neutron streaming
through penetrations, parasitic absorption in structure, and radioactive decay.

In this dataset the TBR ranges from about **1.19** (low Li-6, low BeF₂) to **1.40**
(high Li-6, high BeF₂), so every scanned composition is at least nominally
self-sufficient.

## The two scanned design parameters

| Parameter | Symbol in data | Grid values | Meaning |
|---|---|---|---|
| Li-6 enrichment | `lithium6_enrichments` | 0.07, 0.20, 0.50, 0.80, 1.00 | Atom fraction of ⁶Li in Li. Natural Li is ≈7.5 % ⁶Li. |
| Beryllium multiplier | `beryllium_multipliers` | 0.9, 1.0, 1.1, 1.2, 1.4 | Scale factor on nominal Be atom density (proxy for BeF₂ content). |

## Composition mapping: BeF₂ mol % ↔ beryllium multiplier

The Shift study scanned a dimensionless beryllium multiplier, but the natural design
knob is the salt composition. The skill uses:

```
beryllium_multiplier = (mol% BeF2) / (nominal mol% BeF2),   nominal = 33.33
```

**Why 33.33 mol %?** Standard eutectic FLiBe is the 2:1 molar mixture
`2 LiF · 1 BeF₂`, i.e. ⅓ of the molecules are BeF₂ → 33.33 mol % BeF₂. This is
defined as the nominal composition (beryllium multiplier = 1.0).

**The assumption.** The mapping is exact if the beryllium **atom density** scales
linearly with BeF₂ mole fraction — i.e. the total atom density of the salt is held
constant as the LiF:BeF₂ ratio changes. This is a good first-order approximation and
matches how the multiplier was intended. It neglects the second-order change in the
salt's molar volume (density) as composition shifts. If you have a measured or
modeled density relation `ρ(x_BeF2)`, the more accurate multiplier is

```
beryllium_multiplier = n_Be(x) / n_Be(nominal)
```

where `n_Be` is the actual beryllium atom number density at composition `x`.

**Overriding the nominal.** Both the CLI (`--nominal-bef2 VALUE`) and the library
(`bef2_molpct_to_be_multiplier(molpct, nominal_bef2_molpct=...)`) let you change the
composition that maps to multiplier 1.0, or you can bypass composition entirely and
pass `--be-multiplier` directly.

**Resulting valid range.** The multiplier grid 0.9–1.4 corresponds to BeF₂ contents
of `0.9·33.33 ≈ 30.0` to `1.4·33.33 ≈ 46.67` mol %, bracketing the nominal eutectic.

## How TBR is computed from the table

For each grid point `(Li-6, Be-multiplier)`:

1. The total incident neutrons are `source_strength × irradiation_time`
   (here `7.53e17 n/s × 86400 s = 6.506e22`).
2. The total bred tritium is the volume integral of the tritium (ZAID 1003) number
   density over all mesh cells: `Σ_c n_T(c) · V(c)`.
3. `TBR = tritium_atoms / incident_neutrons`.

This produces a `(num_enrichments × num_multipliers)` table of TBR values that the
skill then interpolates.

## Interpolation method

The skill uses SciPy's `RegularGridInterpolator` with **linear** interpolation over
the 2-D `(Li-6 enrichment, Be multiplier)` grid. At a grid node the interpolant
returns the tabulated value exactly; between nodes it blends the surrounding values
linearly.

Spatial dependence (for flux and number density) is **piecewise constant per mesh
cell**: the code selects the cell containing the requested position, then interpolates
only over the two design parameters. There is no continuous spatial interpolation —
that matches the resolution of the stored data.

## Caveats and limits

- **This is interpolation, not simulation.** Results emulate what Shift would have
  produced; they carry interpolation error between grid points, not Monte Carlo
  statistical uncertainty. Do not quote σ/relative-error figures as if from a fresh
  transport solve.
- **Stay inside the grid.** Outside 30–46.67 mol % BeF₂ or 0.07–1.0 Li-6 enrichment,
  single-point queries are rejected and sweep points return `NaN`. Extrapolation is
  available (`--allow-extrapolation`) but is physically unreliable, especially for
  the threshold `⁷Li` reaction and `(n,2n)` multiplication.
- **Composition mapping is first-order.** See the assumption above; for high-accuracy
  work, supply a density-corrected multiplier.
- **Fixed irradiation scenario.** `source_strength` and `irradiation_time` are baked
  into the precomputed table; this skill does not re-scale to other source histories.
- **TBR here excludes engineering losses.** The values are intrinsic blanket breeding;
  a real design must keep margin for streaming and structural absorption.
