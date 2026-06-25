# Data schema: `neutronics_isotopics.h5`

Layout of the precomputed Shift parameter study consumed by the skill. All datasets
sit at the HDF5 root. Dimensions for the bundled file are shown in brackets.

## Scalars

| Dataset | Type | Value (bundled) | Meaning |
|---|---|---|---|
| `num_cells` | int | 10 | Number of spatial mesh cells. |
| `num_nuclides` | int | 20 | Number of tracked nuclides. |
| `num_enrichments` | int | 5 | Number of Li-6 enrichment cases. |
| `num_multipliers` | int | 5 | Number of beryllium-multiplier cases. |
| `source_strength` | float | 7.53e17 | Neutron source rate (n/s). |
| `irradiation_time` | float | 86400.0 | Irradiation time (s) = 1 day. |

## 1-D arrays

| Dataset | Shape | Value (bundled) | Meaning |
|---|---|---|---|
| `cell_edges` | (num_cells + 1,) = (11,) | 0,10,…,100 | Mesh-cell edges (cm). Cell `k` spans `[edges[k], edges[k+1])`. |
| `cell_volumes` | (num_cells,) = (10,) | all 85000 | Volume of each mesh cell. |
| `lithium6_enrichments` | (num_enrichments,) = (5,) | 0.07, 0.20, 0.50, 0.80, 1.00 | Li-6 enrichment grid (atom fraction). |
| `beryllium_multipliers` | (num_multipliers,) = (5,) | 0.9, 1.0, 1.1, 1.2, 1.4 | Be-multiplier grid. |
| `nuclide_list` | (num_nuclides,) = (20,) | ZAIDs (see below) | Nuclide identifiers, ZAID = 1000·Z + A. |

Bundled `nuclide_list` (ZAIDs): `1001 1002 1003 2003 2004 3006 3007 4009 4010
5010 5011 6012 6013 7014 7015 8016 8017 8018 9019 10020`. Tritium is `1003`.

## Multi-dimensional arrays

| Dataset | Shape | Index order | Meaning |
|---|---|---|---|
| `flux` | (num_cells, num_enrichments, num_multipliers) = (10, 5, 5) | `[cell, li6, be]` | Neutron flux (n/cm²·s). |
| `number_density` | (num_nuclides, num_cells, num_enrichments, num_multipliers) = (20, 10, 5, 5) | `[nuclide, cell, li6, be]` | Per-nuclide number density. |
| `original_number_density` | (num_nuclides, num_enrichments, num_multipliers) = (20, 5, 5) | `[nuclide, li6, be]` | Pre-irradiation (initial) densities. Not used for TBR; retained for reference. |

## Notes

- **Index convention.** The `li6` axis indexes `lithium6_enrichments`; the `be` axis
  indexes `beryllium_multipliers`. Interpolation is performed over these two axes.
- **Cell selection.** `cell_index(position)` uses
  `searchsorted(cell_edges, position, side="right") - 1`, clamped to a valid cell.
  This places position 15 cm in cell 1 (`[10, 20)`). The original code omitted the
  `-1` and selected the cell above; the modernized processor fixes this.
- **TBR derivation.** See `physics.md`. The skill integrates the `1003` row of
  `number_density` over cell volumes and divides by `source_strength·irradiation_time`.
- **Adding data.** A larger study (more grid points or parameters) can be dropped in
  by matching this schema; the processor validates array shapes against the scalar
  dimensions on load and will fail fast on a mismatch.
