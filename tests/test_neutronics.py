"""Tests for the salt-neutronics TBR skill.

Run with ``pytest`` from the skill root. These tests validate physics
consistency (TBR computed two independent ways), interpolation exactness at grid
nodes, the composition mapping, the cell-index bug fix, and bounds handling.
"""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from salt_neutronics import ShiftFlibeProcessor
from salt_neutronics.composition import (
    NOMINAL_BEF2_MOLPCT,
    bef2_molpct_to_be_multiplier,
    be_multiplier_to_bef2_molpct,
)
from salt_neutronics.processor import (
    DEFAULT_DATA_FILE,
    MAGNET_FLUX_PREFERRED_N_PER_CM2_S,
    MAGNET_FLUX_REJECT_N_PER_CM2_S,
    SHIELDING_THICKNESS_CM,
    TRITIUM_ZAID,
    shielding_verdict,
)


@pytest.fixture(scope="module")
def proc() -> ShiftFlibeProcessor:
    return ShiftFlibeProcessor()


def test_default_data_file_exists() -> None:
    assert DEFAULT_DATA_FILE.is_file()


def test_shapes_consistent(proc: ShiftFlibeProcessor) -> None:
    assert proc.flux.shape == (proc.num_cells, proc.num_enrichments, proc.num_multipliers)
    assert proc.number_density.shape == (
        proc.num_nuclides, proc.num_cells, proc.num_enrichments, proc.num_multipliers,
    )
    assert proc.tritium_ratio.shape == (proc.num_enrichments, proc.num_multipliers)


def test_tbr_matches_independent_calculation(proc: ShiftFlibeProcessor) -> None:
    """Recompute TBR straight from the HDF5 and compare to the processor."""
    with h5py.File(DEFAULT_DATA_FILE, "r") as f:
        nd = f["number_density"][:]
        vols = f["cell_volumes"][:]
        nucs = f["nuclide_list"][:].tolist()
        inc = float(f["source_strength"][()]) * float(f["irradiation_time"][()])
    ti = nucs.index(TRITIUM_ZAID)
    expected = np.tensordot(nd[ti], vols, axes=([0], [0])) / inc
    np.testing.assert_allclose(proc.tritium_ratio, expected, rtol=1e-12)


def test_tbr_physically_reasonable(proc: ShiftFlibeProcessor) -> None:
    assert np.all(proc.tritium_ratio > 0.5)
    assert np.all(proc.tritium_ratio < 3.0)


def test_interpolation_exact_at_grid_nodes(proc: ShiftFlibeProcessor) -> None:
    """Linear interpolation must reproduce the table exactly at every node."""
    for i, li6 in enumerate(proc.lithium6_enrichments):
        for j, be in enumerate(proc.beryllium_multipliers):
            got = float(proc.tritium_ratio_interp(li6, be))
            assert got == pytest.approx(proc.tritium_ratio[i, j], rel=1e-12)


def test_interpolation_monotone_midpoint(proc: ShiftFlibeProcessor) -> None:
    """A midpoint value should lie between its bracketing node values."""
    li6 = proc.lithium6_enrichments
    be = proc.beryllium_multipliers
    mid_li6 = 0.5 * (li6[0] + li6[1])
    mid_be = 0.5 * (be[0] + be[1])
    corners = [
        proc.tritium_ratio[0, 0], proc.tritium_ratio[0, 1],
        proc.tritium_ratio[1, 0], proc.tritium_ratio[1, 1],
    ]
    val = float(proc.tritium_ratio_interp(mid_li6, mid_be))
    assert min(corners) <= val <= max(corners)


def test_tbr_array_input_shape(proc: ShiftFlibeProcessor) -> None:
    li6 = np.array([0.1, 0.3, 0.6])
    be = np.array([1.0, 1.1, 1.2])
    out = proc.tritium_ratio_interp(li6, be)
    assert out.shape == (3,)
    assert np.all(np.isfinite(out))


def test_cell_index_fix(proc: ShiftFlibeProcessor) -> None:
    """Position 15 cm lies in cell [10, 20) -> index 1 (was 2 in the old code)."""
    # cell_edges are 0,10,...,100 with 10 cells.
    assert proc.cell_index(15.0) == 1
    assert proc.cell_index(0.0) == 0
    assert proc.cell_index(5.0) == 0
    assert proc.cell_index(95.0) == proc.num_cells - 1
    # Out-of-domain positions clamp into a valid cell.
    assert proc.cell_index(1000.0) == proc.num_cells - 1
    assert proc.cell_index(-5.0) == 0


def test_out_of_bounds_returns_nan(proc: ShiftFlibeProcessor) -> None:
    # Interpolation never raises on OOB; it returns NaN so sweeps degrade gracefully.
    assert np.isnan(float(proc.tritium_ratio_interp(0.075, 5.0)))
    assert proc.in_bounds(0.075, 1.0) is True
    assert proc.in_bounds(0.075, 5.0) is False


def test_check_bounds_raises(proc: ShiftFlibeProcessor) -> None:
    # The explicit guard (used by single-point CLI queries) does raise.
    with pytest.raises(ValueError):
        proc.check_bounds(0.075, 5.0)


def test_extrapolation_allowed() -> None:
    p = ShiftFlibeProcessor(allow_extrapolation=True)
    val = float(p.tritium_ratio_interp(0.075, 1.5))  # just past grid max (1.4)
    assert np.isfinite(val)


def test_unknown_zaid_raises(proc: ShiftFlibeProcessor) -> None:
    with pytest.raises(KeyError):
        proc.nuclide_index(99999)


def test_flux_and_density_finite_at_node(proc: ShiftFlibeProcessor) -> None:
    li6 = proc.lithium6_enrichments[1]
    be = proc.beryllium_multipliers[1]
    assert np.isfinite(float(proc.flux_interp(15.0, li6, be)))
    assert np.isfinite(float(proc.zaid_cell_interp(8016, 15.0, li6, be)))


def test_composition_roundtrip() -> None:
    for bef2 in [25.0, 33.3333, 40.0, 46.6]:
        mult = bef2_molpct_to_be_multiplier(bef2)
        back = be_multiplier_to_bef2_molpct(mult)
        assert float(back) == pytest.approx(bef2, rel=1e-12)


def test_nominal_flibe_maps_to_unity() -> None:
    assert float(bef2_molpct_to_be_multiplier(NOMINAL_BEF2_MOLPCT)) == pytest.approx(1.0)


def test_boundary_composition_not_nan(proc: ShiftFlibeProcessor) -> None:
    """30 mol% BeF2 -> multiplier 0.8999999999999999 must still resolve (FP snap)."""
    mult = float(bef2_molpct_to_be_multiplier(30.0))
    assert mult < 0.9  # genuinely just below the grid edge before snapping
    assert proc.in_bounds(0.075, mult) is True
    assert np.isfinite(float(proc.tritium_ratio_interp(0.075, mult)))


def test_composition_rejects_invalid() -> None:
    with pytest.raises(ValueError):
        bef2_molpct_to_be_multiplier(150.0)


# --- radiation shielding (magnet flux) ------------------------------------- #

def test_shielding_flux_is_flux_at_reference_thickness(proc: ShiftFlibeProcessor) -> None:
    """The magnet shielding flux is the neutron flux at the blanket back face (default 1 m)."""
    assert SHIELDING_THICKNESS_CM == 100.0
    got = float(proc.shielding_flux(0.5, 1.0))
    expected = float(proc.flux_interp(SHIELDING_THICKNESS_CM, 0.5, 1.0))
    assert got == pytest.approx(expected)
    assert got > 0.0


def test_shielding_flux_respects_thickness(proc: ShiftFlibeProcessor) -> None:
    """A thinner blanket shields less, so the magnet sees a higher (or equal) flux."""
    thin = float(proc.shielding_flux(0.5, 1.0, thickness_cm=20.0))
    thick = float(proc.shielding_flux(0.5, 1.0, thickness_cm=100.0))
    assert thin >= thick


def test_shielding_verdict_thresholds() -> None:
    assert shielding_verdict(MAGNET_FLUX_REJECT_N_PER_CM2_S * 10) == "reject"
    assert shielding_verdict(MAGNET_FLUX_PREFERRED_N_PER_CM2_S / 10) == "preferred"
    # Between the two thresholds is acceptable.
    mid = (MAGNET_FLUX_REJECT_N_PER_CM2_S + MAGNET_FLUX_PREFERRED_N_PER_CM2_S) / 2
    assert shielding_verdict(mid) == "acceptable"
