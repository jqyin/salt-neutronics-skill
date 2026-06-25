"""Post-processing of precomputed Shift FLiBe neutronics/depletion results.

This module reads a precomputed parameter study (stored in HDF5) of a molten-salt
(FLiBe: ``LiF`` + ``BeF2``) tritium-breeding blanket and exposes fast interpolation
over the two free design parameters that were scanned:

* **Li-6 enrichment** -- atom fraction of :sup:`6`\\ Li in the lithium.
* **Beryllium multiplier** -- multiplicative factor on the nominal beryllium atom
  density (a proxy for the ``BeF2`` content of the salt; see
  :mod:`salt_neutronics.composition`).

It does **not** run a neutron-transport calculation. It interpolates an existing
table, which is exactly what lets an agent cheaply *mimic* what a full Monte Carlo
(Shift) run on HPC would have produced for an unscanned composition.

The headline quantity is the **tritium breeding ratio (TBR)** -- tritium atoms bred
per source neutron -- which must exceed ~1.0 for a fusion blanket to be
self-sufficient.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import RegularGridInterpolator

__all__ = ["ShiftFlibeProcessor", "GridBounds", "TRITIUM_ZAID", "DEFAULT_DATA_FILE"]

logger = logging.getLogger(__name__)

#: ZAID (1000*Z + A) of tritium / hydrogen-3, the bred nuclide of interest.
TRITIUM_ZAID = 1003

#: Absolute tolerance for snapping a query coordinate onto a grid boundary. This
#: absorbs floating-point noise (e.g. 30 mol% BeF2 -> 0.8999999999999999, which
#: would otherwise fall just outside the multiplier grid minimum of 0.9) without
#: affecting genuinely out-of-grid requests. It is far smaller than any grid
#: spacing in the dataset.
BOUND_ATOL = 1e-9

#: Default precomputed results file, resolved relative to this package's
#: installation so the skill works regardless of the current working directory.
DEFAULT_DATA_FILE = (Path(__file__).resolve().parent.parent / "data" / "neutronics_isotopics.h5")

#: Environment variable that overrides the data file location (e.g. an HPC scratch
#: path), checked when no explicit ``filename`` is passed to the processor.
DATA_FILE_ENV_VAR = "SALT_NEUTRONICS_DATA"


def resolve_data_file(filename: str | Path | None = None) -> Path:
    """Resolve the results-file path: explicit arg > env var > bundled default."""
    if filename is not None:
        return Path(filename)
    env = os.environ.get(DATA_FILE_ENV_VAR)
    return Path(env) if env else DEFAULT_DATA_FILE


@dataclass(frozen=True)
class GridBounds:
    """Inclusive min/max of a scanned parameter, used for bounds checking."""

    name: str
    low: float
    high: float

    def contains(self, value: ArrayLike) -> NDArray[np.bool_]:
        """Return a boolean mask of which ``value`` entries lie within the grid."""
        arr = np.asarray(value, dtype=float)
        return (arr >= self.low) & (arr <= self.high)

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return f"{self.name} in [{self.low:g}, {self.high:g}]"


class ShiftFlibeProcessor:
    """Read and interpolate a precomputed Shift FLiBe neutronics parameter study.

    Parameters
    ----------
    filename:
        Path to the HDF5 results file. Defaults to the data file bundled with
        the skill (:data:`DEFAULT_DATA_FILE`).
    allow_extrapolation:
        If ``False`` (default), requesting a point outside the scanned
        ``(Li-6, Be-multiplier)`` grid raises :class:`ValueError`. If ``True``,
        SciPy linearly extrapolates beyond the grid and a warning is logged --
        extrapolated values are not backed by simulation data and should be
        treated with caution.

    Notes
    -----
    Spatial dependence is treated as *piecewise constant per mesh cell*: flux and
    number density are looked up in the cell that contains the requested
    position, then interpolated only over ``(Li-6, Be-multiplier)``. This matches
    the resolution of the stored data; there is no continuous spatial
    interpolation.
    """

    def __init__(
        self,
        filename: str | Path | None = None,
        *,
        allow_extrapolation: bool = False,
    ) -> None:
        self.filename = resolve_data_file(filename)
        if not self.filename.is_file():
            raise FileNotFoundError(
                f"Neutronics results file not found: {self.filename}. "
                "Pass an explicit path or place 'neutronics_isotopics.h5' in the "
                "skill's data/ directory."
            )
        self.allow_extrapolation = allow_extrapolation

        with h5py.File(self.filename, "r") as h5f:
            # Scalar metadata describing the parameter study.
            self.num_cells = int(h5f["num_cells"][()])
            self.num_nuclides = int(h5f["num_nuclides"][()])
            self.num_enrichments = int(h5f["num_enrichments"][()])
            self.num_multipliers = int(h5f["num_multipliers"][()])
            self.source_strength = float(h5f["source_strength"][()])  # neutrons / s
            self.irradiation_time = float(h5f["irradiation_time"][()])  # s

            # 1-D coordinate / index arrays.
            self.cell_edges = np.asarray(h5f["cell_edges"][:], dtype=float)
            self.beryllium_multipliers = np.asarray(h5f["beryllium_multipliers"][:], dtype=float)
            self.lithium6_enrichments = np.asarray(h5f["lithium6_enrichments"][:], dtype=float)
            self.nuclides = np.asarray(h5f["nuclide_list"][:]).astype(np.int64)
            self.volumes = np.asarray(h5f["cell_volumes"][:], dtype=float)

            # Multi-dimensional result arrays.
            #   flux:           (cell, enrichment, multiplier)
            #   number_density: (nuclide, cell, enrichment, multiplier)
            self.flux = np.asarray(h5f["flux"][:], dtype=float)
            self.number_density = np.asarray(h5f["number_density"][:], dtype=float)

        self._validate_shapes()

        # Map ZAID -> row index once, instead of scanning a list on every call.
        self._zaid_to_index = {int(z): i for i, z in enumerate(self.nuclides)}

        # Convenient bounds objects for validation and error messages.
        self.li6_bounds = GridBounds(
            "Li-6 enrichment",
            float(self.lithium6_enrichments.min()),
            float(self.lithium6_enrichments.max()),
        )
        self.be_bounds = GridBounds(
            "Beryllium multiplier",
            float(self.beryllium_multipliers.min()),
            float(self.beryllium_multipliers.max()),
        )

        self.tritium_ratio = self._compute_tritium_ratio()
        self._tbr_interp = self._make_interpolator(self.tritium_ratio)

    # ------------------------------------------------------------------ #
    # Construction helpers
    # ------------------------------------------------------------------ #
    def _validate_shapes(self) -> None:
        """Fail fast with a clear message if the file is internally inconsistent."""
        expected_flux = (self.num_cells, self.num_enrichments, self.num_multipliers)
        if self.flux.shape != expected_flux:
            raise ValueError(f"flux shape {self.flux.shape} != expected {expected_flux}")

        expected_nd = (self.num_nuclides, self.num_cells, self.num_enrichments, self.num_multipliers)
        if self.number_density.shape != expected_nd:
            raise ValueError(
                f"number_density shape {self.number_density.shape} != expected {expected_nd}"
            )
        if self.volumes.shape != (self.num_cells,):
            raise ValueError(
                f"cell_volumes shape {self.volumes.shape} != expected {(self.num_cells,)}"
            )

    def _compute_tritium_ratio(self) -> NDArray[np.float64]:
        """Tritium atoms bred per source neutron for every grid point.

        Returns
        -------
        ndarray, shape ``(num_enrichments, num_multipliers)``
            ``TBR[i, j]`` for the ``i``-th Li-6 enrichment and ``j``-th Be
            multiplier.

        Notes
        -----
        Total bred tritium in a case is the volume integral of the tritium
        number density over all mesh cells, ``sum_c n_T(c) * V(c)``. Dividing by
        the total incident neutrons (``source_strength * irradiation_time``)
        gives the breeding ratio. The integral is a single vectorized
        ``tensordot`` over the cell axis rather than a Python loop.
        """
        incident_neutrons = self.source_strength * self.irradiation_time
        if incident_neutrons <= 0:
            raise ValueError(
                "Non-positive incident-neutron normalization "
                f"(source_strength={self.source_strength}, "
                f"irradiation_time={self.irradiation_time})."
            )

        tritium_idx = self.nuclide_index(TRITIUM_ZAID)
        # tritium_density: (cell, enrichment, multiplier)
        tritium_density = self.number_density[tritium_idx]
        # Integrate over the cell axis -> (enrichment, multiplier).
        tritium_atoms = np.tensordot(tritium_density, self.volumes, axes=([0], [0]))
        return tritium_atoms / incident_neutrons

    def _make_interpolator(self, values: NDArray[np.float64]) -> RegularGridInterpolator:
        """Build a linear interpolator over the ``(Li-6, Be-multiplier)`` grid."""
        fill_value = None if self.allow_extrapolation else np.nan
        return RegularGridInterpolator(
            (self.lithium6_enrichments, self.beryllium_multipliers),
            values,
            method="linear",
            bounds_error=False,
            fill_value=fill_value,
        )

    # ------------------------------------------------------------------ #
    # Index helpers
    # ------------------------------------------------------------------ #
    def nuclide_index(self, zaid: int) -> int:
        """Row index of ``zaid`` in :attr:`nuclides`.

        Raises
        ------
        KeyError
            If the nuclide is not tracked in this dataset.
        """
        try:
            return self._zaid_to_index[int(zaid)]
        except KeyError:
            raise KeyError(
                f"ZAID {zaid} is not present in this dataset. "
                f"Available ZAIDs: {self.nuclides.tolist()}"
            ) from None

    def cell_index(self, position: float) -> int:
        """Index of the mesh cell that contains ``position``.

        Cells are defined by :attr:`cell_edges`: cell ``k`` spans
        ``[cell_edges[k], cell_edges[k + 1])``. The result is clamped to a valid
        cell so positions on the domain boundary are handled gracefully.

        Notes
        -----
        This corrects an off-by-one error in the original implementation, which
        used ``np.searchsorted`` directly and therefore selected the cell *above*
        the one containing the position.
        """
        idx = int(np.searchsorted(self.cell_edges, position, side="right") - 1)
        return int(np.clip(idx, 0, self.num_cells - 1))

    # ------------------------------------------------------------------ #
    # Bounds checking
    # ------------------------------------------------------------------ #
    def in_bounds(self, li6: ArrayLike, be_mult: ArrayLike) -> bool:
        """``True`` iff every requested point lies inside the scanned grid."""
        li6_s, be_s = self._snap_to_grid(li6, be_mult)
        return bool(
            np.all(self.li6_bounds.contains(li6_s)) and np.all(self.be_bounds.contains(be_s))
        )

    def _snap_to_grid(
        self, li6: ArrayLike, be_mult: ArrayLike
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Snap coordinates within :data:`BOUND_ATOL` of a grid edge onto that edge."""
        return (
            self._snap_axis(li6, self.li6_bounds),
            self._snap_axis(be_mult, self.be_bounds),
        )

    @staticmethod
    def _snap_axis(values: ArrayLike, bounds: GridBounds) -> NDArray[np.float64]:
        arr = np.array(values, dtype=float)  # copy; never mutate caller arrays
        arr = np.where(np.abs(arr - bounds.low) <= BOUND_ATOL, bounds.low, arr)
        arr = np.where(np.abs(arr - bounds.high) <= BOUND_ATOL, bounds.high, arr)
        return arr

    def check_bounds(self, li6: ArrayLike, be_mult: ArrayLike) -> None:
        """Raise :class:`ValueError` if a request is outside the scanned grid.

        This is *not* called automatically by the interpolation methods (which
        return ``NaN`` outside the grid); callers that want a hard failure with a
        descriptive message -- e.g. a single-point query -- invoke it explicitly.
        Honors :attr:`allow_extrapolation`: when extrapolation is permitted, an
        out-of-bounds request only logs a warning instead of raising.
        """
        if self.in_bounds(li6, be_mult):
            return
        msg = (
            "Requested point is outside the simulated grid "
            f"({self.li6_bounds}; {self.be_bounds}). "
            "Interpolation is only reliable inside the grid."
        )
        if self.allow_extrapolation:
            logger.warning("%s Extrapolating anyway because allow_extrapolation=True.", msg)
        else:
            raise ValueError(msg + " Pass allow_extrapolation=True to override.")

    # ------------------------------------------------------------------ #
    # Interpolation API
    # ------------------------------------------------------------------ #
    # These return NaN for points outside the grid (unless allow_extrapolation
    # was set at construction, in which case SciPy linearly extrapolates). They
    # never raise on out-of-bounds input, so a sweep that brushes the grid edge
    # degrades gracefully instead of aborting. Use in_bounds()/check_bounds() to
    # detect or reject extrapolation explicitly.
    def tritium_ratio_interp(self, li6: ArrayLike, be_mult: ArrayLike) -> NDArray[np.float64]:
        """Interpolate the TBR at ``(li6, be_mult)``.

        Accepts scalars or broadcastable arrays and returns an array of the same
        broadcast shape (a 0-d array for scalar inputs). This is the primary
        entry point an agent uses to emulate a Shift TBR result.
        """
        out_shape = np.shape(np.broadcast_arrays(li6, be_mult)[0])
        li6, be_mult = self._snap_to_grid(li6, be_mult)
        points = self._stack_points(li6, be_mult)
        return self._tbr_interp(points).reshape(out_shape)

    def flux_interp(self, x_coord: float, li6: ArrayLike, be_mult: ArrayLike) -> NDArray[np.float64]:
        """Interpolate neutron flux at spatial ``x_coord`` and ``(li6, be_mult)``."""
        out_shape = np.shape(np.broadcast_arrays(li6, be_mult)[0])
        cell = self.cell_index(x_coord)
        interp = self._make_interpolator(self.flux[cell])
        li6, be_mult = self._snap_to_grid(li6, be_mult)
        points = self._stack_points(li6, be_mult)
        return interp(points).reshape(out_shape)

    def zaid_cell_interp(
        self, zaid: int, position: float, li6: ArrayLike, be_mult: ArrayLike
    ) -> NDArray[np.float64]:
        """Interpolate the number density of ``zaid`` at ``position`` and ``(li6, be_mult)``."""
        out_shape = np.shape(np.broadcast_arrays(li6, be_mult)[0])
        zaid_idx = self.nuclide_index(zaid)
        cell = self.cell_index(position)
        interp = self._make_interpolator(self.number_density[zaid_idx, cell])
        li6, be_mult = self._snap_to_grid(li6, be_mult)
        points = self._stack_points(li6, be_mult)
        return interp(points).reshape(out_shape)

    @property
    def cell_centers(self) -> NDArray[np.float64]:
        """Midpoints of each mesh cell, derived from :attr:`cell_edges`."""
        return 0.5 * (self.cell_edges[1:] + self.cell_edges[:-1])

    # ------------------------------------------------------------------ #
    # Internal
    # ------------------------------------------------------------------ #
    @staticmethod
    def _stack_points(li6: ArrayLike, be_mult: ArrayLike) -> NDArray[np.float64]:
        """Pack broadcast ``(li6, be_mult)`` into the ``(N, 2)`` array SciPy wants."""
        li6_b, be_b = np.broadcast_arrays(np.asarray(li6, dtype=float), np.asarray(be_mult, dtype=float))
        return np.column_stack([li6_b.ravel(), be_b.ravel()])
