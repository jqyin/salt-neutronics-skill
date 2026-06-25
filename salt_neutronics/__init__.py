"""Emulate Shift FLiBe neutronics results by interpolating a precomputed MC table.

Public API
----------
- :class:`~salt_neutronics.processor.ShiftFlibeProcessor` -- load and interpolate
  the precomputed parameter study.
- :mod:`salt_neutronics.composition` -- map FLiBe ``BeF2`` mole percent to the
  beryllium multiplier axis.
- :mod:`salt_neutronics.plotting` -- render TBR / flux / density figures.
- :mod:`salt_neutronics.cli` -- command-line interface (``python -m salt_neutronics``).
"""

from __future__ import annotations

from .composition import (
    NOMINAL_BEF2_MOLPCT,
    bef2_molpct_to_be_multiplier,
    be_multiplier_to_bef2_molpct,
)
from .processor import DEFAULT_DATA_FILE, TRITIUM_ZAID, GridBounds, ShiftFlibeProcessor

__version__ = "1.0.0"

__all__ = [
    "ShiftFlibeProcessor",
    "GridBounds",
    "TRITIUM_ZAID",
    "DEFAULT_DATA_FILE",
    "NOMINAL_BEF2_MOLPCT",
    "bef2_molpct_to_be_multiplier",
    "be_multiplier_to_bef2_molpct",
]
