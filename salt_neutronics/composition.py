"""Mapping between FLiBe salt composition and the beryllium multiplier axis.

The Shift parameter study was scanned over a dimensionless *beryllium multiplier*
-- a scale factor on the nominal beryllium atom density. Users and agents,
however, naturally think in terms of **salt composition**, e.g. the mole percent
of ``BeF2`` in the binary ``LiF`` + ``BeF2`` (FLiBe) salt.

Standard "FLiBe" is the eutectic 2 ``LiF`` : 1 ``BeF2`` mixture, i.e.
``2 / 3 = 66.67 mol% LiF`` and ``1 / 3 = 33.33 mol% BeF2``. We define this
nominal composition to correspond to a beryllium multiplier of 1.0.

Modeling assumption (chosen with the user)
------------------------------------------
The beryllium multiplier is taken to be **proportional to the BeF2 mole
fraction**, i.e. proportional to the number of beryllium atoms in the salt::

    be_multiplier = (mol% BeF2) / (nominal mol% BeF2)

with ``nominal mol% BeF2 = 33.33``. This is exact if the beryllium *atom density*
scales linearly with BeF2 content (constant total atom density). It ignores
second-order changes in the salt molar volume as the LiF:BeF2 ratio changes; if a
density-corrected relation is needed, override ``nominal_bef2_molpct`` or replace
this mapping. The scanned multiplier grid (0.9-1.4) then corresponds to roughly
30-46.7 mol% BeF2, i.e. compositions bracketing the nominal eutectic.
"""

from __future__ import annotations

from numpy.typing import ArrayLike, NDArray
import numpy as np

__all__ = [
    "NOMINAL_BEF2_MOLPCT",
    "bef2_molpct_to_be_multiplier",
    "be_multiplier_to_bef2_molpct",
    "lif_molpct_from_bef2",
]

#: Mole percent of BeF2 in nominal (eutectic 2:1) FLiBe, mapped to multiplier 1.0.
NOMINAL_BEF2_MOLPCT = 100.0 / 3.0  # 33.333... mol%


def bef2_molpct_to_be_multiplier(
    bef2_molpct: ArrayLike,
    nominal_bef2_molpct: float = NOMINAL_BEF2_MOLPCT,
) -> NDArray[np.float64]:
    """Convert ``BeF2`` mole percent to the beryllium multiplier used by the data.

    Parameters
    ----------
    bef2_molpct:
        Mole percent of ``BeF2`` in the FLiBe salt (0-100). Scalar or array.
    nominal_bef2_molpct:
        Composition that maps to multiplier 1.0. Defaults to eutectic FLiBe
        (33.33 mol% ``BeF2``).

    Returns
    -------
    ndarray
        Beryllium multiplier(s).
    """
    bef2 = np.asarray(bef2_molpct, dtype=float)
    if np.any((bef2 < 0) | (bef2 > 100)):
        raise ValueError("BeF2 mole percent must be within [0, 100].")
    if nominal_bef2_molpct <= 0:
        raise ValueError("nominal_bef2_molpct must be positive.")
    return bef2 / nominal_bef2_molpct


def be_multiplier_to_bef2_molpct(
    be_multiplier: ArrayLike,
    nominal_bef2_molpct: float = NOMINAL_BEF2_MOLPCT,
) -> NDArray[np.float64]:
    """Inverse of :func:`bef2_molpct_to_be_multiplier`."""
    mult = np.asarray(be_multiplier, dtype=float)
    return mult * nominal_bef2_molpct


def lif_molpct_from_bef2(bef2_molpct: ArrayLike) -> NDArray[np.float64]:
    """Complementary ``LiF`` mole percent for a binary ``LiF`` + ``BeF2`` salt."""
    bef2 = np.asarray(bef2_molpct, dtype=float)
    return 100.0 - bef2
