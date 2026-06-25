"""Plotting helpers for Shift FLiBe neutronics results.

All functions take a :class:`~salt_neutronics.processor.ShiftFlibeProcessor`, write
a PNG to ``output_dir``, and return its path. A non-interactive Matplotlib backend
(``Agg``) is selected on import so the figures render on headless HPC compute
nodes without a display.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe; must precede pyplot import
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .composition import be_multiplier_to_bef2_molpct, NOMINAL_BEF2_MOLPCT  # noqa: E402
from .processor import ShiftFlibeProcessor  # noqa: E402

__all__ = ["plot_tbr_surface", "plot_flux_profiles", "plot_nuclide_density_surface"]


def _ensure_dir(output_dir: str | Path) -> Path:
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    return path


def plot_tbr_surface(
    proc: ShiftFlibeProcessor,
    *,
    nominal_bef2: float = NOMINAL_BEF2_MOLPCT,
    mark: tuple[float, float] | None = None,
    n_grid: int = 80,
    output_dir: str | Path = ".",
    filename: str = "tbr_surface.png",
) -> Path:
    """Heatmap of interpolated TBR over BeF2 mol% and Li-6 enrichment.

    Parameters
    ----------
    mark:
        Optional ``(bef2_mol_percent, li6_enrichment)`` point to annotate -- use
        this to show where a queried composition lands on the surface.
    """
    bef2_lo, bef2_hi = (
        be_multiplier_to_bef2_molpct(proc.be_bounds.low, nominal_bef2),
        be_multiplier_to_bef2_molpct(proc.be_bounds.high, nominal_bef2),
    )
    bef2 = np.linspace(float(bef2_lo), float(bef2_hi), n_grid)
    li6 = np.linspace(proc.li6_bounds.low, proc.li6_bounds.high, n_grid)
    mult = bef2 / nominal_bef2
    BEF2, LI6 = np.meshgrid(bef2, li6, indexing="ij")
    MULT, _ = np.meshgrid(mult, li6, indexing="ij")
    tbr = proc.tritium_ratio_interp(LI6, MULT)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    mesh = ax.pcolormesh(BEF2, LI6, tbr, shading="auto", cmap="viridis")
    cs = ax.contour(BEF2, LI6, tbr, levels=8, colors="white", linewidths=0.6, alpha=0.7)
    ax.clabel(cs, inline=True, fontsize=8, fmt="%.2f")
    fig.colorbar(mesh, ax=ax, label="Tritium breeding ratio (atoms / source n)")
    ax.set_xlabel("BeF$_2$ content (mol%)")
    ax.set_ylabel("Li-6 enrichment (atom fraction)")
    ax.set_title("Tritium breeding ratio vs. FLiBe composition")
    if mark is not None:
        ax.plot(mark[0], mark[1], marker="*", color="red", markersize=16,
                markeredgecolor="white", label="query")
        ax.legend(loc="best")

    out = _ensure_dir(output_dir) / filename
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_flux_profiles(
    proc: ShiftFlibeProcessor,
    *,
    be_multiplier: float | None = None,
    output_dir: str | Path = ".",
    filename: str = "flux_profiles.png",
) -> Path:
    """Neutron flux vs. spatial position, one curve per Li-6 enrichment."""
    if be_multiplier is None:
        # Default to the multiplier closest to nominal (1.0).
        mult_idx = int(np.argmin(np.abs(proc.beryllium_multipliers - 1.0)))
    else:
        mult_idx = int(np.argmin(np.abs(proc.beryllium_multipliers - be_multiplier)))
    be_value = proc.beryllium_multipliers[mult_idx]

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for idx, li_enr in enumerate(proc.lithium6_enrichments):
        ax.plot(proc.cell_centers, proc.flux[:, idx, mult_idx],
                marker="o", markersize=3, label=f"Li-6 = {100 * li_enr:.0f}%")
    ax.set_yscale("log")
    ax.set_xlabel("Spatial position (cm)")
    ax.set_ylabel(r"Neutron flux (n/cm$^2$-s)")
    ax.set_title(f"Neutron flux profile (Be multiplier = {be_value:g})")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)

    out = _ensure_dir(output_dir) / filename
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_nuclide_density_surface(
    proc: ShiftFlibeProcessor,
    zaid: int,
    position: float,
    *,
    nominal_bef2: float = NOMINAL_BEF2_MOLPCT,
    output_dir: str | Path = ".",
    filename: str | None = None,
) -> Path:
    """3-D surface of a nuclide's number density over BeF2 mol% and Li-6."""
    zaid_idx = proc.nuclide_index(zaid)
    cell = proc.cell_index(position)
    nd = proc.number_density[zaid_idx, cell]  # (enrichment, multiplier)

    bef2 = be_multiplier_to_bef2_molpct(proc.beryllium_multipliers, nominal_bef2)
    BEF2, LI6 = np.meshgrid(bef2, proc.lithium6_enrichments)

    fig, ax = plt.subplots(figsize=(7.5, 5.5), subplot_kw={"projection": "3d"})
    ax.plot_surface(BEF2, LI6, nd, cmap="Blues", edgecolor="none")
    ax.set_xlabel("BeF$_2$ content (mol%)")
    ax.set_ylabel("Li-6 enrichment")
    ax.set_zlabel("Number density")
    ax.set_title(f"Number density of ZAID {zaid} at x = {position:g} cm")
    ax.view_init(elev=28, azim=-130)

    out = _ensure_dir(output_dir) / (filename or f"density_zaid_{zaid}.png")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out
