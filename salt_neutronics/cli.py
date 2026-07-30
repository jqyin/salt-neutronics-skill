"""Command-line interface for emulating Shift FLiBe neutronics results.

This CLI is the primary surface an AI agent (or a human) drives to obtain a
tritium breeding ratio (TBR) for an arbitrary FLiBe composition without running a
full Monte Carlo transport calculation. It interpolates the precomputed Shift
parameter study and emits a structured, JSON-friendly "mock simulation report".

Subcommands
-----------
* ``tbr``     -- single-point TBR for one composition (the headline query).
* ``sweep``   -- TBR over a range of compositions; writes a CSV table.
* ``flux``    -- neutron flux at a spatial position for one composition.
* ``density`` -- number density of a nuclide (by ZAID) at a position.

Run ``python -m salt_neutronics.cli <subcommand> --help`` for details.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from .composition import (
    NOMINAL_BEF2_MOLPCT,
    bef2_molpct_to_be_multiplier,
    be_multiplier_to_bef2_molpct,
    lif_molpct_from_bef2,
)
from .processor import (
    MAGNET_FLUX_PREFERRED_N_PER_CM2_S,
    MAGNET_FLUX_REJECT_N_PER_CM2_S,
    SHIELDING_THICKNESS_CM,
    ShiftFlibeProcessor,
    shielding_verdict,
)


# --------------------------------------------------------------------------- #
# Shared helpers
# --------------------------------------------------------------------------- #
def _resolve_multiplier(args: argparse.Namespace) -> tuple[float, float | None]:
    """Return ``(be_multiplier, bef2_molpct)`` from the parsed composition args."""
    if args.be_multiplier is not None and args.bef2 is not None:
        raise SystemExit("Specify only one of --bef2 or --be-multiplier, not both.")
    if args.be_multiplier is not None:
        bef2 = float(be_multiplier_to_bef2_molpct(args.be_multiplier, args.nominal_bef2))
        return float(args.be_multiplier), bef2
    bef2 = NOMINAL_BEF2_MOLPCT if args.bef2 is None else float(args.bef2)
    mult = float(bef2_molpct_to_be_multiplier(bef2, args.nominal_bef2))
    return mult, bef2


def _guard_single_point(
    proc: ShiftFlibeProcessor, li6: float, mult: float, nominal: float, allow: bool
) -> None:
    """Reject an out-of-grid single-point query with an actionable message."""
    if allow or proc.in_bounds(li6, mult):
        return
    bef2_lo = float(be_multiplier_to_bef2_molpct(proc.be_bounds.low, nominal))
    bef2_hi = float(be_multiplier_to_bef2_molpct(proc.be_bounds.high, nominal))
    raise SystemExit(
        "error: requested composition is outside the simulated grid; the result "
        "would not be backed by any Shift run.\n"
        f"  Valid Li-6 enrichment : [{proc.li6_bounds.low:g}, {proc.li6_bounds.high:g}]\n"
        f"  Valid BeF2 content    : [{bef2_lo:.2f}, {bef2_hi:.2f}] mol% "
        f"(Be multiplier [{proc.be_bounds.low:g}, {proc.be_bounds.high:g}])\n"
        "  Re-run within these ranges, or pass --allow-extrapolation to extrapolate (cautiously)."
    )


def _interpretation(tbr: float) -> str:
    """One-line physics reading of a TBR value for the agent to relay."""
    if not np.isfinite(tbr):
        return "TBR is undefined (request outside the simulated grid)."
    if tbr >= 1.1:
        margin = "comfortable"
    elif tbr >= 1.0:
        margin = "marginal"
    else:
        return (
            f"TBR = {tbr:.4f} < 1.0: this composition does NOT breed enough tritium "
            "for self-sufficiency."
        )
    return (
        f"TBR = {tbr:.4f} >= 1.0: tritium self-sufficient with a {margin} margin "
        "(before accounting for losses/streaming)."
    )


def _shielding_block(
    proc: ShiftFlibeProcessor,
    li6: float,
    mult: float,
    thickness_cm: float = SHIELDING_THICKNESS_CM,
) -> dict[str, Any]:
    """Magnet radiation-shielding assessment: flux at the blanket back face + a verdict."""
    flux = float(proc.shielding_flux(li6, mult, thickness_cm))
    verdict = shielding_verdict(flux)
    return {
        "magnet_flux_n_per_cm2_s": flux,
        "thickness_cm": thickness_cm,
        "reject_above_n_per_cm2_s": MAGNET_FLUX_REJECT_N_PER_CM2_S,
        "preferred_at_or_below_n_per_cm2_s": MAGNET_FLUX_PREFERRED_N_PER_CM2_S,
        "verdict": verdict,
        "acceptable": verdict != "reject",
    }


def _build_report(
    proc: ShiftFlibeProcessor, bef2: float, mult: float, li6: float, nominal_bef2: float
) -> dict[str, Any]:
    """Assemble a structured single-point TBR report."""
    in_bounds = proc.in_bounds(li6, mult)
    tbr = float(proc.tritium_ratio_interp(li6, mult))
    return {
        "quantity": "tritium_breeding_ratio",
        "input_composition": {
            "bef2_mol_percent": bef2,
            "lif_mol_percent": float(lif_molpct_from_bef2(bef2)),
            "li6_enrichment": li6,
        },
        "derived": {
            "beryllium_multiplier": mult,
            "nominal_bef2_mol_percent": nominal_bef2,
        },
        "result": {
            "tbr": tbr,
            "interpretation": _interpretation(tbr),
        },
        "shielding": _shielding_block(proc, li6, mult),
        "provenance": {
            "method": "linear interpolation of a precomputed Shift Monte Carlo table",
            "is_interpolated": in_bounds,
            "is_extrapolated": (not in_bounds),
            "data_file": str(proc.filename),
            "grid": {
                "li6_enrichment": proc.lithium6_enrichments.tolist(),
                "beryllium_multiplier": proc.beryllium_multipliers.tolist(),
            },
            "source_strength_n_per_s": proc.source_strength,
            "irradiation_time_s": proc.irradiation_time,
            "incident_neutrons": proc.source_strength * proc.irradiation_time,
        },
    }


def _emit(report: dict[str, Any], as_json: bool) -> None:
    """Print a report either as JSON (machine) or a readable summary (human)."""
    if as_json:
        print(json.dumps(report, indent=2))
        return
    comp = report["input_composition"]
    res = report["result"]
    prov = report["provenance"]
    mode = "INTERPOLATED" if prov["is_interpolated"] else "EXTRAPOLATED (outside grid!)"
    print("=" * 64)
    print("  Mock Shift neutronics result  --  Tritium Breeding Ratio")
    print("=" * 64)
    print(f"  BeF2 content        : {comp['bef2_mol_percent']:.3f} mol%")
    print(f"  LiF  content        : {comp['lif_mol_percent']:.3f} mol%")
    print(f"  Li-6 enrichment     : {comp['li6_enrichment']:.4f}")
    print(f"  Beryllium multiplier: {report['derived']['beryllium_multiplier']:.4f}")
    print("-" * 64)
    print(f"  TBR                 : {res['tbr']:.4f}   [{mode}]")
    print(f"  {res['interpretation']}")
    sh = report.get("shielding")
    if sh:
        print("-" * 64)
        print(
            f"  Magnet shield flux  : {sh['magnet_flux_n_per_cm2_s']:.3e} n/cm^2-s "
            f"@ {sh['thickness_cm']:g} cm   [{sh['verdict'].upper()}]"
        )
        if not sh["acceptable"]:
            print(f"  REJECT: flux exceeds {sh['reject_above_n_per_cm2_s']:.0e} — shielding too weak.")
    print("=" * 64)


# --------------------------------------------------------------------------- #
# Subcommand handlers
# --------------------------------------------------------------------------- #
def cmd_tbr(args: argparse.Namespace) -> int:
    proc = ShiftFlibeProcessor(args.data, allow_extrapolation=args.allow_extrapolation)
    mult, bef2 = _resolve_multiplier(args)
    _guard_single_point(proc, args.li6, mult, args.nominal_bef2, args.allow_extrapolation)
    report = _build_report(proc, bef2, mult, args.li6, args.nominal_bef2)
    _emit(report, args.json)
    if args.output:
        Path(args.output).write_text(json.dumps(report, indent=2))
        if not args.json:
            print(f"\nWrote report to {args.output}")
    return 0


def _axis(value: float | None, rng: list[float] | None, default: float) -> np.ndarray:
    """Build a sweep axis from either a fixed value or a ``LOW HIGH N`` range."""
    if rng is not None:
        low, high, num = rng[0], rng[1], int(rng[2]) if len(rng) > 2 else 50
        return np.linspace(low, high, num)
    return np.array([default if value is None else value], dtype=float)


def cmd_sweep(args: argparse.Namespace) -> int:
    proc = ShiftFlibeProcessor(args.data, allow_extrapolation=args.allow_extrapolation)

    bef2_axis = _axis(args.bef2, args.bef2_range, NOMINAL_BEF2_MOLPCT)
    li6_axis = _axis(args.li6, args.li6_range, 0.075)
    if bef2_axis.size == 1 and li6_axis.size == 1:
        raise SystemExit("sweep needs at least one of --bef2-range or --li6-range.")

    mult_axis = bef2_molpct_to_be_multiplier(bef2_axis, args.nominal_bef2)
    BEF2, LI6 = np.meshgrid(bef2_axis, li6_axis, indexing="ij")
    MULT, _ = np.meshgrid(mult_axis, li6_axis, indexing="ij")
    tbr = proc.tritium_ratio_interp(LI6, MULT)

    out_path = Path(args.output or "tbr_sweep.csv")
    header = "bef2_mol_percent,li6_enrichment,beryllium_multiplier,tbr"
    rows = np.column_stack([BEF2.ravel(), LI6.ravel(), MULT.ravel(), tbr.ravel()])
    np.savetxt(out_path, rows, delimiter=",", header=header, comments="", fmt="%.8e")

    finite = tbr[np.isfinite(tbr)]
    n_oob = int(tbr.size - finite.size)
    print(f"Wrote {rows.shape[0]} rows to {out_path}")
    if n_oob:
        print(f"Note: {n_oob}/{tbr.size} points fell outside the grid and are NaN "
              "(use --allow-extrapolation to extrapolate).")
    if finite.size:
        imax = np.unravel_index(np.nanargmax(tbr), tbr.shape)
        print(
            f"TBR range over sweep: {finite.min():.4f} .. {finite.max():.4f}; "
            f"max at BeF2={BEF2[imax]:.2f} mol%, Li-6={LI6[imax]:.3f}"
        )
    if args.plot:
        from .plotting import plot_tbr_surface

        png = plot_tbr_surface(proc, nominal_bef2=args.nominal_bef2, output_dir=out_path.parent)
        print(f"Wrote plot to {png}")
    return 0


def cmd_flux(args: argparse.Namespace) -> int:
    proc = ShiftFlibeProcessor(args.data, allow_extrapolation=args.allow_extrapolation)
    mult, bef2 = _resolve_multiplier(args)
    _guard_single_point(proc, args.li6, mult, args.nominal_bef2, args.allow_extrapolation)
    value = float(proc.flux_interp(args.position, args.li6, mult))
    report = {
        "quantity": "neutron_flux",
        "position_cm": args.position,
        "cell_index": proc.cell_index(args.position),
        "input_composition": {"bef2_mol_percent": bef2, "li6_enrichment": args.li6},
        "derived": {"beryllium_multiplier": mult},
        "result": {"flux_n_per_cm2_s": value},
    }
    _emit_simple(report, args.json, f"Neutron flux = {value:.6e} n/cm^2-s")
    return 0


def cmd_shielding(args: argparse.Namespace) -> int:
    proc = ShiftFlibeProcessor(args.data, allow_extrapolation=args.allow_extrapolation)
    mult, bef2 = _resolve_multiplier(args)
    _guard_single_point(proc, args.li6, mult, args.nominal_bef2, args.allow_extrapolation)
    block = _shielding_block(proc, args.li6, mult, args.thickness)
    report = {
        "quantity": "magnet_shielding_flux",
        "input_composition": {"bef2_mol_percent": bef2, "li6_enrichment": args.li6},
        "derived": {"beryllium_multiplier": mult},
        "result": block,
    }
    _emit_simple(
        report,
        args.json,
        f"Magnet shielding flux = {block['magnet_flux_n_per_cm2_s']:.6e} n/cm^2-s "
        f"@ {block['thickness_cm']:g} cm blanket -> {block['verdict'].upper()}",
    )
    return 0


def cmd_density(args: argparse.Namespace) -> int:
    proc = ShiftFlibeProcessor(args.data, allow_extrapolation=args.allow_extrapolation)
    mult, bef2 = _resolve_multiplier(args)
    _guard_single_point(proc, args.li6, mult, args.nominal_bef2, args.allow_extrapolation)
    value = float(proc.zaid_cell_interp(args.zaid, args.position, args.li6, mult))
    report = {
        "quantity": "number_density",
        "zaid": args.zaid,
        "position_cm": args.position,
        "cell_index": proc.cell_index(args.position),
        "input_composition": {"bef2_mol_percent": bef2, "li6_enrichment": args.li6},
        "derived": {"beryllium_multiplier": mult},
        "result": {"number_density": value},
    }
    _emit_simple(report, args.json, f"Number density (ZAID {args.zaid}) = {value:.6e}")
    return 0


def _emit_simple(report: dict[str, Any], as_json: bool, human_line: str) -> None:
    if as_json:
        print(json.dumps(report, indent=2))
    else:
        print(human_line)


# --------------------------------------------------------------------------- #
# Argument parser
# --------------------------------------------------------------------------- #
def _add_common(p: argparse.ArgumentParser, *, with_position: bool = False) -> None:
    p.add_argument("--data", type=Path, default=None, help="Path to neutronics_isotopics.h5 (default: bundled).")
    p.add_argument("--bef2", type=float, default=None, help="BeF2 content in mole percent.")
    p.add_argument("--be-multiplier", type=float, default=None, help="Beryllium multiplier (alternative to --bef2).")
    p.add_argument("--li6", type=float, default=0.075, help="Li-6 enrichment atom fraction (default: 0.075, ~natural).")
    p.add_argument("--nominal-bef2", type=float, default=NOMINAL_BEF2_MOLPCT, dest="nominal_bef2",
                   help="BeF2 mol%% that maps to multiplier 1.0 (default: 33.33).")
    p.add_argument("--allow-extrapolation", action="store_true", help="Permit (cautious) extrapolation outside the grid.")
    p.add_argument("--json", action="store_true", help="Emit a JSON report on stdout.")
    if with_position:
        p.add_argument("--position", type=float, required=True, help="Spatial position in cm.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="salt-neutronics",
        description="Emulate Shift FLiBe neutronics results (TBR) by interpolating a precomputed MC table.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_tbr = sub.add_parser("tbr", help="Single-point TBR for one composition.")
    _add_common(p_tbr)
    p_tbr.add_argument("--output", type=Path, default=None, help="Optional path to also write the JSON report.")
    p_tbr.set_defaults(func=cmd_tbr)

    p_sweep = sub.add_parser("sweep", help="Sweep TBR over composition; write CSV table.")
    _add_common(p_sweep)
    p_sweep.add_argument("--bef2-range", nargs="+", type=float, metavar="LOW HIGH [N]",
                         help="Sweep BeF2 mol%% from LOW to HIGH in N steps (default N=50).")
    p_sweep.add_argument("--li6-range", nargs="+", type=float, metavar="LOW HIGH [N]",
                         help="Sweep Li-6 enrichment from LOW to HIGH in N steps (default N=50).")
    p_sweep.add_argument("--output", type=Path, default=None, help="CSV output path (default: tbr_sweep.csv).")
    p_sweep.add_argument("--plot", action="store_true", help="Also render the TBR surface PNG.")
    p_sweep.set_defaults(func=cmd_sweep)

    p_flux = sub.add_parser("flux", help="Neutron flux at a position for one composition.")
    _add_common(p_flux, with_position=True)
    p_flux.set_defaults(func=cmd_flux)

    p_shield = sub.add_parser(
        "shielding", help="Magnet radiation-shielding flux at the blanket back face."
    )
    _add_common(p_shield)
    p_shield.add_argument(
        "--thickness", type=float, default=SHIELDING_THICKNESS_CM,
        help="Blanket thickness in cm at which to evaluate the magnet flux (default: 100 = 1 m).",
    )
    p_shield.set_defaults(func=cmd_shielding)

    p_density = sub.add_parser("density", help="Number density of a nuclide (ZAID) at a position.")
    _add_common(p_density, with_position=True)
    p_density.add_argument("--zaid", type=int, required=True, help="Nuclide ZAID, e.g. 1003 for tritium.")
    p_density.set_defaults(func=cmd_density)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, KeyError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
