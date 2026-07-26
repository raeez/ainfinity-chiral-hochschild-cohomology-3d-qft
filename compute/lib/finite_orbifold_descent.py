"""Finite-orbifold descent profiles used by the Monster/Schellekens surface.

The helpers here are deliberately structural.  They do not compute a
local BV anomaly for a lattice action; they record the descent datum the
manuscript requires: homotopy fixed points, a Dijkgraaf-Witten class in
``H^3(BG; U(1))``, and the associator correction when that class is a
coboundary.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


Scalar = int | float | complex | Fraction


@dataclass(frozen=True)
class HomotopyFixedPointProfile:
    """Cosimplicial totalization profile for finite-group descent."""

    group: str
    algebra: str
    totalization: str
    cosimplicial_levels: tuple[str, ...]
    twisted_sector_required: bool
    requires_trivial_dw_class: bool


@dataclass(frozen=True)
class DijkgraafWittenObstructionProfile:
    """Dijkgraaf-Witten obstruction class for finite-orbifold descent."""

    group: str
    cohomology_group: str
    descent_requires_zero: bool
    corrected_associator_formula: str


def homotopy_fixed_point_profile(
    *,
    group: str = "G",
    algebra: str = "V",
    twisted_sector_required: bool = False,
) -> HomotopyFixedPointProfile:
    """Return the homotopy-fixed-point totalization used for descent."""

    levels = (
        algebra,
        f"Map({group},{algebra})",
        f"Map({group}^2,{algebra})",
        "...",
    )
    return HomotopyFixedPointProfile(
        group=group,
        algebra=algebra,
        totalization=(
            f"{algebra}^h{group} = Tot[{algebra} => "
            f"Map({group},{algebra}) => Map({group}^2,{algebra}) => ...]"
        ),
        cosimplicial_levels=levels,
        twisted_sector_required=twisted_sector_required,
        requires_trivial_dw_class=True,
    )


def corrected_associator_formula() -> str:
    """Return the beta-corrected associator formula as a manuscript string."""

    return (
        "Phi^beta_{g,h,k} = beta(g,h) beta(gh,k) "
        "beta(h,k)^(-1) beta(g,hk)^(-1) Phi_{g,h,k}"
    )


def dw_obstruction_profile(group: str = "G") -> DijkgraafWittenObstructionProfile:
    """Return the obstruction group and the zero-class descent condition."""

    return DijkgraafWittenObstructionProfile(
        group=group,
        cohomology_group=f"H^3(B{group}; U(1))",
        descent_requires_zero=True,
        corrected_associator_formula=corrected_associator_formula(),
    )


def corrected_associator_value(
    *,
    beta_g_h: Scalar,
    beta_gh_k: Scalar,
    beta_h_k: Scalar,
    beta_g_hk: Scalar,
    phi_g_h_k: Scalar,
) -> Scalar:
    """Evaluate ``beta(g,h) beta(gh,k) beta(h,k)^-1 beta(g,hk)^-1 Phi``."""

    if beta_h_k == 0 or beta_g_hk == 0:
        raise ValueError("U(1) beta-values must be non-zero")
    return beta_g_h * beta_gh_k / beta_h_k / beta_g_hk * phi_g_h_k


def schellekens_row_partition() -> dict[str, int]:
    """The Schellekens ``c=24`` partition used by the descent criterion."""

    return {
        "Niemeier": 24,
        "V1=0": 1,
        "nonlat": 46,
    }
