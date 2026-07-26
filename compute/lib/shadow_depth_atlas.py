"""Finite witness for the G/L/C/M shadow-depth atlas and wild boundary."""

from __future__ import annotations

from dataclasses import dataclass
from math import inf
from typing import Optional


@dataclass(frozen=True)
class ShadowDepthClass:
    """One row of the shadow-depth classification."""

    label: str
    name: str
    chirally_koszul: bool
    shadow_depth: Optional[float | int]
    algebraic_depth: Optional[float | int]
    example: str
    role: str


KOSZUL_ATLAS: tuple[ShadowDepthClass, ...] = (
    ShadowDepthClass("G", "Gaussian", True, 2, 0, "Heisenberg", "formal"),
    ShadowDepthClass("L", "Lie-transverse", True, 3, 1, "affine Kac-Moody", "finite Lie-Jacobi layer"),
    ShadowDepthClass("C", "Contact", True, 4, 2, "beta-gamma", "finite quartic contact layer"),
    ShadowDepthClass("M", "Mixed", True, inf, inf, "Virasoro and W_N", "infinite Koszul shadow depth"),
)

WILD_BOUNDARY = ShadowDepthClass(
    "W",
    "Wild",
    False,
    None,
    None,
    "Kronecker K_m, m >= 3",
    "outside the chirally Koszul locus; shadow depth undefined",
)


def koszul_shadow_depth_values() -> tuple[float | int, ...]:
    """Return the only shadow-depth values assigned on the Koszul atlas."""

    return tuple(row.shadow_depth for row in KOSZUL_ATLAS if row.shadow_depth is not None)


def class_profile(label: str) -> ShadowDepthClass:
    """Return the classification row for G/L/C/M/W."""

    normalized = label.strip().upper()
    for row in (*KOSZUL_ATLAS, WILD_BOUNDARY):
        if row.label == normalized:
            return row
    raise ValueError("label must be one of G, L, C, M, W")


def shadow_depth_defined(label: str) -> bool:
    """Whether the class has a shadow-depth value."""

    return class_profile(label).shadow_depth is not None


def wild_boundary_statement() -> dict[str, object]:
    """Return the exact A216 distinction."""

    return {
        "koszul_sequent": "A in Kosz_ch => r_sh(A) in {2,3,4,infty}",
        "wild_condition": "A notin Kosz_ch and E_r(B^ord(A)) does not collapse for all finite r",
        "wild_depth": None,
        "wild_depth_text": "r_sh(A) undefined",
        "not_open_problem": "whether infinite depth exists",
        "open_problem": (
            "whether genuine logarithmic SC^{ch,top}-algebras exist whose "
            "ordered bar spectral sequence refuses the Koszul shadow"
        ),
        "class_m_examples": ("Virasoro", "W_N"),
        "wild_model": WILD_BOUNDARY.example,
    }
