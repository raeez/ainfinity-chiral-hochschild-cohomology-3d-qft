"""Kac-Shapovalov original-complex criterion for Virasoro."""
from __future__ import annotations

from typing import Dict, Tuple


def verma_basis_vector(partition: Tuple[int, ...]) -> Dict[str, object]:
    """Return the Verma basis vector notation for a partition."""
    if any(part <= 0 for part in partition):
        raise ValueError("partition entries must be positive")
    if tuple(sorted(partition, reverse=True)) != partition:
        raise ValueError("partition must be weakly decreasing")

    level = sum(partition)
    modes = " ".join(f"L_-{part}" for part in partition) or "1"
    return {
        "partition": partition,
        "level": level,
        "basis_vector": f"{modes} |h>",
    }


def gram_matrix_profile() -> Dict[str, str]:
    """Return the Kac-Shapovalov Gram matrix formula."""
    return {
        "entry": "<h| L_{mu_m}...L_{mu_1} L_{-lambda_1}...L_{-lambda_l} |h>",
        "basis": "L_-lambda |h>, |lambda|=n",
        "walls": "h=h_{r,s}(c)",
    }


def kac_determinant_profile(level: int) -> Dict[str, object]:
    """Return the Kac determinant shape at a finite level."""
    if level < 0:
        raise ValueError("level must be nonnegative")
    return {
        "level": level,
        "determinant": "det G_n(c,h)=C_n prod_{r,s>=1,rs<=n}(h-h_{r,s}(c))^{p(n-rs)}",
        "finite_level_nonzero_off_walls": True,
        "does_not_bound_inverse_norm": True,
    }


def mk_norm_bound_profile() -> Dict[str, str]:
    """Return the finite-window KS norm estimate for transferred operations."""
    return {
        "bound": "||m_k||_{KS,n} <= C^k max_{|lambda|<=n} ||G_{|lambda|}^{-1}|| P_k(n,c,h)",
        "numerator_source": "Virasoro commutator polynomial envelope",
        "denominator_source": "inverse Kac-Shapovalov propagators",
    }


def finite_propagation_profile() -> Dict[str, str]:
    """Return the raw direct-sum finite-propagation condition."""
    return {
        "condition": "forall n exists K(n): pi_{<=n} m_k|_{A_{<=n}^{otimes k}}=0 for k>K(n)",
        "meaning": "finite conformal-weight support is preserved",
        "failure_mode": "higher transferred operations feed lower finite weights through inverse Shapovalov propagators",
    }


def original_tempered_locus_profile() -> Dict[str, Tuple[str, ...] | str]:
    """Return the KS-original tempered-locus criterion."""
    return {
        "locus": "T_orig",
        "requires": (
            "off Kac walls",
            "sum_{k>=2} ||m_k||_{KS,n}/k! < infinity for every n",
            "finite propagation for every n",
        ),
        "outside": "completed ambient: KS-rho Banach or weight-completed",
    }
