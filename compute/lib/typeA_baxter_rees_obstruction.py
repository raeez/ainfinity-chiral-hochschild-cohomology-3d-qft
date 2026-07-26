"""Type-A Baxter-Rees RTT obstruction tower.

This module records the finite structural formulas used by the
Type-A Baxter-Rees appendix: RTT components, Rees filtration, weightwise
continuation, and the first obstruction equations.
"""
from __future__ import annotations

from typing import Dict, Tuple


def rtt_component_relation() -> Dict[str, str]:
    """Return the RTT matrix and component relations."""
    return {
        "generating_series": (
            "T_ij(u)=delta_ij+sum_{r>=1} T_ij^(r) u^{-r}; "
            "T(u)=sum E_ij tensor T_ij(u)"
        ),
        "r_matrix": "R(u)=1-hbar P/u",
        "matrix_relation": "R(u-v)T_1(u)T_2(v)=T_2(v)T_1(u)R(u-v)",
        "component_relation": (
            "(u-v)[T_ij(u),T_kl(v)]="
            "hbar(T_kj(u)T_il(v)-T_kj(v)T_il(u))"
        ),
    }


def baxter_rees_family_profile() -> Dict[str, str]:
    """Return the Rees algebra and its two fibres."""
    return {
        "weight": "wt(T_ij^(r))=r",
        "filtration": "F_d Y_hbar = span of RTT monomials of total weight <= d",
        "rees_algebra": "R_beta Y_hbar(gl_N)=direct_sum_{d>=0} beta^d F_d Y_hbar(gl_N)",
        "formal_family": "Spf completed R_beta Y_hbar(gl_N) -> Spf C[[beta]]",
        "generic_fiber": "RTT Yangian Y_hbar(gl_N)",
        "special_fiber": "associated graded RTT boundary algebra",
        "boundary_packet": "asymptotic and negative-prefundamental Baxter packet",
    }


def weightwise_continuation_profile(weight: int) -> Dict[str, object]:
    """Return finite beta-polynomial data modulo weights greater than weight."""
    if weight < 0:
        raise ValueError("weight must be nonnegative")

    r_terms = tuple(f"beta^{d} R_{d}^(<= {weight})(u)" for d in range(weight + 1))
    theta_terms = tuple(f"beta^{d} Theta_{d}^(<= {weight})" for d in range(weight + 1))
    return {
        "weight_window": weight,
        "finite_R_degree_bound": weight,
        "finite_Theta_degree_bound": weight,
        "r_terms": r_terms,
        "theta_terms": theta_terms,
        "rtt_scope": "satisfies RTT modulo F_{>w}",
        "reason": "positive RTT weights make each fixed window finite",
    }


def rtt_obstruction_tower_profile() -> Dict[str, object]:
    """Return the first two equations in the RTT deformation complex."""
    return {
        "complex": "C_RTT^bullet(Y_Bax^RTT|_{beta=0})",
        "differential": "d_{Theta_0}=[Theta_0,-]",
        "tangent": "dotTheta = (d/d beta)|_{beta=0} Theta_Phi(beta)",
        "tangent_equation": "d_{Theta_0} dotTheta = 0",
        "second_obstruction": "1/2[dotTheta,dotTheta]+d_{Theta_0}Theta_2",
        "cohomology": "H_RTT^2",
        "tower_meaning": "boundary tensor geometry is the RTT obstruction tower",
    }


def beta_polynomial_terms(symbol: str, degree: int) -> Tuple[str, ...]:
    """Return beta-polynomial terms up to a finite degree."""
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    return tuple(f"beta^{d} {symbol}_{d}" for d in range(degree + 1))
