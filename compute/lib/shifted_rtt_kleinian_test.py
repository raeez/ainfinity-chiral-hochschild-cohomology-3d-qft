"""Shifted RTT pairing-annihilator and rank-one Kleinian test."""
from __future__ import annotations

from math import factorial
from typing import Dict, Tuple


def shifted_generator_profile(mu_i: int, mu_j: int) -> Dict[str, object]:
    """Return the chamber-shifted RTT current exponent."""
    exponent = -mu_i + mu_j
    return {
        "mu_i": mu_i,
        "mu_j": mu_j,
        "u_exponent": exponent,
        "formula": "T_ij^(mu)(u)=u^{-mu_i+mu_j}(delta_ij+sum_{r>=1}T_ij^(r)u^{-r})",
    }


def pairing_annihilator_profile() -> Dict[str, str]:
    """Return the shifted quotient ideal and coideal criterion."""
    return {
        "ideal": "I_mu={x in Y_hbar^RTT : <x,y>=0 for all y in Y_hbar^{RTT,vee}_{>=mu}}",
        "quotient": "Y_hat_{hbar,mu}^RTT=Y_hat_hbar^RTT/I_mu",
        "coideal_condition": "Delta(I_mu) subset I_mu tensor Y_hat + Y_hat tensor I_mu",
        "interpretation": "without the coideal condition the quotient is not a line algebra",
    }


def quantum_determinant_profile(rank: int) -> Dict[str, object]:
    """Return the shifted RTT quantum determinant pattern."""
    if rank <= 0:
        raise ValueError("rank must be positive")
    shifts = tuple(f"u-{a}hbar" for a in range(rank))
    return {
        "rank": rank,
        "term_count": factorial(rank),
        "shifts": shifts,
        "formula": (
            "qdet T(u)=sum_{sigma in S_N}(-1)^sigma "
            "prod_{j=1}^N T_{sigma(j),j}(u-(j-1)hbar)"
        ),
        "boundary_value": "qdet T(u)=P_mu(u)",
    }


def kleinian_boundary_relations(m: int) -> Dict[str, object]:
    """Return the rank-one determinantal/Casimir quotient relations."""
    if m < 2:
        raise ValueError("m must be at least 2 for an A_{m-1} singularity")
    yx_factors = tuple(f"z-{a}hbar" for a in range(m))
    xy_factors = tuple(f"z+{a}hbar" for a in range(1, m + 1))
    return {
        "m": m,
        "commutators": ("[z,x]-hbar x", "[z,y]+hbar y"),
        "yx_relation": "yx-prod_{a=0}^{m-1}(z-a hbar)",
        "xy_relation": "xy-prod_{a=1}^{m}(z+a hbar)",
        "yx_factors": yx_factors,
        "xy_factors": xy_factors,
        "associated_graded": f"C[x,y,z]/(xy-z^{m})",
        "kleinian_type": f"A_{m-1}",
    }


def shifted_rtt_candidate_passes_rank_one_test(has_coideal: bool, gr_relation: str, m: int) -> bool:
    """Return whether a candidate passes the rank-one shifted RTT test."""
    expected = kleinian_boundary_relations(m)["associated_graded"]
    return has_coideal and gr_relation == expected
