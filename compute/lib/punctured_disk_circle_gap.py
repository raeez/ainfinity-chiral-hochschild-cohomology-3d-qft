"""Punctured-disk versus circle chain-model gap.

The angular projection from D^* to S^1 is an E_1 truncation.  It is not
an equivalence of the full configuration-space operadic modules.
"""
from __future__ import annotations

from typing import Dict, Tuple


def rotation_bv_operator(arity: int) -> Dict[str, object]:
    """Return the rotation BV operator on Conf_arity(D^*)."""
    if arity <= 0:
        raise ValueError("arity must be positive")

    contractions = tuple(f"iota_partial_arg_z_{i}" for i in range(1, arity + 1))
    return {
        "operator": "Delta_rot = sum_i iota_{partial_arg_z_i}",
        "arity": arity,
        "contractions": contractions,
        "survives_on": "D^* punctured-disk chain model",
        "circle_image": "ordinary cyclic rotation",
    }


def punctured_disk_circle_gap_profile(arity: int) -> Dict[str, object]:
    """Return the scope of the D^* -> S^1 comparison."""
    if arity <= 0:
        raise ValueError("arity must be positive")

    return {
        "punctured_disk_model": (
            "B^{D^*}_bullet(A) = direct_sum C_log(FM_n(D^*)) tensor A^otimes n"
        ),
        "symmetric_model": (
            "B^{D^*,Sigma}_bullet(A) = direct_sum C(Conf_n(D^*)) "
            "tensor_{Sigma_n} A^otimes n"
        ),
        "circle_model": (
            "B^{S^1}_bullet(A) = direct_sum C(Conf_n(S^1)) "
            "tensor_{Sigma_n} A^otimes n"
        ),
        "quasi_isomorphism_scope": "angular E_1 quotient only",
        "full_operadic_module_equivalence": False,
        "gap_generators": (
            "theta_i = d arg z_i",
            "omega_ij = d log(z_i - z_j)",
            "Delta_rot = sum_i iota_{partial_arg_z_i}",
        ),
        "circle_computes": "annular trace",
        "punctured_disk_computes": "chiral centre with rotation",
        "arity": arity,
    }


def angular_projection_effect() -> Dict[str, Tuple[str, ...]]:
    """Record what radial projection preserves and what it forgets."""
    return {
        "preserves": (
            "cyclic order",
            "wraparound monodromy",
            "E_1 angular quotient",
        ),
        "forgets": (
            "radial collision data",
            "pairwise logarithmic collision forms",
            "rotation BV operator as D^* chain operation",
        ),
    }
