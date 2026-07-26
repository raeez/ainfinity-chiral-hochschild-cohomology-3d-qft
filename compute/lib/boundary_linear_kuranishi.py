"""Boundary-linear Kuranishi reduction formulas.

This module is a small structural oracle for the local exact sector:
integrate out the linearly massive component and record the first
effective Kuranishi couplings and derived-critical bulk algebra.
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial
from typing import Callable, Dict, Tuple

Scalar = Fraction | int
Binary = Callable[[Scalar, Scalar], Scalar]
Ternary = Callable[[Scalar, Scalar, Scalar], Scalar]
Linear = Callable[[Scalar], Scalar]


def _cyclic_triples(t1: Scalar, t2: Scalar, t3: Scalar) -> Tuple[Tuple[Scalar, Scalar, Scalar], ...]:
    return ((t1, t2, t3), (t2, t3, t1), (t3, t1, t2))


def massive_i2(F2B: Binary, A_inv: Linear, t1: Scalar, t2: Scalar) -> Scalar:
    """Second massive-field coefficient i_2 = -A^{-1}F_2^B."""
    return -A_inv(F2B(t1, t2))


def massive_i3(
    F2B: Binary,
    F3B: Ternary,
    A_inv: Linear,
    t1: Scalar,
    t2: Scalar,
    t3: Scalar,
) -> Scalar:
    """Third massive-field coefficient with cyclic tree exchange."""
    exchange = sum(
        F2B(a, massive_i2(F2B, A_inv, b, c))
        for a, b, c in _cyclic_triples(t1, t2, t3)
    )
    return -A_inv(F3B(t1, t2, t3) + exchange)


def kappa2(F2C: Binary, t1: Scalar, t2: Scalar) -> Scalar:
    """First reduced Kuranishi coupling kappa_2 = F_2^C."""
    return F2C(t1, t2)


def kappa3(
    F2B: Binary,
    F2C: Binary,
    F3C: Ternary,
    A_inv: Linear,
    t1: Scalar,
    t2: Scalar,
    t3: Scalar,
) -> Scalar:
    """Cubic reduced coupling: direct cubic minus massive exchange."""
    exchange = sum(
        F2C(a, A_inv(F2B(b, c)))
        for a, b, c in _cyclic_triples(t1, t2, t3)
    )
    return F3C(t1, t2, t3) - exchange


def minimal_line_operation_profile(arity: int) -> Dict[str, object]:
    """Return the minimal pointed line operation normalization."""
    if arity < 2:
        raise ValueError("minimal Kuranishi operations start at arity 2")
    return {
        "arity": arity,
        "factor": Fraction(1, factorial(arity)),
        "has_exterior_product": arity == 2,
        "formula": (
            "m_2(lambda_t1,lambda_t2)=lambda_t1 wedge lambda_t2 "
            "+ 1/2 kappa_2(t1,t2)"
            if arity == 2
            else f"m_{arity}(lambda_t1,...,lambda_t{arity})="
            f"1/{factorial(arity)} kappa_{arity}(t1,...,t{arity})"
        ),
    }


def derived_critical_dga_profile() -> Dict[str, str]:
    """Return the coordinate model for O(dCrit(W_eff))."""
    return {
        "potential": "W_eff(u,gamma)=<gamma,kappa(u)>",
        "algebra": "k[[u,gamma]] tensor Lambda(xi,chi)",
        "differential": (
            "d_W=sum_i dW/du_i d/dxi_i + "
            "sum_a dW/dgamma_a d/dchi_a"
        ),
        "bulk_line_chain": (
            "HH(K_F,p) ~= HH(K_kappa) ~= HH(B_kappa) ~= "
            "O(dCrit(W_eff))"
        ),
        "line_not_bulk": "bulk is Hochschild/derived centre of the line algebra",
    }
