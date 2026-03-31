"""Celestial Boundary Transfer Engine: homotopy transfer and Witt algebra.

Implements the celestial boundary transfer computations from Vol II Part III.

The celestial boundary transfer theorem describes how the bulk A-infinity
structure descends to the boundary via homotopy transfer theory (HTT).
The key objects:

1. **Tree-level transfer**: At arity n, the transferred n-ary operation
   m_n^{tr} receives contributions from binary trees with n leaves.
   The count is the Catalan number C_{n-1}.

2. **Single-particle reduction**: In the filtration quotient F^0/F^1,
   only the linearized (single-particle) contributions survive. All
   higher operations vanish mod F^1, recovering the free-field limit.

3. **Obstruction classes**: The r-th obstruction Ob_r to extending
   the transfer from arity r to arity r+1 lives in H^1(g, delta_0),
   the first cohomology of the linearized complex.

4. **Gauge changes**: At order r, a gauge transformation theta_r
   changes the transferred operation b_r by the coboundary delta_0(theta_r).
   The obstruction [Ob_r] is the residual class in cohomology.

5. **Airy-Witt operators**: D_m = -d/dt * t^{2m+2} on C[t] realize
   the Witt algebra [D_m, D_n] = (m-n) D_{m+n}. These control the
   celestial soft algebra at tree level.

References:
  Vol II: celestial_boundary_transfer_core.tex (Part III)
  Vol I: concordance.tex (Theorem A, homotopy transfer)
  Vol I: higher_genus_modular_koszul.tex (shadow tower)
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any, Dict, List, Optional, Tuple

from sympy import (
    Symbol, Rational, simplify, expand, S, symbols, factorial,
    binomial, Matrix, zeros,
)


# =========================================================================
# 1. HOMOTOPY TRANSFER TREE COUNT
# =========================================================================

def homotopy_transfer_tree_count(n: int) -> int:
    """Count tree-level contributions to the transferred n-ary operation.

    At arity n, the homotopy transfer theorem (HTT) builds the
    transferred A-infinity operation m_n^{tr} as a sum over planar
    binary trees with n leaves. The number of such trees is the
    Catalan number C_{n-1}.

    C_k = (2k)! / ((k+1)! * k!)

    Parameters
    ----------
    n : int
        Arity of the transferred operation. Must be >= 1.

    Returns
    -------
    int
        The Catalan number C_{n-1}.
    """
    if n < 1:
        raise ValueError(f"Arity must be >= 1, got {n}")
    k = n - 1
    # Catalan number C_k = binomial(2k, k) / (k+1)
    return int(binomial(2 * k, k) / (k + 1))


# =========================================================================
# 2. SINGLE-PARTICLE REDUCTION
# =========================================================================

def single_particle_reduction(
    transferred_ops: Dict[int, Any],
    filtration: int = 1,
) -> bool:
    """Single-particle transfer theorem: linearized reduction.

    In the filtration quotient F^0/F^{filtration}, only the linearized
    (single-particle) contributions survive. All higher A-infinity
    operations m_n^{tr} with n >= 3 vanish mod F^1.

    The physical content: at the single-particle level, the celestial
    boundary theory is free (no interactions). All interactions arise
    from multi-particle contributions in F^1 and deeper.

    Parameters
    ----------
    transferred_ops : dict
        Dictionary mapping arity n to the transferred operation value.
        Operations at arity >= 3 should vanish mod F^1.
    filtration : int
        Filtration level. Default 1 (standard single-particle reduction).

    Returns
    -------
    bool
        True if all operations of arity >= 3 vanish in the quotient F^0/F^1.
    """
    if filtration < 1:
        raise ValueError(f"Filtration level must be >= 1, got {filtration}")

    for arity, value in transferred_ops.items():
        if arity >= 3:
            # In the single-particle quotient, higher ops must vanish
            if value != 0:
                return False
    return True


# =========================================================================
# 3. OBSTRUCTION CLASS DEGREE
# =========================================================================

def obstruction_class_degree(r: int) -> Dict[str, Any]:
    """The r-th obstruction class Ob_r lives in H^1(g, delta_0).

    The obstruction to extending the homotopy transfer from arity r
    to arity r+1 is a class [Ob_r] in the first cohomology of the
    linearized differential delta_0. This cohomology group carries:

    - Cohomological degree 1 (one step above the cocycles)
    - Internal degree determined by the arity: the obstruction at
      arity r involves (r+1)-fold compositions, so the internal
      degree is (r+1) - 2 = r - 1 (subtracting 2 for the output
      and the differential).

    Parameters
    ----------
    r : int
        Arity level. Must be >= 2.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'cohomological_degree': always 1 (lives in H^1)
        - 'internal_degree': r - 1
        - 'total_degree': r (sum of cohomological and internal)
        - 'arity': r
        - 'target_arity': r + 1 (the arity being extended to)
    """
    if r < 2:
        raise ValueError(f"Arity level must be >= 2, got {r}")

    return {
        'cohomological_degree': 1,
        'internal_degree': r - 1,
        'total_degree': r,
        'arity': r,
        'target_arity': r + 1,
    }


# =========================================================================
# 4. GAUGE CHANGE AT ORDER r
# =========================================================================

def gauge_change_order_r(r: int) -> Dict[str, Any]:
    """At order r, gauge transformation changes b_r by delta_0(theta_r).

    The gauge freedom at each order r in the homotopy transfer is
    parameterized by a chain homotopy theta_r. The effect on the
    transferred operation b_r is:

        b_r  ->  b_r + delta_0(theta_r)

    where delta_0 is the linearized differential. The obstruction
    [Ob_r] is well-defined in H^1(g, delta_0) precisely because
    gauge changes shift b_r by exact terms.

    Parameters
    ----------
    r : int
        Order of the gauge transformation. Must be >= 2.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'order': r
        - 'gauge_parameter_degree': r - 2 (degree of theta_r)
        - 'shift_by': 'delta_0(theta_r)' (the coboundary)
        - 'residual_class': 'H^1(g, delta_0)' (the obstruction group)
        - 'gauge_freedom_dim': r - 1 (dimension of gauge parameter space
          at arity r, from the space of chain homotopies)
    """
    if r < 2:
        raise ValueError(f"Order must be >= 2, got {r}")

    return {
        'order': r,
        'gauge_parameter_degree': r - 2,
        'shift_by': 'delta_0(theta_r)',
        'residual_class': 'H^1(g, delta_0)',
        'gauge_freedom_dim': r - 1,
    }


# =========================================================================
# 5. AIRY-WITT OPERATORS
# =========================================================================

def airy_witt_operator(m: int, max_degree: int) -> Matrix:
    """Construct the Witt algebra operator L_m on the polynomial ring C[t].

    L_m = -t^{m+1} d/dt

    acting on the monomial basis {1, t, t^2, ..., t^{max_degree}}.

    The Witt algebra relation is:
        [L_m, L_n] = (m - n) L_{m+n}

    These operators control the celestial soft algebra at tree level:
    the soft graviton algebra is a quotient of the Witt algebra, and
    the Airy structure on the moduli space is governed by these operators.

    Parameters
    ----------
    m : int
        Index of the Witt generator. Can be any integer.
    max_degree : int
        Maximum polynomial degree for the truncated representation.
        Must be >= 0.

    Returns
    -------
    sympy.Matrix
        (max_degree+1) x (max_degree+1) matrix representing L_m
        in the basis {1, t, t^2, ..., t^{max_degree}}.

    Notes
    -----
    L_m(t^k) = -t^{m+1} * k * t^{k-1} = -k * t^{m+k}

    So L_m maps t^k to -k * t^{m+k}, which is nonzero in the
    truncation only if m+k <= max_degree and k >= 1.
    """
    if max_degree < 0:
        raise ValueError(f"max_degree must be >= 0, got {max_degree}")

    size = max_degree + 1
    mat = zeros(size, size)

    for k in range(size):
        # L_m(t^k) = -k * t^{m+k}
        output_degree = m + k
        coeff = -k  # vanishes for k=0 (constants are killed)

        if 0 <= output_degree < size and coeff != 0:
            mat[output_degree, k] = S(coeff)

    return mat


# =========================================================================
# 6. WITT COMMUTATION CHECK
# =========================================================================

def witt_commutation_check(
    m1: int,
    m2: int,
    max_degree: int,
) -> Dict[str, Any]:
    """Verify the Witt algebra relation [L_{m1}, L_{m2}] = (m1 - m2) L_{m1+m2}.

    Constructs the three matrices L_{m1}, L_{m2}, L_{m1+m2} in the
    truncated polynomial representation up to degree max_degree, and
    checks whether the commutator equals (m1 - m2) times L_{m1+m2}.

    Parameters
    ----------
    m1, m2 : int
        Indices of the Witt generators.
    max_degree : int
        Maximum polynomial degree for the truncated representation.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'commutator': the matrix [D_{m1}, D_{m2}]
        - 'expected': (m1 - m2) * D_{m1+m2}
        - 'difference': commutator - expected
        - 'is_zero': True if the difference vanishes
        - 'm1': m1, 'm2': m2, 'max_degree': max_degree
    """
    D_m1 = airy_witt_operator(m1, max_degree)
    D_m2 = airy_witt_operator(m2, max_degree)
    D_sum = airy_witt_operator(m1 + m2, max_degree)

    commutator = D_m1 * D_m2 - D_m2 * D_m1
    expected = (m1 - m2) * D_sum
    difference = commutator - expected

    return {
        'commutator': commutator,
        'expected': expected,
        'difference': difference,
        'is_zero': difference.equals(zeros(*difference.shape)),
        'm1': m1,
        'm2': m2,
        'max_degree': max_degree,
    }
