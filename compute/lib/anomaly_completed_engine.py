"""Anomaly-Completed Holography Engine: transgression algebras and Clifford completion.

Implements the anomaly-completed holographic constructions from Vol II Part VII.

The anomaly-completed holographic programme extends the Koszul duality
framework to theories with gravitational anomalies. The key objects:

1. **Transgression algebra** B_Theta: Given a dga B and a Maurer-Cartan
   element Theta in B, the transgression algebra is B_Theta = B * k<eta>
   with eta*b = (-1)^{|b|} b*eta and d(eta) = Theta. This extends B by
   a new generator eta whose differential is the anomaly.

2. **Secondary anomaly** u = eta^2 in B_Theta: The square of the
   transgression generator, which is nonzero because eta is odd-degree
   (in general). The degree of u = 2 * deg(eta).

3. **Neutralization**: A B_Theta-module M is neutralizable if the
   anomaly Theta can be trivialized on M. The obstruction lives in
   Ext^2_B(M, M), and if neutralizable, the moduli of neutralizations
   form an affine space over Ext^1_B(M, M).

4. **Genus-Clifford completion**: At genus g, the gravitational anomaly
   requires g Clifford factors (one per handle), multiplying the
   dimension by 2^g. This is the genus-g line operator algebra.

5. **Holographic encoding**: The transgression algebra packages the
   entire holographic dictionary: bulk = center of B_Theta, boundary
   = B, anomaly = Theta, line operators = modules over B_Theta.

References:
  Vol II: anomaly_completed_core.tex (Part VII)
  Vol II: 3d_gravity.tex (Movement VI, genus-Clifford)
  Vol I: concordance.tex (Theorem C, complementarity)
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any, Dict, List, Optional, Tuple

from sympy import (
    Symbol, Rational, simplify, expand, S, symbols, factorial,
    binomial,
)


# =========================================================================
# 1. TRANSGRESSION ALGEBRA
# =========================================================================

def transgression_algebra(
    B_dim: int,
    theta_degree: int,
) -> Dict[str, Any]:
    """Construct the transgression algebra B_Theta = B * k<eta> / relations.

    Given a dga B of dimension B_dim and a Maurer-Cartan element Theta
    of degree theta_degree, the transgression algebra B_Theta is the
    extension of B by a generator eta with:

    - deg(eta) = theta_degree - 1  (so that d(eta) = Theta has the right degree)
    - eta * b = (-1)^{|b|} * b * eta  (graded commutation)
    - d(eta) = Theta

    As a graded vector space, B_Theta = B tensor k<eta> = B + B*eta,
    so dim(B_Theta) = 2 * dim(B).

    Parameters
    ----------
    B_dim : int
        Dimension of the base dga B. Must be >= 1.
    theta_degree : int
        Degree of the Maurer-Cartan element Theta.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'dim_B_Theta': 2 * B_dim (dimension of transgression algebra)
        - 'eta_degree': theta_degree - 1
        - 'theta_degree': theta_degree
        - 'B_dim': B_dim
        - 'commutation_sign': -1 if eta has odd degree, +1 if even
        - 'is_clifford_type': True if eta^2 != 0 in general
    """
    if B_dim < 1:
        raise ValueError(f"B_dim must be >= 1, got {B_dim}")

    eta_degree = theta_degree - 1
    # Sign in eta*b = (-1)^{|b|} b*eta depends on eta's degree parity
    # but the relation is graded: eta*b = (-1)^{|b|} b*eta
    commutation_sign = (-1) ** eta_degree

    return {
        'dim_B_Theta': 2 * B_dim,
        'eta_degree': eta_degree,
        'theta_degree': theta_degree,
        'B_dim': B_dim,
        'commutation_sign': commutation_sign,
        'is_clifford_type': eta_degree % 2 == 1,  # odd degree => eta^2 can be nonzero
    }


# =========================================================================
# 2. SECONDARY ANOMALY
# =========================================================================

def secondary_anomaly_u(eta_degree: int) -> Dict[str, Any]:
    """Compute the secondary anomaly u = eta^2 in B_Theta.

    The square of the transgression generator eta has degree 2 * eta_degree.
    For odd-degree eta, u = eta^2 is potentially nonzero (it is the
    secondary characteristic class of the anomaly). For even-degree eta,
    eta^2 = 0 by graded commutativity.

    Parameters
    ----------
    eta_degree : int
        Degree of the transgression generator eta.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'u_degree': 2 * eta_degree (degree of u = eta^2)
        - 'eta_degree': eta_degree
        - 'is_nonzero': True if eta has odd degree (u != 0 generically)
        - 'is_central': True (u commutes with everything in B_Theta)
    """
    return {
        'u_degree': 2 * eta_degree,
        'eta_degree': eta_degree,
        'is_nonzero': eta_degree % 2 == 1,
        'is_central': True,  # u = eta^2 is always central in B_Theta
    }


# =========================================================================
# 3. NEUTRALIZATION OBSTRUCTION
# =========================================================================

def neutralization_obstruction_degree() -> int:
    """Obstruction to neutralizing Theta on a module M lives in Ext^2_B(M,M).

    The anomaly Theta acts on a B-module M via the induced B_Theta
    action. Neutralization means finding a trivialization of this
    action, i.e., finding eta_M such that d(eta_M) = Theta|_M.

    The obstruction is a class in Ext^2_B(M, M): the first obstruction
    (Ext^1) parameterizes infinitesimal neutralizations, and the
    second-order obstruction (Ext^2) is the Massey-type product that
    determines whether these can be integrated.

    Returns
    -------
    int
        The Ext degree of the obstruction: always 2.
    """
    return 2


# =========================================================================
# 4. NEUTRALIZATION MODULI DIMENSION
# =========================================================================

def neutralization_moduli_dim(ext1_dim: int) -> Dict[str, Any]:
    """If neutralizable, the moduli of neutralizations = affine space over Ext^1.

    When the obstruction in Ext^2 vanishes (i.e., the module M is
    neutralizable), the space of neutralizations is a torsor over
    Ext^1_B(M, M). In particular, if Ext^1 = 0, the neutralization
    is unique.

    Parameters
    ----------
    ext1_dim : int
        Dimension of Ext^1_B(M, M). Must be >= 0.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'moduli_dim': ext1_dim (dimension of the moduli space)
        - 'ext1_dim': ext1_dim
        - 'is_rigid': ext1_dim == 0 (unique neutralization)
        - 'moduli_type': 'affine' (always an affine space over Ext^1)
    """
    if ext1_dim < 0:
        raise ValueError(f"ext1_dim must be >= 0, got {ext1_dim}")

    return {
        'moduli_dim': ext1_dim,
        'ext1_dim': ext1_dim,
        'is_rigid': ext1_dim == 0,
        'moduli_type': 'affine',
    }


# =========================================================================
# 5. GENUS-CLIFFORD COMPLETION
# =========================================================================

def genus_clifford_completion(
    g: int,
    B_dim: int,
) -> Dict[str, Any]:
    """After g genus-Clifford completions, dimension multiplies by 2^g.

    At genus g, the gravitational anomaly kappa * omega_g introduces
    curvature into the bar complex. Each genus handle contributes one
    Clifford factor to the line operator algebra, multiplying the
    dimension by 2. After g handles, the total Clifford factor is 2^g.

    This is the algebraic shadow of the genus-g partition function:
    the trace over the genus-g line operators factors as
    Z_g = str(Theta^g) where str is the super-trace over the
    Clifford algebra Cl(2g).

    Parameters
    ----------
    g : int
        Genus (number of handles). Must be >= 0.
    B_dim : int
        Dimension of the base algebra B. Must be >= 1.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'completed_dim': B_dim * 2^g
        - 'clifford_factor': 2^g
        - 'genus': g
        - 'B_dim': B_dim
        - 'clifford_rank': 2 * g (rank of Cl(2g))
        - 'supertrace_sign': (-1)^g (sign in the super-trace formula)
    """
    if g < 0:
        raise ValueError(f"Genus must be >= 0, got {g}")
    if B_dim < 1:
        raise ValueError(f"B_dim must be >= 1, got {B_dim}")

    clifford_factor = 2 ** g

    return {
        'completed_dim': B_dim * clifford_factor,
        'clifford_factor': clifford_factor,
        'genus': g,
        'B_dim': B_dim,
        'clifford_rank': 2 * g,
        'supertrace_sign': (-1) ** g,
    }
