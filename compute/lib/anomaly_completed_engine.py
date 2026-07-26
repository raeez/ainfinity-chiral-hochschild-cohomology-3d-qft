"""Anomaly-Completed Holography Engine: transgression algebras and Clifford completion.

Implements the anomaly-completed holographic constructions from Vol II Part VII.

The anomaly-completed holographic programme extends the Koszul duality
framework to theories with gravitational anomalies. The key objects:

1. **Transgression algebra** B_Theta: In a genuinely noncommutative bar
   model the Ore relation is eta*b = (-1)^{|b|} b*eta + iota_Theta(b),
   with d(eta) = Theta. When the anomaly contraction iota_Theta vanishes
   and Theta is central, this specializes to the strict central-shadow
   transgression algebra B_Theta = B * k<eta> with
   eta*b = (-1)^{|b|} b*eta and d(eta) = Theta. This central-shadow
   convention does not twist the differential by d + [eta,-].

   For a genuinely curved complex with d^2 = [Theta,-], the flat
   anomaly completion instead uses the curved Maurer-Cartan equation
   d eta + 1/2[eta,eta] + Theta = 0, equivalently
   d eta + eta^2 + Theta = 0 in the associative convention.

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

def ore_transgression_relation(has_contraction: bool = True) -> Dict[str, Any]:
    """Return the noncommutative Ore transgression relation."""
    relation = "eta b - (-1)^|b| b eta"
    if has_contraction:
        relation += " - iota_Theta(b)"

    return {
        'algebra': 'B_Theta = B<eta>/(relation)',
        'relation': relation,
        'eta_degree': 1,
        'd_eta': 'Theta',
        'ore_automorphism': 'sigma(b)=(-1)^|b| b',
        'ore_sigma_derivation': 'iota_Theta' if has_contraction else '0',
        'central_shadow_specialization': not has_contraction,
    }


def transgression_algebra(
    B_dim: int,
    theta_degree: int,
    eta_power_cutoff: Optional[int] = None,
) -> Dict[str, Any]:
    """Construct the strict central-shadow transgression algebra.

    Given a dga B of dimension B_dim and a closed central element Theta
    of degree theta_degree, the transgression algebra B_Theta is the
    extension of B by a generator eta with:

    - deg(eta) = theta_degree - 1  (so that d(eta) = Theta has the right degree)
    - eta * b = (-1)^{|b|} * b * eta  (graded commutation)
    - d(eta) = Theta

    This is the strict central-shadow convention. For a curved twist
    with d^2 = [Theta,-], use curved_mc_twist_data(): flatness requires
    d eta + 1/2[eta,eta] + Theta = 0.

    As a graded left B-module, B_Theta has basis {eta^n | n >= 0}.
    It is therefore infinite-rank over B unless a finite eta-power
    cutoff is imposed.  The old exterior value 2 * dim(B) is only the
    eta <= 1 truncation.

    Parameters
    ----------
    B_dim : int
        Dimension of the base dga B. Must be >= 1.
    theta_degree : int
        Degree of the closed central element Theta.
    eta_power_cutoff : int | None
        Optional maximum eta power retained.  If None, no finite
        dimension is reported.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'dim_B_Theta': None unless eta_power_cutoff is supplied
        - 'eta_degree': theta_degree - 1
        - 'theta_degree': theta_degree
        - 'B_dim': B_dim
        - 'commutation_sign': -1 if eta has odd degree, +1 if even
        - 'is_clifford_type': True if eta^2 != 0 in general
        - 'is_infinite_rank': True iff eta_power_cutoff is None
        - 'eta_power_basis': 'eta^n, n >= 0'
        - 'dim_eta_le_1': 2 * B_dim, the first eta-truncation only
    """
    if B_dim < 1:
        raise ValueError(f"B_dim must be >= 1, got {B_dim}")
    if eta_power_cutoff is not None and eta_power_cutoff < 0:
        raise ValueError(
            f"eta_power_cutoff must be >= 0 when supplied, got {eta_power_cutoff}"
        )

    eta_degree = theta_degree - 1
    # Sign in eta*b = (-1)^{|b|} b*eta depends on eta's degree parity
    # but the relation is graded: eta*b = (-1)^{|b|} b*eta
    commutation_sign = (-1) ** eta_degree
    finite_dim = (
        None if eta_power_cutoff is None else B_dim * (eta_power_cutoff + 1)
    )

    return {
        'dim_B_Theta': finite_dim,
        'eta_degree': eta_degree,
        'theta_degree': theta_degree,
        'B_dim': B_dim,
        'commutation_sign': commutation_sign,
        'is_clifford_type': eta_degree % 2 == 1,  # odd degree => eta^2 can be nonzero
        'is_infinite_rank': eta_power_cutoff is None,
        'eta_power_basis': 'eta^n, n >= 0',
        'eta_power_cutoff': eta_power_cutoff,
        'dim_eta_le_1': 2 * B_dim,
        'convention': 'strict central-shadow Ore extension',
        'forms_twisted_differential': False,
    }


def su2_anomalous_steinberg_profile(k: Any = Symbol('k')) -> Dict[str, Any]:
    """Return the SU(2) anomaly-completed Steinberg package."""
    kappa = simplify(Rational(3, 4) * (k + 2))
    return {
        'boundary_algebra': 'V_k(sl_2)',
        'koszul_dual_level': '-k-4',
        'modular_characteristic': kappa,
        'anomaly_class': 'Theta = kappa omega_1',
        'geometric_source': 'H^3(SU(2);Z)=Z generated by c_3',
        'level_selects': 'k c_3',
        'ore_relation': ore_transgression_relation(True)['relation'],
        'secondary_anomaly': 'u=eta^2',
        'genus_clifford': (
            'G_g(B_Theta)=B_Theta<alpha_i,beta_i>/'
            '(alpha_i alpha_j+alpha_j alpha_i, '
            'beta_i beta_j+beta_j beta_i, '
            'alpha_i beta_j+beta_j alpha_i-delta_ij u)'
        ),
        'invert_u': 'G_g(B_Theta)[u^-1] ~= Mat_{2^g}(B_Theta[u^-1])',
        'u_zero': 'G_g(B_Theta)/(u) ~= (B_Theta/(u)) tensor exterior(alpha_i,beta_i)',
        'string_witness': 'd eta = Theta',
    }


def curved_mc_twist_data() -> Dict[str, Any]:
    """Record the curved Maurer-Cartan flatness convention.

    For a curved dg algebra with d^2 = [Theta, -] and a degree-one
    element eta, the twisted differential d_eta = d + [eta, -] satisfies

        d_eta^2 = [Theta + d eta + 1/2[eta, eta], -].

    In the associative commutator convention, 1/2[eta, eta] = eta^2.
    """
    return {
        'curvature_before_twist': 'Theta',
        'twisted_differential': 'd_eta = d + [eta,-]',
        'twisted_square': 'd_eta^2 = [Theta + d eta + 1/2[eta,eta], -]',
        'curved_mc_equation': 'd eta + 1/2[eta,eta] + Theta = 0',
        'twisted_curvature': 'Theta + d_eta_generator + eta^2',
        'flatness_equation': 'd_eta_generator + eta^2 + Theta = 0',
        'half_bracket_equals_eta_square_for_odd_eta': True,
        'strict_ore_convention': 'central-shadow: d(eta)=Theta, no twist d+[eta,-]',
    }


def virasoro_curved_class_m_completion() -> Dict[str, Any]:
    r"""Return the Virasoro curved-MC anomaly-completed class-M profile."""
    return {
        'base': 'B_ch(Vir_c)^hat_rho',
        'extension': 'k<eta,u>/(u-eta^2)',
        'completed_complex': 'B_ch(Vir_c)^hat_rho tensor k<eta,u>/(u-eta^2)',
        'curved_mc_equation': 'd eta + eta^2 + Theta_Vir = 0',
        'twisted_differential': 'd_eta = d + [eta,-]',
        'secondary_anomaly': 'u = eta^2',
        'secondary_requires_algebra': True,
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
