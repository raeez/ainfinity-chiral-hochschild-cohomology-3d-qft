"""YM Synthesis Engine: boundary BRST, tangent-to-center, and mixed couplings.

Implements the mathematical content of ym_synthesis_core.tex (Chapter 31,
Part VII) concerning boundary deformation theory in the twisted Yang--Mills
framework.

Key mathematical objects:

1. **Tangent-to-center principle** (thm:twisted-ym-tangent-center):
   For Koszul-admissible boundary data, the first-order deformation space
   is canonically identified with the dual of the Koszul dual center:
     T^1_X(B) = HH^2_ch(A_B) ~ Z(A_B^!)^v tensor omega_X.
   In particular, dim HH^2_ch(A_B) = dim Z(A_B^!).

2. **Center-vanishing rigidity** (cor:twisted-ym-center-rigidity):
   If Z(A^!) = 0, then HH^2 = 0, so B is infinitesimally rigid.

3. **One-parameter criterion** (cor:twisted-ym-one-parameter):
   If dim Z(A^!) = 1, then dim HH^2 = 1: a unique first-order coupling.

4. **Central formality** (thm:twisted-ym-central-formality):
   If HH^1(A_B) = 0 (no infinitesimal automorphisms) and
   HH^{-1}(A_B^!) = 0 (no dual negative-degree classes), then the
   formal boundary moduli problem is smooth (an affine space).

5. **Mixed coupling dimension** (thm:twisted-ym-binary-mixed-couplings):
   For n independent Koszul-admissible boundaries, the mixed coupling
   space factorizes as a product of the individual deformation spaces.
   Genuine mixed couplings arise from tensor products of reduced centers.

6. **Kappa complementarity** (Theorem C, Vol I):
   For free fields / Kac-Moody: kappa(A) + kappa(A!) = 0.
   For Virasoro: kappa + kappa' = 13.

References:
  Vol II: ym_synthesis_core.tex (Chapter 31)
  Vol I: chiral_koszul_pairs.tex, concordance.tex (Theorem C/D)
  Costello-Gwilliam (2017): factorization algebras in QFT
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from sympy import Rational, S, Symbol, simplify


# =========================================================================
# 1. KAPPA DATA (from Vol I Theorem D)
# =========================================================================

# Ground truth kappa values from Vol I landscape_census.tex / Theorem D.
# AP1: never copy between families without recomputing.

def _kappa_data(algebra_type: str, **params) -> Dict[str, Any]:
    """Return kappa(A) and kappa(A!) for a standard family.

    The modular characteristic kappa is the leading coefficient of the
    shadow Postnikov tower.  It controls the genus-1 curvature:
      d^2_B = kappa(A) * omega_g.

    Parameters
    ----------
    algebra_type : str
        One of 'heisenberg', 'affine_sl2', 'virasoro', 'betagamma', 'w3'.
    **params
        Family-specific parameters (k for Heisenberg/affine, c for Virasoro/W3).

    Returns
    -------
    dict with keys 'kappa', 'kappa_dual', 'kappa_sum', 'central_charge',
    'dual_central_charge'.
    """
    if algebra_type == 'heisenberg':
        k = params.get('k', Symbol('k'))
        return {
            'kappa': k,
            'kappa_dual': -k,
            'kappa_sum': S.Zero,
            'central_charge': S.One,
            'dual_central_charge': S.One,
        }
    elif algebra_type == 'affine_sl2':
        # sl_2: dim = 3, h^v = 2.
        # kappa(g_k) = dim(g) * (k + h^v) / (2 * h^v)
        k = params.get('k', Symbol('k'))
        h_dual = S(2)
        dim_g = S(3)
        kappa = dim_g * (k + h_dual) / (2 * h_dual)
        # Feigin-Frenkel dual: k' = -k - 2h^v
        k_prime = -k - 2 * h_dual
        kappa_dual = dim_g * (k_prime + h_dual) / (2 * h_dual)
        return {
            'kappa': kappa,
            'kappa_dual': kappa_dual,
            'kappa_sum': simplify(kappa + kappa_dual),
            'central_charge': 3 * k / (k + h_dual),
            'dual_central_charge': 3 * k_prime / (k_prime + h_dual),
        }
    elif algebra_type == 'virasoro':
        c = params.get('c', Symbol('c'))
        # kappa(Vir_c) = c/2.  Dual: Vir_{26-c}, kappa' = (26-c)/2.
        kappa = c / 2
        kappa_dual = (26 - c) / 2
        return {
            'kappa': kappa,
            'kappa_dual': kappa_dual,
            'kappa_sum': simplify(kappa + kappa_dual),
            'central_charge': c,
            'dual_central_charge': 26 - c,
        }
    elif algebra_type == 'betagamma':
        # betagamma: c = -2, kappa = -1.  Dual: kappa' = 1.
        return {
            'kappa': S(-1),
            'kappa_dual': S(1),
            'kappa_sum': S.Zero,
            'central_charge': S(-2),
            'dual_central_charge': S(2),
        }
    elif algebra_type == 'w3':
        c = params.get('c', Symbol('c'))
        # kappa(W_3) = 5c/6.  Dual: W_3 at c' = 100-c, kappa' = 5(100-c)/6.
        # kappa_sum = 5*100/6 = 250/3.
        kappa = 5 * c / 6
        c_dual = 100 - c
        kappa_dual = 5 * c_dual / 6
        return {
            'kappa': kappa,
            'kappa_dual': kappa_dual,
            'kappa_sum': simplify(kappa + kappa_dual),
            'central_charge': c,
            'dual_central_charge': c_dual,
        }
    else:
        raise ValueError(f"Unknown algebra type: {algebra_type}")


# =========================================================================
# 2. TANGENT-TO-CENTER PRINCIPLE
# =========================================================================

def tangent_to_center_dimension(algebra_type: str, **params) -> int:
    """Compute dim HH^2_ch(A_B) = dim Z(A_B^!) for standard families.

    By the tangent-to-center principle (thm:twisted-ym-tangent-center),
    the first-order deformation space of a Koszul-admissible boundary
    datum is canonically identified with the dual center of A^!.

    For the standard families:
    - Heisenberg H_k: Z(H_k^!) = span{level} => dim = 1
    - Affine sl_2 at generic k: Z(g_k^!) = span{level} => dim = 1
    - Virasoro Vir_c: Z(Vir_{26-c}^!) = span{central charge} => dim = 1

    The single deformation parameter is the level (for KM), the level
    (for Heisenberg), or the central charge (for Virasoro).

    Parameters
    ----------
    algebra_type : str
        One of 'heisenberg', 'affine_sl2', 'virasoro'.

    Returns
    -------
    int
        Dimension of the deformation space.
    """
    # For all three standard families, the center of the Koszul dual
    # is one-dimensional, generated by the single parameter (level or c).
    standard_dims = {
        'heisenberg': 1,
        'affine_sl2': 1,
        'virasoro': 1,
    }
    if algebra_type not in standard_dims:
        raise ValueError(
            f"Unknown algebra type '{algebra_type}'. "
            f"Known types: {list(standard_dims.keys())}"
        )
    return standard_dims[algebra_type]


# =========================================================================
# 3. CENTER-VANISHING RIGIDITY
# =========================================================================

def center_vanishing_rigidity(dual_center_dim: int) -> bool:
    """Check if Z(A^!) = 0 implies HH^2 = 0 (infinitesimal rigidity).

    By cor:twisted-ym-center-rigidity: if the center of the Koszul dual
    is zero, then there are no first-order deformations.

    Parameters
    ----------
    dual_center_dim : int
        Dimension of Z(A^!).

    Returns
    -------
    bool
        True if the boundary condition is infinitesimally rigid.
    """
    return dual_center_dim == 0


# =========================================================================
# 4. ONE-PARAMETER CRITERION
# =========================================================================

def one_parameter_criterion(dual_center_dim: int) -> bool:
    """Check if dim Z(A^!) = 1 implies a unique first-order coupling.

    By cor:twisted-ym-one-parameter: if dim Z(A^!) = 1, then
    dim HH^2_ch(A_B) = 1, so there is exactly one coupling up to scale.

    Parameters
    ----------
    dual_center_dim : int
        Dimension of Z(A^!).

    Returns
    -------
    bool
        True if there is exactly one first-order coupling.
    """
    return dual_center_dim == 1


# =========================================================================
# 5. CENTRAL FORMALITY
# =========================================================================

def central_formality_check(hh1: int, hh_minus1_dual: int) -> bool:
    """Check the central formality criterion.

    By thm:twisted-ym-central-formality: the formal boundary moduli
    problem is smooth if HH^1(A_B) = 0 (no infinitesimal automorphisms)
    AND HH^{-1}(A_B^!) = 0 (no dual negative-degree classes).

    When both vanish, the Schlessinger-Stasheff criterion gives a
    smooth formal moduli problem: the deformation space is an affine
    space on HH^2_ch(A_B) ~ Z(A_B^!)^v.

    Parameters
    ----------
    hh1 : int
        Dimension of HH^1_ch(A_B) (infinitesimal automorphisms).
    hh_minus1_dual : int
        Dimension of HH^{-1}(A_B^!) (dual negative-degree classes).

    Returns
    -------
    bool
        True if the formal moduli problem is smooth.
    """
    return hh1 == 0 and hh_minus1_dual == 0


# =========================================================================
# 6. MIXED COUPLING DIMENSION
# =========================================================================

def mixed_coupling_dimension(
    individual_dims: List[int],
    independent: bool = True,
) -> Dict[str, Any]:
    """Compute the mixed coupling space dimension for n boundaries.

    By thm:twisted-ym-binary-mixed-couplings and its n-ary extension:

    For independent boundaries (no mixed OPE between different A_i):
      dim T^1_X(B_1,...,B_n) = sum_i dim(Z_i) + sum_{i<j} dim(Z~_i) * dim(Z~_j)
    where Z~_i is the reduced center (center modulo the unit).

    The first term is the product of individual deformation spaces.
    The second term gives genuine mixed couplings.

    For the standard families (Heisenberg, affine, Virasoro) with
    dim Z = 1: the reduced center Z~ = 0, so there are NO mixed
    couplings. The total dimension is just the sum of individual dims.

    Parameters
    ----------
    individual_dims : list of int
        Dimensions dim Z(A_i^!) for each boundary.
    independent : bool
        Whether the boundary algebras are independent (no mixed OPE).

    Returns
    -------
    dict with keys:
        'total_dim': total dimension of the coupling space
        'individual_sum': sum of individual dims
        'mixed_dim': dimension of genuinely mixed couplings
        'n_boundaries': number of boundaries
    """
    n = len(individual_dims)
    individual_sum = sum(individual_dims)

    if not independent:
        # Non-independent case: we cannot decompose further without
        # knowing the mixed OPE structure.
        return {
            'total_dim': None,
            'individual_sum': individual_sum,
            'mixed_dim': None,
            'n_boundaries': n,
        }

    # For independent boundaries, the reduced center is Z~_i = Z_i / (unit).
    # If dim Z_i = 1, then Z~_i = 0 (the single generator IS the unit parameter).
    # Mixed couplings = sum_{i<j} dim(Z~_i) * dim(Z~_j).
    # For standard families with dim = 1: Z~ = 0, no mixed couplings.
    reduced_dims = [max(0, d - 1) for d in individual_dims]
    mixed_dim = 0
    for i in range(n):
        for j in range(i + 1, n):
            mixed_dim += reduced_dims[i] * reduced_dims[j]

    total_dim = individual_sum + mixed_dim
    return {
        'total_dim': total_dim,
        'individual_sum': individual_sum,
        'mixed_dim': mixed_dim,
        'n_boundaries': n,
    }


# =========================================================================
# 7. COMPLEMENTARITY CHECK (cross-engine consistency with Vol I)
# =========================================================================

def complementarity_check(algebra_type: str, **params) -> Dict[str, Any]:
    """Verify kappa + kappa' satisfies the expected complementarity sum.

    From Vol I Theorem C:
    - Free fields (Heisenberg, betagamma) and Kac-Moody: kappa + kappa' = 0.
    - Virasoro: kappa + kappa' = 13.
    - W_3: kappa + kappa' = 250/3.

    Parameters
    ----------
    algebra_type : str
        The algebra family.
    **params
        Family-specific parameters.

    Returns
    -------
    dict with keys 'kappa', 'kappa_dual', 'kappa_sum', 'expected_sum', 'match'.
    """
    data = _kappa_data(algebra_type, **params)

    expected_sums = {
        'heisenberg': S.Zero,
        'affine_sl2': S.Zero,
        'virasoro': S(13),
        'betagamma': S.Zero,
        'w3': Rational(250, 3),
    }
    expected = expected_sums[algebra_type]

    kappa_sum = simplify(data['kappa'] + data['kappa_dual'])
    match = simplify(kappa_sum - expected) == 0

    return {
        'kappa': data['kappa'],
        'kappa_dual': data['kappa_dual'],
        'kappa_sum': kappa_sum,
        'expected_sum': expected,
        'match': match,
    }


# =========================================================================
# 8. BOUNDARY BRST DATA
# =========================================================================

@dataclass
class BoundaryBRSTData:
    """Data for a boundary BRST package.

    Packages the three columns of the boundary complex:
      Open_B = B(A_B)  (the bar complex = open boundary)
      Bulk^v_B = Omega(A_B^!)  (the cobar of the dual = dual bulk)
      T^1_X(B) = HH^2_ch(A_B)  (the tangent space)

    The tangent-to-center principle identifies T^1 ~ Z(A^!)^v.
    """
    algebra_type: str
    open_name: str        # B(A_B)
    bulk_dual_name: str   # Omega(A_B^!)
    tangent_dim: int      # dim HH^2_ch = dim Z(A^!)
    kappa: Any
    kappa_dual: Any
    is_rigid: bool        # Z(A^!) = 0 => rigid
    is_one_parameter: bool  # dim Z(A^!) = 1
    is_formal: bool       # formality criterion satisfied


def boundary_brst_data(algebra_type: str, **params) -> BoundaryBRSTData:
    """Construct the full boundary BRST data for a standard family.

    Parameters
    ----------
    algebra_type : str
        One of 'heisenberg', 'affine_sl2', 'virasoro'.
    **params
        Family-specific parameters.

    Returns
    -------
    BoundaryBRSTData
    """
    dim = tangent_to_center_dimension(algebra_type, **params)
    kdata = _kappa_data(algebra_type, **params)

    names = {
        'heisenberg': ('B(H_k)', 'Omega(Sym^ch(V*))'),
        'affine_sl2': ('B(g_k)', 'Omega(g_{k\'}^!)'),
        'virasoro': ('B(Vir_c)', 'Omega(Vir_{26-c})'),
    }
    open_n, bulk_n = names.get(algebra_type, ('B(A)', 'Omega(A^!)'))

    # For standard families: HH^1 = 0 (no automorphisms beyond overall scale),
    # HH^{-1}(A^!) = 0 (no negative-degree Hochschild classes).
    # Hence central formality holds.
    return BoundaryBRSTData(
        algebra_type=algebra_type,
        open_name=open_n,
        bulk_dual_name=bulk_n,
        tangent_dim=dim,
        kappa=kdata['kappa'],
        kappa_dual=kdata['kappa_dual'],
        is_rigid=center_vanishing_rigidity(dim),
        is_one_parameter=one_parameter_criterion(dim),
        is_formal=True,  # all standard families satisfy formality
    )
