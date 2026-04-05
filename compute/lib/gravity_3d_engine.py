"""3D Gravity Compute Engine: Virasoro-Koszul duality and genus expansion.

Implements the gravitational sector of the holomorphic-topological programme.

The Virasoro algebra Vir_c governs 3d gravity via the holomorphic-topological
correspondence: the boundary chiral algebra is Vir_c, the bulk is the
3d gravitational theory on C x R, and the Koszul dual Vir_{26-c} controls
the dual gravitational algebra (c_dual = 26 - c, self-dual at c = 13).

Key mathematical objects:

1. **Virasoro lambda-bracket**: {T_lam T} = (c/12) lam^3 + 2T lam + dT.
   This encodes the OPE T(z)T(w) ~ (c/2)/(z-w)^4 + 2T(w)/(z-w)^2 + dT(w)/(z-w).

2. **Ternary associator**: A_3(T,T,T; lam12, lam23) from the failure of
   associativity of the lambda-bracket. m_3 = -A_3 is the ternary A-infinity
   operation from homotopy transfer.

3. **Quartic contact invariant**: Q^contact_Vir = 10/(c(5c+22)), the first
   nonlinear modular shadow coefficient beyond kappa.

4. **Gravitational kappa**: kappa(Vir_c) = c/2 (Theorem D, modular
   characteristic for 3d gravity).

5. **Genus generating function**: F_g = kappa_eff * lambda_g^FP via the
   A-hat genus formula, where kappa_eff = (c-26)/2 accounts for the
   Koszul dual shift.

6. **R-matrix pole structure**: The classical r-matrix r(z) has order-4
   leading pole with residue c/2, reflecting the quartic pole in the
   Virasoro OPE.

References:
  Vol II: 3d_gravity.tex (Movement VI)
  Vol I: concordance.tex (Theorem D), higher_genus_modular_koszul.tex
  Vol I: nonlinear_modular_shadows.tex (quartic contact)
  Vol I: genus_generating_function.py (A-hat coefficients)
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any, Dict, List, Optional, Tuple

from sympy import (
    Symbol, Rational, simplify, expand, S, symbols, factorial,
    bernoulli, sinh, series, Abs, collect, Poly,
)


# =========================================================================
# 1. VIRASORO LAMBDA-BRACKET
# =========================================================================

def virasoro_lambda_bracket(c=None):
    """Return the Virasoro lambda-bracket {T_lam T}.

    {T_lam T} = (c/12) lam^3 + 2T lam + dT

    This encodes the OPE:
      T(z)T(w) ~ (c/2)/(z-w)^4 + 2T(w)/(z-w)^2 + dT(w)/(z-w)

    via the state-field dictionary: T_{(n)} T corresponds to the
    (n+1)-th pole in the OPE, and {T_lam T} = sum_n T_{(n)}T lam^n/n!.

    Pole dictionary:
      (z-w)^{-4}: c/2    -> T_{(3)}T = c/2   -> coeff of lam^3/3! = c/12
      (z-w)^{-2}: 2T     -> T_{(1)}T = 2T    -> coeff of lam^1/1! = 2T
      (z-w)^{-1}: dT     -> T_{(0)}T = dT    -> coeff of lam^0/0! = dT

    Returns a dict with keys 'lam3', 'lam1_T', 'lam0_dT' for the
    coefficients of lam^3, T*lam, and dT respectively.

    Parameters
    ----------
    c : optional
        Central charge. If None, uses symbolic Symbol('c').
    """
    c_val = Symbol('c') if c is None else S(c)

    return {
        'lam3': c_val / 12,        # coefficient of lam^3
        'lam1_T': S(2),            # coefficient of T * lam
        'lam0_dT': S(1),           # coefficient of dT (= partial T)
        'central_charge': c_val,
    }


# =========================================================================
# 2. VIRASORO ASSOCIATOR A_3
# =========================================================================

def virasoro_associator(c=None, lam12=None, lam23=None):
    r"""Compute the ternary associator A_3(T, T, T; lam12, lam23).

    The associator measures the failure of the lambda-bracket to be
    associative (it satisfies the Jacobi identity instead). For the
    Virasoro algebra with a single generator T of conformal weight 2:

    A_3(T, T, T; lam12, lam23) =
        -d^2 T
        - (2 lam12 + 3 lam23) dT
        - 2 lam23 (2 lam12 + lam23) T
        - (c/12) lam23^3 (2 lam12 + lam23)

    This is computed from the Jacobi identity residue:
        {T_{lam12} {T_{lam23} T}} - {T_{lam23} {T_{lam12} T}}
        - {{T_{lam12} T}_{lam12+lam23} T}

    Returns a dict of coefficients for {d^2T, dT, T, scalar}
    as polynomials in (lam12, lam23, c).
    """
    c_val = Symbol('c') if c is None else S(c)
    l12 = Symbol('lambda_12') if lam12 is None else S(lam12)
    l23 = Symbol('lambda_23') if lam23 is None else S(lam23)

    # Coefficient of d^2 T
    coeff_d2T = S(-1)

    # Coefficient of dT
    coeff_dT = -(2 * l12 + 3 * l23)

    # Coefficient of T
    coeff_T = -2 * l23 * (2 * l12 + l23)

    # Scalar term (no field, pure central extension contribution)
    coeff_scalar = -c_val / 12 * l23**3 * (2 * l12 + l23)

    return {
        'd2T': coeff_d2T,
        'dT': expand(coeff_dT),
        'T': expand(coeff_T),
        'scalar': expand(coeff_scalar),
        'c': c_val,
        'lam12': l12,
        'lam23': l23,
    }


# =========================================================================
# 3. m_3 = -A_3
# =========================================================================

def virasoro_m3_coefficients(c=None, lam12=None, lam23=None):
    r"""Compute m_3(T, T, T; lam12, lam23) = -A_3(T, T, T; lam12, lam23).

    The ternary A-infinity operation m_3 is the negation of the associator.
    This follows from the general homotopy transfer formula: the A-infinity
    operations are defined such that the failure of associativity at order n
    is controlled by m_{n+1}.

    Returns a dict with the same keys as virasoro_associator, but negated.
    """
    A3 = virasoro_associator(c=c, lam12=lam12, lam23=lam23)

    return {
        'd2T': -A3['d2T'],
        'dT': expand(-A3['dT']),
        'T': expand(-A3['T']),
        'scalar': expand(-A3['scalar']),
        'c': A3['c'],
        'lam12': A3['lam12'],
        'lam23': A3['lam23'],
    }


# =========================================================================
# 4. QUARTIC CONTACT INVARIANT
# =========================================================================

def quartic_contact_virasoro(c=None):
    r"""Quartic contact invariant Q^contact_Vir = 10/(c(5c+22)).

    This is the first nonlinear modular shadow coefficient beyond kappa.
    It measures the quartic resonance obstruction in the shadow Postnikov
    tower for the Virasoro algebra.

    Poles: c = 0 and c = -22/5.
    Positive for c > 0.
    The quartic contact invariant controls the genus-2 shadow correction.

    Ground truth: Vol I, thm:nms-virasoro-quartic-explicit.
    """
    c_val = Symbol('c') if c is None else S(c)
    return 10 / (c_val * (5 * c_val + 22))


def quartic_contact_virasoro_exact(c_num, c_den=1):
    """Exact rational quartic contact invariant.

    Uses fractions.Fraction for exact arithmetic.
    """
    c_val = Fraction(c_num, c_den)
    return Fraction(10, 1) / (c_val * (5 * c_val + 22))


# =========================================================================
# 5. GRAVITATIONAL KAPPA
# =========================================================================

def gravity_kappa(c=None):
    r"""Modular characteristic kappa(Vir_c) = c/2.

    This is the genus-1 curvature: d_B^2 = kappa(A) * omega_1.
    For the Virasoro algebra, kappa = c/2.

    Ground truth: Vol I Theorem D, landscape_census.tex.
    """
    c_val = Symbol('c') if c is None else S(c)
    return c_val / 2


def gravity_kappa_exact(c_num, c_den=1):
    """Exact rational kappa using fractions.Fraction."""
    return Fraction(c_num, c_den) / 2


# =========================================================================
# 6. KOSZUL DUAL CENTRAL CHARGE
# =========================================================================

def koszul_dual_central_charge(c=None):
    r"""Koszul dual central charge: c_dual = 26 - c.

    Vir_c^! = Vir_{26-c} as chiral algebras.

    CRITICAL: Self-dual at c = 13, NOT c = 26.
    This is the Feigin-Frenkel involution for Virasoro.
    """
    c_val = Symbol('c') if c is None else S(c)
    return 26 - c_val


# =========================================================================
# 7. R-MATRIX POLE STRUCTURE
# =========================================================================

def gravity_r_matrix_poles(c=None):
    r"""Classical r-matrix pole structure for the Virasoro algebra.

    The r-matrix r(z) = Res^coll_{0,2}(Theta_A) from the binary
    genus-0 shadow of the universal MC element.

    For Virasoro, the r-matrix has the pole structure:
      r(z) = (c/2)/z^4 + 2T/z^2 + dT/z

    This directly mirrors the OPE pole structure of T(z)T(w).

    Returns a dict {pole_order: residue_description}.
    """
    c_val = Symbol('c') if c is None else S(c)

    return {
        4: {'residue': c_val / 2, 'description': 'c/2 (quartic Casimir)'},
        2: {'residue': '2T', 'description': '2T (energy-momentum)'},
        1: {'residue': 'dT', 'description': 'dT (derivative)'},
    }


def gravity_r_matrix_leading_residue(c=None):
    """Leading residue of the r-matrix: c/2 at the order-4 pole."""
    c_val = Symbol('c') if c is None else S(c)
    return c_val / 2


# =========================================================================
# 8. COMPLEMENTARITY CONSTANT
# =========================================================================

def complementarity_constant_virasoro():
    r"""Complementarity constant for Virasoro: kappa + kappa' = 13.

    kappa(Vir_c) + kappa(Vir_{26-c}) = c/2 + (26-c)/2 = 13.

    This is the Virasoro instance of Theorem C (complementarity).
    The constant 13 = kappa(Vir_13) is the self-dual value.

    For comparison:
      - Heisenberg: kappa + kappa' = 0
      - Affine KM: kappa + kappa' = 0
      - W_3: kappa + kappa' = 250/3
    """
    return S(13)


def verify_complementarity(c=None):
    """Verify kappa(Vir_c) + kappa(Vir_{26-c}) = 13."""
    c_val = Symbol('c') if c is None else S(c)
    kappa = gravity_kappa(c_val)
    kappa_dual = gravity_kappa(26 - c_val)
    total = simplify(kappa + kappa_dual)
    return {
        'kappa': kappa,
        'kappa_dual': kappa_dual,
        'sum': total,
        'equals_13': simplify(total - 13) == 0,
    }


# =========================================================================
# 9. GENUS GENERATING FUNCTION (A-HAT GENUS)
# =========================================================================

def _lambda_fp(g):
    r"""Faber-Pandharipande intersection number lambda_g^FP.

    lambda_g^FP = (2^{2g-1} - 1) / 2^{2g-1} * |B_{2g}| / (2g)!

    where B_{2g} is the 2g-th Bernoulli number.

    These are the coefficients in:
      A-hat(x) = (x/2)/sinh(x/2) = sum_{g>=0} (-1)^g lambda_g^FP x^{2g}
    with lambda_0^FP := 1.
    """
    if g < 1:
        raise ValueError(f"Genus must be >= 1, got {g}")
    B_2g = bernoulli(2 * g)
    numerator = (2**(2*g - 1) - 1) * Abs(B_2g)
    denominator = 2**(2*g - 1) * factorial(2 * g)
    return Rational(numerator, denominator)


def _ahat_coefficient(g):
    r"""Coefficient of x^{2g} in A-hat(x) = (x/2)/sinh(x/2).

    Returns (-1)^g * lambda_g^FP for g >= 1, and 1 for g = 0.
    """
    if g == 0:
        return Rational(1)
    return (-1)**g * _lambda_fp(g)


def genus_generating_function_coefficients(c=None, max_genus=5):
    r"""Genus expansion free energy coefficients F_g for Virasoro.

    The genus-g free energy from the A-hat genus formula:
      F_g = kappa_eff * lambda_g^FP

    where kappa_eff = (c - 26)/2 is the EFFECTIVE curvature that
    accounts for the Koszul dual shift. The physical free energy
    uses kappa_eff rather than kappa = c/2 because the partition
    function involves the relative curvature.

    IMPORTANT: The F_g are POSITIVE when c < 26 (Bernoulli sign
    analysis: lambda_g^FP > 0 always, and kappa_eff < 0 for c < 26,
    but the physical free energy has an additional (-1)^g from
    the A-hat series, so F_g = (-1)^{g+1} * kappa_eff * |lambda_g^FP|
    ... actually let's be precise:

    The scalar free energy at genus g is:
      F_g = kappa * lambda_g^FP
    where kappa = c/2 and lambda_g^FP > 0.

    So F_g > 0 iff kappa > 0 iff c > 0.

    The EFFECTIVE (shifted) version uses kappa_eff = (c-26)/2.

    We return BOTH: the raw F_g = kappa * lambda_g^FP and the
    shifted F_g^eff = kappa_eff * lambda_g^FP.

    Explicit values:
      lambda_1 = 1/24
      lambda_2 = 7/5760
      lambda_3 = 31/967680
      lambda_4 = 127/154828800

    Returns
    -------
    dict with keys 'raw' (using kappa=c/2) and 'effective' (using
    kappa_eff=(c-26)/2), each mapping genus g to F_g.
    """
    c_val = Symbol('c') if c is None else S(c)
    kappa = c_val / 2
    kappa_eff = (c_val - 26) / 2

    raw = {}
    effective = {}
    lambda_fp_values = {}

    for g in range(1, max_genus + 1):
        lfp = _lambda_fp(g)
        lambda_fp_values[g] = lfp
        raw[g] = expand(kappa * lfp)
        effective[g] = expand(kappa_eff * lfp)

    return {
        'raw': raw,
        'effective': effective,
        'lambda_fp': lambda_fp_values,
        'kappa': kappa,
        'kappa_eff': kappa_eff,
    }


def ahat_series_coefficients(max_genus=10):
    r"""Return coefficients of x^{2g} in A-hat(x) = (x/2)/sinh(x/2).

    Standard expansion:
      A-hat(x) = 1 - x^2/24 + 7x^4/5760 - 31x^6/967680 + ...

    The sign pattern is (-1)^g for the genus-g coefficient.
    """
    result = {0: Rational(1)}
    for g in range(1, max_genus + 1):
        result[g] = _ahat_coefficient(g)
    return result


def verify_ahat_series(max_genus=6):
    r"""Verify A-hat coefficients against direct series expansion of (x/2)/sinh(x/2)."""
    x = Symbol('x')
    s = series(x / 2 / sinh(x / 2), x, 0, 2 * max_genus + 2)

    results = {}
    all_match = True

    for g in range(max_genus + 1):
        coeff_from_series = Rational(s.coeff(x, 2 * g))
        coeff_from_formula = _ahat_coefficient(g)
        match = simplify(coeff_from_series - coeff_from_formula) == 0
        all_match = all_match and match
        results[g] = {
            'series': coeff_from_series,
            'formula': coeff_from_formula,
            'match': match,
        }

    return {
        'all_match': all_match,
        'by_genus': results,
    }


# =========================================================================
# 10. SHADOW DEPTH AND CLASSIFICATION
# =========================================================================

def virasoro_shadow_class():
    """Shadow depth class for Virasoro: M (mixed, r_max = infinity).

    The Virasoro algebra has infinite shadow depth because the quintic
    obstruction o^(5)_Vir != 0 (the shadow obstruction tower never terminates).

    Contrast with:
      Heisenberg: G (Gaussian, r_max = 2)
      Affine KM:  L (Lie/tree, r_max = 3)
      Beta-gamma: C (contact/quartic, r_max = 4)
    """
    return 'M'


def virasoro_shadow_depth():
    """Shadow depth r_max for Virasoro: infinity."""
    return float('inf')


# =========================================================================
# 11. GENUS-1 HESSIAN CORRECTION
# =========================================================================

def genus1_hessian_correction_virasoro(c=None):
    r"""Genus-1 Hessian correction delta_H^(1)_Vir.

    delta_H^(1)_Vir = 120 / [c^2 (5c + 22)] * x^2

    This is the genus-1 correction to the shadow obstruction tower from the
    quartic contact invariant. It lives in the weight-4 sector.

    Returns the coefficient (without the x^2 factor).
    """
    c_val = Symbol('c') if c is None else S(c)
    return 120 / (c_val**2 * (5 * c_val + 22))


# =========================================================================
# 12. SUMMARY / DIAGNOSTIC
# =========================================================================

def gravity_diagnostic(c_val):
    """Comprehensive diagnostic for 3d gravity at central charge c.

    Returns a dict summarizing all gravitational invariants.
    """
    c = S(c_val)
    c_dual = koszul_dual_central_charge(c)
    kappa = gravity_kappa(c)
    kappa_dual = gravity_kappa(c_dual)
    q_contact = quartic_contact_virasoro(c)
    genus_data = genus_generating_function_coefficients(c, max_genus=4)

    return {
        'c': c,
        'c_dual': c_dual,
        'self_dual': simplify(c - 13) == 0,
        'kappa': kappa,
        'kappa_dual': kappa_dual,
        'kappa_sum': simplify(kappa + kappa_dual),
        'Q_contact': simplify(q_contact),
        'shadow_class': virasoro_shadow_class(),
        'shadow_depth': virasoro_shadow_depth(),
        'F_g_raw': {g: simplify(v) for g, v in genus_data['raw'].items()},
        'delta_H1': simplify(genus1_hessian_correction_virasoro(c)),
    }
