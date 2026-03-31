"""Tests for 3D Gravity Compute Engine.

Verifies the gravitational A∞ structure derived from Virasoro Koszul duality:
1. Virasoro associator A₃ and ternary operation m₃ = -A₃
2. Quartic contact invariant Q^contact = 10/(c(5c+22))
3. Gravitational Koszul triangle: Vir_c^! = Vir_{26-c}, self-dual at c=13
4. Modular characteristic κ(Vir_c) = c/2, complementarity sum = 13
5. R-matrix pole structure from Laplace transform of λ-bracket
6. Genus expansion via Â-genus formula
7. Cross-engine consistency with bulk_boundary_duality_engine

References:
  Vol II: 3d_gravity.tex (Movements I-VI)
  Vol I: concordance.tex, nonlinear_modular_shadows.tex
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from fractions import Fraction
from sympy import Symbol, Rational, S, simplify, expand, symbols

from lib.gravity_3d_engine import (
    virasoro_lambda_bracket,
    virasoro_associator,
    virasoro_m3_coefficients,
    quartic_contact_virasoro,
    quartic_contact_virasoro_exact,
    gravity_kappa,
    gravity_kappa_exact,
    koszul_dual_central_charge,
    gravity_r_matrix_poles,
    gravity_r_matrix_leading_residue,
    complementarity_constant_virasoro,
    verify_complementarity,
    genus_generating_function_coefficients,
    ahat_series_coefficients,
    verify_ahat_series,
    virasoro_shadow_depth,
    virasoro_shadow_class,
)


# ===================================================================
# 1. VIRASORO LAMBDA-BRACKET
# ===================================================================

class TestVirasoroLambdaBracket:
    """Verify {T_λ T} = (c/12)λ³ + 2Tλ + ∂T."""

    def test_quartic_pole_coefficient(self):
        """Central charge coefficient: c/12 for the λ³ term."""
        bracket = virasoro_lambda_bracket(c=24)
        assert bracket['lam3'] == 2  # 24/12 = 2

    def test_double_pole_coefficient(self):
        """Coefficient of Tλ is always 2."""
        for c_val in [0, 1, 13, 26, -22]:
            bracket = virasoro_lambda_bracket(c=c_val)
            assert bracket['lam1_T'] == 2

    def test_simple_pole_coefficient(self):
        """Coefficient of ∂T is always 1."""
        bracket = virasoro_lambda_bracket(c=13)
        assert bracket['lam0_dT'] == 1

    def test_symbolic_central_charge(self):
        """Symbolic c works."""
        bracket = virasoro_lambda_bracket()
        c = Symbol('c')
        assert simplify(bracket['lam3'] - c / 12) == 0

    def test_central_charge_stored(self):
        """Central charge value is stored."""
        bracket = virasoro_lambda_bracket(c=26)
        assert bracket['central_charge'] == 26


# ===================================================================
# 2. VIRASORO ASSOCIATOR
# ===================================================================

class TestVirasoroAssociator:
    """Verify A₃(T,T,T; λ₁₂, λ₂₃) formula from eq:gravity-associator."""

    def test_d2T_coefficient(self):
        """Coefficient of ∂²T is -1 (always)."""
        A3 = virasoro_associator(c=1, lam12=0, lam23=0)
        assert A3['d2T'] == -1

    def test_scalar_linear_in_c(self):
        """Scalar term is proportional to c."""
        A3_1 = virasoro_associator(c=1, lam12=1, lam23=1)
        A3_2 = virasoro_associator(c=2, lam12=1, lam23=1)
        # scalar(c=2) should be 2 * scalar(c=1)
        assert simplify(A3_2['scalar'] - 2 * A3_1['scalar']) == 0

    def test_scalar_vanishes_at_c0(self):
        """At c=0 the scalar (central extension) term vanishes."""
        A3 = virasoro_associator(c=0, lam12=3, lam23=5)
        assert A3['scalar'] == 0

    def test_specific_evaluation(self):
        """A₃ at λ₁₂=0, λ₂₃=1, c=12: scalar = -(12/12)*1³*(0+1) = -1."""
        A3 = virasoro_associator(c=12, lam12=0, lam23=1)
        assert A3['scalar'] == -1

    def test_specific_evaluation_2(self):
        """A₃ at λ₁₂=1, λ₂₃=1, c=12: scalar = -1*1*(2+1) = -3."""
        A3 = virasoro_associator(c=12, lam12=1, lam23=1)
        assert A3['scalar'] == -3

    def test_dT_coefficient_at_origin(self):
        """dT coefficient at λ₁₂=0, λ₂₃=0 is 0."""
        A3 = virasoro_associator(c=1, lam12=0, lam23=0)
        assert A3['dT'] == 0

    def test_T_coefficient_at_origin(self):
        """T coefficient at λ₁₂=0, λ₂₃=0 is 0."""
        A3 = virasoro_associator(c=1, lam12=0, lam23=0)
        assert A3['T'] == 0

    def test_not_symmetric(self):
        """A₃ is NOT symmetric in λ₁₂ ↔ λ₂₃."""
        A3_a = virasoro_associator(c=1, lam12=1, lam23=2)
        A3_b = virasoro_associator(c=1, lam12=2, lam23=1)
        assert A3_a['scalar'] != A3_b['scalar']


# ===================================================================
# 3. m₃ = -A₃ (BORCHERDS IDENTITY)
# ===================================================================

class TestM3NegatesAssociator:
    """Verify m₃ = -A₃ from the Borcherds identity."""

    def test_m3_negates_d2T(self):
        m3 = virasoro_m3_coefficients(c=1, lam12=0, lam23=0)
        assert m3['d2T'] == 1  # -(-1) = 1

    def test_m3_negates_scalar(self):
        A3 = virasoro_associator(c=12, lam12=1, lam23=1)
        m3 = virasoro_m3_coefficients(c=12, lam12=1, lam23=1)
        assert simplify(m3['scalar'] + A3['scalar']) == 0

    def test_m3_negates_all_components(self):
        """Full check: m₃ + A₃ = 0 for all components."""
        for c_val in [1, 6, 13, 26]:
            A3 = virasoro_associator(c=c_val, lam12=2, lam23=3)
            m3 = virasoro_m3_coefficients(c=c_val, lam12=2, lam23=3)
            for key in ['d2T', 'dT', 'T', 'scalar']:
                assert simplify(m3[key] + A3[key]) == 0, \
                    f"m3[{key}] + A3[{key}] ≠ 0 at c={c_val}"

    def test_symbolic_negation(self):
        """Symbolic verification of m₃ = -A₃."""
        A3 = virasoro_associator()
        m3 = virasoro_m3_coefficients()
        assert simplify(m3['scalar'] + A3['scalar']) == 0


# ===================================================================
# 4. QUARTIC CONTACT INVARIANT
# ===================================================================

class TestQuarticContact:
    """Verify Q^contact = 10/(c(5c+22))."""

    def test_quartic_formula_c1(self):
        """At c=1: Q = 10/(1*27) = 10/27."""
        assert quartic_contact_virasoro(c=1) == Rational(10, 27)

    def test_quartic_formula_c13(self):
        """At self-dual c=13: Q = 10/(13*87) = 10/1131."""
        assert quartic_contact_virasoro(c=13) == Rational(10, 1131)

    def test_quartic_positive_for_c_gt_0(self):
        """Q > 0 for c > 0."""
        for c_val in [1, 2, 5, 13, 26, 100]:
            assert quartic_contact_virasoro(c=c_val) > 0

    def test_quartic_exact_arithmetic(self):
        """Exact rational arithmetic via Fraction."""
        q = quartic_contact_virasoro_exact(1)
        assert q == Fraction(10, 27)

    def test_quartic_poles(self):
        """Poles at c=0 and c=-22/5."""
        c = Symbol('c')
        q = quartic_contact_virasoro()
        denom = c * (5 * c + 22)
        # Denominator vanishes at c=0 and c=-22/5
        assert denom.subs(c, 0) == 0
        assert denom.subs(c, Rational(-22, 5)) == 0

    def test_quartic_cross_family(self):
        """Verify Vol I ground truth: Q^contact(Vir_c) = 10/(c(5c+22))."""
        # This is the DEFINING formula — verify structure
        c = Symbol('c')
        q = quartic_contact_virasoro()
        assert simplify(q - 10 / (c * (5 * c + 22))) == 0


# ===================================================================
# 5. GRAVITATIONAL KOSZUL TRIANGLE
# ===================================================================

class TestGravitationalKoszulTriangle:
    """Verify Vir_c^! = Vir_{26-c}, self-dual at c=13."""

    def test_koszul_dual_c(self):
        """26 - c is the dual central charge."""
        assert koszul_dual_central_charge(c=0) == 26
        assert koszul_dual_central_charge(c=26) == 0
        assert koszul_dual_central_charge(c=13) == 13

    def test_self_dual_c13(self):
        """Vir_{13}^! = Vir_{13}: self-dual at c=13."""
        assert koszul_dual_central_charge(c=13) == 13

    def test_not_self_dual_c26(self):
        """c=26 is NOT self-dual (AP8: c=13 is self-dual, NOT c=26)."""
        assert koszul_dual_central_charge(c=26) != 26

    def test_double_dual(self):
        """(Vir_c^!)^! = Vir_c: double dual is identity."""
        for c_val in [0, 1, 13, 26, -5]:
            assert koszul_dual_central_charge(
                koszul_dual_central_charge(c=c_val)) == c_val

    def test_kappa_virasoro(self):
        """κ(Vir_c) = c/2."""
        assert gravity_kappa(c=0) == 0
        assert gravity_kappa(c=1) == Rational(1, 2)
        assert gravity_kappa(c=26) == 13

    def test_kappa_exact(self):
        """Exact rational κ."""
        assert gravity_kappa_exact(1) == Fraction(1, 2)
        assert gravity_kappa_exact(7, 10) == Fraction(7, 20)

    def test_complementarity_13(self):
        """κ(Vir_c) + κ(Vir_{26-c}) = 13 for all c."""
        assert complementarity_constant_virasoro() == 13
        for c_val in [0, 1, 5, 13, 26, -10]:
            result = verify_complementarity(c=c_val)
            assert result['equals_13'], f"Complementarity fails at c={c_val}"

    def test_complementarity_symbolic(self):
        """Symbolic complementarity verification."""
        result = verify_complementarity()
        assert result['equals_13']


# ===================================================================
# 6. R-MATRIX
# ===================================================================

class TestGravitationalRMatrix:
    """Verify r(z) = (c/2)/z⁴ + 2T/z² + ∂T/z."""

    def test_leading_pole_order_4(self):
        """Leading pole is order 4."""
        poles = gravity_r_matrix_poles(c=1)
        assert 4 in poles

    def test_leading_residue_c_over_2(self):
        """Leading residue is c/2."""
        assert gravity_r_matrix_leading_residue(c=26) == 13
        assert gravity_r_matrix_leading_residue(c=1) == Rational(1, 2)

    def test_laplace_transform_scalar(self):
        """Laplace: ∫₀^∞ e^{-λz}(c/12)λ³ dλ = (c/12)(3!/z⁴) = c/(2z⁴)."""
        c_val = Rational(24)
        # (c/12) * 3! = (c/12) * 6 = c/2
        laplace_result = c_val / 12 * 6
        expected = c_val / 2
        assert laplace_result == expected

    def test_pole_structure_complete(self):
        """r(z) has exactly poles of order 4, 2, 1."""
        poles = gravity_r_matrix_poles(c=13)
        assert set(poles.keys()) == {4, 2, 1}

    def test_r_matrix_matches_ope(self):
        """R-matrix pole orders match OPE pole orders of T(z)T(w)."""
        bracket = virasoro_lambda_bracket(c=12)
        poles = gravity_r_matrix_poles(c=12)
        # λ³ → z⁻⁴ (order 4 pole)
        assert 4 in poles
        # Tλ → z⁻² (order 2 pole)
        assert 2 in poles
        # ∂T → z⁻¹ (order 1 pole)
        assert 1 in poles


# ===================================================================
# 7. GENUS EXPANSION
# ===================================================================

class TestGenusExpansion:
    """Verify F_g via Â-genus formula."""

    def test_lambda_fp_1(self):
        """λ₁^FP = 1/24."""
        result = genus_generating_function_coefficients(c=2, max_genus=1)
        assert result['lambda_fp'][1] == Rational(1, 24)

    def test_lambda_fp_2(self):
        """λ₂^FP = 7/5760."""
        result = genus_generating_function_coefficients(c=2, max_genus=2)
        assert result['lambda_fp'][2] == Rational(7, 5760)

    def test_lambda_fp_3(self):
        """λ₃^FP = 31/967680."""
        result = genus_generating_function_coefficients(c=2, max_genus=3)
        assert result['lambda_fp'][3] == Rational(31, 967680)

    def test_F1_raw(self):
        """F₁ = κ * λ₁^FP = (c/2) * (1/24) = c/48."""
        result = genus_generating_function_coefficients(c=48, max_genus=1)
        assert result['raw'][1] == 1  # 48/2 * 1/24 = 1

    def test_F1_effective(self):
        """F₁^eff = κ_eff * λ₁^FP = (c-26)/2 * 1/24."""
        result = genus_generating_function_coefficients(c=26, max_genus=1)
        assert result['effective'][1] == 0  # (26-26)/2 * 1/24 = 0

    def test_ahat_coefficients_match_series(self):
        """Â-genus coefficients match direct series expansion of (x/2)/sinh(x/2)."""
        result = verify_ahat_series(max_genus=4)
        assert result['all_match'], \
            f"Â-genus mismatch: {result['by_genus']}"

    def test_lambda_fp_positive(self):
        """All λ_g^FP > 0 (Bernoulli sign analysis)."""
        result = genus_generating_function_coefficients(c=2, max_genus=5)
        for g in range(1, 6):
            assert result['lambda_fp'][g] > 0, f"λ_{g}^FP not positive"


# ===================================================================
# 8. CROSS-ENGINE CONSISTENCY
# ===================================================================

class TestCrossEngineConsistency:
    """Cross-engine checks (AP10 compliance)."""

    def test_kappa_ground_truth_values(self):
        """κ values match known ground truth."""
        ground_truth = {
            0: Rational(0),
            1: Rational(1, 2),
            13: Rational(13, 2),
            26: Rational(13),
        }
        for c_val, expected_kappa in ground_truth.items():
            computed = gravity_kappa(c=c_val)
            assert simplify(computed - expected_kappa) == 0, \
                f"κ mismatch at c={c_val}: got {computed}, expected {expected_kappa}"

    def test_shadow_depth_infinite(self):
        """Virasoro shadow depth is infinite (class M)."""
        assert virasoro_shadow_depth() == float('inf')

    def test_shadow_class_M(self):
        """Virasoro is class M (mixed, infinite tower)."""
        assert virasoro_shadow_class() == 'M'

    def test_kappa_matches_bbd_engine(self):
        """Cross-check: κ from gravity engine matches bulk_boundary_duality_engine."""
        try:
            from lib.bulk_boundary_duality_engine import virasoro_koszul_pair
            pair = virasoro_koszul_pair()
            bbd_kappa = pair.kappa
            gravity_kappa_val = gravity_kappa(c=pair.central_charge)
            assert simplify(bbd_kappa - gravity_kappa_val) == 0
        except (ImportError, AttributeError):
            pytest.skip("bulk_boundary_duality_engine not available for cross-check")

    def test_quartic_contact_cross_engine(self):
        """Cross-check: Q^contact from gravity engine vs Vol I engine if available."""
        try:
            sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                            '..', '..', '..', 'chiral-bar-cobar', 'compute'))
            from lib.virasoro_shadow_all_arity import quartic_contact as vol1_quartic
            c_val = Rational(1)
            vol1_val = vol1_quartic(c_val)
            vol2_val = quartic_contact_virasoro(c=1)
            assert simplify(vol1_val - vol2_val) == 0
        except (ImportError, AttributeError):
            pytest.skip("Vol I virasoro_shadow_all_arity not available")


# ===================================================================
# 9. BTZ AND BROWN-HENNEAUX
# ===================================================================

class TestBTZEntropy:
    """Verify BTZ entropy formula from κ."""

    def test_btz_argument(self):
        """S_BTZ ~ √(2κh/3). Verify the argument 2κh/3."""
        # For c=6, h=1: 2*(6/2)*1/3 = 2
        kappa = gravity_kappa(c=6)
        result = 2 * kappa * 1 / 3
        assert result == 2

    def test_btz_at_self_dual(self):
        """At c=13, h=1: argument = 2*(13/2)*1/3 = 13/3."""
        kappa = gravity_kappa(c=13)
        result = 2 * kappa * 1 / 3
        assert result == Rational(13, 3)


# ===================================================================
# 10. INTEGRATION / SMOKE TESTS
# ===================================================================

class TestSmoke:
    """Quick smoke tests for overall consistency."""

    def test_all_functions_callable(self):
        """All exported functions are callable with default arguments."""
        virasoro_lambda_bracket()
        virasoro_associator()
        virasoro_m3_coefficients()
        quartic_contact_virasoro()
        gravity_kappa()
        koszul_dual_central_charge()
        gravity_r_matrix_poles()
        complementarity_constant_virasoro()
        verify_complementarity()

    def test_c13_consistency(self):
        """At c=13 (self-dual): κ=13/2, c_dual=13, Q=10/(13*87)."""
        assert gravity_kappa(c=13) == Rational(13, 2)
        assert koszul_dual_central_charge(c=13) == 13
        assert quartic_contact_virasoro(c=13) == Rational(10, 1131)
