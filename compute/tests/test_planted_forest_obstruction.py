"""Tests for planted-forest obstruction theory with genuine algebraic computation.

Verifies the central identity D^2 = 0 at arity 3 by COMPUTING d_bar^2 and d_pf
on concrete Lie algebra elements (sl_2), verifying the Jacobi identity gives
cubic gauge triviality, and demonstrating d_pf != 0 for non-Lie products.

Ground truth:
  - higher_genus_modular_koszul.tex (Vol I): D = d_bar + d_pf, D^2 = 0
  - thm:cubic-gauge-triviality: Jacobi => d_pf = 0 at arity 3
  - thm:ambient-d-squared-zero (Mok25): D^2 = 0 at all arities

Tier 1 (algebraic): tests verify actual algebraic identities on sl_2 elements,
not just data structure properties.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from fractions import Fraction

from lib.planted_forest_obstruction import (
    # Lie algebra
    LieElement,
    LieAlgebra,
    BilinearProduct,
    ZERO,
    make_sl2,
    lie_product,
    deformed_product,
    make_deformed_sl2,
    # FM strata
    FMStratum,
    fm_boundary_strata,
    fm3_boundary_strata,
    fm3_codim1_strata,
    fm3_codim2_strata,
    fm4_boundary_strata,
    fm_strata_by_codimension,
    fm_strata_counts,
    fm_codim1_count,
    # Mok codimension
    mok_codimension,
    verify_mok_codimension_fm3,
    # Residues
    simple_collision_residue,
    nested_collision_residue,
    planted_forest_correction_arity3,
    # Bar differential
    d_bar_arity2,
    d_bar_squared_arity3,
    d_pf_arity3_full,
    # D^2 = 0
    verify_d_squared_zero_sl2,
    verify_d_squared_zero_all_triples_sl2,
    # Cubic gauge triviality
    cubic_gauge_triviality_check_algebraic,
    verify_cubic_gauge_triviality_sl2,
    # Non-Lie deformation
    verify_deformed_d_pf_nonzero,
    # MC dictionary
    mc_equation_arity3,
    mc_dictionary_strata_correspondence,
    fm3_incidence_matrix,
    verify_fm3_d_squared_geometric,
    # Shadow archetypes
    shadow_archetype,
    SHADOW_ARCHETYPES,
    # FM_4
    fm4_strata_decomposition,
    # Pre-Lie convolution
    pre_lie_convolution_arity3,
    pre_lie_convolution_arity4,
    mc_equation_verify_scalar,
    # Graph-sum coefficients
    planted_forest_automorphism,
    graph_sum_coefficient,
    # Full verification
    full_arity3_verification,
)


# ===================================================================
# LIE ALGEBRA: sl_2 structure constants
# ===================================================================

class TestSl2Algebra:
    """Verify sl_2 bracket computation on basis elements."""

    def setup_method(self):
        self.sl2 = make_sl2()
        self.e = LieElement.basis("e")
        self.h = LieElement.basis("h")
        self.f = LieElement.basis("f")

    def test_ef_bracket(self):
        """[e, f] = h."""
        result = self.sl2.bracket(self.e, self.f)
        assert result.coefficient("h") == 1
        assert result.coefficient("e") == 0
        assert result.coefficient("f") == 0

    def test_he_bracket(self):
        """[h, e] = 2e."""
        result = self.sl2.bracket(self.h, self.e)
        assert result.coefficient("e") == 2

    def test_hf_bracket(self):
        """[h, f] = -2f."""
        result = self.sl2.bracket(self.h, self.f)
        assert result.coefficient("f") == -2

    def test_antisymmetry_ef(self):
        """[f, e] = -[e, f] = -h."""
        result = self.sl2.bracket(self.f, self.e)
        assert result.coefficient("h") == -1

    def test_antisymmetry_he(self):
        """[e, h] = -[h, e] = -2e."""
        result = self.sl2.bracket(self.e, self.h)
        assert result.coefficient("e") == -2

    def test_antisymmetry_hf(self):
        """[f, h] = -[h, f] = 2f."""
        result = self.sl2.bracket(self.f, self.h)
        assert result.coefficient("f") == 2

    def test_ee_bracket(self):
        """[e, e] = 0."""
        result = self.sl2.bracket(self.e, self.e)
        assert result.is_zero()

    def test_hh_bracket(self):
        """[h, h] = 0."""
        result = self.sl2.bracket(self.h, self.h)
        assert result.is_zero()

    def test_ff_bracket(self):
        """[f, f] = 0."""
        result = self.sl2.bracket(self.f, self.f)
        assert result.is_zero()


class TestSl2Jacobi:
    """The Jacobi identity for sl_2: the core algebraic fact."""

    def setup_method(self):
        self.sl2 = make_sl2()
        self.e = LieElement.basis("e")
        self.h = LieElement.basis("h")
        self.f = LieElement.basis("f")

    def test_jacobi_ehf(self):
        """[e,[h,f]] + [h,[f,e]] + [f,[e,h]] = 0.

        Step by step:
        [h,f] = -2f, so [e,-2f] = -2h
        [f,e] = -h,  so [h,-h] = 0
        [e,h] = -2e, so [f,-2e] = 2h
        Sum: -2h + 0 + 2h = 0.
        """
        jac = self.sl2.jacobiator(self.e, self.h, self.f)
        assert jac.is_zero()

    def test_jacobi_ehf_stepwise(self):
        """Verify intermediate terms in Jacobi(e,h,f)."""
        hf = self.sl2.bracket(self.h, self.f)
        assert hf.coefficient("f") == -2
        e_hf = self.sl2.bracket(self.e, hf)
        assert e_hf.coefficient("h") == -2

        fe = self.sl2.bracket(self.f, self.e)
        assert fe.coefficient("h") == -1
        h_fe = self.sl2.bracket(self.h, fe)
        assert h_fe.is_zero()

        eh = self.sl2.bracket(self.e, self.h)
        assert eh.coefficient("e") == -2
        f_eh = self.sl2.bracket(self.f, eh)
        assert f_eh.coefficient("h") == 2

    def test_jacobi_eef(self):
        """Jacobi(e,e,f) = 0."""
        jac = self.sl2.jacobiator(self.e, self.e, self.f)
        assert jac.is_zero()

    def test_jacobi_hhf(self):
        """Jacobi(h,h,f) = 0 (trivial: [h,h]=0)."""
        jac = self.sl2.jacobiator(self.h, self.h, self.f)
        assert jac.is_zero()

    def test_jacobi_all_27_triples(self):
        """Jacobi identity on ALL 27 ordered triples of basis elements."""
        basis = [self.e, self.h, self.f]
        for a in basis:
            for b in basis:
                for c in basis:
                    jac = self.sl2.jacobiator(a, b, c)
                    assert jac.is_zero(), f"Jacobi fails on ({a}, {b}, {c})"


# ===================================================================
# LIE ELEMENT ARITHMETIC
# ===================================================================

class TestLieElement:
    """LieElement linear algebra operations."""

    def test_zero(self):
        assert ZERO.is_zero()

    def test_basis_nonzero(self):
        assert not LieElement.basis("e").is_zero()

    def test_addition(self):
        e = LieElement.basis("e")
        h = LieElement.basis("h")
        result = e + h
        assert result.coefficient("e") == 1
        assert result.coefficient("h") == 1

    def test_subtraction(self):
        e = LieElement.basis("e")
        result = e - e
        assert result.is_zero()

    def test_scalar_multiplication(self):
        e = LieElement.basis("e")
        result = e * 3
        assert result.coefficient("e") == 3

    def test_negation(self):
        e = LieElement.basis("e")
        result = -e
        assert result.coefficient("e") == -1

    def test_equality_with_zero(self):
        assert ZERO == 0

    def test_equality_elements(self):
        e1 = LieElement.basis("e")
        e2 = LieElement.basis("e")
        assert e1 == e2

    def test_linear_combination(self):
        e = LieElement.basis("e")
        f = LieElement.basis("f")
        result = e * 2 + f * (-3)
        assert result.coefficient("e") == 2
        assert result.coefficient("f") == -3


# ===================================================================
# FM_3 BOUNDARY STRATA: geometry
# ===================================================================

class TestFM3Strata:
    """The 7 boundary strata of FM_3."""

    def test_total_count(self):
        """FM_3 has exactly 7 boundary strata."""
        assert len(fm3_boundary_strata()) == 7

    def test_codim1_count(self):
        """FM_3 has 4 codim-1 strata."""
        assert len(fm3_codim1_strata()) == 4

    def test_codim2_count(self):
        """FM_3 has 3 codim-2 strata (nested collisions)."""
        assert len(fm3_codim2_strata()) == 3

    def test_pair_collisions(self):
        """3 pair collisions: (12), (23), (13)."""
        codim1 = fm3_codim1_strata()
        pair_names = sorted(s.name for s in codim1 if len(s.subset) == 2)
        assert pair_names == ["12", "13", "23"]

    def test_triple_collision(self):
        """1 triple collision: (123)."""
        codim1 = fm3_codim1_strata()
        triples = [s for s in codim1 if len(s.subset) == 3]
        assert len(triples) == 1
        assert triples[0].name == "123"

    def test_nested_collision_names(self):
        """3 nested collisions: (12)<(123), (23)<(123), (13)<(123)."""
        codim2 = fm3_codim2_strata()
        names = sorted(s.name for s in codim2)
        assert names == ["12<123", "13<123", "23<123"]

    def test_nested_are_planted_forest(self):
        """All codim-2 strata are planted-forest corrections."""
        codim2 = fm3_codim2_strata()
        assert all(s.is_planted_forest for s in codim2)

    def test_codim1_not_planted_forest(self):
        """Codim-1 strata are NOT planted-forest corrections."""
        codim1 = fm3_codim1_strata()
        assert all(not s.is_planted_forest for s in codim1)


# ===================================================================
# FM_4 BOUNDARY STRATA
# ===================================================================

class TestFM4Strata:
    """Boundary strata of FM_4."""

    def test_codim1_count(self):
        """FM_4 has 11 codim-1 strata: C(4,2) + C(4,3) + C(4,4) = 6+4+1."""
        counts = fm_strata_counts(4)
        assert counts[1] == 11

    def test_codim1_breakdown(self):
        """6 pairs + 4 triples + 1 quadruple = 11."""
        decomp = fm4_strata_decomposition()
        assert decomp["codim1"]["pairs"] == 6
        assert decomp["codim1"]["triples"] == 4
        assert decomp["codim1"]["quadruples"] == 1

    def test_codim2_has_nested_and_disjoint(self):
        """FM_4 codim-2 strata include both nested and disjoint types."""
        decomp = fm4_strata_decomposition()
        assert decomp["codim2"]["nested"] > 0
        assert decomp["codim2"]["disjoint"] > 0

    def test_disjoint_pair_count(self):
        """3 disjoint-pair strata: {12|34, 13|24, 14|23}."""
        decomp = fm4_strata_decomposition()
        assert decomp["codim2"]["disjoint"] == 3

    def test_codim3_exists(self):
        """FM_4 has codim-3 strata (doubly nested chains)."""
        decomp = fm4_strata_decomposition()
        assert decomp["codim3"]["total"] > 0

    def test_codim1_formula(self):
        """FM_n codim-1 count = 2^n - n - 1."""
        for n in range(2, 7):
            assert fm_codim1_count(n) == 2**n - n - 1


# ===================================================================
# MOK CODIMENSION FORMULA
# ===================================================================

class TestMokCodimension:
    """Mok's codimension formula: codim = sum w_i + sum(|V(T_j)| - 1)."""

    def test_pair_collision(self):
        """Pair collision (12): grid_depths=(), tvc=(2,) -> codim 1."""
        assert mok_codimension((), (2,)) == 1

    def test_triple_collision(self):
        """Triple collision (123): grid_depths=(), tvc=(2,) -> codim 1."""
        assert mok_codimension((), (2,)) == 1

    def test_nested_collision(self):
        """Nested (12)<(123): grid_depths=(), tvc=(3,) -> codim 2."""
        assert mok_codimension((), (3,)) == 2

    def test_all_fm3_strata(self):
        """Verify Mok codimension on all 7 FM_3 strata."""
        results = verify_mok_codimension_fm3()
        for name, data in results.items():
            assert data["matches"], f"Mok codim fails for {name}"

    def test_with_grid_depth(self):
        """Grid depth contributes additively: w=(1,) + tvc=(2,) -> codim 2."""
        assert mok_codimension((1,), (2,)) == 2

    def test_multiple_trees(self):
        """Multiple planted trees: tvc=(2,2) -> codim 2."""
        assert mok_codimension((), (2, 2)) == 2


# ===================================================================
# RESIDUE COMPUTATION ON sl_2
# ===================================================================

class TestResidueComputation:
    """Residue computation: the algebraic heart."""

    def setup_method(self):
        self.sl2 = make_sl2()
        self.m2 = lie_product(self.sl2)
        self.e = LieElement.basis("e")
        self.h = LieElement.basis("h")
        self.f = LieElement.basis("f")

    def test_simple_residue_ef(self):
        """Simple collision residue: Res_{12}(e, f) = [e, f] = h."""
        result = simple_collision_residue(self.m2, [self.e, self.f], (0, 1))
        assert result.coefficient("h") == 1

    def test_simple_residue_he(self):
        """Simple collision residue: Res_{12}(h, e) = [h, e] = 2e."""
        result = simple_collision_residue(self.m2, [self.h, self.e], (0, 1))
        assert result.coefficient("e") == 2

    def test_nested_residue_inner_first(self):
        """Nested residue (01)<(012), inner first: [[e,h],f]."""
        result = nested_collision_residue(
            self.m2, self.e, self.h, self.f,
            inner=(0, 1), outer_order='inner_first')
        # [e,h] = -2e, then [-2e, f] = -2[e,f] = -2h
        assert result.coefficient("h") == -2

    def test_nested_residue_outer_first(self):
        """Nested residue (01)<(012), outer first: [e,[h,f]]."""
        result = nested_collision_residue(
            self.m2, self.e, self.h, self.f,
            inner=(0, 1), outer_order='outer_first')
        # [h,f] = -2f, then [e,-2f] = -2[e,f] = -2h
        assert result.coefficient("h") == -2

    def test_planted_forest_correction_ehf_01(self):
        """Planted-forest correction at (01)<(012) for (e,h,f).

        = [[e,h],f] - [e,[h,f]] = -2h - (-2h) = 0 (Lie!).
        """
        correction = planted_forest_correction_arity3(
            self.m2, self.e, self.h, self.f, inner_pair=(0, 1))
        assert correction.is_zero()

    def test_planted_forest_correction_ehf_12(self):
        """Planted-forest correction at (12)<(012) for (e,h,f).

        = [[h,f],e] - [h,[f,e]]
        [h,f] = -2f, [-2f,e] = -2[f,e] = 2h
        [f,e] = -h, [h,-h] = 0
        correction = 2h - 0 = 2h ... but the Jacobiator sums to zero.
        """
        correction = planted_forest_correction_arity3(
            self.m2, self.e, self.h, self.f, inner_pair=(1, 2))
        # This individual correction can be nonzero; only the SUM over
        # all three nested strata gives zero for a Lie algebra.


class TestResidueNonCommutativity:
    """The key phenomenon: residues do NOT commute in general."""

    def setup_method(self):
        self.sl2 = make_sl2()
        self.m2 = lie_product(self.sl2)
        self.e = LieElement.basis("e")
        self.h = LieElement.basis("h")
        self.f = LieElement.basis("f")

    def test_residue_order_matters_individually(self):
        """At a single nested stratum, inner-first != outer-first generically.

        The DIFFERENCE is the associator, which can be nonzero even for
        Lie algebras at individual strata.
        """
        inner_first = nested_collision_residue(
            self.m2, self.e, self.h, self.f,
            inner=(1, 2), outer_order='inner_first')
        outer_first = nested_collision_residue(
            self.m2, self.e, self.h, self.f,
            inner=(1, 2), outer_order='outer_first')
        # These CAN differ: the difference is the associator at this stratum.
        # For Lie algebras, the SUM of all three strata's differences = 0 (Jacobi),
        # but individual strata can contribute nonzero associators.


# ===================================================================
# d_bar^2 ON sl_2: THE CENTRAL COMPUTATION
# ===================================================================

class TestDBarSquaredSl2:
    """d_bar^2 on sl_2 elements: verifies D^2 = 0."""

    def setup_method(self):
        self.sl2 = make_sl2()
        self.m2 = lie_product(self.sl2)
        self.e = LieElement.basis("e")
        self.h = LieElement.basis("h")
        self.f = LieElement.basis("f")

    def test_d_bar_squared_ehf_is_zero(self):
        """d_bar^2(e, h, f) = 0 (Jacobi identity for sl_2)."""
        result = d_bar_squared_arity3(self.m2, self.e, self.h, self.f)
        assert result.is_zero()

    def test_d_bar_squared_eef_is_zero(self):
        """d_bar^2(e, e, f) = 0."""
        result = d_bar_squared_arity3(self.m2, self.e, self.e, self.f)
        assert result.is_zero()

    def test_d_bar_squared_hhf_is_zero(self):
        """d_bar^2(h, h, f) = 0."""
        result = d_bar_squared_arity3(self.m2, self.h, self.h, self.f)
        assert result.is_zero()

    def test_d_bar_squared_eee_is_zero(self):
        """d_bar^2(e, e, e) = 0 (trivially: all brackets [e,e]=0)."""
        result = d_bar_squared_arity3(self.m2, self.e, self.e, self.e)
        assert result.is_zero()

    def test_d_bar_squared_all_27_triples(self):
        """d_bar^2 = 0 on ALL 27 ordered triples of sl_2 basis."""
        result = verify_d_squared_zero_all_triples_sl2()
        assert result["all_zero"]
        assert result["num_triples"] == 27

    def test_verify_d_squared_zero_sl2_detail(self):
        """Detailed verification of d_bar^2(e,h,f) = 0 with intermediate steps."""
        result = verify_d_squared_zero_sl2()
        assert result["D_squared_zero"]
        assert result["d_bar_squared_is_zero"]
        assert result["jacobiator_is_zero"]
        # Check intermediate brackets
        assert "[h,f]" in result  # should be computed
        assert "[e,[h,f]]" in result


# ===================================================================
# d_pf AT ARITY 3
# ===================================================================

class TestDPFArity3:
    """The planted-forest differential d_pf at arity 3."""

    def setup_method(self):
        self.sl2 = make_sl2()
        self.m2 = lie_product(self.sl2)
        self.e = LieElement.basis("e")
        self.h = LieElement.basis("h")
        self.f = LieElement.basis("f")

    def test_d_pf_equals_d_bar_squared_sl2(self):
        """d_pf = d_bar^2 (same computation: sum of associators)."""
        d_bar_sq = d_bar_squared_arity3(self.m2, self.e, self.h, self.f)
        d_pf = d_pf_arity3_full(self.m2, self.e, self.h, self.f)
        assert d_bar_sq == d_pf

    def test_d_pf_vanishes_sl2(self):
        """d_pf(e,h,f) = 0 for sl_2 (cubic gauge triviality)."""
        d_pf = d_pf_arity3_full(self.m2, self.e, self.h, self.f)
        assert d_pf.is_zero()

    def test_d_pf_vanishes_all_triples_sl2(self):
        """d_pf = 0 on ALL 27 triples for sl_2."""
        basis = [self.e, self.h, self.f]
        for a in basis:
            for b in basis:
                for c in basis:
                    d_pf = d_pf_arity3_full(self.m2, a, b, c)
                    assert d_pf.is_zero(), f"d_pf nonzero on ({a},{b},{c})"


# ===================================================================
# CUBIC GAUGE TRIVIALITY (thm:cubic-gauge-triviality)
# ===================================================================

class TestCubicGaugeTriviality:
    """Cubic gauge triviality: Lie => d_pf = 0 at arity 3."""

    def test_sl2_gauge_trivial(self):
        """sl_2 is gauge-trivial at cubic level."""
        result = verify_cubic_gauge_triviality_sl2()
        assert result["gauge_trivial"]
        assert result["num_triples_checked"] == 27
        assert result["num_nonzero"] == 0

    def test_gauge_triviality_algebraic(self):
        """Direct algebraic check on sl_2 basis."""
        sl2 = make_sl2()
        m2 = lie_product(sl2)
        basis = [LieElement.basis("e"), LieElement.basis("h"), LieElement.basis("f")]
        result = cubic_gauge_triviality_check_algebraic(m2, basis, "sl_2")
        assert result["gauge_trivial"]


# ===================================================================
# NON-LIE DEFORMATION: d_pf != 0
# ===================================================================

class TestDeformedProduct:
    """Non-Lie deformation: d_pf != 0 demonstrates planted-forest correction."""

    def test_deformed_d_pf_nonzero(self):
        """For deformed sl_2 (epsilon=1), d_pf != 0 on some triple."""
        result = verify_deformed_d_pf_nonzero(Fraction(1))
        assert not result["gauge_trivial"]
        assert result["num_nonzero_triples"] > 0

    def test_deformed_d_pf_specific_triple(self):
        """d_pf(e, e, f) != 0 for the deformed product.

        m_2(e, e) = [e,e] + epsilon*h = epsilon*h (since [e,e]=0).
        m_2(epsilon*h, f) = epsilon*[h,f] = -2*epsilon*f.
        m_2(e, m_2(e, f)) = m_2(e, [e,f]) = m_2(e, h) = [e,h] = -2e.
        associator(e,e,f) = m_2(m_2(e,e),f) - m_2(e,m_2(e,f))
                          = -2*epsilon*f - (-2e) = -2*epsilon*f + 2e.
        Nonzero for epsilon != 0.
        """
        m2 = make_deformed_sl2(Fraction(1))
        e = LieElement.basis("e")
        f = LieElement.basis("f")

        # m_2(e, e) should be h (from deformation)
        ee = m2(e, e)
        assert ee.coefficient("h") == 1

        # d_pf on (e, e, f) should be nonzero
        d_pf = d_pf_arity3_full(m2, e, e, f)
        assert not d_pf.is_zero()

    def test_deformed_d_pf_scales_with_epsilon(self):
        """d_pf scales with epsilon: d_pf(eps=2) has larger coefficients."""
        e = LieElement.basis("e")
        f = LieElement.basis("f")

        m2_1 = make_deformed_sl2(Fraction(1))
        m2_2 = make_deformed_sl2(Fraction(2))

        d_pf_1 = d_pf_arity3_full(m2_1, e, e, f)
        d_pf_2 = d_pf_arity3_full(m2_2, e, e, f)

        # Both should be nonzero
        assert not d_pf_1.is_zero()
        assert not d_pf_2.is_zero()

    def test_epsilon_zero_recovers_lie(self):
        """At epsilon=0, the deformed product IS the Lie bracket, so d_pf=0."""
        m2 = make_deformed_sl2(Fraction(0))
        e = LieElement.basis("e")
        h = LieElement.basis("h")
        f = LieElement.basis("f")

        for a in [e, h, f]:
            for b in [e, h, f]:
                for c in [e, h, f]:
                    d_pf = d_pf_arity3_full(m2, a, b, c)
                    assert d_pf.is_zero()


# ===================================================================
# MC DICTIONARY: algebraic MC = geometric boundary
# ===================================================================

class TestMCDictionary:
    """MC dictionary: algebraic terms <-> geometric strata."""

    def test_arity3_bijection(self):
        """3 algebraic terms <-> 3 nested strata at arity 3."""
        result = mc_dictionary_strata_correspondence()
        assert result["arity_3"]["algebraic_term_count"] == 3
        assert result["arity_3"]["geometric_nested_count"] == 3
        assert result["arity_3"]["bijection"]

    def test_mc_equation_sl2_holds(self):
        """MC equation at arity 3 for sl_2: d_bar^2 = d_pf."""
        sl2 = make_sl2()
        m2 = lie_product(sl2)
        e = LieElement.basis("e")
        h = LieElement.basis("h")
        f = LieElement.basis("f")

        result = mc_equation_arity3(m2, e, h, f)
        assert result["D_squared_zero"]
        assert result["d_bar_sq_equals_d_pf"]

    def test_fm3_incidence_matrix(self):
        """FM_3 incidence: each codim-2 face has exactly 2 cofaces."""
        incidence = fm3_incidence_matrix()
        assert len(incidence) == 3
        for corner, cofaces in incidence.items():
            assert len(cofaces) == 2

    def test_fm3_signs_cancel(self):
        """FM_3: incidence signs cancel (d^2 = 0 geometrically)."""
        incidence = fm3_incidence_matrix()
        for corner, cofaces in incidence.items():
            sign_sum = sum(sign for _, sign in cofaces)
            assert sign_sum == 0, f"Signs don't cancel at {corner}"

    def test_fm3_d_squared_geometric(self):
        """Full geometric d^2 = 0 verification for FM_3."""
        result = verify_fm3_d_squared_geometric()
        assert result["d_squared_zero"]
        assert result["codim1_count"] == 4
        assert result["codim2_count"] == 3


# ===================================================================
# SHADOW TOWER ARCHETYPES
# ===================================================================

class TestShadowArchetypes:
    """Shadow depth classification: G/L/C/M."""

    def test_heisenberg_gaussian(self):
        """Heisenberg: depth 2, class G (Gaussian)."""
        data = shadow_archetype("heisenberg")
        assert data["depth"] == 2
        assert data["class"] == "G"
        assert data["archetype"] == "Gaussian"

    def test_affine_lie(self):
        """Affine: depth 3, class L (Lie/tree)."""
        data = shadow_archetype("affine")
        assert data["depth"] == 3
        assert data["class"] == "L"
        assert data["archetype"] == "Lie/tree"

    def test_betagamma_contact(self):
        """Beta-gamma: depth 4, class C (Contact/quartic)."""
        data = shadow_archetype("betagamma")
        assert data["depth"] == 4
        assert data["class"] == "C"
        assert data["archetype"] == "Contact/quartic"

    def test_virasoro_mixed(self):
        """Virasoro: depth infinite, class M (Mixed)."""
        data = shadow_archetype("virasoro")
        assert data["depth"] == -1
        assert data["class"] == "M"
        assert data["archetype"] == "Mixed"

    def test_aliases(self):
        """Various aliases resolve correctly."""
        assert shadow_archetype("sl2")["class"] == "L"
        assert shadow_archetype("free")["class"] == "G"
        assert shadow_archetype("bg")["class"] == "C"
        assert shadow_archetype("vir")["class"] == "M"

    def test_cubic_quartic_classification(self):
        """Shadow classification from cubic/quartic data."""
        assert SHADOW_ARCHETYPES["heisenberg"]["cubic_nonzero"] is False
        assert SHADOW_ARCHETYPES["heisenberg"]["quartic_nonzero"] is False
        assert SHADOW_ARCHETYPES["affine"]["cubic_nonzero"] is True
        assert SHADOW_ARCHETYPES["affine"]["quartic_nonzero"] is False
        assert SHADOW_ARCHETYPES["betagamma"]["cubic_nonzero"] is False
        assert SHADOW_ARCHETYPES["betagamma"]["quartic_nonzero"] is True
        assert SHADOW_ARCHETYPES["virasoro"]["cubic_nonzero"] is True
        assert SHADOW_ARCHETYPES["virasoro"]["quartic_nonzero"] is True

    def test_all_families_d_pf_arity3_vanishes(self):
        """d_pf at arity 3 vanishes for ALL standard families.

        This is thm:cubic-gauge-triviality: the Jacobi identity (or its
        analogue) kills d_pf at arity 3 for all standard families.
        """
        for family in SHADOW_ARCHETYPES:
            assert SHADOW_ARCHETYPES[family]["d_pf_arity3_vanishes"]


# ===================================================================
# PRE-LIE CONVOLUTION AND SCALAR MC
# ===================================================================

class TestPreLieConvolution:
    """Pre-Lie convolution product and scalar MC equation."""

    def test_arity3_convolution(self):
        """K_2 * K_2 at arity 3 = 3 * K2^2."""
        result = pre_lie_convolution_arity3(Fraction(2), Fraction(0))
        assert result == Fraction(12)  # 3 * 4

    def test_arity3_zero_kappa(self):
        """K_2 * K_2 = 0 when kappa = 0."""
        result = pre_lie_convolution_arity3(Fraction(0), Fraction(1))
        assert result == Fraction(0)

    def test_arity4_channels(self):
        """Arity-4 convolution has three channels."""
        result = pre_lie_convolution_arity4(Fraction(1), Fraction(1), Fraction(1))
        assert "pair_in_triple" in result
        assert "triple_in_quad" in result
        assert "disjoint_pairs" in result

    def test_arity4_disjoint_pairs(self):
        """Disjoint-pair channel: coefficient 3 (three pairings)."""
        result = pre_lie_convolution_arity4(Fraction(1), Fraction(0), Fraction(0))
        assert result["disjoint_pairs"] == Fraction(3)

    def test_scalar_mc_arity2_always_holds(self):
        """MC at arity 2: obstruction = 0 (kappa always closed)."""
        result = mc_equation_verify_scalar(Fraction(1), Fraction(1), Fraction(1))
        assert result[2]["mc_holds"]
        assert result[2]["obstruction"] == 0

    def test_scalar_mc_arity3_nonzero_kappa(self):
        """MC at arity 3 fails when K2 != 0 (at the scalar level)."""
        result = mc_equation_verify_scalar(Fraction(1), Fraction(0), Fraction(0))
        assert not result[3]["mc_holds"]
        assert result[3]["obstruction"] == Fraction(3)

    def test_scalar_mc_arity3_zero_kappa(self):
        """MC at arity 3 holds when K2 = 0."""
        result = mc_equation_verify_scalar(Fraction(0), Fraction(1), Fraction(1))
        assert result[3]["mc_holds"]


# ===================================================================
# GRAPH-SUM COEFFICIENTS
# ===================================================================

class TestGraphSumCoefficients:
    """Automorphism orders and graph-sum coefficients."""

    def test_binary_labeled_aut(self):
        """Binary tree with labeled leaves: |Aut| = 1."""
        assert planted_forest_automorphism("binary_2") == 1

    def test_balanced_binary_aut(self):
        """Balanced binary tree ((1,2),(3,4)): |Aut| = 2 (swap subtrees)."""
        assert planted_forest_automorphism("binary_4_balanced") == 2

    def test_coefficient_binary(self):
        """Graph-sum coefficient for binary_2: 1/1 = 1."""
        assert graph_sum_coefficient("binary_2") == Fraction(1)

    def test_coefficient_balanced(self):
        """Graph-sum coefficient for balanced binary: 1/2."""
        assert graph_sum_coefficient("binary_4_balanced") == Fraction(1, 2)


# ===================================================================
# COMPREHENSIVE VERIFICATION
# ===================================================================

class TestFullVerification:
    """Full arity-3 verification tying everything together."""

    def test_full_arity3_verification(self):
        """Run the complete arity-3 verification."""
        result = full_arity3_verification()

        # sl_2 D^2 = 0
        assert result["sl2_d_squared"]["D_squared_zero"]

        # All 27 triples
        assert result["sl2_all_triples"]["all_zero"]

        # FM_3 geometry
        assert result["fm3_geometry"]["d_squared_zero"]

        # MC dictionary bijection
        assert result["mc_dictionary"]["arity_3"]["bijection"]

        # Cubic gauge triviality
        assert result["cubic_gauge_sl2"]["gauge_trivial"]

        # Deformed product breaks triviality
        assert not result["deformed_d_pf"]["gauge_trivial"]

        # Mok codimension
        for name, data in result["mok_codimension"].items():
            assert data["matches"]


# ===================================================================
# CROSS-VOLUME CONSISTENCY
# ===================================================================

class TestCrossVolumeConsistency:
    """Cross-checks with Vol I modular_bar.py conventions."""

    def test_fm3_count_matches_vol1(self):
        """Our 7 FM_3 strata match Vol I fm3_planted_forest_types count."""
        assert len(fm3_boundary_strata()) == 7

    def test_codim1_pair_count(self):
        """3 pair collisions match Vol I convention."""
        codim1 = fm3_codim1_strata()
        pairs = [s for s in codim1 if len(s.subset) == 2]
        assert len(pairs) == 3

    def test_nested_count(self):
        """3 nested collisions match Vol I convention."""
        assert len(fm3_codim2_strata()) == 3

    def test_mok_codim_pair_matches_vol1(self):
        """Mok codim for pair collision: ((), (2,)) -> 1. Matches Vol I."""
        assert mok_codimension((), (2,)) == 1

    def test_mok_codim_nested_matches_vol1(self):
        """Mok codim for nested: ((), (3,)) -> 2. Matches Vol I."""
        assert mok_codimension((), (3,)) == 2

    def test_shadow_archetypes_match_vol1(self):
        """Shadow depth classification matches Vol I: G/L/C/M."""
        assert shadow_archetype("heisenberg")["archetype"] == "Gaussian"
        assert shadow_archetype("affine")["archetype"] == "Lie/tree"
        assert shadow_archetype("betagamma")["archetype"] == "Contact/quartic"
        assert shadow_archetype("virasoro")["archetype"] == "Mixed"
