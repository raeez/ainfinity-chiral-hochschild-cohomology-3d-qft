r"""Tests for the factorization-modular engine.

Verifies computationally testable claims from three foundational chapters:
  1. factorization_swiss_cheese.tex — Three propagators, Arnold defect, curvature
  2. modular_swiss_cheese_operad.tex — Modular homotopy-Koszulity, product decomposition
  3. relative_feynman_transform.tex — Bicomplex, involutivity, genus spectral sequence

Tests cover:
  1. Three-propagator system (rational, holomorphic, Arakelov)
  2. Arnold defect at genus 0 and genus >= 1
  3. Curvature formula d_fib^2 = kappa * omega_g
  4. Product decomposition of mixed operations (FM_k x E_1)
  5. Swiss-cheese directionality (no open-to-closed)
  6. Modular homotopy-Koszulity consequences (Quillen equivalence, FT^2 ~ id)
  7. Relative Feynman transform bicomplex (D_P^2=0, D_Mod^2=0, anticommutator=0)
  8. Homotopy-involutivity of FT_rel
  9. Genus spectral sequence (E_1 page, d_1 obstruction class)
  10. Three routes relationship (surjection, injection, curvature visibility)
  11. Cross-family kappa consistency
  12. Propagator comparison table

References:
  Vol II: factorization_swiss_cheese.tex, modular_swiss_cheese_operad.tex,
          relative_feynman_transform.tex
  Vol I: concordance.tex (Theorem D), configuration_spaces.tex (FM)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy import Symbol, Rational, simplify, expand, S, symbols


# ===================================================================
# 1. THREE-PROPAGATOR SYSTEM — RATIONAL PROPAGATOR
# ===================================================================

class TestRationalPropagator:
    """Verify rational propagator K^(0)(z,w) = dz/(z-w) at genus 0."""

    def test_residue_is_one(self):
        """Rational propagator has simple pole with residue 1."""
        from lib.factorization_modular_engine import rational_propagator_residue
        result = rational_propagator_residue()
        assert result['residue'] == 1

    def test_simple_pole(self):
        """Rational propagator has pole order 1."""
        from lib.factorization_modular_engine import rational_propagator_residue
        result = rational_propagator_residue()
        assert result['pole_order'] == 1

    def test_genus_zero(self):
        """Rational propagator lives at genus 0."""
        from lib.factorization_modular_engine import rational_propagator_residue
        result = rational_propagator_residue()
        assert result['genus'] == 0

    def test_is_single_valued(self):
        """Rational propagator is single-valued (no B-cycles at genus 0)."""
        from lib.factorization_modular_engine import rational_propagator_residue
        result = rational_propagator_residue()
        assert result['single_valued'] is True

    def test_is_holomorphic(self):
        """Rational propagator is holomorphic (no correction needed)."""
        from lib.factorization_modular_engine import rational_propagator_residue
        result = rational_propagator_residue()
        assert result['holomorphic'] is True

    def test_is_flat(self):
        """Rational propagator is flat (Arnold relation holds at genus 0)."""
        from lib.factorization_modular_engine import rational_propagator_residue
        result = rational_propagator_residue()
        assert result['flat'] is True


# ===================================================================
# 2. ARNOLD DEFECT
# ===================================================================

class TestArnoldDefect:
    """Verify Arnold relation defect at various genera."""

    def test_genus0_defect_is_zero(self):
        """At genus 0, Arnold relation holds: defect = 0."""
        from lib.factorization_modular_engine import arnold_defect
        result = arnold_defect(0)
        assert result['defect'] == 0

    def test_genus0_arnold_holds(self):
        """Arnold relation holds exactly at genus 0."""
        from lib.factorization_modular_engine import arnold_defect
        result = arnold_defect(0)
        assert result['arnold_holds_exactly'] is True

    def test_genus1_defect_nonzero(self):
        """At genus 1, Arnold relation fails: defect = omega_g."""
        from lib.factorization_modular_engine import arnold_defect
        result = arnold_defect(1)
        assert result['defect'] == 'omega_g'
        assert result['arnold_holds_exactly'] is False

    def test_genus2_defect_nonzero(self):
        """At genus 2, Arnold relation fails: defect = omega_g."""
        from lib.factorization_modular_engine import arnold_defect
        result = arnold_defect(2)
        assert result['defect'] == 'omega_g'
        assert result['arnold_holds_exactly'] is False

    def test_genus0_curvature_zero(self):
        """At genus 0: curvature is 0 (flat bar complex)."""
        from lib.factorization_modular_engine import arnold_defect
        result = arnold_defect(0)
        assert result['curvature'] == 0


# ===================================================================
# 3. CURVATURE FORMULA d_fib^2 = kappa * omega_g
# ===================================================================

class TestCurvatureFormula:
    """Verify curvature from propagator: d_fib^2 = kappa * omega_g."""

    def test_genus0_curvature_zero(self):
        """At genus 0, d^2 = 0 regardless of kappa."""
        from lib.factorization_modular_engine import curvature_from_propagator
        c = Symbol('c')
        result = curvature_from_propagator(c / 2, 0)
        assert result['curvature'] == 0
        assert result['d_squared_zero'] is True

    def test_genus1_curvature_virasoro(self):
        """At genus 1, d^2 = kappa * omega_1 for Virasoro (kappa = c/2)."""
        from lib.factorization_modular_engine import curvature_from_propagator
        c = Symbol('c')
        result = curvature_from_propagator(c / 2, 1)
        assert simplify(result['curvature'] - c / 2) == 0

    def test_genus1_curvature_heisenberg(self):
        """At genus 1, d^2 = 1 * omega_1 for Heisenberg (kappa = 1)."""
        from lib.factorization_modular_engine import curvature_from_propagator
        result = curvature_from_propagator(1, 1)
        assert result['curvature'] == 1
        assert result['d_squared_zero'] is False

    def test_zero_kappa_always_flat(self):
        """If kappa = 0, bar complex is flat at all genera."""
        from lib.factorization_modular_engine import curvature_from_propagator
        for g in range(5):
            result = curvature_from_propagator(0, g)
            assert result['d_squared_zero'] is True

    def test_genus2_curvature_w3(self):
        """At genus 2, d^2 coefficient is kappa = 5c/6 for W_3."""
        from lib.factorization_modular_engine import curvature_from_propagator
        c = Symbol('c')
        result = curvature_from_propagator(5 * c / 6, 2)
        assert simplify(result['curvature'] - 5 * c / 6) == 0


# ===================================================================
# 4. PRODUCT DECOMPOSITION OF MIXED OPERATIONS
# ===================================================================

class TestProductDecomposition:
    """Verify FM_k(C) x E_1(m) product decomposition."""

    def test_fm2_dim(self):
        """dim_R FM_2(C) = 2 (one complex parameter)."""
        from lib.factorization_modular_engine import fm_real_dimension
        assert fm_real_dimension(2) == 2

    def test_fm3_dim(self):
        """dim_R FM_3(C) = 4."""
        from lib.factorization_modular_engine import fm_real_dimension
        assert fm_real_dimension(3) == 4

    def test_fm4_dim(self):
        """dim_R FM_4(C) = 6."""
        from lib.factorization_modular_engine import fm_real_dimension
        assert fm_real_dimension(4) == 6

    def test_fm1_dim(self):
        """FM_1(C) is a point, dim = 0."""
        from lib.factorization_modular_engine import fm_real_dimension
        assert fm_real_dimension(1) == 0

    def test_fm0_dim(self):
        """FM_0(C) is a point, dim = 0."""
        from lib.factorization_modular_engine import fm_real_dimension
        assert fm_real_dimension(0) == 0

    def test_e1_contractible(self):
        """E_1(m) is contractible for m >= 2 (homotopy dimension 0)."""
        from lib.factorization_modular_engine import e1_dimension
        for m in range(2, 10):
            assert e1_dimension(m) == 0

    def test_mixed_dim_genus0(self):
        """Mixed operation dim = dim FM_k(C) + 0 at genus 0."""
        from lib.factorization_modular_engine import mixed_operation_dim
        # k=3, m=2 at genus 0
        result = mixed_operation_dim(3, 2, g=0)
        assert result['total_dim'] == 4  # 2*3-2 = 4
        assert result['closed_dim'] == 4
        assert result['open_dim'] == 0

    def test_product_decomposition_flag(self):
        """The operation space always factors as product."""
        from lib.factorization_modular_engine import mixed_operation_dim
        result = mixed_operation_dim(2, 3, g=0)
        assert result['product_decomposition'] is True

    def test_fm_dimension_formula_general(self):
        """dim FM_k(C) = 2k-2 for k = 2..10."""
        from lib.factorization_modular_engine import fm_real_dimension
        for k in range(2, 11):
            assert fm_real_dimension(k) == 2 * k - 2


# ===================================================================
# 5. SWISS-CHEESE DIRECTIONALITY
# ===================================================================

class TestDirectionality:
    """Verify Swiss-cheese directionality: no open-to-closed operations."""

    def test_no_open_to_closed(self):
        """Open-to-closed operations are empty (Axiom 5)."""
        from lib.factorization_modular_engine import swiss_cheese_directionality
        result = swiss_cheese_directionality()
        assert result['open_to_closed'] is False

    def test_closed_to_open_exists(self):
        """Closed-to-open (bulk-to-boundary) operations exist."""
        from lib.factorization_modular_engine import swiss_cheese_directionality
        result = swiss_cheese_directionality()
        assert result['closed_to_open'] is True


# ===================================================================
# 6. MODULAR HOMOTOPY-KOSZULITY CONSEQUENCES
# ===================================================================

class TestModularHomotopyKoszulity:
    """Verify consequences of modular homotopy-Koszulity of SC_mod."""

    def test_bar_cobar_quillen_equivalence(self):
        """Bar-cobar for SC_mod is a Quillen equivalence."""
        from lib.factorization_modular_engine import bar_cobar_is_quillen_equivalence
        assert bar_cobar_is_quillen_equivalence() is True

    def test_feynman_involution_is_identity(self):
        """FT_mod^2 ~ id on closed-colour algebras."""
        from lib.factorization_modular_engine import feynman_involution_squared
        assert feynman_involution_squared() == 'id'

    def test_proof_has_four_steps(self):
        """Modular homotopy-Koszulity proof has 4 steps."""
        from lib.factorization_modular_engine import modular_hkoszul_proof_steps
        steps = modular_hkoszul_proof_steps()
        assert len(steps) == 4

    def test_all_steps_proved(self):
        """All four proof steps are established/proved."""
        from lib.factorization_modular_engine import modular_hkoszul_proof_steps
        steps = modular_hkoszul_proof_steps()
        for key, step in steps.items():
            assert step['status'] in ('established', 'proved'), \
                f"Step {key} has status {step['status']}"


# ===================================================================
# 7. RELATIVE FEYNMAN TRANSFORM — BICOMPLEX
# ===================================================================

class TestBicomplex:
    """Verify the bicomplex structure of the relative Feynman transform."""

    def test_bicomplex_D_P_squared_zero(self):
        """D_P^2 = 0 (genus-preserving differential squares to zero)."""
        from lib.factorization_modular_engine import bicomplex_structure
        bc = bicomplex_structure()
        assert bc['D_P_squared'] == 0

    def test_bicomplex_D_Mod_squared_zero(self):
        """D_Mod^2 = 0 (genus-raising differential squares to zero)."""
        from lib.factorization_modular_engine import bicomplex_structure
        bc = bicomplex_structure()
        assert bc['D_Mod_squared'] == 0

    def test_bicomplex_anticommutator_zero(self):
        """D_P D_Mod + D_Mod D_P = 0 (anticommutator vanishes)."""
        from lib.factorization_modular_engine import bicomplex_structure
        bc = bicomplex_structure()
        assert bc['anticommutator'] == 0

    def test_total_D_squared_zero(self):
        """(D_P + D_Mod)^2 = 0 (total differential squares to zero)."""
        from lib.factorization_modular_engine import bicomplex_structure
        bc = bicomplex_structure()
        assert bc['total_D_squared'] == 0

    def test_verify_bicomplex_with_zeros(self):
        """verify_bicomplex accepts (0, 0, 0) as valid bicomplex."""
        from lib.factorization_modular_engine import verify_bicomplex
        result = verify_bicomplex(0, 0, 0)
        assert result['is_bicomplex'] is True
        assert result['D_P_squared_zero'] is True
        assert result['D_Mod_squared_zero'] is True
        assert result['anticommutator_zero'] is True

    def test_verify_bicomplex_rejects_nonzero(self):
        """verify_bicomplex rejects non-zero D_P^2."""
        from lib.factorization_modular_engine import verify_bicomplex
        result = verify_bicomplex(1, 0, 0)
        assert result['is_bicomplex'] is False
        assert result['D_P_squared_zero'] is False

    def test_D_P_genus_preserving(self):
        """D_P = D_int + D_sep preserves genus."""
        from lib.factorization_modular_engine import bicomplex_structure
        bc = bicomplex_structure()
        assert bc['D_P_genus_shift'] == 0

    def test_D_Mod_genus_raising(self):
        """D_Mod = D_nsep raises genus by 1."""
        from lib.factorization_modular_engine import bicomplex_structure
        bc = bicomplex_structure()
        assert bc['D_Mod_genus_shift'] == 1


# ===================================================================
# 8. HOMOTOPY-INVOLUTIVITY
# ===================================================================

class TestInvolutivity:
    """Verify homotopy-involutivity of the relative Feynman transform."""

    def test_ft_rel_involutive(self):
        """FT_rel(FT_rel(A)) ~ A (homotopy-involutive)."""
        from lib.factorization_modular_engine import relative_ft_involutivity
        assert relative_ft_involutivity() is True


# ===================================================================
# 9. GENUS SPECTRAL SEQUENCE
# ===================================================================

class TestGenusSpectralSequence:
    """Verify genus spectral sequence claims."""

    def test_e1_genus0_tree_level(self):
        """E_1 at genus 0 is tree-level bar cohomology."""
        from lib.factorization_modular_engine import e1_page_at_genus
        result = e1_page_at_genus(0)
        assert result['genus'] == 0
        assert result['obstruction_class'] is None

    def test_e1_genus1_has_obstruction(self):
        """E_1 at genus 1 has obstruction class kappa * omega_1."""
        from lib.factorization_modular_engine import e1_page_at_genus
        result = e1_page_at_genus(1)
        assert result['obstruction_class'] == 'kappa * omega_g'

    def test_d1_genus_shift_is_one(self):
        """d_1 differential raises genus by exactly 1."""
        from lib.factorization_modular_engine import e1_page_at_genus
        result = e1_page_at_genus(0)
        assert result['d1_genus_shift'] == 1

    def test_obstruction_virasoro(self):
        """Obstruction at genus 1 for Virasoro: coefficient c/2."""
        from lib.factorization_modular_engine import d1_obstruction_class
        c = Symbol('c')
        result = d1_obstruction_class(1, c / 2)
        assert simplify(result['obstruction_coefficient'] - c / 2) == 0

    def test_obstruction_vanishes_at_kappa_zero(self):
        """Obstruction vanishes when kappa = 0."""
        from lib.factorization_modular_engine import d1_obstruction_class
        result = d1_obstruction_class(1, 0)
        assert result['vanishes'] is True

    def test_obstruction_heisenberg_genus2(self):
        """Heisenberg (kappa=1) has nonzero obstruction at genus 2."""
        from lib.factorization_modular_engine import d1_obstruction_class
        result = d1_obstruction_class(2, 1)
        assert result['obstruction_coefficient'] == 1
        assert result['vanishes'] is False


# ===================================================================
# 10. THREE ROUTES RELATIONSHIP
# ===================================================================

class TestThreeRoutes:
    """Verify the three-routes relationship."""

    def test_factorization_surjects(self):
        """Factorization (Route B) surjects onto algebraic skeleton (Route C)."""
        from lib.factorization_modular_engine import factorization_surjects_onto_ft
        assert factorization_surjects_onto_ft() is True

    def test_operadic_injects(self):
        """Operadic (Route A) injects into algebraic skeleton (Route C)."""
        from lib.factorization_modular_engine import operadic_injects_into_ft
        assert operadic_injects_into_ft() is True

    def test_ft_does_not_see_curvature(self):
        """Relative FT does NOT distinguish flat from curved models."""
        from lib.factorization_modular_engine import ft_sees_curvature
        assert ft_sees_curvature() is False

    def test_route_B_sees_curvature(self):
        """Only factorization (Route B) sees the curvature."""
        from lib.factorization_modular_engine import three_routes
        routes = three_routes()
        assert routes['route_B']['sees_curvature'] is True
        assert routes['route_A']['sees_curvature'] is False
        assert routes['route_C']['sees_curvature'] is False

    def test_all_routes_see_flat_models(self):
        """All three routes see the flat models."""
        from lib.factorization_modular_engine import three_routes
        routes = three_routes()
        for key in ('route_A', 'route_B', 'route_C'):
            assert routes[key]['sees_flat_models'] is True

    def test_only_factorization_sees_curved(self):
        """Only Route B (factorization) sees the curved model."""
        from lib.factorization_modular_engine import three_routes
        routes = three_routes()
        assert routes['route_B']['sees_curved_model'] is True
        assert routes['route_A']['sees_curved_model'] is False
        assert routes['route_C']['sees_curved_model'] is False


# ===================================================================
# 11. CROSS-FAMILY KAPPA CONSISTENCY
# ===================================================================

class TestKappaConsistency:
    """Verify kappa values across families match Theorem D."""

    def test_heisenberg_kappa(self):
        """kappa(Heisenberg) = 1."""
        from lib.factorization_modular_engine import kappa_for_family
        assert kappa_for_family('heisenberg') == 1

    def test_virasoro_kappa(self):
        """kappa(Virasoro) = c/2."""
        from lib.factorization_modular_engine import kappa_for_family
        c = Symbol('c')
        assert simplify(kappa_for_family('virasoro', c=c) - c / 2) == 0

    def test_w3_kappa(self):
        """kappa(W_3) = 5c/6."""
        from lib.factorization_modular_engine import kappa_for_family
        c = Symbol('c')
        assert simplify(kappa_for_family('w3', c=c) - 5 * c / 6) == 0

    def test_betagamma_kappa(self):
        """kappa(betagamma) = -1."""
        from lib.factorization_modular_engine import kappa_for_family
        assert kappa_for_family('betagamma') == -1

    def test_affine_kappa_sl2(self):
        """kappa(sl_2 affine at level k) = 3(k+2)/4."""
        from lib.factorization_modular_engine import kappa_for_family
        k = Symbol('k')
        # sl_2: dim = 3, h^v = 2
        result = kappa_for_family('affine', k=k, dim_g=3, h_dual=2)
        expected = 3 * (k + 2) / 4
        assert simplify(result - expected) == 0


# ===================================================================
# 12. PROPAGATOR COMPARISON TABLE
# ===================================================================

class TestPropagatorTable:
    """Verify the three-propagator comparison table."""

    def test_three_propagators_exist(self):
        """Table has exactly three entries."""
        from lib.factorization_modular_engine import propagator_genus_table
        table = propagator_genus_table()
        assert len(table) == 3

    def test_rational_in_derived_category(self):
        """Rational propagator lives in derived category."""
        from lib.factorization_modular_engine import propagator_genus_table
        table = propagator_genus_table()
        rational = [p for p in table if p['name'] == 'rational'][0]
        assert rational['category'] == 'derived'

    def test_arakelov_in_coderived_category(self):
        """Arakelov propagator lives in coderived category."""
        from lib.factorization_modular_engine import propagator_genus_table
        table = propagator_genus_table()
        arakelov = [p for p in table if p['name'] == 'arakelov'][0]
        assert arakelov['category'] == 'coderived'

    def test_arakelov_visible_only_to_route_B(self):
        """Curved model (Arakelov) is visible only to Route B (factorization)."""
        from lib.factorization_modular_engine import propagator_genus_table
        table = propagator_genus_table()
        arakelov = [p for p in table if p['name'] == 'arakelov'][0]
        assert arakelov['visible_to'] == ['B']

    def test_holomorphic_propagator_properties(self):
        """Holomorphic propagator at genus >= 1 is multi-valued."""
        from lib.factorization_modular_engine import holomorphic_propagator_properties
        for g in range(1, 5):
            props = holomorphic_propagator_properties(g)
            assert props['single_valued'] is False
            assert props['holomorphic'] is True

    def test_arakelov_propagator_properties(self):
        """Arakelov propagator at genus >= 1 is single-valued but not holomorphic."""
        from lib.factorization_modular_engine import arakelov_propagator_properties
        for g in range(1, 5):
            props = arakelov_propagator_properties(g)
            assert props['single_valued'] is True
            assert props['holomorphic'] is False
            assert props['flat'] is False
