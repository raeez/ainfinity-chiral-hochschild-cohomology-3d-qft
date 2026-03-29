r"""Tests for LG cubic A-infinity chain-level computations.

Verifies the A-infinity structure for the Landau-Ginzburg model with
cubic superpotential W(phi) = g*phi^3/3:

  1. Superpotential derivatives
  2. m_1 = 0 (linearized BRST at phi=0) and m_1^2 = 0
  3. m_2 = free commutative product, associativity
  4. m_3 = 2g (cubic vertex from FM_3(C))
  5. m_{k>=4} = 0 (form-degree vanishing)
  6. A-infinity identities at n=1,2,3,4
  7. Jacobian ring Jac(W) = C[phi]/(phi^2)
  8. Descent to PVA (trivial brackets)
  9. FM_3(C) residue computation

References:
  Vol II: pva-descent.tex, ainfty.py
  Vol I: bar_construction.tex, configuration_spaces.tex
  Costello-Gwilliam (2017): Factorization Algebras in QFT
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy import Symbol, Rational, simplify, expand, S, symbols, diff


# ===================================================================
# 1. SUPERPOTENTIAL AND DERIVATIVES
# ===================================================================

class TestSuperpotential:
    """Verify cubic superpotential and its derivatives."""

    def test_W_at_zero(self):
        """W(0) = 0 (no vacuum energy)."""
        from lib.lg_ainfty_chain_level import cubic_superpotential
        g = Symbol('g')
        assert cubic_superpotential(0, g) == 0

    def test_W_value(self):
        """W(phi) = g*phi^3/3."""
        from lib.lg_ainfty_chain_level import cubic_superpotential
        g, phi = symbols('g phi')
        assert simplify(cubic_superpotential(phi, g) - g * phi**3 / 3) == 0

    def test_W_prime(self):
        """W'(phi) = g*phi^2."""
        from lib.lg_ainfty_chain_level import cubic_W_prime
        g, phi = symbols('g phi')
        assert simplify(cubic_W_prime(phi, g) - g * phi**2) == 0

    def test_W_prime_at_zero(self):
        """W'(0) = 0 (phi=0 is critical)."""
        from lib.lg_ainfty_chain_level import cubic_W_prime
        g = Symbol('g')
        assert cubic_W_prime(0, g) == 0

    def test_W_double_prime(self):
        """W''(phi) = 2g*phi."""
        from lib.lg_ainfty_chain_level import cubic_W_double_prime
        g, phi = symbols('g phi')
        assert simplify(cubic_W_double_prime(phi, g) - 2 * g * phi) == 0

    def test_W_double_prime_at_zero(self):
        """W''(0) = 0 (linearization trivial)."""
        from lib.lg_ainfty_chain_level import cubic_W_double_prime
        g = Symbol('g')
        assert cubic_W_double_prime(0, g) == 0

    def test_W_triple_prime(self):
        """W'''(phi) = 2g (constant)."""
        from lib.lg_ainfty_chain_level import cubic_W_triple_prime
        g = Symbol('g')
        assert cubic_W_triple_prime(g) == 2 * g

    def test_W_triple_prime_numerical(self):
        """W''' = 2g at g=3: W''' = 6."""
        from lib.lg_ainfty_chain_level import cubic_W_triple_prime
        assert cubic_W_triple_prime(3) == 6

    def test_verify_derivatives(self):
        """Cross-check all derivative formulas."""
        from lib.lg_ainfty_chain_level import verify_superpotential_derivatives
        g = Symbol('g')
        result = verify_superpotential_derivatives(g)
        assert result['checks']['W_prime_correct'] is True
        assert result['checks']['W_double_prime_correct'] is True
        assert result['checks']['W_triple_prime_correct'] is True
        assert result['checks']['linearization_trivial'] is True

    def test_sympy_derivative_consistency(self):
        """Verify d^3(g*phi^3/3)/dphi^3 = 2g via sympy diff."""
        g, x = symbols('g x')
        W = g * x**3 / 3
        assert simplify(diff(W, x, 3) - 2 * g) == 0


# ===================================================================
# 2. DIFFERENTIAL m_1
# ===================================================================

class TestM1:
    """Test the linearized BRST differential m_1."""

    def test_m1_phi_zero(self):
        """m_1(phi) = 0 at the linearized level."""
        from lib.lg_ainfty_chain_level import m1_on_polynomial
        g = Symbol('g')
        assert m1_on_polynomial(1, 'phi', g) == 0

    def test_m1_psi_zero(self):
        """m_1(psi) = 0 at the linearized level."""
        from lib.lg_ainfty_chain_level import m1_on_polynomial
        g = Symbol('g')
        assert m1_on_polynomial(1, 'psi', g) == 0

    def test_m1_phi_squared_zero(self):
        """m_1(phi^2) = 0."""
        from lib.lg_ainfty_chain_level import m1_on_polynomial
        g = Symbol('g')
        assert m1_on_polynomial(2, 'phi', g) == 0

    def test_m1_squared_zero(self):
        """m_1^2 = 0 (fundamental identity)."""
        from lib.lg_ainfty_chain_level import verify_m1_squared
        result = verify_m1_squared()
        assert result['is_zero'] is True
        assert result['m1_squared_phi'] == 0
        assert result['m1_squared_psi'] == 0

    def test_m0_zero(self):
        """m_0 = W(0) = 0 (uncurved)."""
        from lib.lg_ainfty_chain_level import verify_m1_squared
        result = verify_m1_squared()
        assert result['curved'] is False
        assert result['m0'] == 0


# ===================================================================
# 3. BINARY OPERATION m_2
# ===================================================================

class TestM2:
    """Test the binary operation m_2."""

    def test_m2_free_product(self):
        """m_2(phi^a, phi^b) = phi^{a+b} (commutative product)."""
        from lib.lg_ainfty_chain_level import m2_free
        assert m2_free(1, 1) == 2
        assert m2_free(1, 2) == 3
        assert m2_free(2, 3) == 5

    def test_m2_commutativity(self):
        """m_2 is commutative: m_2(a,b) = m_2(b,a)."""
        from lib.lg_ainfty_chain_level import m2_free
        assert m2_free(1, 2) == m2_free(2, 1)
        assert m2_free(3, 5) == m2_free(5, 3)

    def test_m2_associativity(self):
        """m_2 is associative: m_2(m_2(a,b),c) = m_2(a,m_2(b,c))."""
        from lib.lg_ainfty_chain_level import m2_free
        a, b, c = 1, 2, 3
        assert m2_free(m2_free(a, b), c) == m2_free(a, m2_free(b, c))

    def test_m2_interaction_zero(self):
        """No tree-level g-correction to m_2."""
        from lib.lg_ainfty_chain_level import m2_interaction
        g = Symbol('g')
        assert m2_interaction(1, 1, g) == 0
        assert m2_interaction(2, 3, g) == 0

    def test_m2_total_equals_free(self):
        """m_2^{total} = m_2^{free} (no g-correction at tree level)."""
        from lib.lg_ainfty_chain_level import m2_total, m2_free
        g = Symbol('g')
        assert m2_total(1, 1, g) == m2_free(1, 1)
        assert m2_total(2, 3, g) == m2_free(2, 3)

    def test_m2_associativity_defect(self):
        """Verify associativity defect is zero."""
        from lib.lg_ainfty_chain_level import verify_m2_associativity_defect
        result = verify_m2_associativity_defect(Symbol('g'))
        assert result['defect'] == 0
        assert result['is_associative'] is True


# ===================================================================
# 4. TERNARY OPERATION m_3
# ===================================================================

class TestM3:
    """Test the cubic interaction m_3."""

    def test_m3_phi_phi_phi(self):
        """m_3(phi, phi, phi) = 2g."""
        from lib.lg_ainfty_chain_level import m3_cubic
        g = Symbol('g')
        assert m3_cubic(1, 1, 1, g) == 2 * g

    def test_m3_numerical(self):
        """m_3(phi, phi, phi) = 6 at g=3."""
        from lib.lg_ainfty_chain_level import m3_cubic
        assert m3_cubic(1, 1, 1, 3) == 6

    def test_m3_higher_powers_zero(self):
        """m_3(phi^2, phi, phi) = 0 (vertex is trilinear)."""
        from lib.lg_ainfty_chain_level import m3_cubic
        g = Symbol('g')
        assert m3_cubic(2, 1, 1, g) == 0
        assert m3_cubic(1, 2, 1, g) == 0
        assert m3_cubic(1, 1, 2, g) == 0

    def test_m3_zero_power_zero(self):
        """m_3(1, phi, phi) = 0."""
        from lib.lg_ainfty_chain_level import m3_cubic
        g = Symbol('g')
        assert m3_cubic(0, 1, 1, g) == 0

    def test_m3_all_zeros(self):
        """m_3(1, 1, 1) = 0."""
        from lib.lg_ainfty_chain_level import m3_cubic
        g = Symbol('g')
        assert m3_cubic(0, 0, 0, g) == 0

    def test_m3_basis_elements(self):
        """m_3 on all basis element triples."""
        from lib.lg_ainfty_chain_level import m3_on_basis_elements
        g = Symbol('g')
        result = m3_on_basis_elements(g)
        assert result[(1, 1, 1)] == 2 * g
        assert result[(0, 0, 0)] == 0
        assert result[(0, 1, 1)] == 0
        assert result[(1, 0, 1)] == 0
        assert result[(1, 1, 0)] == 0

    def test_m3_equals_W_triple_prime(self):
        """m_3(phi,phi,phi) = W''' = 2g (the cubic coupling)."""
        from lib.lg_ainfty_chain_level import m3_cubic, cubic_W_triple_prime
        g = Symbol('g')
        assert m3_cubic(1, 1, 1, g) == cubic_W_triple_prime(g)


# ===================================================================
# 5. HIGHER OPERATIONS m_{k>=4}
# ===================================================================

class TestHigherOperations:
    """Test vanishing of m_{k>=4}."""

    def test_m4_vanishes(self):
        """m_4 = 0 by form-degree counting."""
        from lib.lg_ainfty_chain_level import mk_vanishing_proof
        g = Symbol('g')
        result = mk_vanishing_proof(4, g)
        assert result['vanishes'] is True
        assert result['deficit'] == 4

    def test_m5_vanishes(self):
        """m_5 = 0 by form-degree counting."""
        from lib.lg_ainfty_chain_level import mk_vanishing_proof
        result = mk_vanishing_proof(5, Symbol('g'))
        assert result['vanishes'] is True
        assert result['deficit'] == 4

    def test_m10_vanishes(self):
        """m_10 = 0 by form-degree counting."""
        from lib.lg_ainfty_chain_level import mk_vanishing_proof
        result = mk_vanishing_proof(10, Symbol('g'))
        assert result['vanishes'] is True

    def test_m3_does_not_vanish(self):
        """m_3 does NOT vanish (form degree check confirms)."""
        from lib.lg_ainfty_chain_level import mk_vanishing_proof
        g = Symbol('g')
        result = mk_vanishing_proof(3, g)
        assert result['vanishes'] is False
        assert result['m3_value'] == 2 * g

    def test_m2_no_vertex(self):
        """m_2 comes from free propagator, not vertex."""
        from lib.lg_ainfty_chain_level import mk_vanishing_proof
        result = mk_vanishing_proof(2, Symbol('g'))
        assert result['vanishes'] is False

    def test_deficit_always_4(self):
        """Form degree deficit is always 4 for k >= 4."""
        from lib.lg_ainfty_chain_level import mk_vanishing_proof
        g = Symbol('g')
        for k in range(4, 15):
            result = mk_vanishing_proof(k, g)
            assert result['deficit'] == 4, f"Deficit {result['deficit']} at k={k}"

    def test_vertex_count(self):
        """V = k-2, E = k-3 for tree diagrams."""
        from lib.lg_ainfty_chain_level import mk_vanishing_proof
        g = Symbol('g')
        for k in range(3, 10):
            result = mk_vanishing_proof(k, g)
            assert result['V'] == k - 2
            assert result['E'] == k - 3

    def test_euler_relation(self):
        """V - E = 1 (tree condition)."""
        from lib.lg_ainfty_chain_level import mk_vanishing_proof
        g = Symbol('g')
        for k in range(3, 10):
            result = mk_vanishing_proof(k, g)
            assert result['V'] - result['E'] == 1


# ===================================================================
# 6. A-INFINITY IDENTITIES
# ===================================================================

class TestAInfinityIdentities:
    """Verify A-infinity identities."""

    def test_ainfty_n1(self):
        """n=1 identity: m_1^2 = 0."""
        from lib.lg_ainfty_chain_level import verify_ainfty_n1
        result = verify_ainfty_n1(Symbol('g'))
        assert result['satisfied'] is True

    def test_ainfty_n2(self):
        """n=2 identity: m_1(m_2) + m_2(m_1,-) - m_2(-,m_1) = 0."""
        from lib.lg_ainfty_chain_level import verify_ainfty_n2
        result = verify_ainfty_n2(Symbol('g'))
        assert result['satisfied'] is True
        assert result['total'] == 0

    def test_ainfty_n3(self):
        """n=3 identity: associativity of m_2 (since m_1 = 0)."""
        from lib.lg_ainfty_chain_level import verify_ainfty_n3
        result = verify_ainfty_n3(Symbol('g'))
        assert result['satisfied'] is True
        assert result['associativity_defect'] == 0

    def test_ainfty_n4(self):
        """n=4 identity: first nontrivial, involves m_3."""
        from lib.lg_ainfty_chain_level import verify_ainfty_n4
        g = Symbol('g')
        result = verify_ainfty_n4(g)
        assert result['satisfied'] is True

    def test_ainfty_n4_cancellation(self):
        """n=4: m_2(m_3,id) cancels -m_2(id,m_3)."""
        from lib.lg_ainfty_chain_level import verify_ainfty_n4
        g = Symbol('g')
        result = verify_ainfty_n4(g)
        assert result['term1_m2_m3'] + result['term3_m2_m3'] == 0

    def test_ainfty_n4_m3_m2_terms_zero(self):
        """n=4: all m_3(m_2,id,id) terms vanish."""
        from lib.lg_ainfty_chain_level import verify_ainfty_n4
        g = Symbol('g')
        result = verify_ainfty_n4(g)
        assert result['term2_m3_m2'] == 0
        assert result['term4_m3_m2'] == 0
        assert result['term5_m3_m2'] == 0

    def test_ainfty_n4_m4_zero(self):
        """n=4: m_4 contribution is zero (m_4 = 0)."""
        from lib.lg_ainfty_chain_level import verify_ainfty_n4
        result = verify_ainfty_n4(Symbol('g'))
        assert result['m4_contribution'] == 0

    def test_ainfty_n4_numerical(self):
        """n=4 identity at g=1: 2*1 + 0 - 2*1 + 0 + 0 = 0."""
        from lib.lg_ainfty_chain_level import verify_ainfty_n4
        result = verify_ainfty_n4(1)
        assert result['total'] == 0


# ===================================================================
# 7. JACOBIAN RING
# ===================================================================

class TestJacobianRing:
    """Test Jacobian ring computation."""

    def test_jacobian_dimension(self):
        """dim Jac(g*phi^3/3) = 2 (Milnor number of A_2)."""
        from lib.lg_ainfty_chain_level import jacobian_ring
        result = jacobian_ring(Symbol('g'))
        assert result['dimension'] == 2

    def test_jacobian_milnor(self):
        """Milnor number = 2 for cubic singularity."""
        from lib.lg_ainfty_chain_level import jacobian_ring
        result = jacobian_ring(1)
        assert result['milnor_number'] == 2

    def test_jacobian_basis(self):
        """Basis of Jac = {1, phi}."""
        from lib.lg_ainfty_chain_level import jacobian_ring
        result = jacobian_ring(1)
        assert result['basis'] == ['1', 'phi']

    def test_jacobian_singularity_type(self):
        """Singularity type = A_2."""
        from lib.lg_ainfty_chain_level import jacobian_ring
        result = jacobian_ring(1)
        assert result['singularity_type'] == 'A_2'

    def test_jacobian_multiplication(self):
        """phi * phi = 0 in Jac."""
        from lib.lg_ainfty_chain_level import jacobian_ring
        result = jacobian_ring(1)
        assert result['multiplication_table'][('phi', 'phi')] == '0 (mod phi^2)'

    def test_jacobian_g0_degenerate(self):
        """At g=0: Jac is infinite-dimensional."""
        from lib.lg_ainfty_chain_level import jacobian_ring
        result = jacobian_ring(0)
        assert result['dimension'] == float('inf')

    def test_residue_pairing(self):
        """Residue pairing <1, phi> = 1, <phi, phi> = 0."""
        from lib.lg_ainfty_chain_level import jacobian_ring
        result = jacobian_ring(1)
        assert result['residue_pairing'][('1', 'phi')] == 1
        assert result['residue_pairing'][('phi', 'phi')] == 0


# ===================================================================
# 8. DESCENT TO PVA
# ===================================================================

class TestDescentToPVA:
    """Test PVA descent from the A-infinity structure."""

    def test_pva_brackets_trivial(self):
        """All PVA brackets vanish on Jac(W) for cubic LG."""
        from lib.lg_ainfty_chain_level import descent_to_pva_bracket
        result = descent_to_pva_bracket(Symbol('g'))
        assert result['all_brackets_zero'] is True

    def test_pva_phi_phi_zero(self):
        """{phi_lam phi} = 0 on Jac."""
        from lib.lg_ainfty_chain_level import descent_to_pva_bracket
        result = descent_to_pva_bracket(1)
        assert result['pva_brackets'][('phi', 'phi')] == 0

    def test_chain_level_nontrivial(self):
        """Chain-level A-infinity is nontrivial despite trivial PVA."""
        from lib.lg_ainfty_chain_level import descent_to_pva_bracket
        result = descent_to_pva_bracket(1)
        assert result['chain_level_nontrivial'] is True
        assert result['m3_nonzero'] is True

    def test_homotopy_transfer_m2(self):
        """Transferred m_2 on Jac: phi*phi = 0."""
        from lib.lg_ainfty_chain_level import homotopy_transfer_m2
        result = homotopy_transfer_m2(1)
        assert result[('phi', 'phi')] == '0'
        assert result[('1', 'phi')] == 'phi'


# ===================================================================
# 9. FM_3(C) RESIDUE COMPUTATION
# ===================================================================

class TestFM3Residue:
    """Test FM_3(C) residue computation for m_3."""

    def test_fm3_integral_value(self):
        """FM_3(C) integral = 1 (normalized volume)."""
        from lib.lg_ainfty_chain_level import m3_fm3_integral
        result = m3_fm3_integral()
        assert result['integral_value'] == 1

    def test_fm3_topology(self):
        """FM_3 diagram: V=1, E=0, 3 external legs."""
        from lib.lg_ainfty_chain_level import m3_fm3_integral
        result = m3_fm3_integral()
        assert result['vertices'] == 1
        assert result['propagators'] == 0
        assert result['external_legs'] == 3

    def test_fm3_residue_total(self):
        """Total residue = 2g = m_3 value."""
        from lib.lg_ainfty_chain_level import fm3_residue_computation
        g = Symbol('g')
        result = fm3_residue_computation(g)
        assert result['matches_direct'] is True

    def test_fm3_residue_corolla(self):
        """Corolla D_{123} carries the full residue."""
        from lib.lg_ainfty_chain_level import fm3_residue_computation
        g = Symbol('g')
        result = fm3_residue_computation(g)
        assert result['boundary_strata']['D_123']['residue'] == 2 * g

    def test_fm3_pairwise_zero(self):
        """Pairwise collision strata contribute zero."""
        from lib.lg_ainfty_chain_level import fm3_residue_computation
        g = Symbol('g')
        result = fm3_residue_computation(g)
        assert result['boundary_strata']['D_12']['residue'] == 0
        assert result['boundary_strata']['D_13']['residue'] == 0
        assert result['boundary_strata']['D_23']['residue'] == 0

    def test_fm3_residue_numerical(self):
        """At g=5: total residue = 10."""
        from lib.lg_ainfty_chain_level import fm3_residue_computation
        result = fm3_residue_computation(5)
        assert result['total_residue'] == 10


# ===================================================================
# 10. COMPLETE SUMMARY
# ===================================================================

class TestCompleteSummary:
    """Test the complete A-infinity summary."""

    def test_summary_operations(self):
        """Summary lists all operations correctly."""
        from lib.lg_ainfty_chain_level import lg_cubic_ainfty_summary
        g = Symbol('g')
        result = lg_cubic_ainfty_summary(g)
        assert result['operations']['m0'] == 0
        assert result['operations']['m1'] == 0
        assert result['operations']['m3'] == 2 * g
        assert result['operations']['m4_and_higher'] == 0

    def test_summary_not_curved(self):
        """LG model is uncurved (m_0 = 0)."""
        from lib.lg_ainfty_chain_level import lg_cubic_ainfty_summary
        result = lg_cubic_ainfty_summary(1)
        assert result['is_curved'] is False

    def test_summary_truncation_level(self):
        """Truncation at m_3 (m_{k>=4} = 0)."""
        from lib.lg_ainfty_chain_level import lg_cubic_ainfty_summary
        result = lg_cubic_ainfty_summary(1)
        assert result['truncation_level'] == 3

    def test_summary_jacobian(self):
        """Jacobian ring data in summary."""
        from lib.lg_ainfty_chain_level import lg_cubic_ainfty_summary
        result = lg_cubic_ainfty_summary(1)
        assert result['jacobian_ring']['dimension'] == 2
        assert result['jacobian_ring']['type'] == 'A_2'

    def test_summary_pva_trivial(self):
        """PVA brackets are trivial."""
        from lib.lg_ainfty_chain_level import lg_cubic_ainfty_summary
        result = lg_cubic_ainfty_summary(1)
        assert result['pva_brackets_trivial'] is True

    def test_summary_chain_nontrivial(self):
        """Chain level is nontrivial (m_3 != 0)."""
        from lib.lg_ainfty_chain_level import lg_cubic_ainfty_summary
        result = lg_cubic_ainfty_summary(1)
        assert result['chain_level_nontrivial'] is True

    def test_summary_all_identities(self):
        """All A-infinity identities hold."""
        from lib.lg_ainfty_chain_level import lg_cubic_ainfty_summary
        result = lg_cubic_ainfty_summary(1)
        for n in [1, 2, 3, 4]:
            assert result['ainfty_identities'][f'n{n}'] is True

    def test_summary_numerical_g2(self):
        """At g=2: m_3 = 4, W''' = 4."""
        from lib.lg_ainfty_chain_level import lg_cubic_ainfty_summary
        result = lg_cubic_ainfty_summary(2)
        assert result['operations']['m3'] == 4
        assert result['W_triple_prime'] == 4
