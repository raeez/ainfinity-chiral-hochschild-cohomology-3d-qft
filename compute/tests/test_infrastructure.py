"""
Tests for enhanced infrastructure: sesquilinearity, symmetrize, PVA checker,
Laplace bridge, FM boundary calculus, and Arnold relations.

Test tier: Tier 1 (structural) — verifying mathematical identities.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy import Symbol, Rational, simplify, expand, S, symbols


# ============================================================
# Sesquilinearity verification (ainfty.py)
# ============================================================

class TestSesquilinearityLeft:
    """Verify left sesquilinearity: {∂a_λ b} = -λ{a_λ b}."""

    def test_virasoro_left_sesqui(self):
        """For Virasoro: {∂T_λ T} = -λ{T_λ T}."""
        from lib.ainfty import verify_sesquilinearity_left
        from lib.examples.virasoro import _bracket_TT, _bracket_sesqui_left
        lam = Symbol('lambda')
        c = Symbol('c')

        lhs = _bracket_sesqui_left(lam, c)  # {∂T_λ T} computed from definition
        base = _bracket_TT(lam, c)  # {T_λ T}
        assert verify_sesquilinearity_left(lhs, base, lam) == 0

    def test_virasoro_left_sesqui_c1(self):
        """Left sesquilinearity at c=1 (free boson)."""
        from lib.ainfty import verify_sesquilinearity_left
        from lib.examples.virasoro import _bracket_TT, _bracket_sesqui_left
        lam = Symbol('lambda')
        lhs = _bracket_sesqui_left(lam, 1)
        base = _bracket_TT(lam, 1)
        assert verify_sesquilinearity_left(lhs, base, lam) == 0

    def test_abelian_cs_left_sesqui(self):
        """For abelian CS: {∂J_λ J} = -λ · kλ = -kλ².

        The bracket {J_λ J} = kλ.
        Left sesqui gives: {∂J_λ J} = -λ · kλ = -kλ².
        """
        from lib.ainfty import verify_sesquilinearity_left
        lam = Symbol('lambda')
        k = Symbol('k')
        base = k * lam  # {J_λ J} = kλ
        lhs = -k * lam**2  # {∂J_λ J} = -λ · kλ
        assert verify_sesquilinearity_left(lhs, base, lam) == 0

    def test_su2_left_sesqui_diagonal(self):
        """For su(2): {∂J^1_λ J^1} = -λ · kλ = -kλ²."""
        from lib.ainfty import verify_sesquilinearity_left
        lam = Symbol('lambda')
        k = Symbol('k')
        base = k * lam  # {J^1_λ J^1} = kδ^{11}λ = kλ
        lhs = -k * lam**2  # left sesqui: -λ · base
        assert verify_sesquilinearity_left(lhs, base, lam) == 0


class TestSesquilinearityRight:
    """Verify right sesquilinearity: {a_λ ∂b} = (λ+∂){a_λ b}."""

    def test_virasoro_right_sesqui(self):
        """For Virasoro: {T_λ ∂T} = (λ+∂){T_λ T}."""
        from lib.ainfty import verify_sesquilinearity_right
        from lib.examples.virasoro import (
            _bracket_TT, _bracket_sesqui_right, _apply_partial
        )
        lam = Symbol('lambda')
        c = Symbol('c')

        lhs = _bracket_sesqui_right(lam, c)
        base = _bracket_TT(lam, c)
        assert verify_sesquilinearity_right(lhs, base, lam, _apply_partial) == 0

    def test_virasoro_right_sesqui_c26(self):
        """Right sesquilinearity at c=26 (bosonic string)."""
        from lib.ainfty import verify_sesquilinearity_right
        from lib.examples.virasoro import (
            _bracket_TT, _bracket_sesqui_right, _apply_partial
        )
        lam = Symbol('lambda')
        lhs = _bracket_sesqui_right(lam, 26)
        base = _bracket_TT(lam, 26)
        assert verify_sesquilinearity_right(lhs, base, lam, _apply_partial) == 0


# ============================================================
# LaurentSeries.symmetrize (spectral.py)
# ============================================================

class TestSymmetrize:
    """Test the symmetrize method for extracting commutative products."""

    def test_even_series_unchanged(self):
        """A series with only even powers is its own symmetrization."""
        from lib.spectral import LaurentSeries
        lam = Symbol('lambda')
        f = LaurentSeries({0: 1, 2: 3}, lam)
        sym = f.symmetrize()
        assert sym.coeff(0) == 1
        assert sym.coeff(2) == 3

    def test_odd_series_vanishes(self):
        """A series with only odd powers symmetrizes to zero."""
        from lib.spectral import LaurentSeries
        lam = Symbol('lambda')
        f = LaurentSeries({1: 2, 3: 4}, lam)
        sym = f.symmetrize()
        assert sym.is_zero()

    def test_mixed_series(self):
        """Symmetrize keeps even powers, kills odd powers."""
        from lib.spectral import LaurentSeries
        lam = Symbol('lambda')
        a, b, c = symbols('a b c')
        f = LaurentSeries({0: a, 1: b, 2: c}, lam)
        sym = f.symmetrize()
        assert simplify(sym.coeff(0) - a) == 0
        assert simplify(sym.coeff(1)) == 0  # odd power killed
        assert simplify(sym.coeff(2) - c) == 0

    def test_virasoro_bracket_symmetrize(self):
        """Symmetrize the Virasoro bracket.

        {T_λ T} = ∂T + 2Tλ + (c/12)λ³
        Even part: ∂T (degree 0)
        Odd part: 2Tλ + (c/12)λ³ (degrees 1, 3)

        So symmetrization gives just ∂T.
        """
        from lib.spectral import LaurentSeries
        from lib.examples.virasoro import T, dT
        lam = Symbol('lambda')
        c = Symbol('c')
        bracket = LaurentSeries(
            {0: dT, 1: 2*T, 3: Rational(1, 12)*c}, lam
        )
        sym = bracket.symmetrize()
        assert simplify(sym.coeff(0) - dT) == 0
        assert simplify(sym.coeff(1)) == 0
        assert simplify(sym.coeff(3)) == 0

    def test_singular_symmetrize(self):
        """Symmetrize works on negative powers too."""
        from lib.spectral import LaurentSeries
        lam = Symbol('lambda')
        f = LaurentSeries({-2: 1, -1: 2}, lam)
        sym = f.symmetrize()
        # (-2) is even: survives. (-1) is odd: killed.
        assert sym.coeff(-2) == 1
        assert sym.coeff(-1) == 0


# ============================================================
# PVA Checker with Jacobi and Leibniz (pva.py)
# ============================================================

class TestPVACheckerEnhanced:
    """Test the enhanced PVA checker with Jacobi and Leibniz."""

    def test_free_theory_all_axioms(self):
        """Free theory: bracket = 0, all axioms trivially satisfied."""
        from lib.pva import PVAChecker
        a = Symbol('a')
        partial = Symbol('partial')

        checker = PVAChecker(
            product=lambda x, y: x * y,
            bracket=lambda x, y, l: S.Zero,
            partial=partial,
            generators=[a],
        )
        results = checker.check_all()
        # All axioms should hold
        for axiom, checks in results.items():
            for entry in checks:
                assert entry[-1] == 0, f"Failed: {axiom} on {entry[:-1]}"

    def test_abelian_cs_jacobi(self):
        """Abelian CS: {J_λ J} = kλ. Jacobi is trivially 0=0=0."""
        from lib.pva import PVAChecker
        J = Symbol('J')
        k = Symbol('k')
        partial = Symbol('partial')

        def bracket(a, b, lam):
            if a == J and b == J:
                return k * lam
            return S.Zero

        checker = PVAChecker(
            product=lambda x, y: x * y,
            bracket=bracket,
            partial=partial,
            generators=[J],
        )
        results = checker.check_all()

        # Jacobi: {J_λ {J_μ J}} = {J_λ kμ} = 0, etc.
        for entry in results['jacobi']:
            assert entry[-1] == 0, f"Jacobi failed: {entry}"

    def test_su2_jacobi_via_pva_checker(self):
        """su(2)_k: verify Jacobi through the PVA checker framework.

        This is the first nontrivial use of the enhanced PVA checker.
        The bracket returns linear combinations of generators, which
        _bracket_extended handles by decomposition.
        """
        from lib.pva import PVAChecker
        from lib.examples.nonabelian_cs import su2_lambda_bracket, J

        k = Symbol('k')
        partial = Symbol('partial')

        # Map generator symbols to indices
        idx = {J[1]: 1, J[2]: 2, J[3]: 3}

        def bracket(a, b, lam):
            if a in idx and b in idx:
                return su2_lambda_bracket(idx[a], idx[b], lam, k)
            # Linear extension for composite expressions
            result = S.Zero
            for gen in [J[1], J[2], J[3]]:
                coeff_a = a.coeff(gen) if hasattr(a, 'coeff') else (1 if a == gen else 0)
                if coeff_a != 0:
                    for gen_b in [J[1], J[2], J[3]]:
                        coeff_b = b.coeff(gen_b) if hasattr(b, 'coeff') else (1 if b == gen_b else 0)
                        if coeff_b != 0:
                            result += coeff_a * coeff_b * su2_lambda_bracket(
                                idx[gen], idx[gen_b], lam, k
                            )
            return expand(result)

        checker = PVAChecker(
            product=lambda x, y: x * y,
            bracket=bracket,
            partial=partial,
            generators=[J[1], J[2], J[3]],
        )
        lam = Symbol('lambda')
        mu = Symbol('mu')
        results = checker.check_all(lam=lam, mu=mu)

        # All 27 Jacobi triples should vanish
        for entry in results['jacobi']:
            assert entry[-1] == 0, \
                f"Jacobi failed on ({entry[0]}, {entry[1]}, {entry[2]}): {entry[3]}"

    def test_check_all_returns_all_axioms(self):
        """Verify that check_all returns all four axiom categories."""
        from lib.pva import PVAChecker
        a = Symbol('a')
        partial = Symbol('partial')
        checker = PVAChecker(
            product=lambda x, y: x * y,
            bracket=lambda x, y, l: S.Zero,
            partial=partial,
            generators=[a],
        )
        results = checker.check_all()
        assert 'commutativity' in results
        assert 'skew_symmetry' in results
        assert 'jacobi' in results
        assert 'leibniz' in results


# ============================================================
# Laplace bridge (laplace_bridge.py)
# ============================================================

class TestLaplaceBridge:
    """Test the Laplace transform bridge BR3."""

    def test_constant_bracket(self):
        """Constant bracket c₀ → c₀/z."""
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        z = Symbol('z')
        a = Symbol('a')
        r = lambda_bracket_to_r_matrix({0: a}, z)
        assert simplify(r - a/z) == 0

    def test_linear_bracket(self):
        """Linear bracket c₁λ → c₁/z²."""
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        z = Symbol('z')
        k = Symbol('k')
        r = lambda_bracket_to_r_matrix({1: k}, z)
        assert simplify(r - k/z**2) == 0

    def test_cubic_bracket(self):
        """Cubic bracket c₃λ³ → 6c₃/z⁴."""
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        z = Symbol('z')
        c = Symbol('c')
        r = lambda_bracket_to_r_matrix({3: c}, z)
        assert simplify(r - 6*c/z**4) == 0

    def test_br3_abelian(self):
        """BR3 for abelian CS: {J_λ J}=kλ gives pre-dlog coefficient k/z²."""
        from lib.laplace_bridge import verify_br3_abelian
        k = Symbol('k')
        z = Symbol('z')
        r, expected, diff = verify_br3_abelian(k, z)
        assert diff == 0

    def test_br3_virasoro(self):
        """BR3 for Virasoro: λ-bracket → OPE via Laplace.

        {T_λ T} = ∂T + 2Tλ + (c/12)λ³
        r(z) = ∂T/z + 2T/z² + (c/2)/z⁴
        """
        from lib.laplace_bridge import verify_br3_virasoro
        c = Symbol('c')
        z = Symbol('z')
        r, expected, diff = verify_br3_virasoro(c, z)
        assert diff == 0

    def test_br3_virasoro_fourth_pole(self):
        """The fourth pole in the Virasoro OPE is c/2, not c/12.

        This is the critical coefficient check:
        (c/12) · 3! = c/12 · 6 = c/2.
        """
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        c = Symbol('c')
        z = Symbol('z')
        # λ³ coefficient in bracket is c/12
        r = lambda_bracket_to_r_matrix({3: Rational(1, 12)*c}, z)
        # Should give (c/12) · 6 / z⁴ = c/2 / z⁴
        assert simplify(r - Rational(1, 2)*c/z**4) == 0

    def test_br3_su2_diagonal(self):
        """BR3 for su(2) diagonal: {J^a_λ J^a} = kλ → r^{aa}(z) = k/z²."""
        from lib.laplace_bridge import verify_br3_su2
        k = Symbol('k')
        z = Symbol('z')
        results = verify_br3_su2(k, z)
        # Diagonal entries
        for a in [1, 2, 3]:
            _, _, diff = results[(a, a)]
            assert diff == 0, f"BR3 failed for ({a},{a})"

    def test_br3_su2_offdiagonal(self):
        """BR3 for su(2) off-diagonal: {J^1_λ J^2} = J^3 → r^{12}(z) = J^3/z."""
        from lib.laplace_bridge import verify_br3_su2
        k = Symbol('k')
        z = Symbol('z')
        results = verify_br3_su2(k, z)
        # Off-diagonal: check (1,2) and (2,3)
        for pair in [(1, 2), (2, 3), (1, 3)]:
            _, _, diff = results[pair]
            assert diff == 0, f"BR3 failed for {pair}"

    def test_ope_to_bracket_roundtrip(self):
        """OPE → bracket → r-matrix should recover the OPE."""
        from lib.laplace_bridge import r_matrix_to_lambda_bracket, lambda_bracket_to_r_matrix
        z = Symbol('z')
        # Virasoro OPE: c_0 = dT, c_1 = 2T, c_3 = c/2
        dT = Symbol('dT')
        T = Symbol('T')
        c = Symbol('c')
        ope = {0: dT, 1: 2*T, 3: Rational(1, 2)*c}
        bracket = r_matrix_to_lambda_bracket(ope)
        # bracket should be {0: dT, 1: 2T, 3: c/12}
        assert simplify(bracket[0] - dT) == 0
        assert simplify(bracket[1] - 2*T) == 0
        assert simplify(bracket[3] - Rational(1, 12)*c) == 0

        # Now Laplace back to r-matrix
        r = lambda_bracket_to_r_matrix(bracket, z)
        expected = dT/z + 2*T/z**2 + Rational(1, 2)*c/z**4
        assert simplify(expand(r - expected)) == 0


# ============================================================
# FM Boundary (fm_boundary.py)
# ============================================================

class TestFMBoundary:
    """Test FM compactification boundary strata enumeration."""

    def test_fm2_boundary(self):
        """FM_2(C) has one boundary stratum: D_{12}."""
        from lib.fm_boundary import boundary_strata, count_boundary_strata
        strata = boundary_strata(2)
        assert len(strata) == 1
        assert frozenset({1, 2}) in strata
        assert count_boundary_strata(2) == 1

    def test_fm3_boundary(self):
        """FM_3(C) has 4 boundary strata: D_{12}, D_{13}, D_{23}, D_{123}."""
        from lib.fm_boundary import boundary_strata, count_boundary_strata
        strata = boundary_strata(3)
        assert len(strata) == 4
        assert count_boundary_strata(3) == 4
        assert frozenset({1, 2}) in strata
        assert frozenset({1, 3}) in strata
        assert frozenset({2, 3}) in strata
        assert frozenset({1, 2, 3}) in strata

    def test_fm4_boundary_count(self):
        """FM_4(C) has 2^4 - 4 - 1 = 11 boundary strata."""
        from lib.fm_boundary import count_boundary_strata
        assert count_boundary_strata(4) == 11

    def test_fm5_boundary_count(self):
        """FM_5(C) has 2^5 - 5 - 1 = 26 boundary strata."""
        from lib.fm_boundary import count_boundary_strata
        assert count_boundary_strata(5) == 26

    def test_boundary_decomposition(self):
        """D_{12} in FM_3(C): fiber=FM_2(C), base=FM_2(C)."""
        from lib.fm_boundary import boundary_decomposition
        dec = boundary_decomposition(3, frozenset({1, 2}))
        assert dec['fiber_arity'] == 2
        assert dec['base_arity'] == 2
        assert dec['fiber_dim'] == 2
        assert dec['base_dim'] == 2

    def test_stokes_term_count_3(self):
        """k=3: A∞ identity has 3(3+1)/2 = 6 terms."""
        from lib.fm_boundary import check_stokes_term_count
        result = check_stokes_term_count(3)
        assert result['match']
        assert result['term_count'] == 6

    def test_stokes_term_count_4(self):
        """k=4: A∞ identity has 4(4+1)/2 = 10 terms."""
        from lib.fm_boundary import check_stokes_term_count
        result = check_stokes_term_count(4)
        assert result['match']
        assert result['term_count'] == 10

    def test_lv_sign_from_stokes(self):
        """LV signs from boundary orientation."""
        from lib.fm_boundary import stokes_sign
        # For k=3, S={1,2} (r=0, s=2, t=1): LV sign = (-1)^{0·2+1} = -1
        assert stokes_sign(3, frozenset({1, 2}), convention='lv') == -1
        # S={2,3} (r=1, s=2, t=0): LV sign = (-1)^{1·2+0} = 1
        assert stokes_sign(3, frozenset({2, 3}), convention='lv') == 1

    def test_codim2_corners_fm3(self):
        """FM_3(C) codim-2 corners: D_{ij} ∩ D_{123} (nested)."""
        from lib.fm_boundary import codim2_corners
        corners = codim2_corners(3)
        # Should have corners: each D_{ij} ⊂ D_{123}
        assert len(corners) >= 3

    def test_ainfty_terms_n2(self):
        """n=2 identity has 3 terms: m_1(m_2(a,b)), m_2(m_1(a),b), m_2(a,m_1(b))."""
        from lib.fm_boundary import ainfty_from_stokes
        terms = ainfty_from_stokes(2)
        assert len(terms) == 3  # 2(2+1)/2 = 3


# ============================================================
# Arnold relations (arnold.py)
# ============================================================

class TestArnoldRelations:
    """Test Arnold relations in H*(Conf_n(C)).

    The Arnold relation is a NONTRIVIAL identity: it does NOT hold in
    the free exterior algebra. It holds because of the partial fraction
    identity z_{ij} + z_{jk} + z_{ki} = 0.
    """

    def test_arnold_123_partial_fractions(self):
        """Arnold relation for (1,2,3) via partial fraction identity."""
        from lib.arnold import arnold_relation_partial_fractions
        result = arnold_relation_partial_fractions(1, 2, 3)
        assert result == 0

    def test_arnold_134_partial_fractions(self):
        """Arnold relation for (1,3,4)."""
        from lib.arnold import arnold_relation_partial_fractions
        result = arnold_relation_partial_fractions(1, 3, 4)
        assert result == 0

    def test_arnold_all_triples_4(self):
        """All Arnold relations hold for 4 points."""
        from lib.arnold import verify_arnold_all_triples
        results = verify_arnold_all_triples(4)
        for triple, result in results.items():
            assert result == 0, f"Arnold relation failed for {triple}"

    def test_arnold_all_triples_5(self):
        """All Arnold relations hold for 5 points."""
        from lib.arnold import verify_arnold_all_triples
        results = verify_arnold_all_triples(5)
        for triple, result in results.items():
            assert result == 0, f"Arnold relation failed for {triple}"

    def test_exterior_nonzero(self):
        """Arnold LHS is NONZERO in the free exterior algebra.

        This is the key point: the Arnold relation is not a tautology.
        It requires the partial fraction identity from Conf_n(C).
        """
        from lib.arnold import arnold_relation_exterior
        result = arnold_relation_exterior(1, 2, 3)
        assert not result.is_zero()  # Nonzero in free exterior algebra

    def test_omega_antisymmetry(self):
        """omega_{ij} = -omega_{ji}."""
        from lib.arnold import DifferentialForm
        w12 = DifferentialForm.omega(1, 2)
        w21 = DifferentialForm.omega(2, 1)
        assert (w12 + w21) == 0

    def test_omega_self_wedge(self):
        """omega_{ij} ^ omega_{ij} = 0."""
        from lib.arnold import DifferentialForm
        w12 = DifferentialForm.omega(1, 2)
        result = w12.wedge(w12)
        assert result == 0

    def test_wedge_anticommutativity(self):
        """omega_{12} ^ omega_{34} = -omega_{34} ^ omega_{12}."""
        from lib.arnold import DifferentialForm
        w12 = DifferentialForm.omega(1, 2)
        w34 = DifferentialForm.omega(3, 4)
        result = w12.wedge(w34) + w34.wedge(w12)
        assert result == 0

    def test_orlik_solomon_3(self):
        """Conf_3(C): 3 generators, 1 Arnold relation."""
        from lib.arnold import orlik_solomon_presentation
        os = orlik_solomon_presentation(3)
        assert os['generators'] == 3
        assert os['relations'] == 1

    def test_orlik_solomon_4(self):
        """Conf_4(C): 6 generators, 4 Arnold relations."""
        from lib.arnold import orlik_solomon_presentation
        os = orlik_solomon_presentation(4)
        assert os['generators'] == 6
        assert os['relations'] == 4

    def test_aos_nested_cancellation(self):
        """AOS cancellation at nested corner D_{12} in D_{123}."""
        from lib.arnold import aos_cancellation_at_corner
        result = aos_cancellation_at_corner(
            3, frozenset({1, 2}), frozenset({1, 2, 3})
        )
        assert result['type'] == 'nested'
        assert result['cancels']

    def test_aos_disjoint_cancellation(self):
        """AOS cancellation at disjoint corner D_{12} and D_{34}."""
        from lib.arnold import aos_cancellation_at_corner
        result = aos_cancellation_at_corner(
            4, frozenset({1, 2}), frozenset({3, 4})
        )
        assert result['type'] == 'disjoint'
        assert result['cancels']
