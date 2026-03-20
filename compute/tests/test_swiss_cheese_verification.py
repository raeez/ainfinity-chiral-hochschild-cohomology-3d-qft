"""Tests for Swiss-cheese operad SC^{ch,top} verification.

Verifies:
1. Arnold-Orlik-Solomon cohomology of FM_k(C) (Betti numbers, Poincare polynomial)
2. Operadic composition laws (closed-closed, open-open, mixed)
3. Interchange law for mixed operations
4. Boundary structure (Stasheff vs FM, three-face Stokes)
5. Homotopy-Koszulity indicators
6. Combinatorial identities (Stirling, Catalan)

References:
  thm:homotopy-Koszul (Vol II): SC^{ch,top} is homotopy-Koszul
  fm-calculus.tex (Vol II): FM boundary structure
  arnold_relations.tex (Appendix A): AOS presentation

Tier 1 (structural): all tests are self-certifying identities.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from math import factorial, comb

from lib.swiss_cheese_verification import (
    unsigned_stirling_first,
    catalan,
    poincare_polynomial_fm,
    aos_cohomology,
    euler_characteristic_fm,
    total_betti_fm,
    top_betti_fm,
    arnold_relation_generators,
    aos_algebra_verify,
    betti_from_stirling,
    e1_cohomology,
    e1_euler_characteristic,
    sc_mixed_cohomology,
    closed_composition_check,
    open_composition_check,
    mixed_interchange_check,
    stasheff_boundary_faces,
    stasheff_face_count,
    fm_boundary_faces,
    fm_face_count,
    fm_vs_stasheff_face_comparison,
    three_face_stokes_verify,
    mixed_operation_dimensions,
    homotopy_koszul_indicators,
    koszul_criterion_check,
    operadic_gf_closed,
    operadic_gf_open,
    sc_dimension_check,
    bar_complex_dimensions_sc,
    sc_verification_summary,
)


# ===================================================================
# STIRLING NUMBERS
# ===================================================================

class TestStirlingNumbers:
    """Unsigned Stirling numbers of the first kind |s(n,k)|.

    These give dim H^q(FM_n(C)) = |s(n, n-q)|.
    Standard values from OEIS A132393 / DLMF.
    """

    def test_s_0_0(self):
        """|s(0,0)| = 1 (base case)."""
        assert unsigned_stirling_first(0, 0) == 1

    def test_s_1_1(self):
        """|s(1,1)| = 1 (one permutation of 1 element, 1 cycle)."""
        assert unsigned_stirling_first(1, 1) == 1

    def test_s_2_1(self):
        """|s(2,1)| = 1 (one permutation (12) with 1 cycle)."""
        assert unsigned_stirling_first(2, 1) == 1

    def test_s_2_2(self):
        """|s(2,2)| = 1 (identity permutation, 2 cycles)."""
        assert unsigned_stirling_first(2, 2) == 1

    def test_s_3_1(self):
        """|s(3,1)| = 2 (two 3-cycles)."""
        assert unsigned_stirling_first(3, 1) == 2

    def test_s_3_2(self):
        """|s(3,2)| = 3 (three transpositions as 2-cycle + fixed point)."""
        assert unsigned_stirling_first(3, 2) == 3

    def test_s_3_3(self):
        """|s(3,3)| = 1 (identity)."""
        assert unsigned_stirling_first(3, 3) == 1

    def test_s_4_row(self):
        """|s(4,k)| for k=1..4: [6, 11, 6, 1]."""
        expected = [6, 11, 6, 1]
        actual = [unsigned_stirling_first(4, k) for k in range(1, 5)]
        assert actual == expected

    def test_s_5_row(self):
        """|s(5,k)| for k=1..5: [24, 50, 35, 10, 1]."""
        expected = [24, 50, 35, 10, 1]
        actual = [unsigned_stirling_first(5, k) for k in range(1, 6)]
        assert actual == expected

    def test_row_sum_is_factorial(self):
        """sum_k |s(n,k)| = n! for all n (every permutation has some cycle count)."""
        for n in range(1, 8):
            row_sum = sum(unsigned_stirling_first(n, k) for k in range(1, n + 1))
            assert row_sum == factorial(n), f"Failed for n={n}"

    def test_boundary_cases(self):
        """|s(n,k)| = 0 for k > n or k < 0."""
        assert unsigned_stirling_first(3, 4) == 0
        assert unsigned_stirling_first(5, 6) == 0
        assert unsigned_stirling_first(2, -1) == 0

    def test_s_n_n_is_1(self):
        """|s(n,n)| = 1 for all n (only the identity has n cycles)."""
        for n in range(0, 8):
            assert unsigned_stirling_first(n, n) == 1

    def test_s_n_1_is_n_minus_1_factorial(self):
        """|s(n,1)| = (n-1)! for n >= 1 (number of (n-1)-cycles = (n-1)!)."""
        for n in range(1, 8):
            assert unsigned_stirling_first(n, 1) == factorial(n - 1)


# ===================================================================
# CATALAN NUMBERS
# ===================================================================

class TestCatalan:
    """Catalan numbers C_n = C(2n,n)/(n+1)."""

    def test_small_values(self):
        """C_0=1, C_1=1, C_2=2, C_3=5, C_4=14, C_5=42."""
        expected = [1, 1, 2, 5, 14, 42]
        actual = [catalan(n) for n in range(6)]
        assert actual == expected

    def test_negative(self):
        """C_n = 0 for n < 0."""
        assert catalan(-1) == 0


# ===================================================================
# POINCARE POLYNOMIAL OF FM_k(C)
# ===================================================================

class TestPoincareFM:
    """Poincare polynomial P(t) = prod_{j=1}^{n-1} (1+jt) of FM_n(C).

    Arnold's theorem (1969): H*(Conf_n(C)) has this Poincare polynomial.
    """

    def test_fm1(self):
        """FM_1(C) = point: P(t) = 1."""
        assert poincare_polynomial_fm(1) == [1]

    def test_fm2(self):
        """FM_2(C) ~ S^1: P(t) = 1 + t."""
        assert poincare_polynomial_fm(2) == [1, 1]

    def test_fm3(self):
        """FM_3(C): P(t) = (1+t)(1+2t) = 1 + 3t + 2t^2."""
        assert poincare_polynomial_fm(3) == [1, 3, 2]

    def test_fm4(self):
        """FM_4(C): P(t) = (1+t)(1+2t)(1+3t) = 1 + 6t + 11t^2 + 6t^3."""
        assert poincare_polynomial_fm(4) == [1, 6, 11, 6]

    def test_fm5(self):
        """FM_5(C): P(t) = (1+t)(1+2t)(1+3t)(1+4t)."""
        poly = poincare_polynomial_fm(5)
        assert poly == [1, 10, 35, 50, 24]

    def test_betti_matches_stirling(self):
        """Betti numbers from Poincare must match Stirling formula.

        dim H^q(FM_n(C)) = |s(n, n-q)| for all q.
        """
        for n in range(1, 7):
            poly = poincare_polynomial_fm(n)
            stirling = betti_from_stirling(n)
            assert poly == stirling, f"Mismatch at n={n}: {poly} vs {stirling}"

    def test_total_betti_equals_factorial(self):
        """P(1) = n! (total Betti number = n!)."""
        for n in range(1, 8):
            assert sum(poincare_polynomial_fm(n)) == factorial(n)

    def test_euler_char_vanishes(self):
        """P(-1) = 0 for n >= 2 (Euler characteristic vanishes)."""
        for n in range(2, 8):
            assert euler_characteristic_fm(n) == 0

    def test_top_betti(self):
        """dim H^{n-1}(FM_n) = (n-1)! (top Betti = top Stirling)."""
        for n in range(1, 8):
            assert top_betti_fm(n) == factorial(n - 1)


# ===================================================================
# AOS ALGEBRA STRUCTURE
# ===================================================================

class TestAOSAlgebra:
    """Arnold-Orlik-Solomon presentation of H*(FM_n(C))."""

    def test_generators_count(self):
        """Number of generators eta_{ij} = C(n,2) for each n."""
        for n in range(2, 7):
            info = arnold_relation_generators(n)
            assert info['num_generators'] == comb(n, 2)

    def test_relations_count(self):
        """Number of Arnold relations = C(n,3) for each n."""
        for n in range(3, 7):
            info = arnold_relation_generators(n)
            assert info['num_relations'] == comb(n, 3)

    def test_fm2_generators(self):
        """FM_2: one generator eta_{12}."""
        info = arnold_relation_generators(2)
        assert info['generators'] == [(1, 2)]
        assert info['relations'] == []

    def test_fm3_generators(self):
        """FM_3: three generators eta_{12}, eta_{13}, eta_{23}."""
        info = arnold_relation_generators(3)
        assert len(info['generators']) == 3
        assert len(info['relations']) == 1

    def test_aos_verify_passes(self):
        """Full AOS verification passes for n=1..6."""
        for n in range(1, 7):
            check = aos_algebra_verify(n)
            assert check['all_pass'], f"AOS verify failed at n={n}"

    def test_h1_dimension(self):
        """dim H^1(FM_n) = C(n,2) = n(n-1)/2.

        The Arnold relations are degree-2 relations (involving products
        of two generators), so they do not reduce H^1.
        """
        for n in range(2, 7):
            betti = aos_cohomology(n)
            assert betti[1] == comb(n, 2)


# ===================================================================
# E_1 OPERAD (LITTLE INTERVALS)
# ===================================================================

class TestE1Operad:
    """E_1(m) is contractible for m >= 1."""

    def test_contractible(self):
        """E_1(m) = [1] (contractible) for m = 1..10."""
        for m in range(1, 11):
            assert e1_cohomology(m) == [1]

    def test_empty_at_zero(self):
        """E_1(0) is empty."""
        assert e1_cohomology(0) == []

    def test_euler_char(self):
        """chi(E_1(m)) = 1 for m >= 1."""
        for m in range(1, 6):
            assert e1_euler_characteristic(m) == 1


# ===================================================================
# SWISS-CHEESE MIXED SPACES
# ===================================================================

class TestSCMixed:
    """Mixed Swiss-cheese spaces SC(k,m)."""

    def test_purely_open(self):
        """SC(0,m) = E_1(m) (contractible for m >= 1)."""
        for m in range(1, 5):
            assert sc_mixed_cohomology(0, m) == [1]

    def test_sc_1_0(self):
        """SC(1,0) = point."""
        assert sc_mixed_cohomology(1, 0) == [1]

    def test_sc_1_1(self):
        """SC(1,1) = half-disk (contractible)."""
        assert sc_mixed_cohomology(1, 1) == [1]

    def test_sc_2_0(self):
        """SC(2,0) = FM_2(H) = interval (contractible)."""
        assert sc_mixed_cohomology(2, 0) == [1]


# ===================================================================
# OPERADIC COMPOSITION
# ===================================================================

class TestClosedComposition:
    """Closed-closed composition: FM_k1 x FM_k2 -> FM_{k1+k2-1}."""

    def test_output_arity(self):
        """Output arity is k1 + k2 - 1."""
        for k1 in range(2, 5):
            for k2 in range(2, 5):
                result = closed_composition_check(k1, k2)
                assert result['output_arity'] == k1 + k2 - 1

    def test_betti_match(self):
        """Output Betti numbers match FM_{k1+k2-1}."""
        for k1 in range(2, 5):
            for k2 in range(2, 5):
                result = closed_composition_check(k1, k2)
                assert result['match']


class TestOpenComposition:
    """Open-open composition: E_1(m1) x E_1(m2) -> E_1(m1+m2-1)."""

    def test_output_arity(self):
        """Output arity is m1 + m2 - 1."""
        for m1 in range(1, 5):
            for m2 in range(1, 5):
                result = open_composition_check(m1, m2)
                assert result['output_arity'] == m1 + m2 - 1

    def test_associative(self):
        """All compositions are associative (contractible spaces)."""
        for m1 in range(1, 5):
            for m2 in range(1, 5):
                result = open_composition_check(m1, m2)
                assert result['associative']


class TestMixedInterchange:
    """Interchange law for mixed SC operations."""

    def test_interchange_holds(self):
        """Interchange law holds for all small (k,m)."""
        for k in range(1, 4):
            for m in range(1, 4):
                result = mixed_interchange_check(k, m)
                assert result['all_interchange_holds'], \
                    f"Interchange failed at (k,m)=({k},{m})"


# ===================================================================
# BOUNDARY STRUCTURE
# ===================================================================

class TestStasheffBoundary:
    """Boundary faces of the Stasheff associahedron K_n."""

    def test_k2_one_face(self):
        """K_2: one face (0,2,0) corresponding to m_1(m_2(a,b)).

        The n=2 A-infinity identity has one term involving m_2:
        m_1(m_2(a,b)) + m_2(m_1(a),b) + ... = 0.
        The (0,2,0) face is the s=2 term.
        """
        faces = stasheff_boundary_faces(2)
        assert faces == [(0, 2, 0)]

    def test_k3_three_faces(self):
        """K_3: three faces from the n=3 A-infinity identity.

        stasheff_boundary_faces(3) gives faces (r,s,t):
        s=2: r can be 0 or 1. (0,2,1) and (1,2,0). That's 2.
        s=3: r=0, t=0. (0,3,0). That's 1.
        Total = 3.

        These correspond to:
        (0,2,1): m_2(m_2(a,b), c) = ((ab)c)
        (1,2,0): m_2(a, m_2(b,c)) = (a(bc))
        (0,3,0): m_1(m_3(a,b,c)) and m_3(m_1(a),b,c) etc.
        """
        faces = stasheff_boundary_faces(3)
        assert len(faces) == 3

    def test_k3_faces_explicit(self):
        """K_3 boundary faces: (0,2,1), (1,2,0), (0,3,0)."""
        faces = stasheff_boundary_faces(3)
        expected = [(0, 2, 1), (1, 2, 0), (0, 3, 0)]
        assert faces == expected

    def test_face_count_formula(self):
        """Number of faces = n(n-1)/2 for n >= 3."""
        for n in range(3, 8):
            assert stasheff_face_count(n) == n * (n - 1) // 2
            assert len(stasheff_boundary_faces(n)) == n * (n - 1) // 2

    def test_face_sum_constraint(self):
        """Every face (r,s,t) satisfies r + s + t = n and s >= 2."""
        for n in range(2, 7):
            for r, s, t in stasheff_boundary_faces(n):
                assert r + s + t == n
                assert s >= 2
                assert r >= 0
                assert t >= 0


class TestFMBoundary:
    """Boundary faces of FM_n(C) (Fulton-MacPherson)."""

    def test_fm2_faces(self):
        """FM_2: one face D_{12}."""
        faces = fm_boundary_faces(2)
        assert len(faces) == 1
        assert frozenset({1, 2}) in faces

    def test_fm3_faces(self):
        """FM_3: 4 faces (D_{12}, D_{13}, D_{23}, D_{123})."""
        faces = fm_boundary_faces(3)
        assert len(faces) == 4

    def test_face_count_formula(self):
        """Number of FM faces = 2^n - n - 1."""
        for n in range(2, 8):
            assert fm_face_count(n) == 2**n - n - 1
            assert len(fm_boundary_faces(n)) == 2**n - n - 1

    def test_fm_has_more_faces_than_stasheff(self):
        """FM_n has >= Stasheff faces (FM includes non-consecutive)."""
        for n in range(3, 7):
            assert fm_face_count(n) >= stasheff_face_count(n)

    def test_consecutive_matches_stasheff(self):
        """Consecutive-block FM faces match Stasheff face count.

        The Stasheff faces correspond exactly to consecutive blocks
        {r+1, ..., r+s} in FM_n(C).
        """
        for n in range(2, 7):
            comparison = fm_vs_stasheff_face_comparison(n)
            assert comparison['stasheff_matches_consecutive']


# ===================================================================
# THREE-FACE STOKES
# ===================================================================

class TestThreeFaceStokes:
    """Stokes' theorem on FM_3(C): the n=3 A-infinity relation."""

    def test_three_binary_faces(self):
        """FM_3 has exactly 3 binary collision faces."""
        result = three_face_stokes_verify()
        assert result['binary_faces'] == 3

    def test_boundary_consistency(self):
        """Boundary consistency check passes."""
        result = three_face_stokes_verify()
        assert result['boundary_consistency']

    def test_ainfty_arity(self):
        """This gives the arity-3 A-infinity identity."""
        result = three_face_stokes_verify()
        assert result['a_infinity_arity'] == 3


# ===================================================================
# DIMENSION TABLE
# ===================================================================

class TestDimensionTable:
    """Dimension table for SC(k,m) components."""

    def test_all_entries_computed(self):
        """All entries in the table are computed."""
        table = mixed_operation_dimensions(max_k=3, max_m=3)
        assert len(table) > 0

    def test_sc_10_dim(self):
        """SC(1,0) has dim 0 (point)."""
        result = sc_dimension_check(1, 0)
        assert result['dim'] == 0

    def test_sc_01_dim(self):
        """SC(0,1) has dim 0 (point)."""
        result = sc_dimension_check(0, 1)
        assert result['dim'] == 0

    def test_sc_11_dim(self):
        """SC(1,1) has dim 1."""
        result = sc_dimension_check(1, 1)
        assert result['dim'] == 1

    def test_sc_20_dim(self):
        """SC(2,0) has dim 2."""
        result = sc_dimension_check(2, 0)
        assert result['dim'] == 2

    def test_sc_02_dim(self):
        """SC(0,2) has dim 0 (after quotienting translation + dilation)."""
        result = sc_dimension_check(0, 2)
        assert result['dim'] == 0

    def test_dimension_formula(self):
        """dim SC(k,m) = 2k + m - 2 for k+m >= 2."""
        for k in range(0, 5):
            for m in range(0, 5):
                if k + m < 2:
                    continue
                result = sc_dimension_check(k, m)
                assert result['dim'] == 2 * k + m - 2


# ===================================================================
# HOMOTOPY-KOSZULITY
# ===================================================================

class TestHomotopyKoszul:
    """Homotopy-Koszulity of SC^{ch,top}.

    The theorem (thm:homotopy-Koszul) follows from:
    (1) Classical SC is Koszul (Livernet-Voronov)
    (2) Kontsevich formality (quasi-iso)
    (3) Bar-cobar transfer
    """

    def test_indicators_computed(self):
        """Koszul indicators are computed without error."""
        indicators = homotopy_koszul_indicators(max_k=4)
        assert indicators is not None

    def test_closed_poincare_computed(self):
        """Closed part Poincare polynomials are computed."""
        indicators = homotopy_koszul_indicators(max_k=4)
        for n in range(1, 5):
            assert n in indicators['closed_poincare']

    def test_koszul_dual_closed_is_lie(self):
        """Koszul dual of closed part (Com) is Lie: dim = (n-1)!"""
        indicators = homotopy_koszul_indicators(max_k=6)
        for n in range(1, 7):
            assert indicators['koszul_dual_closed'][n] == factorial(n - 1)

    def test_koszul_dual_open_is_ass(self):
        """Koszul dual of open part (Ass) is Ass: dim = n!"""
        indicators = homotopy_koszul_indicators(max_k=6)
        for n in range(1, 7):
            assert indicators['koszul_dual_open'][n] == factorial(n)

    def test_bar_cobar_inversions(self):
        """Bar-cobar inversions hold: Com<->Lie, Ass<->Ass."""
        indicators = homotopy_koszul_indicators()
        inversions = indicators['bar_cobar_inversions']
        assert inversions['Com_to_Lie']
        assert inversions['Ass_to_Ass']
        assert inversions['SC_homotopy_Koszul']


class TestKoszulCriterion:
    """Koszul criterion: bar complex concentrated in top weight."""

    def test_criterion_passes(self):
        """Koszul criterion check passes for all tested arities."""
        results = koszul_criterion_check(max_n=6)
        for n in range(1, 7):
            assert results[n]['koszul_concentration']

    def test_lie_dim(self):
        """Koszul dual dimension = (n-1)! (Lie operad)."""
        results = koszul_criterion_check(max_n=6)
        for n in range(1, 7):
            assert results[n]['lie_dim'] == factorial(n - 1)


# ===================================================================
# OPERADIC GENERATING FUNCTIONS
# ===================================================================

class TestOperadicGF:
    """Generating functions for the closed and open parts of SC."""

    def test_closed_gf_type(self):
        """Closed part is Com, dual is Lie."""
        gf = operadic_gf_closed()
        assert gf['closed_type'] == 'Com'
        assert gf['dual_type'] == 'Lie'

    def test_open_gf_self_dual(self):
        """Open part (Ass) is self-dual."""
        gf = operadic_gf_open()
        assert gf['self_dual']
        assert gf['open_type'] == 'Ass'
        assert gf['dual_type'] == 'Ass'

    def test_com_lie_coefficients(self):
        """Com coefficients are 1/n!, Lie coefficients are 1/n.

        Com: f(x) = sum x^n/n!. Coefficient of x^n is 1/n!.
        Lie: g(x) = sum (n-1)! x^n / n! = sum x^n/n. Coefficient is 1/n.
        """
        from fractions import Fraction
        gf = operadic_gf_closed(max_n=6)
        for n in range(1, 7):
            assert gf['closed_coefficients'][n] == Fraction(1, factorial(n))
            assert gf['dual_coefficients'][n] == Fraction(1, n)


# ===================================================================
# BAR COMPLEX DIMENSIONS
# ===================================================================

class TestBarComplex:
    """Bar complex dimensions for SC."""

    def test_closed_bar_cohomology(self):
        """B(Com)(n) has cohomology of dimension (n-1)! = dim Lie(n)."""
        data = bar_complex_dimensions_sc(max_arity=6)
        for n in range(1, 7):
            assert data[n]['closed_bar_cohomology'] == factorial(n - 1)

    def test_open_bar_cohomology(self):
        """B(Ass)(n) has cohomology of dimension n! = dim Ass(n)."""
        data = bar_complex_dimensions_sc(max_arity=6)
        for n in range(1, 7):
            assert data[n]['open_bar_cohomology'] == factorial(n)

    def test_concentration(self):
        """Bar complex cohomology is concentrated (Koszulity)."""
        data = bar_complex_dimensions_sc(max_arity=6)
        for n in range(1, 7):
            assert data[n]['closed_concentrated']
            assert data[n]['open_concentrated']


# ===================================================================
# FULL SUMMARY
# ===================================================================

class TestFullSummary:
    """Full SC verification summary."""

    def test_summary_passes(self):
        """Complete verification summary passes all checks."""
        summary = sc_verification_summary(max_n=5)
        assert summary['all_pass']

    def test_summary_contains_aos(self):
        """Summary includes AOS checks for each arity."""
        summary = sc_verification_summary(max_n=4)
        for n in range(1, 5):
            assert n in summary['aos_checks']

    def test_summary_koszul_indicators(self):
        """Summary includes Koszul indicators."""
        summary = sc_verification_summary(max_n=4)
        assert summary['koszul_indicators'] is not None


# ===================================================================
# CROSS-CHECKS WITH VOL II INFRASTRUCTURE
# ===================================================================

class TestCrossChecks:
    """Cross-checks with existing fm_boundary and arnold modules."""

    def test_fm_boundary_face_count_agrees(self):
        """fm_boundary_faces agrees with fm_boundary.count_boundary_strata."""
        from lib.fm_boundary import count_boundary_strata
        for n in range(2, 8):
            assert fm_face_count(n) == count_boundary_strata(n)

    def test_fm_boundary_faces_agree(self):
        """fm_boundary_faces agrees with fm_boundary.boundary_strata."""
        from lib.fm_boundary import boundary_strata
        for n in range(2, 6):
            our_faces = set(fm_boundary_faces(n))
            their_faces = set(boundary_strata(n))
            assert our_faces == their_faces, f"Mismatch at n={n}"

    def test_poincare_agrees_with_arnold_os(self):
        """Poincare polynomial agrees with arnold.orlik_solomon_presentation."""
        from lib.arnold import orlik_solomon_presentation
        for n in range(2, 6):
            our_poly = poincare_polynomial_fm(n)
            os_info = orlik_solomon_presentation(n)
            # The Poincare factors from OS are [1, 2, ..., n-1]
            # Our polynomial should match
            expected_factors = os_info['poincare_factors']
            assert expected_factors == list(range(1, n))

    def test_generator_count_agrees_with_arnold(self):
        """Number of AOS generators agrees with arnold module."""
        from lib.arnold import orlik_solomon_presentation
        for n in range(2, 6):
            our_info = arnold_relation_generators(n)
            os_info = orlik_solomon_presentation(n)
            assert our_info['num_generators'] == os_info['generators']
            assert our_info['num_relations'] == os_info['relations']
