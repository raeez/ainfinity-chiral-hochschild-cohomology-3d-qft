"""
Tests for worked examples: free multiplet, abelian CS, Virasoro, su(2), LG cubic.

Goal: EXPOSE errors in the paper's computations, not just confirm.
Strategy: compare against independently known results from the literature.

Test tiers:
  Tier 1 (structural): self-certifying identities (d²=0, Jacobi, etc.)
  Tier 2 (published): known values from Kac, BPZ, DSK, etc.
  Tier 3 (cross-check): two independent code paths agree
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy import Symbol, Rational, simplify, expand, S, symbols


# ===================================================================
# FREE MULTIPLET
# ===================================================================

class TestFreeMultiplet:
    """Free chiral multiplet: all m_{k>=3} = 0."""

    def test_m2_is_product(self):
        """m_2(a,b) = a*b with no spectral parameter dependence."""
        from lib.examples.free_multiplet import m2
        a, b = Symbol('a'), Symbol('b')
        assert m2(a, b) == a * b

    def test_higher_mk_vanish(self):
        """m_k = 0 for k >= 3."""
        from lib.examples.free_multiplet import m_k
        a, b, c, d = symbols('a b c d')
        assert m_k(3, a, b, c) == 0
        assert m_k(4, a, b, c, d) == 0
        assert m_k(5, a, b, c, d, a) == 0

    def test_ainfty_identity_arity3_odd_elements(self):
        """n=3 identity for degree-1 elements (BV anti-fields).

        In the paper's convention (|m_k|=1-k), for degree-1 elements
        with Q=0 and m_{k>=3}=0:
        epsilon(1,2) = (2-1)*1 = 1, so sign = -1.
        => m_2(m_2(a,b),c) - m_2(a,m_2(b,c)) = 0 (associativity).
        """
        from lib.examples.free_multiplet import m2, m_k
        from lib.ainfty import verify_ainfty_identity
        a, b, c = symbols('a b c')
        # m_1 = Q = 0 for Q-closed elements (cohomology representatives)
        ops = {
            1: lambda x: S.Zero,
            2: lambda x, y: m2(x, y),
            3: lambda x, y, z: m_k(3, x, y, z),
        }
        result = verify_ainfty_identity(ops, [(a, 1), (b, 1), (c, 1)], 3)
        assert simplify(result) == 0

    def test_lambda_bracket_vanishes(self):
        """Lambda-bracket should be zero for free theory."""
        from lib.examples.free_multiplet import lambda_bracket
        J = Symbol('J')
        lam = Symbol('lambda')
        assert lambda_bracket(J, J, lam) == 0


# ===================================================================
# ABELIAN CS
# ===================================================================

class TestAbelianCS:
    """Abelian Chern-Simons boundary current algebra."""

    def test_lambda_bracket_current(self):
        """{J_lambda J} = k*lambda."""
        from lib.examples.abelian_cs import current_lambda_bracket
        J = Symbol('J')
        lam = Symbol('lambda')
        k = Symbol('k')
        result = current_lambda_bracket(J, J, lam, k)
        assert simplify(result - k * lam) == 0

    def test_skew_symmetry(self):
        """Skew-symmetry for abelian current algebra."""
        from lib.examples.abelian_cs import check_skew_symmetry_current
        lam = Symbol('lambda')
        k = Symbol('k')
        assert check_skew_symmetry_current(lam, k) == 0

    def test_jacobi(self):
        """Jacobi for abelian current algebra (trivial)."""
        from lib.examples.abelian_cs import check_jacobi_current
        lam, mu = symbols('lambda mu')
        k = Symbol('k')
        assert check_jacobi_current(lam, mu, k) == 0

    def test_r_matrix_form(self):
        """R-matrix for abelian CS: r(z) = hbar * k*q1*q2 / z."""
        from lib.examples.abelian_cs import abelian_r_matrix
        z = Symbol('z')
        hbar = Symbol('hbar')
        k = Symbol('k')
        q1, q2 = symbols('q1 q2')
        r = abelian_r_matrix(z, hbar, k)
        assert simplify(r - hbar * k * q1 * q2 / z) == 0

    def test_r_matrix_level_zero_collapse(self):
        """At level zero the Heisenberg collision residue vanishes."""
        from lib.examples.abelian_cs import abelian_r_matrix
        z = Symbol('z')
        hbar = Symbol('hbar')
        assert abelian_r_matrix(z, hbar, 0) == 0

    def test_yang_baxter_abelian(self):
        """YBE trivially satisfied for abelian case."""
        from lib.examples.abelian_cs import check_yang_baxter_abelian
        z1, z2 = symbols('z1 z2')
        hbar = Symbol('hbar')
        assert check_yang_baxter_abelian(z1, z2, hbar) == 0


# ===================================================================
# VIRASORO — REAL COMPUTATION
# ===================================================================

class TestVirasoro:
    """Virasoro algebra lambda-bracket checks — ALL GENUINE."""

    def test_ope_to_lambda_coefficient(self):
        """c/2 in OPE fourth pole -> c/12 in lambda-bracket cubic term.

        This conversion factor of 1/3! = 1/6 is easy to get wrong.
        Tier 2 (published): Kac, Vertex Algebras for Beginners.
        """
        from lib.examples.virasoro import ope_to_lambda_bracket_coefficient_check
        c = Symbol('c')
        assert ope_to_lambda_bracket_coefficient_check(c) == 0

    def test_lambda_bracket_form(self):
        """Check {T_lambda T} = dT + 2T*lambda + (c/12)*lambda^3.

        Tier 2 (published): matches standard Virasoro lambda-bracket.
        """
        from lib.examples.virasoro import virasoro_lambda_bracket, T, dT
        lam = Symbol('lambda')
        c = Symbol('c')
        partial = Symbol('partial')

        result = virasoro_lambda_bracket(T, T, lam, partial, c)
        expected = dT + 2 * T * lam + Rational(1, 12) * c * lam**3
        assert simplify(result - expected) == 0

    def test_skew_symmetry(self):
        """Skew-symmetry via Borcherds formula.

        Tier 1 (structural): self-certifying identity.
        """
        from lib.examples.virasoro import check_virasoro_skew_symmetry
        lam = Symbol('lambda')
        partial = Symbol('partial')
        c = Symbol('c')
        assert check_virasoro_skew_symmetry(lam, partial, c) == 0

    def test_jacobi_real(self):
        """REAL Jacobi identity check — NOT `return S.Zero`.

        Tier 1 (structural): self-certifying identity.
        This is the single highest-value test in the codebase.

        Verifies: {T_λ {T_μ T}} − {T_μ {T_λ T}} = {{T_λ T}_{λ+μ} T}
        for symbolic c (arbitrary central charge).
        """
        from lib.examples.virasoro import check_virasoro_jacobi
        lam, mu = symbols('lambda mu')
        c = Symbol('c')
        result = check_virasoro_jacobi(lam, mu, c)
        assert result == 0, f"Virasoro Jacobi FAILED: {result}"

    def test_jacobi_components(self):
        """Component-wise Jacobi check: extract T, dT, d2T coefficients.

        More robust than relying on simplify() for the full expression.
        """
        from lib.examples.virasoro import check_virasoro_jacobi_via_components
        lam, mu = symbols('lambda mu')
        c = Symbol('c')
        results = check_virasoro_jacobi_via_components(lam, mu, c)
        for field_name, coeff in results.items():
            assert coeff == 0, f"Jacobi component '{field_name}' nonzero: {coeff}"

    def test_jacobi_c1(self):
        """Jacobi at c=1 (free boson). Tier 2 (published)."""
        from lib.examples.virasoro import check_virasoro_jacobi
        lam, mu = symbols('lambda mu')
        result = check_virasoro_jacobi(lam, mu, S.One)
        assert result == 0

    def test_jacobi_c26(self):
        """Jacobi at c=26 (bosonic string). Tier 2 (published)."""
        from lib.examples.virasoro import check_virasoro_jacobi
        lam, mu = symbols('lambda mu')
        result = check_virasoro_jacobi(lam, mu, S(26))
        assert result == 0

    def test_jacobi_c_minus2(self):
        """Jacobi at c=-2 (symplectic fermion). Tier 2."""
        from lib.examples.virasoro import check_virasoro_jacobi
        lam, mu = symbols('lambda mu')
        result = check_virasoro_jacobi(lam, mu, S(-2))
        assert result == 0

    def test_sesquilinearity_right(self):
        """Right sesquilinearity: {T_λ ∂T} = (λ+∂){T_λ T}.

        Tier 1 (structural).
        """
        from lib.examples.virasoro import check_virasoro_leibniz
        lam = Symbol('lambda')
        c = Symbol('c')
        assert check_virasoro_leibniz(lam, c) == 0

    def test_central_charge_conventions(self):
        """Cross-check: c is a free parameter in the abstract PVA.

        Tier 2: the Virasoro PVA is defined for all c ∈ C.
        """
        # Verified by the fact that all checks above work for symbolic c
        pass


# ===================================================================
# SU(2) CURRENT ALGEBRA — FIRST NONTRIVIAL JACOBI
# ===================================================================

class TestSU2CurrentAlgebra:
    """su(2)_k current algebra from nonabelian CS.

    FIRST genuinely nontrivial Jacobi check in the codebase.
    """

    def test_bracket_structure_constants(self):
        """Verify {J^1_λ J^2} = J^3 + 0 (no level contribution for a≠b).

        Tier 2 (published): ε^{123} = 1.
        """
        from lib.examples.nonabelian_cs import su2_lambda_bracket, J
        lam = Symbol('lambda')
        k = Symbol('k')
        result = su2_lambda_bracket(1, 2, lam, k)
        assert simplify(result - J[3]) == 0

    def test_bracket_level_term(self):
        """{J^1_λ J^1} = k·λ (diagonal: level contribution only).

        Tier 2 (published).
        """
        from lib.examples.nonabelian_cs import su2_lambda_bracket
        lam = Symbol('lambda')
        k = Symbol('k')
        result = su2_lambda_bracket(1, 1, lam, k)
        assert simplify(result - k * lam) == 0

    def test_bracket_antisymmetry_structure(self):
        """{J^a_λ J^b} = -{J^b_λ J^a} at λ=0 (antisymmetry of structure constants)."""
        from lib.examples.nonabelian_cs import su2_lambda_bracket
        lam = Symbol('lambda')
        k = Symbol('k')
        for a in [1, 2, 3]:
            for b in [1, 2, 3]:
                bracket_ab = su2_lambda_bracket(a, b, S.Zero, k)
                bracket_ba = su2_lambda_bracket(b, a, S.Zero, k)
                assert simplify(bracket_ab + bracket_ba) == 0, \
                    f"Antisymmetry failed for ({a},{b})"

    def test_skew_symmetry_all(self):
        """Skew-symmetry via Borcherds for all 9 pairs.

        Tier 1 (structural).
        """
        from lib.examples.nonabelian_cs import check_su2_skew_symmetry
        lam = Symbol('lambda')
        k = Symbol('k')
        for a in [1, 2, 3]:
            for b in [1, 2, 3]:
                result = check_su2_skew_symmetry(a, b, lam, k)
                assert result == 0, f"Skew-symmetry failed for ({a},{b}): {result}"

    def test_jacobi_123(self):
        """Jacobi for (a,b,c) = (1,2,3): the most informative triple.

        Tier 1 (structural): self-certifying identity.
        """
        from lib.examples.nonabelian_cs import check_su2_jacobi
        lam, mu = symbols('lambda mu')
        k = Symbol('k')
        result = check_su2_jacobi(1, 2, 3, lam, mu, k)
        assert result == 0, f"Jacobi (1,2,3) FAILED: {result}"

    def test_jacobi_112(self):
        """Jacobi for (1,1,2): has both structure and level contributions."""
        from lib.examples.nonabelian_cs import check_su2_jacobi
        lam, mu = symbols('lambda mu')
        k = Symbol('k')
        result = check_su2_jacobi(1, 1, 2, lam, mu, k)
        assert result == 0, f"Jacobi (1,1,2) FAILED: {result}"

    def test_jacobi_all_27(self):
        """Jacobi for ALL 27 triples (a,b,c) ∈ {1,2,3}³.

        Tier 1 (structural): exhaustive check.
        This is the STRONGEST test: all 27 Jacobi identities for su(2).
        """
        from lib.examples.nonabelian_cs import check_su2_jacobi_all
        lam, mu = symbols('lambda mu')
        k = Symbol('k')
        results = check_su2_jacobi_all(lam, mu, k)
        for (a, b, c), val in results.items():
            assert val == 0, f"Jacobi ({a},{b},{c}) FAILED: {val}"

    def test_sugawara_central_charge(self):
        """Sugawara c = 3k/(k+2) for sl_2.

        Tier 2 (published): Kac, Vertex Algebras for Beginners.
        CRITICAL: Sugawara UNDEFINED at k = -2 (critical level).
        """
        from lib.examples.nonabelian_cs import sugawara_central_charge
        # k=1: c = 3/3 = 1
        assert sugawara_central_charge(1) == 1
        # k=2: c = 6/4 = 3/2
        assert sugawara_central_charge(2) == Rational(3, 2)
        # k=10: c = 30/12 = 5/2
        assert sugawara_central_charge(10) == Rational(5, 2)


# ===================================================================
# LG CUBIC
# ===================================================================

class TestLGCubic:
    """Landau-Ginzburg model with W = g*phi^3/3."""

    def test_m1_Q_squared_phi(self):
        """Q²(phi) = 0. Tier 1 (structural: d²=0)."""
        from lib.examples.lg_cubic import check_Q_squared_lg
        phi = Symbol('phi')
        assert check_Q_squared_lg(phi, 'phi') == 0

    def test_m1_Q_squared_psi(self):
        """Q²(psi) = 0. Tier 1."""
        from lib.examples.lg_cubic import check_Q_squared_lg
        psi = Symbol('psi')
        assert check_Q_squared_lg(psi, 'psi') == 0

    def test_m1_is_free(self):
        """Q(phi) = phi_value (represents psi), Q(psi) = 0.

        The LG differential at phi=0 is identical to the free theory.
        """
        from lib.examples.lg_cubic import m1_lg
        phi = Symbol('phi')
        assert m1_lg(phi, 'phi') == phi  # Q maps phi -> psi (value preserved)
        assert m1_lg(phi, 'psi') == 0    # Q(psi) = 0

    def test_m2_free_part_commutative(self):
        """Free part of m_2 is commutative. Tier 1."""
        from lib.examples.lg_cubic import m2_lg_free_part
        a, b = symbols('a b')
        assert simplify(m2_lg_free_part(a, b) - m2_lg_free_part(b, a)) == 0

    def test_m3_proportional_to_g(self):
        """m_3 = 0 when g = 0 (free theory limit). Tier 1."""
        from lib.examples.lg_cubic import m3_lg_vertex
        a, b, c = symbols('a b c')
        assert m3_lg_vertex(a, b, c, 0) == 0

    def test_m3_nonzero(self):
        """m_3 ≠ 0 for g ≠ 0 (the key feature of LG). Tier 1."""
        from lib.examples.lg_cubic import m3_lg_vertex
        a, b, c, g = symbols('a b c g')
        result = m3_lg_vertex(a, b, c, g)
        assert result != 0
        assert simplify(result - 2 * g * a * b * c) == 0

    def test_m3_symmetric(self):
        """m_3 should be graded-symmetric (for degree-0 inputs).

        For three degree-0 fields, the Koszul signs are all +1,
        so m_3(a,b,c) should be symmetric in a,b,c.
        """
        from lib.examples.lg_cubic import m3_lg_vertex
        a, b, c, g = symbols('a b c g')
        assert simplify(m3_lg_vertex(a, b, c, g) - m3_lg_vertex(b, a, c, g)) == 0
        assert simplify(m3_lg_vertex(a, b, c, g) - m3_lg_vertex(a, c, b, g)) == 0
        assert simplify(m3_lg_vertex(a, b, c, g) - m3_lg_vertex(c, b, a, g)) == 0

    def test_m4_vanishes(self):
        """m_4 = 0 by the general formula. Tier 1."""
        from lib.examples.lg_cubic import m_k_lg
        a, b, c, d, g = symbols('a b c d g')
        assert m_k_lg(4, (a, b, c, d), g) == 0

    def test_m5_vanishes(self):
        """m_5 = 0. Tier 1."""
        from lib.examples.lg_cubic import m_k_lg
        a, b, c, d, e, g = symbols('a b c d e g')
        assert m_k_lg(5, (a, b, c, d, e), g) == 0

    def test_truncation_degree_counting_k3(self):
        """k=3: V=1, E=0, m_3 nonzero. Tier 1."""
        from lib.examples.lg_cubic import check_truncation_degree_counting
        info = check_truncation_degree_counting(3)
        assert info['vertices'] == 1
        assert info['internal_edges'] == 0
        assert info['vanishes'] is False

    def test_truncation_degree_counting_k4(self):
        """k=4: V=2, E=1, and W^(4)=0. Tier 1."""
        from lib.examples.lg_cubic import check_truncation_degree_counting
        info = check_truncation_degree_counting(4)
        assert info['vertices'] == 2
        assert info['internal_edges'] == 1
        assert info['vanishes'] is True
        assert info['taylor_coefficient_zero'] is True

    def test_truncation_degree_counting_k5(self):
        """k=5: V=3, E=2, and W^(5)=0. Tier 1."""
        from lib.examples.lg_cubic import check_truncation_degree_counting
        info = check_truncation_degree_counting(5)
        assert info['vertices'] == 3
        assert info['internal_edges'] == 2
        assert info['vanishes'] is True
        assert info['taylor_coefficient_zero'] is True

    def test_truncation_degree_counting_k10(self):
        """k=10: should also vanish. Tier 1."""
        from lib.examples.lg_cubic import check_truncation_degree_counting
        info = check_truncation_degree_counting(10)
        assert info['vanishes'] is True

    def test_m_k_reduces_to_free_at_g0(self):
        """At g=0, LG model = free model: m_2 = product, m_{k>=3} = 0."""
        from lib.examples.lg_cubic import m_k_lg
        a, b, c = symbols('a b c')
        # m_2 at g=0 is still the product
        assert m_k_lg(2, (a, b), 0) == a * b
        # m_3 at g=0 vanishes
        assert m_k_lg(3, (a, b, c), 0) == 0


# ===================================================================
# LAURENT SERIES
# ===================================================================

class TestLaurentSeries:
    """Tests for the spectral parameter algebra."""

    def test_regular_singular_decomposition(self):
        """Decompose m_2(a,b;lambda) into regular and singular parts."""
        from lib.spectral import LaurentSeries, reg_sing_decompose
        lam = Symbol('lambda')
        a, b, c = symbols('a b c')
        series = LaurentSeries({-1: a, 0: b, 1: c}, lam)

        reg, sing = reg_sing_decompose(series)
        assert sing.coeff(-1) == a
        assert sing.coeff(0) == 0
        assert reg.coeff(0) == b
        assert reg.coeff(1) == c
        assert reg.coeff(-1) == 0

    def test_addition(self):
        """Laurent series addition."""
        from lib.spectral import LaurentSeries
        lam = Symbol('lambda')
        s1 = LaurentSeries({-1: 1, 0: 2}, lam)
        s2 = LaurentSeries({-1: 3, 1: 4}, lam)
        s3 = s1 + s2
        assert s3.coeff(-1) == 4
        assert s3.coeff(0) == 2
        assert s3.coeff(1) == 4

    def test_multiplication_pole_times_regular(self):
        """(1/lambda) * lambda = 1."""
        from lib.spectral import LaurentSeries
        lam = Symbol('lambda')
        pole = LaurentSeries({-1: 1}, lam)
        reg = LaurentSeries({1: 1}, lam)
        product = pole * reg
        assert product.coeff(0) == 1
        assert product.coeff(-1) == 0
        assert product.coeff(1) == 0

    def test_zero_series(self):
        """Zero series is zero."""
        from lib.spectral import LaurentSeries
        lam = Symbol('lambda')
        z = LaurentSeries({}, lam)
        assert z.is_zero()


# ===================================================================
# LG TRUNCATION (Taylor coefficient)
# ===================================================================

class TestLGTruncation:
    """Tests related to the LG cubic truncation claim."""

    def test_degree_counting_m4(self):
        """Check the primitive Taylor obstruction for m_4."""
        from lib.examples.lg_cubic import check_truncation_degree_counting
        info = check_truncation_degree_counting(4)
        assert info['vertices'] == 2
        assert info['internal_edges'] == 1
        assert info['primitive_derivative_order'] == 4
        assert info['taylor_coefficient_zero'] is True
        assert info['vanishes'] is True
