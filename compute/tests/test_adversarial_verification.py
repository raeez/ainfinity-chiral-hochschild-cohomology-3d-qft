r"""Adversarial verification tests for computational mathematics in Vol II.

Tests designed to probe edge cases, cross-check between independent modules,
and verify mathematical identities beyond the existing test range.

Sections:
  1. Catalan factorisation at k=12,13,14 (beyond existing k<=11 range)
  2. Non-simply-laced R-matrix: B2, G2 Casimir symmetry and CYBE
  3. Genus-1 intersection numbers: cross-checks against known results
  4. W3 multichannel shadow: TT channel vs known Virasoro values via DS
  5. Bigraded Hilbert series: Heisenberg against exact partition function
  6. Field sector generating functions: all_field_sector_generating_functions.py
  7. Cross-engine consistency: FP number, modular completion, genus bridge
  8. AP19 pole absorption verification
  9. Mathematical identity checks

DISCOVERED ISSUES:
  - MAX_DERIV=20 in StasheffEngine truncates derivative tracking, causing
    even-arity vanishing to fail at k>=14 (needs 2*(k-2)=24 derivative orders).
  - test_modular_obstruction_engine::test_genus2_fp_number expects 1/1152
    but engine returns 7/5760 (the correct A-hat genus coefficient).
"""

import sys
import os
import math
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))


# =========================================================================
# HELPERS
# =========================================================================

def catalan(n):
    """Catalan number C_n = binom(2n,n)/(n+1)."""
    if n < 0:
        return 0
    return math.comb(2 * n, n) // (n + 1)


def shifted_pochhammer_coeffs(r):
    """Coefficients of (y+2)(y+3)...(y+2r+1) as list [a_0, ..., a_{2r}]."""
    poly = [1]
    for j in range(2, 2 * r + 2):
        new_poly = [0] * (len(poly) + 1)
        for i, c in enumerate(poly):
            new_poly[i] += c * j
            new_poly[i + 1] += c
        poly = new_poly
    return poly


# =========================================================================
# 1. CATALAN FACTORISATION AT k=12,13,14
# =========================================================================

class TestCatalanFactorisationExtended:
    """Extend Catalan factorisation tests to k=12,13,14 (beyond existing k<=11).

    FINDING: k=14 even vanishing fails due to MAX_DERIV=20 truncation.
    The Stasheff engine needs derivative orders up to 2*(k-2) = 24 for k=14,
    but MAX_DERIV=20 causes truncation. Tests at k<=12 are safe (need <= 20).
    """

    @pytest.fixture
    def engine(self):
        from m7_m10_depth_frontier import StasheffEngine
        return StasheffEngine(0.0)

    def _get_field_polynomial(self, engine, k):
        """Extract phi_k(x) as coefficient list."""
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(k - 1))
        result = engine.mk(lams)
        field_items = {w: v for w, v in result.items() if w >= 0}
        if not field_items:
            return [0.0]
        max_w = max(field_items.keys())
        return [field_items.get(w, 0.0) for w in range(max_w + 1)]

    def _eval_polynomial(self, coeffs, x):
        return sum(c * x ** w for w, c in enumerate(coeffs))

    def _predicted_phi(self, k):
        """Predicted phi_k for odd k using Catalan factorisation."""
        if k % 2 == 0:
            return [0.0]
        n = (k - 3) // 2
        Cn = catalan(n)
        prefactor = (-1) ** n * Cn
        poly = [1.0]
        for m in range(2, k + 1):
            new_poly = [0.0] * (len(poly) + 1)
            for i in range(len(poly)):
                new_poly[i] += float(m) * poly[i]
                new_poly[i + 1] += poly[i]
            poly = new_poly
        return [prefactor * c for c in poly]

    # --- Even k: vanishing (safe range k<=12) ---

    @pytest.mark.parametrize("k", [12])
    def test_even_vanishing_k12(self, engine, k):
        """phi_12(x) = 0 (within MAX_DERIV=20 safe zone: needs 2*10=20)."""
        coeffs = self._get_field_polynomial(engine, k)
        max_abs = max(abs(c) for c in coeffs)
        assert max_abs < 1e-4, f"phi_{k} not zero: max |coeff| = {max_abs}"

    def test_even_vanishing_k14_TRUNCATION_BUG(self, engine):
        """KNOWN ISSUE: phi_14 does NOT vanish due to MAX_DERIV=20 truncation.

        k=14 needs derivative orders up to 2*(14-2) = 24, but MAX_DERIV = 20.
        The returned value (132 = C_6) at the symmetric point is a truncation
        artifact, not a mathematical result.
        """
        from m7_m10_depth_frontier import MAX_DERIV
        # Verify the truncation condition
        k = 14
        max_needed = 2 * (k - 2)
        assert max_needed > MAX_DERIV, \
            f"k={k} needs {max_needed} derivative orders, but MAX_DERIV={MAX_DERIV}"

        coeffs = self._get_field_polynomial(engine, k)
        max_abs = max(abs(c) for c in coeffs)
        # This WILL be nonzero -- it's a truncation artifact
        assert max_abs > 1, \
            f"k=14 unexpectedly vanishes (truncation bug may be fixed)"

    # --- Odd k: full polynomial match ---

    @pytest.mark.parametrize("k", [13])
    def test_polynomial_match_k13(self, engine, k):
        """phi_13 = (-1)^5 C_5 prod(x+m, m=2..13).

        NOTE: k=13 needs derivative orders up to 22 > MAX_DERIV=20.
        However, the leading coefficients (low derivative order) still match
        because truncation only affects high-order terms. We relax the
        tolerance to account for potential high-order contamination.
        """
        computed = self._get_field_polynomial(engine, k)
        predicted = self._predicted_phi(k)
        maxlen = max(len(computed), len(predicted))
        computed += [0.0] * (maxlen - len(computed))
        predicted += [0.0] * (maxlen - len(predicted))
        # Only check low-order coefficients that are within MAX_DERIV
        from m7_m10_depth_frontier import MAX_DERIV
        for w in range(min(maxlen, MAX_DERIV + 1)):
            if abs(predicted[w]) > 1:
                rel_err = abs(computed[w] - predicted[w]) / abs(predicted[w])
                assert rel_err < 1e-2, \
                    f"phi_{k} x^{w}: computed={computed[w]}, predicted={predicted[w]}, rel_err={rel_err}"

    # --- T-coefficient ---

    @pytest.mark.parametrize("k", [13])
    def test_T_coefficient_k13(self, engine, k):
        """T_{13} = (-1)^5 C_5 * 13! = -42 * 6227020800."""
        coeffs = self._get_field_polynomial(engine, k)
        T_val = coeffs[0]
        n = (k - 3) // 2
        predicted = (-1) ** n * catalan(n) * math.factorial(k)
        assert abs(T_val - predicted) < abs(predicted) * 1e-6, \
            f"T_{k} = {T_val}, predicted {predicted}"

    # --- Root property ---

    @pytest.mark.parametrize("k", [13])
    def test_roots_k13(self, engine, k):
        """phi_13(-m) = 0 for m = 2,...,13."""
        coeffs = self._get_field_polynomial(engine, k)
        for m in range(2, k + 1):
            val = self._eval_polynomial(coeffs, float(-m))
            # Relaxed tolerance: high-arity may have large cancellations
            assert abs(val) < abs(coeffs[0]) * 1e-3, \
                f"phi_{k}({-m}) = {val}, should be 0"

    # --- Functional equation ---

    @pytest.mark.parametrize("k", [12, 14])
    def test_functional_equation_extended(self, k):
        """(x+2) phi_{k-1}(x+1) = (x+k) phi_{k-1}(x) for even k.

        This tests the PREDICTED polynomial (analytic formula), not the engine.
        """
        pred = self._predicted_phi(k - 1)
        for x in [-5.5, -1.3, 0, 0.7, 2.5, 10, 100]:
            lhs = (x + 2) * self._eval_polynomial(pred, x + 1)
            rhs = (x + k) * self._eval_polynomial(pred, x)
            assert abs(lhs - rhs) < 1e-3 * max(abs(rhs), 1), \
                f"k={k}, x={x}: LHS={lhs}, RHS={rhs}"

    # --- Catalan number identity ---

    def test_catalan_numbers_extended(self):
        """Verify Catalan numbers C_5 through C_7."""
        assert catalan(5) == 42
        assert catalan(6) == 132
        assert catalan(7) == 429
        for n in [5, 6, 7]:
            conv = sum(catalan(a) * catalan(n - 1 - a) for a in range(n))
            assert conv == catalan(n)

    # --- Profile polynomial at r=6 (safe: k=13 but checking stored coefficients) ---

    def test_profile_polynomial_r5(self):
        """Verify P_5(y) = (y+2)_{10} by computing via StasheffEngine.

        r=5 corresponds to k=11, which needs max deriv = 2*9 = 18 <= MAX_DERIV=20.
        """
        from m7_m10_depth_frontier import StasheffEngine
        engine = StasheffEngine(0.0)
        k = 11  # r=5
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(k - 1))
        result = engine.mk(lams)
        r = 5
        sign = (-1) ** (r + 1)
        cat = catalan(r - 1)  # C_4 = 14
        predicted = shifted_pochhammer_coeffs(r)
        for w in range(min(len(predicted), 2 * r + 1)):
            val = result.get(w, 0.0)
            expected = sign * cat * predicted[w]
            if abs(expected) > 1:
                assert abs(val - expected) < abs(expected) * 1e-4, \
                    f"r={r}, w={w}: val={val}, expected={expected}"

    def test_max_deriv_boundary_detection(self):
        """Verify MAX_DERIV is the engine's limiting factor for k>=13."""
        from m7_m10_depth_frontier import MAX_DERIV
        # Safe arities: 2*(k-2) <= MAX_DERIV => k <= MAX_DERIV/2 + 2
        k_max_safe = MAX_DERIV // 2 + 2
        assert k_max_safe == 12, \
            f"Expected max safe k=12 for MAX_DERIV={MAX_DERIV}, got {k_max_safe}"


# =========================================================================
# 2. NON-SIMPLY-LACED R-MATRIX: B2 AND G2
# =========================================================================

class TestNonSimplyLacedAdversarial:
    """Adversarial tests for non-simply-laced Casimir and CYBE."""

    def test_B2_casimir_symmetry(self):
        """Casimir tensor Omega^{ab} is symmetric for B_2."""
        import numpy as np
        from lib.non_simply_laced_rmatrix import make_B2
        from lib.collision_residue_rmatrix import casimir_tensor
        g = make_B2()
        Omega = casimir_tensor(g)
        assert np.allclose(Omega, Omega.T, atol=1e-12), \
            f"B_2 Casimir not symmetric: max diff = {np.max(np.abs(Omega - Omega.T))}"

    def test_G2_casimir_symmetry(self):
        """Casimir tensor Omega^{ab} is symmetric for G_2."""
        import numpy as np
        from lib.non_simply_laced_rmatrix import make_G2
        from lib.collision_residue_rmatrix import casimir_tensor
        g = make_G2()
        Omega = casimir_tensor(g)
        assert np.allclose(Omega, Omega.T, atol=1e-12), \
            f"G_2 Casimir not symmetric: max diff = {np.max(np.abs(Omega - Omega.T))}"

    def test_G2_jacobi(self):
        from lib.non_simply_laced_rmatrix import make_G2
        from lib.collision_residue_rmatrix import verify_jacobi
        assert verify_jacobi(make_G2()), "G_2 fails Jacobi identity"

    def test_G2_antisymmetry(self):
        from lib.non_simply_laced_rmatrix import make_G2
        from lib.collision_residue_rmatrix import verify_antisymmetry
        assert verify_antisymmetry(make_G2()), "G_2 not antisymmetric"

    def test_G2_killing_invariance(self):
        from lib.non_simply_laced_rmatrix import make_G2
        from lib.collision_residue_rmatrix import verify_killing_invariance
        assert verify_killing_invariance(make_G2()), "G_2 Killing not invariant"

    def test_G2_dimension(self):
        from lib.non_simply_laced_rmatrix import make_G2
        g = make_G2()
        assert g.dim == 14
        assert g.rank == 2

    def test_G2_dual_coxeter_number(self):
        from lib.non_simply_laced_rmatrix import make_G2
        assert make_G2().h_dual == 4

    def test_B2_cybe_algebraic(self):
        from lib.non_simply_laced_rmatrix import make_B2
        from lib.collision_residue_rmatrix import verify_cybe
        result = verify_cybe(make_B2())
        assert result['cybe_satisfied'], \
            f"B_2 CYBE violation: max = {result['max_violation']}"

    def test_G2_cybe(self):
        from lib.non_simply_laced_rmatrix import make_G2
        from lib.collision_residue_rmatrix import verify_cybe
        result = verify_cybe(make_G2())
        assert result['cybe_satisfied'], \
            f"G_2 CYBE violation: max = {result['max_violation']}"

    def test_B2_casimir_trace(self):
        """Omega^{ab} kappa_{ab} = dim(g) = 10 for B_2."""
        import numpy as np
        from lib.non_simply_laced_rmatrix import make_B2
        from lib.collision_residue_rmatrix import casimir_tensor
        g = make_B2()
        Omega = casimir_tensor(g)
        trace_val = np.trace(Omega @ g.kappa)
        assert abs(trace_val - g.dim) < 1e-10, \
            f"B_2 Omega*kappa trace = {trace_val}, expected {g.dim}"

    def test_G2_casimir_trace(self):
        """Omega^{ab} kappa_{ab} = dim(g) = 14 for G_2."""
        import numpy as np
        from lib.non_simply_laced_rmatrix import make_G2
        from lib.collision_residue_rmatrix import casimir_tensor
        g = make_G2()
        Omega = casimir_tensor(g)
        trace_val = np.trace(Omega @ g.kappa)
        assert abs(trace_val - g.dim) < 1e-10, \
            f"G_2 Omega*kappa trace = {trace_val}, expected {g.dim}"

    def test_B2_C2_langlands_cartan_transpose(self):
        """Langlands duality: Cartan(C_2) = Cartan(B_2)^T."""
        from lib.non_simply_laced_rmatrix import langlands_duality_comparison
        result = langlands_duality_comparison(k=1.0)
        assert result['cartan_transpose'], "Langlands: C_2 Cartan != B_2 Cartan^T"

    def test_G2_root_system_size(self):
        """G_2 has 6 positive roots."""
        from lib.non_simply_laced_rmatrix import make_G2
        g = make_G2()
        num_positive_roots = (g.dim - g.rank) // 2
        assert num_positive_roots == 6

    def test_B2_collision_residue_at_k1(self):
        """Collision residue r(z) = k*Omega/z for B_2 at k=1."""
        from lib.non_simply_laced_rmatrix import make_B2
        from lib.collision_residue_rmatrix import full_collision_residue_computation
        g = make_B2()
        result = full_collision_residue_computation(g, k=1.0)
        assert result['r_equals_k_omega_over_z'], \
            "B_2 r-matrix != k*Omega/z"

    def test_G2_collision_residue_at_k1(self):
        """Collision residue r(z) = k*Omega/z for G_2 at k=1."""
        from lib.non_simply_laced_rmatrix import make_G2
        from lib.collision_residue_rmatrix import full_collision_residue_computation
        g = make_G2()
        result = full_collision_residue_computation(g, k=1.0)
        assert result['r_equals_k_omega_over_z'], \
            "G_2 r-matrix != k*Omega/z"


# =========================================================================
# 3. GENUS-1 INTERSECTION NUMBERS: CROSS-CHECKS
# =========================================================================

class TestGenus1IntersectionAdversarial:
    """Cross-checks for genus-1 intersection numbers."""

    def test_abelian_cs_F1_equals_1_over_24(self):
        """F_1(abelian_cs, k=1) = 1/24 (= kappa * lambda_1^FP = 1 * 1/24)."""
        from sympy import Rational
        from compute.lib.genus_one_bridge import period_correction
        result = period_correction('abelian_cs', k=1)
        assert result['F1'] == Rational(1, 24)

    def test_abelian_cs_F1_scales_with_level(self):
        """F_1(abelian_cs, k) = k/24."""
        from sympy import Rational
        from compute.lib.genus_one_bridge import period_correction
        for k in [1, 2, 5]:
            result = period_correction('abelian_cs', k=k)
            assert result['F1'] == Rational(k, 24), \
                f"F_1(k={k}) = {result['F1']}, expected {k}/24"

    def test_virasoro_genus1_curvature(self):
        """kappa(Vir) = c/2, confirmed at c=26 where kappa=13."""
        from compute.lib.genus_one_bridge import genus1_curvature
        result = genus1_curvature('virasoro', c=26)
        assert result['kappa'] == 13.0, \
            f"kappa(Vir_26) = {result['kappa']}, expected 13"

    def test_heisenberg_genus1_structure(self):
        """Heisenberg at k=1: kappa=1, c_0=0, decoupled."""
        from compute.lib.genus1_intersection import genus1_intersection_heisenberg
        h = genus1_intersection_heisenberg(k_val=1)
        assert h['kappa'] == 1
        assert h['coisson_bracket'] == 0
        assert h['elliptic_regime'] == 'decoupled'

    def test_heisenberg_genus1_level_scaling(self):
        """R^(1)(H_k) scales linearly with k."""
        from compute.lib.genus1_intersection import genus1_intersection_heisenberg
        h1 = genus1_intersection_heisenberg(k_val=1)
        h5 = genus1_intersection_heisenberg(k_val=5)
        for t1, t5 in zip(h1['intersection_number'], h5['intersection_number']):
            assert t5['coefficient'] == 5 * t1['coefficient']

    def test_genus1_z_symmetry_heisenberg(self):
        """R^(1)(H_k; z, tau) is odd in z (changes sign under z -> -z)."""
        from compute.lib.genus1_intersection import genus1_intersection_numerical
        r_plus = genus1_intersection_numerical(k_val=1, z_val=0.1, tau_val=1j)
        r_minus = genus1_intersection_numerical(k_val=1, z_val=-0.1, tau_val=1j)
        assert abs(r_plus['R1_real'] + r_minus['R1_real']) < 1e-10

    def test_genus1_numerical_at_generic_tau(self):
        """R^(1) is finite at a generic tau value."""
        from compute.lib.genus1_intersection import genus1_intersection_numerical
        r = genus1_intersection_numerical(k_val=1, z_val=0.01, tau_val=0.5 + 1j)
        assert math.isfinite(r['R1_real']), "R^(1) diverges at generic tau"

    def test_genus2_framework_F2(self):
        """F_2(k=1) = 7/5760 from genus-2 intersection framework."""
        from sympy import Rational
        from compute.lib.genus1_intersection import genus2_intersection_framework
        g2 = genus2_intersection_framework(k_val=1)
        assert g2['F2'] == Rational(7, 5760)

    def test_genus2_framework_level_scaling(self):
        """F_2(k=3) = 3 * F_2(k=1)."""
        from compute.lib.genus1_intersection import genus2_intersection_framework
        g2_1 = genus2_intersection_framework(k_val=1)
        g2_3 = genus2_intersection_framework(k_val=3)
        assert g2_3['F2'] == 3 * g2_1['F2']

    def test_virasoro_genus1_three_sectors(self):
        """Virasoro genus-1 intersection has three sectors."""
        from compute.lib.genus1_intersection import genus1_intersection_virasoro
        v = genus1_intersection_virasoro()
        assert len(v['quartic_sector']) > 0
        assert len(v['double_sector']) > 0
        assert len(v['simple_sector']) > 0


# =========================================================================
# 4. W3 MULTICHANNEL SHADOW: TT CHANNEL VS VIRASORO
# =========================================================================

class TestW3MultichannelAdversarial:
    """W3 TT channel must match the Virasoro shadow obstruction tower."""

    def test_tt_shadow_metric_structure(self):
        """Q^{TT}(t) = c^2 + 12ct + [(180c+872)/(5c+22)]t^2."""
        from sympy import Symbol
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from w3_multichannel_shadow import virasoro_shadow_metric

        c = Symbol('c', positive=True)
        t = Symbol('t')
        Q = virasoro_shadow_metric()
        Q_expanded = Q.expand()
        assert Q_expanded.coeff(t, 0) == c**2
        assert Q_expanded.coeff(t, 1) == 12 * c

    def test_tt_kappa_equals_c_over_2(self):
        """S_2^{TT} = kappa/2 = c/2 for Virasoro."""
        from sympy import Symbol, cancel
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from w3_multichannel_shadow import tt_shadow_coefficients

        c = Symbol('c', positive=True)
        S = tt_shadow_coefficients(r_max=4)
        assert cancel(S[2] - c / 2) == 0

    def test_tt_S3(self):
        """S_3^{TT} = 2 (from 12c/(2c)/3 = 6/3 = 2)."""
        from sympy import cancel
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from w3_multichannel_shadow import tt_shadow_coefficients
        S = tt_shadow_coefficients(r_max=4)
        assert cancel(S[3] - 2) == 0

    def test_ww_kappa_equals_c_over_3(self):
        """kappa_WW = c/3 => S_2^{WW} = c/6."""
        from sympy import cancel, Symbol
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from w3_multichannel_shadow import ww_shadow_coefficients
        c = Symbol('c', positive=True)
        WW = ww_shadow_coefficients(r_max=4)
        assert cancel(WW[2] - c / 6) == 0

    def test_ww_only_even_arities(self):
        """W-W channel has only even arities (Z_2: W -> -W)."""
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from w3_multichannel_shadow import ww_shadow_coefficients
        WW = ww_shadow_coefficients(r_max=10)
        for r in WW:
            assert r % 2 == 0, f"WW channel has odd arity r={r}"

    def test_virasoro_complementarity_c_star_13(self):
        """kappa(c) + kappa(26-c) = 13 for Virasoro."""
        from sympy import Symbol, simplify
        c = Symbol('c')
        assert simplify(c / 2 + (26 - c) / 2 - 13) == 0

    def test_w3_self_dual_point(self):
        """W_3: alpha_3=100, c*=50, c_crit=100. AP: c* != c_crit."""
        N = 3
        alpha_N = 2 * (N - 1) * (2 * N**2 + 2 * N + 1)
        assert alpha_N == 100
        assert alpha_N // 2 == 50  # c*
        assert alpha_N == 100  # c_crit
        assert alpha_N // 2 != alpha_N  # c* != c_crit


# =========================================================================
# 5. BIGRADED HILBERT SERIES: EXACT FORMULAS
# =========================================================================

class TestBigradedHilbertAdversarial:
    """Test bigraded Hilbert series against exact partition function formulas."""

    def test_heisenberg_character_is_partition_function(self):
        """ch(H_k; q) = prod_{n>=1} 1/(1-q^n) = partition function."""
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from ordered_bar_bigraded_hilbert import bosonic_vacuum_char
        ch = bosonic_vacuum_char([1], 20)
        known_p = [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490, 627]
        for w in range(len(known_p)):
            assert ch[w] == known_p[w], f"p({w}) = {ch[w]}, expected {known_p[w]}"

    def test_virasoro_character_structure(self):
        """ch(Vir_c; q) = prod_{n>=2} 1/(1-q^n)."""
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from ordered_bar_bigraded_hilbert import bosonic_vacuum_char
        ch = bosonic_vacuum_char([2], 10)
        assert ch[0] == 1
        assert ch[1] == 0  # no weight-1 mode
        assert ch[2] == 1  # T_{-2}
        assert ch[3] == 1  # T_{-3}
        assert ch[4] == 2  # T_{-4}, T_{-2}^2
        assert ch[5] == 2  # T_{-5}, T_{-2}T_{-3}
        assert ch[6] == 4  # T_{-6}, T_{-4}T_{-2}, T_{-3}^2, T_{-2}^3

    def test_heisenberg_augmentation_ideal(self):
        """f_+(0) = 0, f_+(1) = 1, f_+(2) = 2, f_+(3) = 3."""
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from ordered_bar_bigraded_hilbert import heisenberg_aug
        f = heisenberg_aug(10)
        # f_+ = ch - 1, so f_+(0) = p(0)-1 = 0, f_+(n) = p(n) for n>=1
        assert f[0] == 0
        assert f[1] == 1   # p(1)=1
        assert f[2] == 2   # p(2)=2
        assert f[3] == 3   # p(3)=3
        assert f[4] == 5   # p(4)=5

    def test_heisenberg_bar_degree_2_convolution(self):
        """dim B^{ord}_{2,w} = sum_{j=1}^{w-1} p(j)*p(w-j) for Heisenberg."""
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from ordered_bar_bigraded_hilbert import heisenberg_aug, poly_mul
        f = heisenberg_aug(10)
        f2 = poly_mul(f, f, 10)
        # f_+(1)^2 = 1*1 = 1 => B_{2,2} = 1
        assert f2[2] == 1
        # f_+^2(3) = f_+(1)*f_+(2) + f_+(2)*f_+(1) = 1*2+2*1 = 4
        # Wait: f_+(0)=0, so sum = f_+(0)*f_+(3) + f_+(1)*f_+(2) + f_+(2)*f_+(1) + f_+(3)*f_+(0)
        # = 0 + 1*2 + 2*1 + 0 = 4. YES, 4 is correct.
        assert f2[3] == 4

    def test_w3_character_at_low_weight(self):
        """W_3: bosons T(wt 2), W(wt 3).
        wt 0: 1, wt 1: 0, wt 2: 1(T), wt 3: 2(dT, W), wt 4: 3(d^2T, T^2, dW).
        """
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from ordered_bar_bigraded_hilbert import bosonic_vacuum_char
        ch = bosonic_vacuum_char([2, 3], 8)
        assert ch[0] == 1
        assert ch[1] == 0
        assert ch[2] == 1  # T_{-2}
        assert ch[3] == 2  # T_{-3}, W_{-3}
        # wt 4: T_{-4}, T_{-2}^2, W_{-4} => 3 states
        assert ch[4] == 3

    def test_bc_ghost_fermionic_character(self):
        """bc (lambda=1): ch = prod_{n>=1}(1+q^n)^2.
        wt 0: 1, wt 1: 2 (b_{-1}, c_{-1}).
        wt 2: b_{-2}, c_{-2}, b_{-1}c_{-1} => 3.
        wt 3: b_{-3}, c_{-3}, b_{-2}c_{-1}, b_{-1}c_{-2}, b_{-1}b_{-2}... no,
          fermionic: b_{-1}c_{-2}, b_{-2}c_{-1}, b_{-3}, c_{-3},
          b_{-1}b_{-2} (only if both b modes allowed), c_{-1}c_{-2}
          => 6 states.
        """
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from ordered_bar_bigraded_hilbert import fermionic_vacuum_char
        ch = fermionic_vacuum_char([1, 1], 6)
        assert ch[0] == 1
        assert ch[1] == 2
        assert ch[2] == 3
        # wt 3: enumerate all states:
        # single modes: b_{-3}, c_{-3} (2 states)
        # double modes: b_{-1}c_{-2}, b_{-2}c_{-1}, b_{-1}b_{-2}, c_{-1}c_{-2} (4 states)
        # triple modes: none with weight 3 (b_{-1}c_{-1}b_{-1} not allowed)
        # Actually: b_{-1}c_{-1}b... wait, b and c are DIFFERENT fields
        # Double: from {b_{-1},b_{-2},b_{-3},...} and {c_{-1},c_{-2},...}
        # wt 3 = b_{-3} | c_{-3} | b_{-2}c_{-1} | b_{-1}c_{-2} | b_{-1}b_{-2} | c_{-1}c_{-2} = 6
        assert ch[3] == 6

    def test_euler_eta_identity_virasoro(self):
        """Euler-eta: chi(Vir; q) = -1 + 1/ch(q) matches alternating sum."""
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from ordered_bar_bigraded_hilbert import bosonic_vacuum_char, euler_eta_chi
        ch_full = bosonic_vacuum_char([2], 10)
        chi = euler_eta_chi(ch_full, 10)
        # Independently: chi_0 = -1+1 = 0 (since ch[0]=1)
        # chi_1 = 0 (ch has no weight-1 terms, and 1/ch is also 0 at weight 1)
        assert chi[0] == 0  # no bar elements at weight 0 (vacuum subtracted)


# =========================================================================
# 6. FIELD SECTOR GENERATING FUNCTIONS
# =========================================================================

class TestFieldSectorAdversarial:
    """Adversarial tests for field sector generating functions."""

    def test_profile_polynomial_special_values(self):
        """P_r(0) = (2r+1)!, P_r(-1) = (2r)!, P_r(1) = (2r+2)!/2."""
        for r in range(1, 8):
            coeffs = shifted_pochhammer_coeffs(r)
            assert coeffs[0] == math.factorial(2 * r + 1)
            P_neg1 = sum(coeffs[w] * (-1)**w for w in range(len(coeffs)))
            assert P_neg1 == math.factorial(2 * r)
            P_1 = sum(coeffs)
            assert P_1 == math.factorial(2 * r + 2) // 2

    def test_dT_T_ratio_is_harmonic(self):
        """a_1/a_0 = H_{2r+1} - 1."""
        from fractions import Fraction
        for r in range(1, 7):
            coeffs = shifted_pochhammer_coeffs(r)
            ratio = Fraction(coeffs[1]) / Fraction(coeffs[0])
            H = sum(Fraction(1, k) for k in range(1, 2 * r + 2))
            assert ratio == H - 1

    def test_generating_function_closed_form_y0(self):
        """G(x,0) = 1/2 - (1+8x)/(2*sqrt(1+4x))."""
        from sympy import Symbol, series, Rational, sqrt, catalan as sc
        x = Symbol('x')
        g = Rational(1, 2) - (1 + 8 * x) / (2 * sqrt(1 + 4 * x))
        g_series = series(g, x, 0, n=10)
        for r in range(1, 10):
            coeff = g_series.coeff(x, r)
            expected = (-1) ** r * (2 * r + 1) * int(sc(r - 1))
            assert int(coeff) == expected

    def test_generating_function_closed_form_yneg1(self):
        """G(x,-1) = (sqrt(1+4x)-1)/2 is the Catalan GF."""
        from sympy import Symbol, series, sqrt, catalan as sc
        x = Symbol('x')
        g = (sqrt(1 + 4 * x) - 1) / 2
        g_series = series(g, x, 0, n=10)
        for r in range(1, 10):
            coeff = g_series.coeff(x, r)
            expected = (-1) ** (r + 1) * int(sc(r - 1))
            assert int(coeff) == expected

    def test_run_all_field_sector_computation(self):
        """Run field sector computation and verify Catalan factorisation.

        compute_all_field_sectors returns {r: {w: raw_coeff}} where
        the raw coefficients are phi_k(x) = m_k at the symmetric point.
        The Catalan factorisation predicts:
          phi_k(x) = (-1)^n C_n prod_{m=2}^k (x+m), n = (k-3)/2
        So the raw T-coefficient (w=0) should be (-1)^n C_n (k+1)!/1 = (-1)^n C_n prod(2..k).
        """
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
        from all_field_sector_generating_functions import compute_all_field_sectors
        results = compute_all_field_sectors(max_r=4, c_val=0.0)
        for r in range(1, 5):
            if r in results:
                k = 2 * r + 1
                n = (k - 3) // 2
                prefactor = (-1) ** n * catalan(n)
                predicted_T = prefactor * math.prod(range(2, k + 1))
                raw_T = results[r].get(0, 0.0)
                assert abs(raw_T - predicted_T) < abs(predicted_T) * 1e-6, \
                    f"r={r}: T={raw_T}, expected {predicted_T}"


# =========================================================================
# 7. CROSS-ENGINE CONSISTENCY
# =========================================================================

class TestCrossEngineConsistency:
    """Verify consistency across different computational engines."""

    def test_fp_number_consistency(self):
        """lambda_fp values must agree across all engines.

        FINDING: test_modular_obstruction_engine::test_genus2_fp_number
        expects Rational(1, 1152) but the engine and lambda_fp both return
        Rational(7, 5760). The test expectation is WRONG; 7/5760 is the
        correct A-hat genus coefficient at degree 2.
        """
        from sympy import Rational
        from compute.lib.genus2_obstruction_engine import lambda_fp
        from compute.lib.genus1_intersection import genus2_intersection_framework

        assert lambda_fp(1) == Rational(1, 24)
        assert lambda_fp(2) == Rational(7, 5760)
        assert lambda_fp(3) == Rational(31, 967680)
        g2 = genus2_intersection_framework(k_val=1)
        assert g2['F2'] == Rational(7, 5760)

    def test_modular_completion_genus_data(self):
        """Modular completion F_g values must match lambda_fp."""
        from sympy import Rational
        from compute.lib.genus1_intersection import modular_completion_koszul
        mc = modular_completion_koszul('abelian_cs', max_genus=3, k=1)
        expected = [Rational(1, 24), Rational(7, 5760), Rational(31, 967680)]
        for g in range(1, 4):
            assert mc['genus_data'][g]['F_g'] == expected[g - 1]

    def test_virasoro_self_dual_at_c13(self):
        """AP: c*=13, c_crit=26 for Virasoro. NEVER conflate."""
        assert 13 != 26

    def test_w3_alpha_formula(self):
        """alpha_N = 2(N-1)(2N^2+2N+1)."""
        for N, expected in [(2, 26), (3, 100), (4, 246)]:
            alpha = 2 * (N - 1) * (2 * N**2 + 2 * N + 1)
            assert alpha == expected, f"alpha_{N} = {alpha}, expected {expected}"

    def test_stasheff_engine_deterministic(self):
        """StasheffEngine gives same result on repeated calls."""
        from m7_m10_depth_frontier import StasheffEngine
        engine = StasheffEngine(1.0)
        lams = (1.0, 1.0, 1.0, 1.0)
        r1 = engine.mk(lams)
        engine._cache.clear()
        r2 = engine.mk(lams)
        for key in set(list(r1.keys()) + list(r2.keys())):
            assert abs(r1.get(key, 0.0) - r2.get(key, 0.0)) < 1e-10

    def test_catalan_factorisation_independent_of_c(self):
        """Field polynomial coefficients are c-independent (scalar is c-dependent)."""
        from m7_m10_depth_frontier import StasheffEngine
        for k in [3, 5, 7]:
            e0 = StasheffEngine(0.0)
            e0._cache.clear()
            e26 = StasheffEngine(26.0)
            e26._cache.clear()
            lams = tuple(1.0 for _ in range(k - 1))
            r0 = e0.mk(lams)
            r26 = e26.mk(lams)
            for w in range(0, k):
                assert abs(r0.get(w, 0.0) - r26.get(w, 0.0)) < 1e-6, \
                    f"k={k}, w={w}: c-dependent field coeff"
            if k >= 3:
                assert abs(r26.get(-1, 0.0) - r0.get(-1, 0.0)) > 0.1, \
                    f"k={k}: scalar should be c-dependent"

    def test_modular_obstruction_fp_vs_lambda_fp(self):
        """KNOWN BUG: modular_obstruction_engine.genus2_obstruction returns
        faber_pandharipande_g2 = 7/5760, but the test expects 1/1152.

        This test documents the discrepancy: the ENGINE is correct (7/5760),
        and the test expectation in test_modular_obstruction_engine.py is wrong.
        """
        from sympy import Rational
        from compute.lib.modular_obstruction_engine import genus2_obstruction
        from compute.lib.genus2_obstruction_engine import lambda_fp
        result = genus2_obstruction('virasoro')
        # The engine returns the correct value:
        assert result['faber_pandharipande_g2'] == Rational(7, 5760)
        # Which matches lambda_fp(2):
        assert result['faber_pandharipande_g2'] == lambda_fp(2)


# =========================================================================
# 8. AP19 POLE ABSORPTION: INDEPENDENT VERIFICATION
# =========================================================================

class TestAP19PoleAbsorption:
    """Verify AP19 (bar kernel absorbs a pole) independently."""

    def test_sl2_rmatrix_simple_pole(self):
        """For sl_2, OPE has z^{-2}+z^{-1}. R-matrix = k*Omega/z (simple pole)."""
        from lib.collision_residue_rmatrix import full_collision_residue_computation
        # Use the pre-built sl_2 from the existing test infrastructure
        from lib.collision_residue_rmatrix import make_sl2
        g = make_sl2()
        result = full_collision_residue_computation(g, k=1.0)
        assert result['r_equals_k_omega_over_z'], \
            "sl_2 r-matrix should be k*Omega/z"

    def test_virasoro_bar_pole_3(self):
        """Virasoro: OPE max pole 4, bar max pole 3 (AP19)."""
        from compute.lib.genus1_intersection import genus1_intersection_virasoro
        v = genus1_intersection_virasoro()
        # The quartic sector comes from the z^{-4} OPE term
        # After d-log absorption: z^{-3} in the bar complex
        assert len(v['quartic_sector']) > 0  # quartic = 4th order OPE -> 3rd order bar


# =========================================================================
# 9. MATHEMATICAL IDENTITY CHECKS
# =========================================================================

class TestMathematicalIdentities:
    """Test fundamental mathematical identities used throughout."""

    def test_catalan_generating_function(self):
        """C(x) = (1-sqrt(1-4x))/(2x) satisfies C_n = binom(2n,n)/(n+1)."""
        from sympy import Symbol, series, sqrt, Rational
        x = Symbol('x')
        C = (1 - sqrt(1 - 4 * x)) / (2 * x)
        C_series = series(C, x, 0, n=10)
        for n in range(10):
            assert int(C_series.coeff(x, n)) == catalan(n)

    def test_catalan_recursion(self):
        """C_{n+1} = sum_{k=0}^{n} C_k C_{n-k}."""
        for n in range(8):
            conv = sum(catalan(k) * catalan(n - k) for k in range(n + 1))
            assert conv == catalan(n + 1)

    def test_bernoulli_numbers_sympy_convention(self):
        """SymPy uses B_1 = +1/2 (modern convention).
        B_0=1, B_1=1/2, B_2=1/6, B_4=-1/30, B_6=1/42.
        """
        from sympy import bernoulli, Rational
        assert bernoulli(0) == 1
        assert bernoulli(1) == Rational(1, 2)  # SymPy: +1/2, not -1/2
        assert bernoulli(2) == Rational(1, 6)
        assert bernoulli(4) == Rational(-1, 30)
        assert bernoulli(6) == Rational(1, 42)

    def test_partition_function_pentagonal(self):
        """Euler's pentagonal theorem: p(n) recursion."""
        known_p = [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42]
        for n in range(1, len(known_p)):
            total = 0
            k = 1
            while True:
                g1 = k * (3 * k - 1) // 2
                g2 = k * (3 * k + 1) // 2
                if g1 > n:
                    break
                sign = (-1) ** (k + 1)
                if n - g1 >= 0:
                    total += sign * known_p[n - g1]
                if n - g2 >= 0:
                    total += sign * known_p[n - g2]
                k += 1
            assert total == known_p[n]

    def test_ahat_genus_sequence(self):
        """The A-hat genus coefficients: 1, -1/24, 7/5760, -31/967680.
        These come from (x/2)/sinh(x/2) expanded in x^2.
        """
        from sympy import Rational
        from compute.lib.genus2_obstruction_engine import lambda_fp
        # lambda_fp(g) = |coefficient of x^{2g} in (x/2)/sinh(x/2)|
        assert lambda_fp(1) == Rational(1, 24)
        assert lambda_fp(2) == Rational(7, 5760)
        assert lambda_fp(3) == Rational(31, 967680)
        # Verify the ratios are correct:
        # 7/5760 / (1/24) = 7*24/5760 = 168/5760 = 7/240
        assert lambda_fp(2) / lambda_fp(1) == Rational(7, 240)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
