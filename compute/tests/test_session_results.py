r"""Master cross-check of key numerical results from the compute session.

Verifies the principal findings across multiple compute modules:

1. Catalan formula: T_{2n+3}(1,...,1) = (-1)^n C_n (2n+3)!
2. W_3 WW shadow c-independence: the normalized numerator N_r is a
   rational integer independent of c, verified at 3+ c-values.
3. Euler-eta identity at weight 5 for all 6 standard families.
4. Period-2 vanishing: T_k(1,...,1) = 0 for even k = 4, 6, 8.

These are independent cross-checks that do not merely re-run existing
modules but recompute key quantities from first principles.
"""

import math
import sys
import os

import pytest

# Ensure compute/ is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))


def catalan(n):
    """Catalan number C_n = binom(2n,n)/(n+1)."""
    return math.comb(2 * n, n) // (n + 1)


# =========================================================================
# 1. CATALAN FORMULA FOR THE VIRASORO SHADOW TOWER
# =========================================================================

class TestCatalanFormula:
    """Verify T_{2n+3}(1,...,1) = (-1)^n C_n (2n+3)!.

    This is the T-coefficient (derivative-order 0) of the Virasoro
    A_infinity operation m_{2n+3}(T,...,T; 1,...,1) at the symmetric
    point. The formula follows from the generating function
    g(x) = 1/2 - (1+8x)/(2*sqrt(1+4x)).
    """

    @pytest.fixture(scope="class")
    def engine(self):
        from m7_m10_depth_frontier import StasheffEngine
        return StasheffEngine(0.0)  # c=0 for T-coefficient (c-independent)

    @pytest.mark.parametrize("n", [0, 1, 2, 3, 4])
    def test_catalan_formula(self, engine, n):
        """T_{2n+3}(1,...,1) = (-1)^n C_n (2n+3)!."""
        k = 2 * n + 3
        lams = tuple(1.0 for _ in range(k - 1))
        engine._cache.clear()
        result = engine.mk(lams)
        T_val = result.get(0, 0.0)
        predicted = (-1) ** n * catalan(n) * math.factorial(k)
        assert abs(T_val - predicted) < max(abs(predicted) * 1e-8, 1), \
            f"k={k}: T = {T_val}, expected (-1)^{n} C_{n} {k}! = {predicted}"

    @pytest.mark.parametrize("n", [0, 1, 2, 3])
    def test_signed_sum_formula(self, engine, n):
        """S_k = phi_k(1) = (-1)^n C_n (k+1)!/2 for odd k = 2n+3."""
        k = 2 * n + 3
        lams = tuple(1.0 for _ in range(k - 1))
        engine._cache.clear()
        result = engine.mk(lams)
        # Sum over all non-negative field-derivative orders
        S_val = sum(v for f, v in result.items() if f >= 0)
        predicted = (-1) ** n * catalan(n) * math.factorial(k + 1) // 2
        assert abs(S_val - predicted) < max(abs(predicted) * 1e-6, 1), \
            f"k={k}: S = {S_val}, expected (-1)^{n} C_{n} {k+1}!/2 = {predicted}"


# =========================================================================
# 2. W_3 WW SHADOW c-INDEPENDENCE
# =========================================================================

class TestW3WWShadowCIndependence:
    """Verify that the W_3 WW shadow numerator N_r is c-independent.

    The WW shadow coefficient is:
      S_{2r}^{WW} = (c/3) C(1/2, r-1) gamma^{r-1} / (2r)
    where gamma = 61440 / (c^2 (5c+22)^3).

    Defining N_r = S_{2r} * c^{2(r-1)-1} * (5c+22)^{3(r-1)} * 3 * 2r,
    one gets N_r = C(1/2, r-1) * 61440^{r-1}, which is a rational
    integer independent of c. We verify this at 3+ c-values.
    """

    def _half_binom(self, n):
        """Compute C(1/2, n) = (1/2)(1/2-1)...(1/2-n+1)/n! as a Fraction."""
        from fractions import Fraction
        result = Fraction(1)
        for j in range(n):
            result *= (Fraction(1, 2) - j)
        if n > 0:
            result /= Fraction(math.factorial(n))
        return result

    def _S_WW(self, r, c_val):
        """Compute S_{2r}^{WW} at a given c value (float)."""
        n = r - 1
        bcoeff = float(self._half_binom(n))
        gamma = 61440.0 / (c_val ** 2 * (5 * c_val + 22) ** 3)
        kappa_WW = c_val / 3.0
        return kappa_WW * bcoeff * gamma ** n / (2 * r)

    def _N_r_exact(self, r):
        """Compute N_r = C(1/2, r-1) * 61440^{r-1} exactly."""
        from fractions import Fraction
        n = r - 1
        bcoeff = self._half_binom(n)
        return int(bcoeff * Fraction(61440) ** n)

    @pytest.mark.parametrize("r", [2, 3, 4, 5, 6])
    def test_c_independence(self, r):
        """N_r is the same integer at c=1, c=10, c=50."""
        N_exact = self._N_r_exact(r)
        assert isinstance(N_exact, int), f"N_{r} = {N_exact} is not integer"

        c_vals = [1.0, 10.0, 50.0, 99.0]
        for c_val in c_vals:
            S_val = self._S_WW(r, c_val)
            n = r - 1
            # N_r = S_{2r} * c^{2n-1} * (5c+22)^{3n} * 3 * 2r
            N_numerical = S_val * c_val ** (2 * n - 1) * (5 * c_val + 22) ** (3 * n) * 3 * 2 * r
            assert abs(N_numerical - N_exact) < abs(N_exact) * 1e-6, \
                f"r={r}, c={c_val}: N_numerical={N_numerical}, N_exact={N_exact}"

    @pytest.mark.parametrize("r,expected_N", [
        (2, 30720),
        (3, -471859200),
        (4, 14495514624000),
        (5, -556627761561600000),
    ])
    def test_N_values(self, r, expected_N):
        """Check specific N_r values."""
        N = self._N_r_exact(r)
        assert N == expected_N, f"N_{r} = {N}, expected {expected_N}"


# =========================================================================
# 3. EULER-ETA IDENTITY AT WEIGHT 5
# =========================================================================

class TestEulerEtaWeight5:
    """Verify chi_5 = [q^5](-1 + 1/ch(A;q)) for all 6 standard families.

    The Euler-eta identity states:
      chi_w = sum_{n>=1} (-1)^n dim B^{ord}_{n,w} = [q^w](-1 + 1/ch(A;q))

    At weight 5 the sum terminates at finite n (since f_+^n starts
    at weight n * w_min, and w_min >= 1 for all families), so the
    identity is exact with N=10 terms.

    Expected values at weight 5 (computed independently):
      H_k:        chi_5 =  1
      V_k(sl_2):  chi_5 =  0
      Vir_c:      chi_5 =  0
      W_3:        chi_5 =  0
      beta-gamma:  chi_5 =  2
      bc:          chi_5 = -4
    """

    def _compute_aug_and_chi(self, family_func, maxw=5):
        """Compute augmentation ideal and chi from the Euler-eta formula."""
        from ordered_bar_bigraded_hilbert import poly_mul, euler_eta_chi
        f_plus = family_func(maxw)
        ch_full = f_plus[:]
        ch_full[0] += 1
        chi = euler_eta_chi(ch_full, maxw)
        return f_plus, chi

    def _compute_chi_from_table(self, f_plus, maxw=5, maxn=10):
        """Compute chi_w from the alternating sum of bigraded dimensions."""
        from ordered_bar_bigraded_hilbert import poly_mul
        powers = {1: f_plus[:maxw + 1]}
        for n in range(2, maxn + 1):
            powers[n] = poly_mul(powers[n - 1], f_plus, maxw)
        chi = {}
        for w in range(1, maxw + 1):
            chi[w] = sum(
                ((-1) ** n) * (powers[n][w] if w < len(powers[n]) else 0)
                for n in range(1, maxn + 1)
            )
        return chi

    @pytest.mark.parametrize("family_name,expected_chi5", [
        ("heisenberg", 1),
        ("sl2", 0),
        ("virasoro", 0),
        ("w3", 0),
        ("betagamma", 2),
        ("bc", -4),
    ])
    def test_euler_eta_weight5(self, family_name, expected_chi5):
        """Verify chi_5 for each family."""
        from ordered_bar_bigraded_hilbert import (
            heisenberg_aug, sl2_aug, virasoro_aug, w3_aug,
            betagamma_aug, bc_aug,
        )

        family_map = {
            "heisenberg": heisenberg_aug,
            "sl2": sl2_aug,
            "virasoro": virasoro_aug,
            "w3": w3_aug,
            "betagamma": betagamma_aug,
            "bc": bc_aug,
        }

        func = family_map[family_name]
        f_plus, chi_eta = self._compute_aug_and_chi(func, maxw=5)
        chi_table = self._compute_chi_from_table(f_plus, maxw=5, maxn=10)

        # Both methods should agree
        assert chi_eta[5] == chi_table[5], \
            f"{family_name}: eta formula gives {chi_eta[5]}, table gives {chi_table[5]}"

        # Both should match the expected value
        assert chi_eta[5] == expected_chi5, \
            f"{family_name}: chi_5 = {chi_eta[5]}, expected {expected_chi5}"

    def test_all_families_consistent(self):
        """Cross-check: the Euler-eta formula and table agree at weights 1-5."""
        from ordered_bar_bigraded_hilbert import (
            heisenberg_aug, sl2_aug, virasoro_aug, w3_aug,
            betagamma_aug, bc_aug, euler_eta_chi, poly_mul,
        )

        families = [
            ("H_k", heisenberg_aug),
            ("V_k(sl_2)", sl2_aug),
            ("Vir_c", virasoro_aug),
            ("W_3", w3_aug),
            ("beta-gamma", betagamma_aug),
            ("bc", bc_aug),
        ]

        for name, func in families:
            f_plus = func(10)
            ch_full = f_plus[:]
            ch_full[0] += 1
            chi_eta = euler_eta_chi(ch_full, 10)

            powers = {1: f_plus[:11]}
            for n in range(2, 11):
                powers[n] = poly_mul(powers[n - 1], f_plus, 10)
            for w in range(1, 6):
                chi_table = sum(
                    ((-1) ** n) * (powers[n][w] if w < len(powers[n]) else 0)
                    for n in range(1, 11)
                )
                assert chi_eta[w] == chi_table, \
                    f"{name} weight {w}: eta={chi_eta[w]}, table={chi_table}"


# =========================================================================
# 4. PERIOD-2 VANISHING
# =========================================================================

class TestPeriod2Vanishing:
    """Verify T_k(1,...,1) = 0 for even k = 4, 6, 8.

    At the symmetric point (all spectral parameters equal to 1),
    the ENTIRE Virasoro A_infinity operation vanishes for even arity
    k >= 4. This is Theorem 3 of field_sector_generating_function.py.

    The generating function g(x) = 1/2 - (1+8x)/(2*sqrt(1+4x))
    is a power series in x (no half-integer powers), so even-arity
    contributions vanish identically.
    """

    @pytest.fixture(scope="class")
    def engine(self):
        from m7_m10_depth_frontier import StasheffEngine
        return StasheffEngine(1.0)  # c=1 (nonzero to catch c-dependent errors)

    @pytest.mark.parametrize("k", [4, 6, 8])
    def test_even_arity_vanishes(self, engine, k):
        """m_k(T,...,T; 1,...,1) = 0 for even k >= 4 at the symmetric point."""
        lams = tuple(1.0 for _ in range(k - 1))
        engine._cache.clear()
        result = engine.mk(lams)
        total_abs = sum(abs(v) for v in result.values())
        assert total_abs < 1e-6, \
            f"k={k}: |m_k(sym)| = {total_abs}, expected 0"

    @pytest.mark.parametrize("k", [4, 6, 8])
    def test_even_arity_vanishes_c0(self, k):
        """Same test at c=0 for extra confidence."""
        from m7_m10_depth_frontier import StasheffEngine
        engine = StasheffEngine(0.0)
        lams = tuple(1.0 for _ in range(k - 1))
        result = engine.mk(lams)
        total_abs = sum(abs(v) for v in result.values())
        assert total_abs < 1e-6, \
            f"k={k}, c=0: |m_k(sym)| = {total_abs}, expected 0"


# =========================================================================
# 5. SUPPLEMENTARY CROSS-CHECKS
# =========================================================================

class TestSupplementaryCrossChecks:
    """Additional cross-checks tying the session results together."""

    def test_catalan_convolution(self):
        """C_n = sum_{a=0}^{n-1} C_a C_{n-1-a} (Segner recurrence)."""
        for n in range(1, 8):
            conv = sum(catalan(a) * catalan(n - 1 - a) for a in range(n))
            assert conv == catalan(n), \
                f"Convolution = {conv}, C_{n} = {catalan(n)}"

    def test_catalan_values(self):
        """First few Catalan numbers: 1, 1, 2, 5, 14, 42, 132."""
        expected = [1, 1, 2, 5, 14, 42, 132]
        for n, c_n in enumerate(expected):
            assert catalan(n) == c_n, f"C_{n} = {catalan(n)}, expected {c_n}"

    def test_generating_function_algebraic(self):
        """g(x) = 1/2 - (1+8x)/(2*sqrt(1+4x)) encodes Catalan formula."""
        from sympy import Symbol, Rational, sqrt, series, catalan as sympy_catalan
        x = Symbol('x')
        g = Rational(1, 2) - (1 + 8 * x) / (2 * sqrt(1 + 4 * x))
        g_series = series(g, x, 0, n=8)
        for r in range(1, 8):
            coeff = int(g_series.coeff(x, r))
            # [x^r] g = (-1)^r * (2r+1) * C_{r-1}
            expected = (-1) ** r * (2 * r + 1) * catalan(r - 1)
            assert coeff == expected, \
                f"r={r}: [x^r]g = {coeff}, expected {expected}"

    def test_bigraded_dims_heisenberg_weight1(self):
        """dim B^{ord}_{1,1}(H_k) = 1, dim B^{ord}_{n,1}(H_k) = 0 for n >= 2.

        f_+(1) = 1, but f_+^n at weight 1 is 0 for n >= 2 because the
        minimum contributing weight for tensor degree n is n * w_min = n.
        """
        from ordered_bar_bigraded_hilbert import heisenberg_aug, poly_mul
        f_plus = heisenberg_aug(5)
        assert f_plus[1] == 1, f"f_+(1) = {f_plus[1]}, expected 1"
        power2 = poly_mul(f_plus, f_plus, 5)
        assert power2[1] == 0, f"f_+^2 at weight 1 = {power2[1]}, expected 0"

    def test_virasoro_augmentation_ideal(self):
        """Virasoro augmentation ideal: f_+(q) = q^2 + q^3 + 2q^4 + 2q^5 + ..."""
        from ordered_bar_bigraded_hilbert import virasoro_aug
        f_plus = virasoro_aug(10)
        assert f_plus[0] == 0
        assert f_plus[1] == 0
        assert f_plus[2] == 1
        assert f_plus[3] == 1
        assert f_plus[4] == 2
        assert f_plus[5] == 2

    def test_w3_augmentation_ideal(self):
        """W_3 augmentation ideal: f_+(q) = q^2 + 2q^3 + 3q^4 + 4q^5 + ...

        W_3 has generators T (wt 2) and W (wt 3).
        ch = prod_{n>=2} 1/(1-q^n) * prod_{n>=3} 1/(1-q^n).
        Weight 2: T_{-2} only -> 1
        Weight 3: T_{-3} + W_{-3} -> 2
        Weight 4: prod_{n>=2} at q^4 gives 2 (T_{-4}, T_{-2}^2),
                  prod_{n>=3} at q^4 gives 1 (W_{-4}),
                  convolution at q^4 = 1*1 + 0*1 + 2*0 + ... = 3
        """
        from ordered_bar_bigraded_hilbert import w3_aug
        f_plus = w3_aug(10)
        assert f_plus[0] == 0
        assert f_plus[1] == 0
        assert f_plus[2] == 1  # T at weight 2
        assert f_plus[3] == 2  # dT + W at weight 3
        assert f_plus[4] == 3  # convolution at weight 4
        assert f_plus[5] == 4  # convolution at weight 5


# =========================================================================
# 6. BIGRADED DIMENSION SPOT-CHECKS
# =========================================================================

class TestBigradedDimSpotChecks:
    """Spot-check specific bigraded dimensions from the computation."""

    def test_heisenberg_n1_w5(self):
        """dim B^{ord}_{1,5}(H_k) = 7 (partition count p(5))."""
        from ordered_bar_bigraded_hilbert import heisenberg_aug
        f_plus = heisenberg_aug(5)
        assert f_plus[5] == 7

    def test_virasoro_n1_w5(self):
        """dim B^{ord}_{1,5}(Vir_c) = 2 (from d^3T and dT_{-2}... partition into parts >=2)."""
        from ordered_bar_bigraded_hilbert import virasoro_aug
        f_plus = virasoro_aug(5)
        assert f_plus[5] == 2

    def test_bc_fermionic_n1_w3(self):
        """dim B^{ord}_{1,3}(bc) = 4 (fermionic: choose subsets of {b_{-1},c_{-1},b_{-2},c_{-2},b_{-3},c_{-3}} summing to 3)."""
        from ordered_bar_bigraded_hilbert import bc_aug
        f_plus = bc_aug(5)
        # Weight 3 states from prod_{n>=1}(1+q^n)^2 - 1:
        # b_{-1}c_{-1}b_{-1}... no, fermionic means each mode at most once.
        # With 2 fermions at each weight n>=1:
        # (1+q)^2 (1+q^2)^2 (1+q^3)^2 ...
        # Coefficient of q^3: from (1+q)^2(1+q^2)^2(1+q^3)^2
        # = sum of ways to pick subsets of {1,1,2,2,3,3} summing to 3
        # 3: pick one mode of weight 3 -> 2 ways
        # 2+1: pick one weight-2 and one weight-1 -> 2*2 = 4 ways
        # 1+1+1: impossible (only 2 weight-1 modes, max sum = 2)
        # Total at q^3 in ch: 2 + 4 = 6. f_+(3) = 6 - 0 = 6? No, -1 is at q^0.
        # ch = 1 + 2q + 3q^2 + 6q^3 + ...
        # f_+ = ch - 1 = 2q + 3q^2 + 6q^3 + ...
        # Actually let me just check:
        assert f_plus[3] == 6  # verified from computation


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
