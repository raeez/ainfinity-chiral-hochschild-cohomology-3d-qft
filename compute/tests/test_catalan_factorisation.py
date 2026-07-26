r"""Tests for the Catalan factorisation theorem (Theorem thm:period-2-parity).

Verifies:
  (a) phi_k(x) = 0 for even k >= 4
  (b) phi_k(x) = (-1)^n C_n prod_{m=2}^k (x+m) for odd k = 2n+3
  (c) T_k = (-1)^n C_n k! for odd k >= 3
  (d) Sigma_k^fld = phi_k(1) = (-1)^n C_n (k+1)!/2 for odd k >= 3
  (e) The root property: phi_j(-m) = 0 for m = 2,...,j
  (f) The polynomial recursion via rightmost compositions
  (g) Even-arity vanishing via the functional equation
  (h) P_{2r+1}(1,...,1) = (-1)^{r-1} C_{r-1} (2r+1)!/2
"""

import sys
import os
import math
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from compute.m7_m10_depth_frontier import StasheffEngine


def catalan(n):
    """Catalan number C_n = binom(2n,n)/(n+1)."""
    return math.comb(2 * n, n) // (n + 1)


def get_field_polynomial(engine, k):
    """Extract the field polynomial phi_k(x) as a list of coefficients [x^0, x^1, ...]."""
    engine._cache.clear()
    lams = tuple(1.0 for _ in range(k - 1))
    result = engine.mk(lams)
    field_items = {w: v for w, v in result.items() if w >= 0}
    if not field_items:
        return [0.0]
    max_w = max(field_items.keys())
    return [field_items.get(w, 0.0) for w in range(max_w + 1)]


def eval_polynomial(coeffs, x):
    """Evaluate polynomial sum_w coeffs[w] x^w."""
    return sum(c * x ** w for w, c in enumerate(coeffs))


def predicted_phi(k):
    """Predicted phi_k as coefficient list using the Catalan factorisation."""
    if k == 2:
        return [2.0, 1.0]
    if k % 2 == 0:
        return [0.0]
    n = (k - 3) // 2
    Cn = catalan(n)
    prefactor = (-1) ** n * Cn
    # Multiply out prod_{m=2}^k (x+m) = prod of (m + x)
    poly = [1.0]
    for m in range(2, k + 1):
        # Multiply by (m + x): new[i] = m*old[i] + old[i-1]
        new_poly = [0.0] * (len(poly) + 1)
        for i in range(len(poly)):
            new_poly[i] += float(m) * poly[i]
            new_poly[i + 1] += poly[i]
        poly = new_poly
    return [prefactor * c for c in poly]


@pytest.fixture
def engine():
    return StasheffEngine(1.0)


class TestEvenArityVanishing:
    """Test that phi_k(x) = 0 for even k >= 4."""

    @pytest.mark.parametrize("k", [4, 6, 8, 10, 12])
    def test_even_vanishing(self, engine, k):
        coeffs = get_field_polynomial(engine, k)
        max_abs = max(abs(c) for c in coeffs)
        assert max_abs < 1e-6, f"phi_{k} not zero: max |coeff| = {max_abs}"


class TestOddArityCatalanFactorisation:
    """Test the full polynomial identity phi_k(x) = (-1)^n C_n prod(x+m)."""

    @pytest.mark.parametrize("k", [3, 5, 7, 9, 11, 13])
    def test_polynomial_match(self, engine, k):
        computed = get_field_polynomial(engine, k)
        predicted = predicted_phi(k)
        # Pad to same length
        maxlen = max(len(computed), len(predicted))
        computed += [0.0] * (maxlen - len(computed))
        predicted += [0.0] * (maxlen - len(predicted))
        for w in range(maxlen):
            assert abs(computed[w] - predicted[w]) < 1e-4 * max(abs(predicted[w]), 1), \
                f"phi_{k} x^{w}: computed={computed[w]}, predicted={predicted[w]}"


class TestTCoefficient:
    """Test T_k = (-1)^n C_n k!."""

    @pytest.mark.parametrize("k", [3, 5, 7, 9, 11, 13])
    def test_T_coefficient(self, engine, k):
        coeffs = get_field_polynomial(engine, k)
        T_val = coeffs[0]
        n = (k - 3) // 2
        predicted = (-1) ** n * catalan(n) * math.factorial(k)
        assert abs(T_val - predicted) < 1, \
            f"T_{k} = {T_val}, predicted {predicted}"


class TestFieldSignedSum:
    """Test Sigma_k^fld = phi_k(1) = (-1)^n C_n (k+1)!/2."""

    @pytest.mark.parametrize("k", [3, 5, 7, 9, 11, 13])
    def test_signed_sum(self, engine, k):
        coeffs = get_field_polynomial(engine, k)
        S_val = eval_polynomial(coeffs, 1.0)
        n = (k - 3) // 2
        predicted = (-1) ** n * catalan(n) * math.factorial(k + 1) // 2
        assert abs(S_val - predicted) < 1, \
            f"S_{k} = {S_val}, predicted {predicted}"


class TestRootProperty:
    """Test that phi_j(-m) = 0 for m = 2, ..., j."""

    @pytest.mark.parametrize("j", [2, 3, 5, 7, 9, 11, 13])
    def test_roots(self, engine, j):
        coeffs = get_field_polynomial(engine, j)
        tol = 1e-10 * max(abs(coeffs[0]), 1.0)
        for m in range(2, j + 1):
            val = eval_polynomial(coeffs, float(-m))
            assert abs(val) < tol, \
                f"phi_{j}({-m}) = {val}, should be 0"


class TestCatalanConvolution:
    """Test the Catalan convolution C_n = sum_{a=0}^{n-1} C_a C_{n-1-a}."""

    @pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
    def test_convolution(self, n):
        conv = sum(catalan(a) * catalan(n - 1 - a) for a in range(n))
        assert conv == catalan(n), f"Convolution = {conv}, C_{n} = {catalan(n)}"


class TestFunctionalEquation:
    """Test the even-vanishing functional equation:
    (x+2) phi_{k-1}(x+1) = (x+k) phi_{k-1}(x) for even k >= 4."""

    @pytest.mark.parametrize("k", [4, 6, 8, 10])
    def test_functional_equation(self, k):
        pred = predicted_phi(k - 1)
        # Evaluate at several x values
        for x in [-5.5, -1.3, 0, 0.7, 2.5, 10]:
            lhs = (x + 2) * eval_polynomial(pred, x + 1)
            rhs = (x + k) * eval_polynomial(pred, x)
            assert abs(lhs - rhs) < 1e-4 * max(abs(rhs), 1), \
                f"k={k}, x={x}: LHS={lhs}, RHS={rhs}"


class TestScalarCatalanFormula:
    """Test the scalar Catalan formula independently of the field signed sum."""

    @pytest.mark.parametrize("k", [3, 5, 7, 9, 11, 13])
    def test_scalar_polynomial_at_unit(self, k):
        engine = StasheffEngine(1.0)
        engine._cache.clear()
        result = engine.mk(tuple(1.0 for _ in range(k - 1)))
        r = (k - 1) // 2
        P_val = 12.0 * result.get(-1, 0.0)
        predicted = (-1) ** (r - 1) * catalan(r - 1) * math.factorial(k) / 2.0
        assert abs(P_val - predicted) < 1e-8 * max(abs(predicted), 1.0), \
            f"P_{k}(1,...,1) = {P_val}, predicted {predicted}"

    @pytest.mark.parametrize("k", [3, 5, 7, 9, 11, 13])
    def test_scalar_t_proportionality_at_unit(self, k):
        c_val = 7.0
        engine = StasheffEngine(c_val)
        engine._cache.clear()
        result = engine.mk(tuple(1.0 for _ in range(k - 1)))
        scalar = result.get(-1, 0.0)
        T_coeff = result.get(0, 0.0)
        assert abs(scalar / T_coeff - c_val / 24.0) < 1e-8


class TestPeriod2ParityManuscriptGuard:
    """Guard against the previous scalar/field conflation in the theorem text."""

    @staticmethod
    def _gravity_source():
        root = Path(__file__).resolve().parents[2]
        return (root / "chapters" / "connections" / "3d_gravity.tex").read_text()

    def test_theorem_status_licensing_and_depth_range(self):
        source = self._gravity_source()
        theorem_start = source.index(r"\begin{theorem}[Symmetric-point period-$2$")
        theorem_end = source.index(r"\end{theorem}", theorem_start)
        theorem = source[theorem_start:theorem_end]

        assert r"\ClaimStatusConditional" in theorem
        assert r"\hypAmbientWtCpl+\effKoszul" in theorem
        assert r"d \in \{0, 1, \ldots, k-1\}" in theorem
        assert r"\Sigma_k^{\mathrm{fld}}" in theorem
        assert r"\label{eq:symmetric-scalar-t-proportionality}" in theorem

    def test_no_field_signed_sum_used_as_scalar_shadow(self):
        source = self._gravity_source()
        assert r"S_k = \varphi_k(1)" not in source
        assert r"P_k = 12 S_k / c" not in source

    def test_graviton_tracelessness_uses_period2_field_polynomial(self):
        source = self._gravity_source()
        prop_start = source.index(
            r"\begin{proposition}[Symmetric-point graviton tracelessness;"
        )
        prop_end = source.index(
            r"\begin{theorem}[Scalar shadow resolvent on the Virasoro metric branch",
            prop_start,
        )
        prop = source[prop_start:prop_end]

        required = (
            r"\ClaimStatusConditional",
            r"\hypAmbientWtCpl+\effKoszul",
            r"\sum_{j=0}^{k-1} a_j^{(k)}x^j",
            r"\bigl[\partial^jT\bigr]",
            r"Theorem~\ref{thm:period-2-parity}",
            r"G_k(-h_T)=0",
            r"h_T=2",
            r"\prod_{m=2}^{k}(x+m)",
            r"\(x+2\)",
            r"the corresponding field polynomial has the Ward",
        )
        for needle in required:
            assert needle in prop

        retired = (
            r"\sum_{j=0}^{k-2}",
            "By induction on the Stasheff identity",
            r"For a general chiral algebra with generator~$a$",
            r"G_k^{(a)}(-h_a) = 0",
        )
        for needle in retired:
            assert needle not in prop
