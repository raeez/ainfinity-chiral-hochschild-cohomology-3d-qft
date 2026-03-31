"""Tests for DG-Shifted Factorization Bridge Engine.

Verifies the dg-shifted factorization bridge computations:
1. Root multiplicity = 1 for all simple Lie algebras
2. BCH coefficients beta_n = 1/n via direct formula and integral
3. Spectral Drinfeld obstruction vanishes for all simple types
4. Spectral Kohno relation reduces to IB for rational r-matrix
5. Jacobi collapse dimension = 1 for all simple types

References:
  Vol II: dg_shifted_factorization_bridge.tex (Part VI)
  Vol II: spectral-braiding-core.tex (Part III)
  Vol I: concordance.tex (MC3, Drinfeld-Kohno)
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from fractions import Fraction
from sympy import Symbol, Rational, S, simplify

from lib.dg_shifted_factorization_engine import (
    root_multiplicity,
    bch_coefficient,
    bch_coefficient_integral,
    spectral_drinfeld_obstruction_vanishes,
    spectral_kohno_check,
    jacobi_collapse_dimension,
)


# ===================================================================
# 1. ROOT MULTIPLICITY
# ===================================================================

class TestRootMultiplicity:
    """Verify root multiplicity = 1 for all simple Lie algebras."""

    def test_sl_n_type_A(self):
        """sl_n (n=2..6): all root multiplicities = 1."""
        for rank in range(1, 6):
            result = root_multiplicity('A', rank)
            assert result['multiplicity'] == 1

    def test_so_odd_type_B(self):
        """so_{2n+1} (n=2..5): root multiplicity = 1."""
        for rank in range(2, 6):
            result = root_multiplicity('B', rank)
            assert result['multiplicity'] == 1

    def test_sp_type_C(self):
        """sp_{2n} (n=2..5): root multiplicity = 1."""
        for rank in range(2, 6):
            result = root_multiplicity('C', rank)
            assert result['multiplicity'] == 1

    def test_so_even_type_D(self):
        """so_{2n} (n=4..7): root multiplicity = 1."""
        for rank in range(4, 8):
            result = root_multiplicity('D', rank)
            assert result['multiplicity'] == 1

    def test_G2(self):
        """G_2: root multiplicity = 1."""
        result = root_multiplicity('G', 2)
        assert result['multiplicity'] == 1
        assert result['dimension'] == 14

    def test_F4(self):
        """F_4: root multiplicity = 1."""
        result = root_multiplicity('F', 4)
        assert result['multiplicity'] == 1
        assert result['dimension'] == 52

    def test_E6(self):
        """E_6: root multiplicity = 1, dim = 78."""
        result = root_multiplicity('E', 6)
        assert result['multiplicity'] == 1
        assert result['dimension'] == 78

    def test_E7(self):
        """E_7: root multiplicity = 1, dim = 133."""
        result = root_multiplicity('E', 7)
        assert result['multiplicity'] == 1
        assert result['dimension'] == 133

    def test_E8(self):
        """E_8: root multiplicity = 1, dim = 248."""
        result = root_multiplicity('E', 8)
        assert result['multiplicity'] == 1
        assert result['dimension'] == 248

    def test_simply_laced_classification(self):
        """ADE are simply-laced; BCF G are not."""
        assert root_multiplicity('A', 3)['is_simply_laced'] is True
        assert root_multiplicity('D', 4)['is_simply_laced'] is True
        assert root_multiplicity('E', 6)['is_simply_laced'] is True
        assert root_multiplicity('B', 3)['is_simply_laced'] is False
        assert root_multiplicity('C', 3)['is_simply_laced'] is False
        assert root_multiplicity('F', 4)['is_simply_laced'] is False
        assert root_multiplicity('G', 2)['is_simply_laced'] is False

    def test_dimension_sl_n(self):
        """dim(sl_n) = n^2 - 1."""
        for n in range(2, 8):
            result = root_multiplicity('A', n - 1)
            assert result['dimension'] == n ** 2 - 1

    def test_invalid_type(self):
        """Unknown Lie type raises ValueError."""
        with pytest.raises(ValueError):
            root_multiplicity('X', 3)

    def test_invalid_rank(self):
        """Invalid rank for type raises ValueError."""
        with pytest.raises(ValueError):
            root_multiplicity('G', 3)  # G only has rank 2


# ===================================================================
# 2. BCH COEFFICIENT
# ===================================================================

class TestBCH:
    """Verify BCH coefficients beta_n = 1/n."""

    def test_beta_1(self):
        """beta_1 = 1."""
        assert bch_coefficient(1) == Fraction(1, 1)

    def test_beta_2(self):
        """beta_2 = 1/2."""
        assert bch_coefficient(2) == Fraction(1, 2)

    def test_beta_3(self):
        """beta_3 = 1/3."""
        assert bch_coefficient(3) == Fraction(1, 3)

    def test_beta_10(self):
        """beta_10 = 1/10."""
        assert bch_coefficient(10) == Fraction(1, 10)

    def test_beta_100(self):
        """beta_100 = 1/100."""
        assert bch_coefficient(100) == Fraction(1, 100)

    def test_integral_formula(self):
        """Verify via integral: int_0^1 (1-t)^{n-1} dt = 1/n."""
        for n in range(1, 12):
            integral_val = bch_coefficient_integral(n)
            assert integral_val == Rational(1, n)

    def test_integral_matches_direct(self):
        """Integral and direct formula agree."""
        for n in range(1, 8):
            direct = bch_coefficient(n)
            integral = bch_coefficient_integral(n)
            assert Fraction(integral.p, integral.q) == direct

    def test_invalid_level(self):
        """Filtration level < 1 raises ValueError."""
        with pytest.raises(ValueError):
            bch_coefficient(0)


# ===================================================================
# 3. SPECTRAL DRINFELD STRICTIFICATION
# ===================================================================

class TestStrictification:
    """Verify spectral Drinfeld obstruction vanishes for all simple types."""

    def test_vanishes_type_A(self):
        """Vanishes for sl_n (n=2..6)."""
        for rank in range(1, 6):
            result = spectral_drinfeld_obstruction_vanishes('A', rank)
            assert result['vanishes'] is True
            assert result['reason'] == 'root_multiplicity_one'

    def test_vanishes_type_B(self):
        """Vanishes for so_{2n+1}."""
        for rank in range(2, 5):
            result = spectral_drinfeld_obstruction_vanishes('B', rank)
            assert result['vanishes'] is True

    def test_vanishes_type_C(self):
        """Vanishes for sp_{2n}."""
        for rank in range(2, 5):
            result = spectral_drinfeld_obstruction_vanishes('C', rank)
            assert result['vanishes'] is True

    def test_vanishes_type_D(self):
        """Vanishes for so_{2n}."""
        for rank in range(4, 7):
            result = spectral_drinfeld_obstruction_vanishes('D', rank)
            assert result['vanishes'] is True

    def test_vanishes_exceptionals(self):
        """Vanishes for all exceptional types."""
        for lie_type, rank in [('G', 2), ('F', 4), ('E', 6), ('E', 7), ('E', 8)]:
            result = spectral_drinfeld_obstruction_vanishes(lie_type, rank)
            assert result['vanishes'] is True


# ===================================================================
# 4. SPECTRAL KOHNO
# ===================================================================

class TestSpectralKohno:
    """Verify spectral Kohno relation for rational r-matrix."""

    def test_ib_holds_sl2(self):
        """IB relation holds for sl_2 (dim=3)."""
        result = spectral_kohno_check(3)
        assert result['ib_relation_holds'] is True

    def test_reduces_to_jacobi(self):
        """Rational r-matrix r(u) = Omega/u reduces spectral Kohno to IB/Jacobi."""
        result = spectral_kohno_check(3)
        assert result['reduces_to_jacobi'] is True

    def test_rational_r_matrix(self):
        """r-matrix type is rational."""
        result = spectral_kohno_check(3)
        assert result['r_matrix_type'] == 'rational'

    def test_simple_pole(self):
        """Rational r-matrix has simple pole at u=0."""
        result = spectral_kohno_check(3)
        assert result['pole_order'] == 1

    def test_ib_holds_all_dims(self):
        """IB relation holds for all dims (Jacobi identity is universal)."""
        for dim in [3, 8, 14, 52, 78, 133, 248]:
            result = spectral_kohno_check(dim)
            assert result['ib_relation_holds'] is True


# ===================================================================
# 5. JACOBI COLLAPSE DIMENSION
# ===================================================================

class TestJacobiCollapse:
    """Verify Jacobi collapse dimension = 1 for all simple types."""

    def test_collapse_dim_type_A(self):
        """Collapse dimension = 1 for sl_n."""
        for rank in range(1, 6):
            result = jacobi_collapse_dimension('A', rank)
            assert result['collapse_dim'] == 1
            assert result['collapses'] is True

    def test_collapse_dim_type_B(self):
        """Collapse dimension = 1 for so_{2n+1}."""
        for rank in range(2, 5):
            result = jacobi_collapse_dimension('B', rank)
            assert result['collapse_dim'] == 1

    def test_collapse_dim_type_C(self):
        """Collapse dimension = 1 for sp_{2n}."""
        for rank in range(2, 5):
            result = jacobi_collapse_dimension('C', rank)
            assert result['collapse_dim'] == 1

    def test_collapse_dim_type_D(self):
        """Collapse dimension = 1 for so_{2n}."""
        for rank in range(4, 7):
            result = jacobi_collapse_dimension('D', rank)
            assert result['collapse_dim'] == 1

    def test_collapse_dim_exceptionals(self):
        """Collapse dimension = 1 for all exceptional types."""
        for lie_type, rank in [('G', 2), ('F', 4), ('E', 6), ('E', 7), ('E', 8)]:
            result = jacobi_collapse_dimension(lie_type, rank)
            assert result['collapse_dim'] == 1

    def test_obstruction_vanishes(self):
        """Jacobi collapse implies obstruction vanishes."""
        for lie_type, rank in [('A', 3), ('B', 3), ('C', 3), ('D', 4), ('E', 6)]:
            result = jacobi_collapse_dimension(lie_type, rank)
            assert result['obstruction_vanishes'] is True

    def test_root_mult_stored(self):
        """Root multiplicity = 1 is stored for reference."""
        result = jacobi_collapse_dimension('A', 4)
        assert result['root_mult'] == 1
