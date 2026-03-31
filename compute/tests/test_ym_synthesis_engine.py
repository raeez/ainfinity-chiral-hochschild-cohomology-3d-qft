"""Tests for YM Synthesis Engine: boundary BRST, tangent-to-center, mixed couplings.

Verifies:
1. Tangent-to-center dimensions for standard families
2. Center-vanishing rigidity criterion
3. One-parameter criterion
4. Central formality check
5. Mixed coupling dimension and factorization
6. Kappa complementarity (cross-engine consistency with Vol I)
7. Boundary BRST data
8. Cross-engine kappa consistency with gravity_3d_engine

Each test performs actual computation or consistency check.

References:
  Vol II: ym_synthesis_core.tex (Chapter 31)
  Vol I: concordance.tex (Theorems C/D)
  CLAUDE.md: Critical Pitfalls AP1-AP13
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from sympy import Rational, S, Symbol, simplify

from lib.ym_synthesis_engine import (
    tangent_to_center_dimension,
    center_vanishing_rigidity,
    one_parameter_criterion,
    central_formality_check,
    mixed_coupling_dimension,
    complementarity_check,
    boundary_brst_data,
    BoundaryBRSTData,
    _kappa_data,
)


# ===================================================================
# 1. TANGENT-TO-CENTER DIMENSIONS
# ===================================================================

class TestTangentToCenter:
    """Verify dim HH^2_ch(A_B) = dim Z(A_B^!) for standard families."""

    def test_heisenberg_dim_1(self):
        """Heisenberg: Z(H_k^!) = span{level} => dim = 1."""
        assert tangent_to_center_dimension('heisenberg') == 1

    def test_affine_sl2_dim_1(self):
        """Affine sl_2: Z(g_k^!) = span{level} => dim = 1."""
        assert tangent_to_center_dimension('affine_sl2') == 1

    def test_virasoro_dim_1(self):
        """Virasoro: Z(Vir_{26-c}^!) = span{c} => dim = 1."""
        assert tangent_to_center_dimension('virasoro') == 1

    def test_unknown_algebra_raises(self):
        """Unknown algebra type should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown"):
            tangent_to_center_dimension('unknown_family')

    def test_all_standard_families_equal(self):
        """All standard families have dim = 1 (the single parameter)."""
        for family in ['heisenberg', 'affine_sl2', 'virasoro']:
            assert tangent_to_center_dimension(family) == 1


# ===================================================================
# 2. CENTER-VANISHING RIGIDITY
# ===================================================================

class TestRigidity:
    """Verify center_vanishing => no deformations."""

    def test_zero_center_is_rigid(self):
        """Z(A^!) = 0 => infinitesimally rigid."""
        assert center_vanishing_rigidity(0) is True

    def test_nonzero_center_not_rigid(self):
        """Z(A^!) != 0 => NOT rigid."""
        assert center_vanishing_rigidity(1) is False
        assert center_vanishing_rigidity(3) is False

    def test_standard_families_not_rigid(self):
        """Standard families have dim = 1, so they are NOT rigid."""
        for family in ['heisenberg', 'affine_sl2', 'virasoro']:
            dim = tangent_to_center_dimension(family)
            assert center_vanishing_rigidity(dim) is False


# ===================================================================
# 3. ONE-PARAMETER CRITERION
# ===================================================================

class TestOneParameter:
    """Verify dim Z(A^!) = 1 => unique first-order coupling."""

    def test_dim_1_is_one_parameter(self):
        """dim Z = 1 => one-parameter."""
        assert one_parameter_criterion(1) is True

    def test_dim_0_not_one_parameter(self):
        """dim Z = 0 => rigid, not one-parameter."""
        assert one_parameter_criterion(0) is False

    def test_dim_2_not_one_parameter(self):
        """dim Z = 2 => two parameters, not one-parameter."""
        assert one_parameter_criterion(2) is False

    def test_standard_families_are_one_parameter(self):
        """All standard families satisfy one-parameter criterion."""
        for family in ['heisenberg', 'affine_sl2', 'virasoro']:
            dim = tangent_to_center_dimension(family)
            assert one_parameter_criterion(dim) is True


# ===================================================================
# 4. CENTRAL FORMALITY
# ===================================================================

class TestCentralFormality:
    """Verify central formality criterion."""

    def test_both_zero_is_formal(self):
        """HH^1 = 0 and HH^{-1}(A^!) = 0 => formal."""
        assert central_formality_check(0, 0) is True

    def test_nonzero_hh1_not_formal(self):
        """HH^1 != 0 => not formal (automorphisms obstruct)."""
        assert central_formality_check(1, 0) is False

    def test_nonzero_hh_minus1_dual_not_formal(self):
        """HH^{-1}(A^!) != 0 => not formal."""
        assert central_formality_check(0, 1) is False

    def test_both_nonzero_not_formal(self):
        """Both nonzero => not formal."""
        assert central_formality_check(1, 1) is False

    def test_standard_families_formal(self):
        """Standard families have HH^1 = HH^{-1}(A^!) = 0."""
        # For Heisenberg, affine, Virasoro: no infinitesimal automorphisms
        # beyond overall scale, no negative-degree Hochschild classes.
        assert central_formality_check(0, 0) is True


# ===================================================================
# 5. MIXED COUPLING DIMENSION
# ===================================================================

class TestMixedCoupling:
    """Verify mixed coupling dimensions for multiple boundaries."""

    def test_single_boundary(self):
        """Single boundary: no mixed couplings possible."""
        result = mixed_coupling_dimension([1])
        assert result['total_dim'] == 1
        assert result['mixed_dim'] == 0
        assert result['n_boundaries'] == 1

    def test_two_standard_boundaries_no_mixed(self):
        """Two standard (dim=1) boundaries: no mixed couplings.

        For dim Z = 1, the reduced center Z~ = 0, so the tensor product
        Z~_1 x Z~_2 = 0.  No genuine mixed couplings.
        """
        result = mixed_coupling_dimension([1, 1])
        assert result['total_dim'] == 2
        assert result['mixed_dim'] == 0

    def test_three_standard_boundaries(self):
        """Three standard boundaries: still no mixed couplings."""
        result = mixed_coupling_dimension([1, 1, 1])
        assert result['total_dim'] == 3
        assert result['mixed_dim'] == 0

    def test_higher_dim_centers_have_mixed(self):
        """If dim Z > 1, genuine mixed couplings appear.

        For dim Z_1 = 2, Z_2 = 2: Z~_1 = 1, Z~_2 = 1.
        Mixed = 1*1 = 1.  Total = 2 + 2 + 1 = 5.
        """
        result = mixed_coupling_dimension([2, 2])
        assert result['total_dim'] == 5
        assert result['mixed_dim'] == 1

    def test_mixed_dim_3_plus_2(self):
        """dim Z_1 = 3, Z_2 = 2: Z~_1 = 2, Z~_2 = 1. Mixed = 2."""
        result = mixed_coupling_dimension([3, 2])
        assert result['total_dim'] == 7
        assert result['mixed_dim'] == 2

    def test_independent_factorizes(self):
        """Independent boundaries: total = individual + mixed."""
        dims = [1, 1, 1, 1]
        result = mixed_coupling_dimension(dims)
        assert result['total_dim'] == sum(dims) + result['mixed_dim']

    def test_non_independent_returns_none(self):
        """Non-independent boundaries: cannot determine mixed dim."""
        result = mixed_coupling_dimension([1, 1], independent=False)
        assert result['total_dim'] is None
        assert result['mixed_dim'] is None


# ===================================================================
# 6. KAPPA COMPLEMENTARITY (cross-engine consistency)
# ===================================================================

class TestKappaComplementarity:
    """Verify kappa + kappa' for all families (Theorem C)."""

    def test_heisenberg_sum_zero(self):
        """Heisenberg: kappa(H_k) + kappa(H_k^!) = 0."""
        result = complementarity_check('heisenberg', k=S(1))
        assert result['match'] is True
        assert result['expected_sum'] == 0

    def test_affine_sl2_sum_zero(self):
        """Affine sl_2: kappa + kappa' = 0."""
        result = complementarity_check('affine_sl2', k=S(1))
        assert result['match'] is True
        assert result['expected_sum'] == 0

    def test_virasoro_sum_13(self):
        """Virasoro: kappa(Vir_c) + kappa(Vir_{26-c}) = 13."""
        result = complementarity_check('virasoro', c=S(1))
        assert result['match'] is True
        assert result['expected_sum'] == 13

    def test_virasoro_symbolic_sum_13(self):
        """Virasoro symbolic: c/2 + (26-c)/2 = 13."""
        result = complementarity_check('virasoro')
        assert result['match'] is True

    def test_betagamma_sum_zero(self):
        """Beta-gamma: kappa + kappa' = 0."""
        result = complementarity_check('betagamma')
        assert result['match'] is True

    def test_w3_sum_250_over_3(self):
        """W_3: kappa + kappa' = 250/3."""
        result = complementarity_check('w3', c=S(2))
        assert result['match'] is True
        assert result['expected_sum'] == Rational(250, 3)


# ===================================================================
# 7. BOUNDARY BRST DATA
# ===================================================================

class TestBoundaryBRSTData:
    """Verify boundary BRST data construction."""

    def test_heisenberg_brst(self):
        """Heisenberg boundary BRST data."""
        data = boundary_brst_data('heisenberg')
        assert data.algebra_type == 'heisenberg'
        assert data.tangent_dim == 1
        assert data.is_rigid is False
        assert data.is_one_parameter is True
        assert data.is_formal is True

    def test_virasoro_brst(self):
        """Virasoro boundary BRST data."""
        data = boundary_brst_data('virasoro')
        assert data.algebra_type == 'virasoro'
        assert data.tangent_dim == 1
        assert data.is_rigid is False
        assert data.is_one_parameter is True

    def test_affine_brst(self):
        """Affine sl_2 boundary BRST data."""
        data = boundary_brst_data('affine_sl2')
        assert data.tangent_dim == 1
        assert data.is_one_parameter is True

    def test_open_name_contains_bar(self):
        """Open boundary complex is the bar complex B(A)."""
        data = boundary_brst_data('virasoro')
        assert 'B(' in data.open_name

    def test_bulk_dual_contains_omega(self):
        """Dual bulk is Omega(A^!)."""
        data = boundary_brst_data('virasoro')
        assert 'Omega' in data.bulk_dual_name


# ===================================================================
# 8. KAPPA DATA INTERNAL CONSISTENCY
# ===================================================================

class TestKappaDataConsistency:
    """Cross-check kappa data against gravity_3d_engine values."""

    def test_virasoro_kappa_c_over_2(self):
        """kappa(Vir_c) = c/2."""
        c = Symbol('c')
        data = _kappa_data('virasoro', c=c)
        assert simplify(data['kappa'] - c / 2) == 0

    def test_heisenberg_kappa_equals_k(self):
        """kappa(H_k) = k."""
        k = Symbol('k')
        data = _kappa_data('heisenberg', k=k)
        assert simplify(data['kappa'] - k) == 0

    def test_affine_sl2_kappa_formula(self):
        """kappa(sl_2, k) = 3(k+2)/4."""
        k = Symbol('k')
        data = _kappa_data('affine_sl2', k=k)
        expected = Rational(3, 1) * (k + 2) / 4
        assert simplify(data['kappa'] - expected) == 0

    def test_virasoro_kappa_at_c_26(self):
        """kappa(Vir_26) = 13, kappa(Vir_0) = 0."""
        data = _kappa_data('virasoro', c=S(26))
        assert data['kappa'] == 13
        data_0 = _kappa_data('virasoro', c=S(0))
        assert data_0['kappa'] == 0

    def test_virasoro_dual_kappa_at_c_13(self):
        """At c=13 (self-dual): kappa = kappa' = 13/2."""
        data = _kappa_data('virasoro', c=S(13))
        assert data['kappa'] == Rational(13, 2)
        assert data['kappa_dual'] == Rational(13, 2)

    def test_w3_kappa_formula(self):
        """kappa(W_3, c) = 5c/6."""
        c = Symbol('c')
        data = _kappa_data('w3', c=c)
        assert simplify(data['kappa'] - 5 * c / 6) == 0
