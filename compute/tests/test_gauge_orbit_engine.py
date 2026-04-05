r"""Tests for gauge orbit unification engine.

Tests the eight numpy-based functions implementing gauge-orbit theorems I-IV
from relative_feynman_transform.tex (subsec:gauge-orbit-canonical), plus
cross-checks against the sympy-based functions in the same engine.

Concrete 4x4 matrix test data:
  D0 = [[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]]  (nilpotent, two Jordan blocks)
  delta chosen to satisfy [D0, delta] + delta^2 = 0 with D0
  Phi = I + small correction (invertible)
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
import pytest

from compute.lib.gauge_orbit_engine import (
    mc_equation_check,
    gauge_conjugation,
    curvature_gauge_invariant,
    bicomplex_from_mc,
    heisenberg_gauge_data,
    virasoro_gauge_data,
    three_models_kappa_agreement,
    clutching_residue_sum,
    # Sympy-based cross-checks
    mc_element_delta,
    gauge_transform_conjugation,
    heisenberg_bar_arity2,
    three_models_same_cohomology,
    standard_kappa_values,
    kappa_duality_check,
)


# =========================================================================
# Standard 4x4 test matrices
# =========================================================================

# D0: nilpotent with two 2x1 Jordan blocks
D0_4x4 = np.array([
    [0, 1, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 0, 0],
], dtype=float)

# delta = 0 trivially satisfies MC equation with any D0
DELTA_ZERO_4x4 = np.zeros((4, 4), dtype=float)

# A nontrivial delta that satisfies MC with D0_4x4.
# We need: D0 @ delta - delta @ D0 + delta^2 = 0.
# Choose delta with entries only in the (0,2) and (2,0) block corners
# so the commutator and square cancel.
# delta = [[0,0,a,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]] => check MC:
#   D0 @ delta = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]] (row 0 of D0 is [0,1,0,0],
#     D0 @ delta: row0 = [0,0,a*0,0] ... let me compute carefully)
# Actually: D0 @ delta - delta @ D0:
#   D0 has (0,1)=1, (2,3)=1. delta has (0,2)=a.
#   (D0 @ delta)[i,j] = sum_k D0[i,k]*delta[k,j]
#   (D0 @ delta)[0,2] = D0[0,1]*delta[1,2] = 0
#   (delta @ D0)[i,j] = sum_k delta[i,k]*D0[k,j]
#   (delta @ D0)[0,3] = delta[0,2]*D0[2,3] = a
#   So commutator has (0,3) = -a. delta^2 = 0 for rank-1. MC residual = -a at (0,3).
# Need something better. Let's use delta with (1,2) entry:
# delta = [[0,0,0,0],[0,0,d,0],[0,0,0,0],[0,0,0,0]]
#   D0 @ delta: (0,j) = D0[0,1]*delta[1,j] => (0,2) = d
#   delta @ D0: (1,j) = delta[1,2]*D0[2,j] => (1,3) = d
#   comm = D0@delta - delta@D0: (0,2) = d, (1,3) = -d.  delta^2 = 0.
#   MC = [[0,0,d,0],[0,0,0,-d],[0,0,0,0],[0,0,0,0]] != 0.
#
# For a true MC solution with D0_4x4, the simplest is delta = 0.
# For a nontrivial example, use a strictly upper-triangular delta
# with matching structure. Take:
#   delta = [[0,0,0,b],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
#   D0 @ delta: all zero (D0[0,1]*delta[1,j]=0 for j=3 since delta[1,3]=0)
#   Wait: delta[0,3]=b. D0@delta[0,3] = D0[0,1]*delta[1,3] = 0.
#   D0@delta = 0 entirely? D0[i,k]*delta[k,j]: only D0[0,1]=1, D0[2,3]=1 nonzero.
#     D0@delta[0,j] = delta[1,j] = 0 for all j. D0@delta[2,j] = delta[3,j] = 0.
#   delta@D0[i,j] = delta[i,k]*D0[k,j]. delta[0,3]*D0[3,j] = 0 (row 3 of D0 is zero).
#   So D0@delta - delta@D0 = 0. And delta^2 = 0. MC holds for any b. Good.

DELTA_NONTRIVIAL_4x4 = np.array([
    [0, 0, 0, 0.5],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
], dtype=float)

# Phi: invertible, identity + small correction mixing both blocks
PHI_4x4 = np.array([
    [1, 0.1, 0.2, 0],
    [0, 1,   0,   0.1],
    [0, 0,   1,   0.1],
    [0, 0,   0,   1],
], dtype=float)

# Bicomplex matrices.
# D_P = D0_4x4 (two Jordan blocks, D_P^2 = 0).
# D_Mod at position (0,3): maps block 1 to block 2, D_Mod^2 = 0,
# and {D_P, D_Mod} = 0 because D_P has no row/column overlap with (0,3).
D_P_4x4 = D0_4x4.copy()

D_Mod_bicomplex_4x4 = np.array([
    [0, 0, 0, 1],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
], dtype=float)


# =========================================================================
# 1. mc_equation_check
# =========================================================================

class TestMCEquationCheck:
    """Tests for mc_equation_check(D0, delta)."""

    def test_zero_delta(self):
        """delta = 0 trivially satisfies MC with any D0."""
        result = mc_equation_check(D0_4x4, DELTA_ZERO_4x4)
        assert result['holds'] is True
        assert result['residual'] < 1e-14

    def test_nontrivial_delta(self):
        """delta with (0,3) entry satisfies MC with block-diagonal D0."""
        result = mc_equation_check(D0_4x4, DELTA_NONTRIVIAL_4x4)
        assert result['holds'] is True
        assert result['residual'] < 1e-12

    def test_random_delta_fails(self):
        """A random delta generically does NOT satisfy MC."""
        rng = np.random.default_rng(42)
        delta_rand = rng.standard_normal((4, 4)) * 0.1
        result = mc_equation_check(D0_4x4, delta_rand)
        assert result['holds'] is False
        assert result['residual'] > 1e-6

    def test_2x2_nilpotent(self):
        """Simple 2x2 nilpotent D0 with zero delta."""
        D0 = np.array([[0, 1], [0, 0]], dtype=float)
        delta = np.zeros((2, 2), dtype=float)
        result = mc_equation_check(D0, delta)
        assert result['holds'] is True

    def test_identity_D0_forces_delta_sq_zero(self):
        """If D0 = 0, MC reduces to delta^2 = 0."""
        D0 = np.zeros((3, 3), dtype=float)
        # Nilpotent delta: strictly upper triangular
        delta = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]], dtype=float)
        # delta^2 = [[0,0,1],[0,0,0],[0,0,0]], nonzero
        result = mc_equation_check(D0, delta)
        assert result['holds'] is False

    def test_D0_zero_with_rank1_delta(self):
        """D0=0, delta nilpotent rank 1 => delta^2=0 => MC holds."""
        D0 = np.zeros((3, 3), dtype=float)
        delta = np.array([[0, 0, 0], [0, 0, 1], [0, 0, 0]], dtype=float)
        result = mc_equation_check(D0, delta)
        assert result['holds'] is True


# =========================================================================
# 2. gauge_conjugation
# =========================================================================

class TestGaugeConjugation:
    """Tests for gauge_conjugation(D_total, Phi)."""

    def test_identity_gauge(self):
        """Phi = I leaves D_total unchanged."""
        D_total = D0_4x4 + DELTA_NONTRIVIAL_4x4
        I = np.eye(4)
        result = gauge_conjugation(D_total, I)
        assert np.allclose(result, D_total)

    def test_conjugated_D_squared_preserved(self):
        """D^2 is preserved under gauge conjugation (similarity)."""
        D_total = D0_4x4.copy()  # D0^2 = 0
        D_conj = gauge_conjugation(D_total, PHI_4x4)
        assert np.allclose(D_conj @ D_conj, np.zeros((4, 4)), atol=1e-12)

    def test_nontrivial_conjugation(self):
        """Conjugation by non-identity Phi changes the matrix."""
        D_total = D0_4x4 + DELTA_NONTRIVIAL_4x4
        D_conj = gauge_conjugation(D_total, PHI_4x4)
        # Should not be equal to D_total generically
        assert not np.allclose(D_conj, D_total)

    def test_double_conjugation_inverse(self):
        """Conjugating by Phi then by Phi^{-1} returns the original."""
        D_total = D0_4x4 + DELTA_NONTRIVIAL_4x4
        D_conj = gauge_conjugation(D_total, PHI_4x4)
        Phi_inv = np.linalg.inv(PHI_4x4)
        D_back = gauge_conjugation(D_conj, Phi_inv)
        assert np.allclose(D_back, D_total, atol=1e-12)


# =========================================================================
# 3. curvature_gauge_invariant
# =========================================================================

class TestCurvatureGaugeInvariant:
    """Tests for curvature_gauge_invariant(d_fib, Phi)."""

    def test_nilpotent_d_fib(self):
        """Nilpotent d_fib => d_fib^2 = 0, invariant under any gauge."""
        assert curvature_gauge_invariant(D0_4x4, PHI_4x4) is True

    def test_identity_gauge(self):
        """Phi = I trivially preserves curvature."""
        d_fib = D0_4x4 + 0.3 * DELTA_NONTRIVIAL_4x4
        assert curvature_gauge_invariant(d_fib, np.eye(4)) is True

    def test_nontrivial_curvature(self):
        """d_fib with nonzero d_fib^2 still satisfies the conjugation identity."""
        # d_fib not nilpotent: use a matrix with d_fib^2 != 0
        d_fib = np.array([
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        ], dtype=float)
        # d_fib^2 = [[0,0,1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]] != 0
        assert curvature_gauge_invariant(d_fib, PHI_4x4) is True

    def test_random_phi(self):
        """Random invertible Phi still preserves the curvature identity."""
        rng = np.random.default_rng(123)
        Phi = np.eye(4) + 0.3 * rng.standard_normal((4, 4))
        # Make sure it's invertible (generically true)
        assert abs(np.linalg.det(Phi)) > 0.01
        d_fib = np.array([
            [0, 2, 0, 0],
            [0, 0, 3, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        ], dtype=float)
        assert curvature_gauge_invariant(d_fib, Phi) is True


# =========================================================================
# 4. bicomplex_from_mc
# =========================================================================

class TestBicomplexFromMC:
    """Tests for bicomplex_from_mc(D_P, D_Mod)."""

    def test_valid_bicomplex(self):
        """D_P on superdiag-1, D_Mod on superdiag-3: valid bicomplex."""
        result = bicomplex_from_mc(D_P_4x4, D_Mod_bicomplex_4x4)
        assert result['D_P_squared_zero'] is True
        assert result['D_Mod_squared_zero'] is True
        assert result['anticommutator_zero'] is True
        assert result['total_squared_zero'] is True
        assert result['all_hold'] is True

    def test_invalid_bicomplex(self):
        """D_Mod at (1,2) with block-diagonal D_P: anticommutator nonzero."""
        # D_Mod with entry at (1,2) connects the two Jordan blocks
        # in a way that makes {D_P, D_Mod} != 0.
        D_Mod_bad = np.array([
            [0, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        ], dtype=float)
        result = bicomplex_from_mc(D_P_4x4, D_Mod_bad)
        assert result['D_P_squared_zero'] is True
        assert result['D_Mod_squared_zero'] is True
        # Anticommutator: D_P @ D_Mod + D_Mod @ D_P should be nonzero
        assert result['anticommutator_zero'] is False
        assert result['all_hold'] is False

    def test_both_zero(self):
        """Zero matrices form a trivial bicomplex."""
        Z = np.zeros((4, 4))
        result = bicomplex_from_mc(Z, Z)
        assert result['all_hold'] is True

    def test_total_squared_iff_bicomplex(self):
        """D_total^2 = 0 iff all three bicomplex conditions hold."""
        result = bicomplex_from_mc(D_P_4x4, D_Mod_bicomplex_4x4)
        assert result['total_squared_zero'] == result['all_hold']


# =========================================================================
# 5. heisenberg_gauge_data
# =========================================================================

class TestHeisenbergGaugeData:
    """Tests for heisenberg_gauge_data(k)."""

    def test_level_1(self):
        data = heisenberg_gauge_data(1)
        assert data['kappa'] == 1
        assert data['shadow_depth'] == 2
        assert data['mc_element_scalar'] == 1
        assert data['shadow_class'] == 'G'

    def test_level_minus_half(self):
        data = heisenberg_gauge_data(-0.5)
        assert data['kappa'] == -0.5
        assert data['shadow_depth'] == 2
        assert data['mc_element_scalar'] == -0.5

    def test_kappa_equals_level(self):
        """kappa(H_k) = k for all k."""
        for k in [0, 1, -1, 0.5, 100, -7.3]:
            data = heisenberg_gauge_data(k)
            assert data['kappa'] == k

    def test_gaussian_class(self):
        """Heisenberg is always Gaussian class (shadow depth 2)."""
        for k in [1, 2, -3, 0.1]:
            data = heisenberg_gauge_data(k)
            assert data['shadow_depth'] == 2
            assert data['shadow_class'] == 'G'


# =========================================================================
# 6. virasoro_gauge_data
# =========================================================================

class TestVirasoroGaugeData:
    """Tests for virasoro_gauge_data(c)."""

    def test_c26(self):
        data = virasoro_gauge_data(26)
        assert data['kappa'] == 13.0
        assert data['shadow_depth'] == float('inf')
        assert data['mc_element_scalar'] == 13.0

    def test_c1(self):
        data = virasoro_gauge_data(1)
        assert data['kappa'] == 0.5

    def test_kappa_formula(self):
        """kappa(Vir_c) = c/2 for all c."""
        for c in [0, 1, 13, 26, -2, 0.5]:
            data = virasoro_gauge_data(c)
            assert abs(data['kappa'] - c / 2.0) < 1e-15

    def test_class_M(self):
        """Virasoro is always class M (infinite shadow obstruction tower)."""
        for c in [1, 13, 26, 100]:
            data = virasoro_gauge_data(c)
            assert data['shadow_class'] == 'M'
            assert data['shadow_depth'] == float('inf')


# =========================================================================
# 7. three_models_kappa_agreement
# =========================================================================

class TestThreeModelsKappaAgreement:
    """Tests for three_models_kappa_agreement(algebra_type, params)."""

    def test_heisenberg(self):
        result = three_models_kappa_agreement('heisenberg', {'k': 3})
        assert result['all_agree'] is True
        assert result['kappa_flat'] == 3.0

    def test_virasoro(self):
        result = three_models_kappa_agreement('virasoro', {'c': 26})
        assert result['all_agree'] is True
        assert result['kappa_flat'] == 13.0

    def test_affine_sl2(self):
        """sl_2: dim=3, h^vee=2. At k=1: kappa = 3*(1+2)/(2*2) = 9/4."""
        result = three_models_kappa_agreement('affine', {
            'k': 1, 'dim_g': 3, 'h_vee': 2,
        })
        assert result['all_agree'] is True
        assert abs(result['kappa_flat'] - 9.0 / 4.0) < 1e-14

    def test_betagamma(self):
        result = three_models_kappa_agreement('betagamma', {})
        assert result['all_agree'] is True
        assert result['kappa_flat'] == -1.0

    def test_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown algebra type"):
            three_models_kappa_agreement('nonexistent', {})


# =========================================================================
# 8. clutching_residue_sum
# =========================================================================

class TestClutchingResidueSum:
    """Tests for clutching_residue_sum(edge_weights)."""

    def test_single_edge(self):
        assert clutching_residue_sum([3.0]) == 3.0

    def test_multiple_edges(self):
        assert abs(clutching_residue_sum([1.0, 2.0, 3.5]) - 6.5) < 1e-15

    def test_empty(self):
        assert clutching_residue_sum([]) == 0.0

    def test_negative_weights(self):
        assert abs(clutching_residue_sum([1.0, -1.0]) - 0.0) < 1e-15

    def test_fraction_weights(self):
        """Exact arithmetic with Fraction inputs."""
        weights = [Fraction(1, 3), Fraction(1, 6), Fraction(1, 2)]
        assert clutching_residue_sum(weights) == 1.0


# =========================================================================
# 9. Cross-checks: numpy vs sympy engines
# =========================================================================

class TestCrossEngineConsistency:
    """Cross-check numpy functions against sympy functions."""

    def test_heisenberg_kappa_sympy_vs_numpy(self):
        """heisenberg_bar_arity2 and heisenberg_gauge_data agree on kappa."""
        sympy_data = heisenberg_bar_arity2(k=7)
        numpy_data = heisenberg_gauge_data(7)
        assert float(sympy_data['kappa']) == numpy_data['kappa']
        assert sympy_data['shadow_depth'] == numpy_data['shadow_depth']

    def test_virasoro_kappa_cross(self):
        """standard_kappa_values for Virasoro matches virasoro_gauge_data."""
        kappas = standard_kappa_values()
        # Symbolic: c/2
        from sympy import Symbol as Sym
        c_sym = Sym('c')
        assert kappas['virasoro']['kappa'] == c_sym / 2
        # Numeric
        data = virasoro_gauge_data(26)
        assert data['kappa'] == 13.0

    def test_kappa_duality_virasoro(self):
        """kappa(Vir_c) + kappa(Vir_{26-c}) = 13."""
        for c_val in [0, 1, 13, 26, 5.5]:
            d1 = virasoro_gauge_data(c_val)
            d2 = virasoro_gauge_data(26.0 - c_val)
            assert abs(d1['kappa'] + d2['kappa'] - 13.0) < 1e-14

    def test_mc_sympy_vs_numpy(self):
        """mc_element_delta (sympy) vs mc_equation_check (numpy) agree."""
        from sympy import Matrix as SMatrix, zeros as szeros
        D0_s = SMatrix([[0, 1, 0, 0], [0, 0, 0, 0],
                        [0, 0, 0, 1], [0, 0, 0, 0]])
        D_s = D0_s  # delta = 0
        sympy_result = mc_element_delta(D_s, D0_s)
        assert sympy_result['mc_satisfied'] is True

        numpy_result = mc_equation_check(D0_4x4, DELTA_ZERO_4x4)
        assert numpy_result['holds'] is True
