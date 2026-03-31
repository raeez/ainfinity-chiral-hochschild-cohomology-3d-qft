"""Tests for Celestial Boundary Transfer Engine.

Verifies the celestial boundary transfer computations:
1. Homotopy transfer tree count (Catalan numbers)
2. Single-particle reduction in filtration quotient
3. Obstruction class degrees in H^1(g, delta_0)
4. Gauge change structure at each order
5. Airy-Witt operators D_m and Witt algebra commutation

References:
  Vol II: celestial_boundary_transfer_core.tex (Part III)
  Vol I: concordance.tex (Theorem A, homotopy transfer)
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from fractions import Fraction
from sympy import Symbol, Rational, S, simplify, expand, symbols, Matrix, zeros

from lib.celestial_boundary_transfer_engine import (
    homotopy_transfer_tree_count,
    single_particle_reduction,
    obstruction_class_degree,
    gauge_change_order_r,
    airy_witt_operator,
    witt_commutation_check,
)


# ===================================================================
# 1. HOMOTOPY TRANSFER TREE COUNT (Catalan numbers)
# ===================================================================

class TestTreeCount:
    """Verify tree counts = Catalan numbers C_{n-1}."""

    def test_arity_1(self):
        """C_0 = 1: a single unary operation."""
        assert homotopy_transfer_tree_count(1) == 1

    def test_arity_2(self):
        """C_1 = 1: one binary tree with 2 leaves."""
        assert homotopy_transfer_tree_count(2) == 1

    def test_arity_3(self):
        """C_2 = 2: two binary trees with 3 leaves."""
        assert homotopy_transfer_tree_count(3) == 2

    def test_arity_4(self):
        """C_3 = 5: five binary trees with 4 leaves."""
        assert homotopy_transfer_tree_count(4) == 5

    def test_arity_5(self):
        """C_4 = 14: fourteen binary trees with 5 leaves."""
        assert homotopy_transfer_tree_count(5) == 14

    def test_arity_6(self):
        """C_5 = 42."""
        assert homotopy_transfer_tree_count(6) == 42

    def test_arity_7(self):
        """C_6 = 132."""
        assert homotopy_transfer_tree_count(7) == 132

    def test_catalan_recurrence(self):
        """Catalan numbers satisfy C_{n+1} = sum_{i=0}^{n} C_i C_{n-i}."""
        for n in range(1, 8):
            c = [homotopy_transfer_tree_count(k + 1) for k in range(n + 2)]
            recurrence_sum = sum(c[i] * c[n - i] for i in range(n + 1))
            assert c[n + 1] == recurrence_sum, f"Catalan recurrence fails at n={n}"

    def test_invalid_arity(self):
        """Arity < 1 raises ValueError."""
        with pytest.raises(ValueError):
            homotopy_transfer_tree_count(0)


# ===================================================================
# 2. SINGLE-PARTICLE REDUCTION
# ===================================================================

class TestSingleParticle:
    """Verify single-particle reduction: higher ops vanish mod F^1."""

    def test_free_field_passes(self):
        """Free field: only m_2 nonzero, higher ops vanish."""
        ops = {2: 1, 3: 0, 4: 0, 5: 0}
        assert single_particle_reduction(ops) is True

    def test_interacting_fails(self):
        """Interacting theory: m_3 nonzero, fails single-particle."""
        ops = {2: 1, 3: 1, 4: 0}
        assert single_particle_reduction(ops) is False

    def test_empty_ops(self):
        """No operations at all: vacuously true."""
        assert single_particle_reduction({}) is True

    def test_only_binary(self):
        """Only binary operation: passes."""
        assert single_particle_reduction({2: 42}) is True

    def test_quartic_nonzero_fails(self):
        """Quartic m_4 nonzero while m_3 = 0: still fails."""
        ops = {2: 1, 3: 0, 4: 1}
        assert single_particle_reduction(ops) is False

    def test_invalid_filtration(self):
        """Filtration < 1 raises ValueError."""
        with pytest.raises(ValueError):
            single_particle_reduction({2: 1}, filtration=0)


# ===================================================================
# 3. OBSTRUCTION CLASS DEGREE
# ===================================================================

class TestObstruction:
    """Verify obstruction class degree structure."""

    def test_cohomological_degree_always_1(self):
        """Obstruction always in H^1."""
        for r in range(2, 10):
            result = obstruction_class_degree(r)
            assert result['cohomological_degree'] == 1

    def test_internal_degree(self):
        """Internal degree = r - 1."""
        for r in range(2, 10):
            result = obstruction_class_degree(r)
            assert result['internal_degree'] == r - 1

    def test_total_degree(self):
        """Total degree = r."""
        for r in range(2, 10):
            result = obstruction_class_degree(r)
            assert result['total_degree'] == r

    def test_target_arity(self):
        """Target arity = r + 1."""
        result = obstruction_class_degree(5)
        assert result['target_arity'] == 6

    def test_invalid_arity(self):
        """Arity < 2 raises ValueError."""
        with pytest.raises(ValueError):
            obstruction_class_degree(1)


# ===================================================================
# 4. GAUGE CHANGE STRUCTURE
# ===================================================================

class TestGaugeChange:
    """Verify gauge change structure at order r."""

    def test_order_stored(self):
        """Order is stored correctly."""
        for r in range(2, 8):
            result = gauge_change_order_r(r)
            assert result['order'] == r

    def test_gauge_parameter_degree(self):
        """Gauge parameter theta_r has degree r - 2."""
        for r in range(2, 8):
            result = gauge_change_order_r(r)
            assert result['gauge_parameter_degree'] == r - 2

    def test_residual_class(self):
        """Residual class is always H^1(g, delta_0)."""
        result = gauge_change_order_r(3)
        assert result['residual_class'] == 'H^1(g, delta_0)'

    def test_gauge_freedom_dim(self):
        """Gauge freedom dimension = r - 1."""
        for r in range(2, 8):
            result = gauge_change_order_r(r)
            assert result['gauge_freedom_dim'] == r - 1

    def test_invalid_order(self):
        """Order < 2 raises ValueError."""
        with pytest.raises(ValueError):
            gauge_change_order_r(1)


# ===================================================================
# 5. AIRY-WITT OPERATORS
# ===================================================================

class TestAiryWitt:
    """Verify Witt algebra L_m = -t^{m+1} d/dt."""

    def test_L0_matrix(self):
        """L_0 = -t d/dt: maps t^k to -k t^k (Euler operator).

        L_0(1) = 0, L_0(t) = -t, L_0(t^2) = -2t^2, ...
        """
        L0 = airy_witt_operator(0, 4)
        # L_0(t^0) = 0 (kills constants)
        assert L0[0, 0] == 0
        # L_0(t^1) = -1 * t^1
        assert L0[1, 1] == -1
        # L_0(t^2) = -2 * t^2
        assert L0[2, 2] == -2
        # L_0(t^3) = -3 * t^3
        assert L0[3, 3] == -3

    def test_L1_matrix(self):
        """L_1 = -t^2 d/dt: maps t^k to -k t^{k+1}.

        L_1(1) = 0, L_1(t) = -t^2, L_1(t^2) = -2t^3, ...
        """
        L1 = airy_witt_operator(1, 5)
        # L_1(t^0) = 0 (kills constants)
        assert L1[1, 0] == 0
        # L_1(t^1) = -1 * t^2
        assert L1[2, 1] == -1
        # L_1(t^2) = -2 * t^3
        assert L1[3, 2] == -2
        # L_1(t^3) = -3 * t^4
        assert L1[4, 3] == -3

    def test_L0_is_diagonal(self):
        """L_0 is diagonal: L_0(t^k) = -k t^k."""
        L0 = airy_witt_operator(0, 5)
        for i in range(6):
            for j in range(6):
                if i == j:
                    assert L0[i, j] == -i  # diagonal entry is -k
                else:
                    assert L0[i, j] == 0

    def test_matrix_size(self):
        """Matrix has correct dimensions."""
        L = airy_witt_operator(0, 7)
        assert L.shape == (8, 8)

    def test_witt_commutation_L0_L1(self):
        """[L_0, L_1] = (0-1) L_1 = -L_1."""
        result = witt_commutation_check(0, 1, 10)
        assert result['is_zero'] is True

    def test_witt_commutation_L0_Lm1(self):
        """[L_0, L_{-1}] = (0-(-1)) L_{-1} = L_{-1}."""
        result = witt_commutation_check(0, -1, 10)
        assert result['is_zero'] is True

    def test_witt_commutation_L1_L2(self):
        """[L_1, L_2] = (1-2) L_3 = -L_3."""
        result = witt_commutation_check(1, 2, 15)
        assert result['is_zero'] is True

    def test_witt_commutation_symmetric(self):
        """[L_m, L_n] = -[L_n, L_m] (antisymmetry of commutator)."""
        for m, n in [(0, 1), (1, 2), (-1, 1), (0, 2)]:
            r1 = witt_commutation_check(m, n, 12)
            r2 = witt_commutation_check(n, m, 12)
            diff = r1['commutator'] + r2['commutator']
            assert diff.equals(zeros(*diff.shape))

    def test_invalid_max_degree(self):
        """max_degree < 0 raises ValueError."""
        with pytest.raises(ValueError):
            airy_witt_operator(0, -1)
