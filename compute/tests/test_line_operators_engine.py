"""Tests for Line Operators Engine.

Verifies SC^{ch,top} operation space dimensions and structure:
1. Operation space dimensions at small arities
2. No-open-to-closed directionality
3. Homotopy-Koszulity indicators (bar concentration)
4. Cross-engine consistency with sc_bar_cobar_engine

Each test performs ACTUAL computation, not lookup.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from math import factorial, comb

from lib.line_operators_engine import (
    fm_dim,
    e1_dim,
    sc_operation_space_dim,
    sc_betti_dim,
    no_open_to_closed_check,
    homotopy_koszul_bar_concentration,
    cross_engine_bridge,
)


# ===================================================================
# 1. OPERATION SPACE DIMENSIONS
# ===================================================================

class TestOperationSpaces:
    """Test SC operation space dimensions at small arities."""

    def test_fm2_dim(self):
        """FM_2(C) has real dim 2(2-1) = 2."""
        assert fm_dim(2) == 2

    def test_fm3_dim(self):
        """FM_3(C) has real dim 2(3-1) = 4."""
        assert fm_dim(3) == 4

    def test_fm1_dim(self):
        """FM_1(C) = point, dim 0."""
        assert fm_dim(1) == 0

    def test_fm0_empty(self):
        """FM_0(C) is empty."""
        assert fm_dim(0) == -1

    def test_e1_contractible(self):
        """E_1(m) is contractible for m >= 1, empty for m=0."""
        assert e1_dim(0) == -1
        assert e1_dim(1) == 0
        assert e1_dim(2) == 0
        assert e1_dim(5) == 0

    def test_sc_closed_only_dim_2(self):
        """SC(2, 0) = FM_2(C), dim 2."""
        assert sc_operation_space_dim(2, 0) == 2

    def test_sc_closed_only_dim_3(self):
        """SC(3, 0) = FM_3(C), dim 4."""
        assert sc_operation_space_dim(3, 0) == 4

    def test_sc_open_only_dim(self):
        """SC(0, 2) = E_1(2), dim 0 (contractible)."""
        assert sc_operation_space_dim(0, 2) == 0

    def test_sc_mixed_1_1(self):
        """SC(1, 1): FM_1(C) x E_1(1) = point x point, dim 0."""
        assert sc_operation_space_dim(1, 1) == 0

    def test_sc_mixed_2_1(self):
        """SC(2, 1): FM_2(C) x E_1(1), dim = 2."""
        assert sc_operation_space_dim(2, 1) == 2

    def test_sc_insufficient_inputs(self):
        """SC(0, 0), SC(1, 0), SC(0, 1) are empty (need >=2 total)."""
        assert sc_operation_space_dim(0, 0) == -1
        assert sc_operation_space_dim(1, 0) == -1
        assert sc_operation_space_dim(0, 1) == -1

    def test_betti_fm2(self):
        """Total Betti of FM_2(C) = 2! = 2."""
        assert sc_betti_dim(2, 0) == 2

    def test_betti_fm3(self):
        """Total Betti of FM_3(C) = 3! = 6."""
        assert sc_betti_dim(3, 0) == 6

    def test_betti_pure_open(self):
        """Total Betti of E_1(m) = 1 (contractible)."""
        assert sc_betti_dim(0, 2) == 1
        assert sc_betti_dim(0, 5) == 1

    def test_betti_empty(self):
        """Total Betti = 0 for empty spaces."""
        assert sc_betti_dim(0, 0) == 0
        assert sc_betti_dim(1, 0) == 0


# ===================================================================
# 2. NO OPEN-TO-CLOSED DIRECTIONALITY
# ===================================================================

class TestDirectionality:
    """Test the no-open-to-closed rule."""

    def test_rule_holds(self):
        """The no-open-to-closed rule must hold."""
        result = no_open_to_closed_check()
        assert result['rule_holds'] is True

    def test_no_violations(self):
        """There should be zero violations."""
        result = no_open_to_closed_check()
        assert len(result['violations']) == 0

    def test_closed_to_closed_exists(self):
        """Closed-to-closed operations DO exist for k >= 2."""
        result = no_open_to_closed_check()
        dims = result['closed_to_closed_dims']
        # (k=2, dim=2), (k=3, dim=4), etc.
        assert len(dims) >= 2
        assert (2, 2) in dims
        assert (3, 4) in dims


# ===================================================================
# 3. HOMOTOPY-KOSZULITY
# ===================================================================

class TestHomotopyKoszul:
    """Test homotopy-Koszulity indicators."""

    def test_concentration_arity_2(self):
        """At arity 2: FM_2(C) has dim 2, Euler = 2!, top Betti = 1! = 1."""
        results = homotopy_koszul_bar_concentration(k_max=2)
        assert 2 in results
        data = results[2]
        assert data['real_dim'] == 2
        assert data['euler_char'] == 2
        assert data['top_betti'] == 1
        assert data['koszul_concentrated'] is True

    def test_concentration_arity_3(self):
        """At arity 3: FM_3(C) has dim 4, Euler = 6, top Betti = 2."""
        results = homotopy_koszul_bar_concentration(k_max=3)
        data = results[3]
        assert data['real_dim'] == 4
        assert data['euler_char'] == 6
        assert data['top_betti'] == 2
        assert data['koszul_concentrated'] is True

    def test_concentration_arity_4(self):
        """At arity 4: FM_4(C) has dim 6, Euler = 24, top Betti = 6."""
        results = homotopy_koszul_bar_concentration(k_max=4)
        data = results[4]
        assert data['real_dim'] == 6
        assert data['euler_char'] == 24
        assert data['top_betti'] == 6
        assert data['koszul_concentrated'] is True

    def test_all_concentrated_up_to_6(self):
        """Koszul concentration holds at all arities 2..6."""
        results = homotopy_koszul_bar_concentration(k_max=6)
        for k in range(2, 7):
            assert results[k]['koszul_concentrated'] is True

    def test_top_betti_is_factorial(self):
        """Arnold: top Betti of FM_k(C) = (k-1)! for all k."""
        results = homotopy_koszul_bar_concentration(k_max=8)
        for k in range(2, 9):
            assert results[k]['top_betti'] == factorial(k - 1)

    def test_euler_is_factorial(self):
        """Euler characteristic of FM_k(C) = k! for all k."""
        results = homotopy_koszul_bar_concentration(k_max=8)
        for k in range(2, 9):
            assert results[k]['euler_char'] == factorial(k)


# ===================================================================
# 4. CROSS-ENGINE CONSISTENCY
# ===================================================================

class TestCrossEngine:
    """Cross-check against sc_bar_cobar_engine if available."""

    def test_bridge_returns_data(self):
        """The cross-engine bridge should return a dict."""
        result = cross_engine_bridge()
        assert isinstance(result, dict)
        assert 'bridge_available' in result

    def test_no_mismatches_closed_only(self):
        """For closed-only operations (m=0), dimensions should agree."""
        result = cross_engine_bridge()
        if result['bridge_available']:
            assert len(result['mismatches']) == 0
