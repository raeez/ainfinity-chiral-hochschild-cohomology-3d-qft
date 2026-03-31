"""Tests for Anomaly-Completed Holography Engine.

Verifies the anomaly-completed holographic constructions:
1. Transgression algebra B_Theta dimension and structure
2. Secondary anomaly u = eta^2 degree computation
3. Neutralization obstruction in Ext^2 and moduli over Ext^1
4. Genus-Clifford completion dimension doubling
5. Cross-engine consistency with gravity_3d_engine genus data

References:
  Vol II: anomaly_completed_core.tex (Part VII)
  Vol II: 3d_gravity.tex (Movement VI)
  Vol I: concordance.tex (Theorem C)
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from fractions import Fraction
from sympy import Symbol, Rational, S, simplify

from lib.anomaly_completed_engine import (
    transgression_algebra,
    secondary_anomaly_u,
    neutralization_obstruction_degree,
    neutralization_moduli_dim,
    genus_clifford_completion,
)


# ===================================================================
# 1. TRANSGRESSION ALGEBRA
# ===================================================================

class TestTransgression:
    """Verify transgression algebra B_Theta = B * k<eta> / relations."""

    def test_dimension_doubling(self):
        """dim(B_Theta) = 2 * dim(B): extension by one generator."""
        for B_dim in [1, 3, 10, 52, 248]:
            result = transgression_algebra(B_dim, theta_degree=2)
            assert result['dim_B_Theta'] == 2 * B_dim

    def test_eta_degree(self):
        """deg(eta) = deg(Theta) - 1."""
        result = transgression_algebra(10, theta_degree=3)
        assert result['eta_degree'] == 2

    def test_eta_degree_for_mc(self):
        """For MC element Theta of degree 2: deg(eta) = 1."""
        result = transgression_algebra(10, theta_degree=2)
        assert result['eta_degree'] == 1

    def test_commutation_sign_odd_eta(self):
        """Odd-degree eta: commutation sign = -1."""
        result = transgression_algebra(10, theta_degree=2)
        # eta has degree 1 (odd), so (-1)^1 = -1
        assert result['commutation_sign'] == -1

    def test_commutation_sign_even_eta(self):
        """Even-degree eta: commutation sign = +1."""
        result = transgression_algebra(10, theta_degree=3)
        # eta has degree 2 (even), so (-1)^2 = +1
        assert result['commutation_sign'] == 1

    def test_clifford_type_odd(self):
        """Odd-degree eta gives Clifford-type algebra (eta^2 can be nonzero)."""
        result = transgression_algebra(10, theta_degree=2)
        assert result['is_clifford_type'] is True

    def test_clifford_type_even(self):
        """Even-degree eta gives exterior-type (eta^2 = 0 by graded commutativity)."""
        result = transgression_algebra(10, theta_degree=3)
        assert result['is_clifford_type'] is False

    def test_invalid_dim(self):
        """B_dim < 1 raises ValueError."""
        with pytest.raises(ValueError):
            transgression_algebra(0, theta_degree=2)

    def test_universal_property_structure(self):
        """B_Theta stores base dimension and theta degree."""
        result = transgression_algebra(14, theta_degree=4)
        assert result['B_dim'] == 14
        assert result['theta_degree'] == 4


# ===================================================================
# 2. SECONDARY ANOMALY
# ===================================================================

class TestSecondaryAnomaly:
    """Verify secondary anomaly u = eta^2."""

    def test_u_degree(self):
        """deg(u) = 2 * deg(eta)."""
        for eta_deg in range(-2, 6):
            result = secondary_anomaly_u(eta_deg)
            assert result['u_degree'] == 2 * eta_deg

    def test_u_nonzero_odd(self):
        """u is generically nonzero for odd-degree eta."""
        for eta_deg in [1, 3, 5, -1]:
            result = secondary_anomaly_u(eta_deg)
            assert result['is_nonzero'] is True

    def test_u_zero_even(self):
        """u vanishes for even-degree eta (by graded commutativity)."""
        for eta_deg in [0, 2, 4, -2]:
            result = secondary_anomaly_u(eta_deg)
            assert result['is_nonzero'] is False

    def test_u_central(self):
        """u = eta^2 is always central in B_Theta."""
        for eta_deg in range(-2, 6):
            result = secondary_anomaly_u(eta_deg)
            assert result['is_central'] is True

    def test_mc_case(self):
        """For MC element (theta_degree=2), eta_degree=1, u_degree=2."""
        result = secondary_anomaly_u(1)
        assert result['u_degree'] == 2
        assert result['is_nonzero'] is True


# ===================================================================
# 3. NEUTRALIZATION
# ===================================================================

class TestNeutralization:
    """Verify neutralization obstruction structure."""

    def test_obstruction_in_ext2(self):
        """Obstruction lives in Ext^2."""
        assert neutralization_obstruction_degree() == 2

    def test_moduli_over_ext1(self):
        """Moduli of neutralizations is affine over Ext^1."""
        result = neutralization_moduli_dim(5)
        assert result['moduli_dim'] == 5
        assert result['moduli_type'] == 'affine'

    def test_rigid_when_ext1_zero(self):
        """Ext^1 = 0 implies unique (rigid) neutralization."""
        result = neutralization_moduli_dim(0)
        assert result['is_rigid'] is True
        assert result['moduli_dim'] == 0

    def test_not_rigid_when_ext1_positive(self):
        """Ext^1 > 0 implies non-rigid neutralization."""
        for d in [1, 3, 10]:
            result = neutralization_moduli_dim(d)
            assert result['is_rigid'] is False

    def test_invalid_ext1(self):
        """Negative Ext^1 dimension raises ValueError."""
        with pytest.raises(ValueError):
            neutralization_moduli_dim(-1)


# ===================================================================
# 4. GENUS-CLIFFORD COMPLETION
# ===================================================================

class TestGenusClifford:
    """Verify genus-Clifford completion dimension doubling."""

    def test_genus_0_no_change(self):
        """Genus 0: no Clifford factors, dimension unchanged."""
        result = genus_clifford_completion(0, 10)
        assert result['completed_dim'] == 10
        assert result['clifford_factor'] == 1

    def test_genus_1_doubles(self):
        """Genus 1: one Clifford factor, dimension doubles."""
        result = genus_clifford_completion(1, 10)
        assert result['completed_dim'] == 20
        assert result['clifford_factor'] == 2

    def test_genus_2_quadruples(self):
        """Genus 2: two Clifford factors, dimension quadruples."""
        result = genus_clifford_completion(2, 10)
        assert result['completed_dim'] == 40
        assert result['clifford_factor'] == 4

    def test_genus_g_general(self):
        """Genus g: 2^g Clifford factor."""
        for g in range(7):
            result = genus_clifford_completion(g, 1)
            assert result['completed_dim'] == 2 ** g

    def test_clifford_rank(self):
        """Clifford rank = 2g."""
        for g in range(6):
            result = genus_clifford_completion(g, 5)
            assert result['clifford_rank'] == 2 * g

    def test_supertrace_sign(self):
        """Super-trace sign = (-1)^g."""
        for g in range(6):
            result = genus_clifford_completion(g, 5)
            assert result['supertrace_sign'] == (-1) ** g

    def test_composition(self):
        """(g1 + g2) Clifford = g1 Clifford * g2 Clifford (multiplicativity)."""
        B_dim = 7
        for g1 in range(4):
            for g2 in range(4):
                r_sum = genus_clifford_completion(g1 + g2, B_dim)
                r1 = genus_clifford_completion(g1, B_dim)
                r2 = genus_clifford_completion(g2, 1)
                # 2^(g1+g2) = 2^g1 * 2^g2
                assert r_sum['clifford_factor'] == r1['clifford_factor'] * r2['clifford_factor']

    def test_invalid_genus(self):
        """Negative genus raises ValueError."""
        with pytest.raises(ValueError):
            genus_clifford_completion(-1, 10)

    def test_invalid_dim(self):
        """B_dim < 1 raises ValueError."""
        with pytest.raises(ValueError):
            genus_clifford_completion(0, 0)


# ===================================================================
# 5. CROSS-ENGINE CONSISTENCY
# ===================================================================

class TestCrossEngine:
    """Consistency with gravity_3d_engine genus data."""

    def test_transgression_with_virasoro(self):
        """Transgression of Vir (infinite-dim) still doubles dimension."""
        # For a finite-dimensional approximation: truncated Virasoro
        # at level N has dim ~ N. The doubling principle still holds.
        for approx_dim in [10, 50, 100]:
            result = transgression_algebra(approx_dim, theta_degree=2)
            assert result['dim_B_Theta'] == 2 * approx_dim

    def test_genus_clifford_base_independence(self):
        """Clifford factor 2^g is independent of the base algebra."""
        for B_dim in [3, 14, 52, 248]:
            for g in [1, 2, 3]:
                result = genus_clifford_completion(g, B_dim)
                assert result['clifford_factor'] == 2 ** g

    def test_anomaly_degree_consistency(self):
        """MC element has degree 2 => eta has degree 1 => u has degree 2.

        This is consistent: the secondary anomaly u = eta^2 has the
        same degree as the primary anomaly Theta, matching the
        gravitational anomaly structure in 3d_gravity.tex.
        """
        trans = transgression_algebra(10, theta_degree=2)
        anom = secondary_anomaly_u(trans['eta_degree'])
        assert anom['u_degree'] == 2  # Same as theta_degree
