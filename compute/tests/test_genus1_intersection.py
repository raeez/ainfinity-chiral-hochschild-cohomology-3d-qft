"""Tests for genus-1 derived intersection numbers.

Verifies R^{(1)}_A(z;τ) for Heisenberg, affine sl₂, and Virasoro,
plus the genus-2 framework and modular completion theorem.
"""
import pytest
import cmath
import math

from sympy import S, Symbol, Rational, simplify

from compute.lib.genus1_intersection import (
    weierstrass_zeta_minus_rational,
    weierstrass_p_minus_rational,
    weierstrass_p2_minus_rational,
    genus1_intersection_heisenberg,
    genus1_intersection_affine_sl2,
    genus1_intersection_virasoro,
    genus1_intersection_numerical,
    genus2_intersection_framework,
    modular_completion_koszul,
)
from compute.lib.genus_one_bridge import genus1_curvature, period_correction
from compute.lib.genus2_obstruction_engine import lambda_fp


class TestWeierstraussExpansions:
    """Test Laurent expansions of Weierstrass functions."""

    def test_zeta_expansion_structure(self):
        """ζ(z)-1/z has only ODD powers of z."""
        terms = weierstrass_zeta_minus_rational(5)
        for t in terms:
            assert t['z_power'] % 2 == 1

    def test_p_expansion_structure(self):
        """℘(z)-1/z² has only EVEN powers of z."""
        terms = weierstrass_p_minus_rational(4)
        for t in terms:
            assert t['z_power'] % 2 == 0

    def test_p2_expansion_leading(self):
        """℘''(z)-6/z⁴ leading term is 6·G₄ (constant in z)."""
        terms = weierstrass_p2_minus_rational(3)
        assert terms[0]['z_power'] == 0
        assert terms[0]['coefficient'] == 6  # 2·1·(2·1-1)·(2·1+1) = 6
        assert terms[0]['G_weight'] == 4

    def test_zeta_coefficients(self):
        """ζ(z)-1/z expansion: all coefficients are -1 × G_{2n}."""
        terms = weierstrass_zeta_minus_rational(5)
        for t in terms:
            assert t['coefficient_factor'] == -1

    def test_p_coefficients(self):
        """℘(z)-1/z² expansion: coefficient of z^{2n} is (2n+1)·G_{2n+2}."""
        terms = weierstrass_p_minus_rational(4)
        for i, t in enumerate(terms):
            n = i + 1
            assert t['coefficient'] == 2 * n + 1


class TestHeisenbergGenus1:
    """Test genus-1 intersection for Heisenberg."""

    def test_structure(self):
        h = genus1_intersection_heisenberg(k_val=1)
        assert h['kappa'] == 1
        assert h['coisson_bracket'] == 0  # c₀ = 0: Heisenberg is abelian
        assert h['elliptic_regime'] == 'decoupled'  # always decoupled since c₀ = 0

    def test_level_scaling(self):
        """R¹ scales linearly with k."""
        h1 = genus1_intersection_heisenberg(k_val=1)
        h3 = genus1_intersection_heisenberg(k_val=3)
        for t1, t3 in zip(h1['intersection_number'], h3['intersection_number']):
            assert t3['coefficient'] == 3 * t1['coefficient']

    def test_leading_quasi_modular(self):
        h = genus1_intersection_heisenberg(k_val=1)
        lead = h['leading_term']
        assert lead['quasi_modular'] is True
        assert lead['z_power'] == 1

    def test_numerical_at_square_torus(self):
        """At τ=i: the leading coefficient is π."""
        r = genus1_intersection_numerical(k_val=1, z_val=0.001, tau_val=1j)
        # R¹ ≈ -G₂(i)·z = -π·z at small z
        expected = -math.pi * 0.001
        assert abs(r['R1_real'] - expected) < 1e-4

    def test_numerical_linearity_in_k(self):
        """R¹(H_k) = k · R¹(H_1) exactly."""
        r1 = genus1_intersection_numerical(k_val=1, z_val=0.1, tau_val=1j)
        r5 = genus1_intersection_numerical(k_val=5, z_val=0.1, tau_val=1j)
        assert abs(r5['R1_real'] / r1['R1_real'] - 5.0) < 1e-10

    def test_real_on_imaginary_tau(self):
        """R¹ is real when τ is purely imaginary and z is real."""
        r = genus1_intersection_numerical(k_val=1, z_val=0.1, tau_val=1j)
        assert abs(r['R1_imag']) < 1e-12

    def test_degeneration_at_cusp(self):
        """As τ→i∞ (cusp): R¹ → constant (rational function of z)."""
        r_near = genus1_intersection_numerical(k_val=1, z_val=0.1, tau_val=5j)
        r_far = genus1_intersection_numerical(k_val=1, z_val=0.1, tau_val=10j)
        # Values should converge as Im(τ) → ∞
        assert abs(r_near['R1_real'] - r_far['R1_real']) < 0.001


class TestVirasoroGenus1:
    """Test genus-1 intersection for Virasoro (three-sector structure)."""

    def test_three_sectors(self):
        v = genus1_intersection_virasoro()
        assert len(v['quartic_sector']) > 0
        assert len(v['double_sector']) > 0
        assert len(v['simple_sector']) > 0

    def test_quartic_leading_is_modular(self):
        """Leading quartic correction is (c/2)·G₄ — weight 4, modular."""
        v = genus1_intersection_virasoro()
        lead = v['quartic_sector'][0]
        assert lead['z_power'] == 0  # constant in z
        assert lead['eisenstein'] == 'G_{4}'

    def test_simple_leading_is_quasi_modular(self):
        """Leading simple correction is -∂T·G₂ — weight 2, quasi-modular."""
        v = genus1_intersection_virasoro()
        lead = v['simple_sector'][0]
        assert lead['z_power'] == 1
        assert lead['eisenstein'] == 'G_{2}'

    def test_c_zero_quartic_vanishes(self):
        """At c=0 (Witt algebra): quartic sector vanishes."""
        v = genus1_intersection_virasoro(c_val=0)
        for t in v['quartic_sector']:
            assert t['numerical_coeff'] == 0


class TestGenus2Framework:
    """Test genus-2 intersection framework."""

    def test_free_energy(self):
        g2 = genus2_intersection_framework(k_val=1)
        assert g2['F2'] == Rational(7, 5760)

    def test_level_scaling(self):
        g2_1 = genus2_intersection_framework(k_val=1)
        g2_3 = genus2_intersection_framework(k_val=3)
        assert g2_3['F2'] == 3 * g2_1['F2']

    def test_r_matrix_tower(self):
        g2 = genus2_intersection_framework(k_val=1)
        assert 0 in g2['r_matrix_genus_tower']
        assert 1 in g2['r_matrix_genus_tower']
        assert 2 in g2['r_matrix_genus_tower']


class TestModularCompletionKoszul:
    """Test modular completion theorem for Koszul algebras."""

    @pytest.mark.parametrize("family,params", [
        ('free_multiplet', {}),
        ('abelian_cs', {'k': 1}),
        ('nonabelian_cs', {'k': 3}),
    ])
    def test_automatic_extension(self, family, params):
        mc = modular_completion_koszul(family, max_genus=3, **params)
        assert mc['is_koszul'] is True
        for g in range(1, 4):
            assert mc['genus_data'][g]['extension'] == \
                'automatic (scalar curvature, nef Hodge bundle)'

    def test_ahat_generating_function(self):
        """F_g = κ·λ_g^FP matches the Â-genus at all computed genera."""
        mc = modular_completion_koszul('abelian_cs', max_genus=3, k=1)
        ahat_coeffs = [Rational(1, 24), Rational(7, 5760), Rational(31, 967680)]
        for g in range(1, 4):
            assert mc['genus_data'][g]['F_g'] == ahat_coeffs[g - 1]

    def test_cross_engine_genus1(self):
        """F₁ from modular completion matches genus_one_bridge."""
        for family, params in [('abelian_cs', {'k': 1}), ('virasoro', {'c': 26})]:
            mc = modular_completion_koszul(family, max_genus=1, **params)
            pc = period_correction(family, **params)
            assert simplify(mc['genus_data'][1]['F_g'] - pc['F1']) == 0

    def test_cross_engine_genus2(self):
        """F₂ from modular completion matches genus-2 obstruction engine."""
        mc = modular_completion_koszul('abelian_cs', max_genus=2, k=1)
        assert mc['genus_data'][2]['F_g'] == lambda_fp(2)
