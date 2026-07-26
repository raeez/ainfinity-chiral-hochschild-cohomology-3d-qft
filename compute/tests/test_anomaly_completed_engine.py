"""Tests for Anomaly-Completed Holography Engine.

Verifies the anomaly-completed holographic constructions:
1. Transgression algebra B_Theta rank and truncations
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
    ore_transgression_relation,
    transgression_algebra,
    curved_mc_twist_data,
    virasoro_curved_class_m_completion,
    su2_anomalous_steinberg_profile,
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

    def test_noncommutative_ore_relation_has_contraction_term(self):
        """The full noncommutative Ore relation includes iota_Theta."""
        result = ore_transgression_relation(has_contraction=True)

        assert result['relation'] == 'eta b - (-1)^|b| b eta - iota_Theta(b)'
        assert result['d_eta'] == 'Theta'
        assert result['ore_automorphism'] == 'sigma(b)=(-1)^|b| b'
        assert result['ore_sigma_derivation'] == 'iota_Theta'
        assert result['central_shadow_specialization'] is False

    def test_central_shadow_is_zero_contraction_specialization(self):
        """The old central-shadow convention is the iota_Theta=0 case."""
        result = ore_transgression_relation(has_contraction=False)

        assert result['relation'] == 'eta b - (-1)^|b| b eta'
        assert result['ore_sigma_derivation'] == '0'
        assert result['central_shadow_specialization'] is True

    def test_infinite_rank_by_default(self):
        """B_Theta is free over B on eta^n, n >= 0."""
        for B_dim in [1, 3, 10, 52, 248]:
            result = transgression_algebra(B_dim, theta_degree=2)
            assert result['dim_B_Theta'] is None
            assert result['is_infinite_rank'] is True
            assert result['eta_power_basis'] == 'eta^n, n >= 0'

    def test_eta_power_truncation_dimension(self):
        """Finite dimensions require an explicit eta-power cutoff."""
        result = transgression_algebra(10, theta_degree=2, eta_power_cutoff=3)
        assert result['dim_B_Theta'] == 40
        assert result['is_infinite_rank'] is False
        assert result['dim_eta_le_1'] == 20

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

    def test_invalid_eta_power_cutoff(self):
        """Negative eta cutoff raises ValueError."""
        with pytest.raises(ValueError):
            transgression_algebra(10, theta_degree=2, eta_power_cutoff=-1)

    def test_universal_property_structure(self):
        """B_Theta stores base dimension and theta degree."""
        result = transgression_algebra(14, theta_degree=4)
        assert result['B_dim'] == 14
        assert result['theta_degree'] == 4
        assert result['convention'] == 'strict central-shadow Ore extension'
        assert result['forms_twisted_differential'] is False

    def test_curved_mc_convention(self):
        """Curved twisting has d eta + eta^2 + Theta = 0."""
        data = curved_mc_twist_data()
        assert data['twisted_differential'] == 'd_eta = d + [eta,-]'
        assert data['twisted_square'] == 'd_eta^2 = [Theta + d eta + 1/2[eta,eta], -]'
        assert data['curved_mc_equation'] == 'd eta + 1/2[eta,eta] + Theta = 0'
        assert data['twisted_curvature'] == 'Theta + d_eta_generator + eta^2'
        assert data['flatness_equation'] == 'd_eta_generator + eta^2 + Theta = 0'
        assert data['half_bracket_equals_eta_square_for_odd_eta'] is True
        assert data['strict_ore_convention'] == 'central-shadow: d(eta)=Theta, no twist d+[eta,-]'

    def test_virasoro_curved_class_m_completion(self):
        """Virasoro class-M anomaly completion uses the curved MC equation."""
        data = virasoro_curved_class_m_completion()
        assert data['base'] == 'B_ch(Vir_c)^hat_rho'
        assert data['extension'] == 'k<eta,u>/(u-eta^2)'
        assert data['curved_mc_equation'] == 'd eta + eta^2 + Theta_Vir = 0'
        assert data['twisted_differential'] == 'd_eta = d + [eta,-]'
        assert data['secondary_anomaly'] == 'u = eta^2'
        assert data['secondary_requires_algebra'] is True

    def test_item19_source_contains_curved_mc_completion(self):
        """Source guard for the curved-MC anomaly-completion formula."""
        root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        source = os.path.join(root, 'chapters', 'connections', 'anomaly_completed_core.tex')
        with open(source, encoding='utf-8') as fh:
            text = fh.read()
        assert 'def:tholog-curved-mc-anomaly-completion' in text
        assert 'd\\eta+\\frac12[\\eta,\\eta]+\\Theta=0' in text
        assert 'd\\eta+\\eta^{2}+\\Theta=0' in text
        assert 'ex:virasoro-curved-mc-anomaly-completion' in text
        assert 'd\\eta+\\eta^{2}+\\Theta_{\\mathrm{Vir}}=0' in text

    def test_item207_source_contains_ore_and_su2_clifford_dichotomy(self):
        """Source guard for the SU(2) anomaly-completed Steinberg package."""
        root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        core = os.path.join(root, 'chapters', 'connections', 'anomaly_completed_core.tex')
        frontier = os.path.join(
            root, 'chapters', 'connections', 'anomaly_completed_frontier.tex'
        )
        with open(core, encoding='utf-8') as fh:
            core_text = fh.read()
        with open(frontier, encoding='utf-8') as fh:
            frontier_text = fh.read()

        assert 'def:tholog-ore-transgression-contraction' in core_text
        assert 'eq:tholog-ore-transgression-contraction' in core_text
        assert '\\eta b-(-1)^{|b|}b\\eta-\\iota_{\\Theta}(b)' in core_text
        assert 'eq:su2-transgression' in frontier_text
        assert '\\eta b-(-1)^{|b|}b\\eta-\\iota_\\Theta(b)' in frontier_text
        assert 'eq:su2-secondary' in frontier_text
        assert 'eq:su2-anomaly-3form' in frontier_text
        assert 'H^3(SU(2);\\mathbb Z)\\cong\\mathbb Z' in frontier_text

    def test_su2_anomalous_steinberg_profile_records_pdf_dichotomy(self):
        """SU(2) profile records kappa, k c_3, u, and the two genus regimes."""
        data = su2_anomalous_steinberg_profile(S(2))

        assert data['boundary_algebra'] == 'V_k(sl_2)'
        assert data['koszul_dual_level'] == '-k-4'
        assert data['modular_characteristic'] == S(3)
        assert data['geometric_source'] == 'H^3(SU(2);Z)=Z generated by c_3'
        assert data['level_selects'] == 'k c_3'
        assert data['ore_relation'] == 'eta b - (-1)^|b| b eta - iota_Theta(b)'
        assert data['secondary_anomaly'] == 'u=eta^2'
        assert 'Mat_{2^g}' in data['invert_u']
        assert 'exterior(alpha_i,beta_i)' in data['u_zero']
        assert data['string_witness'] == 'd eta = Theta'


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
        """Finite Virasoro approximations are infinite-rank until truncated."""
        for approx_dim in [10, 50, 100]:
            result = transgression_algebra(approx_dim, theta_degree=2)
            assert result['is_infinite_rank'] is True
            eta_le_2 = transgression_algebra(
                approx_dim, theta_degree=2, eta_power_cutoff=2
            )
            assert eta_le_2['dim_B_Theta'] == 3 * approx_dim

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
