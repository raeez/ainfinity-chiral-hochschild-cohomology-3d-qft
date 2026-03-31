"""Tests for Celestial Holography Engine: modular obstruction towers and one-wheel.

Verifies:
1. First modular obstruction structure (proportional to kappa)
2. Genus-g obstruction tower (recursive, degree 2)
3. Torsor structure for extensions
4. Uniqueness criterion (H^1 = H^2 = 0)
5. One-wheel Virasoro (W_1 proportional to c/2)
6. Cross-engine kappa consistency
7. Additivity of one-wheel for direct sums
8. Obstruction tower summary

Each test performs actual computation or structural verification.

References:
  Vol II: celestial_holography_core.tex
  Vol I: higher_genus_modular_koszul.tex, concordance.tex (Theorem D)
  CLAUDE.md: Critical Pitfalls AP1-AP13
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from sympy import Rational, S, Symbol, simplify

from lib.celestial_holography_engine import (
    first_modular_obstruction,
    genus_g_obstruction_class,
    torsor_structure,
    uniqueness_criterion,
    one_wheel_virasoro,
    cross_engine_kappa,
    obstruction_tower_summary,
    one_wheel_additivity,
    ObstructionData,
)


# ===================================================================
# 1. FIRST MODULAR OBSTRUCTION
# ===================================================================

class TestFirstObstruction:
    """Verify the first modular obstruction Obs_1."""

    def test_virasoro_obs1_proportional_to_kappa(self):
        """For Virasoro at c: Obs_1 is proportional to kappa = c/2."""
        c = Symbol('c')
        result = first_modular_obstruction(c / 2)
        assert simplify(result['omega_1_leading'] - c / 2) == 0

    def test_heisenberg_obs1_proportional_to_k(self):
        """For Heisenberg at level k: Obs_1 proportional to kappa = k."""
        k = Symbol('k')
        result = first_modular_obstruction(k)
        assert result['omega_1_leading'] == k

    def test_zero_kappa_means_zero_obstruction(self):
        """kappa = 0 => Obs_1 = 0 (no genus-1 obstruction)."""
        result = first_modular_obstruction(0)
        assert result['obs_1_class'] == 'zero'

    def test_nonzero_kappa_means_nonzero_obstruction(self):
        """kappa != 0 => Obs_1 != 0."""
        result = first_modular_obstruction(Rational(1, 2))
        assert result['obs_1_class'] == 'nonzero'

    def test_obs1_always_in_h2(self):
        """Obs_1 always lives in H^2(g^{(1)}, d_pi)."""
        result = first_modular_obstruction(S(5))
        assert result['h2_obstruction'] is True

    def test_virasoro_c_26_kappa_13(self):
        """Virasoro at c=26: kappa = 13, Obs_1 nonzero."""
        result = first_modular_obstruction(S(13))
        assert result['obs_1_class'] == 'nonzero'
        assert result['omega_1_leading'] == 13

    def test_virasoro_c_0_kappa_0(self):
        """Virasoro at c=0: kappa = 0, Obs_1 = 0."""
        result = first_modular_obstruction(S(0))
        assert result['obs_1_class'] == 'zero'


# ===================================================================
# 2. OBSTRUCTION TOWER
# ===================================================================

class TestObstructionTower:
    """Verify the recursive obstruction tower structure."""

    def test_genus_1_obstruction(self):
        """Genus-1 obstruction data has correct structure."""
        obs = genus_g_obstruction_class(1)
        assert obs.genus == 1
        assert obs.obstruction_degree == 2
        assert '(1)' in obs.target_cohomology

    def test_genus_2_obstruction(self):
        """Genus-2 obstruction data."""
        obs = genus_g_obstruction_class(2)
        assert obs.genus == 2
        assert obs.obstruction_degree == 2
        assert '(2)' in obs.target_cohomology

    def test_genus_3_obstruction(self):
        """Genus-3 obstruction data."""
        obs = genus_g_obstruction_class(3)
        assert obs.genus == 3
        assert '(3)' in obs.target_cohomology

    def test_all_obstructions_degree_2(self):
        """Every obstruction class lives in degree 2."""
        for g in range(1, 6):
            obs = genus_g_obstruction_class(g)
            assert obs.obstruction_degree == 2

    def test_torsor_group_is_h1(self):
        """Torsor for extensions is H^1, not H^2."""
        obs = genus_g_obstruction_class(1)
        assert 'H^1' in obs.torsor_group

    def test_genus_0_raises(self):
        """Genus 0 should raise ValueError (no obstruction at genus 0)."""
        with pytest.raises(ValueError, match="Genus must be >= 1"):
            genus_g_obstruction_class(0)

    def test_tower_summary(self):
        """Obstruction tower summary returns correct number of levels."""
        tower = obstruction_tower_summary(max_genus=5)
        assert len(tower) == 5
        assert tower[0].genus == 1
        assert tower[-1].genus == 5


# ===================================================================
# 3. TORSOR STRUCTURE
# ===================================================================

class TestTorsor:
    """Verify the torsor structure for genus-g extensions."""

    def test_h1_zero_trivial_torsor(self):
        """H^1 = 0 => trivial torsor (unique extension)."""
        result = torsor_structure(0)
        assert result['is_trivial_torsor'] is True
        assert result['torsor_exists'] is True

    def test_h1_nonzero_nontrivial_torsor(self):
        """H^1 != 0 => nontrivial torsor (non-unique extension)."""
        result = torsor_structure(1)
        assert result['is_trivial_torsor'] is False
        assert result['torsor_exists'] is True

    def test_h1_large_still_exists(self):
        """Even large H^1 still gives a torsor (if extension exists)."""
        result = torsor_structure(10)
        assert result['torsor_exists'] is True
        assert result['h1_dim'] == 10


# ===================================================================
# 4. UNIQUENESS CRITERION
# ===================================================================

class TestUniqueness:
    """Verify the full uniqueness criterion for modular lifts."""

    def test_both_zero_unique(self):
        """H^1 = H^2 = 0 => unique lift."""
        result = uniqueness_criterion(0, 0)
        assert result['lift_exists'] is True
        assert result['lift_unique'] is True
        assert result['full_uniqueness'] is True

    def test_h2_nonzero_no_lift(self):
        """H^2 != 0 => obstruction, no guaranteed lift."""
        result = uniqueness_criterion(0, 1)
        assert result['lift_exists'] is False
        assert result['lift_unique'] is False

    def test_h1_nonzero_h2_zero_nonunique(self):
        """H^1 != 0 but H^2 = 0 => lift exists but is not unique."""
        result = uniqueness_criterion(1, 0)
        assert result['lift_exists'] is True
        assert result['lift_unique'] is False

    def test_both_nonzero(self):
        """H^1 != 0 and H^2 != 0 => neither exists nor unique."""
        result = uniqueness_criterion(1, 1)
        assert result['lift_exists'] is False
        assert result['lift_unique'] is False


# ===================================================================
# 5. ONE-WHEEL VIRASORO
# ===================================================================

class TestOneWheel:
    """Verify the one-wheel sum for the Virasoro algebra."""

    def test_w1_proportional_to_kappa(self):
        """W_1(Vir_c) is proportional to kappa = c/2."""
        c = Symbol('c')
        result = one_wheel_virasoro(c)
        assert result['w1_proportional_to_kappa'] is True
        assert simplify(result['w1_leading'] - c / 2) == 0

    def test_w1_at_c_1(self):
        """W_1(Vir_1) = 1/2."""
        result = one_wheel_virasoro(S(1))
        assert result['w1_leading'] == Rational(1, 2)

    def test_w1_at_c_26(self):
        """W_1(Vir_26) = 13."""
        result = one_wheel_virasoro(S(26))
        assert result['w1_leading'] == 13

    def test_w1_at_c_0(self):
        """W_1(Vir_0) = 0."""
        result = one_wheel_virasoro(S(0))
        assert result['w1_leading'] == 0

    def test_proportionality_constant_is_1(self):
        """The proportionality constant W_1 / kappa = 1."""
        result = one_wheel_virasoro(Symbol('c'))
        assert result['proportionality_constant'] == 1

    def test_mode_sum_finite(self):
        """Truncated mode sum is finite and positive for c > 0."""
        result = one_wheel_virasoro(S(26), max_modes=5)
        # mode_sum = 5^2 * 26/12 = 25 * 26/12 = 650/12 = 325/6
        assert result['mode_sum_truncated'] == Rational(325, 6)

    def test_kappa_matches_gravity_engine(self):
        """kappa from one-wheel matches gravity_3d_engine convention.

        gravity_3d_engine: kappa(Vir_c) = c/2.
        one_wheel_virasoro: W_1 leading = c/2.
        """
        c = Symbol('c')
        result = one_wheel_virasoro(c)
        assert simplify(result['kappa'] - c / 2) == 0


# ===================================================================
# 6. CROSS-ENGINE KAPPA CONSISTENCY
# ===================================================================

class TestCrossEngine:
    """Verify kappa consistency across engines."""

    def test_virasoro_kappa_c_over_2(self):
        """Virasoro kappa = c/2 matches other engines."""
        c = Symbol('c')
        result = cross_engine_kappa('virasoro', c=c)
        assert simplify(result['kappa'] - c / 2) == 0

    def test_heisenberg_kappa_k(self):
        """Heisenberg kappa = k matches other engines."""
        k = Symbol('k')
        result = cross_engine_kappa('heisenberg', k=k)
        assert simplify(result['kappa'] - k) == 0

    def test_affine_sl2_kappa(self):
        """Affine sl_2: kappa = 3(k+2)/4."""
        k = Symbol('k')
        result = cross_engine_kappa('affine_sl2', k=k)
        expected = Rational(3, 1) * (k + 2) / 4
        assert simplify(result['kappa'] - expected) == 0

    def test_unknown_type_raises(self):
        """Unknown algebra type should raise."""
        with pytest.raises(ValueError, match="Unknown"):
            cross_engine_kappa('unknown')


# ===================================================================
# 7. ADDITIVITY
# ===================================================================

class TestAdditivity:
    """Verify additivity of one-wheel for direct sums."""

    def test_heisenberg_plus_virasoro(self):
        """W(H_k + Vir_c) = W(H_k) + W(Vir_c) = k + c/2."""
        k, c = Symbol('k'), Symbol('c')
        result = one_wheel_additivity(k, c / 2)
        assert result['additive'] is True
        assert simplify(result['kappa_sum'] - k - c / 2) == 0

    def test_two_heisenberg(self):
        """W(H_k1 + H_k2) = k1 + k2."""
        result = one_wheel_additivity(S(1), S(2))
        assert result['kappa_sum'] == 3

    def test_zero_plus_anything(self):
        """W(0 + A) = W(A)."""
        result = one_wheel_additivity(S(0), S(5))
        assert result['kappa_sum'] == 5
