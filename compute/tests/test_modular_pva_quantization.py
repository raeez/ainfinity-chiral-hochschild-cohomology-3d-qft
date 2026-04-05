"""Tests for modular PVA quantization: genus expansion from classical to quantum.

Verifies the genus expansion programme for quantizing PVAs to vertex algebras.
Covers genus-1 obstruction for all families, critical level analysis,
DS identity, quantization tower, and cross-volume consistency.

References:
  modular_pva_quantization.tex (Vol II)
  higher_genus_modular_koszul.tex (Vol I): modular bar, genus spectral sequence
  genus_one_bridge.py (Vol II): kappa, period correction, Arnold defect

Tier 1 (structural): all tests are algebraic identities or consistency checks.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy import Symbol, Rational, simplify, expand, S, symbols, oo


# ===================================================================
# GENUS-0 CLASSICAL DATA
# ===================================================================

class TestGenus0ClassicalData:
    """Verify classical PVA data for all standard families."""

    def test_heisenberg_kappa(self):
        """Heisenberg: kappa = k."""
        from lib.modular_pva_quantization import genus0_classical_data
        k = Symbol('k')
        data = genus0_classical_data('heisenberg', k=k)
        assert simplify(data['kappa'] - k) == 0

    def test_heisenberg_shadow_depth(self):
        """Heisenberg: Gaussian archetype, depth 2."""
        from lib.modular_pva_quantization import genus0_classical_data
        data = genus0_classical_data('heisenberg', k=1)
        assert data['shadow_depth'] == 2
        assert data['shadow_archetype'] == 'G'

    def test_affine_sl2_kappa(self):
        """Affine sl_2: kappa = dim(sl_2)*(k+h^v)/(2*h^v) = 3*(k+2)/4."""
        from lib.modular_pva_quantization import genus0_classical_data
        k = Symbol('k')
        data = genus0_classical_data('affine_sl2', k=k)
        # dim(sl_2) = 3, h^v = 2
        expected_kappa = 3 * (k + 2) / 4
        assert simplify(data['kappa'] - expected_kappa) == 0

    def test_affine_sl3_kappa(self):
        """Affine sl_3: kappa = dim(sl_3)*(k+h^v)/(2*h^v) = 4*(k+3)/3."""
        from lib.modular_pva_quantization import genus0_classical_data
        k = Symbol('k')
        data = genus0_classical_data('affine_sl3', k=k)
        # dim(sl_3) = 8, h^v = 3
        expected_kappa = 8 * (k + 3) / 6
        assert simplify(data['kappa'] - expected_kappa) == 0

    def test_affine_shadow_depth(self):
        """Affine: Lie/tree archetype, depth 3."""
        from lib.modular_pva_quantization import genus0_classical_data
        data = genus0_classical_data('affine_sl2', k=1)
        assert data['shadow_depth'] == 3
        assert data['shadow_archetype'] == 'L'

    def test_virasoro_kappa(self):
        """Virasoro: kappa = c/2."""
        from lib.modular_pva_quantization import genus0_classical_data
        c = Symbol('c')
        data = genus0_classical_data('virasoro', c=c)
        assert simplify(data['kappa'] - c / 2) == 0

    def test_virasoro_shadow_depth_infinite(self):
        """Virasoro: Mixed archetype, infinite depth."""
        from lib.modular_pva_quantization import genus0_classical_data
        data = genus0_classical_data('virasoro', c=1)
        assert data['shadow_depth'] == oo
        assert data['shadow_archetype'] == 'M'

    def test_w3_kappa(self):
        """W_3: kappa = 5c/6."""
        from lib.modular_pva_quantization import genus0_classical_data
        c = Symbol('c')
        data = genus0_classical_data('w3', c=c)
        assert simplify(data['kappa'] - 5 * c / 6) == 0

    def test_w3_resonance_divisor(self):
        """W_3: resonance at c = -22/5."""
        from lib.modular_pva_quantization import genus0_classical_data
        data = genus0_classical_data('w3', c=Symbol('c'))
        assert data['resonance_divisor'] == Rational(-22, 5)

    def test_betagamma_kappa(self):
        """betagamma: kappa = +1 (c = +2)."""
        from lib.modular_pva_quantization import genus0_classical_data
        data = genus0_classical_data('betagamma')
        assert data['kappa'] == 1
        assert data['shadow_depth'] == 4
        assert data['shadow_archetype'] == 'C'

    def test_free_multiplet_kappa(self):
        """Free multiplet: kappa = 1/2."""
        from lib.modular_pva_quantization import genus0_classical_data
        data = genus0_classical_data('free_multiplet')
        assert data['kappa'] == Rational(1, 2)

    def test_genus0_bar_is_flat(self):
        """All families: d_0^2 = 0 at genus 0."""
        from lib.modular_pva_quantization import genus0_classical_data
        for family in ['heisenberg', 'affine_sl2', 'virasoro', 'w3', 'betagamma']:
            data = genus0_classical_data(family, k=1, c=1)
            assert data['bar_differential_squared'] == 0


# ===================================================================
# GENUS-1 OBSTRUCTION
# ===================================================================

class TestGenus1Obstruction:
    """Verify Ob_1 = 0 for all standard families at generic parameters."""

    def test_heisenberg_ob1_vanishes(self):
        """Heisenberg: Ob_1 = 0 (no interaction)."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('heisenberg', k=Symbol('k'))
        assert ob['vanishes'] is True
        assert ob['obstruction'] == 0

    def test_affine_sl2_ob1_vanishes_generic(self):
        """Affine sl_2: Ob_1 = 0 for generic k."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('affine_sl2', k=Symbol('k'))
        assert ob['vanishes'] is True

    def test_affine_sl2_ob1_nonzero_critical(self):
        """Affine sl_2: Ob_1 != 0 at critical level k = -2."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('affine_sl2', k=-2)
        assert ob['vanishes'] is False

    def test_affine_sl3_ob1_vanishes_generic(self):
        """Affine sl_3: Ob_1 = 0 for generic k."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('affine_sl3', k=Symbol('k'))
        assert ob['vanishes'] is True

    def test_affine_sl3_ob1_nonzero_critical(self):
        """Affine sl_3: Ob_1 != 0 at critical level k = -3."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('affine_sl3', k=-3)
        assert ob['vanishes'] is False

    def test_virasoro_ob1_vanishes_generic(self):
        """Virasoro: Ob_1 = 0 for generic c."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('virasoro', c=Symbol('c'))
        assert ob['vanishes'] is True

    def test_virasoro_ob1_nonzero_c0(self):
        """Virasoro: Ob_1 != 0 at c = 0 (degenerate)."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('virasoro', c=0)
        assert ob['vanishes'] is False

    def test_virasoro_ob1_at_c1(self):
        """Virasoro: Ob_1 = 0 at c = 1 (free boson)."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('virasoro', c=1)
        assert ob['vanishes'] is True

    def test_virasoro_ob1_at_c26(self):
        """Virasoro: Ob_1 = 0 at c = 26 (bosonic string)."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('virasoro', c=26)
        assert ob['vanishes'] is True

    def test_w3_ob1_vanishes_generic(self):
        """W_3: Ob_1 = 0 for generic c."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('w3', c=Symbol('c'))
        assert ob['vanishes'] is True

    def test_w3_ob1_nonzero_resonance(self):
        """W_3: Ob_1 != 0 at resonance c = -22/5."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('w3', c=Rational(-22, 5))
        assert ob['vanishes'] is False

    def test_w3_ob1_nonzero_c0(self):
        """W_3: Ob_1 != 0 at c = 0."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('w3', c=0)
        assert ob['vanishes'] is False

    def test_betagamma_ob1_vanishes(self):
        """betagamma: Ob_1 = 0 (quadratic OPE)."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('betagamma')
        assert ob['vanishes'] is True

    def test_free_multiplet_ob1_vanishes(self):
        """Free multiplet: Ob_1 = 0 (trivially)."""
        from lib.modular_pva_quantization import genus1_obstruction
        ob = genus1_obstruction('free_multiplet')
        assert ob['vanishes'] is True


# ===================================================================
# GENUS-1 LOOP EQUATION
# ===================================================================

class TestGenus1LoopEquation:
    """Verify the genus-1 correction Theta_1 = kappa/24."""

    def test_heisenberg_theta1(self):
        """Heisenberg: Theta_1 = k/24 (κ=k)."""
        from lib.modular_pva_quantization import genus1_loop_equation
        k = Symbol('k')
        eq = genus1_loop_equation('heisenberg', k=k)
        assert simplify(eq['theta_1'] - k / 24) == 0

    def test_virasoro_theta1(self):
        """Virasoro: Theta_1 = c/48."""
        from lib.modular_pva_quantization import genus1_loop_equation
        c = Symbol('c')
        eq = genus1_loop_equation('virasoro', c=c)
        assert simplify(eq['theta_1'] - c / 48) == 0

    def test_w3_theta1(self):
        """W_3: Theta_1 = 5c/144."""
        from lib.modular_pva_quantization import genus1_loop_equation
        c = Symbol('c')
        eq = genus1_loop_equation('w3', c=c)
        assert simplify(eq['theta_1'] - 5 * c / 144) == 0

    def test_affine_sl2_theta1(self):
        """Affine sl_2: Theta_1 = kappa/24 = 3*(k+2)/4 / 24 = (k+2)/32."""
        from lib.modular_pva_quantization import genus1_loop_equation
        k = Symbol('k')
        eq = genus1_loop_equation('affine_sl2', k=k)
        expected = 3 * (k + 2) / 4 * Rational(1, 24)
        assert simplify(eq['theta_1'] - expected) == 0

    def test_faber_pandharipande_lambda1(self):
        """lambda_1^FP = 1/24."""
        from lib.modular_pva_quantization import genus1_loop_equation
        eq = genus1_loop_equation('heisenberg', k=1)
        assert eq['faber_pandharipande_lambda1'] == Rational(1, 24)

    def test_D1_squared_restored(self):
        """Period correction restores D_1^2 = 0 for all families."""
        from lib.modular_pva_quantization import genus1_loop_equation
        for family in ['heisenberg', 'virasoro', 'w3', 'betagamma']:
            eq = genus1_loop_equation(family, k=1, c=1)
            assert eq['D1_squared_restored'] is True


# ===================================================================
# CRITICAL LEVEL OBSTRUCTION
# ===================================================================

class TestCriticalLevelObstruction:
    """Verify Ob_1 != 0 at the critical level for affine algebras."""

    def test_sl2_critical_level(self):
        """sl_2: critical level k = -h^v = -2."""
        from lib.modular_pva_quantization import critical_level_obstruction
        data = critical_level_obstruction('sl2')
        assert data['critical_level'] == -2
        assert data['h_dual'] == 2
        assert data['sugawara_defined'] is False
        assert data['obstruction_nonzero'] is True

    def test_sl3_critical_level(self):
        """sl_3: critical level k = -h^v = -3."""
        from lib.modular_pva_quantization import critical_level_obstruction
        data = critical_level_obstruction('sl3')
        assert data['critical_level'] == -3
        assert data['h_dual'] == 3
        assert data['obstruction_nonzero'] is True

    def test_sl4_critical_level(self):
        """sl_4: critical level k = -h^v = -4."""
        from lib.modular_pva_quantization import critical_level_obstruction
        data = critical_level_obstruction('sl4')
        assert data['critical_level'] == -4
        assert data['h_dual'] == 4

    def test_feigin_frenkel_center(self):
        """At critical level: Feigin-Frenkel center exists."""
        from lib.modular_pva_quantization import critical_level_obstruction
        for lie_type in ['sl2', 'sl3']:
            data = critical_level_obstruction(lie_type)
            assert data['feigin_frenkel_center'] is True

    def test_pbw_collapses_at_critical(self):
        """PBW filtration collapses at critical level."""
        from lib.modular_pva_quantization import critical_level_obstruction
        data = critical_level_obstruction('sl2')
        assert data['pbw_intact'] is False


# ===================================================================
# MODULAR BAR CURVATURE
# ===================================================================

class TestModularBarCurvature:
    """Verify d^2_g = kappa * omega_g at each genus."""

    def test_genus0_flat(self):
        """All families: d_0^2 = 0."""
        from lib.modular_pva_quantization import modular_bar_curvature
        for family in ['heisenberg', 'virasoro', 'w3']:
            curv = modular_bar_curvature(family, 0, k=1, c=1)
            assert curv['is_flat'] is True
            assert curv['d_squared'] == 0

    def test_genus1_heisenberg_curvature(self):
        """Heisenberg at genus 1: d^2 = k * omega_1."""
        from lib.modular_pva_quantization import modular_bar_curvature
        k = Symbol('k')
        curv = modular_bar_curvature('heisenberg', 1, k=k)
        assert simplify(curv['d_squared'] - k) == 0

    def test_genus1_virasoro_curvature(self):
        """Virasoro at genus 1: d^2 = c/2 * omega_1."""
        from lib.modular_pva_quantization import modular_bar_curvature
        c = Symbol('c')
        curv = modular_bar_curvature('virasoro', 1, c=c)
        assert simplify(curv['d_squared'] - c / 2) == 0

    def test_genus1_w3_curvature(self):
        """W_3 at genus 1: d^2 = 5c/6 * omega_1."""
        from lib.modular_pva_quantization import modular_bar_curvature
        c = Symbol('c')
        curv = modular_bar_curvature('w3', 1, c=c)
        assert simplify(curv['d_squared'] - 5 * c / 6) == 0

    def test_period_correction_genus1(self):
        """Period correction at genus 1: kappa * lambda_1^FP = kappa/24."""
        from lib.modular_pva_quantization import modular_bar_curvature
        c = Symbol('c')
        curv = modular_bar_curvature('virasoro', 1, c=c)
        assert curv['faber_pandharipande'] == Rational(1, 24)
        assert simplify(curv['period_correction'] - c / 48) == 0

    def test_genus2_faber_pandharipande(self):
        """Genus 2: lambda_2^FP = 7/5760."""
        from lib.modular_pva_quantization import modular_bar_curvature
        curv = modular_bar_curvature('virasoro', 2, c=Symbol('c'))
        assert curv['faber_pandharipande'] == Rational(7, 5760)


# ===================================================================
# DS QUANTIZATION IDENTITY
# ===================================================================

class TestDSQuantizationIdentity:
    """Verify DS(V_k^{quantum}) = W_k^{quantum}."""

    def test_sl2_principal_ds(self):
        """sl_2 principal DS: Virasoro from affine sl_2."""
        from lib.modular_pva_quantization import ds_quantization_identity
        k = Symbol('k')
        result = ds_quantization_identity('sl2', k)
        assert result['identity_verified'] is True
        assert result['status'] == 'PROVED'
        # c_ds should be c_Vir = 1 - 6(k+1)^2/(k+2)
        expected_c = 1 - 6 * (k + 1)**2 / (k + 2)
        assert simplify(result['c_ds'] - expected_c) == 0

    def test_sl3_principal_ds(self):
        """sl_3 principal DS: W_3 from affine sl_3."""
        from lib.modular_pva_quantization import ds_quantization_identity
        k = Symbol('k')
        result = ds_quantization_identity('sl3', k)
        assert result['identity_verified'] is True
        # c_ds should be 2 - 24(k+2)^2/(k+3)
        expected_c = 2 - 24 * (k + 2)**2 / (k + 3)
        assert simplify(result['c_ds'] - expected_c) == 0

    def test_ds_complementarity_sl2(self):
        """sl_2 DS: c(k) + c(k') = 26 (Virasoro complementarity)."""
        from lib.modular_pva_quantization import ds_quantization_identity
        k = Symbol('k')
        result = ds_quantization_identity('sl2', k)
        assert simplify(result['complementarity_sum'] - 26) == 0

    def test_ds_complementarity_sl3(self):
        """sl_3 DS: c(k) + c(k') = 100 (W_3 complementarity)."""
        from lib.modular_pva_quantization import ds_quantization_identity
        k = Symbol('k')
        result = ds_quantization_identity('sl3', k)
        assert simplify(result['complementarity_sum'] - 100) == 0

    def test_non_principal_ds_conjectural(self):
        """Non-principal DS: status CONJECTURAL."""
        from lib.modular_pva_quantization import ds_quantization_identity
        result = ds_quantization_identity('sl3', Symbol('k'), f_nilpotent='subregular')
        assert result['status'] == 'CONJECTURAL'

    def test_ds_affine_obstruction_consistency(self):
        """DS identity: affine Ob_1 status propagates to W-algebra."""
        from lib.modular_pva_quantization import ds_quantization_identity
        k = Symbol('k')
        result = ds_quantization_identity('sl2', k)
        assert result['affine_obstruction_vanishes'] is True


# ===================================================================
# GENUS EXPANSION COEFFICIENTS
# ===================================================================

class TestGenusExpansionCoefficients:
    """Verify Theta_g at each genus."""

    def test_genus0_coefficient(self):
        """Genus 0: classical data."""
        from lib.modular_pva_quantization import genus_expansion_coefficients
        c = Symbol('c')
        coeffs = genus_expansion_coefficients('virasoro', max_g=3, c=c)
        assert coeffs[0]['status'] == 'PROVED'

    def test_genus1_coefficient_virasoro(self):
        """Virasoro genus 1: Theta_1 = c/48."""
        from lib.modular_pva_quantization import genus_expansion_coefficients
        c = Symbol('c')
        coeffs = genus_expansion_coefficients('virasoro', max_g=1, c=c)
        assert simplify(coeffs[1]['theta_g'] - c / 48) == 0

    def test_genus2_coefficient_virasoro(self):
        """Virasoro genus 2: Theta_2 = 7c/11520."""
        from lib.modular_pva_quantization import genus_expansion_coefficients
        c = Symbol('c')
        coeffs = genus_expansion_coefficients('virasoro', max_g=2, c=c)
        expected = c / 2 * Rational(7, 5760)
        assert simplify(coeffs[2]['theta_g'] - expected) == 0

    def test_genus1_coefficient_w3(self):
        """W_3 genus 1: Theta_1 = 5c/144."""
        from lib.modular_pva_quantization import genus_expansion_coefficients
        c = Symbol('c')
        coeffs = genus_expansion_coefficients('w3', max_g=1, c=c)
        assert simplify(coeffs[1]['theta_g'] - 5 * c / 144) == 0


# ===================================================================
# QUANTIZATION OBSTRUCTION TOWER
# ===================================================================

class TestQuantizationObstructionTower:
    """Verify the full obstruction tower Ob_g at each genus."""

    def test_genus0_no_obstruction(self):
        """Genus 0: Ob_0 = 0 (classical axioms)."""
        from lib.modular_pva_quantization import quantization_obstruction_tower
        tower = quantization_obstruction_tower('virasoro', c=Symbol('c'))
        assert tower[0]['vanishes'] is True

    def test_tower_unobstructed_virasoro(self):
        """Virasoro: full tower unobstructed at generic c."""
        from lib.modular_pva_quantization import quantization_obstruction_tower
        tower = quantization_obstruction_tower('virasoro', max_g=3, c=Symbol('c'))
        for g in range(4):
            assert tower[g]['vanishes'] is True

    def test_tower_unobstructed_heisenberg(self):
        """Heisenberg: full tower unobstructed."""
        from lib.modular_pva_quantization import quantization_obstruction_tower
        tower = quantization_obstruction_tower('heisenberg', max_g=3, k=Symbol('k'))
        for g in range(4):
            assert tower[g]['vanishes'] is True

    def test_tower_obstructed_at_critical(self):
        """Affine sl_2 at k=-2: tower obstructed."""
        from lib.modular_pva_quantization import quantization_obstruction_tower
        tower = quantization_obstruction_tower('affine_sl2', max_g=3, k=-2)
        assert tower[1]['vanishes'] is False


# ===================================================================
# VIRASORO QUANTIZATION DATA
# ===================================================================

class TestVirasoroQuantizationData:
    """Full quantization data for Virasoro."""

    def test_virasoro_self_dual_at_13(self):
        """Virasoro self-dual at c = 13 (NOT c = 26)."""
        from lib.modular_pva_quantization import virasoro_quantization_data
        data = virasoro_quantization_data(13)
        assert data['self_dual_point'] == 13
        assert simplify(data['kappa'] - data['kappa_dual']) == 0

    def test_virasoro_complementarity(self):
        """Virasoro: kappa + kappa_dual = 13."""
        from lib.modular_pva_quantization import virasoro_quantization_data
        c = Symbol('c')
        data = virasoro_quantization_data(c)
        assert simplify(data['complementarity_sum'] - 13) == 0

    def test_virasoro_quartic_contact(self):
        """Q^contact_Vir = 10/[c(5c+22)]."""
        from lib.modular_pva_quantization import virasoro_quantization_data
        c = Symbol('c')
        data = virasoro_quantization_data(c)
        expected = 10 / (c * (5 * c + 22))
        assert simplify(data['quartic_contact'] - expected) == 0

    def test_virasoro_quintic_forced(self):
        """Virasoro: shadow obstruction tower is infinite (quintic forced)."""
        from lib.modular_pva_quantization import virasoro_quantization_data
        data = virasoro_quantization_data(1)
        assert data['quintic_forced'] is True
        assert data['shadow_depth'] == oo

    def test_virasoro_quantizes_at_c1(self):
        """Virasoro quantizes at c = 1."""
        from lib.modular_pva_quantization import virasoro_quantization_data
        data = virasoro_quantization_data(1)
        assert data['obstruction_vanishes'] is True


# ===================================================================
# W_3 QUANTIZATION DATA
# ===================================================================

class TestW3QuantizationData:
    """Full quantization data for W_3."""

    def test_w3_complementarity(self):
        """W_3: kappa + kappa_dual = 5*100/6 = 500/6."""
        from lib.modular_pva_quantization import w3_quantization_data
        c = Symbol('c')
        data = w3_quantization_data(c)
        assert simplify(data['complementarity_sum'] - Rational(500, 6)) == 0

    def test_w3_beta_squared(self):
        """W_3: beta^2 = 16/(22 + 5c)."""
        from lib.modular_pva_quantization import w3_quantization_data
        c = Symbol('c')
        data = w3_quantization_data(c)
        expected = 16 / (22 + 5 * c)
        assert simplify(data['beta_squared'] - expected) == 0

    def test_w3_resonant_at_minus22over5(self):
        """W_3: resonant at c = -22/5."""
        from lib.modular_pva_quantization import w3_quantization_data
        data = w3_quantization_data(Rational(-22, 5))
        assert data['is_resonant'] is True
        assert data['obstruction_vanishes'] is False

    def test_w3_not_resonant_generic(self):
        """W_3: not resonant at generic c."""
        from lib.modular_pva_quantization import w3_quantization_data
        data = w3_quantization_data(Symbol('c'))
        assert data['is_resonant'] is False

    def test_w3_quantizes_at_c100(self):
        """W_3 quantizes at c = 100."""
        from lib.modular_pva_quantization import w3_quantization_data
        data = w3_quantization_data(100)
        assert data['obstruction_vanishes'] is True


# ===================================================================
# VERIFY ALL FAMILIES QUANTIZE
# ===================================================================

class TestAllFamiliesQuantize:
    """Verify the landscape-wide quantization result."""

    def test_all_families_ob1_vanishes(self):
        """All standard families have Ob_1 = 0 at generic parameters."""
        from lib.modular_pva_quantization import verify_all_families_quantize
        results = verify_all_families_quantize()
        assert results['_all_quantize_generically'] is True

    def test_each_family_quantizes(self):
        """Each individual family quantizes."""
        from lib.modular_pva_quantization import verify_all_families_quantize
        results = verify_all_families_quantize()
        for family in ['heisenberg', 'affine_sl2', 'affine_sl3',
                        'virasoro', 'w3', 'betagamma', 'free_multiplet']:
            assert results[family]['vanishes'] is True, \
                f"{family} does not quantize generically"

    def test_obstruction_loci_correct(self):
        """Obstruction loci: affine at critical level, Vir/W at c=0."""
        from lib.modular_pva_quantization import verify_all_families_quantize
        results = verify_all_families_quantize()
        assert -2 in results['affine_sl2']['obstruction_locus']
        assert -3 in results['affine_sl3']['obstruction_locus']
        assert 0 in results['virasoro']['obstruction_locus']


# ===================================================================
# QUANTUM CORRECTION FORMULA
# ===================================================================

class TestQuantumCorrectionFormula:
    """Verify explicit Theta_g formulas."""

    def test_genus1_formula_heisenberg(self):
        """Heisenberg: Theta_1 = k/24 (κ=k)."""
        from lib.modular_pva_quantization import quantum_correction_formula
        k = Symbol('k')
        result = quantum_correction_formula('heisenberg', 1, k=k)
        assert simplify(result['theta_g'] - k / 24) == 0
        assert result['lambda_g_FP'] == Rational(1, 24)

    def test_genus1_formula_virasoro(self):
        """Virasoro: Theta_1 = c/48."""
        from lib.modular_pva_quantization import quantum_correction_formula
        c = Symbol('c')
        result = quantum_correction_formula('virasoro', 1, c=c)
        assert simplify(result['theta_g'] - c / 48) == 0

    def test_genus2_formula_virasoro(self):
        """Virasoro: Theta_2 = 7c/11520."""
        from lib.modular_pva_quantization import quantum_correction_formula
        c = Symbol('c')
        result = quantum_correction_formula('virasoro', 2, c=c)
        expected = c / 2 * Rational(7, 5760)
        assert simplify(result['theta_g'] - expected) == 0

    def test_genus3_formula_w3(self):
        """W_3: Theta_3 = (5c/6) * (1/82944)."""
        from lib.modular_pva_quantization import quantum_correction_formula
        c = Symbol('c')
        result = quantum_correction_formula('w3', 3, c=c)
        expected = 5 * c / 6 * Rational(1, 82944)
        assert simplify(result['theta_g'] - expected) == 0

    def test_genus1_proved_status(self):
        """Genus 1: status PROVED."""
        from lib.modular_pva_quantization import quantum_correction_formula
        result = quantum_correction_formula('virasoro', 1, c=1)
        assert result['status'] == 'PROVED'

    def test_genus2_conjectural_status(self):
        """Genus >= 2: status CONJECTURAL (scalar part)."""
        from lib.modular_pva_quantization import quantum_correction_formula
        result = quantum_correction_formula('virasoro', 2, c=1)
        assert 'CONJECTURAL' in result['status']


# ===================================================================
# COMPLEMENTARITY (CROSS-VOLUME CONSISTENCY)
# ===================================================================

class TestComplementarity:
    """Cross-volume consistency of kappa complementarity."""

    def test_heisenberg_complementarity(self):
        """Heisenberg: kappa + kappa_dual = 0."""
        from lib.modular_pva_quantization import complementarity_sum
        k = Symbol('k')
        assert simplify(complementarity_sum('heisenberg', k=k)) == 0

    def test_virasoro_complementarity(self):
        """Virasoro: kappa + kappa_dual = 13."""
        from lib.modular_pva_quantization import complementarity_sum
        c = Symbol('c')
        assert simplify(complementarity_sum('virasoro', c=c) - 13) == 0

    def test_w3_complementarity(self):
        """W_3: kappa + kappa_dual = 500/6."""
        from lib.modular_pva_quantization import complementarity_sum
        c = Symbol('c')
        assert simplify(complementarity_sum('w3', c=c) - Rational(500, 6)) == 0

    def test_affine_sl2_complementarity(self):
        """Affine sl_2: kappa + kappa_dual = 0 (KM anti-symmetry)."""
        from lib.modular_pva_quantization import complementarity_sum
        k = Symbol('k')
        assert simplify(complementarity_sum('affine_sl2', k=k)) == 0


# ===================================================================
# GENUS SPECTRAL SEQUENCE
# ===================================================================

class TestGenusSpectralSequence:
    """Verify genus spectral sequence data."""

    def test_e1_p0_classical(self):
        """E_1 page p=0: classical data."""
        from lib.modular_pva_quantization import genus_spectral_sequence_data
        data = genus_spectral_sequence_data('virasoro', c=Symbol('c'))
        assert data['E1_p0']['d_0_squared'] == 0

    def test_e1_p1_genus1_correction(self):
        """E_1 page p=1: genus-1 correction."""
        from lib.modular_pva_quantization import genus_spectral_sequence_data
        c = Symbol('c')
        data = genus_spectral_sequence_data('virasoro', c=c)
        assert simplify(data['E1_p1']['theta_1'] - c / 48) == 0

    def test_spectral_sequence_distinct_from_pbw(self):
        """Genus spectral sequence is DISTINCT from PBW spectral sequence."""
        from lib.modular_pva_quantization import genus_spectral_sequence_data
        data = genus_spectral_sequence_data('virasoro', c=1)
        assert 'DISTINCT' in data['note']


# ===================================================================
# FULL QUANTIZATION SUMMARY
# ===================================================================

class TestFullQuantizationSummary:
    """Verify the landscape-wide summary table."""

    def test_summary_covers_all_families(self):
        """Summary covers all 7 standard families."""
        from lib.modular_pva_quantization import full_quantization_summary
        summary = full_quantization_summary()
        expected_families = {
            'heisenberg', 'affine_sl2', 'affine_sl3',
            'virasoro', 'w3', 'betagamma', 'free_multiplet'
        }
        assert set(summary.keys()) == expected_families

    def test_summary_all_quantizable(self):
        """All families are quantizable at generic parameters."""
        from lib.modular_pva_quantization import full_quantization_summary
        summary = full_quantization_summary()
        for family, data in summary.items():
            assert data['quantizable'] is True, f"{family} not quantizable"

    def test_summary_shadow_archetypes(self):
        """Shadow archetypes are correctly assigned."""
        from lib.modular_pva_quantization import full_quantization_summary
        summary = full_quantization_summary()
        assert summary['heisenberg']['shadow_archetype'] == 'G'
        assert summary['affine_sl2']['shadow_archetype'] == 'L'
        assert summary['virasoro']['shadow_archetype'] == 'M'
        assert summary['w3']['shadow_archetype'] == 'M'
        assert summary['betagamma']['shadow_archetype'] == 'C'
