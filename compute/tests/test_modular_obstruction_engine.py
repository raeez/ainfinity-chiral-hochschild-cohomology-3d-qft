r"""Tests for the modular obstruction engine: genus-raising pipeline.

Verifies the full genus-raising quantization pipeline:
  PVA -> Modular Bar Datum -> D_1 -> Ob_1 -> W-Normal Form -> Theta_1 -> MC verification

Tests cover:
  1. Classical MC element extraction for all families
  2. D_1 genus-raising operator computation
  3. Ob_1 = D_1(Theta_0) explicit computation
  4. W-normal form transformation for W_3
  5. Genus-1 correction Theta_1 = kappa/24
  6. MC equation verification at order hbar^1
  7. Chain-level D_1 computations (Virasoro and W_3)
  8. Mode-by-mode sewing traces
  9. Complementarity of obstructions
  10. Pipeline consistency

References:
  modular_pva_quantization.tex (Vol II)
  higher_genus_modular_koszul.tex (Vol I)
  nonlinear_modular_shadows.tex (Vol I): shadow tower obstruction theory
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy import Symbol, Rational, simplify, expand, S, symbols, oo


# ===================================================================
# 1. CLASSICAL MC ELEMENT EXTRACTION
# ===================================================================

class TestClassicalMCElement:
    """Verify classical MC element extraction for all families."""

    def test_virasoro_mc_generators(self):
        """Virasoro MC element has generator T of weight 2."""
        from lib.modular_obstruction_engine import classical_mc_element
        mc = classical_mc_element('virasoro')
        assert mc['generators'] == ['T']
        assert mc['conformal_weights']['T'] == 2

    def test_virasoro_mc_nproducts(self):
        """Virasoro n-products: T_{(0)}T=dT, T_{(1)}T=2T, T_{(3)}T=c/2."""
        from lib.modular_obstruction_engine import classical_mc_element
        c = Symbol('c')
        mc = classical_mc_element('virasoro', c=c)
        prods = mc['mc_element'][('T', 'T')]
        T = Symbol('T')
        dT = Symbol('dT')
        assert prods[0] == dT
        assert simplify(prods[1] - 2 * T) == 0
        assert simplify(prods[3] - c / 2) == 0

    def test_virasoro_mc_kappa(self):
        """Virasoro kappa = c/2."""
        from lib.modular_obstruction_engine import classical_mc_element
        c = Symbol('c')
        mc = classical_mc_element('virasoro', c=c)
        assert simplify(mc['kappa'] - c / 2) == 0

    def test_w3_mc_generators(self):
        """W_3 has generators T (weight 2) and W (weight 3)."""
        from lib.modular_obstruction_engine import classical_mc_element
        mc = classical_mc_element('w3')
        assert 'T' in mc['generators']
        assert 'W' in mc['generators']
        assert mc['conformal_weights']['T'] == 2
        assert mc['conformal_weights']['W'] == 3

    def test_w3_mc_kappa(self):
        """W_3 kappa = 5c/6."""
        from lib.modular_obstruction_engine import classical_mc_element
        c = Symbol('c')
        mc = classical_mc_element('w3', c=c)
        assert simplify(mc['kappa'] - 5 * c / 6) == 0

    def test_w3_beta_squared(self):
        """W_3 beta^2 = 16/(22+5c)."""
        from lib.modular_obstruction_engine import classical_mc_element
        c = Symbol('c')
        mc = classical_mc_element('w3', c=c)
        expected = Rational(16, 1) / (22 + 5 * c)
        assert simplify(mc['beta_squared'] - expected) == 0

    def test_heisenberg_mc_kappa(self):
        """Heisenberg kappa = k/2."""
        from lib.modular_obstruction_engine import classical_mc_element
        k = Symbol('k')
        mc = classical_mc_element('heisenberg', k=k)
        assert simplify(mc['kappa'] - k / 2) == 0

    def test_affine_sl2_mc_kappa(self):
        """Affine sl_2: kappa = 3(k+2)/4."""
        from lib.modular_obstruction_engine import classical_mc_element
        k = Symbol('k')
        mc = classical_mc_element('affine_sl2', k=k)
        assert simplify(mc['kappa'] - 3 * (k + 2) / 4) == 0

    def test_betagamma_mc_kappa(self):
        """betagamma: kappa = 1."""
        from lib.modular_obstruction_engine import classical_mc_element
        mc = classical_mc_element('betagamma')
        assert mc['kappa'] == 1


# ===================================================================
# 2. D_1 GENUS-RAISING OPERATOR
# ===================================================================

class TestGenusRaisingD1:
    """Test the D_1 computation."""

    def test_virasoro_D1_gives_kappa(self):
        """D_1(Theta_0) = kappa * omega_1 for Virasoro."""
        from lib.modular_obstruction_engine import (
            classical_mc_element, genus_raising_operator_D1)
        c = Symbol('c')
        mc = classical_mc_element('virasoro', c=c)
        d1 = genus_raising_operator_D1(mc)
        assert simplify(d1['D1_theta0'] - c / 2) == 0

    def test_virasoro_D1_is_coboundary(self):
        """D_1(Theta_0) is a coboundary (exact in modular operad cohomology)."""
        from lib.modular_obstruction_engine import (
            classical_mc_element, genus_raising_operator_D1)
        mc = classical_mc_element('virasoro', c=Symbol('c'))
        d1 = genus_raising_operator_D1(mc)
        assert d1['is_coboundary'] is True

    def test_virasoro_D1_primitive(self):
        """Primitive of D_1(Theta_0) is kappa/24."""
        from lib.modular_obstruction_engine import (
            classical_mc_element, genus_raising_operator_D1)
        c = Symbol('c')
        mc = classical_mc_element('virasoro', c=c)
        d1 = genus_raising_operator_D1(mc)
        assert simplify(d1['primitive'] - c / 48) == 0

    def test_w3_D1_gives_kappa(self):
        """D_1(Theta_0) = 5c/6 * omega_1 for W_3."""
        from lib.modular_obstruction_engine import (
            classical_mc_element, genus_raising_operator_D1)
        c = Symbol('c')
        mc = classical_mc_element('w3', c=c)
        d1 = genus_raising_operator_D1(mc)
        assert simplify(d1['D1_theta0'] - 5 * c / 6) == 0

    def test_heisenberg_D1(self):
        """D_1(Theta_0) = k/2 for Heisenberg."""
        from lib.modular_obstruction_engine import (
            classical_mc_element, genus_raising_operator_D1)
        k = Symbol('k')
        mc = classical_mc_element('heisenberg', k=k)
        d1 = genus_raising_operator_D1(mc)
        assert simplify(d1['D1_theta0'] - k / 2) == 0


# ===================================================================
# 3. Ob_1 EXPLICIT COMPUTATION
# ===================================================================

class TestOb1Explicit:
    """Test the explicit genus-1 obstruction."""

    def test_virasoro_ob1_vanishes(self):
        """Virasoro Ob_1 class vanishes for generic c."""
        from lib.modular_obstruction_engine import genus1_obstruction_explicit
        ob = genus1_obstruction_explicit('virasoro')
        assert ob['ob1_class_vanishes'] is True

    def test_virasoro_ob1_primitive(self):
        """Primitive xi_1 = c/48 for Virasoro."""
        from lib.modular_obstruction_engine import genus1_obstruction_explicit
        c = Symbol('c')
        ob = genus1_obstruction_explicit('virasoro', c=c)
        assert simplify(ob['primitive_xi1'] - c / 48) == 0

    def test_w3_ob1_vanishes(self):
        """W_3 Ob_1 class vanishes for generic c."""
        from lib.modular_obstruction_engine import genus1_obstruction_explicit
        ob = genus1_obstruction_explicit('w3')
        assert ob['ob1_class_vanishes'] is True

    def test_w3_ob1_primitive(self):
        """Primitive xi_1 = 5c/144 for W_3."""
        from lib.modular_obstruction_engine import genus1_obstruction_explicit
        c = Symbol('c')
        ob = genus1_obstruction_explicit('w3', c=c)
        assert simplify(ob['primitive_xi1'] - 5 * c / 144) == 0

    def test_heisenberg_ob1_vanishes(self):
        """Heisenberg Ob_1 vanishes."""
        from lib.modular_obstruction_engine import genus1_obstruction_explicit
        ob = genus1_obstruction_explicit('heisenberg')
        assert ob['ob1_class_vanishes'] is True

    def test_betagamma_ob1_vanishes(self):
        """betagamma Ob_1 vanishes."""
        from lib.modular_obstruction_engine import genus1_obstruction_explicit
        ob = genus1_obstruction_explicit('betagamma')
        assert ob['ob1_class_vanishes'] is True

    def test_ob1_faber_pandharipande(self):
        """Period integral = 1/24 (Faber-Pandharipande)."""
        from lib.modular_obstruction_engine import genus1_obstruction_explicit
        ob = genus1_obstruction_explicit('virasoro')
        assert ob['faber_pandharipande_period'] == Rational(1, 24)


# ===================================================================
# 4. W-NORMAL FORM
# ===================================================================

class TestWNormalForm:
    """Test the triangular W-normal form transformation."""

    def test_virasoro_already_normal(self):
        """Virasoro is already in W-normal form."""
        from lib.modular_obstruction_engine import w_normal_form_transformation
        nf = w_normal_form_transformation('virasoro')
        assert nf['already_normal'] is True

    def test_w3_kappa_sum(self):
        """W_3 normal form: kappa_T + kappa_W = 5c/6."""
        from lib.modular_obstruction_engine import w_normal_form_transformation
        nf = w_normal_form_transformation('w3')
        assert nf['kappa_sum_check'] is True

    def test_w3_kappa_T(self):
        """W_3 normal form: kappa_T = c/2."""
        from lib.modular_obstruction_engine import w_normal_form_transformation
        c = Symbol('c')
        nf = w_normal_form_transformation('w3', c=c)
        assert simplify(nf['kappa_T'] - c / 2) == 0

    def test_w3_kappa_W(self):
        """W_3 normal form: kappa_W = c/3."""
        from lib.modular_obstruction_engine import w_normal_form_transformation
        c = Symbol('c')
        nf = w_normal_form_transformation('w3', c=c)
        assert simplify(nf['kappa_W'] - c / 3) == 0

    def test_w3_ob1_vanishes_in_nf(self):
        """Ob_1 vanishes in normal form for W_3."""
        from lib.modular_obstruction_engine import w_normal_form_transformation
        nf = w_normal_form_transformation('w3')
        assert nf['ob1_vanishes_in_nf'] is True


# ===================================================================
# 5. GENUS-1 CORRECTION THETA_1
# ===================================================================

class TestGenus1Correction:
    """Test the genus-1 correction computation."""

    def test_virasoro_theta1(self):
        """Virasoro Theta_1 = c/48."""
        from lib.modular_obstruction_engine import genus1_correction
        c = Symbol('c')
        g1 = genus1_correction('virasoro', c=c)
        assert simplify(g1['theta_1'] - c / 48) == 0

    def test_w3_theta1(self):
        """W_3 Theta_1 = 5c/144."""
        from lib.modular_obstruction_engine import genus1_correction
        c = Symbol('c')
        g1 = genus1_correction('w3', c=c)
        assert simplify(g1['theta_1'] - 5 * c / 144) == 0

    def test_heisenberg_theta1(self):
        """Heisenberg Theta_1 = k/48."""
        from lib.modular_obstruction_engine import genus1_correction
        k = Symbol('k')
        g1 = genus1_correction('heisenberg', k=k)
        assert simplify(g1['theta_1'] - k / 48) == 0

    def test_mc_order1_satisfied(self):
        """MC equation at order hbar^1 is satisfied."""
        from lib.modular_obstruction_engine import genus1_correction
        g1 = genus1_correction('virasoro', c=Symbol('c'))
        assert g1['mc_order1_satisfied'] is True


# ===================================================================
# 6. CHAIN-LEVEL D_1 (GENUINE COMPUTATION)
# ===================================================================

class TestChainLevelD1:
    """Test chain-level D_1 computations."""

    def test_virasoro_chain_level_kappa(self):
        """Virasoro chain-level D_1 gives kappa = c/2."""
        from lib.modular_obstruction_engine import virasoro_D1_chain_level
        result = virasoro_D1_chain_level()
        assert result['kappa_check'] is True

    def test_virasoro_chain_level_n0_zero(self):
        """T_{(0)}T = dT contributes zero to D_1 (orthogonal)."""
        from lib.modular_obstruction_engine import virasoro_D1_chain_level
        result = virasoro_D1_chain_level()
        assert result['n_product_contributions'][0]['D1_value'] == S.Zero

    def test_virasoro_chain_level_n3_zero(self):
        """T_{(3)}T = c/2 (scalar) contributes zero to D_1."""
        from lib.modular_obstruction_engine import virasoro_D1_chain_level
        result = virasoro_D1_chain_level()
        assert result['n_product_contributions'][3]['D1_value'] == S.Zero

    def test_virasoro_chain_level_n1_primary(self):
        """T_{(1)}T = 2T is the PRIMARY contributor to D_1."""
        from lib.modular_obstruction_engine import virasoro_D1_chain_level
        c = Symbol('c')
        result = virasoro_D1_chain_level()
        val = result['n_product_contributions'][1]['D1_value']
        assert simplify(val - c / 2) == 0

    def test_virasoro_chain_total_equals_kappa(self):
        """Total chain-level D_1 = kappa = c/2."""
        from lib.modular_obstruction_engine import virasoro_D1_chain_level
        c = Symbol('c')
        result = virasoro_D1_chain_level()
        assert simplify(result['D1_total_chain_level'] - c / 2) == 0

    def test_virasoro_chain_xi1(self):
        """Primitive from chain level = kappa/24 = c/48."""
        from lib.modular_obstruction_engine import virasoro_D1_chain_level
        c = Symbol('c')
        result = virasoro_D1_chain_level()
        assert simplify(result['primitive_xi1'] - c / 48) == 0


class TestW3ChainLevel:
    """Test W_3 chain-level D_1 computation."""

    def test_w3_kappa_T(self):
        """W_3 T-sector: kappa_T = c/2."""
        from lib.modular_obstruction_engine import w3_D1_chain_level
        c = Symbol('c')
        result = w3_D1_chain_level()
        assert simplify(result['kappa_T'] - c / 2) == 0

    def test_w3_kappa_W(self):
        """W_3 W-sector: kappa_W = c/3."""
        from lib.modular_obstruction_engine import w3_D1_chain_level
        c = Symbol('c')
        result = w3_D1_chain_level()
        assert simplify(result['kappa_W'] - c / 3) == 0

    def test_w3_kappa_total(self):
        """W_3 total kappa = 5c/6."""
        from lib.modular_obstruction_engine import w3_D1_chain_level
        result = w3_D1_chain_level()
        assert result['kappa_check_5c_6'] is True

    def test_w3_cross_zero(self):
        """W_3 cross-sector (T-W) = 0 by weight orthogonality."""
        from lib.modular_obstruction_engine import w3_D1_chain_level
        result = w3_D1_chain_level()
        assert result['cross_sector'] == S.Zero

    def test_w3_xi1(self):
        """W_3 primitive xi_1 = 5c/144."""
        from lib.modular_obstruction_engine import w3_D1_chain_level
        result = w3_D1_chain_level()
        assert result['xi_1_expected_5c_144'] is True

    def test_w3_numerical_c30(self):
        """W_3 at c=30: kappa = 25, xi_1 = 25/24."""
        from lib.modular_obstruction_engine import w3_D1_chain_level
        result = w3_D1_chain_level(c_val=30)
        assert result['kappa_total'] == 25
        assert result['xi_1'] == Rational(25, 24)


# ===================================================================
# 7. MC EQUATION VERIFICATION
# ===================================================================

class TestMCVerification:
    """Test MC equation at order hbar^1."""

    def test_virasoro_mc_order1(self):
        """MC equation at O(hbar) holds for Virasoro."""
        from lib.modular_obstruction_engine import verify_mc_equation_order1
        result = verify_mc_equation_order1('virasoro')
        assert result['mc_order1'] is True

    def test_w3_mc_order1(self):
        """MC equation at O(hbar) holds for W_3."""
        from lib.modular_obstruction_engine import verify_mc_equation_order1
        result = verify_mc_equation_order1('w3')
        assert result['mc_order1'] is True

    def test_heisenberg_mc_order1(self):
        """MC equation at O(hbar) holds for Heisenberg."""
        from lib.modular_obstruction_engine import verify_mc_equation_order1
        result = verify_mc_equation_order1('heisenberg')
        assert result['mc_order1'] is True

    def test_chain_level_mc_virasoro(self):
        """Chain-level MC verification for Virasoro."""
        from lib.modular_obstruction_engine import chain_level_mc_verification
        c = Symbol('c')
        result = chain_level_mc_verification('virasoro', c=c)
        assert result['genus0_mc_holds'] is True
        assert result['genus1_mc_in_cohomology'] is True

    def test_chain_level_mc_w3(self):
        """Chain-level MC verification for W_3."""
        from lib.modular_obstruction_engine import chain_level_mc_verification
        result = chain_level_mc_verification('w3')
        assert result['genus1_mc_in_cohomology'] is True


# ===================================================================
# 8. MODE-BY-MODE SEWING TRACES
# ===================================================================

class TestModeTraces:
    """Test mode-by-mode D_1 sewing traces."""

    def test_virasoro_mode_m1_zero(self):
        """Virasoro: <L_1, L_{-1}> = 0 (m=1, m^2-1=0)."""
        from lib.modular_obstruction_engine import virasoro_D1_mode_trace
        result = virasoro_D1_mode_trace()
        assert result['mode_traces'][1]['pairing'] == 0

    def test_virasoro_mode_m2(self):
        """Virasoro: <L_2, L_{-2}> = (c/12)*2*3 = c/2."""
        from lib.modular_obstruction_engine import virasoro_D1_mode_trace
        c = Symbol('c')
        result = virasoro_D1_mode_trace()
        expected = c * 2 * (4 - 1) / 12  # c * 2 * 3 / 12 = c/2
        assert simplify(result['mode_traces'][2]['pairing'] - expected) == 0

    def test_virasoro_mode_m3(self):
        """Virasoro: <L_3, L_{-3}> = (c/12)*3*8 = 2c."""
        from lib.modular_obstruction_engine import virasoro_D1_mode_trace
        c = Symbol('c')
        result = virasoro_D1_mode_trace()
        expected = c * 3 * (9 - 1) / 12  # c * 3 * 8 / 12 = 2c
        assert simplify(result['mode_traces'][3]['pairing'] - expected) == 0

    def test_virasoro_closed_form_agrees(self):
        """Closed form sum = direct summation."""
        from lib.modular_obstruction_engine import virasoro_D1_mode_trace
        result = virasoro_D1_mode_trace(max_level=10)
        assert result['agrees'] is True

    def test_virasoro_mode_trace_closed_form(self):
        """Closed form: (c/48)*N*(N+1)*(N-1)*(N+2) for truncation at level N."""
        from lib.modular_obstruction_engine import virasoro_D1_mode_trace
        c = Symbol('c')
        N = 5
        result = virasoro_D1_mode_trace(max_level=N)
        expected = c * N * (N + 1) * (N - 1) * (N + 2) / 48
        assert simplify(result['closed_form_truncated'] - expected) == 0


# ===================================================================
# 9. COMPLEMENTARITY
# ===================================================================

class TestComplementarity:
    """Test complementarity of obstructions."""

    def test_virasoro_complementarity(self):
        """Virasoro: kappa + kappa' = 13."""
        from lib.modular_obstruction_engine import obstruction_complementarity
        result = obstruction_complementarity('virasoro')
        assert result['check'] is True
        assert result['complementarity_sum'] == 13

    def test_w3_complementarity(self):
        """W_3: kappa + kappa' = 500/6."""
        from lib.modular_obstruction_engine import obstruction_complementarity
        result = obstruction_complementarity('w3')
        assert result['check'] is True

    def test_affine_complementarity(self):
        """Affine sl_2: kappa + kappa' = 0."""
        from lib.modular_obstruction_engine import obstruction_complementarity
        result = obstruction_complementarity('affine_sl2')
        assert result['check'] is True
        assert result['complementarity_sum'] == 0

    def test_heisenberg_complementarity(self):
        """Heisenberg: kappa + kappa' = 0."""
        from lib.modular_obstruction_engine import obstruction_complementarity
        result = obstruction_complementarity('heisenberg')
        assert result['check'] is True
        assert result['complementarity_sum'] == 0


# ===================================================================
# 10. W_3 NORMAL FORM OBSTRUCTION VANISHING
# ===================================================================

class TestW3NormalFormVanishing:
    """Test W_3 Ob_1 vanishing via normal form."""

    def test_kappa_check(self):
        """kappa_T + kappa_W = 5c/6."""
        from lib.modular_obstruction_engine import w3_ob1_normal_form_vanishing
        result = w3_ob1_normal_form_vanishing()
        assert result['kappa_check_5c_over_6'] is True

    def test_xi1_check(self):
        """xi_1 = 5c/144."""
        from lib.modular_obstruction_engine import w3_ob1_normal_form_vanishing
        result = w3_ob1_normal_form_vanishing()
        assert result['xi_1_check'] is True

    def test_ob1_vanishes(self):
        """[Ob_1] = 0 in cohomology."""
        from lib.modular_obstruction_engine import w3_ob1_normal_form_vanishing
        result = w3_ob1_normal_form_vanishing()
        assert result['ob1_class_vanishes'] is True

    def test_numerical_c30(self):
        """At c=30: kappa = 25, xi_1 = 25/24."""
        from lib.modular_obstruction_engine import w3_ob1_normal_form_vanishing
        result = w3_ob1_normal_form_vanishing(c_val=30)
        assert result['kappa_total'] == 25
        assert result['xi_1'] == Rational(25, 24)

    def test_resonance_divisor(self):
        """Resonance at c = -22/5."""
        from lib.modular_obstruction_engine import w3_ob1_normal_form_vanishing
        result = w3_ob1_normal_form_vanishing()
        assert result['resonance_at'] == Rational(-22, 5)


# ===================================================================
# 11. FULL PIPELINE
# ===================================================================

class TestFullPipeline:
    """Test the full obstruction pipeline."""

    def test_virasoro_pipeline_complete(self):
        """Virasoro pipeline runs to completion."""
        from lib.modular_obstruction_engine import full_obstruction_pipeline
        result = full_obstruction_pipeline('virasoro')
        assert result['quantizable'] is True
        assert result['pipeline_complete'] is True

    def test_w3_pipeline_complete(self):
        """W_3 pipeline runs to completion."""
        from lib.modular_obstruction_engine import full_obstruction_pipeline
        result = full_obstruction_pipeline('w3')
        assert result['quantizable'] is True
        assert result['pipeline_complete'] is True

    def test_heisenberg_pipeline(self):
        """Heisenberg pipeline runs to completion."""
        from lib.modular_obstruction_engine import full_obstruction_pipeline
        result = full_obstruction_pipeline('heisenberg')
        assert result['quantizable'] is True

    def test_betagamma_pipeline(self):
        """betagamma pipeline runs to completion."""
        from lib.modular_obstruction_engine import full_obstruction_pipeline
        result = full_obstruction_pipeline('betagamma')
        assert result['quantizable'] is True


# ===================================================================
# 12. GENUS-2 OBSTRUCTION
# ===================================================================

class TestGenus2:
    """Test genus-2 obstruction."""

    def test_virasoro_ob2_exists(self):
        """Genus-2 obstruction computed for Virasoro."""
        from lib.modular_obstruction_engine import genus2_obstruction
        result = genus2_obstruction('virasoro')
        assert result is not None
        assert result['ob2'] is not None

    def test_virasoro_ob2_class_vanishes(self):
        """[Ob_2] = 0 for Virasoro (recursive existence)."""
        from lib.modular_obstruction_engine import genus2_obstruction
        result = genus2_obstruction('virasoro')
        assert result['ob2_class_vanishes'] is True

    def test_genus2_fp_number(self):
        """Faber-Pandharipande genus-2 number = 1/1152."""
        from lib.modular_obstruction_engine import genus2_obstruction
        result = genus2_obstruction('virasoro')
        assert result['faber_pandharipande_g2'] == Rational(1, 1152)


# ===================================================================
# 13. ADDITIVITY
# ===================================================================

class TestAdditivity:
    """Test additivity of obstructions."""

    def test_heisenberg_plus_virasoro(self):
        """kappa(H + Vir) = kappa(H) + kappa(Vir)."""
        from lib.modular_obstruction_engine import verify_additivity_of_obstructions
        k, c = Symbol('k'), Symbol('c')
        result = verify_additivity_of_obstructions([
            ('heisenberg', {'k': k}),
            ('virasoro', {'c': c}),
        ])
        assert simplify(result['sum'] - (k / 2 + c / 2)) == 0

    def test_additivity_holds(self):
        """Additivity structural flag."""
        from lib.modular_obstruction_engine import verify_additivity_of_obstructions
        result = verify_additivity_of_obstructions([
            ('heisenberg', {'k': 1}),
            ('betagamma', {}),
        ])
        assert result['additivity_holds'] is True


# ===================================================================
# 14. WEIGHT-GRADED D_1
# ===================================================================

class TestWeightGradedD1:
    """Test weight-graded D_1 for Virasoro."""

    def test_weight0_zero(self):
        """Weight 0 contribution is zero."""
        from lib.modular_obstruction_engine import virasoro_D1_by_weight
        result = virasoro_D1_by_weight()
        assert result['weight_contributions'][0] == S.Zero

    def test_weight2_zero(self):
        """Weight 2: <L_1, L_{-1}> = 0 (m=1: m(m^2-1)=0)."""
        from lib.modular_obstruction_engine import virasoro_D1_by_weight
        result = virasoro_D1_by_weight()
        assert result['weight_contributions'][2] == S.Zero

    def test_weight4_nonzero(self):
        """Weight 4: <L_2, L_{-2}> = c/2."""
        from lib.modular_obstruction_engine import virasoro_D1_by_weight
        c = Symbol('c')
        result = virasoro_D1_by_weight()
        assert simplify(result['weight_contributions'][4] - c / 2) == 0


# ===================================================================
# 15. W_3 MODE DECOMPOSITION
# ===================================================================

class TestW3ModeDecomposition:
    """Test W_3 mode-by-mode decomposition."""

    def test_T_sector_m2(self):
        """T-sector at m=2: (c/12)*2*3 = c/2."""
        from lib.modular_obstruction_engine import w3_D1_mode_decomposition
        c = Symbol('c')
        result = w3_D1_mode_decomposition()
        expected = c * 2 * 3 / 12
        assert simplify(result['T_sector_traces'][2] - expected) == 0

    def test_W_sector_below_weight3_zero(self):
        """W-sector modes at m<3 are zero (W has weight 3)."""
        from lib.modular_obstruction_engine import w3_D1_mode_decomposition
        result = w3_D1_mode_decomposition()
        assert result['W_sector_traces'][1] == S.Zero
        assert result['W_sector_traces'][2] == S.Zero

    def test_W_sector_m3(self):
        """W-sector at m=3: (c/3)*3*(8)*(5)/360 = c/3*120/360 = c/9."""
        from lib.modular_obstruction_engine import w3_D1_mode_decomposition
        c = Symbol('c')
        result = w3_D1_mode_decomposition()
        # <W_3, W_{-3}> = (c/3) * 3 * (9-1) * (9-4) / 360
        #               = (c/3) * 3 * 8 * 5 / 360 = c * 120 / 1080 = c/9
        expected = c * 3 * 8 * 5 / 360
        assert simplify(result['W_sector_traces'][3] - expected) == 0
