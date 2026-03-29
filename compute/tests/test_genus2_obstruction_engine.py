r"""Tests for the genus-2 obstruction engine.

Verifies the full genus-2 MC pipeline:
  1. Faber-Pandharipande numbers lambda_g^FP
  2. Stable graph enumeration at genus 2
  3. Graph amplitudes
  4. Genus-2 obstruction Ob_2 = D_1(Theta_1) + (1/2)*D_2(Theta_0)
  5. Genus-2 free energy F_2 for all families
  6. F_2/F_1 and F_2/F_1^2 universal ratios
  7. Genus-2 MC element assembly
  8. Complementarity at genus 2
  9. Generating function verification
  10. Obstruction tower consistency

References:
  higher_genus_modular_koszul.tex (Vol I): modular bar, genus spectral sequence
  modular_obstruction_engine.py (Vol II): genus-1 pipeline
  genus2_obstruction_engine.py
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy import Symbol, Rational, simplify, expand, S, bernoulli


# ===================================================================
# 1. FABER-PANDHARIPANDE NUMBERS
# ===================================================================

class TestFaberPandharipande:
    """Verify Faber-Pandharipande intersection numbers."""

    def test_lambda_1_is_1_over_24(self):
        """lambda_1^FP = 1/24."""
        from lib.genus2_obstruction_engine import lambda_fp
        assert lambda_fp(1) == Rational(1, 24)

    def test_lambda_2_is_7_over_5760(self):
        """lambda_2^FP = 7/5760."""
        from lib.genus2_obstruction_engine import lambda_fp
        assert lambda_fp(2) == Rational(7, 5760)

    def test_lambda_3(self):
        """lambda_3^FP = 31/967680."""
        from lib.genus2_obstruction_engine import lambda_fp
        assert lambda_fp(3) == Rational(31, 967680)

    def test_lambda_all_positive(self):
        """All lambda_g^FP are positive (Bernoulli sign pattern)."""
        from lib.genus2_obstruction_engine import lambda_fp
        for g in range(1, 8):
            assert lambda_fp(g) > 0, f"lambda_{g} = {lambda_fp(g)} <= 0"

    def test_lambda_decreasing(self):
        """lambda_g^FP decreases super-exponentially."""
        from lib.genus2_obstruction_engine import lambda_fp
        for g in range(1, 6):
            assert lambda_fp(g + 1) < lambda_fp(g), (
                f"lambda_{g+1} >= lambda_{g}")

    def test_lambda_from_bernoulli(self):
        """Verify formula: lambda_g = (2^{2g-1}-1)|B_{2g}| / (2^{2g-1}(2g)!)."""
        from lib.genus2_obstruction_engine import lambda_fp
        from sympy import factorial
        for g in range(1, 5):
            B_2g = bernoulli(2 * g)
            expected = (2**(2*g-1) - 1) * abs(B_2g) / (2**(2*g-1) * factorial(2*g))
            assert lambda_fp(g) == Rational(expected), (
                f"g={g}: {lambda_fp(g)} != {expected}")

    def test_lambda_genus0_raises(self):
        """lambda_fp(0) should raise ValueError."""
        from lib.genus2_obstruction_engine import lambda_fp
        try:
            lambda_fp(0)
            assert False, "Should have raised ValueError"
        except ValueError:
            pass


# ===================================================================
# 2. FREE ENERGY F_g
# ===================================================================

class TestFreeEnergy:
    """Test genus-g free energy F_g = kappa * lambda_g."""

    def test_F1_virasoro(self):
        """F_1(Vir_c) = c/48."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        c = Symbol('c')
        kappa = c / 2
        F1 = F_g_free_energy(kappa, 1)
        assert simplify(F1 - c / 48) == 0

    def test_F2_virasoro(self):
        """F_2(Vir_c) = 7c/11520."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        c = Symbol('c')
        kappa = c / 2
        F2 = F_g_free_energy(kappa, 2)
        assert simplify(F2 - 7 * c / 11520) == 0

    def test_F2_heisenberg(self):
        """F_2(H_k) = 7k/5760 (κ=k)."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        k = Symbol('k')
        kappa = k
        F2 = F_g_free_energy(kappa, 2)
        assert simplify(F2 - 7 * k / 5760) == 0

    def test_F2_w3(self):
        """F_2(W_3) = 7c/6912 (since kappa = 5c/6)."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        c = Symbol('c')
        kappa = 5 * c / 6
        F2 = F_g_free_energy(kappa, 2)
        expected = 5 * c / 6 * Rational(7, 5760)
        assert simplify(F2 - expected) == 0

    def test_F2_betagamma(self):
        """F_2(betagamma) = 7/5760 (kappa = 1)."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        F2 = F_g_free_energy(S.One, 2)
        assert F2 == Rational(7, 5760)


# ===================================================================
# 3. STABLE GRAPHS AT GENUS 2
# ===================================================================

class TestStableGraphs:
    """Test genus-2 stable graph enumeration."""

    def test_six_graphs_at_genus_2(self):
        """There are exactly 6 stable graph types at (g=2, n=0)."""
        from lib.genus2_obstruction_engine import genus2_stable_graphs_n0
        graphs = genus2_stable_graphs_n0()
        assert len(graphs) == 6

    def test_all_genus_2(self):
        """All graphs have total arithmetic genus 2."""
        from lib.genus2_obstruction_engine import verify_genus2_graph_count
        results = verify_genus2_graph_count()
        for name, data in results.items():
            assert data['genus_correct'], f"Graph {name}: genus = {data['total_genus']} != 2"

    def test_h1_consistency(self):
        """h^1 stored in graph matches computed value."""
        from lib.genus2_obstruction_engine import verify_genus2_graph_count
        results = verify_genus2_graph_count()
        for name, data in results.items():
            assert data['h1_consistent'], f"Graph {name}: h1 mismatch"

    def test_smooth_graph_has_no_edges(self):
        """Graph I (smooth) has 0 edges."""
        from lib.genus2_obstruction_engine import genus2_stable_graphs_n0
        graphs = genus2_stable_graphs_n0()
        smooth = graphs[0]
        assert smooth.n_edges == 0
        assert smooth.vertex_genera == (2,)

    def test_separating_graph(self):
        """Graph III (separating) has 2 genus-1 vertices."""
        from lib.genus2_obstruction_engine import genus2_stable_graphs_n0
        graphs = genus2_stable_graphs_n0()
        sep = graphs[2]  # Graph III
        assert sep.vertex_genera == (1, 1)
        assert sep.n_edges == 1
        assert sep.degeneration_type == 'separating'

    def test_three_banana_automorphisms(self):
        """Graph VI (three-banana) has |Aut| = 12."""
        from lib.genus2_obstruction_engine import genus2_stable_graphs_n0
        graphs = genus2_stable_graphs_n0()
        banana = graphs[5]  # Graph VI
        assert banana.aut_order == 12
        assert banana.n_edges == 3

    def test_nonseparating_self_edge(self):
        """Graph II has 1 self-edge on genus-1 vertex."""
        from lib.genus2_obstruction_engine import genus2_stable_graphs_n0
        graphs = genus2_stable_graphs_n0()
        g2 = graphs[1]  # Graph II
        assert g2.n_self_edges == 1
        assert g2.vertex_genera == (1,)
        assert g2.aut_order == 2

    def test_graph_IV_two_self_edges(self):
        """Graph IV has 2 self-edges on genus-0 vertex."""
        from lib.genus2_obstruction_engine import genus2_stable_graphs_n0
        graphs = genus2_stable_graphs_n0()
        g4 = graphs[3]  # Graph IV
        assert g4.n_self_edges == 2
        assert g4.vertex_genera == (0,)
        assert g4.aut_order == 8


# ===================================================================
# 4. GRAPH AMPLITUDES
# ===================================================================

class TestGraphAmplitudes:
    """Test graph amplitude computations."""

    def test_smooth_amplitude_is_F2(self):
        """Graph I amplitude = F_2 = kappa * 7/5760."""
        from lib.genus2_obstruction_engine import (
            genus2_stable_graphs_n0, graph_amplitude_genus2, F_g_free_energy)
        kappa = Symbol('kappa')
        graphs = genus2_stable_graphs_n0()
        amp = graph_amplitude_genus2(graphs[0], kappa)
        expected = F_g_free_energy(kappa, 2)
        assert simplify(amp['amplitude'] - expected) == 0

    def test_smooth_amplitude_weighted(self):
        """Graph I weighted amplitude = F_2 / 1 (aut_order = 1)."""
        from lib.genus2_obstruction_engine import (
            genus2_stable_graphs_n0, graph_amplitude_genus2)
        kappa = Symbol('kappa')
        graphs = genus2_stable_graphs_n0()
        amp = graph_amplitude_genus2(graphs[0], kappa)
        assert amp['aut_order'] == 1
        assert simplify(amp['weighted_amplitude'] - amp['amplitude']) == 0

    def test_graph_II_amplitude(self):
        """Graph II amplitude involves kappa * P = 1."""
        from lib.genus2_obstruction_engine import (
            genus2_stable_graphs_n0, graph_amplitude_genus2)
        kappa = Symbol('kappa')
        graphs = genus2_stable_graphs_n0()
        amp = graph_amplitude_genus2(graphs[1], kappa)
        # amplitude = kappa * (1/kappa) = 1
        assert simplify(amp['amplitude'] - 1) == 0


# ===================================================================
# 5. GENUS-2 OBSTRUCTION Ob_2
# ===================================================================

class TestGenus2Obstruction:
    """Test the genus-2 obstruction computation."""

    def test_D1_theta1_virasoro(self):
        """D_1(Theta_1) for Virasoro."""
        from lib.genus2_obstruction_engine import genus2_obstruction_D1_theta1
        c = Symbol('c')
        kappa = c / 2
        result = genus2_obstruction_D1_theta1(kappa)
        # D_1(Theta_1) = kappa/24 = c/48
        assert simplify(result['D1_theta1'] - c / 48) == 0

    def test_D2_theta0_components(self):
        """D_2(Theta_0) has separating and nonseparating parts."""
        from lib.genus2_obstruction_engine import genus2_obstruction_D2_theta0
        kappa = Symbol('kappa')
        result = genus2_obstruction_D2_theta0(kappa)
        assert result['D2_sep'] is not None
        assert result['D2_nonsep'] is not None

    def test_ob2_virasoro_full(self):
        """Full Ob_2 for Virasoro is computed."""
        from lib.genus2_obstruction_engine import genus2_obstruction_full
        result = genus2_obstruction_full('virasoro')
        assert result['ob2'] is not None
        assert result['ob2_is_exact']

    def test_ob2_heisenberg_full(self):
        """Full Ob_2 for Heisenberg is computed."""
        from lib.genus2_obstruction_engine import genus2_obstruction_full
        result = genus2_obstruction_full('heisenberg')
        assert result['ob2'] is not None
        assert result['ob2_is_exact']

    def test_ob2_w3_full(self):
        """Full Ob_2 for W_3 is computed."""
        from lib.genus2_obstruction_engine import genus2_obstruction_full
        result = genus2_obstruction_full('w3')
        assert result['ob2'] is not None
        assert result['mc_order2_satisfied']


# ===================================================================
# 6. F_2 VERIFICATION FOR ALL FAMILIES
# ===================================================================

class TestF2Verification:
    """Verify F_2 for all standard families."""

    def test_F2_heisenberg_match(self):
        """F_2 for Heisenberg matches expected."""
        from lib.genus2_obstruction_engine import verify_F2_heisenberg
        result = verify_F2_heisenberg()
        assert result['match']

    def test_F2_virasoro_match(self):
        """F_2 for Virasoro matches expected."""
        from lib.genus2_obstruction_engine import verify_F2_virasoro
        result = verify_F2_virasoro()
        assert result['match']

    def test_F2_w3_match(self):
        """F_2 for W_3 matches expected."""
        from lib.genus2_obstruction_engine import verify_F2_w3
        result = verify_F2_w3()
        assert result['match']

    def test_F2_affine_sl2_match(self):
        """F_2 for affine sl_2 matches expected."""
        from lib.genus2_obstruction_engine import verify_F2_affine_sl2
        result = verify_F2_affine_sl2()
        assert result['match']

    def test_F2_heisenberg_k1(self):
        """F_2(H_1) = 7/5760 (κ=1 at k=1)."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        F2 = F_g_free_energy(S.One, 2)  # kappa = 1 for k=1
        assert F2 == Rational(7, 5760)

    def test_F2_virasoro_c26(self):
        """F_2(Vir_26) = 7*13/5760 = 91/5760."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        F2 = F_g_free_energy(S(13), 2)  # kappa = 26/2 = 13
        expected = 13 * Rational(7, 5760)
        assert F2 == expected

    def test_F2_betagamma_numeric(self):
        """F_2(betagamma) = 7/5760 as a float."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        F2 = F_g_free_energy(S.One, 2)
        assert abs(float(F2) - 7/5760) < 1e-15


# ===================================================================
# 7. F_2 / F_1 RATIOS
# ===================================================================

class TestF2Ratios:
    """Test universal ratios involving F_2."""

    def test_F2_over_F1_is_7_over_240(self):
        """F_2/F_1 = 7/240 (kappa-independent)."""
        from lib.genus2_obstruction_engine import F2_over_F1_ratio
        kappa = Symbol('kappa')
        result = F2_over_F1_ratio(kappa)
        assert result['match']
        assert result['kappa_independent']

    def test_F2_over_F1_squared_inversely_proportional(self):
        """F_2/F_1^2 = 7/(10*kappa)."""
        from lib.genus2_obstruction_engine import F2_over_F1_squared
        kappa = Symbol('kappa')
        result = F2_over_F1_squared(kappa)
        assert result['match']

    def test_F2_over_F1_squared_numeric(self):
        """F_2/F_1^2 at kappa=1 is 7/10."""
        from lib.genus2_obstruction_engine import F2_over_F1_squared
        result = F2_over_F1_squared(S.One)
        assert simplify(result['ratio'] - Rational(7, 10)) == 0

    def test_F2_over_F1_ratio_concrete(self):
        """F_2/F_1 = (kappa*7/5760)/(kappa/24) = 7*24/5760 = 7/240."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        kappa = Symbol('kappa')
        F1 = F_g_free_energy(kappa, 1)
        F2 = F_g_free_energy(kappa, 2)
        ratio = simplify(F2 / F1)
        assert ratio == Rational(7, 240)


# ===================================================================
# 8. GENUS-2 MC ELEMENT
# ===================================================================

class TestGenus2MC:
    """Test genus-2 MC element assembly."""

    def test_mc_element_virasoro(self):
        """MC element for Virasoro has three terms."""
        from lib.genus2_obstruction_engine import genus2_mc_element
        result = genus2_mc_element('virasoro')
        assert result['theta_0'] is not None
        assert result['theta_1'] is not None
        assert result['theta_2'] is not None

    def test_mc_element_theta0_is_kappa(self):
        """Theta_0 at scalar level is kappa."""
        from lib.genus2_obstruction_engine import genus2_mc_element
        c = Symbol('c')
        result = genus2_mc_element('virasoro', c=c)
        assert simplify(result['theta_0'] - c / 2) == 0

    def test_mc_element_theta1_is_kappa_over_24(self):
        """Theta_1 = kappa/24."""
        from lib.genus2_obstruction_engine import genus2_mc_element
        c = Symbol('c')
        result = genus2_mc_element('virasoro', c=c)
        assert simplify(result['theta_1'] - c / 48) == 0

    def test_mc_element_theta2_is_F2(self):
        """Theta_2 = F_2 = kappa * 7/5760."""
        from lib.genus2_obstruction_engine import genus2_mc_element
        c = Symbol('c')
        result = genus2_mc_element('virasoro', c=c)
        expected = c / 2 * Rational(7, 5760)
        assert simplify(result['theta_2'] - expected) == 0


# ===================================================================
# 9. COMPLEMENTARITY AT GENUS 2
# ===================================================================

class TestComplementarity:
    """Test Theorem C complementarity at genus 2."""

    def test_F2_complementarity_virasoro(self):
        """F_2(Vir_c) + F_2(Vir_{26-c}) is c-independent."""
        from lib.genus2_obstruction_engine import genus2_complementarity_check
        result = genus2_complementarity_check()
        assert result['F2_complementarity']

    def test_F1_complementarity_virasoro(self):
        """F_1(Vir_c) + F_1(Vir_{26-c}) = 13/24."""
        from lib.genus2_obstruction_engine import genus2_complementarity_check
        result = genus2_complementarity_check()
        assert result['F1_complementarity']

    def test_F2_sum_is_91_over_5760(self):
        """F_2(Vir_c) + F_2(Vir_{26-c}) = 91/5760."""
        from lib.genus2_obstruction_engine import genus2_complementarity_check
        result = genus2_complementarity_check()
        assert result['F2_sum_expected'] == Rational(91, 5760)

    def test_F2_at_self_dual_point(self):
        """At c=13 (self-dual): F_2(Vir_13) = 91/11520."""
        from lib.genus2_obstruction_engine import F_g_free_energy
        kappa_13 = Rational(13, 2)
        F2_13 = F_g_free_energy(kappa_13, 2)
        assert F2_13 == Rational(91, 11520)


# ===================================================================
# 10. GENERATING FUNCTION
# ===================================================================

class TestGeneratingFunction:
    """Test the A-hat generating function verification."""

    def test_ahat_coefficients(self):
        """A-hat generating function coefficients are correct."""
        from lib.genus2_obstruction_engine import ahat_generating_function_coefficients
        result = ahat_generating_function_coefficients(max_genus=4)
        assert result['all_positive']
        assert result['coefficients'][1] == Rational(1, 24)
        assert result['coefficients'][2] == Rational(7, 5760)

    def test_verify_ahat_through_genus4(self):
        """Verify A-hat generating function through genus 4."""
        from lib.genus2_obstruction_engine import verify_ahat_through_genus
        result = verify_ahat_through_genus(max_genus=4)
        for key, val in result['checks'].items():
            assert val, f"Check failed: {key}"

    def test_lambda_2_float(self):
        """lambda_2 = 7/5760 approximately 0.001215..."""
        from lib.genus2_obstruction_engine import lambda_fp
        val = float(lambda_fp(2))
        assert abs(val - 7/5760) < 1e-12


# ===================================================================
# 11. OBSTRUCTION TOWER CONSISTENCY
# ===================================================================

class TestObstructionTower:
    """Test obstruction tower consistency."""

    def test_tower_virasoro(self):
        """Obstruction tower is consistent for Virasoro."""
        from lib.genus2_obstruction_engine import obstruction_tower_consistency
        result = obstruction_tower_consistency('virasoro', max_genus=3)
        assert result['tower_consistent']
        assert result['all_positive']

    def test_tower_heisenberg(self):
        """Obstruction tower is consistent for Heisenberg."""
        from lib.genus2_obstruction_engine import obstruction_tower_consistency
        result = obstruction_tower_consistency('heisenberg', max_genus=3)
        assert result['tower_consistent']

    def test_ratios_kappa_independent(self):
        """F_g/F_1 ratios are independent of kappa."""
        from lib.genus2_obstruction_engine import obstruction_tower_consistency
        result = obstruction_tower_consistency('virasoro', max_genus=4)
        for g, ratio in result['ratios_to_F1'].items():
            # ratio should be a rational number (no c or k dependence)
            assert isinstance(ratio, Rational), (
                f"F_{g}/F_1 = {ratio} is not a Rational")


# ===================================================================
# 12. GRAPH DECOMPOSITION OF F_2
# ===================================================================

class TestGraphDecomposition:
    """Test stable graph decomposition of F_2."""

    def test_heisenberg_only_smooth(self):
        """For Heisenberg, only the smooth graph contributes to F_2."""
        from lib.genus2_obstruction_engine import graph_decomposition_F2_heisenberg
        result = graph_decomposition_F2_heisenberg()
        assert result['match']
        # All non-smooth contributions should be zero
        contribs = result['graph_contributions']
        for name, val in contribs.items():
            if 'smooth' not in name:
                assert val == 0, f"Non-smooth graph {name} has nonzero contribution"

    def test_virasoro_has_all_graph_types(self):
        """For Virasoro, all 6 graph types can contribute."""
        from lib.genus2_obstruction_engine import graph_decomposition_F2_virasoro
        result = graph_decomposition_F2_virasoro()
        assert len(result['graph_contributions']) == 6

    def test_virasoro_universal_formula(self):
        """Virasoro F_2 follows universal formula."""
        from lib.genus2_obstruction_engine import graph_decomposition_F2_virasoro
        result = graph_decomposition_F2_virasoro()
        assert result['lambda_2_FP'] == Rational(7, 5760)
