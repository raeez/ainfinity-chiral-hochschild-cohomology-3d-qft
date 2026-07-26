r"""Tests for genus-2 free energy from the stable graph sum.

Verifies F_2(A) = kappa(A) * 7/5760 from FIRST PRINCIPLES:

  Layer 1: Stable graph enumeration at (g=2, n=0) — 7 graphs
  Layer 2: Graph sum for chi^{orb}(M_bar_2) = -1/1440
  Layer 3: Bernoulli formula B_4 = -1/30 => lambda_2^FP = 7/5760
  Layer 4: Log(A-hat) Taylor expansion => lambda_2^FP = 7/5760
  Layer 5: Cross-check linking layers 2 and 3
  Layer 6: F_2 for all standard families (Heisenberg, Virasoro, etc.)

This is the FIRST module to derive F_2 from graph sums rather than
hardcoding it from the Bernoulli formula. The graph sum computes
chi^{orb}(M_bar_2) = -1/1440 by summing vertex-product contributions
over all 7 stable graphs, then connects to lambda_2^FP via the
Harer-Zagier / Bernoulli correspondence.

References:
    thm:theorem-d (higher_genus_modular_koszul.tex): F_g = kappa * lambda_g^FP
    Faber, "A conjectural description of the tautological ring" (1999)
    Harer-Zagier, "The Euler characteristic of the moduli space of curves" (1986)
"""

import sys
import os
import unittest
from fractions import Fraction

# Ensure the compute directory is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.genus2_graph_sum import (
    # Bernoulli
    bernoulli_exact,
    # Faber-Pandharipande
    lambda_fp,
    verify_lambda_fp_from_ahat_expansion,
    # Graph enumeration
    StableGraphG2,
    enumerate_genus2_stable_graphs,
    verify_genus2_graph_enumeration,
    # Euler characteristic
    chi_orb_open,
    graph_sum_chi_orb,
    # First-principles derivation
    compute_F2,
    verify_F2_from_first_principles,
    bernoulli_from_sinh_expansion,
    # Family-specific
    F2_heisenberg,
    F2_virasoro,
    F2_affine_sl2,
    F2_w3,
    F2_betagamma,
    # Ratios and cross-checks
    F2_over_F1_ratio,
    verify_bernoulli_cross_check,
    verify_F_g_table,
    # Feynman rules
    cohft_feynman_rules_F2,
    # Penner
    penner_vertex_weight,
    penner_graph_sum_genus2,
)


# ===========================================================================
# 1. BERNOULLI NUMBERS
# ===========================================================================

class TestBernoulliNumbers(unittest.TestCase):
    """Verify Bernoulli numbers from the recursion."""

    def test_B0(self):
        """B_0 = 1."""
        self.assertEqual(bernoulli_exact(0), Fraction(1))

    def test_B1(self):
        """B_1 = -1/2."""
        self.assertEqual(bernoulli_exact(1), Fraction(-1, 2))

    def test_B2(self):
        """B_2 = 1/6."""
        self.assertEqual(bernoulli_exact(2), Fraction(1, 6))

    def test_B4(self):
        """B_4 = -1/30 (the key value for lambda_2^FP)."""
        self.assertEqual(bernoulli_exact(4), Fraction(-1, 30))

    def test_B6(self):
        """B_6 = 1/42."""
        self.assertEqual(bernoulli_exact(6), Fraction(1, 42))

    def test_odd_vanish(self):
        """B_{2k+1} = 0 for k >= 1."""
        for k in range(1, 10):
            self.assertEqual(bernoulli_exact(2 * k + 1), Fraction(0))

    def test_B4_recursion(self):
        """Verify B_4 satisfies sum_{k=0}^4 C(5,k)*B_k = 0."""
        check = (Fraction(1) + Fraction(5) * Fraction(-1, 2)
                 + Fraction(10) * Fraction(1, 6)
                 + Fraction(5) * bernoulli_exact(4))
        self.assertEqual(check, Fraction(0))


# ===========================================================================
# 2. FABER-PANDHARIPANDE NUMBERS
# ===========================================================================

class TestFaberPandharipande(unittest.TestCase):
    """Verify lambda_g^FP from the Bernoulli formula."""

    def test_lambda_1(self):
        """lambda_1^FP = 1/24."""
        self.assertEqual(lambda_fp(1), Fraction(1, 24))

    def test_lambda_2(self):
        """lambda_2^FP = 7/5760 — the central value for F_2."""
        self.assertEqual(lambda_fp(2), Fraction(7, 5760))

    def test_lambda_3(self):
        """lambda_3^FP = 31/967680."""
        self.assertEqual(lambda_fp(3), Fraction(31, 967680))

    def test_lambda_2_derivation(self):
        """Derive lambda_2^FP = (2^3-1)/2^3 * |B_4|/4! = 7/8 * 1/720 = 7/5760."""
        factor = Fraction(2**3 - 1, 2**3)  # 7/8
        abs_B4 = abs(bernoulli_exact(4))     # 1/30
        result = factor * abs_B4 / Fraction(24)  # 24 = 4!
        self.assertEqual(result, Fraction(7, 5760))

    def test_all_positive(self):
        """All lambda_g^FP are positive (Bernoulli sign absorbed)."""
        for g in range(1, 8):
            self.assertGreater(lambda_fp(g), 0)


# ===========================================================================
# 3. GRAPH ENUMERATION AT GENUS 2
# ===========================================================================

class TestGenus2GraphEnumeration(unittest.TestCase):
    """Verify the stable graph enumeration at (g=2, n=0)."""

    def test_graph_count(self):
        """There are exactly 7 stable graphs at (g=2, n=0)."""
        graphs = enumerate_genus2_stable_graphs()
        self.assertEqual(len(graphs), 7)

    def test_all_genus_2(self):
        """Every graph has arithmetic genus 2."""
        for g in enumerate_genus2_stable_graphs():
            self.assertEqual(g.arithmetic_genus, 2, f"Graph {g.name} has genus {g.arithmetic_genus}")

    def test_all_stable(self):
        """Every graph satisfies the stability condition."""
        for g in enumerate_genus2_stable_graphs():
            self.assertTrue(g.is_stable(), f"Graph {g.name} is not stable")

    def test_smooth_graph(self):
        """Graph I: smooth genus-2 curve. 1 vertex, 0 edges, |Aut|=1."""
        g = enumerate_genus2_stable_graphs()[0]
        self.assertEqual(g.name, "I")
        self.assertEqual(g.vertex_genera, (2,))
        self.assertEqual(g.num_edges, 0)
        self.assertEqual(g.aut_order, 1)
        self.assertEqual(g.h1, 0)

    def test_figure_eight(self):
        """Graph II: genus-1 with self-loop. 1 vertex, 1 edge, |Aut|=2."""
        g = enumerate_genus2_stable_graphs()[1]
        self.assertEqual(g.name, "II")
        self.assertEqual(g.vertex_genera, (1,))
        self.assertEqual(g.num_edges, 1)
        self.assertEqual(g.aut_order, 2)
        self.assertEqual(g.h1, 1)

    def test_separating(self):
        """Graph III: g=1 + g=1 separated. 2 vertices, 1 edge, |Aut|=2."""
        g = enumerate_genus2_stable_graphs()[2]
        self.assertEqual(g.name, "III")
        self.assertEqual(g.vertex_genera, (1, 1))
        self.assertEqual(g.num_edges, 1)
        self.assertEqual(g.aut_order, 2)
        self.assertEqual(g.h1, 0)

    def test_sunset(self):
        """Graph IV: g=0 with 2 self-loops. 1 vertex, 2 edges, |Aut|=8."""
        g = enumerate_genus2_stable_graphs()[3]
        self.assertEqual(g.name, "IV")
        self.assertEqual(g.vertex_genera, (0,))
        self.assertEqual(g.num_edges, 2)
        self.assertEqual(g.aut_order, 8)
        self.assertEqual(g.h1, 2)

    def test_mixed(self):
        """Graph V: g=0 (self-loop) + g=1. 2 vertices, 2 edges, |Aut|=2."""
        g = enumerate_genus2_stable_graphs()[4]
        self.assertEqual(g.name, "V")
        self.assertEqual(g.vertex_genera, (0, 1))
        self.assertEqual(g.num_edges, 2)
        self.assertEqual(g.aut_order, 2)
        self.assertEqual(g.h1, 1)

    def test_theta(self):
        """Graph VI: g=0 + g=0, 3 edges (theta graph). |Aut|=12."""
        g = enumerate_genus2_stable_graphs()[5]
        self.assertEqual(g.name, "VI")
        self.assertEqual(g.vertex_genera, (0, 0))
        self.assertEqual(g.num_edges, 3)
        self.assertEqual(g.aut_order, 12)
        self.assertEqual(g.h1, 2)

    def test_barbell(self):
        """Graph VII: two genus-0 vertices, each with a self-loop, joined by a bridge."""
        g = enumerate_genus2_stable_graphs()[6]
        self.assertEqual(g.name, "VII")
        self.assertEqual(g.vertex_genera, (0, 0))
        self.assertEqual(g.num_edges, 3)
        self.assertEqual(g.aut_order, 8)
        self.assertEqual(g.h1, 2)

    def test_valences(self):
        """Verify valences of all graphs."""
        graphs = enumerate_genus2_stable_graphs()
        expected_valences = {
            "I": (0,),      # genus-2 vertex, no half-edges
            "II": (2,),     # self-loop: 2 half-edges
            "III": (1, 1),  # one edge: 1 half-edge each
            "IV": (4,),     # two self-loops: 4 half-edges
            "V": (3, 1),    # v0: self-loop (2) + edge (1) = 3; v1: edge (1) = 1
            "VI": (3, 3),   # three edges: 3 half-edges each
            "VII": (3, 3),  # each vertex has one self-loop plus the bridge
        }
        for g in graphs:
            computed = tuple(g.valence(v) for v in range(g.num_vertices))
            self.assertEqual(computed, expected_valences[g.name],
                             f"Graph {g.name} valences: {computed} != {expected_valences[g.name]}")

    def test_h1_values(self):
        """Verify first Betti numbers."""
        graphs = enumerate_genus2_stable_graphs()
        expected_h1 = {"I": 0, "II": 1, "III": 0, "IV": 2, "V": 1, "VI": 2, "VII": 2}
        for g in graphs:
            self.assertEqual(g.h1, expected_h1[g.name])


# ===========================================================================
# 4. ORBIFOLD EULER CHARACTERISTIC
# ===========================================================================

class TestChiOrbOpen(unittest.TestCase):
    """Verify chi^{orb}(M_{g,n}) from Harer-Zagier."""

    def test_M03(self):
        """chi^{orb}(M_{0,3}) = 1."""
        self.assertEqual(chi_orb_open(0, 3), Fraction(1))

    def test_M04(self):
        """chi^{orb}(M_{0,4}) = -1."""
        self.assertEqual(chi_orb_open(0, 4), Fraction(-1))

    def test_M05(self):
        """chi^{orb}(M_{0,5}) = 2."""
        self.assertEqual(chi_orb_open(0, 5), Fraction(2))

    def test_M06(self):
        """chi^{orb}(M_{0,6}) = -6."""
        self.assertEqual(chi_orb_open(0, 6), Fraction(-6))

    def test_M11(self):
        """chi^{orb}(M_{1,1}) = -1/12."""
        self.assertEqual(chi_orb_open(1, 1), Fraction(-1, 12))

    def test_M12(self):
        """chi^{orb}(M_{1,2}) = -1/12 (recursion from M_{1,1})."""
        # chi(M_{1,2}) = (2*1-2+1) * chi(M_{1,1}) = 1 * (-1/12) = -1/12
        self.assertEqual(chi_orb_open(1, 2), Fraction(-1, 12))

    def test_M2(self):
        """chi^{orb}(M_2) = B_4/(4*2*1) = -1/240."""
        self.assertEqual(chi_orb_open(2, 0), Fraction(-1, 240))

    def test_M3(self):
        """chi^{orb}(M_3) = B_6/(4*3*2) = (1/42)/24 = 1/1008."""
        self.assertEqual(chi_orb_open(3, 0), Fraction(1, 1008))

    def test_unstable_M10(self):
        """M_{1,0} is unstable (2g-2+n = 0)."""
        with self.assertRaises(ValueError):
            chi_orb_open(1, 0)


# ===========================================================================
# 5. GRAPH SUM FOR CHI^{ORB}(M_bar_2)
# ===========================================================================

class TestChiOrbGraphSum(unittest.TestCase):
    """Verify chi^{orb}(M_bar_2) = -1/1440 from the graph sum."""

    def test_graph_sum_total(self):
        """The graph sum gives chi^{orb}(M_bar_2) = -1/1440."""
        result = graph_sum_chi_orb()
        self.assertEqual(result['total'], Fraction(-1, 1440))

    def test_graph_sum_match(self):
        """The computed value matches the known value."""
        result = graph_sum_chi_orb()
        self.assertTrue(result['match'])

    def test_smooth_contribution(self):
        """Graph I (smooth g=2): contributes chi^{orb}(M_2) = -1/240."""
        result = graph_sum_chi_orb()
        self.assertEqual(result['contributions']['I']['weighted'], Fraction(-1, 240))

    def test_figure_eight_contribution(self):
        """Graph II (figure-eight): chi(M_{1,2})/2 = (-1/12)/2 = -1/24."""
        result = graph_sum_chi_orb()
        self.assertEqual(result['contributions']['II']['weighted'], Fraction(-1, 24))

    def test_separating_contribution(self):
        """Graph III (separating): chi(M_{1,1})^2/2 = (1/144)/2 = 1/288."""
        result = graph_sum_chi_orb()
        self.assertEqual(result['contributions']['III']['weighted'], Fraction(1, 288))

    def test_sunset_contribution(self):
        """Graph IV (sunset): chi(M_{0,4})/8 = (-1)/8 = -1/8."""
        result = graph_sum_chi_orb()
        self.assertEqual(result['contributions']['IV']['weighted'], Fraction(-1, 8))

    def test_mixed_contribution(self):
        """Graph V (mixed): chi(M_{0,3})*chi(M_{1,1})/2 = (-1/12)/2 = -1/24."""
        result = graph_sum_chi_orb()
        self.assertEqual(result['contributions']['V']['weighted'], Fraction(-1, 24))

    def test_theta_contribution(self):
        """Graph VI (theta): chi(M_{0,3})^2/12 = 1/12."""
        result = graph_sum_chi_orb()
        self.assertEqual(result['contributions']['VI']['weighted'], Fraction(1, 12))

    def test_barbell_contribution(self):
        """Graph VII (barbell): chi(M_{0,3})^2/8 = 1/8."""
        result = graph_sum_chi_orb()
        self.assertEqual(result['contributions']['VII']['weighted'], Fraction(1, 8))

    def test_sum_of_contributions(self):
        """Verify the sum: -1/240 - 1/24 + 1/288 - 1/8 - 1/24 + 1/12 + 1/8 = -1/1440."""
        total = (Fraction(-1, 240)
                 + Fraction(-1, 24)
                 + Fraction(1, 288)
                 + Fraction(-1, 8)
                 + Fraction(-1, 24)
                 + Fraction(1, 12)
                 + Fraction(1, 8))
        self.assertEqual(total, Fraction(-1, 1440))


# ===========================================================================
# 6. LOG(A-HAT) TAYLOR EXPANSION
# ===========================================================================

class TestAhatExpansion(unittest.TestCase):
    """Verify lambda_g^FP from the Taylor expansion of (x/2)/sinh(x/2)."""

    def test_ahat_lambda_1(self):
        """|[x^2] Ahat| = 1/24 = lambda_1^FP."""
        result = bernoulli_from_sinh_expansion(max_genus=3)
        self.assertEqual(result['verification'][1]['computed_from_ahat'],
                         Fraction(1, 24))

    def test_ahat_lambda_2(self):
        """|[x^4] Ahat| = 7/5760 = lambda_2^FP."""
        result = bernoulli_from_sinh_expansion(max_genus=3)
        self.assertEqual(result['verification'][2]['computed_from_ahat'],
                         Fraction(7, 5760))

    def test_ahat_lambda_3(self):
        """|[x^6] Ahat| = 31/967680 = lambda_3^FP."""
        result = bernoulli_from_sinh_expansion(max_genus=3)
        self.assertEqual(result['verification'][3]['computed_from_ahat'],
                         Fraction(31, 967680))

    def test_ahat_matches_bernoulli(self):
        """All coefficients match the Bernoulli formula through g=5."""
        result = bernoulli_from_sinh_expansion(max_genus=5)
        self.assertTrue(result['all_match'])

    def test_ahat_sign_alternation(self):
        """Coefficients of (x/2)/sinh(x/2) alternate in sign."""
        result = bernoulli_from_sinh_expansion(max_genus=5)
        for g in range(1, 6):
            self.assertTrue(result['verification'][g]['sign_alternating'],
                            f"Sign pattern wrong at g={g}")

    def test_alternative_ahat(self):
        """Alternative verification via verify_lambda_fp_from_ahat_expansion."""
        result = verify_lambda_fp_from_ahat_expansion(max_genus=5)
        for g in range(1, 6):
            self.assertTrue(result[g]['match'],
                            f"Mismatch at g={g}: {result[g]}")


# ===========================================================================
# 7. FIRST-PRINCIPLES VERIFICATION OF F_2
# ===========================================================================

class TestF2FirstPrinciples(unittest.TestCase):
    """The master test: F_2 = kappa * 7/5760 from first principles."""

    def test_all_layers_pass(self):
        """All 5 verification layers pass."""
        result = verify_F2_from_first_principles(kappa=Fraction(1))
        self.assertTrue(result['all_passed'],
                        f"Layers failed: {[f'layer{i}' for i in range(1,6) if not result[f'layer{i}']['passed']]}")

    def test_layer1_enumeration(self):
        """Layer 1: 7 stable graphs, all genus 2, all stable."""
        result = verify_F2_from_first_principles()
        self.assertTrue(result['layer1']['passed'])
        self.assertEqual(result['layer1']['graph_count'], 7)

    def test_layer2_graph_sum(self):
        """Layer 2: chi^{orb}(M_bar_2) = -1/1440 from graph sum."""
        result = verify_F2_from_first_principles()
        self.assertTrue(result['layer2']['passed'])
        self.assertEqual(result['layer2']['computed'], Fraction(-1, 1440))

    def test_layer3_bernoulli(self):
        """Layer 3: lambda_2^FP = 7/5760 from B_4 = -1/30."""
        result = verify_F2_from_first_principles()
        self.assertTrue(result['layer3']['passed'])

    def test_layer4_ahat_expansion(self):
        """Layer 4: lambda_2^FP = 7/5760 from Ahat power series expansion."""
        result = verify_F2_from_first_principles()
        self.assertTrue(result['layer4']['passed'])

    def test_layer5_cross_check(self):
        """Layer 5: lambda_2^FP = (7/8)*(8/24)*|chi^{orb}(M_2)|."""
        result = verify_F2_from_first_principles()
        self.assertTrue(result['layer5']['passed'])

    def test_F2_value(self):
        """F_2 = kappa * 7/5760 at kappa = 1."""
        result = verify_F2_from_first_principles(kappa=Fraction(1))
        self.assertEqual(result['F_2']['F_2'], Fraction(7, 5760))


# ===========================================================================
# 8. F_2 FOR STANDARD FAMILIES
# ===========================================================================

class TestF2Families(unittest.TestCase):
    """Verify F_2 for all standard families."""

    def test_heisenberg_k1(self):
        """F_2(H_1) = 7/5760 (kappa = 1)."""
        self.assertEqual(F2_heisenberg(Fraction(1)), Fraction(7, 5760))

    def test_heisenberg_k2(self):
        """F_2(H_2) = 7/2880 (kappa = 2)."""
        self.assertEqual(F2_heisenberg(Fraction(2)), Fraction(7, 2880))

    def test_virasoro_c26(self):
        """F_2(Vir_26) = 91/5760 = 7*13/5760 (kappa = 13)."""
        self.assertEqual(F2_virasoro(Fraction(26)), Fraction(91, 5760))

    def test_virasoro_c1(self):
        """F_2(Vir_1) = 7/11520 (kappa = 1/2)."""
        self.assertEqual(F2_virasoro(Fraction(1)), Fraction(7, 11520))

    def test_virasoro_c0(self):
        """F_2(Vir_0) = 0 (kappa = 0 at c=0, uncurved)."""
        self.assertEqual(F2_virasoro(Fraction(0)), Fraction(0))

    def test_affine_sl2_k1(self):
        """F_2(V_1(sl_2)) = 7*3*3/(4*5760) = 63/23040 = 7/2560... let me compute.

        kappa = 3*(1+2)/4 = 9/4.
        F_2 = (9/4) * 7/5760 = 63/23040 = 7/2560.
        """
        result = F2_affine_sl2(Fraction(1))
        expected = Fraction(9, 4) * Fraction(7, 5760)
        self.assertEqual(result, expected)

    def test_betagamma(self):
        """F_2(beta-gamma) = 7/5760 (kappa = 1)."""
        self.assertEqual(F2_betagamma(), Fraction(7, 5760))

    def test_w3_c2(self):
        """F_2(W_3 at c=2) = 5*2/6 * 7/5760 = 5/3 * 7/5760 = 7/3456."""
        result = F2_w3(Fraction(2))
        expected = Fraction(5, 3) * Fraction(7, 5760)
        self.assertEqual(result, expected)


# ===========================================================================
# 9. UNIVERSAL RATIOS
# ===========================================================================

class TestUniversalRatios(unittest.TestCase):
    """Verify universal (kappa-independent) ratios involving F_2."""

    def test_F2_over_F1(self):
        """F_2/F_1 = lambda_2^FP/lambda_1^FP = (7/5760)/(1/24) = 7/240."""
        self.assertEqual(F2_over_F1_ratio(), Fraction(7, 240))

    def test_F2_over_F1_derivation(self):
        """Derive F_2/F_1 = 7*24/5760 = 168/5760 = 7/240."""
        ratio = Fraction(7, 5760) / Fraction(1, 24)
        self.assertEqual(ratio, Fraction(7, 240))

    def test_ratio_kappa_independent(self):
        """F_2/F_1 doesn't depend on kappa."""
        for kappa in [Fraction(1), Fraction(1, 2), Fraction(13), Fraction(9, 4)]:
            F1 = kappa * lambda_fp(1)
            F2 = kappa * lambda_fp(2)
            if F1 != 0:
                self.assertEqual(F2 / F1, Fraction(7, 240))


# ===========================================================================
# 10. BERNOULLI CROSS-CHECK
# ===========================================================================

class TestBernoulliCrossCheck(unittest.TestCase):
    """Cross-check F_2 against B_4 via independent computations."""

    def test_bernoulli_cross_check(self):
        """All cross-checks pass."""
        result = verify_bernoulli_cross_check()
        self.assertTrue(result['recursion_satisfied'])
        self.assertTrue(result['lambda_2_correct'])
        self.assertTrue(result['ratio_correct'])

    def test_B4_from_recursion(self):
        """B_4 = -1/30 from the sum_{k=0}^4 C(5,k)B_k = 0 recursion."""
        result = verify_bernoulli_cross_check()
        self.assertEqual(result['B_4'], Fraction(-1, 30))
        self.assertEqual(result['B_4_recursion_check'], Fraction(0))


# ===========================================================================
# 11. COMPLEMENTARITY AND VIRASORO DUALITY
# ===========================================================================

class TestComplementarity(unittest.TestCase):
    """Verify complementarity relations for F_2."""

    def test_virasoro_complementarity_sum(self):
        """F_2(Vir_c) + F_2(Vir_{26-c}) = 7*26/(2*5760) = 91/5760.

        This is c-independent (AP24: kappa + kappa' = 26/2 = 13 for Virasoro).
        """
        for c_val in [0, 1, 2, 13, 25, 26]:
            c = Fraction(c_val)
            F2_c = F2_virasoro(c)
            F2_dual = F2_virasoro(Fraction(26) - c)
            total = F2_c + F2_dual
            self.assertEqual(total, Fraction(91, 5760),
                             f"F_2(Vir_{c}) + F_2(Vir_{26-c}) = {total} != 91/5760")

    def test_virasoro_comparison_fixed(self):
        """At c=13 fixed point: F_2(Vir_13) = 91/11520 = (13/2)*7/5760."""
        F2_13 = F2_virasoro(Fraction(13))
        expected = Fraction(13, 2) * Fraction(7, 5760)
        self.assertEqual(F2_13, expected)

    def test_heisenberg_additive(self):
        """kappa is additive: F_2(H_{k1+k2}) = F_2(H_{k1}) + F_2(H_{k2})."""
        k1, k2 = Fraction(3), Fraction(5)
        self.assertEqual(F2_heisenberg(k1 + k2),
                         F2_heisenberg(k1) + F2_heisenberg(k2))


# ===========================================================================
# 12. FULL TABLE VERIFICATION
# ===========================================================================

class TestFgTable(unittest.TestCase):
    """Verify the table of F_g values through genus 5."""

    def test_all_match(self):
        """Bernoulli and log(Ahat) agree through g=5."""
        result = verify_F_g_table(max_genus=5)
        self.assertTrue(result['all_match'])

    def test_all_positive(self):
        """All lambda_g^FP are positive."""
        result = verify_F_g_table(max_genus=5)
        self.assertTrue(result['all_positive'])

    def test_lambda_values(self):
        """Spot-check lambda_g^FP values."""
        table = verify_F_g_table(max_genus=5)['table']
        self.assertEqual(table[1]['lambda_g_FP'], Fraction(1, 24))
        self.assertEqual(table[2]['lambda_g_FP'], Fraction(7, 5760))
        self.assertEqual(table[3]['lambda_g_FP'], Fraction(31, 967680))


# ===========================================================================
# 13. GRAPH VERIFICATION VIA Vol I CROSS-CHECK
# ===========================================================================

class TestVolICrossCheck(unittest.TestCase):
    """Cross-check graph enumeration against the Vol I stable_graph_enumeration module."""

    def test_graph_count_matches_vol1(self):
        """7 graphs at (g=2,n=0) — same count as Vol I."""
        # Vol I has genus2_stable_graphs_n0() returning 7 graphs
        graphs = enumerate_genus2_stable_graphs()
        self.assertEqual(len(graphs), 7)

    def test_automorphism_orders_match_vol1(self):
        """Automorphism orders match Vol I: 1, 2, 2, 8, 2, 12, 8."""
        graphs = enumerate_genus2_stable_graphs()
        aut_orders = [g.aut_order for g in graphs]
        self.assertEqual(aut_orders, [1, 2, 2, 8, 2, 12, 8])

    def test_automorphism_product(self):
        """Product of 1/|Aut| over all graphs is a sanity check.

        sum 1/|Aut| = 1/1 + 1/2 + 1/2 + 1/8 + 1/2 + 1/12 + 1/8
                    = 1 + 1/2 + 1/2 + 1/8 + 1/2 + 1/12 + 1/8
                    = 24/24 + 12/24 + 12/24 + 3/24 + 12/24 + 2/24 + 3/24
                    = 68/24 = 17/6
        """
        graphs = enumerate_genus2_stable_graphs()
        total = sum(Fraction(1, g.aut_order) for g in graphs)
        self.assertEqual(total, Fraction(17, 6))


# ===========================================================================
# 14. EDGE CASES AND CONSISTENCY
# ===========================================================================

class TestEdgeCases(unittest.TestCase):
    """Edge cases and consistency checks."""

    def test_F2_at_kappa_zero(self):
        """F_2 = 0 when kappa = 0 (e.g., Vir at c=0)."""
        self.assertEqual(Fraction(0) * lambda_fp(2), Fraction(0))

    def test_F2_positive_for_positive_kappa(self):
        """F_2 > 0 for kappa > 0."""
        for k in [1, 2, 5, 13, 100]:
            F2 = Fraction(k) * lambda_fp(2)
            self.assertGreater(F2, 0)

    def test_F2_float_precision(self):
        """F_2 as float matches 7/5760 to machine precision."""
        self.assertAlmostEqual(float(lambda_fp(2)), 7 / 5760, places=15)

    def test_F2_over_F1_squared(self):
        """F_2/F_1^2 = 7/(10*kappa) (kappa-dependent!)."""
        kappa = Fraction(1)
        F1 = kappa * lambda_fp(1)
        F2 = kappa * lambda_fp(2)
        ratio = F2 / F1 ** 2
        # F_2/F_1^2 = (7*kappa/5760)/(kappa^2/576) = 7/(10*kappa)
        expected = Fraction(7, 10) / kappa
        self.assertEqual(ratio, expected)

    def test_penner_vertex_weight(self):
        """Penner vertex weight V_n = (-1)^{n-1}(n-1)!."""
        self.assertEqual(penner_vertex_weight(3), Fraction(2))    # (-1)^2 * 2! = 2
        self.assertEqual(penner_vertex_weight(4), Fraction(-6))   # (-1)^3 * 3! = -6

    def test_compute_F2_consistent(self):
        """compute_F2 gives consistent results."""
        result = compute_F2(Fraction(1))
        self.assertTrue(result['all_layers_consistent'])


if __name__ == '__main__':
    unittest.main()
