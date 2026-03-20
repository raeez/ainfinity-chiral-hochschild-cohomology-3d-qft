"""Tests for the cross-volume Hochschild-bulk bridge.

Verifies thm:hochschild-bridge-genus0: ChirHoch*(A) = CHC*(A) at genus 0.

Test tiers:
  Tier 1 (structural): self-certifying identities (Poincare duality, Euler char)
  Tier 2 (published): known dimensions from Kac, FBZ, etc.
  Tier 3 (cross-check): Vol I dimensions = Vol II dimensions

Paper references:
  Vol I: chiral_hochschild_koszul.tex (Theorem H), concordance.tex
  Vol II: bulk_chc.tex, thm:hochschild-bridge-genus0
  Cross-volume: cross_volume_bridge.py (Vol I)
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy import Symbol, Rational, simplify, expand, S, symbols, binomial


# ===================================================================
# ChirHoch DIMENSIONS — Poincare polynomials
# ===================================================================

class TestChirHochDimensions:
    """Verify ChirHoch*(A) Poincare polynomials for standard families.

    For Koszul chiral algebras with n generators:
      ChirHoch^k(A) = binom(n, k).
      Poincare polynomial P(t) = (1 + t)^n.

    This is Theorem H of Vol I.
    """

    def test_heisenberg_poincare(self):
        """ChirHoch*(H) has Poincare polynomial (1+t).

        One generator J, so P(t) = 1 + t.
        Dimensions: [1, 1, 0, 0, ...].
        """
        from lib.hochschild_bulk_bridge import chir_hoch_dimensions
        dims = chir_hoch_dimensions("heisenberg")
        assert dims[0] == 1
        assert dims[1] == 1
        assert dims[2] == 0
        assert dims[3] == 0

    def test_heisenberg_poincare_polynomial(self):
        """P(t) = (1+t) for Heisenberg."""
        from lib.hochschild_bulk_bridge import chir_hoch_poincare
        t = Symbol('t')
        P = chir_hoch_poincare("heisenberg", t)
        assert expand(P - (1 + t)) == 0

    def test_affine_sl2_poincare(self):
        """ChirHoch*(V_k(sl_2)) has Poincare polynomial (1+t)^3.

        Three generators J^1, J^2, J^3.
        Dimensions: [1, 3, 3, 1, 0, 0, ...].
        This is the exterior algebra Lambda*(sl_2).
        """
        from lib.hochschild_bulk_bridge import chir_hoch_dimensions
        dims = chir_hoch_dimensions("affine_sl2")
        assert dims[0] == 1
        assert dims[1] == 3
        assert dims[2] == 3
        assert dims[3] == 1
        assert dims[4] == 0

    def test_affine_sl2_poincare_polynomial(self):
        """P(t) = (1+t)^3 = 1 + 3t + 3t^2 + t^3 for affine sl_2."""
        from lib.hochschild_bulk_bridge import chir_hoch_poincare
        t = Symbol('t')
        P = chir_hoch_poincare("affine_sl2", t)
        expected = expand((1 + t)**3)
        assert expand(P - expected) == 0

    def test_virasoro_poincare(self):
        """ChirHoch*(Vir_c) has Poincare polynomial (1+t).

        One generator T, so P(t) = 1 + t.
        Same Poincare polynomial as Heisenberg (one generator each).
        """
        from lib.hochschild_bulk_bridge import chir_hoch_dimensions
        dims = chir_hoch_dimensions("virasoro")
        assert dims[0] == 1
        assert dims[1] == 1
        assert dims[2] == 0

    def test_betagamma_poincare(self):
        """ChirHoch*(betagamma) has Poincare polynomial (1+t)^2.

        Two generators beta, gamma.
        Dimensions: [1, 2, 1, 0, 0, ...].
        """
        from lib.hochschild_bulk_bridge import chir_hoch_dimensions
        dims = chir_hoch_dimensions("betagamma")
        assert dims[0] == 1
        assert dims[1] == 2
        assert dims[2] == 1
        assert dims[3] == 0

    def test_betagamma_poincare_polynomial(self):
        """P(t) = (1+t)^2 = 1 + 2t + t^2 for betagamma."""
        from lib.hochschild_bulk_bridge import chir_hoch_poincare
        t = Symbol('t')
        P = chir_hoch_poincare("betagamma", t)
        expected = expand((1 + t)**2)
        assert expand(P - expected) == 0

    def test_w3_poincare(self):
        """ChirHoch*(W_3) has Poincare polynomial (1+t)^2.

        Two generators T, W. Same Poincare polynomial as betagamma.
        """
        from lib.hochschild_bulk_bridge import chir_hoch_dimensions
        dims = chir_hoch_dimensions("w3")
        assert dims[0] == 1
        assert dims[1] == 2
        assert dims[2] == 1
        assert dims[3] == 0

    def test_lattice_rank1_poincare(self):
        """ChirHoch*(V_Z) has Poincare polynomial (1+t).

        Rank-1 lattice = Heisenberg. One generator.
        """
        from lib.hochschild_bulk_bridge import chir_hoch_dimensions
        dims = chir_hoch_dimensions("lattice_1")
        assert dims == [1, 1, 0, 0, 0, 0, 0]

    def test_lattice_rank8_poincare(self):
        """ChirHoch*(V_{E_8}) has Poincare polynomial (1+t)^8.

        Eight generators. ChirHoch^4 = binom(8,4) = 70.
        """
        from lib.hochschild_bulk_bridge import chir_hoch_dimensions
        dims = chir_hoch_dimensions("lattice_8")
        assert dims[0] == 1
        assert dims[1] == 8
        assert dims[4] == 70
        assert dims[6] == 28

    def test_lattice_rank24_total_dim(self):
        """Total dim of ChirHoch*(V_Leech) = 2^24 = 16777216."""
        from lib.hochschild_bulk_bridge import total_dimension
        assert total_dimension("lattice_24") == 2**24


# ===================================================================
# BRIDGE: ChirHoch = CHC at genus 0
# ===================================================================

class TestBridge:
    """Verify thm:hochschild-bridge-genus0: ChirHoch^n(A) = CHC^n(A).

    The bridge is a quasi-isomorphism at genus 0. Both sides compute
    the same exterior algebra on the generators.
    """

    def test_bridge_heisenberg(self):
        """Bridge for Heisenberg: ChirHoch = CHC."""
        from lib.hochschild_bulk_bridge import verify_bridge
        result = verify_bridge("heisenberg")
        assert result['match'], f"Mismatch: {result['differences']}"

    def test_bridge_affine_sl2(self):
        """Bridge for affine sl_2: ChirHoch = CHC."""
        from lib.hochschild_bulk_bridge import verify_bridge
        result = verify_bridge("affine_sl2")
        assert result['match'], f"Mismatch: {result['differences']}"

    def test_bridge_virasoro(self):
        """Bridge for Virasoro: ChirHoch = CHC."""
        from lib.hochschild_bulk_bridge import verify_bridge
        result = verify_bridge("virasoro")
        assert result['match'], f"Mismatch: {result['differences']}"

    def test_bridge_betagamma(self):
        """Bridge for betagamma: ChirHoch = CHC."""
        from lib.hochschild_bulk_bridge import verify_bridge
        result = verify_bridge("betagamma")
        assert result['match'], f"Mismatch: {result['differences']}"

    def test_bridge_w3(self):
        """Bridge for W_3: ChirHoch = CHC."""
        from lib.hochschild_bulk_bridge import verify_bridge
        result = verify_bridge("w3")
        assert result['match'], f"Mismatch: {result['differences']}"

    def test_bridge_lattice_rank2(self):
        """Bridge for rank-2 lattice: ChirHoch = CHC."""
        from lib.hochschild_bulk_bridge import verify_bridge
        result = verify_bridge("lattice_2")
        assert result['match'], f"Mismatch: {result['differences']}"

    def test_bridge_all_families(self):
        """Bridge for ALL families simultaneously."""
        from lib.hochschild_bulk_bridge import verify_bridge_all
        results = verify_bridge_all()
        for name, result in results.items():
            assert result['match'], f"{name}: mismatch {result['differences']}"


# ===================================================================
# KAPPA COMPLEMENTARITY (Theorem C)
# ===================================================================

class TestKappaComplementarity:
    """Verify kappa(A) + kappa(A!) is level-independent.

    Theorem C: Q_g(A) + Q_g(A!) = H*(M_g, Z(A)).
    At the scalar level, this means kappa(A) + kappa(A!) is a constant
    depending only on the root datum, not on the level parameter.
    """

    def test_virasoro_kappa_sum(self):
        """kappa(Vir_c) + kappa(Vir_{26-c}) = 13.

        c/2 + (26-c)/2 = 13. Independent of c.
        """
        from lib.hochschild_bulk_bridge import verify_kappa_complementarity
        result = verify_kappa_complementarity("virasoro")
        assert result['match']
        assert result['sum'] == 13

    def test_heisenberg_kappa_sum(self):
        """kappa(H_k) + kappa(H_k!) = 0.

        k + (-k) = 0. Purely anti-symmetric.
        """
        from lib.hochschild_bulk_bridge import verify_kappa_complementarity
        result = verify_kappa_complementarity("heisenberg")
        assert result['match']
        assert result['sum'] == 0

    def test_affine_sl2_kappa_sum(self):
        """kappa(sl_2_k) + kappa(sl_2_{-k-4}) = 0.

        3(k+2)/4 + 3(-k-2)/4 = 0.
        """
        from lib.hochschild_bulk_bridge import verify_kappa_complementarity
        result = verify_kappa_complementarity("affine_sl2")
        assert result['match']
        assert result['sum'] == 0

    def test_betagamma_kappa_sum(self):
        """kappa(betagamma) + kappa(betagamma!) = 0."""
        from lib.hochschild_bulk_bridge import verify_kappa_complementarity
        result = verify_kappa_complementarity("betagamma")
        assert result['match']
        assert result['sum'] == 0

    def test_w3_kappa_sum(self):
        """kappa(W_3_c) + kappa(W_3_{100-c}) = 250/3.

        5c/6 + 5(100-c)/6 = 500/6 = 250/3.
        """
        from lib.hochschild_bulk_bridge import verify_kappa_complementarity
        result = verify_kappa_complementarity("w3")
        assert result['match']
        assert result['sum'] == Rational(250, 3)

    def test_lattice_kappa_sum(self):
        """kappa(V_Lambda) + kappa(V_Lambda!) = 0."""
        from lib.hochschild_bulk_bridge import verify_kappa_complementarity
        result = verify_kappa_complementarity("lattice_1")
        assert result['match']
        assert result['sum'] == 0

    def test_kappa_formula_c_over_2(self):
        """For Virasoro: kappa = c/2, so kappa(Vir_1) = 1/2."""
        from lib.hochschild_bulk_bridge import kappa_value
        c = Symbol('c')
        kap = kappa_value("virasoro")
        assert simplify(kap - c / 2) == 0

    def test_kappa_sl2_formula(self):
        """kappa(sl_2_k) = 3(k+2)/4."""
        from lib.hochschild_bulk_bridge import kappa_value
        k = Symbol('k')
        kap = kappa_value("affine_sl2")
        assert simplify(kap - Rational(3, 4) * (k + 2)) == 0

    def test_kappa_complementarity_all(self):
        """Kappa complementarity holds for all families."""
        from lib.hochschild_bulk_bridge import verify_kappa_complementarity, ALL_FAMILIES
        for name in ALL_FAMILIES:
            result = verify_kappa_complementarity(name)
            assert result['is_level_independent'], \
                f"{name}: kappa sum depends on level: {result['sum']}"


# ===================================================================
# SHADOW ARCHETYPE
# ===================================================================

class TestShadowArchetype:
    """Verify shadow depth classification matches Vol I.

    Shadow depth r_max classifies COMPLEXITY of Koszul algebras:
    G (Gaussian): r_max = 2, L (Lie/tree): r_max = 3,
    C (Contact/quartic): r_max = 4, M (Mixed): r_max = infinity.
    """

    def test_heisenberg_gaussian(self):
        """Heisenberg is Gaussian class (G), depth 2."""
        from lib.hochschild_bulk_bridge import verify_shadow_archetype
        assert verify_shadow_archetype("heisenberg", "G", 2)

    def test_affine_lie_tree(self):
        """Affine sl_2 is Lie/tree class (L), depth 3."""
        from lib.hochschild_bulk_bridge import verify_shadow_archetype
        assert verify_shadow_archetype("affine_sl2", "L", 3)

    def test_betagamma_contact_quartic(self):
        """betagamma is Contact/quartic class (C), depth 4."""
        from lib.hochschild_bulk_bridge import verify_shadow_archetype
        assert verify_shadow_archetype("betagamma", "C", 4)

    def test_virasoro_mixed(self):
        """Virasoro is Mixed class (M), depth infinity."""
        from lib.hochschild_bulk_bridge import verify_shadow_archetype
        assert verify_shadow_archetype("virasoro", "M", float('inf'))

    def test_w3_mixed(self):
        """W_3 is Mixed class (M), depth infinity."""
        from lib.hochschild_bulk_bridge import verify_shadow_archetype
        assert verify_shadow_archetype("w3", "M", float('inf'))

    def test_lattice_gaussian(self):
        """Lattice VOAs are Gaussian class (G), depth 2."""
        from lib.hochschild_bulk_bridge import verify_shadow_archetype
        assert verify_shadow_archetype("lattice_1", "G", 2)
        assert verify_shadow_archetype("lattice_8", "G", 2)
        assert verify_shadow_archetype("lattice_24", "G", 2)

    def test_shadow_depth_does_not_determine_koszulness(self):
        """Both finite and infinite depth algebras are Koszul.

        Heisenberg (depth 2), affine (depth 3), betagamma (depth 4),
        Virasoro (depth inf) are ALL Koszul.
        Shadow depth classifies complexity, not Koszulness.
        """
        from lib.hochschild_bulk_bridge import shadow_archetype
        finite_depth_families = ["heisenberg", "affine_sl2", "betagamma"]
        infinite_depth_families = ["virasoro", "w3"]
        for fam in finite_depth_families:
            arch = shadow_archetype(fam)
            assert arch['shadow_depth'] < float('inf')
        for fam in infinite_depth_families:
            arch = shadow_archetype(fam)
            assert arch['shadow_depth'] == float('inf')


# ===================================================================
# KOSZUL DUAL PAIRING / POINCARE DUALITY
# ===================================================================

class TestKoszulDualPairing:
    """Structural properties of ChirHoch*(A) as an exterior algebra.

    For n generators, ChirHoch*(A) = Lambda*(k^n) has:
    - Euler characteristic 0 (for n >= 1)
    - Poincare duality: dim ChirHoch^k = dim ChirHoch^{n-k}
    - Total dimension 2^n
    """

    def test_euler_char_zero_heisenberg(self):
        """chi(ChirHoch*(H)) = 0."""
        from lib.hochschild_bulk_bridge import chir_hoch_euler_char
        assert chir_hoch_euler_char("heisenberg") == 0

    def test_euler_char_zero_affine(self):
        """chi(ChirHoch*(sl_2)) = 0."""
        from lib.hochschild_bulk_bridge import chir_hoch_euler_char
        assert chir_hoch_euler_char("affine_sl2") == 0

    def test_euler_char_zero_virasoro(self):
        """chi(ChirHoch*(Vir)) = 0."""
        from lib.hochschild_bulk_bridge import chir_hoch_euler_char
        assert chir_hoch_euler_char("virasoro") == 0

    def test_euler_char_zero_all(self):
        """chi(ChirHoch*(A)) = 0 for all nontrivial families."""
        from lib.hochschild_bulk_bridge import chir_hoch_euler_char, ALL_FAMILIES
        for name in ALL_FAMILIES:
            assert chir_hoch_euler_char(name) == 0, \
                f"{name}: Euler char nonzero"

    def test_poincare_duality_heisenberg(self):
        """dim ChirHoch^0(H) = dim ChirHoch^1(H) = 1."""
        from lib.hochschild_bulk_bridge import poincare_duality_check
        result = poincare_duality_check("heisenberg")
        assert result['palindromic']
        assert result['dimensions'] == [1, 1]

    def test_poincare_duality_affine(self):
        """ChirHoch*(sl_2) = [1, 3, 3, 1] is palindromic."""
        from lib.hochschild_bulk_bridge import poincare_duality_check
        result = poincare_duality_check("affine_sl2")
        assert result['palindromic']
        assert result['dimensions'] == [1, 3, 3, 1]

    def test_poincare_duality_betagamma(self):
        """ChirHoch*(betagamma) = [1, 2, 1] is palindromic."""
        from lib.hochschild_bulk_bridge import poincare_duality_check
        result = poincare_duality_check("betagamma")
        assert result['palindromic']
        assert result['dimensions'] == [1, 2, 1]

    def test_poincare_duality_all(self):
        """Poincare duality for all families."""
        from lib.hochschild_bulk_bridge import poincare_duality_check, ALL_FAMILIES
        for name in ALL_FAMILIES:
            result = poincare_duality_check(name)
            assert result['palindromic'], \
                f"{name}: Poincare duality fails {result['dimensions']}"

    def test_total_dim_heisenberg(self):
        """Total dim ChirHoch*(H) = 2^1 = 2."""
        from lib.hochschild_bulk_bridge import total_dimension
        assert total_dimension("heisenberg") == 2

    def test_total_dim_affine(self):
        """Total dim ChirHoch*(sl_2) = 2^3 = 8."""
        from lib.hochschild_bulk_bridge import total_dimension
        assert total_dimension("affine_sl2") == 8

    def test_total_dim_betagamma(self):
        """Total dim ChirHoch*(betagamma) = 2^2 = 4."""
        from lib.hochschild_bulk_bridge import total_dimension
        assert total_dimension("betagamma") == 4


# ===================================================================
# VIRASORO SELF-DUALITY (critical pitfall check)
# ===================================================================

class TestViraseroSelfDuality:
    """Verify Virasoro self-dual at c=13, NOT c=26.

    Vir_c^! = Vir_{26-c}. Self-dual iff c = 13.
    This is a critical pitfall from CLAUDE.md.
    """

    def test_self_dual_at_13(self):
        """Virasoro self-dual at c = 13."""
        from lib.hochschild_bulk_bridge import virasoro_self_dual_point
        data = virasoro_self_dual_point()
        assert data['self_dual_c'] == 13

    def test_not_self_dual_at_26(self):
        """Virasoro is NOT self-dual at c = 26."""
        from lib.hochschild_bulk_bridge import virasoro_self_dual_point
        data = virasoro_self_dual_point()
        assert data['self_dual_c'] != 26
        assert data['NOT_26'] is True

    def test_kappa_at_self_dual_point(self):
        """kappa(Vir_13) = 13/2."""
        from lib.hochschild_bulk_bridge import virasoro_self_dual_point
        data = virasoro_self_dual_point()
        assert data['kappa_at_self_dual'] == Rational(13, 2)


# ===================================================================
# MASTER BRIDGE: all verifications combined
# ===================================================================

class TestMasterBridge:
    """Run the complete bridge verification suite."""

    def test_run_all_bridges(self):
        """All bridges pass for all families."""
        from lib.hochschild_bulk_bridge import run_all_bridges
        results = run_all_bridges()
        for name, fam_data in results.items():
            assert fam_data['bridge']['match'], \
                f"{name}: bridge mismatch"
            assert fam_data['kappa']['match'], \
                f"{name}: kappa complementarity mismatch"
            assert fam_data['poincare_duality']['palindromic'], \
                f"{name}: Poincare duality fails"
            assert fam_data['euler_char'] == 0, \
                f"{name}: Euler char nonzero"
