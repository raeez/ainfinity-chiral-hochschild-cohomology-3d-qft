"""Cross-engine consistency tests: verify invariants agree across independent code paths.

Tests that the SAME mathematical invariant (kappa, dual kappa, shadow class/depth,
complementarity sum) is computed identically by independent engines. This catches
normalization bugs, formula transcription errors, and convention drift.

Ground truth: conclusion.tex:922-930 (the complementarity table).

Engines compared:
  BBD = bulk_boundary_duality_engine (KoszulDualPair dataclass)
  MPQ = modular_pva_quantization (genus0_classical_data)
  HBB = hochschild_bulk_bridge (ChiralAlgebraData dataclass)
  G1  = genus_one_bridge (family registry)

AP10 compliance: no hardcoded expected values in cross-engine tests.
Each test compares two engines against each other.

References:
  conclusion.tex:922-930: ground truth complementarity table
  rosetta_stone.tex:457,763,827: kappa(H_k) = k
  w-algebras-w3.tex:320: W_3 self-dual at c = 50
  ht_bulk_boundary_line.tex:2704: W_3(c)^! = W_3(100-c)
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from sympy import Symbol, Rational, simplify, S, oo


k = Symbol('k')
c = Symbol('c')


# =========================================================================
# Helpers: extract kappa from each engine
# =========================================================================

def _bbd_pair(family):
    from lib.bulk_boundary_duality_engine import (
        heisenberg_koszul_pair, virasoro_koszul_pair,
        affine_sl2_koszul_pair, betagamma_koszul_pair,
        w3_koszul_pair,
    )
    return {
        'heisenberg': heisenberg_koszul_pair,
        'virasoro': virasoro_koszul_pair,
        'affine_sl2': affine_sl2_koszul_pair,
        'betagamma': betagamma_koszul_pair,
        'w3': w3_koszul_pair,
    }[family]()


def _mpq_data(family, **params):
    from lib.modular_pva_quantization import genus0_classical_data
    return genus0_classical_data(family, **params)


def _hbb_data(family):
    from lib.hochschild_bulk_bridge import ALL_FAMILIES
    return ALL_FAMILIES[family]()


def _g1_kappa(g1_family, **params):
    from lib.genus_one_bridge import genus1_curvature
    return genus1_curvature(g1_family, **params)['kappa']


# =========================================================================
# 1. KAPPA AGREEMENT: BBD vs MPQ (both Vol II engines, same convention)
# =========================================================================

class TestKappaBBDvsMPQ:
    """BBD and MPQ should use identical kappa formulas."""

    def test_virasoro(self):
        """kappa(Vir_c): BBD = MPQ = c/2."""
        bbd = _bbd_pair('virasoro').kappa
        mpq = _mpq_data('virasoro', c=c)['kappa']
        assert simplify(bbd - mpq) == 0

    def test_affine_sl2(self):
        """kappa(sl_2_k): BBD = MPQ = 3(k+2)/4."""
        bbd = _bbd_pair('affine_sl2').kappa
        mpq = _mpq_data('affine_sl2', k=k)['kappa']
        assert simplify(bbd - mpq) == 0

    def test_w3(self):
        """kappa(W_3): BBD = MPQ = 5c/6."""
        bbd = _bbd_pair('w3').kappa
        mpq = _mpq_data('w3', c=c)['kappa']
        assert simplify(bbd - mpq) == 0

    def test_heisenberg(self):
        """kappa(H_k): BBD = MPQ."""
        bbd = _bbd_pair('heisenberg').kappa
        mpq = _mpq_data('heisenberg', k=k)['kappa']
        assert simplify(bbd - mpq) == 0

    @pytest.mark.xfail(reason=(
        "Vol I/Vol II sign discrepancy: BBD uses κ(βγ)=-1 (Vol II rosetta_stone.tex:1996), "
        "MPQ uses κ(βγ)=+1 (Vol I preface:3800). Unresolved convention."
    ))
    def test_betagamma(self):
        """kappa(betagamma): BBD vs MPQ — UNRESOLVED SIGN CONVENTION."""
        bbd = _bbd_pair('betagamma').kappa
        mpq = _mpq_data('betagamma')['kappa']
        assert simplify(bbd - mpq) == 0


# =========================================================================
# 2. KAPPA AGREEMENT: BBD vs HBB
# =========================================================================

class TestKappaBBDvsHBB:
    """BBD and HBB kappa comparison."""

    def test_virasoro(self):
        """kappa(Vir_c): BBD = HBB = c/2."""
        bbd = _bbd_pair('virasoro').kappa
        hbb = _hbb_data('virasoro').kappa
        assert simplify(bbd - hbb) == 0

    def test_affine_sl2(self):
        """kappa(sl_2_k): BBD = HBB = 3(k+2)/4."""
        bbd = _bbd_pair('affine_sl2').kappa
        hbb = _hbb_data('affine_sl2').kappa
        assert simplify(bbd - hbb) == 0

    def test_w3(self):
        """kappa(W_3): BBD = HBB = 5c/6."""
        bbd = _bbd_pair('w3').kappa
        hbb = _hbb_data('w3').kappa
        assert simplify(bbd - hbb) == 0

    def test_heisenberg(self):
        """kappa(H_k): BBD = HBB = k."""
        bbd = _bbd_pair('heisenberg').kappa
        hbb = _hbb_data('heisenberg').kappa
        assert simplify(bbd - hbb) == 0

    @pytest.mark.xfail(reason=(
        "Vol I/Vol II sign discrepancy: BBD uses κ(βγ)=-1 (Vol II), "
        "HBB uses κ(βγ)=+1 (Vol I). Unresolved convention."
    ))
    def test_betagamma(self):
        """kappa(betagamma): BBD vs HBB — UNRESOLVED SIGN CONVENTION."""
        bbd = _bbd_pair('betagamma').kappa
        hbb = _hbb_data('betagamma').kappa
        assert simplify(bbd - hbb) == 0


# =========================================================================
# 3. KAPPA AGREEMENT: BBD vs G1 (genus_one_bridge)
# =========================================================================

class TestKappaBBDvsG1:
    """BBD and genus_one_bridge kappa comparison."""

    def test_virasoro(self):
        """kappa(Vir_c): BBD = G1 = c/2."""
        bbd = _bbd_pair('virasoro').kappa
        g1 = _g1_kappa('virasoro', c=c)
        assert simplify(bbd - g1) == 0

    def test_affine_sl2(self):
        """kappa(sl_2_k): BBD = G1 = 3(k+2)/4."""
        bbd = _bbd_pair('affine_sl2').kappa
        g1 = _g1_kappa('nonabelian_cs', k=k)
        assert simplify(bbd - g1) == 0

    def test_w3(self):
        """kappa(W_3): BBD = G1 = 5c/6."""
        bbd = _bbd_pair('w3').kappa
        g1 = _g1_kappa('w3', c=c)
        assert simplify(bbd - g1) == 0

    def test_heisenberg(self):
        """kappa(H_k): BBD = G1 = k."""
        bbd = _bbd_pair('heisenberg').kappa
        g1 = _g1_kappa('abelian_cs', k=k)
        assert simplify(bbd - g1) == 0


# =========================================================================
# 4. SHADOW CLASS/DEPTH AGREEMENT
# =========================================================================

class TestShadowConsistency:
    """Shadow archetype and depth agree across BBD, MPQ, HBB."""

    @pytest.mark.parametrize("family", [
        'heisenberg', 'virasoro', 'affine_sl2', 'betagamma', 'w3',
    ])
    def test_shadow_class_bbd_vs_mpq(self, family):
        """BBD shadow_class = MPQ shadow_archetype for all families."""
        bbd_class = _bbd_pair(family).shadow_class
        mpq_class = _mpq_data(family, k=1, c=1)['shadow_archetype']
        assert bbd_class == mpq_class, (
            f"{family}: BBD={bbd_class}, MPQ={mpq_class}"
        )

    @pytest.mark.parametrize("family", [
        'heisenberg', 'virasoro', 'affine_sl2', 'betagamma', 'w3',
    ])
    def test_shadow_class_bbd_vs_hbb(self, family):
        """BBD shadow_class = HBB shadow_class for all families."""
        bbd_class = _bbd_pair(family).shadow_class
        hbb_class = _hbb_data(family).shadow_class
        assert bbd_class == hbb_class, (
            f"{family}: BBD={bbd_class}, HBB={hbb_class}"
        )

    @pytest.mark.parametrize("family", [
        'heisenberg', 'virasoro', 'affine_sl2', 'betagamma', 'w3',
    ])
    def test_shadow_depth_mpq_vs_hbb(self, family):
        """MPQ shadow_depth = HBB shadow_depth for all families."""
        mpq_depth = _mpq_data(family, k=1, c=1)['shadow_depth']
        hbb_depth = _hbb_data(family).shadow_depth
        # Normalize: sympy oo vs float('inf')
        if mpq_depth == oo:
            mpq_depth = float('inf')
        assert mpq_depth == hbb_depth, (
            f"{family}: MPQ={mpq_depth}, HBB={hbb_depth}"
        )


# =========================================================================
# 5. COMPLEMENTARITY SUM AGREEMENT
# =========================================================================

class TestComplementaritySumConsistency:
    """kappa(A) + kappa(A!) consistent across engines."""

    def test_virasoro_bbd_vs_hbb(self):
        """Virasoro: kappa_sum = 13 in both BBD and HBB."""
        bbd_sum = _bbd_pair('virasoro').kappa_sum
        hbb_sum = _hbb_data('virasoro').kappa_sum
        assert simplify(bbd_sum - hbb_sum) == 0

    def test_affine_sl2_bbd_vs_hbb(self):
        """sl_2: kappa_sum = 0 in both BBD and HBB."""
        bbd_sum = _bbd_pair('affine_sl2').kappa_sum
        hbb_sum = _hbb_data('affine_sl2').kappa_sum
        assert simplify(bbd_sum - hbb_sum) == 0

    def test_w3_bbd_vs_hbb(self):
        """W_3: kappa_sum = 250/3 in both BBD and HBB."""
        bbd_sum = _bbd_pair('w3').kappa_sum
        hbb_sum = _hbb_data('w3').kappa_sum
        assert simplify(bbd_sum - hbb_sum) == 0

    def test_heisenberg_sum_zero_all_engines(self):
        """Heisenberg: kappa_sum = 0 in all engines (convention-independent)."""
        bbd_sum = _bbd_pair('heisenberg').kappa_sum
        hbb_sum = _hbb_data('heisenberg').kappa_sum
        assert simplify(bbd_sum) == 0
        assert simplify(hbb_sum) == 0

    def test_betagamma_sum_zero_all_engines(self):
        """betagamma: kappa_sum = 0 in all engines."""
        bbd_sum = _bbd_pair('betagamma').kappa_sum
        hbb_sum = _hbb_data('betagamma').kappa_sum
        assert simplify(bbd_sum) == 0
        assert simplify(hbb_sum) == 0

    def test_w3_bbd_vs_g1(self):
        """W_3: BBD complementarity sum = G1 complementarity constant."""
        from lib.genus_one_bridge import complementarity_genus1
        bbd_sum = _bbd_pair('w3').kappa_sum
        g1_result = complementarity_genus1('w3', c=c)
        g1_sum = g1_result['kappa_sum']
        assert simplify(bbd_sum - g1_sum) == 0


# =========================================================================
# 6. W_3 DUAL CENTRAL CHARGE (regression test for the 98/5 bug)
# =========================================================================

class TestW3DualCentralCharge:
    """W_3 Koszul dual central charge = 100 - c (not 98/5 - c)."""

    def test_bbd_dual_c(self):
        """BBD: c(W_3^!) = 100 - c."""
        pair = _bbd_pair('w3')
        assert simplify(pair.dual_central_charge - (100 - c)) == 0

    def test_hbb_dual_c(self):
        """HBB: c(W_3^!) = 100 - c."""
        data = _hbb_data('w3')
        assert simplify(data.dual_central_charge - (100 - c)) == 0

    def test_bbd_self_dual_at_50(self):
        """BBD self-dual point = 50."""
        pair = _bbd_pair('w3')
        assert pair.self_dual_point == 50

    def test_kappa_sum_250_3(self):
        """kappa_sum = 250/3 (not 49/3)."""
        pair = _bbd_pair('w3')
        assert pair.kappa_sum == Rational(250, 3)

    def test_bbd_matches_g1(self):
        """BBD W_3 dual kappa matches G1."""
        from lib.genus_one_bridge import complementarity_genus1
        bbd_dual = _bbd_pair('w3').dual_kappa
        g1_result = complementarity_genus1('w3', c=c)
        g1_dual = g1_result['kappa_A_dual']
        assert simplify(bbd_dual - g1_dual) == 0


# =========================================================================
# 7. GROUND TRUTH: conclusion.tex complementarity table
# =========================================================================

class TestGroundTruthTable:
    """Verify engines match conclusion.tex:922-930 ground truth.

    Ground truth (conclusion.tex):
      H_k:   kappa=k,     dual=-k,     sum=0,     class=G
      sl_2:  kappa=3(k+2)/4, dual=3(-k-2)/4, sum=0, class=L
      bg:    kappa=-1,    dual=1,      sum=0,     class=C
      Vir_c: kappa=c/2,   dual=(26-c)/2, sum=13,  class=M
      W_3:   kappa=5c/6,  dual=5(100-c)/6, sum=250/3, class=M
    """

    def test_affine_sl2_kappa(self):
        """Ground truth: kappa(sl_2) = 3(k+2)/4."""
        expected = Rational(3, 4) * (k + 2)
        assert simplify(_bbd_pair('affine_sl2').kappa - expected) == 0
        assert simplify(_hbb_data('affine_sl2').kappa - expected) == 0
        assert simplify(_mpq_data('affine_sl2', k=k)['kappa'] - expected) == 0

    def test_virasoro_kappa(self):
        """Ground truth: kappa(Vir_c) = c/2."""
        expected = c / 2
        assert simplify(_bbd_pair('virasoro').kappa - expected) == 0
        assert simplify(_hbb_data('virasoro').kappa - expected) == 0
        assert simplify(_mpq_data('virasoro', c=c)['kappa'] - expected) == 0

    def test_w3_kappa(self):
        """Ground truth: kappa(W_3) = 5c/6."""
        expected = 5 * c / 6
        assert simplify(_bbd_pair('w3').kappa - expected) == 0
        assert simplify(_hbb_data('w3').kappa - expected) == 0
        assert simplify(_mpq_data('w3', c=c)['kappa'] - expected) == 0

    def test_virasoro_complementarity_sum(self):
        """Ground truth: kappa(Vir_c) + kappa(Vir_{26-c}) = 13."""
        assert _bbd_pair('virasoro').kappa_sum == 13
        assert _hbb_data('virasoro').kappa_sum == 13

    def test_w3_complementarity_sum(self):
        """Ground truth: kappa(W_3_c) + kappa(W_3_{100-c}) = 250/3."""
        assert _bbd_pair('w3').kappa_sum == Rational(250, 3)
        assert _hbb_data('w3').kappa_sum == Rational(250, 3)

    def test_heisenberg_kappa_ground_truth(self):
        """Ground truth: kappa(H_k) = k. All engines now agree."""
        expected = k
        assert simplify(_bbd_pair('heisenberg').kappa - expected) == 0
        assert simplify(_hbb_data('heisenberg').kappa - expected) == 0
        assert simplify(_mpq_data('heisenberg', k=k)['kappa'] - expected) == 0

    def test_betagamma_complementarity_sum_zero(self):
        """All engines agree: kappa(betagamma) + kappa(betagamma!) = 0."""
        assert _bbd_pair('betagamma').kappa_sum == 0
        assert _hbb_data('betagamma').kappa_sum == 0
