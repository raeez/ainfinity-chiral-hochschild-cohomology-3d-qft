"""Tests for Bulk-Boundary-Line Duality Engine (Theorem J).

Verifies:
1. Koszul dual pairs for all standard families
2. kappa complementarity (Theorem C)
3. Line operators as A!-modules (Theorem J)
4. Classical r-matrix and CYBE
5. Cross-volume bridge (Vol I shadow data = Vol II boundary)
6. Self-duality analysis (AP8 compliance)
7. Yangian shadow data
8. Comprehensive duality table
9. Conformal weight spectra
10. Anomaly completion data

Each test performs ACTUAL computation or consistency check.

References:
  Vol II: Theorem J (bulk-boundary-line chapter)
  Vol I: chiral_koszul_pairs.tex, concordance.tex
  CLAUDE.md: Critical Pitfalls AP1-AP13
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from sympy import Symbol, Rational, S, simplify, expand, oo, symbols

from lib.bulk_boundary_duality_engine import (
    # Koszul dual pairs
    KoszulDualPair,
    heisenberg_koszul_pair,
    virasoro_koszul_pair,
    affine_sl2_koszul_pair,
    betagamma_koszul_pair,
    w3_koszul_pair,
    # Complementarity
    verify_kappa_complementarity,
    # Line operators
    LineOperatorData,
    heisenberg_line_operators,
    virasoro_line_operators,
    affine_sl2_line_operators,
    # r-matrix
    classical_r_matrix,
    verify_cybe_abelian,
    verify_cybe_sl2,
    # Cross-volume
    cross_volume_kappa_bridge,
    cross_volume_shadow_class_bridge,
    # Theorem J
    verify_theorem_j_heisenberg,
    verify_theorem_j_virasoro,
    verify_theorem_j_affine,
    # Yangian
    yangian_shadow_from_r_matrix,
    # Self-duality
    self_duality_analysis,
    # Duality table
    comprehensive_duality_table,
    # Spectrum
    conformal_weight_spectrum,
    # Anomaly
    anomaly_completion_data,
)


# ===================================================================
# 1. KOSZUL DUAL PAIRS
# ===================================================================

class TestKoszulDualPairs:
    """Test Koszul dual pair construction for each family."""

    def test_heisenberg_generators(self):
        """Heisenberg: 1 generator J, dual has 1 generator phi."""
        pair = heisenberg_koszul_pair()
        assert pair.algebra_generators == ['J']
        assert pair.dual_generators == ['phi']

    def test_virasoro_generators(self):
        """Virasoro: 1 generator T."""
        pair = virasoro_koszul_pair()
        assert pair.algebra_generators == ['T']

    def test_sl2_generators(self):
        """sl_2: 3 generators e, h, f."""
        pair = affine_sl2_koszul_pair()
        assert len(pair.algebra_generators) == 3

    def test_w3_generators(self):
        """W_3: 2 generators T, W."""
        pair = w3_koszul_pair()
        assert pair.algebra_generators == ['T', 'W']

    def test_betagamma_generators(self):
        """Beta-gamma: 2 generators beta, gamma."""
        pair = betagamma_koszul_pair()
        assert pair.algebra_generators == ['beta', 'gamma']

    def test_heisenberg_not_self_dual(self):
        """Heisenberg is NOT self-dual (CRITICAL PITFALL)."""
        pair = heisenberg_koszul_pair()
        assert pair.self_dual_point is None

    def test_virasoro_self_dual_at_13(self):
        """Virasoro self-dual at c=13 (NOT c=26)."""
        pair = virasoro_koszul_pair()
        assert pair.self_dual_point == S(13)

    def test_w3_self_dual_at_50(self):
        """W_3 self-dual at c = 50."""
        pair = w3_koszul_pair()
        assert pair.self_dual_point == S(50)


# ===================================================================
# 2. KAPPA COMPLEMENTARITY (THEOREM C)
# ===================================================================

class TestKappaComplementarity:
    """Verify kappa(A) + kappa(A!) for all families (Theorem C)."""

    def test_heisenberg_complementarity(self):
        """Heisenberg: kappa + kappa' = 0."""
        pair = heisenberg_koszul_pair()
        result = verify_kappa_complementarity(pair)
        assert result['match']
        assert result['kappa_sum'] == '0'

    def test_virasoro_complementarity(self):
        """Virasoro: kappa + kappa' = 13."""
        pair = virasoro_koszul_pair()
        result = verify_kappa_complementarity(pair)
        assert result['match']

    def test_sl2_complementarity(self):
        """sl_2: kappa + kappa' = 0."""
        pair = affine_sl2_koszul_pair()
        result = verify_kappa_complementarity(pair)
        assert result['match']

    def test_betagamma_complementarity(self):
        """Beta-gamma: kappa + kappa' = 0."""
        pair = betagamma_koszul_pair()
        result = verify_kappa_complementarity(pair)
        assert result['match']

    def test_w3_complementarity(self):
        """W_3: kappa + kappa' = 250/3."""
        pair = w3_koszul_pair()
        result = verify_kappa_complementarity(pair)
        assert result['match']

    def test_heisenberg_level_independent(self):
        """Heisenberg complementarity constant is level-independent."""
        pair = heisenberg_koszul_pair()
        result = verify_kappa_complementarity(pair)
        assert result['level_independent']

    def test_virasoro_level_independent(self):
        """Virasoro complementarity constant is c-independent."""
        pair = virasoro_koszul_pair()
        result = verify_kappa_complementarity(pair)
        assert result['level_independent']

    def test_sl2_level_independent(self):
        """sl_2 complementarity constant is k-independent."""
        pair = affine_sl2_koszul_pair()
        result = verify_kappa_complementarity(pair)
        assert result['level_independent']

    def test_w3_level_independent(self):
        """W_3 complementarity constant is c-independent."""
        pair = w3_koszul_pair()
        result = verify_kappa_complementarity(pair)
        assert result['level_independent']

    def test_km_vs_w_complementarity_constant(self):
        """KM/free fields have constant 0; W-algebras have constant rho*K != 0."""
        heis = verify_kappa_complementarity(heisenberg_koszul_pair())
        sl2 = verify_kappa_complementarity(affine_sl2_koszul_pair())
        vir = verify_kappa_complementarity(virasoro_koszul_pair())
        w3 = verify_kappa_complementarity(w3_koszul_pair())

        # KM/free: sum = 0
        assert heis['kappa_sum'] == '0'
        assert sl2['kappa_sum'] == '0'

        # Virasoro: sum = 13 != 0
        assert vir['kappa_sum'] == '13'

        # W_3: sum = 250/3 != 0
        assert w3['kappa_sum'] == '250/3'


# ===================================================================
# 3. LINE OPERATORS
# ===================================================================

class TestLineOperators:
    """Test line operator data for each family."""

    def test_heisenberg_lines_exist(self):
        """Heisenberg: line operators parametrized by momentum."""
        lines = heisenberg_line_operators()
        assert len(lines) >= 1

    def test_heisenberg_vacuum_weight(self):
        """Heisenberg: p=0 line has weight 0."""
        lines = heisenberg_line_operators(k_val=1)
        vacuum = [l for l in lines if l.label == S.Zero]
        assert len(vacuum) >= 1
        assert vacuum[0].conformal_weight == 0

    def test_heisenberg_p1_weight(self):
        """Heisenberg at k=1: p=1 line has weight 1/2."""
        lines = heisenberg_line_operators(k_val=1)
        p1 = [l for l in lines if l.label == S.One]
        assert len(p1) == 1
        assert p1[0].conformal_weight == Rational(1, 2)

    def test_virasoro_lines_exist(self):
        """Virasoro: line operators are Vir_{26-c} modules."""
        lines = virasoro_line_operators()
        assert len(lines) >= 1

    def test_virasoro_vacuum_line(self):
        """Virasoro: vacuum line has weight 0."""
        lines = virasoro_line_operators()
        vacuum = [l for l in lines if l.conformal_weight == S.Zero]
        assert len(vacuum) >= 1

    def test_sl2_lines_exist(self):
        """sl_2: line operators are V_{-k-4}(sl_2) modules."""
        lines = affine_sl2_line_operators()
        assert len(lines) >= 1

    def test_sl2_trivial_line(self):
        """sl_2: trivial line (j=0) has weight 0."""
        lines = affine_sl2_line_operators()
        trivial = [l for l in lines if l.label == S.Zero]
        assert len(trivial) >= 1
        assert trivial[0].conformal_weight == 0


# ===================================================================
# 4. CLASSICAL r-MATRIX AND CYBE
# ===================================================================

class TestRMatrix:
    """Classical r-matrix and Yang-Baxter equation."""

    def test_heisenberg_r_matrix_pole(self):
        """Heisenberg: r(z) has pole order 1."""
        pair = heisenberg_koszul_pair()
        result = classical_r_matrix(pair)
        assert result['pole_order'] == 1

    def test_virasoro_r_matrix_pole(self):
        """Virasoro: r(z) has pole order 3."""
        pair = virasoro_koszul_pair()
        result = classical_r_matrix(pair)
        assert result['pole_order'] == 3

    def test_sl2_r_matrix_casimir(self):
        """sl_2: r(z) = Omega/z (Casimir r-matrix)."""
        pair = affine_sl2_koszul_pair()
        result = classical_r_matrix(pair)
        assert result['pole_order'] == 1
        assert 'Omega' in result['r_matrix']

    def test_w3_r_matrix_pole(self):
        """W_3: r(z) has pole order 5."""
        pair = w3_koszul_pair()
        result = classical_r_matrix(pair)
        assert result['pole_order'] == 5

    def test_betagamma_r_matrix_pole(self):
        """Beta-gamma: r(z) has pole order 0 (simple pole)."""
        pair = betagamma_koszul_pair()
        result = classical_r_matrix(pair)
        assert result['pole_order'] == 0

    def test_cybe_abelian(self):
        """CYBE for Heisenberg: trivially holds (abelian)."""
        result = verify_cybe_abelian()
        assert result['cybe_holds']

    def test_cybe_sl2(self):
        """CYBE for sl_2: holds by Belavin-Drinfeld."""
        result = verify_cybe_sl2()
        assert result['cybe_holds']

    def test_r_matrix_shadow_interpretation(self):
        """r(z) = Res^coll_{0,2}(Theta_A) for all families."""
        for pair_fn in [heisenberg_koszul_pair, virasoro_koszul_pair,
                        affine_sl2_koszul_pair]:
            pair = pair_fn()
            result = classical_r_matrix(pair)
            assert result['shadow_interpretation'] == 'Res^coll_{0,2}(Theta_A)'


# ===================================================================
# 5. THEOREM J VERIFICATION
# ===================================================================

class TestTheoremJ:
    """Verify Theorem J: lines = A!-modules."""

    def test_theorem_j_heisenberg(self):
        """Theorem J for Heisenberg: lines = H_k^! modules."""
        result = verify_theorem_j_heisenberg()
        assert result['theorem_j_holds']
        assert result['koszul_dual'] == 'Sym^ch(V*)'

    def test_theorem_j_virasoro(self):
        """Theorem J for Virasoro: lines = Vir_{26-c} modules."""
        result = verify_theorem_j_virasoro()
        assert result['theorem_j_holds']
        assert result['koszul_dual'] == 'Vir_{26-c}'

    def test_theorem_j_affine(self):
        """Theorem J for affine sl_2: lines = V_{-k-4}(sl_2) modules."""
        result = verify_theorem_j_affine()
        assert result['theorem_j_holds']
        assert 'Feigin-Frenkel' in result['reason']

    def test_theorem_j_virasoro_self_dual(self):
        """At c=13: lines are self-conjugate."""
        result = verify_theorem_j_virasoro(c_val=13)
        assert result['is_self_dual']

    def test_theorem_j_virasoro_c26(self):
        """At c=26: dual has c'=0."""
        result = verify_theorem_j_virasoro(c_val=26)
        assert simplify(S(result['dual_central_charge'])) == 0


# ===================================================================
# 6. CROSS-VOLUME BRIDGE
# ===================================================================

class TestCrossVolumeBridge:
    """Cross-volume bridge: Vol I data = Vol II data."""

    def test_kappa_bridge_heisenberg(self):
        """Heisenberg: Vol I kappa = Vol II kappa."""
        pair = heisenberg_koszul_pair()
        result = cross_volume_kappa_bridge(pair)
        assert result['match']

    def test_kappa_bridge_virasoro(self):
        """Virasoro: Vol I kappa = Vol II kappa."""
        pair = virasoro_koszul_pair()
        result = cross_volume_kappa_bridge(pair)
        assert result['match']

    def test_shadow_class_heisenberg(self):
        """Heisenberg: shadow class G consistent across volumes."""
        pair = heisenberg_koszul_pair()
        result = cross_volume_shadow_class_bridge(pair)
        assert result['consistent']

    def test_shadow_class_virasoro(self):
        """Virasoro: shadow class M consistent."""
        pair = virasoro_koszul_pair()
        result = cross_volume_shadow_class_bridge(pair)
        assert result['consistent']

    def test_shadow_class_sl2(self):
        """sl_2: shadow class L consistent."""
        pair = affine_sl2_koszul_pair()
        result = cross_volume_shadow_class_bridge(pair)
        assert result['consistent']

    def test_bridge_theorem_reference(self):
        """Bridge references thm:hochschild-bridge-genus0."""
        pair = heisenberg_koszul_pair()
        result = cross_volume_kappa_bridge(pair)
        assert 'hochschild-bridge-genus0' in result['bridge_theorem']


# ===================================================================
# 7. SELF-DUALITY ANALYSIS
# ===================================================================

class TestSelfDuality:
    """Self-duality analysis (AP8 compliance)."""

    def test_virasoro_self_dual_at_13(self):
        """Virasoro: self-dual at c=13."""
        pair = virasoro_koszul_pair()
        result = self_duality_analysis(pair)
        assert result['self_dual_point'] == '13'

    def test_heisenberg_not_self_dual(self):
        """Heisenberg: NEVER self-dual."""
        pair = heisenberg_koszul_pair()
        result = self_duality_analysis(pair)
        assert not result['is_family_self_dual']

    def test_w3_self_dual_at_50(self):
        """W_3: self-dual at c=50."""
        pair = w3_koszul_pair()
        result = self_duality_analysis(pair)
        assert result['self_dual_point'] == '50'

    def test_ap8_warning(self):
        """All self-duality analyses carry AP8 warning."""
        for pair_fn in [heisenberg_koszul_pair, virasoro_koszul_pair,
                        w3_koszul_pair]:
            pair = pair_fn()
            result = self_duality_analysis(pair)
            assert 'AP8' in result['warning']


# ===================================================================
# 8. YANGIAN SHADOW DATA
# ===================================================================

class TestYangianShadow:
    """Yangian shadow data from r-matrices."""

    def test_heisenberg_trivial_yangian(self):
        """Heisenberg: trivial Yangian (abelian)."""
        pair = heisenberg_koszul_pair()
        result = yangian_shadow_from_r_matrix(pair)
        assert result['yangian_type'] == 'trivial'

    def test_sl2_dg_shifted_yangian(self):
        """sl_2: dg-shifted Yangian (non-abelian)."""
        pair = affine_sl2_koszul_pair()
        result = yangian_shadow_from_r_matrix(pair)
        assert result['yangian_type'] == 'dg-shifted'

    def test_virasoro_quantum_corrections(self):
        """Virasoro: has quantum corrections (infinite shadow obstruction tower)."""
        pair = virasoro_koszul_pair()
        result = yangian_shadow_from_r_matrix(pair)
        assert result['quantum_corrections']

    def test_heisenberg_no_quantum_corrections(self):
        """Heisenberg: NO quantum corrections (Gaussian class)."""
        pair = heisenberg_koszul_pair()
        result = yangian_shadow_from_r_matrix(pair)
        assert not result['quantum_corrections']


# ===================================================================
# 9. COMPREHENSIVE DUALITY TABLE
# ===================================================================

class TestDualityTable:
    """Comprehensive duality table across all families."""

    def test_table_has_five_families(self):
        """Table covers all 5 standard families."""
        result = comprehensive_duality_table()
        assert result['num_families'] == 5

    def test_all_complementarity_holds(self):
        """All families satisfy kappa complementarity."""
        result = comprehensive_duality_table()
        assert result['all_complementarity_holds']

    def test_table_entries_complete(self):
        """Each table entry has all required fields."""
        result = comprehensive_duality_table()
        required_fields = ['algebra', 'shadow_class', 'kappa', 'dual_kappa',
                          'complementarity', 'kappa_sum', 'r_matrix_pole',
                          'yangian_type', 'self_dual_point']
        for entry in result['table']:
            for field in required_fields:
                assert field in entry, f"Missing field {field} in {entry['algebra']}"

    def test_shadow_classes_present(self):
        """All four shadow classes (G, L, C, M) are represented."""
        result = comprehensive_duality_table()
        classes = set(entry['shadow_class'] for entry in result['table'])
        assert 'G' in classes
        assert 'L' in classes
        assert 'C' in classes
        assert 'M' in classes


# ===================================================================
# 10. CONFORMAL WEIGHT SPECTRUM
# ===================================================================

class TestConformalWeightSpectrum:
    """Conformal weight spectrum of line operators."""

    def test_heisenberg_spectrum_continuous(self):
        """Heisenberg: continuous spectrum."""
        pair = heisenberg_koszul_pair()
        result = conformal_weight_spectrum(pair)
        assert result['continuous']

    def test_virasoro_spectrum_continuous(self):
        """Virasoro at generic c: continuous spectrum."""
        pair = virasoro_koszul_pair()
        result = conformal_weight_spectrum(pair)
        assert result['continuous']

    def test_sl2_spectrum_discrete(self):
        """sl_2 at integer level: discrete spectrum."""
        pair = affine_sl2_koszul_pair()
        result = conformal_weight_spectrum(pair)
        assert result['discrete_at_integer_level']

    def test_heisenberg_vacuum_weight(self):
        """Heisenberg: p=0 has weight 0."""
        pair = heisenberg_koszul_pair()
        result = conformal_weight_spectrum(pair)
        assert result['spectrum']['p=0'] == 0

    def test_heisenberg_p1_weight(self):
        """Heisenberg at k=1: p=1 has weight 1/2."""
        pair = heisenberg_koszul_pair()
        result = conformal_weight_spectrum(pair)
        assert result['spectrum']['p=1'] == Rational(1, 2)


# ===================================================================
# 11. ANOMALY COMPLETION
# ===================================================================

class TestAnomalyCompletion:
    """Anomaly-completed Koszul duality data."""

    def test_anomaly_data_heisenberg(self):
        """Heisenberg: anomaly completion data exists."""
        pair = heisenberg_koszul_pair()
        result = anomaly_completion_data(pair)
        assert 'holographic_datum_fields' in result
        assert len(result['holographic_datum_fields']) == 6

    def test_anomaly_data_virasoro(self):
        """Virasoro: anomaly completion data."""
        pair = virasoro_koszul_pair()
        result = anomaly_completion_data(pair)
        assert 'genus_1_correction' in result


# ===================================================================
# 12. KOSZUL DUAL CENTRAL CHARGES
# ===================================================================

class TestDualCentralCharges:
    """Verify Koszul dual central charges are correct."""

    def test_virasoro_dual_c(self):
        """Virasoro: c(A!) = 26-c."""
        c = Symbol('c')
        pair = virasoro_koszul_pair()
        assert simplify(pair.dual_central_charge - (26 - c)) == 0

    def test_sl2_dual_level(self):
        """sl_2: dual level k' = -k-4."""
        k = Symbol('k')
        pair = affine_sl2_koszul_pair()
        # dual_central_charge = 3*(-k-4)/(-k-4+2) = 3*(-k-4)/(-k-2)
        expected_c_dual = 3 * (-k - 4) / (-k - 4 + 2)
        assert simplify(pair.dual_central_charge - expected_c_dual) == 0

    def test_w3_dual_c(self):
        """W_3: c(A!) = 100 - c."""
        c = Symbol('c')
        pair = w3_koszul_pair()
        assert simplify(pair.dual_central_charge - (100 - c)) == 0

    def test_virasoro_at_c13_dual_c13(self):
        """At c=13: Vir_13^! = Vir_13 (self-dual)."""
        pair = virasoro_koszul_pair(c_val=13)
        assert pair.dual_central_charge == 13

    def test_virasoro_at_c0_dual_c26(self):
        """At c=0: Vir_0^! = Vir_26."""
        pair = virasoro_koszul_pair(c_val=0)
        assert pair.dual_central_charge == 26

    def test_virasoro_at_c26_dual_c0(self):
        """At c=26: Vir_26^! = Vir_0."""
        pair = virasoro_koszul_pair(c_val=26)
        assert pair.dual_central_charge == 0
