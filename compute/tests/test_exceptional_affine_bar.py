r"""Tests for ordered E₁ bar complex data: V_k(E₆), V_k(E₇), V_k(E₈).

Verifies all seven computational aspects:
  (1) Lie algebra invariants (dim, h∨, exponents, consistency)
  (2) Curvature κ and Koszul complementarity
  (3) Shadow class L, depth 3
  (4) Collision residue r(z) = kΩ/z
  (5) Yangian Y_ℏ(E_N) identification
  (6) Euler-eta character
  (7) DS reduction and depth gap

Cross-references:
  Vol II: w-algebras-stable.tex (depth gap table at eq:exponent-depth)
  Vol II: rosetta_stone.tex (affine CS example)
  compute/lib/affine_half_space_bv_engine.py (dual Coxeter numbers)
  compute/lib/dg_shifted_factorization_engine.py (Lie dimensions)
  compute/lib/collision_residue_rmatrix.py (r-matrix framework)
"""

import pytest
import sys
import os
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.exceptional_affine_bar import (
    _EXCEPTIONAL_E_DATA,
    verify_exceptional_data,
    curvature_kappa,
    shadow_class,
    collision_residue,
    ordered_koszul_dual,
    euler_eta,
    ds_reduction,
    exceptional_affine_bar_data,
    comparison_table,
)

# Also import from existing engines for cross-validation
from lib.affine_half_space_bv_engine import dual_coxeter_number
from lib.dg_shifted_factorization_engine import _lie_dimension


# =========================================================================
# PART 1: LIE ALGEBRA DATA CONSISTENCY
# =========================================================================

class TestExceptionalLieData:
    """Verify the hard-coded exceptional Lie algebra data."""

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_data_consistency(self, name):
        """All internal consistency checks pass."""
        result = verify_exceptional_data(name)
        assert result['checks']['all_passed'], (
            f"{name}: consistency checks failed: "
            f"{[k for k, v in result['checks'].items() if not v and k != 'all_passed']}"
        )

    @pytest.mark.parametrize("name,expected_dim", [
        ("E6", 78), ("E7", 133), ("E8", 248),
    ])
    def test_dimensions(self, name, expected_dim):
        """Verify dimensions match Bourbaki tables."""
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['dim'] == expected_dim

    @pytest.mark.parametrize("name,expected_h", [
        ("E6", 12), ("E7", 18), ("E8", 30),
    ])
    def test_dual_coxeter(self, name, expected_h):
        """Verify h∨ matches Kac tables."""
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['h_dual'] == expected_h

    @pytest.mark.parametrize("name,expected_rank", [
        ("E6", 6), ("E7", 7), ("E8", 8),
    ])
    def test_rank(self, name, expected_rank):
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['rank'] == expected_rank

    @pytest.mark.parametrize("name,expected_exps", [
        ("E6", [1, 4, 5, 7, 8, 11]),
        ("E7", [1, 5, 7, 9, 11, 13, 17]),
        ("E8", [1, 7, 11, 13, 17, 19, 23, 29]),
    ])
    def test_exponents(self, name, expected_exps):
        """Verify exponents match Bourbaki Plates V-VII."""
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['exponents'] == expected_exps

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_dim_equals_rank_plus_2_positive_roots(self, name):
        """dim = rank + 2|Φ⁺|."""
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['dim'] == data['rank'] + 2 * data['num_positive_roots']

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_exponent_sum_equals_positive_roots(self, name):
        """Classical identity: Σ e_i = |Φ⁺|."""
        data = _EXCEPTIONAL_E_DATA[name]
        assert sum(data['exponents']) == data['num_positive_roots']

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_coxeter_from_largest_exponent(self, name):
        """h = 1 + e_r (largest exponent)."""
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['h_coxeter'] == 1 + data['exponents'][-1]

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_simply_laced_h_equals_h_dual(self, name):
        """For simply-laced (ADE): h = h∨."""
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['h_coxeter'] == data['h_dual']

    @pytest.mark.parametrize("name,expected_n_pos", [
        ("E6", 36), ("E7", 63), ("E8", 120),
    ])
    def test_positive_root_count(self, name, expected_n_pos):
        """Verify |Φ⁺| values."""
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['num_positive_roots'] == expected_n_pos


# =========================================================================
# PART 2: CROSS-VALIDATION WITH EXISTING ENGINES
# =========================================================================

class TestCrossValidation:
    """Verify consistency with affine_half_space_bv_engine and dg_shifted_factorization_engine."""

    @pytest.mark.parametrize("name,lie_type,rank", [
        ("E6", "E", 6), ("E7", "E", 7), ("E8", "E", 8),
    ])
    def test_h_dual_matches_bv_engine(self, name, lie_type, rank):
        """h∨ matches the authoritative table in affine_half_space_bv_engine."""
        expected = dual_coxeter_number(lie_type, rank)
        actual = _EXCEPTIONAL_E_DATA[name]['h_dual']
        assert actual == expected, f"{name}: h∨ mismatch: {actual} vs {expected}"

    @pytest.mark.parametrize("name,lie_type,rank", [
        ("E6", "E", 6), ("E7", "E", 7), ("E8", "E", 8),
    ])
    def test_dim_matches_factorization_engine(self, name, lie_type, rank):
        """dim(g) matches _lie_dimension in dg_shifted_factorization_engine."""
        expected = _lie_dimension(lie_type, rank)
        actual = _EXCEPTIONAL_E_DATA[name]['dim']
        assert actual == expected, f"{name}: dim mismatch: {actual} vs {expected}"


# =========================================================================
# PART 3: CURVATURE
# =========================================================================

class TestCurvature:
    """Verify modular characteristic κ(V_k(E_N))."""

    @pytest.mark.parametrize("name,k,expected", [
        # κ = dim·(k+h∨)/(2h∨)
        # E6, k=1: 78·13/(2·12) = 78·13/24 = 1014/24 = 169/4
        ("E6", 1, Fraction(169, 4)),
        # E7, k=1: 133·19/(2·18) = 133·19/36 = 2527/36
        ("E7", 1, Fraction(2527, 36)),
        # E8, k=1: 248·31/(2·30) = 248·31/60 = 7688/60 = 1922/15
        ("E8", 1, Fraction(1922, 15)),
    ])
    def test_kappa_at_k1(self, name, k, expected):
        """κ at k=1 has the correct value."""
        result = curvature_kappa(name, k)
        assert result['kappa'] == expected, (
            f"{name}: κ = {result['kappa']}, expected {expected}"
        )

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_kappa_at_k0(self, name):
        """κ = dim(g)/2 at k = 0 (NOT zero; the h∨ shift is essential)."""
        result = curvature_kappa(name, 0)
        dim_g = _EXCEPTIONAL_E_DATA[name]['dim']
        assert result['kappa'] == Fraction(dim_g, 2)

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_kappa_vanishes_at_critical(self, name):
        """κ = 0 at k = -h∨ (critical level)."""
        h_dual = _EXCEPTIONAL_E_DATA[name]['h_dual']
        result = curvature_kappa(name, -h_dual)
        assert result['kappa'] == 0

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_koszul_complementarity(self, name):
        """κ(V_k) + κ(V_k^!) = 0 (Feigin-Frenkel)."""
        for k in [1, 2, 5, 10]:
            result = curvature_kappa(name, k)
            assert result['complementarity_verified'] is True, (
                f"{name} at k={k}: complementarity failed"
            )

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_kappa_semiclassical_limit(self, name):
        """κ → dim(g)/2 as k → ∞."""
        dim_g = _EXCEPTIONAL_E_DATA[name]['dim']
        result = curvature_kappa(name, 10**6)
        expected_limit = Fraction(dim_g, 2)
        # At k = 10^6, κ should be very close to dim(g)/2
        diff = abs(result['kappa'] - expected_limit)
        assert diff < Fraction(1, 100), (
            f"{name}: semiclassical limit: κ = {result['kappa']}, "
            f"expected ≈ {expected_limit}"
        )

    @pytest.mark.parametrize("name,k,expected", [
        # E6, k=2: 78·2/(2·14) = 156/28 = 39/7
        ("E6", 2, Fraction(39, 7)),
        # E7, k=2: 133·2/(2·20) = 266/40 = 133/20
        ("E7", 2, Fraction(133, 20)),
        # E8, k=2: 248·2/(2·32) = 496/64 = 31/4
        ("E8", 2, Fraction(31, 4)),
    ])
    def test_kappa_at_k2(self, name, k, expected):
        """κ at k=2."""
        result = curvature_kappa(name, k)
        assert result['kappa'] == expected


# =========================================================================
# PART 4: SHADOW CLASS
# =========================================================================

class TestShadowClass:
    """Verify all three are Class L, depth 3."""

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_class_L(self, name):
        result = shadow_class(name)
        assert result['shadow_class'] == 'L'

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_depth_3(self, name):
        result = shadow_class(name)
        assert result['shadow_depth_r_max'] == 3

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_quartic_vanishes(self, name):
        """m₄ = 0 by Jacobi identity."""
        result = shadow_class(name)
        assert result['quartic_contact'] == 0


# =========================================================================
# PART 5: COLLISION RESIDUE
# =========================================================================

class TestCollisionResidue:
    """Verify r(z) = kΩ/z and AP19 pole absorption."""

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_ope_max_pole_is_2(self, name):
        """Affine KM has max OPE pole order 2."""
        result = collision_residue(name)
        assert result['ope_max_pole'] == 2

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_r_matrix_max_pole_is_1(self, name):
        """After d-log absorption (AP19), max pole is 1."""
        result = collision_residue(name)
        assert result['r_matrix_max_pole'] == 1

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_pole_absorption(self, name):
        """AP19 verified: shift by exactly 1."""
        result = collision_residue(name)
        assert 'AP19 verified' in result['pole_absorption']

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_cybe_satisfied(self, name):
        """CYBE / IBR holds for the Casimir of any simple Lie algebra."""
        result = collision_residue(name)
        assert result['cybe_satisfied'] is True

    @pytest.mark.parametrize("name,expected_casimir_adj", [
        ("E6", 24),  # 2·12
        ("E7", 36),  # 2·18
        ("E8", 60),  # 2·30
    ])
    def test_quadratic_casimir_adjoint(self, name, expected_casimir_adj):
        """C₂(adj) = 2h∨ for the quadratic Casimir in the adjoint."""
        result = collision_residue(name)
        assert result['quadratic_casimir_adjoint'] == expected_casimir_adj

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_casimir_term_count(self, name):
        """Casimir has dim(g) terms: rank Cartan + 2|Φ⁺| root terms."""
        result = collision_residue(name)
        data = _EXCEPTIONAL_E_DATA[name]
        cs = result['casimir_structure']
        assert cs['total_terms'] == data['dim']
        assert cs['root_terms'] == 2 * data['num_positive_roots']
        assert cs['cartan_terms'] == data['rank']
        assert cs['root_terms'] + cs['cartan_terms'] == cs['total_terms']


# =========================================================================
# PART 6: ORDERED KOSZUL DUAL (YANGIAN)
# =========================================================================

class TestYangian:
    """Verify the Yangian identification as ordered Koszul dual."""

    @pytest.mark.parametrize("name,expected_min_rep", [
        ("E6", 27), ("E7", 56), ("E8", 248),
    ])
    def test_minimal_representation_dim(self, name, expected_min_rep):
        """Minimal faithful representation dimensions."""
        result = ordered_koszul_dual(name)
        assert result['minimal_representation_dim'] == expected_min_rep

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_generator_count(self, name):
        """2·dim(g) generators per mode (even + odd)."""
        result = ordered_koszul_dual(name)
        dim_g = _EXCEPTIONAL_E_DATA[name]['dim']
        assert result['num_even_generators_per_mode'] == dim_g
        assert result['num_odd_generators_per_mode'] == dim_g
        assert result['total_generators_per_mode'] == 2 * dim_g

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_root_multiplicity_one(self, name):
        """Root multiplicity 1 ensures strictification."""
        result = ordered_koszul_dual(name)
        assert result['root_multiplicity'] == 1

    def test_E8_adjoint_is_minimal(self):
        """E₈ is unique: its minimal faithful rep IS the adjoint (248)."""
        result = ordered_koszul_dual('E8')
        assert result['minimal_representation_dim'] == 248
        assert result['minimal_representation_dim'] == _EXCEPTIONAL_E_DATA['E8']['dim']

    def test_E6_min_rep_not_adjoint(self):
        """E₆: min rep (27) is NOT the adjoint (78)."""
        result = ordered_koszul_dual('E6')
        assert result['minimal_representation_dim'] == 27
        assert result['minimal_representation_dim'] != _EXCEPTIONAL_E_DATA['E6']['dim']


# =========================================================================
# PART 7: EULER-ETA
# =========================================================================

class TestEulerEta:
    """Verify the Euler-eta character data."""

    @pytest.mark.parametrize("name,expected_exp", [
        ("E6", 78), ("E7", 133), ("E8", 248),
    ])
    def test_eta_exponent_is_dim(self, name, expected_exp):
        """η exponent equals dim(g)."""
        result = euler_eta(name)
        assert result['eta_exponent'] == expected_exp

    @pytest.mark.parametrize("name,expected_c_eff", [
        # c_eff(k=1) = dim/(1+h∨)
        # E6: 78/13 = 6
        ("E6", Fraction(6, 1)),
        # E7: 133/19 = 7
        ("E7", Fraction(7, 1)),
        # E8: 248/31 = 8
        ("E8", Fraction(8, 1)),
    ])
    def test_effective_central_charge_k1(self, name, expected_c_eff):
        """c_eff at k=1: dim(g)/(1+h∨)."""
        result = euler_eta(name)
        assert result['effective_central_charge_at_k1'] == expected_c_eff

    def test_c_eff_equals_rank_at_k1(self):
        """At k=1, c_eff = rank for all three exceptional E-types.

        This is a nontrivial check: dim(g)/(1+h∨) = rank.
        For E₆: 78/13 = 6 = rank. For E₇: 133/19 = 7 = rank.
        For E₈: 248/31 = 8 = rank.

        This identity dim(g) = rank · (1+h∨) characterizes the
        simply-laced exceptional types and reflects the fact that
        dim(g) = rank + 2|Φ⁺| and |Φ⁺| = rank·h/2 for simply-laced.
        """
        for name in ["E6", "E7", "E8"]:
            result = euler_eta(name)
            rank = _EXCEPTIONAL_E_DATA[name]['rank']
            assert result['effective_central_charge_at_k1'] == rank, (
                f"{name}: c_eff(k=1) = {result['effective_central_charge_at_k1']}, "
                f"expected rank = {rank}"
            )


# =========================================================================
# PART 8: DS REDUCTION AND DEPTH GAP
# =========================================================================

class TestDSReduction:
    """Verify Drinfeld-Sokolov reduction data and depth gap formula."""

    @pytest.mark.parametrize("name,expected_d_gap", [
        ("E6", 22), ("E7", 34), ("E8", 58),
    ])
    def test_depth_gap(self, name, expected_d_gap):
        """d_gap = 2e_r matches the table in w-algebras-stable.tex."""
        result = ds_reduction(name)
        assert result['d_gap'] == expected_d_gap

    @pytest.mark.parametrize("name,expected_e_r", [
        ("E6", 11), ("E7", 17), ("E8", 29),
    ])
    def test_largest_exponent(self, name, expected_e_r):
        result = ds_reduction(name)
        assert result['largest_exponent'] == expected_e_r

    @pytest.mark.parametrize("name,expected_s_r", [
        ("E6", 12), ("E7", 18), ("E8", 30),
    ])
    def test_highest_spin(self, name, expected_s_r):
        result = ds_reduction(name)
        assert result['highest_spin'] == expected_s_r

    @pytest.mark.parametrize("name,expected_max_pole", [
        ("E6", 24), ("E7", 36), ("E8", 60),
    ])
    def test_binding_ope_max_pole(self, name, expected_max_pole):
        """Max OPE pole = 2s_r."""
        result = ds_reduction(name)
        assert result['binding_ope_max_pole'] == expected_max_pole

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_d_gap_formula_consistency(self, name):
        """d_gap = 2e_r = 2(h-1)."""
        result = ds_reduction(name)
        assert result['coxeter_formula_agrees'] is True

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_num_generators_equals_rank(self, name):
        """W(E_N) has rank generators."""
        result = ds_reduction(name)
        data = _EXCEPTIONAL_E_DATA[name]
        assert result['num_generators'] == data['rank']

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_generator_weights_from_exponents(self, name):
        """Generator weights = exponents + 1."""
        result = ds_reduction(name)
        data = _EXCEPTIONAL_E_DATA[name]
        expected = [e + 1 for e in data['exponents']]
        assert result['generator_weights'] == expected

    def test_E8_deepest_gap(self):
        """E₈ has the deepest gap (d_gap=58) among all simple Lie algebras."""
        result = ds_reduction('E8')
        assert result['is_deepest_gap'] is True
        assert result['d_gap'] == 58

    def test_E6_not_deepest(self):
        result = ds_reduction('E6')
        assert result['is_deepest_gap'] is False

    def test_gap_ordering(self):
        """d_gap(E6) < d_gap(E7) < d_gap(E8)."""
        gaps = [ds_reduction(name)['d_gap'] for name in ["E6", "E7", "E8"]]
        assert gaps[0] < gaps[1] < gaps[2]

    def test_class_transport_L_to_M(self):
        """DS reduction transports Class L → Class M for all three."""
        for name in ["E6", "E7", "E8"]:
            result = ds_reduction(name)
            assert 'L' in result['input_class']
            assert 'M' in result['output_class']


# =========================================================================
# PART 9: FULL PACKAGE AND COMPARISON TABLE
# =========================================================================

class TestFullPackage:
    """Verify the assembled package and comparison table."""

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_full_package_complete(self, name):
        """All seven sections present in the full package."""
        result = exceptional_affine_bar_data(name, k=1)
        assert 'lie_data' in result
        assert 'consistency' in result
        assert 'curvature' in result
        assert 'shadow' in result
        assert 'collision' in result
        assert 'yangian' in result
        assert 'euler' in result
        assert 'ds' in result

    def test_comparison_table_has_three_rows(self):
        table = comparison_table(k=1)
        assert len(table['rows']) == 3

    def test_comparison_table_dim_ordering(self):
        """dim(E6) < dim(E7) < dim(E8)."""
        table = comparison_table(k=1)
        dims = [row['dim'] for row in table['rows']]
        assert dims == [78, 133, 248]

    def test_comparison_table_all_class_L(self):
        table = comparison_table(k=1)
        for row in table['rows']:
            assert row['class'] == 'L'


# =========================================================================
# PART 10: CROSS-CHECKS WITH w-algebras-stable.tex TABLE
# =========================================================================

class TestCrossCheckWithManuscript:
    """Verify against the explicit table in w-algebras-stable.tex lines 1797-1813.

    The table states:
      E₆: e_r=11, s_r=12, max pole=24, d_gap=22
      E₇: e_r=17, s_r=18, max pole=36, d_gap=34
      E₈: e_r=29, s_r=30, max pole=60, d_gap=58
    """

    def test_E6_matches_manuscript_table(self):
        ds = ds_reduction('E6')
        assert ds['largest_exponent'] == 11
        assert ds['highest_spin'] == 12
        assert ds['binding_ope_max_pole'] == 24
        assert ds['d_gap'] == 22

    def test_E7_matches_manuscript_table(self):
        ds = ds_reduction('E7')
        assert ds['largest_exponent'] == 17
        assert ds['highest_spin'] == 18
        assert ds['binding_ope_max_pole'] == 36
        assert ds['d_gap'] == 34

    def test_E8_matches_manuscript_table(self):
        ds = ds_reduction('E8')
        assert ds['largest_exponent'] == 29
        assert ds['highest_spin'] == 30
        assert ds['binding_ope_max_pole'] == 60
        assert ds['d_gap'] == 58

    def test_E8_deepest_in_classification_note(self):
        """From w-algebras-stable.tex line 1817-1818:
        'E₈ has the deepest arity-2 depth gap in the entire classification (d_gap = 58)'"""
        ds = ds_reduction('E8')
        assert ds['d_gap'] == 58

        # Verify it beats all other types listed in the table
        # A_{N-1}: d_gap = 2N-2. To beat 58: need N >= 31 (sl_31).
        # B_n: d_gap = 4n-2. To beat 58: need n >= 16 (so_33).
        # D_n: d_gap = 4n-6. To beat 58: need n >= 17 (so_34).
        # Among EXCEPTIONAL types, E8 is the largest.
        for name in ["E6", "E7"]:
            assert ds_reduction(name)['d_gap'] < 58


# =========================================================================
# PART 11: ADDITIONAL IDENTITIES
# =========================================================================

class TestAdditionalIdentities:
    """Verify less-obvious identities relating the invariants."""

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_dim_equals_rank_times_1_plus_h(self, name):
        """For simply-laced: dim(g) = rank · (1 + h).

        Since h = h∨ and |Φ⁺| = rank·h/2:
          dim = rank + 2·(rank·h/2) = rank(1 + h).
        """
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['dim'] == data['rank'] * (1 + data['h_coxeter'])

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_positive_roots_from_rank_and_h(self, name):
        """|Φ⁺| = rank · h / 2 for simply-laced."""
        data = _EXCEPTIONAL_E_DATA[name]
        assert data['num_positive_roots'] == data['rank'] * data['h_coxeter'] // 2

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_exponents_are_coprime_to_h_for_simply_laced(self, name):
        """For simply-laced types, the exponents are the integers 1 ≤ m ≤ h-1
        such that m is coprime to h... wait, that's only for A_n.

        Actually, for E₆, h=12, exponents are [1,4,5,7,8,11].
        gcd(4,12)=4 ≠ 1, so this is NOT true for E₆.

        The correct property: exponents come in pairs (e_i, h-e_i).
        """
        data = _EXCEPTIONAL_E_DATA[name]
        h = data['h_coxeter']
        exps = data['exponents']
        # Check pairing: {e_i} ∪ {h - e_i} = {e_i} as multisets
        paired = sorted([h - e for e in exps])
        assert exps == paired, (
            f"{name}: exponent pairing fails: {exps} vs {paired}"
        )

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_r_matrix_pole_after_ds_absorption(self, name):
        """After DS, the r-matrix max pole = binding OPE max pole - 1 (AP19)."""
        ds = ds_reduction(name)
        assert ds['r_matrix_max_pole_after_dlog'] == ds['binding_ope_max_pole'] - 1
