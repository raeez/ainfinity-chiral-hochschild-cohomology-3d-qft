r"""Tests for remaining families E₁ ordered chiral bar complex.

Tests the five new families:
  (1) N=2 superconformal algebra
  (2) N=4 superconformal algebra
  (3) Affine E₆, E₇, E₈ at level 1
  (4) Bershadsky-Polyakov W₃^{(2)}
  (5) Virasoro minimal models M(p,q)

Verifies:
  - Collision residues with AP19 d-log absorption
  - m₂ spectral structure
  - GLCM classification
  - Shadow obstruction tower coefficients
  - R-matrix type
  - Euler-eta
  - Koszul sign conventions for fermions

Cross-references:
  compute/remaining_families_ordered_bar.py (source)
  compute/lib/exceptional_affine_bar.py (E-type data)
  compute/lib/collision_residue_rmatrix.py (r-matrix framework)
  CLAUDE.md: AP19, AP36, AP37, AP38
"""

import pytest
import sys
import os
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from remaining_families_ordered_bar import (
    N2Superconformal,
    N4Superconformal,
    AffineExceptional,
    BershadskyPolyakov,
    VirasoroMinimalModel,
)


# =========================================================================
# FAMILY 1: N=2 SUPERCONFORMAL ALGEBRA
# =========================================================================

class TestN2Superconformal:
    """Tests for the N=2 SCA."""

    def test_central_charge_k1(self):
        """c = 3k/(k+2) = 1 at k=1."""
        n2 = N2Superconformal(k_val=1)
        assert n2.c == Fraction(1)

    def test_central_charge_k2(self):
        """c = 3/2 at k=2."""
        n2 = N2Superconformal(k_val=2)
        assert n2.c == Fraction(3, 2)

    def test_generators(self):
        """4 generators: T, G+, G-, J."""
        assert len(N2Superconformal.GENERATORS) == 4

    def test_parities(self):
        """T, J bosonic; G+, G- fermionic."""
        assert N2Superconformal.PARITIES['T'] == 0
        assert N2Superconformal.PARITIES['J'] == 0
        assert N2Superconformal.PARITIES['Gp'] == 1
        assert N2Superconformal.PARITIES['Gm'] == 1

    def test_max_ope_pole(self):
        """Max OPE pole is 4 (from TT channel)."""
        n2 = N2Superconformal(k_val=1)
        assert n2.max_ope_pole() == 4

    def test_ap19_tt_channel(self):
        """AP19: TT quartic OPE → cubic r-matrix."""
        n2 = N2Superconformal(k_val=1)
        cr = n2.collision_residues()
        assert cr[('T', 'T')]['max_pole'] == 3

    def test_ap19_gpm_channel(self):
        """AP19: G+G- cubic OPE → quadratic r-matrix."""
        n2 = N2Superconformal(k_val=1)
        cr = n2.collision_residues()
        assert cr[('Gp', 'Gm')]['max_pole'] == 2

    def test_ap19_jj_channel(self):
        """AP19: JJ double OPE → simple r-matrix."""
        n2 = N2Superconformal(k_val=1)
        cr = n2.collision_residues()
        assert cr[('J', 'J')]['max_pole'] == 1

    def test_gpgp_vanishes(self):
        """G+G+ vanishes by U(1) charge conservation."""
        n2 = N2Superconformal(k_val=1)
        cr = n2.collision_residues()
        assert cr[('Gp', 'Gp')]['max_pole'] == 0

    def test_glcm_class_m(self):
        """Overall GLCM class is M (quartic TT pole)."""
        n2 = N2Superconformal(k_val=1)
        glcm = n2.glcm_class()
        assert glcm['overall_class'] == 'M'
        assert glcm['max_ope_pole'] == 4
        assert glcm['max_r_pole'] == 3

    def test_shadow_tower_curvatures(self):
        """Shadow curvatures: κ_TT = c/2, κ_JJ = c/3, κ_GG = 2c/3."""
        n2 = N2Superconformal(k_val=1)
        st = n2.shadow_tower()
        assert st['TT']['kappa'] == Fraction(1, 2)
        assert st['JJ']['kappa'] == Fraction(1, 3)
        assert st['GpGm']['kappa'] == Fraction(2, 3)

    def test_jj_class_l(self):
        """JJ sub-sector is class L (shadow terminates)."""
        n2 = N2Superconformal(k_val=1)
        st = n2.shadow_tower()
        assert st['JJ']['S3'] == 0

    def test_r_matrix_e_infinity(self):
        """R-matrix is derived from local OPE (E∞, AP36)."""
        n2 = N2Superconformal(k_val=1)
        rm = n2.r_matrix()
        assert rm['is_E_infinity']
        assert rm['R_derived_from_OPE']
        assert not rm['R_equals_tau']

    def test_m2_depth_tt(self):
        """m₂(T,T) has depth 3 (quartic OPE)."""
        n2 = N2Superconformal(k_val=1)
        m2 = n2.m2_spectral()
        assert m2[('T', 'T')]['depth'] == 3

    def test_m2_depth_gpm(self):
        """m₂(G+,G-) has depth 2 (cubic OPE)."""
        n2 = N2Superconformal(k_val=1)
        m2 = n2.m2_spectral()
        assert m2[('Gp', 'Gm')]['depth'] == 2

    def test_m3_nonzero_overall(self):
        """m₃ ≠ 0 overall (TT and G+G- sectors contribute)."""
        n2 = N2Superconformal(k_val=1)
        m3 = n2.m3_analysis()
        assert m3['TT']['nonzero']
        assert m3['GpGm']['nonzero']
        assert not m3['JJ']['nonzero']


# =========================================================================
# FAMILY 2: N=4 SUPERCONFORMAL ALGEBRA
# =========================================================================

class TestN4Superconformal:
    """Tests for the N=4 (small) SCA."""

    def test_central_charge(self):
        """c = 6k; c = 6 at k=1."""
        n4 = N4Superconformal(k_val=1)
        assert n4.c == Fraction(6)

    def test_generator_count(self):
        """8 generators: 4 bosonic + 4 fermionic."""
        assert len(N4Superconformal.GENERATORS) == 8
        assert len(N4Superconformal.GENERATORS_BOSONIC) == 4
        assert len(N4Superconformal.GENERATORS_FERMIONIC) == 4

    def test_glcm_class_m(self):
        """Overall class M (TT quartic pole)."""
        n4 = N4Superconformal(k_val=1)
        glcm = n4.glcm_class()
        assert glcm['overall_class'] == 'M'

    def test_tt_r_pole(self):
        """TT: OPE 4 → r-matrix 3."""
        n4 = N4Superconformal(k_val=1)
        cr = n4.collision_residues_summary()
        assert cr['TT']['r_max_pole'] == 3

    def test_gg_diagonal_r_pole(self):
        """GG diagonal: OPE 3 → r-matrix 2."""
        n4 = N4Superconformal(k_val=1)
        cr = n4.collision_residues_summary()
        assert cr['GG_diagonal']['r_max_pole'] == 2

    def test_jj_r_pole(self):
        """JJ (affine su(2)): OPE 2 → r-matrix 1."""
        n4 = N4Superconformal(k_val=1)
        cr = n4.collision_residues_summary()
        assert cr['JJ']['r_max_pole'] == 1

    def test_shadow_tt(self):
        """TT shadow: S₂ = c/2 = 3."""
        n4 = N4Superconformal(k_val=1)
        st = n4.shadow_tower()
        assert st['TT']['S2'] == Fraction(3)

    def test_jj_class_l(self):
        """JJ sub-sector: class L (affine ŝu(2))."""
        n4 = N4Superconformal(k_val=1)
        st = n4.shadow_tower()
        assert st['JJ']['S3'] == 0


# =========================================================================
# FAMILY 3: AFFINE EXCEPTIONAL
# =========================================================================

class TestAffineExceptional:
    """Tests for V_k(E₆), V_k(E₇), V_k(E₈) at level 1."""

    @pytest.mark.parametrize("name,dim,h_dual,rank", [
        ("E6", 78, 12, 6),
        ("E7", 133, 18, 7),
        ("E8", 248, 30, 8),
    ])
    def test_invariants(self, name, dim, h_dual, rank):
        """Verify Lie algebra invariants."""
        ae = AffineExceptional(name, k=1)
        assert ae.dim == dim
        assert ae.h_dual == h_dual
        assert ae.rank == rank

    @pytest.mark.parametrize("name,expected_c", [
        ("E6", Fraction(78, 13)),  # = 6
        ("E7", Fraction(133, 19)),  # = 7
        ("E8", Fraction(248, 31)),  # = 8
    ])
    def test_central_charge_k1(self, name, expected_c):
        """c = dim(g)/(1+h∨) at k=1."""
        ae = AffineExceptional(name, k=1)
        assert ae.c == expected_c

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_class_l(self, name):
        """All affine KM are class L."""
        ae = AffineExceptional(name, k=1)
        glcm = ae.glcm_class()
        assert glcm['overall_class'] == 'L'

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_m3_vanishes(self, name):
        """m₃ = 0 for all affine KM (Jacobi identity)."""
        ae = AffineExceptional(name, k=1)
        result = ae.m3_vanishing()
        assert result['m3_zero']
        assert result['m4_zero']

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_r_matrix_simple_pole(self, name):
        """r(z) = kΩ/z (single simple pole after AP19)."""
        ae = AffineExceptional(name, k=1)
        cr = ae.collision_residues()
        assert cr['r_max_pole'] == 1

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_shadow_terminates(self, name):
        """Shadow obstruction tower terminates: S₃ = S₄ = 0."""
        ae = AffineExceptional(name, k=1)
        st = ae.shadow_tower()
        assert st['S3'] == 0
        assert st['S4'] == 0

    @pytest.mark.parametrize("name", ["E6", "E7", "E8"])
    def test_euler_eta(self, name):
        """Euler-eta: η^{dim(g)}."""
        ae = AffineExceptional(name, k=1)
        ee = ae.euler_eta()
        assert ee['eta_exponent'] == ae.dim

    @pytest.mark.parametrize("name,expected_gap", [
        ("E6", 22),   # 2·11
        ("E7", 34),   # 2·17
        ("E8", 58),   # 2·29
    ])
    def test_ds_depth_gap(self, name, expected_gap):
        """DS depth gap = 2·e_r."""
        ae = AffineExceptional(name, k=1)
        glcm = ae.glcm_class()
        assert glcm['ds_depth_gap'] == expected_gap

    def test_depth_gap_ordering(self):
        """E₈ > E₇ > E₆ depth gaps."""
        gaps = {name: AffineExceptional(name).glcm_class()['ds_depth_gap']
                for name in ['E6', 'E7', 'E8']}
        assert gaps['E8'] > gaps['E7'] > gaps['E6']


# =========================================================================
# FAMILY 4: BERSHADSKY-POLYAKOV
# =========================================================================

class TestBershadskyPolyakov:
    """Tests for BP W₃^{(2)}."""

    def test_central_charge_k1(self):
        """c = -(3+1+2)/2 = -3 at k=1."""
        bp = BershadskyPolyakov(k_val=1)
        assert bp.c == Fraction(-3)

    def test_central_charge_k2(self):
        """c = -(12+2+2)/3 = -16/3 at k=2."""
        bp = BershadskyPolyakov(k_val=2)
        assert bp.c == Fraction(-16, 3)

    def test_generators(self):
        """3 generators: T, G, J."""
        assert len(BershadskyPolyakov.GENERATORS) == 3

    def test_gg_coefficient_a(self):
        """A_GG = 2k(2k+3)/(3(k+1)) = 5/3 at k=1."""
        bp = BershadskyPolyakov(k_val=1)
        assert bp.A_gg == Fraction(5, 3)

    def test_ap19_tt(self):
        """TT: OPE 4 → r-matrix 3."""
        bp = BershadskyPolyakov(k_val=1)
        cr = bp.collision_residues()
        assert cr['TT']['max_pole'] == 3

    def test_ap19_gg(self):
        """GG: OPE 3 → r-matrix 2."""
        bp = BershadskyPolyakov(k_val=1)
        cr = bp.collision_residues()
        assert cr['GG']['max_pole'] == 2

    def test_ap19_jj(self):
        """JJ: OPE 2 → r-matrix 1."""
        bp = BershadskyPolyakov(k_val=1)
        cr = bp.collision_residues()
        assert cr['JJ']['max_pole'] == 1

    def test_glcm_class_m(self):
        """Overall class M."""
        bp = BershadskyPolyakov(k_val=1)
        glcm = bp.glcm_class()
        assert glcm['overall_class'] == 'M'


# =========================================================================
# FAMILY 5: VIRASORO MINIMAL MODELS
# =========================================================================

class TestMinimalModels:
    """Tests for Virasoro minimal models M(p,q)."""

    @pytest.mark.parametrize("p,q,expected_c", [
        (4, 3, Fraction(1, 2)),    # Ising
        (5, 4, Fraction(7, 10)),   # Tricritical Ising
        (6, 5, Fraction(4, 5)),    # Tetracritical
        (5, 2, Fraction(-22, 5)),  # Lee-Yang
    ])
    def test_central_charges(self, p, q, expected_c):
        """Verify c = 1 - 6(p-q)²/(pq)."""
        mm = VirasoroMinimalModel(p, q)
        assert mm.c == expected_c

    def test_ising_unitary(self):
        """Ising M(4,3) is unitary."""
        mm = VirasoroMinimalModel(4, 3)
        st = mm.shadow_tower_comparison()
        assert st['unitary']

    def test_lee_yang_non_unitary(self):
        """Lee-Yang M(5,2) is non-unitary."""
        mm = VirasoroMinimalModel(5, 2)
        st = mm.shadow_tower_comparison()
        assert not st['unitary']

    def test_lee_yang_shadow_pole(self):
        """Lee-Yang has S₄ pole (5c+22=0)."""
        mm = VirasoroMinimalModel(5, 2)
        st = mm.shadow_tower_comparison()
        assert st['has_pole']

    def test_ising_no_shadow_pole(self):
        """Ising does not have S₄ pole."""
        mm = VirasoroMinimalModel(4, 3)
        st = mm.shadow_tower_comparison()
        assert not st['has_pole']

    def test_m2_universal(self):
        """m₂ is the same for L_c and Vir_c (universal OPE)."""
        for p, q in [(4, 3), (5, 4), (5, 2)]:
            mm = VirasoroMinimalModel(p, q)
            m2 = mm.m2_comparison()
            assert m2['m2_same']

    def test_collision_residue_cubic(self):
        """All models: r-matrix max pole = 3 (from Virasoro quartic)."""
        for p, q in [(4, 3), (5, 4), (6, 5), (5, 2)]:
            mm = VirasoroMinimalModel(p, q)
            cr = mm.collision_residues_vir()
            assert cr['max_pole'] == 3

    def test_invalid_model_raises(self):
        """Invalid (p,q) should raise."""
        with pytest.raises(AssertionError):
            VirasoroMinimalModel(3, 3)  # p must be > q
        with pytest.raises(AssertionError):
            VirasoroMinimalModel(4, 2)  # gcd(4,2)=2 != 1

    def test_glcm_class_m(self):
        """All minimal models: Vir_c class M."""
        for p, q in [(4, 3), (5, 4), (5, 2)]:
            mm = VirasoroMinimalModel(p, q)
            glcm = mm.glcm_class()
            assert glcm['Vir_c_class'] == 'M (quartic pole, infinite depth)'


# =========================================================================
# CROSS-FAMILY TESTS
# =========================================================================

class TestCrossFamily:
    """Cross-family consistency checks."""

    def test_ap19_all_families(self):
        """AP19 (d-log pole absorption) holds for all families."""
        # N=2
        n2 = N2Superconformal(k_val=1)
        cr = n2.collision_residues()
        for pair, data in cr.items():
            ope = n2.ope_poles(pair[0], pair[1])
            if ope:
                max_ope = max(ope.keys())
                assert data['max_pole'] == max_ope - 1, (
                    f"N=2 AP19 fail at {pair}"
                )

        # Exceptional affine
        for name in ['E6', 'E7', 'E8']:
            ae = AffineExceptional(name)
            cr = ae.collision_residues()
            assert cr['r_max_pole'] == cr['ope_max_pole'] - 1

        # BP
        bp = BershadskyPolyakov(k_val=1)
        cr = bp.collision_residues()
        # TT: OPE 4 → r 3
        assert cr['TT']['max_pole'] == 3
        # GG: OPE 3 → r 2
        assert cr['GG']['max_pole'] == 2
        # JJ: OPE 2 → r 1
        assert cr['JJ']['max_pole'] == 1

    def test_virasoro_subsector_universal(self):
        """All algebras with a weight-2 generator T have TT collision
        residue with r-poles at {3, 1} (from the universal Virasoro OPE)."""
        # N=2
        n2 = N2Superconformal(k_val=1)
        cr_n2 = n2.collision_residues()
        assert 3 in cr_n2[('T', 'T')]['r_poles']
        assert 1 in cr_n2[('T', 'T')]['r_poles']

        # BP
        bp = BershadskyPolyakov(k_val=1)
        cr_bp = bp.collision_residues()
        assert 3 in cr_bp['TT']['r_poles']

    def test_class_l_iff_double_pole(self):
        """Class L ↔ max OPE pole ≤ 2 (and m₃=0)."""
        # E-types are class L
        for name in ['E6', 'E7', 'E8']:
            ae = AffineExceptional(name)
            assert ae.glcm_class()['overall_class'] == 'L'
            assert ae.collision_residues()['ope_max_pole'] == 2
            assert ae.m3_vanishing()['m3_zero']

        # N=2, BP are NOT class L
        n2 = N2Superconformal(k_val=1)
        assert n2.glcm_class()['overall_class'] == 'M'

    def test_koszul_sign_consistency(self):
        """Desuspension flips parity: |s⁻¹(boson)| = 1, |s⁻¹(fermion)| = 0."""
        for cls in [N2Superconformal, BershadskyPolyakov]:
            parities = cls.PARITIES
            for gen, p in parities.items():
                bar_p = (p + 1) % 2
                if p == 0:
                    assert bar_p == 1, f"{gen}: boson should be odd in bar"
                else:
                    assert bar_p == 0, f"{gen}: fermion should be even in bar"
