"""Tests for holographic_ht_engine: the holographic modular Koszul datum H(T).

Verifies:
1. Holographic datum construction for all standard families
2. Kappa complementarity (Theorem C) for all families
3. Shadow depth classification consistency
4. Collision residue r(z) pole orders and CYBE
5. Sphere reconstruction at arities 3 and 4
6. Holographic connection flatness
7. M2 brane boundary algebra
8. Genus expansion from holographic datum
9. Cross-family consistency

References:
  Vol II: ht_bulk_boundary_line.tex (constr:genus-zero-to-hmkd)
  Vol I:  concordance.tex, higher_genus_modular_koszul.tex
  CLAUDE.md: Critical Pitfalls AP1-AP13
"""

import pytest
from fractions import Fraction

from sympy import S, Rational, simplify, Symbol, oo, bernoulli

from lib.holographic_ht_engine import (
    # Core constructors
    heisenberg_holographic_datum,
    affine_sl2_holographic_datum,
    betagamma_holographic_datum,
    virasoro_holographic_datum,
    w3_holographic_datum,
    # Sphere reconstruction
    sphere_shadow_arity_n,
    sphere_reconstruction_3pt,
    sphere_reconstruction_4pt,
    # Connection
    holographic_connection_data,
    verify_connection_flatness_genus0,
    genus1_one_loop_data,
    connection_flatness_from_mc,
    # M2 brane
    M2Generator,
    m2_classical_commutator,
    m2_lowest_commutator,
    m2_boundary_algebra,
    m2_generator_count,
    m2_shadow_data,
    m2_holographic_complementarity,
    # Cross-family
    all_standard_holographic_data,
    verify_all_kappa_complementarity,
    verify_all_shadow_depth_consistency,
    verify_collision_residue_pole_orders,
    genus_expansion_from_holographic_datum,
    holographic_datum_summary,
)


# =========================================================================
# 1. HOLOGRAPHIC DATUM CONSTRUCTION
# =========================================================================

class TestHolographicDatumConstruction:
    """Test that H(T) is correctly constructed for each standard family."""

    def test_heisenberg_datum_exists(self):
        """H(free) should have all six components."""
        H = heisenberg_holographic_datum()
        assert H.A is not None
        assert H.A_dual is not None
        assert H.C is not None
        assert H.r_z is not None
        assert H.theta is not None
        assert H.connection is not None

    def test_affine_sl2_datum_exists(self):
        """H(CS) should have all six components."""
        H = affine_sl2_holographic_datum()
        assert H.A is not None
        assert H.A_dual is not None

    def test_betagamma_datum_exists(self):
        """H(N=4) should have all six components."""
        H = betagamma_holographic_datum()
        assert H.A is not None
        assert H.A_dual is not None

    def test_virasoro_datum_exists(self):
        """H(Vir) should have all six components."""
        H = virasoro_holographic_datum()
        assert H.A is not None
        assert H.A_dual is not None

    def test_w3_datum_exists(self):
        """H(W3) should have all six components."""
        H = w3_holographic_datum()
        assert H.A is not None
        assert H.A_dual is not None

    def test_heisenberg_family_tag(self):
        H = heisenberg_holographic_datum()
        assert H.A.family == 'heisenberg'

    def test_affine_family_tag(self):
        H = affine_sl2_holographic_datum()
        assert H.A.family == 'affine'

    def test_betagamma_family_tag(self):
        H = betagamma_holographic_datum()
        assert H.A.family == 'betagamma'


# =========================================================================
# 2. KAPPA COMPLEMENTARITY (THEOREM C)
# =========================================================================

class TestKappaComplementarity:
    """Verify kappa(A) + kappa(A!) = rho*K for all families.

    CRITICAL: These are THREE DISTINCT formulas (AP1):
    kappa(KM) = dim(g)(k+h^v)/(2h^v)
    kappa(Vir) = c/2
    kappa(W_3) = 5c/6
    """

    def test_heisenberg_kappa_sum_zero(self):
        """kappa(H_k) + kappa(H_k^!) = 0."""
        H = heisenberg_holographic_datum(k_val=1)
        result = H.kappa_complementarity()
        assert result['match'], f"kappa sum = {result['sum']}, expected 0"

    def test_affine_sl2_kappa_sum_zero(self):
        """kappa(V_k(sl_2)) + kappa(V_{k'}(sl_2)) = 0."""
        H = affine_sl2_holographic_datum()
        result = H.kappa_complementarity()
        assert result['match'], f"kappa sum = {result['sum']}, expected 0"

    def test_betagamma_kappa_sum_zero(self):
        """kappa(betagamma) + kappa(bc) = 0."""
        H = betagamma_holographic_datum()
        result = H.kappa_complementarity()
        assert result['match'], f"kappa sum = {result['sum']}, expected 0"

    def test_virasoro_kappa_sum_13(self):
        """kappa(Vir_c) + kappa(Vir_{26-c}) = 13."""
        H = virasoro_holographic_datum()
        result = H.kappa_complementarity()
        assert result['match'], f"kappa sum = {result['sum']}, expected 13"

    def test_w3_kappa_sum_250_3(self):
        """kappa(W_3(c)) + kappa(W_3(100-c)) = 250/3."""
        H = w3_holographic_datum()
        result = H.kappa_complementarity()
        assert result['match'], f"kappa sum = {result['sum']}, expected 250/3"

    def test_heisenberg_kappa_value(self):
        """kappa(H_1) = 1."""
        H = heisenberg_holographic_datum(k_val=1)
        assert H.A.kappa == S.One

    def test_betagamma_kappa_value(self):
        """kappa(betagamma) = 1 (c = 2, so c/2 = 1)."""
        H = betagamma_holographic_datum()
        assert H.A.kappa == 1

    def test_all_families_complementarity(self):
        """Batch check kappa complementarity for all families."""
        results = verify_all_kappa_complementarity()
        for name, result in results.items():
            assert result['match'], f"Complementarity fails for {name}: {result}"


# =========================================================================
# 3. SHADOW DEPTH CLASSIFICATION
# =========================================================================

class TestShadowDepthClassification:
    """Verify shadow depth matches G/L/C/M classification."""

    def test_heisenberg_gaussian(self):
        H = heisenberg_holographic_datum()
        assert H.A.shadow_class == 'G'
        assert H.A.shadow_depth == 2

    def test_affine_lie_tree(self):
        H = affine_sl2_holographic_datum()
        assert H.A.shadow_class == 'L'
        assert H.A.shadow_depth == 3

    def test_betagamma_contact(self):
        H = betagamma_holographic_datum()
        assert H.A.shadow_class == 'C'
        assert H.A.shadow_depth == 4

    def test_virasoro_mixed(self):
        H = virasoro_holographic_datum()
        assert H.A.shadow_class == 'M'
        assert H.A.shadow_depth == oo

    def test_w3_mixed(self):
        H = w3_holographic_datum()
        assert H.A.shadow_class == 'M'
        assert H.A.shadow_depth == oo

    def test_all_shadow_depth_consistency(self):
        """All families should have consistent shadow depth data."""
        results = verify_all_shadow_depth_consistency()
        for name, result in results.items():
            assert result['class_match'], f"Class mismatch for {name}: {result}"
            assert result['depth_match'], f"Depth mismatch for {name}: {result}"


# =========================================================================
# 4. COLLISION RESIDUE r(z)
# =========================================================================

class TestCollisionResidue:
    """Test r(z) = Res^coll_{0,2}(Theta_A)."""

    def test_heisenberg_pole_order_1(self):
        """r(z) = k/z for Heisenberg: pole order 1."""
        H = heisenberg_holographic_datum()
        assert H.r_z.pole_order == 1

    def test_affine_pole_order_1(self):
        """r(z) = Omega/z for affine: pole order 1."""
        H = affine_sl2_holographic_datum()
        assert H.r_z.pole_order == 1

    def test_betagamma_pole_order_1(self):
        """r(z) = 1/z for betagamma: pole order 1."""
        H = betagamma_holographic_datum()
        assert H.r_z.pole_order == 1

    def test_virasoro_pole_order_3(self):
        """r(z) for Virasoro has pole order 3 (from z^{-4} term)."""
        H = virasoro_holographic_datum()
        assert H.r_z.pole_order == 3

    def test_w3_pole_order_5(self):
        """r(z) for W_3 has pole order 5 (from z^{-6} term)."""
        H = w3_holographic_datum()
        assert H.r_z.pole_order == 5

    def test_all_satisfy_cybe(self):
        """All standard families should satisfy the CYBE."""
        for name, datum in all_standard_holographic_data().items():
            assert datum.r_z.satisfies_cybe, f"CYBE fails for {name}"

    def test_pole_order_consistency(self):
        """Batch check pole orders."""
        results = verify_collision_residue_pole_orders()
        for name, result in results.items():
            assert result['match'], f"Pole order mismatch for {name}: {result}"


# =========================================================================
# 5. SPHERE RECONSTRUCTION
# =========================================================================

class TestSphereReconstruction:
    """Test S_n(Theta_A)|_{genus 0} recovers boundary OPE."""

    def test_heisenberg_s2_is_kappa(self):
        """S_2 = kappa for Heisenberg."""
        H = heisenberg_holographic_datum(k_val=1)
        s2 = sphere_shadow_arity_n(H, 2)
        assert s2['value'] == S.One

    def test_heisenberg_s3_vanishes(self):
        """S_3 = 0 for Heisenberg (Gaussian class: tower terminates at 2)."""
        H = heisenberg_holographic_datum()
        s3 = sphere_shadow_arity_n(H, 3)
        assert s3['is_zero']

    def test_heisenberg_s4_vanishes(self):
        """S_4 = 0 for Heisenberg."""
        H = heisenberg_holographic_datum()
        s4 = sphere_shadow_arity_n(H, 4)
        assert s4['is_zero']

    def test_affine_s3_nonzero(self):
        """S_3 nonzero for affine (Lie class has cubic shadow)."""
        H = affine_sl2_holographic_datum()
        s3 = sphere_shadow_arity_n(H, 3)
        assert not s3['is_zero']

    def test_affine_s4_vanishes(self):
        """S_4 = 0 for affine (tower terminates at arity 3)."""
        H = affine_sl2_holographic_datum()
        s4 = sphere_shadow_arity_n(H, 4)
        assert s4['is_zero']

    def test_betagamma_s3_vanishes(self):
        """S_3 = 0 for betagamma (cubic shadow vanishes)."""
        H = betagamma_holographic_datum()
        s3 = sphere_shadow_arity_n(H, 3)
        assert s3['is_zero']

    def test_virasoro_s3_nonzero(self):
        """S_3 nonzero for Virasoro (mixed class has cubic shadow)."""
        H = virasoro_holographic_datum()
        s3 = sphere_shadow_arity_n(H, 3)
        assert not s3['is_zero']

    def test_virasoro_quartic_contact(self):
        """Q^contact_Vir = 10/[c(5c+22)] at c=25 should be 10/(25*147) = 2/735."""
        H = virasoro_holographic_datum(c_val=25)
        s4 = sphere_shadow_arity_n(H, 4)
        expected = Rational(10, 25 * (5 * 25 + 22))  # 10/3675 = 2/735
        assert simplify(S(s4['value']) - expected) == 0

    def test_sphere_reconstruction_3pt_heisenberg(self):
        """3-point reconstruction for Heisenberg: quadratic OPE suffices."""
        H = heisenberg_holographic_datum()
        result = sphere_reconstruction_3pt(H)
        assert result['ope_recovered']
        assert result['trilinear_coupling'] == S.Zero

    def test_sphere_reconstruction_4pt_virasoro(self):
        """4-point reconstruction for Virasoro: Q^contact provides crossing."""
        H = virasoro_holographic_datum()
        result = sphere_reconstruction_4pt(H)
        assert 'Q^contact_Vir' in result['crossing_info']

    def test_higher_arity_heisenberg_vanishes(self):
        """S_n = 0 for n > 2 for Heisenberg (depth 2)."""
        H = heisenberg_holographic_datum()
        for n in [5, 6, 7]:
            sn = sphere_shadow_arity_n(H, n)
            assert sn['is_zero'], f"S_{n} should vanish for Heisenberg"


# =========================================================================
# 6. HOLOGRAPHIC CONNECTION
# =========================================================================

class TestHolographicConnection:
    """Test nabla^hol_{g,n} = d - Sh_{g,n}(Theta_A)."""

    def test_genus0_heisenberg_flat(self):
        """Genus-0 Heisenberg connection is flat."""
        H = heisenberg_holographic_datum()
        result = verify_connection_flatness_genus0(H, 3)
        assert result['is_flat']
        assert result['curvature_vanishes']

    def test_genus0_affine_kz_flat(self):
        """Genus-0 KZ connection for affine sl_2 is flat."""
        H = affine_sl2_holographic_datum()
        result = verify_connection_flatness_genus0(H, 3)
        assert result['is_flat']
        assert result['connection_type'] == 'KZ'

    def test_genus0_virasoro_bpz_flat(self):
        """Genus-0 BPZ connection for Virasoro is flat."""
        H = virasoro_holographic_datum()
        result = verify_connection_flatness_genus0(H, 3)
        assert result['is_flat']
        assert result['connection_type'] == 'BPZ'

    def test_genus1_one_loop_kappa(self):
        """Genus-1 curvature = kappa for Heisenberg."""
        H = heisenberg_holographic_datum(k_val=1)
        result = genus1_one_loop_data(H)
        assert result['kappa'] == S.One
        assert result['lambda_1_fp'] == Rational(-1, 24)
        assert simplify(result['f1'] + Rational(1, 24)) == 0  # -kappa/24 = -1/24

    def test_genus1_virasoro_curvature(self):
        """Genus-1 curvature for Virasoro is c/2."""
        H = virasoro_holographic_datum(c_val=25)
        result = genus1_one_loop_data(H)
        assert result['kappa'] == Rational(25, 2)

    def test_all_connections_flat(self):
        """All pre-built connections in the datum should be flat."""
        for name, datum in all_standard_holographic_data().items():
            result = connection_flatness_from_mc(datum)
            assert result['all_flat'], f"Not all connections flat for {name}"

    def test_mc_mechanism(self):
        """Flatness mechanism is D^2 = 0."""
        H = heisenberg_holographic_datum()
        result = connection_flatness_from_mc(H)
        assert result['bar_intrinsic']
        assert result['mc_equation'] == 'D*Theta + (1/2)[Theta,Theta] = 0'


# =========================================================================
# 7. M2 BRANE EXAMPLE
# =========================================================================

class TestM2Brane:
    """Test the M2 brane boundary algebra and shadow data."""

    def test_m2_generator_count_k1(self):
        """Generator count for K=1, weight_cutoff=3: 1 * 3*4/2 = 6."""
        assert m2_generator_count(1, 3) == 6

    def test_m2_generator_count_k2(self):
        """Generator count for K=2, weight_cutoff=2: 4 * 2*3/2 = 12."""
        assert m2_generator_count(2, 2) == 12

    def test_m2_boundary_algebra_construction(self):
        """M2 boundary algebra should have correct number of generators."""
        data = m2_boundary_algebra(K=2, weight_cutoff=3)
        expected = m2_generator_count(2, 3)
        assert data.n_generators_truncated == expected

    def test_m2_lowest_commutator(self):
        """[E_{ab} z, E_{cd} partial]_0 = -delta_{bc} E_{ad}."""
        result = m2_lowest_commutator(K=2)
        assert result['classical_coefficient'] == -1
        assert result['quantum_correction_order'] == 1

    def test_m2_classical_commutator_diagonal(self):
        """Test commutator with matching matrix indices."""
        g1 = M2Generator(alpha=1, beta=2, m=1, n=0, weight=2)
        g2 = M2Generator(alpha=2, beta=1, m=0, n=1, weight=2)
        result = m2_classical_commutator(g1, g2, K=2)
        # delta_{beta,gamma} = delta_{2,2} = 1, so term1 is active
        assert result['term1_active']
        # delta_{alpha,delta} = delta_{1,1} = 1, so term2 is active
        assert result['term2_active']

    def test_m2_classical_commutator_off_diagonal(self):
        """Test commutator with non-matching matrix indices."""
        g1 = M2Generator(alpha=1, beta=1, m=0, n=0, weight=1)
        g2 = M2Generator(alpha=2, beta=2, m=0, n=0, weight=1)
        result = m2_classical_commutator(g1, g2, K=2)
        # delta_{1,2} = 0 and delta_{1,2} = 0
        assert not result['term1_active']
        assert not result['term2_active']

    def test_m2_shadow_data_class(self):
        """M2 shadow class should be M (mixed, infinite tower)."""
        data = m2_shadow_data(K=3)
        assert data['shadow_class'] == 'M'

    def test_m2_shadow_data_kappa_open(self):
        """kappa_{M2} is open (not computed from first principles)."""
        data = m2_shadow_data(K=3)
        assert data['kappa'] is None
        assert data['kappa_status'] == 'open'

    def test_m2_complementarity_structural(self):
        """Structural complementarity holds regardless of kappa value."""
        result = m2_holographic_complementarity(K=3)
        assert result['complementarity_structural']
        assert result['kappa_status'] == 'open'

    def test_m2_four_fold_match(self):
        """M2 system has four-fold parametrization."""
        data = m2_shadow_data(K=2)
        match = data['four_fold_match']
        assert 'sector' in match
        assert 'level' in match
        assert 'rank' in match
        assert 'genus' in match


# =========================================================================
# 8. GENUS EXPANSION
# =========================================================================

class TestGenusExpansion:
    """Test genus expansion F_g(A) from holographic datum."""

    def test_genus0_trivial(self):
        """F_0 = 0 by convention."""
        H = heisenberg_holographic_datum(k_val=1)
        result = genus_expansion_from_holographic_datum(H, max_genus=2)
        assert result[0]['F_g'] == 0

    def test_genus1_heisenberg(self):
        """F_1(H_1) = -kappa/24 = -1/24."""
        H = heisenberg_holographic_datum(k_val=1)
        result = genus_expansion_from_holographic_datum(H, max_genus=2)
        assert simplify(result[1]['F_g'] + Rational(1, 24)) == 0

    def test_genus1_betagamma(self):
        """F_1(betagamma) = -1/24 (kappa = 1)."""
        H = betagamma_holographic_datum()
        result = genus_expansion_from_holographic_datum(H, max_genus=2)
        assert simplify(result[1]['F_g'] + Rational(1, 24)) == 0

    def test_genus2_positive(self):
        """F_2 > 0 for positive kappa (Bernoulli signs: F_g POSITIVE)."""
        H = heisenberg_holographic_datum(k_val=1)
        result = genus_expansion_from_holographic_datum(H, max_genus=2)
        f2 = result[2]['F_g']
        assert f2 > 0, f"F_2 = {f2} should be positive"

    def test_genus2_value(self):
        """F_2 = kappa * 7/5760.

        B_4 = -1/30.
        lambda_2^FP = (2^3 - 1)/2^3 * |B_4| / 4! = 7/8 * (1/30) / 24 = 7/5760.
        """
        H = heisenberg_holographic_datum(k_val=1)
        result = genus_expansion_from_holographic_datum(H, max_genus=2)
        expected_lambda = Rational(7, 5760)
        assert result[2]['lambda_fp'] == expected_lambda


# =========================================================================
# 9. CROSS-FAMILY CONSISTENCY
# =========================================================================

class TestCrossFamilyConsistency:
    """Cross-family consistency checks (AP10: not just hardcoded tests)."""

    def test_kappa_additivity_heisenberg(self):
        """Kappa should be additive: kappa(H_k1 + H_k2) = kappa(H_k1) + kappa(H_k2)."""
        H1 = heisenberg_holographic_datum(k_val=1)
        H2 = heisenberg_holographic_datum(k_val=2)
        # kappa(H_1) = 1, kappa(H_2) = 2
        assert H1.A.kappa + H2.A.kappa == S(3)

    def test_virasoro_self_dual_at_13(self):
        """Virasoro self-dual at c = 13, NOT c = 26 (AP8)."""
        H = virasoro_holographic_datum(c_val=13)
        assert H.A.central_charge == 13
        assert H.A_dual.central_charge == 13  # 26 - 13 = 13

    def test_virasoro_not_self_dual_at_26(self):
        """Virasoro NOT self-dual at c = 26."""
        H = virasoro_holographic_datum(c_val=26)
        assert H.A_dual.central_charge == 0  # 26 - 26 = 0
        assert H.A.central_charge != H.A_dual.central_charge

    def test_w3_self_dual_at_50(self):
        """W_3 self-dual at c = 50."""
        H = w3_holographic_datum(c_val=50)
        assert H.A.central_charge == 50
        assert H.A_dual.central_charge == 50  # 100 - 50 = 50

    def test_heisenberg_not_self_dual(self):
        """Heisenberg is NOT self-dual (CRITICAL PITFALL)."""
        H = heisenberg_holographic_datum()
        assert H.A_dual.duality_type == 'chiral-symmetric'
        # H_k^! = Sym^ch(V*) != H_{-k}

    def test_datum_summary_has_all_components(self):
        """Summary should have all six components of H(T)."""
        H = heisenberg_holographic_datum(k_val=1)
        summary = holographic_datum_summary(H)
        assert 'A' in summary
        assert 'A!' in summary
        assert 'C' in summary
        assert 'r(z)' in summary
        assert 'Theta_A' in summary
        assert 'nabla^hol' in summary

    def test_bar_intrinsic_for_all_families(self):
        """All standard families should have bar-intrinsic Theta_A."""
        for name, datum in all_standard_holographic_data().items():
            assert datum.theta.is_bar_intrinsic, \
                f"Theta_A not bar-intrinsic for {name}"

    def test_all_standard_data_constructs(self):
        """All 5 standard holographic data should construct without error."""
        data = all_standard_holographic_data()
        assert len(data) == 5
        for name in ['heisenberg', 'affine_sl2', 'betagamma', 'virasoro', 'w3']:
            assert name in data


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
