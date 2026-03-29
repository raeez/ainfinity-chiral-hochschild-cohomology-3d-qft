"""Tests for cross-volume deep bridge: Laplace, signs, shadow-boundary, boundary-linear.

Verifies:
1. Laplace bridge for all standard families
2. Roundtrip: bracket -> OPE -> bracket
3. Sign convention compliance (LV/Koszul at arities 2,3,4)
4. Shadow-boundary kappa comparison
5. Kappa complementarity
6. Shadow depth comparison across volumes
7. Boundary-linear example W(x,y) = xy
8. Koszul dual data cross-volume comparison
9. Full bridge summary

Each test performs ACTUAL computation or consistency check.

Ground truth:
  Vol I: cross_volume_bridge.py, e1_shadow_tower.py, CLAUDE.md
  Vol II: laplace_bridge.py, sc_bar_cobar_engine.py
  CLAUDE.md: Critical Pitfalls AP1-AP13 (esp. AP8, AP9, AP10)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from sympy import Symbol, Rational, simplify, expand, S, symbols, factorial

from lib.cross_volume_deep_bridge import (
    # Laplace bridge
    laplace_transform_bracket,
    inverse_laplace_bracket,
    # Family data
    heisenberg_data,
    affine_sl2_data,
    virasoro_data,
    betagamma_data,
    bc_ghost_data,
    w3_data,
    lattice_vz_data,
    ALL_FAMILIES,
    # Laplace verifications
    verify_laplace_heisenberg,
    verify_laplace_virasoro,
    verify_laplace_affine_diagonal,
    verify_laplace_affine_offdiag,
    verify_laplace_betagamma,
    verify_laplace_w3_WW,
    verify_laplace_roundtrip,
    run_all_laplace_verifications,
    # Sign conventions
    ainfty_sign_m1_squared,
    ainfty_sign_m1_m2,
    ainfty_sign_m2_m2_m3,
    ainfty_sign_arity4,
    verify_sign_convention_families,
    # Shadow-boundary comparison
    shadow_boundary_kappa_table,
    shadow_boundary_complementarity_table,
    shadow_depth_comparison,
    # Boundary-linear
    boundary_linear_superpotential,
    verify_boundary_linear_dsquared,
    verify_boundary_linear_koszul,
    # Koszul dual
    koszul_dual_comparison_table,
    # Full bridge
    run_full_cross_volume_bridge,
)


k = Symbol('k')
c = Symbol('c')
z = Symbol('z')


# =========================================================================
# 1. Laplace bridge for standard families
# =========================================================================

class TestLaplaceHeisenberg:
    """Heisenberg: {J_lam J} = k*lam -> r(z) = k/z^2."""

    def test_laplace_match(self):
        result = verify_laplace_heisenberg()
        assert result['match'] is True

    def test_direct_computation(self):
        """Direct: Laplace of lam^1 = 1!/z^2."""
        bracket = {1: k}
        r = laplace_transform_bracket(bracket)
        assert simplify(r - k / z**2) == 0

    def test_ope_coefficient(self):
        """OPE: J(z)J(w) ~ k/(z-w)^2 has pole order 2 coeff k."""
        bracket = {1: k}
        r = laplace_transform_bracket(bracket)
        # r(z) = k/z^2 => OPE is k/(z-w)^2
        assert simplify(r * z**2 - k) == 0


class TestLaplaceVirasoro:
    """Virasoro: {T_lam T} = dT + 2T*lam + (c/12)*lam^3."""

    def test_laplace_match(self):
        result = verify_laplace_virasoro()
        assert result['match'] is True

    def test_c_over_2_fourth_pole(self):
        """The lambda^3 term gives c/2 at z^{-4}."""
        bracket = {3: c / 12}
        r = laplace_transform_bracket(bracket)
        # 3! = 6, so (c/12)*6/z^4 = c/(2z^4)
        expected = c / (2 * z**4)
        assert simplify(r - expected) == 0

    def test_derivative_term(self):
        """The lambda^0 = dT term gives dT/z."""
        dT = Symbol('dT')
        bracket = {0: dT}
        r = laplace_transform_bracket(bracket)
        assert simplify(r - dT / z) == 0


class TestLaplaceAffine:
    """Affine sl_2: diagonal and off-diagonal."""

    def test_diagonal_match(self):
        result = verify_laplace_affine_diagonal()
        assert result['match'] is True

    def test_offdiag_match(self):
        result = verify_laplace_affine_offdiag()
        assert result['match'] is True

    def test_full_sl2_rmatrix_structure(self):
        """Full r-matrix: r^{ab}(z) = eps^{abc}J^c/z + k*delta^{ab}/z^2."""
        J3 = Symbol('J3')
        # Off-diagonal: {J^1_lam J^2} = J^3 -> J^3/z
        bracket_offdiag = {0: J3}
        r_offdiag = laplace_transform_bracket(bracket_offdiag)
        assert simplify(r_offdiag - J3 / z) == 0

        # Diagonal: {J^1_lam J^1} = k*lam -> k/z^2
        bracket_diag = {1: k}
        r_diag = laplace_transform_bracket(bracket_diag)
        assert simplify(r_diag - k / z**2) == 0


class TestLaplaceBetaGamma:
    """Beta-gamma: {beta_lam gamma} = 1 -> r(z) = 1/z."""

    def test_laplace_match(self):
        result = verify_laplace_betagamma()
        assert result['match'] is True

    def test_simple_pole(self):
        """Constant bracket (lam^0) gives simple pole 1/z."""
        bracket = {0: S.One}
        r = laplace_transform_bracket(bracket)
        assert simplify(r - S.One / z) == 0


class TestLaplaceW3:
    """W_3: leading WW bracket gives c/(3z^6)."""

    def test_ww_leading(self):
        result = verify_laplace_w3_WW()
        assert result['match'] is True

    def test_fifth_power_laplace(self):
        """Laplace of lam^5 = 5!/z^6 = 120/z^6."""
        bracket = {5: S.One}
        r = laplace_transform_bracket(bracket)
        assert simplify(r - 120 / z**6) == 0


class TestLaplaceRoundtrip:
    """Roundtrip: bracket -> OPE -> bracket."""

    def test_roundtrip_constant(self):
        assert verify_laplace_roundtrip({0: S.One}) is True

    def test_roundtrip_linear(self):
        assert verify_laplace_roundtrip({1: k}) is True

    def test_roundtrip_cubic(self):
        assert verify_laplace_roundtrip({3: c / 12}) is True

    def test_roundtrip_mixed(self):
        T = Symbol('T')
        dT = Symbol('dT')
        assert verify_laplace_roundtrip({0: dT, 1: 2 * T, 3: c / 12}) is True

    def test_roundtrip_quintic(self):
        assert verify_laplace_roundtrip({5: c / 360}) is True


class TestInverseLaplace:
    """Inverse Laplace: OPE -> bracket."""

    def test_inverse_simple(self):
        """OPE coeff c_0 = 1 -> bracket {0: 1/0! = 1}."""
        result = inverse_laplace_bracket({0: S.One})
        assert result[0] == S.One

    def test_inverse_double_pole(self):
        """OPE coeff c_1 = k -> bracket {1: k/1! = k}."""
        result = inverse_laplace_bracket({1: k})
        assert simplify(result[1] - k) == 0

    def test_inverse_fourth_pole(self):
        """OPE coeff c_3 = c/2 -> bracket {3: (c/2)/3! = c/12}."""
        result = inverse_laplace_bracket({3: c / 2})
        assert simplify(result[3] - c / 12) == 0


# =========================================================================
# 2. Sign convention checks
# =========================================================================

class TestSignConventions:
    """LV/Koszul sign conventions for A_infinity."""

    def test_m1_squared_uncurved(self):
        """m_1^2 = 0 in the uncurved case."""
        data = ainfty_sign_m1_squared()
        assert 'm_1^2 = 0' in data['relation']

    def test_m1_m2_koszul_sign(self):
        """Koszul sign (-1)^{|a|} in m_1 m_2 relation."""
        data = ainfty_sign_m1_m2()
        assert '(-1)^{|a|}' in data['sign']

    def test_m2_m3_relation(self):
        """Arity-3 relation has 6 terms (K_4 pentagon)."""
        data = ainfty_sign_m2_m2_m3()
        assert data['num_terms'] == 6

    def test_arity4_relation(self):
        """Arity-4 relation has 14 terms (K_5 associahedron)."""
        data = ainfty_sign_arity4()
        assert data['num_terms'] == 14

    def test_all_families_uncurved_g0(self):
        """All families are uncurved at genus 0."""
        results = verify_sign_convention_families()
        for name, data in results.items():
            assert data['genus_0_uncurved'] is True

    def test_m1_squared_zero_g0(self):
        """d^2 = 0 at genus 0 for all families."""
        results = verify_sign_convention_families()
        for name, data in results.items():
            assert data['m1_squared_zero_g0'] is True

    def test_curved_families_at_g1(self):
        """Families with nonzero kappa are curved at genus 1."""
        results = verify_sign_convention_families()
        # All families except possibly some edge cases have kappa != 0
        for name, data in results.items():
            if data['kappa'] != S.Zero:
                assert data['m1_squared_curved_g1'] is True


# =========================================================================
# 3. Shadow-boundary kappa comparison
# =========================================================================

class TestShadowBoundaryKappa:
    """Vol I kappa vs Vol II boundary curvature."""

    def test_all_families_match(self):
        """kappa_vol1 = kappa_vol2 for all families."""
        table = shadow_boundary_kappa_table()
        for name, data in table.items():
            assert data['match'] is True, f"kappa mismatch for {name}"

    def test_heisenberg_kappa(self):
        table = shadow_boundary_kappa_table()
        assert simplify(table['Heisenberg H_k']['kappa_vol1'] - k) == 0

    def test_virasoro_kappa(self):
        table = shadow_boundary_kappa_table()
        assert simplify(table['Virasoro Vir_c']['kappa_vol1'] - c / 2) == 0

    def test_sl2_kappa(self):
        table = shadow_boundary_kappa_table()
        expected = Rational(3, 4) * (k + 2)
        assert simplify(table['Affine sl_2 at level k']['kappa_vol1'] - expected) == 0


# =========================================================================
# 4. Kappa complementarity
# =========================================================================

class TestKappaComplementarity:
    """kappa(A) + kappa(A!) is constant (independent of parameters)."""

    def test_all_constant(self):
        """Sum is parameter-independent for all families."""
        table = shadow_boundary_complementarity_table()
        for name, data in table.items():
            assert data['is_constant'] is True, \
                f"kappa sum not constant for {name}: {data['sum']}"

    def test_heisenberg_sum_zero(self):
        """Heisenberg: kappa + kappa' = 0."""
        table = shadow_boundary_complementarity_table()
        assert table['Heisenberg H_k']['sum'] == 0

    def test_virasoro_sum_13(self):
        """Virasoro: kappa + kappa' = 13."""
        table = shadow_boundary_complementarity_table()
        assert table['Virasoro Vir_c']['sum'] == 13

    def test_sl2_sum_zero(self):
        """Affine sl_2: kappa + kappa' = 0."""
        table = shadow_boundary_complementarity_table()
        assert table['Affine sl_2 at level k']['sum'] == 0

    def test_w3_sum_constant(self):
        """W_3: kappa + kappa' = 250/3."""
        table = shadow_boundary_complementarity_table()
        w3_sum = table['W_3 at central charge c']['sum']
        assert simplify(w3_sum - Rational(250, 3)) == 0

    def test_betagamma_bc_complementary(self):
        """Beta-gamma and bc are Koszul dual with kappa sum = 0."""
        table = shadow_boundary_complementarity_table()
        bg_sum = table['Beta-gamma']['sum']
        assert bg_sum == 0


# =========================================================================
# 5. Shadow depth comparison
# =========================================================================

class TestShadowDepthComparison:
    """Shadow depth classification consistent across volumes."""

    def test_all_match(self):
        table = shadow_depth_comparison()
        for name, data in table.items():
            assert data['match'] is True

    def test_heisenberg_gaussian(self):
        table = shadow_depth_comparison()
        assert table['Heisenberg H_k']['vol1_class'] == 'G'
        assert table['Heisenberg H_k']['vol1_depth'] == 2

    def test_sl2_lie(self):
        table = shadow_depth_comparison()
        assert table['Affine sl_2 at level k']['vol1_class'] == 'L'
        assert table['Affine sl_2 at level k']['vol1_depth'] == 3

    def test_virasoro_mixed(self):
        table = shadow_depth_comparison()
        assert table['Virasoro Vir_c']['vol1_class'] == 'M'
        assert table['Virasoro Vir_c']['vol1_depth'] == float('inf')

    def test_betagamma_contact(self):
        table = shadow_depth_comparison()
        assert table['Beta-gamma']['vol1_class'] == 'C'
        assert table['Beta-gamma']['vol1_depth'] == 4

    def test_lattice_gaussian(self):
        table = shadow_depth_comparison()
        assert table['Lattice V_Z']['vol1_class'] == 'G'
        assert table['Lattice V_Z']['vol1_depth'] == 2


# =========================================================================
# 6. Boundary-linear example: W(x,y) = xy
# =========================================================================

class TestBoundaryLinear:
    """Boundary-linear superpotential W = xy."""

    def test_is_quadratic(self):
        data = boundary_linear_superpotential()
        assert data['is_quadratic'] is True

    def test_m3_zero(self):
        """Quadratic W => m_k = 0 for k >= 3."""
        data = boundary_linear_superpotential()
        assert data['m3_zero'] is True
        assert data['m4_zero'] is True

    def test_bar_strict_dg(self):
        """B^{intr} is a strict dg algebra."""
        data = boundary_linear_superpotential()
        assert data['bar_is_strict_dg'] is True

    def test_critical_point_isolated(self):
        """Critical locus is the origin (isolated)."""
        data = boundary_linear_superpotential()
        assert data['critical_dim'] == 0

    def test_jacobi_ring_dim_1(self):
        """Jacobi ring C[x,y]/(y,x) = C has dim 1."""
        data = boundary_linear_superpotential()
        assert data['jacobi_dim'] == 1

    def test_d_squared_zero(self):
        result = verify_boundary_linear_dsquared()
        assert result['d_squared_zero'] is True

    def test_koszulness(self):
        result = verify_boundary_linear_koszul()
        assert result['is_koszul'] is True
        assert result['is_formal'] is True
        assert result['shadow_class'] == 'G'


# =========================================================================
# 7. Koszul dual comparison table
# =========================================================================

class TestKoszulDualTable:
    """Cross-volume Koszul dual data."""

    def test_heisenberg_not_self_dual(self):
        """Heisenberg is NOT self-dual (critical pitfall)."""
        table = koszul_dual_comparison_table()
        assert table['Heisenberg']['self_dual'] is False

    def test_virasoro_self_dual_at_13(self):
        """Virasoro self-dual at c=13, NOT c=26."""
        table = koszul_dual_comparison_table()
        assert table['Virasoro']['self_dual_at'] == 13
        assert table['Virasoro']['NOT_self_dual_at_26'] is True

    def test_sl2_critical_warning(self):
        """sl_2 self-dual at critical level k=-2 (Sugawara undefined)."""
        table = koszul_dual_comparison_table()
        assert table['Affine sl_2']['self_dual_at'] == -2
        assert 'Sugawara' in table['Affine sl_2']['critical_level_warning']

    def test_w3_dual_at_100_minus_c(self):
        """W_3: c + c' = 100."""
        table = koszul_dual_comparison_table()
        assert table['W_3']['c_sum'] == 100

    def test_virasoro_kappa_sum(self):
        """Virasoro: kappa_sum = 13."""
        table = koszul_dual_comparison_table()
        assert table['Virasoro']['kappa_sum'] == 13


# =========================================================================
# 8. Full bridge verification
# =========================================================================

class TestFullBridge:
    """Run the complete cross-volume bridge."""

    def test_all_laplace_pass(self):
        """All Laplace verifications pass."""
        results = run_all_laplace_verifications()
        for name, data in results.items():
            assert data['match'] is True, f"Laplace failed for {name}"

    def test_full_bridge_runs(self):
        """Full bridge runs without error."""
        summary = run_full_cross_volume_bridge()
        assert 'summary' in summary

    def test_full_bridge_laplace_count(self):
        """All Laplace checks pass in full bridge."""
        summary = run_full_cross_volume_bridge()
        s = summary['summary']
        assert s['laplace_passed'] == s['laplace_total']

    def test_full_bridge_kappa_count(self):
        """All kappa checks pass in full bridge."""
        summary = run_full_cross_volume_bridge()
        s = summary['summary']
        assert s['kappa_passed'] == s['kappa_total']

    def test_full_bridge_depth_count(self):
        """All depth checks pass in full bridge."""
        summary = run_full_cross_volume_bridge()
        s = summary['summary']
        assert s['depth_passed'] == s['depth_total']


# =========================================================================
# 9. Family data consistency
# =========================================================================

class TestFamilyDataConsistency:
    """Cross-checks on family data (AP10 compliance: don't just hardcode)."""

    def test_seven_families(self):
        """There are 7 standard families."""
        assert len(ALL_FAMILIES) == 7

    def test_heisenberg_central_charge_1(self):
        """Heisenberg always has c = 1."""
        fam = heisenberg_data()
        assert fam.central_charge == S.One

    def test_virasoro_dual_cc(self):
        """Vir_c^! has central charge 26-c."""
        fam = virasoro_data()
        assert simplify(fam.dual_central_charge - (26 - c)) == 0

    def test_sl2_kappa_formula(self):
        """kappa(sl_2) = dim(g)*(k+h^v)/(2*h^v) = 3*(k+2)/4."""
        fam = affine_sl2_data()
        expected = Rational(3, 4) * (k + 2)
        assert simplify(fam.kappa - expected) == 0

    def test_w3_kappa_formula(self):
        """kappa(W_3) = 5c/6."""
        fam = w3_data()
        assert simplify(fam.kappa - 5 * c / 6) == 0

    def test_bc_kappa(self):
        """bc: kappa = c/2 = -26/2 = -13."""
        fam = bc_ghost_data()
        assert fam.kappa == S(-13)

    def test_kappa_complementarity_all_families(self):
        """kappa + kappa' is parameter-free for all families."""
        for fam_fn in ALL_FAMILIES:
            fam = fam_fn()
            kappa_sum = simplify(expand(fam.kappa + fam.dual_kappa))
            # Must have no free symbols from level parameters
            assert len(kappa_sum.free_symbols) == 0, \
                f"{fam.name}: kappa sum = {kappa_sum} has free symbols"
