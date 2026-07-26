r"""Tests for collision residue r-matrix: r(z) = Res^coll_{0,2}(Theta_A).

Verifies the identification r(z) = k * Omega / z for affine g_k
from first principles, including:
- Lie algebra axiom verification
- Casimir tensor computation
- d-log pole absorption (AP19)
- Classical Yang-Baxter equation (CYBE) via IBR
- Cross-family consistency

References:
  Vol I: higher_genus_modular_koszul.tex (thm:mc2-bar-intrinsic)
  Vol II: ht_bulk_boundary_line_core.tex (holographic datum)
  AP19: The bar kernel absorbs a pole
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.collision_residue_rmatrix import (
    make_sl2, make_sl3, make_u1,
    verify_jacobi, verify_killing_invariance, verify_antisymmetry,
    casimir_tensor, casimir_tensor_explicit,
    AffineOPE, collision_residue_rmatrix,
    affine_trace_form_scale, affine_kz_scale, affine_level_prefixed_kz_scale,
    verify_cybe, verify_pole_absorption,
    virasoro_ope_to_rmatrix,
    full_collision_residue_computation,
)


# =========================================================================
# PART 1: Lie algebra axioms
# =========================================================================

class TestLieAlgebraAxioms:
    """Verify Lie algebra data is correct before using it."""

    def test_sl2_jacobi(self):
        """Jacobi identity for sl_2."""
        g = make_sl2()
        assert verify_jacobi(g), "sl_2 fails Jacobi identity"

    def test_sl2_antisymmetry(self):
        """Antisymmetry f^{ab}_c = -f^{ba}_c for sl_2."""
        g = make_sl2()
        assert verify_antisymmetry(g), "sl_2 structure constants not antisymmetric"

    def test_sl2_killing_invariance(self):
        """Ad-invariance of Killing form for sl_2."""
        g = make_sl2()
        assert verify_killing_invariance(g), "sl_2 Killing form not ad-invariant"

    def test_sl2_dimension(self):
        """sl_2 has dimension 3."""
        g = make_sl2()
        assert g.dim == 3
        assert g.rank == 1
        assert g.h_dual == 2

    def test_sl3_jacobi(self):
        """Jacobi identity for sl_3."""
        g = make_sl3()
        assert verify_jacobi(g), "sl_3 fails Jacobi identity"

    def test_sl3_antisymmetry(self):
        """Antisymmetry for sl_3."""
        g = make_sl3()
        assert verify_antisymmetry(g), "sl_3 structure constants not antisymmetric"

    def test_sl3_killing_invariance(self):
        """Ad-invariance of Killing form for sl_3."""
        g = make_sl3()
        assert verify_killing_invariance(g), "sl_3 Killing form not ad-invariant"

    def test_sl3_dimension(self):
        """sl_3 has dimension 8."""
        g = make_sl3()
        assert g.dim == 8
        assert g.rank == 2
        assert g.h_dual == 3

    def test_u1_jacobi(self):
        """Jacobi identity for u(1) (trivial)."""
        g = make_u1()
        assert verify_jacobi(g)

    def test_u1_abelian(self):
        """u(1) has vanishing structure constants."""
        g = make_u1()
        assert np.max(np.abs(g.f)) < 1e-15


# =========================================================================
# PART 2: Casimir tensor
# =========================================================================

class TestCasimirTensor:
    """Verify Casimir tensor Omega = sum kappa^{ab} t_a tensor t_b."""

    def test_sl2_casimir_explicit(self):
        """Casimir for sl_2: Omega = (1/2)(h tensor h) + e tensor f + f tensor e.

        In our basis {e=0, h=1, f=2}:
          kappa = [[0, 0, 1],
                   [0, 2, 0],
                   [1, 0, 0]]
          kappa^{-1} = [[0, 0, 1],
                        [0, 1/2, 0],
                        [1, 0, 0]]

        So Omega^{00} = 0, Omega^{02} = 1 (e tensor f),
           Omega^{11} = 1/2 (h tensor h),
           Omega^{20} = 1 (f tensor e).
        """
        g = make_sl2()
        omega = casimir_tensor(g)

        # Check specific entries
        assert abs(omega[0, 2] - 1.0) < 1e-12, f"Omega(e,f) = {omega[0,2]}, expected 1"
        assert abs(omega[2, 0] - 1.0) < 1e-12, f"Omega(f,e) = {omega[2,0]}, expected 1"
        assert abs(omega[1, 1] - 0.5) < 1e-12, f"Omega(h,h) = {omega[1,1]}, expected 0.5"
        assert abs(omega[0, 0]) < 1e-12, f"Omega(e,e) = {omega[0,0]}, expected 0"
        assert abs(omega[0, 1]) < 1e-12, f"Omega(e,h) = {omega[0,1]}, expected 0"
        assert abs(omega[1, 0]) < 1e-12, f"Omega(h,e) = {omega[1,0]}, expected 0"
        assert abs(omega[1, 2]) < 1e-12, f"Omega(h,f) = {omega[1,2]}, expected 0"
        assert abs(omega[2, 1]) < 1e-12, f"Omega(f,h) = {omega[2,1]}, expected 0"
        assert abs(omega[2, 2]) < 1e-12, f"Omega(f,f) = {omega[2,2]}, expected 0"

    def test_sl2_casimir_is_inverse_killing(self):
        """Casimir tensor = inverse Killing form."""
        g = make_sl2()
        omega = casimir_tensor(g)
        product = g.kappa @ omega

        # kappa * kappa^{-1} = identity
        assert np.max(np.abs(product - np.eye(3))) < 1e-12

    def test_sl2_casimir_trace(self):
        """Tr(Omega) = sum kappa^{aa} = dim(g) for sl_N.

        For sl_2: Omega^{00} + Omega^{11} + Omega^{22} = 0 + 1/2 + 0 = 1/2.
        Wait: Tr(kappa^{-1}) for sl_2 with our normalization.

        kappa^{-1} diagonal: 0, 1/2, 0 -> Tr = 1/2.
        But there are off-diagonal terms: kappa^{02} = kappa^{20} = 1.
        Trace only counts diagonal: Tr = 1/2.

        Actually, Tr(Omega) = sum_a Omega^{aa} = sum_a kappa^{aa} = 1/2.
        The quadratic Casimir eigenvalue in the adjoint is
        C_2(adj) = 2h^v = 4 for sl_2 (with our normalization).
        """
        g = make_sl2()
        omega = casimir_tensor(g)
        tr = np.trace(omega)
        assert abs(tr - 0.5) < 1e-12, f"Tr(Omega) = {tr}, expected 0.5"

    def test_sl3_casimir_is_inverse_killing(self):
        """Casimir for sl_3 is inverse of Killing form."""
        g = make_sl3()
        omega = casimir_tensor(g)
        product = g.kappa @ omega
        assert np.max(np.abs(product - np.eye(8))) < 1e-12

    def test_sl3_casimir_symmetry(self):
        """Casimir tensor is symmetric: Omega^{ab} = Omega^{ba}."""
        g = make_sl3()
        omega = casimir_tensor(g)
        assert np.max(np.abs(omega - omega.T)) < 1e-12

    def test_u1_casimir(self):
        """Casimir for u(1): Omega = 1 (identity)."""
        g = make_u1()
        omega = casimir_tensor(g)
        assert abs(omega[0, 0] - 1.0) < 1e-12


# =========================================================================
# PART 3: Collision residue r-matrix
# =========================================================================

class TestCollisionResidue:
    """Core computation: r(z) = Res^coll_{0,2}(Theta_A) for affine g_k."""

    def test_sl2_k1_rmatrix(self):
        """r-matrix for sl_2 at k=1: r(z) = Omega/z.

        OPE: J^a(z) J^b(w) ~ kappa^{ab}/(z-w)^2 + f^{ab}_c J^c/(z-w)
        d-log extraction: z^{-2} -> z^{-1}, z^{-1} -> z^{0}.
        Singular part: r(z) = kappa/z.

        The tensor r-matrix is r_{12}(z) = Omega/z where Omega = kappa^{-1}.
        """
        g = make_sl2()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)

        assert result['pole_absorption_verified'], "Pole absorption failed"
        assert result['r_matrix_max_pole'] == 1, f"Expected max pole 1, got {result['r_matrix_max_pole']}"
        assert result['r_equals_k_omega_over_z'], "r(z) != k*Omega/z"

    def test_sl2_k1_pole_coefficient(self):
        """At k=1, the pole-1 coefficient of r(z) is exactly kappa_{ab}."""
        g = make_sl2()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)

        pole_coeff = result['r_pole_coefficients'][1]
        expected = g.kappa.copy()  # k=1, so k*kappa = kappa
        assert np.max(np.abs(pole_coeff - expected)) < 1e-12

    def test_sl2_k_general_rmatrix(self):
        """r-matrix for sl_2 at general k: r(z) = k*Omega/z."""
        g = make_sl2()
        for k in [1, 2, 5, 10, -1, 0.5]:
            ope = AffineOPE(g=g, k=float(k))
            result = collision_residue_rmatrix(ope)
            assert result['r_equals_k_omega_over_z'], f"Failed at k={k}"

            # Check pole coefficient = k * kappa
            pole_coeff = result['r_pole_coefficients'][1]
            expected = k * g.kappa
            assert np.max(np.abs(pole_coeff - expected)) < 1e-12, \
                f"Pole coeff mismatch at k={k}"

    def test_sl2_k0_trace_pole_vanishes_regular_lie_part_survives(self):
        """At k=0 only the singular trace-form pole vanishes.

        The affine simple-pole bracket f^{ab}_c J^c is still present in
        the OPE and becomes the regular z^0 part after d-log extraction.
        Thus k=0 is not a symmetric-category collapse for non-abelian g.
        """
        g = make_sl2()
        result = collision_residue_rmatrix(AffineOPE(g=g, k=0.0))

        pole_coeff = result['r_pole_coefficients'][1]
        assert np.max(np.abs(pole_coeff)) < 1e-12

        regular = result['r_regular_part']
        assert regular is not None
        assert np.max(np.abs(regular)) > 0
        assert abs(regular[0, 2, 1] - 1.0) < 1e-12
        assert abs(regular[1, 0, 0] - 2.0) < 1e-12

    def test_trace_form_kz_normalizations_are_distinct_at_k0(self):
        """Trace-form and KZ affine normalizations have different k=0 tests."""
        h_dual = 2.0

        assert affine_trace_form_scale(0.0) == 0.0
        assert affine_level_prefixed_kz_scale(0.0, h_dual) == 0.0
        assert affine_kz_scale(0.0, h_dual) == 1.0 / h_dual

        with pytest.raises(ZeroDivisionError):
            affine_kz_scale(-h_dual, h_dual)
        with pytest.raises(ZeroDivisionError):
            affine_level_prefixed_kz_scale(-h_dual, h_dual)

    def test_sl3_k1_rmatrix(self):
        """r-matrix for sl_3 at k=1: r(z) = Omega/z."""
        g = make_sl3()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)

        assert result['pole_absorption_verified'], "Pole absorption failed for sl_3"
        assert result['r_matrix_max_pole'] == 1
        assert result['r_equals_k_omega_over_z'], "r(z) != Omega/z for sl_3"

    def test_sl3_pole_coefficient(self):
        """Pole-1 coefficient for sl_3 at k=1 is kappa_{ab}."""
        g = make_sl3()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)

        pole_coeff = result['r_pole_coefficients'][1]
        assert np.max(np.abs(pole_coeff - g.kappa)) < 1e-12

    def test_u1_k1_rmatrix(self):
        """r-matrix coefficient for u(1) at k=1: 1/z (Heisenberg).

        This is the abelian case: Omega_H = 1 in the chosen basis, so
        the tensor kernel k*Omega_H/z is stored as coefficient k/z.
        """
        g = make_u1()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)

        assert result['pole_absorption_verified']
        assert result['r_matrix_max_pole'] == 1
        # k * kappa = 1 * 1 = 1
        pole_coeff = result['r_pole_coefficients'][1]
        assert abs(pole_coeff[0, 0] - 1.0) < 1e-12

    def test_no_higher_poles(self):
        """For KM, the r-matrix has NO poles beyond order 1.

        OPE max pole = 2 -> r-matrix max pole = 1. Nothing at z^{-2} or higher.
        """
        for g in [make_sl2(), make_sl3(), make_u1()]:
            ope = AffineOPE(g=g, k=1.0)
            result = collision_residue_rmatrix(ope)
            poles = result['r_pole_coefficients']
            for order in poles:
                assert order <= 1, f"Unexpected pole at order {order} for {g.name}"

    def test_regular_part_is_structure_constants(self):
        """The regular part of r(z) comes from the simple-pole OPE term.

        OPE at z^{-1}: f^{ab}_c J^c -> d-log extraction -> z^{0} (regular).
        So the regular part is the structure constant tensor.
        """
        g = make_sl2()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)

        reg = result['r_regular_part']
        # reg should be the structure constant tensor f^{ab}_c
        assert reg is not None, "Regular part should be non-None"
        # For sl_2: f is a 3x3x3 tensor
        assert reg.shape == (3, 3, 3), f"Expected shape (3,3,3), got {reg.shape}"
        # Check specific values
        assert abs(reg[0, 2, 1] - 1.0) < 1e-12   # f^{ef}_h = 1
        assert abs(reg[1, 0, 0] - 2.0) < 1e-12   # f^{he}_e = 2


# =========================================================================
# PART 4: Classical Yang-Baxter equation
# =========================================================================

class TestCYBE:
    """Verify the classical Yang-Baxter equation via IBR."""

    def test_sl2_cybe(self):
        """CYBE for sl_2 Casimir: [Omega_{12}+Omega_{13}, Omega_{23}] = 0."""
        g = make_sl2()
        result = verify_cybe(g)

        assert result['ibr_holds'], \
            f"IBR failed for sl_2, max violation = {result['ibr_max_violation']}"
        assert result['centrality_holds'], \
            f"Casimir centrality failed, max violation = {result['centrality_max_violation']}"
        assert result['cybe_satisfied'], "CYBE not satisfied for sl_2"

    def test_sl3_cybe(self):
        """CYBE for sl_3 Casimir."""
        g = make_sl3()
        result = verify_cybe(g)

        assert result['ibr_holds'], \
            f"IBR failed for sl_3, max violation = {result['ibr_max_violation']}"
        assert result['centrality_holds'], \
            f"Casimir centrality failed for sl_3, max violation = {result['centrality_max_violation']}"
        assert result['cybe_satisfied'], "CYBE not satisfied for sl_3"

    def test_ibr_violation_magnitude(self):
        """IBR violation should be machine-epsilon level."""
        for g in [make_sl2(), make_sl3()]:
            result = verify_cybe(g)
            assert result['ibr_max_violation'] < 1e-10, \
                f"IBR violation {result['ibr_max_violation']} too large for {g.name}"

    def test_centrality_violation_magnitude(self):
        """Casimir centrality violation should be machine-epsilon level."""
        for g in [make_sl2(), make_sl3()]:
            result = verify_cybe(g)
            assert result['centrality_max_violation'] < 1e-10, \
                f"Centrality violation {result['centrality_max_violation']} too large for {g.name}"


# =========================================================================
# PART 5: Pole absorption (AP19)
# =========================================================================

class TestPoleAbsorption:
    """Verify AP19: the bar kernel absorbs a pole."""

    def test_km_pole_absorption(self):
        """KM: OPE poles at z^{-2}, z^{-1} -> r-matrix pole at z^{-1} only.

        OPE has max pole 2. r-matrix has max pole 1. Shift = 1.
        """
        ope_poles = {2: 'k*kappa (central term)', 1: 'f^{ab}_c J^c (current)'}
        r_poles = {1: 'k*kappa/z (Casimir)'}
        result = verify_pole_absorption(ope_poles, r_poles)

        assert result['ap19_verified'], "AP19 not verified for KM"
        assert result['ope_max_pole'] == 2
        assert result['r_max_pole'] == 1
        assert result['shift_by_one']

    def test_virasoro_pole_absorption(self):
        """Virasoro: OPE z^{-4}, z^{-2}, z^{-1} -> r-matrix z^{-3}, z^{-1}.

        OPE max pole = 4. r-matrix max pole = 3. Shift = 1.
        The z^{-1} OPE term becomes regular (z^{0}).
        """
        ope_poles = {4: 'c/2', 2: '2T', 1: 'dT'}
        r_poles = {3: '(c/2)/z^3', 1: '2T/z'}
        result = verify_pole_absorption(ope_poles, r_poles)

        assert result['ap19_verified'], "AP19 not verified for Virasoro"
        assert result['ope_max_pole'] == 4
        assert result['r_max_pole'] == 3
        assert result['shift_by_one']

    def test_virasoro_rmatrix_structure(self):
        """Virasoro r-matrix: r(z) = (c/2)/z^3 + 2T/z.

        NOT (c/2)/z^4 + 2T/z^2 + dT/z (that's the OPE, not the r-matrix).
        """
        for c in [1, 26, 13, 0.5]:
            result = virasoro_ope_to_rmatrix(c)

            # r-matrix has poles at z^{-3} and z^{-1}
            assert 3 in result['r_matrix_poles'], "Missing z^{-3} pole"
            assert 1 in result['r_matrix_poles'], "Missing z^{-1} pole"
            assert result['r_max_pole'] == 3, f"Expected max pole 3, got {result['r_max_pole']}"

            # No z^{-4} pole in r-matrix (that's the OPE, not r-matrix)
            assert 4 not in result['r_matrix_poles'], \
                "z^{-4} should NOT appear in r-matrix (AP19)"

            # No z^{-2} pole in r-matrix
            assert 2 not in result['r_matrix_poles'], \
                "z^{-2} should NOT appear in r-matrix"

    def test_heisenberg_pole_absorption(self):
        """Heisenberg: OPE z^{-2} only -> r-matrix z^{-1} only.

        J(z) J(w) ~ k/(z-w)^2. After d-log: tensor kernel
        k*Omega_H/z, stored here as coefficient k/z.
        """
        ope_poles = {2: 'k (central term)'}
        r_poles = {1: 'k/z'}
        result = verify_pole_absorption(ope_poles, r_poles)

        assert result['ap19_verified']
        assert result['ope_max_pole'] == 2
        assert result['r_max_pole'] == 1

    def test_numerical_pole_absorption_sl2(self):
        """Numerical verification of pole absorption for sl_2."""
        g = make_sl2()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)

        # OPE max pole = 2
        assert ope.max_pole_order == 2
        # r-matrix max pole = 1
        assert result['r_matrix_max_pole'] == 1
        # Shift = 1
        assert result['pole_shift'] == 1


# =========================================================================
# PART 6: Full pipeline
# =========================================================================

class TestFullPipeline:
    """End-to-end collision residue computation."""

    def test_sl2_k1_full(self):
        """Full pipeline for sl_2 at k=1."""
        g = make_sl2()
        result = full_collision_residue_computation(g, k=1.0)

        assert result['lie_algebra_valid'], "Lie algebra axioms failed"
        assert result['r_equals_k_omega_over_z'], "r(z) != Omega/z"
        assert result['cybe_satisfied'], "CYBE failed"
        assert result['all_checks_passed'], "Not all checks passed"

        # kappa(sl_2, k=1) = 3*(1+2)/(2*2) = 9/4
        assert abs(result['kappa_A'] - 9.0/4) < 1e-12, \
            f"kappa = {result['kappa_A']}, expected 9/4"

    def test_sl2_k2_full(self):
        """Full pipeline for sl_2 at k=2."""
        g = make_sl2()
        result = full_collision_residue_computation(g, k=2.0)

        assert result['all_checks_passed']

        # kappa(sl_2, k=2) = 3*(2+2)/(2*2) = 3
        assert abs(result['kappa_A'] - 3.0) < 1e-12

    def test_sl3_k1_full(self):
        """Full pipeline for sl_3 at k=1."""
        g = make_sl3()
        result = full_collision_residue_computation(g, k=1.0)

        assert result['lie_algebra_valid'], "sl_3 Lie algebra axioms failed"
        assert result['r_equals_k_omega_over_z'], "r(z) != Omega/z for sl_3"
        assert result['cybe_satisfied'], "CYBE failed for sl_3"
        assert result['all_checks_passed'], "Not all checks passed for sl_3"

        # kappa(sl_3, k=1) = 8*(1+3)/(2*3) = 16/3
        assert abs(result['kappa_A'] - 16.0/3) < 1e-12, \
            f"kappa = {result['kappa_A']}, expected 16/3"

    def test_u1_k1_full(self):
        """Full pipeline for u(1) at k=1 (Heisenberg)."""
        g = make_u1()
        result = full_collision_residue_computation(g, k=1.0)

        assert result['lie_algebra_valid']
        assert result['r_equals_k_omega_over_z']
        assert result['cybe_satisfied']
        assert result['all_checks_passed']

        # kappa(Heisenberg, k=1) = 1 (abelian convention)
        assert abs(result['kappa_A'] - 1.0) < 1e-12

    def test_sl2_critical_level(self):
        """sl_2 at critical level k = -h^v = -2.

        kappa(sl_2, k=-2) = 3*(-2+2)/(2*2) = 0.
        The r-matrix still exists: r(z) = -2*Omega/z.
        The Sugawara construction is UNDEFINED (AP: "Sugawara UNDEFINED
        at critical level"), but the bar complex r-matrix is well-defined.
        """
        g = make_sl2()
        result = full_collision_residue_computation(g, k=-2.0)

        assert result['lie_algebra_valid']
        assert result['r_equals_k_omega_over_z'], \
            "r-matrix should still be well-defined at critical level"
        assert result['cybe_satisfied']

        # kappa = 0 at critical level
        assert abs(result['kappa_A']) < 1e-12, \
            f"kappa should be 0 at critical level, got {result['kappa_A']}"


# =========================================================================
# PART 7: Cross-checks and consistency
# =========================================================================

class TestCrossChecks:
    """Cross-family consistency, known values, and edge cases."""

    def test_kappa_formula_sl2(self):
        """kappa(sl_2, k) = 3(k+2)/4 = dim(g)(k+h^v)/(2h^v).

        AP1 warning: kappa formulas are distinct for each family.
        """
        g = make_sl2()
        for k in [1, 2, 3, -1, 0, 10]:
            result = full_collision_residue_computation(g, k=float(k))
            expected = 3.0 * (k + 2) / 4.0
            assert abs(result['kappa_A'] - expected) < 1e-12, \
                f"kappa(sl_2, k={k}) = {result['kappa_A']}, expected {expected}"

    def test_kappa_formula_sl3(self):
        """kappa(sl_3, k) = 8(k+3)/6 = 4(k+3)/3."""
        g = make_sl3()
        for k in [1, 2, 3, -1, 0]:
            result = full_collision_residue_computation(g, k=float(k))
            expected = 8.0 * (k + 3) / 6.0
            assert abs(result['kappa_A'] - expected) < 1e-12, \
                f"kappa(sl_3, k={k}) = {result['kappa_A']}, expected {expected}"

    def test_feigin_frenkel_duality(self):
        """Feigin-Frenkel: k <-> -k - 2h^v.

        kappa(g, k) + kappa(g, -k-2h^v) = 0 for KM (AP24: complementarity).
        """
        for g, name in [(make_sl2(), 'sl2'), (make_sl3(), 'sl3')]:
            for k in [1, 2, 3, 5]:
                k_dual = -k - 2 * g.h_dual
                kappa_k = g.dim * (k + g.h_dual) / (2 * g.h_dual)
                kappa_dual = g.dim * (k_dual + g.h_dual) / (2 * g.h_dual)
                kappa_sum = kappa_k + kappa_dual
                assert abs(kappa_sum) < 1e-12, \
                    f"kappa + kappa' = {kappa_sum} != 0 for {name} at k={k}"

    def test_casimir_eigenvalue_adjoint(self):
        """Quadratic Casimir eigenvalue in adjoint representation.

        For sl_N with our normalization (Killing form / (2h^v)):
          C_2(adj) = sum_{a,b} kappa^{ab} * ad(t_a) * ad(t_b)

        The adjoint Casimir eigenvalue is 2h^v with our normalization.
        We compute: Tr(C_2(adj)) = sum_{a,b} kappa^{ab} sum_c f^{ac}_d f^{bd}_c
        and check it equals 2h^v * dim.
        """
        for g in [make_sl2(), make_sl3()]:
            omega = casimir_tensor(g)
            d = g.dim

            # C_2(adj)^c_e = sum_{a,b} kappa^{ab} f^{ac}_d f^{bd}_e  ... wait
            # Actually: C_2(adj) = sum_{a,b} kappa^{ab} ad(t_a) ad(t_b)
            # ad(t_a)^c_d = f^{ac}_d
            # C_2(adj)^c_e = sum_{a,b,d} kappa^{ab} f^{ac}_d f^{bd}_e

            casimir_adj = np.zeros((d, d), dtype=float)
            for c in range(d):
                for e in range(d):
                    val = 0.0
                    for a in range(d):
                        for b in range(d):
                            if abs(omega[a, b]) < 1e-15:
                                continue
                            for dd in range(d):
                                val += omega[a, b] * g.f[a, c, dd] * g.f[b, dd, e]
                    casimir_adj[c, e] = val

            # Should be proportional to identity: C_2(adj) = 2h^v * Id
            # with our normalization (Killing / (2h^v))
            expected_eigenvalue = 2 * g.h_dual
            ratio = casimir_adj / expected_eigenvalue

            # Check it's proportional to identity
            assert np.max(np.abs(ratio - np.eye(d))) < 1e-10, \
                f"Casimir not proportional to identity for {g.name}: " \
                f"max deviation = {np.max(np.abs(ratio - np.eye(d)))}"

    def test_ope_max_pole_order(self):
        """KM OPE max pole order is 2 (double pole from central extension)."""
        for g in [make_sl2(), make_sl3(), make_u1()]:
            ope = AffineOPE(g=g, k=1.0)
            assert ope.max_pole_order == 2

    def test_r_matrix_only_simple_pole(self):
        """After d-log extraction, KM r-matrix has ONLY a simple pole.

        The double pole from the central extension becomes a simple pole.
        The simple pole from the structure constants becomes regular.
        No higher poles exist.
        """
        for g in [make_sl2(), make_sl3()]:
            ope = AffineOPE(g=g, k=1.0)
            result = collision_residue_rmatrix(ope)
            poles = result['r_pole_coefficients']

            # Only pole at order 1
            assert set(poles.keys()) == {1}, \
                f"Expected only pole at order 1, got {set(poles.keys())} for {g.name}"

    def test_r_matrix_level_scaling(self):
        """r(z) scales linearly with level k: r_k(z) = k * r_1(z).

        This follows from the OPE: the double pole is k*kappa, so
        after d-log extraction, r(z) = k * kappa / z = k * Omega / z.
        """
        g = make_sl2()
        ope1 = AffineOPE(g=g, k=1.0)
        res1 = collision_residue_rmatrix(ope1)
        coeff1 = res1['r_pole_coefficients'][1]

        for k in [2, 3, 5, -1, 0.5]:
            ope_k = AffineOPE(g=g, k=float(k))
            res_k = collision_residue_rmatrix(ope_k)
            coeff_k = res_k['r_pole_coefficients'][1]

            expected = k * coeff1
            assert np.max(np.abs(coeff_k - expected)) < 1e-12, \
                f"Level scaling failed at k={k}"


# =========================================================================
# PART 8: AP19 explicit numerical examples
# =========================================================================

class TestAP19Explicit:
    """Explicit numerical verification of AP19 pole absorption.

    AP19: "the r-matrix lives one pole order below the OPE."
    """

    def test_km_ope_has_two_poles(self):
        """KM OPE: poles at z^{-2} (central) and z^{-1} (current)."""
        g = make_sl2()
        ope = AffineOPE(g=g, k=1.0)

        # Double pole coefficient (n=2): k * kappa
        c2 = ope.c_n_scalar(2)
        assert c2 is not None
        assert np.max(np.abs(c2 - g.kappa)) < 1e-12

        # Simple pole coefficient (n=1): structure constants
        c1 = ope.c_n_current(1)
        assert c1 is not None
        assert c1.shape == (3, 3, 3)

        # No triple pole
        assert ope.c_n_scalar(3) is None
        assert ope.c_n_current(3) is None

    def test_rmatrix_has_one_pole(self):
        """After d-log: r-matrix has single pole at z^{-1}."""
        g = make_sl2()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)

        # Only one pole
        assert len(result['r_pole_coefficients']) == 1
        assert 1 in result['r_pole_coefficients']

        # The coefficient is k * kappa = kappa (at k=1)
        pole_coeff = result['r_pole_coefficients'][1]
        assert np.max(np.abs(pole_coeff - g.kappa)) < 1e-12

    def test_bosonic_even_pole_rule(self):
        """AP19 remark: for a bosonic algebra, d-log extraction sends
        z^{-2n} to z^{-(2n-1)} (odd order). The r-matrix has
        NO even-order poles for a bosonic algebra with even-order OPE poles.

        KM: OPE has z^{-2} (even). r-matrix has z^{-1} (odd). Correct.
        Virasoro: OPE has z^{-4}, z^{-2} (even). r-matrix has z^{-3}, z^{-1} (odd).
        """
        # KM check
        g = make_sl2()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)
        for order in result['r_pole_coefficients']:
            assert order % 2 == 1, f"Even-order pole {order} in KM r-matrix"

        # Virasoro check
        vir = virasoro_ope_to_rmatrix(26)
        for order in vir['r_matrix_poles']:
            assert order % 2 == 1, f"Even-order pole {order} in Virasoro r-matrix"
