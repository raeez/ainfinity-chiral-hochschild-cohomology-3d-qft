r"""Tests for non-simply-laced collision residue r-matrices: B_2 and C_2.

Verifies:
  1. Lie algebra axioms for so(5) ≅ sp(4) (Jacobi, antisymmetry, Killing invariance)
  2. Casimir tensor computation with root-length dependent coefficients
  3. Collision residue r(z) = k * Omega / z
  4. CYBE via IBR (infinitesimal braid relations)
  5. Langlands duality B_2 <-> C_2
  6. FM-integral coefficients at degree 3
  7. Root decomposition of Casimir (short vs long root contributions)

References:
  Vol II, ordered_associative_chiral_kd_core.tex (non-simply-laced bar complex)
  Vol II, ht_bulk_boundary_line_core.tex (r-matrix identification)
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.non_simply_laced_rmatrix import (
    make_B2, make_C2,
    rep_rmatrix_B2_defining,
    langlands_duality_comparison,
    fm_integral_degree3_nonsimplylaced,
    casimir_root_decomposition,
    run_full_computation,
)
from lib.collision_residue_rmatrix import (
    verify_jacobi, verify_killing_invariance, verify_antisymmetry,
    casimir_tensor, AffineOPE, collision_residue_rmatrix,
    verify_cybe, full_collision_residue_computation,
)


# =========================================================================
# PART 1: Lie algebra axioms for B_2 / C_2
# =========================================================================

class TestLieAlgebraAxiomsB2C2:
    """Verify Lie algebra data is correct for so(5) ≅ sp(4)."""

    def test_B2_jacobi(self):
        """Jacobi identity for B_2 = so(5)."""
        g = make_B2()
        assert verify_jacobi(g), "B_2 fails Jacobi identity"

    def test_B2_antisymmetry(self):
        """Antisymmetry f^{ab}_c = -f^{ba}_c for B_2."""
        g = make_B2()
        assert verify_antisymmetry(g), "B_2 structure constants not antisymmetric"

    def test_B2_killing_invariance(self):
        """Ad-invariance of Killing form for B_2."""
        g = make_B2()
        assert verify_killing_invariance(g), "B_2 Killing form not ad-invariant"

    def test_B2_dimension(self):
        """B_2 = so(5) has dimension 10, rank 2, h^v = 3."""
        g = make_B2()
        assert g.dim == 10
        assert g.rank == 2
        assert g.h_dual == 3

    def test_C2_jacobi(self):
        """Jacobi identity for C_2 = sp(4)."""
        g = make_C2()
        assert verify_jacobi(g), "C_2 fails Jacobi identity"

    def test_C2_antisymmetry(self):
        """Antisymmetry for C_2."""
        g = make_C2()
        assert verify_antisymmetry(g), "C_2 structure constants not antisymmetric"

    def test_C2_killing_invariance(self):
        """Ad-invariance of Killing form for C_2."""
        g = make_C2()
        assert verify_killing_invariance(g), "C_2 Killing form not ad-invariant"

    def test_C2_dimension(self):
        """C_2 = sp(4) has dimension 10, rank 2, h^v = 3."""
        g = make_C2()
        assert g.dim == 10
        assert g.rank == 2
        assert g.h_dual == 3

    def test_isomorphism(self):
        """so(5) ≅ sp(4): same structure constants."""
        g_B2 = make_B2()
        g_C2 = make_C2()
        assert np.allclose(g_B2.f, g_C2.f, atol=1e-12), \
            "B_2 and C_2 should have identical structure constants"
        assert np.allclose(g_B2.kappa, g_C2.kappa, atol=1e-12), \
            "B_2 and C_2 should have identical Killing form"


# =========================================================================
# PART 2: Casimir tensor for B_2 / C_2
# =========================================================================

class TestCasimirB2C2:
    """Verify Casimir tensor for non-simply-laced types."""

    def test_casimir_is_inverse_killing(self):
        """Casimir = inverse of Killing form."""
        g = make_B2()
        omega = casimir_tensor(g)
        product = g.kappa @ omega
        assert np.max(np.abs(product - np.eye(10))) < 1e-10, \
            "Omega is not inverse of kappa"

    def test_casimir_symmetry(self):
        """Casimir tensor is symmetric."""
        g = make_B2()
        omega = casimir_tensor(g)
        assert np.max(np.abs(omega - omega.T)) < 1e-12, \
            "Casimir not symmetric"

    def test_casimir_B2_equals_C2(self):
        """B_2 and C_2 Casimirs are identical (same abstract algebra)."""
        omega_B2 = casimir_tensor(make_B2())
        omega_C2 = casimir_tensor(make_C2())
        assert np.allclose(omega_B2, omega_C2, atol=1e-12)

    def test_root_length_dependence(self):
        """Casimir coefficients differ for short vs long roots.

        For so(5) ≅ sp(4) with basis from C_2 labelling:
          - alpha_1 = e_1-e_2 (short for C_2): kappa(e, f) = 2, Casimir coeff = 1/2
          - alpha_2 = 2e_2 (long for C_2): kappa(e, f) = 1, Casimir coeff = 1
          - alpha_1+alpha_2 (short): kappa(e, f) = 2, Casimir coeff = 1/2
          - 2alpha_1+alpha_2 (long): kappa(e, f) = 1, Casimir coeff = 1

        The Casimir coefficient for a root alpha is 1/kappa(e_alpha, f_alpha).
        For non-simply-laced types, this depends on the root length.
        """
        g = make_B2()
        decomp = casimir_root_decomposition(g)

        # Short roots (C_2 labelling): e_{a1}, e_{a1+a2}
        # kappa(e,f) = 2, Casimir coeff = 1/2
        assert abs(decomp['root_casimir_coefficients']['e_{a1}']['e_x_f'] - 0.5) < 1e-12
        assert abs(decomp['root_casimir_coefficients']['e_{a1+a2}']['e_x_f'] - 0.5) < 1e-12

        # Long roots (C_2 labelling): e_{a2}, e_{2a1+a2}
        # kappa(e,f) = 1, Casimir coeff = 1
        assert abs(decomp['root_casimir_coefficients']['e_{a2}']['e_x_f'] - 1.0) < 1e-12
        assert abs(decomp['root_casimir_coefficients']['e_{2a1+a2}']['e_x_f'] - 1.0) < 1e-12

    def test_casimir_coefficient_is_reciprocal_kappa(self):
        """For each root alpha: Omega(e_alpha, f_alpha) = 1/kappa(e_alpha, f_alpha)."""
        g = make_B2()
        decomp = casimir_root_decomposition(g)

        for label, data in decomp['root_casimir_coefficients'].items():
            assert data['predicted'] is not None
            assert abs(data['e_x_f'] - data['predicted']) < 1e-12, \
                f"Casimir coeff {data['e_x_f']} != 1/kappa = {data['predicted']} for {label}"

    def test_cartan_part_of_casimir(self):
        """Cartan part of Casimir is B^{-1} where B is symmetrised Cartan restricted to h."""
        g = make_B2()
        decomp = casimir_root_decomposition(g)

        # The Cartan Killing form in our basis is:
        kappa_h = decomp['kappa_cartan']
        B_inv = decomp['B_inv_cartan']

        # Verify B_inv is indeed the inverse of kappa_h
        assert np.allclose(kappa_h @ B_inv, np.eye(2), atol=1e-12)

    def test_casimir_adjoint_eigenvalue(self):
        """Quadratic Casimir eigenvalue in the adjoint = 2*h^v = 6."""
        g = make_B2()
        omega = casimir_tensor(g)
        d = g.dim

        # C_2(adj)^c_e = sum_{a,b} kappa^{ab} f^{ac}_d f^{bd}_e
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

        # Should be 2*h^v * Id = 6 * Id
        expected = 2 * g.h_dual  # = 6
        ratio = casimir_adj / expected
        assert np.max(np.abs(ratio - np.eye(d))) < 1e-10, \
            f"Casimir eigenvalue wrong: max deviation {np.max(np.abs(ratio - np.eye(d)))}"


# =========================================================================
# PART 3: Collision residue r-matrix
# =========================================================================

class TestCollisionResidueB2C2:
    """Collision residue r(z) = k * Omega / z for B_2 and C_2."""

    def test_B2_k1_rmatrix(self):
        """r-matrix for B_2 at k=1."""
        g = make_B2()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)
        assert result['pole_absorption_verified']
        assert result['r_matrix_max_pole'] == 1
        assert result['r_equals_k_omega_over_z']

    def test_C2_k1_rmatrix(self):
        """r-matrix for C_2 at k=1."""
        g = make_C2()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)
        assert result['pole_absorption_verified']
        assert result['r_matrix_max_pole'] == 1
        assert result['r_equals_k_omega_over_z']

    def test_B2_k_general(self):
        """r-matrix for B_2 at general levels."""
        g = make_B2()
        for k in [1, 2, 5, -1, -3, 0.5]:
            ope = AffineOPE(g=g, k=float(k))
            result = collision_residue_rmatrix(ope)
            assert result['r_equals_k_omega_over_z'], f"Failed at k={k}"

    def test_pole_coefficient_is_k_kappa(self):
        """Pole-1 coefficient of r(z) = k * kappa."""
        g = make_B2()
        for k in [1, 2, -1]:
            ope = AffineOPE(g=g, k=float(k))
            result = collision_residue_rmatrix(ope)
            pole = result['r_pole_coefficients'][1]
            expected = k * g.kappa
            assert np.max(np.abs(pole - expected)) < 1e-12, f"Mismatch at k={k}"

    def test_no_higher_poles(self):
        """r-matrix has only simple pole (from double-pole OPE)."""
        for g in [make_B2(), make_C2()]:
            ope = AffineOPE(g=g, k=1.0)
            result = collision_residue_rmatrix(ope)
            for order in result['r_pole_coefficients']:
                assert order <= 1, f"Unexpected pole at order {order}"

    def test_level_scaling(self):
        """r(z) scales linearly with level: r_k = k * r_1."""
        g = make_B2()
        res1 = collision_residue_rmatrix(AffineOPE(g=g, k=1.0))
        coeff1 = res1['r_pole_coefficients'][1]

        for k in [2, 3, -1, 0.5]:
            res_k = collision_residue_rmatrix(AffineOPE(g=g, k=float(k)))
            coeff_k = res_k['r_pole_coefficients'][1]
            assert np.max(np.abs(coeff_k - k * coeff1)) < 1e-12


# =========================================================================
# PART 4: CYBE verification
# =========================================================================

class TestCYBE_B2C2:
    """Classical Yang-Baxter equation for B_2 and C_2."""

    def test_B2_ibr(self):
        """IBR [Omega_{12}+Omega_{13}, Omega_{23}] = 0 for B_2."""
        g = make_B2()
        result = verify_cybe(g)
        assert result['ibr_holds'], \
            f"IBR failed, max violation = {result['ibr_max_violation']}"

    def test_B2_centrality(self):
        """Casimir centrality [Omega, Delta(x)] = 0 for B_2."""
        g = make_B2()
        result = verify_cybe(g)
        assert result['centrality_holds'], \
            f"Centrality failed, max violation = {result['centrality_max_violation']}"

    def test_B2_cybe(self):
        """Full CYBE for B_2."""
        g = make_B2()
        result = verify_cybe(g)
        assert result['cybe_satisfied']

    def test_C2_cybe(self):
        """Full CYBE for C_2."""
        g = make_C2()
        result = verify_cybe(g)
        assert result['cybe_satisfied']

    def test_violation_is_machine_epsilon(self):
        """IBR and centrality violations should be at machine precision."""
        for g in [make_B2(), make_C2()]:
            result = verify_cybe(g)
            assert result['ibr_max_violation'] < 1e-10
            assert result['centrality_max_violation'] < 1e-10


# =========================================================================
# PART 5: Langlands duality
# =========================================================================

class TestLanglandsDuality:
    """Langlands duality B_2^L = C_2."""

    def test_cartan_transpose(self):
        """A_{C_2} = A_{B_2}^T."""
        result = langlands_duality_comparison(1.0)
        assert result['cartan_transpose']

    def test_abstract_casimir_match(self):
        """Abstract Casimir is the same (same Lie algebra)."""
        result = langlands_duality_comparison(1.0)
        assert result['casimir_match']

    def test_rep_rmatrix_match(self):
        """4-dim rep R-matrices match (same algebra, same rep)."""
        result = langlands_duality_comparison(1.0)
        assert result['rep_rmatrix_match']

    def test_symmetriser_swap(self):
        """D_{B_2} and D_{C_2} are related by swapping entries."""
        result = langlands_duality_comparison(1.0)
        assert result['symmetriser_swap']

    def test_symmetrised_cartan_matrices(self):
        """B_{B_2} = [[4,-2],[-2,2]] and B_{C_2} = [[2,-2],[-2,4]].

        These are related by swapping diagonal entries: Langlands duality
        exchanges long and short root lengths.
        """
        result = langlands_duality_comparison(1.0)
        B_B2 = result['B_B2']
        B_C2 = result['B_C2']

        assert np.allclose(B_B2, [[4, -2], [-2, 2]], atol=1e-12)
        assert np.allclose(B_C2, [[2, -2], [-2, 4]], atol=1e-12)

        # B_C2 = flip of B_B2 (swap rows and columns)
        P = np.array([[0, 1], [1, 0]])
        assert np.allclose(B_C2, P @ B_B2 @ P, atol=1e-12)

    def test_dual_level(self):
        """Langlands dual level: k^L = r^2 * k where r^2 = long^2/short^2 = 2."""
        result = langlands_duality_comparison(3.0)
        assert result['k_dual'] == 6.0  # 2 * 3

    def test_dual_coxeter_number_match(self):
        """Both B_2 and C_2 have h^v = 3."""
        g_B2 = make_B2()
        g_C2 = make_C2()
        assert g_B2.h_dual == 3
        assert g_C2.h_dual == 3


# =========================================================================
# PART 6: FM-integral coefficients
# =========================================================================

class TestFMIntegral:
    """FM-integral coefficients at degree 3 for non-simply-laced types."""

    def test_simply_laced_reference(self):
        """Simply-laced: B(1,1) = 1."""
        result = fm_integral_degree3_nonsimplylaced([1.0, 2.0])
        assert abs(result['simply_laced_reference'] - 1.0) < 1e-12

    def test_short_short(self):
        """Two short roots: B(1,1) = 1 (same as simply-laced)."""
        result = fm_integral_degree3_nonsimplylaced([1.0, 2.0])
        data = result['beta_integrals'][(1.0, 1.0)]
        assert abs(data['beta'] - 1.0) < 1e-12

    def test_short_long(self):
        """Short + long roots: B(1,2) = 1/2."""
        result = fm_integral_degree3_nonsimplylaced([1.0, 2.0])
        data = result['beta_integrals'][(1.0, 2.0)]
        assert abs(data['beta'] - 0.5) < 1e-12

    def test_long_long(self):
        """Two long roots: B(2,2) = 1/6."""
        result = fm_integral_degree3_nonsimplylaced([1.0, 2.0])
        data = result['beta_integrals'][(2.0, 2.0)]
        assert abs(data['beta'] - 1.0/6) < 1e-12

    def test_symmetry(self):
        """B(a,b) = B(b,a)."""
        result = fm_integral_degree3_nonsimplylaced([1.0, 2.0])
        sl = result['beta_integrals'][(1.0, 2.0)]['beta']
        ls = result['beta_integrals'][(2.0, 1.0)]['beta']
        assert abs(sl - ls) < 1e-12


# =========================================================================
# PART 7: Full pipeline
# =========================================================================

class TestFullPipeline:
    """End-to-end computation for B_2 and C_2."""

    def test_B2_full_k1(self):
        """Full pipeline for B_2 at k=1."""
        g = make_B2()
        result = full_collision_residue_computation(g, k=1.0)
        assert result['all_checks_passed']

    def test_C2_full_k1(self):
        """Full pipeline for C_2 at k=1."""
        g = make_C2()
        result = full_collision_residue_computation(g, k=1.0)
        assert result['all_checks_passed']

    def test_kappa_formula(self):
        """kappa(g_k) = dim(g) * (k + h^v) / (2 * h^v) for both B_2 and C_2.

        dim = 10, h^v = 3. kappa(k=1) = 10 * 4 / 6 = 20/3.
        """
        for g in [make_B2(), make_C2()]:
            result = full_collision_residue_computation(g, k=1.0)
            expected = 10.0 * (1 + 3) / (2 * 3)  # = 20/3
            assert abs(result['kappa_A'] - expected) < 1e-12, \
                f"kappa = {result['kappa_A']}, expected {expected}"

    def test_critical_level(self):
        """Critical level k = -h^v = -3 gives kappa = 0."""
        for g in [make_B2(), make_C2()]:
            result = full_collision_residue_computation(g, k=-3.0)
            assert abs(result['kappa_A']) < 1e-12

    def test_feigin_frenkel_duality(self):
        """FF duality: kappa(k) + kappa(-k-2h^v) = 0."""
        g = make_B2()
        for k in [1, 2, 5]:
            k_dual = -k - 2 * g.h_dual  # = -k - 6
            kappa_k = g.dim * (k + g.h_dual) / (2 * g.h_dual)
            kappa_dual = g.dim * (k_dual + g.h_dual) / (2 * g.h_dual)
            assert abs(kappa_k + kappa_dual) < 1e-12, \
                f"FF failed at k={k}: sum = {kappa_k + kappa_dual}"

    def test_full_run(self):
        """Run the complete computation and verify all checks pass."""
        results = run_full_computation(k=1.0, verbose=False)
        assert results['all_pass']


# =========================================================================
# PART 8: Representation-theoretic R-matrix
# =========================================================================

class TestRepRMatrix:
    """R-matrix in representations."""

    def test_4dim_rmatrix_is_symmetric(self):
        """The Casimir element in rep x rep is symmetric under swap."""
        g = make_B2()
        R = rep_rmatrix_B2_defining(g, k=1.0)

        # R lives in V x V where V = C^4.
        # Omega_{12} should satisfy: P * Omega_{12} * P = Omega_{21} = Omega_{12}
        # where P is the swap operator.
        n = 4
        P = np.zeros((n**2, n**2))
        for i in range(n):
            for j in range(n):
                P[i * n + j, j * n + i] = 1.0

        # Omega_{21} = P Omega_{12} P
        R_swapped = P @ R @ P
        # Should equal R (Casimir is symmetric element of g x g)
        assert np.allclose(R, R_swapped, atol=1e-12), \
            "R-matrix should be symmetric under particle swap"

    def test_4dim_rmatrix_traceless(self):
        """Trace of Casimir in rep x rep vanishes for sp(4).

        Tr(Omega) = Tr(sum kappa^{ab} rho(t_a) rho(t_b)) involves the
        quadratic Casimir, which has trace = C_2(V) * dim(V).
        But the TENSOR Omega has Tr = sum_a kappa^{aa} = Tr(kappa^{-1}).
        For sp(4) this is generally nonzero (it's the sum of the diagonal
        entries of the inverse Killing form).

        Actually, the trace of the 16x16 Omega should be:
        Tr(Omega) = sum_{a,b} kappa^{ab} Tr(rho(t_a)) Tr(rho(t_b)) = 0
        since sp(4) elements are traceless in the fundamental rep.
        """
        g = make_B2()
        R = rep_rmatrix_B2_defining(g, k=1.0)
        assert abs(np.trace(R)) < 1e-12, \
            f"Trace of R = {np.trace(R)}, expected 0"

    def test_4dim_rmatrix_scaling(self):
        """R scales linearly with k."""
        g = make_B2()
        R1 = rep_rmatrix_B2_defining(g, k=1.0)
        for k in [2, 3, -1, 0.5]:
            Rk = rep_rmatrix_B2_defining(g, k=float(k))
            assert np.allclose(Rk, k * R1, atol=1e-12)
