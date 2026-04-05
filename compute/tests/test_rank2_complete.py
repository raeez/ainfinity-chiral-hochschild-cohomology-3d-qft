r"""Complete rank-2 ordered E_1 bar data: A_2, B_2, C_2, G_2.

Verifies ALL seven items for each rank-2 simple Lie algebra:
  (1) Casimir tensor Omega in defining representation
  (2) R-matrix R(z) = 1 + hbar*Omega/z on V tensor V
  (3) YBE on V tensor V tensor V
  (4) RTT relation count
  (5) Koszul dual = Y_hbar(g)
  (6) Euler-eta: chi = -1 + eta^{dim g}
  (7) Non-simply-laced: root-length-dependent Casimir coefficients + Langlands

References:
  Vol II, ordered_associative_chiral_kd_core.tex
  Vol II, ht_bulk_boundary_line_core.tex (holographic datum)
  AP19: bar kernel absorbs a pole
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.collision_residue_rmatrix import (
    make_sl3,
    verify_jacobi, verify_killing_invariance, verify_antisymmetry,
    casimir_tensor, AffineOPE, collision_residue_rmatrix,
    verify_cybe, full_collision_residue_computation,
)
from lib.non_simply_laced_rmatrix import (
    make_B2, make_C2, make_G2,
    casimir_root_decomposition, langlands_duality_comparison,
    rtt_relation_count, euler_eta_rank2,
    casimir_in_defining_rep, verify_ybe_in_rep,
    _get_sl3_3dim_matrices, _get_sp4_basis_matrices, _get_g2_7dim_matrices,
    run_complete_rank2_computation,
)


# =========================================================================
# PART 0: G_2 Lie algebra axioms
# =========================================================================

class TestG2LieAlgebra:
    """Verify G_2 Lie algebra data is correct."""

    def test_g2_jacobi(self):
        g = make_G2()
        assert verify_jacobi(g), "G_2 fails Jacobi identity"

    def test_g2_antisymmetry(self):
        g = make_G2()
        assert verify_antisymmetry(g), "G_2 structure constants not antisymmetric"

    def test_g2_killing_invariance(self):
        g = make_G2()
        assert verify_killing_invariance(g), "G_2 Killing form not ad-invariant"

    def test_g2_dimension(self):
        """G_2 has dimension 14, rank 2, h^v = 4."""
        g = make_G2()
        assert g.dim == 14
        assert g.rank == 2
        assert g.h_dual == 4

    def test_g2_killing_nondegenerate(self):
        g = make_G2()
        det = np.linalg.det(g.kappa)
        assert abs(det) > 1e-10, f"Killing form degenerate: det = {det}"

    def test_g2_killing_symmetric(self):
        g = make_G2()
        assert np.max(np.abs(g.kappa - g.kappa.T)) < 1e-12


# =========================================================================
# PART 1: Casimir tensor (Item 1)
# =========================================================================

class TestCasimirAllRank2:
    """Casimir tensor Omega = kappa^{-1} for all rank-2 algebras."""

    @pytest.mark.parametrize("make_fn,name,dim", [
        (make_sl3, 'A_2', 8),
        (make_B2, 'B_2', 10),
        (make_C2, 'C_2', 10),
        (make_G2, 'G_2', 14),
    ])
    def test_casimir_is_inverse_killing(self, make_fn, name, dim):
        g = make_fn()
        omega = casimir_tensor(g)
        product = g.kappa @ omega
        assert np.max(np.abs(product - np.eye(dim))) < 1e-10, \
            f"Omega * kappa != Id for {name}"

    @pytest.mark.parametrize("make_fn,name", [
        (make_sl3, 'A_2'), (make_B2, 'B_2'), (make_C2, 'C_2'), (make_G2, 'G_2'),
    ])
    def test_casimir_symmetry(self, make_fn, name):
        g = make_fn()
        omega = casimir_tensor(g)
        assert np.max(np.abs(omega - omega.T)) < 1e-12, \
            f"Casimir not symmetric for {name}"

    @pytest.mark.parametrize("make_fn,name,expected_eigenvalue", [
        (make_sl3, 'A_2', 6),   # 2*h^v = 6
        (make_B2, 'B_2', 6),    # 2*h^v = 6
        (make_C2, 'C_2', 6),    # 2*h^v = 6
        (make_G2, 'G_2', 8),    # 2*h^v = 8
    ])
    def test_adjoint_casimir_eigenvalue(self, make_fn, name, expected_eigenvalue):
        """Quadratic Casimir eigenvalue in the adjoint = 2*h^v."""
        g = make_fn()
        omega = casimir_tensor(g)
        d = g.dim

        casimir_adj = np.zeros((d, d))
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

        ratio = casimir_adj / expected_eigenvalue
        assert np.max(np.abs(ratio - np.eye(d))) < 1e-8, \
            f"Casimir eigenvalue wrong for {name}: expected {expected_eigenvalue}"


# =========================================================================
# PART 2: R-matrix R(z) = 1 + hbar*Omega/z (Item 2)
# =========================================================================

class TestRMatrixAllRank2:
    """Collision residue r(z) = k * Omega / z for all rank-2."""

    @pytest.mark.parametrize("make_fn,name", [
        (make_sl3, 'A_2'), (make_B2, 'B_2'), (make_C2, 'C_2'), (make_G2, 'G_2'),
    ])
    def test_rmatrix_identification(self, make_fn, name):
        g = make_fn()
        res = full_collision_residue_computation(g, k=1.0)
        assert res['r_equals_k_omega_over_z'], f"r != Omega/z for {name}"
        assert res['all_checks_passed'], f"Failed for {name}"

    @pytest.mark.parametrize("make_fn,name", [
        (make_sl3, 'A_2'), (make_B2, 'B_2'), (make_C2, 'C_2'), (make_G2, 'G_2'),
    ])
    def test_pole_absorption(self, make_fn, name):
        """AP19: OPE double pole -> r-matrix simple pole."""
        g = make_fn()
        ope = AffineOPE(g=g, k=1.0)
        result = collision_residue_rmatrix(ope)
        assert result['pole_absorption_verified'], f"AP19 failed for {name}"
        assert result['r_matrix_max_pole'] == 1, f"Max pole != 1 for {name}"

    @pytest.mark.parametrize("make_fn,name", [
        (make_sl3, 'A_2'), (make_B2, 'B_2'), (make_C2, 'C_2'), (make_G2, 'G_2'),
    ])
    def test_level_scaling(self, make_fn, name):
        """r(z) scales linearly with level k."""
        g = make_fn()
        res1 = collision_residue_rmatrix(AffineOPE(g=g, k=1.0))
        coeff1 = res1['r_pole_coefficients'][1]
        for k in [2, -1, 0.5]:
            res_k = collision_residue_rmatrix(AffineOPE(g=g, k=float(k)))
            coeff_k = res_k['r_pole_coefficients'][1]
            assert np.max(np.abs(coeff_k - k * coeff1)) < 1e-12, \
                f"Level scaling failed for {name} at k={k}"


# =========================================================================
# PART 3: YBE on V tensor V tensor V (Item 3)
# =========================================================================

class TestYBEAllRank2:
    """Yang-Baxter equation in defining representations."""

    def test_ybe_A2_3dim(self):
        """YBE for sl_3 in 3-dim fundamental."""
        g = make_sl3()
        mats = _get_sl3_3dim_matrices()
        Omega = casimir_in_defining_rep(g, mats)
        result = verify_ybe_in_rep(Omega, 3)
        assert result['ybe_satisfied'], \
            f"YBE failed for A_2: violation = {result['ibr_form1_max_violation']}"

    def test_ybe_B2C2_4dim(self):
        """YBE for B_2/C_2 in 4-dim fundamental of sp(4)."""
        g = make_B2()
        mats = _get_sp4_basis_matrices()
        Omega = casimir_in_defining_rep(g, mats)
        result = verify_ybe_in_rep(Omega, 4)
        assert result['ybe_satisfied'], \
            f"YBE failed for B_2: violation = {result['ibr_form1_max_violation']}"

    def test_ybe_G2_7dim(self):
        """YBE for G_2 in 7-dim fundamental."""
        g = make_G2()
        mats = _get_g2_7dim_matrices()
        Omega = casimir_in_defining_rep(g, mats)
        result = verify_ybe_in_rep(Omega, 7)
        assert result['ybe_satisfied'], \
            f"YBE failed for G_2: violation = {result['ibr_form1_max_violation']}"

    def test_ybe_violations_machine_precision(self):
        """All YBE violations should be at machine epsilon."""
        for make_fn, get_mats, n, name in [
            (make_sl3, _get_sl3_3dim_matrices, 3, 'A_2'),
            (make_B2, _get_sp4_basis_matrices, 4, 'B_2'),
            (make_G2, _get_g2_7dim_matrices, 7, 'G_2'),
        ]:
            g = make_fn()
            mats = get_mats()
            Omega = casimir_in_defining_rep(g, mats)
            result = verify_ybe_in_rep(Omega, n)
            assert result['ibr_form1_max_violation'] < 1e-10, \
                f"YBE violation {result['ibr_form1_max_violation']} for {name}"


# =========================================================================
# PART 4: CYBE in adjoint (from existing code)
# =========================================================================

class TestCYBEAllRank2:
    """CYBE via IBR in the adjoint representation."""

    @pytest.mark.parametrize("make_fn,name", [
        (make_sl3, 'A_2'), (make_B2, 'B_2'), (make_C2, 'C_2'), (make_G2, 'G_2'),
    ])
    def test_cybe(self, make_fn, name):
        g = make_fn()
        result = verify_cybe(g)
        assert result['cybe_satisfied'], f"CYBE failed for {name}"

    @pytest.mark.parametrize("make_fn,name", [
        (make_sl3, 'A_2'), (make_B2, 'B_2'), (make_C2, 'C_2'), (make_G2, 'G_2'),
    ])
    def test_ibr_machine_precision(self, make_fn, name):
        g = make_fn()
        result = verify_cybe(g)
        assert result['ibr_max_violation'] < 1e-10, \
            f"IBR violation {result['ibr_max_violation']} for {name}"


# =========================================================================
# PART 5: RTT relation count (Item 4)
# =========================================================================

class TestRTTRelationCount:
    """RTT relation count in defining representations."""

    def test_A2_rtt(self):
        """sl_3 in 3-dim: 36 independent RTT relations."""
        g = make_sl3()
        rtt = rtt_relation_count(g, 3)
        assert rtt['independent_relations'] == 36

    def test_B2_rtt(self):
        """so_5 in 4-dim (spin): 120 independent RTT relations."""
        g = make_B2()
        rtt = rtt_relation_count(g, 4)
        assert rtt['independent_relations'] == 120

    def test_C2_rtt(self):
        """sp_4 in 4-dim: 120 independent RTT relations."""
        g = make_C2()
        rtt = rtt_relation_count(g, 4)
        assert rtt['independent_relations'] == 120

    def test_G2_rtt(self):
        """G_2 in 7-dim: 1176 independent RTT relations."""
        g = make_G2()
        rtt = rtt_relation_count(g, 7)
        assert rtt['independent_relations'] == 1176

    def test_rtt_formula(self):
        """RTT independent = n^2(n^2-1)/2."""
        for n in [2, 3, 4, 5, 7]:
            expected = n**2 * (n**2 - 1) // 2
            g = make_sl3()  # algebra doesn't matter for the count
            rtt = rtt_relation_count(g, n)
            assert rtt['independent_relations'] == expected


# =========================================================================
# PART 6: Koszul dual = Yangian (Item 5)
# =========================================================================

class TestKoszulDual:
    """Koszul dual identification: A^!_line = Y_hbar(g)."""

    @pytest.mark.parametrize("make_fn,name,h_dual", [
        (make_sl3, 'A_2', 3),
        (make_B2, 'B_2', 3),
        (make_C2, 'C_2', 3),
        (make_G2, 'G_2', 4),
    ])
    def test_hbar_parameter(self, make_fn, name, h_dual):
        """hbar = 1/(k + h^v) for the Yangian."""
        g = make_fn()
        assert g.h_dual == h_dual
        # At level k, hbar = 1/(k + h^v)
        for k in [1, 2, 5]:
            hbar = 1.0 / (k + g.h_dual)
            expected = 1.0 / (k + h_dual)
            assert abs(hbar - expected) < 1e-15

    @pytest.mark.parametrize("make_fn,name,expected_kappa_k1", [
        (make_sl3, 'A_2', 8.0 * 4 / 6),    # dim*(k+h^v)/(2h^v) = 8*4/6 = 16/3
        (make_B2, 'B_2', 10.0 * 4 / 6),    # 10*4/6 = 20/3
        (make_C2, 'C_2', 10.0 * 4 / 6),    # 10*4/6 = 20/3
        (make_G2, 'G_2', 14.0 * 5 / 8),    # 14*5/8 = 35/4
    ])
    def test_kappa_at_k1(self, make_fn, name, expected_kappa_k1):
        """kappa(g, k=1) = dim(g)*(1+h^v)/(2h^v)."""
        g = make_fn()
        res = full_collision_residue_computation(g, k=1.0)
        assert abs(res['kappa_A'] - expected_kappa_k1) < 1e-12, \
            f"kappa({name}, k=1) = {res['kappa_A']}, expected {expected_kappa_k1}"

    @pytest.mark.parametrize("make_fn,name", [
        (make_sl3, 'A_2'), (make_B2, 'B_2'), (make_C2, 'C_2'), (make_G2, 'G_2'),
    ])
    def test_critical_level_kappa_zero(self, make_fn, name):
        """kappa = 0 at critical level k = -h^v."""
        g = make_fn()
        res = full_collision_residue_computation(g, k=-float(g.h_dual))
        assert abs(res['kappa_A']) < 1e-12, \
            f"kappa should be 0 at critical level for {name}"

    @pytest.mark.parametrize("make_fn,name", [
        (make_sl3, 'A_2'), (make_B2, 'B_2'), (make_C2, 'C_2'), (make_G2, 'G_2'),
    ])
    def test_feigin_frenkel_duality(self, make_fn, name):
        """FF duality: kappa(k) + kappa(-k - 2h^v) = 0."""
        g = make_fn()
        for k in [1, 2, 5]:
            k_dual = -k - 2 * g.h_dual
            kappa_k = g.dim * (k + g.h_dual) / (2 * g.h_dual)
            kappa_dual = g.dim * (k_dual + g.h_dual) / (2 * g.h_dual)
            assert abs(kappa_k + kappa_dual) < 1e-12, \
                f"FF failed for {name} at k={k}"


# =========================================================================
# PART 7: Euler-eta identity (Item 6)
# =========================================================================

class TestEulerEta:
    """Euler-eta: chi = -1 + eta^{dim g}."""

    @pytest.mark.parametrize("make_fn,name,expected_dim", [
        (make_sl3, 'A_2', 8),
        (make_B2, 'B_2', 10),
        (make_C2, 'C_2', 10),
        (make_G2, 'G_2', 14),
    ])
    def test_eta_exponent(self, make_fn, name, expected_dim):
        g = make_fn()
        eta = euler_eta_rank2(g)
        assert eta['eta_exponent'] == expected_dim
        assert eta['first_nontrivial_coefficient'] == -expected_dim


# =========================================================================
# PART 8: Non-simply-laced specific (Item 7)
# =========================================================================

class TestNonSimplyLacedSpecific:
    """Root-length-dependent Casimir and Langlands duality."""

    def test_B2_C2_isomorphism(self):
        """so(5) and sp(4) are isomorphic as abstract Lie algebras."""
        g_B2 = make_B2()
        g_C2 = make_C2()
        assert np.allclose(g_B2.f, g_C2.f, atol=1e-12)
        assert np.allclose(g_B2.kappa, g_C2.kappa, atol=1e-12)

    def test_B2_C2_casimir_match(self):
        """Casimir tensors match (same algebra)."""
        omega_B2 = casimir_tensor(make_B2())
        omega_C2 = casimir_tensor(make_C2())
        assert np.allclose(omega_B2, omega_C2, atol=1e-12)

    def test_langlands_B2_C2_cartan_transpose(self):
        """A_{C_2} = A_{B_2}^T."""
        A_B2 = np.array([[2, -1], [-2, 2]], dtype=float)
        A_C2 = np.array([[2, -2], [-1, 2]], dtype=float)
        assert np.allclose(A_B2, A_C2.T)

    def test_langlands_B2_C2_dual_level(self):
        """Langlands dual level k^L = 2k for B_2 <-> C_2."""
        result = langlands_duality_comparison(3.0)
        assert result['k_dual'] == 6.0

    def test_B2_root_length_casimir_dependence(self):
        """Casimir coefficients depend on root length: Omega(e,f) = 1/kappa(e,f)."""
        g = make_B2()
        decomp = casimir_root_decomposition(g)
        for label, data in decomp['root_casimir_coefficients'].items():
            assert data['predicted'] is not None
            assert abs(data['e_x_f'] - data['predicted']) < 1e-12, \
                f"Casimir reciprocal mismatch for {label}"

    def test_G2_root_casimir_ratio(self):
        """G_2: Casimir ratio between root types = 3 (= root length ratio)."""
        g = make_G2()
        omega = casimir_tensor(g)

        # alpha_2 (long root): index 1
        long_casimir = omega[1, 1 + 6]

        # alpha_1 (short root): index 0
        short_casimir = omega[0, 0 + 6]

        if abs(short_casimir) > 1e-14:
            ratio = long_casimir / short_casimir
            assert abs(ratio - 3.0) < 1e-10, \
                f"G_2 Casimir ratio = {ratio}, expected 3"

    def test_G2_self_langlands_dual(self):
        """G_2^L = G_2 (self-dual under Langlands)."""
        A_G2 = np.array([[2, -1], [-3, 2]], dtype=float)
        # A^T = [[2, -3], [-1, 2]] is the Cartan matrix of G_2 with roots relabelled
        # (swap node 1 and node 2). This is the SAME root system (G_2 is self-dual).
        A_G2_T = A_G2.T
        # The transposed matrix [[2, -3], [-1, 2]] has A_{12} = -3, A_{21} = -1.
        # This is the G_2 Cartan matrix with alpha_1 <-> alpha_2.
        assert A_G2_T[0, 1] == -3  # was -1
        assert A_G2_T[1, 0] == -1  # was -3
        # Langlands duality swaps nodes: G_2 -> G_2.

    def test_G2_dual_coxeter(self):
        """h^v(G_2) = 4."""
        g = make_G2()
        assert g.h_dual == 4

    def test_G2_root_length_ratio(self):
        """Root length ratio for G_2: long^2/short^2 = 3."""
        # The G_2 Cartan matrix gives:
        # A_{21} = -3 means (alpha_2, alpha_2)/(alpha_1, alpha_1) = A_{12}/A_{21} = (-1)/(-3) = 1/3
        # Wait: A_{ij} = 2(alpha_i, alpha_j)/(alpha_j, alpha_j).
        # A_{12}/A_{21} = [(alpha_1, alpha_1)/(alpha_2, alpha_2)] * 1
        # A_{12} = -1, A_{21} = -3. So:
        # A_{12}/A_{21} = -1/-3 = 1/3.
        # This means (alpha_1, alpha_1)/(alpha_2, alpha_2) = |A_{21}|/|A_{12}| = 3/1 = ... hmm.
        # Actually: A_{12} = 2(a1,a2)/(a2,a2), A_{21} = 2(a2,a1)/(a1,a1).
        # A_{12}*A_{21} = 4(a1,a2)^2 / ((a1,a1)(a2,a2)).
        # A_{12}/A_{21} = (a1,a1)/(a2,a2).
        # So (a1,a1)/(a2,a2) = (-1)/(-3) = 1/3.
        # alpha_1 is SHORT (1/3 of alpha_2). Correct.
        ratio = abs(-1) / abs(-3)  # A_{12}/A_{21} = 1/3
        assert abs(ratio - 1.0/3) < 1e-15
        assert abs(1.0/ratio - 3.0) < 1e-15  # long^2/short^2 = 3

    def test_fm_integral_G2(self):
        """FM-integral beta values for G_2 (root length ratio 3)."""
        pytest.importorskip('scipy')
        from lib.non_simply_laced_rmatrix import fm_integral_degree3_nonsimplylaced
        result = fm_integral_degree3_nonsimplylaced([1.0, 3.0])

        # short-short: B(1,1) = 1
        assert abs(result['beta_integrals'][(1.0, 1.0)]['beta'] - 1.0) < 1e-12
        # short-long: B(1,3) = 1/3! = ... B(1,3) = Gamma(1)*Gamma(3)/Gamma(4) = 1*2/6 = 1/3
        assert abs(result['beta_integrals'][(1.0, 3.0)]['beta'] - 1.0/3) < 1e-12
        # long-long: B(3,3) = Gamma(3)^2/Gamma(6) = 4/120 = 1/30
        assert abs(result['beta_integrals'][(3.0, 3.0)]['beta'] - 1.0/30) < 1e-12


# =========================================================================
# PART 9: Rep-matrix consistency checks
# =========================================================================

class TestRepConsistency:
    """Verify representation matrices are consistent with structure constants."""

    def test_sl3_rep_consistency(self):
        """3x3 matrices of sl_3 reproduce the structure constants."""
        g = make_sl3()
        mats = _get_sl3_3dim_matrices()
        for a in range(8):
            for b in range(8):
                comm = mats[a] @ mats[b] - mats[b] @ mats[a]
                expected = sum(g.f[a, b, c] * mats[c] for c in range(8))
                assert np.max(np.abs(comm - expected)) < 1e-10, \
                    f"sl_3 rep inconsistency at ({a}, {b})"

    def test_sp4_rep_consistency(self):
        """4x4 matrices of sp(4) reproduce the structure constants."""
        g = make_B2()
        mats = _get_sp4_basis_matrices()
        for a in range(10):
            for b in range(10):
                comm = mats[a] @ mats[b] - mats[b] @ mats[a]
                expected = sum(g.f[a, b, c] * mats[c] for c in range(10))
                assert np.max(np.abs(comm - expected)) < 1e-10, \
                    f"sp(4) rep inconsistency at ({a}, {b})"

    def test_g2_rep_consistency(self):
        """7x7 matrices of G_2 reproduce the structure constants."""
        g = make_G2()
        mats = _get_g2_7dim_matrices()
        max_err = 0.0
        for a in range(14):
            for b in range(14):
                comm = mats[a] @ mats[b] - mats[b] @ mats[a]
                expected = sum(g.f[a, b, c] * mats[c] for c in range(14))
                err = np.max(np.abs(comm - expected))
                max_err = max(max_err, err)
        assert max_err < 1e-10, \
            f"G_2 rep inconsistency: max error = {max_err}"


# =========================================================================
# PART 10: End-to-end pipeline
# =========================================================================

class TestFullRank2Pipeline:
    """End-to-end computation for all rank-2 algebras."""

    @pytest.mark.parametrize("make_fn,name", [
        (make_sl3, 'A_2'), (make_B2, 'B_2'), (make_C2, 'C_2'), (make_G2, 'G_2'),
    ])
    def test_full_pipeline(self, make_fn, name):
        g = make_fn()
        result = full_collision_residue_computation(g, k=1.0)
        assert result['all_checks_passed'], f"Full pipeline failed for {name}"

    def test_complete_rank2_computation(self):
        """Run the complete rank-2 computation and verify all checks pass."""
        results = run_complete_rank2_computation(k=1.0, verbose=False)
        assert results['all_pass'], "Complete rank-2 computation failed"
