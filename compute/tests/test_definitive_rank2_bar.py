r"""Definitive tests for ordered E_1 bar complex data: ALL rank-2 simple Lie algebras.

Tests the complete package for {A_2 (sl_3), B_2 (so_5), C_2 (sp_4), G_2}:
  (1) Casimir Omega in the defining representation.
  (2) R(z) = 1 + k*Omega/z on V tensor V.
  (3) CYBE verification on V^{x3} (numerical at 5 values of k).
  (4) RTT relation count.
  (5) Euler-eta: chi = -1 + eta^d with d = 8, 10, 10, 14.
  (6) Non-simply-laced: root-length Casimir coefficients, Langlands duality.
  (7) DS reduction to W-algebras. Depth gaps: 4, 6, 6, 10.

References:
  Vol II, ordered_associative_chiral_kd_core.tex (ordered bar complex)
  Vol II, rosetta_stone.tex (affine CS example)
  Vol II, w-algebras-stable.tex (exponent-depth correspondence)
  AP19: The bar kernel absorbs a pole.
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.non_simply_laced_rmatrix import (
    make_B2, make_C2, make_G2,
    casimir_root_decomposition,
    langlands_duality_comparison,
    ds_reduction_rank2,
    verify_ybe_numerical_multi_k,
    verify_ybe_in_rep,
    rtt_relation_count,
    euler_eta_rank2,
    casimir_in_defining_rep,
    _get_sl3_3dim_matrices,
    _get_sp4_basis_matrices,
    _get_g2_7dim_matrices,
    run_definitive_rank2,
)
from lib.collision_residue_rmatrix import (
    make_sl3,
    verify_jacobi, verify_killing_invariance, verify_antisymmetry,
    casimir_tensor,
    verify_cybe, full_collision_residue_computation,
)


# ==========================================================================
# PART 1: Lie algebra axioms for all rank-2 types
# ==========================================================================

class TestLieAlgebraAxioms:
    """Verify Lie algebra data for all rank-2 simple Lie algebras."""

    @pytest.fixture(params=[
        ('A_2', make_sl3, 8, 2, 3),
        ('B_2', make_B2, 10, 2, 3),
        ('C_2', make_C2, 10, 2, 3),
        ('G_2', make_G2, 14, 2, 4),
    ], ids=['A2', 'B2', 'C2', 'G2'])
    def algebra(self, request):
        name, maker, dim, rank, h_dual = request.param
        g = maker()
        return name, g, dim, rank, h_dual

    def test_jacobi(self, algebra):
        name, g, _, _, _ = algebra
        assert verify_jacobi(g), f"{name} fails Jacobi identity"

    def test_antisymmetry(self, algebra):
        name, g, _, _, _ = algebra
        assert verify_antisymmetry(g), f"{name} structure constants not antisymmetric"

    def test_killing_invariance(self, algebra):
        name, g, _, _, _ = algebra
        assert verify_killing_invariance(g), f"{name} Killing form not ad-invariant"

    def test_dimensions(self, algebra):
        name, g, dim, rank, h_dual = algebra
        assert g.dim == dim, f"{name}: dim = {g.dim} != {dim}"
        assert g.rank == rank, f"{name}: rank = {g.rank} != {rank}"
        assert g.h_dual == h_dual, f"{name}: h^v = {g.h_dual} != {h_dual}"


# ==========================================================================
# PART 2: Casimir tensor and CYBE in defining representation
# ==========================================================================

class TestCasimirAndCYBE:
    """Verify Casimir tensors and classical Yang-Baxter equation."""

    def _get_rep_data(self, name):
        """Return (g, rep_dim, rep_mats) for each algebra."""
        if name == 'A_2':
            return make_sl3(), 3, _get_sl3_3dim_matrices()
        elif name in ('B_2', 'C_2'):
            g = make_B2() if name == 'B_2' else make_C2()
            return g, 4, _get_sp4_basis_matrices()
        elif name == 'G_2':
            return make_G2(), 7, _get_g2_7dim_matrices()

    @pytest.mark.parametrize('name', ['A_2', 'B_2', 'C_2', 'G_2'])
    def test_casimir_nondegenerate(self, name):
        """Casimir tensor kappa^{-1} is non-degenerate."""
        g, _, _ = self._get_rep_data(name)
        omega = casimir_tensor(g)
        assert abs(np.linalg.det(omega)) > 1e-10, f"{name}: Casimir tensor is degenerate"

    @pytest.mark.parametrize('name', ['A_2', 'B_2', 'C_2', 'G_2'])
    def test_cybe_abstract(self, name):
        """CYBE verified in abstract (Lie algebra) form."""
        g, _, _ = self._get_rep_data(name)
        res = full_collision_residue_computation(g, 1.0)
        assert res['cybe_satisfied'], f"{name}: abstract CYBE failed"

    @pytest.mark.parametrize('name', ['A_2', 'B_2', 'C_2', 'G_2'])
    def test_ibr_in_rep(self, name):
        """Infinitesimal braid relation [O12, O13 + O23] = 0 in representation."""
        g, rep_dim, rep_mats = self._get_rep_data(name)
        Omega_rep = casimir_in_defining_rep(g, rep_mats)
        ibr = verify_ybe_in_rep(Omega_rep, rep_dim)
        assert ibr['ybe_satisfied'], \
            f"{name}: IBR fails in {rep_dim}-dim rep, max violation = {ibr['ibr_form1_max_violation']:.2e}"

    @pytest.mark.parametrize('name', ['A_2', 'B_2', 'C_2', 'G_2'])
    def test_cybe_numerical_multi_k(self, name):
        """CYBE verified numerically at 5 values of k in representation."""
        g, rep_dim, rep_mats = self._get_rep_data(name)
        Omega_rep = casimir_in_defining_rep(g, rep_mats)
        k_values = [0.5, 1.0, 2.0, 3.0, 5.0]
        result = verify_ybe_numerical_multi_k(Omega_rep, rep_dim, k_values)
        assert result['ibr_satisfied'], f"{name}: IBR failed"
        assert result['all_passed'], \
            f"{name}: CYBE failed at some k-value: {result['results_per_k']}"


# ==========================================================================
# PART 3: RTT relation counts
# ==========================================================================

class TestRTTRelations:
    """Verify RTT relation counts match expected values."""

    @pytest.mark.parametrize('name,rep_dim,expected', [
        ('A_2', 3, 36),     # 3^2 * (3^2-1)/2 = 9*8/2 = 36
        ('B_2', 4, 120),    # 4^2 * (4^2-1)/2 = 16*15/2 = 120
        ('C_2', 4, 120),    # same as B_2
        ('G_2', 7, 1176),   # 7^2 * (7^2-1)/2 = 49*48/2 = 1176
    ])
    def test_rtt_count(self, name, rep_dim, expected):
        if name == 'A_2':
            g = make_sl3()
        elif name == 'B_2':
            g = make_B2()
        elif name == 'C_2':
            g = make_C2()
        elif name == 'G_2':
            g = make_G2()
        rtt = rtt_relation_count(g, rep_dim)
        assert rtt['independent_relations'] == expected, \
            f"{name}: RTT count = {rtt['independent_relations']} != {expected}"


# ==========================================================================
# PART 4: Euler-eta invariants
# ==========================================================================

class TestEulerEta:
    """Verify Euler-eta chi = -1 + eta^{dim g}."""

    @pytest.mark.parametrize('name,expected_dim', [
        ('A_2', 8),    # dim(sl_3) = 8
        ('B_2', 10),   # dim(so_5) = 10
        ('C_2', 10),   # dim(sp_4) = 10
        ('G_2', 14),   # dim(G_2) = 14
    ])
    def test_eta_exponent(self, name, expected_dim):
        if name == 'A_2':
            g = make_sl3()
        elif name == 'B_2':
            g = make_B2()
        elif name == 'C_2':
            g = make_C2()
        elif name == 'G_2':
            g = make_G2()
        eta = euler_eta_rank2(g)
        assert eta['eta_exponent'] == expected_dim, \
            f"{name}: eta exponent = {eta['eta_exponent']} != {expected_dim}"
        assert eta['euler_eta'] == f'-1 + eta^{{{expected_dim}}}', \
            f"{name}: wrong formula: {eta['euler_eta']}"


# ==========================================================================
# PART 5: DS reduction and depth gaps
# ==========================================================================

class TestDSReduction:
    """Verify DS reduction data for all rank-2 types."""

    @pytest.fixture
    def ds_data(self):
        return ds_reduction_rank2()

    @pytest.mark.parametrize('name,exponents,w_name,d_gap', [
        ('A_2', [1, 2], 'W_3', 4),
        ('B_2', [1, 3], 'W(B_2)', 6),
        ('C_2', [1, 3], 'W(C_2)', 6),
        ('G_2', [1, 5], 'W(G_2)', 10),
    ])
    def test_ds_reduction(self, ds_data, name, exponents, w_name, d_gap):
        entry = ds_data[name]
        assert entry['exponents'] == exponents, \
            f"{name}: exponents = {entry['exponents']} != {exponents}"
        assert entry['w_algebra'] == w_name, \
            f"{name}: W-algebra = {entry['w_algebra']} != {w_name}"
        assert entry['d_gap'] == d_gap, \
            f"{name}: d_gap = {entry['d_gap']} != {d_gap}"

    @pytest.mark.parametrize('name', ['A_2', 'B_2', 'C_2', 'G_2'])
    def test_coxeter_formula(self, ds_data, name):
        """d_gap = 2*(h-1) where h is the Coxeter number."""
        entry = ds_data[name]
        assert entry['coxeter_formula_agrees'], \
            f"{name}: d_gap = {entry['d_gap']} != 2*(h-1) = {entry['d_gap_from_coxeter']}"

    @pytest.mark.parametrize('name', ['A_2', 'B_2', 'C_2', 'G_2'])
    def test_exponent_sum(self, ds_data, name):
        """sum(exponents) = |Phi+|."""
        data = ds_data[name]['lie_data']
        assert sum(data['exponents']) == data['num_positive_roots'], \
            f"{name}: sum(exp) = {sum(data['exponents'])} != |Phi+| = {data['num_positive_roots']}"


# ==========================================================================
# PART 6: Non-simply-laced structure
# ==========================================================================

class TestNonSimplyLaced:
    """Verify root-length dependent Casimir coefficients and Langlands duality."""

    def test_b2_c2_casimir_match(self):
        """B_2 and C_2 have the same abstract Casimir (isomorphic Lie algebras)."""
        g_B2 = make_B2()
        g_C2 = make_C2()
        omega_B2 = casimir_tensor(g_B2)
        omega_C2 = casimir_tensor(g_C2)
        assert np.allclose(omega_B2, omega_C2, atol=1e-12), \
            "B_2 and C_2 should have identical abstract Casimir tensors"

    def test_langlands_cartan_transpose(self):
        """A_{C_2} = A_{B_2}^T (Langlands duality)."""
        lang = langlands_duality_comparison(1.0)
        assert lang['cartan_transpose'], "Cartan matrices not transposes"

    def test_langlands_symmetriser_swap(self):
        """D_{C_2} swaps entries of D_{B_2}."""
        lang = langlands_duality_comparison(1.0)
        assert lang['symmetriser_swap'], "Symmetrisers not swapped"

    def test_b2_root_length_coefficients(self):
        """B_2 Casimir has root-length dependent coefficients."""
        g = make_B2()
        decomp = casimir_root_decomposition(g)
        # Long roots (alpha_1, alpha_1+2*alpha_2): kappa(e,f) should differ from short
        # In the sp(4) representation:
        # alpha_1 (C_2 short/B_2 long): kappa = 2
        # alpha_2 (C_2 long/B_2 short): kappa = 1
        for label, data in decomp['root_casimir_coefficients'].items():
            assert abs(data['kappa_ef']) > 0.5, \
                f"Root {label}: kappa(e,f) = {data['kappa_ef']} too small"
            # Casimir coefficient = 1/kappa(e,f) (predicted)
            if data['predicted'] is not None:
                assert abs(data['e_x_f'] - data['predicted']) < 1e-6, \
                    f"Root {label}: Omega coeff {data['e_x_f']} != 1/kappa = {data['predicted']}"

    def test_g2_root_length_ratio(self):
        """G_2 has root length ratio long^2/short^2 = 3."""
        g = make_G2()
        omega = casimir_tensor(g)
        # alpha_2 (long root): kappa(e,f) = 1 (normalised)
        # alpha_1 (short root): kappa(e,f) = 3 (from 2/(alpha,alpha) = 2/(2/3) = 3)
        kappa_long = g.kappa[1, 7]   # e_{a2}, f_{a2}
        kappa_short = g.kappa[0, 6]  # e_{a1}, f_{a1}
        # Ratio should be 3 (short/long in kappa = 3, because kappa ~ 1/|alpha|^2)
        ratio = kappa_short / kappa_long
        assert abs(ratio - 3.0) < 0.1, \
            f"G_2 root length ratio: kappa(short)/kappa(long) = {ratio} != 3"

    def test_g2_self_dual(self):
        """G_2 is Langlands self-dual."""
        ds = ds_reduction_rank2()
        assert ds['G_2']['lie_data']['langlands_dual'] == 'G_2', \
            "G_2 should be Langlands self-dual"


# ==========================================================================
# PART 7: Full pipeline integration test
# ==========================================================================

class TestFullPipeline:
    """Integration test: run the complete definitive computation."""

    def test_all_rank2_pass(self):
        """The complete rank-2 definitive computation passes all checks."""
        results = run_definitive_rank2(verbose=False)
        assert results['all_pass'], "Not all rank-2 checks passed"

    def test_all_collision_residues(self):
        """Full collision residue pipeline passes for all algebras."""
        results = run_definitive_rank2(verbose=False)
        for name in ['A_2', 'B_2', 'C_2', 'G_2']:
            r = results[name]['full_collision_residue']
            assert r['all_checks_passed'], f"{name}: collision residue pipeline failed"
            assert r['cybe_satisfied'], f"{name}: abstract CYBE failed"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
