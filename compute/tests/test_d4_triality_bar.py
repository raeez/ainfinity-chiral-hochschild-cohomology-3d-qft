r"""Tests for the D_4 triality-equivariant ordered bar complex.

Verifies all six computational aspects of V_k(so_8):
  (1) D_4 Lie algebra data consistency (dim, h^v, exponents, repeated exponent)
  (2) so(8) structure constants: Jacobi, antisymmetry, ad-invariance
  (3) Three 8-dim representations: brackets, Casimir, triality equality
  (4) Split Casimir: spectra match across triality triple
  (5) Triality-invariant bar complex: Poincare series
  (6) DS reduction: depth gap = 10 (e_r = 5, not 3)

Cross-references:
  Vol II: rosetta_stone.tex, comp:lattice-voa-D4-ordered-bar
  Vol II: ordered_associative_chiral_kd_core.tex
  Vol II: dg_shifted_factorization_bridge.tex (D_4 star discussion)
"""

import pytest
import numpy as np
import sys
import os
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.d4_triality_bar import (
    _D4_DATA,
    verify_d4_data,
    make_so8,
    _vector_rep_so8,
    _spinor_reps_so8,
    _casimir_in_rep,
    _split_casimir_in_rep,
    compute_casimir_all_reps,
    compute_split_casimir_all_reps,
    compute_r_matrix,
    verify_triality_invariance,
    compute_koszul_duals,
    tensor_product_decompositions,
    triality_invariant_bar_complex,
    ds_reduction_d4,
    verify_representations,
    complete_d4_triality_bar,
)

from lib.collision_residue_rmatrix import (
    verify_jacobi, verify_killing_invariance, verify_antisymmetry,
)


# =========================================================================
# PART 1: D_4 LIE ALGEBRA DATA
# =========================================================================

class TestD4Data:
    """Verify the D_4 Lie algebra data."""

    def test_dimension(self):
        assert _D4_DATA['dim'] == 28

    def test_rank(self):
        assert _D4_DATA['rank'] == 4

    def test_dual_coxeter(self):
        assert _D4_DATA['h_dual'] == 6

    def test_exponents(self):
        assert _D4_DATA['exponents'] == [1, 3, 3, 5]

    def test_repeated_exponent(self):
        """D_4 is the UNIQUE simple Lie algebra with a repeated exponent."""
        exps = _D4_DATA['exponents']
        assert exps.count(3) == 2, "D_4 must have exponent 3 with multiplicity 2"

    def test_exponent_sum(self):
        """sum(exponents) = |Phi^+|."""
        assert sum(_D4_DATA['exponents']) == _D4_DATA['num_positive_roots']

    def test_dim_root_decomposition(self):
        """dim = rank + 2|Phi^+|."""
        assert _D4_DATA['dim'] == _D4_DATA['rank'] + 2 * _D4_DATA['num_positive_roots']

    def test_coxeter_from_exponent(self):
        """h = 1 + e_r."""
        assert _D4_DATA['h_coxeter'] == 1 + max(_D4_DATA['exponents'])

    def test_all_consistency_checks(self):
        result = verify_d4_data()
        assert result['checks']['all_passed'], (
            f"D_4 consistency checks failed: "
            f"{[k for k, v in result['checks'].items() if not v and k != 'all_passed']}"
        )


# =========================================================================
# PART 2: so(8) STRUCTURE CONSTANTS
# =========================================================================

class TestSO8LieAlgebra:
    """Verify the so(8) Lie algebra construction."""

    def test_dimension(self):
        g = make_so8()
        assert g.dim == 28
        assert g.rank == 4
        assert g.h_dual == 6

    def test_jacobi_identity(self):
        g = make_so8()
        assert verify_jacobi(g), "so(8) fails Jacobi identity"

    def test_antisymmetry(self):
        g = make_so8()
        assert verify_antisymmetry(g), "so(8) structure constants not antisymmetric"

    def test_killing_invariance(self):
        g = make_so8()
        assert verify_killing_invariance(g), "so(8) Killing form not ad-invariant"

    def test_killing_form_positive_definite(self):
        """Our normalised Killing form should be the identity (orthonormal basis)."""
        g = make_so8()
        assert np.allclose(g.kappa, np.eye(28)), "Killing form should be identity"


# =========================================================================
# PART 3: THREE 8-DIM REPRESENTATIONS
# =========================================================================

class TestRepresentations:
    """Verify the three 8-dimensional representations of so(8)."""

    def test_vector_rep_bracket(self):
        """8_v preserves the Lie bracket."""
        g = make_so8()
        rho = _vector_rep_so8()
        max_err = 0
        for a in range(28):
            for b in range(a+1, 28):
                comm = rho[a] @ rho[b] - rho[b] @ rho[a]
                expected = np.zeros((8, 8), dtype=rho.dtype)
                for c in range(28):
                    if abs(g.f[a, b, c]) > 1e-15:
                        expected += g.f[a, b, c] * rho[c]
                max_err = max(max_err, np.max(np.abs(comm - expected)))
        assert max_err < 1e-12, f"8_v bracket error: {max_err}"

    def test_spinor_rep_bracket(self):
        """8_s preserves the Lie bracket."""
        g = make_so8()
        rho_s, _ = _spinor_reps_so8()
        max_err = 0
        for a in range(28):
            for b in range(a+1, 28):
                comm = rho_s[a] @ rho_s[b] - rho_s[b] @ rho_s[a]
                expected = np.zeros((8, 8), dtype=rho_s.dtype)
                for c in range(28):
                    if abs(g.f[a, b, c]) > 1e-15:
                        expected += g.f[a, b, c] * rho_s[c]
                max_err = max(max_err, np.max(np.abs(comm - expected)))
        assert max_err < 1e-12, f"8_s bracket error: {max_err}"

    def test_conjugate_spinor_rep_bracket(self):
        """8_c preserves the Lie bracket."""
        g = make_so8()
        _, rho_c = _spinor_reps_so8()
        max_err = 0
        for a in range(28):
            for b in range(a+1, 28):
                comm = rho_c[a] @ rho_c[b] - rho_c[b] @ rho_c[a]
                expected = np.zeros((8, 8), dtype=rho_c.dtype)
                for c in range(28):
                    if abs(g.f[a, b, c]) > 1e-15:
                        expected += g.f[a, b, c] * rho_c[c]
                max_err = max(max_err, np.max(np.abs(comm - expected)))
        assert max_err < 1e-12, f"8_c bracket error: {max_err}"

    def test_casimir_triality_equality(self):
        """C_2(8_v) = C_2(8_s) = C_2(8_c): the triality miracle."""
        result = compute_casimir_all_reps()
        assert result['triality_casimir_equality'], (
            f"Casimir eigenvalues differ: "
            f"8_v={result['eigenvalue_vector']}, "
            f"8_s={result['eigenvalue_spinor']}, "
            f"8_c={result['eigenvalue_conjugate_spinor']}"
        )

    def test_casimir_is_scalar(self):
        """Casimir is proportional to identity (irreducibility)."""
        result = compute_casimir_all_reps()
        assert result['is_scalar_vector'], "C_2(8_v) not scalar"
        assert result['is_scalar_spinor'], "C_2(8_s) not scalar"
        assert result['is_scalar_conjugate_spinor'], "C_2(8_c) not scalar"

    def test_casimir_value(self):
        """C_2 = -7 for all three reps (with our normalisation)."""
        result = compute_casimir_all_reps()
        assert abs(result['eigenvalue_vector'] - (-7.0)) < 1e-10


# =========================================================================
# PART 4: SPLIT CASIMIR AND TRIALITY
# =========================================================================

class TestSplitCasimir:
    """Verify the split Casimir and triality at the representation level."""

    def test_spectra_match(self):
        """Split Casimir spectra are identical across all three reps."""
        data = compute_split_casimir_all_reps()
        assert data['triality_spectra_all_match'], (
            "Split Casimir spectra differ across triality triple"
        )

    def test_spectrum_decomposition(self):
        """Spectrum matches tensor product decomposition 8 x 8 = 1 + 28 + 35."""
        data = compute_split_casimir_all_reps()
        eigs = data['eigenvalues_vector']
        # Count multiplicities
        mult_neg1 = np.sum(np.abs(eigs - (-1.0)) < 1e-6)
        mult_pos1 = np.sum(np.abs(eigs - 1.0) < 1e-6)
        mult_7 = np.sum(np.abs(eigs - 7.0) < 1e-6)
        assert mult_neg1 == 35, f"Expected mult(-1) = 35 (35_v), got {mult_neg1}"
        assert mult_pos1 == 28, f"Expected mult(+1) = 28 (adjoint), got {mult_pos1}"
        assert mult_7 == 1, f"Expected mult(+7) = 1 (trivial), got {mult_7}"

    def test_trace_zero(self):
        """Tr(Omega) = 0 (from tracelessness of generators)."""
        data = compute_split_casimir_all_reps()
        assert abs(data['trace_vector']) < 1e-10
        assert abs(data['trace_spinor']) < 1e-10
        assert abs(data['trace_conjugate_spinor']) < 1e-10


# =========================================================================
# PART 5: TRIALITY INVARIANCE
# =========================================================================

class TestTriality:
    """Verify triality invariance of the bar complex."""

    def test_abstract_casimir_invariance(self):
        result = verify_triality_invariance()
        assert result['abstract_casimir_invariance']

    def test_representation_spectra_match(self):
        result = verify_triality_invariance()
        assert result['representation_spectra_match']

    def test_R_sigma_equals_R(self):
        """R^sigma = R for all sigma in S_3."""
        result = verify_triality_invariance()
        assert result['R_sigma_equals_R']


# =========================================================================
# PART 6: BAR COMPLEX AND POINCARE SERIES
# =========================================================================

class TestBarComplex:
    """Verify the triality-invariant ordered bar complex."""

    def test_shadow_class_L(self):
        bc = triality_invariant_bar_complex(k=1)
        assert bc['shadow_class'] == 'L'
        assert bc['shadow_depth'] == 3

    def test_curvature_complementarity(self):
        bc = triality_invariant_bar_complex(k=1)
        assert bc['complementarity_verified']

    def test_curvature_value(self):
        """kappa(V_1(so_8)) = 28*(1+6)/(2*6) = 196/12 = 49/3."""
        bc = triality_invariant_bar_complex(k=1)
        assert bc['kappa'] == Fraction(49, 3)

    def test_poincare_n0(self):
        bc = triality_invariant_bar_complex(k=1)
        assert bc['poincare_S3_invariant'][0] == 1

    def test_poincare_n1(self):
        """At arity 1, dim = (28 + 3*14 + 2*7)/6 = 84/6 = 14."""
        bc = triality_invariant_bar_complex(k=1)
        assert bc['poincare_S3_invariant'][1] == 14

    def test_poincare_n2(self):
        """At arity 2, dim = (784 + 588 + 98)/6 = 1470/6 = 245."""
        bc = triality_invariant_bar_complex(k=1)
        assert bc['poincare_S3_invariant'][2] == 245

    def test_fixed_subalgebra_transposition(self):
        """Fixed subalgebra of a transposition is so(7), dim 21."""
        bc = triality_invariant_bar_complex(k=1)
        name, dim = bc['fixed_subalgebras']['transposition']
        assert dim == 21
        assert 'so(7)' in name

    def test_fixed_subalgebra_3cycle(self):
        """Fixed subalgebra of a 3-cycle is G_2, dim 14."""
        bc = triality_invariant_bar_complex(k=1)
        name, dim = bc['fixed_subalgebras']['3-cycle']
        assert dim == 14
        assert 'G_2' in name


# =========================================================================
# PART 7: DS REDUCTION AND DEPTH
# =========================================================================

class TestDSReduction:
    """Verify DS reduction to W(D_4)."""

    def test_exponents(self):
        ds = ds_reduction_d4()
        assert ds['exponents'] == [1, 3, 3, 5]

    def test_largest_exponent(self):
        ds = ds_reduction_d4()
        assert ds['largest_exponent'] == 5

    def test_generator_weights(self):
        ds = ds_reduction_d4()
        assert ds['generator_weights'] == [2, 4, 4, 6]

    def test_depth_gap(self):
        """d_gap = 2*e_r = 2*5 = 10."""
        ds = ds_reduction_d4()
        assert ds['d_gap'] == 10

    def test_highest_spin(self):
        ds = ds_reduction_d4()
        assert ds['highest_spin'] == 6

    def test_binding_ope_max_pole(self):
        ds = ds_reduction_d4()
        assert ds['binding_ope_max_pole'] == 12

    def test_repeated_exponent_multiplicity(self):
        """The exponent 3 has multiplicity 2 (connected to triality)."""
        ds = ds_reduction_d4()
        assert ds['exponent_multiplicity_3'] == 2


# =========================================================================
# PART 8: KOSZUL DUALS
# =========================================================================

class TestKoszulDuals:
    """Verify the three Koszul dual presentations."""

    def test_koszul_dual_identification(self):
        kd = compute_koszul_duals()
        assert kd['koszul_dual'] == 'Y_hbar^{dg}(so_8)'
        assert kd['cohomology'] == 'Y_hbar(so_8)'

    def test_triality_equivalence(self):
        kd = compute_koszul_duals()
        assert kd['triality_equivalence']

    def test_root_multiplicity_one(self):
        """Root multiplicity 1 ensures complete strictification."""
        kd = compute_koszul_duals()
        assert kd['root_multiplicity'] == 1


# =========================================================================
# PART 9: COMPLETE PACKAGE
# =========================================================================

class TestCompletePackage:
    """Verify the complete D_4 triality bar complex package."""

    def test_complete_runs(self):
        """The complete computation finishes without error."""
        data = complete_d4_triality_bar(k=1)
        assert data is not None
        assert 'triality' in data
        assert 'bar_complex' in data

    def test_k_minus_h_dual(self):
        """At the critical level k = -h^v = -6, kappa should vanish
        (actually kappa = dim(g)*(k+h^v)/(2h^v) = 28*0/12 = 0)."""
        bc = triality_invariant_bar_complex(k=-6)
        assert bc['kappa'] == 0

    def test_koszul_complementarity_general(self):
        """kappa + kappa' = 0 for arbitrary level."""
        for k in [1, 2, -3, Fraction(1, 2)]:
            bc = triality_invariant_bar_complex(k=k)
            assert bc['complementarity_verified'], f"Failed at k={k}"
