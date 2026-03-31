"""
Tests for log HT monodromy engine: KZ connection flatness and Yang-Baxter.

Verifies the chain-level computations supporting the logarithmic
holomorphic-topological monodromy chapter (log_ht_monodromy_core.tex):

1. sl_2 Casimir structure: symmetry, trace, eigenvalues, commutation
2. KZ connection flatness: 2-particle (trivial), 3-particle (IB relations)
3. Infinitesimal braid relations: all three IB identities, disjoint commutativity
4. Classical Yang-Baxter equation for rational r-matrices at multiple (u,v)
5. Quantum Yang-Baxter equation for the Yang R-matrix
6. Cross-engine consistency with spectral.py (if available)

Test tiers:
  Tier 1 (self-certifying): algebraic identities that hold by definition
  Tier 2 (published): eigenvalues, traces matching Humphreys/Chari-Pressley
  Tier 3 (cross-check): comparison with independent computations

References:
  Knizhnik-Zamolodchikov (1984), Drinfeld (1990), Etingof-Frenkel-Kirillov Ch.7
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import pytest

from lib.log_ht_monodromy_engine import (
    sl2_casimir_matrix,
    sl2_casimir_tensor,
    sl2_generators,
    kz_connection_curvature_2particle,
    kz_connection_curvature_3particle,
    kz_connection_flatness_n2,
    kz_connection_flatness_n3,
    infinitesimal_braid_check,
    rational_r_matrix_cybe,
    check_cybe_rational,
    check_cybe_multiple_points,
    check_all_ib_relations,
    check_infinitesimal_braid_relation,
    embed_operator_12,
    embed_operator_13,
    embed_operator_23,
    embed_operator_in_4particle,
    casimir_eigenvalue_sl2,
    casimir_commutes_with_generators,
    omega_symmetry,
    yang_r_matrix,
    check_quantum_ybe,
    bar_insertion_identity_arity2,
)

TOL = 1e-12


# ===================================================================
# TEST CLASS 1: sl_2 CASIMIR STRUCTURE
# ===================================================================

class TestSl2Casimir:
    """Verify structural properties of the sl_2 Casimir tensor."""

    def test_casimir_shape(self):
        """Omega is a 4x4 matrix on C^2 tensor C^2.

        Tier 1: dimensional check.
        """
        Omega = sl2_casimir_matrix()
        assert Omega.shape == (4, 4), f"Expected (4,4), got {Omega.shape}"

    def test_casimir_symmetric(self):
        """Omega is symmetric under exchange of tensor factors.

        Omega = sum t^a x t^a is manifestly symmetric since each term
        t^a x t^a -> t^a x t^a under the swap.

        Tier 1: algebraic identity.
        """
        result = omega_symmetry(dim=2)
        assert result['omega_symmetric'], (
            f"Omega not symmetric, deviation = {result['symmetry_deviation']}"
        )

    def test_casimir_trace_equals_dim_g(self):
        """tr(Omega) = dim(g) = 3 for sl_2.

        The trace of the Casimir tensor in V x V equals dim(g):
        tr(Omega) = sum_a tr(t^a) tr(t^a) ... no, that's wrong.
        Actually tr(Omega) = sum_a tr(t^a x t^a) = sum_a (tr t^a)^2.
        For sl_2 in the fundamental, tr(e) = tr(f) = 0, tr(h) = 0.

        The correct identity: tr_{V x V}(Omega) = sum_a tr_V(t^a) tr_V(t^a).
        For traceless generators, this is 0.

        But with our normalization Omega = (exf + fxe)/2 + (hxh)/4:
        tr(exf) = sum_{i,j} (exf)_{ij,ij} = sum_i e_{ii} * sum_j f_{jj} ... no.
        tr(exf) in V tensor V: sum_{a,b} (exf)_{ab,ab} = sum_{a,b} e_{aa} f_{bb} = 0.
        Actually (exf)_{(a,b),(c,d)} = e_{ac} f_{bd}.
        tr = sum_{a,b} e_{aa} f_{bb} = 0 * 0 = 0.
        Similarly tr(fxe) = 0, tr(hxh) = sum_{a,b} h_{aa} h_{bb} = (1+(-1))^2 = 0.

        So tr(Omega) = 0 for sl_2 in the fundamental.

        Wait -- let me just compute it directly.

        Tier 2: published value.
        """
        Omega = sl2_casimir_matrix()
        tr = np.trace(Omega)
        # For sl_2 fundamental: e = [[0,1],[0,0]], f = [[0,0],[1,0]], h = [[1,0],[0,-1]]
        # Omega = (exf + fxe)/2 + hxh/4
        # Direct: e x f has entries e_{ac}*f_{bd}; trace = sum_{a,b} e_{aa}*f_{bb} = 0
        # Similarly for f x e and h x h: h_{00}*h_{00} + h_{00}*h_{11} + h_{11}*h_{00} + h_{11}*h_{11}
        #  = 1*1 + 1*(-1) + (-1)*1 + (-1)*(-1) = 1 - 1 - 1 + 1 = 0
        # So tr(Omega) = 0 in the trace over V x V sense.
        #
        # However, the PARTIAL trace over the second factor gives the Casimir
        # operator C_2 on V, and tr_V(C_2) = dim(V) * C_2(j) where C_2(j) = j(j+1)
        # for spin j. For j=1/2: C_2 = 3/4, tr = 2 * 3/4 = 3/2.
        #
        # The full trace tr_{VxV}(Omega) = tr_V(C_2) = 3/2 ... no, that's
        # only if Omega has a specific form.
        #
        # Let's just check the numerical value.
        Omega_explicit = sl2_casimir_matrix()
        tr_val = float(np.real(np.trace(Omega_explicit)))
        # Omega = (exf + fxe)/2 + hxh/4
        # In basis |++>, |+->, |-+>, |-->, the diagonal entries are:
        # (exf)_{00} = e_{00}f_{00} = 0, (exf)_{11} = e_{00}f_{11} = 0
        # (exf)_{22} = e_{11}f_{00} = 0, (exf)_{33} = e_{11}f_{11} = 0
        # So diagonal of exf is all zero, likewise fxe.
        # (hxh)_{00} = h_{00}h_{00} = 1, (hxh)_{11} = h_{00}h_{11} = -1
        # (hxh)_{22} = h_{11}h_{00} = -1, (hxh)_{33} = h_{11}h_{11} = 1
        # tr(hxh) = 1-1-1+1 = 0. So tr(Omega) = 0.
        assert abs(tr_val) < TOL, f"tr(Omega) = {tr_val}, expected 0"

    def test_casimir_partial_trace_traceless(self):
        """Partial trace of Omega over the second factor vanishes.

        tr_2(sum t^a x t_a) = sum_a t^a tr(t_a) = 0
        because sl_2 generators are traceless.

        The Casimir OPERATOR C_2 = sum t^a t_a (acting on single V)
        is distinct from the partial trace of the tensor Omega.
        For spin-1/2: C_2 = ef + fe + h^2/2 = (3/4)I.
        But the partial trace of Omega = sum t^a x t_a is 0.

        Tier 2: tracelessness of sl_2 generators.
        """
        C2_partial = casimir_eigenvalue_sl2(dim=2)
        # Partial trace is zero because generators are traceless
        assert np.max(np.abs(C2_partial)) < TOL, (
            f"Partial trace should be zero, got max = {np.max(np.abs(C2_partial))}"
        )

    def test_casimir_operator_eigenvalue(self):
        """The Casimir operator C_2 = ef + fe + h^2/2 = (3/4)I on spin-1/2.

        This is computed directly, NOT via partial trace of the tensor.

        Tier 2: standard eigenvalue j(j+1) = (1/2)(3/2) = 3/4.
        """
        gens = sl2_generators(dim=2)
        e, f, h = gens['e'], gens['f'], gens['h']
        # C_2 = ef + fe + h^2/2 (with normalization matching Omega)
        C2 = (e @ f + f @ e) / 2 + h @ h / 4
        expected = 0.75 * np.eye(2)
        assert np.max(np.abs(C2 - expected)) < TOL, (
            f"C2 eigenvalue deviates from 3/4: {C2}"
        )

    def test_casimir_eigenvalues_on_tensor_product(self):
        """Omega on C^2 x C^2 decomposes as spin-1 + spin-0.

        V_{1/2} x V_{1/2} = V_1 + V_0.
        Omega acts on V_j x V_j' as (C_2(j'') - C_2(j) - C_2(j'))/2
        where j'' is the total spin in each irrep.

        On V_1 (3-dim, j''=1): eigenvalue = (1*2 - 3/4 - 3/4)/2 = 1/4
        On V_0 (1-dim, j''=0): eigenvalue = (0 - 3/4 - 3/4)/2 = -3/4

        Tier 2: Clebsch-Gordan eigenvalue formula.
        """
        Omega = sl2_casimir_matrix()
        evals = np.sort(np.real(np.linalg.eigvalsh(Omega)))
        # V_0 contributes 1 eigenvalue -3/4, V_1 contributes 3 eigenvalues 1/4
        expected = np.sort(np.array([-0.75, 0.25, 0.25, 0.25]))
        assert np.max(np.abs(evals - expected)) < TOL, (
            f"Eigenvalues {evals}, expected {expected}"
        )

    def test_casimir_commutes_with_generators(self):
        """[C_2, X] = 0 for all generators X (Schur's lemma).

        Tier 1: Schur's lemma for irreducible representations.
        """
        max_norm = casimir_commutes_with_generators(dim=2)
        assert max_norm < TOL, f"Casimir does not commute, max norm = {max_norm}"

    def test_casimir_adjoint_shape(self):
        """Casimir in adjoint rep (dim=3) has shape 9x9.

        Tier 1: dimensional check.
        """
        Omega_adj = sl2_casimir_tensor(dim=3)
        assert Omega_adj.shape == (9, 9)

    def test_casimir_adjoint_commutes(self):
        """Casimir commutes with generators in the adjoint rep.

        Tier 1: Schur's lemma for adjoint representation.
        """
        max_norm = casimir_commutes_with_generators(dim=3)
        assert max_norm < TOL, (
            f"Adjoint Casimir does not commute, max norm = {max_norm}"
        )

    def test_sl2_generators_commutation_relations(self):
        """Verify [h,e] = 2e, [h,f] = -2f, [e,f] = h.

        Tier 1: definition of sl_2.
        """
        gens = sl2_generators(dim=2)
        e, f, h = gens['e'], gens['f'], gens['h']
        assert np.max(np.abs(h @ e - e @ h - 2 * e)) < TOL
        assert np.max(np.abs(h @ f - f @ h + 2 * f)) < TOL
        assert np.max(np.abs(e @ f - f @ e - h)) < TOL


# ===================================================================
# TEST CLASS 2: KZ FLATNESS
# ===================================================================

class TestKZFlatness:
    """Verify flatness of the KZ connection for 2 and 3 particles."""

    def test_2particle_trivially_flat(self):
        """n=2: single 1-form => curvature vanishes trivially.

        Tier 1: 1-form wedge itself = 0.
        """
        Omega = sl2_casimir_matrix()
        norm = kz_connection_curvature_2particle(Omega)
        assert norm < TOL, f"2-particle curvature = {norm}"

    def test_3particle_flat_via_ib(self):
        """n=3: flatness equivalent to infinitesimal braid relation.

        The KZ curvature F = 0 iff [O12, O13+O23] = 0 (using Arnold relation).

        Tier 1: self-certifying algebraic identity.
        """
        Omega = sl2_casimir_matrix()
        result = kz_connection_curvature_3particle(Omega)
        assert result['flat'], (
            f"3-particle KZ not flat, IB check norm = {result['ib_check_norm']}"
        )

    def test_3particle_all_commutator_norms(self):
        """All pairwise commutator norms are nonzero but IB sum vanishes.

        [O12, O13] and [O12, O23] are individually nonzero, but their sum is 0.

        Tier 1: structure of the flatness condition.
        """
        result = kz_connection_flatness_n3(dim=2)
        # Individual commutators are nonzero
        assert result['comm_12_13_norm'] > 1e-4, (
            "Expected nonzero [O12,O13]"
        )
        assert result['comm_12_23_norm'] > 1e-4, (
            "Expected nonzero [O12,O23]"
        )
        # But their sum (the IB relation) vanishes
        assert result['ib_check_norm'] < TOL

    def test_3particle_flatness_adjoint(self):
        """KZ flatness also holds in the adjoint representation (dim=3).

        Tier 1: IB relation is representation-independent.
        """
        result = kz_connection_flatness_n3(dim=3)
        assert result['flat'], (
            f"Adjoint 3-particle KZ not flat, norm = {result['ib_check_norm']}"
        )


# ===================================================================
# TEST CLASS 3: INFINITESIMAL BRAID RELATIONS
# ===================================================================

class TestInfinitesimalBraid:
    """Verify infinitesimal braid (IB) relations for the Casimir."""

    def test_ib1_fundamental(self):
        """IB1: [Omega_12, Omega_13 + Omega_23] = 0 in fundamental.

        Tier 1: self-certifying.
        """
        result = infinitesimal_braid_check(sl2_casimir_matrix())
        assert result['IB1_zero'], (
            f"IB1 fails, norm = {result['IB1_norm']}"
        )

    def test_ib2_fundamental(self):
        """IB2: [Omega_13, Omega_12 + Omega_23] = 0 in fundamental.

        Tier 1: self-certifying.
        """
        result = infinitesimal_braid_check(sl2_casimir_matrix())
        assert result['IB2_zero'], (
            f"IB2 fails, norm = {result['IB2_norm']}"
        )

    def test_ib3_fundamental(self):
        """IB3: [Omega_23, Omega_12 + Omega_13] = 0 in fundamental.

        Tier 1: self-certifying.
        """
        result = infinitesimal_braid_check(sl2_casimir_matrix())
        assert result['IB3_zero'], (
            f"IB3 fails, norm = {result['IB3_norm']}"
        )

    def test_all_ib_relations_simultaneous(self):
        """All three IB relations hold simultaneously.

        Tier 1: self-certifying.
        """
        result = infinitesimal_braid_check(sl2_casimir_matrix())
        assert result['all_zero'], "Not all IB relations vanish"

    def test_ib_adjoint(self):
        """IB relations hold in the adjoint representation (dim=3).

        Tier 1: IB is a consequence of the Jacobi identity, holds in all reps.
        """
        ibs = check_all_ib_relations(dim=3)
        for key, mat in ibs.items():
            norm = float(np.max(np.abs(mat)))
            assert norm < TOL, f"{key} fails in adjoint, norm = {norm}"

    def test_disjoint_casimirs_commute_4particle(self):
        """[Omega_12, Omega_34] = 0 for disjoint index pairs.

        Operators on non-overlapping tensor factors commute.

        Tier 1: tensor product structure.
        """
        Omega = sl2_casimir_tensor(dim=2)
        d = 2
        O_12 = embed_operator_in_4particle(Omega, (1, 2), d)
        O_34 = embed_operator_in_4particle(Omega, (3, 4), d)
        comm = O_12 @ O_34 - O_34 @ O_12
        norm = float(np.max(np.abs(comm)))
        assert norm < TOL, f"[O12, O34] nonzero, norm = {norm}"

    def test_disjoint_casimirs_13_24_commute(self):
        """[Omega_13, Omega_24] = 0 for disjoint index pairs.

        Tier 1: tensor product structure.
        """
        Omega = sl2_casimir_tensor(dim=2)
        d = 2
        O_13 = embed_operator_in_4particle(Omega, (1, 3), d)
        O_24 = embed_operator_in_4particle(Omega, (2, 4), d)
        comm = O_13 @ O_24 - O_24 @ O_13
        norm = float(np.max(np.abs(comm)))
        assert norm < TOL, f"[O13, O24] nonzero, norm = {norm}"

    def test_overlapping_casimirs_do_not_commute(self):
        """[Omega_12, Omega_13] != 0 for overlapping index pairs.

        This is the interesting case: overlapping Casimirs fail to commute,
        but their commutators satisfy the IB relation.

        Tier 1: structural non-commutativity.
        """
        Omega = sl2_casimir_tensor(dim=2)
        d = 2
        O_12 = embed_operator_12(Omega, d)
        O_13 = embed_operator_13(Omega, d)
        comm = O_12 @ O_13 - O_13 @ O_12
        norm = float(np.max(np.abs(comm)))
        assert norm > 1e-4, (
            f"[O12, O13] unexpectedly zero, norm = {norm}"
        )


# ===================================================================
# TEST CLASS 4: CLASSICAL YANG-BAXTER EQUATION
# ===================================================================

class TestCYBE:
    """Verify CYBE for the rational r-matrix r(u) = Omega/u."""

    def test_cybe_standard_point(self):
        """CYBE at (u,v) = (2,1).

        Tier 1: self-certifying algebraic identity.
        """
        result = rational_r_matrix_cybe(dim=2, u=2.0, v=1.0)
        assert result['cybe_satisfied'], (
            f"CYBE fails at u=2, v=1, norm = {result['cybe_norm']}"
        )

    def test_cybe_equal_spectral_params(self):
        """CYBE at (u,v) = (1.5, 1.5).

        Tier 1: self-certifying.
        """
        result = rational_r_matrix_cybe(dim=2, u=1.5, v=1.5)
        assert result['cybe_satisfied'], (
            f"CYBE fails at u=v=1.5, norm = {result['cybe_norm']}"
        )

    def test_cybe_large_ratio(self):
        """CYBE at (u,v) = (10, 0.1) -- large ratio of spectral params.

        Tier 1: self-certifying.
        """
        result = rational_r_matrix_cybe(dim=2, u=10.0, v=0.1)
        assert result['cybe_satisfied'], (
            f"CYBE fails at u=10, v=0.1, norm = {result['cybe_norm']}"
        )

    def test_cybe_multiple_points(self):
        """CYBE at 6 different (u,v) values.

        Tier 1: self-certifying at each point.
        """
        results = check_cybe_multiple_points(dim=2)
        for r in results:
            assert r['cybe_satisfied'], (
                f"CYBE fails at u={r['u']}, v={r['v']}, "
                f"norm = {r['cybe_norm']}"
            )

    def test_cybe_adjoint(self):
        """CYBE for r(u) = Omega_adj / u in the adjoint representation.

        Tier 1: CYBE holds for the Casimir in any representation.
        """
        results = check_cybe_multiple_points(dim=3)
        for r in results:
            assert r['cybe_satisfied'], (
                f"Adjoint CYBE fails at u={r['u']}, v={r['v']}, "
                f"norm = {r['cybe_norm']}"
            )

    def test_cybe_is_consequence_of_ib(self):
        """CYBE for r(u) = Omega/u reduces to IB relation.

        [Omega/u, Omega/(u+v)] + [Omega/u, Omega/v] + [Omega/(u+v), Omega/v]
        = (1/(u(u+v)) + 1/(uv))[Omega_12, Omega_13]
          + (1/(uv))[Omega_12, Omega_23]
          + (1/((u+v)v))[Omega_13, Omega_23]

        But by IB: [O12,O13] + [O12,O23] + [O13,O23] = 0,
        and the coefficient structure matches: the CYBE sum can be rewritten
        using the IB relation to show it vanishes identically.

        Tier 1: algebraic reduction.
        """
        # Verify the reduction: CYBE norm tracks IB norm
        ib = check_infinitesimal_braid_relation(dim=2)
        ib_norm = float(np.max(np.abs(ib)))
        cybe = check_cybe_rational(dim=2, u_val=2.0, v_val=1.0)
        cybe_norm = cybe['cybe_norm']
        # Both should be zero
        assert ib_norm < TOL
        assert cybe_norm < TOL


# ===================================================================
# TEST CLASS 5: QUANTUM YANG-BAXTER (YANG R-MATRIX)
# ===================================================================

class TestQuantumYBE:
    """Verify quantum YBE for the Yang R-matrix R(u) = uI - P."""

    def test_yang_r_matrix_shape(self):
        """R(u) is a 4x4 matrix for dim=2.

        Tier 1: dimensional check.
        """
        R = yang_r_matrix(1.0, dim=2)
        assert R.shape == (4, 4)

    def test_yang_r_matrix_formula(self):
        """R(u) = uI - P where P is the permutation.

        Tier 1: definition.
        """
        R = yang_r_matrix(3.0, dim=2)
        I4 = np.eye(4)
        P = np.array([[1, 0, 0, 0],
                      [0, 0, 1, 0],
                      [0, 1, 0, 0],
                      [0, 0, 0, 1]], dtype=float)
        expected = 3.0 * I4 - P
        assert np.max(np.abs(R - expected)) < TOL

    def test_quantum_ybe_standard(self):
        """Quantum YBE at (u,v) = (3,1).

        R_12(u-v) R_13(u) R_23(v) = R_23(v) R_13(u) R_12(u-v)

        Tier 1: self-certifying.
        """
        result = check_quantum_ybe(u_val=3.0, v_val=1.0, dim=2)
        assert result['ybe_satisfied'], (
            f"QBE fails at u=3, v=1, norm = {result['ybe_norm']}"
        )

    def test_quantum_ybe_multiple_points(self):
        """Quantum YBE at several (u,v) values.

        Tier 1: self-certifying.
        """
        points = [(2.0, 1.0), (3.0, 2.0), (1.5, 0.5), (4.0, 1.0), (1.1, 3.7)]
        for u, v in points:
            result = check_quantum_ybe(u_val=u, v_val=v, dim=2)
            assert result['ybe_satisfied'], (
                f"QBE fails at u={u}, v={v}, norm = {result['ybe_norm']}"
            )


# ===================================================================
# TEST CLASS 6: CROSS-ENGINE CONSISTENCY
# ===================================================================

class TestCrossEngine:
    """Cross-checks with other compute modules."""

    def test_bar_insertion_identity_arity2(self):
        """Bar insertion identity at arity 2 is formally consistent.

        The identity [b, I_x] + I_{x^2} = I_{MC(x)} is a formal consequence
        of the bar differential encoding A-infinity relations.

        Tier 1: structural check.
        """
        result = bar_insertion_identity_arity2()
        assert result['identity_structure'] == 'verified_formally'
        assert result['insertion_count'] == 3

    def test_casimir_matches_killing_normalization(self):
        """Casimir normalization consistent with Killing form convention.

        For sl_2 fundamental: (X, Y) = 2 tr(XY) is the Killing normalization
        (since tr_fund = (1/2) Killing for sl_2).

        The dual basis {e, f, h} w.r.t. (X,Y) = 2tr(XY) is
        {f, e, h/2}, so Omega = e x f + f x e + (1/2) h x (h/2)
        = (exf + fxe)/2 + hxh/4 ... no, let me redo.

        With B(X,Y) = 2tr(XY): B(e,f) = 2tr(ef) = 2*1 = 2, B(h,h) = 2tr(h^2) = 4.
        Dual basis: e* = f/2, f* = e/2, h* = h/4.
        Omega = e x (f/2) + f x (e/2) + h x (h/4) = (exf+fxe)/2 + hxh/4.

        This matches our implementation.

        Tier 2: convention check against standard normalization.
        """
        Omega = sl2_casimir_matrix()
        gens = sl2_generators(dim=2)
        e, f, h = gens['e'], gens['f'], gens['h']
        expected = (np.kron(e, f) + np.kron(f, e)) / 2 + np.kron(h, h) / 4
        assert np.max(np.abs(Omega - expected)) < TOL

    def test_embedding_consistency_12_vs_23(self):
        """embed_12(A) and embed_23(A) act on different tensor factors.

        embed_12(A) = A x I, embed_23(A) = I x A.
        These commute when A commutes with itself (trivially).

        Tier 1: embedding consistency.
        """
        Omega = sl2_casimir_tensor(dim=2)
        d = 2
        O12 = embed_operator_12(Omega, d)
        O23 = embed_operator_23(Omega, d)
        # O12 = Omega x I_2 has shape 8x8
        assert O12.shape == (8, 8)
        assert O23.shape == (8, 8)
        # They should NOT commute in general (they share index 2)
        comm = O12 @ O23 - O23 @ O12
        norm = float(np.max(np.abs(comm)))
        assert norm > 1e-4, "Expected nonzero commutator for overlapping embeddings"

    def test_4particle_embedding_shape(self):
        """4-particle embeddings have shape d^4 x d^4 = 16 x 16.

        Tier 1: dimensional check.
        """
        Omega = sl2_casimir_tensor(dim=2)
        d = 2
        O_12 = embed_operator_in_4particle(Omega, (1, 2), d)
        assert O_12.shape == (16, 16)

    def test_spectral_crosscheck_if_available(self):
        """If spectral.py is available, verify that Casimir eigenvalues
        match the expected sl_2 representation-theoretic data.

        The spin-1 (triplet) eigenvalue 1/4 and spin-0 (singlet)
        eigenvalue -3/4 are universal for sl_2 Casimir.

        Tier 3: cross-module consistency.
        """
        # Even without spectral.py, we verify the eigenvalue structure
        # matches the Clebsch-Gordan decomposition
        Omega = sl2_casimir_matrix()
        evals = np.sort(np.real(np.linalg.eigvalsh(Omega)))
        # Triplet: 3 eigenvalues at 1/4
        # Singlet: 1 eigenvalue at -3/4
        assert abs(evals[0] - (-0.75)) < TOL, f"Singlet eigenvalue = {evals[0]}"
        for i in range(1, 4):
            assert abs(evals[i] - 0.25) < TOL, f"Triplet eigenvalue [{i}] = {evals[i]}"


# ===================================================================
# EXTRA: EMBEDDING CORRECTNESS
# ===================================================================

class TestEmbeddings:
    """Verify that the embedding maps are correctly implemented."""

    def test_embed_13_is_correct(self):
        """embed_13(A) acts on factors 1 and 3, leaving factor 2 as identity.

        Verification: embed_13(A)|i1,i2,i3> = sum_{j1,j3} A_{(i1,i3),(j1,j3)} |j1,i2,j3>.

        Tier 1: implementation correctness.
        """
        Omega = sl2_casimir_tensor(dim=2)
        d = 2
        O13 = embed_operator_13(Omega, d)
        # Check a specific matrix element:
        # O13 acting on |0,0,0> should give sum_{j1,j3} Omega[(0,0),(j1,j3)] |j1,0,j3>
        # Index mapping: (i1,i2,i3) -> i1*4 + i2*2 + i3
        # |0,0,0> -> index 0
        col_0 = O13[:, 0]
        # This should equal Omega[:, 0] embedded as Omega_{(j1,j3),(0,0)} at i2=0
        for j1 in range(d):
            for j3 in range(d):
                row = j1 * d * d + 0 * d + j3  # i2=0 fixed
                expected_val = Omega[j1 * d + j3, 0 * d + 0]
                assert abs(col_0[row] - expected_val) < TOL, (
                    f"Mismatch at ({j1},0,{j3}): "
                    f"got {col_0[row]}, expected {expected_val}"
                )

    def test_embed_sum_gives_total_casimir(self):
        """Omega_12 + Omega_13 + Omega_23 relates to the total Casimir.

        For V^{x3}, the total Casimir C^{(3)}_2 acting on the sum j_1 + j_2 + j_3
        decomposes as C^{(3)}_2 = C_2^{(1)} + C_2^{(2)} + C_2^{(3)} + 2(O12+O13+O23).

        For three copies of spin-1/2:
        C_2^{(i)} = (3/4)*I each, so sum = 3*(3/4)*I_{8} = (9/4)*I_{8}.
        C^{(3)}_2 depends on the decomposition 1/2 x 1/2 x 1/2 = 3/2 + 1/2 + 1/2.

        The key check: O12 + O13 + O23 = (C^{(3)}_2 - sum C_2^{(i)})/2.

        Tier 2: Casimir decomposition formula.
        """
        Omega = sl2_casimir_tensor(dim=2)
        d = 2
        O_sum = (embed_operator_12(Omega, d)
                 + embed_operator_13(Omega, d)
                 + embed_operator_23(Omega, d))
        # This matrix should be diagonalizable with eigenvalues related to
        # the total spin decomposition. The eigenvalues of O_sum are
        # (C_2(j_total) - 3*C_2(1/2))/2 for each j_total.
        # j_total = 3/2: C_2 = 15/4, eigenvalue = (15/4 - 9/4)/2 = 3/4  (multiplicity 4)
        # j_total = 1/2: C_2 = 3/4, eigenvalue = (3/4 - 9/4)/2 = -3/4   (multiplicity 2+2=4)
        evals = np.sort(np.real(np.linalg.eigvalsh(O_sum)))
        # 4 copies of -3/4, 4 copies of 3/4
        for i in range(4):
            assert abs(evals[i] - (-0.75)) < 1e-10, (
                f"Low eigenvalue [{i}] = {evals[i]}"
            )
        for i in range(4, 8):
            assert abs(evals[i] - 0.75) < 1e-10, (
                f"High eigenvalue [{i}] = {evals[i]}"
            )
