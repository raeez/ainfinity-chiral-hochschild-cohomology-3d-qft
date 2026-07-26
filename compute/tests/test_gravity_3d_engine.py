"""Tests for 3D Gravity Compute Engine.

Verifies the gravitational A∞ structure derived from Virasoro line comparison:
1. Virasoro associator A₃ and ternary operation m₃ = -A₃
2. Quartic contact invariant Q^contact = 10/(c(5c+22))
3. Gravitational Koszul triangle: same-family line representative Vir_{26-c}
4. Modular characteristic κ(Vir_c) = c/2, complementarity sum = 13
5. Laplace kernel and bar collision r-matrix pole structures
6. Genus expansion via Â-genus formula
7. Cross-engine consistency with bulk_boundary_duality_engine

References:
  Vol II: 3d_gravity.tex (Movements I-VI)
  Vol I: concordance.tex, nonlinear_modular_shadows.tex
"""
import sys
import os
import math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from fractions import Fraction
from pathlib import Path
from sympy import Symbol, Rational, S, simplify, expand, symbols, diff, sqrt, series, binomial, log, pi, I, sinh, exp

from lib.gravity_3d_engine import (
    virasoro_lambda_bracket,
    virasoro_associator,
    virasoro_m3_coefficients,
    hpl_planar_binary_tree_count,
    virasoro_hpl_tree_profile,
    virasoro_hpl_sdr_data,
    virasoro_ds_linear_sdr_identities,
    virasoro_s4_completed_ambient_requirement,
    virasoro_exact_gravity_scope_profile,
    brown_henneaux_chiral_chart_scope,
    gravitational_mc_bridge_scope,
    brown_henneaux_line_test_package_scope,
    virasoro_bar_intrinsic_mc_shadow_scope,
    genus0_directed_product_decomposition_scope,
    quantum_group_clutching_equivariance_scope,
    modular_operad_unitality_scope,
    modular_bar_reduction_scope,
    affine_modular_bar_datum_scope,
    affine_e3_topological_km_scope,
    principal_ds_e3_topological_scope,
    good_graded_ds_e3_topological_scope,
    virasoro_scalar_bar_trace_profile,
    genus1_virasoro_mc_scope,
    genus1_virasoro_amplitudes_scope,
    genus1_modular_anomaly_scope,
    genus1_virasoro_kzb_shadow_kernel_scope,
    gravity_weinberg_ward_residue_scope,
    gravity_cachazo_strominger_ward_package_scope,
    gravity_chy_quartic_contact_scope,
    virasoro_shadow_metric_coefficients,
    gravity_infinite_soft_shadow_hierarchy_scope,
    virasoro_shadow_metric_pole_profile,
    stasheff_scalar_pole_cancellation_scope,
    virasoro_catalan_shape_factor,
    virasoro_shadow_closed_form_coefficient,
    shadow_closed_form_scope,
    catalan_dynkin_field_polynomial,
    catalan_dynkin_parity_scope,
    crossing_stasheff_scope,
    shapovalov_channel_norm_squared,
    shapovalov_projected_channel_norm_squared,
    shapovalov_bootstrap_scope,
    virasoro_large_c_shadow_asymptotics,
    large_c_bootstrap_scope,
    otoc_braiding_phase,
    otoc_r_matrix_scope,
    mss_bound_value,
    mss_annular_bar_scope,
    scrambling_time_from_amplitude,
    scrambling_time_scope,
    w3_w_line_shadow_coefficient,
    w3_w_line_closed_form_scope,
    maloney_witten_scope_profile,
    page_curve_profile,
    page_curve_scope,
    desitter_central_charge,
    desitter_horizon_entropy_from_radius,
    desitter_horizon_entropy_from_c,
    desitter_shadow_profile,
    jt_wp_spectral_curve_y,
    jt_disk_density,
    jt_wp_to_density_balance,
    jt_schwarzian_scope,
    colored_partition_number,
    k3_borcherds_scalar_bridge_scope,
    k3_scalar_no_promotion_scope,
    heptagon_growth_bound_constant,
    heptagon_growth_bound_scope,
    factorization_l0_gevrey_constant,
    factorization_l0_gevrey_scope,
    factorization_bcov_candidate_location,
    factorization_bcov_candidate_residue,
    factorization_bcov_instanton_scope,
    heisenberg_zwegers_shadow_scope,
    partition_number,
    virasoro_vacuum_verma_coefficient,
    virasoro_vacuum_verma_asymptotic,
    virasoro_hardy_ramanujan_cardy_scope,
    quartic_contact_virasoro,
    quartic_contact_virasoro_exact,
    gravity_kappa,
    gravity_kappa_exact,
    koszul_dual_central_charge,
    gravity_laplace_kernel_poles,
    gravity_r_matrix_poles,
    gravity_r_matrix_leading_residue,
    virasoro_collision_cybe_scope,
    virasoro_ds_hpl_transfer_scope,
    ds_ordered_bar_intertwining_scope,
    principal_ds_coproduct_primitivity_scope,
    virasoro_primary_ward_connection,
    virasoro_primary_ward_log,
    virasoro_primary_ward_even_coefficients,
    virasoro_primary_ward_even_exponents,
    complementarity_constant_virasoro,
    verify_complementarity,
    verdier_modular_s_shadow_scope,
    genus_generating_function_coefficients,
    ahat_series_coefficients,
    verify_ahat_series,
    virasoro_shadow_depth,
    virasoro_shadow_class,
)


# ===================================================================
# 1. VIRASORO LAMBDA-BRACKET
# ===================================================================

class TestVirasoroLambdaBracket:
    """Verify {T_λ T} = (c/12)λ³ + 2Tλ + ∂T."""

    def test_quartic_pole_coefficient(self):
        """Central charge coefficient: c/12 for the λ³ term."""
        bracket = virasoro_lambda_bracket(c=24)
        assert bracket['lam3'] == 2  # 24/12 = 2

    def test_double_pole_coefficient(self):
        """Coefficient of Tλ is always 2."""
        for c_val in [0, 1, 13, 26, -22]:
            bracket = virasoro_lambda_bracket(c=c_val)
            assert bracket['lam1_T'] == 2

    def test_simple_pole_coefficient(self):
        """Coefficient of ∂T is always 1."""
        bracket = virasoro_lambda_bracket(c=13)
        assert bracket['lam0_dT'] == 1

    def test_symbolic_central_charge(self):
        """Symbolic c works."""
        bracket = virasoro_lambda_bracket()
        c = Symbol('c')
        assert simplify(bracket['lam3'] - c / 12) == 0

    def test_central_charge_stored(self):
        """Central charge value is stored."""
        bracket = virasoro_lambda_bracket(c=26)
        assert bracket['central_charge'] == 26


# ===================================================================
# 2. VIRASORO ASSOCIATOR
# ===================================================================

class TestVirasoroAssociator:
    """Verify A₃(T,T,T; λ₁₂, λ₂₃) formula from eq:gravity-associator."""

    def test_d2T_coefficient(self):
        """Coefficient of ∂²T is -1 (always)."""
        A3 = virasoro_associator(c=1, lam12=0, lam23=0)
        assert A3['d2T'] == -1

    def test_scalar_linear_in_c(self):
        """Scalar term is proportional to c."""
        A3_1 = virasoro_associator(c=1, lam12=1, lam23=1)
        A3_2 = virasoro_associator(c=2, lam12=1, lam23=1)
        # scalar(c=2) should be 2 * scalar(c=1)
        assert simplify(A3_2['scalar'] - 2 * A3_1['scalar']) == 0

    def test_scalar_vanishes_at_c0(self):
        """At c=0 the scalar (central extension) term vanishes."""
        A3 = virasoro_associator(c=0, lam12=3, lam23=5)
        assert A3['scalar'] == 0

    def test_specific_evaluation(self):
        """A₃ at λ₁₂=0, λ₂₃=1, c=12: scalar = -(12/12)*1³*(0+1) = -1."""
        A3 = virasoro_associator(c=12, lam12=0, lam23=1)
        assert A3['scalar'] == -1

    def test_specific_evaluation_2(self):
        """A₃ at λ₁₂=1, λ₂₃=1, c=12: scalar = -1*1*(2+1) = -3."""
        A3 = virasoro_associator(c=12, lam12=1, lam23=1)
        assert A3['scalar'] == -3

    def test_dT_coefficient_at_origin(self):
        """dT coefficient at λ₁₂=0, λ₂₃=0 is 0."""
        A3 = virasoro_associator(c=1, lam12=0, lam23=0)
        assert A3['dT'] == 0

    def test_T_coefficient_at_origin(self):
        """T coefficient at λ₁₂=0, λ₂₃=0 is 0."""
        A3 = virasoro_associator(c=1, lam12=0, lam23=0)
        assert A3['T'] == 0

    def test_not_symmetric(self):
        """A₃ is NOT symmetric in λ₁₂ ↔ λ₂₃."""
        A3_a = virasoro_associator(c=1, lam12=1, lam23=2)
        A3_b = virasoro_associator(c=1, lam12=2, lam23=1)
        assert A3_a['scalar'] != A3_b['scalar']


# ===================================================================
# 3. m₃ = -A₃ (BORCHERDS IDENTITY)
# ===================================================================

class TestM3NegatesAssociator:
    """Verify m₃ = -A₃ from the Borcherds identity."""

    def test_m3_negates_d2T(self):
        m3 = virasoro_m3_coefficients(c=1, lam12=0, lam23=0)
        assert m3['d2T'] == 1  # -(-1) = 1

    def test_m3_negates_scalar(self):
        A3 = virasoro_associator(c=12, lam12=1, lam23=1)
        m3 = virasoro_m3_coefficients(c=12, lam12=1, lam23=1)
        assert simplify(m3['scalar'] + A3['scalar']) == 0

    def test_m3_negates_all_components(self):
        """Full check: m₃ + A₃ = 0 for all components."""
        for c_val in [1, 6, 13, 26]:
            A3 = virasoro_associator(c=c_val, lam12=2, lam23=3)
            m3 = virasoro_m3_coefficients(c=c_val, lam12=2, lam23=3)
            for key in ['d2T', 'dT', 'T', 'scalar']:
                assert simplify(m3[key] + A3[key]) == 0, \
                    f"m3[{key}] + A3[{key}] ≠ 0 at c={c_val}"

    def test_symbolic_negation(self):
        """Symbolic verification of m₃ = -A₃."""
        A3 = virasoro_associator()
        m3 = virasoro_m3_coefficients()
        assert simplify(m3['scalar'] + A3['scalar']) == 0


# ===================================================================
# 4. HPL MINIMAL TRANSFER DATA
# ===================================================================

class TestVirasoroHPLTransfer:
    """Verify the explicit HPL transfer data used for Virasoro class M."""

    def test_planar_binary_tree_counts(self):
        """|PRT_n| is the Catalan number C_{n-1}."""
        assert [hpl_planar_binary_tree_count(n) for n in range(2, 7)] == [
            1, 2, 5, 14, 42
        ]

    def test_tree_profile_records_transfer_labels(self):
        """An n-ary transfer summand has p at the root, i at leaves, and h internally."""
        profile = virasoro_hpl_tree_profile(5)
        assert profile['tree_set'] == 'PRT_5'
        assert profile['summands'] == 14
        assert profile['leaves'] == 5
        assert profile['binary_vertices_mu'] == 4
        assert profile['internal_h_edges'] == 3
        assert profile['root_label'] == 'p'
        assert profile['leaf_label'] == 'i'
        assert profile['sign_convention'] == 'suspended-bar Koszul sign'
        assert profile['internal_edge_order'] == 'left-to-right depth-first planar order'
        assert profile['orientation_line'] == 'det(k^{E_int(T)})'
        assert profile['free_pm_choice'] is False
        assert profile['formula_label'] == 'eq:gravity-suspended-hpl-transfer'

    def test_sdr_identity_and_antighost_weight_normalization(self):
        """The antighost homotopy is h|A_w = G_0/w on positive weights."""
        data = virasoro_hpl_sdr_data(weight=4)
        assert data['homotopy_identity'] == 'Qh + hQ = id - i∘p'
        assert 'pi = 1_H' in data['side_conditions']
        assert data['finite_hpl_datum'] == (
            'Q_DS',
            'mu=m_2',
            'delta_DS',
            'i',
            'p',
            'h_DS',
            'complete filtration',
        )
        assert data['positive_weight_homotopy'] == 'h|A_w = G_0/w'
        assert data['requires_antighost'] == 'G_0 with [Q,G_0] = L_0'
        assert data['homotopy_coefficient'] == Rational(1, 4)

    def test_ds_linear_sdr_identities_generatorwise(self):
        """The displayed DS maps contract only the linear d0 complex."""
        data = virasoro_ds_linear_sdr_identities(k=0)
        assert data['complex'] == 'C_DS_lin over C[partial]'
        assert data['level_condition'] == 'k != -2'
        assert data['h2_zero']
        assert data['ph_zero']
        assert data['hi_zero']
        assert data['homotopy_identity']
        assert data['full_brst_retract'] is False
        assert data['perturbation_needed_for_full_brst'] == 'delta = d_BRST - d0'

    def test_ds_linear_sdr_rejects_critical_level(self):
        """The formulas use W=(k+2)T and are undefined at k=-2."""
        with pytest.raises(ValueError):
            virasoro_ds_linear_sdr_identities(k=-2)

    def test_ds_hpl_transfer_scope_is_homotopy_coherent_not_strict_yangian(self):
        """The DS-HPL theorem constructs a spectral package, not a strict Yangian."""
        scope = virasoro_ds_hpl_transfer_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['constructed_object'] == 'homotopy-coherent Virasoro spectral package'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'delta')
        assert 'hypAmbientWtCpl' in scope['requires']
        assert 'hypKZSDR' in scope['requires']
        assert 'HPL-transferred A_infinity morphism Delta_z' in scope['proved_core']
        assert 'homotopy-coherent A_infinity Yang-Baxter relation' in scope['proved_core']
        assert scope['collision_kernel'] == 'r_coll_Vir(z) = (c/2)/z^3 + 2T/z'
        assert 'strict dg-shifted Yangian presentation of H_Vir' in scope['not_proved_here']
        assert scope['is_strict_or_dg_shifted_yangian_assertion'] is False
        assert scope['promotes_after_open_axioms'] is True

    def test_ds_sdr_proposition_is_linear_and_licensed(self):
        """The manuscript must not state the linear SDR as a full BRST retract."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Linear DS deformation retract;"
        )
        end = source.index(r"\end{proposition}", start)
        prop = source[start:end]
        local = source[source.rfind("Assume", 0, start):end]

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\alpha+\gamma$",
            r"\hypAmbientWtCpl",
            r"k\ne -2",
            r"C^{\mathrm{DS}}_{\mathrm{lin}}",
            r"\C[\partial]",
            r"\operatorname{id} - ip = d_0 h + h d_0",
            r"d_{\mathrm{BRST}}",
            "perturbed object",
        )
        for needle in required:
            assert needle in prop or needle in local

        retired = (
            r"\begin{proposition}[SDR for the DS complex]",
            r"retracts $C^{\mathrm{DS}}$ onto the",
            "Direct verification on each generator.",
        )
        for needle in retired:
            assert needle not in local

    def test_ds_hpl_transfer_theorem_does_not_overclaim_yangian(self):
        """The theorem statement must separate HPL transfer from Yangian promotion."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[DS-HPL Virasoro spectral transfer;"
        )
        end = source.index(r"\end{theorem}", start)
        theorem = source[start:end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma+\delta$",
            r"k\ne -2",
            r"\hypAmbientWtCpl",
            r"\hypKZSDR",
            r"\mathsf Y_{\mathrm{Vir}}^{\mathrm{HPL}}",
            r"\widetilde i=s_C i s_H^{-1}",
            r"E_{\mathrm{int}}(T)=\{e_1<\cdots<e_N\}",
            r"\det(\Bbbk^{E_{\mathrm{int}}(T)})",
            r"\label{eq:gravity-suspended-hpl-transfer}",
            r"\sum_{T\in\operatorname{PRT}_n}\widetilde m_T",
            r"there is no free \(\pm\)-choice",
            r"r_{\mathrm{Vir}}^{\mathrm{coll}}(z)",
            r"\frac{c/2}{z^3}+\frac{2T}{z}",
            r"homotopy-coherent \(\Ainf\) Yang--Baxter relation",
            "does not assert a strict or dg-shifted Yangian",
            "rationality, translation-compatibility, and sign-normalised",
        )
        for needle in required:
            assert needle in theorem

        retired = (
            r"\begin{theorem}[DS-HPL transfer]",
            "is a dg-shifted Yangian on the Virasoro algebra",
            r"r^{\mathrm{Vir}}(z) = (c/2)/z^3 + 2T/z",
            "The transferred data",
            r"\pm\,p",
            "higher-order HPL corrections",
        )
        for needle in retired:
            assert needle not in theorem

    def test_ds_hpl_honest_gaps_match_the_theorem_scope(self):
        """The honest-gaps remark must not contradict the theorem statement."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{remark}[Honest gaps]")
        end = source.index(r"\end{remark}", start)
        remark = source[start:end]
        flat = " ".join(remark.split())

        assert "Three verifications remain beyond the HPL transfer core" in flat
        assert "translation compatibility" in remark
        assert "rational strictification" in remark
        assert "dg-shifted Yangian rather than only the collision-residue" in flat
        assert "The weight-filtration convergence and the morphism-transfer compatibility are the HPL output" in flat
        assert "Two verifications remain beyond the proved core" not in remark

    def test_ds_ordered_bar_scope_is_filtered_not_unconditional(self):
        """The ordered DS comparison is filtered and conditional."""
        scope = ds_ordered_bar_intertwining_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'delta')
        assert 'principal DS chart at non-critical level' in scope['requires']
        assert 'hypAmbientWtCpl' in scope['requires']
        assert 'hypKZSDR' in scope['requires']
        assert 'finite-weight BRST concentration' in scope['requires']
        assert 'ordered BRST-bar bicomplex' in scope['requires']
        assert 'HPL-transferred R-descent datum' in scope['requires']
        assert scope['object'] == 'filtered ordered-bar comparison in the homotopy category'
        assert scope['comparison'] == (
            'B^ord(W_k(g)) -> H^0_DS(B^ord(V_k(g))) filtered quasi-isomorphism'
        )
        assert scope['base_space'] == 'Conf^<_bullet(R) x FM_bullet(C)'
        assert scope['descent_square'] == 'commutes up to transferred R-descent homotopy'
        assert 'unconditional ordered DS functor' in scope['not_asserted']
        assert 'strict equality of ordered bar dg-coalgebras' in scope['not_asserted']
        assert 'R-descent unchanged by DS reduction' in scope['not_asserted']
        assert 'strict dg-shifted Yangian presentation' in scope['not_asserted']
        assert scope['h0_requires_concentration'] is True
        assert scope['r_matrix_transferred_not_copied'] is True
        assert scope['yangian_promotion_separate'] is True

    def test_ds_ordered_bar_theorem_and_remarks_are_scoped(self):
        """The manuscript must not state unconditional ordered DS functoriality."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Filtered DS comparison for the ordered bar"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        yangian_start = source.index(r"\begin{remark}[Consequences for the gravitational Yangian]", proof_end)
        yangian_end = source.index(r"\end{remark}", yangian_start)
        yangian = source[yangian_start:yangian_end]
        table_start = source.index(r"\begin{remark}[The gauge-gravity-matter trichotomy:", yangian_end)
        table_end = source.index(r"\end{remark}", table_start)
        table = source[table_start:table_end]
        combined = " ".join((block + yangian + table).split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma+\delta$",
            r"\hypAmbientWtCpl",
            r"\hypKZSDR",
            r"finite-weight BRST concentration",
            r"HPL-transferred \(R\)-descent",
            r"filtered quasi-isomorphism",
            r"\simeq_{\mathrm{filt}}",
            r"homotopy category of ordered bar dg-coalgebras",
            r"does not construct an unconditional ordered DS functor",
            r"nor a strict dg-shifted Yangian presentation",
            r"right-hand side would be \(H^\bullet_{\mathrm{DS}}\)",
            r"not copied unchanged from \(V_k(\fg)\)",
            r"transferred datum is the one used",
            r"separate Yangian-promotion hypotheses",
            r"\simeq_{\mathrm{filt}}",
            r"compatible with DS reduction on this filtered principal surface",
            r"ordered-bar shadow under the hypotheses",
            r"not a primitive classification before the \(R\)-descent datum is fixed",
            r"primitive only in the scoped generator sense",
        )
        for needle in required:
            assert " ".join(needle.split()) in combined

        retired = (
            r"\begin{theorem}[DS--ordered bar intertwining",
            r"commutes with the ordered bar construction",
            r"as dg-coalgebras over",
            r"commutes up to canonical homotopy",
            r"The proof extends the symmetric bar argument",
            r"Three inputs are required",
            r"without $\Sigma_n$-symmetrisation",
            r"$[d_{\mathrm{BRST}}, d_{\mathrm{bar}}^{\mathrm{ord}}] = 0$",
            r"The PBW filtration ensures a comparison",
            r"acts on the algebra without changing the descent map",
            r"the $R$-matrix is a property of the OPE, not of the ghost sector",
            r"implies that the dg-shifted Yangian of $\cW_k(\fg)$ is obtained",
            r"two-colour architecture is fully compatible with DS reduction",
            r"extends to the ordered bar complex via",
            r"coincides with the symmetric bar class because DS reduction preserves both",
            r"The gauge sector is the unique case where the coproduct is non-primitive",
            r"\begin{remark}[Consequences for the gravitational Yangian] \begin{remark}",
            r"\begin{center} \begin{center}",
        )
        for needle in retired:
            assert " ".join(needle.split()) not in combined

    def test_principal_ds_coproduct_primitivity_scope_is_generator_level(self):
        """The coproduct theorem is generator-level and conditional on the tree count."""
        scope = principal_ds_coproduct_primitivity_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['base_case'] == 'Delta_{z,2}^Vir(T,T) = 0'
        assert scope['all_degree_statement'] == 'generator-level principal DS primitivity'
        assert scope['licensing_tags'] == ('alpha', 'gamma', 'delta')
        assert 'signed ghost-defect lemma for HPL morphism trees' in scope['requires']
        assert 'degree-two Virasoro source-tree ghost defect' in scope['proved_here']
        assert 'ordinary Hopf primitivity for arbitrary composites' in scope['not_asserted']
        assert scope['degree_two_source_tree_ghost_shift'] == -1
        assert scope['degree_two_target_tree_ghost_shift'] == -1
        assert scope['projection_accepts_only_bidegree'] == (0, 0)

    def test_coproduct_primitivity_theorems_are_scoped_and_licensed(self):
        """The manuscript must not state unconditional all-field primitivity."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()

        prop_start = source.index(
            r"\begin{proposition}[Degree-two Virasoro generator coproduct;"
        )
        prop_end = source.index(r"\end{proposition}", prop_start)
        prop = source[prop_start:prop_end]

        theorem_start = source.index(
            r"\begin{theorem}[Principal-DS generator coproduct primitivity;"
        )
        theorem_end = source.index(r"\end{theorem}", theorem_start)
        theorem = source[theorem_start:theorem_end]

        assert r"\ClaimStatusConditional" in prop
        assert r"licensing $\alpha+\gamma+\delta$" in prop
        assert r"\hypAmbientWtCpl" in prop
        assert r"\hypKZSDR" in prop
        assert "This is a generator statement" in prop

        required = (
            r"\ClaimStatusConditional",
            r"ghost-zero BRST projection",
            r"\hypAmbientWtCpl",
            r"\hypKZSDR",
            "signed ghost-defect",
            r"\Delta_{z,n}^W(W_s)",
            "does not assert ordinary Hopf-primitivity",
        )
        for needle in required:
            assert needle in theorem

        retired = (
            r"\begin{proposition}[Coproduct primitivity at degree~$2$]",
            r"\begin{theorem}[Gravitational coproduct primitivity]",
            r"\Delta_{z,n}^W \;=\; 0",
            r"\text{for every } x \in W_k(\mathfrak{g})",
            "strict primitivity holds at all degrees for all principal DS reductions",
            "The two families exhaust the HPL expansion at each degree. Source trees have ghost~$+1$",
        )
        block = prop + theorem + source[
            source.index(r"\begin{remark}[Scope and the infinite tower]"):
            source.index(r"\end{remark}", source.index(r"\begin{remark}[Scope and the infinite tower]"))
        ]
        for needle in retired:
            assert needle not in block

    def test_good_graded_ds_e3_topological_scope_is_conditional(self):
        """The non-principal DS theorem is a good-graded package theorem."""
        scope = good_graded_ds_e3_topological_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'epsilon')
        assert scope['hypothesis_package'] == 'hypDSBRST'
        assert 'nilpotent f with good DS grading' in scope['requires']
        assert 'total DS differential Q_DS,f = Q_CS + Q_red,f' in scope['requires']
        assert scope['brst_identity_scope'] == 'on Q_DS,f-cohomology'
        assert scope['topological_structure_target'] == (
            'H^bullet_{Q_DS,f} Zder^ch(W^k(g,f))'
        )
        assert scope['strict_raw_chain_level_statement'] is False
        assert scope['unreduced_affine_current_qcs_exactness'] is False
        assert scope['normalisation_dependent_improvement_coefficients'] is True

    def test_good_graded_ds_manuscript_does_not_use_unreduced_current_proof(self):
        """The general DS theorem must not reuse the old Q_CS-current proof."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Good-graded DS cohomological $\Ethree$-topologisation"
        )
        end = source.index(r"\end{theorem}", start)
        theorem = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]
        bp_remark_start = source.index(r"\begin{remark}[Specialisation to Bershadsky", proof_end)
        bp_remark_end = source.index(r"\end{remark}", bp_remark_start)
        bp_remark = source[bp_remark_start:bp_remark_end]

        required_statement = (
            r"\ClaimStatusConditional",
            r"\hypDSBRST",
            r"Q_{\mathrm{DS},f}=Q_{\mathrm{CS}}+Q_{\mathrm{red},f}",
            r"H_{Q_{\mathrm{DS},f}}^\bullet\Zder^{\mathrm{ch}}(\cW)",
            r"\text{on $Q_{\mathrm{DS},f}$-cohomology.}",
        )
        for needle in required_statement:
            assert needle in theorem

        required_proof = (
            r"Q_{\mathrm{DS},f}\)-closed observable",
            r"Q_{\mathrm{DS},f}\text{-cohomology}",
            r"E_2^{\mathrm{top}}\otimes E_1^{\mathrm{top}}",
            r"does not use",
            r"unreduced bulk",
        )
        for needle in required_proof:
            assert needle in block

        assert r"Q_{\mathrm{DS},f_{\min}}" in bp_remark
        assert "normalisation-dependent" in bp_remark

        retired = (
            r"\label{eq:current-Q-exact-general}",
            r"\label{eq:T-imp-general}",
            r"\label{eq:G-prime-general}",
            r"J_a \;=\; [Q_{\mathrm{CS}},\, \bar c_a]",
            r"T_{\mathrm{DS}}(f) \;=\; [Q_{\mathrm{CS}},\, G'_f]",
            r"G'_{f_{\min}}(z) \;=\;",
        )
        for needle in retired:
            assert needle not in block
            assert needle not in bp_remark

    def test_gravity_finite_presentation_uses_full_hpl_datum(self):
        """The theorem must not compress the HPL transfer to bare (m2,h)."""
        root = Path(__file__).resolve().parents[2]
        gravity = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = gravity.index(r"\begin{theorem}[Finite presentation of the gravitational $\Ainf$")
        end = gravity.index(
            r"\begin{proposition}[Symmetric-point graviton tracelessness;", start
        )
        block = gravity[start:end]

        required = (
            r"\ClaimStatusConditional",
            r"\hypAmbientWtCpl",
            r"\hypKZSDR",
            r"\effKoszul",
            r"k_{\mathrm{aff}}",
            r"Q_{\mathrm{DS}}h_{\mathrm{DS}}+h_{\mathrm{DS}}Q_{\mathrm{DS}}",
            r"\mathrm{id}-i\circ p",
            r"\delta_{\mathrm{DS}}",
            "the finite DS/KZ HPL datum produces the entire infinite tower",
        )
        for needle in required:
            assert needle in block

        retired = (
            "uniquely determined by the two data",
            "one binary operation and one homotopy produce",
            r"level~$k$",
            r"presented by $(m_2, h)$",
        )
        for needle in retired:
            assert needle not in block

    def test_antighost_homotopy_rejects_nonpositive_weight(self):
        """The formula G_0/w is not a weight-zero homotopy."""
        with pytest.raises(ValueError):
            virasoro_hpl_sdr_data(weight=0)

    def test_s4_forces_completed_ambient(self):
        """Nonzero bar weight 4 is a completed/pro class-M statement."""
        for c_val in [1, 13, 26]:
            requirement = virasoro_s4_completed_ambient_requirement(c=c_val)
            assert requirement['S4'] != 0
            assert requirement['bar_weight'] == 4
            assert 'weight-completed/pro' in requirement['ambient']
            assert requirement['raw_direct_sum_chain_statement'] is False
            assert requirement['completed_statement'] is True
            assert requirement['singular_values'] == (S(0), Rational(-22, 5))


class TestExactGravityScope:
    """Verify the exact algebraic gravity statement is scoped correctly."""

    def test_completed_derived_centre_statement(self):
        """The rigorous statement is completed E_3-topological algebra."""
        profile = virasoro_exact_gravity_scope_profile(c=12)
        assert profile['boundary_algebra'] == 'Vir_c'
        assert profile['exact_algebraic_statement'] == (
            'Z_der^ch(Vir_c)^hat_rho in E_3^top-Alg'
        )
        assert profile['rho_condition'] == '0 < rho < |c|/6'
        assert profile['rho_bound'] == 2
        assert profile['completion'] == 'weight-completed/pro Banach completion'

    def test_brown_henneaux_is_dictionary_not_path_integral(self):
        """Brown-Henneaux supplies c, not a constructed metric path integral."""
        profile = virasoro_exact_gravity_scope_profile(c=13)
        assert profile['brown_henneaux_dictionary'] == 'c = 3 ell / (2 G_N)'
        assert profile['physical_reading'] == 'boundary-CFT holographic interpretation'
        assert profile['dynamical_metric_path_integral_constructed'] is False
        assert profile['required_physical_hypotheses'] == (
            'BHdict',
            'ModInv',
            'VacDom',
            'SadDom',
            'Borel/Stokes',
        )

    def test_brown_henneaux_chiral_chart_scope_profile(self):
        """The Brown-Henneaux theorem is an external chiral chart."""
        profile = brown_henneaux_chiral_chart_scope()
        assert profile['claim_status'] == 'ProvedElsewhere'
        assert profile['boundary_chart'] == 'b_BH'
        assert profile['licensing_tags'] == ('alpha', 'beta', 'gamma')
        assert profile['level'] == 'k = ell / (4 G_N)'
        assert profile['chiral_central_charge'] == (
            'c_ch = 6k = 3 ell / (2 G_N)'
        )
        assert profile['full_ads3_asymptotic_symmetry'] == (
            'Vir_c_ch direct_sum Vir_c_ch'
        )
        assert 'hypBHdict' in profile['requires']
        assert 'one oriented SL(2,R) Chern-Simons factor' in profile['requires']
        assert profile['used_as'] == 'external boundary chart A_{b_BH} = Vir_c_ch'
        assert profile['derived_from_bar_complex'] is False
        assert profile['dynamical_metric_path_integral_constructed'] is False

    def test_brown_henneaux_theorem_is_licensed_external_chart(self):
        """The manuscript theorem must not present BH as bar-complex output."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Brown--Henneaux chiral boundary chart;"
        )
        end = source.index(r"\end{theorem}", start)
        theorem = source[start:end]

        required = (
            r"\ClaimStatusProvedElsewhere",
            r"licensing $\alpha+\beta+\gamma$",
            r"\hypBHdict",
            r"\cite{BrownHenneaux,Witten88}",
            r"one oriented $SL(2,\R)$ Chern--Simons",
            r"c_{\mathrm{ch}} = 6k = \frac{3\ell}{2G_N}",
            "two commuting Virasoro",
            "not a theorem of the bar complex",
        )
        for needle in required:
            assert needle in theorem

        retired = (
            r"\begin{theorem}[Brown--Henneaux]",
            "The asymptotic symmetry algebra of one $SL(2,\\R)$ Chern--Simons",
            "$k = \\ell/(4G)$, $c = 6k = 3\\ell/(2G)$",
        )
        for needle in retired:
            assert needle not in theorem

    def test_gravitational_mc_bridge_scope_profile(self):
        """The MC package is conditional and separates boundary from bulk."""
        profile = gravitational_mc_bridge_scope()
        assert profile['claim_status'] == 'Conditional'
        assert profile['licensing_tags'] == (
            'alpha',
            'beta',
            'gamma',
            'delta',
            'epsilon',
        )
        assert profile['boundary_face'] == 'A_{b_BH} = Vir_{c_ch}'
        assert profile['bulk_after_comparison'] == 'Z_der^ch(Vir_{c_ch})'
        assert profile['closed_face_is_boundary_virasoro'] is False
        assert 'hypBHdict' in profile['requires']
        assert 'physics-bridge BV datum' in profile['requires']
        assert 'hypAmbientWtCpl' in profile['requires']
        assert 'effKoszul' in profile['requires']
        assert 'abstract spectral line braiding' in profile['licensed_projections']
        assert profile['collision_kernel'] == (
            'r_coll_Vir(z) = (c_ch/2)/z^3 + 2T/z'
        )
        assert 'closed face equals Vir_c' in profile['not_asserted']
        assert 'BTZ black holes are proved MC deformations' in profile['not_asserted']

    def test_gravitational_mc_theorem_is_scoped_and_not_closed_virasoro(self):
        """The manuscript must not identify the closed/bulk face with Vir_c."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Brown--Henneaux Virasoro MC bridge;"
        )
        end = source.index(r"\end{theorem}", start)
        theorem = source[start:end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing",
            r"$\alpha+\beta+\gamma+\delta+\varepsilon$",
            r"\hypBHdict+\hypAmbientWtCpl+\effKoszul",
            r"A_{b_{\mathrm{BH}}}=\mathrm{Vir}_{c_{\mathrm{ch}}}",
            r"\alpha_{\mathrm{grav}}\in\mc(\gSC_T)",
            "boundary/chiral face",
            "completed perturbative line category",
            "abstract open-colour line braiding",
            r"\Theta_{\mathrm{Vir}_{c_{\mathrm{ch}}}}",
            r"\Zder^{\mathrm{ch}}(\mathrm{Vir}_{c_{\mathrm{ch}}})",
            "not the boundary Virasoro algebra",
            "not part of this theorem",
        )
        for needle in required:
            assert needle in theorem

        retired = (
            r"\begin{theorem}[Gravitational MC element]",
            r"\item closed face $= \mathrm{Vir}_c$",
            r"Items \textup{(i)--(v)} are proved unconditionally",
            "The gauge-fixed BV-BRST complex of Chern--Simons theory is cubic",
            "so the one-loop finiteness input",
        )
        for needle in retired:
            assert needle not in theorem

    def test_btz_deformation_is_candidate_not_proved_geometry(self):
        """BTZ geometry is only represented after comparison data are fixed."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{definition}[Candidate BTZ deformation]")
        end = source.index(r"\end{definition}", start)
        definition = source[start:end]

        assert "After the Brown--Henneaux, Cardy, and BTZ comparison data" in definition
        assert "candidate MC deformation" in definition
        assert "represented in the algebraic model" in definition
        assert "A BTZ black hole of mass~$M$ and angular momentum~$J$ is a MC" not in definition

    def test_brown_henneaux_newton_constant_normalization_is_gn(self):
        """Active Vol II Brown-Henneaux surfaces use G_N, not bare G."""
        root = Path(__file__).resolve().parents[2]
        rel_paths = (
            "chapters/connections/3d_gravity.tex",
            "chapters/connections/thqg_perturbative_finiteness.tex",
            "chapters/connections/universal_holography_functor.tex",
            "chapters/connections/ym_synthesis_core.tex",
            "chapters/connections/celestial_holography_core.tex",
            "chapters/frame/preface.tex",
            "chapters/theory/introduction.tex",
        )
        retired = (
            r"3\ell/(2G)",
            r"\frac{3\ell}{2G}",
            r"\ell/(4G)",
            r"\frac{3\ell}{4G}",
            r"\frac{\ell}{32G}",
            r"8Gh",
            r"8G\ell",
            r"G \to 0",
            r"G/\ell",
            r"G^2/\ell^2",
        )
        for rel_path in rel_paths:
            text = (root / rel_path).read_text()
            for needle in retired:
                assert needle not in text, f"{needle} remains in {rel_path}"

    def test_scalar_bar_trace_is_not_metric_path_integral(self):
        """Z_bar^grav is a scalar boundary trace, not the physical path integral."""
        profile = virasoro_scalar_bar_trace_profile(c=18)
        assert 'Tr_{B_ch(Vir_c)^hat_rho}' in profile['formula']
        assert 'Theta_Vir_c^(g)' in profile['formula']
        assert profile['trace_space'] == 'B_ch(Vir_c)^hat_rho'
        assert profile['rho_bound'] == 3
        assert profile['object_type'] == 'scalar boundary trace'
        assert profile['is_metric_path_integral'] is False

    def test_maloney_witten_sum_is_not_phi_hol_output(self):
        """Phi_hol produces the scalar seed, not the MW modular orbit sum."""
        profile = maloney_witten_scope_profile()
        assert profile['bar_trace'] == 'Tr_Bord(Vir_c) exp(Theta)'
        assert profile['bar_trace_role'] == 'perturbative thermal-AdS/BTZ seed'
        assert 'SL2Z/Gamma_inf' in profile['maloney_witten_sum']
        assert profile['phi_hol_outputs_maloney_witten_sum'] is False
        assert profile['bar_trace_equals_maloney_witten_sum'] is False
        assert 'saddle set' in profile['requires_for_maloney_witten_sum']
        assert 'ensemble prescription' in profile['requires_for_maloney_witten_sum']

    def test_scalar_genus_tower_has_no_page_or_btz_stokes_jump(self):
        """The FP scalar tower is convergent and has no Page/BTZ Stokes jump."""
        profile = maloney_witten_scope_profile()
        assert profile['scalar_series_variable'] == 'hbar^2'
        assert profile['scalar_series_radius'] == '4*pi^2'
        assert profile['ordinary_gevrey1_borel_transform'] == 'entire'
        assert profile['page_or_btz_stokes_from_scalar_tower'] is False
        assert 'raw two-sector transseries' in profile['page_stokes_requires']
        assert 'saddle extraction' in profile['page_stokes_requires']

    def test_main_synopsis_blocks_mw_page_stokes_overreach(self):
        """Source guard for the A215 top-level synopsis repair."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "main.tex").read_text()
        flat = " ".join(source.split())
        assert "Maloney--Witten" in flat
        assert "not the Maloney--Witten modular orbit sum" in flat
        assert "no Page or BTZ Stokes jump follows from the scalar tower" in flat
        assert "Borel resummation in the $\\kappa_{\\mathrm{ch}} < 0$ sector" not in source


class TestPageCurveProfile:
    """Guard the conditional two-sector Page-profile model."""

    def test_page_curve_crossing_and_real_window(self):
        """The branch crossing gives t_P=3S_BH/13 in 0<c<26."""
        profile = page_curve_profile(13, 26)
        assert profile['hawking_rate'] == Rational(13, 6)
        assert profile['island_rate'] == Rational(13, 6)
        assert profile['t_page'] == 6
        assert profile['s_page'] == 13
        assert profile['t_evap'] == 12
        assert profile['crossing_balance'] == 0

        asymmetric = page_curve_profile(6, 26)
        assert asymmetric['hawking_rate'] == 1
        assert asymmetric['island_rate'] == Rational(10, 3)
        assert asymmetric['t_page'] == 6
        assert asymmetric['s_page'] == 6
        assert asymmetric['t_evap'] == Rational(39, 5)
        assert asymmetric['crossing_balance'] == 0

        c, s_bh = symbols('c S_BH')
        symbolic = page_curve_profile(c, s_bh)
        assert symbolic['t_page'] == 3 * s_bh / 13
        assert symbolic['s_page'] == c * s_bh / 26
        assert simplify(symbolic['crossing_balance']) == 0

    def test_page_curve_rejects_degenerate_real_window(self):
        """The displayed decreasing branch is not a real Page curve outside 0<c<26."""
        for bad_c in (0, 26, -1, 30):
            with pytest.raises(ValueError, match="0 < c < 26"):
                page_curve_profile(bad_c, 26)

        with pytest.raises(ValueError, match="S_BH > 0"):
            page_curve_profile(13, 0)

    def test_page_curve_scope_separates_transseries_from_scalar_tower(self):
        """The Page profile is conditional on the raw transseries and entropy functional."""
        scope = page_curve_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta', 'epsilon')
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'effScalarShadowProj' in scope['hypotheses']
        assert 'hypModularCardy' in scope['hypotheses']
        assert 'real same-family window 0<c<26' in scope['hypotheses']
        assert scope['t_page'] == '3*S_BH/13'
        assert scope['scalar_tower_alone_derives_page_curve'] is False
        assert scope['valid_for_c_ge_26_without_extra_hypothesis'] is False
        assert scope['post_evaporation_profile_asserted'] is False
        assert scope['borel_singularity_gives_exp_page_time_formula'] is False
        assert scope['line_verdier_sector_not_naive_koszul_slogan'] is True

    def test_page_curve_proposition_is_conditional_and_stokes_scoped(self):
        """Source guard for the conditional Page-profile repair."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\subsection{Conditional Page profile")
        end = source.index(r"\subsection{de~Sitter entropy", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\begin{proposition}[Conditional Page profile from the two-sector raw",
            r"\ClaimStatusConditional",
            r"$\hypAmbientWtCpl+\effScalarShadowProj+\hypModularCardy$",
            r"$0<c<26$",
            r"t_{\mathrm{evap}}",
            r"S_{\mathrm{rad}}^{\mathrm{model}}",
            r"$0\le t\le t_{\mathrm{evap}}$",
            "same-family line--Verdier comparison rate",
            r"\label{eq:page-stokes-wall}",
            "equal-real-part condition",
            r"\frac{I_{\mathrm{I}}(t)-I_{\mathrm{H}}(t)}{\hbar}",
            r"B(\cA_{\mathrm{line}})",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "not asserted outside the real window",
            "not a consequence of scalar complementarity alone",
            "Borel singularity is not a second formula for the Page time",
        )
        for needle in required_flat:
            assert needle in flat

        assert "$t=3S_{\\mathrm{BH}}/13$ in the real two-sector model" in flat
        assert "scalar shadow" in block

        retired = (
            r"\subsection{Page curve from gravitational Koszul duality}",
            r"\begin{proposition}[Gravitational Page curve]",
            "The island phase uses the Koszul dual rate",
            r"the Koszul dual saddle $B(\cA^!)$ dominates",
            r"dominance to $B(\cA^!)$",
            r"where the Koszul dual sector $B(\cA^!)$ provides",
            r"\label{eq:page-borel-scale}",
            r"\frac{1}{\hbar}\,e^{S_{\mathrm{BH}}}",
            r"\frac{1}{\hbar}\,e^{\zeta_P/\hbar}",
            "Page time is the corresponding non-perturbative scale",
            "Page curve from gravitational Koszul duality",
        )
        for needle in retired:
            assert needle not in block


class TestDesitterShadowScope:
    """Guard the de Sitter scalar shadow and entropy normalization."""

    def test_desitter_entropy_formula_three_ways(self):
        """Area law and central-charge normalization give the same entropy."""
        assert desitter_central_charge(2, 3) == 1
        assert desitter_horizon_entropy_from_radius(2, 3) == pi / 3
        assert desitter_horizon_entropy_from_c(1) == pi / 3

        ell, g_n = symbols('ell_dS G_N')
        c_expr = desitter_central_charge(ell, g_n)
        assert c_expr == 3 * ell / (2 * g_n)
        assert simplify(
            desitter_horizon_entropy_from_radius(ell, g_n)
            - desitter_horizon_entropy_from_c(c_expr)
        ) == 0

    def test_desitter_profile_records_scope_and_fixed_point(self):
        """The dS statement is a scalar real-section statement, not dS/CFT."""
        profile = desitter_shadow_profile(26, 3)
        assert profile['claim_status'] == 'Conditional'
        assert profile['licensing_tags'] == ('alpha', 'gamma', 'epsilon')
        assert 'de Sitter real-section metric normalization' in profile['hypotheses']
        assert 'hypAmbientWtCpl' in profile['hypotheses']
        assert 'effScalarShadowProj' in profile['hypotheses']
        assert profile['c_dS'] == 13
        assert profile['kappa_dS'] == Rational(13, 2)
        assert profile['entropy_from_radius'] == 13 * pi / 3
        assert profile['entropy_from_c'] == 13 * pi / 3
        assert profile['entropy_balance'] == 0
        assert profile['fixed_point_c'] == 13
        assert profile['fixed_point_kappa'] == Rational(13, 2)
        assert profile['fixed_point_entropy'] == 13 * pi / 3
        assert profile['nariai_geometry_constructed'] is False
        assert profile['desitter_hilbert_space_constructed'] is False
        assert profile['dscft_correlator_functor_constructed'] is False
        assert profile['banks_dimension_theorem'] is False
        assert profile['literal_complex_wick_rotation_gives_real_coefficients'] is False

    def test_desitter_profile_rejects_nonpositive_metric_data(self):
        """The real-section comparison requires positive radius and G_N."""
        with pytest.raises(ValueError, match="radius"):
            desitter_central_charge(0, 1)
        with pytest.raises(ValueError, match="Newton"):
            desitter_horizon_entropy_from_radius(1, 0)
        with pytest.raises(ValueError, match="central charge"):
            desitter_horizon_entropy_from_c(-1)

    def test_desitter_proposition_is_real_section_scoped(self):
        """Source guard for the A396 de Sitter repair."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\subsection{de~Sitter entropy")
        end = source.index(r"\subsection{JT gravity", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\begin{proposition}[de~Sitter scalar shadow and horizon-normalized",
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\gamma+\varepsilon$",
            r"$\hypAmbientWtCpl+\effScalarShadowProj$",
            r"c_{\mathrm{dS}}:=\frac{3\ell_{\mathrm{dS}}}{2G_N}",
            r"S_{\mathrm{dS}}",
            r"\frac{A_{\mathrm{hor}}}{4G_N}",
            r"\frac{2\pi\ell_{\mathrm{dS}}}{4G_N}",
            r"\frac{\pi\ell_{\mathrm{dS}}}{2G_N}",
            r"\frac{\pi c_{\mathrm{dS}}}{3}",
            r"F_g^{\mathrm{dS},\mathrm{sc}}",
            r"\kappa_{\mathrm{dS}}=\frac{13}{2}",
            r"S_{\mathrm{dS}}=\frac{13\pi}{3}",
            "No Nariai geometry",
            "no de~Sitter Hilbert space",
            "no dS/CFT",
            r"\ClaimStatusHeuristic",
            r"\dim(\mathcal H_{\mathrm{dS}})\stackrel{\mathrm{heur}}{=}",
        )
        for needle in required:
            assert needle in block

        assert "not a theorem of the Virasoro scalar shadow tower" in flat
        assert "literal complex Wick rotation without this real-section choice" in flat

        retired = (
            r"\begin{proposition}[de~Sitter shadow obstruction tower]",
            r"Under the Wick rotation $\ell \mapsto i\ell$ mapping",
            r"\pi\ell_{\mathrm{dS}}/(2G)",
            r"c_{\mathrm{dS}} = 3\ell_{\mathrm{dS}}/(2G)",
            r"F_g = (c_{\mathrm{dS}}/2)\lambda_g^{\mathrm{FP}}",
            r"The Nariai limit $c_{\mathrm{dS}} = 13$",
            "The Banks conjecture",
            "is consistent with\n the convergent shadow obstruction tower",
        )
        for needle in retired:
            assert needle not in block


class TestJTSchwarzianScope:
    """Guard the conditional Schwarzian/JT scalar comparison."""

    def test_wp_curve_density_conversion(self):
        """The compact sine curve gives the JT density only after continuation."""
        z, E = symbols('z E')
        y = jt_wp_spectral_curve_y(z)
        assert y.subs(z, 0) == 0
        assert diff(y, z).subs(z, 0) == Rational(1, 2)
        assert series(y, z, 0, 5).removeO() == z / 2 - pi**2 * z**3 / 3

        assert jt_disk_density(0) == 0
        assert jt_disk_density(1) == sinh(2 * pi) / (4 * pi**2)
        assert jt_wp_to_density_balance(E) == 0

    def test_jt_density_rejects_negative_energy(self):
        """The physical disk density is stated on the nonnegative energy ray."""
        with pytest.raises(ValueError, match="nonnegative"):
            jt_disk_density(-1)

    def test_jt_scope_blocks_scalar_overclaims(self):
        """The scope object separates comparison data from Vol II proof."""
        scope = jt_schwarzian_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta', 'epsilon')
        assert 'Schwarzian comparison datum' in scope['hypotheses']
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'effScalarShadowProj' in scope['hypotheses']
        assert 'JT spectral-curve normalization' in scope['hypotheses']
        assert 'energy-contour convention' in scope['hypotheses']
        assert 'matrix-integral or Stokes completion datum' in scope['hypotheses']
        assert scope['compact_curve'] == 'x=z^2, y_WP(z)=sin(2*pi*z)/(4*pi)'
        assert scope['physical_density'] == 'rho_0(E)=sinh(2*pi*sqrt(E))/(4*pi^2)'
        assert scope['branch_of_sqrt_x_required'] is True
        assert scope['physical_density_requires_contour'] is True
        assert scope['scalar_shadow_alone_proves_wp_volumes'] is False
        assert scope['scalar_shadow_alone_completes_jt_series'] is False
        assert scope['physical_jt_density_from_sine_curve_without_contour'] is False
        assert scope['bernoulli_decay_is_jt_borel_completion'] is False
        assert scope['schur_q_series_supplies_uv_completion'] is False
        assert scope['miki_stress_sheaf_derives_jt_measure'] is False

    def test_jt_proposition_is_conditional_and_contour_scoped(self):
        """Source guard for the Schwarzian/JT repair."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\subsection{JT gravity")
        end = source.index(r"\begin{remark}[Class-$\mathcal{S}$", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\begin{proposition}[Conditional Schwarzian/JT scalar limit;",
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma+\delta+\varepsilon$",
            r"$\hypAmbientWtCpl+\effScalarShadowProj$",
            r"\mathfrak s_{\mathrm{JT}}",
            r"x=z^2",
            r"y_{\mathrm{WP}}(z)=\frac{\sin(2\pi z)}{4\pi}",
            r"z=i\sqrt{E}",
            r"\rho_0(E)=\frac{1}{i\pi}\,y_{\mathrm{WP}}(i\sqrt{E})",
            r"\frac{\sinh(2\pi\sqrt{E})}{4\pi^2}",
            r"\hypStokes",
            r"\begin{remark}[Normalization check]",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "not a proof of the Weil--Petersson volume theorem",
            "not a consequence of Theorem",
            "not an analytic completion of the divergent JT expansion",
            "requires the matrix-integral or Stokes datum",
            "does not construct the non-perturbative JT measure",
        )
        for needle in required_flat:
            assert needle in flat

        retired = (
            r"\begin{proposition}[JT shadow correspondence]",
            r"In the Schwarzian limit $c \to \infty$:",
            "The shadow metric degenerates to the JT spectral curve",
            r"y = \sin(2\pi\sqrt{x})/(4\pi)",
            r"The Weil--Petersson volumes $V_{g,1}(b)$ are reproduced",
            "providing an\n analytic completion",
            "The statement is the Schwarzian limit of the displayed shadow",
        )
        for needle in retired:
            assert needle not in block

    def test_class_s_jt_remark_no_longer_promotes_scalar_lane(self):
        """Downstream bridge must not repeat the old JT completion claim."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{remark}[Schwarzian/JT comparison shadow")
        end = source.index(r"\begin{remark}[What \(K^{\kappaChHodge}=8\)", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            "conditional Schwarzian comparison",
            r"\mathfrak s_{\mathrm{JT}}",
            r"x=z^2,\qquad y_{\mathrm{WP}}(z)=\frac{\sin(2\pi z)}{4\pi}",
            r"\rho_0(E)=\frac{\sinh(2\pi\sqrt E)}{4\pi^2}",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "not a derivation of the JT curve or the JT measure",
            "does not converge the asymptotically divergent",
            "requires the matrix-integral or Stokes datum",
        )
        for needle in required_flat:
            assert needle in flat

        retired = (
            "the limit reproduces the JT spectral curve",
            r"y = \sin(2\pi\sqrt{x})/(4\pi)",
            "Cardy asymptotic of the Schur",
            "supplies the UV completion",
        )
        for needle in retired:
            assert needle not in flat


class TestK3BorcherdsScalarBridge:
    """Guard the K3 x E Borcherds product as a scalar theorem."""

    def test_colored_partition_number_gottsche_values(self):
        """Göttsche specialization gives p_24(5)=176256 directly."""
        expected = [1, 24, 324, 3200, 25650, 176256]
        got = [colored_partition_number(24, n) for n in range(6)]
        assert got == expected

        with pytest.raises(ValueError, match="colors"):
            colored_partition_number(-1, 5)
        with pytest.raises(ValueError, match="degree"):
            colored_partition_number(24, -1)

    def test_borcherds_scalar_scope_records_normalizations(self):
        """The scalar bridge separates BPS, OP, Hilbert, and Bruinier data."""
        scope = k3_borcherds_scalar_bridge_scope()
        assert scope['claim_status'] == 'ProvedElsewhere'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'epsilon')
        assert 'Gritsenko-Nikulin half-K3 Jacobi normalization' in scope['hypotheses']
        assert 'Borcherds singular-theta lift' in scope['hypotheses']
        assert 'DMVV second-quantized K3 elliptic-genus scalar' in scope['hypotheses']
        assert 'Gottsche Hilbert-scheme specialization' in scope['hypotheses']
        assert 'Bruinier Heegner projection in the half-K3 convention' in scope['hypotheses']
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'effScalarShadowProj' in scope['hypotheses']
        assert scope['phi10_un'] == 'Delta_5^2'
        assert scope['d5_factor'] == Rational(1, 64)
        assert scope['phi10_op_factor'] == Rational(1, 4096)
        assert scope['bps_scalar_prefactor_delta5_inverse_square'] == 1
        assert scope['op_dt_scalar_prefactor_delta5_inverse_square'] == -4096
        assert scope['p24_5'] == 176256
        assert scope['hilb5_k3_euler'] == 176256
        assert scope['bruinier_reduced_c2_triangle'] == 0
        assert scope['bruinier_reduced_c3_half_k3'] == -64
        assert scope['phi_minus2_1_reduced_c3_other_convention'] == -8
        assert scope['uses_phi_minus2_1_reduced_convention'] is False
        assert scope['p24_5_is_bruinier_obstruction_coefficient'] is False
        assert scope['scalar_identity_promotes_to_gravity_line_trace'] is False
        assert scope['scalar_identity_proves_maloney_witten_equality'] is False
        assert scope['one_variable_hilbert_series_is_full_siegel_form'] is False
        assert scope['op_normalization_changes_bkm_denominator_algebra'] is False

    def test_borcherds_scalar_theorem_is_licensed_and_no_promotion(self):
        """Source guard for the Borcherds scalar bridge theorem."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\section{The Borcherds-product scalar bridge}")
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        flat = " ".join(block.split())

        required = (
            r"\begin{theorem}[K3$\times E$ Borcherds scalar shadow;",
            r"\ClaimStatusProvedElsewhere",
            r"licensing $\beta+\gamma+\varepsilon$",
            "Gritsenko--Nikulin half-K3 normalisation",
            r"$\hypAmbientWtCpl+\effScalarShadowProj$",
            r"\Phi_{10}^{\mathrm{un}}(Z)",
            r"\Delta_5(Z)^2",
            r"D_5=64^{-1}\Delta_5",
            r"\Phi_{10}^{\mathrm{OP}}=D_5^2=4096^{-1}\Delta_5^2",
            r"Z_{\mathrm{OP/DT}}=-(\Phi_{10}^{\mathrm{OP}})^{-1}",
            r"Z^{K3\times E}_{\mathrm{BPS}}(Z)",
            r"\Delta_5(Z)^{-2}",
            r"p_{24}(5)",
            "176256",
            r"c_2^{\triangle}=0",
            r"c_3^{\mathrm{Br}}=-64[H_3]",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "not a BV or Bruinier obstruction coefficient",
            "half-K3 Gritsenko--Nikulin scalar convention",
            "not equal to the full three-variable Siegel modular form",
        )
        for needle in required_flat:
            assert needle in flat

        retired = (
            r"\begin{theorem}[K3$\times E$ Borcherds scalar shadow]",
            "This coefficient is a BV",
            "This coefficient is a Bruinier obstruction coefficient",
            "one-variable Hilbert-scheme series is equal to the full",
        )
        for needle in retired:
            assert needle not in flat

    def test_scalar_bridge_residual_blocks_operator_promotion(self):
        """The following caveat must continue to block scalar-to-operator promotion."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{remark}[Scalar bridge residual]")
        end = source.index(r"\begin{proposition}[K3 scalar no-promotion;", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            "is a scalar theorem",
            "proves neither",
            r"\tr_{B^{\mathrm{ord}}(\Vir_c)}",
            r"\bigl(\Phi_{10}^{\mathrm{un}}\bigr)^{-1}",
            "nor a Maloney--Witten equality",
            "gravity-line-hall-borcherds-comparison",
            "does not prove the square",
        )
        for needle in required:
            assert needle in flat

    def test_k3_scalar_no_promotion_scope_records_missing_comparisons(self):
        """The no-promotion theorem is scalar non-faithfulness, not a slogan."""
        scope = k3_scalar_no_promotion_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'epsilon')
        assert scope['ambient'] == 'hypAmbientWtCpl'
        assert scope['effectiveness'] == 'effScalarShadowProj'
        assert scope['uses_scalar_non_faithfulness_lemma'] is True
        assert scope['scalar_character_functor_faithful'] is False
        assert scope['acyclic_filtered_summand_preserves_scalar_character'] is True
        assert scope['acyclic_filtered_summand_preserves_chain_object'] is False
        assert scope['scalar_data']['phi10_un'] == 'Delta_5^2'
        assert scope['scalar_data']['z_bps'] == '(Phi_10^{un})^{-1}'
        assert scope['scalar_data']['kappa_bkm'] == 5
        assert scope['scalar_data']['kappa_tuple_fourfold'] == (0, 3, 5, 24)
        assert scope['missing_beta_comparison_data'] == (
            'positive-half Hall-Borcherds E1-chiral bialgebra morphism',
            'completed Drinfeld-double extension through the current envelope on E',
            'filtered SC^{ch,top} morphism to the Virasoro gravity-line boundary algebra',
            'derived-centre trace compatibility with the BPS scalar character',
        )
        assert scope['determines_filtered_sc_chtop_morphism'] is False
        assert scope['determines_ordered_virasoro_bar_trace'] is False
        assert scope['determines_maloney_witten_sum'] is False
        assert scope['scalar_shadow_is_object_equivalence'] is False

    def test_k3_scalar_no_promotion_proposition_is_licensed_and_uses_lemma(self):
        """Source guard for prop:3dg-k3-scalar-no-promotion."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{proposition}[K3 scalar no-promotion;")
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing \(\beta+\gamma+\varepsilon\)",
            r"\hypAmbientWtCpl",
            r"\effScalarShadowProj",
            r"\label{prop:3dg-k3-scalar-no-promotion}",
            r"\Phi_{10}^{\mathrm{un}}=\Delta_5^2",
            r"Z_{\mathrm{BPS}}^{K3\times E}=\bigl(\Phi_{10}^{\mathrm{un}}\bigr)^{-1}",
            r"c_{\phi_{0,1}^{K3}}(0)/2=5",
            r"\{\kappa_{\mathrm{cat}},\kappaChHodge^{\mathrm{Heis}},",
            r"=\{0,3,5,24\}",
            r"\cref{lem:scalar-non-faithfulness}",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "do not determine a filtered",
            "nor a chain-level identification with the ordered Virasoro bar trace",
            "at least the following \\(\\beta\\)-comparison data",
            "positive-half Hall--Borcherds comparison",
            "completed Drinfeld double",
            "filtered \\(\\SCchtop\\)-morphism",
            "derived chiral centres",
            "acyclic filtered summand",
            "Items \\textup{(i)}--\\textup{(iv)} are exactly those missing",
        )
        for needle in required_flat:
            assert needle in flat

        retired = (
            r"\begin{proposition}[K3 scalar no-promotion]",
            "These scalar data determine a filtered",
            "scalar shadow to chain-level gravity requires no",
        )
        for needle in retired:
            assert needle not in flat


class TestHeptagonGrowthBound:
    """Guard the finite-type heptagon Gevrey-1 coefficient bound."""

    def test_heptagon_growth_constant_depends_on_generator_rank(self):
        """The constant has an essential R factor."""
        assert heptagon_growth_bound_constant(3, Rational(2, 5)) == (
            12 * exp(pi * sqrt(Rational(2, 3)))
        )
        assert heptagon_growth_bound_constant(2, 7) == (
            56 * exp(pi * sqrt(Rational(2, 3)))
        )

        with pytest.raises(ValueError, match="rank"):
            heptagon_growth_bound_constant(0, 1)
        with pytest.raises(ValueError, match="norm"):
            heptagon_growth_bound_constant(1, -1)

    def test_heptagon_growth_scope_is_local_not_borel_global(self):
        scope = heptagon_growth_bound_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('gamma', 'epsilon')
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'effKoszul' in scope['hypotheses']
        assert 'finite strong-generator profile (W,R)' in scope['hypotheses']
        assert scope['constant_parameters'] == ('W', 'R', 'M')
        assert scope['constant_formula'] == '4*R*max(1,M)*exp(pi*sqrt(2/3))'
        assert scope['depends_only_on_W_M'] is False
        assert scope['requires_generator_rank'] is True
        assert scope['tree_shape_bound'] == 'Catalan_n <= 4^n'
        assert scope['weight_decomposition_bound'] == 'p(n) <= exp(pi*sqrt(2*n/3))'
        assert scope['partition_absorbed_exponentially'] is True
        assert scope['pbw_symmetrization_bound'] == 'n!'
        assert scope['borel_transform_local_radius'] == '>= 1/C(W,R,M)'
        assert scope['proves_sectorial_continuation'] is False
        assert scope['proves_borel_singularity_location'] is False
        assert scope['proves_maloney_witten_interpretation'] is False

    def test_heptagon_growth_bound_source_is_finite_type_and_licensed(self):
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{proposition}[Heptagon growth bound")
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"\hypAmbientWtCpl+\effKoszul",
            "finite strong-generator profile",
            r"G=\bigoplus_{1\le j\le W}G_j",
            r"R=\dim G<\infty",
            "Lyndon--PBW basis",
            r"C=C(W,R,M)",
            r"C(W,R,M)=4R\max(1,M)e^{\pi\sqrt{2/3}}",
            r"4^n R^n \max(1,M)^n p(n)\,n!",
            "The dependence on \\(R\\) is essential",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "weight decompositions are bounded by the partition function",
            "ordered Lie--PBW symmetrisation contributes at most \\(n!\\)",
            "without a finite strong-generator profile the number of weight-\\(\\le W\\) labels is not controlled",
        )
        for needle in required_flat:
            assert needle in flat

        retired = (
            r"\begin{proposition}[Heptagon growth bound on ordered bar coefficients]",
            r"C(W, M)",
            "depending only on the generation weight and OPE norm",
            "The factor \\(n!\\) cannot be improved",
        )
        for needle in retired:
            assert needle not in block

    def test_heptagon_growth_bound_downstream_does_not_reuse_two_parameter_constant(self):
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{proposition}[Factorisation $L_0$-spectrum")
        end = source.index(r"\end{proof}", start)
        block = source[start:end]

        assert "finite strong-generator profile" in block
        assert "C_{\\mathrm{fact}}(W,U,R,D,M,P)^{n}\\, n!" in block
        assert "C(W,R,M)^{n}\\, n!" not in block
        assert "same \\(C(W,R,M)" not in block
        assert "C(W, M)" not in block


class TestFactorizationL0Gevrey:
    """Guard the finite-profile factorisation-BV Gevrey-1 bound."""

    def test_factorization_constant_has_bv_rank_and_propagator_norm(self):
        assert factorization_l0_gevrey_constant(3, 5, Rational(2, 5), 7) == (
            840 * exp(pi * sqrt(Rational(2, 3)))
        )
        assert factorization_l0_gevrey_constant(2, 4, 3, Rational(1, 2)) == (
            192 * exp(pi * sqrt(Rational(2, 3)))
        )
        with pytest.raises(ValueError, match="bar generator"):
            factorization_l0_gevrey_constant(0, 1)
        with pytest.raises(ValueError, match="BV generator"):
            factorization_l0_gevrey_constant(1, 0)
        with pytest.raises(ValueError, match="OPE"):
            factorization_l0_gevrey_constant(1, 1, -1, 1)
        with pytest.raises(ValueError, match="propagator"):
            factorization_l0_gevrey_constant(1, 1, 1, -1)

    def test_factorization_scope_is_local_not_resurgent(self):
        scope = factorization_l0_gevrey_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('gamma', 'epsilon')
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'effKoszul' in scope['hypotheses']
        assert 'finite factorisation-BV graph profile (U,D)' in scope['hypotheses']
        assert scope['constant_parameters'] == ('W', 'U', 'R', 'D', 'M', 'P')
        assert scope['constant_formula'] == (
            '8*R*D*max(1,M,P)*exp(pi*sqrt(2/3))'
        )
        assert scope['same_as_heptagon_constant'] is False
        assert scope['splitting_bound'] == 'n+1 <= 2^n'
        assert scope['bar_tree_shape_bound'] == '4^j'
        assert scope['bv_graph_skeleton_bound'] == '4^i'
        assert scope['factorial_absorption'] == 'i!*j! <= n!'
        assert scope['borel_transform_local_radius'] == (
            '>= 1/C_fact(W,U,R,D,M,P)'
        )
        assert scope['proves_sectorial_continuation'] is False
        assert scope['proves_borel_singularity_location'] is False
        assert scope['proves_bcov_instanton'] is False
        assert scope['proves_maloney_witten_interpretation'] is False

    def test_factorization_l0_proposition_is_finite_profile_and_licensed(self):
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{proposition}[Factorisation $L_0$-spectrum")
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing \(\gamma+\varepsilon\)",
            r"\hypAmbientWtCpl+\effKoszul",
            "finite factorisation-BV graph profile",
            r"H_{\mathrm{BV}}=\bigoplus_{1\le j\le U}H_{\mathrm{BV},j}",
            r"D=\dim H_{\mathrm{BV}}<\infty",
            "regularised vertex/propagator kernel norm bounded by \\(P\\)",
            r"C_{\mathrm{fact}}(W,U,R,D,M,P)^{n}\, n!",
            r"C_{\mathrm{fact}}(W,U,R,D,M,P)",
            r"=8RD\,\max(1,M,P)\,e^{\pi\sqrt{2/3}}",
            "This is a local Gevrey-\\(1\\) coefficient bound",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "There are \\(n+1\\le2^n\\) such splittings",
            "finite graph profile gives at most \\(D^i\\)",
            "regularised propagator/vertex norm contributes",
            "loop-order automorphism denominators only reduce",
            "\\(i!\\,j!\\le n!\\)",
            "did not construct sectorial Borel continuation",
        )
        for needle in required_flat:
            assert needle in flat

        retired = (
            r"\begin{proposition}[Factorisation $L_0$-spectrum Gevrey-$1$ bound]",
            "with the same \\(C(W,R,M)=4R\\max(1,M)e^{\\pi\\sqrt{2/3}}\\)",
            "the three factors combine exactly as in the heptagon proof",
            "yielding the same\n\\(C(W,R,M)\\) constant",
            "factorisation BV presentation and the\nheptagon chain-level presentation are two faces",
            "equivalently\n\\autoref{prop:3dg-fact-bv-l0-gevrey}",
        )
        for needle in retired:
            assert needle not in block
            assert needle not in source


class TestFactorizationBcovInstanton:
    """Guard the conditional BCOV saddle datum and normalization."""

    def test_location_and_residue_are_separate_data(self):
        """The Borel location is an action difference, not the beta coefficient."""
        c = Symbol('c')
        assert factorization_bcov_candidate_location() == S.One / (2 * pi)
        assert simplify(
            factorization_bcov_candidate_residue(c=c) - c * (c - 25) / 24
        ) == 0
        assert factorization_bcov_candidate_residue(c=24) == -1
        assert factorization_bcov_candidate_residue(c=25) == 0
        assert factorization_bcov_candidate_residue(c=0) == 0
        assert factorization_bcov_candidate_residue(c=13) == Rational(-13, 2)

    def test_scope_is_conditional_stokes_data_not_local_qme(self):
        """Local BV/QME data do not supply the saddle theorem."""
        scope = factorization_bcov_instanton_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'gamma', 'delta')
        assert 'Borel-coordinate normalisation' in scope['hypotheses']
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'hypStokes' in scope['hypotheses']
        assert 'nondegenerate factorisation-BV saddle' in scope['hypotheses']
        assert 'saddle-to-ordered-bar comparison' in scope['hypotheses']
        assert scope['action_difference'] == (
            'A_BV(Phi_*)=S_eff(Phi_*)-S_eff(Phi_0)'
        )
        assert scope['candidate_location'] == '1/(2*pi)'
        assert scope['candidate_residue'] == 'c*(c-25)/24'
        assert scope['action_location_and_residue_are_separate'] is True
        assert scope['action_equals_beta_residue'] is False
        assert scope['old_action_formula_asserted'] is False
        assert scope['local_qme_proves_singularity'] is False
        assert scope['gevrey_bound_proves_singularity'] is False
        assert scope['requires_sectorial_comparison'] is True
        assert scope['requires_one_loop_determinant'] is True
        assert scope['proves_maloney_witten_interpretation'] is False

    def test_bcov_proposition_is_statused_and_separates_action_from_residue(self):
        """The manuscript must not identify the beta coefficient with the action."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Factorisation-BV BCOV saddle datum"
        )
        prop_end = source.index(r"\end{proposition}", start)
        proof_end = source.index(r"\end{proof}", prop_end)
        block = source[start:proof_end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing",
            r"$\alpha+\gamma+\delta$",
            r"$\hypAmbientWtCpl+\hypStokes$",
            r"A_{\mathrm{BV}}(\Phi_\ast)",
            r"S^{\mathrm{eff}}[\Phi_\ast]-S^{\mathrm{eff}}[\Phi_0]",
            r"A_{\mathrm{BV}}(\Phi_\ast)=\frac{1}{2\pi}",
            r"\rho_c=\frac{c(c-25)}{24}",
            r"(\zeta_\ast,\operatorname{Res}_{\zeta_\ast})",
            "action location and the one-loop residue are separate data",
            "\\(c(c-25)/24\\) is not the instanton action",
            "saddle-to-ordered-bar comparison",
            "Maloney--Witten interpretation",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "After the Borel coordinate is fixed",
            "an exponential sector \\(e^{-A/\\hbar}\\)",
            "singularity at \\(\\zeta=A\\)",
            "separate determinant computation fixes the residue",
            "do not prove the existence of \\(\\Phi_\\ast\\)",
            "neither the nondegenerate solution, nor sectorial continuation",
        )
        for needle in required_flat:
            assert needle in flat

        retired = (
            r"S^{\mathrm{eff}}[\Phi_\ast]=2\pi\cdot c(c-25)/24",
            "then the corresponding candidate Borel-plane action is",
            "The candidate residue is \\(\\rho_c=c(c-25)/24\\)",
            "Under the hypotheses, the location and residue are the standard",
            "singularities of the Borel transform occur at saddle actions",
        )
        for needle in retired:
            assert needle not in block

    def test_alien_derivative_theorem_uses_saddle_datum(self):
        """The downstream alien derivative theorem must invoke the datum."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Alien-derivative closure on factorisation"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        assert r"\hypStokes" in block
        assert r"\autoref{prop:3dg-fact-bv-bcov-instanton}" in block
        assert r"\zeta_{\ast} = 1/(2\pi)" in block
        assert r"\rho_{c}\cdot\widehat{\mathcal F}_{A^{!}}" in block


class TestHeisenbergZwegersShadow:
    """Guard the rank-one inverse-eta Heisenberg shadow theorem."""

    def test_heisenberg_shadow_scope_is_fock_not_eta_power(self):
        scope = heisenberg_zwegers_shadow_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('gamma', 'epsilon')
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'effScalarShadowProj' in scope['hypotheses']
        assert 'rank-one standard Fock vacuum trace' in scope['hypotheses']
        assert 'nonzero Heisenberg level' in scope['hypotheses']
        assert scope['raw_character'] == 'eta(tau)^(-1)'
        assert scope['coefficient_sequence'] == 'p(n)'
        assert scope['modular_weight'] == Rational(-1, 2)
        assert scope['zwegers_shadow_vanishes'] is True
        assert scope['mock_completion_required'] is False
        assert scope['borel_transform'] == 'entire of exponential type zero'
        assert scope['finite_borel_singularities'] is False
        assert scope['level_changes_oscillator_count'] is False
        assert scope['unbranched_real_eta_power_modular_asserted'] is False
        assert scope['determinant_line_power_requires_branch_multiplier'] is True
        assert scope['determinant_power_weight_formula'] == (
            '-r/2 after branch/multiplier datum'
        )
        assert scope['proves_maloney_witten_path_integral'] is False

    def test_heisenberg_shadow_source_is_licensed_and_scoped(self):
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\subsection{Zwegers shadow for the Heisenberg chiral algebra}"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing \(\gamma+\varepsilon\)",
            r"\hypAmbientWtCpl+\effScalarShadowProj",
            r"k\in\R^{\times}",
            "standard Fock vacuum trace",
            "inverse eta multiplier",
            r"\mathcal F^{\mathrm{raw}}_{\mathcal H_k}(\tau)=\eta(\tau)^{-1}",
            "does not assert modularity for unbranched real powers",
            "rescaling the generator",
            "partitions of \\(n\\)",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "A determinant-line power \\(\\eta(\\tau)^{-r}\\) is a different scalar shadow",
            "requires a branch and multiplier datum before it has weight \\(-r/2\\)",
            "completed orbit object",
            "Rademacher expansion of a weight-\\(-1/2\\) modular form",
            "entire of exponential type zero",
            "no finite Borel-plane singularity or Stokes residue",
            "zero image under the Zwegers shadow operator",
        )
        for needle in required_flat:
            assert needle in flat

        retired = (
            r"\begin{proposition}[Zwegers shadow for free-boson Heisenberg]",
            r"For $A = \mathcal H_{k}$ at any level $k \in \R$",
            r"The determinant-line shadow with exponent $k$ is $\eta(\tau)^{-k}$.",
            r"The Borel transform is entire, the Laplace integral is absolutely",
            r"convergent on $\{\mathrm{Re}\,\zeta > 0\}$",
            "Classical modularity implies the shadow vanishes",
        )
        for needle in retired:
            assert needle not in block

    def test_heisenberg_shadow_downstream_bridge_uses_fock_scope(self):
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{conjecture}[Off-Koszul Maloney--Witten bridge")
        end = source.index(r"\end{conjecture}", start)
        block = source[start:end]
        assert "standard rank-one\nFock trace and inverse eta multiplier" in block
        assert "For $A = \\mathcal H_{k}$ at any level" not in block

        assembly_start = source.index(r"\begin{remark}[Conditional assembly]")
        assembly_end = source.index(r"\end{remark}", assembly_start)
        assembly = source[assembly_start:assembly_end]
        assert "rank-one inverse-eta Fock trace" in assembly
        assert "exponential type zero" in assembly
        assert "Laplace integral converges in every direction" not in assembly


class TestVirasoroHardyRamanujanCardy:
    """Guard the Virasoro vacuum coefficient formula and Cardy separation."""

    def test_vacuum_verma_coefficients_are_partitions_without_part_one(self):
        expected = [1, 0, 1, 1, 2, 2, 4, 4, 7, 8, 12, 14]
        got = [virasoro_vacuum_verma_coefficient(n) for n in range(12)]
        assert got == expected
        for n in range(1, 18):
            assert virasoro_vacuum_verma_coefficient(n) == (
                partition_number(n) - partition_number(n - 1)
            )
        with pytest.raises(ValueError, match="degree"):
            partition_number(-1)
        with pytest.raises(ValueError, match="degree"):
            virasoro_vacuum_verma_coefficient(-1)

    def test_hardy_ramanujan_prefactor_is_partition_difference_prefactor(self):
        """The ratio is below 1 at finite n and increases toward the leading term."""
        ratios = []
        for n in (50, 100, 200):
            exact = virasoro_vacuum_verma_coefficient(n)
            asym = float(virasoro_vacuum_verma_asymptotic(n).evalf())
            ratios.append(exact / asym)
        assert ratios[0] < ratios[1] < ratios[2] < 1
        assert ratios[0] > 0.75
        assert ratios[2] > 0.85
        with pytest.raises(ValueError, match="positive"):
            virasoro_vacuum_verma_asymptotic(0)

    def test_scope_separates_verma_coefficients_from_cardy_density(self):
        scope = virasoro_hardy_ramanujan_cardy_scope()
        assert scope['claim_status_vacuum_coefficients'] == 'ProvedHere'
        assert scope['claim_status_cardy_density'] == 'ProvedElsewhere'
        assert scope['licensing_tags'] == ('gamma', 'delta')
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'hypModularCardy for physical CFT comparison' in scope['hypotheses']
        assert scope['vacuum_character'] == 'prod_{m>=2}(1-q^m)^(-1)'
        assert scope['coefficient_formula'] == 'p(n)-p(n-1)'
        assert scope['asymptotic_prefactor'] == 'pi/(12*sqrt(2))'
        assert scope['asymptotic_exponential'] == 'exp(pi*sqrt(2*n/3))'
        assert scope['independent_of_c'] is True
        assert scope['cardy_is_verma_coefficient_identity'] is False
        assert scope['bare_verma_module_has_modular_invariance'] is False
        assert scope['vacuum_pqw_equals_physical_density'] is False
        assert scope['zwegers_shadow_from_verma_coefficients_alone'] is False

    def test_hardy_ramanujan_cardy_proposition_is_licensed_and_proved(self):
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Virasoro vacuum Verma coefficient growth;"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"\ClaimStatusProvedElsewhere",
            r"\rho_{\mathrm{vac}}(n)",
            r"p(-1)=0",
            r"p(n)-p(n-1)",
            r"\frac{\pi}{12\sqrt2}",
            r"\prod_{m\ge2}(1-q^m)^{-1}",
            r"\cite{HardyRamanujan1918}",
            r"\cite{Cardy86}",
        )
        for needle in required:
            assert needle in block

        required_flat = (
            "only vacuum null is the translation null",
            "partitions of \\(n\\) containing a part~\\(1\\)",
            "Cardy statement is not a Verma-module coefficient identity",
            "does not follow from \\eqref{eq:3dg-hardy-ramanujan-cardy}",
            "modular covariance, vacuum dominance, and the physical effective central charge",
            "not present in the bare universal Virasoro vacuum module",
            "licensing \\(\\gamma+\\delta\\)",
            "\\(\\hypAmbientWtCpl+\\hypModularCardy\\)",
        )
        for needle in required_flat:
            assert needle in flat

        retired = (
            r"\begin{proposition}[Virasoro vacuum-module coefficient growth and",
            "by the Hardy--Ramanujan asymptotic for the partition function. If one",
            "The second statement is not a Verma-module coefficient identity.",
            "matching\n\\autoref{prop:3dg-hardy-ramanujan-cardy}",
        )
        for needle in retired:
            assert needle not in block

    def test_downstream_zwegers_text_uses_cardy_clause_not_verma_formula(self):
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        theorem_start = source.index(
            r"\begin{theorem}[Explicit Zwegers $\mu$-shadow for $\Vir_c$"
        )
        theorem_end = source.index(r"\end{proof}", theorem_start)
        theorem_block = source[theorem_start:theorem_end]
        assert "the \\(\\hypModularCardy\\) clause" in theorem_block
        assert "not the bare vacuum Verma\ncoefficient formula" in theorem_block
        assert "Eichler-derivative recovers the\nphysical Cardy density" in theorem_block

        remark_start = source.index(r"\begin{remark}[What remains conditional]")
        remark_end = source.index(r"\end{remark}", remark_start)
        remark = source[remark_start:remark_end]
        flat = " ".join(remark.split())
        assert "not a consequence of the vacuum Verma coefficient formula alone" in flat
        assert "PBW--Hardy--Ramanujan computation" in flat
        assert "\\(\\hypModularCardy\\) supplies the physical Cardy tail" in flat
        assert "Hardy--Ramanujan--Cardy density" not in source


# ===================================================================
# 5. QUARTIC CONTACT INVARIANT
# ===================================================================

class TestQuarticContact:
    """Verify Q^contact = 10/(c(5c+22))."""

    def test_quartic_formula_c1(self):
        """At c=1: Q = 10/(1*27) = 10/27."""
        assert quartic_contact_virasoro(c=1) == Rational(10, 27)

    def test_quartic_formula_c13(self):
        """At comparison-fixed c=13: Q = 10/(13*87) = 10/1131."""
        assert quartic_contact_virasoro(c=13) == Rational(10, 1131)

    def test_quartic_positive_for_c_gt_0(self):
        """Q > 0 for c > 0."""
        for c_val in [1, 2, 5, 13, 26, 100]:
            assert quartic_contact_virasoro(c=c_val) > 0

    def test_quartic_exact_arithmetic(self):
        """Exact rational arithmetic via Fraction."""
        q = quartic_contact_virasoro_exact(1)
        assert q == Fraction(10, 27)

    def test_quartic_poles(self):
        """Poles at c=0 and c=-22/5."""
        c = Symbol('c')
        q = quartic_contact_virasoro()
        denom = c * (5 * c + 22)
        # Denominator vanishes at c=0 and c=-22/5
        assert denom.subs(c, 0) == 0
        assert denom.subs(c, Rational(-22, 5)) == 0

    def test_quartic_cross_family(self):
        """Verify Vol I ground truth: Q^contact(Vir_c) = 10/(c(5c+22))."""
        # This is the DEFINING formula — verify structure
        c = Symbol('c')
        q = quartic_contact_virasoro()
        assert simplify(q - 10 / (c * (5 * c + 22))) == 0


# ===================================================================
# 6. GRAVITATIONAL KOSZUL TRIANGLE
# ===================================================================

class TestGravitationalKoszulTriangle:
    """Verify the Virasoro line-side comparison representative."""

    def test_koszul_dual_c(self):
        """26 - c is the same-family comparison central charge."""
        assert koszul_dual_central_charge(c=0) == 26
        assert koszul_dual_central_charge(c=26) == 0
        assert koszul_dual_central_charge(c=13) == 13

    def test_line_comparison_fixed_c13(self):
        """The same-family Virasoro representative is fixed at c=13."""
        assert koszul_dual_central_charge(c=13) == 13

    def test_line_comparison_not_fixed_c26(self):
        """c=26 is not fixed by the same-family comparison."""
        assert koszul_dual_central_charge(c=26) != 26

    def test_double_dual(self):
        """The same-family central-charge involution squares to the identity."""
        for c_val in [0, 1, 13, 26, -5]:
            assert koszul_dual_central_charge(
                koszul_dual_central_charge(c=c_val)) == c_val

    def test_kappa_virasoro(self):
        """κ(Vir_c) = c/2."""
        assert gravity_kappa(c=0) == 0
        assert gravity_kappa(c=1) == Rational(1, 2)
        assert gravity_kappa(c=26) == 13

    def test_kappa_exact(self):
        """Exact rational κ."""
        assert gravity_kappa_exact(1) == Fraction(1, 2)
        assert gravity_kappa_exact(7, 10) == Fraction(7, 20)

    def test_complementarity_13(self):
        """κ(Vir_c) + κ(Vir_{26-c}) = 13 for all c."""
        assert complementarity_constant_virasoro() == 13
        for c_val in [0, 1, 5, 13, 26, -10]:
            result = verify_complementarity(c=c_val)
            assert result['equals_13'], f"Complementarity fails at c={c_val}"

    def test_complementarity_symbolic(self):
        """Symbolic complementarity verification."""
        result = verify_complementarity()
        assert result['equals_13']

    def test_verdier_modular_s_shadow_scope_is_not_partition_function(self):
        """The Verdier-S square is an algebraic shadow statement."""
        scope = verdier_modular_s_shadow_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma')
        assert scope['verdier_action'] == 'c -> 26 - c'
        assert scope['degree_zero_input'] == 'Theta^(1)_0(Vir_c) = (c/2) omega_1'
        assert scope['degree_zero_output'] == '((26-c)/2) S_sh(omega_1)'
        assert scope['commuting_square'] is True
        assert scope['physical_torus_partition_function_statement'] is False
        assert scope['dedekind_multiplier_required_for_partition_function'] is True
        assert scope['e2_anomaly_required_for_partition_function'] is True
        assert '0 <= c <= 26' in scope['cardy_expression_scope']
        assert 'Delta >= 0' in scope['cardy_expression_scope']
        assert 'Dedekind eta multiplier is trivial' in scope['not_asserted']
        assert 'E2 is a modular form without completion' in scope['not_asserted']

    def test_verdier_s_proposition_is_shadow_scoped(self):
        """The manuscript must not state D and S commute for the physical partition function."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Formal Verdier--$S$ square for the genus-one shadow;"
        )
        end = source.index(r"\end{proposition}", start)
        prop = source[start:end]
        post = source[end:source.index(r"\subsubsection*{Anomaly completion}", end)]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma$",
            r"\hypAmbientWtCpl",
            r"S_{\mathrm{sh}}",
            "Hodge line and the Dedekind multiplier line",
            r"\Theta_{\mathrm{Vir}_c}^{(1),\mathrm{sh}}",
            r"\frac{26-c}{2}\,S_{\mathrm{sh}}(\omega_1)",
            "not about",
            "physical torus partition function",
        )
        for needle in required:
            assert needle in prop

        assert r"0\le c\le 26" in post
        assert r"\Delta\ge0" in post
        assert "physical modular-invariance and vacuum-dominance hypotheses" in post
        assert r"2\pi\sqrt{\frac{c\Delta}{6}}" in post

        retired = (
            r"\begin{proposition}[Commutativity of $\mathbb{D}$ and $S$]",
            r"S \circ \mathbb{D}(\Theta_{\mathrm{Vir}_c}^{(1)}(\tau))",
            "commute on genus-$1$ MC elements",
            "The commutativity produces a paired Cardy formula",
            r"= 2\pi\sqrt{c\Delta/6} + 2\pi\sqrt{(26-c)\Delta/6}$",
        )
        block = prop + post
        for needle in retired:
            assert needle not in block

    def test_btz_complementarity_remark_has_physical_and_real_scope(self):
        """The later BTZ entropy copy must carry the same Cardy and real-domain scope."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{remark}[BTZ complementarity and the Hawking--Page transition]")
        end = source.index(r"\end{remark}", start)
        remark = source[start:end]

        required = (
            "modular invariance, vacuum",
            "saddle comparison hypotheses",
            r"\(0\le c\le26\)",
            r"shifted energy \(h\ge0\)",
            "paired scalar entropy",
            "same-family",
            r"2\pi\sqrt{h/6}\bigl(\sqrt{c} + \sqrt{26-c}\,\bigr)",
            r"maximized on \(0\le c\le26\)",
        )
        for needle in required:
            assert needle in remark

        retired = (
            "After the Cardy/BTZ entropy dictionary is imposed, the Koszul pair",
            "suggests a paired entropy",
        )
        for needle in retired:
            assert needle not in remark

    def test_gravity_koszul_triangle_is_central_projection_guard(self):
        """The manuscript must not identify C[[c]] with the whole bulk."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Central-sector gravitational Koszul triangle;"
        )
        end = source.index(r"\end{theorem}", start)
        theorem = source[start:end]

        required = (
            r"\ClaimStatusConditional",
            r"\hypBHdict+\hypAmbientWtCpl+\effKoszul",
            r"b_{\mathrm{BH}}",
            r"A_{b_{\mathrm{BH}}}=\mathrm{Vir}_c",
            r"c = 6k = 3\ell/(2G_N)",
            r"\widehat{\mathcal C}^{\,b_{\mathrm{BH}},\mathrm{pert}}_{\mathrm{line}}",
            r"\Zder^{\mathrm{ch}}(\mathrm{Vir}_c)",
            r"\pi_{\mathrm{cent}}",
            r"\HH^0_{\mathrm{GF}}(\mathrm{Vir})\oplus",
            r"\C[\![c]\!]",
            "it is not the full derived chiral",
        )
        for needle in required:
            assert needle in theorem

        retired = (
            r"\begin{theorem}[Gravitational Koszul triangle]",
            r"\textbf{Bulk central-extension projection} $\simeq \C[\![c]\!]$",
            "The three vertices of the bulk-boundary-line triangle are:",
            "from which three functors extract boundary, lines, and bulk",
            "No additional input beyond homotopy-Koszulity",
        )
        for needle in retired:
            assert needle not in theorem

    def test_steinberg_remark_keeps_bulk_out_of_bar_extraction(self):
        """The Steinberg analogy must keep the bulk as chiral Hochschild."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\begin{remark}[Steinberg analogy with central projection]")
        end = source.index(r"\end{remark}", start)
        remark = source[start:end]

        assert "bulk vertex is instead computed by chiral Hochschild cochains" in remark
        assert r"\xrightarrow{\pi_{\mathrm{cent}}}\C[\![c]\!]" in remark
        assert "universal central-extension projection" in remark
        assert "three functors extract boundary, lines, and bulk" not in remark

    def test_gravity_koszul_dual_branch_names_comparison_surface(self):
        """The manuscript must not collapse line duality to raw Verdier duality."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Same-family Virasoro representative of the line-side dual;"
        )
        end = source.index(r"\end{proposition}", start)
        prop = source[start:end]

        required = (
            r"\ClaimStatusConditional",
            r"\hypAmbientWtCpl+\effKoszul",
            r"A^{i,\mathrm{ord}}_{\mathrm{Vir}}",
            r"H^\bullet\!\left(\barB^{\mathrm{ord}}(\mathrm{Vir}_c)\right)^\vee",
            r"(\mathrm{Vir}_c)^{!_{\mathrm{line}},\infty}",
            r"\chi_{\mathrm{Vir}}^{\mathrm{line}}",
            r"\mathrm{Vir}_{26-c}",
            r"\kappaChHodge(\mathrm{Vir}_c) + \kappaChHodge(\mathrm{Vir}_{26-c})",
            "does not identify the Verdier dual",
            r"Conjecture~\textup{\ref{conj:gravity-line-identification}}",
        )
        for needle in required:
            assert needle in prop

        retired = (
            r"\begin{proposition}[Virasoro Koszul-dual branch]",
            r"The completed Koszul-dual branch of $\mathrm{Vir}_c$ is",
            r"$(\mathrm{Vir}_c)^!_\infty$",
            "represented by\n$\\mathrm{Vir}_{26-c}$",
        )
        for needle in retired:
            assert needle not in prop

    def test_brown_henneaux_line_test_package_scope(self):
        """The six-component package is a conditional test package, not the full open sector."""
        scope = brown_henneaux_line_test_package_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'delta', 'epsilon')
        assert scope['package'] == 'C_grav^test'
        assert scope['ambient'] == 'hypAmbientWtCpl'
        assert scope['boundary_algebra'] == 'A_{b_BH} = Vir_c'
        assert scope['bulk_object'] == 'Zder^ch(Vir_c)'
        assert scope['central_projection'] == 'C[[c]]'
        assert scope['central_projection_is_full_bulk'] is False
        assert scope['primitive_coproduct_scope'] == 'generator-level Virasoro T'
        assert scope['all_degree_primitivity_status'] == 'conditional on signed ghost-defect lemma'
        assert scope['line_model_status'] == 'Conjectural via gravity-line-identification'
        assert scope['annulus_trace_scope'] == 'vacuum one-boundary scalar character'
        assert scope['physical_torus_partition_function_statement'] is False
        assert scope['constructs_level_A_gravity_line_operator_algebra'] is False
        assert 'C[[c]] is the full derived chiral centre' in scope['not_asserted']
        assert 'vacuum character is the physical torus partition function' in scope['not_asserted']

    def test_grav_open_six_requirements_are_scoped(self):
        """The manuscript must keep the Brown-Henneaux line package conditional and projected."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Brown--Henneaux algebraic line test package;"
        )
        end = source.index(r"\end{proposition}", start)
        prop = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma+\delta+\varepsilon$",
            r"\hypBHdict+\hypAmbientWtCpl+\effKoszul",
            r"A_{b_{\mathrm{BH}}}=\End(b_{\mathrm{BH}})\simeq\mathrm{Vir}_c",
            r"\Zder^{\mathrm{ch}}(\mathrm{Vir}_c)",
            "only its central-extension projection",
            r"\Delta_z^{\mathrm{Vir}}(T)",
            "signed ghost-defect lemma",
            "is not proved by the",
            r"Conjecture~\textup{\ref{conj:gravity-line-identification}}",
            "not the full physical torus partition function",
            r"F_g(\mathrm{Vir}_c)=\kappaChHodge(\mathrm{Vir}_c)",
            "not a construction of the level-",
            r"\Phi_{10}^{\mathrm{un}}",
        )
        for needle in required:
            assert needle in prop

        assert "not the full derived chiral centre" in block
        assert "metric path integral" in block
        assert "BTZ saddle sum" in block

        retired = (
            "The open-sector factorization dg-category of 3d gravity is",
            r"\begin{proposition}[Gravitational primitive package]",
            r"Definition~\textup{\ref{def:grav-open-algebraic}} satisfies",
            "Primitive coproduct by construction",
            r"\cC_{\mathrm{line}} \simeq",
            "the genus-one closed-string partition function restricted",
            "the tensor product is trivial",
            r"\Delta_z(x) = \tau_z(x) \otimes 1 + 1 \otimes x",
        )
        for needle in retired:
            assert needle not in block

    def test_virasoro_bar_intrinsic_mc_shadow_scope(self):
        """The Virasoro MC package is the positive-genus bar correction."""
        scope = virasoro_bar_intrinsic_mc_shadow_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'epsilon')
        assert 'hypBHdict' in scope['requires']
        assert 'hypAmbientWtCpl' in scope['requires']
        assert 'effKoszul' in scope['requires']
        assert scope['mc_element'] == 'Theta_grav = Theta_Vir_c = D_Vir_c - d0'
        assert scope['genus_zero_operations_in_theta'] is False
        assert 'Stasheff' in scope['genus_zero_role']
        assert scope['f1'] == 'F_1(Vir_c) = c/48'
        assert 'lambda_g^FP' in scope['uniform_weight_scalar_lane']
        assert scope['non_scalar_higher_genus_component'] == 'stable-graph coderivation d_Vir_c^(g)'
        assert 'Q_contact_Vir=10/(c(5c+22))' in scope['shadow_projections']
        assert scope['central_projection_is_full_bulk'] is False
        assert scope['physical_partition_function_statement'] is False
        assert 'Theta_grav includes sum_{k>=2} alpha_k' in scope['not_asserted']
        assert 'C[[c]] is equivalent to the full derived chiral centre' in scope['not_asserted']

    def test_genus0_directed_product_scope(self):
        """The genus-zero product is directed, not a disjoint-colour product."""
        scope = genus0_directed_product_decomposition_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('alpha', 'gamma')
        assert 'hypAmbientWtCpl' in scope['requires']
        assert scope['product_kind'] == 'directed colour-filtered product'
        assert scope['ordinary_disjoint_product'] is False
        assert scope['color_order'] == ('ch', 'bdy', 'tr')
        assert scope['operation_spaces']['tr_output'] == 'FM_k(C) x E1(m) x E1(p)'
        assert scope['mixed_inputs_to_tr_present'] is True
        assert scope['strict_after_chart_choice'] is True
        assert 'module category over the SCchtop algebra' in scope['algebra_scope']
        assert 'ordinary product of two coloured operads on disjoint colour sets' in scope['not_asserted']

    def test_genus0_product_decomposition_is_not_external_product_in_manuscript(self):
        """The proposition must name the directed product and retain mixed tr inputs."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Directed genus-$0$ product decomposition;"
        )
        end = source.index(r"\end{proposition}", start)
        prop = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\alpha+\gamma$",
            r"\hypAmbientWtCpl",
            r"\SCchtop \;\boxtimes_{\mathrm{dir}}\; \Eone^{\mathsf{tr}}",
            r"\mathsf{ch}\leq\mathsf{bdy}\leq\mathsf{tr}",
            r"\FM_k(\C)\times\Eone(m)\times\Eone(p)",
            "mixed",
            r"\mathsf{ch}$/$\mathsf{bdy}",
            "ordinary product of",
            "disjoint colour sets",
        )
        for needle in required:
            assert needle in prop

        assert "strict directed product" in block
        assert "module category over the" in source

        retired = (
            r"\begin{proposition}[Product decomposition at genus~$0$]",
            "isomorphic as a three-coloured topological operad to the\nexternal product",
            r"\SCchtop \;\times\; \Eone^{\mathsf{tr}}",
            r"Therefore the operad structure is the external product",
            r"underlying object carries an $\SCchtop$-algebra",
        )
        for needle in retired:
            assert needle not in block

    def test_qt_equivariance_scope_names_all_annular_data(self):
        """Graph equivariance needs braid descent and inverse orientation."""
        scope = quantum_group_clutching_equivariance_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma')
        assert 'hypAmbientWtCpl' in scope['requires']
        assert 'pure-braid descent to Aut(Gamma)' in scope['requires']
        assert 'inverse-orientation relation R_21(-z) R_12(z) = id' in scope['requires']
        assert scope['endpoint_exchange_mechanism'] == 'quasi-triangular coproduct comparison'
        assert 'positive braid lifts' in scope['edge_permutation_mechanism']
        assert scope['literal_commutativity_of_monodromy_products'] is False
        assert scope['inverse_orientation_from_quasitriangularity'] is False
        assert scope['braiding_is_symmetric_group_action_before_descent'] is False
        assert 'quasi-triangularity plus YBE alone' in scope['not_asserted'][0]

    def test_qt_equivariance_manuscript_rejects_false_qt_ybe_shortcut(self):
        """The proposition must not derive symmetric graph equivariance from QT+YBE alone."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Braided annular clutching is graph-equivariant]"
        )
        end = source.index(r"\end{proposition}", start)
        prop = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]

        required_prop = (
            r"\ClaimStatusProvedHere",
            r"licensing $\alpha+\beta+\gamma$",
            r"\hypAmbientWtCpl",
            "pure-braid",
            r"R_{21}(-z)R_{12}(z)=\id",
            "not a consequence of",
            "quasi-triangularity alone",
        )
        for needle in required_prop:
            assert needle in prop

        required_proof = (
            "No commutativity of different monodromy",
            "Matsumoto",
            "This is braid coherence",
            "not the ordinary equality",
            "setting $\\Delta=\\id$",
        )
        for needle in required_proof:
            assert needle in block

        retired_block = (
            r"\begin{proposition}[Quasi-triangularity implies clutching equivariance]",
            r"\prod_{i=1}^{k} \mathrm{Mon}(R)_{e_{\tau(i)}}",
            "as operators on the sewing, because the monodromy insertions act",
            r"quasi-triangularity with $\Delta = \mathrm{id}$",
        )
        for needle in retired_block:
            assert needle not in block

        retired_source = (
            "proved unconditionally for $E_1$-chiral quantum groups",
            "quasi-triangularity + YBE pair is \\emph{sufficient}",
        )
        for needle in retired_source:
            assert needle not in source

    def test_qt_equivariance_scope_propagates_to_active_summaries(self):
        """Introduction and modular-operad summary must carry the repaired scope."""
        root = Path(__file__).resolve().parents[2]
        intro = (root / "chapters" / "theory" / "introduction.tex").read_text()
        modop = (
            root / "chapters" / "theory" / "modular_swiss_cheese_operad.tex"
        ).read_text()
        standalone = (root / "standalone" / "stokes_gap_kzb_regularity.tex").read_text()

        for text in (intro, modop, standalone):
            assert "pure-braid descent" in text
            assert "inverse" in text

        assert "braided-annular equivariance" in intro
        assert "braided annular datum" in modop
        assert "Four are proved unconditionally" not in standalone
        assert "quasi-triangularity plus\nYang--Baxter" not in standalone

    def test_modular_operad_unitality_scope_is_counit_normalised(self):
        """Unitality uses the diagonal bicomodule, not a disk monodromy slogan."""
        scope = modular_operad_unitality_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma')
        assert 'unit-extended stable-graph chart' in scope['requires']
        assert 'counit-normalised annular sewing' in scope['requires']
        assert 'diagonal bicomodule C_Delta' in scope['requires']
        assert 'hypAmbientWtCpl' in scope['requires']
        assert scope['exceptional_vertex_stable'] is False
        assert scope['exceptional_vertex_kind'] == 'unstable genus-0 two-flag identity vertex'
        assert scope['annular_unit_seam'] == 'C_Delta identity bicomodule after unit/counit'
        assert 'square_C C_Delta' in scope['cotensor_identity']
        assert scope['uses_contractible_disk_argument'] is False
        assert scope['proves_general_clutching_existence'] is False
        assert scope['proves_composition_associativity'] is False
        assert scope['independent_of_genus_and_shadow_class_after_normalisation'] is True

    def test_modular_operad_unitality_manuscript_uses_unit_extended_graphs(self):
        """The unitality proposition must not use the simply-connected unit argument."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Unit-normalised annular sewing is unital at every genus]"
        )
        end = source.index(r"\end{proposition}", start)
        prop = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]

        required_prop = (
            r"\ClaimStatusProvedHere",
            r"licensing $\alpha+\beta+\gamma$",
            r"\hypAmbientWtCpl",
            "unit-extended stable-graph chart",
            "counit-normalised annular sewing",
            "unstable genus-$0$ two-flag identity vertex",
            r"(\varepsilon\otimes\id)\Delta_z=\id",
            r"\mathrm{Mon}(R)_{\eta,e}=\id",
            "does not assert the existence",
        )
        for needle in required_prop:
            assert needle in prop

        required_proof = (
            "not an additional stable vertex",
            "diagonal bicomodule is the",
            r"M \,\square_C\, C_\Delta",
            "not the assertion that a\nunit component bounds a simply connected disk",
            "strictly unital",
            r"\xi_{\Gamma^\eta}=\xi_{\Gamma/e}",
        )
        for needle in required_proof:
            assert needle in block

        retired = (
            r"\begin{proposition}[Unitality of $\mathcal{O}^{\Ainf\text{-ch}}$",
            "unit vertex~$u$ is a punctured sphere",
            "simply\nconnected",
            "parallel transport of any flat\nconnection around a contractible loop",
            "monodromy around a\ncontractible cycle",
            "has insufficient\ninputs",
        )
        for needle in retired:
            assert needle not in block

    def test_modular_operad_unitality_scope_propagates_to_active_summaries(self):
        """Summaries must name unit/counit normalisation, not only a vacuum slogan."""
        root = Path(__file__).resolve().parents[2]
        intro = (root / "chapters" / "theory" / "introduction.tex").read_text()
        modop = (
            root / "chapters" / "theory" / "modular_swiss_cheese_operad.tex"
        ).read_text()
        standalone = (root / "standalone" / "stokes_gap_kzb_regularity.tex").read_text()

        assert "unit-normalised annular unitality" in intro
        assert "counit-normalised diagonal bicomodule" in modop
        assert "unit/counit normalisation of the annular seam" in modop
        assert "unit-normalised annular sewing" in standalone
        assert "unitality (vacuum axiom)" not in standalone

    def test_modular_bar_reduction_scope_requires_square_zero_internal_data(self):
        """The uncurved modular-bar theorem cannot consume raw curvature."""
        scope = modular_bar_reduction_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'epsilon')
        assert 'hypAmbientWtCpl' in scope['requires']
        assert 'square-zero internal differential or S-tail curvature datum' in scope['requires']
        assert scope['abstract_theorem'] == 'thm:modular-bar'
        assert scope['abstract_theorem_requires_internal_square_zero'] is True
        assert scope['positive_genus_curvature_is_modular_bar_datum_by_itself'] is False
        assert scope['curvature_term'] == 'kappaChHodge(A) omega_g'
        assert 'd_S^2=0' in scope['curvature_absorption']
        assert scope['one_edge_maps_degree'] == 1
        assert 'genus 0 affine non-critical' in scope['covered_loci']
        assert 'general A_infinity-chiral genus >= 1 one-edge construction' in scope['open_loci']
        assert scope['proves_general_positive_genus_clutching_maps'] is False
        assert scope['proves_concrete_operadic_associativity_by_itself'] is False

    def test_modular_bar_reduction_manuscript_is_a_criterion_not_a_construction(self):
        """The proposition must not feed curved d^2 directly into thm:modular-bar."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Modular-bar criterion for annular clutching data]"
        )
        end = source.index(r"\end{proposition}", start)
        prop = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]

        required_prop = (
            r"\ClaimStatusProvedHere",
            r"licensing",
            r"\alpha+\beta+\gamma+\varepsilon",
            r"\hypAmbientWtCpl",
            r"\(S\)-tail curvature datum",
            "Internal square-zero datum",
            r"d_S:=d+D_S",
            "curvature data,\n not yet a modular bar datum",
            r"d_S\circ\Delta_\epsilon+\Delta_\epsilon\circ d_S=0",
            "does not by itself construct",
            "or prove concrete operadic associativity",
        )
        for needle in required_prop:
            assert needle in prop

        required_proof = (
            "applies only to a\nmodular bar datum",
            r"(D_{\mathrm{int}}+D_{\mathrm{exp}})^2",
            "requires a square-zero internal differential",
            "converted into an uncurved modular-bar input",
            "not a consequence of the abstract\ntheorem",
            "if the positive-genus curvature has not\nbeen absorbed",
        )
        for needle in required_proof:
            assert needle in block

        retired = (
            r"\begin{proposition}[Reduction to modular bar datum axioms]",
            "The four axioms decompose as follows",
            "The vertex bar differential satisfies $d^2 = 0$\n at genus~$0$ and",
            "The modular bar differential $D^2 = 0$\n\\textup{(}Theorem~\\textup{\\ref{thm:modular-bar}}",
            "Parts~(a) and~(c) are immediate",
            "Proposition~\\ref{prop:ainf-chiral-modular-bar-reduction}\nverifies all four axioms",
        )
        for needle in retired:
            assert needle not in block

    def test_affine_modular_bar_datum_scope_is_exact_kz_kzb_locus(self):
        """The affine corollary is not an all-genera generic-level theorem."""
        scope = affine_modular_bar_datum_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'epsilon')
        assert 'hypAmbientWtCpl' in scope['requires']
        assert scope['criterion_used'] == 'prop:ainf-chiral-modular-bar-reduction'
        assert scope['abstract_theorem'] == 'thm:modular-bar'
        assert 'genus 0 at every non-critical level' in scope['covered_loci']
        assert any('integrable level' in locus for locus in scope['covered_loci'])
        assert 'KZ logarithmic regular singular' in scope['genus_zero_mechanism']
        assert 'regular-singular KZB extension' in scope['integrable_all_genera_mechanism']
        assert 'curvature compensation' in scope['integrable_all_genera_mechanism']
        assert scope['positive_genus_generic_nonintegral_claimed'] is False
        assert scope['stokes_gap_closed'] is False
        assert scope['proves_modular_operad_axioms_by_itself'] is False

    def test_affine_modular_bar_datum_manuscript_excludes_generic_positive_genus(self):
        """The corollary must name the covered loci and the Stokes exclusion."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{corollary}[Affine modular-bar datum on the KZ/KZB covered loci]"
        )
        end = source.index(r"\end{corollary}", start)
        cor = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]
        remark_end = source.index(r"\end{remark}", proof_end)
        remark = source[proof_end:remark_end]

        required_corollary = (
            r"\ClaimStatusProvedHere",
            r"\alpha+\beta+\gamma+\varepsilon",
            r"\hypAmbientWtCpl",
            r"genus~\(0\), at every non-critical level",
            r"arbitrary genus, at integrable level",
            "regular-singular KZB extension",
            "curvature\n compensation datum",
            "asserts no modular bar datum",
        )
        for needle in required_corollary:
            assert needle in cor

        required_proof = (
            "vertex bar differential is uncurved",
            "Drinfeld pentagon",
            "hexagon",
            "semisimple integrable\npositive-energy module category",
            r"\(d_S^2=0\)",
            "without Stokes sector ambiguity",
            "the criterion therefore does not\napply",
        )
        for needle in required_proof:
            assert needle in block

        assert "The corollary alone supplies only the\nmodular bar datum" in remark

        retired = (
            r"\begin{corollary}[Modular bar datum for affine chiral",
            "with KZ\\slash KZB monodromy forms a modular bar datum",
            "At genus~$0$, the same holds at every non-critical level",
            "single analytic question: does the KZB connection at\ngeneric level have regular singularities",
        )
        for needle in retired:
            assert needle not in block

    def test_affine_e3_topological_km_scope_is_cohomological(self):
        """The affine E3 theorem is cohomological, not raw-chain-level."""
        scope = affine_e3_topological_km_scope()
        assert scope['claim_status'] == 'ProvedHere'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'epsilon')
        assert 'non-critical level k != -h^vee' in scope['requires']
        assert 'Q_CS-cohomology ambient' in scope['requires']
        assert scope['sugawara_denominator'] == '2(k + h^vee)'
        assert scope['critical_level_excluded'] is True
        assert scope['brst_identity_scope'] == 'on Q_CS-cohomology'
        assert scope['topological_structure_target'] == 'H^bullet_{Q_CS} Zder^ch(V_k(g))'
        assert scope['applies_dunn_to_bicoloured_scchtop'] is False
        assert scope['strict_raw_chain_level_statement'] is False
        assert scope['strict_raw_chain_level_frontier'] == 'rem:frontier-class-L-strict-chain-level'

    def test_affine_e3_topological_km_manuscript_names_cohomology_target(self):
        """The theorem must not claim raw-chain-level E3 topologisation."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Cohomological $\Ethree$-topological structure for affine"
        )
        end = source.index(r"\end{theorem}", start)
        theorem = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]

        required_statement = (
            r"\ClaimStatusProvedHere",
            r"\alpha+\beta+\gamma+\varepsilon",
            r"\hypKMHTBV",
            r"k \ne -h^\vee",
            r"H^\bullet_{Q_{\mathrm{CS}}}\Zder^{\mathrm{ch}}(V_k(\fg))",
            r"\(Q_{\mathrm{CS}}\)-cohomological topologisation statement",
            "does\nnot assert a strict raw-chain-level",
        )
        for needle in required_statement:
            assert needle in theorem

        required_proof = (
            r"\quad\text{on }Q_{\mathrm{CS}}\text{-cohomology}",
            r"T_{\mathrm{Sug},(0)}=L_{-1}",
            "holomorphic translations act trivially in cohomology",
            "single-colour external\nproduct",
            "Dunn additivity gives",
            r"rem:frontier-class-L-strict-chain-level",
        )
        for needle in required_proof:
            assert needle in block

        retired = (
            r"\begin{theorem}[$\Ethree$-topological for affine Kac--Moody]",
            "The derived chiral center\n$\\Zder^{\\mathrm{ch}}(V_k(\\fg))$ carries",
            "The theorem establishes\n$T_{\\mathrm{Sug}} = [Q_{\\mathrm{CS}}, G_{\\mathrm{Sug}}]$ on\n$H^\\bullet",
        )
        for needle in retired:
            assert needle not in block

    def test_principal_ds_e3_topological_scope_is_conditional(self):
        """The DS E3 theorem is over Q_DS and the DS primitive package."""
        scope = principal_ds_e3_topological_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'epsilon')
        assert scope['hypothesis_package'] == 'hypDSBRST'
        assert 'principal nilpotent f_prin' in scope['requires']
        assert 'total DS differential Q_DS = Q_CS + Q_red' in scope['requires']
        assert scope['brst_identity_scope'] == 'on Q_DS-cohomology'
        assert scope['topological_structure_target'] == (
            'H^bullet_{Q_DS} Zder^ch(W^k(g,f_prin))'
        )
        assert scope['applies_dunn_to_bicoloured_scchtop'] is False
        assert scope['strict_raw_chain_level_statement'] is False
        assert scope['unreduced_cartan_current_qcs_exactness'] is False
        assert scope['normalisation_dependent_cartan_antighost_coefficient'] is True

    def test_principal_ds_e3_manuscript_uses_qds_not_unreduced_qcs(self):
        """The DS theorem must name the DS differential and exclude old overclaims."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Principal DS cohomological $\Ethree$-topologisation"
        )
        end = source.index(r"\end{theorem}", start)
        theorem = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]
        intro = (root / "chapters" / "theory" / "introduction.tex").read_text()

        required_statement = (
            r"\ClaimStatusConditional",
            r"\alpha+\beta+\gamma+\varepsilon",
            r"\hypDSBRST",
            r"Q_{\mathrm{DS}}=Q_{\mathrm{CS}}+Q_{\mathrm{red}}",
            r"H_{Q_{\mathrm{DS}}}^\bullet\Zder^{\mathrm{ch}}(\cW)",
            r"does not assert that Cartan currents are",
            r"does not assert a strict raw-chain-level",
        )
        for needle in required_statement:
            assert needle in theorem

        required_proof = (
            r"T_{\mathrm{DS},(0)}=L_{-1}",
            r"Q_{\mathrm{DS}}\text{-cohomology}",
            "single-colour external product",
            r"E_2^{\mathrm{top}}\otimes E_1^{\mathrm{top}}",
            r"does not identify affine Cartan currents with \(Q_{\mathrm{CS}}\)-",
        )
        for needle in required_proof:
            assert needle in block

        assert r"\hypDSBRST" in intro
        assert r"T_{\mathrm{DS}}=[Q_{\mathrm{DS}},G'_{\mathrm{DS}}]" in intro
        assert (
            "not the false statement that Cartan currents are\n"
            "\\(Q_{\\mathrm{CS}}\\)-exact"
        ) in intro

        retired = (
            r"\begin{theorem}[Cohomological $\Ethree$-topological structure via DS reduction]",
            r"T_{\mathrm{DS}} \;=\; [Q_{\mathrm{CS}},\, G']",
            "The key observation is that each Cartan current",
            r"\label{eq:current-Q-exact}",
            r"\label{eq:G-imp-formula}",
            r"\tfrac{1}{4}\,\partial\bar c_h",
        )
        for needle in retired:
            assert needle not in block

    def test_gravity_mc_package_is_positive_genus_and_bulk_projected(self):
        """The manuscript must not fold genus zero into Theta_grav or collapse the bulk."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        bulk_start = source.index(r"\subsubsection*{Bulk (derived center)}")
        bulk_end = source.index(r"\subsubsection*{Classical $r$-matrix and coproduct}", bulk_start)
        bulk = source[bulk_start:bulk_end]
        start = source.index(
            r"\begin{theorem}[Bar-intrinsic Virasoro MC shadow;"
        )
        end = source.index(r"\end{theorem}", start)
        theorem = source[start:end]
        proof_end = source.index(r"\end{proof}", end)
        block = source[start:proof_end]

        bulk_required = (
            r"\pi_{\mathrm{cent}}\colon",
            r"\Zder^{\mathrm{ch}}(\mathrm{Vir}_c)",
            r"\longrightarrow",
            r"\C[\![c]\!]",
            "not an equivalence with the full derived chiral",
            "central projection",
        )
        for needle in bulk_required:
            assert needle in bulk

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma+\varepsilon$",
            r"\hypBHdict+\hypAmbientWtCpl+\effKoszul",
            r"D_{\mathrm{Vir}_c}",
            r"\Theta_{\mathrm{grav}}",
            r"\Theta_{\mathrm{Vir}_c}",
            r"D_{\mathrm{Vir}_c}-d_0",
            r"\sum_{g\ge1}\hbar^g d_{\mathrm{Vir}_c}^{(g)}",
            r"[d_0,\Theta_{\mathrm{grav}}]",
            "not part of the positive-genus",
            "uniform-weight scalar lane",
            "full non-scalar genus component",
            r"\mathfrak{Q}^{\mathrm{contact}}_{\mathrm{Vir}}",
            "not replacements for the full MC element",
        )
        for needle in required:
            assert needle in theorem

        assert r"\Defcyc(\cA)\widehat{\otimes}\Gmod" in block
        assert "physical torus partition function" not in block

        retired = (
            r"\begin{theorem}[Gravitational MC element: primitive package]",
            "The full modular MC element of $\\mathrm{Vir}_c$ is",
            r"\sum_{k \ge 2} \alpha_k",
            r"\sum_{g \ge 1}\,",
            r"\Theta_{\mathrm{grav}} = D_A - d_0",
            "guarantees that $\\Theta_{\\mathrm{grav}}$ satisfies the MC",
            r"$\kappa = c/2$ \textup{(degree~$2$)}",
            r"\mathcal{Z}^{\mathrm{der}}_{\mathrm{ch}}(\mathrm{Vir}_c)"
            "\n\\;\\simeq\\;",
        )
        joined = bulk + block
        for needle in retired:
            assert needle not in joined


def test_e3_scope_map_keeps_logarithmic_wp_frontier():
    """Logarithmic triplets are not inherited E3-topological rows."""
    root = Path(__file__).resolve().parents[2]
    source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
    flat = " ".join(source.split())

    required = (
        "non-unitary logarithmic triplet algebras $\\cW(p)$ at $p\\ge2$ "
        "are not absorbed into case~\\ref{item:scope-orbifold}",
        "requires a separate logarithmic HT BV datum",
        "regular/logarithmic amplitude estimates of "
        "Conjecture~\\ref{conj:tempered-stratum-contains-wp}",
        "frontier row rather than an inherited $\\Ethree$-topological theorem",
        "The logarithmic triplet family $\\cW(p)$ is excluded from those two "
        "proved lanes at the current scope",
    )
    for needle in required:
        assert needle in flat

    retired = (
        "No example in the standard landscape is currently known to fall here",
        "If no member of the standard landscape populates case~\\ref{item:scope-frontier}",
        "Conjecture~\\ref{conj:E3-topological-general} is vacuously a theorem",
        "logarithmic triplet algebras $\\cW(p)$ lie in case~\\ref{item:scope-orbifold}",
    )
    for needle in retired:
        assert needle not in flat


# ===================================================================
# 7. R-MATRIX
# ===================================================================

class TestGravitationalRMatrix:
    """Verify r^L(z) and the bar collision r-matrix are distinct."""

    def test_laplace_kernel_leading_pole_order_4(self):
        """The pre-dlog Laplace/OPE kernel has leading order 4."""
        poles = gravity_laplace_kernel_poles(c=1)
        assert 4 in poles

    def test_collision_leading_pole_order_3(self):
        """The bar collision r-matrix has leading order 3."""
        poles = gravity_r_matrix_poles(c=1)
        assert 3 in poles
        assert 4 not in poles

    def test_leading_residue_c_over_2(self):
        """Leading collision residue is c/2."""
        assert gravity_r_matrix_leading_residue(c=26) == 13
        assert gravity_r_matrix_leading_residue(c=1) == Rational(1, 2)

    def test_laplace_transform_scalar(self):
        """Laplace: ∫₀^∞ e^{-λz}(c/12)λ³ dλ = (c/12)(3!/z⁴) = c/(2z⁴)."""
        c_val = Rational(24)
        # (c/12) * 3! = (c/12) * 6 = c/2
        laplace_result = c_val / 12 * 6
        expected = c_val / 2
        assert laplace_result == expected

    def test_laplace_pole_structure_complete(self):
        """r^L(z) has exactly poles of order 4, 2, 1."""
        poles = gravity_laplace_kernel_poles(c=13)
        assert set(poles.keys()) == {4, 2, 1}

    def test_collision_pole_structure_complete(self):
        """r^{coll}(z) has exactly poles of order 3 and 1."""
        poles = gravity_r_matrix_poles(c=13)
        assert set(poles.keys()) == {3, 1}

    def test_collision_r_matrix_absorbs_one_pole(self):
        """Collision pole orders are OPE pole orders shifted down by one."""
        bracket = virasoro_lambda_bracket(c=12)
        laplace = gravity_laplace_kernel_poles(c=12)
        collision = gravity_r_matrix_poles(c=12)

        assert bracket['lam3'] == 1
        assert 4 in laplace and 2 in laplace and 1 in laplace
        assert 3 in collision and 1 in collision
        assert 4 not in collision
        assert 2 not in collision

    def test_collision_cybe_scope_profile(self):
        """The gravitational CYBE is a collision-residue MC statement."""
        scope = virasoro_collision_cybe_scope()
        assert scope['kernel'] == 'r_coll_Vir(z) = (c/2)/z^3 + 2T/z'
        assert scope['pre_comparison_form'] == 'arity-three MC equation for Theta_Vir_c'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma')
        assert 'Arnold relation on FM_3(C)' in scope['proof_mechanism']
        assert scope['is_finite_dimensional_casimir_cybe'] is False
        assert scope['is_strict_quantum_ybe_for_fusion_kernel'] is False

    def test_gravity_ybe_proposition_is_collision_residue_scoped(self):
        """The manuscript must not state an unqualified Virasoro YBE."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Collision-residue gravitational CYBE;"
        )
        end = source.index(r"\end{proposition}", start)
        prop = source[start:end]
        flat = " ".join(prop.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma$",
            r"\hypBHdict+\hypAmbientWtCpl",
            r"r^{\mathrm{coll}}_{\mathrm{Vir}}",
            r"\Rescoll_{0,2}(\Theta_{\mathrm{Vir}_c})",
            r"\frac{c/2}{z^3}+\frac{2T}{z}",
            r"Theorem~\ref{thm:thqg-V-cybe-from-arnold}",
        )
        for needle in required:
            assert needle in prop

        assert "ordinary finite-dimensional Casimir identity" in flat
        assert "not the strict quantum YBE for the fusion kernel" in flat

        retired = (
            r"\begin{proposition}[Gravitational YBE]",
            "The CYBE",
            r"is the residue of the Arnold relation on $\FM_3(\C)$",
        )
        for needle in retired:
            assert needle not in prop


# ===================================================================
# 8. PRIMARY-STATE WARD FACTOR
# ===================================================================

class TestVirasoroPrimaryWardFactor:
    """Verify the scalar primary-state Ward factor is only a projection."""

    def test_connection_is_log_derivative(self):
        """d(2h log z - c/(4z^2)) = 2h/z + c/(2z^3)."""
        c, h, z = symbols('c h z')
        assert simplify(
            diff(virasoro_primary_ward_log(c=c, h=h, z=z), z)
            - virasoro_primary_ward_connection(c=c, h=h, z=z)
        ) == 0

    def test_connection_uses_collision_c_over_2_not_lambda_c_over_12(self):
        """The z^-3 central coefficient is c/2 after collision extraction."""
        z = Symbol('z')
        connection = virasoro_primary_ward_connection(c=24, h=0, z=z)
        assert simplify(connection - 12 / z**3) == 0

    def test_first_four_even_coefficients(self):
        """exp(-c/(4z^2)) has coefficients 1, -c/4, c^2/32, -c^3/384."""
        c = Symbol('c')
        coeffs = virasoro_primary_ward_even_coefficients(c=c, max_k=3)
        assert coeffs[0] == 1
        assert coeffs[1] == -c / 4
        assert coeffs[2] == c**2 / 32
        assert coeffs[3] == -c**3 / 384

    def test_no_odd_inverse_powers_after_primary_factor(self):
        """After removing z^(2h), every exponent is even."""
        exponents = virasoro_primary_ward_even_exponents(max_k=8)
        assert exponents == [0, -2, -4, -6, -8, -10, -12, -14, -16]
        assert all(exp % 2 == 0 for exp in exponents)


# ===================================================================
# 9. GENUS EXPANSION
# ===================================================================

class TestGenusExpansion:
    """Verify F_g via Â-genus formula."""

    def test_lambda_fp_1(self):
        """λ₁^FP = 1/24."""
        result = genus_generating_function_coefficients(c=2, max_genus=1)
        assert result['lambda_fp'][1] == Rational(1, 24)

    def test_lambda_fp_2(self):
        """λ₂^FP = 7/5760."""
        result = genus_generating_function_coefficients(c=2, max_genus=2)
        assert result['lambda_fp'][2] == Rational(7, 5760)

    def test_lambda_fp_3(self):
        """λ₃^FP = 31/967680."""
        result = genus_generating_function_coefficients(c=2, max_genus=3)
        assert result['lambda_fp'][3] == Rational(31, 967680)

    def test_F1_raw(self):
        """F₁ = κ * λ₁^FP = (c/2) * (1/24) = c/48."""
        result = genus_generating_function_coefficients(c=48, max_genus=1)
        assert result['raw'][1] == 1  # 48/2 * 1/24 = 1

    def test_F1_effective(self):
        """F₁^eff = κ_eff * λ₁^FP = (c-26)/2 * 1/24."""
        result = genus_generating_function_coefficients(c=26, max_genus=1)
        assert result['effective'][1] == 0  # (26-26)/2 * 1/24 = 0

    def test_ahat_coefficients_match_series(self):
        """Â-genus coefficients match direct series expansion of (x/2)/sinh(x/2)."""
        result = verify_ahat_series(max_genus=4)
        assert result['all_match'], \
            f"Â-genus mismatch: {result['by_genus']}"

    def test_lambda_fp_positive(self):
        """All λ_g^FP > 0 (Bernoulli sign analysis)."""
        result = genus_generating_function_coefficients(c=2, max_genus=5)
        for g in range(1, 6):
            assert result['lambda_fp'][g] > 0, f"λ_{g}^FP not positive"


class TestGenusOneVirasoroMCScope:
    """Guard the genus-1 Virasoro MC scalar/Ward normalization split."""

    def test_scope_separates_hodge_trace_from_vacuum_casimir(self):
        """At c=24, κ=12, F₁=1/2, and the Ward constant is -1."""
        scope = genus1_virasoro_mc_scope(c=24)
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'hypVirTorusBlock')
        assert scope['kappa'] == 12
        assert scope['hodge_trace'] == Rational(1, 2)
        assert scope['vacuum_one_point_constant'] == -1
        assert scope['two_point_central_kernel'] == '(c/2) P4(z,tau)'
        assert scope['P2_definition'] == 'P2=wp+E2/12'
        assert scope['P4_definition'] == 'P4=(1/6) d_z^2 P2'
        assert scope['full_all_n_stable_graph_series_evaluated'] is False
        assert scope['regular_torus_blocks_determined_by_kappa'] is False
        assert scope['central_kernel_is_wp'] is False
        assert scope['casimir_equals_scalar_hodge_trace'] is False

    def test_genus1_mc_theorem_is_conditionally_scoped(self):
        """The theorem must not sell scalar data as full torus blocks."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Genus-$1$ Virasoro MC scalar and Ward package;"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma+\delta$",
            r"\hypAmbientWtCpl+\hypVirTorusBlock",
            r"\Theta^{(1)}_{n=0}",
            r"\kappaChHodge(\mathrm{Vir}_c)",
            r"\frac{c}{2}\,\omega_1",
            r"T_{(3)}T=c/2",
            r"\operatorname{CT}_{q}",
            r"Z^{\mathrm{Vir}}_1",
            r"q^{-c/24}\prod_{m\ge2}(1-q^m)^{-1}",
            r"\frac{c}{2}\,P_4(z_{12},\tau)",
            r"P_4(z,\tau)=\frac16\partial_z^2P_2(z,\tau)",
            "regular torus blocks and convergence of the full",
            "regular part is a torus block depending on the chosen module",
            "module-dependent torus conformal block",
        )
        for needle in required:
            assert needle in block

        retired = (
            r"\frac{c}{2}\,\wp(z_{12};\tau)",
            r"\sum_{n=1}^{\infty}",
            r"\frac{E_2(\tau)}{24}\cdot T",
            "genus-$1$\nfree energy density",
            "where $p(n)$ counts partitions",
            r"E_2 \cdot T \otimes T",
            "The full all-$n$ stable graph series is evaluated by",
            "The regular part is determined by",
        )
        for needle in retired:
            assert needle not in block


class TestGenusOneVirasoroAmplitudesScope:
    """Guard the genus-1 Virasoro Ward-amplitude theorem."""

    def test_scope_includes_p1_derivative_and_regular_block_hypothesis(self):
        """At c=24, det-shadow coefficient is -1/2 and Ward constant is -1."""
        scope = genus1_virasoro_amplitudes_scope(c=24)
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'hypVirTorusBlock')
        assert scope['scalar_hodge_trace'] == Rational(1, 2)
        assert scope['determinant_line_derivative_coefficient'] == Rational(-1, 2)
        assert scope['vacuum_one_point_constant'] == -1
        assert scope['generic_vacuum_character_product_starts_at'] == 2
        assert scope['P1_local_expansion'] == 'P1=z^(-1)+O(z)'
        assert scope['two_point_simple_pole_derivative_vanishes'] is True
        assert scope['two_point_regular_block_from_hypothesis'] is True
        assert scope['three_point_requires_p1_derivative'] is True
        assert 'P1(z_ij,tau) partial_z_j A2' in scope['three_point_pairwise_singular_terms']
        assert scope['regular_blocks_determined_by_kappa'] is False
        assert scope['minimal_model_quotients_included'] is False
        assert scope['physical_torus_partition_function_statement'] is False

    def test_genus1_amplitudes_theorem_is_conditional_and_has_derivative_term(self):
        """The three-point Ward recursion must include the simple-pole term."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Genus-$1$ Virasoro Ward amplitudes;"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma+\delta$",
            r"\hypAmbientWtCpl+\hypVirTorusBlock",
            r"\eta(\tau)^{-\kappaChHodge(\mathrm{Vir}_c)}",
            r"-\frac{c}{48}E_2(\tau)",
            r"q^{-c/24}\prod_{n\ge2}(1-q^n)^{-1}",
            "extra minimal-model singular-vector quotients",
            r"P_1(z,\tau)=z^{-1}+O(z)",
            r"P_4(z,\tau)=\frac16\partial_z^2P_2(z,\tau)",
            r"A^{(1),\mathrm{reg}}_2(\tau)",
            r"P_1(z_{ij},\tau)\,\partial_{z_j}",
            r"A^{(1)}_2(T,T;\,z_j,z_k;\tau)",
            "not a full non-chiral physical torus partition",
            r"P_1(z_{12},\tau)\partial_{z_2}A^{(1)}_1(T;\tau)=0",
            r"P_1\partial A^{(1)}_2",
        )
        for needle in required:
            assert needle in block

        retired = (
            r"\begin{theorem}[Genus-$1$ Virasoro amplitudes: scalar trace and OPE part]",
            r"A^{(1)}_1(T;\tau)=-(c/24)\,E_2(\tau)",
            "The remaining term is a\nholomorphic torus block depending on the chosen Virasoro module; it is\nnot determined by the scalar Hodge trace.\n\\end{theorem}",
        )
        for needle in retired:
            assert needle not in block


class TestGenusOneModularAnomalyScope:
    """Guard the E2 modular anomaly and bar-transgression comparison."""

    def test_scope_separates_universal_E2_anomaly_from_kappa_weight(self):
        """At c=24, κ=12 and the determinant E2 coefficient is -1/2."""
        scope = genus1_modular_anomaly_scope(c=24)
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'hypEisComp')
        assert scope['E2_anomaly_coefficient'] == 12
        assert scope['E2_anomaly_is_universal'] is True
        assert scope['E2_anomaly_is_kappa_inverse'] is False
        assert scope['E2hat_is_modular_weight_2'] is True
        assert scope['q_dq_log_eta'] == 'E2/24'
        assert scope['bar_curvature_coefficient'] == 12
        assert scope['determinant_line_derivative_coefficient'] == Rational(-1, 2)
        assert scope['nonholomorphic_correction_requires_comparison'] is True
        assert scope['nonholomorphic_correction_is_transgression_generator'] is False
        assert scope['eta_square_equals_bar_curvature_asserted'] is False

    def test_modular_anomaly_proposition_is_conditionally_scoped(self):
        """The proposition must not claim a kappa^{-1} anomaly or eta^2 equality."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Modular anomaly of $E_2$ and the genus-$1$ bar"
        )
        prop_end = source.index(r"\end{proposition}", start)
        proof_end = source.index(r"\end{proof}", prop_end)
        block = source[start:proof_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma+\delta$",
            r"\hypAmbientWtCpl+\hypEisComp",
            r"d^2_{\bar B}=\kappaChHodge(\mathrm{Vir}_c)\,\omega_1",
            r"E_2(-1/\tau) = \tau^2 E_2(\tau) + 12\tau/(2\pi i)",
            r"The anomalous term \(12\tau/(2\pi i)\) is universal",
            r"-\kappaChHodge(\mathrm{Vir}_c)/24",
            r"\widehat{E}_2(\tau) := E_2(\tau) - 3/(\pi\,\mathrm{Im}\,\tau)",
            r"Under \(\hypEisComp\)",
            r"does not identify this scalar function with the algebra generator",
            r"does not assert the unlicensed equality",
            r"q\,\partial_q\log \eta(\tau)^{-\kappaChHodge(\mathrm{Vir}_c)}",
            "not a chain-level\nconstruction of the bar-transgression complex",
        )
        for needle in required:
            assert needle in block

        retired = (
            r"$\kappa^{-1}$ times the non-separating degeneration class",
            r"element $\eta$ of the bar-transgression complex",
            r"whose square $\eta^2 = u =",
            r"\kappa \cdot \omega_1$ kills the curvature",
            "is the transgression\n element",
        )
        for needle in retired:
            assert needle not in block


class TestGenusOneKZBShadowKernelScope:
    """Guard the genus-1 Virasoro KZB shadow-kernel scope."""

    def test_scope_is_connection_shadow_not_full_quantum_r_matrix(self):
        """At c=24, the scalar P2 coefficient is 12 and the rest is conditional."""
        scope = genus1_virasoro_kzb_shadow_kernel_scope(c=24)
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'hypKZBHeat')
        assert scope['object'] == 'connection-level KZB shadow kernel'
        assert scope['elliptic_scalar_kernel'] == '(c/2) P2(z,tau)'
        assert scope['scalar_kernel_coefficient'] == 12
        assert scope['P2_definition'] == 'P2=wp+E2/12'
        assert scope['contact_term'] == 'C_T^KZB(tau) T'
        assert scope['chosen_contact_normalization'] == 'C_T^KZB(tau)=2 E2(tau)'
        assert scope['contact_term_universal_without_kzb_heat'] is False
        assert scope['linearized_kzb_flatness_only'] is True
        assert scope['full_quantum_R_matrix_constructed'] is False
        assert scope['nonlinear_qybe_proved'] is False
        assert scope['higher_hbar_corrections_determined'] is False
        assert scope['central_tt_kernel_is_P4_not_P2'] is True

    def test_genus1_r_matrix_theorem_is_kzb_shadow_scoped(self):
        """The theorem must not construct a full modular quantum R-matrix."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Genus-$1$ Virasoro KZB shadow kernel;"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma+\delta$",
            r"\hypAmbientWtCpl+\hypKZBHeat",
            r"\nabla^{\mathrm{KZB},\mathrm{sh}}",
            r"P_2(z,\tau)=\wp(z;\tau)+\frac{E_2(\tau)}{12}",
            r"r^{\mathrm{KZB},\mathrm{sh}}_1(z;\tau)",
            r"\frac{c}{2}P_2(z,\tau)",
            r"\mathcal C_T^{\mathrm{KZB}}(\tau)\,T",
            r"\mathcal C_T^{\mathrm{KZB}}(\tau)=2E_2(\tau)",
            r"\((c/2)P_4\)",
            r"\nabla^{\mathrm{KZB}}_\tau r_0",
            "not construct the full quantum",
            "does not prove the\nnonlinear quantum Yang--Baxter equation",
            "does not determine the\nhigher",
            "theta/Kronecker heat equation",
            r"\cite{Bernard1988,CalaqueEnriquezEtingof2009}",
        )
        for needle in required:
            assert needle in block

        retired = (
            r"\begin{theorem}[Genus-$1$ spectral $R$-matrix]",
            r"R^{\mathrm{mod}}(z;\hbar,\tau)",
            r"R_0(z) + \hbar^2\, r_1(z;\tau) + O(\hbar^4)",
            r"2T \cdot E_2(\tau)",
            "one-loop graviton self-energy corrections",
            "integrability condition for the $\\tau$-dependent family of\n$R$-matrices",
        )
        for needle in retired:
            assert needle not in block


class TestGravityWeinbergWardResidueScope:
    """Guard the degree-2 Virasoro Ward residue and Weinberg comparison."""

    def test_scope_separates_ward_residue_from_physical_soft_factor(self):
        """The algebraic residue is W^(0), while Weinberg needs CelSoft data."""
        scope = gravity_weinberg_ward_residue_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'hypCelSoft')
        assert scope['algebraic_object'] == 'degree-2 Virasoro Ward residue'
        assert scope['local_residue_dictionary'] == (
            'W_i^(j)=Res (z-z_i)^j T(z) O_i dz = T_(j) O_i'
        )
        assert scope['simple_pole_residue'] == 'W_i^(0)=partial_{z_i}'
        assert scope['double_pole_residue'] == 'W_i^(1)=h_i'
        assert scope['physical_weinberg_factor'] == (
            'sum_i epsilon_{mu nu} p_i^mu p_i^nu/(q dot p_i)'
        )
        assert scope['physical_factor_is_spin_two'] is True
        assert scope['requires_celestial_bms_comparison'] is True
        assert scope['algebraic_residue_equals_physical_weinberg_without_comparison'] is False
        assert scope['translation_invariance_derives_spin_two_factor'] is False
        assert scope['old_mode_shift_T_p_plus_1_retired'] is True

    def test_weinberg_theorem_is_comparison_scoped(self):
        """The theorem must not identify a raw chiral residue with the 4d theorem."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        dict_start = source.index(r"\subsubsection*{The dictionary: pole order vs.\ soft order}")
        dict_end = source.index(r"\subsubsection*{Leading soft theorem from $m_2$}", dict_start)
        dictionary = source[dict_start:dict_end]
        start = source.index(
            r"\begin{theorem}[Degree-$2$ Virasoro Ward residue and Weinberg comparison;"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]

        dictionary_required = (
            r"W^{(j)}_i",
            r"T_{(j)}\cO_i(z_i)",
            r"(z - z_i)^j\,T(z)\,\cO_i(z_i)\,dz",
            r"not this index before the",
            r"\hypCelSoft",
        )
        for needle in dictionary_required:
            assert needle in dictionary

        theorem_required = (
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma+\delta$",
            r"\hypAmbientWtCpl+\hypCelSoft",
            r"S_2(z_0)",
            r"\frac{h_i}{(z_0-z_i)^2}",
            r"\frac{\partial_{z_i}}{z_0-z_i}",
            r"W_i^{(0)}=\Res_{z_0=z_i}S_2(z_0)\,dz_0=\partial_{z_i}",
            r"W_i^{(1)}=\Res_{z_0=z_i}(z_0-z_i)S_2(z_0)\,dz_0=h_i",
            r"\mathsf{CelSoft}_{\hypCelSoft}",
            r"\varepsilon_{\mu\nu}p_i^\mu p_i^\nu",
            r"q\cdot p_i",
            "Without\n\\(\\hypCelSoft\\)",
            "Virasoro Ward-residue identity, not a derivation",
            r"\cite{Weinberg65,HMPS15,Strominger18}",
        )
        for needle in theorem_required:
            assert needle in block

        retired = (
            r"T_{(p+1)}\cO_i(z_i)",
            r"(z - z_i)^{p+1}",
            r"\begin{theorem}[Weinberg from the degree-$2$ shadow]",
            r"S^{(0)}_i",
            r"\frac{\partial_{z_i}}{z_0 - z_i}",
            "translation invariance\n$\\sum_i \\partial_{z_i}",
            r"(p_i\cdot\epsilon)/(p_i\cdot q)",
            "forces an infinite sequence of new soft\nWard identities",
        )
        for needle in retired:
            assert needle not in dictionary
            assert needle not in block

    def test_genus1_summary_uses_eiscomp_for_e2_completion(self):
        """The E2 summary must not identify the scalar completion with a generator."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(r"\emph{1.\ Eisenstein series $E_2(\tau)$.}")
        end = source.index(r"\emph{2.\ Dedekind $\eta$-function.}", start)
        block = source[start:end]
        assert r"\hypEisComp" in block
        assert "analytic image of the completed transgression" in block
        assert "transgression element killing the bar curvature" not in block


class TestGravityCachazoStromingerWardPackageScope:
    """Guard the degree-2 conformal Ward package entering subleading soft."""

    def test_scope_separates_double_pole_from_physical_cs_factor(self):
        """The double pole is a weight term; physical CS needs CelSoft and degree 3."""
        scope = gravity_cachazo_strominger_ward_package_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'hypCelSoft')
        assert scope['algebraic_object'] == (
            'degree-2 holomorphic stress-tensor Ward package'
        )
        assert scope['ward_operator'] == (
            'sum_i (Y(z_i) partial_{z_i} + h_i partial Y(z_i))'
        )
        assert scope['double_pole_residue'] == 'W_i^(1)=h_i'
        assert scope['simple_pole_residue_required'] == 'W_i^(0)=partial_{z_i}'
        assert scope['double_pole_role'] == 'conformal-weight term only'
        assert scope['physical_cachazo_strominger_factor'] == (
            'sum_i epsilon_{mu nu} p_i^mu q_rho J_i^{rho nu}/(q dot p_i)'
        )
        assert scope['requires_antiholomorphic_completion'] is True
        assert scope['requires_degree_three_shadow_channel'] is True
        assert scope['double_pole_alone_equals_physical_cs'] is False
        assert scope['global_conformal_ward_identity_equals_physical_cs_without_comparison'] is False

    def test_cs_theorem_is_comparison_boundary_not_raw_double_pole(self):
        """The theorem must not identify h_i/(z-z_i)^2 with Cachazo-Strominger."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Degree-$2$ conformal Ward package and the"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        table_start = source.index(r"\begin{center}", source.index(r"\begin{remark}[Novelty of the higher soft shadows]"))
        table_end = source.index(r"\end{center}", table_start)
        table = source[table_start:table_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma+\delta$",
            r"\hypAmbientWtCpl+\hypCelSoft",
            r"\mathcal W_Y^{(2),\mathrm{hol}}",
            r"Y(z_i)\partial_{z_i}",
            r"h_i\,\partial Y(z_i)",
            r"W_i^{(1)}=h_i",
            r"W_i^{(0)}=\partial_{z_i}",
            r"\mathsf{CelSoft}_{\hypCelSoft}",
            r"degree-$3$ shadow channel",
            r"S_{\mathrm{CS}}^{(1)}",
            r"\varepsilon_{\mu\nu}p_i^\mu q_\rho J_i^{\rho\nu}",
            r"q\cdot p_i",
            "The double-pole residue alone is therefore not the momentum-space\nsubleading theorem",
            r"\cite{CachazoStrominger14,PateRaclariuStromingerYuan19,Strominger18}",
        )
        for needle in required:
            assert needle in block

        retired = (
            r"\begin{theorem}[Cachazo--Strominger from the degree-$2$",
            "The double-pole residue of $S_2$ gives the subleading\nsoft graviton theorem",
            r"S^{(1)}_i",
            r"\frac{h_i}{(z_0 - z_i)^2}",
            "gives the global conformal Ward identity, which is the\nsubleading soft graviton theorem",
            r"L_0 + \bar{L}_0",
            r"$S^{(1)}$ & $S_2|_{\text{double pole}}$ & $m_2$",
        )
        for needle in retired:
            assert needle not in block
            assert needle not in table

        assert r"$S^{(1)}$ & $S_3$; $S_2|_{\text{double pole}}$ only the weight term" in table

    def test_split_movement_summary_uses_repaired_soft_hierarchy(self):
        """The split Movement VI-X file must not sell the old soft dictionary."""
        root = Path(__file__).resolve().parents[2]
        source = (
            root / "chapters" / "connections" / "thqg_3d_gravity_movements_vi_x.tex"
        ).read_text()
        start = source.index("The detailed derivation of the soft graviton comparison")
        end = source.index(r"\S\ref{subsec:gravity-soft-graviton-ainf}", start)
        block = source[start:end]

        required = (
            "the degree-$2$ Virasoro Ward package and Weinberg comparison",
            "Cachazo--Strominger superrotation channel from the degree-$3$\nshadow",
            "together with the $S_2$ weight term",
            "the quartic\nsub-subleading contact channel",
        )
        for needle in required:
            assert needle in block

        retired = (
            "Weinberg from the simple-pole residue of $S_2$",
            "Cachazo--Strominger from the double-pole residue",
            "Cachazo--He--Yuan from the ternary operation $m_3$",
        )
        for needle in retired:
            assert needle not in block


class TestGravityCHYQuarticContactScope:
    """Guard the quartic-contact scope of the sub-subleading soft channel."""

    def test_scope_separates_m3_precursor_from_quartic_contact_channel(self):
        """The physical sub-subleading comparison is not the raw ternary operation."""
        scope = gravity_chy_quartic_contact_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'hypCelSoft')
        assert scope['algebraic_object'] == 'degree-4 quartic contact sub-subleading shadow'
        assert scope['binary_primary_T2_residue'] == 'T_(2) O_i = 0'
        assert scope['ternary_role'] == 'm3 is cubic precursor entering [C,C] at degree 4'
        assert scope['virasoro_cubic_shadow_gauge_trivial'] is True
        assert scope['first_nontrivial_virasoro_graviton_channel'] == 'S4 quartic contact'
        assert scope['quartic_contact_coefficient'] == '10/(c(5c+22))'
        assert scope['large_c_asymptotic'] == '2/c^2'
        assert scope['degree_four_shadow_formula'] == 'Sh(Q)|soft + [Sh(C),Sh(C)]|soft'
        assert scope['physical_chy_hsl_requires_celsoft'] is True
        assert scope['physical_cs_hamada_shiu_requires_celsoft'] is True
        assert scope['m3_alone_equals_physical_chy'] is False
        assert scope['m3_alone_equals_physical_subsubleading'] is False
        assert scope['binary_residue_equals_subsubleading'] is False

    def test_chy_theorem_is_quartic_contact_scoped(self):
        """The theorem must not advertise m3 itself as the soft theorem."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Quartic contact sub-subleading shadow and"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        table_start = source.index(r"\begin{center}", source.index(r"\begin{remark}[Novelty of the higher soft shadows]"))
        table_end = source.index(r"\end{center}", table_start)
        table = source[table_start:table_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing",
            r"\hypAmbientWtCpl+\hypCelSoft",
            "not a binary Ward residue and is not the ternary operation by itself",
            r"S^{(2)}_i|_{\text{binary}} = 0",
            r"\FM_3(\C)",
            r"[\Sh(\mathfrak C),\Sh(\mathfrak C)]",
            r"S^{(2),\mathrm{alg}}_{\mathrm{Vir}}",
            r"\mathfrak{Q}^{\mathrm{contact}}_{\mathrm{Vir}}",
            r"\frac{10}{c(5c+22)}",
            r"S^{(2)}",
            r"\Sh_{0,n}(\mathfrak Q(\mathrm{Vir}_c))",
            r"\bigl[\Sh_{0,n}(\mathfrak C(\mathrm{Vir}_c)),",
            r"\mathsf{CelSoft}_{\hypCelSoft}",
            "Cachazo--Strominger / Hamada--Shiu sub-subleading\nsoft-graviton package",
            r"\cite{CachazoStrominger14,HamadaShiuSubsubSoft}",
        )
        for needle in required:
            assert needle in block

        table_required = (
            r"$S^{(2)}$ & \(S_4\) plus \([S_3,S_3]\); for Virasoro \(S_3\) is exact",
            r"\(m_4\) plus \([m_3,m_3]\)",
            "Cubic exact / quartic non-trivial",
        )
        for needle in table_required:
            assert needle in table

        retired = (
            r"\begin{theorem}[Conditional sub-subleading soft shadow from the ternary operation]",
            "It requires the\nternary $\\Ainf$ operation $m_3$",
            "quadratic in the translation generator, the algebraic incarnation\nof the $J^{\\mu\\nu}J^{\\rho\\sigma}$ structure",
            r"$S^{(2)}$ & $S_3$ (primary), $S_4$ (graviton)",
            r"$m_3$, $m_4$ & Gauge-trivial / non-trivial",
            "Cachazo--He--Yuan from the ternary operation $m_3$",
            "Hamada--Shiu--Li--Strominger",
            r"\cite{HamadaShiuSubsubSoft,LiStromingerHigherSoft}",
        )
        for needle in retired:
            assert needle not in block
            assert needle not in table


class TestGravityInfiniteSoftShadowHierarchyScope:
    """Guard the infinite Virasoro soft-shadow hierarchy scope."""

    def test_shadow_metric_coefficients_compute_s5(self):
        """The normalized scalar S5 comes from the metric branch, not raw m5."""
        c = Symbol('c')
        coeffs = virasoro_shadow_metric_coefficients(max_r=5)
        assert simplify(coeffs[2] - c / 2) == 0
        assert simplify(coeffs[3] - 2) == 0
        assert simplify(coeffs[4] - S(10) / (c * (5 * c + 22))) == 0
        assert simplify(coeffs[5] + S(48) / (c**2 * (5 * c + 22))) == 0

    def test_scope_separates_support_shadow_from_physical_soft_tower(self):
        """The class-M hierarchy is algebraic until CelSoft is supplied."""
        scope = gravity_infinite_soft_shadow_hierarchy_scope()
        c = Symbol('c')
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'hypCelSoft')
        assert scope['ambient_object'] == 'completed class-M Virasoro support-shadow packet'
        assert scope['soft_order_to_arity'] == 'p maps to r=p+2 after CelSoft comparison'
        assert scope['degree_r_channel'] == 'normalized Sh(Theta_r) plus lower MC brackets'
        assert scope['raw_operation_alone_controls_soft_order'] is False
        assert scope['higher_physical_soft_theorem_asserted_without_celsoft'] is False
        assert scope['independent_primitive_classes_asserted'] is False
        assert simplify(scope['S5'] + S(48) / (c**2 * (5 * c + 22))) == 0
        assert scope['quintic_bracket_subchannel'] == '20/(c^2(5c+22))'
        assert scope['table_source_language'] == (
            'Theta_(p+2) projection with lower MC bracket terms, not m_(p+2) alone'
        )

    def test_infinite_soft_tower_theorem_is_shadow_projection_scoped(self):
        """The theorem must not identify raw m_(p+2) with physical soft order."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Infinite Virasoro soft-shadow hierarchy"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        table_start = source.index(
            r"\begin{center}",
            source.index(r"\begin{remark}[Novelty of the higher soft shadows]"),
        )
        table_end = source.index(r"\end{center}", table_start)
        table = source[table_start:table_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing",
            r"\hypAmbientWtCpl+\hypCelSoft",
            r"Q_{\mathrm{Vir}}(t)=(c+6t)^2+\frac{80t^2}{5c+22}",
            r"\sqrt{Q_{\mathrm{Vir}}(0)}=c",
            r"it is not the raw",
            r"transferred operation \(m_r\) by itself",
            r"S_r=[t^r]\,H(t)/r",
            "does not assert\n infinitely many independent primitive generators",
            r"S_5=-\frac{48}{c^2(5c+22)}",
            r"20/[c^2(5c+22)]",
            r"\mathsf{CelSoft}_{\hypCelSoft}",
            "not a physical higher soft theorem",
            r"c=3\ell/(2G_N)",
            r"Computation",
            r"\ref{comp:S5-verification}",
        )
        for needle in required:
            assert needle in block

        table_required = (
            r"$S^{(3)}$ & \(S_5\) plus degree-$5$ bracket sub-channels",
            r"\(\Theta_5\) after shadow projection; raw \(m_5\) is only a",
            r"degree-\(p+2\) normalized shadow projection plus lower brackets",
            r"\(\Theta_{p+2}\) in the degree-\(p+2\) MC equation, not raw",
            r"No finite class-\(\mathbf M\) support cutoff",
        )
        for needle in table_required:
            assert needle in table

        retired = (
            r"\begin{theorem}[Conditional soft-shadow hierarchy from the",
            r"Soft order $p$ is controlled by the shadow degree",
            r"via the $\Ainf$ operation $m_{p+2}$",
            "independent algebraic shadow classes",
            r"$S^{(3)}$ & $S_5$ & $m_5$ via $o_5 \neq 0$",
            r"$S^{(p)}$ & $S_{p+2}$ & $m_{p+2}$ via $o_{p+2}$",
            r"Non-trivial for $p$ even",
        )
        for needle in retired:
            assert needle not in block
            assert needle not in table

    def test_split_soft_graviton_higher_hierarchy_is_support_scoped(self):
        """The split soft-graviton file must not count physical Ward identities."""
        root = Path(__file__).resolve().parents[2]
        source = (
            root / "chapters" / "connections" / "thqg_soft_graviton_theorems.tex"
        ).read_text()
        intro_start = source.index(
            r"\subsection{Higher celestial soft factors and the infinite tower}"
        )
        intro_end = source.index(
            r"\begin{proposition}[Inductive non-termination for class $\mathbf{M}$",
            intro_start,
        )
        intro = source[intro_start:intro_end]
        prop_start = intro_end
        prop_end = source.index(r"\end{proposition}", prop_start)
        prop = source[prop_start:prop_end]

        intro_required = (
            "the algebraic support-shadow hierarchy has\nno finite arity cutoff",
            "After a celestial comparison datum is fixed",
            "not a count of\nindependent physical Ward identities",
        )
        for needle in intro_required:
            assert needle in intro

        prop_required = (
            r"\ClaimStatusConditional",
            r"licensing",
            r"\hypAmbientWtCpl+\hypCelSoft",
            "detects\n infinitely many support channels",
            "not an\n independent tower of physical Ward identities",
        )
        for needle in prop_required:
            assert needle in prop

        retired = (
            "the celestial soft hierarchy is an\ninfinite tower of coupled Ward identities",
            "this gives the corresponding\n infinite hierarchy of soft Ward identities",
        )
        for needle in retired:
            assert needle not in intro
            assert needle not in prop


class TestShadowClosedFormScope:
    """Guard the closed-form scalar shadow formula and Catalan identity."""

    def test_closed_form_matches_metric_branch_through_degree_14(self):
        """The Catalan formula is the scalar metric branch, not raw m_r."""
        c = Symbol('c')
        coeffs = virasoro_shadow_metric_coefficients(max_r=14, c=c)
        for r in range(4, 15):
            closed = virasoro_shadow_closed_form_coefficient(r, c=c)
            assert simplify(coeffs[r] - closed) == 0

    def test_shape_factor_generating_function(self):
        """The displayed algebraic generating function has the stated coefficients."""
        x, t = symbols('x t')
        gf = (
            sqrt((1 - t)**2 + 4 * x * t**2) - (1 - t)
        ) / (2 * x * t**2)
        expansion = series(gf, t, 0, 8).removeO().expand()
        for n in range(8):
            expected = virasoro_catalan_shape_factor(n + 4, x=x)
            assert simplify(expansion.coeff(t, n) - expected) == 0

    def test_half_binomial_identity_is_correct(self):
        """The Catalan half-binomial identity has denominator 2^(2m-1)."""
        for m in range(1, 8):
            correct = (-1)**(m - 1) * S(math.factorial(2 * m - 2)) / (
                S(2)**(2 * m - 1)
                * S(math.factorial(m - 1))
                * S(math.factorial(m))
            )
            old = (-1)**(m - 1) * S(math.factorial(2 * m - 2)) / (
                S(4)**m
                * S(m)
                * S(math.factorial(m - 1))
                * S(math.factorial(m))
            )
            assert simplify(binomial(S(1) / 2, m) - correct) == 0
            assert simplify(binomial(S(1) / 2, m) - old) != 0

    def test_scope_is_normalized_scalar_t_line_branch(self):
        """The helper records the hypotheses and excludes raw operations."""
        scope = shadow_closed_form_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('gamma', 'epsilon')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'effScalarShadowProj')
        assert scope['object'] == 'normalized scalar T-line shadow coefficient'
        assert scope['branch'] == 'sqrt(Q_Vir(0))=c'
        assert scope['excluded_central_charges'] == (S.Zero, -S(22) / 5)
        assert scope['Q_Vir'] == '(c+6t)^2 + 80 t^2/(5c+22)'
        assert scope['coefficient_rule'] == 'S_r=[t^r]H(t)/r'
        assert scope['raw_transferred_mr_formula'] is False
        assert scope['full_ordered_bar_invariant_before_projection'] is False
        assert scope['nonconstant_pole_divisor'] == (
            'c^(r-3)(5c+22)^floor((r-2)/2)'
        )
        assert scope['rational_scalar_denominators_are_poles'] is False

    def test_shadow_closed_form_theorem_is_scalar_branch_scoped(self):
        """The theorem must not advertise an unlicensed full shadow invariant."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Closed-form normalized scalar shadow coefficient formula;"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        remark_start = source.index(
            r"\begin{remark}[Structure of the closed form]", proof_end
        )
        remark_end = source.index(r"\end{remark}", remark_start)
        remark = source[remark_start:remark_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\gamma+\varepsilon$",
            r"\hypAmbientWtCpl+\effScalarShadowProj",
            r"normalized scalar \(T\)-line",
            r"\sqrt{Q_{\mathrm{Vir}}(0)}=c",
            r"c\ne0,-22/5",
            r"S_r=[t^r]H(t)/r",
            r"not a formula for the raw transferred operation \(m_r\)",
            r"not a formula for the full ordered-bar invariant",
            "nonconstant pole divisor",
            "rational scalar",
            r"u=-6t/c",
            r"\binom{1/2}{m}=(-1)^{m-1}\frac{C_{m-1}}{2^{2m-1}}",
        )
        for needle in required:
            assert needle in block

        remark_required = (
            "normalized\nscalar shadow coefficient",
            "nonconstant pole divisor",
            "source of the $(5c+22)$ pole-divisor factors",
        )
        for needle in remark_required:
            assert needle in remark

        retired = (
            r"\begin{theorem}[Closed-form shadow coefficient formula]",
            r"4^m \cdot m",
            "ordered bar complex\ninvariants",
            "the denominator $c^{r-3}$",
        )
        for needle in retired:
            assert needle not in block
            assert needle not in remark

    def test_advertised_catalan_copies_do_not_use_retired_identity(self):
        """Related Catalan scripts must not preserve the false denominator."""
        root = Path(__file__).resolve().parents[2]
        gravity = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        ordered = (root / "compute" / "ordered_e1_shadow_catalan.py").read_text()
        wn = (root / "compute" / "wn_universal_catalan.py").read_text()
        texts = [gravity, ordered, wn]
        retired = (
            "4^m \\cdot m",
            "4^n * n * (2n-3)!!/n!",
            "4^k * k * (2k-1)",
            "The denominator is exactly",
            "E_1 ordered bar complex\ninvariants",
            "MAIN THEOREM (Proved here computationally",
            "for ALL principal W_N",
            "This gets complicated",
            "Wait,",
        )
        for text in texts:
            for needle in retired:
                assert needle not in text
        assert "conditional Catalan pattern" in wn
        assert "not a proof of universal\nprincipal-W_N closed forms" in wn


class TestCatalanDynkinParityScope:
    """Guard Catalan-Dynkin parity as a scoped symmetric-point recursion."""

    def test_catalan_dynkin_polynomials_satisfy_rightmost_recursion(self):
        """The h=2 Dynkin string is the formal solution of the recursion."""
        x = Symbol('x')
        for arity in range(4, 12):
            rhs = S.Zero
            for j in range(2, arity):
                rhs += (
                    (-1)**(arity - j)
                    * catalan_dynkin_field_polynomial(j, x + arity - j)
                    * catalan_dynkin_field_polynomial(arity + 1 - j, x)
                )
            rhs = simplify(-rhs)
            assert simplify(rhs - catalan_dynkin_field_polynomial(arity, x)) == 0

    def test_catalan_dynkin_even_vanishes_and_odd_roots_are_fixed(self):
        """Even arities vanish; odd arities have roots -2,...,-k."""
        x = Symbol('x')
        for arity in (4, 6, 8, 10, 12):
            assert catalan_dynkin_field_polynomial(arity, x) == 0

        for arity in (3, 5, 7, 9, 11):
            poly = catalan_dynkin_field_polynomial(arity, x)
            n = (arity - 3) // 2
            cat = math.comb(2 * n, n) // (n + 1)
            expected = S((-1)**n * cat)
            for m in range(2, arity + 1):
                expected *= x + m
                assert simplify(poly.subs(x, -m)) == 0
            assert simplify(poly - expected) == 0

    def test_catalan_dynkin_scope_records_non_universal_status(self):
        """The scope guard must reject the old arbitrary class-M reading."""
        scope = catalan_dynkin_parity_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'gamma', 'epsilon')
        assert 'chosen one-generator branch' in scope['hypotheses']
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'effKoszul' in scope['hypotheses']
        assert scope['object'] == 'symmetric-point field polynomial'
        assert scope['full_ordered_spectral_polynomial_statement'] is False
        assert scope['all_one_generator_class_M_statement'] is False
        assert scope['arbitrary_conformal_weight_statement'] is False
        assert scope['diagonal_reflection_quotient_only'] is True

    def test_catalan_dynkin_theorem_is_scoped(self):
        """The theorem must not assert parity for every one-generator class-M algebra."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Catalan--Dynkin parity for a Virasoro-normalized"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\gamma+\varepsilon$",
            r"\hypAmbientWtCpl+\effKoszul",
            r"chosen one-generator branch",
            r"\varphi_2^a(x)=x+2",
            r"\varphi_3^a(x)=(x+2)(x+3)",
            r"\varphi_{2n+3}^a(x)=(-1)^nC_n\prod_{m=2}^{2n+3}(x+m)",
            r"\varphi_{2q}^a(x)=0",
            r"This is not a theorem about every one-generator class-$\mathbf{M}$ algebra",
            r"diagonal quotient",
            r"not a statement about the full ordered spectral polynomial",
        )
        for needle in required:
            assert needle in block

        retired = (
            r"Let $\cA$ be any one-generator $\SCchtop$-algebra of class",
            r"generator of conformal weight~$s$",
            r"holds for every one-generator class-$\mathbf{M}$ algebra",
            r"at even degrees this involution acts by $(-1)$, killing the field polynomial",
            r"at odd degrees it acts trivially, preserving it",
        )
        for needle in retired:
            assert needle not in block


class TestCrossingStasheffScope:
    """Guard Stasheff crossing as a chiral channel identity, not full bootstrap."""

    def test_scope_separates_chiral_channel_from_full_bootstrap(self):
        """The scope guard must keep analytic bootstrap data outside d_B^2=0."""
        scope = crossing_stasheff_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma')
        assert 'ordered FM_3 collision chart' in scope['hypotheses']
        assert 'chiral OPE-channel comparison' in scope['hypotheses']
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert scope['object'] == 'singular-OPE channel identity'
        assert scope['requires_fourth_state_pairing'] is True
        assert scope['full_conformal_bootstrap_statement'] is False
        assert scope['physical_crossing_requires_delta_data'] is True
        assert scope['minimal_m2_nonassociativity_statement'] is False
        assert scope['positivity_or_unitarity_statement'] is False

    def test_virasoro_contact_cancels_associator(self):
        """The targeted Virasoro reading is m3=-A3 in the chosen gauge."""
        samples = (
            (26, 2, 3),
            (Rational(5, 2), Rational(1, 3), Rational(2, 5)),
            (Symbol('c'), Symbol('lambda_12'), Symbol('lambda_23')),
        )
        for c_val, lam12, lam23 in samples:
            associator = virasoro_associator(c=c_val, lam12=lam12, lam23=lam23)
            contact = virasoro_m3_coefficients(c=c_val, lam12=lam12, lam23=lam23)
            for key in ('d2T', 'dT', 'T', 'scalar'):
                assert simplify(associator[key] + contact[key]) == 0

    def test_crossing_theorem_is_scoped(self):
        """The theorem must not identify d_B^2=0 with full bootstrap crossing."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Chiral OPE-channel crossing from the Stasheff relation"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        remark_end = source.index(r"\end{remark}", proof_end)
        block = source[start:remark_end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma$",
            r"\hypAmbientWtCpl",
            r"ordered \(\FM_3(\C)\) collision chart",
            r"chiral OPE-channel comparison",
            r"singular-OPE channel identity",
            r"after pairing the output with a fourth state",
            r"not the full conformal bootstrap crossing equation",
            r"analytic conformal-block convergence",
            r"positivity",
            r"not a statement that a strict minimal \(m_2\) fails to be associative",
            r"m_3=-A_3",
        )
        for needle in required:
            assert " ".join(needle.split()) in flat

        retired = (
            r"\begin{theorem}[Crossing = $d_{\barB}^2 = 0$]",
            r"is the four-point crossing equation of the conformal bootstrap",
            r"On cohomology, the relation becomes exactly the four-point",
            r"m_3^H \ne 0 and the crossing equation acquires",
            r"The $\Ainf$ tower solves the full bootstrap hierarchy",
            r"the bar complex gives exact values",
        )
        for needle in retired:
            assert " ".join(needle.split()) not in flat


class TestShapovalovBootstrapScope:
    """Guard Shapovalov bootstrap input as channel positivity, not coordinate bounds."""

    def test_projected_channel_norm_is_contractive(self):
        """Orthogonal quotient projection decreases the Shapovalov channel norm."""
        coefficients = (3, 4, 12)
        gram = (
            (1, 0, 0),
            (0, 2, 0),
            (0, 0, 3),
        )
        full = shapovalov_channel_norm_squared(coefficients, gram=gram)
        projected = shapovalov_projected_channel_norm_squared(
            coefficients,
            kept_indices=(0, 2),
            gram=gram,
        )
        assert full == S(473)
        assert projected == S(441)
        assert simplify(full - projected) == S(32)
        assert projected <= full

    def test_scope_records_channel_norm_not_coordinatewise_bound(self):
        """The scope guard must reject a pointwise OPE-coefficient bound."""
        scope = shapovalov_bootstrap_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'beta', 'gamma', 'delta')
        assert 'finite-level Shapovalov normalisation' in scope['hypotheses']
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'unitary Hilbert quotient' in scope['hypotheses']
        assert scope['object'] == 'finite-level channel Gram norm'
        assert scope['propagator_source'] == (
            'inverse finite-level Shapovalov Gram matrix'
        )
        assert scope['coordinatewise_ope_bound_asserted'] is False
        assert scope['channel_norm_contraction_asserted'] is True
        assert scope['unitary_minimal_models_only'] == (
            'c=1-6/(m(m+1)), m>=2'
        )
        assert scope['generic_minimal_model_unitarity_asserted'] is False
        assert scope['raw_verma_bar_implies_bpz_constraints'] is False
        assert scope['raw_hpl_summands_may_have_kac_poles'] is True
        assert scope['scalar_shadow_projection_is_separate'] is True

    def test_bootstrap_shapovalov_theorem_is_channel_scoped(self):
        """The theorem must state Gram contraction, not coefficient-wise bounds."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Shapovalov Gram positivity and channel contraction"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma+\delta$",
            r"finite-level Shapovalov normalisation",
            r"\hypAmbientWtCpl",
            r"unitary Hilbert quotient",
            r"chiral OPE-channel vector",
            r"\bigl\|P_N v_{ij}^{(N)}\bigr\|_{G_N}^2",
            r"\bigl\|v_{ij}^{(N)}\bigr\|_{G_N}^2",
            r"does not imply a coefficient-wise inequality",
            r"unitary minimal-model",
            r"m\ge2",
            r"raw Verma bar complex alone",
            r"Only after the scalar shadow projection",
        )
        for needle in required:
            assert " ".join(needle.split()) in flat

        retired = (
            r"\begin{theorem}[Bootstrap bounds from the Shapovalov form]",
            r"|C^k_{ij,\,\mathrm{unitary}}|^2",
            r"|C^k_{ij,\,\mathrm{Verma}}|^2",
            r"coefficients are bounded above by the Verma-module values",
            r"by projecting out null states",
            r"or $c = c_{p,q}$ with $h = h_{r,s}$",
        )
        for needle in retired:
            assert " ".join(needle.split()) not in flat


class TestLargeCBootstrapScope:
    """Guard the large-c statement as a scalar shadow/contact branch."""

    def test_large_c_asymptotics_from_scalar_metric_branch(self):
        """Leading constants come from S_r=[t^r]H/r, not a prose table."""
        asymptotics = virasoro_large_c_shadow_asymptotics(max_r=9)
        expected = {
            4: S(2),
            5: -S(48) / 5,
            6: S(48),
            7: -S(1728) / 7,
            8: S(1296),
            9: -S(6912),
        }
        for r, constant in expected.items():
            row = asymptotics[r]
            assert simplify(row['leading_constant'] - constant) == 0
            assert simplify(row['expected_leading_constant'] - constant) == 0
            assert row['decay_power'] == 2 - r
        assert asymptotics[4]['error_order'] == 'O(c^(-3))'
        assert asymptotics[9]['error_order'] == 'O(c^(-8))'

    def test_scope_separates_scalar_branch_from_full_large_c_bootstrap(self):
        """The scalar identity/contact lane is not the full analytic bootstrap."""
        scope = large_c_bootstrap_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('gamma', 'epsilon', 'delta')
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'effScalarShadowProj' in scope['hypotheses']
        assert 'identity-block contact comparison' in scope['hypotheses']
        assert scope['object'] == 'scalar Virasoro T-line shadow branch'
        assert scope['large_c_asymptotic'] == (
            'S_r=8*(-6)^(r-4)*c^(2-r)/r+O(c^(1-r))'
        )
        assert scope['S4_leading'] == '2/c^2'
        assert scope['S5_leading'] == '-48/(5c^3)'
        assert scope['scalar_radius'] == 'c*sqrt((5c+22)/(180c+872)) ~ c/6'
        assert scope['full_large_c_bootstrap_statement'] is False
        assert scope['non_vacuum_blocks_controlled'] is False
        assert scope['single_valued_crossing_asserted'] is False
        assert scope['positivity_or_spectral_density_asserted'] is False
        assert scope['identity_contact_lane_only_after_comparison'] is True
        assert scope['old_coefficient_exponent'] == 'not [t^r]H=O(c^(3-r))'

    def test_large_c_bootstrap_proposition_is_scalar_scoped(self):
        """The proposition must not identify the scalar tower with the full bootstrap."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Scalar large-$c$ Virasoro shadow branch"
        )
        theorem_end = source.index(r"\end{proposition}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\gamma+\varepsilon+\delta$",
            r"\hypAmbientWtCpl+\effScalarShadowProj",
            r"identity-block contact comparison",
            r"H(t)=t^2\sqrt{Q_{\mathrm{Vir}}(t)}",
            r"S_r=[t^r]H(t)/r",
            r"S_r",
            r"\frac{8(-6)^{r-4}}{r}\,c^{2-r}",
            r"O(c^{1-r})",
            r"R_{\mathrm{scal}}",
            r"\sim c/6",
            r"not a statement that the full large-\(c\)",
            r"Non-vacuum blocks",
            r"positivity",
            r"analytic convergence",
        )
        for needle in required:
            assert " ".join(needle.split()) in flat

        retired = (
            r"\begin{proposition}[Large-$c$ expansion from the shadow obstruction tower]",
            r"matches the large-$c$ conformal bootstrap",
            r"The $n$-th order correction to the bootstrap is controlled by $S_{n+2}$",
            r"the shadow obstruction tower is the $1/c$ expansion of the bootstrap",
            r"[t^r]H = O(c^{3-r})",
        )
        for needle in retired:
            assert " ".join(needle.split()) not in flat


class TestOTOCRMatrixScope:
    """Guard OTOC/R-matrix input as block monodromy, not full chaos."""

    def test_diagonal_braiding_phase_and_inverse_orientation(self):
        """The diagonal phase is only an oriented braid specialization."""
        neutral = otoc_braiding_phase(
            h_p=Rational(4, 3),
            h_v=Rational(1, 2),
            h_w=Rational(5, 6),
        )
        assert simplify(neutral - 1) == 0

        positive = otoc_braiding_phase(
            h_p=Rational(5, 3),
            h_v=Rational(1, 2),
            h_w=Rational(5, 6),
        )
        inverse = otoc_braiding_phase(
            h_p=Rational(5, 3),
            h_v=Rational(1, 2),
            h_w=Rational(5, 6),
            orientation='inverse',
        )
        assert simplify(positive * inverse - 1) == 0

        with pytest.raises(ValueError):
            otoc_braiding_phase(1, 0, 0, orientation='clockwise')

    def test_scope_is_block_monodromy_not_full_otoc(self):
        """The R-matrix supplies block monodromy, not thermal chaos by itself."""
        scope = otoc_r_matrix_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'hypVirTorusBlock' in scope['hypotheses']
        assert 'conformal-block local system' in scope['hypotheses']
        assert 'chosen OTO continuation path' in scope['hypotheses']
        assert 'thermal trace coefficients' in scope['hypotheses']
        assert 'spectral-line R-matrix comparison' in scope['hypotheses']
        assert scope['object'] == 'chiral conformal-block monodromy operator'
        assert scope['diagonal_phase'] == (
            'exp(2*pi*i*(h_p-h_V-h_W)) for positive orientation'
        )
        assert scope['full_normalized_otoc_asserted'] is False
        assert scope['thermal_coefficients_computed_by_r_matrix'] is False
        assert scope['anti_holomorphic_sector_included'] is False
        assert scope['lyapunov_exponent_determined'] is False
        assert scope['scalar_phase_without_diagonalization'] is False
        assert scope['boltzmann_weight_formula_asserted'] is False
        assert scope['monodromy_matrix_required_generically'] is True

    def test_otoc_r_matrix_theorem_is_monodromy_scoped(self):
        """The theorem must not sell R-matrix monodromy as the full OTOC."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[OTOC block monodromy from the \(R\)-matrix"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma+\delta$",
            r"\hypAmbientWtCpl+\hypVirTorusBlock",
            r"conformal-block local system",
            r"torus Virasoro block/thermal trace normalisation",
            r"OTO continuation path",
            r"thermal trace data",
            r"F_{\mathrm{OTOC}}^{\chi,\mathrm{mon}}",
            r"\rho_{\mathrm{blk}}(\gamma_{\mathrm{OTO}})_{p}^{\,q}",
            r"multiplicity-free and diagonalizes the braiding",
            r"\exp\!\bigl(2\pi i(h_p-h_V-h_W)\bigr)",
            r"opposite OTO orientation gives the inverse phase",
            r"does not determine the thermal coefficients",
            r"anti-holomorphic sector",
            r"Lyapunov exponent",
            r"monodromy is a matrix, not a scalar phase",
        )
        for needle in required:
            assert " ".join(needle.split()) in flat

        retired = (
            r"\begin{theorem}[Monodromy contribution to the OTOC block expansion]",
            r"F_{\mathrm{OTOC}}^{\mathrm{mon}}(t)",
            r"\sum_p\, |C_{VWp}|^2",
            r"e^{-2\pi i(h_p - h_V - h_W)}",
            r"e^{-\beta h_p}\, / Z(\beta)",
            r"the OTOC receives this $R$-matrix monodromy contribution inside the thermal conformal-block sum",
            r"the phase becomes",
            r"with thermal weights $e^{-\beta h_p}$",
        )
        for needle in retired:
            assert " ".join(needle.split()) not in flat

        mss_start = source.index(r"\noindent\textbf{Vacuum block and the \(R\)-matrix trace.}")
        mss_end = source.index(r"For the identity exchange", mss_start)
        mss = " ".join(source[mss_start:mss_end].split())
        assert "does not by itself determine the Lorentzian growth" in mss
        assert "that growth comes from the analytic continuation of the Virasoro block" in mss
        assert "the phase becomes" not in mss


class TestMSSAnnularBarScope:
    """Guard the MSS theorem as an analytic strip theorem with annular input."""

    def test_mss_bound_value_scales_as_inverse_beta(self):
        """The analytic bound is exactly 2*pi/beta."""
        assert simplify(mss_bound_value(2) - pi) == 0
        assert simplify(mss_bound_value(Rational(1, 2)) - 4 * pi) == 0
        with pytest.raises(ValueError):
            mss_bound_value(0)

    def test_scope_separates_analytic_mss_from_annular_bar_input(self):
        """The bar complex supplies periodicity data, not MSS boundedness."""
        scope = mss_annular_bar_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('beta', 'gamma', 'delta')
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'hypVirTorusBlock' in scope['hypotheses']
        assert 'hypModularCardy' in scope['hypotheses']
        assert 'MSS normalized thermal OTOC strip hypotheses' in scope['hypotheses']
        assert 'boundedness on the half-strip' in scope['hypotheses']
        assert 'HHLL identity-block dominance for saturation' in scope['hypotheses']
        assert scope['object'] == 'normalized physical thermal OTOC strip bound'
        assert scope['analytic_theorem'] == 'Maldacena-Shenker-Stanford'
        assert scope['bound'] == 'lambda_L <= 2*pi/beta'
        assert scope['identity_block_scale'] == 'exp(2*pi*t/beta)'
        assert scope['bar_complex_role'] == (
            'genus-1 annular curvature and wrap-around periodicity datum'
        )
        assert scope['bar_curvature_alone_proves_bound'] is False
        assert scope['boundedness_or_positivity_from_bar_complex'] is False
        assert scope['thermal_coefficients_computed_by_bar_complex'] is False
        assert scope['wraparound_mode_alone_equals_lyapunov'] is False
        assert scope['full_normalized_otoc_from_annular_bar_alone'] is False
        assert scope['identity_block_saturation_requires_HHLL'] is True
        assert scope['identity_block_saturation_requires_vacuum_dominance'] is True

    def test_mss_theorem_is_analytic_and_annular_scoped(self):
        """The theorem must not sell annular curvature as a proof of MSS."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[MSS analytic strip bound with annular-bar input"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\beta+\gamma+\delta$",
            r"\hypAmbientWtCpl+\hypVirTorusBlock+\hypModularCardy",
            r"normalized thermal OTOC strip hypotheses",
            r"Maldacena--Shenker--Stanford hypotheses",
            r"boundedness on the half-strip",
            r"thermal analyticity of width",
            r"\lambda_L \;\le\; \frac{2\pi}{\beta}",
            r"does not supply boundedness on the half-strip",
            r"does not identify the wrap-around mode with the Lyapunov exponent",
            r"HHLL identity-block dominance",
            r"e^{2\pi t/\beta}",
            r"pre-scrambling range",
            r"controlled by the identity-block exchange",
            r"not a statement about an arbitrary analytic function",
            r"does not derive the MSS hypotheses from the bar complex",
        )
        for needle in required:
            assert " ".join(needle.split()) in flat

        table_start = source.index(r"\subsubsection*{Summary: the chaos dictionary}")
        table_end = source.index(r"%====================================================================", table_start)
        table = " ".join(source[table_start:table_end].split())
        table_required = (
            r"MSS Lyapunov bound \(\lambda_L\le2\pi/\beta\)",
            r"MSS strip theorem plus annular periodicity datum",
            r"The bar complex records the algebraic inputs",
        )
        for needle in table_required:
            assert " ".join(needle.split()) in table

        retired = (
            r"\begin{theorem}[MSS bound and annular bar curvature]",
            r"The analytic bound is the MSS theorem.",
            r"a function analytic in a strip of width~$\beta$",
            r"must satisfy $\lambda \le 2\pi/\beta$",
            r"The number \(2\pi/\beta\) is the fundamental wrap-around mode",
            r"the inequality is the analytic strip theorem",
            r"contains the complete quantum chaotic content of the theory",
            r"Lyapunov $\lambda_L = 2\pi/\beta$",
            r"KMS periodicity of $d_{\mathrm{wrap}}$",
        )
        for needle in retired:
            assert " ".join(needle.split()) not in flat
            assert " ".join(needle.split()) not in table


class TestScramblingTimeScope:
    """Guard the scrambling scale as conditional HHLL physics."""

    def test_scrambling_time_keeps_probe_amplitude(self):
        """The threshold solves A*c^(-1)*exp(2*pi*t/beta)=1."""
        assert simplify(
            scrambling_time_from_amplitude(beta=2, c=100, amplitude=4)
            - log(25) / pi
        ) == 0
        assert simplify(
            scrambling_time_from_amplitude(beta=2 * pi, c=64, amplitude=1)
            - log(64)
        ) == 0
        with pytest.raises(ValueError):
            scrambling_time_from_amplitude(beta=0, c=10, amplitude=1)
        with pytest.raises(ValueError):
            scrambling_time_from_amplitude(beta=1, c=10, amplitude=0)

    def test_scope_separates_shadow_truncation_from_full_scrambling(self):
        """The shadow tower gives a branch scale, not full thermal chaos."""
        scope = scrambling_time_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('gamma', 'epsilon', 'delta')
        assert 'hypAmbientWtCpl' in scope['hypotheses']
        assert 'effScalarShadowProj' in scope['hypotheses']
        assert 'hypVirTorusBlock' in scope['hypotheses']
        assert 'hypModularCardy' in scope['hypotheses']
        assert 'HHLL identity-block dominance' in scope['hypotheses']
        assert 'nonzero O(1) probe normalization' in scope['hypotheses']
        assert scope['object'] == 'HHLL identity-block scrambling scale'
        assert scope['expansion_parameter'] == (
            'u(t)=A_id*c^(-1)*exp(2*pi*t/beta)'
        )
        assert scope['threshold_equation'] == 'u(t_*) asymp 1'
        assert scope['scale_with_amplitude'] == (
            't_*=(beta/(2*pi))*log(c/A_id)+O_beta(1)'
        )
        assert scope['leading_scale'] == 't_*=(beta/(2*pi))*log(c)+O_beta(1)'
        assert scope['exact_log_c_without_amplitude_asserted'] is False
        assert scope['shadow_tower_alone_proves_scrambling'] is False
        assert scope['full_thermal_otoc_from_scalar_shadow_asserted'] is False
        assert scope['nth_physical_otoc_correction_from_Sr_asserted'] is False
        assert scope['all_fm4_strata_equal_at_tstar_asserted'] is False
        assert scope['post_tstar_interior_dominance_asserted'] is False
        assert scope['fully_scrambled_from_bar_complex_alone'] is False
        assert scope['finite_degree_failure_scale_on_selected_branch'] is True

    def test_scrambling_proposition_is_conditional_and_amplitude_scoped(self):
        """The proposition must not state exact log(c) scrambling from shadows alone."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Conditional scrambling scale on the HHLL"
        )
        prop_end = source.index(r"\end{proposition}", start)
        proof_end = source.index(r"\end{proof}", prop_end)
        block = source[start:proof_end]
        remark_start = source.index(
            r"\begin{remark}[The conditional bar-complex interpretation", proof_end
        )
        remark_end = source.index(r"\end{remark}", remark_start)
        remark = source[remark_start:remark_end]
        table_start = source.index(r"\subsubsection*{Summary: the chaos dictionary}")
        table_end = source.index(r"%====================================================================", table_start)
        table = source[table_start:table_end]
        combined = " ".join((block + remark + table).split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\gamma+\varepsilon+\delta$",
            r"\hypAmbientWtCpl+\effScalarShadowProj+\hypVirTorusBlock+\hypModularCardy",
            r"A_{\mathrm{id}}\in\mathbb R_{>0}",
            r"F_{\mathrm{OTOC}}(t)=1-u(t)+O(u(t)^2)",
            r"\frac{\beta}{2\pi}\log\frac{c}{A_{\mathrm{id}}}+O_\beta(1)",
            r"\frac{\beta}{2\pi}\log c+O_\beta(1)",
            r"without the additive normalization constant is not asserted",
            r"does not determine the full thermal OTOC expansion",
            r"S_r=\frac{8(-6)^{r-4}}{r}c^{2-r}+O(c^{1-r})",
            r"does not identify the scalar shadow tower with the full thermal OTOC expansion",
            r"does not claim that the bar complex alone proves scrambling",
            r"not a proof that every",
            r"Fulton--MacPherson stratum contributes with equal weight",
            r"not a consequence of the ordered bar complex alone",
            r"Scrambling scale \(t_*=(\beta/2\pi)\log c+O_\beta(1)\)",
            r"HHLL identity block plus scalar shadow truncation scale",
        )
        for needle in required:
            assert " ".join(needle.split()) in combined

        retired = (
            r"\begin{proposition}[Scrambling time from the $1/c$ expansion]",
            r"The scrambling time of the Virasoro system at large~$c$ is",
            r"t_* \;=\; \frac{\beta}{2\pi}\,\log c,",
            r"The $n$-th correction to the four-point function is $O(c^{-n})",
            r"the accumulated perturbative series $\sum_{n=1}^N c^{-n}\,e^{n\lambda_L t}$",
            r"Before $t_*$: the bar differential is well approximated by $m_2$",
            r"At $t = t_*$: all boundary strata",
            r"For $t > t_*$: the interior of $\FM_4(\C)$ dominates",
            r"The system is fully scrambled",
            r"bar complex requires all operations simultaneously",
            r"The scrambling time is therefore the \emph{degree transition",
            r"Scrambling $t_* = (\beta/2\pi)\log c$",
            r"Shadow obstruction tower breakdown: $S_r \sim c^{2-r}$",
        )
        for needle in retired:
            assert " ".join(needle.split()) not in combined


class TestW3WLineClosedFormScope:
    """Guard the W3 W-line scalar closed form and complementarity scope."""

    def test_w3_w_line_formula_matches_metric_branch_through_degree_14(self):
        """The W3 formula is the scalar W-line branch with sqrt(Q_W(0))=2c/3."""
        c = Symbol('c')
        t = Symbol('t')
        delta = S(122880) / (c**2 * (5 * c + 22)**3)
        h_series = series(
            t**2 * (2 * c / 3) * sqrt(1 + delta * t**2),
            t,
            0,
            15,
        ).removeO()
        for arity in range(2, 15):
            expected = simplify(expand(h_series).coeff(t, arity) / arity)
            actual = w3_w_line_shadow_coefficient(arity, c=c)
            assert simplify(actual - expected) == 0

    def test_w3_w_line_special_values_and_pole_cleared_constants(self):
        """S4, S6, and pole-cleared constants fix the normalization."""
        c = Symbol('c')
        assert simplify(w3_w_line_shadow_coefficient(2, c=c) - c / 3) == 0
        assert w3_w_line_shadow_coefficient(3, c=c) == 0
        assert simplify(
            w3_w_line_shadow_coefficient(4, c=c)
            - S(10240) / (c * (5 * c + 22)**3)
        ) == 0
        assert simplify(
            w3_w_line_shadow_coefficient(6, c=c)
            + S(209715200) / (c**3 * (5 * c + 22)**6)
        ) == 0

        for arity in (4, 6, 8, 10):
            r = arity // 2
            coeff = w3_w_line_shadow_coefficient(arity, c=c)
            cleared = simplify(
                coeff * c**(2 * r - 3) * (5 * c + 22)**(3 * (r - 1))
            )
            assert c not in cleared.free_symbols

    def test_w3_raw_complementarity_is_not_asserted(self):
        """The raw rational coefficient is not invariant under c -> 100-c."""
        c = Symbol('c')
        s4 = w3_w_line_shadow_coefficient(4, c=c)
        raw_difference = simplify(s4 - s4.subs(c, 100 - c))
        assert raw_difference != 0
        scope = w3_w_line_closed_form_scope()
        assert scope['raw_complementarity_for_S'] is False
        assert scope['pole_cleared_normalization_is_constant'] is True

    def test_w3_scope_is_scalar_w_line_branch(self):
        """The helper records W-line, projection, and excluded charges."""
        scope = w3_w_line_closed_form_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('alpha', 'gamma', 'epsilon')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'effScalarShadowProj')
        assert scope['object'] == 'normalized scalar W-line shadow coefficient'
        assert scope['line_branch'] == 'x_T=0, x_W != 0'
        assert scope['branch'] == 'sqrt(Q_W(0))=2c/3'
        assert scope['excluded_central_charges'] == (S.Zero, -S(22) / 5)
        assert scope['odd_coefficients_vanish'] is True
        assert scope['raw_two_variable_shadow_tensor_formula'] is False

    def test_w3_w_line_theorem_is_scoped_and_corrected(self):
        """The theorem must not preserve the old false binomial/complementarity text."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[$\cW_3$ $W$-line normalized scalar shadow coefficients;"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        remark_start = source.index(r"\begin{remark}[The number $30720$]", proof_end)
        remark_end = source.index(r"\end{remark}", remark_start)
        remark = source[remark_start:remark_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\gamma+\varepsilon$",
            r"\hypAmbientWtCpl+\effScalarShadowProj",
            r"\sqrt{Q_W(0)}=2c/3",
            r"c\ne0,-22/5",
            r"S_2^W=c/3",
            r"S_{2r+1}^W=0",
            r"not a formula for the full two-variable",
            r"\binom{1/2}{n}=(-1)^{n-1}\frac{C_{n-1}}{2^{2n-1}}",
            r"C_{r-1}=\frac{2(2r-3)}{r}C_{r-2}",
            r"N_{2r}^W",
            "The raw coefficient is not invariant",
        )
        for needle in required:
            assert needle in block

        remark_required = (
            "pole-cleared numerator",
            "not an additional central-charge dependence",
            r"\delta_3/4",
        )
        for needle in remark_required:
            assert needle in remark

        retired = (
            r"\begin{theorem}[$\cW_3$ $W$-line shadow coefficients]",
            r"C_{j-1}/(4^j \cdot j)",
            r"Verified through degree~$14$",
            r"leaves the normalised",
            r"$c$-independent",
            r"c \to 100 - c$ sends $S_{2r}^W(c)",
        )
        for needle in retired:
            assert needle not in block
            assert needle not in remark

    def test_w3_compute_copies_use_correct_w_line_normalization(self):
        """W3 scripts must not preserve the old c/3 square-root normalization."""
        root = Path(__file__).resolve().parents[2]
        texts = [
            (root / "compute" / "w3_shadow_coefficients.py").read_text(),
            (root / "compute" / "w3_shadow_closed_form.py").read_text(),
            (root / "compute" / "w3_multichannel_shadow.py").read_text(),
        ]
        retired = (
            "Complementarity: S_{2r}^W(c) + S_{2r}^W(100-c) is constant in c",
            "/ (3 * r * (2r-3) * c^{2r-3}",
            "Q^{WW}(u) = (c/3)²",
            "G^{WW}(u) = (c/3)√",
            "γ = 61440/(c²(5c+22)³). Complementarity",
            "S_2^{WW} = c/6",
        )
        for text in texts:
            for needle in retired:
                assert needle not in text
        assert "Q^{WW}(u) = (2c/3)^2" in texts[2]
        assert "delta = 122880/(c²(5c+22)³)" in texts[2]


class TestStasheffScalarPoleCancellationScope:
    """Guard Stasheff cancellation as scalar shadow pole cancellation."""

    def test_shadow_metric_pole_profile_ignores_rational_scalars(self):
        """The pole divisor is in Q[c]; integer denominators are separate."""
        c = Symbol('c')
        profile = virasoro_shadow_metric_pole_profile(max_r=12)
        for r, row in profile.items():
            expected = c**(r - 3) * (5 * c + 22)**((r - 2) // 2)
            assert simplify(row['pole_divisor'] - expected) == 0
            assert row['pole_divisor_matches'] is True
            assert row['rational_scalar_denominator'].is_number
        assert profile[6]['rational_scalar_denominator'] == 3
        assert profile[7]['rational_scalar_denominator'] == 7
        assert profile[11]['rational_scalar_denominator'] == 11

    def test_scope_is_scalar_projection_not_full_mr_cancellation(self):
        """Raw HPL/Kac denominators are not erased before scalar projection."""
        scope = stasheff_scalar_pole_cancellation_scope()
        assert scope['claim_status'] == 'Conditional'
        assert scope['licensing_tags'] == ('gamma', 'epsilon')
        assert scope['hypotheses'] == ('hypAmbientWtCpl', 'effScalarShadowProj')
        assert scope['object'] == 'normalized scalar shadow pole divisor'
        assert scope['nonconstant_pole_divisor'] == (
            'c^(r-3)(5c+22)^floor((r-2)/2)'
        )
        assert scope['rational_scalar_denominators_ignored'] is True
        assert scope['raw_hpl_summands_may_have_kac_poles'] is True
        assert scope['full_transferred_mr_pole_cancellation_asserted'] is False
        assert scope['kac_table_poles_cancel_after_scalar_projection'] is True
        assert scope['planarity_literal_nonplanar_tree_claim'] is False
        assert scope['operator_to_scalar_projection_constructed_all_arities'] is False

    def test_stasheff_cancellation_theorem_is_scalar_projection_scoped(self):
        """The theorem must not claim cancellation in the full A-infinity operation."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "3d_gravity.tex").read_text()
        start = source.index(
            r"\begin{theorem}[Scalar shadow pole cancellation"
        )
        theorem_end = source.index(r"\end{theorem}", start)
        proof_end = source.index(r"\end{proof}", theorem_end)
        block = source[start:proof_end]
        remark_start = source.index(r"\begin{remark}[Projected Stasheff cancellation", proof_end)
        remark_end = source.index(r"\end{remark}", remark_start)
        remark = source[remark_start:remark_end]
        prop_start = source.index(r"\begin{proposition}[Pole-divisor formula", remark_end)
        prop_end = source.index(r"\end{proof}", prop_start)
        prop = source[prop_start:prop_end]

        required = (
            r"\ClaimStatusConditional",
            r"licensing",
            r"\hypAmbientWtCpl+\effScalarShadowProj",
            "normalized scalar shadow projection",
            "nonconstant pole divisor",
            r"c^{r-3}(5c+22)^{\lfloor(r-2)/2\rfloor}",
            "Rational scalar denominators",
            "raw HPL summands and raw transferred",
            "operations may contain",
            "does not construct the all-arity",
            r"\Pi_{\mathrm{sh}}^{(r)}",
        )
        for needle in required:
            assert needle in block

        remark_required = (
            "projected singularities it excludes",
            "not a literal statement that every tree amplitude",
            "with a\nKac-table denominator is non-planar",
            "after scalar shadow projection",
        )
        for needle in remark_required:
            assert needle in remark

        prop_required = (
            r"\begin{proposition}[Pole-divisor formula",
            "nonconstant pole divisor",
            "Rational scalar factors are not included",
            "No Kac-table central-charge poles appear in the scalar shadow",
        )
        for needle in prop_required:
            assert needle in prop

        retired = (
            r"\begin{theorem}[Stasheff cancellation]",
            r"\textbf{cancel in the full $\Ainf$ computation}",
            "cancel in the \\emph{total} $m_k$",
            "appears in the individual HPL tree\namplitudes, yet cancels in their sum",
            "every tree amplitude whose\ndenominator contains a Kac-table factor is a non-planar",
            r"\begin{proposition}[Denominator formula]",
            r"\operatorname{denom}(S_r)",
        )
        for needle in retired:
            assert needle not in block
            assert needle not in remark
            assert needle not in prop


# ===================================================================
# 10. CROSS-ENGINE CONSISTENCY
# ===================================================================

class TestCrossEngineConsistency:
    """Cross-engine checks (AP10 compliance)."""

    def test_kappa_ground_truth_values(self):
        """κ values match known ground truth."""
        ground_truth = {
            0: Rational(0),
            1: Rational(1, 2),
            13: Rational(13, 2),
            26: Rational(13),
        }
        for c_val, expected_kappa in ground_truth.items():
            computed = gravity_kappa(c=c_val)
            assert simplify(computed - expected_kappa) == 0, \
                f"κ mismatch at c={c_val}: got {computed}, expected {expected_kappa}"

    def test_shadow_depth_infinite(self):
        """Virasoro shadow depth is infinite (class M)."""
        assert virasoro_shadow_depth() == float('inf')

    def test_shadow_class_M(self):
        """Virasoro is class M (mixed, infinite tower)."""
        assert virasoro_shadow_class() == 'M'

    def test_kappa_matches_bbd_engine(self):
        """Cross-check: κ from gravity engine matches bulk_boundary_duality_engine."""
        try:
            from lib.bulk_boundary_duality_engine import virasoro_koszul_pair
            pair = virasoro_koszul_pair()
            bbd_kappa = pair.kappa
            gravity_kappa_val = gravity_kappa(c=pair.central_charge)
            assert simplify(bbd_kappa - gravity_kappa_val) == 0
        except (ImportError, AttributeError):
            pytest.skip("bulk_boundary_duality_engine not available for cross-check")

    def test_quartic_contact_cross_engine(self):
        """Cross-check: Q^contact from gravity engine vs Vol I engine if available."""
        try:
            sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                            '..', '..', '..', 'chiral-bar-cobar', 'compute'))
            from lib.virasoro_shadow_all_arity import quartic_contact as vol1_quartic
            c_val = Rational(1)
            vol1_val = vol1_quartic(c_val)
            vol2_val = quartic_contact_virasoro(c=1)
            assert simplify(vol1_val - vol2_val) == 0
        except (ImportError, AttributeError):
            pytest.skip("Vol I virasoro_shadow_all_arity not available")


# ===================================================================
# 11. BTZ AND BROWN-HENNEAUX
# ===================================================================

class TestBTZEntropy:
    """Verify BTZ entropy formula from κ."""

    def test_btz_argument(self):
        """S_BTZ ~ √(2κh/3). Verify the argument 2κh/3."""
        # For c=6, h=1: 2*(6/2)*1/3 = 2
        kappa = gravity_kappa(c=6)
        result = 2 * kappa * 1 / 3
        assert result == 2

    def test_btz_at_comparison_fixed_point(self):
        """At c=13, h=1: argument = 2*(13/2)*1/3 = 13/3."""
        kappa = gravity_kappa(c=13)
        result = 2 * kappa * 1 / 3
        assert result == Rational(13, 3)


# ===================================================================
# 12. INTEGRATION / SMOKE TESTS
# ===================================================================

class TestSmoke:
    """Quick smoke tests for overall consistency."""

    def test_all_functions_callable(self):
        """All exported functions are callable with default arguments."""
        virasoro_lambda_bracket()
        virasoro_associator()
        virasoro_m3_coefficients()
        quartic_contact_virasoro()
        gravity_kappa()
        koszul_dual_central_charge()
        gravity_laplace_kernel_poles()
        gravity_r_matrix_poles()
        virasoro_exact_gravity_scope_profile()
        virasoro_scalar_bar_trace_profile()
        virasoro_primary_ward_connection()
        virasoro_primary_ward_log()
        virasoro_primary_ward_even_coefficients()
        virasoro_primary_ward_even_exponents()
        complementarity_constant_virasoro()
        verify_complementarity()

    def test_c13_fixed_point_consistency(self):
        """At c=13 fixed point: κ=13/2, c_dual=13, Q=10/(13*87)."""
        assert gravity_kappa(c=13) == Rational(13, 2)
        assert koszul_dual_central_charge(c=13) == 13
        assert quartic_contact_virasoro(c=13) == Rational(10, 1131)
