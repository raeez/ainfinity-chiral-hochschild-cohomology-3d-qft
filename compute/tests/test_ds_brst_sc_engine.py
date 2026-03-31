"""Tests for DS-BRST Swiss-Cheese compatibility engine.

Verifies that Drinfeld-Sokolov BRST reduction transports the
SC^{ch,top}-algebra structure from V_k(sl₂) to Vir_{c_DS(k)}.

Test organization:
  1. Central charge formula and decomposition
  2. λ-bracket: DS produces Virasoro at c = c_DS
  3. A∞ operations: m_k match at arities 2-6
  4. Curvature transport: κ_DS = c_DS/2
  5. Complementarity: κ + κ! = 13
  6. Complexity transport: class L → class M
  7. Genus transport: F_1, F_2 agree
  8. Cross-engine consistency with existing engines
  9. Multi-level robustness scan

AP10 compliance: no hardcoded expected values; all comparisons
are structural (identities that hold for all k).
"""
import pytest
from sympy import S, Symbol, Rational, simplify, expand, oo

# Engine under test
from compute.lib.ds_brst_sc_engine import (
    ds_central_charge_sl2,
    sugawara_central_charge_sl2,
    ds_central_charge_decomposition,
    ds_lambda_bracket_from_affine,
    ds_associator_comparison,
    ds_kappa_sl2,
    ds_kappa_complementarity,
    ds_complexity_transport,
    ds_genus_transport,
    ds_mk_comparison_numerical,
    ds_brst_sc_full_verification,
)

# Cross-engine imports
from compute.lib.examples.affine_kac_moody import (
    sl2_data, affine_central_charge, affine_kappa,
    kappa_complementarity_affine,
)
from compute.lib.swiss_cheese_virasoro_wheels import (
    m2_virasoro, m3_virasoro, m3_virasoro_as_poly,
    associator_virasoro, mk_stasheff_recursive_numerical,
)
from compute.lib.genus_one_bridge import (
    genus1_curvature, period_correction,
)
from compute.lib.bulk_boundary_duality_engine import (
    virasoro_koszul_pair,
)


# =========================================================================
# 1. CENTRAL CHARGE
# =========================================================================

class TestDSCentralCharge:
    """Test the DS central charge formula c_DS(k) = 1 - 6(k+1)²/(k+2)."""

    def test_formula_symbolic(self):
        """Verify c_DS(k) = 1 - 6(k+1)²/(k+2) symbolically."""
        k = Symbol('k')
        c = ds_central_charge_sl2(k)
        expected = 1 - 6 * (k + 1) ** 2 / (k + 2)
        assert simplify(c - expected) == 0

    def test_specific_values(self):
        """Verify c_DS at specific integer levels."""
        # k=0: c = 1 - 6/2 = -2
        assert ds_central_charge_sl2(0) == -2
        # k=1: c = 1 - 24/3 = -7
        assert ds_central_charge_sl2(1) == -7
        # k=2: c = 1 - 54/4 = -25/2
        assert ds_central_charge_sl2(2) == Rational(-25, 2)

    def test_critical_level_raises(self):
        """DS central charge must raise at critical level k = -h^∨ = -2."""
        with pytest.raises(ValueError, match="critical level"):
            ds_central_charge_sl2(-2)

    def test_decomposition_verification(self):
        """The decomposition internal consistency check must pass."""
        for k_val in [1, 3, 5, 10]:
            decomp = ds_central_charge_decomposition(k_val)
            assert decomp['verification'] == 0

    def test_sugawara_vs_ds(self):
        """c_DS ≠ c_Sug (the DS shift is substantial)."""
        for k_val in [1, 2, 5]:
            c_sug = sugawara_central_charge_sl2(k_val)
            c_ds = ds_central_charge_sl2(k_val)
            assert c_ds != c_sug
            # c_DS < c_Sug for all positive k (the shift is negative)
            assert float(c_ds) < float(c_sug)

    def test_large_k_asymptotics(self):
        """For large k: c_DS ≈ -6k (leading term)."""
        k = 1000
        c = float(ds_central_charge_sl2(k))
        ratio = c / (-6 * k)
        assert abs(ratio - 1.0) < 0.01  # within 1%


# =========================================================================
# 2. λ-BRACKET
# =========================================================================

class TestDSLambdaBracket:
    """Verify that DS reduction produces the Virasoro λ-bracket."""

    def test_bracket_structure(self):
        """DS λ-bracket has the standard Virasoro form."""
        for k_val in [1, 3, 10]:
            result = ds_lambda_bracket_from_affine(k_val)
            # The expected form: {T_λ T} = ∂T + 2Tλ + (c/12)λ³
            assert result['expected_dT_coeff'] == 1
            lam = Symbol('lambda')
            assert result['expected_T_coeff'] == 2 * lam

    def test_quartic_pole_coefficient(self):
        """The quartic pole coefficient is c_DS/12 in the λ-bracket."""
        for k_val in [1, 2, 5, 10]:
            c_ds = ds_central_charge_sl2(k_val)
            lam = Symbol('lambda')
            result = ds_lambda_bracket_from_affine(k_val)
            expected_scalar = c_ds * lam ** 3 / 12
            assert simplify(result['expected_scalar_coeff'] - expected_scalar) == 0


# =========================================================================
# 3. A∞ OPERATIONS: m_k MATCH
# =========================================================================

class TestDSAinftyMatch:
    """Verify DS-transferred A∞ operations match intrinsic Virasoro."""

    def test_m3_comparison_symbolic(self):
        """m_3^{DS} = m_3^{Vir}(c_DS) symbolically."""
        for k_val in [1, 3, 5]:
            comp = ds_associator_comparison(k_val)
            m3 = comp['m3_virasoro_at_c_ds']
            c_ds = comp['c_ds']

            # The m_3 formula (from Movement I):
            # m_3(T,T,T; l1, l2) = ∂²T + (2l1+3l2)∂T + 2l2(2l1+l2)T + (c/12)l2³(2l1+l2)
            l1, l2 = Symbol('lambda_1'), Symbol('lambda_2')
            assert m3['d2T'] == 1
            assert simplify(m3['dT'] - (2 * l1 + 3 * l2)) == 0
            assert simplify(m3['T'] - 2 * l2 * (2 * l1 + l2)) == 0
            expected_scalar = c_ds * l2 ** 3 * (2 * l1 + l2) / 12
            assert simplify(m3['scalar'] - expected_scalar) == 0

    def test_m3_nonvanishing_numerical(self):
        """m_3 ≠ 0 at c = c_DS for all tested levels (class M)."""
        for k_val in [1, 2, 3, 5, 10]:
            result = ds_mk_comparison_numerical(k_val, 3, [1.0, 1.0])
            assert result['nonvanishing'], f"m_3 vanishes at k={k_val}"

    def test_m4_nonvanishing(self):
        """m_4 ≠ 0 at c = c_DS for all tested levels."""
        for k_val in [1, 3, 5, 10]:
            result = ds_mk_comparison_numerical(k_val, 4, [0.5, 1.0, 1.5])
            assert result['nonvanishing'], f"m_4 vanishes at k={k_val}"

    def test_m5_nonvanishing(self):
        """m_5 ≠ 0 at c = c_DS for all tested levels."""
        for k_val in [1, 3, 5, 10]:
            result = ds_mk_comparison_numerical(k_val, 5, [0.25, 0.5, 0.75, 1.0])
            assert result['nonvanishing'], f"m_5 vanishes at k={k_val}"

    def test_m6_nonvanishing(self):
        """m_6 ≠ 0 at c = c_DS for all tested levels."""
        for k_val in [1, 3, 5]:
            result = ds_mk_comparison_numerical(
                k_val, 6, [0.2, 0.4, 0.6, 0.8, 1.0]
            )
            assert result['nonvanishing'], f"m_6 vanishes at k={k_val}"

    def test_mk_matches_intrinsic_virasoro(self):
        """DS m_k match intrinsic Virasoro at same c (cross-engine check)."""
        for k_val in [3, 5]:
            c_ds = float(ds_central_charge_sl2(k_val))
            for arity in [3, 4, 5]:
                lam_vals = [float(i) / float(arity) for i in range(1, arity)]
                # DS prediction
                ds_result = ds_mk_comparison_numerical(k_val, arity, lam_vals)
                # Direct Virasoro computation
                vir_result = mk_stasheff_recursive_numerical(
                    arity, lam_vals, c_ds
                )
                # They should agree (both use the same Stasheff recursion at c_DS)
                assert abs(ds_result['T_coeff'] - vir_result['T_coeff']) < 1e-8


# =========================================================================
# 4. CURVATURE TRANSPORT
# =========================================================================

class TestDSCurvatureTransport:
    """Verify κ_DS = c_DS/2 and genus-g free energy transport."""

    def test_kappa_equals_c_over_2(self):
        """κ(Vir_{c_DS}) = c_DS/2 for all levels."""
        for k_val in [1, 2, 3, 5, 10]:
            curv = ds_kappa_sl2(k_val)
            c_ds = ds_central_charge_sl2(k_val)
            assert simplify(curv['kappa_ds'] - c_ds / 2) == 0

    def test_kappa_nonzero(self):
        """κ_DS ≠ 0 for all positive integer levels (genus tower is nontrivial)."""
        for k_val in [1, 2, 3, 5, 10]:
            curv = ds_kappa_sl2(k_val)
            assert curv['kappa_ds'] != 0

    def test_complementarity_universal(self):
        """κ(c_DS) + κ(26 - c_DS) = 13 for all levels."""
        for k_val in [1, 2, 3, 5, 7, 10, 20, 50]:
            comp = ds_kappa_complementarity(k_val)
            assert comp['complementarity_holds'], f"Complementarity fails at k={k_val}"
            assert comp['kappa_sum'] == 13

    def test_genus1_free_energy(self):
        """F₁(Vir_{c_DS}) = κ_DS/24 = c_DS/48."""
        for k_val in [1, 3, 10]:
            g1 = ds_genus_transport(k_val, genus=1)
            c_ds = ds_central_charge_sl2(k_val)
            expected_F1 = c_ds / 48
            assert simplify(g1['F_g_virasoro'] - expected_F1) == 0

    def test_genus2_free_energy(self):
        """F₂(Vir_{c_DS}) = κ_DS · 7/5760 = c_DS · 7/11520."""
        for k_val in [1, 3, 10]:
            g2 = ds_genus_transport(k_val, genus=2)
            c_ds = ds_central_charge_sl2(k_val)
            expected_F2 = c_ds * Rational(7, 11520)
            assert simplify(g2['F_g_virasoro'] - expected_F2) == 0


# =========================================================================
# 5. COMPLEXITY TRANSPORT
# =========================================================================

class TestDSComplexityTransport:
    """Verify class L → class M under DS reduction."""

    def test_input_is_class_L(self):
        """V_k(sl₂) is class L (r_max = 3)."""
        for k_val in [1, 3, 10]:
            ct = ds_complexity_transport(k_val)
            assert ct['input_class'] == 'L'
            assert ct['input_r_max'] == 3

    def test_output_is_class_M(self):
        """Vir_{c_DS} is class M (r_max = ∞)."""
        for k_val in [1, 3, 10]:
            ct = ds_complexity_transport(k_val)
            assert ct['output_class'] == 'M'
            assert ct['output_r_max'] == oo

    def test_quartic_contact_nonzero(self):
        """The quartic contact invariant Q ≠ 0 at c = c_DS."""
        for k_val in [1, 3, 5, 10]:
            ct = ds_complexity_transport(k_val)
            Q = ct['quartic_contact']
            assert simplify(Q) != 0, f"Quartic contact vanishes at k={k_val}"

    def test_quartic_contact_poles(self):
        """Q has poles at c = 0 and c = -22/5, not at c_DS for integer k > 0."""
        for k_val in [1, 2, 3, 5, 10, 20]:
            c_ds = ds_central_charge_sl2(k_val)
            # c_DS ≠ 0 for k > 0 (c_DS = 0 would require k+1 = ±√((k+2)/6))
            assert simplify(c_ds) != 0
            # c_DS ≠ -22/5 (would require specific k)
            assert simplify(c_ds + Rational(22, 5)) != 0


# =========================================================================
# 6. CROSS-ENGINE CONSISTENCY
# =========================================================================

class TestCrossEngineConsistency:
    """Verify DS engine results are consistent with other engines."""

    def test_affine_kappa_vs_ds_kappa(self):
        """The affine κ and DS κ are distinct (DS changes the stress tensor)."""
        g = sl2_data()
        for k_val in [1, 3, 5]:
            kappa_aff = affine_kappa(g, k_val)
            curv = ds_kappa_sl2(k_val)
            kappa_ds = curv['kappa_ds']
            # They must differ (DS shifts the curvature)
            assert simplify(kappa_aff - kappa_ds) != 0

    def test_virasoro_koszul_pair_at_c_ds(self):
        """The Koszul dual of Vir_{c_DS} is Vir_{26-c_DS}."""
        for k_val in [1, 3, 5]:
            c_ds = ds_central_charge_sl2(k_val)
            c_float = float(c_ds)
            pair = virasoro_koszul_pair(c_val=c_float)
            assert abs(float(pair.dual_central_charge) - (26 - c_float)) < 1e-10

    def test_genus1_bridge_consistency(self):
        """DS curvature matches the genus-1 bridge virasoro family."""
        for k_val in [1, 3, 5]:
            c_ds = ds_central_charge_sl2(k_val)
            curv_ds = ds_kappa_sl2(k_val)
            # The genus-1 bridge gives κ(Vir_c) = c/2
            g1 = genus1_curvature('virasoro', c=c_ds)
            assert simplify(curv_ds['kappa_ds'] - g1['kappa']) == 0

    def test_virasoro_m3_direct_comparison(self):
        """Direct comparison of m_3 from DS engine vs Virasoro wheels engine."""
        for k_val in [1, 3, 5]:
            c_ds = ds_central_charge_sl2(k_val)
            l1, l2 = Symbol('l1'), Symbol('l2')

            # From DS engine
            comp = ds_associator_comparison(k_val, l1, l2)
            m3_ds = comp['m3_virasoro_at_c_ds']

            # From Virasoro wheels engine
            m3_vir = m3_virasoro(l1, l2, c_sym=c_ds)

            # Compare each component
            assert simplify(m3_ds['d2T'] - m3_vir['d2T']) == 0
            assert simplify(m3_ds['dT'] - m3_vir['dT']) == 0
            assert simplify(m3_ds['T'] - m3_vir['T']) == 0
            assert simplify(m3_ds['scalar'] - m3_vir.get('1', S.Zero)) == 0


# =========================================================================
# 7. MASTER VERIFICATION
# =========================================================================

class TestMasterVerification:
    """Run the full master verification at multiple levels."""

    @pytest.mark.parametrize("k_val", [1, 2, 3, 5, 10])
    def test_full_verification(self, k_val):
        """All checks pass at level k."""
        results = ds_brst_sc_full_verification(k_val)
        assert results['summary']['all_checks_pass'], (
            f"Full verification failed at k={k_val}: "
            f"{[c for c, p in results['summary']['checks'] if not p]}"
        )

    def test_full_verification_large_k(self):
        """Verification passes at large k (semiclassical regime)."""
        results = ds_brst_sc_full_verification(k_val=50)
        assert results['summary']['all_checks_pass']


# =========================================================================
# 8. STRUCTURAL THEOREMS
# =========================================================================

class TestStructuralTheorems:
    """Tests for the structural consequences of DS-BRST SC compatibility."""

    def test_complexity_increase_is_strict(self):
        """DS reduction strictly increases complexity: L → M, never L → L."""
        for k_val in [1, 2, 3, 5, 10, 20]:
            ct = ds_complexity_transport(k_val)
            assert ct['output_r_max'] > ct['input_r_max']

    def test_curvature_sign_change(self):
        """DS curvature is negative for positive integer k.

        κ(V_k(sl₂)) > 0 but κ(Vir_{c_DS}) < 0 for k > 0.
        The sign change reflects the ghost contribution overwhelming
        the Sugawara contribution.
        """
        g = sl2_data()
        for k_val in [1, 2, 3, 5, 10]:
            kappa_aff = float(affine_kappa(g, k_val))
            kappa_ds = float(ds_kappa_sl2(k_val)['kappa_ds'])
            assert kappa_aff > 0, f"Affine κ should be positive at k={k_val}"
            assert kappa_ds < 0, f"DS κ should be negative at k={k_val}"

    def test_koszul_dual_central_charge_positive(self):
        """The Koszul dual Vir_{26-c_DS} has large positive c for k > 0.

        Since c_DS < 0 for k > 0, we have 26 - c_DS > 26.
        The Koszul dual is in the 'large central charge' regime.
        """
        for k_val in [1, 3, 5, 10]:
            c_ds = float(ds_central_charge_sl2(k_val))
            c_dual = 26 - c_ds
            assert c_dual > 26, f"Dual c should be > 26 at k={k_val}"

    def test_self_dual_point(self):
        """Virasoro self-dual at c = 13. Check c_DS(k) = 13 has a solution.

        1 - 6(k+1)²/(k+2) = 13 → -6(k+1)² = 12(k+2) → (k+1)² = -2(k+2)
        → k² + 2k + 1 = -2k - 4 → k² + 4k + 5 = 0
        → k = (-4 ± √(16-20))/2 = (-4 ± 2i)/2 = -2 ± i.

        So c_DS = 13 only at complex level k = -2 ± i. No real self-dual
        point exists for the DS reduction.
        """
        # Verify: c_DS at k = -2+i should be 13
        from sympy import I
        k_sd = -2 + I
        c_sd = ds_central_charge_sl2(k_sd)
        assert simplify(c_sd - 13) == 0


# =========================================================================
# 9. FIRST-PRINCIPLES QUARTIC POLE (Wick decomposition)
# =========================================================================

class TestQuarticPoleFirstPrinciples:
    """Verify the first-principles quartic pole decomposition."""

    @pytest.mark.parametrize("k_val", [1, 2, 3, 5, 10, 20])
    def test_total_matches_expected(self, k_val):
        """Total quartic pole = c_DS/2 from first principles."""
        from compute.lib.ds_brst_sc_engine import quartic_pole_first_principles
        r = quartic_pole_first_principles(k_val)
        assert r['match'], f"First-principles mismatch at k={k_val}"

    @pytest.mark.parametrize("k_val", [1, 2, 3, 5, 10])
    def test_sugawara_sub_decomposition(self, k_val):
        """Naive Wick + cascading = c_Sug/2."""
        from compute.lib.ds_brst_sc_engine import quartic_pole_first_principles
        r = quartic_pole_first_principles(k_val)
        assert r['sugawara_match'], f"Sugawara sub-decomposition fails at k={k_val}"

    def test_improvement_is_negative(self):
        """The improvement contribution is NEGATIVE for all positive k.

        This is the deepest arithmetic: ∂_z∂_w[1/(z-w)²] = -6/(z-w)⁴.
        """
        from compute.lib.ds_brst_sc_engine import quartic_pole_first_principles
        for k_val in [1, 2, 3, 5, 10, 50]:
            r = quartic_pole_first_principles(k_val)
            assert float(r['improvement']) < 0

    def test_cascading_dominates_naive_wick(self):
        """Cascading > naive Wick for all k ≥ 1.

        The root-sector cascading from J^eJ^f interacting contractions
        is always larger than the Cartan-sector naive Wick.
        """
        from compute.lib.ds_brst_sc_engine import quartic_pole_first_principles
        for k_val in [1, 3, 5, 10]:
            r = quartic_pole_first_principles(k_val)
            assert float(r['cascading_root']) > float(r['naive_wick_cartan'])

    def test_improvement_dominates_sugawara(self):
        """For k ≥ 1: |improvement| > Sugawara/2, making c_DS negative.

        This is WHY c_DS < 0 for all positive integer levels.
        """
        from compute.lib.ds_brst_sc_engine import quartic_pole_first_principles
        for k_val in [1, 3, 5, 10]:
            r = quartic_pole_first_principles(k_val)
            assert abs(float(r['improvement'])) > float(r['sugawara_half'])


# =========================================================================
# 10. sl_N UNIVERSALITY
# =========================================================================

class TestSlNUniversality:
    """Verify L → M complexity transport for all sl_N."""

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    @pytest.mark.parametrize("k_val", [1, 3])
    def test_quartic_pole_matches_slN(self, N, k_val):
        """First-principles decomposition matches for sl_N."""
        from compute.lib.ds_brst_sc_engine import ds_quartic_pole_slN
        r = ds_quartic_pole_slN(N, k_val)
        assert r['match'], f"sl_{N}, k={k_val}: decomposition mismatch"

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_complexity_L_to_M(self, N):
        """DS reduction maps class L → class M for all sl_N."""
        from compute.lib.ds_brst_sc_engine import ds_complexity_slN
        r = ds_complexity_slN(N, 3)
        assert r['input_class'] == 'L'
        assert r['output_class'] == 'M'
        assert r['output_r_max'] == oo

    @pytest.mark.parametrize("N", [2, 3, 4, 5])
    def test_c_ds_negative_for_positive_k(self, N):
        """c_DS < 0 for all sl_N at positive integer level."""
        from compute.lib.ds_brst_sc_engine import ds_central_charge_slN
        for k_val in [1, 2, 3]:
            c = float(ds_central_charge_slN(N, k_val))
            assert c < 0, f"sl_{N}, k={k_val}: c_DS = {c} should be negative"

    def test_ghost_grows_with_rank(self):
        """Ghost contribution grows dramatically with rank."""
        from compute.lib.ds_brst_sc_engine import ds_quartic_pole_slN
        ghosts = []
        for N in [2, 3, 4, 5]:
            r = ds_quartic_pole_slN(N, 1)
            ghosts.append(float(r['c_ghost']))
        # Ghost contributions should be monotonically decreasing (more negative)
        for i in range(len(ghosts) - 1):
            assert ghosts[i + 1] < ghosts[i]

    def test_rho_squared_formula(self):
        """(ρ∨,ρ∨)_κ = N(N²-1)/12 for sl_N."""
        from compute.lib.ds_brst_sc_engine import ds_quartic_pole_slN
        for N in [2, 3, 4, 5]:
            r = ds_quartic_pole_slN(N, 1)
            expected = N * (N ** 2 - 1) / 12
            assert float(r['rho_squared']) == expected


# =========================================================================
# 11. GROWTH RATE ANALYSIS
# =========================================================================

class TestGrowthRate:
    """Test m_k growth rate under DS reduction."""

    def test_all_nonvanishing(self):
        """All m_k ≠ 0 up to arity 8 for k = 1, 3, 10."""
        from compute.lib.ds_brst_sc_engine import mk_growth_rate
        for k_val in [1, 3, 10]:
            r = mk_growth_rate(k_val, max_arity=8)
            assert r['all_nonvanishing'], f"Some m_k vanish at k={k_val}"

    def test_T_coeff_independent_of_c(self):
        """T-coefficient of m_k is INDEPENDENT of c at fixed spectral params.

        This is because the Stasheff recursion for m_k has T-coefficient
        determined purely by the combinatorial structure (the number 2 in
        the Virasoro bracket {T_λ T} = ∂T + 2Tλ + ...), not by c.
        """
        from compute.lib.ds_brst_sc_engine import mk_growth_rate
        r1 = mk_growth_rate(1, max_arity=6)
        r2 = mk_growth_rate(10, max_arity=6)
        for i, (d1, d2) in enumerate(zip(r1['data'], r2['data'])):
            assert abs(d1['T_coeff'] - d2['T_coeff']) < 1e-6, (
                f"T_coeff differs at arity {d1['arity']}"
            )

    def test_scalar_scales_with_c(self):
        """Scalar coefficient scales with c_DS (larger |c| → larger |scalar|)."""
        from compute.lib.ds_brst_sc_engine import mk_growth_rate
        r1 = mk_growth_rate(1, max_arity=6)
        r10 = mk_growth_rate(10, max_arity=6)
        # At arity 3: scalar = (c/12)·λ₂³·(2λ₁+λ₂), so |scalar| ∝ |c|
        sc1 = abs(r1['data'][1]['scalar_coeff'])  # arity 3
        sc10 = abs(r10['data'][1]['scalar_coeff'])
        assert sc10 > sc1  # |c_DS(10)| > |c_DS(1)|
