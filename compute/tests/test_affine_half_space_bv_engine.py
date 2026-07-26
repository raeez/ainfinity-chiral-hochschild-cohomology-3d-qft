"""Tests for the affine half-space BV engine.

Verifies dual Coxeter numbers, effective level shifts, one-loop exactness,
two-loop vanishing, and DS-compatibility of the level shift for the
affine half-space BV chapter (affine_half_space_bv.tex).

Test organization:
  1. TestDualCoxeter: h^vee for all classical and exceptional families
  2. TestLevelShift: K_eff = k + hbar*(-h^vee) at various parameters
  3. TestOneLoopGraphCount: graph counts at n=1,2 vertices
  4. TestTwoLoopVanishing: structural vanishing argument
  5. TestDSCompatibility: level shift survives DS for sl_2->Vir, sl_3->W_3
  6. TestKappaConsistency: cross-checks between kappa and central charge
  7. TestEdgeCases: critical level, exceptional algebras, boundary cases

AP10 compliance: structural identities (additivity, DS functoriality)
used as cross-checks, not just hardcoded expected values.
"""
import pytest
from pathlib import Path
from sympy import Rational, S, simplify

from compute.lib.affine_half_space_bv_engine import (
    dual_coxeter_number,
    effective_level_shift,
    image_charge_propagator_profile,
    one_loop_graph_count,
    reflected_arnold_relation_profile,
    reflected_obstruction_profile,
    reflected_weight_composition,
    two_loop_vanishing_reason,
    verify_level_shift_ds_compatible,
    kappa_kac_moody,
    central_charge_sugawara,
)


# =========================================================================
# 1. DUAL COXETER NUMBERS
# =========================================================================

class TestDualCoxeter:
    """Verify h^vee for all simple Lie algebras against Kac tables."""

    def test_sl2(self):
        """h^vee(sl_2) = h^vee(A_1) = 2."""
        assert dual_coxeter_number('A', 1) == 2

    def test_sl3(self):
        """h^vee(sl_3) = h^vee(A_2) = 3."""
        assert dual_coxeter_number('A', 2) == 3

    def test_sl4(self):
        """h^vee(sl_4) = h^vee(A_3) = 4."""
        assert dual_coxeter_number('A', 3) == 4

    def test_so5(self):
        """h^vee(so_5) = h^vee(B_2) = 2*2-1 = 3."""
        assert dual_coxeter_number('B', 2) == 3

    def test_sp4(self):
        """h^vee(sp_4) = h^vee(C_2) = 2+1 = 3."""
        assert dual_coxeter_number('C', 2) == 3

    def test_B2_equals_C2(self):
        """B_2 ~ C_2 isomorphism: h^vee must agree."""
        assert dual_coxeter_number('B', 2) == dual_coxeter_number('C', 2)

    def test_G2(self):
        """h^vee(G_2) = 4."""
        assert dual_coxeter_number('G', 2) == 4

    def test_E8(self):
        """h^vee(E_8) = 30."""
        assert dual_coxeter_number('E', 8) == 30

    def test_F4(self):
        """h^vee(F_4) = 9."""
        assert dual_coxeter_number('F', 4) == 9

    def test_E6(self):
        """h^vee(E_6) = 12."""
        assert dual_coxeter_number('E', 6) == 12

    def test_E7(self):
        """h^vee(E_7) = 18."""
        assert dual_coxeter_number('E', 7) == 18

    def test_D4(self):
        """h^vee(D_4) = 2*4-2 = 6."""
        assert dual_coxeter_number('D', 4) == 6

    def test_unknown_type_raises(self):
        """Unknown Lie type should raise ValueError."""
        with pytest.raises(ValueError):
            dual_coxeter_number('X', 3)


# =========================================================================
# 2. EFFECTIVE LEVEL SHIFT
# =========================================================================

class TestLevelShift:
    """Verify K_eff = k + hbar*(-h^vee) = k - hbar*h^vee."""

    def test_sl2_k10_hbar1(self):
        """sl_2, k=10, hbar=1: K_eff = 10 - 1*2 = 8."""
        result = effective_level_shift(10, 1, 'A', 1)
        assert result == Rational(8)

    def test_sl2_k5_hbar1(self):
        """sl_2, k=5, hbar=1: K_eff = 5 - 2 = 3."""
        assert effective_level_shift(5, 1, 'A', 1) == Rational(3)

    def test_sl2_shifted_zero(self):
        """sl_2, k=2, hbar=1: K_eff = 2-2 = 0.

        This is the zero of the shifted effective level in the
        half-space BV convention. It is not the affine critical level,
        which is k=-h^vee=-2 for sl_2.
        """
        assert effective_level_shift(2, 1, 'A', 1) == Rational(0)

    def test_sl3_k10_hbar1(self):
        """sl_3, k=10, hbar=1: K_eff = 10 - 3 = 7."""
        assert effective_level_shift(10, 1, 'A', 2) == Rational(7)

    def test_hbar_zero_no_shift(self):
        """At hbar=0, K_eff = k (no quantum correction)."""
        assert effective_level_shift(10, 0, 'A', 1) == Rational(10)

    def test_hbar_half(self):
        """sl_2, k=10, hbar=1/2: K_eff = 10 - (1/2)*2 = 9."""
        result = effective_level_shift(10, Rational(1, 2), 'A', 1)
        assert result == Rational(9)

    def test_E8_k30_hbar1(self):
        """E_8, k=30, hbar=1: K_eff = 30 - 30 = 0.

        This is shifted-effective-level zero, not the affine critical
        level k=-h^vee=-30.
        """
        assert effective_level_shift(30, 1, 'E', 8) == Rational(0)

    def test_rational_level(self):
        """Fractional level k=5/2 for sl_2: K_eff = 5/2 - 2 = 1/2."""
        result = effective_level_shift(Rational(5, 2), 1, 'A', 1)
        assert result == Rational(1, 2)


# =========================================================================
# 3. ONE-LOOP GRAPH COUNTS
# =========================================================================

class TestOneLoopGraphCount:
    """Verify 1-loop Feynman graph counts."""

    def test_n1_tadpole(self):
        """n=1: exactly 1 graph (tadpole)."""
        assert one_loop_graph_count(1) == 1

    def test_n2_bubble(self):
        """n=2: exactly 1 graph (bubble)."""
        assert one_loop_graph_count(2) == 1

    def test_n0_raises(self):
        """n=0 is invalid (must have at least 1 vertex)."""
        with pytest.raises(ValueError):
            one_loop_graph_count(0)


# =========================================================================
# 4. TWO-LOOP VANISHING
# =========================================================================

class TestTwoLoopVanishing:
    """Verify the structural argument for 2-loop vanishing."""

    def test_reason_nonempty(self):
        """The vanishing reason should be a non-empty string."""
        reason = two_loop_vanishing_reason()
        assert isinstance(reason, str)
        assert len(reason) > 50

    def test_reason_mentions_beta_leg(self):
        """The argument should mention the beta-leg constraint."""
        reason = two_loop_vanishing_reason()
        assert 'beta' in reason.lower()

    def test_reason_mentions_figure_eight(self):
        """The argument should identify the figure-8 graph."""
        reason = two_loop_vanishing_reason()
        assert 'figure-8' in reason.lower() or 'figure_8' in reason.lower()

    def test_reason_mentions_propagator(self):
        """The argument should mention propagators."""
        reason = two_loop_vanishing_reason()
        assert 'propagator' in reason.lower()


# =========================================================================
# 4b. REFLECTED HALF-SPACE PROPAGATOR
# =========================================================================

class TestReflectedHalfSpacePropagator:
    """Verify the image-charge and reflected Arnold package."""

    def test_image_charge_profile(self):
        profile = image_charge_propagator_profile()

        assert profile["half_space"] == "H={t>=0}"
        assert "P_H(x,y)" in profile["propagator"]
        assert "-sigma P" in profile["propagator"]
        assert profile["reflection"] == "sigma(phi)=(-1)^{w_phi} phi"
        assert "OPE_boundary" in profile["boundary_ope"]
        assert "-(-1)^{w_j}" in profile["boundary_ope"]

    def test_reflected_weight_composition_mod_two(self):
        assert reflected_weight_composition(1, 0, 1) == 0
        assert reflected_weight_composition(1, 1, 1) == 1
        assert reflected_weight_composition(0, 0, 2) == 0

    def test_half_space_reduction_theorem_is_conditional_reduction(self):
        """The reflected identity verifies supplied BV data; it does not construct them."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "connections" / "affine_half_space_bv.tex").read_text()
        start = source.index(r"\begin{theorem}[Reduction from the reflected weight identity;")
        end = source.index(r"\end{theorem}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing \(\gamma+\delta+\varepsilon\)",
            r"Assumption~\ref{assum:kz-analytic-sdr-package}",
            "renormalized logarithmic half-space BV theory",
            "reflected weight identity",
            "This theorem is a reduction theorem, not a construction",
            "Costello effective action, counterterms, or all-loop QME",
            "After the renormalized BV data have been constructed",
        )
        for needle in required:
            assert " ".join(needle.split()) in flat

        retired = (
            r"\begin{theorem}[Reduction from the reflected weight identity]%",
            "the reflected weight identity constructs the Costello effective action",
            "the reflected weight identity proves the all-loop QME",
        )
        for needle in retired:
            assert needle not in block

    def test_reflected_obstruction_profile(self):
        profile = reflected_obstruction_profile()

        assert "Res_{D12}" in profile["quadratic"]
        assert "P(x1,x2)-sigma P(x1,bar x2)" in profile["quadratic"]
        assert "sum_cyc" in profile["cubic"]
        assert "[pi_2,pi_2]" in profile["cubic"]
        assert profile["vanishing"] == "O_3=0 by reflected Arnold relation"

    def test_reflected_arnold_relation_profile(self):
        profile = reflected_arnold_relation_profile()

        assert "omega_12^H wedge omega_23^H" in profile["relation"]
        assert profile["cancels_cubic_obstruction"] is True
        assert "doubled FM compactification" in profile["mechanism"]

    def test_affine_half_space_source_contains_image_charge_theorem(self):
        import os

        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        source = os.path.join(repo_root, 'chapters', 'connections', 'affine_half_space_bv.tex')
        with open(source, encoding='utf-8') as handle:
            text = handle.read()

        assert 'thm:image-charge-reflected-arnold-cancellation' in text
        assert 'eq:half-space-image-charge-propagator' in text
        assert 'eq:half-space-reflected-field-weight' in text
        assert 'eq:half-space-boundary-ope-coefficient' in text
        assert 'eq:half-space-reflected-weight-composition' in text
        assert 'eq:half-space-quadratic-reflected-obstruction' in text
        assert 'eq:half-space-cubic-reflected-obstruction' in text
        assert 'eq:half-space-reflected-arnold-relation' in text


# =========================================================================
# 5. DS COMPATIBILITY
# =========================================================================

class TestDSCompatibility:
    """Verify that the level shift is compatible with DS reduction."""

    def test_sl2_to_virasoro(self):
        """DS reduction sl_2 -> Virasoro: level shift is compatible."""
        result = verify_level_shift_ds_compatible('A', 1, 'principal')
        assert result['compatible'] is True
        assert result['h_vee'] == 2
        assert result['ds_target'] == 'Virasoro'
        assert result['dim_g'] == 3

    def test_sl3_to_W3(self):
        """DS reduction sl_3 -> W_3: level shift is compatible."""
        result = verify_level_shift_ds_compatible('A', 2, 'principal')
        assert result['compatible'] is True
        assert result['h_vee'] == 3
        assert result['ds_target'] == 'W_3'
        assert result['dim_g'] == 8

    def test_sl2_central_charge_at_shifted_level(self):
        """For sl_2 at k=10, DS at shifted level k_eff=8 gives
        c_DS = 1 - 6*(8+1)^2/(8+2) = 1 - 6*81/10 = 1 - 486/10 = -238/5."""
        result = verify_level_shift_ds_compatible('A', 1, 'principal')
        # k_test = 10, k_eff_test = 10 - 2 = 8
        assert result['k_eff_test'] == Rational(8)
        # c_DS(8) = 1 - 6*(8+1)^2/(8+2) = 1 - 6*81/10 = 1 - 486/10
        expected_c = 1 - Rational(6) * Rational(81) / Rational(10)
        assert result['c_ds_shifted'] == expected_c

    def test_nonprincipal_raises(self):
        """Non-principal DS reduction is not supported."""
        with pytest.raises(ValueError, match="principal"):
            verify_level_shift_ds_compatible('A', 1, 'subregular')


# =========================================================================
# 6. KAPPA CONSISTENCY
# =========================================================================

class TestKappaConsistency:
    """Cross-check kappa = (k+h^vee)*dim(g)/(2*h^vee) for Kac-Moody algebras.

    CORRECTED: kappa(KM) != c/2.  The correct formula is
    kappa = t*d/(2*h^vee) where t = k + h^vee.
    This was an AP10 violation (wrong code + wrong test = passing test).
    """

    def test_kappa_correct_formula_sl2(self):
        """kappa(sl_2, k) = 3(k+2)/4 for generic k (NOT c/2)."""
        k = Rational(7)
        kap = kappa_kac_moody(k, 'A', 1)
        # sl_2: d=3, h^vee=2, so kappa = (k+2)*3/4
        assert kap == Rational(3, 4) * (k + 2)
        # Verify NOT equal to c/2
        c = central_charge_sugawara(k, 'A', 1)
        assert kap != c / 2

    def test_kappa_correct_formula_sl3(self):
        """kappa(sl_3, k) = 4(k+3)/3 for generic k (NOT c/2)."""
        k = Rational(5)
        kap = kappa_kac_moody(k, 'A', 2)
        # sl_3: d=8, h^vee=3, so kappa = (k+3)*8/6 = 4(k+3)/3
        assert kap == Rational(4, 3) * (k + 3)
        # Verify NOT equal to c/2
        c = central_charge_sugawara(k, 'A', 2)
        assert kap != c / 2

    def test_kappa_antisymmetry_sl2(self):
        """kappa(sl_2, k) + kappa(sl_2, -k-4) = 0."""
        k = Rational(7)
        kap = kappa_kac_moody(k, 'A', 1)
        kap_dual = kappa_kac_moody(-k - 4, 'A', 1)
        assert kap + kap_dual == 0

    def test_kappa_critical_vanishes(self):
        """At the critical level k = -h^vee, kappa_ch vanishes."""
        # sl_2: critical at k = -2
        assert kappa_kac_moody(-2, 'A', 1) == 0

    def test_kappa_E8_positive(self):
        """kappa(E_8, k=31) should be positive (unitary regime)."""
        kap = kappa_kac_moody(31, 'E', 8)
        assert kap is not None
        assert kap > 0

    def test_kappa_E8_formula(self):
        """kappa(E_8, k) = 62(k+30)/15."""
        k = Rational(31)
        kap = kappa_kac_moody(k, 'E', 8)
        assert kap == Rational(62, 15) * (k + 30)
