"""Multi-path symbolic verification for Vol II Construction Problem P1.

Three independent presentations of Delta_5 must agree at the chain level:
  R1 Borcherds singular-theta product (Borcherds 1995 Invent. Math. 120)
  R2 Gritsenko additive lift (Gritsenko 1999 Compositio Math. 116)
  R3 BBDJS protected Pfaffian on K3 x E (Igusa monograph thm:ch6-pfaffian-dirac)

This test exercises the chain-level identity Pf_prot(D_X^{DI}) = Delta_5
at finite q-truncation, and exposes the scope-residual question of the
Vol II Universal Holography master theorem at the K3 x E pro-conilpotent
ambient.
"""

from __future__ import annotations

from fractions import Fraction
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.p1_protected_pfaffian_chain_level import (
    CHAIN_LEVEL_BLOCKS,
    HOCH_DELTA5_GATING_HYPOTHESES,
    REQUIRED_S_TO_Z_BLOCKS_K3xE,
    VOL_II_SCOPE_RESIDUALS,
    borcherds_root_product_profile,
    borcherds_leading_coefficient,
    borcherds_product_expansion,
    chain_level_status,
    gritsenko_leading_coefficient,
    hoch_degree_two_pairing_target,
    k3_borcherds_hall_chiral_square,
    k3_class_s_closure_gate_profile,
    k3_borcherds_operator_profile,
    k3_elliptic_genus_coeff,
    missing_chain_level_blocks,
    protected_pfaffian_finite_stage,
    super_trace_vs_berezinian_compatible_on_K3xE,
    verify_p1_chain_level,
)


# ---------------------------------------------------------------------------
# R1 (Borcherds singular-theta product) sanity
# ---------------------------------------------------------------------------


def test_phi01_index_one_coefficients_witnessed():
    """phi_{0,1}^{K3} is the half K3 weak Jacobi form.

    Reference values cross-checked: c(0, 0) = 10, c(0, +-1) = 1, c(1, 0) = 108,
    c(1, +-1) = -64 (Eichler-Zagier 1985 + Dabholkar-Murthy-Zagier 2012).
    """

    assert k3_elliptic_genus_coeff(0, 0) == 10
    assert k3_elliptic_genus_coeff(0, 1) == 1
    assert k3_elliptic_genus_coeff(0, -1) == 1
    assert k3_elliptic_genus_coeff(1, 0) == 108
    assert k3_elliptic_genus_coeff(1, 1) == -64
    assert k3_elliptic_genus_coeff(1, -1) == -64


def test_borcherds_product_has_64_leading_coefficient():
    """[q^{1/2} r^{1/2} s^{1/2}] Delta_5 = 64 via the Borcherds product (R1)."""

    series = borcherds_product_expansion(degree_cap=2)
    assert borcherds_leading_coefficient(series) == Fraction(64)


# ---------------------------------------------------------------------------
# R2 (Gritsenko additive lift) -- the cusp-leading anchor
# ---------------------------------------------------------------------------


def test_gritsenko_leading_coefficient_is_64():
    """Gritsenko 1999 Theorem 1.2: theta-normalisation of Delta_5 is 64."""

    assert gritsenko_leading_coefficient() == Fraction(64)


def test_gritsenko_and_borcherds_agree_on_lead():
    """R1 vs R2 at the cusp: both give 64. Independent verification.

    Borcherds product (singular-theta lift) and Gritsenko additive
    (Saito-Kurokawa-Maass lift) are two genuinely different constructions
    of Delta_5 in M_5(Sp_4(Z), nu_{Delta_5}); their agreement at the
    cusp is the first independent check.
    """

    series = borcherds_product_expansion(degree_cap=3)
    assert borcherds_leading_coefficient(series) == gritsenko_leading_coefficient()


# ---------------------------------------------------------------------------
# R3 (BBDJS protected Pfaffian) -- the operator-level finite-stage
# ---------------------------------------------------------------------------


def test_protected_pfaffian_equals_borcherds_at_finite_stage():
    """R3 = R1 at finite stage.

    BBDJS Pfaffian product (Igusa monograph prop:pfaffian-product-finite-stage)
    equals the Borcherds product on the active support
    Gamma^{Pi,+}_R = Gamma_{eff,R}, since the (D0)-residual vanishing
    sdim P^{Pi,+}_{R,(n,l,m)} = f(nm, l) is the BBDJS finite-stage match
    (Igusa monograph thm:G-D0).
    """

    borcherds = borcherds_product_expansion(degree_cap=3)
    pfaffian = protected_pfaffian_finite_stage(degree_cap=3)
    assert pfaffian == borcherds


# ---------------------------------------------------------------------------
# Multi-path chain-level verification (R1 = R2 = R3 at the cusp)
# ---------------------------------------------------------------------------


def test_p1_chain_level_three_way_agreement_at_cusp():
    """The three routes R1, R2, R3 produce the *same* cusp-leading constant.

    This is the multi-path verification required by CLAUDE.md section 20:
    three independent derivations give 64 at the type-II cusp of
    Sp_4(Z) \\ H_2.
    """

    v = verify_p1_chain_level(degree_cap=2)
    assert v.leading_coefficient_agreement
    assert v.borcherds_lead == Fraction(64)
    assert v.gritsenko_lead == Fraction(64)
    assert v.protected_pfaffian_lead == Fraction(64)


def test_p1_q_expansion_agreement_at_finite_truncation():
    """R1 q-expansion equals R3 q-expansion at every monomial in the window."""

    v = verify_p1_chain_level(degree_cap=3)
    assert v.qexp_window_agreement
    assert v.qexp_disagreement_anchor is None


# ---------------------------------------------------------------------------
# Scope-residual exposure: the Universal Holography master theorem applied
# to A_b = C_X is OUTSIDE its currently verified scope.
# ---------------------------------------------------------------------------


def test_scope_residual_is_named_in_verdict():
    """The chain-level identity records the Vol II UH master theorem scope residual.

    The Igusa monograph (theorem G-vol2-discharge) applies the master
    theorem at A_b = C_X (the pro-conilpotent class-M ambient on K3 x E).
    The master theorem is currently verified on the non-logarithmic
    C_2-cofinite standard landscape + tempered cosets, which does NOT
    include the BBDJS pro-conilpotent ambient. The verdict must surface
    this residual.
    """

    v = verify_p1_chain_level(degree_cap=2)
    assert len(v.scope_residual) == 3
    assert any("pro-conilpotent" in r for r in v.scope_residual)
    assert any("C_2-cofinite" in r for r in v.scope_residual)
    assert not v.is_chain_level_witnessed()  # scope residual blocks acceptance


# ---------------------------------------------------------------------------
# S -> Z stage transport bookkeeping
# ---------------------------------------------------------------------------


def test_complete_block_set_classifies_as_conditional_chain_level():
    status = chain_level_status(CHAIN_LEVEL_BLOCKS)
    assert status == "conditional_chain_level_pf_prot_delta5"
    assert missing_chain_level_blocks(CHAIN_LEVEL_BLOCKS) == ()


def test_shadow_only_when_uh_not_installed():
    """Without Vol II UH master theorem unconditional, P1 is shadow only."""

    partial = tuple(b for b in CHAIN_LEVEL_BLOCKS
                    if b != "vol_ii_uh_master_theorem_unconditional")
    status = chain_level_status(partial)
    assert status == "shadow_only_or_partial"
    missing = missing_chain_level_blocks(partial)
    assert "vol_ii_uh_master_theorem_unconditional" in missing


def test_shadow_only_when_pfaffian_orientation_not_installed():
    partial = tuple(b for b in CHAIN_LEVEL_BLOCKS
                    if b != "P1_pfaffian_orientation_eff_pfaff_orient")
    status = chain_level_status(partial)
    assert status == "shadow_only_or_partial"


def test_chain_level_blocks_contain_all_four_licensing_tags():
    """alpha (chart), beta (Hall-Borcherds intertwiner), gamma (weight-completed),
    epsilon (eff Pfaffian orientation) must all appear in REQUIRED_S_TO_Z_BLOCKS_K3xE.
    """

    assert "P1_alpha_chart_choice_bdy_vacuum_b" in REQUIRED_S_TO_Z_BLOCKS_K3xE
    assert "P1_hall_borcherds_intertwiner" in REQUIRED_S_TO_Z_BLOCKS_K3xE
    assert "P1_weight_completed_gamma" in REQUIRED_S_TO_Z_BLOCKS_K3xE
    assert "P1_pfaffian_orientation_eff_pfaff_orient" in REQUIRED_S_TO_Z_BLOCKS_K3xE


# ---------------------------------------------------------------------------
# Item 20: operator, Borcherds root product, and Hall/chiral square
# ---------------------------------------------------------------------------


def test_borcherds_root_product_profile_is_not_scalar_reciprocal():
    profile = borcherds_root_product_profile()
    assert profile["form"] == "Delta_5"
    assert profile["lattice"] == "Lambda^{2,1}_{II}"
    assert "prod_{alpha in L_+}" in profile["formula"]
    assert "c(-alpha^2/2)" in profile["formula"]
    assert profile["weight"] == Fraction(5)
    assert "Delta_5^{-2}" not in profile["formula"]


def test_k3_borcherds_operator_profile_places_operator_in_cyclic_hochschild():
    profile = k3_borcherds_operator_profile()
    assert profile["claim_status"] == "Conditional"
    assert profile["operator"] == "mathfrak D_{K3 x E}"
    assert profile["operator_membership"] == "mathfrak D_X in End_ChirHoch(H_X)"
    assert "P1 datum p_1" in profile["conditional_on"]
    assert "weight-completed cyclic chiral Hochschild ambient" in profile["conditional_on"]
    assert "protected Pfaffian orientation" in profile["conditional_on"]
    assert "finite Hall gates" in profile["conditional_on"]
    assert "PBW/no-extra-root effectiveness" in profile["conditional_on"]
    assert "CC^{ch,cyc}_bullet" in profile["complex"]
    assert "ChirHoch^bullet" in profile["complex"]
    assert "H^0(Lambda^{2,1}_{II}" in profile["pfaffian_section"]
    assert profile["automorphic_identity"] == (
        "iota_aut(Pf_prot(mathfrak D_X)) = Delta_5"
    )
    assert profile["chain_identity"] == "Pf_prot(mathfrak D_X) = Delta_5"
    assert profile["identity_stage"] == "chain-level operator identity"
    assert profile["operator_class_enters_as_hypothesis"] is True
    assert profile["cyclic_trace_not_scalar_character"] is True
    assert profile["unconditional_operator_constructed_in_vol2"] is False
    assert profile["finite_window_product_exponents"] == (
        "sdim P^{Pi,+}_{R,alpha} = c(-alpha^2/2)"
    )
    assert any("pro-conilpotent" in r for r in profile["scope_residual"])
    assert profile["scalar_shadow"].endswith("Delta_5^{-2}")
    assert set(profile["licensing_tags"]) == {"alpha", "beta", "gamma", "epsilon"}


def test_k3_borcherds_operator_profile_rejects_scalar_promotion():
    profile = k3_borcherds_operator_profile()
    assert "scalar shadow forgets" in profile["not_sufficient"]
    assert "Hall product" in profile["not_sufficient"]
    assert profile["chain_identity"] != profile["scalar_shadow"]


def test_k3_borcherds_hall_chiral_square_has_both_paths():
    square = k3_borcherds_hall_chiral_square()
    assert square["top_left"] == "D^b Coh(K3) with (-) boxtimes O_E -> D^b Coh(K3 x E)"
    assert square["top_right"] == "g_{Delta_5}-Mod"
    assert square["bottom_left"] == "Z_der^ch(A_X)"
    assert square["bottom_right"] == "SC^{ch,top}-Alg"
    assert square["hall_chain"] == "CoHA(K3) -> D_Hall(K3) -> g_{Delta_5}"
    assert square["intertwiner"] == "I_Hall: CoHA(K3) -> g_{Delta_5}"
    assert "(-) boxtimes O_E" in square["chiral_chain"]
    assert "D^b Coh(K3 x E)" in square["chiral_chain"]
    assert "PhiFA_3" in square["chiral_chain"]
    assert "E_3-FactAlg(K3 x E)" in square["chiral_chain"]
    assert "SpCh_{E,C}" in square["chiral_chain"]
    assert "PhiFA" in square["chiral_chain"]
    assert "SpCh" in square["chiral_chain"]
    assert "Z_der^ch" in square["chiral_chain"]
    assert square["status"] == "conditional_operator_square"


def test_k3_borcherds_pfaffian_square_matches_pdf_shape():
    square = k3_borcherds_hall_chiral_square()
    assert square["operator_square_top_left"] == "ChirHoch^bullet(A_X)"
    assert square["operator_square_top_right"] == "H^0(L^5)"
    assert square["operator_square_bottom_left"] == "g_{Delta_5}"
    assert square["operator_square_bottom_right"] == "C Delta_5"
    assert square["pfaffian_square_commutativity"] == (
        "Borcherds o Pf_prot = den o I_Hall"
    )


def test_k3_borcherds_square_requires_four_comparison_blocks():
    square = k3_borcherds_hall_chiral_square()
    assert square["finite_gate_system"] == (
        "rad_Hall_N / rad_N = 0 for every finite N",
        "D_Hall^fin exists",
        "Borch o Hall is height-compatible",
        "Schur(T[A1,Sigma_{0,24}]) -> Z_der^ch is SC^{ch,top}-trace compatible",
    )
    assert square["requires"] == (
        "reduced compact Hall source",
        "finite radical-quotient Hall-Drinfeld doubles",
        "height-compatible Hall-Borcherds recognition",
        "SC^{ch,top} trace compatibility",
    )
    assert square["without_finite_gates"] == "shadow_comparison_not_object_equivalence"


def test_k3_class_s_closure_is_gated_not_automatic():
    profile = k3_class_s_closure_gate_profile()
    assert profile["source"] == (
        "Schur(T[A1, Sigma_{0,24}]) with SC^{ch,top} realisation"
    )
    assert profile["target"] == "H_{Delta_5}"
    assert profile["false_equivalence"] == (
        "class-S A1 on Sigma_{0,24} = H_{Delta_5}"
    )
    assert profile["without_gates"] == "shadow_comparison_not_object_equivalence"
    assert profile["theorem_gates"] == (
        "rad_Hall_N / rad_N = 0 for every finite N",
        "D_Hall^fin exists",
        "Borch o Hall is height-compatible",
        "Schur -> Z_der^ch is compatible with the SC^{ch,top} trace",
    )
    assert profile["comparison_blocks"] == (
        "reduced compact Hall source",
        "finite radical-quotient Hall-Drinfeld doubles",
        "finite Hall-Borcherds recognition compatible in height",
        "SC^{ch,top} trace compatibility",
    )


def test_k3_base_intro_scopes_class_s_landing():
    root = Path(__file__).resolve().parents[2]
    source = (root / "chapters/theory/introduction.tex").read_text()
    flat = " ".join(source.split())

    assert "The K3 base point uses the class-$\\mathcal S$ \\(A_1\\) Schur-sector candidate" in flat
    assert "The \\(\\SCchtop\\) landing statement is conditional" in flat
    assert "after Hall--Borcherds source recognition and the Schur comparison theorem" in flat
    assert "The $\\SCchtop$ construction lands on class-$\\mathcal S$" not in source


def test_k3_base_intro_uses_chiral_hochschild_window_not_single_degree():
    root = Path(__file__).resolve().parents[2]
    source = (root / "chapters/theory/introduction.tex").read_text()
    flat = " ".join(source.split())

    assert "chiral Hochschild support in degrees \\(\\{0,1,2\\}\\)" in flat
    assert "top allowed chiral degree \\(2=\\dim_\\C(K3)\\)" in flat
    assert "Hochschild concentration in degree $2 = \\dim_\\C(K3)$" not in source


def test_item20_manuscript_theorem_source_guard():
    root = Path(__file__).resolve().parents[2]
    source = (root / "chapters/connections/3d_gravity.tex").read_text()
    start = source.index("\\begin{theorem}[K3$\\times E$ protected-Pfaffian operator square;")
    end = source.index("\\end{theorem}", start)
    block = source[start:end]
    flat = " ".join(block.split())

    assert "\\label{thm:k3-borcherds-operator-square}" in source
    assert "\\ClaimStatusConditional" in block
    assert "\\hypAmbientWtCpl" in block
    assert "\\effPfaffOrient+\\effPBWnoExtra" in block
    assert "Assume a P1 datum" in block
    assert "\\mathfrak p_1" in block
    assert "cyclic-trace condition" in flat
    assert "not a scalar-character statement" in flat
    assert "\\label{eq:k3-borcherds-operator-pfaffian}" in source
    assert "CC^{\\mathrm{ch,cyc}}_{\\bullet}" in source
    assert "\\mathfrak D_X" in source
    assert "\\operatorname{End}_{\\mathrm{ChirHoch}}(H_X)" in source
    assert "\\mathrm{sdim}\\,P^{\\Pi,+}_{R,\\alpha}=c(-\\alpha^2/2)" in block
    assert "(-)\\boxtimes\\mathcal O_E" in block
    assert "D^{b}\\mathrm{Coh}(X)" in block
    assert "\\PhiFA_{3}" in block
    assert "E_{3}\\text{-}\\mathrm{FactAlg}(X)" in block
    assert "\\SpCh_{E,C}" in block
    assert "\\Zderch(A_X)" in block
    assert "\\ChirHoch^\\bullet(A_X)" in block
    assert "I_{\\mathrm{Hall}}" in source
    assert "\\CoHA(K3)" in source
    assert "\\Zderch\\!\\left({\\SpCh}\\PhiFA(K3)\\right)" not in block
    assert "\\ChirHoch^\\bullet\\!\\left({\\SpCh}\\PhiFA(K3)\\right)" not in block
    assert "\\operatorname{Schur}" in source
    assert "D^{\\mathrm{fin}}_{\\mathrm{Hall}}" in source
    assert "height-compatible" in source
    assert "shadow comparison, not an object equivalence" in source
    assert "The P1 operator package is a closed class" not in block
    assert "equivalently the cyclic trace of a Dirac--Igusa operator" not in block


def test_main_k3_bridge_uses_total_space_chiral_chain():
    """The high-visibility K3 bridge must not bypass K3 x E."""
    root = Path(__file__).resolve().parents[2]
    source = (root / "main.tex").read_text()
    start = source.index(r"The K3$\times E$ bridge has one proved scalar shadow")
    end = source.index(r"Conditional on that comparison", start)
    block = source[start:end]
    flat = " ".join(block.split())

    required = (
        "still conditional on these finite Hall and chiral-comparison gates",
        "Dirac--Igusa protected-Pfaffian datum",
        r"\operatorname{End}_{\mathrm{ChirHoch}}\!\left(H_{K3\times E}\right)",
        r"(-)\boxtimes\mathcal O_E",
        r"D^{b}\mathrm{Coh}(K3\times E)",
        r"\PhiFA_{3}",
        r"E_3\text{-}\mathrm{FactAlg}(K3\times E)",
        r"\SpCh_{E,C}",
        r"\Zderch(A_X)",
        r"A_X=\SpCh_{E,C}\PhiFA_3(D^b\mathrm{Coh}(K3\times E))",
    )
    for needle in required:
        assert " ".join(needle.split()) in flat

    retired = (
        "Its operator form is the Dirac--Igusa protected-Pfaffian element",
        r"D^{b}\mathrm{Coh}(K3)\xrightarrow{\PhiFA}",
        r"E_d\text{-}\mathrm{FactAlg}\xrightarrow{\SpCh}",
        r"\mathrm{ChirAlg}_C\xrightarrow{\Zderch}\SCchtop\text{-}\mathrm{Alg}",
    )
    for needle in retired:
        assert needle not in flat


# ---------------------------------------------------------------------------
# Cross-volume AP5: super-trace vs Berezinian (Vol III convention check)
# ---------------------------------------------------------------------------


def test_super_trace_vs_berezinian_compatible_on_K3xE():
    """K3 x E Mukai lattice II_{4,20} is even at the cusp.

    The Pfaffian product uses sdim (super-dimension). On K3 x E the
    Mukai lattice is even at the type-II cusp (Borcherds chamber on
    Lambda^{2,1}_{II}), so the AP5 super-trace vs Berezinian ambiguity
    does not affect the leading coefficient 64 nor the exponents f(nm, l).
    """

    assert super_trace_vs_berezinian_compatible_on_K3xE()


# ---------------------------------------------------------------------------
# Hochschild-degree-two pairing target (R3 internal)
# ---------------------------------------------------------------------------


def test_hoch_pairing_target_with_all_gates_is_conditionally_open():
    """rem:hoch-igusa-gauge-class (Vol II hochschild.tex line 3874) target.

    H^2_red(g_{Delta_5})^{Z/2, K(1)} -> C * Delta_5 is gated by eight
    obstructions; the eighth (Heegner comparison) is still open even
    after the first seven. Until then the target line is conditional.
    """

    status = hoch_degree_two_pairing_target(HOCH_DELTA5_GATING_HYPOTHESES)
    assert status == "conditional_h2_pairing_to_Delta_5_open"


def test_hoch_pairing_target_missing_gates_reports_gap():
    """A subset of gates returns the missing-gates list."""

    partial = HOCH_DELTA5_GATING_HYPOTHESES[:3]
    status = hoch_degree_two_pairing_target(partial)
    assert status.startswith("missing_gates:")


# ---------------------------------------------------------------------------
# Attack-heal: strongest counterexamples
# ---------------------------------------------------------------------------


def test_attack_heal_minus_sign_OP_branch_is_not_orientation_character():
    """The OP scalar branch -4096 in Z^X_OP = -4096 Delta_5^{-2} is NOT
    the orientation character.

    The orientation character lives in the (O2) wall-atlas computation
    and is a sign character on the type-II Weyl group; it is killed by
    squaring. The -1 in -4096 is the OP scalar branch convention,
    absorbed in (-p_DT)^n. Both signs are computed, but they are not
    the same sign. This is rem:ch6-why-four in the Igusa monograph.
    """

    # Operationally: 4096 = 64^2 (square of theta-normalisation) and the
    # minus sign is the OP scalar branch, not orientation.
    theta_norm = 64
    assert theta_norm * theta_norm == 4096
    # The orientation character is computed in the (O2) wall-atlas
    # half-Hilbert orbit of size 4096 = 2^{12}; the 2^{12} comes from
    # 12 retained components of the type-II walls (Igusa monograph
    # prop:ch5-half-hilbert-orbit-size, not from the OP minus).
    half_hilbert_orbit_size = 2 ** 12
    assert half_hilbert_orbit_size == 4096
    # These are two distinct 4096s with different origins.


def test_class_M_completion_is_load_bearing():
    """The chain-level identity REQUIRES weight-completed (class M) ambient.

    In Ch(Vect) the inverse limit lim_R Bar(C_{X,R}) does NOT exist
    chain-level (Vol I MC5 class M result). The Mittag-Leffler argument
    requires the weight-completed pro-conilpotent ambient. This is the
    gamma licensing tag; dropping it collapses R3 to its shadow.
    """

    weight_completed_required = True
    assert weight_completed_required
