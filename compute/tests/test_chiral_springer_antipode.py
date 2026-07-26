from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.chiral_springer_antipode import (
    HOPF_ANTIPODE_IDENTITIES,
    MC_ORIENTATION_REVERSAL_EQUATION,
    antipode_identity_set_complete,
    chiral_springer_antipode_profile,
    mc_orientation_reversal_status,
)


def test_chiral_springer_carrier_uses_complementary_lagrangian():
    profile = chiral_springer_antipode_profile()

    assert profile.steinberg_object == "St^ch(A) = L_A x^h_{M_vac} L_A^vee"
    assert profile.lagrangian_left == "L_A"
    assert profile.lagrangian_right == "L_A^vee (Verdier/Koszul-complementary)"


def test_convolution_and_antipode_formulas_match_contract():
    profile = chiral_springer_antipode_profile()

    assert profile.convolution_space == "L_A x^h_M L_A^vee x^h_M L_A"
    assert profile.convolution_formula == (
        "F star G = p_13,*(p_12^*F tensor p_23^*G)"
    )
    assert profile.point_involution == "iota(x, y) = (y, x)"
    assert profile.sheaf_antipode == "S(F) = D_Verdier(iota^*F)"


def test_four_hopf_antipode_identities_are_exact():
    profile = chiral_springer_antipode_profile()

    assert profile.hopf_identities == HOPF_ANTIPODE_IDENTITIES
    assert antipode_identity_set_complete(profile.hopf_identities)
    assert profile.hopf_identities == (
        "m(S tensor 1) Delta = eta epsilon",
        "m(1 tensor S) Delta = eta epsilon",
        "Delta S = (S tensor S) Delta^op",
        "epsilon S = epsilon",
    )


def test_mc_homotopy_gate_distinguishes_full_antipode_from_shadow():
    blocked = mc_orientation_reversal_status(class_label="M", has_homotopy=False)
    supplied = mc_orientation_reversal_status(class_label="M", has_homotopy=True)

    assert blocked["condition"] == MC_ORIENTATION_REVERSAL_EQUATION
    assert blocked["status"] == "associated_graded_anti_equivalence_only"
    assert blocked["full_antipode"] is False
    assert supplied["status"] == "conditional_full_antipode"
    assert supplied["full_antipode"] is True


def test_class_g_has_vacuous_higher_mc_obstruction():
    status = mc_orientation_reversal_status(class_label="G", has_homotopy=False)

    assert status["homotopy_required"] is False
    assert status["status"] == "verified_gaussian_antipode"
    assert status["full_antipode"] is True


def test_chiral_springer_antipode_source_guard():
    root = Path(__file__).resolve().parents[2]
    source = (root / "chapters/connections/spectral-braiding-frontier.tex").read_text()

    assert "\\label{constr:chiral-springer-antipode-carrier}" in source
    assert "\\Stch(\\cA)" in source
    assert "\\cL_{\\cA}^{\\vee}" in source
    assert "F\\star G" in source
    assert "p_{13,*}" in source
    assert "\\mathbb{D}_{\\mathrm{Verdier}}(\\iota^*F)" in source
    assert "S_{\\cA}(\\Theta_{\\cA})" in source
    assert "-\\Theta_{\\cA}+d_{\\Theta_{\\cA}}(\\Xi)" in source
    assert "associated-graded anti-equivalence" in source
