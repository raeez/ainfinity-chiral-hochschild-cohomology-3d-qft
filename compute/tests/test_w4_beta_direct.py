"""Direct W4 finite-rank attack on the beta_4 discriminator."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path

from compute.lib.w4_beta_direct import (
    direct_w4_attack_report,
    harmonic_beta4_candidate,
    quadratic_beta4_candidate,
    triangular_beta4_candidate,
    w4_attack_discriminates_polynomial_interpolants,
    w4_conditional_A4_from_scaling,
    w4_kappa_over_c_from_lanes,
    w4_kappa_ratio_to_virasoro,
    w4_lane_sum_beta_candidate,
    w4_required_A5_from_bridge,
    w4_required_riccati_ratio,
    w4_riccati_bridge_requirement,
    w4_spin_lanes,
)


VOL2_ROOT = Path(__file__).resolve().parents[2]


def test_w4_direct_spin_lane_norms():
    lanes = w4_spin_lanes()
    assert lanes[2].top_ope_norm_over_c == Fraction(1, 2)
    assert lanes[3].top_ope_norm_over_c == Fraction(1, 3)
    assert lanes[4].top_ope_norm_over_c == Fraction(1, 4)
    assert w4_kappa_over_c_from_lanes() == Fraction(13, 12)
    assert w4_kappa_ratio_to_virasoro() == Fraction(13, 6)


def test_w4_direct_lane_sum_selects_thirteen():
    assert w4_lane_sum_beta_candidate() == Fraction(13)
    assert harmonic_beta4_candidate() == Fraction(13)
    assert triangular_beta4_candidate() == Fraction(15)
    assert quadratic_beta4_candidate() == Fraction(16)
    assert w4_attack_discriminates_polynomial_interpolants()


def test_w4_direct_attack_keeps_beta_law_conjectural():
    report = direct_w4_attack_report()
    assert report.lane_sum == Fraction(13)
    assert report.required_bridge_ratio == Fraction(-52, 5)
    assert report.conditional_A4_from_scaling == Fraction(13, 3)
    assert report.required_A5_from_bridge == Fraction(-676, 15)
    assert report.observed_full_miura_A5 is None
    assert report.proves_beta4 is False
    assert "A_5(W4)" in report.missing_bridge
    assert "-52/5" in report.missing_bridge


def test_w4_required_riccati_bridge_is_explicit_and_absent():
    bridge = w4_riccati_bridge_requirement()
    assert w4_required_riccati_ratio() == Fraction(-52, 5)
    assert w4_conditional_A4_from_scaling() == Fraction(13, 3)
    assert w4_required_A5_from_bridge() == Fraction(-676, 15)
    assert bridge.required_ratio == Fraction(-52, 5)
    assert bridge.required_A5 == Fraction(-676, 15)
    assert bridge.observed_full_miura_A5 is None
    assert bridge.bridge_verified is False
    assert "No full W4 Miura/OPE computation" in bridge.obstruction


def test_w4_missing_bridge_status_is_pinned_in_manuscript():
    beta_tex = (
        VOL2_ROOT / "chapters/theory/beta_N_closed_form_all_platonic.tex"
    ).read_text()
    tempered_tex = (
        VOL2_ROOT / "chapters/theory/wn_tempered_closure_platonic.tex"
    ).read_text()

    beta_label = beta_tex.index(r"\label{conj:beta-N-harmonic-closed-form}")
    beta_status_window = beta_tex[beta_label: beta_label + 220]
    beta_flat = " ".join(beta_tex.split())
    assert r"\ClaimStatusConjectured" in beta_status_window
    assert "test the harmonic scaling hypothesis rather than prove it" in beta_flat
    assert "This is a spin-lane witness, not a Riccati bridge" in beta_flat
    assert "full Miura/OPE computation" in beta_tex
    assert "is absent" in beta_tex

    tempered_label = tempered_tex.index(r"\label{thm:wn-tempered-all-N}")
    tempered_status_window = tempered_tex[tempered_label: tempered_label + 220]
    tempered_flat = " ".join(tempered_tex.split())
    assert r"\ClaimStatusConditional" in tempered_status_window
    assert "finite envelope as an explicit assumption" in tempered_flat
    assert "does not derive it from Fateev--Lukyanov" in tempered_flat
