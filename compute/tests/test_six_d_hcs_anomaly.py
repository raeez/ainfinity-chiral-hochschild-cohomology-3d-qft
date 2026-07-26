"""Checks for the six-dimensional hCS one-loop anomaly package."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path

from compute.lib.six_d_hcs_anomaly import (
    k3_euler_anomaly_factor,
    quartic_trace_identity_status,
    shadow_vs_gauge_cancellation_profile,
    six_d_hcs_anomaly_profile,
)


def test_six_d_hcs_local_anomaly_profile_records_pdf_formulas() -> None:
    """The CY3 hCS anomaly lives in local BV cohomology."""
    profile = six_d_hcs_anomaly_profile()

    assert "Omega^{0,*}(Y,g)[1]" in profile.field_complex
    assert "1/2 <A,dbar A>" in profile.classical_action
    assert "1/6 <A,[A,A]>" in profile.classical_action
    assert "Tr_g(A partial A partial A partial A)" in (
        profile.local_anomaly_functional
    )
    assert profile.local_cohomology_class == "H^1_loc(Obs_cl(Y), Q+{S,-})"
    assert profile.quartic_trace_identity == (
        "Tr_g(X^4)=lambda_g (Tr_g(X^2))^2"
    )


def test_k3_euler_factor_cancels_one_loop_denominator() -> None:
    """Integration over K3 contributes ``int c_2(TK3)=24``."""
    factor = k3_euler_anomaly_factor()

    assert factor["integral_c2_TK3"] == 24
    assert factor["integral_e_TK3"] == 24
    assert factor["heat_kernel_denominator"] == Fraction(1, 24)
    assert factor["product"] == 1
    assert "int_{K3} c_2(TK3)" in factor["formula"]
    assert factor["heat_kernel_source"] == (
        "analytic heat-kernel / Bernoulli regularisation"
    )
    assert factor["euler_source"] == "Euler class / Chern-Gauss-Bonnet"
    assert "does not produce the heat-kernel denominator" in factor["hrr_role"]
    assert factor["not_euler_source"] == "Hirzebruch signature formula"
    assert factor["integral_p1_TK3"] == -48
    assert factor["signature"] == -16
    assert "sigma(K3)" in factor["signature_formula_check"]


def test_k3_euler_source_is_not_signature_formula_in_manuscript() -> None:
    root = Path(__file__).resolve().parents[2]
    source = (root / "chapters/connections/bv_brst.tex").read_text()
    flat = " ".join(source.split())

    assert "The Euler characteristic of $K3$ is $\\chi(K3)=24$ by Euler-class integration" in flat
    assert "\\chi(K3)=\\int_{K3}e(T_{K3})" in flat
    assert "\\sigma(K3)=-16=\\frac13\\int_{K3}p_1(T_{K3})" in flat
    assert "It does not produce the heat-kernel denominator" in flat
    assert "\\left(\\frac{|B_2|}{4}\\right)\\int_{K3}e(T_{K3})" in flat
    assert "which follows from the Euler-class term in Riemann--Roch" not in flat
    assert "The Euler characteristic of $K3$ is $\\chi(K3)=24$ by the Hirzebruch signature formula" not in flat


def test_deligne_quartic_identity_locus_is_not_all_simple_lie() -> None:
    """The rigid quartic identity is a Deligne-locus condition."""
    a1 = quartic_trace_identity_status("A1")
    e6 = quartic_trace_identity_status("E6")
    a2_refined = quartic_trace_identity_status("A2-refined")

    assert a1["quartic_factors_through_quadratic_square"] is True
    assert a1["requires_extra_cubic_inflow"] is False
    assert e6["quartic_factors_through_quadratic_square"] is False
    assert e6["requires_extra_cubic_inflow"] is True
    assert a2_refined["quartic_factors_through_quadratic_square"] is True
    assert a2_refined["requires_extra_cubic_inflow"] is True


def test_shadow_tower_does_not_construct_green_schwarz_cancellation() -> None:
    """The algebraic shadow tower is a detector, not the BV cancellation."""
    profile = shadow_vs_gauge_cancellation_profile()

    assert profile["shadow_detects_residual_quartic_class"] is True
    assert profile["shadow_constructs_green_schwarz_cancellation"] is False
    assert profile["green_schwarz_or_axion_counterterm_is_extra_data"] is True
