"""Structural checks for the six-dimensional hCS one-loop anomaly."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


DELIGNE_RIGID_QUARTIC = frozenset({"A1", "G2", "D4", "F4", "E7", "E8"})


@dataclass(frozen=True)
class SixDHCSAnomalyProfile:
    """Local BV anomaly profile for holomorphic Chern-Simons on a CY3."""

    field_complex: str
    classical_action: str
    local_anomaly_functional: str
    local_cohomology_class: str
    quartic_trace_identity: str


def six_d_hcs_anomaly_profile() -> SixDHCSAnomalyProfile:
    """Return the formula package for the CY3 hCS one-loop anomaly."""

    return SixDHCSAnomalyProfile(
        field_complex="A in Omega^{0,*}(Y,g)[1]",
        classical_action=(
            "int_Y Omega_Y wedge (1/2 <A,dbar A> + 1/6 <A,[A,A]>)"
        ),
        local_anomaly_functional=(
            "A^(1) = int_Y Omega_Y wedge Tr_g(A partial A partial A partial A)"
        ),
        local_cohomology_class="H^1_loc(Obs_cl(Y), Q+{S,-})",
        quartic_trace_identity="Tr_g(X^4)=lambda_g (Tr_g(X^2))^2",
    )


def k3_euler_anomaly_factor() -> dict[str, Fraction | int | str]:
    """Return the K3 Euler factor multiplying the one-loop coefficient."""

    heat_kernel_denominator = Fraction(1, 24)
    euler_number = 24
    p1_integral = -48
    signature = -16
    return {
        "integral_c2_TK3": euler_number,
        "integral_e_TK3": euler_number,
        "heat_kernel_denominator": heat_kernel_denominator,
        "product": heat_kernel_denominator * euler_number,
        "formula": "(1/24) int_{K3} c_2(TK3) = 1",
        "heat_kernel_source": "analytic heat-kernel / Bernoulli regularisation",
        "euler_source": "Euler class / Chern-Gauss-Bonnet",
        "hrr_role": (
            "checks chi(O_K3)=2 and chi(T_K3)=-20; "
            "does not produce the heat-kernel denominator"
        ),
        "not_euler_source": "Hirzebruch signature formula",
        "signature_formula_check": (
            "sigma(K3)=(1/3) int_{K3} p_1(TK3)=-16"
        ),
        "integral_p1_TK3": p1_integral,
        "signature": signature,
    }


def quartic_trace_identity_status(lie_label: str) -> dict[str, bool | str]:
    """Classify whether the simple label is in the rigid Deligne locus."""

    if lie_label in DELIGNE_RIGID_QUARTIC:
        return {
            "quartic_factors_through_quadratic_square": True,
            "requires_extra_cubic_inflow": False,
            "status": "rigid Deligne quartic identity",
        }
    if lie_label == "A2-refined":
        return {
            "quartic_factors_through_quadratic_square": True,
            "requires_extra_cubic_inflow": True,
            "status": "requires critical twist plus cubic inflow",
        }
    return {
        "quartic_factors_through_quadratic_square": False,
        "requires_extra_cubic_inflow": lie_label in {"A2", "E6"},
        "status": "not in the rigid cancellation locus",
    }


def shadow_vs_gauge_cancellation_profile() -> dict[str, bool | str]:
    """Separate algebraic shadow detection from gauge-theoretic cancellation."""

    return {
        "shadow_detects_residual_quartic_class": True,
        "shadow_constructs_green_schwarz_cancellation": False,
        "green_schwarz_or_axion_counterterm_is_extra_data": True,
        "distinction": (
            "The shadow tower records the residual quartic class; it does not "
            "construct the BV counterterm that cancels the hCS anomaly."
        ),
    }
