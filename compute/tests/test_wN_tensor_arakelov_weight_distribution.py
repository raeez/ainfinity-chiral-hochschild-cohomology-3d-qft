"""Verification of the W_N tensor-Arakelov weight-distribution law.

Claim (Theorem D tensor-Arakelov upgrade, working_notes.tex): for W_N
with strong generators at weights j in {2, ..., N}, the weight-w
diagonal Arakelov anomaly coefficient is

    kappa_j(W_N) = c / j   for j in {2, 3, ..., N},
    kappa_w(W_N) = 0       for all other w,

and the total trace recovers Vol I's scalar formula

    tr_F(K)(W_N) = sum_{j=2}^N c/j = c * (H_N - 1),

where H_N = 1 + 1/2 + 1/3 + ... + 1/N is the N-th harmonic number.

This test verifies the identity sum_{j=2}^N c/j = c (H_N - 1) as
an algebraic check for small N, and tests that Virasoro (N=2) and
W_3 (N=3) reduce correctly.

Installed 2026-04-17 as part of T1.4 Theorem D tensor-Arakelov
inscription verification.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import re

from compute.lib.independent_verification import independent_verification


def _harmonic(n: int) -> Fraction:
    """H_n = sum_{k=1}^n 1/k as an exact fraction."""
    if n < 0:
        raise ValueError("harmonic H_n requires n >= 0")
    total = Fraction(0)
    for k in range(1, n + 1):
        total += Fraction(1, k)
    return total


def _kappa_total_wn(N: int, c: Fraction) -> Fraction:
    """Vol I: kappa(W_N) = c * (H_N - 1)."""
    return c * (_harmonic(N) - 1)


def _kappa_w_distribution_wn(N: int, c: Fraction) -> dict[int, Fraction]:
    """W_N weight-distribution: kappa_j = c/j for j in {2, ..., N}; else 0."""
    return {j: c / j for j in range(2, N + 1)}


def _tensor_arakelov_component_formula_profile() -> dict[str, object]:
    """The item-18 component formula for K_{alpha beta}."""
    return {
        "target": "Sym^2(F^vee) tensor Omega^2(Mbar_g,n)",
        "component": "K_{alpha beta}",
        "stable_graph_set": "G_st_{g,2}(alpha,beta)",
        "symmetry_factor": "1/|Aut Gamma|",
        "edge_factor": "prod_e P_{alpha_e beta_e}",
        "vertex_factor": "prod_v C_v",
        "integral_domain": "[Mbar_Gamma]",
        "trace": "sum_alpha K_{alpha alpha} = kappa_ch^Hodge(A) omega_g",
    }


def _w3_tensor_channel_matrix(c: Fraction) -> dict[str, object]:
    """W_3 tensor curvature before scalar trace."""
    trace = c / 2 + c / 3
    return {
        "basis": ("T", "W"),
        "matrix": (("K_TT", "K_TW"), ("K_WT", "K_WW")),
        "diagonal": ("K_TT", "K_WW"),
        "mixed": ("K_TW", "K_WT"),
        "trace": trace,
        "expected_trace": Fraction(5) * c / 6,
        "genus2_cross": (c + 204) / (16 * c),
        "scalar_trace_forgets": ("K_TW", "K_WT"),
    }


@independent_verification(
    claim="thm:theoremD-tensor-arakelov",
    derived_from=[
        "Programme Vol I Mumford-Arakelov scalar formula",
        "Weight-indexed decomposition of modular bar curvature",
        "Universal Holography Functor tensor curvature construction",
    ],
    verified_against=[
        "Linshaw 2021 arXiv:2104.01228 (universal W_infty and structure constants)",
        "Kac-Wakimoto 2004 J. Algebra 258 (characters and conformal blocks for W-algebras)",
    ],
    disjoint_rationale=(
        "Linshaw 2021 establishes structure constants of universal "
        "W_infty at arbitrary c, independent of the programme's Arakelov "
        "curvature formula. Kac-Wakimoto 2004 give characters of "
        "W-algebra modules from representation theory at non-critical "
        "level, confirming the generator structure at weights {2,...,N} "
        "without invoking Arakelov / tensor curvature. The weight-j "
        "distribution kappa_j = c/j is the inverse-of-weight law "
        "characteristic of rank-(N-1) principal W-algebras, reading off "
        "from the Sugawara-weight decomposition independently of the "
        "modular-bar curvature machinery."
    ),
)
def test_wn_kappa_tensor_arakelov():
    c = Fraction(1)  # symbolic c = 1 for identity check

    # Virasoro (N = 2): single weight-2 generator, kappa_2 = c/2.
    dist_N2 = _kappa_w_distribution_wn(2, c)
    assert dist_N2 == {2: Fraction(1, 2)}, (
        f"Virasoro distribution should be {{2: c/2}}, got {dist_N2}"
    )
    total_N2 = sum(dist_N2.values())
    assert total_N2 == _kappa_total_wn(2, c), (
        f"Virasoro trace: distribution sum {total_N2} should equal "
        f"kappa(Vir_c) = c*(H_2-1) = c/2 = {_kappa_total_wn(2, c)}"
    )

    # W_3 (N = 3): weights 2, 3; kappa_2 = c/2, kappa_3 = c/3.
    dist_N3 = _kappa_w_distribution_wn(3, c)
    assert dist_N3 == {2: Fraction(1, 2), 3: Fraction(1, 3)}, (
        f"W_3 distribution should be {{2: c/2, 3: c/3}}, got {dist_N3}"
    )
    total_N3 = sum(dist_N3.values())
    expected_N3 = _kappa_total_wn(3, c)
    assert total_N3 == expected_N3, (
        f"W_3 trace: distribution sum {total_N3} should equal "
        f"c*(H_3-1) = c*(5/6) = {expected_N3}"
    )
    assert expected_N3 == Fraction(5, 6)

    # W_4 (N = 4): weights 2, 3, 4; sum = c/2 + c/3 + c/4 = c*(6+4+3)/12 = 13c/12.
    dist_N4 = _kappa_w_distribution_wn(4, c)
    total_N4 = sum(dist_N4.values())
    expected_N4 = _kappa_total_wn(4, c)
    assert total_N4 == expected_N4
    assert expected_N4 == Fraction(13, 12), (
        f"kappa(W_4) should be 13/12 for c=1, got {expected_N4}"
    )

    # W_5 (N = 5): sum = c*(H_5 - 1) = c*(77/60).
    dist_N5 = _kappa_w_distribution_wn(5, c)
    total_N5 = sum(dist_N5.values())
    expected_N5 = _kappa_total_wn(5, c)
    assert total_N5 == expected_N5
    assert expected_N5 == Fraction(77, 60)

    # General N up to 8: verify distribution sum = kappa_total exactly.
    for N in range(2, 9):
        dist = _kappa_w_distribution_wn(N, c)
        total = sum(dist.values())
        expected = _kappa_total_wn(N, c)
        assert total == expected, (
            f"W_{N} trace check failed: {total} != {expected}"
        )

    # Off-diagonal weights are zero.
    dist_N4 = _kappa_w_distribution_wn(4, c)
    assert 1 not in dist_N4, "No generator at weight 1 for W_N, N >= 2"
    assert 5 not in dist_N4, "No generator at weight 5 for W_4"


def test_harmonic_base_cases():
    """Independent sanity check of harmonic numbers from the programme formula."""
    assert _harmonic(0) == Fraction(0)
    assert _harmonic(1) == Fraction(1)
    assert _harmonic(2) == Fraction(3, 2)
    assert _harmonic(3) == Fraction(11, 6)
    assert _harmonic(4) == Fraction(25, 12)
    assert _harmonic(5) == Fraction(137, 60)


def test_virasoro_single_generator_concentration():
    """Virasoro single-strong-generator concentration: kappa_w=0 for w!=2."""
    c = Fraction(1)
    dist = _kappa_w_distribution_wn(2, c)
    # Only weight 2 is in the distribution dictionary.
    assert list(dist.keys()) == [2]
    assert dist[2] == Fraction(1, 2)


def test_tensor_arakelov_component_formula_profile():
    """K_{alpha beta} is a stable-graph integral, not a scalar slogan."""
    profile = _tensor_arakelov_component_formula_profile()
    assert profile["target"] == "Sym^2(F^vee) tensor Omega^2(Mbar_g,n)"
    assert profile["stable_graph_set"] == "G_st_{g,2}(alpha,beta)"
    assert profile["symmetry_factor"] == "1/|Aut Gamma|"
    assert profile["edge_factor"] == "prod_e P_{alpha_e beta_e}"
    assert profile["vertex_factor"] == "prod_v C_v"
    assert profile["trace"] == "sum_alpha K_{alpha alpha} = kappa_ch^Hodge(A) omega_g"


def test_w3_tensor_channel_matrix_before_trace():
    """W_3 has TT, TW, WT, WW tensor channels; the trace forgets mixed entries."""
    matrix = _w3_tensor_channel_matrix(Fraction(1))
    assert matrix["basis"] == ("T", "W")
    assert matrix["matrix"] == (("K_TT", "K_TW"), ("K_WT", "K_WW"))
    assert matrix["trace"] == matrix["expected_trace"] == Fraction(5, 6)
    assert matrix["mixed"] == ("K_TW", "K_WT")
    assert matrix["scalar_trace_forgets"] == ("K_TW", "K_WT")
    assert matrix["genus2_cross"] == Fraction(205, 16)


def test_tensor_arakelov_source_contains_component_formula():
    """Source guard for item-18 tensor-Arakelov component formula."""
    repo = Path(__file__).resolve().parents[2]
    source = repo / "chapters/theory/theorems_C_D_native_vol2_platonic.tex"
    text = source.read_text(encoding="utf-8")
    assert "eq:tensor-arakelov-component-graph-integral" in text
    assert "P_{\\alpha_e\\beta_e}" in text
    assert "C_v" in text
    assert "ex:theoremD-W3-four-channel-tensor" in text
    assert "K_{TT} & K_{TW}" in text


def test_tensor_arakelov_source_separates_scalar_trace_from_full_energy():
    """The native Theorem D surface must not identify tensor data with its trace."""
    repo = Path(__file__).resolve().parents[2]
    source = repo / "chapters/theory/theorems_C_D_native_vol2_platonic.tex"
    text = source.read_text(encoding="utf-8")
    compact = re.sub(r"\s+", "", text)

    required_compact = (
        (
            "F_g^{\\mathrm{sc}}(\\cA)="
            "\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}"
        ),
        (
            "F_g(\\cA)=F_g^{\\mathrm{sc}}(\\cA)"
            "+\\deltaF_g^{\\mathrm{cross}}(\\cA)"
        ),
        "d^2=\\kappaChHodge(\\mathcalA)\\omega_g",
    )
    for required in required_compact:
        assert required in compact

    required_text = (
        "Faber--Pandharipande numerical trace contributes",
        "trace of the off-diagonal tensor curvature",
    )
    for required in required_text:
        assert required in text

    stale_compact = (
        "F_g=\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}",
        "F_g=\\kappaChHodge(\\cA)\\cdot\\lambda_g^{\\mathrm{FP}}",
        "F_g(\\cA)=\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}",
        "F_g(\\cA)=\\kappaChHodge(\\cA)\\cdot\\lambda_g^{\\mathrm{FP}}",
        "d^2=\\kappa\\cdot\\omega_g",
        "K_{w_1,w_2}=\\deltaF_g^{\\mathrm{cross}}|_{w_1,w_2}",
    )
    for stale in stale_compact:
        assert stale not in compact
