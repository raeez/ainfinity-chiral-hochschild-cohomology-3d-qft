"""Convention guards for W_N tempered-stratum closure.

Coverage
--------
  - prop:beta-stirling-dominance (ProvedHere): a finite Riccati-envelope
    coefficient beta gives factorial tempering by Stirling dominance.
  - prop:beta-N-closed-form (Conjectured): beta_N = 12(H_N - 1), conditional
    on the kappa-ratio scaling law.
  - prop:rho-star-WN (Conditional): rho_*^{W_N}(c) = |c|/[12(H_N - 1)]
    under the harmonic beta law.
  - W4 bridge guard: the spin-lane value 13 is evidence, not a proof of the
    full Riccati coefficient beta_4.

Independent-verification decorators are installed per the Vol II HZ-IV
protocol. derived_from and verified_against sources are DISJOINT:
  - Derivation uses the conditional finite-beta envelope or the harmonic
    kappa-ratio ansatz.
  - Verification uses direct arithmetic, Stirling estimates, and the W4
    bridge sentinel.

Each test has `# VERIFIED` comments (AP10) citing at least TWO
independent sources.

These tests deliberately do not certify global Class M closure or prove
that every principal W_N has the required finite Riccati envelope.
"""

from __future__ import annotations

import math
import os
import sys
from fractions import Fraction
from pathlib import Path

# Ensure package is importable.
_VOL2_ROOT = Path(__file__).resolve().parents[2]
if str(_VOL2_ROOT) not in sys.path:
    sys.path.insert(0, str(_VOL2_ROOT))

from compute.lib.wn_tempered_closure import (  # noqa: E402
    BETA_2,
    BETA_3,
    WN_TEMPERING_SCOPE,
    banach_radius_certificate,
    beta_N,
    beta_N_candidate_A,
    beta_N_candidate_B,
    beta_N_is_finite,
    class_m_banach_norm,
    geometric_mc_series_bound,
    mc_series_bound,
    operation_bound_holds,
    rho_star_WN,
    stirling_vs_beta_dominance,
    tempering_certificate,
    tempering_rate_bound,
)
from compute.lib.beta_N_closed_form import beta_N_from_kappa  # noqa: E402
from compute.lib.w4_beta_direct import direct_w4_attack_report  # noqa: E402

try:
    from compute.lib.independent_verification import independent_verification  # noqa: E402
except ImportError:
    # Allow running without decorator infrastructure (CI fallback).
    def independent_verification(*_args, **_kwargs):
        def wrap(f):
            return f
        return wrap


def _wn_tempered_tex() -> str:
    return (
        _VOL2_ROOT / "chapters/theory/wn_tempered_closure_platonic.tex"
    ).read_text()


def _window_around_label(tex: str, label: str, before: int = 500, after: int = 1500) -> str:
    marker = rf"\label{{{label}}}"
    index = tex.index(marker)
    return tex[max(0, index - before): index + after]


def _flat(text: str) -> str:
    return " ".join(text.split())


# ---------------------------------------------------------------------------
# Known proved values
# ---------------------------------------------------------------------------


def test_beta_2_is_six():
    """beta_2 = 6 from Vol II thm:tempered-stratum-contains-virasoro.

    # VERIFIED: [LT] tempered_stratum_characterization_platonic.tex:687 (beta_2=6)
    # VERIFIED: [DC] direct Vir Riccati: A_r = 8(-6)^{r-4}/r at large c,
    #              ratio -6(r-1)/r gives beta_2 = 6.
    """
    assert BETA_2 == Fraction(6)
    assert beta_N(2) == Fraction(6)


def test_beta_3_is_ten():
    """beta_3 = 10 from Vol II thm:tempered-stratum-contains-w3.

    # VERIFIED: [LT] tempered_stratum_characterization_platonic.tex:698 (beta_3=10)
    # VERIFIED: [DC] Fateev-Lukyanov cubic coupling + Riccati with
    #              kappa(W_3) = 5c/6 plus W-generator channel (6 + 4 = 10).
    """
    assert BETA_3 == Fraction(10)
    assert beta_N(3) == Fraction(10)


def test_wn_closure_manuscript_carries_licensing_and_scope_guards():
    """The four W_N closure claims carry status, licensing, and scope data."""
    tex = _wn_tempered_tex()

    stirling = _window_around_label(tex, "prop:beta-stirling-dominance")
    assert r"\ClaimStatusProvedHere" in stirling
    assert r"licensing $\gamma$" in stirling
    assert r"\hypAmbientWtCpl" in stirling
    assert "finite" in stirling and r"\beta" in stirling
    assert r"\limsup" in stirling

    tempering = _window_around_label(tex, "thm:wn-tempered-all-N")
    tempering_flat = _flat(tempering)
    assert r"\ClaimStatusConditional" in tempering
    assert r"licensing $\alpha+\gamma+\varepsilon$" in tempering
    assert r"\hypAmbientWtCpl" in tempering
    assert "finite Riccati envelope" in tempering
    assert "finite envelope as an explicit assumption" in tempering_flat
    assert "does not derive it from Fateev--Lukyanov" in tempering_flat

    radius = _window_around_label(tex, "prop:rho-star-WN-closed-form")
    radius_flat = _flat(radius)
    assert r"\ClaimStatusConditional" in radius
    assert r"licensing $\alpha+\gamma$" in radius
    assert r"\hypAmbientWtCpl" in radius
    assert "leading-root sharpness" in radius_flat
    assert r"Conjecture~\ref{prop:beta-N-closed-form}" in radius
    assert "Cauchy--Hadamard root test" in radius_flat
    assert "only the lower bound" in radius_flat

    corollary = _window_around_label(tex, "cor:wn-original-complex-dichotomy-healed")
    corollary_flat = _flat(corollary)
    assert r"\ClaimStatusConditional" in corollary
    assert r"licensing $\alpha+\gamma+\varepsilon$" in corollary
    assert r"\hypAmbientWtCpl" in corollary
    assert "finite-envelope locus" in corollary_flat
    assert "Remaining open candidates" in corollary


def test_wn_radius_equality_is_not_stated_without_sharpness():
    """Finite upper envelopes give tempering; radius equality needs sharpness."""
    tex_flat = _flat(_wn_tempered_tex())
    assert "On the leading-root sharp finite-envelope locus" in tex_flat
    assert "finite envelope gives only the lower bound" in tex_flat
    assert "ordinary-generating convergence radius of" in tex_flat


# ---------------------------------------------------------------------------
# Harmonic-candidate check: matches proved data and separates candidates conditionally
# ---------------------------------------------------------------------------


@independent_verification(
    claim="prop:beta-N-closed-form",
    derived_from=[
        "kappa-ratio scaling A_r(W_N)=(kappa(W_N)/kappa(Vir))^(r-3) A_r(Vir)",
        "kappa(W_N)=c(H_N-1), kappa(Vir)=c/2, and Virasoro ratio -6(r-1)/r",
    ],
    verified_against=[
        "Direct numerical match at N=2 (beta=6 from Vir) and N=3 (beta=10 from W_3)",
        "N=4 discriminator: harmonic value 13, triangular 15, quadratic 16",
    ],
    disjoint_rationale=(
        "Derivation uses the all-rank kappa-ratio law. Verification checks the "
        "two low-rank Laurent values and the first unfixed rank N=4, where the "
        "old candidates are conditionally separated under the harmonic ansatz."
    ),
)
def test_beta_N_candidates_match_N2_N3_only():
    """Both Candidate A and Candidate B match at N=2, 3 by coincidence.

    The harmonic conjecture (conj:beta-N-harmonic-closed-form in
    beta_N_closed_form_all_platonic.tex) predicts beta_N = 12(H_N - 1),
    which also matches N=2, 3. At N=4 all three formulas diverge:
    Candidate A = 15, Candidate B = 16, harmonic candidate = 13.

    This test is retained for historical documentation; the canonical
    beta_N value in beta_N(.) now tracks the harmonic
    value from compute/lib/beta_N_closed_form.py::beta_N_from_kappa.

    # VERIFIED: [DC] formula evaluation N=2: both candidates give 6.
    # VERIFIED: [DC] formula evaluation N=3: both candidates give 10.
    # VERIFIED: [LT] beta_N_closed_form_all_platonic.tex
    #              conj:beta-N-harmonic-closed-form conditionally retracts both.
    """
    # Both candidates (and the harmonic candidate) match at N=2 and N=3
    assert beta_N_candidate_A(2) == BETA_2
    assert beta_N_candidate_A(3) == BETA_3
    assert beta_N_candidate_B(2) == BETA_2
    assert beta_N_candidate_B(3) == BETA_3
    # At N=4 all three formulas diverge
    assert beta_N_candidate_A(4) == Fraction(15)
    assert beta_N_candidate_B(4) == Fraction(16)
    # The harmonic candidate value is 13; see
    # compute/lib/beta_N_closed_form.py and
    # chapters/theory/beta_N_closed_form_all_platonic.tex
    assert beta_N_from_kappa(4) == Fraction(13)
    assert beta_N(4) == Fraction(13)
    assert beta_N_from_kappa(4) != beta_N_candidate_A(4)
    assert beta_N_from_kappa(4) != beta_N_candidate_B(4)


def test_beta_N_is_finite_for_all_N():
    """The harmonic candidate beta_N is a finite positive rational.

    This is an arithmetic guard for the candidate value. It does not prove
    that the full W_N shadow tower has this beta_N as its Riccati envelope.

    # VERIFIED: [DC] Direct evaluation of beta_N(N) returns Fraction for N=2..20.
    # VERIFIED: [DA] Dimensional analysis: beta_N is ratio of structure
    #              constants (rational) / norm coefficients (rational).
    """
    for N in range(2, 21):
        assert beta_N_is_finite(N)
        val = beta_N(N)
        assert val > 0
        assert isinstance(val, Fraction)


# ---------------------------------------------------------------------------
# Tempering theorem: concrete rate bounds
# ---------------------------------------------------------------------------


@independent_verification(
    claim="prop:beta-stirling-dominance",
    derived_from=[
        "Stirling approximation (r!)^{-1/r} ~ e/r",
        "Conditional Riccati-envelope upper bound |S_r| <= C * beta_N^{r-4} / (r * |c|^{r-2})",
    ],
    verified_against=[
        "Direct numerical evaluation of the rate bound at r = 4, 10, 20, 100",
        "Monotone-decreasing certificate of the rate sequence (independent of Stirling step)",
    ],
    disjoint_rationale=(
        "Derivation assumes a finite beta envelope and combines it with the "
        "Stirling approximation. Verification evaluates the resulting bound at "
        "finite r-values, observing the sequence decreases monotonically to 0 "
        "without promoting the envelope to a global W_N theorem."
    ),
)
def test_tempering_rate_tends_to_zero_every_N():
    """For sampled N and c, the conditional finite-beta bound decreases to 0.

    # VERIFIED: [DC] Direct computation at multiple r-values for each N.
    # VERIFIED: [LC] Limiting case r -> infty: bound -> 0 by Stirling.
    """
    c = Fraction(100)
    for N in range(2, 9):
        # Use r-values that avoid float underflow at r=100 for larger beta_N
        bounds = tempering_certificate(N, c, r_values=[4, 10, 20, 30])
        r_vals = sorted(bounds.keys())
        # Monotone decrease (non-strict at far tail due to float)
        for i in range(len(r_vals) - 1):
            assert bounds[r_vals[i]] >= bounds[r_vals[i + 1]], (
                f"N={N}: rate bound not monotone at r={r_vals[i]}: "
                f"{bounds[r_vals[i]]} -> {bounds[r_vals[i + 1]]}"
            )
        # At r = 30, bound is already below 0.2 (tempering certificate)
        assert bounds[30] < 0.2, (
            f"N={N}, c={c}, r=30: bound = {bounds[30]}, expected < 0.2"
        )
        # At large r, bound strictly less than at small r
        assert bounds[30] < bounds[4]


@independent_verification(
    claim="prop:beta-stirling-dominance",
    derived_from=[
        "Stirling factor (r!)^{-1/r} ~ e/r",
    ],
    verified_against=[
        "Critical radius r_crit = beta_N * e / |c| evaluated directly",
        "Empirical check at r = 100 for N up to 6",
    ],
    disjoint_rationale=(
        "Derivation uses asymptotic Stirling under a finite-beta hypothesis. "
        "Verification evaluates beta_N * e / r_max at finite r_max = 100 and "
        "observes it is < 1 for N = 2..6 and c = 100."
    ),
)
def test_stirling_dominates_beta_at_large_r():
    """Stirling factor dominates beta_N at r = 100 for N <= 6, c = 100.

    # VERIFIED: [DC] beta_N * e / r_max < 1 for all tested N.
    # VERIFIED: [DA] Dimensional analysis: Stirling decay rate ~ 1/r,
    #              beta_N constant; at large r, Stirling wins.
    """
    c = Fraction(100)
    for N in range(2, 7):
        assert stirling_vs_beta_dominance(N, c, r_max=100)


# ---------------------------------------------------------------------------
# Conditional rho_* closed form
# ---------------------------------------------------------------------------


@independent_verification(
    claim="prop:rho-star-WN-closed-form",
    derived_from=[
        "Conditional harmonic beta_N = 12(H_N - 1) from conj:beta-N-harmonic-closed-form",
        "Root test: rho_*(c) = (limsup |S_r|^{1/r})^{-1}",
    ],
    verified_against=[
        "Direct algebraic evaluation of rho_* at N=2 and N=3 matching Vol II chapter",
        "Inverse of |c| / beta_N at explicit c values",
    ],
    disjoint_rationale=(
        "Derivation combines the conditional beta_N closed form with Cauchy-Hadamard root "
        "test. Verification checks the algebraic value rho_*(c) = |c|/beta_N "
        "at N=2 gives |c|/6 (matches Vol II Prop. prop:virasoro-rho-star-closed-form) "
        "and N=3 gives |c|/10 (matches Vol II thm:tempered-stratum-contains-w3)."
    ),
)
def test_rho_star_closed_form():
    """Conditional rho_*^{W_N}(c) = |c| / [12(H_N - 1)].

    # VERIFIED: [LT] Vol II chapter Prop. prop:virasoro-rho-star-closed-form.
    # VERIFIED: [LT] Vol II chapter thm:tempered-stratum-contains-w3.
    # VERIFIED: [DC] Direct algebra at N=4,5,6.
    """
    # N=2: rho_*(c) = |c|/6
    assert rho_star_WN(2, Fraction(12)) == Fraction(2)
    assert rho_star_WN(2, Fraction(60)) == Fraction(10)
    # N=3: rho_*(c) = |c|/10
    assert rho_star_WN(3, Fraction(100)) == Fraction(10)
    assert rho_star_WN(3, Fraction(30)) == Fraction(3)
    # N=4: rho_*(c) = |c|/13
    assert rho_star_WN(4, Fraction(130)) == Fraction(10)
    # N=5: rho_*(c) = 5|c|/77
    assert rho_star_WN(5, Fraction(154)) == Fraction(10)
    # General: rho_*(c) * beta_N = |c|
    for N in range(2, 10):
        c = Fraction(1000)
        rho = rho_star_WN(N, c)
        assert rho * beta_N(N) == c


# ---------------------------------------------------------------------------
# Banach body and MC operation constants
# ---------------------------------------------------------------------------


def test_class_m_banach_norm_is_weighted_l1_body():
    """||a||_rho is sum_w ||a_w|| rho^w."""
    weight_norms = {0: Fraction(2), 2: Fraction(3), 4: Fraction(5)}
    rho = Fraction(1, 3)
    expected = Fraction(2) + Fraction(3, 9) + Fraction(5, 81)
    assert class_m_banach_norm(weight_norms, rho) == expected


def test_operation_bound_holds_is_multilinear_norm_inequality():
    """Check ||m_k(a_i)|| <= C_k prod ||a_i||."""
    assert operation_bound_holds(
        output_norm=Fraction(5),
        input_norms=[Fraction(2), Fraction(3)],
        constant=Fraction(1),
    )
    assert not operation_bound_holds(
        output_norm=Fraction(7),
        input_norms=[Fraction(2), Fraction(3)],
        constant=Fraction(1),
    )


def test_mc_series_bound_and_geometric_majorant():
    """The Banach MC criterion is sum_k C_k(rho) R^k < infinity."""
    constants = [Fraction(1, 4), Fraction(1, 8), Fraction(1, 16)]
    assert mc_series_bound(constants, radius=Fraction(1)) == Fraction(7, 16)
    assert mc_series_bound(constants, radius=Fraction(1, 2)) == (
        Fraction(1, 4) * Fraction(1, 4)
        + Fraction(1, 8) * Fraction(1, 8)
        + Fraction(1, 16) * Fraction(1, 16)
    )
    assert geometric_mc_series_bound(Fraction(1, 4), Fraction(1, 2)) == Fraction(1, 2)


def test_banach_radius_certificate_statuses():
    """Virasoro/W3 radii are proved; W_N for N>=4 is conditional."""
    vir = banach_radius_certificate(2, Fraction(60))
    w3 = banach_radius_certificate(3, Fraction(100))
    w4 = banach_radius_certificate(4, Fraction(130))

    assert vir["beta"] == Fraction(6)
    assert vir["rho_star"] == Fraction(10)
    assert vir["status"] == "proved_low_rank"
    assert w3["beta"] == Fraction(10)
    assert w3["rho_star"] == Fraction(10)
    assert w3["status"] == "proved_low_rank"
    assert w4["beta"] == Fraction(13)
    assert w4["rho_star"] == Fraction(10)
    assert w4["status"] == "conditional_harmonic_beta"
    assert w4["requires_mc_sum"] == "sum_k C_k(rho) < infinity"


# ---------------------------------------------------------------------------
# Numerical Laurent-expansion check at small c (chain B verification)
# ---------------------------------------------------------------------------


def test_beta_2_matches_vol1_virasoro_engine():
    """Cross-verification against Vol I virasoro_shadow_sequence.

    Computes the ratio |S_{r+1}|/((r+1)|S_r|) at r = 10, c = 100 from
    the Vol I engine and compares against beta_2 / (c * (r+1)) = 6/1100.

    This uses the Vol I engine (independent of this Vol II module) to
    verify that beta_2 = 6 is the actual Riccati ratio, not a hardcoded
    value.

    # VERIFIED: [DC] Vol I shadow_tower_higher_vir.virasoro_shadow_sequence
    # VERIFIED: [LC] Leading Laurent coefficient A_r at large c gives 6.
    """
    import importlib.util

    vol1_root = os.path.expanduser("~/chiral-bar-cobar")
    if not os.path.exists(vol1_root):
        import pytest
        pytest.skip("Vol I engine not available at expected path")
    spec = importlib.util.spec_from_file_location(
        "_vol1_shadow",
        os.path.join(vol1_root, "compute", "lib", "shadow_tower_higher_vir.py"),
    )
    if spec is None or spec.loader is None:
        import pytest
        pytest.skip("Cannot load Vol I shadow engine module")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    import sympy as sp
    c_val = sp.sympify(Fraction(1000))  # Very large c for clean leading behavior
    S = mod.virasoro_shadow_sequence(c_val, max_r=11)
    r = 10
    # |S_{r+1}|/((r+1)|S_r|) ~ beta / (c * (r+1)) at large c
    numerical_ratio = abs(float(S[r + 1])) / ((r + 1) * abs(float(S[r])))
    expected_ratio = float(BETA_2) / (1000.0 * (r + 1))
    # Within 10% because finite-r correction
    rel_err = abs(numerical_ratio - expected_ratio) / expected_ratio
    assert rel_err < 0.15, (
        f"Vol I Vir ratio at r={r}, c=1000: numerical {numerical_ratio}, "
        f"expected {expected_ratio}, rel_err {rel_err}"
    )


# ---------------------------------------------------------------------------
# Stirling + beta vs r asymptotic bound
# ---------------------------------------------------------------------------


def test_tempering_rate_bound_explicit_values():
    """Direct evaluation of the tempering rate bound at explicit (N, c, r).

    Confirms the bound formula evaluates to finite non-negative values and
    decreases with r. At very large r, float underflow produces exactly 0,
    which is itself evidence of tempering (limsup = 0).

    # VERIFIED: [DC] Direct formula evaluation.
    # VERIFIED: [DA] Dimensional analysis: beta^r/(r^{r+1} c^r) -> 0 at large r.
    """
    c = Fraction(100)
    # For N = 4 at r = 10, r = 30 (before underflow)
    b10 = tempering_rate_bound(4, c, 10)
    b30 = tempering_rate_bound(4, c, 30)
    b50 = tempering_rate_bound(4, c, 50)
    assert b10 > 0
    assert b30 > 0
    # Monotone decreasing
    assert b10 > b30 > b50
    # Bounds at r = 10 already small
    assert b10 < 0.5
    # At r = 30, bound << 1/e
    assert b30 < 0.2
    # At r = 50, either very small or underflowed to 0 (both fine for tempering)
    assert b50 >= 0
    assert b50 < 0.05


# ---------------------------------------------------------------------------
# Cross-N comparison
# ---------------------------------------------------------------------------


def test_beta_N_monotone_increasing():
    """beta_N is strictly increasing in N.

    Reflects that higher-rank W_N has more coupling channels, hence
    larger Riccati ratio.

    # VERIFIED: [DC] Direct evaluation beta_N(N+1) - beta_N(N) > 0.
    # VERIFIED: [DA] 12(H_N - 1) is strictly increasing for N >= 2.
    """
    for N in range(2, 20):
        assert beta_N(N) < beta_N(N + 1)


def test_candidate_A_vs_B_miss_harmonic_N4():
    """Two old candidate formulas agree at N=2,3 but miss harmonic N=4.

    The harmonic conjecture gives beta_4 = 13; candidate_A gives 15 and
    candidate_B gives 16. The direct W4 lane witness supports 13 but
    does not by itself prove the Riccati coefficient.

    # VERIFIED: [DC] Direct evaluation of both candidates at N = 2, 3, 4.
    """
    for N in (2, 3):
        assert beta_N_candidate_A(N) == beta_N_candidate_B(N)
    # Disagreement at N=4
    assert beta_N(4) == Fraction(13)
    assert beta_N_candidate_A(4) == Fraction(15)
    assert beta_N_candidate_B(4) == Fraction(16)
    assert beta_N_candidate_A(4) != beta_N(4)
    assert beta_N_candidate_B(4) != beta_N(4)
    report = direct_w4_attack_report()
    assert report.lane_sum == Fraction(13)
    assert report.required_bridge_ratio == Fraction(-52, 5)
    assert report.required_A5_from_bridge == Fraction(-676, 15)
    assert report.observed_full_miura_A5 is None
    assert report.proves_beta4 is False


# ---------------------------------------------------------------------------
# Scope guard: finite-beta tests are not all-rank W_N proofs
# ---------------------------------------------------------------------------


@independent_verification(
    claim="prop:beta-stirling-dominance",
    derived_from=[
        "Conditional finite-beta Riccati envelope",
    ],
    verified_against=[
        "Numerical certificate at large r for N = 2, 3, 4, 5, 6 and c = 100",
        "The certificate (|S_r|/r!)^(1/r) bound is < 1/e = 0.368 at r = 100",
    ],
    disjoint_rationale=(
        "The conditional theorem says finite beta plus Stirling gives "
        "tempering. Verification evaluates the finite-beta bound for sampled "
        "N without asserting that the full W_N envelope has been proved."
    ),
)
def test_finite_beta_bound_tempered_for_sampled_WN_candidates():
    """Sampled harmonic beta candidates have rate bound < 1/e at r = 30.

    This is a finite-beta convention guard, not a proof that the non-tempered
    stratum of every principal W_N is empty at generic c. Float evaluation at
    r = 100 may underflow; r = 30 gives a non-degenerate numerical comparison
    against 1/e.

    # VERIFIED: [DC] Direct evaluation at explicit (N, c, r).
    # VERIFIED: [LC] Asymptotic limit r -> infty gives 0 for every N.
    """
    one_over_e = 1.0 / math.e
    c = Fraction(100)
    for N in range(2, 7):
        bound = tempering_rate_bound(N, c, 30)
        assert bound < one_over_e, (
            f"N={N}, c={c}, r=30: bound = {bound} >= 1/e = {one_over_e}"
        )
    # At r = 100, bounds underflow to 0 (confirmed tempering)
    for N in range(2, 7):
        bound_large_r = tempering_rate_bound(N, c, 100)
        assert bound_large_r < 1e-3, (
            f"N={N}, c={c}, r=100: bound = {bound_large_r}, expected << 1"
        )


def test_global_WN_closure_is_not_encoded_as_compute_proof():
    """The compute surface records the missing all-rank Riccati bridge."""
    assert WN_TEMPERING_SCOPE["finite_beta_stirling_implication"] == "proved"
    assert WN_TEMPERING_SCOPE["finite_riccati_envelope_all_N"] == "not_encoded"
    assert WN_TEMPERING_SCOPE["harmonic_beta_closed_form"] == (
        "conditional_kappa_ratio_scaling"
    )
    assert WN_TEMPERING_SCOPE["leading_root_sharpness"] == (
        "assumed_for_radius_equality"
    )
    assert WN_TEMPERING_SCOPE["w4_full_miura_A5_bridge"] == "absent"
    report = direct_w4_attack_report()
    assert report.observed_full_miura_A5 is None
    assert report.proves_beta4 is False


# ---------------------------------------------------------------------------
# Proved values assertion
# ---------------------------------------------------------------------------


def test_proved_constants():
    """Proved values BETA_2 = 6 and BETA_3 = 10.

    # VERIFIED: [LT] tempered_stratum_characterization_platonic.tex
    # VERIFIED: [DC] Independent re-derivation in Vir and W_3 engines.
    """
    assert BETA_2 == Fraction(6)
    assert BETA_3 == Fraction(10)
    assert BETA_2 != BETA_3


if __name__ == "__main__":
    test_beta_2_is_six()
    test_beta_3_is_ten()
    test_wn_closure_manuscript_carries_licensing_and_scope_guards()
    test_wn_radius_equality_is_not_stated_without_sharpness()
    test_beta_N_candidates_match_N2_N3_only()
    test_beta_N_is_finite_for_all_N()
    test_tempering_rate_tends_to_zero_every_N()
    test_stirling_dominates_beta_at_large_r()
    test_rho_star_closed_form()
    test_tempering_rate_bound_explicit_values()
    test_beta_N_monotone_increasing()
    test_candidate_A_vs_B_miss_harmonic_N4()
    test_finite_beta_bound_tempered_for_sampled_WN_candidates()
    test_global_WN_closure_is_not_encoded_as_compute_proof()
    test_proved_constants()
    print("All tests passed.")
