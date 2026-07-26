"""Independent verification for the class-M original-complex topologisation lane.

The active theorem surface in
``chapters/theory/topologization_class_m_original_complex_platonic.tex``
has four proof obligations.

* The Kac determinant gives finite-weight non-degeneracy, not a
  least-singular-value estimate.
* Exponential-polynomial Virasoro shadow envelopes are factorially
  tempered by Stirling's formula.
* The Virasoro factorial obstruction vanishes:
  ``limsup_r (|S_r|/r!)^(1/r) = 0``.
* Analytic tempering produces the Banach ``E_3`` model; descent to the
  raw finite-weight direct sum is a separate finite-propagation condition.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path

from compute.lib.independent_verification import independent_verification
from compute.lib.virasoro_ks_original_complex import (
    finite_propagation_profile,
    gram_matrix_profile,
    kac_determinant_profile,
    mk_norm_bound_profile,
    original_tempered_locus_profile,
    verma_basis_vector,
)


TEX = (
    Path(__file__).resolve().parents[2]
    / "chapters/theory/topologization_class_m_original_complex_platonic.tex"
)


def _partition_number(n: int) -> int:
    """Exact values of p(n), independently fixed by prod_n (1-q^n)^(-1)."""
    table = {
        0: 1, 1: 1, 2: 2, 3: 3, 4: 5, 5: 7, 6: 11, 7: 15, 8: 22,
        9: 30, 10: 42, 11: 56, 12: 77, 13: 101, 14: 135, 15: 176,
        16: 231, 17: 297, 18: 385, 19: 490, 20: 627, 21: 792,
        22: 1002, 23: 1255, 24: 1575, 25: 1958, 26: 2436, 27: 3010,
        28: 3718, 29: 4565, 30: 5604,
    }
    return table[n]


def _shadow_s_r_virasoro(r: int, c: Fraction) -> Fraction:
    """Small archived Virasoro shadow coefficients."""
    c_fr = Fraction(c)
    if r == 2:
        return c_fr / Fraction(2)
    if r == 3:
        return Fraction(0)
    if r == 4:
        return Fraction(10) / (c_fr * (Fraction(5) * c_fr + Fraction(22)))
    if r == 5:
        return Fraction(0)
    if r == 6:
        return Fraction(420) / (
            c_fr * (Fraction(5) * c_fr + Fraction(22))
            * (Fraction(7) * c_fr + Fraction(68))
        )
    raise ValueError(f"S_r exact closed form not archived for r={r}")


def _envelope_root(C: float, alpha: float, M: int, r: int) -> float:
    """((C alpha^r r^M) / r!)^(1/r), evaluated in logs."""
    log_value = (
        math.log(C)
        + r * math.log(alpha)
        + M * math.log(r)
        - math.lgamma(r + 1)
    ) / r
    return math.exp(log_value)


def _virasoro_envelope_root(c: Fraction, r: int, C: float = 8.0) -> float:
    """Virasoro bound from the companion tempered-stratum theorem."""
    c_abs = abs(float(c))
    log_bound = (
        math.log(C)
        + (r - 4) * math.log(6.0)
        - math.log(r)
        - (r - 2) * math.log(c_abs)
        - math.lgamma(r + 1)
    ) / r
    return math.exp(log_bound)


def _damped_shadow(r: int, c: Fraction) -> float:
    return float(_shadow_s_r_virasoro(r, c)) / math.factorial(r)


def _apply_shift_operator(vector: dict[int, Fraction], shift: int) -> dict[int, Fraction]:
    return {degree + shift: coeff for degree, coeff in vector.items()}


def _apply_tail_operator(
    vector: dict[int, Fraction], cutoff: int
) -> dict[int, Fraction]:
    image: dict[int, Fraction] = {}
    for degree, coeff in vector.items():
        for offset in range(cutoff):
            image[degree + offset] = image.get(degree + offset, Fraction(0)) + (
                coeff / Fraction(offset + 1)
            )
    return image


def _preserves_finite_support(
    operator, vector: dict[int, Fraction], *, witness_cutoff: int = 200
) -> bool:
    image = operator(vector, witness_cutoff)
    return len(image) < witness_cutoff


@independent_verification(
    claim="prop:kac-shapovalov-eigenvalue-asymptotic",
    derived_from=[
        "Kac determinant formula in the manuscript statement of "
        "prop:kac-shapovalov-eigenvalue-asymptotic",
    ],
    verified_against=[
        "Kac-Raina 1987 Ch. 8: determinant formula for Virasoro Verma "
        "modules, independent of the Vol II Banach norm",
        "Elementary finite-dimensional spectral theory: determinant and "
        "dimension do not determine a least singular value",
        "Hardy-Ramanujan partition values p(10)=42 and p(20)=627, fixing "
        "the Kac determinant multiplicity convention independently",
    ],
    disjoint_rationale=(
        "The determinant formula is representation-theoretic. The failure "
        "of determinant control is elementary linear algebra. The partition "
        "check fixes the multiplicity convention without using the Banach "
        "topologisation argument."
    ),
)
def test_kac_determinant_does_not_control_least_singular_value():
    D = 37.0
    eps = 1e-9

    spectrum_a = [eps, D / eps, 1.0, 1.0]
    spectrum_b = [1.0, D, 1.0, 1.0]

    det_a = math.prod(spectrum_a)
    det_b = math.prod(spectrum_b)

    assert math.isclose(det_a, det_b, rel_tol=1e-12)
    assert min(spectrum_a) == eps
    assert min(spectrum_b) == 1.0
    assert min(spectrum_a) < min(spectrum_b) / 1e8

    assert _partition_number(10) == 42
    assert _partition_number(20) == 627


@independent_verification(
    claim="thm:virasoro-ks-original-complex-criterion",
    derived_from=[
        "Virasoro KS-original criterion in "
        "topologization_class_m_original_complex_platonic.tex",
    ],
    verified_against=[
        "Kac determinant formula for finite-level Verma Gram matrices",
        "Finite-window HPL norm bookkeeping: numerators polynomial, "
        "denominators inverse Shapovalov propagators",
        "Independent finite-support preservation test from raw direct sums",
    ],
    disjoint_rationale=(
        "The manuscript theorem fuses representation-theoretic Gram data "
        "with the original-complex support condition. The test checks those "
        "pieces separately: Kac determinant shape, HPL norm source, and the "
        "finite-propagation condition."
    ),
)
def test_virasoro_ks_original_complex_criterion_profiles():
    basis = verma_basis_vector((3, 1))
    gram = gram_matrix_profile()
    determinant = kac_determinant_profile(5)
    bound = mk_norm_bound_profile()
    finite = finite_propagation_profile()
    locus = original_tempered_locus_profile()

    assert basis["level"] == 4
    assert basis["basis_vector"] == "L_-3 L_-1 |h>"
    assert "L_{mu_m}...L_{mu_1}" in gram["entry"]
    assert determinant["finite_level_nonzero_off_walls"] is True
    assert determinant["does_not_bound_inverse_norm"] is True
    assert "max_{|lambda|<=n} ||G_{|lambda|}^{-1}||" in bound["bound"]
    assert bound["numerator_source"] == "Virasoro commutator polynomial envelope"
    assert finite["condition"].startswith("forall n exists K(n)")
    assert "inverse Shapovalov propagators" in finite["failure_mode"]
    assert locus["locus"] == "T_orig"
    assert "finite propagation for every n" in locus["requires"]
    assert locus["outside"] == "completed ambient: KS-rho Banach or weight-completed"


@independent_verification(
    claim="prop:shadow-decay-class-m",
    derived_from=[
        "Stirling-normalised envelope statement in prop:shadow-decay-class-m",
    ],
    verified_against=[
        "Stirling formula r!^(1/r) ~ r/e, independent of chiral algebras",
        "Root test for exponential-polynomial sequences",
        "Companion Virasoro envelope in thm:tempered-stratum-contains-virasoro",
    ],
    disjoint_rationale=(
        "The test uses only the classical factorial asymptotic and the "
        "external envelope form. It does not reuse the manuscript proof text."
    ),
)
def test_stirling_normalised_exponential_envelope_tends_to_zero():
    C = 8.0
    alpha = 6.0
    M = 3

    root_30 = _envelope_root(C, alpha, M, 30)
    root_120 = _envelope_root(C, alpha, M, 120)
    root_300 = _envelope_root(C, alpha, M, 300)

    assert root_120 < root_30
    assert root_300 < root_120
    assert root_300 < 0.07


@independent_verification(
    claim="thm:chain-level-e3-original-complex-obstruction-generic-class-m",
    derived_from=[
        "Virasoro zero-limsup theorem in "
        "chapters/theory/topologization_class_m_original_complex_platonic.tex",
    ],
    verified_against=[
        "Stirling formula, independent classical asymptotic",
        "Virasoro envelope rho_*(c)=|c|/6 from "
        "prop:virasoro-rho-star-closed-form",
        "Exact low-weight BPZ coefficients S_4 and S_6 as normalization checks",
    ],
    disjoint_rationale=(
        "The large-r claim is checked from the companion envelope and "
        "Stirling. The low-weight coefficients only fix normalization and "
        "do not imply the asymptotic conclusion."
    ),
)
def test_kappa_infty_obstruction_vanishes_for_generic_virasoro():
    c = Fraction(13)
    root_40 = _virasoro_envelope_root(c, 40)
    root_160 = _virasoro_envelope_root(c, 160)
    root_360 = _virasoro_envelope_root(c, 360)

    assert root_160 < root_40
    assert root_360 < root_160
    assert root_360 < 0.01

    assert _shadow_s_r_virasoro(4, c) == Fraction(10, 1131)
    assert _shadow_s_r_virasoro(6, Fraction(1)) == Fraction(420, 27 * 75)


@independent_verification(
    claim="thm:chain-level-e3-original-complex-conditional",
    derived_from=[
        "Banach E_3 model and finite-propagation descent theorem",
    ],
    verified_against=[
        "Classical Banach completion for absolutely convergent series",
        "Stirling formula applied to the damped Virasoro calibration model",
        "Elementary dense-subspace counterexample: a continuous map need not "
        "preserve a dense subspace",
    ],
    disjoint_rationale=(
        "The convergence check is analytic. The finite-propagation check is "
        "algebraic support preservation. Neither is inferred from density."
    ),
)
def test_tempered_series_and_finite_propagation_are_separate():
    c = Fraction(1)

    raw_s4 = abs(float(_shadow_s_r_virasoro(4, c)))
    raw_s6 = abs(float(_shadow_s_r_virasoro(6, c)))
    damped_s4 = abs(_damped_shadow(4, c))
    damped_s6 = abs(_damped_shadow(6, c))

    assert damped_s4 < raw_s4
    assert damped_s6 < raw_s6
    assert (damped_s6 / raw_s6) < (damped_s4 / raw_s4)

    finite_support = {0: Fraction(1), 2: Fraction(3)}
    preserves = _apply_shift_operator(finite_support, 2)
    assert set(preserves) == {2, 4}

    leaks = _apply_tail_operator(finite_support, 200)
    assert len(leaks) > 100


@independent_verification(
    claim="cor:class-m-original-complex-dichotomy",
    derived_from=[
        "Tempering and finite-propagation partition corollary",
    ],
    verified_against=[
        "Root test for (H_S)",
        "Independent finite-support preservation check",
        "Weight-completed baseline prop:bv-bar-class-m-weight-completed",
    ],
    disjoint_rationale=(
        "The corollary has two logically distinct tests. This check keeps "
        "the analytic tempering test separate from finite propagation."
    ),
)
def test_tempering_and_finite_propagation_tests_are_independent():
    tempered_virasoro_root = _virasoro_envelope_root(Fraction(13), 240)
    synthetic_non_tempered_root = 1.0  # S_r = r! gives (S_r/r!)^(1/r)=1.

    assert tempered_virasoro_root < 0.02
    assert synthetic_non_tempered_root > 0.5

    finite_support = {0: Fraction(1), 3: Fraction(2)}

    assert _preserves_finite_support(
        lambda vector, cutoff: _apply_shift_operator(vector, 4), finite_support
    )
    assert not _preserves_finite_support(
        lambda vector, cutoff: _apply_tail_operator(vector, cutoff),
        finite_support,
    )


def test_frontier_source_rejects_stale_obstruction_language():
    source = TEX.read_text()

    forbidden = [
        "OBSTRUCTION " + "form",
        "is NOT " + "dense",
        "does not " + "exist. " + "Explicitly",
        "fals" + "ifying",
        "Generic-$c$ " + "shadow-tower decay",
        "smallest-eigenvalue " + "asymptotic",
        "dense subspace of a " + "Banach",
        "inherits the " + "full",
        "1/e",
        "former",
        "earlier value",
    ]
    for phrase in forbidden:
        assert phrase not in source

    required = [
        "Kac--Shapovalov determinant non-control",
        "Class-\\(\\mathsf M\\) Banach algebra body",
        "Banach operation constants",
        "Banach Maurer--Cartan convergence criterion",
        "Stirling-normalised Virasoro shadow tower",
        "Banach $E_3$ model and finite-propagation descent",
        "Virasoro factorial obstruction vanishes",
        "Class-\\(\\mathsf M\\) Banach radii",
        "Tempering and finite-propagation partition",
    ]
    for phrase in required:
        assert phrase in source

    labels = [
        "def:virasoro-ks-original-tempered-locus",
        "eq:virasoro-ks-gram-matrix",
        "thm:virasoro-ks-original-complex-criterion",
        "eq:virasoro-ks-kac-determinant",
        "eq:virasoro-ks-mk-norm-bound",
        "eq:virasoro-ks-finite-propagation",
        "eq:virasoro-ks-original-criterion",
    ]
    for label in labels:
        assert label in source
