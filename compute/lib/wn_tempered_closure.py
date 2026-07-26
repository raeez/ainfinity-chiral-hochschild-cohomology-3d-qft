r"""W_N tempered-stratum convention guards.

Context
-------
Vol II chapter chapters/theory/tempered_stratum_characterization_platonic.tex
proved that generic Vir_c (N=2) and generic W_{3,c} (N=3) are analytically
tempered with ordinary-generating convergence radius rho_*(c) = |c| / beta_N
for beta_2 = 6, beta_3 = 10. Extension to W_N for N >= 4 was stated as
conj:wn-tempered-general.

This module does not prove the all-rank W_N Riccati envelope. It records
the conditional arithmetic used by the manuscript:

* if a finite Riccati-envelope coefficient beta_N is available for a fixed
  principal W_N, Stirling dominance gives analytic tempering;
* the harmonic value

    beta_N = 12(H_N - 1) = sum_{s=2}^{N} 12/s,

  is the kappa-ratio candidate, not an independently encoded proof of the
  full W_N shadow recursion;
* the W4 bridge from full Miura/OPE data to A_5(W4) is tracked separately
  in compute/lib/w4_beta_direct.py and remains absent;
* the equality rho_* = |c| / beta_N requires leading-root sharpness of the
  ordinary shadow sequence. A finite upper envelope alone gives only the
  lower bound rho_* >= |c| / beta_N.

The finite-beta Stirling implication is:

    limsup_r (|S_r(W_{N,c})|/r!)^{1/r}
      <= limsup_r (C(c) * beta_N^{r-4} / (r * |c|^{r-2} * r!))^{1/r}
      = beta_N/|c| * 1 * 0  =  0.

The beta-free universal factor is (r!)^{-1/r} ~ e/r -> 0; every finite
beta_N is absorbed. This implication is a convention guard for the tests,
not a proof that every principal W_N has the required envelope.

Conditional harmonic candidate beta_N = 12(H_N - 1)
--------------------------------------------------
The leading-Laurent coefficient obeys the kappa-ratio scaling law

    A_r(W_N) = [kappa(W_N) / kappa(Vir)]^{r-3} A_r(Vir).

Since kappa(W_N) = c(H_N - 1), kappa(Vir) = c/2, and
A_r(Vir)/A_{r-1}(Vir) = -6(r-1)/r, one obtains

    A_r(W_N)/A_{r-1}(W_N)
      = -12(H_N - 1)(r-1)/r.

Thus beta_N = 12(H_N - 1). This gives beta_2 = 6, beta_3 = 10,
beta_4 = 13, beta_5 = 77/5, beta_6 = 87/5. The earlier triangular
candidate (N+1)(N+2)/2 and quadratic candidate N^2 - N + 4 are separated
at N = 4 conditional on the kappa-ratio scaling law.

Engine outputs
--------------
This module provides:
  - beta_N(N) : returns the harmonic candidate beta_N.
  - beta_N_candidate_A(N) : triangular comparison candidate.
  - beta_N_candidate_B(N) : quadratic comparison candidate.
  - rho_star_WN(N, c) : |c| / beta_N for the chosen candidate.
  - tempering_rate_bound(N, c, r) : upper bound on (|S_r|/r!)^{1/r}
    that tends to 0 regardless of beta_N's exact form.
  - stirling_vs_beta_dominance(N, c, r) : certificate that beta_N is
    asymptotically subdominant to Stirling, hence tempering.

Convention
----------
beta_N is the harmonic candidate used for finite-beta convention tests.
It is proved at N=2 and N=3, conditional for N>=4. beta_N_candidate_A and
beta_N_candidate_B are retained as comparison formulas separated by the
N=4 harmonic candidate.

Dependencies
------------
    sympy, math, fractions
    compute/lib/independent_verification.py (decorator, not used here)
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from compute.lib.beta_N_closed_form import beta_N_from_kappa, harmonic_number


WN_TEMPERING_SCOPE = {
    "finite_beta_stirling_implication": "proved",
    "finite_riccati_envelope_all_N": "not_encoded",
    "harmonic_beta_closed_form": "conditional_kappa_ratio_scaling",
    "leading_root_sharpness": "assumed_for_radius_equality",
    "w4_full_miura_A5_bridge": "absent",
}


# ---------------------------------------------------------------------------
# Known data from Vol II tempered-stratum chapter
# ---------------------------------------------------------------------------

# Proved closed forms (Vir and W_3)
BETA_2 = Fraction(6)  # Virasoro: proved in thm:tempered-stratum-contains-virasoro
BETA_3 = Fraction(10)  # W_3: proved in thm:tempered-stratum-contains-w3


def beta_N_candidate_A(N: int) -> Fraction:
    r"""Comparison candidate A: beta_N = (N+1)(N+2)/2.

    Matches N=2 and N=3 by coincidence, but predicts beta_4 = 15 while
    the harmonic candidate gives beta_4 = 13.
    """
    if N < 2:
        raise ValueError(f"beta_N requires N >= 2, got N = {N}")
    return Fraction((N + 1) * (N + 2), 2)


def beta_N_candidate_B(N: int) -> Fraction:
    r"""Comparison candidate B: beta_N = N^2 - N + 4.

    Matches N=2 and N=3 by coincidence, but predicts beta_4 = 16 while
    the harmonic candidate gives beta_4 = 13.
    """
    if N < 2:
        raise ValueError(f"beta_N requires N >= 2, got N = {N}")
    return Fraction(N * N - N + 4)


def beta_N(N: int) -> Fraction:
    r"""Harmonic beta_N candidate: beta_N = 12(H_N - 1).

    Closed form:
        beta_N = 12 * (H_N - 1) = sum_{s=2}^{N} 12/s

    This is a proved value at N=2,3 and a conditional kappa-ratio
    candidate for N>=4.
    """
    return beta_N_from_kappa(N)


def beta_N_is_finite(N: int) -> bool:
    r"""The harmonic candidate beta_N is finite for every N >= 2.

    This checks the candidate arithmetic only. It does not prove that the
    full W_N shadow tower has this Riccati coefficient or any other finite
    envelope.
    """
    if N < 2:
        return False
    b = beta_N(N)
    return b > 0 and isinstance(b, Fraction)


# ---------------------------------------------------------------------------
# Tempering rate bound and Stirling dominance certificate
# ---------------------------------------------------------------------------


def tempering_rate_bound(N: int, c: Fraction, r: int) -> float:
    r"""Conditional upper bound on (|S_r(W_N,c)|/r!)^{1/r}.

    Assuming a finite Riccati-style envelope

        |S_r(W_N,c)| <= C(N,c) * beta_N^{r-4} / (r * |c|^{r-2})

    with Stirling's approximation (r!)^{1/r} ~ r/e, we obtain

        (|S_r|/r!)^{1/r}
          <= (C(N,c))^{1/r} * beta_N^{(r-4)/r} * r^{-1/r}
             * |c|^{-(r-2)/r} * (r!)^{-1/r}
          ~ 1 * beta_N * 1 * 1/|c| * e/r
          = beta_N * e / (r * |c|).

    The ratio converges to 0 regardless of beta_N because (r!)^{-1/r} ~ e/r
    decays faster than any polynomial in r. Hence tempering follows from
    the finite-beta envelope. The envelope itself is not proved by this
    function.

    Returns the upper-bound value at finite r.
    """
    if r < 4:
        raise ValueError("Bound requires r >= 4")
    if c == 0:
        raise ValueError("Bound excludes c = 0")
    beta = float(beta_N(N))
    c_abs = abs(float(c))
    C_envelope = 8.0  # universal envelope constant from Laurent tier structure
    # Direct formula for upper bound
    numerator = C_envelope * (beta ** (r - 4))
    denominator = r * (c_abs ** (r - 2)) * math.factorial(r)
    bound = (numerator / denominator) ** (1.0 / r)
    return bound


def stirling_vs_beta_dominance(N: int, c: Fraction, r_max: int = 20) -> bool:
    r"""Certificate: Stirling (r!)^{-1/r} dominates beta_N^{(r-4)/r} at large r.

    For tempering (limsup_r (|S_r|/r!)^{1/r} = 0) to hold, it suffices that
    the Stirling factor (r!)^{-1/r} ~ e/r decays faster than the beta factor
    beta_N^{(r-4)/r} ~ beta_N stays bounded.

    At r_max, check that
        beta_N * e / r_max  <  1
    which is equivalent to r_max > beta_N * e / |c|. For the harmonic
    values through N <= 6 and c = 100 this is immediate.

    Returns True if Stirling dominates.
    """
    beta = float(beta_N(N))
    c_abs = abs(float(c))
    critical_r = beta * math.e / c_abs
    return r_max > critical_r


# ---------------------------------------------------------------------------
# rho_* closed form for W_N
# ---------------------------------------------------------------------------


def rho_star_WN(N: int, c: Fraction) -> Fraction:
    r"""Conditional ordinary-generating convergence radius.

        rho_*^{W_N}(c) = |c| / beta_N = |c| / [12(H_N - 1)]

    At N=2: rho_*(c) = |c|/6 (Virasoro). At N=3:
    rho_*(c) = |c|/10 (W_3). At N=4: rho_*(c) = |c|/13.

    For N>=4 this is conditional on the harmonic kappa-ratio law, the full
    Riccati bridge, and leading-root sharpness. Without sharpness, a finite
    envelope only certifies the lower bound rho_* >= |c| / beta_N.
    """
    if c == 0:
        raise ValueError("rho_* excludes c = 0")
    return abs(c) / beta_N(N)


# ---------------------------------------------------------------------------
# Banach completion and MC convergence guards
# ---------------------------------------------------------------------------


def class_m_banach_norm(weight_norms: Mapping[int, Fraction],
                        rho: Fraction) -> Fraction:
    r"""Weighted Banach norm for the class-M completion.

        ||a||_rho = sum_{w>=0} ||a_w||_w * rho^w

    The input is a finite dictionary of weight norms.  Infinite series are
    handled in the manuscript by the Cauchy criterion; this helper is the
    finite truncation used by the tests.
    """
    if rho <= 0:
        raise ValueError("rho must be positive")
    total = Fraction(0)
    for weight, norm in weight_norms.items():
        if weight < 0:
            raise ValueError("weights must be non-negative")
        if norm < 0:
            raise ValueError("weight norms must be non-negative")
        total += norm * (rho ** weight)
    return total


def operation_bound_holds(output_norm: Fraction,
                          input_norms: Sequence[Fraction],
                          constant: Fraction) -> bool:
    r"""Check ||m_k(a_1,...,a_k)|| <= C_k prod_i ||a_i||."""
    if constant < 0:
        raise ValueError("operation constant must be non-negative")
    product = Fraction(1)
    for norm in input_norms:
        if norm < 0:
            raise ValueError("input norms must be non-negative")
        product *= norm
    return output_norm <= constant * product


def mc_series_bound(operation_constants: Sequence[Fraction],
                    radius: Fraction = Fraction(1)) -> Fraction:
    r"""Finite truncation of sum_{k>=2} C_k(rho) R^k."""
    if radius < 0:
        raise ValueError("radius must be non-negative")
    total = Fraction(0)
    for index, constant in enumerate(operation_constants, start=2):
        if constant < 0:
            raise ValueError("operation constants must be non-negative")
        total += constant * (radius ** index)
    return total


def geometric_mc_series_bound(first_constant: Fraction,
                              ratio: Fraction,
                              radius: Fraction = Fraction(1)) -> Fraction:
    r"""Closed majorant for C_k <= first_constant * ratio^(k-2)."""
    if first_constant < 0:
        raise ValueError("first_constant must be non-negative")
    if not (0 <= ratio * radius < 1):
        raise ValueError("geometric MC majorant requires ratio * radius < 1")
    return first_constant * (radius ** 2) / (1 - ratio * radius)


def banach_radius_certificate(N: int, c: Fraction) -> Dict[str, object]:
    r"""Radius status for the class-M Banach completion.

    N=2 (Virasoro) and N=3 (W_3) are proved rows.  N>=4 records the
    conditional harmonic beta law plus finite-envelope and leading-root
    sharpness hypotheses.
    """
    status = "proved_low_rank" if N in (2, 3) else "conditional_harmonic_beta"
    return {
        "N": N,
        "beta": beta_N(N),
        "rho_star": rho_star_WN(N, c),
        "status": status,
        "requires_mc_sum": "sum_k C_k(rho) < infinity",
    }


# ---------------------------------------------------------------------------
# Check against proved values
# ---------------------------------------------------------------------------


def sanity_check_known_values() -> Dict[str, bool]:
    r"""Verify beta_N formula against the two proved values.

    Expected:
        beta_N(2) = 6 (from Vir)
        beta_N(3) = 10 (from W_3)

    Returns a dict mapping "N=k" to True if the formula matches.
    """
    results = {}
    results["N=2 candidate_A"] = beta_N_candidate_A(2) == BETA_2
    results["N=3 candidate_A"] = beta_N_candidate_A(3) == BETA_3
    results["N=2 candidate_B"] = beta_N_candidate_B(2) == BETA_2
    results["N=3 candidate_B"] = beta_N_candidate_B(3) == BETA_3
    return results


def discriminating_N4_prediction() -> Dict[str, Fraction]:
    r"""Discriminating values at N=4.

    The harmonic candidate value is 13. Candidate_A predicts 15 and
    Candidate_B predicts 16; the three values are distinct.
    """
    return {
        "candidate_A (triangular, (N+1)(N+2)/2)": beta_N_candidate_A(4),
        "candidate_B (quadratic, N^2 - N + 4)": beta_N_candidate_B(4),
        "canonical (harmonic, 12(H_N-1))": beta_N(4),
    }


# ---------------------------------------------------------------------------
# Tempering theorem: direct certificate
# ---------------------------------------------------------------------------


def tempering_certificate(N: int, c: Fraction, r_values: Optional[List[int]] = None) -> Dict[int, float]:
    r"""Return the conditional tempering rate bound at a sequence of r-values.

    Conditional statement: if the fixed W_N shadow tower obeys a finite
    Riccati envelope outside the severe Kac-zero locus D_Kac^sev(W_N), then

        limsup_r (|S_r|/r!)^{1/r} = 0.

    This function computes the finite-r upper bound
        rate(N, c, r) = (C(N,c) * beta_N^{r-4} / (r * |c|^{r-2} * r!))^{1/r}
    which tends to 0. Returns a dict {r: bound}.

    The values monotonically decrease with r, demonstrating convergence
    to 0. This is the concrete numerical content of tempering.
    """
    if r_values is None:
        r_values = [4, 6, 8, 10, 15, 20, 30, 50, 100]
    return {r: tempering_rate_bound(N, c, r) for r in r_values}


# ---------------------------------------------------------------------------
# Self-test on import
# ---------------------------------------------------------------------------


def _self_test():
    """Basic self-test to catch regressions at import time."""
    # Canonical form matches both proved values
    assert beta_N(2) == Fraction(6), f"beta_N(2) = {beta_N(2)}, expected 6"
    assert beta_N(3) == Fraction(10), f"beta_N(3) = {beta_N(3)}, expected 10"
    # Canonical harmonic form is monotone increasing
    for N in range(2, 10):
        assert beta_N(N) < beta_N(N + 1), \
            f"beta_N monotonicity failed at N={N}"
    # Tempering rate decreases
    bounds = tempering_certificate(N=3, c=Fraction(100))
    r_vals = sorted(bounds.keys())
    for i in range(len(r_vals) - 1):
        assert bounds[r_vals[i]] > bounds[r_vals[i + 1]], \
            f"Tempering rate not monotone-decreasing at r={r_vals[i]}"
    # Stirling dominates at large enough r
    assert stirling_vs_beta_dominance(N=6, c=Fraction(100), r_max=100)


_self_test()
