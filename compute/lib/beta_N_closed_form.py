r"""beta_N harmonic candidate under the kappa-ratio scaling ansatz.

The Vol II tempered-stratum chapters
(tempered_stratum_characterization_platonic.tex and
wn_tempered_closure_platonic.tex) prove the finite-beta Stirling
implication for W_N after a finite Riccati envelope is supplied. The
candidate beta_N is defined via the leading-Laurent ratio identity

    A_{r}(W_N) / A_{r-1}(W_N)  =  -beta_N * (r-1) / r    (r >= 5).

The proved data were:
    beta_2 = 6    (Virasoro)
    beta_3 = 10   (W_3)
and two extrapolation candidates were stated:
    Candidate A:   beta_N = (N+1)(N+2)/2           -> beta_4 = 15
    Candidate B:   beta_N = N^2 - N + 4            -> beta_4 = 16

CONDITIONAL SCALING ANSATZ (this module).
----------------------------------------
The leading-Laurent coefficient A_r(W_N) is assumed to satisfy the
kappa-ratio scaling

    A_r(W_N) = [kappa(W_N) / kappa(Vir)]^{r-3} . A_r(Vir) (T-line),

which is the UNIQUE scaling consistent with Vol II's A_4^{W_3} = 10/3 datum
(line 626 of tempered_stratum_characterization_platonic.tex):

    A_4(W_3) = (kappa(W_3)/kappa(Vir)) . A_4(Vir) = (5/3) . 2 = 10/3.

Under this kappa-ratio scaling ansatz, the leading ratio becomes

    A_r(W_N) / A_{r-1}(W_N)  =  -6 . (r-1) / r . [kappa(W_N)/kappa(Vir)]
                             =  -6 . (r-1) / r . [c(H_N-1) / (c/2)]
                             =  -12 (H_N - 1) . (r-1) / r.

Therefore the harmonic candidate is

    beta_N  =  12 * (H_N - 1)                                          (C)

where H_N = sum_{j=1}^{N} 1/j is the N-th harmonic number.

Arithmetic checks
-----------------
    N=2:  beta_2  =  12 * (1/2)     =  6      PROVED DATUM
    N=3:  beta_3  =  12 * (5/6)     =  10     PROVED DATUM
    N=4:  beta_4  =  12 * (13/12)   =  13     CONDITIONAL CANDIDATE
    N=5:  beta_5  =  12 * (77/60)   =  77/5   CONDITIONAL CANDIDATE
    N=6:  beta_6  =  12 * (29/20)   =  87/5   CONDITIONAL CANDIDATE
    N=7:  beta_7  =  12 * (223/140) =  669/35 CONDITIONAL CANDIDATE

beta_N is rational-valued.  Under the harmonic ansatz it is integral at
N=2,3,4 and non-integral for every N >= 5 by the Bertrand-prime
denominator obstruction.  This separates both Candidate A and Candidate B,
which demand integer beta_N for all N.

Relation to prior candidates
---------------------------
Candidate A (beta_N = (N+1)(N+2)/2) predicts integer values, matches N=2,3
by coincidence, predicts beta_4 = 15. CONDITIONALLY EXCLUDED by the
harmonic candidate value 13.

Candidate B (beta_N = N^2-N+4) predicts integer values, matches N=2,3
by coincidence, predicts beta_4 = 16. CONDITIONALLY EXCLUDED by the
harmonic candidate value 13.

Candidate C (this module, harmonic ansatz) gives beta_N = 12(H_N-1),
a rational-valued function of N through the harmonic number H_N. Matches
N=2,3 by the proved low-rank data and predicts beta_4 = 13 if the
scaling ansatz holds.

Evidence and missing bridge
---------------------------
The kappa-ratio scaling is checked via:
    (a) direct master equation evaluation for S_5(W_3) leading-c coefficient:
        A_5(W_3) = -80/3 = (5/3)^2 . (-48/5) = (kappa_ratio)^2 . A_5(Vir).
        Consistent with A_r(W_N) = (kappa_ratio)^{r-3} . A_r(Vir).
    (b) Vol I thm:universal-asymptotic-factor: C_{A,ell} on the T-line is
        6 uniformly across class M; the scaling A_r = (kappa_ratio)^{r-3}
        reflects the absorption of the W-generator channel into a kappa
        rescaling at the full-shadow level.
    (c) Consistency with Vol II A_4^{W_3} = 10/3 EXACTLY equals (5/3) * 2.

These checks do not compute A_5(W4) from the full W4 Miura/OPE data.
That bridge is recorded as absent in compute/lib/w4_beta_direct.py.

Files
-----
    compute/lib/beta_N_closed_form.py  (this file)
    compute/tests/test_beta_N_closed_form.py
    chapters/theory/beta_N_closed_form_all_platonic.tex (Vol II inscription)

Independent verification decorator
----------------------------------
See @independent_verification use in test_beta_N_closed_form.py.

Module API
----------
    harmonic_number(N)      -- H_N = sum_{j=1}^{N} 1/j
    beta_N_from_kappa(N)    -- 12 * (H_N - 1)  (conditional candidate)
    beta_N_candidate_A(N)   -- (N+1)(N+2)/2    (triangular comparison)
    beta_N_candidate_B(N)   -- N^2 - N + 4     (quadratic comparison)
    A_r_W_N(N, r)           -- leading Laurent coefficient at the T-line level
    verify_scaling_law(N, r)        -- check A_r = kappa_ratio^{r-3} * A_r^Vir
    verify_ratio_identity(N, r)     -- check A_r/A_{r-1} = -beta_N (r-1)/r
"""

from __future__ import annotations

from fractions import Fraction
from typing import List


# -----------------------------------------------------------------------------
# Harmonic number and kappa
# -----------------------------------------------------------------------------


def harmonic_number(N: int) -> Fraction:
    """H_N = sum_{j=1}^{N} 1/j.

    Examples:
        H_1 = 1
        H_2 = 3/2
        H_3 = 11/6
        H_4 = 25/12
        H_5 = 137/60

    Raises ValueError for N < 1.
    """
    if N < 1:
        raise ValueError(f"H_N requires N >= 1, got N = {N}")
    total = Fraction(0)
    for j in range(1, N + 1):
        total += Fraction(1, j)
    return total


def kappa_WN_ratio(N: int) -> Fraction:
    """kappa(W_N) / kappa(Vir) = 2 * (H_N - 1).

    kappa(W_N) = c (H_N - 1).
    kappa(Vir) = c / 2.
    Ratio = 2(H_N - 1).

    Examples:
        N=2: 2 * (1/2) = 1  (W_2 = Virasoro)
        N=3: 2 * (5/6) = 5/3
        N=4: 2 * (13/12) = 13/6
        N=5: 2 * (77/60) = 77/30

    Raises ValueError for N < 2.
    """
    if N < 2:
        raise ValueError(f"kappa_WN_ratio requires N >= 2, got N = {N}")
    return 2 * (harmonic_number(N) - Fraction(1))


# -----------------------------------------------------------------------------
# beta_N harmonic candidate (Candidate C)
# -----------------------------------------------------------------------------


def beta_N_from_kappa(N: int) -> Fraction:
    r"""Harmonic candidate beta_N = 12 * (H_N - 1).

    Conditional on the kappa-ratio scaling of the leading-Laurent
    coefficient A_r(W_N) = (kappa_ratio)^{r-3} * A_r(Vir), and checked
    against Vol II's A_4^{W_3} = 10/3 datum.

    Values:
        beta_2 = 6
        beta_3 = 10
        beta_4 = 13
        beta_5 = 77/5
        beta_6 = 87/5

    Separates Candidate A (15 at N=4) and Candidate B (16 at N=4) under
    the harmonic ansatz.

    Raises ValueError for N < 2.
    """
    if N < 2:
        raise ValueError(f"beta_N requires N >= 2, got N = {N}")
    return 12 * (harmonic_number(N) - Fraction(1))


def beta_N_candidate_A(N: int) -> Fraction:
    """(N+1)(N+2)/2 = T_{N+1} (triangular comparison formula)."""
    if N < 2:
        raise ValueError(f"beta_N requires N >= 2, got N = {N}")
    return Fraction((N + 1) * (N + 2), 2)


def beta_N_candidate_B(N: int) -> Fraction:
    """N^2 - N + 4 (quadratic comparison formula)."""
    if N < 2:
        raise ValueError(f"beta_N requires N >= 2, got N = {N}")
    return Fraction(N * N - N + 4)


# -----------------------------------------------------------------------------
# Leading Laurent coefficients
# -----------------------------------------------------------------------------


def A_r_Vir(r: int) -> Fraction:
    r"""Virasoro leading-Laurent coefficient: A_r(Vir) = 8 * (-6)^{r-4} / r.

    From Vol I thm:shadow-tower-asymptotic-closed-form.

    Values:
        A_2 = 8 * (-6)^{-2} / 2 = 8 / (36*2) = 1/9
            (formal; only r >= 4 is the intended domain of the closed form)
        A_4 = 8 / 4 = 2
        A_5 = 8 * (-6) / 5 = -48/5
        A_6 = 8 * 36 / 6 = 48
        A_7 = 8 * (-216) / 7 = -1728/7
    """
    if r < 2:
        raise ValueError(f"A_r requires r >= 2, got r = {r}")
    return Fraction(8) * Fraction(-6) ** (r - 4) / Fraction(r)


def A_r_WN(N: int, r: int) -> Fraction:
    r"""Conditional leading-Laurent coefficient A_r(W_N).

    Under the kappa-ratio scaling ansatz:

        A_r(W_N) = [kappa(W_N) / kappa(Vir)]^{r-3} * A_r(Vir)
                  = [2 * (H_N - 1)]^{r-3} * 8 * (-6)^{r-4} / r

    Specializations:
        N=2: A_r(W_2) = A_r(Vir)  [scaling factor (1)^{r-3} = 1]
        N=3: A_4 = (5/3) * 2 = 10/3 (matches Vol II prop:w3-shadow-leading-asymptotic(3))
        N=4: A_4 = (13/6) * 2 = 13/3
        N=4: A_5 = (13/6)^2 * (-48/5) = -676/15
    """
    if N < 2:
        raise ValueError(f"A_r_WN requires N >= 2, got N = {N}")
    if r < 2:
        raise ValueError(f"A_r requires r >= 2, got r = {r}")
    kr = kappa_WN_ratio(N)
    return kr ** (r - 3) * A_r_Vir(r)


# -----------------------------------------------------------------------------
# Verification routines
# -----------------------------------------------------------------------------


def verify_ratio_identity(N: int, r: int) -> bool:
    r"""Verify the candidate identity A_r(W_N) / A_{r-1}(W_N).

    Returns True iff the harmonic-ansatz coefficients satisfy
    A_r(W_N) / A_{r-1}(W_N) = -beta_N * (r-1) / r. Raises ValueError
    for r < 4.
    """
    if r < 4:
        raise ValueError(f"ratio identity requires r >= 4, got r = {r}")
    ar = A_r_WN(N, r)
    arm1 = A_r_WN(N, r - 1)
    ratio = ar / arm1
    expected = -beta_N_from_kappa(N) * Fraction(r - 1, r)
    return ratio == expected


def verify_known_data() -> List[str]:
    r"""Check that beta_N_from_kappa matches the proved data:
        beta_2 = 6, beta_3 = 10.

    Returns a list of failure strings; empty list means all checks passed.
    """
    failures = []
    b2 = beta_N_from_kappa(2)
    if b2 != Fraction(6):
        failures.append(f"beta_2 = {b2}, expected 6")
    b3 = beta_N_from_kappa(3)
    if b3 != Fraction(10):
        failures.append(f"beta_3 = {b3}, expected 10")
    return failures


def verify_a4_w3_matches_vol2() -> bool:
    r"""Vol II prop:w3-shadow-leading-asymptotic (3) states A_4^{W_3} = 10/3.

    Verify this against the kappa-ratio scaling:
        A_4(W_3) = (kappa(W_3)/kappa(Vir))^{4-3} * A_4(Vir)
                 = (5/3) * 2 = 10/3.
    """
    return A_r_WN(3, 4) == Fraction(10, 3)


def beta_4_predictions() -> dict:
    """Discriminating predictions at N=4.

    Returns a dict of candidate -> predicted beta_4.
    """
    return {
        "candidate_A (triangular)": beta_N_candidate_A(4),
        "candidate_B (quadratic)": beta_N_candidate_B(4),
        "candidate_C (harmonic kappa-ratio ansatz)": beta_N_from_kappa(4),
    }


def beta_table(N_max: int = 10) -> List[tuple]:
    """Produce a table of beta_N values for N = 2, ..., N_max.

    Returns a list of (N, beta_A, beta_B, beta_C) tuples.
    """
    rows = []
    for N in range(2, N_max + 1):
        rows.append((
            N,
            beta_N_candidate_A(N),
            beta_N_candidate_B(N),
            beta_N_from_kappa(N),
        ))
    return rows


# -----------------------------------------------------------------------------
# Main (demonstration)
# -----------------------------------------------------------------------------


def main():
    print("=" * 70)
    print("beta_N HARMONIC CANDIDATE: KAPPA-RATIO ANSATZ")
    print("=" * 70)
    print()
    print("Candidate: beta_N = 12 * (H_N - 1)")
    print("    where H_N = sum_{j=1}^{N} 1/j")
    print()

    # Table
    print("  N |  Cand A  |  Cand B  |  Harmonic ansatz (C)")
    print("----+----------+----------+----------------------")
    for N, bA, bB, bC in beta_table(10):
        print(f" {N:2d} | {str(bA):8s} | {str(bB):8s} | {str(bC)}")
    print()

    # Known data verification
    failures = verify_known_data()
    if failures:
        print("!! KNOWN-DATA FAILURES:")
        for f in failures:
            print(f"   {f}")
    else:
        print("Known data verified: beta_2 = 6, beta_3 = 10.")
    print()

    # A_4(W_3) consistency
    a4_w3_match = verify_a4_w3_matches_vol2()
    print(f"Vol II A_4^{{W_3}} = 10/3 consistency: {a4_w3_match}")
    print()

    # Ratio identity verification across r, N
    print("Ratio identity verification A_r/A_{r-1} = -beta_N*(r-1)/r:")
    print("  N | r | A_r       | A_{r-1}  | ratio     | expected | match")
    for N in range(2, 6):
        for r in range(5, 8):
            ar = A_r_WN(N, r)
            arm1 = A_r_WN(N, r - 1)
            ratio = ar / arm1
            expected = -beta_N_from_kappa(N) * Fraction(r - 1, r)
            match = ratio == expected
            print(f"  {N} | {r} | {str(ar):9s} | {str(arm1):9s} | "
                  f"{str(ratio):9s} | {str(expected):9s} | {match}")
    print()

    # Discrimination at N=4
    print("Discrimination at N=4:")
    preds = beta_4_predictions()
    for cand, val in preds.items():
        print(f"  {cand}: beta_4 = {val}")
    print()
    print("Harmonic ansatz gives beta_4 = 13.")
    print("Candidate A prediction (15) conditionally excluded.")
    print("Candidate B prediction (16) conditionally excluded.")


if __name__ == "__main__":
    main()
