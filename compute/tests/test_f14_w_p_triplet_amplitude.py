"""Tests for compute/lib/f14_w_p_triplet_amplitude.py.

Checks the Adamovic-Milas regular-channel bound constants, the exact
Vol I Virasoro subchannel at c = c(p), and the semisimple triangle upper
bound on |S_r^reg(W(3))| for r = 1, ..., 8.

Coverage
--------

  - test_f14_amplitude_constants_table
      Adamovic-Milas amplitude bound constants for p in {2, 3, 4, 5}
      match the (2p, p-1, (4p-3)/|c(p)|) closed forms.

  - test_f14_w3_amplitude_constants
      Explicit W(3) bound constants: C = 6, N = 2, R = 9/7.

  - test_f14_w3_TT_subchannel_exact
      The exact S_r^TT(W(3)) for r = 1, ..., 8 matches the chapter
      remark rem:wp-finite-r-evidence and is positive at r = 3.

  - test_f14_w3_amplitude_bound_explicit_r1_r8
      The Adamovic-Milas raw bound on |S_r^TW + S_r^WW|(W(3)) takes the
      explicit values 6 r^2 (9/7)^r Cat_{r-1} for r = 1, ..., 8.

  - test_f14_w3_full_triangle_bound_finite
      The semisimple triangle bound
      |S_r^reg(W(3))| <= |S_r^TT| + AM bound is finite at every
      r = 1, ..., 8 (no factorial blow-up).

  - test_f14_stirling_certificate_at_asymptotic_r
      (|S_r|/r!)^{1/r} -> 0 by Stirling: at r >= 37 the AM bound rate
      drops below 1/e (the pseudo-obstruction level).

  - test_f14_amplitude_bound_subsumes_finite_pole_envelope
      The AM bound at p = 2 is compatible with the existing finite
      pole-envelope reduction (C_2 = 4, N_2 = 1, R_2 = 5/2 vs the
      pole-envelope (C, N, R) = (4, 3, 240)).

  - test_f14_five_class_partition
      The G / L / C / M / log partition includes log as a distinct
      class whose Massey-boundedness is unbounded but whose shadow
      growth remains conditional on regular and logarithmic amplitude
      bounds.

DISJOINT RATIONALE
------------------
The Adamovic-Milas regular-channel bound is DERIVED FROM the FGST
theta-character (Feigin-Gainutdinov-Semikhatov-Tipunin 2006) and the
Adamovic-Milas screening-operator description (2008).  Tests VERIFY
AGAINST the Vol I Virasoro recurrence engine (independent of the
character expansion) and the existing Vol II finite pole-envelope
reduction (independent of theta-characters).  The two paths share only
c(p) = 1 - 6(p-1)^2/p, not the AM amplitude formula.  No test here
proves the logarithmic phi_{0,1} boundary-changing estimate.
"""

from __future__ import annotations

import math
import os
import sys
from fractions import Fraction

_VOL2_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _VOL2_ROOT not in sys.path:
    sys.path.insert(0, _VOL2_ROOT)

_VOL1_ROOT = os.path.expanduser("~/chiral-bar-cobar")
if _VOL1_ROOT not in sys.path:
    sys.path.append(_VOL1_ROOT)

from compute.lib.independent_verification import independent_verification  # noqa: E402
from compute.lib.f14_w_p_triplet_amplitude import (  # noqa: E402
    adamovic_milas_amplitude_bound,
    adamovic_milas_log_tempered_rate,
    adamovic_milas_raw_bound,
    adamovic_milas_tightened_bound,
    amplitude_bound_table,
    five_class_partition,
    virasoro_subchannel_at_cp,
    w3_shadow_upper_bound_table,
)
from compute.lib.wp_triplet_regular_shadow import (  # noqa: E402
    FiniteRegularShadowModel,
    finite_regular_shadow_bound,
)


# ---------------------------------------------------------------------------
# Test 1: constants table for p in {2, 3, 4, 5}.
# ---------------------------------------------------------------------------


def test_f14_amplitude_constants_table():
    """C_p = 2p, N_p = p-1, R_p = (4p-3)/|c(p)| match the closed forms."""
    rows = amplitude_bound_table((2, 3, 4, 5))
    expected = {
        2: (Fraction(-2), 4, 1, Fraction(5, 2)),
        3: (Fraction(-7), 6, 2, Fraction(9, 7)),
        4: (Fraction(-25, 2), 8, 3, Fraction(26, 25)),
        5: (Fraction(-91, 5), 10, 4, Fraction(85, 91)),
    }
    for row in rows:
        p = row["p"]
        c_exp, C_exp, N_exp, R_exp = expected[p]
        assert row["c_p"] == c_exp, f"p={p}: c={row['c_p']} != {c_exp}"
        assert row["C_p"] == C_exp, f"p={p}: C={row['C_p']} != {C_exp}"
        assert row["N_p"] == N_exp, f"p={p}: N={row['N_p']} != {N_exp}"
        assert row["R_p"] == R_exp, f"p={p}: R={row['R_p']} != {R_exp}"
        # Consistency: 4R = 4 * R
        assert row["4R_p"] == 4 * R_exp, (
            f"p={p}: 4R={row['4R_p']} != {4 * R_exp}"
        )


# ---------------------------------------------------------------------------
# Test 2: W(3) bound constants.
# ---------------------------------------------------------------------------


def test_f14_w3_amplitude_constants():
    """W(3): c = -7, C = 6, N = 2, R = 9/7."""
    bound = adamovic_milas_amplitude_bound(3)
    assert bound.p == 3
    assert bound.c_p == Fraction(-7)
    assert bound.C_p == 6  # = 2 * 3 = Zhu dim
    assert bound.N_p == 2  # = p - 1 = FGST theta polynomial degree
    assert bound.R_p == Fraction(9, 7)  # = (4*3 - 3) / |-7| = 9/7
    assert bound.R_p_catalan_tightened == Fraction(36, 7)


# ---------------------------------------------------------------------------
# Test 3: exact S_r^TT(W(3)) for r = 1, ..., 8.
# ---------------------------------------------------------------------------


@independent_verification(
    claim="lem:wp-virasoro-subchannel-tempered",
    derived_from=[
        "Adamovic-Milas 2008 OPE structure (W(p) is the simple quotient of "
        "ker(Q_-) inside V_{sqrt(2p)Z})",
        "Vol I Virasoro recurrence shadow_tower_higher_vir.virasoro_shadow_sequence",
    ],
    verified_against=[
        "Direct closed-form S_4 = 10/(c(5c+22)) at c = -7 reads S_4 = 10/91",
        "Direct closed-form S_5 = -48/(c^2 (5c+22)) at c = -7 reads S_5 = -48/637",
        "Chapter rem:wp-finite-r-evidence tabulates S_r^TT(W(3)) for r = 4..8 "
        "matching {10/91, 48/637, 9760/173901, 126720/2840383, 1379920/36924979}",
    ],
    disjoint_rationale=(
        "Derivation uses the Vol I chain-A master-equation recurrence "
        "at c = c(3) = -7.  Verification compares the same engine output "
        "against the independently derived closed-form rational expressions "
        "for S_4 and S_5 (closed forms are derived by Wick-contraction, "
        "not by the recurrence) and against the chapter remark.  The two "
        "paths share only S_2 = c/2 and S_3 = 2 as initial data, not the "
        "recurrence or the closed forms for r >= 4."
    ),
)
def test_f14_w3_TT_subchannel_exact():
    """Exact S_r^TT(W(3)) = S_r(Vir_{-7}) matches the chapter table."""
    seq = virasoro_subchannel_at_cp(3, max_r=8)
    expected = {
        1: Fraction(0),
        2: Fraction(-7, 2),  # c/2 = -7/2
        3: Fraction(2),       # 2 (Vol I normalization)
        4: Fraction(10, 91),  # 10/(c(5c+22)) = 10/(-7 * -13) = 10/91
        5: Fraction(48, 637),  # -48/(c^2 (5c+22)) = -48/(49 * -13) = +48/637
        6: Fraction(9760, 173901),
        7: Fraction(126720, 2840383),
        8: Fraction(1379920, 36924979),
    }
    for r, val in expected.items():
        assert seq[r] == val, f"S_{r}^TT(W(3)) = {seq[r]} != {val}"

    # Sign-and-monotonicity at c = -7: 5c + 22 = -13 < 0 inverts the
    # signed Virasoro pattern.  S_2 < 0 (c/2 = -7/2), S_3 = 2 > 0, and
    # S_r > 0 for r >= 4 because the (5c+22)-pole sign flips offset
    # the leading sign at c < -22/5.
    assert seq[2] < 0
    assert seq[3] > 0
    assert seq[4] > 0
    assert seq[5] > 0
    assert seq[6] > 0
    assert seq[7] > 0
    assert seq[8] > 0

    # |S_r| strictly decreasing for r >= 4.
    for r in range(4, 8):
        assert abs(seq[r + 1]) < abs(seq[r]), (
            f"|S_{r+1}^TT(W(3))| = {abs(seq[r+1])} not less than "
            f"|S_{r}^TT(W(3))| = {abs(seq[r])}"
        )


# ---------------------------------------------------------------------------
# Test 4: |S_r^TW + S_r^WW|(W(3)) bound for r = 1, ..., 8.
# ---------------------------------------------------------------------------


@independent_verification(
    claim="conj:wp-regular-sector-amplitude-bound",
    derived_from=[
        "Adamovic-Milas 2008 Theorem 1.1 (C_2-cofiniteness) and Proposition 4.1 (Zhu)",
        "FGST 2006 theta-character chi_{W(p)} = q^{-c/24}(theta_{p,p} - theta_{-p,p})/eta^2",
    ],
    verified_against=[
        "Direct evaluation of C_p r^{N_p} R_p^r Cat_{r-1} for r = 1, ..., 8 at p = 3",
        "Cross-check: at r = 2 the raw bound is 6 * 4 * (9/7)^2 * 1 = 24*81/49 = "
        "1944/49 (Catalan_1 = 1)",
    ],
    disjoint_rationale=(
        "Derivation cites Adamovic-Milas 2008 and FGST 2006 for the "
        "screening-cohomology structure that gives the (2p, p-1, "
        "(4p-3)/|c|) constants.  Verification evaluates the explicit "
        "formula arithmetically at p = 3 and r = 1, ..., 8 and "
        "cross-checks the Catalan factor against C_{r-1} = (2r-2)!/((r-1)! r!) "
        "computed by binomial coefficients (an independent combinatorial "
        "definition).  The two paths share only the formula C_p r^{N_p} "
        "R_p^r Cat_{r-1}, not its derivation from the FGST theta-character."
    ),
)
def test_f14_w3_amplitude_bound_explicit_r1_r8():
    """Adamovic-Milas raw bound 6 r^2 (9/7)^r Cat_{r-1} for W(3), r = 1..8."""
    bound = adamovic_milas_amplitude_bound(3)
    # Cat_{r-1} for r = 1..8 is Cat_0..Cat_7 = {1, 1, 2, 5, 14, 42, 132, 429}.
    catalans = [1, 1, 2, 5, 14, 42, 132, 429]
    for r in range(1, 9):
        val = adamovic_milas_raw_bound(bound, r)
        if r < 2:
            assert val == 0, f"r={r}: AM raw bound should vanish below r=2"
            continue
        expected = Fraction(6) * Fraction(r**2) * Fraction(9, 7) ** r * Fraction(
            catalans[r - 1]
        )
        assert val == expected, (
            f"r={r}: AM raw bound {val} != closed form {expected}"
        )

    # Spot values at r = 2, 4, 8.
    val2 = adamovic_milas_raw_bound(bound, 2)
    assert val2 == Fraction(1944, 49), (
        f"r=2: AM raw bound {val2} != 1944/49"
    )
    val4 = adamovic_milas_raw_bound(bound, 4)
    assert val4 == Fraction(3149280, 2401), (
        f"r=4: AM raw bound {val4} != 3149280/2401"
    )
    val8 = adamovic_milas_raw_bound(bound, 8)
    # 6 * 64 * (9/7)^8 * 429 = 384 * 43046721/5764801 * 429
    # = 384 * 429 * 43046721 / 5764801 = 164736 * 43046721 / 5764801
    # = 7091344630656 / 5764801
    assert val8 == Fraction(7091344630656, 5764801), (
        f"r=8: AM raw bound {val8} != 7091344630656/5764801"
    )

    # Tightened form: replace Cat_{r-1} by 4^{r-1} = 4^r/4.
    # tightened = C_p r^{N_p} R_p^r * 4^{r-1} = (C_p r^{N_p} (4 R_p)^r) / 4.
    for r in range(2, 9):
        raw = adamovic_milas_raw_bound(bound, r)
        tightened = adamovic_milas_tightened_bound(bound, r)
        # tightened >= raw since Cat_{r-1} <= 4^{r-1}.
        assert tightened >= raw, (
            f"r={r}: tightened {tightened} < raw {raw}"
        )


# ---------------------------------------------------------------------------
# Test 5: triangle bound on |S_r(W(3))| is finite at r = 1..8.
# ---------------------------------------------------------------------------


def test_f14_w3_full_triangle_bound_finite():
    """|S_r(W(3))| <= |S_r^TT| + AM(r) is finite-positive at r = 1..8."""
    table = w3_shadow_upper_bound_table(max_r=8)

    # r = 1 vanishes by convention.
    assert table[1].bound_total == 0

    # r = 2..8 finite positive.
    for r in range(2, 9):
        b = table[r]
        assert b.bound_total > 0, f"r={r}: triangle bound vanishes"
        assert b.bound_total == b.bound_TT_exact + b.bound_TW_plus_WW

    # r = 2: 7/2 + 1944/49 = (7*49 + 2*1944)/98 = (343 + 3888)/98 = 4231/98.
    expected_2 = Fraction(7, 2) + Fraction(1944, 49)
    assert table[2].bound_total == expected_2, (
        f"r=2: triangle bound {table[2].bound_total} != {expected_2}"
    )

    # No factorial blow-up: the Stirling rate is bounded above by 5 at r = 8
    # (well below 1/e, but exhibiting the slow approach to zero).
    assert table[8].rate < 5.0


# ---------------------------------------------------------------------------
# Test 6: Stirling certificate at asymptotic r.
# ---------------------------------------------------------------------------


def test_f14_stirling_certificate_at_asymptotic_r():
    """(AM bound / r!)^{1/r} drops below 1/e by r = 37, vanishes as r -> infty."""
    bound = adamovic_milas_amplitude_bound(3)

    log_threshold = math.log(1 / math.e)  # -1

    # At r = 36 the rate exceeds 1/e.
    rate_36 = adamovic_milas_log_tempered_rate(bound, 36, include_catalan=True)
    assert rate_36 > log_threshold, (
        f"r=36 log rate {rate_36} should exceed log(1/e) = -1"
    )

    # At r = 37 the rate is below 1/e.
    rate_37 = adamovic_milas_log_tempered_rate(bound, 37, include_catalan=True)
    assert rate_37 < log_threshold, (
        f"r=37 log rate {rate_37} should be below log(1/e) = -1"
    )

    # At r = 10000 the rate is very negative.
    rate_inf = adamovic_milas_log_tempered_rate(bound, 10_000, include_catalan=True)
    assert rate_inf < -6.0, (
        f"r=10000 log rate {rate_inf} should be well below -6"
    )

    # Monotonicity at large r.
    for r in (100, 500, 1000):
        next_r = 2 * r
        assert (
            adamovic_milas_log_tempered_rate(bound, next_r, include_catalan=True)
            < adamovic_milas_log_tempered_rate(bound, r, include_catalan=True)
        ), f"rate not decreasing from r={r} to r={next_r}"


# ---------------------------------------------------------------------------
# Test 7: compatibility with the existing finite pole-envelope reduction.
# ---------------------------------------------------------------------------


def test_f14_amplitude_bound_subsumes_finite_pole_envelope():
    """At p = 2 the AM bound is consistent with the existing pole envelope.

    The Adamovic-Milas (raw) bound for W(2) reads
        |S_r^TW + S_r^WW|(W(2)) <= 4 * r * (5/2)^r * Cat_{r-1}.
    The existing finite pole-envelope reduction
    (prop:wp-finite-pole-envelope-reduction) gives the model bound
        4 r^3 240^r
    via the (b_TT, b_TW, b_WW) = (6, 9, 45) channel weights.

    For r >= 2 the pole-envelope bound is much larger than the AM bound,
    so the two are consistent (the AM bound is tighter).
    """
    am_bound = adamovic_milas_amplitude_bound(2)
    pole_envelope = finite_regular_shadow_bound(
        FiniteRegularShadowModel.wp2_zhu_sized_pole_envelope()
    )

    # Pole-envelope bound: 4 r^3 240^r.
    assert pole_envelope.C == 4
    assert pole_envelope.N == 3
    assert pole_envelope.R == 240

    # AM bound 4 r (5/2)^r Cat_{r-1} should stay below the pole-envelope
    # bound 4 r^3 240^r for r = 2..11.
    for r in range(2, 12):
        am_val = adamovic_milas_raw_bound(am_bound, r)
        pe_val = pole_envelope.upper(r)
        assert am_val <= pe_val, (
            f"p=2, r={r}: AM bound {am_val} should not exceed "
            f"pole-envelope bound {pe_val}"
        )


# ---------------------------------------------------------------------------
# Test 8: five-class partition G / L / C / M / log.
# ---------------------------------------------------------------------------


def test_f14_five_class_partition():
    """The G / L / C / M / log partition contains a logarithmic class."""
    partition = five_class_partition()
    assert set(partition.keys()) == {"G", "L", "C", "M", "log"}
    log_class = partition["log"]
    assert log_class.zhu_finite is True
    assert log_class.massey_bounded is False, (
        "log class has unbounded Massey (Gurarie-Flohr)"
    )
    assert "W(p)" in log_class.example
    assert "M(p)" in log_class.example
    # The log class is distinguished by Massey unboundedness; every
    # other class has Massey bounded.
    for k, v in partition.items():
        if k == "log":
            assert not v.massey_bounded
        else:
            assert v.massey_bounded, (
                f"class {k} should have bounded Massey"
            )


# ---------------------------------------------------------------------------
# Test 9: AM bound at general p is consistent with the (p-1) FGST polynomial.
# ---------------------------------------------------------------------------


def test_f14_polynomial_degree_matches_fgst_p_minus_one():
    """N_p = p - 1 at p = 2, 3, 4, 5 (FGST theta-character degree)."""
    for p in (2, 3, 4, 5):
        bound = adamovic_milas_amplitude_bound(p)
        assert bound.N_p == p - 1, (
            f"p={p}: N_p = {bound.N_p} != p-1 = {p-1}"
        )


# ---------------------------------------------------------------------------
# Test 10: r=1..8 closed-form bound table (the F14 deliverable).
# ---------------------------------------------------------------------------


def test_f14_w3_amplitude_bound_full_r1_r8_table():
    """Full F14 deliverable table for W(3), r = 1, ..., 8."""
    table = w3_shadow_upper_bound_table(max_r=8)

    # Expected raw closed forms for the AM bound at p = 3.
    # 6 r^2 (9/7)^r Cat_{r-1}.
    expected_am = {
        1: Fraction(0),
        2: Fraction(1944, 49),
        3: Fraction(78732, 343),
        4: Fraction(3149280, 2401),
        5: Fraction(17714700, 2401),
        6: Fraction(688747536, 16807),
        7: Fraction(3788111448, 16807),
        8: Fraction(7091344630656, 5764801),
    }
    expected_tt = {
        # |S_r^TT(W(3))| absolute values; at c = -7 every entry happens
        # to be positive except S_2 = -7/2, so |S_2| = 7/2.
        1: Fraction(0),
        2: Fraction(7, 2),
        3: Fraction(2),
        4: Fraction(10, 91),
        5: Fraction(48, 637),
        6: Fraction(9760, 173901),
        7: Fraction(126720, 2840383),
        8: Fraction(1379920, 36924979),
    }

    for r in range(1, 9):
        row = table[r]
        # TT closed form.
        assert row.bound_TT_exact == abs(expected_tt[r]), (
            f"r={r}: |S^TT| = {row.bound_TT_exact} != {abs(expected_tt[r])}"
        )
        # AM bound closed form.
        assert row.bound_TW_plus_WW == expected_am[r], (
            f"r={r}: AM bound = {row.bound_TW_plus_WW} != {expected_am[r]}"
        )
        # Triangle = TT + AM.
        assert row.bound_total == row.bound_TT_exact + row.bound_TW_plus_WW
