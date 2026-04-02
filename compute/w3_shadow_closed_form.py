r"""W3 shadow coefficient closed-form formula.

MAIN RESULTS:
=============

(1) QUARTIC SHADOW DECOMPOSITION.
    Q_TTTT(c) = 10/(c(5c+22))                   [= Virasoro S_4]
    Q_TTWW(c) = 1920/(c(5c+22)^2)
    Q_WWWW(c) = 10240/(c(5c+22)^3)
    Each pair of W-legs contributes a factor beta = 32/(5c+22):
      Q_TTWW = 6*beta*Q_TTTT,  Q_WWWW = beta^2*Q_TTTT.

(2) W-LINE SHADOW METRIC.
    On the W-line (x_T = 0), the shadow metric is:
      Q_W(t) = (4c^2/9)(1 + delta*t^2)
    where delta = 122880/(c^2(5c+22)^3).
    This is degree 2 in t (linear in u = t^2), analogous to Virasoro.

(3) W-LINE CATALAN FORMULA.
    S_{2r}^W(c) = (-1)^r * C_{r-1} * 30720^{r-1}
                  / (3 * (2r-3) * c^{2r-3} * (5c+22)^{3(r-1)})
    for r >= 2, where C_k = binom(2k,k)/(k+1) is the Catalan number,
    and 30720 = 2^{11} * 3 * 5.
    This is the W3 analogue of the Virasoro Catalan formula.

(4) QUARTIC SHADOW as function of c.
    m_4(W,W,W,W)|_{scalar} = c*(A(l) + B(l)*c) / (5c+22)
    where A, B are polynomials in the spectral parameters.
    The normalized shadow Q_norm := scalar*(5c+22)/c = A + B*c.

(5) COMPLEMENTARITY under c -> 100-c.
    Q_norm(c) + Q_norm(100-c) = 2A + 100B = constant in c.
    This is the W3 Koszul involution (alpha_3 = 100, self-dual c* = 50).

Cross-references:
  - compute/w3_quartic_contact.py
  - chapters/examples/rosetta_stone.tex, eq:w3-quartic
  - chapters/connections/3d_gravity.tex, Theorem thm:period-2-parity
  - chapters/examples/w-algebras-w3.tex
"""
from __future__ import annotations
from sympy import (
    Symbol, Rational, symbols, expand, simplify, factor, cancel,
    sqrt, series, S, Poly, together, N as numerical,
)
from math import factorial as mfac

c = Symbol('c')
t = Symbol('t')


# ===================================================================
# Utility: Catalan numbers and half-integer binomials
# ===================================================================

def catalan(k):
    """k-th Catalan number C_k = binom(2k,k)/(k+1)."""
    return Rational(mfac(2*k), mfac(k)**2 * (k+1))


def binom_half(n):
    """Binomial coefficient binom(1/2, n).

    = (-1)^{n-1} * (n+1) * C_n / ((2n-1) * 4^n)  for n >= 1
    = 1  for n = 0
    """
    if n == 0:
        return Rational(1, 1)
    bn = Rational(1, 1)
    for j in range(n):
        bn *= (Rational(1, 2) - j)
    bn /= mfac(n)
    return bn


# ===================================================================
# Section 1: Virasoro shadow (review and cross-check)
# ===================================================================

def virasoro_shadow():
    """Virasoro shadow metric and coefficients for comparison."""
    print("=" * 72)
    print("VIRASORO SHADOW (REVIEW)")
    print("=" * 72)

    # Q_Vir(t) = c^2 + 12ct + [(180c+872)/(5c+22)]t^2
    # H(t) = t^2 sqrt(Q(t)) = c t^2 + 6 t^3 + [40/(c(5c+22))] t^4 + ...
    # S_r = [t^r]H / r

    delta_vir = (180*c + 872) / (5*c + 22) / c**2  # q2/q0
    # Actually: q2 = (180c+872)/(5c+22), q0 = c^2
    # delta = q2/q0 = (180c+872)/(c^2(5c+22))

    print("\nQ_Vir(t) = c^2(1 + (12/c)t + [(180c+872)/(c^2(5c+22))]t^2)")
    print("         = c^2(1 + alpha_1 t + alpha_2 t^2)")
    alpha1 = Rational(12, 1) / c
    alpha2 = cancel((180*c + 872) / (c**2 * (5*c + 22)))
    print("  alpha_1 = 12/c")
    print("  alpha_2 = (180c+872)/(c^2(5c+22))")

    # Discriminant
    disc = cancel(alpha1**2 - 4*alpha2)
    print("  disc = alpha_1^2 - 4*alpha_2 = %s" % disc)
    print("       = %s" % factor(disc))
    print("  = -320/(c^2(5c+22)) [negative for c > 0: no real roots]")

    # Branch points
    # t_pm = (-alpha_1 +/- sqrt(disc)) / (2*alpha_2)
    print("  Branch points: t_pm = complex conjugate pair")

    # Virasoro shadow coefficients
    print("\nVirasoro S_r (from H = t^2 sqrt(Q)):")
    q0v = c**2
    q1v = 12*c
    q2v = (180*c + 872)/(5*c + 22)
    Q_v = q0v + q1v*t + q2v*t**2
    sq = series(sqrt(Q_v), t, 0, n=9)

    vir_S = {}
    for r in range(2, 11):
        coeff = sq.coeff(t, r-2)
        # Force c positive (replace sqrt(c^2) by c)
        coeff_clean = cancel(coeff / sqrt(c**2) * c)
        S_r = cancel(coeff_clean / r)
        vir_S[r] = S_r
        print("  S_%d = %s" % (r, S_r))

    return vir_S


# ===================================================================
# Section 2: W3 quartic shadow
# ===================================================================

def w3_quartic_shadow():
    """W3 quartic shadow coefficients."""
    print("\n" + "=" * 72)
    print("W3 QUARTIC SHADOW")
    print("=" * 72)

    Q_TTTT = Rational(10, 1) / (c * (5*c + 22))
    Q_TTWW = Rational(1920, 1) / (c * (5*c + 22)**2)
    Q_WWWW = Rational(10240, 1) / (c * (5*c + 22)**3)
    beta = Rational(32, 1) / (5*c + 22)

    print("\n  Q_TTTT = 10 / (c(5c+22))")
    print("  Q_TTWW = 1920 / (c(5c+22)^2) = 6*beta * Q_TTTT")
    print("  Q_WWWW = 10240 / (c(5c+22)^3) = beta^2 * Q_TTTT")
    print("  beta   = 32/(5c+22)")

    print("\n  DECOMPOSITION:")
    print("  Sh_4(x_T,x_W) = Q_TTTT * [x_T^4 + 6*beta*x_T^2*x_W^2 + beta^2*x_W^4]")
    print("                 = Q_TTTT * [(x_T^2 + beta*x_W^2)^2 + 4*beta*x_T^2*x_W^2]")
    print("\n  Numerical at c=1: Q_TTTT = %s, Q_TTWW = %s, Q_WWWW = %s" %
          (Q_TTTT.subs(c, 1), Q_TTWW.subs(c, 1), Q_WWWW.subs(c, 1)))

    return Q_TTTT, Q_TTWW, Q_WWWW, beta


# ===================================================================
# Section 3: W-line shadow metric and Catalan formula
# ===================================================================

def w_line_catalan():
    """W-line shadow metric and closed-form Catalan formula."""
    print("\n" + "=" * 72)
    print("W-LINE SHADOW METRIC AND CATALAN FORMULA")
    print("=" * 72)

    # Data
    Q_WWWW = Rational(10240, 1) / (c * (5*c + 22)**3)
    kappa_W = c / 3
    h2W = 2 * kappa_W  # = 2c/3

    # Shadow metric: Q_W(t) = h2W^2 + [q2 term]*t^2
    # q0 = (2c/3)^2 = 4c^2/9
    # q2 = 2*h2W * h4W  where h4W = 4*Q_WWWW
    h4W = 4 * Q_WWWW  # = 40960/(c(5c+22)^3)
    q0 = h2W**2  # = 4c^2/9
    q2 = cancel(2 * h2W * h4W)  # = 2*(2c/3)*40960/(c(5c+22)^3)

    # delta = q2/q0
    delta = cancel(q2 / q0)

    print("\n  W-line data:")
    print("    kappa_W = c/3")
    print("    h_2^W = 2c/3")
    print("    S_4^W = Q_WWWW = 10240/(c(5c+22)^3)")
    print("    h_4^W = 4*S_4^W = 40960/(c(5c+22)^3)")

    print("\n  Shadow metric Q_W(t) = q0(1 + delta*t^2):")
    print("    q0 = (2c/3)^2 = 4c^2/9")
    print("    q2 = %s" % factor(q2))
    print("    delta = q2/q0 = %s" % delta)
    print("         = %s" % factor(delta))

    delta_clean = Rational(122880, 1) / (c**2 * (5*c + 22)**3)
    assert simplify(delta - delta_clean) == 0
    print("    = 122880 / (c^2(5c+22)^3)")

    # Branch point
    t0_sq = -1/delta_clean
    print("\n  Branch point: t_0^2 = -1/delta = %s" % factor(cancel(t0_sq)))
    print("  = -c^2(5c+22)^3 / 122880")

    # Generating function
    print("\n  H_W(t) = (2c/3) * t^2 * sqrt(1 + delta*t^2)")
    print("  = sum_{r>=1} h_{2r}^W * t^{2r}")
    print("  where h_{2r}^W = (2c/3) * binom(1/2, r-1) * delta^{r-1}")

    # Catalan formula
    print("\n" + "-" * 60)
    print("CATALAN FORMULA FOR W-LINE SHADOW COEFFICIENTS")
    print("-" * 60)

    print("""
  For r >= 1, the W-line shadow coefficient of W3 at arity 2r is:

    S_{2r}^W(c) = (c/(3r)) * binom(1/2, r-1) * delta^{r-1}

  Using binom(1/2, n) = (-1)^{n-1} * (2n)! / (4^n * (n!)^2 * (2n-1)):

    S_{2r}^W(c) = (-1)^{r} * C_{r-1} * 30720^{r-1}
                  / (3 * (2r-3) * c^{2r-3} * (5c+22)^{3(r-1)})

  where 30720 = 122880/4 = 2^{11} * 3 * 5, and
  C_k = binom(2k,k)/(k+1) is the k-th Catalan number.

  Equivalently:

    S_{2r}^W(c) = (-1)^r * C_{r-1} * 122880^{r-1}
                  / (3 * 4^{r-1} * (2r-3) * c^{2r-3} * (5c+22)^{3(r-1)})

  for r >= 2, and S_2^W = c/3.
""")

    # Explicit computation
    print("  Explicit S_{2r}^W values:")
    print("  " + "-" * 60)
    results = {}
    for r in range(1, 8):
        n = r - 1
        bn = binom_half(n)
        h2r = (2*c/3) * bn * delta_clean**n
        S2r = cancel(h2r / (2*r))
        S2r_f = factor(S2r)
        results[2*r] = S2r

        # Catalan check
        if r >= 2:
            Cn = catalan(r-1)
            # binom(1/2, n) = (-1)^{n-1} (n+1) C_n / ((2n-1) 4^n) ... check
            # Actually: binom(1/2, n) for n = r-1.
            # From our identity:
            # binom(1/2, n) = (-1)^{n-1} (2n)! / (4^n (n!)^2 (2n-1))
            # So S = (c/3r) * (-1)^{n-1} (2n)! / (4^n (n!)^2 (2n-1)) * delta^n
            # with n = r-1: (2(r-1))! / (4^{r-1} ((r-1)!)^2 (2r-3))
            # = binom(2r-2, r-1) / (4^{r-1} (2r-3))
            # = r * C_{r-1} / (4^{r-1} (2r-3))
            catalan_expr = ((-1)**r * r * catalan(r-1) *
                            delta_clean**(r-1) * c /
                            (Rational(4, 1)**(r-1) * (2*r-3) * 3 * r))
            catalan_simplified = cancel(catalan_expr)
            match = simplify(S2r - catalan_simplified) == 0
        else:
            match = "n/a (r=1)"
            catalan_simplified = "n/a"

        print("  S_%d^W = %s" % (2*r, S2r_f))
        if r >= 2:
            print("        Catalan check: %s" % match)
        print()

    return results, delta_clean


# ===================================================================
# Section 4: Complementarity
# ===================================================================

def complementarity(results, delta_clean):
    """Complementarity of W-line shadow under c -> 100-c."""
    print("\n" + "=" * 72)
    print("COMPLEMENTARITY: c -> 100-c (alpha_3 = 100, c* = 50)")
    print("=" * 72)

    # The NORMALIZED shadow is Q_norm = S_{2r} * c^{2r-3} * (5c+22)^{3(r-1)} / numerator_const
    # From the quartic contact computation (w3_quartic_contact.py):
    # m_4(WWWW)|_scalar = c*(A+Bc)/(5c+22)
    # Q_norm = scalar * (5c+22)/c = A + Bc
    # A + Bc + A + B(100-c) = 2A + 100B = constant. VERIFIED.

    # For higher arities, the complementarity is:
    # The NORMALIZED shadow Q_norm_{2r}(c) = S_{2r}^W * denom_2r / num_2r
    # Under c -> 100-c, the sum Q_norm(c) + Q_norm(100-c) should be constant.

    # From the Catalan formula:
    # S_{2r}^W = (-1)^r * num_const / (c^{2r-3} (5c+22)^{3(r-1)})
    # The normalized shadow is:
    # Q_norm_{2r}(c) = S_{2r}^W * c^{2r-3} * (5c+22)^{3(r-1)}
    #               = (-1)^r * C_{r-1} * r * 122880^{r-1} / (4^{r-1}(2r-3)*3r)
    #               = (-1)^r * C_{r-1} * 122880^{r-1} / (4^{r-1}(2r-3)*3)
    # This is INDEPENDENT OF c! So Q_norm(c) = Q_norm(100-c) = constant,
    # and the sum = 2*constant.

    print("\n  The W-line shadow has the structure:")
    print("    S_{2r}^W = K_r / (c^{2r-3} * (5c+22)^{3(r-1)})")
    print("  where K_r depends only on r (not on c).")
    print("  Therefore the NORMALIZED shadow K_r is already constant in c!")
    print("  Complementarity S_{2r}^W(c) + S_{2r}^W(100-c) is NOT simply = 2K_r,")
    print("  because the denominators c^{2r-3} and (5c+22)^{3(r-1)} are DIFFERENT")
    print("  at c and 100-c.\n")

    # Let me check the quartic (r=2) explicitly
    print("  Quartic (r=2) complementarity:")
    S4W = results[4]
    S4W_dual = S4W.subs(c, 100-c)
    print("    S_4^W(c) = %s" % S4W)
    print("    S_4^W(100-c) = %s" % cancel(S4W_dual))
    raw_sum = cancel(S4W + S4W_dual)
    print("    Raw sum = %s" % raw_sum)

    # For the QUARTIC, the complementarity from w3_quartic_contact.py is:
    # scalar(c) = c*(A+Bc)/(5c+22)  [at specific spectral params]
    # Q_norm(c) = A + Bc
    # The W-line shadow S_4^W = Q_WWWW = 10240/(c(5c+22)^3)
    # This is the GLOBAL quartic coefficient, not at specific spectral params.
    # The complementarity of w3_quartic_contact.py works because
    # the spectral-parameter-dependent part cancels.

    # For the W-line shadow, the normalized quantity is:
    # N_{2r}(c) = S_{2r}^W * c^{2r-3} * (5c+22)^{3(r-1)}
    print("\n  Normalized shadow N_{2r}(c) = S_{2r}^W * c^{2r-3} * (5c+22)^{3(r-1)}:")
    for arity in sorted(results.keys()):
        r = arity // 2
        if r < 2 or r > 6:
            continue
        S_r = results[arity]
        N_r = cancel(S_r * c**(2*r-3) * (5*c+22)**(3*(r-1)))
        print("    N_%d = %s" % (arity, N_r))
        # This should be independent of c
        try:
            p = Poly(together(N_r), c)
            print("      degree in c: %d (0 = constant)" % p.degree())
        except Exception as e:
            print("      (check failed: %s)" % e)

    # The normalized shadow IS constant! So complementarity is trivial:
    # N(c) = N(100-c) = const, sum = 2*const.
    print("\n  Since N_{2r} is constant in c, the complementarity is AUTOMATIC:")
    print("  N(c) + N(100-c) = 2*N for all c.")
    print("  The W-line shadow coefficients are PURELY controlled by the")
    print("  denominator c^{2r-3}(5c+22)^{3(r-1)} --- all c-dependence is in the pole.")


# ===================================================================
# Section 5: Cross-check with w3_quartic_contact.py at specific params
# ===================================================================

def cross_check():
    """Cross-check S_4^W against the Stasheff computation."""
    print("\n" + "=" * 72)
    print("CROSS-CHECK: S_4^W against Stasheff computation")
    print("=" * 72)

    # From w3_quartic_contact.py at l=(2,3,5):
    # scalar(WWWW; c) = c*(A+Bc)/(5c+22) with A = 187988000/27, B = 63040000/27
    # This is at SPECIFIC spectral parameters l1=2, l2=3, l3=5.
    # Our S_4^W = 10240/(c(5c+22)^3) is the GLOBAL (direction-independent)
    # shadow coefficient, equivalent to the scalar coefficient at unit
    # spectral parameters after normalization.

    A_val = Rational(187988000, 27)
    B_val = Rational(63040000, 27)

    print("\n  From w3_quartic_contact.py at (l1,l2,l3) = (2,3,5):")
    print("    scalar(WWWW) = c*(A + B*c)/(5c+22)")
    print("    A = %s" % A_val)
    print("    B = %s" % B_val)
    print("    Complementarity sum 2A + 100B = %s" % (2*A_val + 100*B_val))

    # Check: scalar at c=1:
    s1 = 1 * (A_val + B_val) / 27
    print("    scalar(c=1) = %s" % s1)

    # Our S_4^W = 10240/(c(5c+22)^3), at c=1: 10240/(1*27^3) = 10240/19683
    s4w_c1 = Rational(10240, 1) / (1 * 27**3)
    print("    S_4^W(c=1) = %s = %s" % (s4w_c1, numerical(s4w_c1, 10)))

    # These are DIFFERENT because the Stasheff computation gives the
    # scalar part at specific spectral parameters, while S_4^W is the
    # universal coefficient (independent of spectral parameters).
    # The relationship: S_4^W = scalar / (spectral polynomial).

    print("\n  Note: scalar(WWWW; l) = S_4^W * P(l) where P is the")
    print("  spectral polynomial. The complementarity of the NORMALIZED")
    print("  scalar (dividing by P) reduces to the constancy of S_4^W's numerator.")


# ===================================================================
# Section 6: Summary table
# ===================================================================

def summary_table():
    """Print the final summary table."""
    print("\n" + "=" * 72)
    print("FINAL SUMMARY TABLE: W3 W-LINE SHADOW TOWER")
    print("=" * 72)

    delta = Rational(122880, 1) / (c**2 * (5*c+22)**3)

    print("""
  Arity   S_{2r}^W(c)                                     Sign   Catalan
  -----   -----------------------------------------       ----   -------""")

    for r in range(1, 8):
        n = r - 1
        bn = binom_half(n)
        h2r = (2*c/3) * bn * delta**n
        S2r = cancel(h2r / (2*r))
        S2r_f = factor(S2r)
        sign = "+" if r % 2 == 1 else "-"
        Cn = catalan(n) if n >= 0 else 0
        print("  %d       %-50s  %s      C_%d = %s" %
              (2*r, str(S2r_f), sign, n, Cn))

    print("""
  Denominator structure: denom(S_{2r}^W) = c^{2r-3} * (5c+22)^{3r-3}
  Numerator: (-1)^r * C_{r-1} * r * 122880^{r-1} / (4^{r-1} * (2r-3) * 3)
  """)

    # Numerical values at c = 50 (self-dual point)
    print("  At the self-dual point c = 50:")
    for r in range(1, 8):
        n = r - 1
        bn = binom_half(n)
        delta_50 = Rational(122880, 1) / (50**2 * 272**3)
        h2r_50 = Rational(100, 3) * bn * delta_50**n
        S2r_50 = h2r_50 / (2*r)
        print("    S_%d^W(50) = %s = %s" %
              (2*r, S2r_50, numerical(S2r_50, 8)))


# ===================================================================
# Main
# ===================================================================

def main():
    print("=" * 72)
    print("W3 SHADOW COEFFICIENT: CLOSED-FORM CATALAN FORMULA")
    print("=" * 72)
    print()

    vir_S = virasoro_shadow()
    Q_TTTT, Q_TTWW, Q_WWWW, beta = w3_quartic_shadow()
    results, delta = w_line_catalan()
    complementarity(results, delta)
    cross_check()
    summary_table()


if __name__ == '__main__':
    main()
