r"""W₃ E₁ ordered multi-channel shadow coefficients to order 20.

The W₃ algebra has generators T (weight 2) and W (weight 3).
The ordered E₁ shadow obstruction tower decomposes into THREE channels:

  S_r^{TT}  — T-T channel (Virasoro sector, decouples)
  S_r^{TW}  — T-W cross channel on the mixed ray
  S_r^{WW}  — W-W channel (sextic pole from {W_λ W})

CHANNEL 1 (TT):
  Shadow metric Q^{TT}(t) = c² + 12ct + [(180c+872)/(5c+22)]t² (degree 2).
  H^{TT}(t) = t² √Q^{TT}(t). S_r = [t^r]H/r. Complementarity c → 26-c.

CHANNEL 2 (WW):
  Z₂ symmetry (W → -W) forces odd arities to vanish.
  κ_WW = c/3 (Shapovalov norm from W_{(5)}W = c/3).
  Q_{WWWW} = 10240/(c(5c+22)³) (quartic contact from w3_quartic_contact.py).
  Shadow metric on the W-line: Q^{WW}(u) = (c/3)² + [20480/(3(5c+22)³)]u
  where u = t². Generating function G_WW(u) = (c/3)√(1 + γu) with
  γ = 61440/(c²(5c+22)³). Complementarity c → 100-c.

CHANNEL 3 (mixed TᵖW² ray):
  The simplest mixed sector has p T-inputs and 2 W-inputs (n_W even by Z₂).
  Shadow metric on the (1,ε) ray:
  Q_mix(t;ε) = Q^{TT}(t) + ε² · δQ(t) + O(ε⁴)
  where δQ encodes the WW perturbation.

The TT shadow metric Q^{TT} has degree 2 (3 parameters: κ, S_3, S_4).
The WW shadow metric Q^{WW}(u) has degree 1 in u = t² (2 parameters: κ_WW, S_4^{WW}).
The higher degree would come from independent S_6^{WW} data; the assumption
of degree-1 Q^{WW} (verified in w3_shadow_coefficients.py) means S_6^{WW}
is determined. All Sh_{2r}^{WW} follow from the binomial expansion.

Cross-references:
  w3_quartic_contact.py: quartic contact Q_{TTTT}, Q_{TTWW}, Q_{WWWW}
  w3_shadow_coefficients.py: preliminary analysis
  ordered_e1_shadow_catalan.py: Virasoro closed-form Catalan formula
  3d_gravity.tex §shadow-table: Virasoro shadow obstruction tower
  w-algebras-w3.tex: W₃ algebra structure and A∞ operations
"""

from __future__ import annotations

import json
import math
import os
import sys
from fractions import Fraction
from typing import Dict, Tuple

from sympy import (
    S, Symbol, Rational, cancel, expand, factor, series, simplify, sqrt, symbols,
    Poly, together,
)


c = Symbol('c', positive=True)
t = Symbol('t')
u = Symbol('u')  # u = t² for the WW channel


# =====================================================================
# 1. TT CHANNEL (VIRASORO)
# =====================================================================

def virasoro_shadow_metric():
    """Q^{TT}(t) = c² + 12ct + [(180c+872)/(5c+22)]t²."""
    return c**2 + 12*c*t + (180*c + 872) / (5*c + 22) * t**2


def tt_shadow_coefficients(r_max: int = 20) -> Dict[int, object]:
    """Compute S_r^{TT} for r = 2,...,r_max.

    H^{TT}(t) = t² √Q^{TT}(t). S_r = [t^r]H / r = [t^{r-2}]√Q / r.
    """
    Q = virasoro_shadow_metric()
    sqQ = series(sqrt(Q), t, 0, n=r_max - 1)
    results = {}
    for r in range(2, r_max + 1):
        Sh_r = cancel(sqQ.coeff(t, r - 2))
        results[r] = cancel(Sh_r / r)
    return results


# =====================================================================
# 2. WW CHANNEL (PURE W SECTOR)
# =====================================================================

def ww_shadow_metric():
    r"""Q^{WW}(u) with u = t² (only even arities).

    κ_WW = c/3. S_4^{WW} = 10240/(c(5c+22)³).
    G_WW(u) = √Q^{WW}(u) = (c/3)·√(1 + γ·u).
    γ = 2·S_4^{WW} / κ_WW = 61440/(c²(5c+22)³).
    Q^{WW}(u) = (c/3)²·(1 + γ·u) = c²/9 + [20480/(3(5c+22)³)]·u.
    """
    kappa_WW = c / 3
    gamma = Rational(61440) / (c**2 * (5*c + 22)**3)
    Q_WW = kappa_WW**2 * (1 + gamma * u)
    return Q_WW, gamma


def ww_shadow_coefficients(r_max: int = 20) -> Dict[int, object]:
    """Compute S_{2r}^{WW} for r = 1,...,r_max//2.

    H^{WW}(t) = t² √Q^{WW}(t²). For even arity 2r:
    S_{2r}^{WW} = [t^{2r}] H^{WW}(t) / (2r)
                = [u^{r-1}] (c/3)·√(1+γu) / (2r)
                = (c/3)·C(1/2, r-1)·γ^{r-1} / (2r)

    where C(1/2, n) = binom(1/2, n) = (1/2)(1/2-1)...(1/2-n+1)/n!
    """
    kappa_WW = c / 3
    gamma = Rational(61440) / (c**2 * (5*c + 22)**3)

    results = {}
    for r in range(1, r_max // 2 + 1):
        arity = 2 * r
        n = r - 1  # power of γ

        # Binomial coefficient binom(1/2, n)
        bcoeff = Rational(1)
        for j in range(n):
            bcoeff *= (Rational(1, 2) - j)
        if n > 0:
            bcoeff /= Rational(math.factorial(n))

        # [t^{2r}]H = [u^{r-1}] G_WW(u) = κ_WW · bcoeff · γ^n
        Sh_2r = kappa_WW * bcoeff * gamma**n

        # S_{2r} = Sh_{2r} / (2r)
        S_2r = cancel(Sh_2r / (2 * r))
        results[arity] = S_2r

    return results


# =====================================================================
# 3. MIXED (Tᵖ W²) CHANNEL
# =====================================================================

def mixed_shadow_metric():
    r"""Shadow metric perturbation for the mixed (T^p, W^2) sector.

    On the ray η = T + ε·W (with α=1, β=ε), the shadow data are:
      Sh_2(1,ε) = c/2 + (c/3)ε²
      Sh_3(1,ε) = 2 + B₃·ε²           [B₃ from m_3 mixed]
      Sh_4(1,ε) = S_4^{TT} + Q_{TTWW}·ε² + O(ε⁴)

    where S_4^{TT} = 10/(c(5c+22)) and Q_{TTWW} = 1920/(c(5c+22)²).

    Define the generating function on this ray:
      √Q_ray(t;ε) = [c/2 + (c/3)ε²] + [2 + B₃ε²]t + [S_4^{TT} + Q_{TTWW}ε²]t² + ...
    where we use √Q_ray meaning the square root of the shadow metric.

    For degree-2 shadow metric: √Q_ray = a₀ + a₁t + a₂t² exactly, and
      Q_ray = a₀² + 2a₀a₁t + (a₁² + 2a₀a₂)t².
    Each a_i depends on ε:
      a₀ = c/2 + (c/3)ε²
      a₁ = 2 + B₃·ε²
      a₂ = 10/(c(5c+22)) + Q_{TTWW}·ε²

    The mixed shadow coefficients S_r^{mix} are the ε² coefficients of S_r(1,ε):
      S_r^{mix} = ∂²S_r/∂ε² |_{ε=0} / 2 = [ε²] S_r(1,ε)

    From the generating function:
      H_ray(t;ε) = t² √Q_ray(t;ε)
    Expanding in ε:
      √Q_ray(t;ε) = √Q^{TT}(t) + ε² · δ(√Q)/δε² + O(ε⁴)

    The perturbation δ(√Q) comes from the variation of (a₀, a₁, a₂):
      δa₀ = c/3,  δa₁ = B₃,  δa₂ = Q_{TTWW}
    and
      δ(√Q) = [√Q is a₀ + a₁t + a₂t² when Q is degree 2]
      δ(√Q) = δa₀ + δa₁·t + δa₂·t²  (direct perturbation of the sqrt series)

    WAIT: this is only correct if Q is EXACTLY degree 2 on the perturbed ray.
    For the full W₃ shadow metric on a general ray, Q could have higher degree
    (from the interaction of TT and WW poles). But to leading order in ε²,
    the perturbation of Q is:
      δQ(t) = 2(a₀·δa₀) + 2(a₀·δa₁ + a₁·δa₀)t + 2(a₁·δa₁ + a₀·δa₂ + a₂·δa₀)t²
    which IS degree 2 (matching the degree of Q^{TT}).

    So: δ(√Q)(t) = δQ(t) / (2√Q^{TT}(t)).

    The mixed shadow coefficients are then:
      Sh_r^{mix} = [t^{r-2}] δ(√Q)(t)
      S_r^{mix} = Sh_r^{mix} / r

    DETERMINING B₃:
    B₃ = ∂²[Sh_3(1,ε)]/∂ε²|₀ / 2 = the cubic shadow in the TW² sector.
    This comes from m_3 with one T-input and two W-inputs.
    By conformal weight: m_3(T,W,W) output has weight ≤ h_T+h_W+h_W - |m_3| = 2+3+3-1 = 7.
    For the SCALAR output: weight = 0, degree = 7. But the maximum spectral degree
    from the OPE poles of {T_λ W} (order 2) composed with {W_μ W} (order 6) is bounded.
    m_3 has 2 spectral params. The scalar term requires the spectral polynomial
    to have total degree ≥ 7. From the λ-bracket compositions, the maximum degree
    from {T_{l1} {W_{l2} W}} is deg({T_λ·}) + deg({W_λ W}) = 3 + 5 = 8, but the
    scalar extraction at weight 0 requires specific cancellations.

    From the COMPUTATION in w3_quartic_contact.py (m3_compute function):
    m_3(T,W,W; l1, l2) = -[{m_2(T,W;l1)}_{l1+l2} W] + [{T_{l1}} m_2(W,W;l2)]

    Let me compute this symbolically.
    """
    from sympy import symbols as sym_symbols

    l1, l2 = sym_symbols('l1 l2')
    T_s = Symbol('T'); W_s = Symbol('W')
    dT_s = Symbol('dT'); dW_s = Symbol('dW')
    d2T_s = Symbol('d2T'); d2W_s = Symbol('d2W')
    d3T_s = Symbol('d3T')
    L_s = Symbol('Lambda'); dL_s = Symbol('dLambda')

    beta_sq = Rational(16) / (22 + 5*c)

    # m_2(T,W;l1) = {T_{l1} W} = ∂W + 3W·l1
    m2_TW = dW_s + 3*W_s*l1

    # m_2(W,W;l2) = {W_{l2} W}
    m2_WW = (Rational(1, 360)*c*l2**5
             + Rational(1, 3)*T_s*l2**3
             + Rational(1, 2)*dT_s*l2**2
             + (beta_sq*L_s + Rational(3, 10)*d2T_s)*l2
             + beta_sq/2*dL_s + Rational(1, 15)*d3T_s)

    # Term 1: {m_2(T,W;l1)}_{l1+l2} W}
    # m_2(T,W;l1) = ∂W + 3W·l1
    # Left sesquilinearity: {∂W_{L} W} = -L·{W_L W} (∂W → -L factor)
    #                       {W_{L} W}·3l1 = 3l1·{W_L W}
    # where L = l1+l2.
    L = l1 + l2
    # {W_L W}:
    WLW = (Rational(1, 360)*c*L**5
           + Rational(1, 3)*T_s*L**3
           + Rational(1, 2)*dT_s*L**2
           + (beta_sq*L_s + Rational(3, 10)*d2T_s)*L
           + beta_sq/2*dL_s + Rational(1, 15)*d3T_s)

    # {∂W_L W} = -L · {W_L W} (by left sesquilinearity: ∂a → -λ)
    term1_dW = -L * WLW  # from the ∂W part of m_2(T,W)
    term1_W = 3*l1 * WLW  # from the 3W·l1 part
    term1 = expand(term1_dW + term1_W)

    # Term 2: {T_{l1} m_2(W,W;l2)}
    # m_2(W,W;l2) = scalar + T terms + Λ terms + W terms... wait, no W in m_2(W,W).
    # Actually m_2(W,W) has outputs: scalar (c/360 l2^5), T, ∂T, ∂²T, ∂³T, Λ, ∂Λ.
    # {T_{l1} scalar} = 0
    # {T_{l1} T} = ∂T + 2T·l1 + (c/12)l1³
    # {T_{l1} ∂T} = (l1+∂)(∂T + 2T·l1 + (c/12)l1³)
    #             = l1·(∂T + 2T·l1 + (c/12)l1³) + (∂²T + 2∂T·l1)
    #             = l1·∂T + 2T·l1² + (c/12)l1⁴ + ∂²T + 2∂T·l1
    #             = ∂²T + 3l1·∂T + 2l1²·T + (c/12)l1⁴
    # {T_{l1} ∂²T} = (l1+∂)²(∂T + 2T·l1 + (c/12)l1³)
    # {T_{l1} ∂³T} = (l1+∂)³(∂T + 2T·l1 + (c/12)l1³)
    # {T_{l1} Λ} = ∂Λ + 4Λ·l1
    # {T_{l1} ∂Λ} = (l1+∂)(∂Λ + 4Λ·l1) = l1·∂Λ + 4Λ·l1² + d2Λ + 4∂Λ·l1
    #              ... but we don't track d2Λ. For the SCALAR extraction, we
    #              only need the scalar part, which comes from {T_{l1} T-field} → scalar.
    # The scalar of {T_{l1} T} = (c/12)l1³.
    # The scalar of {T_{l1} ∂^n T} = (c/12) · (l1+∂)^n applied to l1³ ... hmm,
    # but (l1+∂)^n acts on the OUTPUT fields. Let me think again.

    # Right sesquilinearity: {a_λ ∂^n b} = (λ+∂)^n {a_λ b}
    # So {T_{l1} ∂^n T} = (l1+∂)^n {T_{l1} T} = (l1+∂)^n (∂T + 2T·l1 + (c/12)l1³)

    # For the SCALAR part: we need the coefficient of the identity field (1).
    # In {T_{l1} T} = ∂T + 2T·l1 + (c/12)l1³, the scalar is (c/12)l1³.
    # (l1+∂)^n of the scalar part: (l1+∂)^n ((c/12)l1³) = (c/12)·(l1+∂)^n(l1³).
    # But (l1+∂) acts on the OUTPUT. The scalar l1³ is a spectral polynomial,
    # NOT a field. So ∂ of a scalar = 0. Therefore (l1+∂)^n(scalar) = l1^n · scalar.
    # So: scalar of {T_{l1} ∂^n T} = l1^n · (c/12)l1³ = (c/12)l1^{n+3}.

    # Collecting term 2 scalar:
    # m_2(W,W;l2) = (c/360)l2^5 · (scalar)
    #             + (l2³/3)·T + (l2²/2)·∂T + (β²L_s + 3/10·∂²T)·l2 + (β²/2·∂Λ + 1/15·∂³T)
    # The T, ∂T, ∂²T, ∂³T coefficients:
    coeff_T = l2**3 / 3
    coeff_dT = l2**2 / 2
    coeff_d2T = Rational(3, 10) * l2
    coeff_d3T = Rational(1, 15)

    # {T_{l1} m_2(W,W;l2)}|_scalar = Σ_n coeff_{∂^n T} · scalar of {T_{l1} ∂^n T}
    # = Σ_n coeff_{∂^n T} · (c/12) · l1^{n+3}

    term2_scalar = (c / 12) * (
        coeff_T * l1**3     # n=0: T
        + coeff_dT * l1**4  # n=1: ∂T
        + coeff_d2T * l1**5 # n=2: ∂²T
        + coeff_d3T * l1**6 # n=3: ∂³T
    )
    term2_scalar = expand(term2_scalar)

    # The Λ and ∂Λ terms: {T_{l1} Λ} has NO scalar part (Λ has weight 4,
    # and {T_λ Λ} = ∂Λ + 4Λλ which has no scalar). Similarly for ∂Λ.
    # So the Λ terms contribute zero to the scalar. ✓

    # Term 1 scalar: extract scalar from term1
    # term1 = expand(term1_dW + term1_W) where both involve WLW.
    # The scalar part of WLW = (c/360)L^5 (the only scalar term in {W_L W}).
    # term1_dW = -L · WLW: scalar part = -L · (c/360)L^5 = -(c/360)L^6
    # term1_W = 3l1 · WLW: scalar part = 3l1 · (c/360)L^5 = (c/120)l1·L^5
    term1_scalar = expand(-(c/360)*L**6 + (c/120)*l1*L**5)

    # m_3(T,W,W; l1, l2)|_scalar = -(term1_scalar - term2_scalar)
    m3_TWW_scalar = expand(-(term1_scalar - term2_scalar))

    print(f"  m_3(T,W,W; l1, l2)|_scalar = {m3_TWW_scalar}")
    print(f"  Simplified: {factor(m3_TWW_scalar)}")

    # At symmetric point l1 = l2 = 1:
    m3_TWW_sym = m3_TWW_scalar.subs([(l1, 1), (l2, 1)])
    m3_TWW_sym = cancel(m3_TWW_sym)
    print(f"  At symmetric point (1,1): {m3_TWW_sym}")
    print(f"  = {factor(m3_TWW_sym)}")

    return m3_TWW_scalar


# =====================================================================
# 4. FULL COMPUTATION
# =====================================================================

def compute_all_channels(r_max: int = 20):
    """Compute all three W₃ shadow channels to order r_max."""

    results = {'TT': {}, 'WW': {}, 'mixed_B3': None}

    # ─────────────────────────────────────────────────────────────
    print("═" * 78)
    print(f"  W₃ E₁ ORDERED MULTI-CHANNEL SHADOW COEFFICIENTS TO ORDER {r_max}")
    print("═" * 78)

    # ═════════════════════════════════════════════════════════════
    # CHANNEL 1: TT
    # ═════════════════════════════════════════════════════════════
    print("\n" + "═" * 78)
    print("  CHANNEL 1: S_r^{TT}  [Virasoro, T-sector decouples]")
    print("  Q^{TT}(t) = c² + 12ct + [(180c+872)/(5c+22)]t²")
    print("  H^{TT}(t) = t² √Q^{TT}(t),  S_r = [t^r]H/r")
    print("═" * 78)

    S_TT = tt_shadow_coefficients(r_max)
    results['TT'] = S_TT
    for r in range(2, r_max + 1):
        S_r = S_TT[r]
        # Factor to get a clean form
        try:
            S_r_clean = factor(S_r)
        except Exception:
            S_r_clean = S_r
        print(f"  r={r:2d}:  S_{r}^{{TT}} = {S_r_clean}")

    # ═════════════════════════════════════════════════════════════
    # CHANNEL 2: WW
    # ═════════════════════════════════════════════════════════════
    print("\n" + "═" * 78)
    print("  CHANNEL 2: S_r^{WW}  [W-W channel, even arities only]")
    print("  Q^{WW}(u) = c²/9 + [20480/(3(5c+22)³)]u,  u = t²")
    print("  G^{WW}(u) = (c/3)√(1+γu),  γ = 61440/(c²(5c+22)³)")
    print("  S_{2r}^{WW} = (c/3)·C(1/2,r-1)·γ^{r-1} / (2r)")
    print("═" * 78)

    S_WW = ww_shadow_coefficients(r_max)
    results['WW'] = S_WW
    for arity in sorted(S_WW.keys()):
        S_r = S_WW[arity]
        try:
            S_r_clean = factor(S_r)
        except Exception:
            S_r_clean = S_r
        print(f"  r={arity:2d}:  S_{arity}^{{WW}} = {S_r_clean}")

    # ═════════════════════════════════════════════════════════════
    # CHANNEL 3: Mixed (cubic data)
    # ═════════════════════════════════════════════════════════════
    print("\n" + "═" * 78)
    print("  CHANNEL 3: MIXED TW² (cubic shadow B₃)")
    print("═" * 78)

    m3_TWW_scalar = mixed_shadow_metric()
    results['mixed_B3'] = m3_TWW_scalar

    # ═════════════════════════════════════════════════════════════
    # COMPLEMENTARITY
    # ═════════════════════════════════════════════════════════════
    print("\n" + "═" * 78)
    print("  COMPLEMENTARITY CHECKS")
    print("═" * 78)

    # TT: c → 26-c (Virasoro self-dual at c=13)
    print("\n  TT channel: c → 26-c")
    for r in [4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20]:
        if r not in S_TT:
            continue
        S_r = S_TT[r]
        S_r_dual = S_r.subs(c, 26 - c)
        S_sum = cancel(S_r + S_r_dual)
        # Check if S_sum is c-independent
        p = Poly(together(S_sum) * (5*c + 22)**(r-2) * c**(r-3), c)
        is_const = p.degree() == 0 if p.total_degree() >= 0 else False
        # For Virasoro, the complementarity is:
        # S_r^{norm}(c) + S_r^{norm}(26-c) = const
        # where S_r^{norm} = S_r · c^{r-3} · (5c+22)^{floor((r-2)/2)}
        #
        # Actually the Virasoro complementarity is more subtle.
        # Let me just compute the sum at specific c values.
        s13 = S_r.subs(c, 13)
        s13_dual = S_r_dual.subs(c, 13)
        print(f"    r={r:2d}: S(c) + S(26-c) at c=13: {cancel(s13 + s13_dual)} (should be 2·S(13) by symmetry)")

    # WW: c → 100-c (W₃ self-dual at c=50)
    print("\n  WW channel: c → 100-c")
    print("  The WW channel has complementarity S_{2r}^{WW,norm}(c) + S_{2r}^{WW,norm}(100-c) = const")
    print("  where the normalization removes the common denominator factors.")
    for arity in [2, 4, 6, 8, 10]:
        if arity not in S_WW:
            continue
        S_r = S_WW[arity]
        S_r_dual = S_r.subs(c, 100 - c)
        # At c=50 (self-dual point):
        s50 = S_r.subs(c, 50)
        s50_dual = S_r_dual.subs(c, 50)
        print(f"    r={arity:2d}: S(50) = {s50}, S(100-50) = {s50_dual}, match: {cancel(s50 - s50_dual) == 0}")

    # ═════════════════════════════════════════════════════════════
    # CLOSED-FORM FORMULAS
    # ═════════════════════════════════════════════════════════════
    print("\n" + "═" * 78)
    print("  CLOSED-FORM FORMULAS")
    print("═" * 78)

    print("\n  TT channel (Virasoro Catalan formula):")
    print("    S_2^{TT} = c/2")
    print("    S_3^{TT} = 2")
    print("    S_r^{TT} = (-6)^{r-4} · D / (2r · c^{r-3}) · F_r(D/144)")
    print("    where D = 80/(5c+22), F_r(x) = Σ_{j=0}^{⌊(r-4)/2⌋} (-1)^j Cat_j C(r-4,2j) x^j")

    print("\n  WW channel (generalized Catalan for the sextic pole):")
    print("    S_2^{WW} = c/6")
    print("    S_{2r}^{WW} = (c/3) · C(1/2, r-1) · γ^{r-1} / (2r)")
    print("    where γ = 61440/(c²(5c+22)³)")
    print("    = (-1)^{r-2} · (2r-4)!! · c · 61440^{r-1} / (3 · (2r)!! · (2r) · c^{2r-2} · (5c+22)^{3r-3})")

    # Simplify the WW formula
    print("\n  WW closed-form via generalized Catalan:")
    print("    C(1/2, n) = (-1)^{n+1} · Cat_{n-1} / (2^{2n-1})")
    print("    where Cat_m = C(2m,m)/(m+1) is the m-th Catalan number.")
    print("    So for n ≥ 1:")
    print("    S_{2r}^{WW} = (-1)^r · Cat_{r-2} · c · (61440)^{r-1} / (3 · 2^{2r-1} · (2r) · c^{2r-2} · (5c+22)^{3r-3})")

    # Verify the Catalan connection for WW
    print("\n  Verification of Catalan connection for WW:")
    for r in range(2, 6):
        n = r - 1
        bcoeff = Rational(1)
        for j in range(n):
            bcoeff *= (Rational(1, 2) - j)
        if n > 0:
            bcoeff /= Rational(math.factorial(n))

        # Check against Catalan formula: C(1/2, n) = (-1)^{n+1} Cat_{n-1} / 2^{2n-1}
        if n >= 1:
            cat_nm1 = math.comb(2*(n-1), n-1) // n  # Cat_{n-1}
            catalan_formula = Rational((-1)**(n+1) * cat_nm1, 2**(2*n-1))
            match = bcoeff == catalan_formula
            print(f"    n={n}: C(1/2,{n}) = {bcoeff}, Catalan formula = {catalan_formula}, match = {match}")
        else:
            print(f"    n=0: C(1/2,0) = 1 (base case)")

    return results


# =====================================================================
# 5. NUMERICAL VERIFICATION
# =====================================================================

def numerical_check():
    """Verify the symbolic shadow coefficients at specific c values."""
    print("\n" + "═" * 78)
    print("  NUMERICAL VERIFICATION AT SPECIFIC c VALUES")
    print("═" * 78)

    S_TT = tt_shadow_coefficients(20)
    S_WW = ww_shadow_coefficients(20)

    for cv in [Rational(1), Rational(10), Rational(25), Rational(50), Rational(99)]:
        print(f"\n  c = {cv}:")

        # TT channel at this c
        print("    TT channel:")
        for r in [4, 6, 8, 10, 12, 14, 16, 18, 20]:
            val = S_TT[r].subs(c, cv)
            print(f"      S_{r:2d}^{{TT}} = {float(val):+.15e}")

        # WW channel at this c
        print("    WW channel:")
        for r in [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]:
            if r in S_WW:
                val = S_WW[r].subs(c, cv)
                print(f"      S_{r:2d}^{{WW}} = {float(val):+.15e}")

    # Verify TT matches the Virasoro Catalan formula
    print("\n  Cross-check: TT vs Virasoro Catalan formula")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from ordered_e1_shadow_catalan import shadow_catalan_exact
        cv_frac = Fraction(25)
        catalan_results = shadow_catalan_exact(cv_frac, r_max=20)
        for r in range(4, 21):
            tt_val = float(S_TT[r].subs(c, 25))
            cat_val = float(catalan_results[r])
            match = abs(tt_val - cat_val) < 1e-12 * max(abs(tt_val), abs(cat_val), 1e-300)
            if not match:
                print(f"    r={r}: TT={tt_val:.15e}, Catalan={cat_val:.15e} *** MISMATCH ***")
            else:
                print(f"    r={r}: OK ({tt_val:.10e})")
    except ImportError:
        print("    (ordered_e1_shadow_catalan.py not available for cross-check)")


# =====================================================================
# 6. COMPREHENSIVE COMPLEMENTARITY ANALYSIS
# =====================================================================

def complementarity_analysis():
    """Full complementarity analysis for all three channels."""
    print("\n" + "═" * 78)
    print("  COMPREHENSIVE COMPLEMENTARITY ANALYSIS")
    print("═" * 78)

    S_TT = tt_shadow_coefficients(20)
    S_WW = ww_shadow_coefficients(20)

    # TT: Virasoro complementarity c → 26-c
    print("\n  TT channel: Virasoro complementarity c → 26-c")
    print("  Normalized: S_r^{norm}(c) = S_r(c) · c^{r-3} · (5c+22)^{⌊(r-2)/2⌋}")
    for r in range(4, 21):
        S_r = S_TT[r]
        S_r_dual = S_r.subs(c, 26 - c)
        # Compute the normalized sum
        p = (r - 2) // 2
        norm = c**(r-3) * (5*c + 22)**p
        S_norm = cancel(S_r * norm)
        S_norm_dual = cancel(S_r_dual.subs(c, 26-c) * (26-c)**(r-3) * (5*(26-c)+22)**p)
        # Actually, let me just evaluate the sum at several c
        s_vals = []
        for cv in [1, 5, 10, 13, 20, 25]:
            val = float(S_r.subs(c, cv)) + float(S_r.subs(c, 26-cv))
            s_vals.append(val)
        # Check if all values are the same (constant sum)
        if len(set(f"{v:.10f}" for v in s_vals)) == 1:
            const = s_vals[0]
            print(f"    r={r:2d}: S(c) + S(26-c) = {const:+.10e}  [CONSTANT] ✓")
        else:
            # Not constant — the normalized version might be
            norm_vals = []
            for cv in [1, 5, 10, 13, 20, 25]:
                n1 = float(S_r.subs(c, cv)) * cv**(r-3) * (5*cv+22)**p
                n2 = float(S_r.subs(c, 26-cv)) * (26-cv)**(r-3) * (5*(26-cv)+22)**p
                norm_vals.append(n1 + n2)
            if max(norm_vals) - min(norm_vals) < 1e-6 * max(abs(v) for v in norm_vals):
                print(f"    r={r:2d}: S_norm(c) + S_norm(26-c) = {norm_vals[0]:+.10e}  [NORM CONST] ✓")
            else:
                print(f"    r={r:2d}: NOT constant. Values: {[f'{v:.6e}' for v in norm_vals[:4]]}")

    # WW: c → 100-c (W₃ complementarity)
    print("\n  WW channel: c → 100-c")
    print("  For Q^{WW}(u) = (c/3)²(1+γu), the complementarity structure is:")
    print("  S_{2r}^{WW}(c) involves (c/3)·γ^{r-1} = c/(3)·(61440)^{r-1}/(c²(5c+22)³)^{r-1}")
    for arity in [4, 6, 8, 10, 12, 14, 16, 18, 20]:
        if arity not in S_WW:
            continue
        S_r = S_WW[arity]
        S_r_dual = S_r.subs(c, 100 - c)
        # Check at several c values
        s_vals = []
        for cv in [1, 10, 25, 50, 75, 90, 99]:
            val = float(S_r.subs(c, cv)) + float(S_r_dual.subs(c, cv))
            s_vals.append(val)
        if max(s_vals) - min(s_vals) < 1e-8 * max(abs(v) for v in s_vals if v != 0):
            const = s_vals[0]
            print(f"    r={arity:2d}: S(c) + S(100-c) = {const:+.10e}  [CONSTANT] ✓")
        else:
            # Try normalized
            r_half = arity // 2
            p = 3 * (r_half - 1)
            norm_vals = []
            for cv in [1, 10, 25, 50, 75, 90, 99]:
                n1 = float(S_r.subs(c, cv)) * cv**(2*r_half-2) * (5*cv+22)**p
                n2 = float(S_r.subs(c, 100-cv)) * (100-cv)**(2*r_half-2) * (5*(100-cv)+22)**p
                norm_vals.append(n1 + n2)
            spread = max(norm_vals) - min(norm_vals)
            scale = max(abs(v) for v in norm_vals) if any(v != 0 for v in norm_vals) else 1
            if spread < 1e-6 * scale:
                print(f"    r={arity:2d}: S_norm(c) + S_norm(100-c) = {norm_vals[0]:+.10e}  [NORM CONST] ✓")
            else:
                print(f"    r={arity:2d}: Values: {[f'{v:.4e}' for v in norm_vals[:4]]}")


# =====================================================================
# 7. SAVE RESULTS
# =====================================================================

def save_results(r_max: int = 20):
    """Save all shadow coefficients to a JSON file."""
    S_TT = tt_shadow_coefficients(r_max)
    S_WW = ww_shadow_coefficients(r_max)

    data = {
        'description': f'W₃ E₁ ordered multi-channel shadow coefficients to order {r_max}',
        'channels': {
            'TT': {
                'description': 'Virasoro (T-sector decouples). Complementarity c → 26-c.',
                'shadow_metric': 'Q^{TT}(t) = c² + 12ct + [(180c+872)/(5c+22)]t²',
                'coefficients': {},
            },
            'WW': {
                'description': 'Pure W sector, even arities only. Complementarity c → 100-c.',
                'shadow_metric': 'Q^{WW}(u) = c²/9 + [20480/(3(5c+22)³)]u, u=t²',
                'gamma': '61440/(c²(5c+22)³)',
                'formula': 'S_{2r}^{WW} = (c/3)·C(1/2,r-1)·γ^{r-1}/(2r)',
                'coefficients': {},
            },
        },
        'complementarity': {
            'TT': 'c → 26-c (Virasoro, c*=13)',
            'WW': 'c → 100-c (W₃, c*=50, c_crit=100)',
        },
    }

    for r in range(2, r_max + 1):
        S_r = S_TT[r]
        data['channels']['TT']['coefficients'][str(r)] = str(S_r)

    for arity in sorted(S_WW.keys()):
        S_r = S_WW[arity]
        data['channels']['WW']['coefficients'][str(arity)] = str(S_r)

    # Numerical values at c = 25 and c = 50
    for cv_name, cv_val in [('c=25', 25), ('c=50', 50)]:
        data[f'numerical_{cv_name}'] = {
            'TT': {str(r): float(S_TT[r].subs(c, cv_val)) for r in range(2, r_max+1)},
            'WW': {str(a): float(S_WW[a].subs(c, cv_val)) for a in sorted(S_WW.keys())},
        }

    outpath = os.path.join(os.path.dirname(__file__), 'w3_multichannel_shadow_results.json')
    with open(outpath, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"\n  Results saved to {outpath}")

    return data


# =====================================================================
# MAIN
# =====================================================================

if __name__ == '__main__':
    print("╔" + "═" * 76 + "╗")
    print("║  W₃ E₁ ORDERED MULTI-CHANNEL SHADOW COMPUTATION                         ║")
    print("╚" + "═" * 76 + "╝")

    results = compute_all_channels(r_max=20)
    complementarity_analysis()
    numerical_check()
    data = save_results(r_max=20)
