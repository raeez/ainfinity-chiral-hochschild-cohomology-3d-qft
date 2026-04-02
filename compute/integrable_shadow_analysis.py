#!/usr/bin/env python3
r"""Shadow tower analysis at special integrable values of the central charge.

Computes shadow invariants S_r, convergence radii R(c), Darboux asymptotics,
and m_k operations at the following physically distinguished values:

    c = 1/2  (Ising model, minimal model M(3,4))
    c = 1    (free boson at self-dual radius)
    c = 13   (Koszul self-dual point: Vir_13^! = Vir_13)
    c = 25   (just below critical string c=26)
    c = 26   (critical string: kappa_eff = 0)

The shadow metric is Q_L(t) = c^2 + 12c*t + ((180c+872)/(5c+22))*t^2
and the shadow coefficients S_r = a_{r-2}/r where a_n are the Taylor
coefficients of sqrt(Q_L(t)).

The convergence radius is R(c) = |t_+| where t_+/- are the zeros of Q_L.
"""
from __future__ import annotations

import sys
import os
import math
import cmath
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'lib'))

from symbolic_stasheff import (
    mk_exact_numerical,
    m4_virasoro_symbolic,
    m5_virasoro_symbolic,
)
from shadow_borel_resurgence import (
    VirasoroShadowData,
    shadow_coefficients,
    shadow_coefficients_fraction,
    darboux_coefficients,
    borel_singularities,
    asymptotic_prediction,
)


# =====================================================================
# Section 1: Shadow data at special central charges
# =====================================================================

SPECIAL_CHARGES = {
    'Ising (c=1/2)':   0.5,
    'Free boson (c=1)': 1.0,
    'Self-dual (c=13)': 13.0,
    'Near-critical (c=25)': 25.0,
    'Critical string (c=26)': 26.0,
}


def print_header(title: str):
    print()
    print('=' * 72)
    print(f'  {title}')
    print('=' * 72)


def analyze_shadow_data(label: str, c_val: float, r_max: int = 30):
    """Complete shadow analysis at a given central charge."""
    print_header(f'{label}  (c = {c_val})')

    d = VirasoroShadowData(c_val)

    print(f'\n--- Basic shadow invariants ---')
    print(f'  kappa(Vir_c) = c/2 = {d.kappa}')
    print(f'  kappa(Vir_c^!) = (26-c)/2 = {(26 - c_val)/2}')
    print(f'  kappa_eff = kappa + kappa^! = 26/2 - c = {13 - c_val}')
    print(f'    (genus tower vanishes iff kappa_eff = 0, i.e. c = 26)')
    print(f'  S_2 = kappa = c/2 = {c_val/2}')
    print(f'  S_3 = alpha = 2  (universal, c-independent)')
    print(f'  S_4 = 10 / (c*(5c+22)) = {d.S4:.10f}')
    print(f'  Delta = 40/(5c+22) = {d.Delta:.10f}')

    print(f'\n--- Shadow metric Q_L(t) = q0 + q1*t + q2*t^2 ---')
    print(f'  q0 = c^2 = {d.q0}')
    print(f'  q1 = 12c = {d.q1}')
    print(f'  q2 = (180c+872)/(5c+22) = {d.q2:.10f}')
    disc = d.q1**2 - 4 * d.q0 * d.q2
    print(f'  discriminant(Q_L) = q1^2 - 4*q0*q2 = {disc:.10f}')

    print(f'\n--- Branch points of sqrt(Q_L) ---')
    print(f'  t_+ = {d.t_plus}')
    print(f'  t_- = {d.t_minus}')
    print(f'  |t_+| = |t_-| = R(c) = {d.R:.10f}')
    print(f'  arg(t_+) = {d.theta:.10f}  ({d.theta * 180 / math.pi:.4f} degrees)')

    print(f'\n--- Convergence analysis ---')
    print(f'  R(c) = {d.R:.10f}')
    print(f'  rho = 1/R = {d.rho:.10f}')
    if d.R > 1:
        print(f'  *** CONVERGENT: R > 1 — shadow tower converges at t=1 ***')
        print(f'  Rate: geometric convergence with ratio 1/R = {1/d.R:.10f}')
    else:
        print(f'  *** DIVERGENT: R < 1 — shadow tower DIVERGES at t=1 ***')
        print(f'  Rate: geometric divergence with ratio 1/R = {1/d.R:.10f}')
        print(f'  Borel resummation required for physical interpretation')

    # Borel singularities
    bs = borel_singularities(c_val)
    print(f'\n--- Borel singularities (instanton actions) ---')
    print(f'  A_+ = {bs["A_plus"]}')
    print(f'  A_- = {bs["A_minus"]}')
    print(f'  |A_+| = {bs["A_plus_mod"]:.10f}')
    print(f'  Stokes direction = {bs["stokes_direction"]:.10f} rad '
          f'({bs["stokes_direction"] * 180 / math.pi:.4f} deg)')

    # Darboux asymptotics
    dd = darboux_coefficients(c_val)
    print(f'\n--- Darboux asymptotics: S_r ~ -A * rho^r * r^(-5/2) * cos(r*omega + phi) ---')
    print(f'  A (amplitude) = {dd.amplitude:.10f}')
    print(f'  rho (growth rate) = {dd.rho:.10f}')
    print(f'  omega (oscillation) = {dd.omega:.10f}')
    print(f'  phi (phase) = {dd.phase:.10f}')

    # Shadow coefficients
    coeffs = shadow_coefficients(c_val, r_max)
    print(f'\n--- Shadow coefficients S_2, ..., S_{min(r_max,15)} ---')
    for r in range(2, min(r_max + 1, 16)):
        sr = coeffs[r]
        asym = asymptotic_prediction(c_val, r) if r >= 6 else None
        if asym is not None:
            print(f'  S_{r:2d} = {sr:20.12f}    (Darboux: {asym:20.12f},'
                  f' ratio: {sr/asym:.6f})')
        else:
            print(f'  S_{r:2d} = {sr:20.12f}')

    # Convergence of partial sums at t=1
    partial = 0.0
    print(f'\n--- Partial sums G_N(1) = sum_{{r=2}}^N S_r ---')
    for N in [5, 10, 15, 20, 25, 30]:
        if N > r_max:
            break
        partial = sum(coeffs[r] for r in range(2, N + 1))
        print(f'  G_{N:2d}(1) = {partial:20.12f}')

    return d, coeffs


def analyze_m_operations(label: str, c_val: float):
    """Compute m_2, m_3, m_4 at specific spectral parameter values."""
    print(f'\n--- A_infinity operations at {label} ---')

    # m_2 at lambda = 1
    m2 = mk_exact_numerical(2, [1.0], c_val)
    print(f'\n  m_2(T,T; lambda=1):')
    print(f'    dT coeff  = {m2["dT_coeff"]}')
    print(f'    T coeff   = {m2["T_coeff"]}')
    print(f'    scalar    = {m2["scalar_coeff"]}')

    # m_3 at lambda = (1,1)
    m3 = mk_exact_numerical(3, [1.0, 1.0], c_val)
    print(f'\n  m_3(T,T,T; lambda=(1,1)):')
    print(f'    d2T coeff = {m3["d2T_coeff"]}')
    print(f'    dT coeff  = {m3["dT_coeff"]}')
    print(f'    T coeff   = {m3["T_coeff"]}')
    print(f'    scalar    = {m3["scalar_coeff"]}')

    # m_4 at lambda = (1,1,1)
    m4 = mk_exact_numerical(4, [1.0, 1.0, 1.0], c_val)
    print(f'\n  m_4(T,T,T,T; lambda=(1,1,1)):')
    for key in ['d2T_coeff', 'dT_coeff', 'T_coeff', 'scalar_coeff']:
        if key in m4:
            print(f'    {key:12s} = {m4[key]}')
    print(f'    nonvanishing = {m4["nonvanishing"]}')

    return m2, m3, m4


# =====================================================================
# Section 2: Exact rational computation at c = 1/2 and c = 1
# =====================================================================

def exact_rational_analysis():
    """Compute exact rational shadow coefficients at c=1/2 and c=1."""
    print_header('Exact rational shadow coefficients')

    for c_num, c_den, label in [(1, 2, 'c = 1/2 (Ising)'),
                                  (1, 1, 'c = 1 (free boson)')]:
        print(f'\n--- {label} ---')
        try:
            exact = shadow_coefficients_fraction(c_num, c_den, r_max=20)
            for r in range(2, 16):
                if r in exact:
                    print(f'  S_{r:2d} = {exact[r]}  = {float(exact[r]):.15f}')
        except Exception as e:
            print(f'  (Exact computation failed: {e})')
            print(f'  Falling back to float:')
            coeffs = shadow_coefficients(c_num / c_den, r_max=20)
            for r in range(2, 16):
                print(f'  S_{r:2d} = {coeffs[r]:.15f}')


# =====================================================================
# Section 3: Convergence threshold analysis
# =====================================================================

def convergence_threshold():
    """Find the exact central charge c* where R(c) = 1."""
    print_header('Convergence threshold: R(c*) = 1')
    print('Scanning c from 0.1 to 30...')

    # Binary search for R(c) = 1
    c_lo, c_hi = 0.1, 30.0
    for _ in range(100):
        c_mid = (c_lo + c_hi) / 2
        d = VirasoroShadowData(c_mid)
        if d.R > 1:
            c_hi = c_mid
        else:
            c_lo = c_mid

    c_star = (c_lo + c_hi) / 2
    d_star = VirasoroShadowData(c_star)
    print(f'\n  c* = {c_star:.10f}')
    print(f'  R(c*) = {d_star.R:.10f}')

    # Table of R(c) across the range
    print(f'\n  Table: R(c) at selected values')
    print(f'  {"c":>8s}  {"R(c)":>12s}  {"rho":>12s}  {"convergent?":>12s}')
    for c_val in [0.5, 1.0, 2.0, 4.0, 6.0, c_star, 7.0, 10.0, 13.0, 25.0, 26.0]:
        d = VirasoroShadowData(c_val)
        conv = 'YES' if d.R > 1 else 'NO'
        label = f' <-- c*' if abs(c_val - c_star) < 0.001 else ''
        print(f'  {c_val:8.4f}  {d.R:12.6f}  {d.rho:12.6f}  {conv:>12s}{label}')

    return c_star


# =====================================================================
# Section 4: Ising model bar complex analysis
# =====================================================================

def ising_analysis():
    """Special analysis of the Ising model (c=1/2).

    The Ising model is the first minimal model M(3,4) with exactly 3
    primaries: 1 (h=0), sigma (h=1/16), epsilon (h=1/2).
    The bar complex on the Verma module vs the irreducible module
    should differ: the irreducible has null vectors removed.
    """
    print_header('Ising model c=1/2: finite representation theory')

    c = 0.5
    print(f'\n  Central charge: c = 1/2')
    print(f'  Kac table: M(3,4) minimal model')
    print(f'  Primaries: 1 (h=0), sigma (h=1/16), epsilon (h=1/2)')
    print(f'  Fusion rules: sigma x sigma = 1 + epsilon')
    print(f'                sigma x epsilon = sigma')
    print(f'                epsilon x epsilon = 1')

    print(f'\n  The Virasoro bar complex at c=1/2 is computed on the')
    print(f'  stress-tensor T (the generator of Vir_{c}). The shadow tower')
    print(f'  S_r encodes the higher A_infinity operations m_k via the')
    print(f'  shadow generating function.')

    # Null vector at level 2: at c=1/2, the Verma module V_{1/2,0}
    # has a null vector at level 2
    print(f'\n  Null vector structure at c=1/2:')
    print(f'    Level 1: L_{-1}|0> (always null in vacuum module)')
    print(f'    Level 2: (L_{-2} - 3/(2(2*0+1)) L_{-1}^2)|0> ... ')
    print(f'    The Kac determinant at level 2: det_2 = c/2 * (5c+22)/10')
    det2 = c / 2 * (5 * c + 22) / 10
    print(f'    det_2(c=1/2) = {det2} {"(nondegenerate)" if abs(det2) > 1e-10 else "(DEGENERATE)"}')

    print(f'\n  The irreducible quotient of the vacuum Verma module at c=1/2')
    print(f'  has finite-dimensional weight spaces but is still infinite-dimensional.')
    print(f'  The bar complex B(Vir_{{1/2}}) is infinite, but the null vectors')
    print(f'  create additional relations that could (in principle) truncate')
    print(f'  the shadow tower at finite order.')

    # What the shadow tower actually sees
    print(f'\n  What the shadow tower computes at c=1/2:')
    print(f'  The m_k operations are defined on the VERMA module (universal envelope).')
    print(f'  On the IRREDUCIBLE module, the null vectors impose additional')
    print(f'  constraints. The bar complex of the irreducible module is a')
    print(f'  QUOTIENT of the bar complex of the Verma module by the null ideal.')
    print(f'  For c=1/2, the quotient has fewer generators but the shadow tower')
    print(f'  (computed here on the Verma/universal level) is an upper bound.')


# =====================================================================
# Section 5: Self-dual analysis at c=13
# =====================================================================

def self_dual_analysis():
    """Analyze the self-dual structure at c=13."""
    print_header('Self-dual point c=13: Koszul involution = identity')

    c = 13.0
    d = VirasoroShadowData(c)

    print(f'\n  Vir_13^! = Vir_{{26-13}} = Vir_13  (SELF-DUAL)')
    print(f'  kappa = c/2 = {c/2}')
    print(f'  kappa^! = (26-c)/2 = {(26-c)/2}')
    print(f'  kappa + kappa^! = {c/2 + (26-c)/2}  (cf. AP24: NOT zero for Virasoro)')
    print(f'  kappa_eff = 13 - 13 = 0  (NO! kappa_eff = kappa(matter) + kappa(ghost))')
    print(f'    Correction: kappa_eff depends on the FULL theory, not Vir alone.')
    print(f'    For pure Virasoro at c=13: the genus tower has kappa_eff = 13 - 26/2 = 0?')
    print(f'    No: kappa_eff = c_total - 26 = 13 - 26 = -13 \\neq 0.')
    print(f'    The genus tower does NOT vanish at c=13.')

    print(f'\n  Self-duality consequences for the bar complex:')
    print(f'  B(Vir_13) is self-dual: there exists an isomorphism')
    print(f'    B(Vir_13) \\simeq B(Vir_13)^\\vee')
    print(f'  This gives a nondegenerate pairing on the bar complex,')
    print(f'  which constrains the m_k operations to satisfy')
    print(f'  cyclic A_infinity relations (Kontsevich-Soibelman).')

    print(f'\n  The shadow tower at c=13 is:')
    coeffs = shadow_coefficients(c, r_max=20)
    for r in range(2, 16):
        print(f'    S_{r:2d} = {coeffs[r]:.15f}')

    # Check: does self-duality impose S_r = f(S_r) for some involution?
    # The Koszul involution c -> 26-c sends S_4(c) -> S_4(26-c)
    # At c=13 this is a fixed point.
    c_dual = 26 - c
    d_dual = VirasoroShadowData(c_dual)
    print(f'\n  Self-duality check: S_4(c) vs S_4(26-c)')
    print(f'    S_4(13)   = {d.S4:.15f}')
    print(f'    S_4(13)   = {d_dual.S4:.15f}  (should match)')
    print(f'    Difference = {abs(d.S4 - d_dual.S4):.2e}')

    print(f'\n  The shadow metric at c=13:')
    print(f'    Q_L(t) = {d.q0} + {d.q1}*t + {d.q2:.10f}*t^2')
    print(f'    Discriminant = {d.q1**2 - 4*d.q0*d.q2:.10f}')
    print(f'    Branch points: t_+ = {d.t_plus}')
    print(f'                   t_- = {d.t_minus}')
    print(f'    R(13) = {d.R:.10f}')
    print(f'    The tower CONVERGES (R > 1).')


# =====================================================================
# Section 6: Critical string analysis at c=25 and c=26
# =====================================================================

def critical_string_analysis():
    """Analyze c=25 (near-critical) and c=26 (critical string)."""
    print_header('Near-critical and critical string: c=25, c=26')

    for c in [25.0, 26.0]:
        d = VirasoroShadowData(c)
        print(f'\n  --- c = {c} ---')
        print(f'  kappa = {d.kappa}')
        print(f'  kappa^! = {(26-c)/2}')

        if c == 26:
            print(f'  *** CRITICAL STRING: kappa^! = 0 ***')
            print(f'  Vir_26^! = Vir_0 (the trivial Virasoro)')
            print(f'  The genus tower F_g vanishes: kappa_eff = 0')
            print(f'  BUT the shadow tower S_r persists (S_r != 0)')
            print(f'  The shadow tower encodes the PERTURBATIVE STRING AMPLITUDES')
            print(f'  decoupled from the genus expansion.')
        if c == 25:
            print(f'  Just below critical: c = 26 - 1')
            print(f'  The genus tower has small kappa_eff = -1/2')
            print(f'  The shadow tower and genus tower are both present.')

        print(f'  R(c) = {d.R:.10f}')
        print(f'  rho = {d.rho:.10f}')
        print(f'  Convergent: {"YES" if d.R > 1 else "NO"}')

        if d.R > 1:
            coeffs = shadow_coefficients(c, r_max=20)
            partial = sum(coeffs[r] for r in range(2, 21))
            print(f'  G_20(1) = {partial:.15f}')
            print(f'  Rapid convergence: |S_20| = {abs(coeffs[20]):.2e}')


# =====================================================================
# Section 7: Comparative table
# =====================================================================

def comparative_table():
    """Print a comprehensive comparison table."""
    print_header('COMPARATIVE TABLE: Shadow invariants at special central charges')

    c_vals = [0.5, 1.0, 13.0, 25.0, 26.0]
    labels = ['Ising', 'Free boson', 'Self-dual', 'Near-crit', 'Crit string']

    # Header
    print(f'\n  {"":>12s}', end='')
    for lab in labels:
        print(f'  {lab:>14s}', end='')
    print()
    print(f'  {"":>12s}', end='')
    for c in c_vals:
        print(f'  {"c="+str(c):>14s}', end='')
    print()
    print('  ' + '-' * (12 + 16 * len(c_vals)))

    rows = []
    for c in c_vals:
        d = VirasoroShadowData(c)
        coeffs = shadow_coefficients(c, r_max=10)
        rows.append({
            'kappa': d.kappa,
            'kappa!': (26 - c) / 2,
            'S_2': c / 2,
            'S_3': 2.0,
            'S_4': d.S4,
            'S_5': coeffs.get(5, 0.0),
            'R': d.R,
            'rho': d.rho,
            'convergent': 'YES' if d.R > 1 else 'NO',
            'q0': d.q0,
            'q1': d.q1,
            'q2': d.q2,
        })

    for key in ['kappa', 'kappa!', 'S_2', 'S_3', 'S_4', 'S_5',
                'R', 'rho', 'convergent', 'q0', 'q1', 'q2']:
        print(f'  {key:>12s}', end='')
        for row in rows:
            val = row[key]
            if isinstance(val, str):
                print(f'  {val:>14s}', end='')
            elif isinstance(val, float):
                if abs(val) > 100:
                    print(f'  {val:>14.4f}', end='')
                elif abs(val) > 0.01:
                    print(f'  {val:>14.6f}', end='')
                else:
                    print(f'  {val:>14.8f}', end='')
            print()


# =====================================================================
# Section 8: Convergence radius formula verification
# =====================================================================

def radius_formula_verification():
    """Verify the convergence radius formula.

    The user's formula: R = c * sqrt((5c+22)/(180c+872))

    Let's verify this against the actual R = |t_+| computation.
    """
    print_header('Convergence radius formula verification')

    print(f'\n  Claimed formula: R(c) = c * sqrt((5c+22)/(180c+872))')
    print(f'  Q_L zeros: t_+/- = (-q1 +/- sqrt(q1^2 - 4*q0*q2)) / (2*q2)')
    print(f'  q0=c^2, q1=12c, q2=(180c+872)/(5c+22)')

    print(f'\n  {"c":>6s}  {"R(actual)":>12s}  {"R(formula)":>12s}  {"match?":>8s}')
    for c in [0.5, 1.0, 2.0, 5.0, 10.0, 13.0, 25.0, 26.0, 50.0]:
        d = VirasoroShadowData(c)
        R_actual = d.R
        # User formula
        R_formula = c * math.sqrt((5*c + 22) / (180*c + 872))
        match = 'YES' if abs(R_actual - R_formula) / R_actual < 0.01 else 'NO'
        print(f'  {c:6.1f}  {R_actual:12.6f}  {R_formula:12.6f}  {match:>8s}')

    # What IS the correct closed-form?
    # R = |t_+| where t_+ = (-12c + sqrt(144c^2 - 4c^2 * (180c+872)/(5c+22))) / (2*(180c+872)/(5c+22))
    # disc = 144c^2 - 4c^2*(180c+872)/(5c+22) = c^2[144 - 4(180c+872)/(5c+22)]
    # = c^2 * [144(5c+22) - 4(180c+872)] / (5c+22)
    # = c^2 * [720c + 3168 - 720c - 3488] / (5c+22)
    # = c^2 * [-320] / (5c+22)
    # = -320*c^2 / (5c+22)
    # So disc < 0 always (for c > 0): branch points are ALWAYS complex conjugate.
    # sqrt(disc) = i * c * sqrt(320/(5c+22)) = i * c * 8*sqrt(5/(5c+22))
    #
    # t_+/- = (-12c +/- i*c*8*sqrt(5/(5c+22))) / (2*(180c+872)/(5c+22))
    #        = c * (-12 +/- 8i*sqrt(5/(5c+22))) * (5c+22) / (2*(180c+872))
    #
    # |t_+|^2 = c^2 * (5c+22)^2 / (4*(180c+872)^2) * (144 + 320/(5c+22))
    #         = c^2 * (5c+22)^2 / (4*(180c+872)^2) * (144(5c+22) + 320) / (5c+22)
    #         = c^2 * (5c+22) * (720c + 3168 + 320) / (4*(180c+872)^2)
    #         = c^2 * (5c+22) * (720c + 3488) / (4*(180c+872)^2)
    #         = c^2 * (5c+22) * 16*(45c + 218) / (4*(180c+872)^2)
    #         Hmm, 720c + 3488 = 4*(180c + 872)
    #         So |t_+|^2 = c^2 * (5c+22) * 4*(180c+872) / (4*(180c+872)^2)
    #                    = c^2 * (5c+22) / (180c+872)
    #
    # THEREFORE: R = c * sqrt((5c+22)/(180c+872))  ✓
    #
    # The formula IS correct!

    print(f'\n  Derivation:')
    print(f'  disc(Q_L) = -320*c^2/(5c+22) < 0 always (for c>0)')
    print(f'  So branch points are ALWAYS complex conjugate.')
    print(f'  |t_+|^2 = c^2 * (5c+22)/(180c+872)')
    print(f'  R(c) = c * sqrt((5c+22)/(180c+872))  [EXACT]')

    # Now recheck user values
    print(f'\n  User-provided values check:')
    for c, claimed_R in [(1.0, 0.160), (0.5, 0.119)]:
        R = c * math.sqrt((5*c+22)/(180*c+872))
        print(f'  c={c}: R = {R:.6f} (user claimed ~{claimed_R})')

    # At c=25
    c = 25.0
    R_25 = c * math.sqrt((5*c+22)/(180*c+872))
    print(f'  c=25: R = {R_25:.6f} (user claimed ~4.14)')

    # c=26
    c = 26.0
    R_26 = c * math.sqrt((5*c+22)/(180*c+872))
    print(f'  c=26: R = {R_26:.6f}')


# =====================================================================
# Main
# =====================================================================

if __name__ == '__main__':
    # Section 1: Full analysis at each special charge
    for label, c_val in SPECIAL_CHARGES.items():
        analyze_shadow_data(label, c_val, r_max=30)
        analyze_m_operations(label, c_val)

    # Section 2: Exact rational computation
    exact_rational_analysis()

    # Section 3: Convergence threshold
    c_star = convergence_threshold()

    # Section 4: Ising special analysis
    ising_analysis()

    # Section 5: Self-dual analysis
    self_dual_analysis()

    # Section 6: Critical string
    critical_string_analysis()

    # Section 7: Comparative table
    comparative_table()

    # Section 8: Radius formula
    radius_formula_verification()

    print('\n\n' + '=' * 72)
    print('  COMPUTATION COMPLETE')
    print('=' * 72)
