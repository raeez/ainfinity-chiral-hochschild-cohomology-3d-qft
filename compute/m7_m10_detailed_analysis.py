r"""Detailed analysis of the period-2 pattern and secondary vanishing in m_2 through m_10.

Key findings from m7_m10_depth_frontier.py:
  - STRUCTURAL GAP at d=k confirmed for ALL k=2,...,10 (Theorem thm:gap-migration(ii))
  - PERIOD-2 PATTERN: Even arities (4,6,8,10) have secondary vanishing at depths 0 and 1.
    Odd arities (3,5,7,9) fill all depths 0,...,k-2.
  - The n=4 anomaly is NOT isolated: it repeats at EVERY even arity.
  - At symmetric point λ_i = λ: m_k|_T VANISHES IDENTICALLY for even k.
  - Scalar P_k(1,...,1) vanishes for even k >= 4.

This script performs deeper investigation:
  1. Exact depth spectra at more values of c
  2. Verify the even-arity palindrome vanishing at the symmetric point
  3. Test whether m_{2r}|_T factors through (l_1 - l_{2r-1}) generically
  4. Compute exact scalar shadow values and compare with shadow_borel_resurgence
  5. L^1 norm growth analysis for Gevrey class determination
"""

from __future__ import annotations
import sys
import os
import time
import random
import math
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from compute.m7_m10_depth_frontier import StasheffEngine, extract_depth_spectrum


def verify_even_palindrome(engine: StasheffEngine, k: int, n_tests: int = 30) -> Dict:
    """Test whether m_k|_T vanishes identically at the symmetric point l_i = l.

    At the symmetric point all l_i = l, the ENTIRE field-dependent sector
    should vanish for even k (if the palindrome pattern holds).
    """
    rng = random.Random(12345 + k)
    max_T_total = 0.0
    max_scalar = 0.0

    for _ in range(n_tests):
        lam = rng.uniform(-5.0, 5.0)
        lams = tuple(lam for _ in range(k - 1))
        result = engine.mk(lams)

        # Sum of all T-dependent field contributions
        T_total = sum(abs(v) for f, v in result.items() if f >= 0)
        max_T_total = max(max_T_total, T_total)
        max_scalar = max(max_scalar, abs(result.get(-1, 0.0)))

    return {
        'k': k,
        'max_T_at_symmetric': max_T_total,
        'T_vanishes_at_symmetric': max_T_total < 1e-8,
        'max_scalar_at_symmetric': max_scalar,
        'scalar_vanishes_at_symmetric': max_scalar < 1e-8,
    }


def test_antisymmetric_factor(engine: StasheffEngine, k: int, n_tests: int = 50) -> Dict:
    """Test whether m_k|_T factors through (l_1 - l_{k-1}).

    For even k: check if m_k(l_1,...,l_{k-1}) vanishes when l_1 = l_{k-1}.
    """
    rng = random.Random(54321 + k)
    max_T_at_constraint = 0.0
    max_T_generic = 0.0

    for _ in range(n_tests):
        # Generic point
        lams_generic = tuple(rng.uniform(0.1, 5.0) for _ in range(k - 1))
        result_gen = engine.mk(lams_generic)
        T_gen = sum(abs(v) for f, v in result_gen.items() if f >= 0)
        max_T_generic = max(max_T_generic, T_gen)

        # Constrained: l_1 = l_{k-1}
        lams_list = list(lams_generic)
        lams_list[-1] = lams_list[0]
        lams_constrained = tuple(lams_list)
        result_con = engine.mk(lams_constrained)
        T_con = sum(abs(v) for f, v in result_con.items() if f >= 0)
        max_T_at_constraint = max(max_T_at_constraint, T_con)

    return {
        'k': k,
        'max_T_generic': max_T_generic,
        'max_T_at_l1_eq_lk': max_T_at_constraint,
        'factors_through_l1_minus_lk': max_T_at_constraint < 1e-6 * max_T_generic,
    }


def compute_L1_norms_detailed(c_val: float, max_arity: int = 10,
                               n_samples: int = 500) -> Dict[int, float]:
    """Compute L^1 norms with more samples for better Gevrey analysis."""
    engine = StasheffEngine(c_val)
    rng = random.Random(99999)
    norms = {}

    for k in range(2, max_arity + 1):
        engine._cache.clear()
        total = 0.0
        for _ in range(n_samples):
            lams = tuple(rng.uniform(-1.0, 1.0) for _ in range(k - 1))
            result = engine.mk(lams)
            T_coeff = abs(result.get(0, 0.0))
            total += T_coeff
        norms[k] = total / n_samples
    return norms


def scalar_polynomial_at_unit(c_val: float, max_arity: int = 10) -> Dict[int, float]:
    """Compute P_k(1,...,1) for all arities (the shadow polynomial at unit point)."""
    engine = StasheffEngine(c_val)
    results = {}
    for k in range(2, max_arity + 1):
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(k - 1))
        result = engine.mk(lams)
        scalar = result.get(-1, 0.0)
        P_k = scalar / (c_val / 12.0) if abs(c_val) > 1e-14 else 0.0
        results[k] = P_k
    return results


def shadow_comparison(c_val: float, max_arity: int = 10):
    """Compare scalar tower with shadow coefficients from shadow_borel_resurgence."""
    from compute.lib.shadow_borel_resurgence import shadow_coefficients

    shadow = shadow_coefficients(c_val, r_max=max_arity + 2)
    unit_polys = scalar_polynomial_at_unit(c_val, max_arity)

    print(f"\n  Shadow comparison at c = {c_val}:")
    print(f"  {'k':>4} {'P_k(1,...,1)':>18} {'S_{k-1}':>18} {'S_k':>18}")
    for k in range(2, max_arity + 1):
        P_k = unit_polys.get(k, 0.0)
        S_k = shadow.get(k, 0.0)
        S_km1 = shadow.get(k - 1, 0.0)
        print(f"  {k:>4} {P_k:>18.6f} {S_km1:>18.10f} {S_k:>18.10f}")


def test_depth_0_1_even(engine: StasheffEngine, k: int, n_tests: int = 100) -> Dict:
    """More thorough test of depth 0 and depth 1 vanishing at even arities.

    Depth 0 = weight k-1 = d^{k-2} T coefficient
    Depth 1 = weight k-2 = d^{k-3} T coefficient
    """
    rng = random.Random(11111 + k)
    max_d0 = 0.0  # depth 0 = deriv order k-2
    max_d1 = 0.0  # depth 1 = deriv order k-3
    max_d2 = 0.0  # depth 2 = deriv order k-4

    d0_order = k - 2  # d^{k-2} T
    d1_order = k - 3  # d^{k-3} T
    d2_order = k - 4  # d^{k-4} T

    for _ in range(n_tests):
        lams = tuple(rng.uniform(-3.0, 3.0) for _ in range(k - 1))
        result = engine.mk(lams)
        max_d0 = max(max_d0, abs(result.get(d0_order, 0.0)))
        if d1_order >= 0:
            max_d1 = max(max_d1, abs(result.get(d1_order, 0.0)))
        if d2_order >= 0:
            max_d2 = max(max_d2, abs(result.get(d2_order, 0.0)))

    return {
        'k': k,
        'depth_0_max': max_d0,
        'depth_0_vanishes': max_d0 < 1e-8,
        'depth_1_max': max_d1,
        'depth_1_vanishes': max_d1 < 1e-8,
        'depth_2_max': max_d2,
        'depth_2_vanishes': max_d2 < 1e-8,
    }


def run_detailed_analysis():
    print("=" * 80)
    print("DETAILED ANALYSIS: PERIOD-2 PATTERN AND SECONDARY VANISHING")
    print("=" * 80)

    c_val = 1.0
    engine = StasheffEngine(c_val)

    # 1. Verify even-arity palindrome vanishing at symmetric point
    print("\n--- 1. Even-arity palindrome vanishing at symmetric point ---")
    for k in range(2, 11):
        engine._cache.clear()
        pal = verify_even_palindrome(engine, k)
        label = "EVEN" if k % 2 == 0 else "odd"
        T_status = "VANISHES" if pal['T_vanishes_at_symmetric'] else f"nonzero ({pal['max_T_at_symmetric']:.3e})"
        sc_status = "VANISHES" if pal['scalar_vanishes_at_symmetric'] else f"nonzero ({pal['max_scalar_at_symmetric']:.3e})"
        print(f"  m_{k:2d} ({label:4s}): T|_sym {T_status:40s}  scalar|_sym {sc_status}")

    # 2. Test (l_1 - l_{k-1}) factor for even k
    print("\n--- 2. Antisymmetric factor (l_1 - l_{k-1}) test ---")
    for k in [4, 6, 8, 10]:
        engine._cache.clear()
        af = test_antisymmetric_factor(engine, k, n_tests=30)
        factors = "YES" if af['factors_through_l1_minus_lk'] else "NO"
        print(f"  m_{k:2d}: factors through (l_1 - l_{k-1})? {factors}  "
              f"(T|_constrained = {af['max_T_at_l1_eq_lk']:.3e}, "
              f"T|_generic = {af['max_T_generic']:.3e})")

    # 3. Detailed depth 0/1 vanishing at even arities
    print("\n--- 3. Depth 0 and 1 vanishing at even arities ---")
    for k in range(2, 11):
        engine._cache.clear()
        d01 = test_depth_0_1_even(engine, k)
        label = "EVEN" if k % 2 == 0 else "odd"
        d0_val = 'ZERO' if d01['depth_0_vanishes'] else f"{d01['depth_0_max']:.3e}"
        d1_val = 'ZERO' if d01['depth_1_vanishes'] else f"{d01['depth_1_max']:.3e}"
        d2_val = 'ZERO' if d01['depth_2_vanishes'] else f"{d01['depth_2_max']:.3e}"
        d0 = f"depth_0={d0_val}"
        d1 = f"depth_1={d1_val}"
        d2 = f"depth_2={d2_val}"
        print(f"  m_{k:2d} ({label:4s}): {d0:25s}  {d1:25s}  {d2}")

    # 4. L^1 norm growth for Gevrey analysis
    print("\n--- 4. L^1 norm growth (Gevrey-1 analysis) ---")
    for c_val in [1.0, 13.0, 26.0]:
        print(f"\n  c = {c_val}:")
        norms = compute_L1_norms_detailed(c_val, max_arity=10, n_samples=300)
        print(f"  {'k':>4} {'||m_k|_T||':>16} {'ratio to (k-1)!':>20} {'ratio_{k/k-1}':>16}")
        for k in range(2, 11):
            norm = norms[k]
            factorial_ratio = norm / math.factorial(k - 1) if math.factorial(k - 1) > 0 else 0.0
            if k > 2 and norms.get(k-1, 0.0) > 1e-300:
                growth_ratio = norm / norms[k-1]
            else:
                growth_ratio = float('nan')
            print(f"  {k:>4} {norm:>16.6e} {factorial_ratio:>20.8f} {growth_ratio:>16.4f}")

    # 5. Scalar polynomial P_k(1,...,1) pattern
    print("\n--- 5. Scalar polynomial P_k(1,...,1) ---")
    for c_val in [1.0, 13.0, 26.0]:
        P = scalar_polynomial_at_unit(c_val, 10)
        print(f"\n  c = {c_val}:")
        for k in range(2, 11):
            v = P.get(k, 0.0)
            label = "ZERO" if abs(v) < 1e-6 else f"{v:.6f}"
            # Check if P_k(1,...,1) = (-1)^{k/2} * (k-1)!! or similar
            if abs(v) > 1e-6:
                fac_ratio = v / math.factorial(k)
                print(f"    P_{k}(1,...,1) = {label:>20s}   ratio to k! = {fac_ratio:.8f}")
            else:
                print(f"    P_{k}(1,...,1) = {label:>20s}")

    # 6. Shadow coefficient comparison
    print("\n--- 6. Shadow coefficient comparison ---")
    try:
        for c_val in [1.0, 13.0, 26.0]:
            shadow_comparison(c_val, 10)
    except Exception as e:
        print(f"  Shadow comparison failed: {e}")

    # 7. c-independence verification at specific depths
    print("\n--- 7. c-independence of T-dependent sector ---")
    print("  Checking that T-dependent coefficients are the same at c=1, c=13, c=26, c=100:")
    for k in range(2, 9):
        rng = random.Random(77777 + k)
        lams = tuple(rng.uniform(0.5, 2.0) for _ in range(k - 1))

        T_coeffs = {}
        for c_val in [1.0, 13.0, 26.0, 100.0]:
            eng = StasheffEngine(c_val)
            result = eng.mk(lams)
            T_coeffs[c_val] = result.get(0, 0.0)

        # Check if all T-coefficients are the same
        vals = list(T_coeffs.values())
        spread = max(abs(v - vals[0]) for v in vals) if vals else 0.0
        max_val = max(abs(v) for v in vals) if vals else 0.0
        rel_spread = spread / max_val if max_val > 1e-300 else 0.0
        status = "c-INDEPENDENT" if rel_spread < 1e-10 else f"c-DEPENDENT (spread={rel_spread:.3e})"
        print(f"  m_{k:2d} T-coeff: {status}")

    # 8. FINAL CONSOLIDATED TABLE
    print(f"\n\n{'=' * 80}")
    print("CONSOLIDATED DEPTH SPECTRUM TABLE: m_2 through m_10")
    print(f"{'=' * 80}")
    print()
    print(f"{'k':>3} {'parity':>7} {'Spec(m_k|_T)':>35} {'gap@k':>7} {'d_min':>6} "
          f"{'0,1 vanish':>11} {'scalar':>8} {'sym vanish':>11}")
    print("-" * 100)

    engine = StasheffEngine(1.0)
    for k in range(2, 11):
        engine._cache.clear()
        spec = extract_depth_spectrum(engine, k, n_samples=30, seed=42+k)
        d01 = test_depth_0_1_even(engine, k, n_tests=30)
        pal = verify_even_palindrome(engine, k, n_tests=20)

        parity = "even" if k % 2 == 0 else "odd"
        gap = "YES" if spec['gap_at_k'] else "NO"
        d01_str = "YES" if (d01['depth_0_vanishes'] and d01['depth_1_vanishes']) else "no"
        sc = f"d={spec['scalar_depth']}" if spec['scalar_present'] else "-"
        sym = "YES" if pal['T_vanishes_at_symmetric'] else "no"

        print(f"{k:>3} {parity:>7} {str(spec['depths_T']):>35} {gap:>7} "
              f"{spec['min_depth_T']:>6} {d01_str:>11} {sc:>8} {sym:>11}")


if __name__ == '__main__':
    run_detailed_analysis()
