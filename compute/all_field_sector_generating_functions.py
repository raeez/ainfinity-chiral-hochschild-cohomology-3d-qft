r"""All-field-sector generating functions for the Virasoro A_infinity shadow tower.

GOAL: Extend the ALGEBRAIC generating function
  g(x) = 1/2 - (1+8x)/(2*sqrt(1+4x))
from the T-coefficient (w=0) to ALL field sectors (w=0,1,2,...).

COMPUTATION PLAN:
  (1) dT_k(1,...,1) for odd k=3,...,25 — generating function for w=1 sector
  (2) d^2 T_k(1,...,1) for odd k=5,...,25 — generating function for w=2 sector
  (3) Ratio dT/T at the symmetric point: simple function of k?
  (4) Normalized coefficients f_w(r) = [d^w T coeff at sym pt] / ((2r)! * (-1)^{r+1} * C_{r-1})
  (5) All-field generating function G(x,y) = sum_{k,w} [coeff] * x^{(k-3)/2} * y^w

Dependencies: m7_m10_depth_frontier.py (StasheffEngine)
"""

from __future__ import annotations

import sys
import os
import math
from fractions import Fraction
from typing import Dict, List, Tuple, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))

from m7_m10_depth_frontier import StasheffEngine


def catalan(n: int) -> int:
    """Catalan number C_n = binom(2n,n)/(n+1)."""
    if n < 0:
        return 0
    return math.comb(2 * n, n) // (n + 1)


def compute_all_field_sectors(max_r: int = 12, c_val: float = 0.0) -> Dict:
    r"""Compute all field sectors at the symmetric point for odd arities.

    For each r = 1, 2, ..., max_r (arity k = 2r+1):
      Compute m_{2r+1}(T,...,T; 1,...,1) and extract coefficient of d^w T
      for each w = 0, 1, 2, ...

    Returns dict[r] -> dict[w] -> coefficient
    """
    engine = StasheffEngine(c_val)
    results = {}

    for r in range(1, max_r + 1):
        k = 2 * r + 1
        n_lams = k - 1
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(n_lams))
        result = engine.mk(lams)

        # Extract coefficients by derivative order w
        sector = {}
        for w in range(-1, k + 1):
            v = result.get(w, 0.0)
            if abs(v) > 1e-6:
                sector[w] = v
        results[r] = sector

    return results


def compute_even_field_sectors(max_r: int = 12, c_val: float = 0.0) -> Dict:
    r"""Compute all field sectors at the symmetric point for EVEN arities.

    For each r = 2, 3, ..., max_r (arity k = 2r):
    """
    engine = StasheffEngine(c_val)
    results = {}

    for r in range(2, max_r + 1):
        k = 2 * r
        n_lams = k - 1
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(n_lams))
        result = engine.mk(lams)

        sector = {}
        for w in range(-1, k + 1):
            v = result.get(w, 0.0)
            if abs(v) > 1e-6:
                sector[w] = v
        results[r] = sector

    return results


def main():
    print("=" * 80)
    print("ALL-FIELD-SECTOR GENERATING FUNCTIONS")
    print("=" * 80)

    # ===================================================================
    # Part 0: Compute raw data at c=0 (no scalar sector, pure field)
    # ===================================================================
    max_r = 8  # r=8 => k=17; higher arities get very slow
    print(f"\nComputing field sectors for r=1..{max_r} (odd arities k=3..{2*max_r+1}), c=0...")
    data = compute_all_field_sectors(max_r=max_r, c_val=0.0)

    print("\n" + "=" * 80)
    print("PART 1: RAW COEFFICIENTS AT SYMMETRIC POINT (c=0)")
    print("=" * 80)

    for r in sorted(data):
        k = 2 * r + 1
        sector = data[r]
        nonzero = {w: v for w, v in sector.items() if w >= 0}
        print(f"\n  r={r:2d} (k={k:2d}):")
        for w in sorted(nonzero):
            print(f"    d^{w}T coeff = {nonzero[w]:30.6f}")

    # ===================================================================
    # Part 1: d^w T coefficients — extract and display
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 1a: TABLE OF d^w T COEFFICIENTS BY SECTOR")
    print("=" * 80)

    # Find max w that appears
    max_w_seen = 0
    for r in data:
        for w in data[r]:
            if w >= 0 and abs(data[r][w]) > 1e-6:
                max_w_seen = max(max_w_seen, w)

    header = f"{'r':>3s} {'k':>3s}"
    for w in range(max_w_seen + 1):
        header += f" {'d^'+str(w)+'T':>22s}"
    print(header)
    print("-" * len(header))

    for r in sorted(data):
        k = 2 * r + 1
        line = f"{r:3d} {k:3d}"
        for w in range(max_w_seen + 1):
            val = data[r].get(w, 0.0)
            if abs(val) < 1e-6:
                line += f" {'0':>22s}"
            else:
                line += f" {val:22.4f}"
        print(line)

    # ===================================================================
    # Part 2: Known T-coefficient (w=0) verification
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 2: T-COEFFICIENT (w=0) VERIFICATION")
    print("=" * 80)

    print(f"\n  {'r':>3s} {'numerical':>22s} {'formula':>22s} {'match':>8s}")
    for r in sorted(data):
        T_num = data[r].get(0, 0.0)
        sign = (-1) ** (r + 1)
        cat_val = catalan(r - 1)
        fact_val = math.factorial(2 * r)
        T_formula = sign * (2 * r + 1) * cat_val * fact_val
        match = abs(T_num - T_formula) < max(1.0, abs(T_formula) * 1e-8)
        print(f"  {r:3d} {T_num:22.4f} {T_formula:22d} {str(match):>8s}")

    # ===================================================================
    # Part 3: dT-coefficient (w=1) — search for generating function
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 3: dT-COEFFICIENT (w=1) — SEQUENCE ANALYSIS")
    print("=" * 80)

    print("\n  Sequence of dT coefficients at symmetric point:")
    dT_seq = []
    for r in sorted(data):
        val = data[r].get(1, 0.0)
        dT_seq.append((r, val))
        k = 2 * r + 1
        print(f"  r={r:2d} (k={k:2d}): dT coeff = {val:30.6f}")

    # Try normalizing by (2r)! * (-1)^{r+1} * C_{r-1}
    print("\n  Normalized dT coefficients: dT / [(-1)^{r+1} * C_{r-1} * (2r)!]")
    for r, val in dT_seq:
        sign = (-1) ** (r + 1)
        cat_val = catalan(r - 1)
        fact_val = math.factorial(2 * r)
        norm = sign * cat_val * fact_val
        if abs(norm) > 1e-10:
            ratio = val / norm
            print(f"  r={r:2d}: dT / norm = {ratio:22.8f}")
        else:
            print(f"  r={r:2d}: norm = 0 (skipped)")

    # ===================================================================
    # Part 3b: dT/T ratio at symmetric point
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 3b: dT/T RATIO AT SYMMETRIC POINT")
    print("=" * 80)

    print(f"\n  {'r':>3s} {'dT':>22s} {'T':>22s} {'dT/T':>16s}")
    for r in sorted(data):
        T_val = data[r].get(0, 0.0)
        dT_val = data[r].get(1, 0.0)
        if abs(T_val) > 1e-10:
            ratio = dT_val / T_val
            print(f"  {r:3d} {dT_val:22.4f} {T_val:22.4f} {ratio:16.8f}")
        else:
            print(f"  {r:3d} {dT_val:22.4f} {T_val:22.4f} {'---':>16s}")

    # Check if dT/T = a*r + b
    print("\n  Testing dT/T = linear in r or k:")
    ratios = []
    for r in sorted(data):
        T_val = data[r].get(0, 0.0)
        dT_val = data[r].get(1, 0.0)
        if abs(T_val) > 1e-10:
            ratios.append((r, dT_val / T_val))

    if len(ratios) >= 2:
        # Check successive differences
        print(f"  {'r':>3s} {'dT/T':>16s} {'diff':>16s} {'2nd diff':>16s}")
        for i, (r, ratio) in enumerate(ratios):
            diff_str = ""
            diff2_str = ""
            if i > 0:
                diff = ratio - ratios[i-1][1]
                diff_str = f"{diff:16.8f}"
            if i > 1:
                d1 = ratios[i][1] - ratios[i-1][1]
                d0 = ratios[i-1][1] - ratios[i-2][1]
                diff2_str = f"{d1 - d0:16.8f}"
            print(f"  {r:3d} {ratio:16.8f} {diff_str:>16s} {diff2_str:>16s}")

    # ===================================================================
    # Part 4: Normalized f_w(r) = [d^w T coeff] / [(2r)! * (-1)^{r+1} * C_{r-1}]
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 4: NORMALIZED COEFFICIENTS f_w(r)")
    print("  f_w(r) = [d^w T coeff at sym pt] / [(2r)! * (-1)^{r+1} * C_{r-1}]")
    print("=" * 80)

    # The w=0 case gives f_0(r) = (2r+1) by Theorem 1
    for w in range(max_w_seen + 1):
        print(f"\n  --- w={w} sector ---")
        fw_seq = []
        for r in sorted(data):
            val = data[r].get(w, 0.0)
            sign = (-1) ** (r + 1)
            cat_val = catalan(r - 1)
            fact_val = math.factorial(2 * r)
            norm = sign * cat_val * fact_val
            if abs(norm) > 1e-10:
                fw = val / norm
                fw_seq.append((r, fw))
                print(f"    r={r:2d}: f_{w}(r) = {fw:22.8f}")
            else:
                print(f"    r={r:2d}: (norm=0, skip)")

        # Try to recognize f_w(r) as a polynomial in r
        if len(fw_seq) >= 3:
            print(f"\n    Checking if f_{w}(r) is polynomial in r:")
            # Compute successive differences to determine degree
            vals = [v for _, v in fw_seq]
            rs = [r for r, _ in fw_seq]
            for deg in range(6):
                diffs = vals[:]
                for d in range(deg + 1):
                    diffs = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
                if len(diffs) >= 2:
                    # Check if constant
                    spread = max(abs(d - diffs[0]) for d in diffs)
                    if spread < 1e-4:
                        print(f"    f_{w}(r) has finite differences of degree {deg+1} ~ {diffs[0]:.6f}")
                        if spread < 1e-4:
                            print(f"    => f_{w}(r) is a polynomial of degree <= {deg+1} in r")
                            break
            else:
                print(f"    Not detected as polynomial of degree <= 6")

            # Attempt polynomial fit
            _try_polynomial_fit(w, fw_seq)

    # ===================================================================
    # Part 5: All-field generating function G(x,y)
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 5: ALL-FIELD GENERATING FUNCTION G(x,y)")
    print("  G(x,y) = sum_{r>=1, w>=0} f_w(r) * (-1)^{r+1} * C_{r-1} * x^r * y^w")
    print("  so that at y=0, G(x,0) = g(x) = 1/2 - (1+8x)/(2*sqrt(1+4x))")
    print("=" * 80)

    # The key question: can we write G(x,y) in closed form?
    # We know g_0(x) = sum_r (-1)^{r+1} * (2r+1) * C_{r-1} * x^r
    #                = 1/2 - (1+8x)/(2*sqrt(1+4x))
    # If f_w(r) = P_w(r) * (2r+1) or some polynomial, then
    # g_w(x) = sum_r (-1)^{r+1} * P_w(r) * (2r+1) * C_{r-1} * x^r
    # can be obtained by differential operators on g_0.

    # First: compute g_w(x) = sum_r (-1)^{r+1} * f_w(r) * C_{r-1} * x^r
    print("\n  Computing sector generating functions g_w(x):")
    for w in range(min(max_w_seen + 1, 6)):
        print(f"\n  --- g_{w}(x) coefficients ---")
        print(f"  {'r':>3s} {'(-1)^{r+1}*f_w*C_{r-1}':>28s}")
        for r in sorted(data):
            val = data[r].get(w, 0.0)
            sign = (-1) ** (r + 1)
            cat_val = catalan(r - 1)
            fact_val = math.factorial(2 * r)
            norm = sign * cat_val * fact_val
            if abs(norm) > 1e-10:
                fw = val / norm
                gw_coeff = (-1)**(r+1) * fw * catalan(r-1)
                print(f"  {r:3d} {gw_coeff:28.8f}")

    # ===================================================================
    # Part 5b: Try to identify G(x,y) algebraically
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 5b: IDENTIFYING THE ALGEBRAIC FORM OF G(x,y)")
    print("=" * 80)

    # Key insight: the generating function for sum_r (-1)^{r+1} * p(r) * C_{r-1} * x^r
    # where p(r) is polynomial in r can be computed by the identity:
    #   sum_r r^j * C_{r-1} * (-x)^r  is related to derivatives of C(x) = (1-sqrt(1-4x))/(2x)
    # Let's use the fact that the Catalan generating function is C(x) = (1-sqrt(1-4x))/(2x)
    # and sum C_{n-1} * x^n = x*C(x).

    # Actually, let's use: sum_{r>=1} C_{r-1} * z^r = z * C(z) where C(z) = (1-sqrt(1-4z))/(2z)
    # So sum C_{r-1} * z^r = (1-sqrt(1-4z))/2
    # And sum (2r+1) * C_{r-1} * z^r = (2*z*d/dz + 1) * (1-sqrt(1-4z))/2

    # For the generating function with ALTERNATING signs:
    # sum (-1)^{r+1} * C_{r-1} * x^r = -(sum C_{r-1} * (-x)^r)
    #   = -(-x) * C(-x) = x * C(-x)
    #   = x * (1 - sqrt(1+4x)) / (-2x) = -(1-sqrt(1+4x))/2 = (sqrt(1+4x)-1)/2

    # So: H(x) = sum_{r>=1} (-1)^{r+1} * C_{r-1} * x^r = (sqrt(1+4x)-1)/2

    # And: sum (-1)^{r+1} * (2r+1) * C_{r-1} * x^r = (2x*d/dx + 1) * H(x)
    #     = (2x * 4/(2*sqrt(1+4x))) + H(x)
    #     Wait, let me be careful.

    # H(x) = (sqrt(1+4x)-1)/2
    # H'(x) = 2/sqrt(1+4x)
    # (2x*d/dx + 1)*H = 2x*H'(x) + H(x)
    #   = 4x/sqrt(1+4x) + (sqrt(1+4x)-1)/2
    #   = (8x + (1+4x) - sqrt(1+4x)) / (2*sqrt(1+4x))
    #   Hmm, let me just check: this should equal g(x) = 1/2 - (1+8x)/(2*sqrt(1+4x))

    # g(x) = 1/2 - (1+8x)/(2*sqrt(1+4x))
    #       = 1/2 - 1/(2*sqrt(1+4x)) - 4x/sqrt(1+4x)
    # 2x*H' + H = 4x/sqrt(1+4x) + sqrt(1+4x)/2 - 1/2
    # Not the same. Let me just compute numerically and match.

    # Actually the known formula is: g(x) = sum_{r>=1} (-1)^{r+1} * (2r+1) * C_{r-1} * x^r
    # Let's define F(x) = sum_{r>=0} C_r * x^{r+1} = x * C(x) = (1-sqrt(1-4x))/2
    # Then F(-x) = sum_{r>=0} C_r * (-x)^{r+1} = (-x) * C(-x)
    #   = (-x)*(1-sqrt(1+4x))/(-2x) = (1-sqrt(1+4x))/2
    # So sum_{r>=1} (-1)^r * C_{r-1} * x^r = F(-x)
    #   = (1-sqrt(1+4x))/2
    # And sum_{r>=1} (-1)^{r+1} * C_{r-1} * x^r = -F(-x) = (sqrt(1+4x)-1)/2

    # Define H_0(x) = (sqrt(1+4x)-1)/2 = sum_{r>=1} (-1)^{r+1} * C_{r-1} * x^r

    # We want: g_w(x) = sum_{r>=1} (-1)^{r+1} * f_w(r) * C_{r-1} * x^r
    # where f_0(r) = 2r+1.

    # For w=0: g_0(x) = sum (-1)^{r+1} * (2r+1) * C_{r-1} * x^r
    #   = sum (-1)^{r+1} * C_{r-1} * x^r + 2 * sum (-1)^{r+1} * r * C_{r-1} * x^r
    #   = H_0(x) + 2 * x * H_0'(x)
    #   Since x*d/dx (x^r) = r*x^r, so sum r * a_r * x^r = x*d/dx(sum a_r * x^r)

    # H_0'(x) = 2/sqrt(1+4x)
    # x*H_0'(x) = 2x/sqrt(1+4x)
    # g_0(x) = H_0(x) + 2*x*H_0'(x) = (sqrt(1+4x)-1)/2 + 4x/sqrt(1+4x)
    #   = [sqrt(1+4x)*(sqrt(1+4x)-1) + 8x] / (2*sqrt(1+4x))
    #   Hmm, let me just verify numerically.

    # Verification
    import sympy
    x_sym = sympy.Symbol('x')
    H0 = (sympy.sqrt(1 + 4*x_sym) - 1) / 2
    g0_test = H0 + 2 * x_sym * sympy.diff(H0, x_sym)
    g0_test_simplified = sympy.simplify(g0_test)
    g0_known = sympy.Rational(1,2) - (1 + 8*x_sym)/(2*sympy.sqrt(1+4*x_sym))
    diff_check = sympy.simplify(g0_test - g0_known)
    print(f"\n  Verification: g_0(x) = H_0 + 2*x*H_0' ?")
    print(f"    H_0 + 2*x*H_0' = {g0_test_simplified}")
    print(f"    g_0(known)      = {g0_known}")
    print(f"    difference      = {diff_check}")

    # ===================================================================
    # Part 6: Identify f_w(r) as polynomial and compute sector gen func
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 6: POLYNOMIAL IDENTIFICATION OF f_w(r)")
    print("=" * 80)

    # Collect f_w(r) data
    fw_data = {}  # fw_data[w] = [(r, f_w(r)), ...]
    for w in range(max_w_seen + 1):
        fw_data[w] = []
        for r in sorted(data):
            val = data[r].get(w, 0.0)
            sign = (-1) ** (r + 1)
            cat_val = catalan(r - 1)
            fact_val = math.factorial(2 * r)
            norm = sign * cat_val * fact_val
            if abs(norm) > 1e-10:
                fw = val / norm
                fw_data[w].append((r, fw))

    # For each w, try to fit f_w(r) as a polynomial in r using Lagrange interpolation
    import numpy as np

    for w in range(min(max_w_seen + 1, 8)):
        if not fw_data.get(w):
            continue
        rs = np.array([r for r, _ in fw_data[w]], dtype=float)
        fws = np.array([f for _, f in fw_data[w]], dtype=float)

        if len(rs) < 2:
            continue

        # Try degrees 1, 2, 3, ... until residual is tiny
        best_deg = None
        for deg in range(1, min(8, len(rs))):
            coeffs = np.polyfit(rs, fws, deg)
            residual = np.max(np.abs(np.polyval(coeffs, rs) - fws))
            if residual < 1e-4:
                best_deg = deg
                best_coeffs = coeffs
                break

        if best_deg is not None:
            print(f"\n  w={w}: f_{w}(r) is polynomial of degree {best_deg}")
            # Convert to exact rational form
            poly_str = "  f_{}(r) = ".format(w)
            terms = []
            for i, c in enumerate(best_coeffs):
                power = best_deg - i
                c_round = round(c, 1)
                if abs(c_round) < 1e-8:
                    continue
                # Try to identify as simple fraction
                frac = Fraction(c_round).limit_denominator(1000)
                if power == 0:
                    terms.append(f"{frac}")
                elif power == 1:
                    terms.append(f"{frac}*r")
                else:
                    terms.append(f"{frac}*r^{power}")
            print(f"    {' + '.join(terms)}")
            print(f"    Coefficients (high to low): {[round(c, 6) for c in best_coeffs]}")
            print(f"    Max residual: {np.max(np.abs(np.polyval(best_coeffs, rs) - fws)):.2e}")

            # Verify at all data points
            print(f"    Verification:")
            for r_val, fw_val in fw_data[w]:
                fitted = np.polyval(best_coeffs, r_val)
                print(f"      r={int(r_val):2d}: data={fw_val:18.6f}, fit={fitted:18.6f}, "
                      f"err={abs(fw_val-fitted):.2e}")
        else:
            print(f"\n  w={w}: f_{w}(r) is NOT polynomial of degree <= 7")
            print(f"    Values: {[(int(r), round(f, 6)) for r, f in fw_data[w]]}")

    # ===================================================================
    # Part 7: EXACT rational identification via Fraction arithmetic
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 7: EXACT RATIONAL f_w(r) VALUES")
    print("=" * 80)

    for w in range(min(max_w_seen + 1, 8)):
        if not fw_data.get(w):
            continue
        print(f"\n  w={w}:")
        for r, fw in fw_data[w]:
            # Try to express as simple fraction
            frac = Fraction(fw).limit_denominator(10**8)
            # Also try expressing as ratio involving known sequences
            r_int = int(r)
            print(f"    r={r_int:2d}: f_{w} = {fw:24.10f}  ~ {frac}")

    # ===================================================================
    # Part 8: EXACT computation at high precision
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 8: HIGH-PRECISION RATIONAL RECONSTRUCTION")
    print("=" * 80)

    # Compute f_w(r) * (2r)! to get integers (if they are)
    for w in range(min(max_w_seen + 1, 8)):
        if not fw_data.get(w):
            continue
        print(f"\n  w={w}: checking if d^{w}T coeff / [(-1)^{{r+1}} * C_{{r-1}}] is integer:")
        for r in sorted(data):
            val = data[r].get(w, 0.0)
            sign = (-1) ** (int(r) + 1)
            cat_val = catalan(int(r) - 1)
            if abs(cat_val) > 0:
                reduced = val / (sign * cat_val)
                # Check if integer
                nearest_int = round(reduced)
                is_int = abs(reduced - nearest_int) < max(1.0, abs(reduced) * 1e-8)
                if is_int and abs(nearest_int) > 0:
                    print(f"    r={int(r):2d}: coeff/(sign*Cat) = {nearest_int:20d}  (integer: {is_int})")

    # ===================================================================
    # Part 9: Ratio analysis — d^w T / d^0 T
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 9: RATIO d^w T / T AT SYMMETRIC POINT")
    print("=" * 80)

    for w in range(1, min(max_w_seen + 1, 8)):
        print(f"\n  w={w}: d^{w}T / T ratio:")
        ratios_w = []
        for r in sorted(data):
            T_val = data[r].get(0, 0.0)
            dw_val = data[r].get(w, 0.0)
            if abs(T_val) > 1e-10 and abs(dw_val) > 1e-10:
                ratio = dw_val / T_val
                ratios_w.append((int(r), ratio))
                # Express as fraction
                frac = Fraction(ratio).limit_denominator(10**6)
                print(f"    r={int(r):2d}: ratio = {ratio:24.10f}  ~ {frac}")

        # Check if ratio is polynomial in r
        if len(ratios_w) >= 3:
            rs_arr = np.array([r for r, _ in ratios_w])
            rats_arr = np.array([v for _, v in ratios_w])
            for deg in range(1, 6):
                c_fit = np.polyfit(rs_arr, rats_arr, deg)
                res = np.max(np.abs(np.polyval(c_fit, rs_arr) - rats_arr))
                if res < 1e-6:
                    print(f"    => Polynomial of degree {deg}: coeffs = {[round(c, 8) for c in c_fit]}")
                    break

    # ===================================================================
    # Part 10: The master two-variable generating function
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 10: MASTER TWO-VARIABLE GENERATING FUNCTION G(x,y)")
    print("=" * 80)

    # If f_w(r) = P_w(r) (polynomial in r), then
    # G(x,y) = sum_w y^w * g_w(x)
    # where g_w(x) = sum_r (-1)^{r+1} * P_w(r) * C_{r-1} * x^r
    # Each g_w can be computed from H_0(x) and its derivatives using
    # the identity: sum_r r^j * (-1)^{r+1} * C_{r-1} * x^r = (x*d/dx)^j H_0(x)

    # Build: g_w(x) from polynomial P_w(r) using (x*d/dx) operators on H_0

    import sympy
    x_sym = sympy.Symbol('x')
    y_sym = sympy.Symbol('y')
    r_sym = sympy.Symbol('r')

    H0_sym = (sympy.sqrt(1 + 4*x_sym) - 1) / 2

    # Compute (x*d/dx)^n H_0 for n = 0, 1, 2, ...
    xD_powers = [H0_sym]
    for n in range(1, 8):
        prev = xD_powers[-1]
        xD_powers.append(sympy.simplify(x_sym * sympy.diff(prev, x_sym)))

    print("\n  (x*d/dx)^n H_0(x):")
    for n, expr in enumerate(xD_powers):
        print(f"    n={n}: {sympy.simplify(expr)}")

    # Now for each w, if f_w(r) = a_0 + a_1*r + a_2*r^2 + ...,
    # then g_w(x) = a_0*H_0 + a_1*(xD)H_0 + a_2*(xD)^2 H_0 + ...
    print("\n  Sector generating functions g_w(x):")
    G_total = sympy.S.Zero
    for w in range(min(max_w_seen + 1, 6)):
        if not fw_data.get(w) or len(fw_data[w]) < 3:
            continue
        rs_arr = np.array([r for r, _ in fw_data[w]], dtype=float)
        fws_arr = np.array([f for _, f in fw_data[w]], dtype=float)

        # Find polynomial degree
        best_deg_w = None
        best_coeffs_w = None
        for deg in range(1, min(8, len(rs_arr))):
            coeffs = np.polyfit(rs_arr, fws_arr, deg)
            residual = np.max(np.abs(np.polyval(coeffs, rs_arr) - fws_arr))
            if residual < 1e-3:
                best_deg_w = deg
                best_coeffs_w = coeffs
                break

        if best_deg_w is None:
            print(f"\n  w={w}: Cannot identify polynomial f_{w}(r)")
            continue

        # Convert numpy coeffs to exact rationals
        # Coeffs are [c_deg, c_{deg-1}, ..., c_1, c_0] (high to low power)
        poly_coeffs = {}  # power -> rational coefficient
        for i, c in enumerate(best_coeffs_w):
            power = best_deg_w - i
            frac = Fraction(c).limit_denominator(1000)
            if abs(float(frac) - c) < 1e-4:
                poly_coeffs[power] = frac
            else:
                poly_coeffs[power] = c

        # Build g_w(x) = sum_j poly_coeffs[j] * (xD)^j H_0
        gw = sympy.S.Zero
        for j, coeff in poly_coeffs.items():
            if j < len(xD_powers):
                gw += sympy.Rational(coeff.numerator, coeff.denominator) * xD_powers[j] \
                    if isinstance(coeff, Fraction) else float(coeff) * xD_powers[j]

        gw_simplified = sympy.simplify(gw)
        print(f"\n  w={w}: f_{w}(r) ~ {dict(sorted(poly_coeffs.items(), reverse=True))}")
        print(f"         g_{w}(x) = {gw_simplified}")

        # Verify: expand and check coefficients match
        gw_series = sympy.series(gw_simplified, x_sym, 0, n=min(max_r, 8) + 1)
        print(f"         Series: {gw_series}")

        G_total += gw_simplified * y_sym**w

    G_simplified = sympy.simplify(G_total)
    print(f"\n  G(x,y) = {G_simplified}")

    # Try to simplify G(x,y) as algebraic expression
    # Factor out common denominators
    G_collected = sympy.collect(sympy.expand(G_total), y_sym)
    print(f"\n  G(x,y) collected: {G_collected}")

    # ===================================================================
    # Part 11: Even arities — verify all vanish at symmetric point
    # ===================================================================
    print("\n" + "=" * 80)
    print("PART 11: EVEN ARITIES AT SYMMETRIC POINT")
    print("=" * 80)

    even_data = compute_even_field_sectors(max_r=8, c_val=0.0)
    for r in sorted(even_data):
        k = 2 * r
        sector = even_data[r]
        total = sum(abs(v) for v in sector.values())
        print(f"  k={k:2d}: total |m_k(sym)| = {total:.6e}, vanishes = {total < 1e-6}")

    print("\n" + "=" * 80)
    print("COMPUTATION COMPLETE")
    print("=" * 80)


def _try_polynomial_fit(w, fw_seq):
    """Try to fit f_w(r) as exact polynomial using Lagrange interpolation with rationals."""
    if len(fw_seq) < 2:
        return

    try:
        import numpy as np
    except ImportError:
        return

    rs = np.array([r for r, _ in fw_seq], dtype=float)
    fws = np.array([f for _, f in fw_seq], dtype=float)

    # Try degrees 1 through 6
    for deg in range(1, min(7, len(rs))):
        coeffs = np.polyfit(rs, fws, deg)
        residual = np.max(np.abs(np.polyval(coeffs, rs) - fws))
        if residual < 1e-4:
            print(f"    Polynomial degree {deg} fit: residual = {residual:.2e}")
            print(f"    Coefficients: {[round(c, 6) for c in coeffs]}")
            return True
    return False


if __name__ == '__main__':
    main()
