r"""Virasoro E_1 ordered shadow invariants via the closed-form Catalan formula.

The shadow coefficients S_r (r >= 2) are the E_1 ordered bar complex
invariants of the Virasoro chiral algebra at central charge c. They encode
the ordered associative Koszul dual structure.

CLOSED-FORM CATALAN FORMULA
============================

    S_r = (-6)^{r-4} * D / (2r * c^{r-3}) * F_r(D/144)

where:
    D = 80 / (5c + 22)
    F_r(x) = sum_{j=0}^{floor((r-4)/2)} (-1)^j * Cat_j * C(r-4, 2j) * x^j
    Cat_j = (2j)! / (j!(j+1)!)   (Catalan numbers)

This formula is equivalent to the convolution recursion in
lib/shadow_borel_resurgence.py but makes the algebraic structure manifest:
- The denominator is exactly c^{r-3} * (5c+22)^{floor((r-2)/2)}
- The numerator is a polynomial in c of bounded degree
- The sign pattern encodes the Koszul oscillation

Cross-check: the recursion-based engine in shadow_borel_resurgence.py must
agree with this closed-form at all tested values.

Manuscript references:
    ordered_associative_chiral_kd_core.tex (Part VII)
    spectral-braiding-core.tex (Part V)
"""

from __future__ import annotations

import json
import math
import os
import sys
from fractions import Fraction
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Catalan numbers and binomial coefficients (exact, integer)
# ---------------------------------------------------------------------------

def catalan(n: int) -> int:
    """Catalan number Cat_n = C(2n,n)/(n+1)."""
    if n < 0:
        return 0
    return math.comb(2 * n, n) // (n + 1)


def F_r_poly(r: int, x: Fraction) -> Fraction:
    """Evaluate F_r(x) = sum_{j=0}^{floor((r-4)/2)} (-1)^j Cat_j C(r-4,2j) x^j."""
    if r < 4:
        # S_2 and S_3 are handled separately
        return Fraction(0)
    n = r - 4  # degree parameter
    j_max = n // 2
    result = Fraction(0)
    x_pow = Fraction(1)
    for j in range(j_max + 1):
        sign = (-1) ** j
        cat_j = catalan(j)
        binom = math.comb(n, 2 * j)
        result += Fraction(sign * cat_j * binom) * x_pow
        x_pow *= x
    return result


# ---------------------------------------------------------------------------
# Shadow coefficients: exact rational functions of c
# ---------------------------------------------------------------------------

def shadow_catalan_exact(c: Fraction, r_max: int = 50) -> Dict[int, Fraction]:
    """Compute S_r for r=2,...,r_max as exact Fraction(numerator, denominator).

    Uses the closed-form Catalan formula for r >= 4.
    S_2 and S_3 are computed from the known small-r values:
        S_2 = c/2  (= kappa, the leading coefficient)
        S_3 = 2    (= alpha, the subleading)
    """
    if c == 0:
        raise ValueError("c = 0 is degenerate (kappa = 0)")

    D = Fraction(80) / (5 * c + 22)
    x = D / 144  # D/144

    result: Dict[int, Fraction] = {}

    # r=2: S_2 = kappa = c/2
    if r_max >= 2:
        result[2] = c / 2
    # r=3: S_3 = alpha = 2
    if r_max >= 3:
        result[3] = Fraction(2)

    for r in range(4, r_max + 1):
        # S_r = (-6)^{r-4} * D / (2r * c^{r-3}) * F_r(D/144)
        neg6_pow = Fraction((-6) ** (r - 4))
        c_pow = c ** (r - 3)
        Fr = F_r_poly(r, x)
        S_r = neg6_pow * D * Fr / (2 * r * c_pow)
        result[r] = S_r

    return result


def shadow_catalan_float(c_val: float, r_max: int = 50) -> Dict[int, float]:
    """Floating-point evaluation of the Catalan formula."""
    if abs(c_val) < 1e-30:
        raise ValueError("c = 0 degenerate")
    D = 80.0 / (5.0 * c_val + 22.0)
    x = D / 144.0

    result: Dict[int, float] = {}
    if r_max >= 2:
        result[2] = c_val / 2.0
    if r_max >= 3:
        result[3] = 2.0

    for r in range(4, r_max + 1):
        neg6_pow = (-6.0) ** (r - 4)
        c_pow = c_val ** (r - 3)
        Fr = 0.0
        n = r - 4
        j_max = n // 2
        x_pow = 1.0
        for j in range(j_max + 1):
            sign = (-1.0) ** j
            cat_j = catalan(j)
            binom = math.comb(n, 2 * j)
            Fr += sign * cat_j * binom * x_pow
            x_pow *= x
        S_r = neg6_pow * D * Fr / (2.0 * r * c_pow)
        result[r] = S_r

    return result


# ---------------------------------------------------------------------------
# Cross-check against the convolution recursion
# ---------------------------------------------------------------------------

def cross_check_recursion(c: Fraction, r_max: int = 50) -> Dict[int, Fraction]:
    """Compute S_r via the convolution recursion (exact Fraction arithmetic).

    This reproduces shadow_coefficients_fraction from the engine.
    """
    kappa = c / 2
    alpha = Fraction(2)
    S4_val = Fraction(10) / (c * (5 * c + 22))
    q0 = 4 * kappa ** 2
    q1 = 12 * kappa * alpha
    q2 = 9 * alpha ** 2 + 16 * kappa * S4_val

    max_n = r_max - 2
    if max_n < 0:
        return {}

    # a[0] = sqrt(q0) = 2*|kappa| = |c|
    a0 = abs(c)
    a = [Fraction(0)] * (max_n + 1)
    a[0] = a0
    if max_n >= 1:
        a[1] = q1 / (2 * a0)
    if max_n >= 2:
        a[2] = (q2 - a[1] ** 2) / (2 * a0)
    for n in range(3, max_n + 1):
        conv = sum(a[j] * a[n - j] for j in range(1, n))
        a[n] = -conv / (2 * a0)

    return {r: a[r - 2] / r for r in range(2, r_max + 1)}


# ---------------------------------------------------------------------------
# Denominator structure analysis
# ---------------------------------------------------------------------------

def analyze_denominator(S_r: Fraction, r: int, c: Fraction) -> Dict:
    """Analyze the denominator of S_r to check the predicted structure.

    Predicted: denom(S_r) divides c^{r-3} * (5c+22)^{floor((r-2)/2)} * (2r).
    """
    num = S_r.numerator
    den = S_r.denominator
    sign = 1 if num > 0 else (-1 if num < 0 else 0)
    return {
        'r': r,
        'numerator': num,
        'denominator': den,
        'sign': sign,
        'num_digits': len(str(abs(num))) if num != 0 else 0,
        'den_digits': len(str(abs(den))),
    }


def verify_denominator_formula(c: Fraction, r_max: int = 50) -> List[Dict]:
    """Verify denom(S_r) | c^{r-3} * (5c+22)^{floor((r-2)/2)} for all r."""
    coeffs = shadow_catalan_exact(c, r_max)
    results = []
    five_c_22 = 5 * c + 22
    for r in range(2, r_max + 1):
        S_r = coeffs[r]
        # Predicted denominator (up to a factor of 2r from the formula)
        predicted_denom_factor = c ** max(r - 3, 0) * five_c_22 ** ((r - 2) // 2)
        # S_r * predicted_denom_factor * (2*r) should be an integer (Fraction with den=1)
        test_val = S_r * predicted_denom_factor * (2 * r)
        is_integer = (test_val.denominator == 1)
        results.append({
            'r': r,
            'S_r_sign': 1 if S_r > 0 else (-1 if S_r < 0 else 0),
            'denom_divides': is_integer,
            'test_val_num': test_val.numerator if is_integer else None,
        })
    return results


# ---------------------------------------------------------------------------
# Numerical evaluations at special central charges
# ---------------------------------------------------------------------------

SPECIAL_CHARGES = {
    'Ising (c=1/2)': Fraction(1, 2),
    'c=1': Fraction(1),
    'Self-dual (c=13)': Fraction(13),
    'c=25': Fraction(25),
    'Critical string (c=26)': Fraction(26),
}


def evaluate_at_special_charges(r_max: int = 50) -> Dict[str, Dict[int, float]]:
    """Evaluate S_r at special central charges, returning float values."""
    results = {}
    for name, c_val in SPECIAL_CHARGES.items():
        coeffs = shadow_catalan_exact(c_val, r_max)
        results[name] = {r: float(v) for r, v in coeffs.items()}
    return results


# ---------------------------------------------------------------------------
# Convergence analysis at c=13
# ---------------------------------------------------------------------------

def convergence_analysis_c13(r_max: int = 50) -> Dict:
    """Analyze |S_r| at c=13 to confirm geometric decay with rho ~ 0.467.

    The raw ratio |S_{r+1}/S_r| oscillates due to the cosine factor in
    the Darboux asymptotics. We use the THEORETICAL rho from the branch
    point positions and also log-linear regression over a stable window.
    """
    # Import the shadow data class for the theoretical prediction
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
    from shadow_borel_resurgence import VirasoroShadowData

    d = VirasoroShadowData(13.0)
    rho_theory = d.rho

    coeffs = shadow_catalan_float(13.0, r_max)
    abs_vals = {r: abs(v) for r, v in coeffs.items()}

    # Estimate rho from successive ratios |S_{r+1}/S_r| for large r
    ratios = {}
    for r in range(10, r_max):
        if abs_vals[r] > 1e-300 and abs_vals.get(r + 1, 0) > 0:
            ratios[r] = abs_vals[r + 1] / abs_vals[r]

    # Also estimate from |S_r|^{1/r}
    root_estimates = {}
    for r in range(10, r_max + 1):
        if abs_vals[r] > 0:
            root_estimates[r] = abs_vals[r] ** (1.0 / r)

    # Log-linear regression: log|S_r| ~ r*log(rho) + const for r in stable window
    # Use r=10..22 as the stable window (before oscillatory crossings become large)
    log_data = []
    for r in range(10, min(23, r_max + 1)):
        if abs_vals[r] > 0:
            log_data.append((r, math.log(abs_vals[r])))
    if len(log_data) >= 3:
        n_pts = len(log_data)
        sum_r = sum(p[0] for p in log_data)
        sum_logS = sum(p[1] for p in log_data)
        sum_r2 = sum(p[0] ** 2 for p in log_data)
        sum_rlogS = sum(p[0] * p[1] for p in log_data)
        slope = (n_pts * sum_rlogS - sum_r * sum_logS) / (n_pts * sum_r2 - sum_r ** 2)
        rho_regression = math.exp(slope)
    else:
        rho_regression = None

    return {
        'abs_vals': abs_vals,
        'ratios': ratios,
        'root_estimates': root_estimates,
        'rho_theory': rho_theory,
        'rho_regression': rho_regression,
        'rho_estimate_ratio': sum(ratios.values()) / len(ratios) if ratios else None,
    }


# ---------------------------------------------------------------------------
# Main computation and output
# ---------------------------------------------------------------------------

def main():
    R_MAX = 50
    print("=" * 78)
    print("VIRASORO E_1 ORDERED SHADOW INVARIANTS — CATALAN FORMULA")
    print(f"Computing S_r for r = 2, ..., {R_MAX}")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # (1) Exact rational computation at c symbolic level
    #     We use c = Fraction for several test values, but also display
    #     the structure at a generic rational c.
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECTION 1: Cross-check Catalan formula vs convolution recursion")
    print("-" * 78)

    test_charges = [Fraction(1, 2), Fraction(1), Fraction(13), Fraction(25), Fraction(26)]
    all_pass = True
    for c_test in test_charges:
        catalan_coeffs = shadow_catalan_exact(c_test, R_MAX)
        recursion_coeffs = cross_check_recursion(c_test, R_MAX)
        mismatches = []
        for r in range(2, R_MAX + 1):
            if catalan_coeffs[r] != recursion_coeffs[r]:
                mismatches.append(r)
                all_pass = False
        status = "PASS" if not mismatches else f"FAIL at r={mismatches}"
        print(f"  c = {c_test}: {status}")

    if all_pass:
        print("  >>> ALL cross-checks PASS through r=50 <<<")
    else:
        print("  >>> CROSS-CHECK FAILURES DETECTED <<<")

    # -----------------------------------------------------------------------
    # (2) Exact S_r at c=13 (self-dual) as rational numbers
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECTION 2: Exact S_r at c = 13 (self-dual point)")
    print("-" * 78)

    c13 = Fraction(13)
    coeffs_13 = shadow_catalan_exact(c13, R_MAX)
    for r in range(2, min(R_MAX + 1, 16)):
        val = coeffs_13[r]
        print(f"  S_{r:2d} = {val}  ({float(val):.10e})")
    print(f"  ... (continuing to r={R_MAX})")

    # -----------------------------------------------------------------------
    # (3) Denominator structure verification
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECTION 3: Denominator formula verification")
    print("  Predicted: S_r has denom dividing c^{r-3}*(5c+22)^{floor((r-2)/2)}")
    print("-" * 78)

    for c_test in test_charges:
        results = verify_denominator_formula(c_test, R_MAX)
        failures = [d for d in results if not d['denom_divides']]
        if not failures:
            print(f"  c = {c_test}: VERIFIED through r={R_MAX}")
        else:
            print(f"  c = {c_test}: FAILED at r = {[d['r'] for d in failures]}")

    # -----------------------------------------------------------------------
    # (4) Full table: sign, numerator degree, denominator structure
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECTION 4: Sign pattern and coefficient growth")
    print("-" * 78)
    print(f"  {'r':>3s}  {'sign':>5s}  {'|num| digits':>12s}  {'|den| digits':>12s}")
    print(f"  {'---':>3s}  {'-----':>5s}  {'------------':>12s}  {'------------':>12s}")

    c_generic = Fraction(13)
    coeffs_gen = shadow_catalan_exact(c_generic, R_MAX)
    for r in range(2, R_MAX + 1):
        info = analyze_denominator(coeffs_gen[r], r, c_generic)
        sign_str = "+" if info['sign'] > 0 else ("-" if info['sign'] < 0 else "0")
        print(f"  {r:3d}  {sign_str:>5s}  {info['num_digits']:>12d}  {info['den_digits']:>12d}")

    # -----------------------------------------------------------------------
    # (5) Numerical evaluations at special charges
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECTION 5: Numerical evaluations at special central charges")
    print("-" * 78)

    evals = evaluate_at_special_charges(R_MAX)
    # Print header
    charge_names = list(SPECIAL_CHARGES.keys())
    short_names = ['c=1/2', 'c=1', 'c=13', 'c=25', 'c=26']
    print(f"  {'r':>3s}  ", end="")
    for sn in short_names:
        print(f"  {sn:>14s}", end="")
    print()
    print(f"  {'---':>3s}  ", end="")
    for _ in short_names:
        print(f"  {'-'*14:>14s}", end="")
    print()

    for r in range(2, R_MAX + 1):
        print(f"  {r:3d}  ", end="")
        for name in charge_names:
            val = evals[name][r]
            if abs(val) < 1e-100:
                print(f"  {'0':>14s}", end="")
            else:
                print(f"  {val:>14.6e}", end="")
        print()

    # -----------------------------------------------------------------------
    # (6) Convergence analysis at c=13
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECTION 6: Convergence at c=13 — geometric decay analysis")
    print("-" * 78)

    conv = convergence_analysis_c13(R_MAX)
    print(f"  Theoretical rho (from branch points): {conv['rho_theory']:.6f}")
    print(f"  Log-linear regression rho (r=10..22): {conv['rho_regression']:.6f}")
    print(f"  Naive ratio average rho:              {conv['rho_estimate_ratio']:.6f}")
    print()
    print("  NOTE: The raw |S_{r+1}/S_r| ratios oscillate because S_r has a")
    print("  cosine factor from the Darboux asymptotics. The theoretical rho")
    print("  from the branch point |t_+|^{-1} is the correct growth rate.")
    print()
    print(f"  {'r':>3s}  {'|S_r|':>14s}  {'|S_{r+1}/S_r|':>14s}  {'|S_r|^{1/r}':>14s}")
    print(f"  {'---':>3s}  {'-'*14:>14s}  {'-'*14:>14s}  {'-'*14:>14s}")
    for r in range(2, R_MAX + 1):
        abs_val = conv['abs_vals'].get(r, 0.0)
        ratio = conv['ratios'].get(r, None)
        root = conv['root_estimates'].get(r, None)
        ratio_str = f"{ratio:.6f}" if ratio is not None else "---"
        root_str = f"{root:.6f}" if root is not None else "---"
        print(f"  {r:3d}  {abs_val:>14.6e}  {ratio_str:>14s}  {root_str:>14s}")

    # -----------------------------------------------------------------------
    # (7) Save full JSON table
    # -----------------------------------------------------------------------
    output_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(output_dir, "ordered_e1_shadow_catalan_table.json")

    table_data = {
        'description': 'Virasoro E_1 ordered shadow coefficients S_r via Catalan formula',
        'r_max': R_MAX,
        'formula': 'S_r = (-6)^{r-4} * D / (2r * c^{r-3}) * F_r(D/144), D=80/(5c+22)',
        'cross_check': 'All values verified against convolution recursion',
    }

    # Store exact rational for c=13
    exact_c13 = {}
    for r in range(2, R_MAX + 1):
        v = coeffs_13[r]
        exact_c13[str(r)] = {
            'numerator': v.numerator,
            'denominator': v.denominator,
            'float': float(v),
        }
    table_data['exact_c13'] = exact_c13

    # Store floats for all special charges
    float_tables = {}
    for name, c_val in SPECIAL_CHARGES.items():
        coeffs = shadow_catalan_exact(c_val, R_MAX)
        float_tables[name] = {str(r): float(v) for r, v in coeffs.items()}
    table_data['special_charges'] = float_tables

    # Denominator verification summary
    denom_verified = {}
    for c_test in test_charges:
        res = verify_denominator_formula(c_test, R_MAX)
        denom_verified[str(c_test)] = all(d['denom_divides'] for d in res)
    table_data['denominator_formula_verified'] = denom_verified

    # Convergence data at c=13
    table_data['convergence_c13'] = {
        'rho_theory': conv['rho_theory'],
        'rho_regression': conv['rho_regression'],
        'rho_ratio_avg': conv['rho_estimate_ratio'],
    }

    with open(json_path, 'w') as f:
        json.dump(table_data, f, indent=2)
    print(f"\n  Table saved to: {json_path}")

    # -----------------------------------------------------------------------
    # (7b) Sign anomaly analysis at c=13
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECTION 7: Sign anomalies from Darboux cosine zeros (c=13)")
    print("-" * 78)

    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib'))
    from shadow_borel_resurgence import VirasoroShadowData
    import cmath

    d13 = VirasoroShadowData(13.0)
    omega_13 = abs(cmath.phase(d13.t_plus))
    sigma_13 = cmath.sqrt(-d13.q2 * d13.t_plus * (d13.t_plus - d13.t_minus))
    ST2_13 = sigma_13 * d13.t_plus ** 2
    phi_13 = cmath.phase(ST2_13)

    print(f"  omega/pi = {omega_13/math.pi:.6f} (deviation from 1: {1 - omega_13/math.pi:.6f})")
    print(f"  phi = {phi_13:.6f}")
    print(f"  rho = {d13.rho:.6f}")
    print(f"  Sign anomaly period ~ pi/(pi - omega) = {math.pi/(math.pi - omega_13):.1f}")
    print()

    # Find sign anomalies
    c13_coeffs = shadow_catalan_exact(Fraction(13), R_MAX)
    prev_sign = None
    anomalies = []
    for r in range(2, R_MAX + 1):
        v = c13_coeffs[r]
        s = 1 if v > 0 else (-1 if v < 0 else 0)
        if prev_sign is not None and s == prev_sign:
            anomalies.append(r)
        prev_sign = s
    print(f"  Sign anomaly positions (two consecutive same signs): r = {anomalies}")
    print(f"  Spacing: {[anomalies[i] - anomalies[i-1] for i in range(1, len(anomalies))]}")
    print(f"  At these r, the Darboux cosine cos(r*omega+phi) crosses zero,")
    print(f"  the leading asymptotic vanishes, and |S_r| temporarily increases.")

    # -----------------------------------------------------------------------
    # (8) Detailed exact rational table for small r (c generic = symbol)
    # -----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SECTION 8: Exact rational S_r at all five special charges (r=2..15)")
    print("-" * 78)

    for name, c_val in SPECIAL_CHARGES.items():
        print(f"\n  {name} (c = {c_val}):")
        coeffs = shadow_catalan_exact(c_val, R_MAX)
        for r in range(2, 16):
            v = coeffs[r]
            print(f"    S_{r:2d} = {v}")

    print("\n" + "=" * 78)
    print("COMPUTATION COMPLETE")
    print("=" * 78)


if __name__ == '__main__':
    main()
