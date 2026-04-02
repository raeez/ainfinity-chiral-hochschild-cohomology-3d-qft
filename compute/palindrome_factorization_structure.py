r"""Virasoro palindrome / factorization structure of m_k|_T for k=2,...,20.

Comprehensive numerical investigation at c=1:

(1) Even k: does m_k|_T vanish at the symmetric point lambda_i = 1? (Predicted: YES)
(2) Even k: what is the minimum depth? (Predicted: 2)
(3) Even k: does m_k|_T have any linear factors (lambda_i - lambda_j) or (lambda_i + lambda_j)?
(4) Odd k: verify T_k(1,...,1) = (-1)^n * C_n * k! with n = (k-3)/2
(5) Leading quadratic form at depth 2 for even k: compute the matrix M_k (tridiagonal, anti-palindromic)
(6) Null vectors of M_k: are they palindromic? Is there a pattern?

Additionally: test factorization for k=4 (known: triple factorization),
k=6, k=8 (known: irreducible), k=10, 12, 14, 16, 18, 20.
"""

from __future__ import annotations
import sys
import os
import math
import time
import random
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from compute.m7_m10_depth_frontier import StasheffEngine


# ====================================================================
# QUESTION 1: Even-arity symmetric point vanishing
# ====================================================================

def test_symmetric_point_vanishing(max_k=20, c_val=1.0):
    """Test whether m_k|_T vanishes at lambda_i = 1 for all i, for even k."""
    print("=" * 100)
    print(f"QUESTION 1: Even-arity symmetric point vanishing (c={c_val})")
    print("=" * 100)
    engine = StasheffEngine(c_val)

    print(f"\n  {'k':>3} {'parity':>6} {'|T-sector|':>15} {'T-coeff':>15} "
          f"{'signed_sum':>15} {'scalar':>15} {'vanishes?':>10}")
    print("  " + "-" * 87)

    for k in range(2, max_k + 1):
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(k - 1))
        t0 = time.time()
        result = engine.mk(lams)
        elapsed = time.time() - t0

        fields = {f: v for f, v in result.items() if f >= 0}
        T_coeff = result.get(0, 0.0)
        signed_sum = sum(v for v in fields.values())
        abs_sum = sum(abs(v) for v in fields.values())
        scalar = result.get(-1, 0.0)
        parity = "even" if k % 2 == 0 else "odd"
        vanishes = "YES" if abs_sum < 1e-6 else "NO"

        print(f"  {k:>3} {parity:>6} {abs_sum:>15.6e} {T_coeff:>15.4f} "
              f"{signed_sum:>15.4f} {scalar:>15.6e} {vanishes:>10}"
              + (f"  [{elapsed:.1f}s]" if elapsed > 1.0 else ""))


# ====================================================================
# QUESTION 2: Minimum depth for even k
# ====================================================================

def test_minimum_depth(max_k=20, c_val=1.0, n_samples=100, seed=99):
    """For each even k, find the minimum depth at which m_k|_T is nonzero.

    Depth of d^w T is (k-1-w). Minimum depth = max populated derivative order + 1
    subtracted from k-1.
    """
    print("\n" + "=" * 100)
    print(f"QUESTION 2: Minimum depth for even k (c={c_val})")
    print("=" * 100)

    engine = StasheffEngine(c_val)
    rng = random.Random(seed)

    print(f"\n  {'k':>3} {'min_depth':>10} {'max_deriv_order':>16} "
          f"{'populated_fields':>50}")
    print("  " + "-" * 85)

    for k in range(4, max_k + 1, 2):
        populated = set()
        for _ in range(n_samples):
            engine._cache.clear()
            lams = tuple(rng.uniform(0.1, 3.0) for _ in range(k - 1))
            result = engine.mk(lams)
            for f, v in result.items():
                if f >= 0 and abs(v) > 1e-10:
                    populated.add(f)

        if populated:
            max_w = max(populated)
            min_depth = k - 1 - max_w
        else:
            max_w = -1
            min_depth = k + 1  # only scalar

        pop_str = ", ".join(f"d^{w}T" if w > 0 else "T" for w in sorted(populated))
        if len(pop_str) > 48:
            pop_str = pop_str[:45] + "..."
        print(f"  {k:>3} {min_depth:>10} {max_w:>16} {pop_str:>50}")


# ====================================================================
# QUESTION 3: Linear factor test for even k
# ====================================================================

def test_linear_factors(max_k=16, c_val=1.0, n_samples=200, seed=42):
    """For each even k, test whether m_k|_T has linear factors
    (lambda_i - lambda_j) or (lambda_i + lambda_j).

    Method: set lambda_j = lambda_i (for minus) or lambda_j = -lambda_i (for plus)
    and check if ALL T-sector coefficients vanish.
    """
    print("\n" + "=" * 100)
    print(f"QUESTION 3: Linear factor test for even k (c={c_val})")
    print("=" * 100)

    engine = StasheffEngine(c_val)
    rng = random.Random(seed)

    for k in range(4, max_k + 1, 2):
        n_lams = k - 1
        print(f"\n  k={k}: testing (lambda_i - lambda_j) and (lambda_i + lambda_j) factors")

        # Test (lambda_i - lambda_j) for all pairs
        minus_factors = []
        for i in range(n_lams):
            for j in range(i + 1, n_lams):
                max_T = 0.0
                max_T_gen = 0.0
                for _ in range(n_samples):
                    engine._cache.clear()
                    base_lams = [rng.uniform(0.2, 3.0) for _ in range(n_lams)]
                    # Generic
                    result_gen = engine.mk(tuple(base_lams))
                    T_gen = sum(abs(v) for f, v in result_gen.items() if f >= 0)
                    max_T_gen = max(max_T_gen, T_gen)

                    # Constrained: lambda_j = lambda_i
                    constrained = list(base_lams)
                    constrained[j] = constrained[i]
                    engine._cache.clear()
                    result_con = engine.mk(tuple(constrained))
                    T_con = sum(abs(v) for f, v in result_con.items() if f >= 0)
                    max_T = max(max_T, T_con)

                ratio = max_T / max_T_gen if max_T_gen > 1e-10 else 0.0
                if ratio < 1e-8:
                    minus_factors.append((i + 1, j + 1))

        # Test (lambda_i + lambda_j) for all pairs
        plus_factors = []
        for i in range(n_lams):
            for j in range(i + 1, n_lams):
                max_T = 0.0
                max_T_gen = 0.0
                for _ in range(n_samples):
                    engine._cache.clear()
                    base_lams = [rng.uniform(0.2, 3.0) for _ in range(n_lams)]
                    # Generic
                    result_gen = engine.mk(tuple(base_lams))
                    T_gen = sum(abs(v) for f, v in result_gen.items() if f >= 0)
                    max_T_gen = max(max_T_gen, T_gen)

                    # Constrained: lambda_j = -lambda_i
                    constrained = list(base_lams)
                    constrained[j] = -constrained[i]
                    engine._cache.clear()
                    result_con = engine.mk(tuple(constrained))
                    T_con = sum(abs(v) for f, v in result_con.items() if f >= 0)
                    max_T = max(max_T, T_con)

                ratio = max_T / max_T_gen if max_T_gen > 1e-10 else 0.0
                if ratio < 1e-8:
                    plus_factors.append((i + 1, j + 1))

        if minus_factors:
            print(f"    (lambda_i - lambda_j) factors: {minus_factors}")
        else:
            print(f"    NO (lambda_i - lambda_j) factors found")

        if plus_factors:
            print(f"    (lambda_i + lambda_j) factors: {plus_factors}")
        else:
            print(f"    NO (lambda_i + lambda_j) factors found")

        # Also test palindrome factor (lambda_i - lambda_{k-i})
        print(f"    Palindrome pairs (i, k-i):")
        for i in range(n_lams):
            j = n_lams - 1 - i
            if j <= i:
                break
            is_minus = (i + 1, j + 1) in minus_factors
            is_plus = (i + 1, j + 1) in plus_factors
            print(f"      (lambda_{i+1} - lambda_{j+1}): {'FACTOR' if is_minus else 'no'}"
                  f"    (lambda_{i+1} + lambda_{j+1}): {'FACTOR' if is_plus else 'no'}")


# ====================================================================
# QUESTION 4: Catalan verification for odd k
# ====================================================================

def verify_catalan_pattern(max_k=19, c_val=1.0):
    """Verify T_k(1,...,1) = (-1)^n * C_n * k! for odd k, with n = (k-3)/2.

    Also verify S_k(1,...,1) = (-1)^n * C_n * (k+1)!/2.
    """
    print("\n" + "=" * 100)
    print(f"QUESTION 4: Catalan pattern for odd k (c={c_val})")
    print("=" * 100)

    engine = StasheffEngine(c_val)

    print(f"\n  {'k':>3} {'n':>3} {'C_n':>8} {'T_predicted':>20} {'T_actual':>20} "
          f"{'match':>6} {'S_predicted':>20} {'S_actual':>20} {'match':>6}")
    print("  " + "-" * 112)

    for k in range(3, max_k + 1, 2):
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(k - 1))
        t0 = time.time()
        result = engine.mk(lams)
        elapsed = time.time() - t0

        fields = {f: v for f, v in result.items() if f >= 0}
        T_val = result.get(0, 0.0)
        S_val = sum(v for v in fields.values())

        n = (k - 3) // 2
        C_n = math.comb(2 * n, n) // (n + 1) if n >= 0 else 1
        T_pred = (-1)**n * C_n * math.factorial(k)
        S_pred = (-1)**n * C_n * math.factorial(k + 1) // 2

        T_match = "YES" if abs(T_val - T_pred) < max(1, abs(T_pred) * 1e-8) else "NO"
        S_match = "YES" if abs(S_val - S_pred) < max(1, abs(S_pred) * 1e-8) else "NO"

        suffix = f"  [{elapsed:.1f}s]" if elapsed > 1.0 else ""
        print(f"  {k:>3} {n:>3} {C_n:>8} {T_pred:>20} {T_val:>20.1f} {T_match:>6}"
              f" {S_pred:>20} {S_val:>20.1f} {S_match:>6}{suffix}")


# ====================================================================
# QUESTION 5: Leading quadratic form at depth 2 for even k
# ====================================================================

def compute_depth2_quadratic_form(k, c_val=1.0, epsilon=1e-6):
    """Compute the matrix M_k of the leading quadratic form at depth 2 for even k.

    At depth 2, the field is d^{k-3}T. The coefficient of d^{k-3}T in m_k|_T
    is a polynomial of degree 2 in (lambda_1, ..., lambda_{k-1}). We extract
    the quadratic form by numerical differentiation.

    M_k[i,j] = d^2 / (d lambda_i d lambda_j) [coeff of d^{k-3}T at lambda=0]
    """
    engine = StasheffEngine(c_val)
    n = k - 1  # number of spectral parameters
    target_field = k - 3  # derivative order for depth 2

    # First check that d^{k-3}T is the shallowest field (highest deriv order populated)
    # at generic lambda
    rng = random.Random(777)
    test_lams = tuple(rng.uniform(0.5, 2.0) for _ in range(n))
    engine._cache.clear()
    test_result = engine.mk(test_lams)
    max_populated = max((f for f, v in test_result.items() if f >= 0 and abs(v) > 1e-10),
                        default=-1)

    # The quadratic form matrix
    M = [[0.0] * n for _ in range(n)]

    # Use central differences around lambda = 0
    base = [0.0] * n

    def eval_field(lams_list):
        engine._cache.clear()
        result = engine.mk(tuple(lams_list))
        return result.get(target_field, 0.0)

    # Diagonal: M[i,i] = (f(+e_i) + f(-e_i) - 2*f(0)) / eps^2
    f0 = eval_field(base)

    for i in range(n):
        lp = list(base)
        lp[i] = epsilon
        lm = list(base)
        lm[i] = -epsilon
        fp = eval_field(lp)
        fm = eval_field(lm)
        M[i][i] = (fp + fm - 2 * f0) / (epsilon ** 2)

    # Off-diagonal: M[i,j] = (f(+e_i+e_j) - f(+e_i-e_j) - f(-e_i+e_j) + f(-e_i-e_j)) / (4*eps^2)
    for i in range(n):
        for j in range(i + 1, n):
            lpp = list(base); lpp[i] = epsilon; lpp[j] = epsilon
            lpm = list(base); lpm[i] = epsilon; lpm[j] = -epsilon
            lmp = list(base); lmp[i] = -epsilon; lmp[j] = epsilon
            lmm = list(base); lmm[i] = -epsilon; lmm[j] = -epsilon
            fpp = eval_field(lpp)
            fpm = eval_field(lpm)
            fmp = eval_field(lmp)
            fmm = eval_field(lmm)
            M[i][j] = (fpp - fpm - fmp + fmm) / (4 * epsilon ** 2)
            M[j][i] = M[i][j]

    return M, max_populated, target_field


def analyze_quadratic_forms(even_ks=None, c_val=1.0):
    """Compute and analyze M_k for several even k values."""
    if even_ks is None:
        even_ks = [4, 6, 8, 10, 12]
    print("\n" + "=" * 100)
    print(f"QUESTION 5: Depth-2 quadratic form M_k for even k (c={c_val})")
    print("=" * 100)

    for k in even_ks:
        print(f"\n  k={k}: target field d^{k-3}T (depth 2)")
        t0 = time.time()
        M, max_pop, target = compute_depth2_quadratic_form(k, c_val)
        elapsed = time.time() - t0
        n = k - 1

        print(f"    Max populated field: d^{max_pop}T, target: d^{target}T")
        print(f"    Computation time: {elapsed:.1f}s")

        # Print the matrix
        print(f"    M_{k} ({n}x{n}):")
        for i in range(n):
            row_str = "    "
            for j in range(n):
                val = M[i][j]
                if abs(val) < 1e-6:
                    row_str += f"{'0':>10}"
                else:
                    # Round to nearest integer if close
                    ival = round(val)
                    if abs(val - ival) < max(1e-4, abs(val) * 1e-6):
                        row_str += f"{ival:>10}"
                    else:
                        row_str += f"{val:>10.3f}"
            print(row_str)

        # Check if tridiagonal
        is_tridiagonal = True
        for i in range(n):
            for j in range(n):
                if abs(i - j) > 1 and abs(M[i][j]) > 1e-4:
                    is_tridiagonal = False
                    break

        print(f"    Tridiagonal: {'YES' if is_tridiagonal else 'NO'}")

        # Check anti-palindrome: M[i,j] = -M[n-1-i, n-1-j]?
        is_anti_palindrome = True
        for i in range(n):
            for j in range(n):
                if abs(M[i][j] + M[n-1-i][n-1-j]) > 1e-4 * max(1, abs(M[i][j])):
                    is_anti_palindrome = False
                    break

        print(f"    Anti-palindrome (M[i,j] = -M[n-1-i,n-1-j]): {'YES' if is_anti_palindrome else 'NO'}")

        # Check palindrome: M[i,j] = M[n-1-i, n-1-j]?
        is_palindrome = True
        for i in range(n):
            for j in range(n):
                if abs(M[i][j] - M[n-1-i][n-1-j]) > 1e-4 * max(1, abs(M[i][j])):
                    is_palindrome = False
                    break

        print(f"    Palindrome (M[i,j] = M[n-1-i,n-1-j]): {'YES' if is_palindrome else 'NO'}")

        # Diagonal and sub/super-diagonal values
        diag = [round(M[i][i]) for i in range(n)]
        super_diag = [round(M[i][i+1]) for i in range(n-1)] if n > 1 else []
        sub_diag = [round(M[i+1][i]) for i in range(n-1)] if n > 1 else []
        print(f"    Diagonal: {diag}")
        print(f"    Super-diagonal: {super_diag}")
        print(f"    Sub-diagonal: {sub_diag}")

        # Eigenvalues (via numpy if available, else manual for small matrices)
        try:
            import numpy as np
            evals = np.linalg.eigvalsh(np.array(M, dtype=float))
            evals_sorted = sorted(evals)
            print(f"    Eigenvalues: {[f'{e:.4f}' for e in evals_sorted]}")

            # Null space
            null_evals = [i for i, e in enumerate(evals_sorted) if abs(e) < 1e-4]
            if null_evals:
                print(f"    Null eigenvalues at indices: {null_evals}")
                # Compute eigenvectors
                evals_full, evecs = np.linalg.eigh(np.array(M, dtype=float))
                for idx in null_evals:
                    real_idx = list(sorted(range(len(evals_full)),
                                           key=lambda i: evals_full[i]))[idx]
                    vec = evecs[:, real_idx]
                    # Normalize so max component is 1
                    vec = vec / vec[abs(vec).argmax()]
                    print(f"    Null vector {idx}: [{', '.join(f'{v:.4f}' for v in vec)}]")
                    # Check palindrome of null vector
                    rev = vec[::-1]
                    is_pal = np.allclose(vec, rev, atol=1e-4)
                    is_anti_pal = np.allclose(vec, -rev, atol=1e-4)
                    print(f"      Palindromic: {'YES' if is_pal else 'NO'}, "
                          f"Anti-palindromic: {'YES' if is_anti_pal else 'NO'}")
            else:
                print(f"    No null eigenvalues (M_k is nondegenerate)")

            # Determinant
            det = np.linalg.det(np.array(M, dtype=float))
            print(f"    det(M_{k}) = {det:.6e}")

            # Rank
            rank = np.linalg.matrix_rank(np.array(M, dtype=float), tol=1e-4)
            print(f"    rank(M_{k}) = {rank}")

        except ImportError:
            print(f"    (numpy not available for eigenvalue computation)")


# ====================================================================
# QUESTION 6: Null vectors of M_k — palindromic structure
# ====================================================================

def analyze_null_vectors(even_ks=None, c_val=1.0):
    """Detailed analysis of null vectors of M_k."""
    if even_ks is None:
        even_ks = [4, 6, 8, 10, 12]
    print("\n" + "=" * 100)
    print(f"QUESTION 6: Null vector palindrome analysis (c={c_val})")
    print("=" * 100)

    try:
        import numpy as np
    except ImportError:
        print("  numpy required for this analysis")
        return

    for k in even_ks:
        print(f"\n  k={k}:")
        M, _, _ = compute_depth2_quadratic_form(k, c_val)
        M_np = np.array(M, dtype=float)
        n = k - 1

        evals, evecs = np.linalg.eigh(M_np)
        sorted_idx = np.argsort(evals)
        evals = evals[sorted_idx]
        evecs = evecs[:, sorted_idx]

        # All eigenvectors
        print(f"    Eigenvalue spectrum:")
        for i in range(n):
            vec = evecs[:, i]
            vec = vec / vec[np.abs(vec).argmax()]
            rev = vec[::-1]
            pal = "PAL" if np.allclose(vec, rev, atol=1e-3) else ""
            anti = "ANTI-PAL" if np.allclose(vec, -rev, atol=1e-3) else ""
            sym_type = pal or anti or "MIXED"
            print(f"      eval={evals[i]:>12.4f}  vec=[{', '.join(f'{v:>8.4f}' for v in vec)}]  {sym_type}")


# ====================================================================
# COMPREHENSIVE FACTORIZATION TEST (k=4 through k=20)
# ====================================================================

def comprehensive_factorization_test(max_k=20, c_val=1.0, n_samples=80, seed=314):
    """For each even k from 4 to max_k, test:
    - Palindrome factors (lambda_i - lambda_{k-i})
    - Complementary sum factors (lambda_i + lambda_{k-i})
    - Linear factors involving adjacent pairs
    - Whether the T-sector polynomial is irreducible

    Uses the ratio test: if setting a constraint makes the T-sector vanish
    relative to the generic case, the constraint defines a factor.
    """
    print("\n" + "=" * 100)
    print(f"COMPREHENSIVE FACTORIZATION TEST (c={c_val})")
    print("=" * 100)

    engine = StasheffEngine(c_val)
    rng = random.Random(seed)

    for k in range(4, max_k + 1, 2):
        n_lams = k - 1
        print(f"\n  k={k} (n_lams={n_lams}):")
        t0 = time.time()

        # First compute generic magnitude
        max_generic = 0.0
        for _ in range(n_samples):
            engine._cache.clear()
            lams = tuple(rng.uniform(0.3, 2.5) for _ in range(n_lams))
            result = engine.mk(lams)
            T_val = sum(abs(v) for f, v in result.items() if f >= 0)
            max_generic = max(max_generic, T_val)

        # Test palindrome factors (lambda_i - lambda_{n-i})
        pal_factors = []
        for i in range(n_lams // 2):
            j = n_lams - 1 - i  # palindrome partner
            max_constrained = 0.0
            for _ in range(n_samples):
                engine._cache.clear()
                lams = [rng.uniform(0.3, 2.5) for _ in range(n_lams)]
                lams[j] = lams[i]  # set lambda_{j+1} = lambda_{i+1}
                result = engine.mk(tuple(lams))
                T_val = sum(abs(v) for f, v in result.items() if f >= 0)
                max_constrained = max(max_constrained, T_val)
            ratio = max_constrained / max_generic if max_generic > 0 else 0
            is_factor = ratio < 1e-7
            pal_factors.append((i + 1, j + 1, ratio, is_factor))

        # Test complementary sum factors (lambda_i + lambda_{n-i})
        sum_factors = []
        for i in range(n_lams // 2):
            j = n_lams - 1 - i
            max_constrained = 0.0
            for _ in range(n_samples):
                engine._cache.clear()
                lams = [rng.uniform(0.3, 2.5) for _ in range(n_lams)]
                lams[j] = -lams[i]  # set lambda_{j+1} = -lambda_{i+1}
                result = engine.mk(tuple(lams))
                T_val = sum(abs(v) for f, v in result.items() if f >= 0)
                max_constrained = max(max_constrained, T_val)
            ratio = max_constrained / max_generic if max_generic > 0 else 0
            is_factor = ratio < 1e-7
            sum_factors.append((i + 1, j + 1, ratio, is_factor))

        # Test adjacent difference factors (lambda_i - lambda_{i+1})
        adj_factors = []
        for i in range(n_lams - 1):
            j = i + 1
            max_constrained = 0.0
            for _ in range(n_samples):
                engine._cache.clear()
                lams = [rng.uniform(0.3, 2.5) for _ in range(n_lams)]
                lams[j] = lams[i]
                result = engine.mk(tuple(lams))
                T_val = sum(abs(v) for f, v in result.items() if f >= 0)
                max_constrained = max(max_constrained, T_val)
            ratio = max_constrained / max_generic if max_generic > 0 else 0
            is_factor = ratio < 1e-7
            adj_factors.append((i + 1, j + 1, ratio, is_factor))

        # Test skip-1 factors (lambda_i - lambda_{i+2})
        skip_factors = []
        for i in range(n_lams - 2):
            j = i + 2
            max_constrained = 0.0
            for _ in range(n_samples):
                engine._cache.clear()
                lams = [rng.uniform(0.3, 2.5) for _ in range(n_lams)]
                lams[j] = lams[i]
                result = engine.mk(tuple(lams))
                T_val = sum(abs(v) for f, v in result.items() if f >= 0)
                max_constrained = max(max_constrained, T_val)
            ratio = max_constrained / max_generic if max_generic > 0 else 0
            is_factor = ratio < 1e-7
            skip_factors.append((i + 1, j + 1, ratio, is_factor))

        elapsed = time.time() - t0

        # Report
        n_pal = sum(1 for _, _, _, f in pal_factors if f)
        n_sum = sum(1 for _, _, _, f in sum_factors if f)
        n_adj = sum(1 for _, _, _, f in adj_factors if f)
        n_skip = sum(1 for _, _, _, f in skip_factors if f)
        total_linear = n_pal + n_sum + n_adj + n_skip

        print(f"    Palindrome factors: {n_pal}/{len(pal_factors)}")
        for i, j, ratio, is_f in pal_factors:
            marker = "*** FACTOR ***" if is_f else f"ratio={ratio:.2e}"
            print(f"      (lam_{i} - lam_{j}): {marker}")

        print(f"    Complementary sum factors: {n_sum}/{len(sum_factors)}")
        for i, j, ratio, is_f in sum_factors:
            marker = "*** FACTOR ***" if is_f else f"ratio={ratio:.2e}"
            print(f"      (lam_{i} + lam_{j}): {marker}")

        if n_adj > 0:
            print(f"    Adjacent diff factors: {n_adj}/{len(adj_factors)}")
            for i, j, ratio, is_f in adj_factors:
                if is_f:
                    print(f"      (lam_{i} - lam_{j}): *** FACTOR ***")

        if n_skip > 0:
            print(f"    Skip-1 factors: {n_skip}/{len(skip_factors)}")
            for i, j, ratio, is_f in skip_factors:
                if is_f:
                    print(f"      (lam_{i} - lam_{j}): *** FACTOR ***")

        degree_str = f"degree >= {k-1 - total_linear}" if total_linear > 0 else f"degree {k-1}"
        irred_str = "IRREDUCIBLE" if total_linear == 0 else f"{total_linear} linear factors"
        print(f"    => {irred_str}, remaining polynomial {degree_str}")
        print(f"    [{elapsed:.1f}s]")


# ====================================================================
# EVEN-ODD COMPREHENSIVE SUMMARY
# ====================================================================

def comprehensive_summary(max_k=20, c_val=1.0):
    """Single-pass: compute all m_k at symmetric point for k=2,...,max_k.
    Report T-coefficient, signed sum, abs sum, scalar, match to Catalan (odd k).
    """
    print("\n" + "=" * 100)
    print(f"COMPREHENSIVE SUMMARY: m_k at symmetric point (c={c_val})")
    print("=" * 100)

    engine = StasheffEngine(c_val)

    print(f"\n  {'k':>3} {'parity':>6} {'T-coeff':>20} {'signed_sum':>20} "
          f"{'|T-sector|':>20} {'Catalan_pred':>20} {'match':>6}")
    print("  " + "-" * 100)

    for k in range(2, max_k + 1):
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(k - 1))
        t0 = time.time()
        result = engine.mk(lams)
        elapsed = time.time() - t0

        fields = {f: v for f, v in result.items() if f >= 0}
        T_val = result.get(0, 0.0)
        signed = sum(v for v in fields.values())
        absval = sum(abs(v) for v in fields.values())
        parity = "even" if k % 2 == 0 else "odd"

        if k % 2 == 1 and k >= 3:
            n = (k - 3) // 2
            C_n = math.comb(2 * n, n) // (n + 1) if n >= 0 else 1
            pred = (-1)**n * C_n * math.factorial(k)
            match = "YES" if abs(T_val - pred) < max(1, abs(pred) * 1e-6) else "NO"
            pred_str = f"{pred}"
        else:
            pred_str = "N/A"
            match = "vanish" if absval < 1e-6 else "NONZERO"

        suffix = f"  [{elapsed:.1f}s]" if elapsed > 1.0 else ""
        print(f"  {k:>3} {parity:>6} {T_val:>20.1f} {signed:>20.1f} "
              f"{absval:>20.6e} {pred_str:>20} {match:>6}{suffix}")


# ====================================================================
# MAIN
# ====================================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Virasoro palindrome/factorization structure")
    parser.add_argument('--max-k', type=int, default=14,
                        help='Maximum arity (default 14; 20 is slow)')
    parser.add_argument('--c', type=float, default=1.0, help='Central charge')
    parser.add_argument('--question', type=int, default=0,
                        help='Run only question N (0=all)')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: fewer samples, smaller max_k')
    args = parser.parse_args()

    max_k = args.max_k
    c_val = args.c
    if args.quick:
        max_k = min(max_k, 10)

    if args.question in (0, 1):
        test_symmetric_point_vanishing(max_k=max_k, c_val=c_val)

    if args.question in (0, 2):
        test_minimum_depth(max_k=max_k, c_val=c_val)

    if args.question in (0, 4):
        verify_catalan_pattern(max_k=max_k if max_k % 2 == 1 else max_k - 1, c_val=c_val)

    if args.question in (0, 3):
        comprehensive_factorization_test(max_k=min(max_k, 16), c_val=c_val)

    if args.question in (0, 5):
        even_ks = [k for k in range(4, min(max_k, 14) + 1, 2)]
        analyze_quadratic_forms(even_ks=even_ks, c_val=c_val)

    if args.question in (0, 6):
        even_ks = [k for k in range(4, min(max_k, 14) + 1, 2)]
        analyze_null_vectors(even_ks=even_ks, c_val=c_val)

    if args.question == 0:
        comprehensive_summary(max_k=max_k, c_val=c_val)


if __name__ == '__main__':
    main()
