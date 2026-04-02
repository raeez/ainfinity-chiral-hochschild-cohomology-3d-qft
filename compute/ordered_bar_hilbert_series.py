r"""Ordered bar complex: full dimension sequences and Hilbert series.

Computes dim(B^{ord}_n), Hilbert series, bar cohomology, and Euler
characteristics for all standard chiral algebra families.

MATHEMATICAL FRAMEWORK:

There are THREE distinct bar complexes (AP37):

(a) FG bar B^{FG}(A): uses ONLY the zeroth product a_{(0)}b.
    This is the bar complex of A as a chiral Lie algebra.
    For A = V_k(g), this is the Chevalley-Eilenberg complex C^*(g).

(b) Full symmetric bar B^{Sigma}(A): uses ALL OPE products with
    Sigma_n-coinvariants. This is Vol I Theorem A.

(c) Ordered bar B^{ord}(A): uses ALL OPE products, retains ordering.
    No Sigma_n quotient. This is the object of Part VII.

For computation, the DEPTH FILTRATION separates:
  - Depth 0: the zeroth product {a_{(0)} b} (Lie bracket for KM)
  - Depth 1: the first product {a_{(1)} b} (central term k*kappa for KM)
  - Depth p: the p-th product {a_{(p)} b}

At depth 0, B^{ord}(A) reduces to the bar complex of the Lie algebra g.
This is NOT the associative bar complex of U(g) restricted to generators.
Rather, it is the KOSZUL COMPLEX / CE complex of g, which uses the
exterior coalgebra structure on g (not the tensor coalgebra).

CRITICAL DISTINCTION:
  - The ORDERED bar complex B^{ord}(A) is the bar complex of A as an
    ASSOCIATIVE chiral algebra, using the tensor coalgebra T^c(sA).
  - At depth 0, the differential uses only m_2 = Lie bracket.
  - But m_2 = Lie bracket is NOT associative, so d^2 != 0 on T^c(sg).
  - The resolution: the ordered bar complex of V_k(g) uses the
    FULL associative product in U(g), not just the Lie bracket.
  - The KOSZUL COMPLEX of g (exterior coalgebra) is the SYMMETRIC bar
    complex, not the ordered one.

For the ORDERED bar complex of V_k(g) at depth 0:
  This is the bar complex B(U(g)) of the universal enveloping algebra,
  which is infinite-dimensional. The standard resolution:
  B_n(U(g)) = U(g)^{otimes n} (n copies of U(g)).

WHAT WE COMPUTE:
  1. Generator-level dimensions (counting tensor monomials)
  2. The Koszul / CE complex (the symmetric / antisymmetric version)
  3. Hilbert series (closed form)
  4. Graded characters with conformal weight
  5. Comparison of ordered vs symmetric bar

References:
  ordered_associative_chiral_kd_core.tex
  ordered_e1_shadow_sl2.py
  collision_residue_rmatrix.py
  Loday-Vallette: Algebraic Operads (2012), Chapter 2
"""

from fractions import Fraction
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
from math import comb, factorial
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from lib.ordered_chiral_kd_engine import (
    bar_differential_ordered,
    d_squared_check,
    shuffle_product,
)


# =========================================================================
# UTILITY
# =========================================================================

def enumerate_words(generators: List[str], length: int) -> List[Tuple[str, ...]]:
    """All ordered words of given length."""
    if length == 0:
        return [()]
    if length == 1:
        return [(g,) for g in generators]
    shorter = enumerate_words(generators, length - 1)
    return [w + (g,) for w in shorter for g in generators]


def matrix_rank_fraction(matrix: List[List[Fraction]]) -> int:
    """Rank of a Fraction matrix via row echelon."""
    if not matrix or not matrix[0]:
        return 0
    nrows = len(matrix)
    ncols = len(matrix[0])
    M = [row[:] for row in matrix]
    rank = 0
    for col in range(ncols):
        pivot_row = None
        for row in range(rank, nrows):
            if M[row][col] != Fraction(0):
                pivot_row = row
                break
        if pivot_row is None:
            continue
        M[rank], M[pivot_row] = M[pivot_row], M[rank]
        pivot_val = M[rank][col]
        for row in range(nrows):
            if row != rank and M[row][col] != Fraction(0):
                factor = M[row][col] / pivot_val
                for c in range(ncols):
                    M[row][c] -= factor * M[rank][c]
        rank += 1
    return rank


def compute_tensor_dim(dims: Dict[int, int], n: int, total_w: int) -> int:
    """Dimension of n-fold tensor at total weight total_w."""
    if n == 0:
        return 1 if total_w == 0 else 0
    if n == 1:
        return dims.get(total_w, 0)
    result = 0
    min_w = min((w for w in dims if dims[w] > 0), default=0)
    for w1 in range(min_w, total_w - (n - 1) * min_w + 1):
        if dims.get(w1, 0) > 0:
            sub = compute_tensor_dim(dims, n - 1, total_w - w1)
            result += dims[w1] * sub
    return result


def partitions_min_part(w, min_part=2):
    """Number of partitions of w with all parts >= min_part."""
    if w == 0:
        return 1
    if w < min_part:
        return 0
    dp = [0] * (w + 1)
    dp[0] = 1
    for part in range(min_part, w + 1):
        for j in range(part, w + 1):
            dp[j] += dp[j - part]
    return dp[w]


# =========================================================================
# CE COMPLEX COMPUTATION
# =========================================================================

def build_ce_differential_matrix(
    generators: List[str],
    struct: Dict,
    p: int,
) -> List[List[Fraction]]:
    """Build the matrix of the CE differential d: Lambda^p(g) -> Lambda^{p+1}(g).

    The CE differential on the exterior algebra is:
      d(x_1 ^ ... ^ x_p) = sum_{i<j} (-1)^{i+j} [x_i, x_j] ^ x_1 ^ ... hat_i ... hat_j ... ^ x_p

    We represent Lambda^p(g) by sorted p-subsets of generators (with a fixed ordering).
    """
    n = len(generators)
    gen_idx = {g: i for i, g in enumerate(generators)}

    # Basis for Lambda^p: sorted p-element subsets
    from itertools import combinations
    source_basis = list(combinations(range(n), p))
    target_basis = list(combinations(range(n), p + 1))

    source_idx = {b: i for i, b in enumerate(source_basis)}
    target_idx = {b: i for i, b in enumerate(target_basis)}

    nrows = len(target_basis)
    ncols = len(source_basis)
    matrix = [[Fraction(0)] * ncols for _ in range(nrows)]

    for col_j, subset in enumerate(source_basis):
        # Apply CE differential
        subset_list = list(subset)
        for ii in range(p):
            for jj in range(ii + 1, p):
                a = generators[subset_list[ii]]
                b = generators[subset_list[jj]]
                bracket = struct.get((a, b), {})

                # Sign: (-1)^{ii + jj}
                sign = (-1) ** (ii + jj)

                for result_gen, coeff in bracket.items():
                    if coeff == Fraction(0):
                        continue
                    result_idx = gen_idx[result_gen]

                    # Build the new (p-1)-subset: remove positions ii, jj
                    remaining = [subset_list[k] for k in range(p) if k != ii and k != jj]

                    # Insert result_idx and sort to get a (p)-subset
                    # Wait: we removed 2 and add 1, so we go from p to p-1.
                    # CE differential maps Lambda^p -> Lambda^{p+1}? No!
                    # CE cohomological differential: d: C^p -> C^{p+1}
                    # CE homological differential: d: Lambda^p -> Lambda^{p-1}
                    # For CE HOMOLOGY: d maps Lambda^p(g) -> Lambda^{p-1}(g)

                    # Let me use the HOMOLOGICAL convention:
                    # d(x_1 ^ ... ^ x_p) = sum_{i<j} (-1)^{i+j+1} [x_i,x_j] ^ hat_i ^ hat_j
                    # Result is in Lambda^{p-1}
                    pass

        # Restart with correct convention for HOMOLOGY
        pass

    # Better approach: use homological CE differential
    # d: Lambda^p(g) -> Lambda^{p-1}(g)
    # d(x_{i_1} ^ ... ^ x_{i_p}) = sum_{s<t} (-1)^{s+t} [x_{i_s}, x_{i_t}] ^ x_{i_1} ^ ... hat_s ... hat_t ... ^ x_{i_p}
    # Result is in Lambda^{p-1}(g)

    return None  # Will implement below


def compute_ce_homology(generators: List[str], struct: Dict, max_p: int):
    """Compute CE homology of a Lie algebra.

    Uses the Koszul / CE complex: Lambda^p(g) with differential
    d: Lambda^p -> Lambda^{p-1}.

    d(e_{i_1} ^ ... ^ e_{i_p}) =
      sum_{1<=s<t<=p} (-1)^{s+t} [e_{i_s}, e_{i_t}] ^ e_{i_1} ^ ... hat_s ... hat_t ... ^ e_{i_p}
    """
    from itertools import combinations

    n = len(generators)
    gen_idx = {g: i for i, g in enumerate(generators)}

    results = {}

    # For each p, compute the differential matrix d_p: Lambda^p -> Lambda^{p-1}
    # Basis: sorted index tuples

    # First compute dimensions
    for p in range(0, min(max_p + 2, n + 2)):
        dim = comb(n, p)
        results[p] = {'dim': dim}

    # Compute differential matrices and their ranks
    ranks = {}
    for p in range(1, min(max_p + 2, n + 1)):
        source_basis = list(combinations(range(n), p))
        target_basis = list(combinations(range(n), p - 1))

        source_map = {b: i for i, b in enumerate(source_basis)}
        target_map = {b: i for i, b in enumerate(target_basis)}

        nrows = len(target_basis)
        ncols = len(source_basis)

        if nrows == 0 or ncols == 0:
            ranks[p] = 0
            continue

        matrix = [[Fraction(0)] * ncols for _ in range(nrows)]

        for col_j, subset in enumerate(source_basis):
            subset_list = list(subset)
            for s in range(p):
                for t in range(s + 1, p):
                    a = generators[subset_list[s]]
                    b = generators[subset_list[t]]
                    bracket = struct.get((a, b), {})

                    sign_st = (-1) ** (s + t)

                    for result_gen, coeff in bracket.items():
                        if coeff == Fraction(0):
                            continue
                        result_gidx = gen_idx[result_gen]

                        # Build remaining after removing positions s and t
                        remaining = [subset_list[k] for k in range(p) if k != s and k != t]

                        # Insert result_gidx into remaining and sort
                        new_set = sorted(remaining + [result_gidx])
                        new_tuple = tuple(new_set)

                        # Check for duplicates (result_gidx already in remaining)
                        if len(set(new_set)) < len(new_set):
                            continue  # Exterior algebra: duplicate -> 0

                        if new_tuple not in target_map:
                            continue

                        # Compute the sign from reordering
                        # We need the sign of the permutation that puts
                        # [result_gidx] + remaining into sorted order
                        # The unsorted list is [result_gidx] + remaining
                        unsorted = [result_gidx] + remaining
                        # Count inversions to get the sign
                        reorder_sign = 1
                        for ii in range(len(unsorted)):
                            for jj in range(ii + 1, len(unsorted)):
                                if unsorted[ii] > unsorted[jj]:
                                    reorder_sign *= -1

                        row_i = target_map[new_tuple]
                        matrix[row_i][col_j] += Fraction(sign_st * reorder_sign) * coeff

        ranks[p] = matrix_rank_fraction(matrix)

    # Compute homology: H_p = ker(d_p) / im(d_{p+1})
    # dim H_p = dim Lambda^p - rank(d_p) - rank(d_{p+1})
    for p in range(0, min(max_p + 1, n + 1)):
        dim_p = comb(n, p)
        rank_d_p = ranks.get(p, 0)  # rank of d_p: Lambda^p -> Lambda^{p-1}
        rank_d_above = ranks.get(p + 1, 0)  # rank of d_{p+1}: Lambda^{p+1} -> Lambda^p
        dim_ker = dim_p - rank_d_p
        dim_H = dim_ker - rank_d_above
        results[p] = {
            'dim': dim_p,
            'rank_d': rank_d_p,
            'dim_ker': dim_ker,
            'rank_d_above': rank_d_above,
            'dim_H': dim_H,
        }

    return results


# =========================================================================
# LIE ALGEBRA DATA
# =========================================================================

def build_sln_data(N: int) -> Tuple[List[str], Dict]:
    """Build generators and structure constants for sl_N."""
    generators = []
    gen_to_matrix = {}

    for i in range(N):
        for j in range(N):
            if i != j:
                label = f"E{i+1}{j+1}"
                generators.append(label)
                mat = [[Fraction(0)] * N for _ in range(N)]
                mat[i][j] = Fraction(1)
                gen_to_matrix[label] = mat

    for i in range(N - 1):
        label = f"H{i+1}"
        generators.append(label)
        mat = [[Fraction(0)] * N for _ in range(N)]
        mat[i][i] = Fraction(1)
        mat[i+1][i+1] = Fraction(-1)
        gen_to_matrix[label] = mat

    def mat_mul(A, B):
        n = len(A)
        C = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    C[i][j] += A[i][k] * B[k][j]
        return C

    def mat_sub(A, B):
        n = len(A)
        return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]

    def decompose(mat, N):
        result = {}
        for i in range(N):
            for j in range(N):
                if i != j and mat[i][j] != Fraction(0):
                    label = f"E{i+1}{j+1}"
                    result[label] = mat[i][j]
        diag = [mat[i][i] for i in range(N)]
        for k in range(1, N):
            c_k = sum(diag[:k])
            if c_k != Fraction(0):
                result[f"H{k}"] = c_k
        return result

    struct = {}
    for g1 in generators:
        for g2 in generators:
            m1 = gen_to_matrix[g1]
            m2_mat = gen_to_matrix[g2]
            comm = mat_sub(mat_mul(m1, m2_mat), mat_mul(m2_mat, m1))
            dec = decompose(comm, N)
            struct[(g1, g2)] = {k: Fraction(v) for k, v in dec.items() if v != Fraction(0)}

    return generators, struct


SL2_GENS = ['e', 'h', 'f']
SL2_STRUCT = {}
SL2_STRUCT[('e', 'f')] = {'h': Fraction(1)}
SL2_STRUCT[('f', 'e')] = {'h': Fraction(-1)}
SL2_STRUCT[('h', 'e')] = {'e': Fraction(2)}
SL2_STRUCT[('e', 'h')] = {'e': Fraction(-2)}
SL2_STRUCT[('h', 'f')] = {'f': Fraction(-2)}
SL2_STRUCT[('f', 'h')] = {'f': Fraction(2)}
for x in SL2_GENS:
    for y in SL2_GENS:
        if (x, y) not in SL2_STRUCT:
            SL2_STRUCT[(x, y)] = {}


# =========================================================================
# ORDERED BAR OF U(g): DIMENSION ANALYSIS
# =========================================================================

def ug_hilbert_series_coeffs(N: int, max_degree: int) -> Dict[int, int]:
    r"""Dimension of U(sl_N) at PBW degree d.

    By PBW, U(g) = Sym(g) as graded vector space.
    dim Sym^d(g) = C(dim_g + d - 1, d) where dim_g = N^2 - 1.

    For the ordered bar complex B_n(U(g)):
      B_n = (augmentation ideal of U(g))^{tensor n}
    where the augmentation ideal = U(g)_+ = bigoplus_{d >= 1} Sym^d(g).
    """
    dim_g = N * N - 1
    dims = {}
    for d in range(0, max_degree + 1):
        dims[d] = comb(dim_g + d - 1, d)
    return dims


# =========================================================================
# FAMILY ANALYSES
# =========================================================================

def analyze_heisenberg():
    """Heisenberg H_k: 1 generator J (weight 1), c_0 = 0."""
    print("=" * 72)
    print("FAMILY 1: HEISENBERG H_k")
    print("  Generator: J (weight 1)")
    print("  Lambda-bracket: {J_lambda J} = k*lambda")
    print("  Zeroth product: {J_{(0)} J} = 0  (c_0 = 0)")
    print("  Bar differential d = 0 (zeroth product vanishes)")
    print("=" * 72)

    # The Lie algebra is 1-dimensional abelian: g = k*J, [J,J] = 0.
    # CE complex: Lambda^0(g) = k, Lambda^1(g) = k*J, Lambda^p = 0 for p >= 2.
    # CE differential d = 0 (abelian).
    # CE homology: H_0 = k, H_1 = k, H_p = 0 for p >= 2.

    print("\n--- CE complex (= symmetric/Koszul bar, depth 0) ---")
    print("  Lambda^0(g) = k (dim 1)")
    print("  Lambda^1(g) = k*J (dim 1)")
    print("  Lambda^p(g) = 0 for p >= 2")
    print("  d = 0 (abelian Lie algebra)")
    print("  H_0 = k, H_1 = k")
    print("  Poincare polynomial: 1 + t")

    print("\n--- Ordered bar complex B^{ord}(U(g)) ---")
    print("  U(g) = k[J] (polynomial ring in J).")
    print("  Augmentation ideal: J*k[J] = span{J, J^2, J^3, ...}")
    print("  dim at PBW degree d: 1 for all d >= 1.")

    # The ordered bar complex of U(g) = k[J]:
    # B_n = (J*k[J])^{tensor n}
    # At PBW total degree d, bar degree n:
    # compositions of d into n positive parts = C(d-1, n-1)
    print("\n  dim B^{ord}_n at PBW degree d:")
    print(f"  {'d \\ n':>8}", end="")
    for n in range(1, 8):
        print(f" | {'n='+str(n):>6}", end="")
    print()
    print(f"  {'-'*65}")
    for d in range(1, 11):
        print(f"  {'d='+str(d):>8}", end="")
        for n in range(1, 8):
            if n <= d:
                dim = comb(d - 1, n - 1)
            else:
                dim = 0
            print(f" | {dim:>6}", end="")
        print()

    print(f"\n  Bigraded Hilbert series:")
    print(f"    h(t,s) = sum_{{n>=1, d>=n}} C(d-1,n-1) t^n s^d")
    print(f"           = ts/(1-s)  *  1/(1-ts/(1-s))")
    print(f"           = ts / (1 - s - ts)")
    print(f"    At s=1 (total Hilbert): diverges (U(g) is infinite-dimensional)")
    print(f"    At s=q (conformal weight): h(t,q) = tq / (1 - q - tq)")

    # The bar differential for U(k[J]) = k[J]:
    # d[J^a | J^b] = J^{a+b} (concatenation in the polynomial ring)
    # This is the bar complex of a commutative algebra, which is acyclic
    # in positive degrees: H_0 = k, H_n = 0 for n >= 1.
    # (The normalized bar complex of a polynomial ring is contractible.)
    print(f"\n  Bar cohomology of U(g) = k[J]:")
    print(f"    H_n(B(k[J])) = k for n=0, 0 for n >= 1")
    print(f"    (Polynomial ring is smooth, bar complex is a resolution of k)")
    print(f"    Alternatively: Tor^{{k[J]}}_n(k, k) = Lambda^n_k(dJ) = k for n=0,1, 0 for n>=2")
    print(f"    (Koszul resolution: k[J] tensor Lambda(dJ) with d(dJ) = J)")

    print(f"\n  At the GENERATOR level (only using J, not J^2, J^3,...):")
    print(f"    B^{{ord,gen}}_n = k (1-dimensional for all n)")
    print(f"    The differential d[J|...|J] uses {'{J_{(0)} J}'} = 0")
    print(f"    So d = 0 and H_n = k for all n >= 1")

    print(f"\n  SUMMARY:")
    dim_seq = [1] * 10
    coh_seq = [1] * 10
    print(f"    Generator-level dim B_n: {dim_seq}")
    print(f"    Generator-level dim H_n: {coh_seq}")
    print(f"    Gen-level Hilbert:  h_B(t) = t/(1-t)")
    print(f"    Gen-level h_coh(t): t/(1-t)  (d = 0)")
    print(f"    CE Poincare:  1 + t")
    print(f"    Graded char:  h(t,q) = tq/(1-tq)")

    return {'family': 'H_k', 'gen_dim': [1]*10, 'ce_poincare': [1, 1]}


def analyze_sl2():
    """Affine V_k(sl_2): 3 generators, depth-0 = CE(sl_2)."""
    print("\n" + "=" * 72)
    print("FAMILY 2: AFFINE V_k(sl_2)")
    print("  Generators: e, h, f (weight 1, dim g = 3)")
    print("  Zeroth product: Lie bracket [e,f]=h, [h,e]=2e, [h,f]=-2f")
    print("=" * 72)

    # CE complex
    print("\n--- CE complex (exterior algebra, depth 0) ---")
    ce = compute_ce_homology(SL2_GENS, SL2_STRUCT, 6)

    print(f"\n  {'p':>4} | {'dim Lambda^p':>14} | {'rank d_p':>10} | {'rank d_{p+1}':>14} | {'dim H_p':>10}")
    print(f"  {'-'*60}")
    ce_coh = []
    for p in range(0, 4):
        if p in ce:
            r = ce[p]
            print(f"  {p:>4} | {r['dim']:>14} | {r['rank_d']:>10} | {r['rank_d_above']:>14} | {r['dim_H']:>10}")
            ce_coh.append(r['dim_H'])

    print(f"\n  CE homology: {ce_coh}")
    print(f"  Expected: H_0=1, H_1=0, H_2=0, H_3=1")
    print(f"  Poincare polynomial: 1 + t^3")

    # Ordered bar complex B^{ord}(U(sl_2))
    print(f"\n--- Ordered bar B^{{ord}}(U(sl_2)) ---")
    print(f"  U(sl_2) = k<e,h,f>/(Serre relations)")
    print(f"  By PBW: dim U(sl_2) at degree d = C(3+d-1, d) = C(d+2, 2)")

    print(f"\n  PBW dimensions of U(sl_2):")
    for d in range(0, 8):
        dim_d = comb(d + 2, 2)
        print(f"    degree {d}: dim = {dim_d}")

    # dim of augmentation ideal at degree d (d >= 1)
    print(f"\n  dim B^{{ord}}_n at PBW total degree d:")
    print(f"  {'d \\ n':>8}", end="")
    for n in range(1, 7):
        print(f" | {'n='+str(n):>8}", end="")
    print()
    print(f"  {'-'*65}")

    # Augmentation ideal dims: degree d >= 1, dim = C(d+2,2)
    aug_dims = {d: comb(d + 2, 2) for d in range(1, 20)}

    for d in range(1, 13):
        print(f"  {'d='+str(d):>8}", end="")
        for n in range(1, 7):
            dim = compute_tensor_dim(aug_dims, n, d)
            print(f" | {dim:>8}", end="")
        print()

    # Generator-level ordered bar
    print(f"\n--- Generator-level ordered bar ---")
    print(f"  B^{{ord,gen}}_n = span of all ordered words [g_1|...|g_n]")
    print(f"  with g_i in {{e, h, f}}, dim = 3^n")
    print(f"  The differential d uses m_2 = Lie bracket, which maps to generators.")
    print(f"  BUT d^2 != 0 with the Lie bracket (Lie bracket is not associative).")
    print(f"  The generator-level space is NOT a sub-complex of B^{{ord}}(U(g)).")
    print(f"  It IS a sub-complex of the CE complex (exterior coalgebra).")

    gen_dim_seq = [3**n for n in range(1, 11)]
    print(f"\n  dim B^{{ord,gen}}_n = 3^n: {gen_dim_seq[:6]}")
    print(f"  Hilbert series: h(t) = 3t/(1-3t)")
    print(f"  Graded char: h(t,q) = 3tq/(1-3tq) (all weight 1)")

    # Full ordered bar cohomology = Tor^{U(sl_2)}_*(k, k) = H_CE(sl_2)
    # by the standard comparison: the bar complex of U(g) computes
    # Tor^{U(g)}(k,k) which equals CE homology H_*(g, k).
    print(f"\n  Ordered bar cohomology H_*(B^{{ord}}(U(sl_2))) = H^CE_*(sl_2):")
    print(f"    H_0 = k (1-dim)")
    print(f"    H_1 = 0")
    print(f"    H_2 = 0")
    print(f"    H_3 = k (1-dim)")
    print(f"    H_n = 0 for n >= 4")
    print(f"    Poincare polynomial: 1 + t^3")
    print(f"    Euler characteristic: 1 - 1 = 0")

    return {
        'family': 'V_k(sl_2)',
        'generators': 3,
        'gen_dim': gen_dim_seq,
        'ce_homology': ce_coh,
        'ce_poincare': '1 + t^3',
    }


def analyze_slN(N: int, max_ce_p: int = None):
    """Affine V_k(sl_N) for general N."""
    dim_g = N * N - 1
    if max_ce_p is None:
        max_ce_p = dim_g

    print(f"\n{'=' * 72}")
    print(f"FAMILY 3.{N}: AFFINE V_k(sl_{N})")
    print(f"  dim g = {dim_g}, rank = {N-1}")
    print(f"{'=' * 72}")

    generators, struct = build_sln_data(N)
    assert len(generators) == dim_g, f"Expected {dim_g} generators, got {len(generators)}"

    # CE complex
    feasible_p = max_ce_p
    # Check if computation is feasible
    total_basis = sum(comb(dim_g, p) for p in range(feasible_p + 2))
    while total_basis > 50000 and feasible_p > 3:
        feasible_p -= 1
        total_basis = sum(comb(dim_g, p) for p in range(feasible_p + 2))

    print(f"\n--- CE complex (computing up to degree {feasible_p}) ---")
    ce = compute_ce_homology(generators, struct, feasible_p)

    print(f"\n  {'p':>4} | {'dim Lambda^p':>14} | {'rank d_p':>10} | {'rank d_{p+1}':>14} | {'dim H_p':>10}")
    print(f"  {'-'*60}")
    ce_coh = []
    for p in range(0, feasible_p + 1):
        if p in ce:
            r = ce[p]
            print(f"  {p:>4} | {r['dim']:>14} | {r['rank_d']:>10} | {r['rank_d_above']:>14} | {r['dim_H']:>10}")
            ce_coh.append(r['dim_H'])

    # Expected CE cohomology
    exponents = list(range(1, N))
    degrees = [2*e + 1 for e in exponents]
    print(f"\n  Exponents of sl_{N}: {exponents}")
    print(f"  CE generator degrees: {degrees}")
    print(f"  Expected Poincare: prod (1 + t^d) for d in {degrees}")

    # Compute expected Betti numbers
    max_deg = sum(degrees)
    poincare = [0] * (max_deg + 1)
    poincare[0] = 1
    for d in degrees:
        new_p = [0] * (max_deg + 1)
        for i in range(max_deg + 1):
            if poincare[i] != 0:
                new_p[i] += poincare[i]
                if i + d <= max_deg:
                    new_p[i + d] += poincare[i]
        poincare = new_p

    expected_betti = {}
    for i, c in enumerate(poincare):
        if c > 0:
            expected_betti[i] = c

    print(f"  Expected Betti numbers: {expected_betti}")
    print(f"  Total dim H* = {sum(poincare)} = 2^{N-1}")

    # Verify computed vs expected
    match = True
    for p in range(min(len(ce_coh), max_deg + 1)):
        expected = poincare[p] if p < len(poincare) else 0
        if ce_coh[p] != expected:
            print(f"  MISMATCH at p={p}: computed {ce_coh[p]}, expected {expected}")
            match = False
    if match and len(ce_coh) > 0:
        print(f"  Computed CE homology MATCHES expected up to degree {len(ce_coh)-1}.")

    # Generator-level bar
    gen_dim_seq = [dim_g**n for n in range(1, 7)]
    print(f"\n  Generator-level dim B_n = {dim_g}^n: {gen_dim_seq}")
    print(f"  Hilbert series: h(t) = {dim_g}t/(1-{dim_g}t)")
    print(f"  Graded char: h(t,q) = {dim_g}tq/(1-{dim_g}tq) (all weight 1)")

    return {
        'family': f'V_k(sl_{N})',
        'dim_g': dim_g,
        'gen_dim': gen_dim_seq,
        'ce_coh': ce_coh,
        'expected_poincare': poincare,
    }


def analyze_virasoro():
    """Virasoro Vir_c: 1 generator T (weight 2)."""
    print(f"\n{'=' * 72}")
    print(f"FAMILY 4: VIRASORO Vir_c")
    print(f"  Generator: T (weight 2)")
    print(f"  OPE: T(z)T(w) ~ (c/2)/(z-w)^4 + 2T(w)/(z-w)^2 + dT(w)/(z-w)")
    print(f"  Lambda-bracket: {{T_lambda T}} = dT + 2*lambda*T + (c/12)*lambda^3")
    print(f"  Zeroth product: {{T_{{(0)}} T}} = dT")
    print(f"{'=' * 72}")

    print(f"\n  CRITICAL: The zeroth product dT is a DESCENDANT (weight 3),")
    print(f"  not a generator. The generator-level bar complex is NOT closed")
    print(f"  under d for non-abelian Virasoro.")

    # Generator-level analysis
    print(f"\n--- Generator-level ordered bar ---")
    print(f"  B^{{ord,gen}}_n = k for all n (only one generator T)")
    print(f"  d[T|T] = dT (exits the generator space)")
    print(f"  So: the generator-level complex is trivially acyclic")
    print(f"  (everything is cocycles, nothing is boundaries within generators)")

    gen_dim = [1] * 10
    print(f"  dim B^{{ord,gen}}_n: {gen_dim}")
    print(f"  Hilbert series: h(t) = t/(1-t)")

    # Full bar complex
    print(f"\n--- Full bar complex (all descendants) ---")
    print(f"  The augmentation ideal of V(c,0) (vacuum Verma) is spanned by")
    print(f"  all states of weight >= 2: L_{{-n1}} ... L_{{-nk}} |0> with sum n_i >= 2.")

    aug_dims = {}
    for w in range(0, 20):
        if w >= 2:
            aug_dims[w] = partitions_min_part(w, 2)
        else:
            aug_dims[w] = 0  # Augmentation ideal: no vacuum, no weight 1

    print(f"\n  Augmentation ideal dimensions by weight:")
    for w in range(2, 16):
        print(f"    weight {w}: dim = {aug_dims[w]}")

    print(f"\n  dim B^{{ord,full}}_{'{n,w}'} (bar degree n, conformal weight w):")
    print(f"  {'w \\ n':>8}", end="")
    for n in range(1, 7):
        print(f" | {'n='+str(n):>8}", end="")
    print()
    print(f"  {'-'*65}")

    total_by_n = defaultdict(int)
    for w in range(2, 16):
        print(f"  {'w='+str(w):>8}", end="")
        for n in range(1, 7):
            d = compute_tensor_dim(aug_dims, n, w)
            total_by_n[n] += d
            print(f" | {d:>8}", end="")
        print()

    print(f"\n  Total (up to weight 15):")
    for n in range(1, 7):
        print(f"    B^{{ord,full}}_{n} (wt<=15): {total_by_n[n]}")

    print(f"\n  Graded Hilbert series:")
    print(f"    h(t,q) = chi_aug(q)*t / (1 - chi_aug(q)*t)")
    print(f"    chi_aug(q) = prod_{{n>=2}} 1/(1-q^n) - 1")
    print(f"               = q^2 + q^3 + 2q^4 + 2q^5 + 4q^6 + 4q^7 + 7q^8 + ...")

    # Bar cohomology
    print(f"\n--- Bar cohomology ---")
    print(f"  H_*(B^{{ord}}(U(Vir))) = Tor^{{U(Vir)}}_*(k, k) = H^CE_*(Vir)")
    print(f"  The Virasoro algebra is infinite-dimensional: dim Vir = infinity.")
    print(f"  CE homology of Vir is known (Feigin-Fuchs, Goncharova):")
    print(f"    H_0 = k")
    print(f"    H_1 = 0 (Vir has no 1-cohomology)")
    print(f"    H_2 = k (the Virasoro cocycle: c/12 * omega)")
    print(f"    H_n: computed by Goncharova (1973), Feigin-Fuchs (1984)")
    print(f"    Poincare series: prod_{{n>=1}} (1 + t^{{2n-1}}) (for Lie_+(Vir))")
    print(f"    = (1+t)(1+t^3)(1+t^5)... = 1 + t + t^3 + t^4 + t^5 + t^6 + ...")

    print(f"\n  For the FULL Virasoro algebra Vir = Lie_+(Vir) + k*C:")
    print(f"    The central element C contributes an extra factor (1+t) to")
    print(f"    the Poincare series.")

    return {
        'family': 'Virasoro',
        'gen_dim': gen_dim,
        'aug_dims': aug_dims,
    }


def analyze_w3():
    """W_3 algebra: generators T (wt 2), W (wt 3)."""
    print(f"\n{'=' * 72}")
    print(f"FAMILY 5: W_3 ALGEBRA")
    print(f"  Generators: T (weight 2), W (weight 3)")
    print(f"{'=' * 72}")

    print(f"\n--- Generator-level ordered bar ---")
    print(f"  B^{{ord,gen}}_n has 2^n elements (2 generators)")
    gen_dim = [2**n for n in range(1, 11)]
    print(f"  dim B_n: {gen_dim[:8]}")
    print(f"  Hilbert series: h(t) = 2t/(1-2t)")

    print(f"\n  Graded by conformal weight:")
    print(f"  At degree n, weight w: count words with k W's and (n-k) T's")
    print(f"  where 2(n-k) + 3k = w, so k = w - 2n. Need 0 <= k <= n.")
    print(f"  dim = C(n, w-2n) if 0 <= w-2n <= n, else 0.")

    print(f"\n  {'w \\ n':>8}", end="")
    for n in range(1, 9):
        print(f" | {'n='+str(n):>6}", end="")
    print()
    print(f"  {'-'*75}")

    for w in range(2, 25):
        row_nonzero = False
        for n in range(1, 9):
            k = w - 2 * n
            if 0 <= k <= n:
                row_nonzero = True
                break
        if not row_nonzero:
            continue
        print(f"  {'w='+str(w):>8}", end="")
        for n in range(1, 9):
            k = w - 2 * n
            if 0 <= k <= n:
                dim = comb(n, k)
            else:
                dim = 0
            print(f" | {dim:>6}", end="")
        print()

    print(f"\n  Graded Hilbert series:")
    print(f"    h(t,q) = t(q^2+q^3) / (1 - t(q^2+q^3))")

    # Full bar complex
    print(f"\n--- Full bar complex ---")
    print(f"  W_3 vacuum module character (augmentation ideal):")
    print(f"    chi_aug(q) = prod_{{n>=2}} 1/(1-q^n) * prod_{{n>=3}} 1/(1-q^n) - 1")

    # Compute W_3 augmentation ideal dims
    # The W_3 has two families of descendants:
    # From T: L_{-n} for n >= 2 (Virasoro modes)
    # From W: W_{-n} for n >= 3 (W-modes)
    # The vacuum module is the free algebra generated by these, modulo null vectors.
    # At generic c, there are no null vectors, so:
    # chi_vac(q) = prod_{n>=2} 1/(1-q^n) * prod_{n>=3} 1/(1-q^n)

    def w3_aug_dim(max_w):
        """Dimensions of W_3 vacuum module (augmentation ideal)."""
        dims = [0] * (max_w + 1)
        dims[0] = 1

        # Multiply by prod_{n>=2} 1/(1-q^n) (Virasoro part)
        for n in range(2, max_w + 1):
            for w in range(n, max_w + 1):
                dims[w] += dims[w - n]

        # Multiply by prod_{n>=3} 1/(1-q^n) (W part)
        for n in range(3, max_w + 1):
            for w in range(n, max_w + 1):
                dims[w] += dims[w - n]

        # Wait -- this double-counts. Let me use proper convolution.
        # Actually, the correct character is:
        # chi(q) = prod_{n>=2} 1/(1-q^n) * prod_{n>=3} 1/(1-q^n)
        # = [prod_{n>=2} 1/(1-q^n)]  *  [prod_{n>=3} 1/(1-q^n)]
        # This is a PRODUCT of two partition generating functions.

        # Compute each factor separately, then convolve.
        vir = [0] * (max_w + 1)
        vir[0] = 1
        for n in range(2, max_w + 1):
            for w in range(n, max_w + 1):
                vir[w] += vir[w - n]

        w_part = [0] * (max_w + 1)
        w_part[0] = 1
        for n in range(3, max_w + 1):
            for w in range(n, max_w + 1):
                w_part[w] += w_part[w - n]

        # Convolve
        result = [0] * (max_w + 1)
        for i in range(max_w + 1):
            for j in range(max_w + 1 - i):
                result[i + j] += vir[i] * w_part[j]

        # Augmentation ideal: subtract the vacuum (weight 0)
        result[0] -= 1
        return {w: result[w] for w in range(max_w + 1)}

    w3_aug = w3_aug_dim(15)
    print(f"\n  W_3 augmentation ideal dims:")
    for w in range(2, 16):
        if w3_aug.get(w, 0) > 0:
            print(f"    weight {w}: dim = {w3_aug[w]}")

    print(f"\n  dim B^{{ord,full}}_{'{n,w}'} for W_3:")
    print(f"  {'w \\ n':>8}", end="")
    for n in range(1, 6):
        print(f" | {'n='+str(n):>8}", end="")
    print()
    print(f"  {'-'*55}")

    for w in range(2, 14):
        if w3_aug.get(w, 0) > 0 or w <= 5:
            print(f"  {'w='+str(w):>8}", end="")
            for n in range(1, 6):
                d = compute_tensor_dim(w3_aug, n, w)
                print(f" | {d:>8}", end="")
            print()

    return {
        'family': 'W_3',
        'gen_dim': gen_dim,
        'weights': [2, 3],
    }


def analyze_betagamma():
    """Beta-gamma: generators beta (wt 1), gamma (wt 0)."""
    print(f"\n{'=' * 72}")
    print(f"FAMILY 6: BETA-GAMMA SYSTEM")
    print(f"  Generators: beta (weight 1), gamma (weight 0)")
    print(f"  OPE: beta(z)gamma(w) ~ 1/(z-w)")
    print(f"  Lambda-bracket: {{beta_lambda gamma}} = 1")
    print(f"  Zeroth product: {{beta_{{(0)}} gamma}} = 1 (scalar)")
    print(f"{'=' * 72}")

    print(f"\n--- Generator-level ordered bar ---")
    gen_dim = [2**n for n in range(1, 11)]
    print(f"  dim B_n = 2^n: {gen_dim[:8]}")
    print(f"  Hilbert series: h(t) = 2t/(1-2t)")

    print(f"\n  The zeroth product gives the Heisenberg Lie algebra h_1:")
    print(f"    [beta, gamma] = 1  (central)")
    print(f"  This is a 2-dim non-abelian Lie algebra (plus the central element 1).")

    # The Lie algebra structure: h_1 = span(beta, gamma, C)
    # with [beta, gamma] = C, [beta, C] = [gamma, C] = 0.
    # If we work with just {beta, gamma} and the bracket landing in scalars,
    # the CE complex is the exterior algebra on {beta, gamma} with
    # d(beta ^ gamma) = [beta, gamma] = 1 (a scalar).

    # CE of h_1 (3-dim Heisenberg):
    # Lambda^0 = k, Lambda^1 = k^3, Lambda^2 = k^3, Lambda^3 = k
    # Poincare: (1+t)^3 if abelian, but h_1 is not abelian.
    # The CE differential: d: Lambda^2 -> Lambda^1 sends
    # beta^gamma -> [beta,gamma] = C, beta^C -> 0, gamma^C -> 0
    # d: Lambda^3 -> Lambda^2 sends
    # beta^gamma^C -> [beta,gamma]^C - [beta,C]^gamma + [gamma,C]^beta = C^C = 0

    h1_gens = ['beta', 'gamma', 'C']
    h1_struct = {}
    h1_struct[('beta', 'gamma')] = {'C': Fraction(1)}
    h1_struct[('gamma', 'beta')] = {'C': Fraction(-1)}
    for x in h1_gens:
        for y in h1_gens:
            if (x, y) not in h1_struct:
                h1_struct[(x, y)] = {}

    print(f"\n--- CE complex of h_1 (3-dim Heisenberg Lie algebra) ---")
    ce = compute_ce_homology(h1_gens, h1_struct, 3)
    ce_coh = []
    print(f"\n  {'p':>4} | {'dim Lambda^p':>14} | {'rank d_p':>10} | {'rank d_{p+1}':>14} | {'dim H_p':>10}")
    print(f"  {'-'*60}")
    for p in range(0, 4):
        if p in ce:
            r = ce[p]
            print(f"  {p:>4} | {r['dim']:>14} | {r['rank_d']:>10} | {r['rank_d_above']:>14} | {r['dim_H']:>10}")
            ce_coh.append(r['dim_H'])

    print(f"\n  CE homology of h_1: {ce_coh}")
    print(f"  Expected: H_0=1, H_1=2, H_2=2, H_3=1  (6 total)")
    print(f"  Note: for h_1, H*(h_1) is NOT the exterior algebra.")
    print(f"  It depends on whether we use trivial or adjoint coefficients.")

    print(f"\n  Graded by conformal weight (beta=wt 1, gamma=wt 0):")
    print(f"  {'w \\ n':>8}", end="")
    for n in range(1, 9):
        print(f" | {'n='+str(n):>6}", end="")
    print()
    print(f"  {'-'*75}")
    for w in range(0, 10):
        print(f"  {'w='+str(w):>8}", end="")
        for n in range(1, 9):
            # w = number of beta's, C(n,w) if 0<=w<=n
            dim = comb(n, w) if 0 <= w <= n else 0
            print(f" | {dim:>6}", end="")
        print()

    print(f"\n  Graded Hilbert series: h(t,q) = t(1+q) / (1 - t(1+q))")

    return {
        'family': 'beta-gamma',
        'gen_dim': gen_dim,
        'weights': [1, 0],
        'ce_coh': ce_coh,
    }


# =========================================================================
# LATEX TABLES
# =========================================================================

def generate_latex_tables():
    """Generate LaTeX-ready tables."""
    print(f"\n{'=' * 72}")
    print(f"LATEX-READY TABLES")
    print(f"{'=' * 72}")

    print(r"""
% TABLE 1: Generator-level bar complex dimensions
\begin{table}[ht]
\centering
\caption{Ordered bar complex dimensions $\dim B^{\mathrm{ord,gen}}_n$ at the generator level.}
\label{tab:ordered-bar-gen-dim}
\begin{tabular}{l c r r r r r r c}
\toprule
Family & $g$ & $n{=}1$ & $n{=}2$ & $n{=}3$ & $n{=}4$ & $n{=}5$ & $n{=}6$ & $h_B(t)$ \\
\midrule
$H_k$ & $1$ & $1$ & $1$ & $1$ & $1$ & $1$ & $1$ & $\frac{t}{1-t}$ \\
$V_k(\mathfrak{sl}_2)$ & $3$ & $3$ & $9$ & $27$ & $81$ & $243$ & $729$ & $\frac{3t}{1-3t}$ \\
$V_k(\mathfrak{sl}_3)$ & $8$ & $8$ & $64$ & $512$ & $4096$ & $32768$ & $262144$ & $\frac{8t}{1-8t}$ \\
$V_k(\mathfrak{sl}_4)$ & $15$ & $15$ & $225$ & $3375$ & $50625$ & $759375$ & $\cdots$ & $\frac{15t}{1-15t}$ \\
$V_k(\mathfrak{sl}_N)$ & $N^2{-}1$ & $N^2{-}1$ & $(N^2{-}1)^2$ & $\cdots$ & $\cdots$ & $\cdots$ & $\cdots$ & $\frac{(N^2{-}1)t}{1-(N^2{-}1)t}$ \\
$\mathrm{Vir}_c$ & $1$ & $1$ & $1$ & $1$ & $1$ & $1$ & $1$ & $\frac{t}{1-t}$ \\
$\mathcal{W}_3$ & $2$ & $2$ & $4$ & $8$ & $16$ & $32$ & $64$ & $\frac{2t}{1-2t}$ \\
$\beta\gamma$ & $2$ & $2$ & $4$ & $8$ & $16$ & $32$ & $64$ & $\frac{2t}{1-2t}$ \\
\bottomrule
\end{tabular}
\end{table}
""")

    print(r"""
% TABLE 2: Chevalley-Eilenberg cohomology (= bar cohomology at depth 0)
\begin{table}[ht]
\centering
\caption{Chevalley--Eilenberg homology $H^{\mathrm{CE}}_p(\mathfrak{g})$ (ordered bar cohomology at depth~$0$).
For $\mathfrak{sl}_N$, the Poincar\'e polynomial is $\prod_{i=1}^{N-1}(1+t^{2i+1})$.}
\label{tab:ce-homology}
\begin{tabular}{l r r r r r r r r r l}
\toprule
$\mathfrak{g}$ & $H_0$ & $H_1$ & $H_2$ & $H_3$ & $H_4$ & $H_5$ & $H_6$ & $H_7$ & $H_8$ & Poincar\'e \\
\midrule
$k$ (abelian) & $1$ & $1$ & $0$ & $0$ & $0$ & $0$ & $0$ & $0$ & $0$ & $1+t$ \\
$\mathfrak{sl}_2$ & $1$ & $0$ & $0$ & $1$ & $0$ & $0$ & $0$ & $0$ & $0$ & $1+t^3$ \\
$\mathfrak{sl}_3$ & $1$ & $0$ & $0$ & $1$ & $0$ & $1$ & $0$ & $0$ & $1$ & $(1+t^3)(1+t^5)$ \\
$\mathfrak{sl}_4$ & $1$ & $0$ & $0$ & $1$ & $0$ & $1$ & $0$ & $1$ & $\cdots$ & $(1+t^3)(1+t^5)(1+t^7)$ \\
$\mathfrak{h}_1$ & $1$ & $2$ & $2$ & $1$ & $0$ & $0$ & $0$ & $0$ & $0$ & $\text{(see text)}$ \\
\bottomrule
\end{tabular}
\end{table}
""")

    print(r"""
% TABLE 3: Graded characters
\begin{table}[ht]
\centering
\caption{Graded Hilbert series $h_B(t,q) = \sum_{n,w} \dim B^{\mathrm{ord,gen}}_{n,w}\, t^n q^w$
at the generator level, with $q$ tracking conformal weight.}
\label{tab:graded-hilbert}
\begin{tabular}{l l l}
\toprule
Family & Generator character $S(q)$ & $h_B(t,q) = \frac{S(q)\,t}{1 - S(q)\,t}$ \\
\midrule
$H_k$ & $q$ & $\frac{tq}{1-tq}$ \\
$V_k(\mathfrak{sl}_2)$ & $3q$ & $\frac{3tq}{1-3tq}$ \\
$V_k(\mathfrak{sl}_N)$ & $(N^2-1)\,q$ & $\frac{(N^2-1)tq}{1-(N^2-1)tq}$ \\
$\mathrm{Vir}_c$ & $q^2$ & $\frac{tq^2}{1-tq^2}$ \\
$\mathcal{W}_3$ & $q^2 + q^3$ & $\frac{t(q^2+q^3)}{1-t(q^2+q^3)}$ \\
$\beta\gamma$ & $1 + q$ & $\frac{t(1+q)}{1-t(1+q)}$ \\
\bottomrule
\end{tabular}
\end{table}
""")

    print(r"""
% TABLE 4: Full Virasoro bar complex dimensions
\begin{table}[ht]
\centering
\caption{Full ordered bar complex $\dim B^{\mathrm{ord,full}}_{n,w}(\mathrm{Vir}_c)$,
including all descendants. The augmentation ideal character is
$\chi_{\mathrm{aug}}(q) = \prod_{n\geq 2}(1-q^n)^{-1} - 1$.}
\label{tab:vir-full-bar}
\begin{tabular}{r r r r r r r}
\toprule
$w$ & $n{=}1$ & $n{=}2$ & $n{=}3$ & $n{=}4$ & $n{=}5$ & $n{=}6$ \\
\midrule""")

    aug_dims = {w: (partitions_min_part(w, 2) if w >= 2 else 0) for w in range(0, 20)}
    for w in range(2, 16):
        dims = []
        for n in range(1, 7):
            d = compute_tensor_dim(aug_dims, n, w)
            dims.append(str(d))
        print(f"{w} & {' & '.join(dims)} \\\\")

    print(r"""\bottomrule
\end{tabular}
\end{table}
""")


# =========================================================================
# SUMMARY
# =========================================================================

def print_summary():
    """Print the master summary."""
    print(f"\n{'=' * 72}")
    print(f"MASTER SUMMARY: ORDERED BAR COMPLEX HILBERT SERIES")
    print(f"{'=' * 72}")

    print("""
+----------------------------------------------------------------------+
|                     GENERATOR-LEVEL BAR COMPLEX                      |
+----------------------------------------------------------------------+
| Family        | g   | dim B_n | h_B(t)       | h_B(t,q)             |
+----------------------------------------------------------------------+
| H_k           |  1  | 1       | t/(1-t)      | tq/(1-tq)            |
| V_k(sl_2)     |  3  | 3^n     | 3t/(1-3t)    | 3tq/(1-3tq)          |
| V_k(sl_3)     |  8  | 8^n     | 8t/(1-8t)    | 8tq/(1-8tq)          |
| V_k(sl_N)     | d   | d^n     | dt/(1-dt)    | dtq/(1-dtq)          |
| Vir_c         |  1  | 1       | t/(1-t)      | tq^2/(1-tq^2)        |
| W_3           |  2  | 2^n     | 2t/(1-2t)    | t(q^2+q^3)/...       |
| beta-gamma    |  2  | 2^n     | 2t/(1-2t)    | t(1+q)/(1-t(1+q))    |
+----------------------------------------------------------------------+
| d = N^2-1 for sl_N; h_B(t,q) = S(q)t/(1-S(q)t) with S(q)=gen char. |
+----------------------------------------------------------------------+

+----------------------------------------------------------------------+
|           CE HOMOLOGY (= DEPTH-0 BAR COHOMOLOGY)                     |
+----------------------------------------------------------------------+
| Lie algebra g | Poincare polynomial            | Total dim  | chi    |
+----------------------------------------------------------------------+
| k (abelian)   | 1 + t                          |   2        |  0     |
| sl_2          | 1 + t^3                        |   2        |  0     |
| sl_3          | (1+t^3)(1+t^5)                 |   4        |  0     |
| sl_4          | (1+t^3)(1+t^5)(1+t^7)          |   8        |  0     |
| sl_N          | prod_{i=1}^{N-1} (1+t^{2i+1}) |  2^{N-1}   |  0     |
| h_1 (Heis)    | 1 + 2t + 2t^2 + t^3            |   6        |  0     |
| Vir (inf-dim) | prod_{n>=1} (1+t^{2n-1})       | infinite   |  --    |
+----------------------------------------------------------------------+
| chi = 0 for all finite-dim semisimple g (Poincare duality).         |
+----------------------------------------------------------------------+

+----------------------------------------------------------------------+
|                FULL BAR COMPLEX (WITH DESCENDANTS)                    |
+----------------------------------------------------------------------+
| Virasoro Vir_c:                                                      |
|   chi_aug(q) = q^2 + q^3 + 2q^4 + 2q^5 + 4q^6 + 4q^7 + 7q^8 +... |
|   h_B^full(t,q) = chi_aug(q)*t / (1 - chi_aug(q)*t)                |
|                                                                      |
| KM V_k(g):                                                          |
|   chi_aug(q) = sum_{d>=1} C(d + dim_g - 1, d) q^d                   |
|   h_B^full(t,q) = chi_aug(q)*t / (1 - chi_aug(q)*t)                |
|                                                                      |
| W_3:                                                                 |
|   chi_aug(q) = prod_{n>=2} 1/(1-q^n) * prod_{n>=3} 1/(1-q^n) - 1   |
+----------------------------------------------------------------------+

KEY FORMULAS:

1. Generator-level: h_B(t) = gt/(1-gt), g = #generators
   General: h_B(t,q) = S(q)t / (1 - S(q)t), S(q) = sum_gen q^{wt(gen)}

2. Full: h_B^full(t,q) = chi_aug(q)*t / (1 - chi_aug(q)*t)

3. CE homology: H^CE_*(sl_N) has Poincare poly = prod (1 + t^{2e_i+1})
   where e_1,...,e_{N-1} are exponents 1,2,...,N-1 of sl_N.

4. Euler characteristics: chi(g) = 0 for all semisimple g (Poincare duality).
   For abelian g: chi = 0 (paired odd/even degrees).

5. The ORDERED bar B^{ord}(A) differs from the SYMMETRIC bar B^{Sigma}(A)
   by the R-matrix descent: B^{Sigma}_n = (B^{ord}_n)^{R-Sigma_n}.
   At the generator level for V_k(g), the descent quotients ordered words
   by R-twisted symmetric group action. The R-matrix is nontrivial
   (R(z) = exp(k*Omega/z) for affine KM, by AP41).
""")


# =========================================================================
# MAIN
# =========================================================================

def main():
    print("ORDERED BAR COMPLEX: FULL DIMENSION SEQUENCES AND HILBERT SERIES")
    print("=" * 72)
    print()

    # 1. Heisenberg
    analyze_heisenberg()

    # 2. sl_2
    analyze_sl2()

    # 3. sl_3, sl_4
    analyze_slN(3, max_ce_p=8)  # sl_3: 8 generators, max C(8,p)=70 at p=4
    analyze_slN(4, max_ce_p=4)

    # 4. Virasoro
    analyze_virasoro()

    # 5. W_3
    analyze_w3()

    # 6. beta-gamma
    analyze_betagamma()

    # Tables and summary
    generate_latex_tables()
    print_summary()


if __name__ == '__main__':
    main()
