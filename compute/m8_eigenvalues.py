r"""Compute exact eigenvalues of the depth-2 quadratic form matrix M_8.

The matrix is:
M = ( 4  -2   0   0   0   0   0 )
    (-2  -4   4   0   0   0   0 )
    ( 0   4   2  -5   0   0   0 )
    ( 0   0  -5   0   5   0   0 )
    ( 0   0   0   5  -2  -4   0 )
    ( 0   0   0   0  -4   4   2 )
    ( 0   0   0   0   0   2  -4 )

The characteristic polynomial is:
λ^7 - 126λ^5 + 4512λ^3 - 46816λ = λ(λ^6 - 126λ^4 + 4512λ^2 - 46816)

So λ=0 is an eigenvalue. Setting μ=λ^2:
μ^3 - 126μ^2 + 4512μ - 46816 = 0

This is a cubic in μ. Let's solve it.
"""
import math


def solve_cubic(a, b, c, d):
    """Solve ax^3 + bx^2 + cx + d = 0 using Cardano's formula."""
    # Normalize to x^3 + px + q = 0 after depressed cubic transform
    b, c, d = b/a, c/a, d/a
    # x = t - b/3
    p = c - b**2/3
    q = 2*b**3/27 - b*c/3 + d

    discriminant = -(4*p**3 + 27*q**2)
    print(f"  Depressed cubic: t^3 + {p:.4f}t + {q:.4f}")
    print(f"  Discriminant: {discriminant:.4f}")

    if discriminant > 0:
        # Three real roots
        m = 2*math.sqrt(-p/3)
        theta = math.acos(3*q/(p*m)) / 3
        roots = [
            m * math.cos(theta - 2*math.pi*k/3) - b/3
            for k in range(3)
        ]
    else:
        # One real root (shouldn't happen here)
        D = -discriminant/108
        sqD = math.sqrt(D) if D >= 0 else 0
        u = (-q/2 + sqD)**(1/3) if -q/2+sqD >= 0 else -(q/2-sqD)**(1/3)
        v = (-q/2 - sqD)**(1/3) if -q/2-sqD >= 0 else -(q/2+sqD)**(1/3)
        roots = [u + v - b/3]
        print(f"  Only one real root: {roots[0]:.6f}")

    return sorted(roots)


def main():
    print("=" * 90)
    print("EIGENVALUES OF THE DEPTH-2 QUADRATIC FORM M_8")
    print("=" * 90)

    print("\nCharacteristic polynomial: λ(λ^6 - 126λ^4 + 4512λ^2 - 46816) = 0")
    print("Setting μ = λ^2: μ^3 - 126μ^2 + 4512μ - 46816 = 0")

    # Solve μ^3 - 126μ^2 + 4512μ - 46816 = 0
    mu_roots = solve_cubic(1, -126, 4512, -46816)

    print(f"\n  μ roots: {[f'{r:.6f}' for r in mu_roots]}")

    all_eigenvalues = [0.0]
    for mu in mu_roots:
        if mu >= 0:
            all_eigenvalues.append(math.sqrt(mu))
            all_eigenvalues.append(-math.sqrt(mu))
        else:
            print(f"  WARNING: negative μ = {mu:.6f} => complex eigenvalues")

    all_eigenvalues.sort()
    print(f"\n  All eigenvalues of M_8:")
    for ev in all_eigenvalues:
        print(f"    λ = {ev:>12.6f}")

    # Verify: product of nonzero eigenvalues
    nonzero_ev = [ev for ev in all_eigenvalues if abs(ev) > 1e-8]
    product = 1.0
    for ev in nonzero_ev:
        product *= ev
    print(f"\n  Product of nonzero eigenvalues: {product:.4f}")
    print(f"  (Should relate to cofactor of zero eigenvalue)")

    # Sum of squares
    sum_sq = sum(ev**2 for ev in all_eigenvalues)
    print(f"  Sum of eigenvalue squares: {sum_sq:.4f} (should be tr(M^2) = 252)")

    # Check: are the μ-roots nice?
    print(f"\n  Checking if μ-roots are rational or have simple form:")
    for mu in mu_roots:
        # Is it an integer?
        if abs(mu - round(mu)) < 1e-4:
            print(f"    μ = {int(round(mu))} (integer)")
        else:
            # Is it a simple fraction?
            for den in range(2, 20):
                if abs(mu * den - round(mu * den)) < 1e-3:
                    print(f"    μ = {int(round(mu*den))}/{den}")
                    break
            else:
                print(f"    μ = {mu:.6f} (not a simple rational)")

    # The anti-palindromic property means eigenvalues come in (λ, -λ) pairs plus 0
    print(f"\n  Anti-palindromic check:")
    print(f"  Eigenvalues should come in ±λ pairs (plus 0):")
    pos = sorted([ev for ev in all_eigenvalues if ev > 1e-8])
    neg = sorted([ev for ev in all_eigenvalues if ev < -1e-8], reverse=True)
    for p, n in zip(pos, neg):
        print(f"    +{p:.6f}, {n:.6f}, sum = {p+n:.6e}")

    # RANK: the matrix has rank 6 (one zero eigenvalue)
    print(f"\n  RANK OF M_8: 6 (one zero eigenvalue)")
    print(f"  The null vector corresponds to the symmetric-point direction")

    # Find the null vector by solving M v = 0
    # M is tridiagonal: a_i v_{i-1} + d_i v_i + a_{i+1} v_{i+1} = 0
    # With a = [-2, 4, -5, 5, -4, 2], d = [4, -4, 2, 0, -2, 4, -4]
    M = [
        [4, -2, 0, 0, 0, 0, 0],
        [-2, -4, 4, 0, 0, 0, 0],
        [0, 4, 2, -5, 0, 0, 0],
        [0, 0, -5, 0, 5, 0, 0],
        [0, 0, 0, 5, -2, -4, 0],
        [0, 0, 0, 0, -4, 4, 2],
        [0, 0, 0, 0, 0, 2, -4],
    ]

    # Forward elimination for null space
    # Row 1: 4v1 - 2v2 = 0 => v2 = 2v1
    # Row 2: -2v1 - 4v2 + 4v3 = 0 => -2v1 - 8v1 + 4v3 = 0 => v3 = 10v1/4 = 5v1/2
    # Row 3: 4v2 + 2v3 - 5v4 = 0 => 8v1 + 5v1 - 5v4 = 0 => v4 = 13v1/5
    # Row 4: -5v3 + 5v5 = 0 => v5 = v3 = 5v1/2
    # Row 5: 5v4 - 2v5 - 4v6 = 0 => 13v1 - 5v1 - 4v6 = 0 => v6 = 8v1/4 = 2v1
    # Row 6: -4v5 + 4v6 + 2v7 = 0 => -10v1 + 8v1 + 2v7 = 0 => v7 = v1
    # Row 7: 2v6 - 4v7 = 0 => 4v1 - 4v1 = 0. CHECK!

    v1 = 10.0  # scale to clear denominators
    v = [v1, 2*v1, 5*v1/2, 13*v1/5, 5*v1/2, 2*v1, v1]
    # With v1=10: [10, 20, 25, 26, 25, 20, 10]
    print(f"\n  Null vector (scaled to integers): {[int(x) for x in v]}")

    # Verify
    print(f"  Verification Mv = 0:")
    for i in range(7):
        val = sum(M[i][j] * v[j] for j in range(7))
        print(f"    Row {i+1}: {val:.4f}")

    # Interpretation
    print(f"\n  INTERPRETATION:")
    print(f"  The null vector [10, 20, 25, 26, 25, 20, 10] is PALINDROMIC!")
    print(f"  It corresponds to the 'balanced' spectral configuration.")
    print(f"  At this ratio, the depth-2 coefficient d^5T vanishes.")
    print(f"  The SYMMETRIC point (1,...,1) is NOT this ratio [10,20,25,26,25,20,10]")
    print(f"  but Q(1,...,1) = 0 by the row-sum condition (not the null vector).")

    # Wait: row sums are NOT all zero. Let me recheck.
    print(f"\n  Row sums of M:")
    for i in range(7):
        rs = sum(M[i][j] for j in range(7))
        print(f"    Row {i+1}: {rs}")

    # Q(1,...,1) = sum of all entries
    total = sum(M[i][j] for i in range(7) for j in range(7))
    print(f"  Q(1,...,1) = sum of all M_{'{ij}'} = {total}")
    print(f"  This is 0, confirming symmetric-point vanishing.")
    print(f"  But this is (1,...,1)^T M (1,...,1) = sum of row sums = {total}")
    row_sums = [sum(M[i][j] for j in range(7)) for i in range(7)]
    print(f"  Individual row sums: {row_sums}")
    print(f"  Sum of row sums: {sum(row_sums)}")

    # (1,...,1) is NOT the null vector. Q(1,...,1) = 0 because
    # sum of row sums = 0, but individual row sums are nonzero.
    # The null VECTOR is [10,20,25,26,25,20,10].

    # At the null vector, Q vanishes because Mv = 0.
    # At (1,...,1), Q vanishes because (1,...,1) is orthogonal to
    # M(1,...,1) = (row sums), and (1,...,1)·(row sums) = 0.
    # Actually Q(1,...,1) = (1,...,1)^T M (1,...,1) = sum of all M entries = 0.

    print(f"\n  KEY DISTINCTION:")
    print(f"  Null vector: Mv = 0 => Q vanishes TO ALL ORDERS in direction v")
    print(f"  Symmetric (1,...,1): Q(1,...,1) = 0 but M(1,...,1) ≠ 0")
    print(f"  The depth-2 vanishing at the symmetric point is a CODIMENSION-1 accident,")
    print(f"  not a rank deficiency. The rank deficiency is along the palindromic")
    print(f"  direction [10,20,25,26,25,20,10].")


if __name__ == '__main__':
    main()
