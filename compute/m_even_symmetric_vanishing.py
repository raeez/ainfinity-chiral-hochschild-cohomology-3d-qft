r"""Investigate the symmetric-point vanishing pattern for even arities.

At the symmetric point l_1 = l_2 = ... = l_{k-1} = lambda:
  m_2: NOT zero (just m_2(T,T;lambda) = dT + 2*lambda*T + c/12*lambda^3)
  m_4: ZERO for ALL lambda (confirmed numerically and symbolically)
  m_6: T-sector ZERO to ~1e-6, scalar near zero. Is this genuine?
  m_8: T-sector ZERO to ~1e-4, scalar nonzero
  m_10: T-sector nonzero at large lambda

The question: is the T-sector of m_{2r}(T^{2r}; lambda,...,lambda)
identically zero for all r >= 2?

The depth 0 vanishing (d^{k-2}T coefficient = 0) for even k is CLEAN:
it's a genuine structural vanishing. But the full symmetric-point behavior
requires further investigation.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from compute.m7_m10_depth_frontier import StasheffEngine


def symmetric_scan(c_val, k, lam_values):
    """Scan m_k at the symmetric point for many lambda values."""
    engine = StasheffEngine(c_val)

    print(f"\n  m_{k} at c={c_val}, symmetric point l_i = lambda:")
    print(f"  {'lambda':>10} {'|T-sector|':>14} {'d0':>14} {'d1':>14} {'d2':>14} {'scalar':>14}")

    for lam in lam_values:
        engine._cache.clear()
        lams = tuple(lam for _ in range(k - 1))
        result = engine.mk(lams)

        T_total = sum(abs(result.get(d, 0.0)) for d in range(k))
        d0 = abs(result.get(k - 2, 0.0))  # depth 0 = deriv order k-2
        d1 = abs(result.get(k - 3, 0.0))  # depth 1 = deriv order k-3
        d2 = abs(result.get(k - 4, 0.0)) if k >= 4 else 0.0
        sc = abs(result.get(-1, 0.0))

        print(f"  {lam:>10.4f} {T_total:>14.6e} {d0:>14.6e} {d1:>14.6e} {d2:>14.6e} {sc:>14.6e}")


def investigate_depth_0_mechanism():
    """Investigate WHY depth 0 vanishes for even k.

    Depth 0 = the coefficient of d^{k-2}T.
    This is the LEADING field term (highest derivative).

    For k=4: depth 0 would be d^2 T coefficient.
    The known m_4 has d^2 T = -(A1+A2+B1+B2+B3)|_{d3T} from the Stasheff RHS.
    At k=4, this requires a d^3T term in the RHS, which after h gives d^2T.

    Actually for the DIRECT formula (no homotopy): the d^2T term of m_4 = -RHS|_{d^2T}.
    This vanishes because of cancellation among the 5 Stasheff compositions.
    Is this the palindrome mechanism?

    Test: at even k, does the d^{k-2}T coefficient have a structural reason to vanish,
    or is it specific to Virasoro?
    """
    print("\n--- Depth 0 vanishing mechanism investigation ---")
    print("  At each even arity, d^{k-2}T coefficient of m_k should be the")
    print("  sum of d^{k-2}T contributions from all Stasheff compositions.")
    print("  If this is a STRUCTURAL vanishing (not accident), it should hold")
    print("  for ALL c values, which we already verified (c-independence).")
    print()
    print("  The depth-0 coefficient has weight w=k-2, so its spectral-param")
    print("  coefficient is LINEAR in the l_i (depth d=1, but wait d=k-1-w=1...")
    print("  No: depth 0 has w=k-2, d=k-1-w=1. That contradicts! Let me recheck.")
    print()
    print("  Actually: the FIELD d^wT has weight w. Its coefficient has spectral degree d.")
    print("  The weight-depth identity says w + d = k-1.")
    print("  So d^{k-2}T has w=k-2, d=k-1-(k-2)=1.")
    print("  Wait, but we're calling this 'depth 0' -- there's a convention issue.")
    print()
    print("  RESOLUTION: The 'depth' in our numerical code is d = k-1-w where w is")
    print("  the derivative order. For depth 0: d=0 means w=k-1, which would be")
    print("  d^{k-1}T. But in the depth spectrum table, the LEADING field has")
    print("  the HIGHEST derivative order (= deepest penetration into the tower).")
    print()
    print("  Looking at the data:")
    print("  m_4: max fields are d^2T and d^1T. d^2T has w=2, d=k-1-2=1. d^1T has w=1, d=2.")
    print("  So Spec(m_4|_T) = {1, 2} (depth 1 and 2), mapping to d^2T and dT.")
    print("  But our table says Spec(m_4|_T) = {2, 3}.")
    print()
    print("  I think there's a depth convention issue. Let me check the manuscript.")
    print("  From the theorem: 'Every T-dependent monomial in m_n satisfies w + d = n - 1'")
    print("  At n=4, m_4 has d^2T (w=2, so d=1) and dT (w=1, so d=2).")
    print("  But the manuscript table says Spec(m_4|_T) = {2, 3}.")
    print()
    print("  Ah! The 'depth d' in the manuscript is the SPECTRAL DEGREE of the coefficient,")
    print("  not what I'm computing. The spectral degree of the coefficient POLYNOMIAL.")
    print("  d^2T has coefficient (l_2 + 2*l_3) which is degree 1 in the l's.")
    print("  Wait, that would be depth 1, not depth 2.")
    print()
    print("  Let me look at what the m6_depth_spectrum code actually does.")


def recheck_depth_convention():
    """Recheck what 'depth' means in the codebase.

    From m6_depth_spectrum.py depth_spectrum():
    p = Poly(coeff, *lam_symbols)
    degs = set()
    for monom, _ in p.as_dict().items():
        degs.add(sum(monom))

    So depth = total degree of the coefficient POLYNOMIAL in the spectral params.
    This is the spectral degree.

    For m_4(T,T,T,T; l1,l2,l3):
    - d^2T coefficient: some polynomial in l1,l2,l3 of certain degree
    - dT coefficient: polynomial of certain degree
    - T coefficient: polynomial of certain degree

    The weight-depth identity says the coefficient of d^wT has spectral degree
    d = k-1-w. So:
    d^2T: w=2, d = 4-1-2 = 1 -> coeff is degree 1 in l's
    dT: w=1, d = 4-1-1 = 2 -> coeff is degree 2
    T: w=0, d = 4-1-0 = 3 -> coeff is degree 3
    scalar: d = k+1 = 5

    But our numerical code uses derivative_order as the key, with:
    depth = k - 1 - derivative_order
    For d^2T: derivative_order = 2, depth = 4-1-2 = 1
    For dT: derivative_order = 1, depth = 4-1-1 = 2
    For T: derivative_order = 0, depth = 4-1-0 = 3

    So in our code, the depths should be {1, 2, 3} for m_4.
    But the manuscript says m_4|_T has Spec = {2, 3}, which means
    depth 1 (= d^2T) vanishes.

    And our earlier numerical output showed:
    m_4: Spec(m_k|_T) = [2, 3]  -> depths 2 and 3 populated
    This means d^2T coefficient (depth 1) VANISHES, and dT and T are nonzero.

    But wait, our precision check said:
    m_4: depth_0 = ZERO
    at derivative order k-2 = 2, which IS d^2T.
    So "depth 0" in our precision check is wrong -- it should be depth 1.
    The bug: in precise_depth_check, I'm mapping d (depth) correctly as k-1-deriv_order.

    Actually wait. In the data table output:
    Row k=4, column d=0: ZERO, d=1: ZERO, d=2: 3.27e+01, d=3: 1.64e+02

    These columns are labeled d=0, d=1, etc. But what does d mean?
    Looking at precise_depth_check: I store depth = k-1-deriv_order.
    For k=4: deriv_order=2 -> depth=1, deriv_order=1 -> depth=2, deriv_order=0 -> depth=3.
    So the data at depth=0 and depth=1:
    depth=0 corresponds to deriv_order=3 = d^3T. d^3T is NOT a valid field for m_4
    (max derivative is k-2=2). So depth 0 should indeed be "absent" structurally.
    And depth=1 corresponds to d^2T, which we see is zero for m_4.

    But the initial summary from the first script showed m_4: Spec = [2, 3], d_min=2.
    The "SECONDARY VANISHING at depths: [1]" message means depth 1 vanishes (= d^2T = 0).

    OK so the picture is now clear:
    - Depth 0 for m_k requires d^{k-1}T. But the Stasheff recursion starting from
      m_2 (which has max derivative 1) and m_3 (max derivative 2) can produce at
      most d^{k-2}T at arity k. So depth 0 is STRUCTURALLY impossible for k >= 3.
      (Actually depth 0 exists for k=2 via dT, and for k=3 via d^2T.)
      Wait: for k=5, the output shows depth 0 is present. Let me re-examine.

      m_5 has fields up to d^3T (k-2=3). d^3T has depth = k-1-3 = 1.
      Then depth 0 = d^4T (k-1-0 = 4)? No, k-1=4, so depth 0 would need w=4,
      but max w for m_5 is k-2=3. So depth 0 requires d^{k-1}T = d^4T which
      should not be present.

      But the output says m_5 has depth 0! Let me recheck.
    """
    print("\n--- Depth convention recheck ---")

    engine = StasheffEngine(1.0)

    for k in range(2, 11):
        engine._cache.clear()
        lams = tuple(1.0 + 0.1 * i for i in range(k - 1))
        result = engine.mk(lams)

        print(f"\n  m_{k} at lams = {lams[:min(5, len(lams))]}{'...' if len(lams) > 5 else ''}:")
        for d_order in sorted([d for d in result if d >= 0], reverse=True):
            depth = k - 1 - d_order
            print(f"    d^{d_order}T (weight={d_order}, depth={depth}): coeff = {result[d_order]:.6e}")
        sc = result.get(-1, 0.0)
        if abs(sc) > 1e-14:
            print(f"    scalar (depth={k+1}): coeff = {sc:.6e}")


def main():
    print("=" * 80)
    print("SYMMETRIC POINT AND DEPTH CONVENTION INVESTIGATION")
    print("=" * 80)

    # 1. Recheck depth convention
    recheck_depth_convention()

    # 2. Symmetric point scan for even arities
    lam_values = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 0.001]
    for k in [4, 6, 8]:
        symmetric_scan(1.0, k, lam_values)

    # 3. Verify even/odd pattern at the correct depth convention
    print(f"\n\n{'=' * 80}")
    print("CORRECTED DEPTH SPECTRUM WITH PROPER CONVENTION")
    print(f"{'=' * 80}")
    print()
    print("Convention: depth d = spectral degree of coefficient polynomial")
    print("           = k - 1 - w where w = derivative order")
    print("           d^wT has depth d = k-1-w")
    print()

    import random
    rng = random.Random(54321)

    for k in range(2, 11):
        engine = StasheffEngine(1.0)
        # Sample many points and track which depths are populated
        depths_present = set()
        for _ in range(200):
            lams = tuple(rng.uniform(-3.0, 3.0) for _ in range(k - 1))
            result = engine.mk(lams)
            for d_order, coeff in result.items():
                if d_order == -1:
                    if abs(coeff) > 1e-8:
                        depths_present.add(k + 1)
                else:
                    if abs(coeff) > 1e-8:
                        depth = k - 1 - d_order
                        depths_present.add(depth)

        # What's the maximum derivative order in the result?
        max_d_order = k - 2  # Expected: max field is d^{k-2}T

        parity = "even" if k % 2 == 0 else "odd"
        depths_T = sorted(d for d in depths_present if d <= k - 1)
        scalar = (k + 1) in depths_present
        gap_at_k = k not in depths_present

        # Theoretical depths: d^{k-2}T (depth 1), d^{k-3}T (depth 2), ..., T (depth k-1)
        # Plus possible d^{k-1}T (depth 0) if the recursion generates it
        # Actually: m_k can have fields d^0T through d^{k-2}T
        # giving depths k-1 down to 1.
        # BUT wait: looking at m_2, m_3:
        # m_2 = dT + 2*l*T + c/12*l^3: fields dT (depth 0) and T (depth 1)
        # m_3 = d^2T + ...*dT + ...*T + scalar: fields d^2T (depth 0), dT (depth 1), T (depth 2)
        # So depth 0 = d^{k-2}T for k=2,3. These exist.
        # For k=4: d^2T is the highest field (k-2=2), depth = 4-1-2 = 1.
        # So depth 0 would be d^3T, which is NOT present in m_4.

        # Hmm, but looking at m_5: the code said Spec includes depth 0.
        # m_5 has d^3T as highest field (k-2=3), depth = 5-1-3 = 1.
        # So depth 0 = d^4T... which m_5 should NOT have.

        # Unless the recursion generates d^4T through high-derivative compositions?
        # Let me check directly.

        print(f"  m_{k:2d} ({parity:4s}): Spec|_T = {depths_T}, "
              f"gap at k={k}: {gap_at_k}, scalar: {'d=' + str(k+1) if scalar else '-'}")


if __name__ == '__main__':
    main()
