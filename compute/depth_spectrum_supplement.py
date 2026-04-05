r"""Supplement: L^1 norms for k=12,13 with 500 samples + scalar polynomial.

Runs faster than the 2000-sample main computation.
"""
import sys, os, math, random, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from compute.m7_m10_depth_frontier import StasheffEngine


def main():
    engine = StasheffEngine(1.0)

    # L^1 norms (500 samples, [-1,1]^{k-1})
    print("L^1 norms (500 samples each):")
    print(f"{'k':>3} {'||m_k|_T||':>16} {'||m_k||_all':>16} {'||m_k|_sc||':>16} {'|T|/k!':>16} {'time':>8}")
    print("-" * 85)
    rng = random.Random(99999)
    norms = {}
    for k in range(2, 14):
        t0 = time.time()
        total_T = 0.0
        total_all = 0.0
        total_sc = 0.0
        n_samples = 500
        for _ in range(n_samples):
            engine._cache.clear()
            lams = tuple(rng.uniform(-1.0, 1.0) for _ in range(k - 1))
            result = engine.mk(lams)
            total_T += abs(result.get(0, 0.0))
            total_all += sum(abs(v) for f, v in result.items() if f >= 0)
            total_sc += abs(result.get(-1, 0.0))
        elapsed = time.time() - t0
        L1_T = total_T / n_samples
        L1_all = total_all / n_samples
        L1_sc = total_sc / n_samples
        norms[k] = {'L1_T': L1_T, 'L1_all': L1_all, 'L1_scalar': L1_sc}
        ratio_fac = L1_T / math.factorial(k)
        print(f"{k:>3} {L1_T:>16.8e} {L1_all:>16.8e} {L1_sc:>16.8e} {ratio_fac:>16.10e} {elapsed:>7.1f}s")

    # Gevrey ratios
    print("\nGevrey-1 normalized ratios:")
    for k in range(3, 14):
        r_k = norms[k]['L1_T'] / math.factorial(k)
        r_km1 = norms[k-1]['L1_T'] / math.factorial(k-1)
        if r_km1 > 0:
            ratio = r_k / r_km1
            print(f"  k={k:>2}: {ratio:.8f}")

    # Scalar polynomial
    print("\nScalar polynomial P_k(1,...,1):")
    print(f"{'k':>3} {'P_k':>24} {'P_k/k!':>18}")
    print("-" * 50)
    for k in range(2, 14):
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(k - 1))
        result = engine.mk(lams)
        sc = result.get(-1, 0.0)
        P_k = sc * 12.0  # P_k = sc / (c/12) with c=1
        ratio = P_k / math.factorial(k) if abs(P_k) > 1e-6 else 0.0
        if abs(P_k) < 1e-6:
            print(f"{k:>3} {'0':>24} {'0':>18}")
        else:
            print(f"{k:>3} {P_k:>24.4f} {ratio:>18.10f}")

    # Successive odd ratios
    print("\nOdd P_k ratios:")
    prev_k = None
    prev_P = None
    for k in range(2, 14):
        engine._cache.clear()
        lams = tuple(1.0 for _ in range(k - 1))
        result = engine.mk(lams)
        P_k = result.get(-1, 0.0) * 12.0
        if k % 2 == 1 and abs(P_k) > 1e-6:
            if prev_P is not None:
                print(f"  P_{k}/P_{prev_k} = {P_k/prev_P:.4f}")
            prev_k = k
            prev_P = P_k


if __name__ == '__main__':
    main()
