"""L^1 norms ||m_k|_T|| for k=2,...,15 using double precision.

For the L^1 norm (average of |T-coefficient| over random spectral parameters),
floating-point cancellation is much less severe than at the symmetric point
because we average over many different parameter configurations.
"""
import sys, os, time, math, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

from compute.m7_m10_depth_frontier import StasheffEngine

c_val = 1.0
n_samples = 2000

print("=" * 95)
print("L^1 NORMS ||m_k|_T|| for k=2,...,15")
print(f"Monte Carlo: {n_samples} samples over [-1,1]^{{k-1}}")
print("=" * 95)

print(f"\n{'k':>3} {'||m_k|_T||':>16} {'k!':>18} {'||/k!':>18} {'succ ratio':>14} {'time':>8}")
print("-" * 85)

norms = {}

for k in range(2, 16):
    engine = StasheffEngine(c_val)
    rng = random.Random(77777 + k)
    total_T = 0.0

    t0 = time.time()
    for trial in range(n_samples):
        engine._cache.clear()
        lams = tuple(rng.uniform(-1.0, 1.0) for _ in range(k - 1))
        result = engine.mk(lams)
        T_coeff = abs(result.get(0, 0.0))
        total_T += T_coeff
    dt = time.time() - t0

    L1 = total_T / n_samples
    norms[k] = L1
    kfact = math.factorial(k)
    gev = L1 / kfact

    if k > 2 and norms.get(k-1, 0) > 1e-300:
        prev_gev = norms[k-1] / math.factorial(k-1)
        ratio = gev / prev_gev if prev_gev > 1e-300 else float('nan')
    else:
        ratio = float('nan')

    ratio_str = f"{ratio:.6f}" if not math.isnan(ratio) else "-"
    print(f"{k:>3} {L1:>16.6e} {kfact:>18} {gev:>18.8e} {ratio_str:>14} {dt:>7.1f}s")

# Gevrey analysis
print("\n" + "=" * 95)
print("GEVREY-1 ANALYSIS")
print("=" * 95)
print()
print("If ||m_k|_T|| ~ C * A^k * k!, then ||m_k|_T||/k! ~ C * A^k")
print("and successive ratios (||m_k||/k!) / (||m_{k-1}||/(k-1)!) -> A.")
print()

print(f"  {'k':>3} {'||m_k|_T||/k!':>18} {'ratio':>14} {'log10(||/k!)':>16}")
print("  " + "-" * 55)

ratios = []
for k in range(2, 16):
    gev = norms[k] / math.factorial(k)
    log_gev = math.log10(gev) if gev > 0 else float('nan')

    if k > 2:
        prev_gev = norms[k-1] / math.factorial(k-1)
        if prev_gev > 1e-300:
            r = gev / prev_gev
            ratios.append((k, r))
            print(f"  {k:>3} {gev:>18.8e} {r:>14.6f} {log_gev:>16.4f}")
        else:
            print(f"  {k:>3} {gev:>18.8e} {'N/A':>14} {log_gev:>16.4f}")
    else:
        print(f"  {k:>3} {gev:>18.8e} {'-':>14} {log_gev:>16.4f}")

if ratios:
    all_r = [r for _, r in ratios]
    late_r = [r for k, r in ratios if k >= 10]
    print(f"\n  Overall average ratio: {sum(all_r)/len(all_r):.6f}")
    if late_r:
        print(f"  Late-stage (k>=10) average: {sum(late_r)/len(late_r):.6f}")
    print(f"  The ratio ||m_k|_T||/k! at k=15: {norms[15]/math.factorial(15):.4e}")
    print(f"  Approximate order: ~10^{math.log10(norms[15]/math.factorial(15)):.1f}")

# Even/odd split
print("\n" + "=" * 95)
print("EVEN/ODD NORM DICHOTOMY")
print("=" * 95)
print()
print("Due to palindromic cancellation, even-k T-coefficients are suppressed.")
print()
print(f"  {'k':>3} {'||m_k|_T||':>16} {'||m_k|_T||/k!':>18} {'parity':>8}")
print("  " + "-" * 50)
for k in range(2, 16):
    gev = norms[k] / math.factorial(k)
    parity = "ODD" if k % 2 == 1 else "even"
    print(f"  {k:>3} {norms[k]:>16.6e} {gev:>18.8e} {parity:>8}")
