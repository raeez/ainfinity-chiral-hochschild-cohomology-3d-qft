r"""High-precision Virasoro symmetric-point T_k(1,...,1) using mpmath.

Uses arbitrary-precision arithmetic to eliminate floating-point cancellation errors
that become severe at k >= 14 in double precision.
"""
import sys, os, time, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

import mpmath
mpmath.mp.dps = 50  # 50 decimal digits of precision

mpf = mpmath.mpf
ZERO = mpf(0)
ONE = mpf(1)
TWO = mpf(2)
TWELVE = mpf(12)

MAX_DERIV = 30

# ===== Field-dict operations with mpmath =====

def fd_add(*dicts, signs=None):
    if signs is None:
        signs = [ONE] * len(dicts)
    result = {}
    for d, s in zip(dicts, signs):
        for f, c in d.items():
            result[f] = result.get(f, ZERO) + s * c
    return {k: v for k, v in result.items() if abs(v) > mpf('1e-45')}

def fd_scale(d, factor):
    return {f: factor * c for f, c in d.items() if abs(factor * c) > mpf('1e-45')}

def fd_apply_partial(d):
    result = {}
    for f, c in d.items():
        if f == -1: continue
        if f < 0: continue
        if f + 1 <= MAX_DERIV:
            result[f + 1] = result.get(f + 1, ZERO) + c
    return {k: v for k, v in result.items() if abs(v) > mpf('1e-45')}

def fd_apply_shift_partial(d, shift):
    shifted = fd_scale(d, shift)
    partial = fd_apply_partial(d)
    return fd_add(shifted, partial)

def fd_apply_shift_partial_n(d, shift, n):
    result = dict(d)
    for _ in range(n):
        result = fd_apply_shift_partial(result, shift)
    return result

# ===== Core operations =====

def m2_num(lam, c):
    return {1: ONE, 0: TWO * lam, -1: c * lam**3 / TWELVE}

def m3_num(l1, l2, c):
    return {
        2: ONE,
        1: TWO * l1 + mpf(3) * l2,
        0: mpf(4) * l1 * l2 + TWO * l2**2,
        -1: c * l2**3 * (TWO * l1 + l2) / TWELVE,
    }

def compose_into_mk_slot(mk_func, inner, slot, outer_lams, n_slots, c):
    base = mk_func(*outer_lams, c)
    result = {}
    for field, coeff in inner.items():
        if field == -1: continue
        n = field
        if n < 0: continue
        if slot < n_slots - 1:
            slot_lam = outer_lams[slot]
            factor = (-slot_lam)**n
            for bf, bc in base.items():
                result[bf] = result.get(bf, ZERO) + coeff * factor * bc
        else:
            shift = sum(outer_lams)
            shifted_base = fd_apply_shift_partial_n(base, shift, n)
            for bf, bc in shifted_base.items():
                result[bf] = result.get(bf, ZERO) + coeff * bc
    return {k: v for k, v in result.items() if abs(v) > mpf('1e-45')}


class MPStasheffEngine:
    """High-precision Stasheff recursion engine using mpmath."""

    def __init__(self, c_val):
        self.c = mpf(c_val)
        self._cache = {}

    def mk(self, lams):
        k = len(lams) + 1
        if k < 2:
            raise ValueError(f"Arity must be >= 2, got {k}")

        cache_key = lams
        if cache_key in self._cache:
            return self._cache[cache_key]

        if k == 2:
            result = m2_num(lams[0], self.c)
        elif k == 3:
            result = m3_num(lams[0], lams[1], self.c)
        else:
            result = self._stasheff_rhs(k, lams)
            result = fd_scale(result, mpf(-1))

        self._cache[cache_key] = result
        return result

    def mk_func_factory(self, arity):
        def func(*args):
            lams = args[:-1]
            return self.mk(tuple(lams))
        return func

    def _stasheff_rhs(self, k, lams):
        total = {}
        lam_list = list(lams)

        for j in range(2, k):
            i = k + 1 - j
            if i < 2: continue

            for s in range(k - j + 1):
                inner_lams = tuple(lam_list[s:s+j-1])
                inner_result = self.mk(inner_lams)
                outer_lams = self._compute_outer_lams(k, lam_list, s, j)

                comp = compose_into_mk_slot(
                    self.mk_func_factory(i),
                    inner_result, s, list(outer_lams), i, self.c
                )

                sign = mpf(-1)**s
                for f, v in comp.items():
                    total[f] = total.get(f, ZERO) + sign * v

        return {k_: v for k_, v in total.items() if abs(v) > mpf('1e-45')}

    def _compute_outer_lams(self, k, lam_list, s, j):
        n = k
        outer_params = []

        for p in range(s):
            if p < len(lam_list):
                outer_params.append(lam_list[p])

        block_end = s + j - 1
        if block_end < n - 1:
            merged = sum(lam_list[s:s+j], ZERO)
            outer_params.append(merged)

        for p in range(s + j, len(lam_list)):
            outer_params.append(lam_list[p])

        expected = k - j
        if len(outer_params) != expected:
            raise ValueError(
                f"Outer params mismatch: got {len(outer_params)}, "
                f"expected {expected} (k={k}, j={j}, s={s})"
            )

        return tuple(outer_params)


def catalan(n):
    if n < 0: return 0
    return math.comb(2*n, n) // (n+1)


def main():
    c_val = 1
    engine = MPStasheffEngine(c_val)
    max_k = 25

    print("=" * 120)
    print("VIRASORO SYMMETRIC-POINT T-COEFFICIENT T_k(1,...,1)  [50-digit precision]")
    print("=" * 120)

    print(f"\n{'k':>3} {'r':>3} {'T_k numerical':>35} {'Catalan formula':>35} {'match':>6} {'dt':>8} {'cache':>8}")
    print("-" * 105)

    T_values = {}
    P_values = {}
    all_T_match = True

    for k in range(2, max_k + 1):
        lams = tuple(ONE for _ in range(k - 1))
        t0 = time.time()
        result = engine.mk(lams)
        dt = time.time() - t0

        T_num = result.get(0, ZERO)
        scalar = result.get(-1, ZERO)
        P_num = scalar * TWELVE / mpf(c_val) if abs(scalar) > mpf('1e-45') else ZERO

        T_values[k] = T_num
        P_values[k] = P_num

        # Catalan formula
        if k == 2:
            T_cat = mpf(2)
        elif k % 2 == 0:
            T_cat = ZERO
        else:
            r = (k - 1) // 2
            T_cat = mpf((-1)**(r+1)) * mpf(k) * mpf(catalan(r-1)) * mpf(math.factorial(2*r))

        if abs(T_cat) > mpf('1e-10'):
            rel_err = abs(T_num - T_cat) / abs(T_cat)
            match = rel_err < mpf('1e-30')
        else:
            rel_err = abs(T_num)
            match = rel_err < mpf('1e-30')

        if not match:
            all_T_match = False

        r_str = str((k-1)//2) if k % 2 == 1 else "-"

        def fmt(x):
            x = float(x)
            if abs(x) < 1e-30:
                return "0"
            if abs(x) < 1e10:
                return f"{x:.1f}"
            return f"{x:.6e}"

        match_str = "OK" if match else "FAIL"
        print(f"{k:>3} {r_str:>3} {fmt(T_num):>35} {fmt(T_cat):>35} {match_str:>6} {dt:>7.2f}s {len(engine._cache):>8}")

    print(f"\nCATALAN THEOREM VERIFIED through k={max_k}: {'YES -- ALL MATCH' if all_T_match else 'NO'}")

    # ========== Even-k vanishing ==========
    print("\n" + "-" * 80)
    print("EVEN-k VANISHING DETAIL:")
    for k in range(4, max_k + 1, 2):
        T = T_values[k]
        print(f"  k={k:>2}: T = {mpmath.nstr(T, 6, strip_zeros=False)}")

    # ========== Scalar P_k ==========
    print("\n" + "=" * 120)
    print("SCALAR P_k(1,...,1)")
    print("=" * 120)
    print(f"\n{'k':>3} {'r':>3} {'P_k numerical':>35} {'formula':>35} {'match':>6}")
    print("-" * 90)

    all_P_match = True
    for k in range(2, max_k + 1):
        P_num = P_values[k]

        if k == 2:
            P_form = ONE
        elif k % 2 == 0:
            P_form = ZERO
        else:
            r = (k - 1) // 2
            P_form = mpf((-1)**(r-1)) * mpf(catalan(r-1)) * mpf(math.factorial(2*r+1)) / TWO

        if abs(P_form) > mpf('1e-10'):
            rel = abs(P_num - P_form) / abs(P_form)
            match = rel < mpf('1e-20')
        else:
            match = abs(P_num) < mpf('1e-20')

        if not match:
            all_P_match = False

        def fmt(x):
            x = float(x)
            if abs(x) < 1e-20:
                return "0"
            if abs(x) < 1e10:
                return f"{x:.1f}"
            return f"{x:.6e}"

        r_str = str((k-1)//2) if k % 2 == 1 else "-"
        print(f"{k:>3} {r_str:>3} {fmt(P_num):>35} {fmt(P_form):>35} {'OK' if match else 'FAIL':>6}")

    print(f"\nSCALAR FORMULA VERIFIED: {'YES' if all_P_match else 'NEEDS INVESTIGATION'}")

    # T/P ratio
    print("\n  T_k / P_k for odd k:")
    for k in range(3, max_k + 1, 2):
        if abs(P_values[k]) > mpf('1e-10'):
            ratio = T_values[k] / P_values[k]
            expected = TWO / mpf(k + 1)
            print(f"  k={k:>2}: T/P = {float(ratio):.10f},  2/(k+1) = {float(expected):.10f},  match: {abs(ratio - expected) < mpf('1e-20')}")

    # ========== Reference table ==========
    print("\n" + "=" * 120)
    print("REFERENCE: Catalan formula components")
    print("=" * 120)
    print(f"\n{'k':>3} {'r':>3} {'C_{r-1}':>12} {'(2r)!':>22} {'(-1)^{r+1}':>12} {'T_k formula':>30}")
    print("-" * 90)
    for k in range(3, max_k + 1, 2):
        r = (k - 1) // 2
        C = catalan(r - 1)
        sgn = (-1)**(r+1)
        fac = math.factorial(2*r)
        val = sgn * k * C * fac
        print(f"{k:>3} {r:>3} {C:>12} {fac:>22} {sgn:>12} {val:>30}")

    print("\n" + "=" * 120)
    print("TOTAL COMPUTATION TIME: will print after L^1 norms")
    print("=" * 120)


if __name__ == '__main__':
    t_global = time.time()
    main()
    print(f"\n*** Total wall time: {time.time() - t_global:.1f}s ***")
