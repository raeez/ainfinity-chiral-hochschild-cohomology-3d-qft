r"""Extend Virasoro A∞ depth spectra to arities 7-10 via numerical Stasheff recursion.

Strategy: Full symbolic computation at arity 7 is prohibitively expensive (tracking
d5T and higher, with hundreds of composition terms). Instead, we use a NUMERICAL
engine at specific c values and many random spectral parameter samples to:

  (a) Extract depth spectra (which derivative orders are populated)
  (b) Compute scalar (shadow) coefficients S_k
  (c) Compute L^1 norms for Gevrey-1 analysis
  (d) Test for secondary vanishing (n=4 anomaly pattern at n=8?)
  (e) Test palindrome factorization patterns

The numerical engine represents m_k(T,...,T; λ₁,...,λ_{k-1}) as a dictionary
{derivative_order: polynomial_coefficient} where the coefficient is a NUMBER
(for fixed spectral parameters and c).

Depth spectrum extraction: run at many random λ values and record which
derivative orders are nonzero. A depth d is populated if the coefficient of
∂^{k-2-d}T (weight w = k-2-d, depth d = k-1-w) is generically nonzero.

By the weight-depth identity: w + d = k-1, so field ∂^w T has depth d = k-1-w.
The scalar contact term has depth k+1.
"""

from __future__ import annotations
import sys
import os
import time
import random
import math
from typing import Dict, List, Tuple, Optional

# Ensure compute/ is importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ===== Numerical field-coefficient representation =====
# A "field dict" maps derivative_order -> coefficient (float)
# derivative_order = -1 means scalar (the '1' field)
# derivative_order = 0 means T
# derivative_order = 1 means dT
# derivative_order = n means d^n T

MAX_DERIV = 10  # track up to d^10 T

def fd_add(*dicts: Dict[int, float], signs=None) -> Dict[int, float]:
    """Add field-coeff dicts with optional signs."""
    if signs is None:
        signs = [1.0] * len(dicts)
    result: Dict[int, float] = {}
    for d, s in zip(dicts, signs):
        for f, c in d.items():
            result[f] = result.get(f, 0.0) + s * c
    return {k: v for k, v in result.items() if abs(v) > 1e-300}


def fd_scale(d: Dict[int, float], factor: float) -> Dict[int, float]:
    """Scale all coefficients."""
    return {f: factor * c for f, c in d.items() if abs(factor * c) > 1e-300}


def fd_apply_partial(d: Dict[int, float]) -> Dict[int, float]:
    """Apply ∂ to a field-coeff dict: ∂(d^n T) = d^{n+1} T, ∂(scalar) = 0."""
    result: Dict[int, float] = {}
    for f, c in d.items():
        if f == -1:  # scalar
            continue
        if f + 1 <= MAX_DERIV:
            result[f + 1] = result.get(f + 1, 0.0) + c
    return {k: v for k, v in result.items() if abs(v) > 1e-300}


def fd_apply_shift_partial(d: Dict[int, float], shift: float) -> Dict[int, float]:
    """Apply (shift + ∂) to a field-coeff dict."""
    shifted = fd_scale(d, shift)
    partial = fd_apply_partial(d)
    return fd_add(shifted, partial)


def fd_apply_shift_partial_n(d: Dict[int, float], shift: float, n: int) -> Dict[int, float]:
    """Apply (shift + ∂)^n."""
    result = dict(d)
    for _ in range(n):
        result = fd_apply_shift_partial(result, shift)
    return result


# ===== Core operations (numerical) =====

def m2_num(lam: float, c: float) -> Dict[int, float]:
    """m₂(T,T; λ) = dT + 2Tλ + (c/12)λ³"""
    return {1: 1.0, 0: 2.0 * lam, -1: c * lam**3 / 12.0}


def m3_num(l1: float, l2: float, c: float) -> Dict[int, float]:
    """m₃(T,T,T; l₁, l₂)"""
    return {
        2: 1.0,
        1: 2.0 * l1 + 3.0 * l2,
        0: 4.0 * l1 * l2 + 2.0 * l2**2,
        -1: c * l2**3 * (2.0 * l1 + l2) / 12.0,
    }


# ===== Sesquilinearity composition engine =====

def compose_into_m2_left_num(inner: Dict[int, float], lam_outer: float, c: float) -> Dict[int, float]:
    """m₂(inner, T; λ_outer). Left sesquilinearity: d^n T -> (-λ)^n."""
    base = m2_num(lam_outer, c)
    result: Dict[int, float] = {}
    for field, coeff in inner.items():
        if field == -1:  # scalar maps to 0
            continue
        n = field  # derivative order
        if n < 0:
            continue
        factor = (-lam_outer)**n
        for bf, bc in base.items():
            result[bf] = result.get(bf, 0.0) + coeff * factor * bc
    return {k: v for k, v in result.items() if abs(v) > 1e-300}


def compose_into_m2_right_num(inner: Dict[int, float], lam_outer: float, c: float) -> Dict[int, float]:
    """m₂(T, inner; λ_outer). Right sesquilinearity: d^n T -> (λ+∂)^n."""
    base = m2_num(lam_outer, c)
    result: Dict[int, float] = {}
    for field, coeff in inner.items():
        if field == -1:
            continue
        n = field
        if n < 0:
            continue
        shifted_base = fd_apply_shift_partial_n(base, lam_outer, n)
        for bf, bc in shifted_base.items():
            result[bf] = result.get(bf, 0.0) + coeff * bc
    return {k: v for k, v in result.items() if abs(v) > 1e-300}


def compose_into_mk_slot_num(mk_func, inner: Dict[int, float], slot: int,
                              outer_lams: List[float], n_slots: int, c: float) -> Dict[int, float]:
    """Compose inner output into slot `slot` of m_k evaluated at outer_lams.

    Sesquilinearity rules:
      slot < n_slots-1: d^n T -> (-outer_lams[slot])^n * m_k(T,...,T)
      slot = n_slots-1: d^n T -> (sum(outer_lams) + ∂)^n * m_k(T,...,T)
    """
    base = mk_func(*outer_lams, c)
    result: Dict[int, float] = {}

    for field, coeff in inner.items():
        if field == -1:
            continue
        n = field
        if n < 0:
            continue

        if slot < n_slots - 1:
            slot_lam = outer_lams[slot]
            factor = (-slot_lam)**n
            for bf, bc in base.items():
                result[bf] = result.get(bf, 0.0) + coeff * factor * bc
        else:
            shift = sum(outer_lams)
            shifted_base = fd_apply_shift_partial_n(base, shift, n)
            for bf, bc in shifted_base.items():
                result[bf] = result.get(bf, 0.0) + coeff * bc

    return {k: v for k, v in result.items() if abs(v) > 1e-300}


# ===== Memoized recursive m_k computation =====

class StasheffEngine:
    """Numerical Stasheff recursion engine.

    Computes m_k(T,...,T; λ₁,...,λ_{k-1}) for Virasoro at fixed c.

    Uses memoization keyed by (arity, spectral_params_tuple).
    """

    def __init__(self, c_val: float):
        self.c = c_val
        self._cache: Dict[Tuple, Dict[int, float]] = {}

    def mk(self, lams: Tuple[float, ...]) -> Dict[int, float]:
        """Compute m_k(T,...,T; lams) where k = len(lams) + 1."""
        k = len(lams) + 1
        if k < 2:
            raise ValueError(f"Arity must be >= 2, got {k}")

        # Check cache
        cache_key = lams
        if cache_key in self._cache:
            return self._cache[cache_key]

        if k == 2:
            result = m2_num(lams[0], self.c)
        elif k == 3:
            result = m3_num(lams[0], lams[1], self.c)
        else:
            # General Stasheff recursion:
            # m_k = - sum_{i+j=k+1, i,j>=2} sum_s (-1)^s m_i(..., m_j(...), ...)
            result = self._stasheff_rhs(k, lams)
            result = fd_scale(result, -1.0)

        self._cache[cache_key] = result
        return result

    def mk_func_factory(self, arity: int):
        """Return a function (*lams, c) -> field_dict for compose_into_mk_slot_num."""
        def func(*args):
            # Last arg is c (ignored, we use self.c)
            lams = args[:-1]
            return self.mk(tuple(lams))
        return func

    def _stasheff_rhs(self, k: int, lams: Tuple[float, ...]) -> Dict[int, float]:
        """Compute the RHS of the arity-k Stasheff identity.

        RHS = sum over all (i,j) with i+j=k+1, i,j>=2:
              sum over s=0,...,k-j:
                (-1)^s * m_i(..., m_j(inputs[s:s+j]; inner_lams), ...; outer_lams)
        """
        total: Dict[int, float] = {}
        lam_list = list(lams)  # k-1 spectral params

        # For each partition i+j = k+1 with i,j >= 2
        for j in range(2, k):  # inner arity j
            i = k + 1 - j      # outer arity i
            if i < 2:
                continue

            # For each insertion position s = 0, ..., k-j
            for s in range(k - j + 1):
                # Inner m_j takes inputs at positions s, s+1, ..., s+j-1
                # Inner spectral params: lams[s], ..., lams[s+j-2] (j-1 params)
                inner_lams = tuple(lam_list[s:s+j-1])

                # Compute inner result
                inner_result = self.mk(inner_lams)

                # Compute outer spectral params
                # The outer m_i has i inputs and i-1 spectral params
                # Outer params: original params with inner block contracted
                outer_lams = self._compute_outer_lams(k, lam_list, s, j)

                # Compose inner output into slot s of outer m_i
                inner_slot = s  # which slot of the outer the inner output goes into
                comp = compose_into_mk_slot_num(
                    self.mk_func_factory(i),
                    inner_result, inner_slot, list(outer_lams), i, self.c
                )

                sign = (-1.0)**s
                for f, v in comp.items():
                    total[f] = total.get(f, 0.0) + sign * v

        return {k_: v for k_, v in total.items() if abs(v) > 1e-300}

    def _compute_outer_lams(self, k: int, lam_list: List[float],
                             s: int, j: int) -> Tuple[float, ...]:
        """Compute the outer spectral parameters after inner m_j insertion at position s.

        The inner block takes positions [s, s+j-1].
        Inner params: lam_list[s:s+j-1] (between inputs s and s+j-1).
        After contraction: the inner block becomes a single node.

        The outer params are formed by:
        - params before the block: lam_list[0:s] (unchanged, s params)
          Wait, not s params. The params between positions 0..s-1 are
          lam_list[0:s-1], giving s-1 params.
          BUT position s-1 is adjacent to the block start, so lam_list[s-1]
          is the param between position s-1 and position s (= start of block).
          After contraction, this becomes the param between position s-1 and
          the inner output. So it survives as-is.

        Actually, the convention from the arity-5/6 code is:
        - The gap between the inner output and the next input after the block
          absorbs the sum of the inner params plus the boundary param.

        Let me follow the pattern from the existing code:
        - For inner at position s taking j inputs with inner params lam_list[s:s+j-1]:
          - Params before: lam_list[0:s]  (these are params between positions 0..s-1 and s)
            Actually wait: lam_list[0] is between position 0 and 1, etc.
            lam_list[p] is between position p and p+1.
          - The inner block spans positions s to s+j-1.
          - Internal params (consumed): lam_list[s] to lam_list[s+j-2]
            (params between consecutive inputs within the block)
          - After contraction: the block becomes a single position.
          - The params in the merged sequence:
            * lam_list[0] through lam_list[s-2]: between original positions before block (s-1 params)
            * If s > 0: lam_list[s-1] is between position s-1 and the block start
              After contraction: between position s-1 and the inner output. Survives.
            * The MERGED param between inner output and next position after block:
              Sum of consumed params + boundary.
              From the arity-5 code pattern:
              C0: inner at pos 0-2, outer params [l1+l2+l3, l4]
                merged param = l1+l2+l3 = sum(lam_list[0:3])
              C1: inner at pos 1-3, outer params [l1, l2+l3+l4]
                merged param = l2+l3+l4 = sum(lam_list[1:4])

              So the merged param at the position of the inner output is:
              sum(lam_list[s:s+j-1]) = sum of all params from start of block to
              one past the end of the block.
              Wait: for C0 with j=3 inner, sum(lam_list[0:3]) = l1+l2+l3,
              but lam_list[0:3] = [l1,l2,l3]. The inner only consumed l1,l2
              (j-1=2 params). The extra l3 is the boundary param after the block.

              So: merged param = sum(lam_list[s:s+j-1]) when s+j-1 < k-1
                  (i.e., when the block doesn't end at the rightmost position).
              When the block IS rightmost: no param after the block. The merged
              param is just sum(lam_list[s:s+j-1]).

              Actually let me re-examine. For C2 in arity-5: inner at pos 2-4
              (j=3, s=2), outer params [l1, l2]. Here:
              lam_list[2:4] = [l3,l4], sum = l3+l4.
              But the outer params are [l1, l2], NOT [l1, l2, l3+l4].
              Because the inner output is in the RIGHTMOST slot of the outer m3,
              there's no param after it.

              OK so the pattern is:
              Build outer params from left to right:
              1. All params strictly before the block: lam_list[0:s-1] if s > 0 -> s-1 params
                 Wait no: lam_list[0] through lam_list[s-2] gives s-1 params.
                 But if s=1, that's lam_list[0:0] = [] (0 params before).
                 Hmm, for C1 (s=1, j=3), outer params = [l1, ...].
                 l1 = lam_list[0], which is between pos 0 and pos 1 = start of block.
                 After contraction: between pos 0 and inner output. So it survives.

              Let me think of it as: after contraction, the i outer inputs are at
              positions [0, ..., s-1, OUTPUT, s+j, ..., k-1].
              The outer spectral params are between consecutive outer positions:
              - Between outer pos p and p+1 (for p = 0, ..., i-2):
                If neither p nor p+1 is the OUTPUT position:
                  the param is the original lam_list[original_pos_of(p)].
                If p+1 = OUTPUT position: the param is lam_list[s-1] (the original
                  param between position s-1 and position s).
                If p = OUTPUT position: the param is sum of lam_list[s:s+j-1]
                  (params from start of block up to but not including the param
                  between block end and next position). Wait, that would be
                  sum(lam_list[s:s+j-1]) = sum of j-1 internal params.
                  But for C0 in arity-5: output at slot 0, sum should be...
                  After contraction: [OUTPUT, 3, 4]. Param between OUTPUT and 3
                  is the combined param spanning from block-start(0) to 3.
                  Original params in that span: lam_list[0], lam_list[1], lam_list[2]
                  = l1, l2, l3. Sum = l1+l2+l3. But internal params are only
                  lam_list[0], lam_list[1] (j-1=2 params). The extra lam_list[2]=l3
                  is the param between pos 2 (block end) and pos 3 (first after block).

                  So: param between OUTPUT and first-after-block =
                  sum(internal params) + boundary_after = sum(lam_list[s:s+j-1]).

                  For s=0, j=3: sum(lam_list[0:3]) = l1+l2+l3. Here s+j-1 = 2,
                  so lam_list[0:3] actually is lam_list[0], lam_list[1], lam_list[2].
                  That's 3 params, which is j=3 params... but j-1=2 internal params.
                  Wait: lam_list[s:s+j-1] = lam_list[0:2] = [l1, l2] (2 params).
                  sum = l1+l2 != l1+l2+l3.

                  Hmm, there's an off-by-one. Let me recount.
                  For C0: outer params = [l1+l2+l3, l4] from the existing code.
                  Inner at positions 0,1,2 (j=3 inputs, s=0).
                  Boundary between block-end (pos 2) and first-after (pos 3) = lam_list[2] = l3.
                  Internal params: lam_list[0]=l1, lam_list[1]=l2 (j-1=2 params).
                  So merged param = l1 + l2 + l3 = sum(internal) + boundary = sum(lam_list[s:s+j]).
                  That's lam_list[0:3] = [l1,l2,l3], sum = l1+l2+l3.

                  For C2: inner at positions 2,3,4 (j=3, s=2). Output is rightmost.
                  Internal params: lam_list[2]=l3, lam_list[3]=l4 (j-1=2 params).
                  No boundary after (block ends at last position).
                  So there is no merged param on the right. Outer params = [l1, l2].
                  Before-block params: lam_list[0]=l1, lam_list[1]=l2 (s-1=1 params
                  between positions before block, but wait s=2, so positions before
                  block are 0,1. Params between them: lam_list[0]=l1. Plus param
                  between pos 1 and block start (pos 2): lam_list[1]=l2.
                  Total: [l1, l2]. Correct!

        OK so the general rule:
        Outer params = params_before + [merged_param] + params_after
        where:
          params_before = lam_list[0:s-1] if s >= 1, else []
            Wait: lam_list[0:s-1] gives params from slot 0 to slot s-2.
            But we also need the param between position s-1 and the block start.
            Hmm no: if s=2, params_before should be [l1, l2]:
            l1 between pos 0 and 1, l2 between pos 1 and 2=block start.
            That's lam_list[0:2] = [l1, l2] = lam_list[0:s]. So s params.
            But wait: for s=0 there are no params before.
            lam_list[0:0] = []. Correct.
            For s=1: lam_list[0:1] = [l1]. That's the param between pos 0
            and pos 1 = block start. After contraction: between pos 0 and OUTPUT.
            Correct.

          Actually: params_before = lam_list[0:s] gives s params... but that
          CAN'T be right because if s > 0, lam_list[s-1] is the param between
          position s-1 and position s = block start. After contraction, this
          becomes the param between position s-1 and OUTPUT. It does survive.
          And lam_list[0:s] = [lam_list[0], ..., lam_list[s-1]].
          For C1 (s=1): lam_list[0:1] = [l1]. Outer params = [l1, merged].
          merged = sum(lam_list[1:4]) = l2+l3+l4. Outer = [l1, l2+l3+l4]. Correct!

          For C0 (s=0): lam_list[0:0] = []. Outer params = [merged, ...].
          merged = sum(lam_list[0:3]) = l1+l2+l3. params_after = lam_list[3:4] = [l4].
          Outer = [l1+l2+l3, l4]. Correct!

          For C2 (s=2): lam_list[0:2] = [l1, l2]. No merged param (block is rightmost).
          Outer = [l1, l2]. Correct!

        So: params_before = list(lam_list[0:s])  (s items, can be empty)
            merged_param = sum(lam_list[s:s+j-1])  if s+j-1 < k-1 (not rightmost)
            Wait: s+j-1 is the index of the last inner position.
            The block spans positions s to s+j-1.
            If s+j-1 < k-1 (not the last position), then there's a param after:
            lam_list[s+j-1] between block-end and first-after.
            merged_param = sum(lam_list[s:s+j])... no.

            Let me just directly compute:
            merged_param = sum(lam_list[s:s+j-1])  -- this gives internal params
            PLUS the boundary: if s+j <= k-1, add lam_list[s+j-1].
            Hmm, for C0: s=0, j=3: sum(lam_list[0:2]) = l1+l2, boundary = lam_list[2]=l3.
            Total = l1+l2+l3. But sum(lam_list[0:3]) = l1+l2+l3 directly.

            So merged_param = sum(lam_list[s:s+j-1+1]) = sum(lam_list[s:s+j])
            when the block is not rightmost?
            s+j-1 < k-1 means s+j < k, so s+j <= k-1, and lam_list has k-1 entries.
            lam_list[s:s+j] is valid (j entries).

            For C0: lam_list[0:3] = [l1,l2,l3]. Sum = l1+l2+l3. Correct!
            For C1: lam_list[1:4] = [l2,l3,l4]. Sum = l2+l3+l4. Correct!

            When the block IS rightmost (s+j-1 = k-1, i.e. s+j = k):
            No merged param needed (output is at the rightmost slot of outer).
            Just: params_before = lam_list[0:s].

            params_after = lam_list[s+j:k-1] if s+j < k-1 else []
            (params after the merged gap)

        Let me just implement this directly:
        """
        n = k  # total arity
        outer_params = []

        # Params before block: lam_list[0:s] (between positions 0..s-1 and block start)
        for p in range(s):
            if p < len(lam_list):
                outer_params.append(lam_list[p])

        # Merged param (if block is not at the rightmost position)
        block_end = s + j - 1  # last position in block
        if block_end < n - 1:
            # Block is not rightmost. Merged param = sum of params spanning the gap.
            # The gap covers: from block start to the first position after block.
            # Params: lam_list[s], lam_list[s+1], ..., lam_list[s+j-1]
            # That's j params.
            merged = sum(lam_list[s:s+j])
            outer_params.append(merged)

        # Params after the block+gap: lam_list[s+j:k-2]
        # These are params between positions after the gap.
        for p in range(s + j, len(lam_list)):
            outer_params.append(lam_list[p])

        # The outer should have i-1 = k-j params
        expected = k - j
        if len(outer_params) != expected:
            # Debugging: this shouldn't happen
            raise ValueError(
                f"Outer params mismatch: got {len(outer_params)}, "
                f"expected {expected} (k={k}, j={j}, s={s})"
            )

        return tuple(outer_params)


def extract_depth_spectrum(engine: StasheffEngine, k: int,
                            n_samples: int = 50,
                            seed: int = 42) -> Dict[str, object]:
    """Extract the depth spectrum of m_k by sampling random spectral params.

    Returns:
        dict with:
          'populated_depths': set of depths that are nonzero for some sample
          'field_populated': {deriv_order: True} for populated fields
          'scalar_populated': True/False
          'min_depth': minimum populated depth in T-dependent sector
          'max_depth_T': maximum populated T-dependent depth
          'gap_at_k': whether depth k is absent (structural gap)
          'scalar_depth': depth of scalar = k+1
    """
    rng = random.Random(seed)
    n_lams = k - 1

    # Track which derivative orders are populated
    populated: Dict[int, bool] = {}
    scalar_values = []

    for _ in range(n_samples):
        lams = tuple(rng.uniform(0.1, 5.0) for _ in range(n_lams))
        result = engine.mk(lams)

        for deriv_order, coeff in result.items():
            if abs(coeff) > 1e-10:
                populated[deriv_order] = True

        scalar_values.append(result.get(-1, 0.0))

    # Convert to depth spectrum
    # Weight w = deriv_order, depth d = k-1-w for T-dependent
    # Scalar: depth = k+1
    depths_T = set()
    for deriv_order in populated:
        if deriv_order == -1:
            continue
        d = k - 1 - deriv_order
        depths_T.add(d)

    scalar_present = -1 in populated

    all_depths = depths_T.copy()
    if scalar_present:
        all_depths.add(k + 1)

    return {
        'populated_depths': sorted(all_depths),
        'depths_T': sorted(depths_T),
        'scalar_present': scalar_present,
        'min_depth_T': min(depths_T) if depths_T else None,
        'max_depth_T': max(depths_T) if depths_T else None,
        'gap_at_k': k not in all_depths,
        'scalar_depth': k + 1,
        'k': k,
    }


def compute_L1_norm_T(engine: StasheffEngine, k: int,
                       n_samples: int = 200, seed: int = 123) -> float:
    """Estimate ||m_k|_T||_{L^1} by averaging |T-coefficient| over random λ.

    The T-coefficient of m_k has depth k-1 (the deepest T-dependent depth).
    We average |coeff of T| over random unit-hypercube spectral params.
    """
    rng = random.Random(seed)
    n_lams = k - 1
    total = 0.0

    for _ in range(n_samples):
        lams = tuple(rng.uniform(-1.0, 1.0) for _ in range(n_lams))
        result = engine.mk(lams)
        T_coeff = abs(result.get(0, 0.0))
        total += T_coeff

    return total / n_samples


def compute_shadow_coefficient(engine: StasheffEngine, k: int,
                                n_samples: int = 200, seed: int = 456) -> Tuple[float, float]:
    """Estimate the shadow coefficient S_k from the scalar part of m_k.

    The scalar part is (c/12) * P_k(λ) where P_k is homogeneous of degree k+1.
    S_k is related to the coefficient of the leading shadow tower term.

    Returns (mean_scalar, ratio_to_c12) where ratio = mean_scalar / (c/12).
    """
    c = engine.c
    rng = random.Random(seed)
    n_lams = k - 1
    total = 0.0

    for _ in range(n_samples):
        # Use unit spectral params to get a normalizable quantity
        lams = tuple(1.0 for _ in range(n_lams))
        result = engine.mk(lams)
        scalar = result.get(-1, 0.0)
        total += scalar

    mean_scalar = total / n_samples
    ratio = mean_scalar / (c / 12.0) if abs(c) > 1e-14 else 0.0
    return mean_scalar, ratio


def compute_scalar_at_unit(engine: StasheffEngine, k: int) -> float:
    """Compute scalar part of m_k at λ₁=...=λ_{k-1}=1.

    This gives (c/12) * P_k(1,...,1) which is the "total shadow weight".
    """
    lams = tuple(1.0 for _ in range(k - 1))
    result = engine.mk(lams)
    return result.get(-1, 0.0)


def test_secondary_vanishing(engine: StasheffEngine, k: int,
                              n_samples: int = 100, seed: int = 789) -> Dict[str, object]:
    """Test for secondary vanishing at arity k.

    At k=4, depths 0 and 1 vanish (the n=4 anomaly).
    Check whether similar vanishing occurs at k=8.
    """
    rng = random.Random(seed)
    n_lams = k - 1

    # For each possible derivative order, collect values
    max_deriv = k - 2  # highest possible derivative order
    field_magnitudes: Dict[int, List[float]] = {d: [] for d in range(max_deriv + 1)}
    field_magnitudes[-1] = []  # scalar

    for _ in range(n_samples):
        lams = tuple(rng.uniform(0.5, 3.0) for _ in range(n_lams))
        result = engine.mk(lams)

        for d in range(max_deriv + 1):
            field_magnitudes[d].append(abs(result.get(d, 0.0)))
        field_magnitudes[-1].append(abs(result.get(-1, 0.0)))

    # Determine vanishing: max absolute value < threshold
    threshold = 1e-8
    vanishing = {}
    for d in range(max_deriv + 1):
        max_val = max(field_magnitudes[d]) if field_magnitudes[d] else 0.0
        depth = k - 1 - d
        vanishing[depth] = max_val < threshold

    return {
        'k': k,
        'vanishing_depths': {d: v for d, v in vanishing.items() if v},
        'nonvanishing_depths': {d: v for d, v in vanishing.items() if not v},
        'scalar_max': max(field_magnitudes[-1]) if field_magnitudes[-1] else 0.0,
    }


def test_palindrome_at_symmetric(engine: StasheffEngine, k: int) -> Dict[str, float]:
    """Evaluate m_k|_T at the symmetric point λ₁=...=λ_{k-1}=λ for several λ values.

    Returns the T-coefficient values at different λ, to check factorization.
    """
    results = {}
    for lam in [0.0, 0.5, 1.0, 2.0, -1.0, 0.1, -0.5]:
        lams = tuple(lam for _ in range(k - 1))
        result = engine.mk(lams)
        results[lam] = {
            'T_coeff': result.get(0, 0.0),
            'dT_coeff': result.get(1, 0.0),
            'scalar': result.get(-1, 0.0),
        }
    return results


# ===== Main computation =====

def run_full_analysis():
    """Run the complete analysis for arities 2-10."""

    print("=" * 80)
    print("VIRASORO A∞ DEPTH SPECTRUM: ARITIES 2-10")
    print("Numerical Stasheff recursion engine")
    print("=" * 80)
    print()

    # Test at multiple c values
    c_values = [1.0, 13.0, 26.0, 100.0]

    for c_val in c_values:
        print(f"\n{'=' * 80}")
        print(f"CENTRAL CHARGE c = {c_val}")
        print(f"{'=' * 80}")

        engine = StasheffEngine(c_val)

        # === Verify arities 2-6 ===
        print("\n--- Verification: arities 2-6 ---")
        for k in range(2, 7):
            t0 = time.time()
            spectrum = extract_depth_spectrum(engine, k, n_samples=30)
            elapsed = time.time() - t0
            print(f"  m_{k}: Spec = {spectrum['populated_depths']}, "
                  f"gap at {k}: {spectrum['gap_at_k']}, "
                  f"min_T = {spectrum['min_depth_T']}, "
                  f"time = {elapsed:.2f}s")

        # === Compute arities 7-10 ===
        print("\n--- NEW: arities 7-10 ---")
        for k in range(7, 11):
            print(f"\n  Computing m_{k}...")
            engine._cache.clear()  # Clear cache to avoid memory buildup
            t0 = time.time()

            try:
                spectrum = extract_depth_spectrum(engine, k, n_samples=20, seed=42+k)
                elapsed = time.time() - t0
                print(f"  m_{k}: Spec = {spectrum['populated_depths']}")
                print(f"         Spec|_T = {spectrum['depths_T']}")
                print(f"         gap at d={k}: {spectrum['gap_at_k']}")
                print(f"         min_depth_T = {spectrum['min_depth_T']}")
                print(f"         scalar present: {spectrum['scalar_present']} (depth {spectrum['scalar_depth']})")
                print(f"         time = {elapsed:.2f}s")

                # Secondary vanishing test
                vanish = test_secondary_vanishing(engine, k, n_samples=20, seed=42+k)
                if vanish['vanishing_depths']:
                    print(f"         SECONDARY VANISHING at depths: {list(vanish['vanishing_depths'].keys())}")
                else:
                    print(f"         No secondary vanishing detected")

            except Exception as e:
                elapsed = time.time() - t0
                print(f"  m_{k}: FAILED after {elapsed:.2f}s: {e}")
                import traceback
                traceback.print_exc()
                break

    # === Detailed analysis at c = 1 (fast) ===
    print(f"\n\n{'=' * 80}")
    print("DETAILED ANALYSIS AT c = 1")
    print(f"{'=' * 80}")

    engine = StasheffEngine(1.0)

    # Scalar coefficients at unit spectral params
    print("\n--- Scalar coefficients at λ₁=...=λ_{k-1}=1 ---")
    scalar_data = {}
    for k in range(2, 11):
        try:
            engine._cache.clear()
            sc = compute_scalar_at_unit(engine, k)
            ratio = sc / (1.0 / 12.0) if abs(sc) > 1e-300 else 0.0
            scalar_data[k] = (sc, ratio)
            print(f"  m_{k}: scalar = {sc:.10e}, P_{k}(1,...,1) = {ratio:.6f}")
        except Exception as e:
            print(f"  m_{k}: FAILED: {e}")
            break

    # L^1 norms for Gevrey analysis
    print("\n--- L^1 norms ||m_k|_T|| (averaged over random λ in [-1,1]^{k-1}) ---")
    l1_data = {}
    for k in range(2, 11):
        try:
            engine._cache.clear()
            l1 = compute_L1_norm_T(engine, k, n_samples=100)
            l1_data[k] = l1
            # Gevrey-1: ||m_k|_T|| ~ C * A^k * k!
            # Ratio test: ||m_{k+1}|| / (k * ||m_k||)
            if k > 2 and k - 1 in l1_data and l1_data[k - 1] > 1e-300:
                ratio = l1 / ((k - 1) * l1_data[k - 1])
                print(f"  m_{k}: ||m_k|_T|| = {l1:.6e}, ratio = {ratio:.4f}")
            else:
                print(f"  m_{k}: ||m_k|_T|| = {l1:.6e}")
        except Exception as e:
            print(f"  m_{k}: FAILED: {e}")
            break

    # Palindrome test at c = 1, symmetric point
    print("\n--- Palindrome test at symmetric point λ_i = λ ---")
    for k in [4, 8]:
        try:
            engine._cache.clear()
            pal = test_palindrome_at_symmetric(engine, k)
            print(f"\n  m_{k} at symmetric point:")
            for lam, data in sorted(pal.items()):
                print(f"    λ={lam:5.1f}: T={data['T_coeff']:12.6e}, "
                      f"dT={data['dT_coeff']:12.6e}, "
                      f"sc={data['scalar']:12.6e}")
            # Check if T_coeff vanishes at λ=0
            if abs(pal[0.0]['T_coeff']) < 1e-10:
                print(f"    -> T-coeff vanishes at λ=0 (palindrome factor λ present)")
        except Exception as e:
            print(f"  m_{k}: FAILED: {e}")


    # === Final summary table ===
    print(f"\n\n{'=' * 80}")
    print("SUMMARY: DEPTH SPECTRA m_2 through m_10")
    print(f"{'=' * 80}")
    print()

    engine = StasheffEngine(1.0)
    print(f"{'Arity':>6} {'Spec(m_k|_T)':>30} {'Gap':>12} {'Scalar':>10} {'Min_T':>8} {'n=4 anom':>10}")
    print("-" * 80)

    for k in range(2, 11):
        try:
            engine._cache.clear()
            spec = extract_depth_spectrum(engine, k, n_samples=30, seed=42+k)
            vanish = test_secondary_vanishing(engine, k, n_samples=30, seed=42+k)

            depths_str = str(spec['depths_T'])
            gap_str = f"d={k} {'absent' if spec['gap_at_k'] else 'PRESENT'}"
            sc_str = f"d={spec['scalar_depth']}" if spec['scalar_present'] else "absent"
            min_str = str(spec['min_depth_T'])
            anom_str = ""
            if vanish['vanishing_depths']:
                anom_depths = sorted(vanish['vanishing_depths'].keys())
                anom_str = f"gap@{anom_depths}"

            print(f"{k:>6} {depths_str:>30} {gap_str:>12} {sc_str:>10} {min_str:>8} {anom_str:>10}")
        except Exception as e:
            print(f"{k:>6} FAILED: {e}")
            break


if __name__ == '__main__':
    run_full_analysis()
