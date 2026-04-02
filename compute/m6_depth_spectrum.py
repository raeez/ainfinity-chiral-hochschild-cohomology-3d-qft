r"""Compute m₆(T,T,T,T,T,T; λ₁,...,λ₅) via the Stasheff recursion.

Extends the symbolic_stasheff engine to arity 6, adding d4T field tracking.

The arity-6 Stasheff identity (with m₁ = 0):
  m₆ = -[m₂∘m₅ + m₅∘m₂ + m₃∘m₄ + m₄∘m₃]

where each term is a sum over all insertion positions with alternating signs.

Goal: extract the depth spectrum and verify gap migration at d=6.
"""
from __future__ import annotations

from sympy import Symbol, S, expand, symbols, Poly, collect, Rational

# ===== Extended field system =====
T_sym = Symbol('T')
dT_sym = Symbol('dT')
d2T_sym = Symbol('d2T')
d3T_sym = Symbol('d3T')
d4T_sym = Symbol('d4T')
c_sym = Symbol('c')

FIELD_ORDER = ['d4T', 'd3T', 'd2T', 'dT', 'T', '1']
DERIV_SHIFT = {'T': 'dT', 'dT': 'd2T', 'd2T': 'd3T', 'd3T': 'd4T'}
# Derivative order of each field symbol
DERIV_ORDER = {'T': 0, 'dT': 1, 'd2T': 2, 'd3T': 3, 'd4T': 4, '1': -1}


def _add(*dicts, signs=None):
    """Add field-coeff dicts with optional signs."""
    if signs is None:
        signs = [1] * len(dicts)
    result = {}
    for d, s in zip(dicts, signs):
        for f, c in d.items():
            val = result.get(f, S.Zero)
            result[f] = expand(val + s * c)
    return {k: v for k, v in result.items() if v != 0}


def _scale(d, factor):
    """Scale all coefficients in a field dict."""
    return {f: expand(factor * v) for f, v in d.items() if v != 0}


def _apply_partial(d):
    """Apply ∂ to a field-coeff dict."""
    result = {}
    for f, c in d.items():
        if f == '1':
            continue
        new_f = DERIV_SHIFT.get(f)
        if new_f:
            val = result.get(new_f, S.Zero)
            result[new_f] = expand(val + c)
    return {k: v for k, v in result.items() if v != 0}


def _apply_shift_partial(d, shift):
    """Apply (shift + ∂) to a field-coeff dict."""
    shifted = _scale(d, shift)
    partial = _apply_partial(d)
    return _add(shifted, partial)


def _apply_shift_partial_n(d, shift, n):
    """Apply (shift + ∂)^n to a field-coeff dict via recursion."""
    if n == 0:
        return dict(d)
    result = d
    for _ in range(n):
        result = _apply_shift_partial(result, shift)
    return result


# ===== Core operations =====

def m2_at(lam, c=None):
    """m₂(T,T; λ) = dT + 2Tλ + (c/12)λ³"""
    cc = c if c is not None else c_sym
    return {'dT': S.One, 'T': 2 * lam, '1': cc * lam**3 / 12}


def m3_at(l1, l2, c=None):
    """m₃(T,T,T; l₁, l₂)"""
    cc = c if c is not None else c_sym
    return {
        'd2T': S.One,
        'dT': expand(2*l1 + 3*l2),
        'T': expand(4*l1*l2 + 2*l2**2),
        '1': expand(cc * l2**3 * (2*l1 + l2) / 12),
    }


def m4_at(l1, l2, l3, c=None):
    """m₄(T,T,T,T; l₁,l₂,l₃) — from symbolic_stasheff engine."""
    from compute.lib.symbolic_stasheff import m4_virasoro_symbolic
    cc = c if c is not None else c_sym
    res = m4_virasoro_symbolic(l1, l2, l3, cc)
    return res['m4']


def m5_at(l1, l2, l3, l4, c=None):
    """m₅(T,T,T,T,T; l₁,l₂,l₃,l₄) — from symbolic_stasheff engine."""
    from compute.lib.symbolic_stasheff import m5_virasoro_symbolic
    cc = c if c is not None else c_sym
    res = m5_virasoro_symbolic(l1, l2, l3, l4, cc)
    return res['m5']


# ===== Sesquilinearity composition =====

def compose_into_m2_left(inner, lam_outer, c=None):
    """m₂(inner, T; λ_outer).

    Left sesquilinearity: m₂(d^n T · f, T; λ) = f · (-λ)^n · m₂(T,T;λ).
    Scalars map to 0.
    """
    cc = c if c is not None else c_sym
    lo = lam_outer
    base = m2_at(lo, cc)
    result = {}
    for field, coeff in inner.items():
        if coeff == 0 or field == '1':
            continue
        n = DERIV_ORDER.get(field, -1)
        if n < 0:
            continue
        factor = (-lo)**n
        for bf, bc in base.items():
            val = result.get(bf, S.Zero)
            result[bf] = expand(val + coeff * factor * bc)
    return {k: v for k, v in result.items() if v != 0}


def compose_into_m2_right(inner, lam_outer, c=None):
    """m₂(T, inner; λ_outer).

    Right sesquilinearity: m₂(T, d^n T · f; λ) = f · (λ+∂)^n · m₂(T,T;λ).
    Scalars map to 0.
    """
    cc = c if c is not None else c_sym
    lo = lam_outer
    base = m2_at(lo, cc)
    result = {}
    for field, coeff in inner.items():
        if coeff == 0 or field == '1':
            continue
        n = DERIV_ORDER.get(field, -1)
        if n < 0:
            continue
        # Apply (lo + ∂)^n to base
        shifted_base = _apply_shift_partial_n(base, lo, n)
        for bf, bc in shifted_base.items():
            val = result.get(bf, S.Zero)
            result[bf] = expand(val + coeff * bc)
    return {k: v for k, v in result.items() if v != 0}


def compose_into_mk_slot(mk_base_func, inner, slot, outer_lams, n_slots, c=None):
    """Compose inner output into slot `slot` of m_k evaluated at outer_lams.

    mk_base_func: function(outer_lams..., c) -> field dict for m_k(T,...,T; outer_lams)
    inner: field dict (output of inner operation)
    slot: which slot to insert into (0-indexed)
    outer_lams: list of k-1 spectral parameters for the outer m_k
    n_slots: number of input slots of outer m_k (= k)

    Sesquilinearity rules:
      slot 0 (leftmost): d^n T -> (-outer_lams[0])^n * m_k(T,...,T)
      slot i (0 < i < n_slots-1): d^n T -> (-outer_lams[i])^n * m_k(T,...,T)
      slot n_slots-1 (rightmost): d^n T -> (sum(outer_lams) + ∂)^n * m_k(T,...,T)
    """
    cc = c if c is not None else c_sym
    base = mk_base_func(*outer_lams, cc)

    result = {}
    for field, coeff in inner.items():
        if coeff == 0 or field == '1':
            continue
        n = DERIV_ORDER.get(field, -1)
        if n < 0:
            continue

        if slot < n_slots - 1:
            # Non-rightmost: factor = (-outer_lams[slot])^n
            slot_lam = outer_lams[slot]
            factor = (-slot_lam)**n
            for bf, bc in base.items():
                val = result.get(bf, S.Zero)
                result[bf] = expand(val + coeff * factor * bc)
        else:
            # Rightmost: (sum(outer_lams) + ∂)^n applied to base
            shift = sum(outer_lams)
            shifted_base = _apply_shift_partial_n(base, shift, n)
            for bf, bc in shifted_base.items():
                val = result.get(bf, S.Zero)
                result[bf] = expand(val + coeff * bc)

    return {k: v for k, v in result.items() if v != 0}


def stasheff_rhs_arity6(l1, l2, l3, l4, l5, c=None):
    r"""Compute the RHS of the arity-6 Stasheff identity.

    m₆ = -RHS where RHS = sum of all compositions m_i ∘ m_j with i+j=7, i,j≥2.

    Compositions at arity 6 (6 inputs, partition into inner j and outer i):

    (i=2, j=5): 2 terms
      s=0: m₂(m₅(12345), 6; Σ₁₅+l₅)  -- wait, let me re-index.

    With 6 inputs and 5 spectral params l₁,...,l₅:

    (i=2, j=5): inner takes 5 consecutive inputs
      s=0: +m₂(m₅(T^5; l₁,l₂,l₃,l₄), T; l₁+l₂+l₃+l₄+l₅)
      s=1: -m₂(T, m₅(T^5; l₂,l₃,l₄,l₅); l₁)

    (i=5, j=2): inner takes 2 consecutive inputs
      s=0: +m₅(m₂(T,T;l₁), T, T, T, T; l₁+l₂, l₃, l₄, l₅)
      s=1: -m₅(T, m₂(T,T;l₂), T, T, T; l₁, l₂+l₃, l₄, l₅)
      s=2: +m₅(T, T, m₂(T,T;l₃), T, T; l₁, l₂, l₃+l₄, l₅)
      s=3: -m₅(T, T, T, m₂(T,T;l₄), T; l₁, l₂, l₃, l₄+l₅)
      s=4: +m₅(T, T, T, T, m₂(T,T;l₅); l₁, l₂, l₃, l₄)

    (i=3, j=4): inner takes 4 consecutive inputs
      s=0: +m₃(m₄(T^4; l₁,l₂,l₃), T, T; l₁+l₂+l₃+l₄, l₅)
      s=1: -m₃(T, m₄(T^4; l₂,l₃,l₄), T; l₁, l₂+l₃+l₄+l₅)
      s=2: +m₃(T, T, m₄(T^4; l₃,l₄,l₅); l₁, l₂)

    (i=4, j=3): inner takes 3 consecutive inputs
      s=0: +m₄(m₃(T^3; l₁,l₂), T, T, T; l₁+l₂+l₃, l₄, l₅)  -- wait, spectral params?

    Actually let me be more careful about spectral parameter merging.

    For m_i(... m_j(inputs[s:s+j]; inner_lams) ...; outer_lams):
    - The inner m_j takes inputs at positions s, s+1, ..., s+j-1
    - The inner spectral params are l_{s+1}, ..., l_{s+j-1} (the j-1 params between consecutive inputs)
    - After composition: i inputs remain, with i-1 spectral params.
    - The outer spectral params: take the original l₁,...,l₅, remove l_{s+1},...,l_{s+j-1},
      and merge l_s and l_{s+j} into a single parameter l_s + l_{s+j}...

    Actually no. Let me think about this differently.

    With 6 inputs labeled a₁,...,a₆ and 5 spectral params l₁,...,l₅ (where l_i sits
    between a_i and a_{i+1}):

    For inner m_j at position s (starting at input a_{s+1}):
    - Inner takes inputs a_{s+1}, ..., a_{s+j}
    - Inner spectral params: l_{s+1}, ..., l_{s+j-1}
    - After composition: remaining inputs are a₁,...,a_s, [inner output], a_{s+j+1},...,a₆
    - Remaining spectral params: l₁,...,l_s, l_{s+j},...,l₅
      BUT the spectral param between a_s and [inner output] is l_s → becomes l_{combined}
      Wait, actually the outer spectral params are:
      l₁,...,l_{s-1}, [l_s if s>0 else nothing], [l_{s+j} if s+j<6 else nothing], l_{s+j+1},...,l₅

    Let me use a concrete indexing. Inputs: positions 0,1,2,3,4,5.
    Spectral params: l[0]=l₁ between pos 0-1, l[1]=l₂ between 1-2, ..., l[4]=l₅ between 4-5.

    For inner m_j taking positions [s, s+1, ..., s+j-1]:
    - Inner spectral params: l[s], l[s+1], ..., l[s+j-2]  (j-1 params)
    - After merging, the outer sees i inputs at positions:
      0,...,s-1, [output], s+j,...,5
    - Outer spectral params (i-1 params):
      For position p in the merged sequence:
        l'[0] = l[0], ..., l'[s-2] = l[s-2] (params before the inner block)
        l'[s-1] = l[s-1] IF s > 0, merged with nothing special
        ...

    This indexing is getting complicated. Let me just hardcode each composition carefully.
    """
    cc = c if c is not None else c_sym
    lams = [l1, l2, l3, l4, l5]

    terms = []  # list of (field_dict, sign)

    # ========================================
    # (i=2, j=5): m₂ ∘ m₅
    # ========================================

    # s=0: m₂(m₅(T^5; l₁,l₂,l₃,l₄), T; l₁+l₂+l₃+l₄+l₅)
    # Inner: positions 0-4, params l₁,l₂,l₃,l₄
    # Outer: 2 inputs, param = sum of all
    inner_05 = m5_at(l1, l2, l3, l4, cc)
    comp_25_s0 = compose_into_m2_left(inner_05, l1+l2+l3+l4+l5, cc)
    terms.append((comp_25_s0, +1))
    print("  Done: m₂∘m₅ s=0")

    # s=1: m₂(T, m₅(T^5; l₂,l₃,l₄,l₅); l₁)
    # Inner: positions 1-5, params l₂,l₃,l₄,l₅
    inner_15 = m5_at(l2, l3, l4, l5, cc)
    comp_25_s1 = compose_into_m2_right(inner_15, l1, cc)
    terms.append((comp_25_s1, -1))
    print("  Done: m₂∘m₅ s=1")

    # ========================================
    # (i=5, j=2): m₅ ∘ m₂
    # ========================================

    for s in range(5):
        inner_lam = lams[s]
        m2_out = m2_at(inner_lam, cc)

        # Merge spectral params: remove l[s], replace by merged params
        # Outer m₅ has 5 inputs, 4 spectral params
        # If inner is at position s (between inputs s and s+1):
        #   Merged outer params: l₁,...,l_{s-1}, l_s + l_{s+1}, l_{s+2},...,l₅
        #   Wait no. Inner m₂ takes inputs at positions s and s+1.
        #   After merging: positions 0,...,s-1,[output],s+2,...,5 → 5 positions
        #   Outer params: l₁,...,l_{s-1}, (if s>0: l_s is between pos s-1 and output),
        #                 (l_{s+1} is between output and pos s+2), l_{s+2},...,l₅
        #   But inner m₂ consumed l_{s+1} (the param between its two inputs).
        #   Wait, inner m₂ between positions s and s+1 uses the spectral param l_{s+1}?
        #   No: l_i is between input i and input i+1. So l_s is between input s and s+1.
        #   Inner m₂ at positions (s, s+1) uses spectral param l_s (= lams[s]).
        #
        #   After merging: 5 inputs. Outer spectral params (4 params):
        #   If s=0: merged inputs are [m₂ result], 2, 3, 4, 5
        #     Outer params: l₁+l₂, l₃, l₄, l₅  wait...
        #     No: the original l₁ is consumed by inner m₂.
        #     Between [output] and input 2: this is l₂ (originally between 1 and 2).
        #     So outer params: l₂, l₃, l₄, l₅  with output in slot 0.
        #     But wait: the output has spectral weight l₁ embedded. So the effective
        #     outer param between merged-input-0 and merged-input-1 should account for this.
        #
        #   Actually in the Stasheff recursion formalism:
        #   m₅(T,...,m₂(T,T;l_{s+1}),...,T) where m₂ replaces inputs at positions s, s+1.
        #   The outer m₅ sees 5 inputs with 4 spectral params.
        #   The spectral params for the outer m₅ are obtained by:
        #     - Taking the original 5 params l₁,...,l₅
        #     - Removing l_{s+1} (consumed by inner)
        #     - Merging l_s with l_{s+1}: the param l_s (before inner block) remains,
        #       and l_{s+2} (after inner block) remains. The consumed param is l_{s+1}.
        #
        # Let me use the standard convention from the existing code (stasheff_rhs_arity5):
        # Inner m₂ at position s takes the spectral param lams[s].
        # After composition, the outer m₅ has spectral params obtained by:
        #   removing lams[s] and merging: lams[:s] + [lams[s] + lams[s+1]] + lams[s+2:]
        #   if s < 4 (not the last position)
        #   For s=4 (last): lams[:4] (just drop the last one)
        #
        # Wait, that's what the existing arity-5 code does. Let me recheck.

        if s < 4:
            outer_lams = list(lams[:s]) + [lams[s] + lams[s+1]] + list(lams[s+2:])
        else:
            outer_lams = list(lams[:s]) + list(lams[s+1:])

        # outer_lams should have 4 elements for m₅
        assert len(outer_lams) == 4, f"s={s}: outer_lams has {len(outer_lams)} elements"

        # Now compose m₂ output into slot s of m₅
        comp = compose_into_mk_slot(m5_at, m2_out, s, outer_lams, 5, cc)
        sign = (-1)**s
        terms.append((comp, sign))
        print(f"  Done: m₅∘m₂ s={s}")

    # ========================================
    # (i=3, j=4): m₃ ∘ m₄
    # ========================================

    for s in range(3):
        # Inner m₄ at positions [s, s+1, s+2, s+3]
        # Inner spectral params: lams[s], lams[s+1], lams[s+2]
        inner_m4 = m4_at(lams[s], lams[s+1], lams[s+2], cc)

        # Outer m₃ has 3 inputs, 2 spectral params
        # Original params: l₁,...,l₅. Inner consumes lams[s], lams[s+1], lams[s+2].
        # Remaining params for outer:
        # Before inner: lams[:s-1] (params between inputs before the block)
        #   If s=0: no params before
        #   If s=1: lams[0] (between input 0 and [output])
        #   If s=2: lams[0], lams[1] (between inputs 0,1 and [output])
        # After inner: lams[s+3:] (params between [output] and inputs after)
        #   If s=0: lams[3], lams[4]
        #   If s=1: lams[4]
        #   If s=2: (none)

        # For m₃ with 3 inputs, outer needs 2 spectral params.
        # Input arrangement after merge:
        #   Inputs before inner: positions 0,...,s-1
        #   [inner output]
        #   Inputs after inner: positions s+4,...,5
        # Total: s + 1 + (6-s-4) = 3 inputs. Good.

        # Outer spectral params:
        # Between input(s-1) and [output]: this is the param that was originally
        # between input(s-1) and input(s) = lams[s-1] if s > 0.
        # But wait, inner consumed positions s through s+3. The param lams[s-1]
        # is between position s-1 and position s. Position s is the first inner input.
        # After merging, lams[s-1] becomes the param between last-before-input and output.

        # Between [output] and input(s+4): this is lams[s+3] (between position s+3 and s+4).

        # So outer params:
        #   slot_before = lams[:s] (s params before the inner block)
        #   But we only have s inputs before the block, which means s-1 params between them
        #   plus 1 param between last-before and output.
        #   Wait, for 3 outer inputs we need exactly 2 outer params.

        # Let me think slot by slot:
        #   s=0: [m₄(0123)], 4, 5. Inner consumed lams[0,1,2].
        #     Outer param 1: between [output] and input 4 = lams[3]
        #     Outer param 2: between input 4 and input 5 = lams[4]
        #     Output in slot 0. Outer params: [lams[3], lams[4]]
        #     But outer sum = l₁+l₂+l₃ (from inner) → first outer param should be l₁+l₂+l₃+l₄? No.
        #     Actually: the outer m₃ just sees 3 inputs with params [l₄, l₅].
        #     The inner output carries spectral content, but the outer params are just [l₄, l₅].
        #     Wait, but what about the sum? The param between inputs in m₃ is just the
        #     spectral variable of that slot. Let me re-examine.

        #   Actually, from the existing arity-5 code:
        #   C0 = _compose_m3_slot(m3_012_full, 0, l1 + l2 + l3, l4, cc)
        #   So when inner m₃(l₁,l₂) output goes into slot 0 of outer m₃,
        #   the outer params are (l₁+l₂+l₃, l₄) — the first outer param is the SUM
        #   of the inner params plus the next original param.

        #   General rule: for inner at position s taking k consecutive inputs with
        #   params lams[s]...lams[s+j-2], the outer params merge as:
        #   ...lams[s-1], SUM(lams[s]...lams[s+j-1]), lams[s+j]...
        #   Wait no. Let me look at this from the composition formula perspective.

        #   In the A∞ composition: m_i(..., m_j(a_{s+1},...,a_{s+j}; λ_inner), ...; λ_outer)
        #   where the outer m_i has i=n-j+1 inputs and i-1=n-j spectral params.
        #   The outer spectral params are: all original params EXCEPT those consumed by inner,
        #   with the params at the boundary of the inner block merged.

        #   For arity 6, inner m₄ at position s (inputs s,s+1,s+2,s+3):
        #   Inner params: lams[s], lams[s+1], lams[s+2]
        #   Outer m₃ params (2 params):
        #     If s=0: inner takes inputs 0-3, outer sees [output],4,5
        #       Outer params: [lams[3], lams[4]]
        #       But with spectral sum from inner: the param between [output] and input 4
        #       carries the total spectral content. In the lambda-bracket formalism,
        #       the outer param is lams[0]+lams[1]+lams[2]+lams[3] for the first slot?
        #       NO. Looking at the arity-5 code more carefully:

        #   In stasheff_rhs_arity5:
        #   C0: m₃(m₃(T^3;l₁,l₂), T, T; l₁+l₂+l₃, l₄)
        #   So outer params are [l₁+l₂+l₃, l₄]. The first outer param is the SUM of the
        #   inner params (l₁,l₂) plus the NEXT original param (l₃).

        #   C1: m₃(T, m₃(T^3;l₂,l₃), T; l₁, l₂+l₃+l₄)
        #   Outer params: [l₁, l₂+l₃+l₄]. Second outer param = sum of inner params + next.

        #   C2: m₃(T, T, m₃(T^3;l₃,l₄); l₁, l₂)
        #   Outer params: [l₁, l₂]. Just the remaining original params.

        #   Pattern: the merged param at the inner's position is the sum of all
        #   original params from the start of the inner block to the end of the
        #   gap it creates. Specifically:
        #   If inner takes positions [s...s+j-1] with params lams[s]...lams[s+j-2],
        #   the merged param at that position is lams[s-1] + ... + lams[s+j-2] + ???

        #   Actually, simpler: the outer params are the original params with the inner
        #   block contracted. The inner block of j inputs has j-1 internal params that
        #   are consumed. The boundary params survive but get merged:
        #   - The param before the block (lams[s-1]) and the param after (lams[s+j-1])
        #     are NOT consumed. The j-1 internal params (lams[s]...lams[s+j-2]) ARE consumed.
        #   - But the slot in the outer operation where the inner output sits needs a
        #     spectral param on each side:
        #     Left side: original lams[s-1] (if s > 0)
        #     Right side: original lams[s+j-1] (if s+j-1 < 5)
        #   - These ARE the outer params. But wait, that doesn't match the arity-5 code.

        #   In the arity-5 code for C0: m₃(m₃(T,T,T;l₁,l₂), T, T; L, M)
        #   Inner at positions 0,1,2 with params l₁,l₂. The inner output is at slot 0 of outer.
        #   Outer m₃ has params L = l₁+l₂+l₃, M = l₃... no, L = l₁+l₂+l₃ and M = l₄.
        #   Hmm. So L = sum of lams[0:3] = l₁+l₂+l₃.

        #   This suggests: in the lambda-bracket formalism, when inner m_j output
        #   is placed in slot p of outer m_i, the spectral parameter at slot p
        #   of the outer is the SUM of all original params that correspond to
        #   that gap. The gap between [output] and the next input spans from
        #   input s to input s+j, so the params are lams[s-1], lams[s], ..., lams[s+j-1]?
        #
        #   For C0: outer params are [l₁+l₂+l₃, l₄].
        #   Inner at pos 0-2. The gap between [output] and input 3 spans
        #   params l₁, l₂, l₃ (lams[0], lams[1], lams[2]). Sum = l₁+l₂+l₃. Correct!
        #   The gap between input 3 and input 4 is just l₄ (lams[3]). Correct!

        #   For C1: outer params are [l₁, l₂+l₃+l₄].
        #   Inner at pos 1-3. The gap between input 0 and [output] is l₁ (lams[0]).
        #   The gap between [output] and input 4 spans lams[1]+lams[2]+lams[3] = l₂+l₃+l₄. Correct!

        #   For C2: outer params are [l₁, l₂].
        #   Inner at pos 2-4. The gap between input 0 and input 1 is l₁.
        #   The gap between input 1 and [output] is l₂ (lams[1]).
        #   Inner output is rightmost, no gap after. Correct!

        #   General rule: outer params are formed by taking the original 5 params
        #   and replacing the inner block with a single gap whose param is the
        #   sum of all original params that fall within that block.
        #   Block at positions [s, s+j-1]: the internal params are lams[s]...lams[s+j-2].
        #   PLUS the boundary param after the block: lams[s+j-1] IF s+j-1 < 5.
        #   Wait no, that gives l₁+l₂+l₃ for C0 which includes lams[2]=l₃ but
        #   lams[2] is between input 2 and input 3, and input 3 is the NEXT input after the block.
        #   So the "boundary after" param IS included.

        #   Let me just use: the merged param at slot p of outer = sum of all original
        #   params in the range that collapses.

        # OK let me just build the outer params directly for each case.

        if s == 0:
            # Inner m₄ at pos 0-3, params l₁,l₂,l₃. Outer m₃: [output],4,5
            # Outer param 1 (between [output] and 4): l₁+l₂+l₃+l₄
            # Wait, from the pattern above: the gap between [output] and input 4
            # spans from input 0 to input 4, so params l₁,l₂,l₃,l₄... no.
            # From the arity-5 pattern: inner at pos 0-2 gives outer param l₁+l₂+l₃.
            # That's 3 params for a block of 3 inputs (j=3 means j-1=2 internal params
            # plus 1 boundary param = 3 total summed).
            # For j=4: j-1=3 internal params (l₁,l₂,l₃) plus 1 boundary param (l₄) = 4?
            # But wait: internal params of block at positions 0-3 are lams[0],lams[1],lams[2] (between 0-1, 1-2, 2-3).
            # Boundary param after: lams[3] (between 3 and 4).
            # Sum = l₁+l₂+l₃+l₄.
            outer_p = [l1+l2+l3+l4, l5]
            inner_slot = 0
        elif s == 1:
            # Inner m₄ at pos 1-4, params l₂,l₃,l₄. Outer m₃: 0,[output],5
            # Gap between 0 and [output]: just l₁
            # Gap between [output] and 5: internal params l₂,l₃,l₄ + boundary l₅ = l₂+l₃+l₄+l₅
            outer_p = [l1, l2+l3+l4+l5]
            inner_slot = 1
        elif s == 2:
            # Inner m₄ at pos 2-5, params l₃,l₄,l₅. Outer m₃: 0,1,[output]
            # Gap between 0 and 1: l₁
            # Gap between 1 and [output]: l₂
            # [output] is rightmost, no gap after.
            outer_p = [l1, l2]
            inner_slot = 2

        comp = compose_into_mk_slot(m3_at, inner_m4, inner_slot, outer_p, 3, cc)
        sign = (-1)**s
        terms.append((comp, sign))
        print(f"  Done: m₃∘m₄ s={s}")

    # ========================================
    # (i=4, j=3): m₄ ∘ m₃
    # ========================================

    for s in range(4):
        # Inner m₃ at positions [s, s+1, s+2]
        # Inner params: lams[s], lams[s+1]
        inner_m3 = m3_at(lams[s], lams[s+1], cc)

        # Outer m₄ has 4 inputs, 3 spectral params.
        # Block at positions [s, s+2]:
        #   Internal params: lams[s], lams[s+1] (between s-s+1 and s+1-s+2)
        #   Boundary after: lams[s+2] (between s+2 and s+3) if s+2 < 5

        if s == 0:
            # Inner at pos 0-2, params l₁,l₂. Outer m₄: [output],3,4,5
            # Gap: [output]-3: l₁+l₂+l₃. 3-4: l₄. 4-5: l₅.
            outer_p = [l1+l2+l3, l4, l5]
            inner_slot = 0
        elif s == 1:
            # Inner at pos 1-3, params l₂,l₃. Outer m₄: 0,[output],4,5
            # Gap: 0-[output]: l₁. [output]-4: l₂+l₃+l₄. 4-5: l₅.
            outer_p = [l1, l2+l3+l4, l5]
            inner_slot = 1
        elif s == 2:
            # Inner at pos 2-4, params l₃,l₄. Outer m₄: 0,1,[output],5
            # Gap: 0-1: l₁. 1-[output]: l₂. [output]-5: l₃+l₄+l₅.
            outer_p = [l1, l2, l3+l4+l5]
            inner_slot = 2
        elif s == 3:
            # Inner at pos 3-5, params l₄,l₅. Outer m₄: 0,1,2,[output]
            # Gap: 0-1: l₁. 1-2: l₂. 2-[output]: l₃.
            outer_p = [l1, l2, l3]
            inner_slot = 3

        comp = compose_into_mk_slot(m4_at, inner_m3, inner_slot, outer_p, 4, cc)
        sign = (-1)**s
        terms.append((comp, sign))
        print(f"  Done: m₄∘m₃ s={s}")

    # ========================================
    # Sum all terms
    # ========================================
    total = {}
    for term_dict, sign in terms:
        for f, v in term_dict.items():
            val = total.get(f, S.Zero)
            total[f] = expand(val + sign * v)

    return {k: v for k, v in total.items() if v != 0}


def m6_virasoro_symbolic(l1, l2, l3, l4, l5, c=None):
    """Compute m₆(T^6; l₁,...,l₅) = -RHS."""
    cc = c if c is not None else c_sym
    rhs = stasheff_rhs_arity6(l1, l2, l3, l4, l5, cc)
    m6 = {f: expand(-v) for f, v in rhs.items()}
    return {'m6': m6, 'rhs': rhs}


def depth_spectrum(field_dict, lam_symbols):
    """Extract the depth spectrum from a field-coeff dict.

    Returns dict {field: set of total degrees in lambda symbols}.
    """
    result = {}
    all_degs = set()
    for f in FIELD_ORDER:
        coeff = field_dict.get(f, S.Zero)
        if coeff == 0:
            continue
        coeff = expand(coeff)
        p = Poly(coeff, *lam_symbols)
        degs = set()
        for monom, _ in p.as_dict().items():
            degs.add(sum(monom))
        result[f] = sorted(degs)
        all_degs |= degs
    return result, sorted(all_degs)


if __name__ == '__main__':
    import sys

    l1, l2, l3, l4, l5 = symbols('l1 l2 l3 l4 l5')
    c = c_sym

    print("=" * 70)
    print("COMPUTING m₆(T,T,T,T,T,T; l₁,...,l₅) via Stasheff recursion")
    print("=" * 70)
    print()

    # First verify arities 2-5
    print("--- Verification: arities 2-5 ---")

    m2 = m2_at(l1, c)
    ds2, spec2 = depth_spectrum(m2, [l1])
    print(f"m₂: fields={ds2}, Spec={spec2}")

    m3 = m3_at(l1, l2, c)
    ds3, spec3 = depth_spectrum(m3, [l1, l2])
    print(f"m₃: fields={ds3}, Spec={spec3}")

    m4 = m4_at(l1, l2, l3, c)
    ds4, spec4 = depth_spectrum(m4, [l1, l2, l3])
    print(f"m₄: fields={ds4}, Spec={spec4}")

    print()
    print("Computing m₅...")
    m5 = m5_at(l1, l2, l3, l4, c)
    ds5, spec5 = depth_spectrum(m5, [l1, l2, l3, l4])
    print(f"m₅: fields={ds5}, Spec={spec5}")

    print()
    print("=" * 70)
    print("Computing m₆ (this may take a few minutes)...")
    print("=" * 70)
    sys.stdout.flush()

    result = m6_virasoro_symbolic(l1, l2, l3, l4, l5, c)
    m6 = result['m6']

    print()
    print("--- m₆ field content ---")
    for f in FIELD_ORDER:
        v = m6.get(f, S.Zero)
        print(f"  {f}: {'NONZERO' if v != 0 else 'ZERO'}")

    ds6, spec6 = depth_spectrum(m6, [l1, l2, l3, l4, l5])
    print()
    print(f"m₆ depth spectrum: {ds6}")
    print(f"m₆ Spec = {spec6}")
    print()

    # Gap migration check
    print("=" * 70)
    print("GAP MIGRATION THEOREM VERIFICATION")
    print("=" * 70)
    gap_at_6 = 6 not in spec6
    scalar_at_7 = 7 in spec6 and '1' in ds6 and 7 in ds6['1']
    max_field_depth = max(d for f, degs in ds6.items() if f != '1' for d in degs) if any(f != '1' for f in ds6) else -1

    print(f"  Gap at d=6: {'YES ✓' if gap_at_6 else 'NO ✗'}")
    print(f"  Scalar at d=7: {'YES ✓' if scalar_at_7 else 'NO ✗'}")
    print(f"  Max field depth: {max_field_depth} (should be ≤ 5)")

    min_depth = min(spec6) if spec6 else None
    print(f"  Min depth: {min_depth}")
    print(f"  n=4 anomaly pattern (min depth > 0): {'YES' if min_depth and min_depth > 0 else 'NO'}")

    print()
    print("=" * 70)
    print("COMPARISON TABLE")
    print("=" * 70)
    print(f"  n=2: Spec={spec2}, gap at 2: {2 not in spec2}")
    print(f"  n=3: Spec={spec3}, gap at 3: {3 not in spec3}")
    print(f"  n=4: Spec={spec4}, gap at 4: {4 not in spec4}")
    print(f"  n=5: Spec={spec5}, gap at 5: {5 not in spec5}")
    print(f"  n=6: Spec={spec6}, gap at 6: {6 not in spec6}")

    # Extract leading scalar coefficient at maximum depth
    print()
    print("=" * 70)
    print("SCALAR SECTOR ANALYSIS")
    print("=" * 70)
    scalar_coeff = expand(m6.get('1', S.Zero))
    if scalar_coeff != 0:
        p = Poly(scalar_coeff, l1, l2, l3, l4, l5)
        max_deg = max(sum(m) for m in p.as_dict().keys())
        print(f"  Scalar max degree: {max_deg}")
        print(f"  Scalar c-dependence: {'c' in str(scalar_coeff)}")
        # Extract just the max-degree terms
        max_terms = {m: v for m, v in p.as_dict().items() if sum(m) == max_deg}
        print(f"  Number of max-degree monomials: {len(max_terms)}")
        # Check if all max-degree coefficients are proportional to c
        all_prop_c = all('c' in str(v) for v in max_terms.values())
        print(f"  All max-degree terms proportional to c: {all_prop_c}")
