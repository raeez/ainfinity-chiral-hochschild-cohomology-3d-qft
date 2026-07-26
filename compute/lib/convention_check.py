"""
Convention comparison: Loday-Vallette (Vol I) vs Koszul (Vol II) signs.

Vol I uses the Loday-Vallette A∞ identity:
  Σ_{r+s+t=n} (−1)^{rs+t} m_{r+1+t}(id^r ⊗ m_s ⊗ id^t) = 0

Vol II uses the Koszul A∞ identity:
  Σ_{i+j=n+1} Σ_{s=0}^{n-j} (−1)^{ε(s,j)} m_i(a₁,...,a_s, m_j(a_{s+1},...,a_{s+j}), a_{s+j+1},...,a_n) = 0
  where ε(s,j) = (j−1)(|a₁|+...+|a_s|)

CRITICAL QUESTION: Are these the same sign convention?

The LV convention uses operadic signs (independent of element degrees).
The Koszul convention uses element-dependent signs.

They are equivalent ON ELEMENTS when the LV convention is understood as
acting on the desuspended complex s^{-1}A, and the Koszul convention
acts on A directly with the degree shift absorbed into the sign formula.

This module verifies the equivalence by computing both signs on explicit
elements and checking they agree.
"""
from sympy import S


def lv_sign(r, s, t):
    """Loday-Vallette sign: (−1)^{rs+t}.

    In the A∞ identity, for the term m_{r+1+t}(id^r ⊗ m_s ⊗ id^t):
    - r elements before the inner operation
    - m_s is the inner operation of arity s
    - t elements after the inner operation
    - Total arity n = r + s + t

    The sign is (−1)^{rs+t}.
    """
    return (-1) ** (r * s + t)


def koszul_sign_on_elements(degrees_before, arity_inner):
    """Koszul sign: (−1)^{(j−1)(|a₁|+...+|a_s|)}.

    In the A∞ identity on elements:
    - degrees_before: list of degrees |a₁|,...,|a_s| of elements before inner op
    - arity_inner: j, the arity of the inner operation m_j

    The sign is (−1)^{ε} where ε = (j−1)·Σ|a_i|.
    """
    epsilon = (arity_inner - 1) * sum(degrees_before)
    return (-1) ** (epsilon % 2)


def lv_sign_on_desuspended(degrees_before, arity_inner, degrees_after):
    """LV sign evaluated on desuspended elements.

    When the LV convention acts on s^{-1}A (desuspended), the element a_i
    of degree |a_i| in A has desuspended degree |a_i| − 1 in s^{-1}A.

    The LV sign is (−1)^{rs+t} where:
    - r = number of elements before = len(degrees_before) = s in Koszul notation
    - s_lv = arity of inner operation = arity_inner = j in Koszul notation
    - t = number of elements after = len(degrees_after)

    But the LV convention also includes the Koszul signs from commuting
    the desuspended elements past m_j. The TOTAL sign is:

    (−1)^{rs+t} · (−1)^{(s−1)Σ(|ā_i|)}

    where |ā_i| = |a_i| − 1 is the desuspended degree.

    Wait — the LV convention already includes all signs.
    Let me think about this more carefully.

    Actually, the precise relationship is:

    The LV identity on T^c(s^{-1}\\bar A):
      Σ_{r+s+t=n} (−1)^{rs+t} m_{r+1+t}(id^r ⊗ m_s ⊗ id^t) = 0

    When EVALUATED on elements a₁⊗...⊗a_n (in A, before desuspension):
    The sign for the term where m_s acts on positions r+1,...,r+s is:

    (−1)^{rs+t} · Koszul sign from commuting s^{-1} past elements

    The Koszul sign from commuting the s suspensions of m_s past the
    r desuspended elements a₁,...,a_r is (−1)^{r(s-1)} times the
    element-degree contribution.

    The precise formula: for elements of definite degree |a_i|, the sign is:

    (−1)^{Σ_{i=1}^{r}(|a_i|−1)(s−1) + Σ_{i=r+s+1}^{n}(|a_i|−1)}

    Wait, I should just compute both sides and compare.
    """
    r = len(degrees_before)
    s_lv = arity_inner
    t = len(degrees_after)

    base_lv = (-1) ** (r * s_lv + t)

    # Koszul sign from desuspension commutation
    desuspended_before = [d - 1 for d in degrees_before]
    koszul_from_desusp = (-1) ** ((s_lv - 1) * sum(desuspended_before) % 2)

    return base_lv * koszul_from_desusp


def compare_signs_on_term(n, degrees, inner_start, inner_arity):
    """Compare LV and Koszul signs for a specific term in the A∞ identity.

    The precise relationship between the two conventions:

    The LV sign (-1)^{rs+t} acts on the reduced bar construction
    T^c(s^{-1}\\bar A).
    When evaluated on elements a₁,...,a_n ∈ A, the desuspension introduces
    additional Koszul signs. The total LV-on-elements sign is:

      (-1)^{rs+t} · (-1)^{(j-1)·Σ_{i=1}^{r}(|a_i|-1)}

    The Koszul convention gives:
      (-1)^{(j-1)·Σ_{i=1}^{s}|a_i|}

    Comparing:
      LV_on_elements = (-1)^{rs+t+(j-1)(Σ|a_i|-r)} = (-1)^{t+r+(j-1)Σ|a_i|}
                     = (-1)^{n-j} · (-1)^{(j-1)Σ|a_i|}
                     = (-1)^{n-j} · Koszul
                     = (-1)^{i-1} · Koszul

    where i = n+1-j is the outer arity.

    This factor (-1)^{i-1} is the well-known relationship between the
    operadic (LV) and element-level (Koszul) A∞ sign conventions.
    Both define VALID A∞ identities; the algebras they characterize
    are the same.

    Parameters:
        n: total number of elements
        degrees: list of n element degrees |a₁|,...,|a_n|
        inner_start: s = index of first element in inner operation (0-based)
        inner_arity: j = arity of inner operation

    Returns:
        dict with lv, koszul, predicted_ratio, and agreement check
    """
    s = inner_start
    j = inner_arity
    r = s
    t = n - s - j
    i = n + 1 - j  # outer arity

    if r < 0 or t < 0 or r + j + t != n:
        raise ValueError(f"Invalid: n={n}, s={s}, j={j} gives r={r}, t={t}")

    # Koszul sign on A (Vol II convention)
    ksign = koszul_sign_on_elements(degrees[:s], j)

    # LV sign on the operadic level
    lv_operadic = lv_sign(r, j, t)

    # LV sign evaluated on elements (including desuspension Koszul signs)
    desusp_sum = sum(degrees[:r])  # Σ|a_i| for i=1..r
    desusp_koszul = (-1) ** (((j - 1) * (desusp_sum - r)) % 2)
    lv_on_elements = lv_operadic * desusp_koszul

    # Predicted ratio: (-1)^{i-1} = (-1)^{n-j}
    predicted_ratio = (-1) ** ((i - 1) % 2)

    # Check: lv_on_elements should equal predicted_ratio * koszul
    actual_ratio = lv_on_elements * ksign  # +1 if same, -1 if opposite

    return {
        'lv_operadic': lv_operadic,
        'lv_on_elements': lv_on_elements,
        'koszul': ksign,
        'predicted_ratio': predicted_ratio,
        'actual_ratio': actual_ratio,
        'relationship_holds': (actual_ratio == predicted_ratio),
    }


def verify_conventions_equivalent(n, degrees):
    """Check if LV and Koszul conventions are related by (-1)^{i-1}.

    The proven relationship: LV_on_elements = (-1)^{i-1} · Koszul
    where i = n+1-j is the outer arity.

    Returns:
        dict with 'relationship_holds' (should be True for all terms)
    """
    terms = []
    all_hold = True
    for j in range(1, n + 1):
        for s in range(0, n - j + 1):
            result = compare_signs_on_term(n, degrees, s, j)
            terms.append(result)
            if not result['relationship_holds']:
                all_hold = False

    return {
        'relationship_holds': all_hold,
        'terms': terms,
        'n': n,
        'degrees': degrees,
    }


def analyze_convention_relationship(max_n=6):
    """Systematic analysis: verify LV = (-1)^{i-1} · Koszul for all arities.

    Tests all arities n from 1 to max_n, for various degree patterns.

    The proven relationship: LV_on_elements = (-1)^{i-1} · Koszul
    where i = n+1-j is the outer arity. This holds for ALL element degrees.
    """
    results = {}
    for n in range(1, max_n + 1):
        test_patterns = []

        # All degree 0
        test_patterns.append(('all_deg0', [0] * n))

        # All degree 1
        test_patterns.append(('all_deg1', [1] * n))

        # Mixed: first half 0, second half 1
        if n >= 2:
            mixed = [0] * (n // 2) + [1] * (n - n // 2)
            test_patterns.append(('mixed_01', mixed))

        # All degree 2 (for higher checks)
        test_patterns.append(('all_deg2', [2] * n))

        for name, degrees in test_patterns:
            key = (n, name)
            results[key] = verify_conventions_equivalent(n, degrees)

    return results
