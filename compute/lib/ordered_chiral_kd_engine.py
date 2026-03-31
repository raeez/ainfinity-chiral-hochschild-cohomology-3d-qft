"""Ordered Associative Chiral KD Engine: ordered bar complex computations.

Implements the ordered (E_1 / associative) sector of chiral Koszul duality.
The bar complex for an associative algebra uses ordered tensor products
(words) with the bar differential collapsing consecutive pairs via m_2.

The key mathematical objects:

1. **Deconcatenation coproduct**: For a word [a1, ..., an], all ordered
   splits (left, right) where left + right = word and both nonempty.
   This is the R-direction (topological, E_1) coproduct.

2. **Bar differential (ordered)**: d[a1|...|an] = sum +/- [a1|...|m2(ai,ai+1)|...|an]
   for consecutive pairs. This is the C-direction differential restricted
   to the ordered sector.

3. **Shuffle product**: All interleavings of two words preserving internal
   order. This gives the commutative product on the bar complex.

4. **Opposite involution**: Reversal w -> w^op. An involution on the
   bar complex interchanging A and A^op.

References:
  Vol II: chapters/connections/ordered_associative_chiral_kd_core.tex
  Vol I: bar_construction.tex, cobar_construction.tex
  Loday-Vallette (2012): Algebraic Operads, Chapter 2
"""
from __future__ import annotations

from itertools import combinations
from math import comb, factorial
from typing import Callable, Dict, List, Optional, Tuple

from fractions import Fraction


# =========================================================================
# 1. DECONCATENATION COPRODUCT
# =========================================================================

def deconcatenation_coproduct(word: List[str]) -> List[Tuple[List[str], List[str]]]:
    """Compute all ordered splits of a word into two nonempty parts.

    For a word [a1, a2, ..., an], returns all pairs (left, right) where
    left = [a1, ..., ai] and right = [a_{i+1}, ..., an] for 1 <= i < n.

    This is the deconcatenation (or deshuffle) coproduct of the tensor
    coalgebra, which is the R-direction (topological) coproduct of the
    bar complex.

    Args:
        word: A list of generator labels.

    Returns:
        List of (left, right) pairs. Empty if len(word) < 2.
    """
    n = len(word)
    if n < 2:
        return []
    splits = []
    for i in range(1, n):
        left = word[:i]
        right = word[i:]
        splits.append((left, right))
    return splits


# =========================================================================
# 2. BAR DIFFERENTIAL (ORDERED / CONSECUTIVE COLLAPSE)
# =========================================================================

def bar_differential_ordered(
    word: List[str],
    m2_func: Callable[[str, str], Dict[str, Fraction]],
) -> Dict[Tuple[str, ...], Fraction]:
    """Apply the ordered bar differential to a word.

    The bar differential for the ordered (associative) bar complex is:
        d[a1|...|an] = sum_{i=1}^{n-1} (-1)^{i-1} [a1|...|m2(ai,ai+1)|...|an]

    where m2 is the binary product. The sign (-1)^{i-1} is the Koszul sign
    from the desuspension convention (cohomological grading, |s^{-1}| = -1).

    Args:
        word: List of generator labels [a1, ..., an].
        m2_func: Binary product m2(a, b) -> {result_label: coefficient}.
            Returns a dict representing a linear combination.

    Returns:
        Dict mapping result words (as tuples) to coefficients (Fraction).
        Represents a linear combination of bar elements.
    """
    n = len(word)
    if n < 2:
        return {}

    result: Dict[Tuple[str, ...], Fraction] = {}

    for i in range(n - 1):
        # Collapse word[i] and word[i+1] via m2
        sign = Fraction((-1) ** i)
        products = m2_func(word[i], word[i + 1])

        for prod_label, prod_coeff in products.items():
            # Build the new word: word[:i] + [prod_label] + word[i+2:]
            new_word = tuple(word[:i]) + (prod_label,) + tuple(word[i + 2:])
            coeff = sign * prod_coeff
            result[new_word] = result.get(new_word, Fraction(0)) + coeff

    # Remove zero entries
    return {k: v for k, v in result.items() if v != Fraction(0)}


def _apply_differential(
    element: Dict[Tuple[str, ...], Fraction],
    m2_func: Callable[[str, str], Dict[str, Fraction]],
) -> Dict[Tuple[str, ...], Fraction]:
    """Apply the bar differential to a linear combination of words.

    Args:
        element: Dict mapping words (tuples) to coefficients.
        m2_func: Binary product function.

    Returns:
        Dict mapping result words to coefficients.
    """
    result: Dict[Tuple[str, ...], Fraction] = {}
    for word_tuple, word_coeff in element.items():
        d_word = bar_differential_ordered(list(word_tuple), m2_func)
        for new_word, new_coeff in d_word.items():
            total = word_coeff * new_coeff
            result[new_word] = result.get(new_word, Fraction(0)) + total
    return {k: v for k, v in result.items() if v != Fraction(0)}


# =========================================================================
# 3. D-SQUARED CHECK
# =========================================================================

def d_squared_check(
    word: List[str],
    m2_func: Callable[[str, str], Dict[str, Fraction]],
) -> Dict[Tuple[str, ...], Fraction]:
    """Verify d^2 = 0 on a given word.

    Computes d(d(word)) and returns the result. For an associative m2,
    the result should be identically zero (d^2 = 0 is equivalent to
    associativity of m2).

    Args:
        word: List of generator labels.
        m2_func: Associative binary product.

    Returns:
        Dict mapping result words to coefficients. Should be empty ({})
        if m2 is associative.
    """
    # First application: d(word)
    d_once = bar_differential_ordered(word, m2_func)

    if not d_once:
        return {}

    # Second application: d(d(word))
    return _apply_differential(d_once, m2_func)


# =========================================================================
# 4. SHUFFLE PRODUCT
# =========================================================================

def shuffle_product(w1: List[str], w2: List[str]) -> List[Tuple[str, ...]]:
    """Compute all (p,q)-shuffles of two words.

    A (p,q)-shuffle of words w1 = [a1,...,ap] and w2 = [b1,...,bq] is a
    permutation sigma of [a1,...,ap,b1,...,bq] such that the relative
    order of the a's and the b's is preserved.

    The number of shuffles is C(p+q, p) = (p+q)! / (p! q!).

    Args:
        w1: First word (list of labels).
        w2: Second word (list of labels).

    Returns:
        List of all shuffled words (as tuples), preserving internal order
        of both w1 and w2.
    """
    p = len(w1)
    q = len(w2)
    n = p + q

    if p == 0:
        return [tuple(w2)] if q > 0 else [()]
    if q == 0:
        return [tuple(w1)]

    result = []

    # Choose which positions in the merged word of length n hold elements of w1
    for positions in combinations(range(n), p):
        merged = [''] * n
        pos_set = set(positions)

        # Place w1 elements at chosen positions
        w1_idx = 0
        for pos in positions:
            merged[pos] = w1[w1_idx]
            w1_idx += 1

        # Place w2 elements in remaining positions
        w2_idx = 0
        for pos in range(n):
            if pos not in pos_set:
                merged[pos] = w2[w2_idx]
                w2_idx += 1

        result.append(tuple(merged))

    return result


# =========================================================================
# 5. OPPOSITE INVOLUTION
# =========================================================================

def opposite_involution(word: List[str]) -> List[str]:
    """Apply the opposite involution: reverse the word.

    For an associative algebra A, the opposite algebra A^op has the
    same underlying vector space but reversed multiplication:
        m2^{op}(a, b) = m2(b, a).

    On the bar complex, the opposite involution acts by reversing the
    order of the tensor factors:
        [a1|a2|...|an]^op = [an|...|a2|a1].

    This is an involution: (w^op)^op = w.

    Args:
        word: List of generator labels.

    Returns:
        Reversed word.
    """
    return list(reversed(word))


# =========================================================================
# 6. STANDARD EXAMPLE: COMMUTATIVE m2 (POLYNOMIAL RING)
# =========================================================================

def commutative_m2(a: str, b: str) -> Dict[str, Fraction]:
    """Commutative multiplication: m2(a,b) = a*b (concatenated label).

    This is the simplest associative product for testing d^2=0.
    The product is commutative: m2(a,b) = m2(b,a), and associative:
    m2(m2(a,b),c) = m2(a,m2(b,c)).

    Args:
        a: First generator label.
        b: Second generator label.

    Returns:
        Dict with single entry: concatenated label -> 1.
    """
    # Canonical ordering for commutativity
    label = ''.join(sorted([a, b]))
    return {label: Fraction(1)}


def free_associative_m2(a: str, b: str) -> Dict[str, Fraction]:
    """Free associative multiplication: m2(a,b) = ab (ordered concatenation).

    Non-commutative but associative. For testing d^2=0 in the
    non-commutative setting.

    Args:
        a: First generator label.
        b: Second generator label.

    Returns:
        Dict with single entry: ordered concatenation -> 1.
    """
    return {a + b: Fraction(1)}
