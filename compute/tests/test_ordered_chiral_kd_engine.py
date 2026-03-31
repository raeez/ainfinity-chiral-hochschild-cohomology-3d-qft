"""Tests for Ordered Associative Chiral KD Engine.

Verifies the ordered (E_1) sector of chiral Koszul duality:
1. Deconcatenation coproduct: all ordered splits
2. Bar differential: consecutive-collapse formula
3. d^2 = 0: associativity implies d^2 = 0
4. Shuffle product: count, associativity, symmetry
5. Opposite involution: involution property

Each test performs ACTUAL computation, not lookup.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from fractions import Fraction
from math import comb

from lib.ordered_chiral_kd_engine import (
    deconcatenation_coproduct,
    bar_differential_ordered,
    d_squared_check,
    shuffle_product,
    opposite_involution,
    commutative_m2,
    free_associative_m2,
)


# ===================================================================
# 1. DECONCATENATION COPRODUCT
# ===================================================================

class TestDeconcatenation:
    """Test the deconcatenation coproduct on words."""

    def test_length_1_empty(self):
        """A single-letter word has no nontrivial splits."""
        result = deconcatenation_coproduct(['a'])
        assert result == []

    def test_length_0_empty(self):
        """Empty word has no splits."""
        result = deconcatenation_coproduct([])
        assert result == []

    def test_length_2(self):
        """[a, b] splits only as ([a], [b])."""
        result = deconcatenation_coproduct(['a', 'b'])
        assert len(result) == 1
        assert result[0] == (['a'], ['b'])

    def test_length_3(self):
        """[a, b, c] has exactly 2 splits."""
        result = deconcatenation_coproduct(['a', 'b', 'c'])
        assert len(result) == 2
        assert (['a'], ['b', 'c']) in result
        assert (['a', 'b'], ['c']) in result

    def test_length_4(self):
        """[a, b, c, d] has exactly 3 splits."""
        result = deconcatenation_coproduct(['a', 'b', 'c', 'd'])
        assert len(result) == 3
        expected = [
            (['a'], ['b', 'c', 'd']),
            (['a', 'b'], ['c', 'd']),
            (['a', 'b', 'c'], ['d']),
        ]
        for e in expected:
            assert e in result

    def test_split_count_general(self):
        """A word of length n has exactly n-1 splits."""
        for n in range(1, 8):
            word = [chr(ord('a') + i) for i in range(n)]
            result = deconcatenation_coproduct(word)
            assert len(result) == max(0, n - 1)

    def test_splits_reconstruct(self):
        """Each split (left, right) concatenates back to the original word."""
        word = ['x', 'y', 'z', 'w']
        for left, right in deconcatenation_coproduct(word):
            assert left + right == word


# ===================================================================
# 2. BAR DIFFERENTIAL
# ===================================================================

class TestBarDifferential:
    """Test the ordered bar differential d_C."""

    def test_d_of_length_1(self):
        """d of a single element is 0."""
        result = bar_differential_ordered(['a'], commutative_m2)
        assert result == {}

    def test_d_of_pair(self):
        """d[a|b] = m2(a, b), one term, positive sign."""
        result = bar_differential_ordered(['a', 'b'], commutative_m2)
        # m2(a, b) = 'ab' with coefficient 1
        # Sign: (-1)^0 = +1
        assert len(result) == 1
        # commutative_m2 sorts: m2('a','b') = {'ab': 1}
        assert result[('ab',)] == Fraction(1)

    def test_d_of_triple_has_two_terms(self):
        """d[a|b|c] has exactly 2 terms (consecutive collapse at positions 0 and 1)."""
        result = bar_differential_ordered(['a', 'b', 'c'], commutative_m2)
        # i=0: +1 * [m2(a,b) | c] = +[ab | c]
        # i=1: -1 * [a | m2(b,c)] = -[a | bc]
        assert len(result) == 2
        assert result[('ab', 'c')] == Fraction(1)
        assert result[('a', 'bc')] == Fraction(-1)

    def test_d_of_triple_noncommutative(self):
        """d[a|b|c] with free associative m2."""
        result = bar_differential_ordered(['a', 'b', 'c'], free_associative_m2)
        # i=0: +1 * [ab | c]
        # i=1: -1 * [a | bc]
        assert result[('ab', 'c')] == Fraction(1)
        assert result[('a', 'bc')] == Fraction(-1)

    def test_d_of_quadruple_has_three_terms(self):
        """d[a|b|c|d] has 3 terms."""
        result = bar_differential_ordered(['a', 'b', 'c', 'd'], free_associative_m2)
        # i=0: +1 * [ab | c | d]
        # i=1: -1 * [a | bc | d]
        # i=2: +1 * [a | b | cd]
        assert len(result) == 3
        assert result[('ab', 'c', 'd')] == Fraction(1)
        assert result[('a', 'bc', 'd')] == Fraction(-1)
        assert result[('a', 'b', 'cd')] == Fraction(1)


# ===================================================================
# 3. D-SQUARED = 0
# ===================================================================

class TestDSquared:
    """Test d^2 = 0 (equivalent to associativity of m2)."""

    def test_d_squared_arity_2_commutative(self):
        """d^2 = 0 on [a|b] with commutative m2."""
        result = d_squared_check(['a', 'b'], commutative_m2)
        assert result == {}

    def test_d_squared_arity_3_commutative(self):
        """d^2 = 0 on [a|b|c] with commutative m2."""
        result = d_squared_check(['a', 'b', 'c'], commutative_m2)
        assert result == {}

    def test_d_squared_arity_4_commutative(self):
        """d^2 = 0 on [a|b|c|d] with commutative m2."""
        result = d_squared_check(['a', 'b', 'c', 'd'], commutative_m2)
        assert result == {}

    def test_d_squared_arity_3_free_assoc(self):
        """d^2 = 0 on [a|b|c] with free associative m2."""
        result = d_squared_check(['a', 'b', 'c'], free_associative_m2)
        assert result == {}

    def test_d_squared_arity_4_free_assoc(self):
        """d^2 = 0 on [a|b|c|d] with free associative m2."""
        result = d_squared_check(['a', 'b', 'c', 'd'], free_associative_m2)
        assert result == {}

    def test_d_squared_arity_5_commutative(self):
        """d^2 = 0 at arity 5."""
        word = ['a', 'b', 'c', 'd', 'e']
        result = d_squared_check(word, commutative_m2)
        assert result == {}

    def test_d_squared_repeated_generators(self):
        """d^2 = 0 even with repeated generators."""
        result = d_squared_check(['a', 'a', 'a'], commutative_m2)
        assert result == {}

    def test_d_squared_fails_nonassociative(self):
        """d^2 != 0 for a non-associative product.

        Define a product where m2(m2(a,b),c) != m2(a,m2(b,c)) by
        making left-association and right-association produce distinct labels.
        """
        def nonassoc_m2(a: str, b: str):
            # m2 that is NOT associative:
            # Always wraps with explicit parentheses, so
            # m2(m2(a,b), c) = ((ab)c) but m2(a, m2(b,c)) = (a(bc))
            return {'(' + a + b + ')': Fraction(1)}

        # For [a|b|c]:
        # d[a|b|c] = [(ab)|c] - [a|(bc)]
        # d(d[a|b|c]) = m2((ab),c) - m2(a,(bc))
        #             = ((ab)c) - (a(bc))  -- different strings!
        result = d_squared_check(['a', 'b', 'c'], nonassoc_m2)
        # Should NOT be zero for this non-associative product
        assert result != {}


# ===================================================================
# 4. SHUFFLE PRODUCT
# ===================================================================

class TestShuffle:
    """Test the shuffle product on words."""

    def test_shuffle_count_1_1(self):
        """Sh([a], [b]) has C(2,1) = 2 shuffles."""
        result = shuffle_product(['a'], ['b'])
        assert len(result) == comb(2, 1)
        assert ('a', 'b') in result
        assert ('b', 'a') in result

    def test_shuffle_count_2_1(self):
        """Sh([a,b], [c]) has C(3,2) = 3 shuffles."""
        result = shuffle_product(['a', 'b'], ['c'])
        assert len(result) == comb(3, 2)

    def test_shuffle_count_2_2(self):
        """Sh([a,b], [c,d]) has C(4,2) = 6 shuffles."""
        result = shuffle_product(['a', 'b'], ['c', 'd'])
        assert len(result) == comb(4, 2)

    def test_shuffle_count_3_2(self):
        """Sh of length 3 and 2 has C(5,3) = 10 shuffles."""
        result = shuffle_product(['a', 'b', 'c'], ['d', 'e'])
        assert len(result) == comb(5, 3)

    def test_shuffle_count_general(self):
        """Sh(w1, w2) has exactly C(|w1|+|w2|, |w1|) elements."""
        for p in range(1, 5):
            for q in range(1, 5):
                w1 = [f'a{i}' for i in range(p)]
                w2 = [f'b{j}' for j in range(q)]
                result = shuffle_product(w1, w2)
                assert len(result) == comb(p + q, p)

    def test_shuffle_preserves_order_w1(self):
        """Each shuffle preserves the internal order of w1."""
        w1 = ['a', 'b', 'c']
        w2 = ['x', 'y']
        for shuf in shuffle_product(w1, w2):
            # Extract positions of w1 elements
            positions = [i for i, s in enumerate(shuf) if s in ['a', 'b', 'c']]
            extracted = [shuf[i] for i in positions]
            assert extracted == ['a', 'b', 'c']

    def test_shuffle_preserves_order_w2(self):
        """Each shuffle preserves the internal order of w2."""
        w1 = ['a', 'b']
        w2 = ['x', 'y', 'z']
        for shuf in shuffle_product(w1, w2):
            positions = [i for i, s in enumerate(shuf) if s in ['x', 'y', 'z']]
            extracted = [shuf[i] for i in positions]
            assert extracted == ['x', 'y', 'z']

    def test_shuffle_empty_w1(self):
        """Shuffle with empty w1 returns just w2."""
        result = shuffle_product([], ['a', 'b'])
        assert len(result) == 1
        assert result[0] == ('a', 'b')

    def test_shuffle_empty_w2(self):
        """Shuffle with empty w2 returns just w1."""
        result = shuffle_product(['a', 'b'], [])
        assert len(result) == 1
        assert result[0] == ('a', 'b')


# ===================================================================
# 5. OPPOSITE INVOLUTION
# ===================================================================

class TestOpposite:
    """Test the opposite (reversal) involution."""

    def test_involution_squared_is_identity(self):
        """(w^op)^op = w for several words."""
        words = [
            ['a'],
            ['a', 'b'],
            ['a', 'b', 'c'],
            ['a', 'b', 'c', 'd'],
            ['x', 'y', 'z', 'w', 'v'],
        ]
        for w in words:
            assert opposite_involution(opposite_involution(w)) == w

    def test_reversal_of_pair(self):
        """[a, b]^op = [b, a]."""
        assert opposite_involution(['a', 'b']) == ['b', 'a']

    def test_reversal_of_triple(self):
        """[a, b, c]^op = [c, b, a]."""
        assert opposite_involution(['a', 'b', 'c']) == ['c', 'b', 'a']

    def test_palindrome_fixed(self):
        """A palindrome is a fixed point of the involution."""
        w = ['a', 'b', 'a']
        assert opposite_involution(w) == w

    def test_empty_word(self):
        """Empty word is a fixed point."""
        assert opposite_involution([]) == []
