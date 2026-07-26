"""BCOV degree-zero Calabi--Yau threefold constants.

This test separates the degree-zero BCOV/Gromov--Witten invariant from
the marked one-point Faber--Pandharipande scalar used in the shadow
obstruction tower.
"""

from __future__ import annotations

import unittest
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial


@lru_cache(maxsize=None)
def bernoulli(n: int) -> Fraction:
    if n == 0:
        return Fraction(1)
    if n == 1:
        return Fraction(-1, 2)
    if n > 1 and n % 2 == 1:
        return Fraction(0)

    total = Fraction(0)
    for k in range(n):
        total += Fraction(comb(n + 1, k)) * bernoulli(k)
    return -total / Fraction(n + 1)


def unmarked_cy3_hodge_integral(g: int) -> Fraction:
    if g < 2:
        raise ValueError("the unmarked CY3 Hodge integral requires g >= 2")
    return (
        abs(bernoulli(2 * g))
        / Fraction(2 * g)
        * abs(bernoulli(2 * g - 2))
        / Fraction(2 * g - 2)
        / Fraction(factorial(2 * g - 2))
    )


def bcov_degree_zero_free_energy(chi: int, g: int) -> Fraction:
    if g == 1:
        return Fraction(chi, 576)
    if g < 1:
        raise ValueError("BCOV genus must be positive")
    return Fraction((-1) ** g) * Fraction(chi, 2) * unmarked_cy3_hodge_integral(g)


def marked_lambda_fp(g: int) -> Fraction:
    if g < 1:
        raise ValueError("lambda_g^FP requires g >= 1")
    return (
        Fraction(2 ** (2 * g - 1) - 1, 2 ** (2 * g - 1))
        * abs(bernoulli(2 * g))
        / Fraction(factorial(2 * g))
    )


class TestBCOVDegreeZeroConstants(unittest.TestCase):
    def test_quintic_genus_one_constant(self):
        self.assertEqual(bcov_degree_zero_free_energy(-200, 1), Fraction(-25, 72))

    def test_quintic_genus_two_constant(self):
        self.assertEqual(unmarked_cy3_hodge_integral(2), Fraction(1, 2880))
        self.assertEqual(bcov_degree_zero_free_energy(-200, 2), Fraction(-5, 144))

    def test_marked_fp_scalar_is_different_from_bcov_degree_zero(self):
        kappa_bcov_quintic = Fraction(-200, 24)
        collapsed = kappa_bcov_quintic * marked_lambda_fp(2)
        self.assertEqual(collapsed, Fraction(-35, 3456))
        self.assertNotEqual(collapsed, bcov_degree_zero_free_energy(-200, 2))


if __name__ == "__main__":
    unittest.main()
