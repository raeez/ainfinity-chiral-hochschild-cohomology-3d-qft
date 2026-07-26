"""Tests for the boundary-linear Kuranishi reduction formulas."""
import os
import sys
from fractions import Fraction

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.boundary_linear_kuranishi import (
    derived_critical_dga_profile,
    kappa2,
    kappa3,
    massive_i2,
    massive_i3,
    minimal_line_operation_profile,
)


def _mul(scale):
    return lambda x, y: Fraction(scale) * Fraction(x) * Fraction(y)


def _tri(scale):
    return lambda x, y, z: Fraction(scale) * Fraction(x) * Fraction(y) * Fraction(z)


def test_massive_i2_is_negative_schur_complement():
    F2B = _mul(6)
    A_inv = lambda x: Fraction(x, 3)

    assert massive_i2(F2B, A_inv, 2, 5) == Fraction(-20)


def test_massive_i3_includes_cyclic_exchange():
    F2B = _mul(2)
    F3B = _tri(5)
    A_inv = lambda x: Fraction(x, 4)

    t1, t2, t3 = 2, 3, 5
    i2_35 = massive_i2(F2B, A_inv, 3, 5)
    i2_52 = massive_i2(F2B, A_inv, 5, 2)
    i2_23 = massive_i2(F2B, A_inv, 2, 3)
    expected_exchange = (
        F2B(2, i2_35)
        + F2B(3, i2_52)
        + F2B(5, i2_23)
    )
    expected = -A_inv(F3B(t1, t2, t3) + expected_exchange)

    assert massive_i3(F2B, F3B, A_inv, t1, t2, t3) == expected


def test_kappa2_and_kappa3_direct_minus_exchange():
    F2B = _mul(3)
    F2C = _mul(7)
    F3C = _tri(11)
    A_inv = lambda x: Fraction(x, 5)

    t1, t2, t3 = 2, 3, 4
    assert kappa2(F2C, t1, t2) == Fraction(42)

    exchange = (
        F2C(2, A_inv(F2B(3, 4)))
        + F2C(3, A_inv(F2B(4, 2)))
        + F2C(4, A_inv(F2B(2, 3)))
    )
    assert kappa3(F2B, F2C, F3C, A_inv, t1, t2, t3) == F3C(t1, t2, t3) - exchange


def test_minimal_line_operation_normalizations():
    m2 = minimal_line_operation_profile(2)
    m3 = minimal_line_operation_profile(3)
    m5 = minimal_line_operation_profile(5)

    assert m2["factor"] == Fraction(1, 2)
    assert m2["has_exterior_product"] is True
    assert "lambda_t1 wedge lambda_t2" in m2["formula"]

    assert m3["factor"] == Fraction(1, 6)
    assert m3["has_exterior_product"] is False
    assert "1/6 kappa_3" in m3["formula"]

    assert m5["factor"] == Fraction(1, 120)
    assert "1/120 kappa_5" in m5["formula"]


def test_minimal_line_operation_rejects_unary():
    with pytest.raises(ValueError):
        minimal_line_operation_profile(1)


def test_derived_critical_dga_profile_is_bulk_not_line():
    profile = derived_critical_dga_profile()

    assert profile["potential"] == "W_eff(u,gamma)=<gamma,kappa(u)>"
    assert profile["algebra"] == "k[[u,gamma]] tensor Lambda(xi,chi)"
    assert "dW/du_i d/dxi_i" in profile["differential"]
    assert "dW/dgamma_a d/dchi_a" in profile["differential"]
    assert "HH(K_F,p)" in profile["bulk_line_chain"]
    assert "O(dCrit(W_eff))" in profile["bulk_line_chain"]
    assert profile["line_not_bulk"] == "bulk is Hochschild/derived centre of the line algebra"


def test_active_core_source_contains_kuranishi_skeleton():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    source = os.path.join(repo_root, 'chapters', 'connections', 'ht_bulk_boundary_line_core.tex')
    with open(source, encoding='utf-8') as handle:
        text = handle.read()

    assert r'\frac{1}{q!}F_q^I' in text
    assert 'eq:massive-field-i2' in text
    assert 'eq:massive-field-i3' in text
    assert r'F_3^C(t_1,t_2,t_3)' in text
    assert r'F_2^C\Bigl(t_1,\,A^{-1}F_2^I(t_2,t_3)\Bigr)' in text
    assert 'eq:weff-derived-critical-dga' in text
    assert 'eq:weff-derived-critical-differential' in text
