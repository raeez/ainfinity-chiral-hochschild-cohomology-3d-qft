"""Exact finite checks for Dolbeault residue traces and current coefficients."""

from itertools import combinations
from math import factorial
import json
import platform
import sympy as sp

z = sp.symbols("z1 z2 z3")
zb = sp.symbols("zb1 zb2 zb3")
x, s, w, xb, sb, wb = sp.symbols("x s w xb sb wb")
lam, t = sp.symbols("lambda t")


def add_term(out, mask, value):
    out[mask] = sp.expand(out.get(mask, 0) + value)
    if out[mask] == 0:
        del out[mask]


def contract(poly, index):
    out = {}
    bit = 1 << index
    for mask, coeff in poly.items():
        if mask & bit:
            sign = (-1) ** ((mask & (bit - 1)).bit_count())
            add_term(out, mask ^ bit, sign * coeff)
    return out


def delta(poly, barred):
    out = {}
    for index, var in enumerate(barred):
        for mask, coeff in contract(poly, index).items():
            add_term(out, mask, -sp.diff(coeff, var))
    return out


def wedge(a, b):
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            if ma & mb:
                continue
            inversions = sum(
                1 for i in range(3) for j in range(3)
                if ma & (1 << i) and mb & (1 << j) and i > j
            )
            add_term(out, ma | mb, (-1) ** inversions * ca * cb)
    return out


def residue(poly, pair, jet_orders=(0, 0)):
    i, j = pair
    k = next(v for v in range(3) if v not in pair)
    substitutions = {
        z[i]: x + s, z[j]: s, z[k]: w,
        zb[i]: xb + sb, zb[j]: sb, zb[k]: wb,
    }
    vector_image = {i: {1: 1}, j: {1: -1, 2: 1}, k: {4: 1}}
    transformed = {}
    for mask, coeff in poly.items():
        term = {0: coeff.subs(substitutions, simultaneous=True)}
        for index in range(3):
            if mask & (1 << index):
                term = wedge(term, vector_image[index])
        for new_mask, value in term.items():
            add_term(transformed, new_mask, value)
    ni, nj = jet_orders
    out = {}
    for mask, coeff in contract(transformed, 0).items():
        value = (-1) ** ni * sp.diff(coeff, x, ni + nj + 1)
        value = value.subs({x: 0, xb: 0})
        add_term(out, mask >> 1, value)
    return out


def plus(a, b):
    out = dict(a)
    for mask, coeff in b.items():
        add_term(out, mask, coeff)
    return out


def run():
    polynomials = [
        1 + z[0] * zb[0] + z[1] ** 2 * zb[2],
        (z[0] - z[1]) ** 3 * (zb[0] + 2 * zb[1]) * zb[2],
        z[0] ** 4 * z[1] * (zb[0] * zb[1] + zb[2] ** 2),
        (z[1] - z[2]) ** 5 * (zb[0] + zb[1] + zb[2]) ** 2,
    ]
    checks = 0
    square_checks = 0
    for coeff in polynomials:
        for mask in range(8):
            alpha = {mask: coeff}
            assert not delta(delta(alpha, zb), zb)
            square_checks += 1
            for pair in combinations(range(3), 2):
                for jets in [(0, 0), (1, 0), (0, 1), (1, 2), (2, 1)]:
                    residual = plus(
                        delta(residue(alpha, pair, jets), (sb, wb)),
                        residue(delta(alpha, zb), pair, jets),
                    )
                    assert not residual, (pair, jets, mask, residual)
                    checks += 1

    jet_normalizations = 0
    for ni in range(5):
        for nj in range(5):
            n = ni + nj + 2
            analytic = sp.diff((z[0] - z[1]) ** -2, z[0], ni, z[1], nj)
            coefficient = (-1) ** ni * factorial(n - 1)
            assert sp.simplify(analytic - coefficient * (z[0] - z[1]) ** -n) == 0
            residue_coefficient = coefficient * (-1) ** (n - 1) / sp.Integer(factorial(n - 1))
            transpose_sign = (-1) ** (n - 1)
            assert residue_coefficient * transpose_sign == (-1) ** ni
            jet_normalizations += 1

    # A normal test germ x has residue pi. All bar derivatives have zero residue.
    alpha = {1: (z[0] - z[1]) * (1 + z[2] * zb[2])}
    normal_trace = residue(alpha, (0, 1))
    assert normal_trace == {0: 1 + w * wb}

    q1, q2, q3, q4, zz = sp.symbols("q1 q2 q3 q4 zz")
    creation = q1 + q2 * zz + q3 * zz**2 + q4 * zz**3
    left = sp.expand(creation**2 * q1 + 2 * lam * zz**-2 * creation).coeff(zz, 0)
    right = sp.expand(creation * q1**2 + 2 * lam * zz**-2 * q1).coeff(zz, 0)
    associator = sp.expand(left - right)
    assert associator == 2 * lam * q3

    # Endpoint evaluation of a total derivative kills the interval one-form.
    c = (t - t**2) / 2
    g = t**3 / 6 - t**2 / 4 + t / 12
    assert sp.diff(g, t) + c == sp.Rational(1, 12)
    assert all(f.subs(t, a) == 0 for f in (c, g) for a in (0, 1))
    return {
        "status": "PASS",
        "python": platform.python_version(), "sympy": sp.__version__,
        "test_chain_square_checks": square_checks,
        "residue_anticommutation_checks": checks,
        "jet_normalizations": jet_normalizations,
        "regular_associator": str(associator),
        "normal_test_residue_without_pi": str(normal_trace),
        "scope": "Exact finite polynomial-jet and exterior-sign checks. The manuscript proves the smooth-chain statements. No composite collision-stratum algebra is certified.",
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
