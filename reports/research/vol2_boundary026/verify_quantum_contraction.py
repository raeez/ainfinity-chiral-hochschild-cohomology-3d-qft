"""Exact finite checks for the compact quantum current complex."""

from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json
import platform

N = 3
ZERO = (0,) * N


def tidy(p):
    return {k: v for k, v in p.items() if v}


def add(*polys):
    out = defaultdict(F)
    for p in polys:
        for key, value in p.items():
            out[key] += value
    return tidy(out)


def scale(p, coefficient):
    return tidy({key: coefficient * value for key, value in p.items()})


def basis(delta=ZERO, y=ZERO, mask=0, hbar=0):
    return {(tuple(delta), tuple(y), mask, hbar): F(1)}


def delta_derivative(p, i):
    out = defaultdict(F)
    for (delta, y, mask, hbar), c in p.items():
        if delta[i]:
            new = list(delta)
            new[i] -= 1
            out[(tuple(new), y, mask, hbar)] += delta[i] * c
    return tidy(out)


def odd_derivative(p, i):
    out = defaultdict(F)
    for (delta, y, mask, hbar), c in p.items():
        if mask & (1 << i):
            sign = (-1) ** ((mask & ((1 << i) - 1)).bit_count())
            out[(delta, y, mask ^ (1 << i), hbar)] += sign * c
    return tidy(out)


def multiply_delta(p, i):
    out = {}
    for (delta, y, mask, hbar), c in p.items():
        new = list(delta)
        new[i] += 1
        out[(tuple(new), y, mask, hbar)] = c
    return out


def multiply_odd(p, i):
    out = {}
    for (delta, y, mask, hbar), c in p.items():
        if not mask & (1 << i):
            sign = (-1) ** ((mask & ((1 << i) - 1)).bit_count())
            out[(delta, y, mask | (1 << i), hbar)] = sign * c
    return out


def hbar_times(p):
    return {(d, y, m, h + 1): c for (d, y, m, h), c in p.items()}


def d0(p):
    return add(*(multiply_delta(odd_derivative(p, i), i) for i in range(N)))


def h0(p):
    out = {}
    for (delta, y, mask, hbar), c in p.items():
        q = sum(delta) + mask.bit_count()
        if q:
            one = {(delta, y, mask, hbar): c / q}
            out = add(out, *(multiply_odd(delta_derivative(one, i), i)
                             for i in range(N)))
    return out


def pi(p):
    return {key: c for key, c in p.items() if sum(key[0]) == 0 and key[2] == 0}


def operators(covariance):
    def K(p):
        return add(*(scale(delta_derivative(delta_derivative(p, i), j), c)
                     for (i, j), c in covariance.items()))

    def D(p):
        contractions = add(*(scale(add(
            delta_derivative(odd_derivative(p, i), j),
            delta_derivative(odd_derivative(p, j), i)), c)
            for (i, j), c in covariance.items()))
        return add(d0(p), hbar_times(contractions))

    def exponential(p, sign):
        out, term, order = p, p, 0
        while term:
            order += 1
            term = scale(hbar_times(K(term)), F(sign, order))
            out = add(out, term)
        return out

    def projection(p):
        return pi(exponential(p, -1))

    def H(p):
        return exponential(h0(exponential(p, -1)), 1)

    return K, D, exponential, projection, H


def multiply(p, q):
    out = defaultdict(F)
    for (da, ya, ma, ha), ca in p.items():
        for (db, yb, mb, hb), cb in q.items():
            if ma & mb:
                continue
            swaps = sum((mb & ((1 << i) - 1)).bit_count()
                        for i in range(N) if ma & (1 << i))
            key = (tuple(a + b for a, b in zip(da, db)),
                   tuple(a + b for a, b in zip(ya, yb)), ma | mb, ha + hb)
            out[key] += (-1) ** swaps * ca * cb
    return tidy(out)


def variable(i, kind):
    exponents = tuple(int(i == j) for j in range(N))
    return basis(delta=exponents) if kind == 'delta' else basis(y=exponents)


counts = []
for level in (F(0), F(1), F(-2)):
    # Undifferentiated jets at -1, 0, 2 with the displayed bulk level.
    covariance = {(0, 1): level, (0, 2): level / 9, (1, 2): level / 4}
    K, D, exponential, projection, H = operators(covariance)
    checked = 0
    for delta in product(range(7), repeat=N):
        if sum(delta) > 6:
            continue
        for mask in range(1 << N):
            p = basis(delta=delta, mask=mask)
            assert not D(D(p)), ('D square', level, delta, mask)
            assert exponential(D(p), -1) == d0(exponential(p, -1))
            assert add(D(H(p)), H(D(p))) == add(p, scale(projection(p), -1))
            assert not H(H(p)), ('H square', level, delta, mask)
            assert not projection(H(p))
            assert not projection(D(p))
            assert not H(projection(p))
            assert projection(projection(p)) == projection(p)
            checked += 1

    x = [add(variable(i, 'delta'), variable(i, 'y')) for i in range(N)]
    y = [variable(i, 'y') for i in range(N)]
    delta = [variable(i, 'delta') for i in range(N)]
    theta = [basis(mask=1 << i) for i in range(N)]
    binary = multiply(x[0], x[1])
    binary_H = scale(add(multiply(theta[0], add(x[1], y[1])),
                         multiply(theta[1], add(x[0], y[0]))), F(1, 2))
    assert H(binary) == binary_H
    cubic = multiply(binary, x[2])
    cubic_H = {}
    correction = {}
    for i in range(N):
        j, k = [a for a in range(N) if a != i]
        core = add(multiply(y[j], y[k]),
                   scale(add(multiply(delta[j], y[k]),
                             multiply(y[j], delta[k])), F(1, 2)),
                   scale(multiply(delta[j], delta[k]), F(1, 3)))
        correction = add(correction,
                         scale(hbar_times(theta[i]), -F(2, 3) * covariance[j, k]))
        cubic_H = add(cubic_H, multiply(theta[i], core))
    assert H(cubic) == add(cubic_H, correction)
    defect = add(D(cubic_H), scale(add(cubic, scale(projection(cubic), -1)), -1))
    assert defect == scale(D(correction), -1)
    if level:
        assert defect
    counts.append({'bulk_level': str(level), 'basis_monomials': checked,
                   'maximum_delta_degree': 6, 'odd_masks': 8})

print(json.dumps({'python': platform.python_version(), 'arithmetic': 'fractions.Fraction',
                  'result': 'PASS', 'cases': counts,
                  'identities': ['D^2=0', 'ND=d0N', 'DH+HD=1-iPi', 'H^2=0',
                                 'PiH=PiD=Hi=0', 'Pi^2=Pi',
                                 'printed binary and cubic formulas',
                                 'nonzero defect without cubic quantum correction'],
                  'scope': 'Three even pairs, three odd generators, exact finite degrees. '
                           'The manuscript proof establishes all finite lists and degrees.'}, indent=2))
