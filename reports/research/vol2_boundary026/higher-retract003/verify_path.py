"""Exact new path and higher-coefficient identities, with all tensor signs."""
from pathlib import Path
import json
import hashlib

helper = Path('reports/research/vol2_boundary026/verify_quantum_contraction.py')
exec(helper.read_text().split('counts = []')[0])

# Index zero is the path direction: y_0=tau and theta_0=d tau.
# Indices 1,2,3 are the physical compact-current generators.
N = 4
ZERO = (0,) * N


def basis(delta=ZERO, y=ZERO, mask=0, hbar=0):
    return {(tuple(delta), tuple(y), mask, hbar): F(1)}


def d0(p):
    return add(*(multiply_delta(odd_derivative(p, i), i) for i in (1, 2, 3)))


def even_y_derivative(p, i):
    out = defaultdict(F)
    for (delta, y, mask, hbar), coefficient in p.items():
        if y[i]:
            new = list(y)
            new[i] -= 1
            out[(delta, tuple(new), mask, hbar)] += y[i] * coefficient
    return tidy(out)


covariance = {(1, 2): F(1), (1, 3): F(1, 9), (2, 3): F(1, 4)}
K, D, exponential, projection, H = operators(covariance)


def total(p):
    return add(D(p), multiply_odd(even_y_derivative(p, 0), 0))


tau = variable(0, 'y')
dtau = basis(mask=1)
theta = {i: basis(mask=1 << i) for i in (1, 2, 3)}
delta = {i: variable(i, 'delta') for i in (1, 2, 3)}
y = {i: variable(i, 'y') for i in (1, 2, 3)}
x = {i: add(delta[i], y[i]) for i in (1, 2, 3)}


def P(p):
    out = {}
    for (da, ya, mask, hbar), coefficient in p.items():
        assert da[0] == ya[0] == 0 and not mask & 1
        term = scale(basis(y=ya, hbar=hbar), coefficient)
        for i in (1, 2, 3):
            replacement = add(multiply(tau, delta[i]), multiply(dtau, theta[i]))
            for _ in range(da[i]):
                term = multiply(term, replacement)
        for i in (1, 2, 3):
            if mask & (1 << i):
                term = multiply(term, multiply(tau, theta[i]))
        out = add(out, term)
    return out


def T(p):
    return exponential(P(exponential(p, -1)), 1)


def parity(p):
    parities = {mask.bit_count() % 2 for _, _, mask, _ in p}
    assert len(parities) <= 1
    return next(iter(parities), 0)


def bracket(p, q, differential):
    return add(differential(multiply(p, q)),
               scale(multiply(differential(p), q), -1),
               scale(multiply(p, differential(q)), -(-1) ** parity(p)))


def G2(a, b):
    return add(T(multiply(a, b)), scale(multiply(T(a), T(b)), -1))


def G3(a, b, c):
    return add(T(multiply(multiply(a, b), c)),
               scale(multiply(T(multiply(a, b)), T(c)), -1),
               scale(multiply(T(multiply(a, c)), T(b)),
                     -(-1) ** (parity(b) * parity(c))),
               scale(multiply(T(a), T(multiply(b, c))), -1),
               scale(multiply(multiply(T(a), T(b)), T(c)), 2))


def evaluate(p, endpoint):
    out = defaultdict(F)
    for (da, ya, mask, hbar), coefficient in p.items():
        if mask & 1:
            continue
        powers = list(ya)
        powers[0] = 0
        out[(da, tuple(powers), mask, hbar)] += coefficient * endpoint ** ya[0]
    return tidy(out)


def integrate(p):
    out = defaultdict(F)
    for (da, ya, mask, hbar), coefficient in p.items():
        if mask & 1:
            powers = list(ya)
            powers[0] = 0
            out[(da, tuple(powers), mask ^ 1, hbar)] += coefficient / (ya[0] + 1)
    return tidy(out)


samples = list(x.values()) + list(theta.values()) + [
    multiply(x[1], x[1]), multiply(theta[1], x[2]),
    multiply(theta[1], theta[2])]
binary_count = ternary_count = 0
for a in samples:
    assert total(T(a)) == T(D(a))
    assert evaluate(T(a), 0) == projection(a)
    assert evaluate(T(a), 1) == a
    assert integrate(T(a)) == H(a)
    for b in samples:
        lhs = add(total(G2(a, b)), bracket(T(a), T(b), total))
        rhs = add(T(bracket(a, b, D)), G2(D(a), b),
                  scale(G2(a, D(b)), (-1) ** parity(a)))
        assert lhs == rhs, ('binary', a, b)
        assert not evaluate(G2(a, b), 1)
        binary_count += 1
        for c in samples:
            lhs = add(total(G3(a, b, c)),
                      bracket(G2(a, b), T(c), total),
                      scale(bracket(G2(a, c), T(b), total),
                            (-1) ** (parity(b) * parity(c))),
                      bracket(T(a), G2(b, c), total))
            rhs = add(G3(D(a), b, c),
                      scale(G3(a, D(b), c), (-1) ** parity(a)),
                      scale(G3(a, b, D(c)), (-1) ** (parity(a) + parity(b))),
                      G2(bracket(a, b, D), c),
                      scale(G2(bracket(a, c, D), b),
                            (-1) ** (parity(b) * parity(c))),
                      scale(G2(bracket(b, c, D), a),
                            (-1) ** (parity(a) * (parity(b) + parity(c)))))
            assert lhs == rhs, ('ternary', a, b, c)
            ternary_count += 1

one_minus_tau2 = add(basis(), scale(multiply(tau, tau), -1))
expected2 = scale(hbar_times(one_minus_tau2), -covariance[1, 2])
assert G2(x[1], x[2]) == expected2
expected3 = scale(hbar_times(hbar_times(multiply(one_minus_tau2, one_minus_tau2))),
                  2 * covariance[1, 2] * covariance[1, 3])
assert G3(multiply(x[1], x[1]), x[2], x[3]) == expected3
assert total(expected2) == scale(multiply(hbar_times(tau), dtau), 2 * covariance[1, 2])
assert bracket(T(x[1]), T(x[2]), total) == scale(total(expected2), -1)

print(json.dumps({'result': 'PASS', 'python': platform.python_version(),
                  'arithmetic': 'exact rational polynomial including tau and odd d-tau',
                  'unary_examples': len(samples), 'binary_identities': binary_count,
                  'ternary_identities': ternary_count,
                  'extra_checks': ['endpoint maps', 'integrated unary path equals H',
                                   'printed binary and cubic path coefficients',
                                   'target bracket cancels path coefficient differential'],
                  'helper_sha256': hashlib.sha256(helper.read_bytes()).hexdigest(),
                  'scope': 'New path identities and higher coefficients only.'}, indent=2))
