"""
Laplace transform bridge between λ-bracket and classical r-matrix.

The bridge relation (BR3 in the paper):
  r(z) = Res_{λ=0} e^{-λz} {a_λ b}

For POLYNOMIAL λ-brackets, this is a finite sum using:
  ∫_0^∞ λ^n e^{-λz} dλ = n! / z^{n+1}   (Re(z) > 0)

So if {a_λ b} = Σ_{n≥0} c_n λ^n, then:
  r(z) = Σ_{n≥0} c_n · n! / z^{n+1}

Conversely, given an OPE a(z)b(w) ~ Σ c_n / (z-w)^{n+1}:
  n-products: a_n b = c_n
  λ-bracket: {a_λ b} = Σ c_n λ^n / n!

The Laplace bridge connects the z-space (OPE/r-matrix) and
λ-space (λ-bracket) descriptions.

Paper references: Section 18 (spectral-braiding.tex), BR3 axiom.
Cross-volume: Vol I concordance.tex, r-matrix conventions.
"""
from sympy import (Symbol, Rational, factorial, expand, simplify, S,
                   symbols, Poly, oo)


def lambda_bracket_to_r_matrix(bracket_coeffs, z):
    """Convert polynomial λ-bracket to r-matrix via Laplace transform.

    If {a_λ b} = Σ_{n≥0} c_n λ^n, then:
      r(z) = Σ_{n≥0} c_n · n! / z^{n+1}

    Parameters:
        bracket_coeffs: dict {power: coefficient} for the λ-bracket
                       (non-negative powers only; negative powers ignored)
        z: symbol for the z-variable

    Returns:
        Symbolic expression for r(z).
    """
    result = S.Zero
    for n, c_n in bracket_coeffs.items():
        if n < 0:
            continue  # Only regular (non-negative) powers contribute
        result += c_n * factorial(n) / z**(n + 1)
    return expand(result)


def r_matrix_to_lambda_bracket(ope_coeffs):
    """Convert OPE coefficients to λ-bracket.

    Given OPE a(z)b(w) ~ Σ_{n≥0} c_n / (z-w)^{n+1}:
      {a_λ b} = Σ_{n≥0} c_n · λ^n / n!

    Parameters:
        ope_coeffs: dict {pole_order: coefficient}
                   where pole_order n means c_n / (z-w)^{n+1}

    Returns:
        dict {power: coefficient} for the λ-bracket
    """
    result = {}
    for n, c_n in ope_coeffs.items():
        if n < 0:
            continue
        result[n] = c_n / factorial(n)
    return result


def verify_br3_abelian(k, z):
    """Verify BR3 for abelian current algebra.

    {J_λ J} = kλ, so bracket_coeffs = {1: k}.
    The pre-dlog Laplace/OPE kernel is r^L(z) = k · 1! / z^2 = k/z².

    The collision kernel used for Yangian-style spectral braiding is obtained
    after the bar-residue extraction and has one pole fewer.

    Returns:
        (laplace_kernel, expected, difference) — difference should be 0.
    """
    bracket_coeffs = {1: k}  # {J_λ J} = kλ
    r = lambda_bracket_to_r_matrix(bracket_coeffs, z)
    expected = k / z**2
    return r, expected, simplify(r - expected)


def verify_br3_virasoro(c, z):
    """Verify BR3 for Virasoro algebra.

    {T_λ T} = ∂T + 2Tλ + (c/12)λ³
    bracket_coeffs = {0: ∂T, 1: 2T, 3: c/12}

    The pre-dlog Laplace/OPE kernel is
    r^L(z) = ∂T · 0!/z¹ + 2T · 1!/z² + (c/12) · 3!/z⁴
           = ∂T/z + 2T/z² + c/2 · 1/z⁴

    This matches the Virasoro OPE:
      T(z)T(w) ~ (c/2)/(z-w)⁴ + 2T/(z-w)² + ∂T/(z-w)

    Returns:
        (laplace_kernel, expected, difference) — difference should be 0.
    """
    dT = Symbol('dT')
    T = Symbol('T')
    bracket_coeffs = {0: dT, 1: 2*T, 3: Rational(1, 12)*c}
    r = lambda_bracket_to_r_matrix(bracket_coeffs, z)
    expected = dT/z + 2*T/z**2 + Rational(1, 2)*c/z**4
    return r, expected, simplify(expand(r - expected))


def verify_br3_su2(k, z):
    """Verify BR3 for su(2) current algebra.

    {J^a_λ J^b} = ε^{abc}J^c + kδ^{ab}λ

    For a=b (diagonal): {J^a_λ J^a} = kλ
      pre-dlog coefficient r^{aa}(z) = k/z²

    For a≠b (off-diagonal): {J^a_λ J^b} = ε^{abc}J^c
      r^{ab}(z) = ε^{abc}J^c/z

    The full pre-dlog Laplace/OPE kernel:
      r^L(z) = Σ_{a,b} r^{ab}(z) · e_a ⊗ e_b
    where e_a are the basis elements.

    Trace-form collision r-kernel for ĝ_k after the bar-residue extraction:
      r_k(z) = kΩ/z  where Ω = Σ J^a ⊗ J^a (Casimir)

    In unit normalization this gives r^{ab} = δ^{ab}/z, but the pre-dlog diagonal
    coefficient above is k/z²
    and ε^{abc}J^c/z for off-diagonal. The discrepancy is because the
    full r-matrix involves BOTH the structure constant and level terms.

    The pre-dlog Laplace kernel for ĝ_k from the λ-bracket:
      r^{ab}(z) = ε^{abc}J^c · 0!/z + kδ^{ab} · 1!/z²
                = ε^{abc}J^c/z + kδ^{ab}/z²

    This encodes the FULL KM OPE:
      J^a(z)J^b(w) ~ kδ^{ab}/(z-w)² + ε^{abc}J^c/(z-w)

    Returns dict (a,b) -> (laplace_kernel, expected, difference).
    """
    from .examples.nonabelian_cs import levi_civita, killing_form

    J = [None, Symbol('J1'), Symbol('J2'), Symbol('J3')]
    results = {}

    for a in [1, 2, 3]:
        for b in [1, 2, 3]:
            # λ-bracket coefficients for {J^a_λ J^b}
            struct_part = sum(levi_civita(a, b, c) * J[c] for c in [1, 2, 3])
            level_part = k * killing_form(a, b)
            bracket_coeffs = {}
            if struct_part != 0:
                bracket_coeffs[0] = struct_part
            if level_part != 0:
                bracket_coeffs[1] = level_part

            r = lambda_bracket_to_r_matrix(bracket_coeffs, z)

            # Expected from OPE
            expected = S.Zero
            if struct_part != 0:
                expected += struct_part / z
            if level_part != 0:
                expected += level_part / z**2

            diff = simplify(expand(r - expected))
            results[(a, b)] = (r, expected, diff)

    return results
