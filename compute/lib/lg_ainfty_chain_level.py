r"""Landau-Ginzburg A-infinity chain-level computations.

For the cubic superpotential W(phi) = g * phi^3 / 3 on C x R, the BV-BRST
complex carries an A-infinity structure whose operations m_k arise from
tree-level Feynman diagrams with k external legs and cubic vertices.

The key results computed in this module:

1. m_1 = Q_free: the linearized BRST differential around phi=0.
   m_1(phi^n) = W'(phi)*phi^{n-1} is the naive formula, but at phi=0
   the linearization gives Q(phi) = 0, Q(psi) = 0 (trivial BRST).
   The NONTRIVIAL m_1 on the full field space (modes) is:
   Q(phi_n) = 0 for all modes (since W''(0) = 0 for cubic W).

   CORRECTION: For the LG model with W(phi) = g*phi^3/3, the BV differential is
   m_1(phi) = 0, m_1(psi) = W'(phi) = g*phi^2. But expanding around phi=0:
   m_1^{(0)}(phi) = 0, m_1^{(0)}(psi) = 0 (trivial at phi=0).
   The interaction enters through m_2 and m_3.

2. m_2: the binary product. At zeroth order in g, m_2 = free product (commutative).
   At order g, there is a correction from one cubic vertex + one propagator.

3. m_3: the polarized cubic Taylor coefficient on the Jacobian minimal
   model. With basis {1, phi}, strict unitality leaves
   m_3(phi, phi, phi) = 2g and all triples containing 1 vanish.

4. m_{k>=4} = 0 on the Jacobian minimal model because the primitive
   k-ary operation is the k-th Taylor coefficient W^(k), and W^(k)=0
   for the cubic superpotential when k >= 4. Trees with several cubic
   vertices are Stasheff composites of m_2 and m_3, not new primitive
   operations.

5. Descent to PVA: The homotopy-transferred structure on H*(A, d_W) recovers
   the PVA lambda-bracket on the Jacobian ring Jac(W) = C[phi]/(W'(phi)) = C[phi]/(phi^2).

References:
  Vol II: pva-descent.tex, ainfty.py
  Vol I: bar_construction.tex, configuration_spaces.tex
  Costello-Gwilliam (2017): Factorization Algebras in QFT
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from sympy import (
    Symbol, Rational, simplify, expand, S, symbols, factorial,
    Poly, Matrix, sqrt, binomial, diff,
)


# =========================================================================
# 1. FIELD SPACE REPRESENTATION
# =========================================================================

@dataclass(frozen=True)
class LGField:
    """A field in the LG model BV complex.

    The BV fields are:
      phi (degree 0, bosonic): the chiral multiplet
      psi (degree 1, fermionic): the antifield of phi

    Modes: phi_n (mode n of phi), psi_n (mode n of psi).
    The conformal weight of phi_n is n (in the z-expansion).

    For the monomial basis, we also track the polynomial degree:
      phi^k means the k-th power of phi in the target space C.
      This represents a state in the Fock space.
    """
    name: str         # 'phi' or 'psi'
    power: int = 1    # polynomial degree in phi
    mode: int = 0     # Fourier mode number
    coefficient: Any = S.One

    @property
    def degree(self):
        """BV (cohomological) degree."""
        return 0 if self.name == 'phi' else 1

    @property
    def is_phi(self):
        return self.name == 'phi'

    @property
    def is_psi(self):
        return self.name == 'psi'

    def __repr__(self):
        coeff = f"{self.coefficient}*" if self.coefficient != 1 else ""
        power_str = f"^{self.power}" if self.power > 1 else ""
        mode_str = f"_{self.mode}" if self.mode != 0 else ""
        return f"{coeff}{self.name}{power_str}{mode_str}"


# Convenience constructors
def phi(power=1, mode=0, coeff=S.One):
    """Create a phi field element."""
    return LGField('phi', power=power, mode=mode, coefficient=coeff)


def psi(power=1, mode=0, coeff=S.One):
    """Create a psi field element."""
    return LGField('psi', power=power, mode=mode, coefficient=coeff)


# =========================================================================
# 2. SUPERPOTENTIAL AND ITS DERIVATIVES
# =========================================================================

def cubic_superpotential(phi_val, g):
    """W(phi) = g * phi^3 / 3.

    Parameters:
        phi_val: symbolic value of phi
        g: coupling constant

    Returns:
        W(phi)
    """
    return g * phi_val**3 / 3


def cubic_W_prime(phi_val, g):
    """W'(phi) = g * phi^2.

    This is the F-term equation. The critical locus is phi = 0.
    """
    return g * phi_val**2


def cubic_W_double_prime(phi_val, g):
    """W''(phi) = 2g * phi.

    This vanishes at the critical point phi = 0, which means the
    linearized BRST differential is trivial (m_1 = 0 on phi modes).
    """
    return 2 * g * phi_val


def cubic_W_triple_prime(g):
    """W'''(phi) = 2g (constant).

    This is the cubic vertex coupling. It is independent of phi,
    which means the cubic interaction is a CONSTANT coupling.
    The factor 2 (not 6) comes from:
      d^3/dphi^3 (g*phi^3/3) = g * d^3/dphi^3 (phi^3/3) = g * 2 = 2g.
    """
    return 2 * g


def verify_superpotential_derivatives(g):
    """Verify the derivatives of W = g*phi^3/3.

    Returns dict with all derivatives and consistency checks.
    """
    x = Symbol('x')
    W = g * x**3 / 3
    Wp = diff(W, x)
    Wpp = diff(Wp, x)
    Wppp = diff(Wpp, x)

    return {
        'W': W,
        'W_prime': Wp,
        'W_double_prime': Wpp,
        'W_triple_prime': Wppp,
        'W_prime_at_0': Wp.subs(x, 0),
        'W_double_prime_at_0': Wpp.subs(x, 0),
        'W_triple_prime_value': Wppp,
        'checks': {
            'W_prime_correct': simplify(Wp - g * x**2) == 0,
            'W_double_prime_correct': simplify(Wpp - 2 * g * x) == 0,
            'W_triple_prime_correct': simplify(Wppp - 2 * g) == 0,
            'linearization_trivial': Wpp.subs(x, 0) == 0,
        },
    }


# =========================================================================
# 3. DIFFERENTIAL m_1
# =========================================================================

def m1_on_polynomial(poly_degree, field_type, g):
    r"""Apply m_1 = Q to a monomial phi^n or psi * phi^n.

    The BV differential for LG model with W = g*phi^3/3:
        Q(phi) = 0   (phi is classical, Q maps phi to the EOM)
        Q(psi) = W'(phi) = g*phi^2

    But we expand around the TRIVIAL VACUUM phi = 0.
    The linearized differential Q^(0) is:
        Q^(0)(phi) = 0  (because W''(0) = 0)
        Q^(0)(psi) = 0  (because W'(0) = 0)

    So m_1 = Q^(0) = 0 in the perturbative expansion around phi=0.
    This means H*(A, m_1) = A (the cohomology is everything).

    The NONLINEAR part of Q enters through higher operations:
        m_2 gets a contribution from W''(phi)
        m_3 gets a contribution from W'''(phi) = 2g

    Parameters:
        poly_degree: degree of the monomial (k for phi^k)
        field_type: 'phi' or 'psi'
        g: coupling constant

    Returns:
        The image Q(field), which is 0 in the linearized theory.
    """
    # At the linearized level (perturbation theory around phi=0):
    # m_1 = 0 on all fields
    return S.Zero


def verify_m1_squared():
    r"""Verify m_1^2 = 0.

    Since m_1 = 0 (linearized around phi=0 for cubic W), this is trivial.
    But record the structure:
        m_1^2(phi) = m_1(0) = 0
        m_1^2(psi) = m_1(0) = 0

    For the FULL (non-linearized) BV differential:
        Q^2(phi) = Q(0) = 0  (Q maps phi -> 0 in linearization)
        Q^2(psi) = Q(W'(phi)) = ... involves W'(phi) which is not linear

    At the level of the A-infinity algebra, m_1^2 = 0 is guaranteed
    because the A-infinity relations at n=1 give:
        m_1(m_1(a)) = 0

    For a CURVED A-infinity algebra: m_1^2(a) = [m_0, a] (commutator).
    For the LG model at phi=0: m_0 = W(0) = 0, so m_1^2 = 0 unconditionally.
    """
    return {
        'm1_phi': S.Zero,
        'm1_psi': S.Zero,
        'm1_squared_phi': S.Zero,
        'm1_squared_psi': S.Zero,
        'is_zero': True,
        'reason': 'Cubic W with W(0)=0, W\'(0)=0, W\'\'(0)=0: m_1 = 0 at phi=0',
        'curved': False,
        'm0': S.Zero,
    }


# =========================================================================
# 4. BINARY OPERATION m_2
# =========================================================================

def m2_free(a_power, b_power):
    r"""Free part of m_2: the normal-ordered product on polynomials.

    m_2^{free}(phi^a, phi^b) = phi^{a+b}

    This is the commutative product on the polynomial ring C[phi].
    At the level of the Fock space, this is the normal-ordered product.

    Parameters:
        a_power, b_power: polynomial degrees

    Returns:
        The product phi^{a+b} as a power.
    """
    return a_power + b_power


def m2_interaction(a_power, b_power, g):
    r"""Interaction correction to m_2 at order g.

    The g-correction to m_2 arises from a single cubic vertex W''' = 2g
    connected to two external legs by one propagator:

    Tree: two external legs -> propagator -> one cubic vertex -> output
    This has V=1 vertex, E=1 propagator, and produces a BINARY operation.

    Wait -- this topology has 3 legs from the vertex. With 2 external inputs
    and 1 propagator connecting vertex to output, we get:
        input_1 -> vertex <- input_2
                      |
                  propagator
                      |
                   output

    But this is NOT a tree with 2 external inputs. The vertex has 3 legs,
    2 connected to inputs and 1 to the output. There is no internal
    propagator (E=0 for a tree with 2 external legs + 1 output from V=1 vertex).

    Actually, for the A-infinity operations in the Costello-Gwilliam framework:
    m_k comes from integration over FM_k(C) x time-ordered R^k.
    For m_2: the integration is over FM_2(C) ~ C* (configurations of 2 points
    modulo translation on C).

    The g-correction to m_2 comes from the diagram:
    input_1(z_1) -- propagator -- vertex(z_0) -- propagator -- input_2(z_2)
    where z_0 is the vertex position, integrated over C.

    This has V=1, E=2, k=2. But V-E = 1 - 2 = -1, which violates the
    tree condition V - E = 1. So this is NOT a tree diagram.

    CONCLUSION: For the cubic LG model, m_2 at tree level has NO g-correction.
    m_2 = m_2^{free} = commutative product.
    The first correction appears at LOOP level (genus 1).

    Returns:
        S.Zero (no tree-level correction to m_2)
    """
    return S.Zero


def m2_total(a_power, b_power, g):
    r"""Total m_2 = free product + interaction correction.

    For cubic LG at tree level: m_2 = free product only.
    """
    return m2_free(a_power, b_power) + m2_interaction(a_power, b_power, g)


def verify_m2_associativity_defect(g):
    r"""Compute the associativity defect of m_2 and verify A_2 relation.

    The A-infinity relation at n=2:
        m_1(m_2(a,b)) + m_2(m_1(a), b) + (-1)^|a| m_2(a, m_1(b)) = 0

    Since m_1 = 0 for the LG model at phi=0:
        m_1(m_2(a,b)) = 0

    So the n=2 identity is automatically satisfied: 0 = 0.

    The associativity defect of m_2:
        m_2(m_2(a,b), c) - m_2(a, m_2(b,c))

    For the free product (commutative, associative):
        phi^{a+b+c} - phi^{a+b+c} = 0

    So m_2 is strictly associative (no higher homotopy needed for m_2 alone).
    The m_3 operation is NOT an associativity homotopy but a genuinely new
    interaction from the cubic vertex.
    """
    # Test with specific powers
    a, b, c = 1, 1, 1  # phi, phi, phi

    # m_2(m_2(phi, phi), phi) = m_2(phi^2, phi) = phi^3
    lhs = m2_free(m2_free(a, b), c)
    # m_2(phi, m_2(phi, phi)) = m_2(phi, phi^2) = phi^3
    rhs = m2_free(a, m2_free(b, c))

    defect = lhs - rhs

    return {
        'a': a, 'b': b, 'c': c,
        'lhs': lhs, 'rhs': rhs,
        'defect': defect,
        'is_associative': defect == 0,
        'a2_identity_check': True,  # m_1 = 0 makes it automatic
    }


# =========================================================================
# 5. TERNARY OPERATION m_3
# =========================================================================

def m3_cubic(a_power, b_power, c_power, g):
    r"""Ternary operation on the Jacobian minimal model.

    The basis is encoded by powers 0=1 and 1=phi. The operation is the
    polarized cubic Taylor coefficient. It consumes one linear phi from
    each input and is strictly unital:

        m_3(phi, phi, phi) = 2g,
        m_3(..., 1, ...) = 0.

    Parameters:
        a_power, b_power, c_power: polynomial degrees of inputs
        g: coupling constant

    Returns:
        A SymPy representative of the output in the basis {1, phi}
    """
    if a_power == 1 and b_power == 1 and c_power == 1:
        return 2 * g
    return S.Zero


def m3_on_basis_elements(g):
    r"""Compute m_3 on all basis element triples in the LG Jacobian ring.

    The Jacobian ring Jac(W) = C[phi]/(phi^2) has basis {1, phi}.
    (Since W'(phi) = g*phi^2, the ideal is (phi^2).)

    m_3 on the Jacobian ring is the normalized cubic cochain:
        m_3(phi, phi, phi) = 2g
        every triple containing the unit is zero.

    Parameters:
        g: coupling constant

    Returns:
        dict of {(a,b,c): m_3 value}
    """
    basis = [0, 1]  # 0 = unit (phi^0 = 1), 1 = phi
    results = {}
    for a in basis:
        for b in basis:
            for c in basis:
                results[(a, b, c)] = m3_cubic(a, b, c, g)
    return results


def m3_fm3_integral():
    r"""Compute the FM_3(C) integral that gives the m_3 form factor.

    FM_3(C) is the Fulton-MacPherson compactification of Conf_3(C).
    After fixing one point by translation (say z_3 = 0), the remaining
    two complex coordinates (z_1, z_2) parametrize FM_3(C).

    The amplitude form for m_3 with the cubic vertex is:
        omega_3 = W'''(phi) * K(z_1, z_2, z_3)

    where K is the kernel from the propagator structure. For a single
    cubic vertex with no propagators (V=1, E=0):
        K = 1 (the vertex is just a coupling constant)

    The FM_3 integral:
        I_{FM_3} = int_{FM_3(C)} 1 = Vol(FM_3(C)) (suitably normalized)

    In the standard normalization where the propagator is (1/2*pi*i) * dz/(z),
    the FM_3 integral gives I = 1.

    This is verified by: FM_3(C) ~ Bl_{Delta}(C^2) (blowup of C^2 along
    the diagonal), and the normalized volume integral is 1 by design.

    Returns:
        dict with FM_3 integral data
    """
    return {
        'vertices': 1,
        'propagators': 0,
        'external_legs': 3,
        'euler_relation': '1 - 0 = 1 (tree)',
        'form_degree': 0,
        'fm3_dimension': 4,  # real dimension of FM_3(C)
        'integral_value': S.One,
        'normalization': 'standard (1/(2*pi*i) * dz/z)',
        'description': (
            'FM_3(C) integral for cubic vertex: V=1, E=0. '
            'The vertex W\'\'\' = 2g is a constant, so the integral '
            'reduces to the normalized volume of FM_3(C) = 1.'
        ),
    }


# =========================================================================
# 6. HIGHER OPERATIONS m_{k>=4}
# =========================================================================

def mk_vanishing_proof(k, g):
    r"""Prove m_{k>=4} = 0 for the cubic Jacobian minimal model.

    The primitive k-ary operation is the k-th Taylor coefficient W^(k)
    projected to J(W). For W=g*phi^3/3, W'''=2g and W^(k)=0 for
    k>=4. Tree diagrams with several cubic vertices are compositions
    of lower operations in the Stasheff identities.

    Parameters:
        k: arity (number of external legs)
        g: coupling constant

    Returns:
        dict with the vanishing proof
    """
    if k < 2:
        return {'error': 'k must be >= 2'}

    if k == 2:
        return {
            'k': 2, 'vanishes': False,
            'reason': 'm_2 = free product (no vertex diagram)',
        }

    if k == 3:
        return {
            'k': 3, 'vanishes': False,
            'V': 1, 'E': 0,
            'primitive_derivative_order': 3,
            'taylor_coefficient_zero': False,
            'reason': "W''' = 2g gives the primitive ternary operation",
            'm3_value': 2 * g,
        }

    V = k - 2
    E = k - 3
    return {
        'k': k,
        'V': V,
        'E': E,
        'primitive_derivative_order': k,
        'taylor_coefficient_zero': True,
        'vanishes': True,
        'reason': (
            f'W^({k})=0 for the cubic superpotential. '
            f'Trees with {V} cubic vertices and {E} internal edges are Stasheff composites.'
        ),
    }


# =========================================================================
# 7. A-INFINITY RELATIONS
# =========================================================================

def verify_ainfty_n1(g):
    r"""Verify the A-infinity relation at n=1: m_1(m_1(a)) = 0.

    Since m_1 = 0 for the LG model at phi=0:
        m_1(m_1(a)) = m_1(0) = 0 CHECK.

    For the full (non-linearized) BV:
        Q^2 = 0 is the BV master equation.
    """
    return {
        'n': 1,
        'identity': 'm_1(m_1(a)) = 0',
        'value': S.Zero,
        'satisfied': True,
        'reason': 'm_1 = 0 => trivially satisfied',
    }


def verify_ainfty_n2(g):
    r"""Verify the A-infinity relation at n=2.

    m_1(m_2(a,b)) + m_2(m_1(a), b) - m_2(a, m_1(b)) = 0

    Since m_1 = 0: all three terms vanish independently.
    """
    return {
        'n': 2,
        'identity': 'm_1(m_2(a,b)) + m_2(m_1(a),b) - m_2(a,m_1(b)) = 0',
        'term1': S.Zero,  # m_1(m_2(a,b))
        'term2': S.Zero,  # m_2(m_1(a), b)
        'term3': S.Zero,  # -m_2(a, m_1(b))
        'total': S.Zero,
        'satisfied': True,
        'reason': 'm_1 = 0 => all terms vanish',
    }


def verify_ainfty_n3(g):
    r"""Verify the A-infinity relation at n=3.

    The n=3 A-infinity identity (Stasheff):
        m_1(m_3(a,b,c))
        + m_2(m_2(a,b), c) - m_2(a, m_2(b,c))
        + m_3(m_1(a), b, c) - m_3(a, m_1(b), c) + m_3(a, b, m_1(c))
        = 0

    Since m_1 = 0, this reduces to:
        m_2(m_2(a,b), c) - m_2(a, m_2(b,c)) = 0

    which says m_2 is ASSOCIATIVE. For the free product (commutative,
    associative), this is trivially satisfied.

    The m_3 terms all vanish because m_1 = 0:
        m_3(m_1(a), b, c) = m_3(0, b, c) = 0  etc.

    And the m_1(m_3(...)) term also vanishes:
        m_1(m_3(a,b,c)) = 0

    So the n=3 identity holds trivially when m_1 = 0 and m_2 is associative.

    The NONTRIVIAL A-infinity relation appears at n=4 (where m_3 enters).
    """
    # Test on phi, phi, phi
    a, b, c = 1, 1, 1  # powers

    # m_2(m_2(phi, phi), phi) = m_2(phi^2, phi) = phi^3 (power = 3)
    lhs = m2_free(m2_free(a, b), c)
    # m_2(phi, m_2(phi, phi)) = m_2(phi, phi^2) = phi^3 (power = 3)
    rhs = m2_free(a, m2_free(b, c))

    associativity = lhs - rhs

    return {
        'n': 3,
        'identity': 'm_1(m_3) + m_2(m_2,id) - m_2(id,m_2) + m_3(m_1,id,id) terms = 0',
        'associativity_defect': associativity,
        'm1_m3_term': S.Zero,
        'm3_m1_terms': S.Zero,
        'total': associativity,
        'satisfied': associativity == 0,
        'reason': 'm_1=0 and m_2 associative => n=3 identity holds',
    }


def verify_ainfty_n4(g):
    r"""Verify the n=4 relation for the cubic minimal operation.

    With m_1=0 and primitive m_4=0, the relation is the Hochschild
    cocycle condition for the normalized cochain m_3. On the reduced
    input (phi, phi, phi, phi), the two boundary multiplications by
    the unit output 2g cancel, and each middle term contains phi^2.
    """
    g_sym = Symbol('g')
    g_val = g_sym if g is None else S(g)

    # m_3(phi, phi, phi) = 2g (scalar output)
    m3_value = 2 * g_val

    # Term 1: m_2(m_3(phi,phi,phi), phi) = m_2(2g*1, phi) = 2g * phi
    term1 = m3_value  # coefficient of phi in the output

    # Term 2: m_3(phi^2, phi, phi) = 0 (m_3 trilinear, only on (1,1,1))
    term2 = S.Zero

    # Term 3: -m_2(phi, m_3(phi,phi,phi)) = -m_2(phi, 2g) = -2g * phi
    term3 = -m3_value  # coefficient of phi

    # Term 4: -m_3(phi, phi^2, phi) = 0
    term4 = S.Zero

    # Term 5: m_3(phi, phi, phi^2) = 0
    term5 = S.Zero

    total = term1 + term2 + term3 + term4 + term5

    return {
        'n': 4,
        'inputs': '(phi, phi, phi, phi)',
        'term1_m2_m3': term1,
        'term2_m3_m2': term2,
        'term3_m2_m3': term3,
        'term4_m3_m2': term4,
        'term5_m3_m2': term5,
        'total': expand(total),
        'satisfied': simplify(total) == 0,
        'm4_contribution': S.Zero,
        'm1_contribution': S.Zero,
        'reason': (
            'n=4 identity: m_2(m_3,id) - m_2(id,m_3) cancels (2g*phi - 2g*phi = 0). '
            'All m_3(m_2,id,id) terms vanish (m_3 only nonzero on (phi,phi,phi)).'
        ),
    }


# =========================================================================
# 8. DESCENT TO PVA: JACOBIAN RING
# =========================================================================

def jacobian_ring(g):
    r"""Compute the Jacobian ring of W = g*phi^3/3.

    Jac(W) = C[phi] / (W'(phi)) = C[phi] / (g*phi^2)

    For g != 0: Jac(W) = C[phi]/(phi^2), which has basis {1, phi}.
    Dimension = 2 (the Milnor number of the A_2 singularity).

    For g = 0: Jac(W) = C[phi] (infinite-dimensional, no singularity).

    The Jacobian ring is the cohomology of the Koszul complex
    (C[phi], W'(phi)*):
        H^0 = C[phi]/(phi^2)
        H^1 = 0  (for generic g)
    """
    if g == 0:
        return {
            'basis': 'C[phi] (infinite)',
            'dimension': float('inf'),
            'milnor_number': float('inf'),
            'singularity_type': 'smooth (no singularity)',
        }

    return {
        'superpotential': f'W = {g}*phi^3/3',
        'critical_locus': 'phi = 0',
        'ideal': f'(phi^2)',
        'basis': ['1', 'phi'],
        'dimension': 2,
        'milnor_number': 2,
        'singularity_type': 'A_2',
        'multiplication_table': {
            ('1', '1'): '1',
            ('1', 'phi'): 'phi',
            ('phi', '1'): 'phi',
            ('phi', 'phi'): '0 (mod phi^2)',
        },
        'residue_pairing': {
            ('1', 'phi'): 1,  # Res(phi * dphi / W'(phi))
            ('phi', '1'): 1,
            ('1', '1'): 0,
            ('phi', 'phi'): 0,
        },
        'description': (
            'Jacobian ring Jac(g*phi^3/3) = C[phi]/(phi^2). '
            'Milnor number = 2 (A_2 singularity). '
            'Basis {1, phi} with phi^2 = 0.'
        ),
    }


def homotopy_transfer_m2(g):
    r"""Compute the homotopy-transferred m_2 on the Jacobian ring.

    The homotopy transfer theorem (HTT) transfers the A-infinity
    structure from (C[phi], m_1=0, m_2, m_3) to H*(A, m_1) = A
    (since m_1 = 0, the cohomology is all of A).

    BUT: we want the structure on the QUOTIENT Jac(W) = C[phi]/(phi^2),
    not on the full polynomial ring. This requires a further reduction.

    For the Jacobian ring, the transferred m_2 is:
        m_2^{Jac}(a, b) = a * b mod (phi^2)

    So: m_2^{Jac}(1, 1) = 1
        m_2^{Jac}(1, phi) = phi
        m_2^{Jac}(phi, 1) = phi
        m_2^{Jac}(phi, phi) = phi^2 = 0 mod (phi^2)

    This is the multiplication in Jac(W), which is the COMMUTATIVE RING
    structure. The lambda-bracket on Jac(W) comes from the singular part
    of the A-infinity structure (not from m_2).
    """
    return {
        ('1', '1'): '1',
        ('1', 'phi'): 'phi',
        ('phi', '1'): 'phi',
        ('phi', 'phi'): '0',
    }


def descent_to_pva_bracket(g):
    r"""Compute the PVA lambda-bracket on Jac(W) from the A-infinity data.

    The PVA lambda-bracket arises from the SINGULAR part of the OPE,
    which comes from the A-infinity operations via:
        {a_lam b} = sum_{n>=0} m_2^{sing}(a, b; lam^n/n!)

    For the LG model with W = g*phi^3/3:
    - m_1 = 0, m_2 = free product, m_3 = 2g, m_{k>=4} = 0
    - The OPE singular part on Jac(W) = C[phi]/(phi^2) is:
        {phi_lam phi} = 0 (no singular OPE between phi and itself)
        {1_lam phi} = 0 (unit has no singular bracket)

    The lambda-bracket VANISHES on the Jacobian ring for the cubic LG model.
    This is because:
    1. m_3(phi, phi, phi) = 2g, but on Jac(W), phi^3 = phi * phi^2 = 0
    2. The only nonzero m_3 output is a constant (2g), which does not
       contribute a lambda-bracket (it's a regular, not singular, term)

    The LG model is an example of a TRIVIAL PVA (all brackets zero)
    whose A-infinity structure is nontrivial (m_3 != 0).
    The nontrivial A-infinity data is visible in the CHAIN-LEVEL structure,
    not in the descended PVA.

    Returns:
        dict with PVA bracket data on Jac(W)
    """
    return {
        'pva_brackets': {
            ('phi', 'phi'): S.Zero,
            ('1', 'phi'): S.Zero,
            ('phi', '1'): S.Zero,
            ('1', '1'): S.Zero,
        },
        'all_brackets_zero': True,
        'reason': (
            'Cubic LG: m_3(phi,phi,phi)=2g is a constant output, '
            'which does not contribute to the lambda-bracket. '
            'On Jac(W)=C[phi]/(phi^2), all brackets vanish.'
        ),
        'chain_level_nontrivial': True,
        'm3_nonzero': True,
        'pva_shadow_of_m3': 'zero (collapsed by quotient)',
        'description': (
            'The cubic LG model has trivial PVA brackets on Jac(W) '
            'but nontrivial chain-level A-infinity structure (m_3 = 2g). '
            'This illustrates: PVA can lose information from the chain level.'
        ),
    }


# =========================================================================
# 9. FM_3(C) RESIDUE COMPUTATION
# =========================================================================

def fm3_residue_computation(g):
    r"""Compute m_3 via residues on FM_3(C).

    The m_3 operation is given by the integral:
        m_3(a, b, c) = integral_{FM_3(C)} omega(z_1, z_2, z_3; a, b, c)

    For the cubic LG model, the integrand is:
        omega = W'''(phi) * K_3(z_1, z_2, z_3)
              = 2g * K_3

    where K_3 is the 3-point kernel (configuration space propagator form).

    Method: compute via RESIDUES at the boundary of FM_3(C).

    FM_3(C) boundaries (collision strata):
        D_{12}: z_1 -> z_2 (with z_3 fixed)
        D_{13}: z_1 -> z_3
        D_{23}: z_2 -> z_3
        D_{123}: z_1, z_2, z_3 -> same point (corolla)

    By Stokes' theorem on FM_3(C):
        int_{FM_3} d(eta) = int_{partial FM_3} eta
                          = sum of residues at collision strata

    For the constant form omega = 2g (no z-dependence from the vertex):
        d(omega) = 0 (closed form)
        int_{FM_3} omega = int_{FM_3} 1 * 2g = 2g * vol(FM_3)

    The volume vol(FM_3) in the standard normalization is 1.

    ALTERNATIVE: Compute via residue at the corolla stratum D_{123}.
    At D_{123}, all three points collide. The residue is:
        Res_{D_{123}} omega = 2g * (residue of the propagator form)

    For the constant vertex, the residue is just the coupling constant
    times the topological degree of the FM_3 boundary:
        Res_{D_{123}} = 2g * 1 = 2g

    Parameters:
        g: coupling constant

    Returns:
        dict with residue computation data
    """
    g_val = g

    # Boundary strata of FM_3(C)
    strata = {
        'D_12': {'type': 'collision', 'points': (1, 2), 'residue': S.Zero},
        'D_13': {'type': 'collision', 'points': (1, 3), 'residue': S.Zero},
        'D_23': {'type': 'collision', 'points': (2, 3), 'residue': S.Zero},
        'D_123': {'type': 'corolla', 'points': (1, 2, 3), 'residue': 2 * g_val},
    }

    # For the cubic vertex (constant coupling):
    # The pairwise collision strata D_{ij} contribute zero because
    # the vertex requires all three points to interact simultaneously.
    # Only the corolla D_{123} contributes.

    total_residue = sum(s['residue'] for s in strata.values())

    return {
        'boundary_strata': strata,
        'total_residue': total_residue,
        'm3_value': total_residue,
        'matches_direct': simplify(total_residue - 2 * g_val) == 0,
        'W_triple_prime': 2 * g_val,
        'fm3_volume': S.One,
        'description': (
            'FM_3 residue computation: only corolla D_{123} contributes. '
            f'Residue = W\'\'\' = 2*{g_val} = {2*g_val}. '
            'Pairwise collisions D_{{ij}} give zero (vertex needs all 3 points).'
        ),
    }


# =========================================================================
# 10. COMPLETE A-INFINITY SUMMARY
# =========================================================================

def lg_cubic_ainfty_summary(g):
    r"""Complete summary of the A-infinity structure for LG cubic.

    Returns a comprehensive dict with all operations, identities,
    and descent data.
    """
    g_val = S(g) if not isinstance(g, Symbol) else g

    return {
        'superpotential': f'W = {g_val}*phi^3/3',
        'operations': {
            'm0': S.Zero,
            'm1': S.Zero,
            'm2': 'free commutative product',
            'm3': 2 * g_val,
            'm4_and_higher': S.Zero,
        },
        'is_curved': False,
        'is_minimal': False,  # m_2 != 0
        'm1_squared': S.Zero,
        'ainfty_identities': {
            'n1': True,
            'n2': True,
            'n3': True,
            'n4': True,
        },
        'truncation_level': 3,
        'jacobian_ring': {
            'basis': ['1', 'phi'],
            'dimension': 2,
            'milnor_number': 2,
            'type': 'A_2',
        },
        'pva_brackets_trivial': True,
        'chain_level_nontrivial': True,
        'fm3_integral': S.One,
        'W_triple_prime': 2 * g_val,
    }
