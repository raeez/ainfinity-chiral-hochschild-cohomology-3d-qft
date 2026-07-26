"""
Landau-Ginzburg model with cubic superpotential.

Setup: Chiral multiplet phi on C x R with superpotential W(phi) = g*phi^3/3.

BV fields: phi (degree 0), psi (antifield, degree 1)

CRITICAL DISTINCTION (from CLAUDE.md: curved A-infinity):
The A-infinity structure from BV-BRST around the TRIVIAL vacuum phi=0 has:
  m_1 = Q_free: Q(phi) = psi, Q(psi) = 0, Q^2 = 0
  m_0 = 0 (uncurved: no vacuum energy for polynomial W)

The superpotential W enters through HIGHER operations, not m_1:
  m_2: free product + O(g) correction from one vertex + one propagator
  m_3: nonzero, proportional to g (from cubic vertex W''' = 2g)
  m_{k>=4} = 0 as primitive operations because W^(k)=0 for k>=4

Physical reasoning:
- m_1 = Q is the linearized BRST differential around phi=0
- W'(0) = 0, W''(0) = 0 for cubic W, so the linearization is free
- The cubic vertex W'''=2g is the interaction vertex
- higher cubic trees are Stasheff composites of m_2 and m_3

Cohomology: H*(A, Q) = C (same as free theory — each (phi_n, psi_n)
is an acyclic doublet under Q).

The Jacobian ring C[phi]/(phi^2) emerges via A-infinity homotopy transfer,
NOT from H*(A,Q) directly. The homotopy-transferred m_2 on the subspace
of phi^0 and phi^1 modes encodes the Jacobian ring structure.

Paper reference: Section 11 (examples-computing.tex), second example.
"""
from sympy import Symbol, Rational, S, symbols, expand


# ---------------------------------------------------------------------------
# BV Field content
# ---------------------------------------------------------------------------
# Working with mode expansion: phi_n (degree 0), psi_n (degree 1).
# The coupling constant g parametrizes the superpotential.

def m1_lg(field_value, field_type):
    """Differential Q = Q_free for the LG model expanded around phi=0.

    Q(phi_n) = psi_n,   Q(psi_n) = 0.

    This satisfies Q^2 = 0:
      Q^2(phi_n) = Q(psi_n) = 0.
      Q^2(psi_n) = Q(0) = 0.

    The interaction W = g*phi^3/3 does NOT modify m_1 when expanding
    around the trivial vacuum phi=0:
      W'(phi) = g*phi^2, so W'(0) = 0
      W''(phi) = 2g*phi, so W''(0) = 0
    The linearized BRST differential is therefore identical to free theory.

    Parameters:
        field_value: symbolic value of the field
        field_type: 'phi' or 'psi'

    Returns:
        Q(field): psi_value if phi, 0 if psi
    """
    if field_type == 'phi':
        return field_value  # Q(phi_n) = psi_n (value carries over)
    elif field_type == 'psi':
        return S.Zero  # Q(psi_n) = 0
    else:
        raise ValueError(f"Unknown field type: {field_type}")


def check_Q_squared_lg(field_value, field_type):
    """Verify Q^2 = 0 for the LG model.

    Q^2(phi) = Q(Q(phi)) = Q(psi) = 0 ✓
    Q^2(psi) = Q(0) = 0 ✓
    """
    Qf = m1_lg(field_value, field_type)
    if field_type == 'phi':
        Q2f = m1_lg(Qf, 'psi')  # Q(psi) = 0
    else:
        Q2f = S.Zero
    return Q2f


def m2_lg_free_part(a, b):
    """Free part of m_2: the normal-ordered product.

    At g=0, m_2 is the free product (same as free_multiplet.m2).
    This is the regular part — no spectral parameter dependence.
    """
    return a * b


def m3_lg_vertex(a, b, c, g):
    """Ternary operation from the cubic vertex.

    m_3(a, b, c) is proportional to g * W''' = 2g (the cubic coupling).

    For the tree diagram with 3 inputs and 1 cubic vertex:
    - V = 1 vertex (the W''' = 2g vertex)
    - E = 0 internal edges (propagators)
    - The vertex directly connects all 3 inputs

    The m_3 operation is:
      m_3(a, b, c) = 2g · (FM_3 integral of the tree amplitude)

    For degree-0 inputs (all phi's):
      m_3(phi, phi, phi) ~ 2g · (form factor from FM_3(C))

    For the SYMBOLIC computation, we record:
      m_3 = 2g · a · b · c · (geometric factor from FM_3)

    The geometric factor is an integral over FM_3(C) which is a point
    (after collapsing the three points to the vertex). So:
      m_3(a, b, c) = 2g · a · b · c  (up to the FM_3 normalization)

    CAUTION: The FM_3 normalization depends on conventions.
    The factor 2 comes from W'''(phi) = d^3(phi^3/3)/dphi^3 = 2.
    """
    return 2 * g * a * b * c


def m_k_lg(k, args, g):
    """General m_k for the LG cubic model.

    k=1: Q_free (differential)
    k=2: commutative product (free part; g-correction not yet implemented)
    k=3: 2g * a * b * c (cubic vertex)
    k>=4: 0 as a primitive operation because W^(k)=0 for cubic W

    Parameters:
        k: arity
        args: tuple of k field values
        g: coupling constant
    """
    if k <= 0:
        raise ValueError(f"Arity must be positive, got {k}")
    if k == 1:
        raise ValueError("Use m1_lg() directly (needs field_type)")
    if k == 2:
        return m2_lg_free_part(args[0], args[1])
    if k == 3:
        return m3_lg_vertex(args[0], args[1], args[2], g)
    return S.Zero  # m_{k>=4} = 0


def check_truncation_degree_counting(k):
    """Legacy helper for the cubic truncation statement.

    The primitive k-ary operation is the k-th Taylor coefficient W^(k).
    For W=g*phi^3/3, W'''=2g and W^(k)=0 for k>=4. Trees with several
    cubic vertices are Stasheff composites, not new primitive operations.

    Returns:
        dict with tree topology data and Taylor-coefficient analysis
    """
    if k < 2:
        return {'arity': k, 'error': 'arity must be >= 2'}

    if k == 2:
        return {
            'arity': 2,
            'vertices': 0,
            'internal_edges': 0,
            'ghost_degree': 0,
            'expected_mk_degree': -1,
            'conclusion': 'm_2 comes from free propagator, not vertex diagram',
            'vanishes': False,
        }

    V = k - 2   # vertices for tree diagrams
    E = k - 3   # internal edges (propagators)

    ghost_degree = V - E  # = 1 for all k >= 3

    # Form degree analysis
    vanishes = (k >= 4)

    return {
        'arity': k,
        'vertices': V,
        'internal_edges': E,
        'ghost_degree': ghost_degree,
        'expected_mk_degree': 1 - k,
        'primitive_derivative_order': k,
        'taylor_coefficient_zero': vanishes,
        'vanishes': vanishes,
        'reason': "W^(k)=0 for cubic W" if vanishes else "W''' gives the cubic operation" if k == 3 else 'free product',
    }


def verify_ainfty_n3_lg(a, b, c, g):
    """Verify the n=3 A∞ identity for the LG model.

    For degree-0 elements (all phi) with m_1 = Q_free:

    Σ terms = m_1(m_3(a,b,c)) ± m_3(m_1(a),b,c) ± m_3(a,m_1(b),c) ± m_3(a,b,m_1(c))
            + m_2(m_2(a,b),c) ± m_2(a,m_2(b,c))

    For Q-closed elements (a,b,c ∈ ker Q, i.e., representing cohomology classes):
      m_1(anything) = Q(anything) is zero on the LHS since m_1(m_3) is Q-exact
      m_3(Q(a),...) = 0 since Q(a) = 0

    So the identity reduces to:
      Q(m_3(a,b,c)) + m_2(m_2(a,b),c) ± m_2(a,m_2(b,c)) = 0

    For the free m_2 (at g=0 order):
      m_2(m_2(a,b),c) ± m_2(a,m_2(b,c)) = 2abc (anti-associativity for deg 0)

    And Q(m_3(a,b,c)) = Q(2g·abc) which depends on the degree of abc.

    This is a consistency check, not a full proof of the A∞ identity
    (which requires the spectral parameter structure).
    """
    # m_2 at zeroth order = free product
    term_m2m2 = m2_lg_free_part(m2_lg_free_part(a, b), c)
    term_m2m2_2 = m2_lg_free_part(a, m2_lg_free_part(b, c))

    # For degree-0 elements with Q=0, the n=3 identity (at order g) gives:
    # Q(m_3(a,b,c)) + (m_2 anti-assoc terms at order g^0) = 0
    # This is a statement about how the g-correction to m_2 compensates m_3.

    # At zeroth order in g: anti-associativity of free m_2
    free_n3 = expand(term_m2m2 + term_m2m2_2)  # Should be 2abc for deg 0

    # At first order in g: m_3 term
    m3_term = m3_lg_vertex(a, b, c, g)

    return {
        'free_antiassoc': free_n3,
        'm3_contribution': m3_term,
        'note': 'Full A∞ identity requires spectral parameters and g-corrected m_2',
    }
