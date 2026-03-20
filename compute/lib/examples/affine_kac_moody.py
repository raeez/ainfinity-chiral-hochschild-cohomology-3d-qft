"""Affine Kac-Moody PVA: boundary algebra from Chern-Simons theory.

For a simple Lie algebra g with level k, the lambda-bracket is:
  {J^a_lambda J^b} = Sigma_c f^{ab}_c J^c + k*kappa(a,b)*lambda

This is the critical missing example family in Vol II, bridging the
gap between Vol I (which has extensive affine theory) and Vol II
(which only had free/abelian/nonabelian CS at generic level).

References:
  - Vol I: kac_moody.tex (affine foundations), affine_foundations.tex
  - Vol II: pva-descent.tex (D2-D6), spectral-braiding.tex (r-matrix)
  - De Sole-Kac (2006): PVA formalism
"""
from sympy import Symbol, Rational, expand, simplify, S, symbols, Matrix
from dataclasses import dataclass, field
from typing import Dict, Tuple, List, Optional
from fractions import Fraction


# =========================================================================
# Lie algebra data
# =========================================================================

@dataclass
class LieAlgebraData:
    """Structure constants, Killing form, and invariants of a simple Lie algebra.

    Attributes:
        name: identifier (e.g., 'sl2', 'sl3', 'u1')
        dim: dimension of g
        rank: rank of g
        h_dual: dual Coxeter number h^v
        basis_labels: names for basis elements
        structure_constants: dict (a, b, c) -> f^{ab}_c
        killing_form: dict (a, b) -> kappa(a, b)
        is_abelian: True for u(1)
    """
    name: str
    dim: int
    rank: int
    h_dual: int
    basis_labels: List[str]
    structure_constants: Dict[Tuple[int, int, int], int]
    killing_form: Dict[Tuple[int, int], int]
    is_abelian: bool = False

    def f(self, a, b, c):
        """Structure constant f^{ab}_c."""
        return self.structure_constants.get((a, b, c), 0)

    def kappa(self, a, b):
        """Killing form kappa(a, b) (normalized by 1/(2h^v))."""
        return self.killing_form.get((a, b), 0)


def sl2_data():
    """sl_2 Lie algebra data.

    Basis: e (index 1), h (index 2), f (index 3).
    Bracket: [e,f] = h, [h,e] = 2e, [h,f] = -2f.
    Killing form (normalized by 1/(2h^v) = 1/4):
      kappa(e,f) = kappa(f,e) = 1, kappa(h,h) = 2.
    dim = 3, rank = 1, h^v = 2.
    """
    # Structure constants: [J^a, J^b] = sum_c f^{ab}_c J^c
    # Convention: 1=e, 2=h, 3=f
    sc = {}

    # [e, f] = h  =>  f^{1,3,2} = 1
    sc[(1, 3, 2)] = 1
    # [f, e] = -h  =>  f^{3,1,2} = -1
    sc[(3, 1, 2)] = -1

    # [h, e] = 2e  =>  f^{2,1,1} = 2
    sc[(2, 1, 1)] = 2
    # [e, h] = -2e  =>  f^{1,2,1} = -2
    sc[(1, 2, 1)] = -2

    # [h, f] = -2f  =>  f^{2,3,3} = -2
    sc[(2, 3, 3)] = -2
    # [f, h] = 2f  =>  f^{3,2,3} = 2
    sc[(3, 2, 3)] = 2

    # Killing form normalized by 1/(2h^v): kappa(e,f) = kappa(f,e) = 1, kappa(h,h) = 2
    kf = {}
    kf[(1, 3)] = 1   # kappa(e, f)
    kf[(3, 1)] = 1   # kappa(f, e)
    kf[(2, 2)] = 2   # kappa(h, h)

    return LieAlgebraData(
        name='sl2',
        dim=3,
        rank=1,
        h_dual=2,
        basis_labels=['e', 'h', 'f'],
        structure_constants=sc,
        killing_form=kf,
    )


def sl3_data():
    """sl_3 Lie algebra data.

    Chevalley basis: e_1, e_2, e_{12}, f_1, f_2, f_{12}, h_1, h_2.
    dim = 8, rank = 2, h^v = 3.

    Indices: 1=e_1, 2=e_2, 3=e_{12}, 4=f_1, 5=f_2, 6=f_{12}, 7=h_1, 8=h_2.

    Structure constants from the standard Chevalley relations.
    Killing form normalized by 1/(2h^v) = 1/6.
    """
    sc = {}

    # Simple root commutators
    # [e_1, f_1] = h_1
    sc[(1, 4, 7)] = 1; sc[(4, 1, 7)] = -1
    # [e_2, f_2] = h_2
    sc[(2, 5, 8)] = 1; sc[(5, 2, 8)] = -1
    # [e_{12}, f_{12}] = h_1 + h_2
    sc[(3, 6, 7)] = 1; sc[(3, 6, 8)] = 1
    sc[(6, 3, 7)] = -1; sc[(6, 3, 8)] = -1

    # Cartan action on positive roots
    # [h_1, e_1] = 2*e_1
    sc[(7, 1, 1)] = 2; sc[(1, 7, 1)] = -2
    # [h_1, e_2] = -e_2
    sc[(7, 2, 2)] = -1; sc[(2, 7, 2)] = 1
    # [h_1, e_{12}] = e_{12}
    sc[(7, 3, 3)] = 1; sc[(3, 7, 3)] = -1

    # [h_2, e_1] = -e_1
    sc[(8, 1, 1)] = -1; sc[(1, 8, 1)] = 1
    # [h_2, e_2] = 2*e_2
    sc[(8, 2, 2)] = 2; sc[(2, 8, 2)] = -2
    # [h_2, e_{12}] = e_{12}
    sc[(8, 3, 3)] = 1; sc[(3, 8, 3)] = -1

    # Cartan action on negative roots
    # [h_1, f_1] = -2*f_1
    sc[(7, 4, 4)] = -2; sc[(4, 7, 4)] = 2
    # [h_1, f_2] = f_2
    sc[(7, 5, 5)] = 1; sc[(5, 7, 5)] = -1
    # [h_1, f_{12}] = -f_{12}
    sc[(7, 6, 6)] = -1; sc[(6, 7, 6)] = 1

    # [h_2, f_1] = f_1
    sc[(8, 4, 4)] = 1; sc[(4, 8, 4)] = -1
    # [h_2, f_2] = -2*f_2
    sc[(8, 5, 5)] = -2; sc[(5, 8, 5)] = 2
    # [h_2, f_{12}] = -f_{12}
    sc[(8, 6, 6)] = -1; sc[(6, 8, 6)] = 1

    # Root raising/lowering
    # [e_1, e_2] = e_{12}
    sc[(1, 2, 3)] = 1; sc[(2, 1, 3)] = -1
    # [f_1, f_2] = -f_{12}  (sign from Serre relations)
    sc[(4, 5, 6)] = -1; sc[(5, 4, 6)] = 1

    # Mixed positive-negative (non-simple roots)
    # Computed from matrix representation: e_i = E_{row,col} in 3x3 matrices
    # [e_1, f_{12}] = -f_2
    sc[(1, 6, 5)] = -1; sc[(6, 1, 5)] = 1
    # [e_2, f_{12}] = f_1
    sc[(2, 6, 4)] = 1; sc[(6, 2, 4)] = -1
    # [e_{12}, f_1] = -e_2  (from [E_{13}, E_{21}] = -E_{23})
    sc[(3, 4, 2)] = -1; sc[(4, 3, 2)] = 1
    # [e_{12}, f_2] = e_1  (from [E_{13}, E_{32}] = E_{12})
    sc[(3, 5, 1)] = 1; sc[(5, 3, 1)] = -1

    # Killing form (normalized by 1/(2h^v) = 1/6)
    # kappa(e_i, f_i) = 1 for simple roots, kappa(e_{12}, f_{12}) = 1
    # kappa(h_1, h_1) = 2, kappa(h_2, h_2) = 2, kappa(h_1, h_2) = -1
    kf = {}
    kf[(1, 4)] = 1; kf[(4, 1)] = 1   # kappa(e_1, f_1)
    kf[(2, 5)] = 1; kf[(5, 2)] = 1   # kappa(e_2, f_2)
    kf[(3, 6)] = 1; kf[(6, 3)] = 1   # kappa(e_{12}, f_{12})
    kf[(7, 7)] = 2                     # kappa(h_1, h_1)
    kf[(8, 8)] = 2                     # kappa(h_2, h_2)
    kf[(7, 8)] = -1; kf[(8, 7)] = -1  # kappa(h_1, h_2)

    return LieAlgebraData(
        name='sl3',
        dim=8,
        rank=2,
        h_dual=3,
        basis_labels=['e_1', 'e_2', 'e_{12}', 'f_1', 'f_2', 'f_{12}', 'h_1', 'h_2'],
        structure_constants=sc,
        killing_form=kf,
    )


def u1_data():
    """u(1) = abelian Lie algebra, dim = 1.

    The affine u(1)_k IS the Heisenberg algebra.
    h^v = 0 (abelian), but we set h_dual = 0 for consistency.
    Killing form: kappa(1,1) = 1 (single generator).
    """
    return LieAlgebraData(
        name='u1',
        dim=1,
        rank=1,
        h_dual=0,
        basis_labels=['J'],
        structure_constants={},
        killing_form={(1, 1): 1},
        is_abelian=True,
    )


# =========================================================================
# Symbolic generators
# =========================================================================

def _make_generators(g):
    """Create symbolic generator symbols J_1, ..., J_{dim(g)} and their derivatives."""
    n = g.dim
    gens = [None]  # 1-indexed
    dgens = [None]
    for i in range(1, n + 1):
        label = g.basis_labels[i - 1] if i - 1 < len(g.basis_labels) else f'J{i}'
        gens.append(Symbol(f'J_{label}'))
        dgens.append(Symbol(f'dJ_{label}'))
    return gens, dgens


# =========================================================================
# Lambda-bracket
# =========================================================================

def affine_lambda_bracket(g, k, a, b, lam):
    """Compute the lambda-bracket {J^a_lambda J^b} for affine g at level k.

    {J^a_lambda J^b} = Sigma_c f^{ab}_c J^c + k * kappa(a, b) * lambda

    Parameters:
        g: LieAlgebraData instance
        k: level (Symbol or numeric)
        a, b: generator indices (1-indexed)
        lam: spectral parameter (Symbol or numeric)

    Returns:
        Symbolic expression: linear combination of generator symbols + constant * lambda.
    """
    gens, _ = _make_generators(g)

    # Structure constant part: sum_c f^{ab}_c J^c
    struct_part = S.Zero
    for c in range(1, g.dim + 1):
        fc = g.f(a, b, c)
        if fc != 0:
            struct_part += fc * gens[c]

    # Level part: k * kappa(a, b) * lambda
    kab = g.kappa(a, b)
    level_part = k * kab * lam

    return expand(struct_part + level_part)


def affine_lambda_bracket_numeric(g, k_val, a, b, lam_val):
    """Evaluate the lambda-bracket at numeric k and lambda.

    Returns a dict {generator_index: coefficient, 'const': scalar}.
    Useful for testing without symbolic overhead.
    """
    result = {}
    for c in range(1, g.dim + 1):
        fc = g.f(a, b, c)
        if fc != 0:
            result[c] = fc
    kab = g.kappa(a, b)
    if kab != 0:
        result['lambda_coeff'] = k_val * kab
    return result


# =========================================================================
# PVA axiom verification
# =========================================================================

def verify_pva_jacobi_affine(g, k):
    """Check Jacobi identity for all triples (a, b, c) of generators.

    {J^a_lam {J^b_mu J^c}} - {J^b_mu {J^a_lam J^c}} = {{J^a_lam J^b}_{lam+mu} J^c}

    For a current algebra, the Jacobi identity reduces to the Lie algebra
    Jacobi identity for structure constants plus invariance of the Killing form.

    Returns:
        dict mapping (a, b, c) -> symbolic result (all should be 0).
    """
    lam, mu = symbols('lambda mu')
    gens, _ = _make_generators(g)
    nu = lam + mu
    n = g.dim

    results = {}
    for a in range(1, n + 1):
        for b in range(1, n + 1):
            for c in range(1, n + 1):
                # Term 1: {J^a_lam {J^b_mu J^c}}
                # Inner: {J^b_mu J^c} = sum_d f^{bc}_d J^d + k*kappa(b,c)*mu
                # Outer on J^d part: {J^a_lam J^d} = sum_e f^{ad}_e J^e + k*kappa(a,d)*lam
                # Outer on constant part: {J^a_lam const} = 0
                term1 = S.Zero
                for d in range(1, n + 1):
                    fbc_d = g.f(b, c, d)
                    if fbc_d == 0:
                        continue
                    bracket_a_d = affine_lambda_bracket(g, k, a, d, lam)
                    term1 += fbc_d * bracket_a_d

                # Term 2: -{J^b_mu {J^a_lam J^c}}
                term2 = S.Zero
                for d in range(1, n + 1):
                    fac_d = g.f(a, c, d)
                    if fac_d == 0:
                        continue
                    bracket_b_d = affine_lambda_bracket(g, k, b, d, mu)
                    term2 += fac_d * bracket_b_d

                lhs = expand(term1 - term2)

                # RHS: {{J^a_lam J^b}_{lam+mu} J^c}
                # {J^a_lam J^b} = sum_d f^{ab}_d J^d + k*kappa(a,b)*lam
                # Bracket the J^d part: {J^d_nu J^c}
                # The constant k*kappa(a,b)*lam part: {const_nu J^c} = 0
                rhs = S.Zero
                for d in range(1, n + 1):
                    fab_d = g.f(a, b, d)
                    if fab_d == 0:
                        continue
                    bracket_d_c = affine_lambda_bracket(g, k, d, c, nu)
                    rhs += fab_d * bracket_d_c

                rhs = expand(rhs)
                results[(a, b, c)] = simplify(expand(lhs - rhs))

    return results


def verify_skew_symmetry_affine(g, k):
    """Check skew-symmetry {J^a_lam J^b} = -{J^b_{-lam-partial} J^a} for all pairs.

    For generators of a current algebra, skew-symmetry follows from:
    - Antisymmetry of structure constants: f^{ab}_c = -f^{ba}_c
    - Symmetry of Killing form: kappa(a,b) = kappa(b,a)

    The Borcherds formula gives:
    -{J^b_{-lam-partial} J^a} = -(sum_c f^{ba}_c J^c + k*kappa(b,a)*(-lam-partial))
    partial acts on everything to the right, but J^c and k*kappa are field-independent,
    so partial(J^c) = dJ^c and partial(const) = 0.

    At the n=0 level (no derivatives): -{sum_c f^{ba}_c J^c - k*kappa(b,a)*lam}
    = sum_c f^{ab}_c J^c + k*kappa(a,b)*lam  (using antisymmetry + symmetry of kappa)

    Returns:
        dict mapping (a, b) -> symbolic result (all should be 0).
    """
    lam = Symbol('lambda')
    n = g.dim
    results = {}

    for a in range(1, n + 1):
        for b in range(1, n + 1):
            lhs = affine_lambda_bracket(g, k, a, b, lam)

            # RHS via Borcherds: at generator level, the relevant n-products are
            # (J^b)_0 (J^a) = sum_c f^{ba}_c J^c  and  (J^b)_1 (J^a) = k*kappa(b,a)
            # Borcherds skew: {J^a_lam J^b} = -sum_n (-lam-partial)^n/n! (J^b_n J^a)
            # n=0: -(J^b_0 J^a) = -sum_c f^{ba}_c J^c = sum_c f^{ab}_c J^c  [antisymmetry]
            # n=1: -(-lam)(J^b_1 J^a) = lam * k*kappa(b,a) = k*kappa(a,b)*lam  [symmetry]
            # (partial on constant = 0, higher n-products vanish)
            rhs_n0 = S.Zero
            gens, _ = _make_generators(g)
            for c in range(1, n + 1):
                fba_c = g.f(b, a, c)
                if fba_c != 0:
                    rhs_n0 += fba_c * gens[c]
            rhs_n0 = -rhs_n0  # the minus from Borcherds

            rhs_n1 = k * g.kappa(b, a) * lam  # -(-lam)*kappa = lam*kappa

            rhs = expand(rhs_n0 + rhs_n1)
            results[(a, b)] = simplify(expand(lhs - rhs))

    return results


# =========================================================================
# Central charge, kappa, and Feigin-Frenkel duality
# =========================================================================

def affine_central_charge(g, k):
    """Sugawara central charge: c = k * dim(g) / (k + h^v).

    For sl_2: c = 3k/(k+2).
    For sl_N: c = k(N^2-1)/(k+N).

    CRITICAL: UNDEFINED at k = -h^v (critical level).
    At the critical level, the Sugawara construction fails because
    the denominator vanishes. This is NOT 'c diverges' but rather
    the Sugawara tensor is not defined.

    Parameters:
        g: LieAlgebraData
        k: level (sympy expression or numeric)

    Returns:
        c as sympy expression. Raises ValueError at critical level.
    """
    if g.is_abelian:
        # For u(1): c = 1 per free boson (no Sugawara denominator issue)
        return S.One

    k_sym = S(k)
    h = S(g.h_dual)
    if simplify(k_sym + h) == 0:
        raise ValueError(
            f"Sugawara central charge is UNDEFINED at critical level k = -h^v = {-g.h_dual}. "
            f"The Sugawara tensor does not exist at the critical level (Feigin-Frenkel)."
        )
    return k_sym * g.dim / (k_sym + h)


def sugawara_central_charge(g, k):
    """Alias for affine_central_charge. Sugawara c = k*dim(g)/(k + h^v).

    UNDEFINED at critical level k = -h^v.
    """
    return affine_central_charge(g, k)


def affine_kappa(g, k):
    """Modular Koszul curvature: kappa = c/2.

    For affine g at level k: kappa = k*dim(g) / (2*(k + h^v)).

    Vol I convention: kappa(A) = c(A)/2 for all chiral algebras.
    """
    if g.is_abelian:
        return Rational(1, 2)
    return affine_central_charge(g, k) / 2


def ff_dual_level(g, k):
    """Feigin-Frenkel dual level: k' = -k - 2*h^v.

    This is the involution of the affine Weyl group that exchanges
    V_k(g) and V_{k'}(g) at the level of Koszul duality.

    For sl_2: k' = -k - 4.
    For sl_3: k' = -k - 6.
    For sl_N: k' = -k - 2N.
    """
    if g.is_abelian:
        # u(1): self-dual, k' = -k
        return -S(k)
    return -S(k) - 2 * g.h_dual


def kappa_complementarity_affine(g, k):
    """Verify kappa-complementarity: kappa(V_k) + kappa(V_{k'}) = dim(g).

    For sl_N at any non-critical k:
      c(k) + c(k') = 2*dim(g)
    so kappa(k) + kappa(k') = dim(g).

    Returns:
        The sum kappa(k) + kappa(k') - dim(g), which should be 0.
    """
    if g.is_abelian:
        # u(1): kappa(k) = 1/2, kappa(-k) = 1/2, sum = 1 = dim
        k_prime = ff_dual_level(g, k)
        return affine_kappa(g, k) + affine_kappa(g, k_prime) - g.dim

    k_sym = S(k)
    k_prime = ff_dual_level(g, k)
    kap = affine_kappa(g, k_sym)
    kap_prime = affine_kappa(g, k_prime)
    return simplify(kap + kap_prime - g.dim)


# =========================================================================
# Classical r-matrix
# =========================================================================

def _inverse_killing_form(g):
    """Compute the inverse Killing form kappa^{ab}.

    The Casimir element uses kappa^{ab} (the inverse matrix of kappa_{ab}),
    NOT kappa_{ab} itself. For the Chevalley basis, the Killing form matrix
    is not diagonal, so kappa^{ab} != kappa_{ab} in general.

    Returns a dict (a, b) -> kappa^{ab} for nonzero entries.
    """
    n = g.dim
    # Build the Killing form matrix
    K = Matrix.zeros(n, n)
    for a in range(1, n + 1):
        for b in range(1, n + 1):
            K[a - 1, b - 1] = g.kappa(a, b)

    # Invert
    K_inv = K.inv()

    result = {}
    for a in range(n):
        for b in range(n):
            val = K_inv[a, b]
            if val != 0:
                result[(a + 1, b + 1)] = val
    return result


def classical_r_matrix(g):
    """Classical r-matrix: r(z) = Omega / z where Omega is the Casimir.

    Omega = sum_{a,b} kappa^{ab} J_a tensor J_b

    where kappa^{ab} is the INVERSE Killing form (not kappa_{ab}).

    For the current algebra, the Laplace transform of the lambda-bracket
    gives r(z) = Omega/z. This satisfies the classical Yang-Baxter equation.

    Returns:
        (casimir_pairs, pole_order) where casimir_pairs is a list of
        (a, b, kappa^{ab}) entries and pole_order = 1.
    """
    kappa_inv = _inverse_killing_form(g)
    casimir_pairs = [(a, b, c) for (a, b), c in kappa_inv.items()]
    return casimir_pairs, 1  # Simple pole at z=0


def _casimir_tensor(g):
    """Compute the Casimir tensor Omega = sum kappa^{ab} e_a tensor e_b.

    Returns a dict (a, b) -> coefficient for nonzero entries.
    Uses the INVERSE Killing form kappa^{ab}, which is required for
    the Casimir to be ad-invariant.
    """
    return _inverse_killing_form(g)


def verify_cybe(g):
    """Verify the classical Yang-Baxter equation for r(u) = Omega/u.

    The spectral-parameter CYBE:
      [r_12(u), r_13(u+v)] + [r_12(u), r_23(v)] + [r_13(u+v), r_23(v)] = 0

    For r(u) = Omega/u, substituting and clearing denominators gives:
      v[Omega_12, Omega_13] + (u+v)[Omega_12, Omega_23] + u[Omega_13, Omega_23] = 0

    Matching coefficients of u and v independently, both conditions reduce to
    the ad-invariance of the Casimir:

      [Omega_12, Omega_13 + Omega_23] = 0

    This is equivalent to: for all x in g,
      sum_{a,b} kappa^{ab} ([x, e_a] tensor e_b + e_a tensor [x, e_b]) = 0

    or in components: for each (m, j, l),
      sum_i kappa^{ij} f^{mi}_a kappa^{al} + sum_i kappa^{im} f^{ij}_a kappa^{al} = 0
      [first term: commutator in slot 2; second: in slot 3]

    Actually we verify the simplest form: ad-invariance of the Casimir,
    which in components means:
      sum_b kappa^{ab} f^{xb}_c + sum_a kappa^{ac} f^{xa}_b  -- NO

    The cleanest test: [Omega_12, Omega_13 + Omega_23] = 0 as a tensor in g^{tensor 3}.

    Returns:
        dict (a, b, c) -> coefficient. All should be 0.
    """
    n = g.dim
    omega = _casimir_tensor(g)

    # [Omega_12, Omega_13 + Omega_23] in g^{tensor 3}
    #
    # Omega_12 = sum_{i,j} omega(i,j) e_i (x) e_j (x) 1
    # Omega_13 = sum_{k,l} omega(k,l) e_k (x) 1 (x) e_l
    # Omega_23 = sum_{k,l} omega(k,l) 1 (x) e_k (x) e_l
    #
    # [Omega_12, Omega_13]: commutator acts in slot 1
    # = sum_{i,j,k,l} omega(i,j)*omega(k,l) [e_i, e_k] (x) e_j (x) e_l
    # = sum_{i,j,k,l,m} omega(i,j)*omega(k,l)*f^{ik}_m  e_m (x) e_j (x) e_l
    #
    # [Omega_12, Omega_23]: commutator acts in slot 2
    # = sum_{i,j,k,l} omega(i,j)*omega(k,l) e_i (x) [e_j, e_k] (x) e_l
    # = sum_{i,j,k,l,m} omega(i,j)*omega(k,l)*f^{jk}_m  e_i (x) e_m (x) e_l

    result = {}

    # Term 1: [Omega_12, Omega_13]^{m,j,l}
    for (i, j), oij in omega.items():
        for (k, l), okl in omega.items():
            for m in range(1, n + 1):
                fikm = g.f(i, k, m)
                if fikm != 0:
                    key = (m, j, l)
                    result[key] = result.get(key, S.Zero) + oij * okl * fikm

    # Term 2: [Omega_12, Omega_23]^{i,m,l}
    for (i_idx, j), oij in omega.items():
        for (k, l), okl in omega.items():
            for m in range(1, n + 1):
                fjkm = g.f(j, k, m)
                if fjkm != 0:
                    key = (i_idx, m, l)
                    result[key] = result.get(key, S.Zero) + oij * okl * fjkm

    # Simplify and filter
    simplified = {}
    for key, val in result.items():
        v = simplify(val)
        if v != 0:
            simplified[key] = v
    return simplified


# =========================================================================
# Shadow data
# =========================================================================

def shadow_data_affine(g, k):
    """Shadow archetype and depth data for affine g at level k.

    For any simple g:
    - Shadow archetype: L (Lie/tree, terminates at depth 3)
    - Shadow depth: 3 (cubic shadow from Lie bracket, no quartic)
    - kappa: c/2
    - Cubic shadow C_3: nonzero (from structure constants f^{ab}_c)
    - Quartic contact Q_4: zero (no self-referential OPE generates quartic)

    For u(1):
    - Shadow archetype: G (Gaussian, terminates at depth 2)
    - Shadow depth: 2

    Vol I: affine_foundations.tex, concordance.tex sec:concordance-shadow-archetypes.
    """
    if g.is_abelian:
        return {
            'archetype': 'G',     # Gaussian
            'depth': 2,
            'kappa': affine_kappa(g, k),
            'cubic_shadow': S.Zero,
            'quartic_contact': S.Zero,
            'description': 'Gaussian (Heisenberg): terminates at depth 2',
        }

    return {
        'archetype': 'L',     # Lie/tree
        'depth': 3,
        'kappa': affine_kappa(g, k),
        'cubic_shadow': 'nonzero',  # Proportional to structure constants
        'quartic_contact': S.Zero,
        'description': 'Lie/tree: terminates at depth 3 (cubic from bracket)',
    }


# =========================================================================
# Lie algebra Jacobi identity (underlying finite-dimensional algebra)
# =========================================================================

def verify_lie_jacobi(g):
    """Verify the Jacobi identity [a,[b,c]] + [b,[c,a]] + [c,[a,b]] = 0.

    This is the FOUNDATION for the PVA Jacobi identity: the affine
    PVA Jacobi reduces to the finite-dimensional Lie Jacobi plus
    invariance of the Killing form.

    Returns:
        dict (a, b, c) -> result for all triples. All should be 0.
    """
    n = g.dim
    results = {}
    for a in range(1, n + 1):
        for b in range(1, n + 1):
            for c in range(1, n + 1):
                # [a, [b, c]] = sum_d f^{bc}_d * f^{ad}_e for each e
                val = {}
                for d in range(1, n + 1):
                    fbc_d = g.f(b, c, d)
                    if fbc_d == 0:
                        continue
                    for e in range(1, n + 1):
                        fad_e = g.f(a, d, e)
                        if fad_e != 0:
                            val[e] = val.get(e, 0) + fbc_d * fad_e

                # + [b, [c, a]]
                for d in range(1, n + 1):
                    fca_d = g.f(c, a, d)
                    if fca_d == 0:
                        continue
                    for e in range(1, n + 1):
                        fbd_e = g.f(b, d, e)
                        if fbd_e != 0:
                            val[e] = val.get(e, 0) + fca_d * fbd_e

                # + [c, [a, b]]
                for d in range(1, n + 1):
                    fab_d = g.f(a, b, d)
                    if fab_d == 0:
                        continue
                    for e in range(1, n + 1):
                        fcd_e = g.f(c, d, e)
                        if fcd_e != 0:
                            val[e] = val.get(e, 0) + fab_d * fcd_e

                # All coefficients should be 0
                nonzero = {e: v for e, v in val.items() if v != 0}
                results[(a, b, c)] = nonzero  # empty dict = identity holds

    return results


def verify_killing_invariance(g):
    """Verify ad-invariance of the Killing form:

    kappa([a,b], c) + kappa(b, [a,c]) = 0 for all a, b, c.

    Equivalently: sum_d f^{ab}_d kappa(d, c) + sum_d f^{ac}_d kappa(b, d) = 0.

    This is the second ingredient (along with Lie Jacobi) that guarantees
    the PVA Jacobi identity for a current algebra.

    Returns:
        dict (a, b, c) -> result for all triples. All should be 0.
    """
    n = g.dim
    results = {}
    for a in range(1, n + 1):
        for b in range(1, n + 1):
            for c in range(1, n + 1):
                # sum_d f^{ab}_d kappa(d, c)
                term1 = 0
                for d in range(1, n + 1):
                    fab_d = g.f(a, b, d)
                    if fab_d != 0:
                        kdc = g.kappa(d, c)
                        term1 += fab_d * kdc

                # sum_d f^{ac}_d kappa(b, d)
                term2 = 0
                for d in range(1, n + 1):
                    fac_d = g.f(a, c, d)
                    if fac_d != 0:
                        kbd = g.kappa(b, d)
                        term2 += fac_d * kbd

                results[(a, b, c)] = term1 + term2

    return results


# =========================================================================
# Cross-volume bridge helpers
# =========================================================================

def kappa_from_vol1_formula(g, k):
    """Compute kappa using the Vol I formula: kappa = c/2.

    This is the authoritative formula from Vol I concordance.tex.
    Should agree exactly with affine_kappa().
    """
    return affine_kappa(g, k)


def shadow_depth_from_vol1(g):
    """Shadow depth from Vol I classification.

    Vol I concordance.tex sec:concordance-shadow-archetypes:
    - G (Gaussian): depth 2 (Heisenberg, u(1))
    - L (Lie/tree): depth 3 (affine, simple g)
    - C (contact): depth 4 (beta-gamma)
    - M (mixed): depth infinity (Virasoro, W_N)
    """
    if g.is_abelian:
        return 2
    return 3
