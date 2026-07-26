"""3D Gravity Compute Engine: Virasoro line comparison and genus expansion.

Implements the gravitational sector of the holomorphic-topological programme.

The Virasoro algebra Vir_c governs the boundary-CFT reading of 3d
gravity via the holomorphic-topological correspondence: the boundary
chiral algebra is Vir_c, the exact algebraic bulk is the completed
derived chiral centre on C x R, and the same-family Virasoro comparison
representative Vir_{26-c} models the line-side dual only after the
ordered-bar/Verdier comparison datum is installed
(c_dual = 26 - c, comparison-fixed at c = 13).
This engine does not assert a dynamical-metric path integral.

Key mathematical objects:

1. **Virasoro lambda-bracket**: {T_lam T} = (c/12) lam^3 + 2T lam + dT.
   This encodes the OPE T(z)T(w) ~ (c/2)/(z-w)^4 + 2T(w)/(z-w)^2 + dT(w)/(z-w).

2. **Ternary associator**: A_3(T,T,T; lam12, lam23) from the failure of
   associativity of the lambda-bracket. m_3 = -A_3 is the ternary A-infinity
   operation from homotopy transfer.

3. **Quartic contact invariant**: Q^contact_Vir = 10/(c(5c+22)), the first
   nonlinear modular shadow coefficient beyond kappa.

4. **Gravitational kappa**: kappa(Vir_c) = c/2 (Theorem D, modular
   characteristic for 3d gravity).

5. **Genus generating function**: F_g = kappa_eff * lambda_g^FP via the
   A-hat genus formula, where kappa_eff = (c-26)/2 accounts for the
   Koszul dual shift.

6. **R-matrix pole structure**: The pre-dlog Laplace/OPE kernel has an
   order-4 leading pole with residue c/2. The bar collision r-matrix
   absorbs one pole and has leading order 3:
   r^{coll}(z) = (c/2)/z^3 + 2T/z.

References:
  Vol II: 3d_gravity.tex (Movement VI)
  Vol I: concordance.tex (Theorem D), higher_genus_modular_koszul.tex
  Vol I: nonlinear_modular_shadows.tex (quartic contact)
  Vol I: genus_generating_function.py (A-hat coefficients)
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any, Dict, List, Optional, Tuple

from sympy import (
    Symbol, Rational, simplify, expand, S, symbols, factorial,
    bernoulli, sin, sinh, series, sqrt, Abs, collect, Poly, log,
    fraction, binomial, limit, oo, exp, I, pi,
)


# =========================================================================
# 1. VIRASORO LAMBDA-BRACKET
# =========================================================================

def virasoro_lambda_bracket(c=None):
    """Return the Virasoro lambda-bracket {T_lam T}.

    {T_lam T} = (c/12) lam^3 + 2T lam + dT

    This encodes the OPE:
      T(z)T(w) ~ (c/2)/(z-w)^4 + 2T(w)/(z-w)^2 + dT(w)/(z-w)

    via the state-field dictionary: T_{(n)} T corresponds to the
    (n+1)-th pole in the OPE, and {T_lam T} = sum_n T_{(n)}T lam^n/n!.

    Pole dictionary:
      (z-w)^{-4}: c/2    -> T_{(3)}T = c/2   -> coeff of lam^3/3! = c/12
      (z-w)^{-2}: 2T     -> T_{(1)}T = 2T    -> coeff of lam^1/1! = 2T
      (z-w)^{-1}: dT     -> T_{(0)}T = dT    -> coeff of lam^0/0! = dT

    Returns a dict with keys 'lam3', 'lam1_T', 'lam0_dT' for the
    coefficients of lam^3, T*lam, and dT respectively.

    Parameters
    ----------
    c : optional
        Central charge. If None, uses symbolic Symbol('c').
    """
    c_val = Symbol('c') if c is None else S(c)

    return {
        'lam3': c_val / 12,        # coefficient of lam^3
        'lam1_T': S(2),            # coefficient of T * lam
        'lam0_dT': S(1),           # coefficient of dT (= partial T)
        'central_charge': c_val,
    }


# =========================================================================
# 2. VIRASORO ASSOCIATOR A_3
# =========================================================================

def virasoro_associator(c=None, lam12=None, lam23=None):
    r"""Compute the ternary associator A_3(T, T, T; lam12, lam23).

    The associator measures the failure of the lambda-bracket to be
    associative (it satisfies the Jacobi identity instead). For the
    Virasoro algebra with a single generator T of conformal weight 2:

    A_3(T, T, T; lam12, lam23) =
        -d^2 T
        - (2 lam12 + 3 lam23) dT
        - 2 lam23 (2 lam12 + lam23) T
        - (c/12) lam23^3 (2 lam12 + lam23)

    This is computed from the Jacobi identity residue:
        {T_{lam12} {T_{lam23} T}} - {T_{lam23} {T_{lam12} T}}
        - {{T_{lam12} T}_{lam12+lam23} T}

    Returns a dict of coefficients for {d^2T, dT, T, scalar}
    as polynomials in (lam12, lam23, c).
    """
    c_val = Symbol('c') if c is None else S(c)
    l12 = Symbol('lambda_12') if lam12 is None else S(lam12)
    l23 = Symbol('lambda_23') if lam23 is None else S(lam23)

    # Coefficient of d^2 T
    coeff_d2T = S(-1)

    # Coefficient of dT
    coeff_dT = -(2 * l12 + 3 * l23)

    # Coefficient of T
    coeff_T = -2 * l23 * (2 * l12 + l23)

    # Scalar term (no field, pure central extension contribution)
    coeff_scalar = -c_val / 12 * l23**3 * (2 * l12 + l23)

    return {
        'd2T': coeff_d2T,
        'dT': expand(coeff_dT),
        'T': expand(coeff_T),
        'scalar': expand(coeff_scalar),
        'c': c_val,
        'lam12': l12,
        'lam23': l23,
    }


# =========================================================================
# 3. m_3 = -A_3
# =========================================================================

def virasoro_m3_coefficients(c=None, lam12=None, lam23=None):
    r"""Compute m_3(T, T, T; lam12, lam23) = -A_3(T, T, T; lam12, lam23).

    The ternary A-infinity operation m_3 is the negation of the associator.
    This follows from the general homotopy transfer formula: the A-infinity
    operations are defined such that the failure of associativity at order n
    is controlled by m_{n+1}.

    Returns a dict with the same keys as virasoro_associator, but negated.
    """
    A3 = virasoro_associator(c=c, lam12=lam12, lam23=lam23)

    return {
        'd2T': -A3['d2T'],
        'dT': expand(-A3['dT']),
        'T': expand(-A3['T']),
        'scalar': expand(-A3['scalar']),
        'c': A3['c'],
        'lam12': A3['lam12'],
        'lam23': A3['lam23'],
    }


# =========================================================================
# 4. HPL MINIMAL TRANSFER DATA
# =========================================================================

def hpl_planar_binary_tree_count(n: int):
    r"""Number of planar rooted binary trees with n labelled leaves.

    The homological perturbation transfer of a binary operation mu=m_2
    sums over PRT_n.  Its cardinality is the Catalan number C_{n-1}.
    """
    if n < 2:
        raise ValueError("PRT_n is used for n >= 2 in the binary HPL transfer")
    return factorial(2 * n - 2) / (factorial(n) * factorial(n - 1))


def virasoro_hpl_tree_profile(n: int) -> Dict[str, Any]:
    r"""Combinatorial profile of an n-ary HPL summand family for Virasoro.

    Each tree is evaluated in the suspended-bar convention: leaves are
    labelled by \tilde i, internal homotopy edges by \tilde h, internal
    vertices by suspended products, and the root by \tilde p.  The
    left-to-right depth-first order on internal edges fixes the
    orientation line, so there is no independent plus/minus choice.
    """
    return {
        'arity': n,
        'tree_set': f'PRT_{n}',
        'summands': hpl_planar_binary_tree_count(n),
        'leaves': n,
        'binary_vertices_mu': n - 1,
        'internal_h_edges': n - 2,
        'root_label': 'p',
        'leaf_label': 'i',
        'sign_convention': 'suspended-bar Koszul sign',
        'internal_edge_order': 'left-to-right depth-first planar order',
        'orientation_line': 'det(k^{E_int(T)})',
        'free_pm_choice': False,
        'formula_label': 'eq:gravity-suspended-hpl-transfer',
    }


def virasoro_hpl_sdr_data(weight=None) -> Dict[str, Any]:
    r"""SDR identities used by the Virasoro minimal transfer.

    On a positive conformal-weight complement, an antighost G_0 satisfying
    [Q,G_0]=L_0 gives the normalized homotopy h|A_w = G_0/w.
    """
    data: Dict[str, Any] = {
        'homotopy_identity': 'Qh + hQ = id - i∘p',
        'side_conditions': ('pi = 1_H', 'ph = 0', 'hi = 0', 'h^2 = 0'),
        'finite_hpl_datum': (
            'Q_DS',
            'mu=m_2',
            'delta_DS',
            'i',
            'p',
            'h_DS',
            'complete filtration',
        ),
        'positive_weight_homotopy': 'h|A_w = G_0/w',
        'requires_antighost': 'G_0 with [Q,G_0] = L_0',
    }
    if weight is not None:
        weight_val = S(weight)
        if weight_val <= 0:
            raise ValueError("the G_0/w homotopy is only defined for positive weight")
        data['weight'] = weight_val
        data['homotopy_coefficient'] = 1 / weight_val
    return data


def virasoro_ds_linear_sdr_identities(k=None) -> Dict[str, Any]:
    r"""Verify the linear DS SDR identities on generators.

    The associated-graded linear DS complex has generators W,V,b,U,c with

      d0(b)=V, d0(U)=2c, d0(W)=d0(V)=d0(c)=0,
      h(V)=b, h(c)=U/2, h(W)=h(b)=h(U)=0,
      ip(W)=W and ip=0 on V,b,U,c.

    This verifies id - ip = d0 h + h d0 on the C[partial]-generators.
    It is not a contraction of d_BRST = d0 + delta before HPL perturbation.
    """
    k_val = Symbol('k') if k is None else S(k)
    if simplify(k_val + 2) == 0:
        raise ValueError("the DS linear SDR uses W=(k+2)T and requires k != -2")

    basis = ('W', 'V', 'b', 'U', 'c')
    zero = {x: S.Zero for x in basis}

    def vec(**coeffs):
        out = dict(zero)
        for key, value in coeffs.items():
            out[key] = S(value)
        return out

    d0 = {
        'W': vec(),
        'V': vec(),
        'b': vec(V=1),
        'U': vec(c=2),
        'c': vec(),
    }
    h = {
        'W': vec(),
        'V': vec(b=1),
        'b': vec(),
        'U': vec(),
        'c': vec(U=Rational(1, 2)),
    }
    ip = {
        'W': vec(W=1),
        'V': vec(),
        'b': vec(),
        'U': vec(),
        'c': vec(),
    }

    def apply(linear_map, input_vec):
        out = vec()
        for source, coeff in input_vec.items():
            if coeff == 0:
                continue
            image = linear_map[source]
            for target, image_coeff in image.items():
                out[target] += coeff * image_coeff
        return {key: simplify(value) for key, value in out.items()}

    def add(v1, v2):
        return {key: simplify(v1[key] + v2[key]) for key in basis}

    def sub(v1, v2):
        return {key: simplify(v1[key] - v2[key]) for key in basis}

    identity = {x: vec(**{x: 1}) for x in basis}
    lhs = {x: sub(identity[x], ip[x]) for x in basis}
    rhs = {x: add(apply(d0, h[x]), apply(h, d0[x])) for x in basis}
    homotopy_identity = all(lhs[x] == rhs[x] for x in basis)

    h_squared = {
        x: apply(h, h[x])
        for x in basis
    }
    h2_zero = all(h_squared[x] == zero for x in basis)
    hi_zero = h['W'] == zero
    ph_zero = all(target != 'W' or coeff == 0 for x in basis for target, coeff in h[x].items())

    return {
        'complex': 'C_DS_lin over C[partial]',
        'level_condition': 'k != -2',
        'k': k_val,
        'basis': basis,
        'differential': 'd0(b)=V, d0(U)=2c, d0(W)=d0(V)=d0(c)=0',
        'homotopy': 'h(V)=b, h(c)=U/2, h(W)=h(b)=h(U)=0',
        'homotopy_degree': -1,
        'pi_identity': True,
        'h2_zero': h2_zero,
        'ph_zero': ph_zero,
        'hi_zero': hi_zero,
        'homotopy_identity': homotopy_identity,
        'full_brst_retract': False,
        'perturbation_needed_for_full_brst': 'delta = d_BRST - d0',
    }


def virasoro_s4_completed_ambient_requirement(c=None) -> Dict[str, Any]:
    r"""Record the class-M ambient forced by the Virasoro quartic shadow."""
    c_val = Symbol('c') if c is None else S(c)
    s4 = quartic_contact_virasoro(c_val)
    return {
        'S4': s4,
        'bar_weight': 4,
        'ambient': 'weight-completed/pro (hypAmbientWtCpl)',
        'raw_direct_sum_chain_statement': False,
        'completed_statement': True,
        'singular_values': (S(0), Rational(-22, 5)),
    }


def virasoro_exact_gravity_scope_profile(c=None) -> Dict[str, Any]:
    r"""Exact algebraic scope of the Virasoro 3d-gravity rung.

    The rigorous statement is the completed derived-centre claim

      Z_der^ch(Vir_c)^hat_rho in E_3^top-Alg,  0 < rho < |c|/6.

    Brown-Henneaux supplies the physical dictionary c=3 ell/(2G_N);
    it does not turn this algebraic statement into a dynamical-metric
    path integral without additional modular and saddle hypotheses.
    """
    c_val = Symbol('c') if c is None else S(c)
    return {
        'boundary_algebra': 'Vir_c',
        'exact_algebraic_statement': (
            'Z_der^ch(Vir_c)^hat_rho in E_3^top-Alg'
        ),
        'rho_condition': '0 < rho < |c|/6',
        'rho_bound': Abs(c_val) / 6,
        'completion': 'weight-completed/pro Banach completion',
        'brown_henneaux_dictionary': 'c = 3 ell / (2 G_N)',
        'physical_reading': 'boundary-CFT holographic interpretation',
        'dynamical_metric_path_integral_constructed': False,
        'required_physical_hypotheses': (
            'BHdict',
            'ModInv',
            'VacDom',
            'SadDom',
            'Borel/Stokes',
        ),
    }


def brown_henneaux_chiral_chart_scope() -> Dict[str, Any]:
    r"""Scope of the Brown-Henneaux chart used by the gravity chapter.

    Brown-Henneaux supplies an external asymptotic-symmetry theorem and
    Witten supplies the Chern-Simons comparison.  The bar complex uses
    the resulting Virasoro boundary algebra after this chart is fixed;
    it does not derive Newton's constant or a metric path integral.
    """
    return {
        'claim_status': 'ProvedElsewhere',
        'boundary_chart': 'b_BH',
        'source_theorem': 'Brown-Henneaux asymptotic charge algebra',
        'cs_comparison': 'one oriented SL(2,R) Chern-Simons factor',
        'level': 'k = ell / (4 G_N)',
        'chiral_central_charge': 'c_ch = 6k = 3 ell / (2 G_N)',
        'full_ads3_asymptotic_symmetry': 'Vir_c_ch direct_sum Vir_c_ch',
        'licensing_tags': ('alpha', 'beta', 'gamma'),
        'requires': (
            'hypBHdict',
            'Brown-Henneaux asymptotic boundary conditions',
            'one oriented SL(2,R) Chern-Simons factor',
            'Chern-Simons/Virasoro boundary-charge comparison',
            'orientation convention for the two chiral factors',
        ),
        'used_as': 'external boundary chart A_{b_BH} = Vir_c_ch',
        'derived_from_bar_complex': False,
        'dynamical_metric_path_integral_constructed': False,
    }


def gravitational_mc_bridge_scope() -> Dict[str, Any]:
    r"""Scope of the Brown-Henneaux Virasoro MC bridge.

    The MC package is internal only after the external Brown-Henneaux
    chart, the physics-bridge BV/parametrix hypotheses, the completed
    ambient, and the line-side Koszul-effectiveness hypothesis have all
    been fixed.  Its boundary face is Virasoro; the closed/bulk object
    is the derived chiral centre after the Hochschild comparison.
    """
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'delta', 'epsilon'),
        'boundary_chart': 'b_BH',
        'boundary_face': 'A_{b_BH} = Vir_{c_ch}',
        'chiral_central_charge': 'c_ch = 6k = 3 ell / (2 G_N)',
        'mc_package': 'alpha_grav in MC(g_SC_T)',
        'requires': (
            'hypBHdict',
            'physics-bridge BV datum',
            'factorized logarithmic parametrix',
            'one-loop finiteness',
            'product-polynomial interactions',
            'hypAmbientWtCpl',
            'effKoszul',
            'completed perturbative line category near b_BH',
        ),
        'licensed_projections': (
            'boundary/chiral Virasoro face',
            'completed open-colour line module category',
            'abstract spectral line braiding',
            'bar-intrinsic Virasoro genus-shadow MC element',
            'Virasoro PVA lambda-bracket shadow',
        ),
        'collision_kernel': 'r_coll_Vir(z) = (c_ch/2)/z^3 + 2T/z',
        'bulk_after_comparison': 'Z_der^ch(Vir_{c_ch})',
        'closed_face_is_boundary_virasoro': False,
        'not_asserted': (
            'unconditional Chern-Simons-to-SC bridge without BV/parametrix data',
            'closed face equals Vir_c',
            'line category is Vir_{26-c}-modules',
            'spectral braiding is the Virasoro fusion kernel',
            'BTZ black holes are proved MC deformations',
            'physical torus modular invariance or vacuum dominance',
        ),
    }


def brown_henneaux_line_test_package_scope() -> Dict[str, Any]:
    r"""Scope of the Brown-Henneaux algebraic line test package.

    The package Perf(Vir_c) with the transferred generator coproduct and
    collision-residue kernel is a conditional comparison model for the
    Brown-Henneaux boundary chart.  It is not the full gravitational
    open-sector factorization category and it does not construct the
    level-A gravity-line operator algebra.
    """
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'delta', 'epsilon'),
        'package': 'C_grav^test',
        'ambient': 'hypAmbientWtCpl',
        'boundary_chart': 'b_BH',
        'boundary_algebra': 'A_{b_BH} = Vir_c',
        'bulk_object': 'Zder^ch(Vir_c)',
        'central_projection': 'C[[c]]',
        'central_projection_is_full_bulk': False,
        'primitive_coproduct_scope': 'generator-level Virasoro T',
        'primitive_coproduct': 'Delta_z^Vir(T) = tau_z(T) tensor 1 + 1 tensor T',
        'all_degree_primitivity_status': 'conditional on signed ghost-defect lemma',
        'line_model': 'highest-weight Vir_{26-c}-modules',
        'line_model_status': 'Conjectural via gravity-line-identification',
        'annulus_trace': 'q^(-c/24) product_{n>=2} (1-q^n)^(-1)',
        'annulus_trace_scope': 'vacuum one-boundary scalar character',
        'physical_torus_partition_function_statement': False,
        'genus_shadow': 'F_g(Vir_c) = kappaChHodge(Vir_c) lambda_g^FP',
        'physical_gravity_requires': (
            'ModInv',
            'VacDom',
            'SadDom',
            'Cardy/BTZ comparison',
        ),
        'constructs_level_A_gravity_line_operator_algebra': False,
        'not_asserted': (
            'Perf(Vir_c) is the full gravitational open-sector factorization category',
            'C[[c]] is the full derived chiral centre',
            'all elements have a proved primitive coproduct',
            'Vir_{26-c}-modules are the proved line category',
            'vacuum character is the physical torus partition function',
            'Pentagon-face scalar trace Phi_10^un is constructed',
        ),
    }


def virasoro_bar_intrinsic_mc_shadow_scope() -> Dict[str, Any]:
    r"""Scope of the bar-intrinsic Virasoro MC shadow.

    Vol I constructs the positive-genus bar-intrinsic MC element
    Theta_A = D_A - d0.  For the Virasoro boundary chart this is the
    completed stable-graph correction.  Genus-zero A-infinity
    operations belong to d0; finite shadow coefficients are projections
    of Theta, not replacements for it.
    """
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'epsilon'),
        'requires': (
            'hypBHdict',
            'hypAmbientWtCpl',
            'effKoszul',
            'Vol I bar-intrinsic MC2 theorem',
            'non-degenerate invariant form',
        ),
        'total_coderivation': 'D_Vir_c = d0 + sum_{g>=1} hbar^g d_Vir_c^(g)',
        'mc_element': 'Theta_grav = Theta_Vir_c = D_Vir_c - d0',
        'mc_equation': '[d0,Theta_grav] + 1/2[Theta_grav,Theta_grav] = 0',
        'completed_target': 'Def_cyc(Vir_c) completed_tensor G_mod',
        'genus_zero_operations_in_theta': False,
        'genus_zero_role': 'm_k are part of d0 and satisfy Stasheff identities',
        'genus_one_scalar_trace': 'Theta_grav^(1,sc) = kappaChHodge(Vir_c) eta tensor Lambda',
        'f1': 'F_1(Vir_c) = c/48',
        'uniform_weight_scalar_lane': (
            'F_g(Vir_c) = kappaChHodge(Vir_c) lambda_g^FP'
        ),
        'non_scalar_higher_genus_component': 'stable-graph coderivation d_Vir_c^(g)',
        'shadow_projections': (
            'kappaChHodge(Vir_c)=c/2',
            'C_Vir=-c',
            'Q_contact_Vir=10/(c(5c+22))',
            'o5(Vir_c) nonzero generically on residue-survival surface',
        ),
        'central_projection_is_full_bulk': False,
        'physical_partition_function_statement': False,
        'not_asserted': (
            'Theta_grav includes sum_{k>=2} alpha_k',
            'finite shadow projections replace the full MC element',
            'scalar genus formula is the full non-scalar genus component',
            'C[[c]] is equivalent to the full derived chiral centre',
            'physical torus partition function follows from MC2',
        ),
    }


def genus0_directed_product_decomposition_scope() -> Dict[str, Any]:
    r"""Scope of the genus-zero directed product decomposition.

    The genus-zero A_infinity-chiral operad is a strict directed product
    of operation spaces after the local HT slab chart is fixed.  The product
    is directed by the colour order ch <= bdy <= tr; it is not the ordinary
    product of two coloured operads on disjoint colour sets.
    """
    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('alpha', 'gamma'),
        'requires': (
            'genus-0 HT slab chart',
            'hypAmbientWtCpl',
        ),
        'product_kind': 'directed colour-filtered product',
        'product_notation': 'SCchtop boxtimes_dir E1^tr',
        'ordinary_disjoint_product': False,
        'color_order': ('ch', 'bdy', 'tr'),
        'nonempty_rule': 'every input colour is <= output colour',
        'operation_spaces': {
            'ch_output': 'FM_k(C)',
            'bdy_output': 'FM_k(C) x E1(m)',
            'tr_output': 'FM_k(C) x E1(m) x E1(p)',
        },
        'mixed_inputs_to_tr_present': True,
        'strict_after_chart_choice': True,
        'composition': 'componentwise FM and Stasheff insertions',
        'algebra_scope': (
            'A_infinity algebra object in the module category over the '
            'SCchtop algebra (A_ch,A_bdy)'
        ),
        'not_asserted': (
            'ordinary product of two coloured operads on disjoint colour sets',
            'transverse colour isolated from chiral and boundary inputs',
            'higher-genus product decomposition',
        ),
    }


def quantum_group_clutching_equivariance_scope() -> Dict[str, Any]:
    r"""Scope of the braided-annular graph-equivariance proposition.

    Quasi-triangularity alone intertwines a coproduct with the opposite
    coproduct.  It does not make a braided category symmetric, and it does
    not imply the inverse-orientation relation for a reversed annular seam.
    Exact stable-graph equivariance therefore uses the full annular sewing
    datum recorded here.
    """
    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('alpha', 'beta', 'gamma'),
        'requires': (
            'stable-graph sewing orientation',
            'quasi-triangular comparison R Delta = Delta^op R',
            'Yang-Baxter braid coherence',
            'pure-braid descent to Aut(Gamma)',
            'inverse-orientation relation R_21(-z) R_12(z) = id',
            'hypAmbientWtCpl',
        ),
        'endpoint_exchange_mechanism': 'quasi-triangular coproduct comparison',
        'edge_permutation_mechanism': 'positive braid lifts; YBE gives reduced-word coherence',
        'aut_descent_mechanism': 'pure-braid ambiguity acts trivially after sewing',
        'self_edge_flip_mechanism': 'orientation reversal uses inverse braiding',
        'literal_commutativity_of_monodromy_products': False,
        'inverse_orientation_from_quasitriangularity': False,
        'braiding_is_symmetric_group_action_before_descent': False,
        'general_ainf_chiral_scope': (
            'hypothesis unless the modular bar datum supplies the same '
            'braided annular compatibility'
        ),
        'not_asserted': (
            'quasi-triangularity plus YBE alone imply Aut(Gamma)-equivariance',
            'monodromy factors can be freely permuted as commuting products',
            'R_21(z) = R(z)^-1 follows from quasi-triangularity',
            'a braided category has an ordinary symmetric-group action before descent',
        ),
    }


def modular_operad_unitality_scope() -> Dict[str, Any]:
    r"""Scope of the unit-normalised annular sewing theorem.

    The modular unit is an unstable genus-zero identity component adjoined
    to the stable-graph category.  Unitality is proved by the strict unit,
    the counit-normalised annular seam, and the diagonal bicomodule identity;
    it is not proved by treating the unit component as a simply-connected
    stable vertex.
    """
    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('alpha', 'beta', 'gamma'),
        'requires': (
            'unit-extended stable-graph chart',
            'strict colour identities in SCchtop boxtimes_dir E1^tr',
            'counit-normalised annular sewing',
            'diagonal bicomodule C_Delta',
            'R-matrix unit normalisation (epsilon tensor id)R=1=(id tensor epsilon)R',
            'hypAmbientWtCpl',
        ),
        'exceptional_vertex_stable': False,
        'exceptional_vertex_kind': 'unstable genus-0 two-flag identity vertex',
        'unit_colours': {
            'ch': 'eta_ch in FM_1(C)',
            'bdy': 'eta_bdy in FM_0(C) x E1(1)',
            'tr': 'eta_tr in FM_0(C) x E1(0) x E1(1)',
        },
        'annular_unit_seam': 'C_Delta identity bicomodule after unit/counit',
        'cotensor_identity': (
            'M square_C C_Delta = M = C_Delta square_C M'
        ),
        'coproduct_counit_axioms': (
            '(epsilon tensor id) Delta_z = id',
            '(id tensor epsilon) Delta_z ~= tau_z',
            'tau_0 = id on the exceptional-edge chart',
        ),
        'monodromy_on_unit': 'Mon(R)_{eta,e}=id by counit normalisation',
        'uses_contractible_disk_argument': False,
        'proves_general_clutching_existence': False,
        'proves_composition_associativity': False,
        'independent_of_genus_and_shadow_class_after_normalisation': True,
        'not_asserted': (
            'the unit vertex is a stable vertex',
            'the unit seam is trivial because a punctured sphere is simply connected',
            'unitality constructs the general positive-genus clutching maps',
            'composition associativity follows from unitality',
        ),
    }


def modular_bar_reduction_scope() -> Dict[str, Any]:
    r"""Scope of the annular-data to modular-bar-datum reduction.

    The stable-graph modular bar theorem is an uncurved square-zero
    statement.  A positive-genus Arakelov curvature computation is not an
    internal differential for that theorem until an S-tail/twisting datum has
    absorbed the curvature.
    """
    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'epsilon'),
        'requires': (
            'hypAmbientWtCpl',
            'nodal bar complex data',
            'annular one-edge expansion maps',
            'square-zero internal differential or S-tail curvature datum',
            'anticommutation with one-edge expansions',
            'codimension-two cancellation',
        ),
        'abstract_theorem': 'thm:modular-bar',
        'abstract_theorem_requires_internal_square_zero': True,
        'positive_genus_curvature_is_modular_bar_datum_by_itself': False,
        'curvature_term': 'kappaChHodge(A) omega_g',
        'curvature_absorption': 'S-tail/twisting maps D_S with d_S^2=0',
        'one_edge_maps_degree': 1,
        'one_edge_maps_source': 'annular bar complex with R-matrix monodromy',
        'covered_loci': (
            'genus 0 affine non-critical',
            'all genera affine integrable with KZB regular singularity',
            'Heisenberg scalar monodromy branch',
        ),
        'open_loci': (
            'general A_infinity-chiral genus >= 1 one-edge construction',
            'generic non-integral affine genus >= 1 Stokes compatibility',
            'positive-genus curvature absorption without supplied S-tail data',
        ),
        'codimension_two_cases': (
            'disjoint collars: Koszul sign cancellation',
            'two separating edges: KZ/KZB pentagon on covered loci',
            'one nonseparating edge: KZ/KZB hexagon on covered loci',
        ),
        'proves_general_positive_genus_clutching_maps': False,
        'proves_concrete_operadic_associativity_by_itself': False,
        'not_asserted': (
            'd^2 = kappa omega_g is an internal differential for thm:modular-bar',
            'Theorem modular-bar alone proves concrete clutching associativity',
            'annular maps exist for every A_infinity-chiral algebra at genus >= 1',
            'generic non-integral KZB Stokes compatibility is proved here',
        ),
    }


def affine_modular_bar_datum_scope() -> Dict[str, Any]:
    r"""Scope of the affine modular-bar datum corollary.

    The affine KZ/KZB input gives a modular bar datum only on the loci where
    the annular one-edge maps satisfy the four modular-bar criterion axioms.
    Genus zero uses the logarithmic KZ connection at every non-critical level;
    positive genus uses the regular-singular KZB extension at integrable level.
    Generic non-integral positive genus still has the Stokes compatibility gap.
    """
    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'epsilon'),
        'requires': (
            'hypAmbientWtCpl',
            'finite-dimensional simple Lie algebra g',
            'non-critical level k != -h^vee',
            'affine nodal bar complex data',
            'annular KZ/KZB one-edge maps',
            'square-zero internal differential after curvature compensation',
            'KZ/KZB codimension-two pentagon and hexagon',
        ),
        'criterion_used': 'prop:ainf-chiral-modular-bar-reduction',
        'abstract_theorem': 'thm:modular-bar',
        'covered_loci': (
            'genus 0 at every non-critical level',
            (
                'all genera at integrable level in the semisimple integrable '
                'positive-energy module category'
            ),
        ),
        'genus_zero_mechanism': (
            'uncurved internal differential; KZ logarithmic regular singular '
            'connection; Drinfeld pentagon and hexagon'
        ),
        'integrable_all_genera_mechanism': (
            'regular-singular KZB extension plus S-tail/twisting curvature '
            'compensation'
        ),
        'positive_genus_generic_nonintegral_claimed': False,
        'stokes_gap_closed': False,
        'proves_modular_operad_axioms_by_itself': False,
        'corollary_supplies': (
            'modular bar datum on covered loci',
            'D^2=0 on Bmod(V_k(g)) on covered loci',
            'complete conilpotent modular dg coalgebra on covered loci',
        ),
        'not_asserted': (
            'generic non-integral all-genera affine modular bar datum',
            'KZB Stokes compatibility at genus >= 1',
            'codimension-two cancellation at generic non-integral positive genus',
            'full modular operad axioms without pure-braid descent and inverse orientation',
        ),
    }


def affine_e3_topological_km_scope() -> Dict[str, Any]:
    r"""Scope of the affine Kac-Moody E3-topological theorem.

    The proved affine statement is the cohomological topologisation of the
    derived chiral centre.  Non-criticality makes the Sugawara antighost
    denominator invertible.  It gives holomorphic translations that vanish
    on Q_CS-cohomology, hence local constancy in the holomorphic plane and
    an E_3^top structure after the single-colour external-product/Dunn step.
    A strict raw-chain-level identity on the original BV complex is a
    separate frontier problem.
    """
    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'epsilon'),
        'requires': (
            'finite-dimensional simple Lie algebra g',
            'non-critical level k != -h^vee',
            'holomorphic Chern-Simons BV factorisation algebra on X x R',
            'boundary chart V_k(g)',
            'Q_CS-cohomology ambient',
            'Sugawara antighost primitive on cohomology',
            'local constancy recognition for factorisation algebras',
            'single-colour external product before Dunn additivity',
        ),
        'sugawara_denominator': '2(k + h^vee)',
        'critical_level_excluded': True,
        'brst_identity_scope': 'on Q_CS-cohomology',
        'translation_class': 'partial_z is zero in Q_CS-cohomology',
        'topological_structure_target': 'H^bullet_{Q_CS} Zder^ch(V_k(g))',
        'e2_upgrade': 'E2^hol becomes E2^top on Q_CS-cohomology',
        'e3_mechanism': 'E2^top tensor E1^top -> E3^top after external product',
        'applies_dunn_to_bicoloured_scchtop': False,
        'strict_raw_chain_level_statement': False,
        'strict_raw_chain_level_frontier': 'rem:frontier-class-L-strict-chain-level',
        'shifted_cs_level': 'k + h^vee',
        'not_asserted': (
            'strict raw-chain-level E3^top structure before Q_CS-cohomology',
            'critical-level Sugawara topologisation',
            'Dunn additivity applied directly to the bicoloured SC^{ch,top} operad',
            'general conformal chiral algebra topologisation',
        ),
    }


def principal_ds_e3_topological_scope() -> Dict[str, Any]:
    r"""Scope of the principal DS E3-topological theorem.

    The principal DS statement is a cohomological topologisation theorem
    over the total DS BRST differential.  The formal topologisation follows
    from a DS primitive for the stress tensor, but that primitive is data of
    the DS reduction package; it is not the assertion that Cartan currents
    are Q_CS-boundaries in the unreduced holomorphic Chern-Simons bulk.
    """
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'epsilon'),
        'hypothesis_package': 'hypDSBRST',
        'requires': (
            'finite-dimensional simple Lie algebra g',
            'principal nilpotent f_prin',
            'non-critical level k != -h^vee',
            'Costello-Gaiotto DS boundary condition',
            'boundary factorisation algebra identified with W^k(g,f_prin)',
            'complete filtered DS ghost ambient',
            'total DS differential Q_DS = Q_CS + Q_red',
            'DS antighost primitive G_DS prime on Q_DS-cohomology',
            'single-colour external product before Dunn additivity',
        ),
        'brst_identity_scope': 'on Q_DS-cohomology',
        'stress_identity': 'T_DS = [Q_DS, G_DS_prime] on Q_DS-cohomology',
        'topological_structure_target': 'H^bullet_{Q_DS} Zder^ch(W^k(g,f_prin))',
        'e2_upgrade': 'E2^hol becomes E2^top on Q_DS-cohomology',
        'e3_mechanism': 'E2^top tensor E1^top -> E3^top after external product',
        'applies_dunn_to_bicoloured_scchtop': False,
        'strict_raw_chain_level_statement': False,
        'unreduced_cartan_current_qcs_exactness': False,
        'normalisation_dependent_cartan_antighost_coefficient': True,
        'not_asserted': (
            'Cartan currents are Q_CS-exact in the unreduced hCS bulk',
            'a universal invariant 1/4 Cartan-antighost coefficient',
            'strict raw-chain-level E3^top structure before Q_DS-cohomology',
            'non-principal DS families without their separate good-grading package',
        ),
    }


def good_graded_ds_e3_topological_scope() -> Dict[str, Any]:
    r"""Scope of the good-graded/non-principal DS E3-topological theorem."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'epsilon'),
        'hypothesis_package': 'hypDSBRST',
        'requires': (
            'finite-dimensional simple Lie algebra g',
            'nilpotent f with good DS grading',
            'non-critical level k != -h^vee',
            'Costello-Gaiotto DS boundary condition for f',
            'boundary factorisation algebra identified with W^k(g,f)',
            'complete filtered or branched-cover DS ghost ambient',
            'total DS differential Q_DS,f = Q_CS + Q_red,f',
            'DS antighost primitive G_DS,f prime on Q_DS,f-cohomology',
        ),
        'brst_identity_scope': 'on Q_DS,f-cohomology',
        'stress_identity': 'T_DS(f) = [Q_DS,f, G_DS,f_prime] on Q_DS,f-cohomology',
        'topological_structure_target': 'H^bullet_{Q_DS,f} Zder^ch(W^k(g,f))',
        'strict_raw_chain_level_statement': False,
        'unreduced_affine_current_qcs_exactness': False,
        'normalisation_dependent_improvement_coefficients': True,
        'not_asserted': (
            'Q_CS-exactness of all affine currents in the unreduced hCS bulk',
            'a universal invariant Cartan-antighost coefficient for BP or subregular DS',
            'generic DS-Hochschild chain-level transport',
            'non-good-graded nilpotent topologisation',
        ),
    }


def virasoro_scalar_bar_trace_profile(c=None) -> Dict[str, Any]:
    r"""Scalar boundary trace associated to the completed Virasoro bar.

    This is the algebraic scalar trace

      Z_bar^grav(hbar;rho)
        = Tr_{B_ch(Vir_c)^hat_rho} exp(sum_{g>=0} hbar^g Theta^{(g)}).

    It is not automatically the gravitational path integral.
    """
    c_val = Symbol('c') if c is None else S(c)
    return {
        'formula': (
            'Z_bar^grav(hbar;rho) = '
            'Tr_{B_ch(Vir_c)^hat_rho} '
            'exp(sum_{g>=0} hbar^g Theta_Vir_c^(g))'
        ),
        'trace_space': 'B_ch(Vir_c)^hat_rho',
        'rho_bound': Abs(c_val) / 6,
        'object_type': 'scalar boundary trace',
        'is_metric_path_integral': False,
        'needs_for_btz_or_desitter_saddle': (
            'BHdict',
            'ModInv',
            'VacDom',
            'SadDom',
            'Borel/Stokes',
        ),
    }


def genus1_virasoro_mc_scope(c=None) -> Dict[str, Any]:
    r"""Scope guard for the genus-1 Virasoro MC scalar and Ward package.

    The completed modular-bar scalar lane fixes the self-sewing
    coefficient

        Theta_0^(1) = kappaChHodge(Vir_c) omega_1 = (c/2) omega_1

    and hence the Hodge trace F_1 = c/48.  The torus Ward-block
    normalization is separate: the generic vacuum character gives
    constant one-point term -c/24, and the two-point central singular
    kernel is (c/2) P_4, not (c/2) wp.
    """
    c_val = Symbol('c') if c is None else S(c)
    kappa = c_val / 2
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': ('hypAmbientWtCpl', 'hypVirTorusBlock'),
        'stable_graph_expression': 'completed genus-1 modular bar expansion',
        'scalar_self_sewing': (
            'Theta_0^(1)=kappaChHodge(Vir_c) omega_1=(c/2) omega_1'
        ),
        'mode_coefficient': 'T_(3)T=c/2',
        'kappa': simplify(kappa),
        'hodge_trace': simplify(kappa / 24),
        'vacuum_one_point_constant': simplify(-c_val / 24),
        'vacuum_character': 'q^(-c/24) prod_{m>=2}(1-q^m)^(-1)',
        'one_point_formula': (
            'A_1^(1)(T;tau)=q d_q log Z_1^Vir(tau)'
        ),
        'two_point_singular_part': (
            '(c/2) P4(z,tau) + 2 P2(z,tau) A_1^(1)(T;tau)'
        ),
        'two_point_central_kernel': '(c/2) P4(z,tau)',
        'P2_definition': 'P2=wp+E2/12',
        'P4_definition': 'P4=(1/6) d_z^2 P2',
        'P4_local_expansion': 'P4=z^(-4)+O(1)',
        'full_all_n_stable_graph_series_evaluated': False,
        'regular_torus_blocks_determined_by_kappa': False,
        'central_kernel_is_wp': False,
        'casimir_equals_scalar_hodge_trace': False,
        'not_asserted': (
            'the full all-n stable graph series is evaluated by kappa',
            'regular torus blocks are determined by the scalar Hodge trace',
            'the two-point central kernel is (c/2) wp',
            'the Casimir term -c/24 equals the Hodge trace c/48',
            'sum_n (c/12)(n^3-n) q^n derives the scalar genus-1 trace',
            'E2/24 times T is the genus-1 one-point Ward block',
        ),
    }


def genus1_virasoro_amplitudes_scope(c=None) -> Dict[str, Any]:
    r"""Scope guard for the genus-1 Virasoro Ward-amplitude theorem.

    The theorem separates the scalar Hodge trace, determinant-line
    shadow, and chiral Virasoro Ward block.  The universal data are the
    Ward-recursive singular terms.  Regular torus blocks, minimal-model
    quotient effects, and physical non-chiral partition functions are
    not determined by the scalar Hodge trace.
    """
    c_val = Symbol('c') if c is None else S(c)
    kappa = c_val / 2
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': ('hypAmbientWtCpl', 'hypVirTorusBlock'),
        'scalar_hodge_trace': simplify(kappa / 24),
        'determinant_line_derivative_coefficient': simplify(-kappa / 24),
        'vacuum_one_point_constant': simplify(-c_val / 24),
        'generic_vacuum_character': 'q^(-c/24) prod_{n>=2}(1-q^n)^(-1)',
        'generic_vacuum_character_product_starts_at': 2,
        'P1_local_expansion': 'P1=z^(-1)+O(z)',
        'P2_definition': 'P2=wp+E2/12',
        'P4_definition': 'P4=(1/6) d_z^2 P2',
        'two_point_simple_pole_derivative_vanishes': True,
        'two_point_regular_block_from_hypothesis': True,
        'three_point_requires_p1_derivative': True,
        'three_point_pairwise_singular_terms': (
            '(c/2) P4(z_ij,tau) A1 + '
            '2 P2(z_ij,tau) A2(z_j,z_k) + '
            'P1(z_ij,tau) partial_z_j A2(z_j,z_k)'
        ),
        'regular_blocks_determined_by_kappa': False,
        'minimal_model_quotients_included': False,
        'physical_torus_partition_function_statement': False,
        'not_asserted': (
            'A1(T;tau)=-(c/24) E2(tau) for the generic Virasoro vacuum module',
            'the three-point singular recursion has no simple-pole derivative term',
            'regular torus blocks are determined by kappaChHodge',
            'minimal-model singular-vector quotients have the generic vacuum character',
            'the chiral Ward block is the full non-chiral torus partition function',
        ),
    }


def genus1_modular_anomaly_scope(c=None) -> Dict[str, Any]:
    r"""Scope guard for the genus-1 E2 modular-anomaly proposition.

    The standard modular-form statement is unconditional:

        E2(-1/tau) = tau^2 E2(tau) + 12 tau/(2*pi*i),
        E2hat = E2 - 3/(pi Im tau) is modular of weight 2.

    The bar-theoretic identification of the non-holomorphic correction
    with an algebraic transgression neutralizing kappa*omega_1 is a
    comparison datum.  The E2 anomaly is not kappa^{-1}; the kappa
    factor enters by multiplying the determinant-line derivative.
    """
    c_val = Symbol('c') if c is None else S(c)
    kappa = c_val / 2
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': ('hypAmbientWtCpl', 'hypEisComp'),
        'E2_transformation': 'E2(-1/tau)=tau^2 E2(tau)+12 tau/(2 pi i)',
        'E2_anomaly_coefficient': S(12),
        'E2_anomaly_is_universal': True,
        'E2_anomaly_is_kappa_inverse': False,
        'E2hat': 'E2 - 3/(pi Im tau)',
        'E2hat_is_modular_weight_2': True,
        'q_dq_log_eta': 'E2/24',
        'determinant_line_derivative_coefficient': simplify(-kappa / 24),
        'bar_curvature_class': 'kappaChHodge(Vir_c) omega_1',
        'bar_curvature_coefficient': simplify(kappa),
        'nonholomorphic_correction_requires_comparison': True,
        'nonholomorphic_correction_is_transgression_generator': False,
        'eta_square_equals_bar_curvature_asserted': False,
        'not_asserted': (
            'the E2 anomalous term is kappa^{-1} times the nonseparating class',
            'the scalar function -3/(pi Im tau) is literally the algebra generator eta',
            'eta^2 equals kappaChHodge(Vir_c) omega_1 without comparison data',
            'the non-holomorphic completion constructs the chain-level bar transgression',
        ),
    }


def genus1_virasoro_kzb_shadow_kernel_scope(c=None) -> Dict[str, Any]:
    r"""Scope guard for the genus-1 Virasoro KZB shadow kernel.

    The theorem records the connection-level KZB shadow

        r_1^KZB,sh = (c/2) P2(z,tau) + C_T^KZB(tau) T,

    with P2=wp+E2/12 and with the T-contact coefficient fixed only after
    the KZB heat-operator normalization has been chosen.  It is not the
    construction of a full quantum R-matrix, and it is not the torus
    stress-tensor central singular kernel, which uses P4.
    """
    c_val = Symbol('c') if c is None else S(c)
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': ('hypAmbientWtCpl', 'hypKZBHeat'),
        'object': 'connection-level KZB shadow kernel',
        'elliptic_scalar_kernel': '(c/2) P2(z,tau)',
        'scalar_kernel_coefficient': simplify(c_val / 2),
        'P2_definition': 'P2=wp+E2/12',
        'contact_term': 'C_T^KZB(tau) T',
        'chosen_contact_normalization': 'C_T^KZB(tau)=2 E2(tau)',
        'contact_term_universal_without_kzb_heat': False,
        'linearized_kzb_flatness_only': True,
        'full_quantum_R_matrix_constructed': False,
        'nonlinear_qybe_proved': False,
        'higher_hbar_corrections_determined': False,
        'central_tt_kernel_is_P4_not_P2': True,
        'not_asserted': (
            'a full quantum R^mod(z;hbar,tau) is constructed',
            'the nonlinear quantum Yang-Baxter equation is proved',
            'the contact term 2T E2 is universal without hypKZBHeat',
            'the torus stress-tensor central singular kernel is (c/2)P2',
            'higher hbar^(2g) corrections are determined',
        ),
    }


def gravity_weinberg_ward_residue_scope() -> Dict[str, Any]:
    r"""Scope guard for the degree-2 Ward residue and Weinberg comparison.

    The algebraic calculation is the Virasoro primary Ward residue

        W_i^(j) = Res (z-z_i)^j T(z) O_i dz = T_(j) O_i,

    so the simple-pole residue is W_i^(0)=partial_i.  The physical
    Weinberg factor is a spin-two momentum-space image only after the
    celestial soft comparison datum has installed the Mellin residue,
    asymptotic states, polarization data, and BMS normalization.
    """
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': ('hypAmbientWtCpl', 'hypCelSoft'),
        'algebraic_object': 'degree-2 Virasoro Ward residue',
        'shadow_connection_kernel': (
            'sum_i (h_i/(z0-z_i)^2 + partial_{z_i}/(z0-z_i))'
        ),
        'local_residue_dictionary': (
            'W_i^(j)=Res (z-z_i)^j T(z) O_i dz = T_(j) O_i'
        ),
        'simple_pole_residue': 'W_i^(0)=partial_{z_i}',
        'double_pole_residue': 'W_i^(1)=h_i',
        'physical_weinberg_factor': (
            'sum_i epsilon_{mu nu} p_i^mu p_i^nu/(q dot p_i)'
        ),
        'physical_factor_is_spin_two': True,
        'requires_celestial_bms_comparison': True,
        'algebraic_residue_equals_physical_weinberg_without_comparison': False,
        'translation_invariance_derives_spin_two_factor': False,
        'old_mode_shift_T_p_plus_1_retired': True,
        'not_asserted': (
            'T_(p+1) is the simple-pole Ward residue',
            'translation invariance derives the spin-two Weinberg factor',
            'the algebraic residue alone contains the soft energy denominator',
            'the physical factor is the photon-like p.epsilon/(p.q)',
        ),
    }


def gravity_cachazo_strominger_ward_package_scope() -> Dict[str, Any]:
    r"""Scope guard for the degree-2 Ward package entering subleading soft.

    The degree-2 Virasoro calculation gives the holomorphic
    stress-tensor Ward operator

        W_Y = sum_i (Y_i partial_i + h_i partial Y_i).

    The double-pole residue is only the conformal-weight term.  The
    physical Cachazo-Strominger factor requires the celestial/BMS
    comparison, the antiholomorphic/polarization data, and the
    degree-3 shadow channel of the soft hierarchy.
    """
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': ('hypAmbientWtCpl', 'hypCelSoft'),
        'algebraic_object': 'degree-2 holomorphic stress-tensor Ward package',
        'ward_operator': 'sum_i (Y(z_i) partial_{z_i} + h_i partial Y(z_i))',
        'double_pole_residue': 'W_i^(1)=h_i',
        'simple_pole_residue_required': 'W_i^(0)=partial_{z_i}',
        'double_pole_role': 'conformal-weight term only',
        'physical_cachazo_strominger_factor': (
            'sum_i epsilon_{mu nu} p_i^mu q_rho J_i^{rho nu}/(q dot p_i)'
        ),
        'requires_antiholomorphic_completion': True,
        'requires_degree_three_shadow_channel': True,
        'double_pole_alone_equals_physical_cs': False,
        'global_conformal_ward_identity_equals_physical_cs_without_comparison': False,
        'not_asserted': (
            'h_i/(z0-z_i)^2 is the Cachazo-Strominger soft factor',
            'the double-pole residue alone gives the subleading soft theorem',
            'L0 eigenvalue alone is the angular momentum operator',
            'the degree-2 Ward package fixes the momentum-space denominator',
        ),
    }


def gravity_chy_quartic_contact_scope() -> Dict[str, Any]:
    r"""Scope guard for the sub-subleading Virasoro quartic contact channel.

    The binary Ward kernel has no T_(2) residue on a primary.  The
    ternary operation m3 is the cubic precursor, but the Virasoro
    sub-subleading graviton shadow is the degree-4 contact channel
    together with the cubic self-interaction in the degree-4 flatness
    equation.  The physical Cachazo--Strominger/Hamada--Shiu reading
    requires the celestial soft comparison datum.
    """
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': ('hypAmbientWtCpl', 'hypCelSoft'),
        'algebraic_object': 'degree-4 quartic contact sub-subleading shadow',
        'binary_primary_T2_residue': 'T_(2) O_i = 0',
        'ternary_role': 'm3 is cubic precursor entering [C,C] at degree 4',
        'virasoro_cubic_shadow_gauge_trivial': True,
        'first_nontrivial_virasoro_graviton_channel': 'S4 quartic contact',
        'quartic_contact_coefficient': '10/(c(5c+22))',
        'large_c_asymptotic': '2/c^2',
        'degree_four_shadow_formula': 'Sh(Q)|soft + [Sh(C),Sh(C)]|soft',
        'physical_chy_hsl_requires_celsoft': True,
        'physical_cs_hamada_shiu_requires_celsoft': True,
        'm3_alone_equals_physical_chy': False,
        'm3_alone_equals_physical_subsubleading': False,
        'binary_residue_equals_subsubleading': False,
        'not_asserted': (
            'm3 alone is the sub-subleading soft theorem',
            'the ternary operation is an independent physical sub-subleading theorem',
            'the binary Ward kernel contributes T_(2) on primaries',
            'quartic contact is a physical 4d soft theorem without hypCelSoft',
        ),
    }


def virasoro_shadow_metric_coefficients(max_r: int = 5, c=None) -> Dict[int, Any]:
    r"""Normalized Virasoro scalar shadow coefficients from the metric branch.

    The branch is fixed by ``sqrt(Q_Vir(0)) = c``:

        H(t) = t^2 c sqrt(1 + 12 t/c
               + (36 + 80/(5c+22)) t^2/c^2),
        S_r = [t^r]H(t)/r.

    This computes the scalar/support-shadow normalization.  It is not
    the raw transferred operation m_r.
    """
    if max_r < 2:
        raise ValueError("shadow coefficients start at arity r=2")
    c_val = Symbol('c') if c is None else S(c)
    t = Symbol('t')
    q_branch = 1 + (S(12) / c_val) * t + (
        S(36) + S(80) / (5 * c_val + 22)
    ) * t**2 / c_val**2
    h_series = series(t**2 * c_val * sqrt(q_branch), t, 0, max_r + 1).removeO()
    return {
        r: simplify(expand(h_series).coeff(t, r) / r)
        for r in range(2, max_r + 1)
    }


def virasoro_catalan_shape_factor(r: int, x=None):
    r"""Catalan shape factor in the normalized Virasoro scalar branch.

    For \(r\ge4\),

        F_r(x) = sum_j (-1)^j C_j binom(r-4, 2j) x^j.

    This is the shape factor in the scalar \(T\)-line coefficient after
    scalar shadow projection; it is not an unprojected ordered-bar
    invariant.
    """
    if r < 4:
        raise ValueError("the Catalan shape factor is defined for r >= 4")
    x_val = Symbol('x') if x is None else S(x)
    total = S.Zero
    for j in range((r - 4) // 2 + 1):
        catalan_j = factorial(2 * j) / (factorial(j) * factorial(j + 1))
        total += (-1)**j * catalan_j * binomial(r - 4, 2 * j) * x_val**j
    return simplify(total)


def virasoro_shadow_closed_form_coefficient(r: int, c=None):
    r"""Closed-form normalized scalar \(T\)-line shadow coefficient.

    The branch is

        Q_Vir(t) = (c+6t)^2 + 80 t^2/(5c+22),
        H(t) = t^2 sqrt(Q_Vir(t)),  sqrt(Q_Vir(0)) = c,
        S_r = [t^r]H(t)/r.

    The output is the scalar-projected coefficient, not the raw
    transferred operation \(m_r\).
    """
    if r < 4:
        raise ValueError("the closed form starts at arity r=4")
    c_val = Symbol('c') if c is None else S(c)
    D = S(80) / (5 * c_val + 22)
    shape = virasoro_catalan_shape_factor(r, D / 144)
    return simplify((-6)**(r - 4) * D * shape / (2 * r * c_val**(r - 3)))


def shadow_closed_form_scope() -> Dict[str, Any]:
    r"""Scope guard for the closed-form Virasoro scalar shadow formula."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('gamma', 'epsilon'),
        'hypotheses': ('hypAmbientWtCpl', 'effScalarShadowProj'),
        'object': 'normalized scalar T-line shadow coefficient',
        'branch': 'sqrt(Q_Vir(0))=c',
        'excluded_central_charges': (S.Zero, -S(22) / 5),
        'Q_Vir': '(c+6t)^2 + 80 t^2/(5c+22)',
        'H': 't^2 sqrt(Q_Vir(t))',
        'coefficient_rule': 'S_r=[t^r]H(t)/r',
        'shape_factor': (
            'F_r(x)=sum_j (-1)^j C_j binom(r-4,2j) x^j'
        ),
        'closed_form': (
            'S_r=(-6)^(r-4) D F_r(D/144)/(2 r c^(r-3))'
        ),
        'D': '80/(5c+22)',
        'half_binomial_identity': (
            'binom(1/2,m)=(-1)^(m-1) C_(m-1)/2^(2m-1)'
        ),
        'raw_transferred_mr_formula': False,
        'full_ordered_bar_invariant_before_projection': False,
        'nonconstant_pole_divisor': 'c^(r-3)(5c+22)^floor((r-2)/2)',
        'rational_scalar_denominators_are_poles': False,
        'verified_degrees': tuple(range(4, 15)),
    }


def catalan_dynkin_field_polynomial(arity: int, x=None):
    r"""Symmetric-point field polynomial for the Catalan-Dynkin branch.

    This is the formal consequence of the Virasoro-normalized rightmost
    recursion

      phi_2=x+2,  phi_3=(x+2)(x+3).

    It is a statement about the symmetric-point field polynomial, not the full
    ordered spectral polynomial.
    """
    if arity < 2:
        raise ValueError("field polynomials start at arity 2")
    x_val = Symbol('x') if x is None else S(x)
    if arity == 2:
        return simplify(x_val + 2)
    if arity % 2 == 0:
        return S.Zero
    n = (arity - 3) // 2
    cat = factorial(2 * n) / (factorial(n) * factorial(n + 1))
    product = S.One
    for m in range(2, arity + 1):
        product *= x_val + m
    return simplify((-1)**n * cat * product)


def catalan_dynkin_parity_scope() -> Dict[str, Any]:
    r"""Scope guard for the Catalan-Dynkin parity theorem."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'gamma', 'epsilon'),
        'hypotheses': (
            'chosen one-generator branch',
            'hypAmbientWtCpl',
            'effKoszul',
            'Virasoro-normalized Dynkin string phi_2=x+2, phi_3=(x+2)(x+3)',
            'rightmost-reduction identity at the symmetric point',
        ),
        'object': 'symmetric-point field polynomial',
        'even_arity': 'phi_k=0 for even k>=4',
        'odd_arity': 'phi_(2n+3)=(-1)^n C_n prod_(m=2)^(2n+3)(x+m)',
        'root_killing_points': 'x=-j for the inner arity j',
        'full_ordered_spectral_polynomial_statement': False,
        'all_one_generator_class_M_statement': False,
        'arbitrary_conformal_weight_statement': False,
        'diagonal_reflection_quotient_only': True,
    }


def crossing_stasheff_scope() -> Dict[str, Any]:
    r"""Scope guard for Stasheff crossing versus physical bootstrap crossing."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'beta', 'gamma'),
        'hypotheses': (
            'ordered FM_3 collision chart',
            'chiral OPE-channel comparison',
            'hypAmbientWtCpl',
        ),
        'object': 'singular-OPE channel identity',
        'bar_identity': 'arity-3 Stasheff/Borcherds relation d_B^2=0',
        'requires_fourth_state_pairing': True,
        's_channel': 'left collision (i,j) first',
        't_channel': 'right collision (j,k) first',
        'contact_term': 'chain-level ternary homotopy/contact distribution',
        'virasoro_contact': 'm3=-A3 in the chosen Virasoro Stasheff gauge',
        'affine_case': 'Lie-Jacobi binary channel, no Virasoro wheel',
        'full_conformal_bootstrap_statement': False,
        'physical_crossing_requires_delta_data': True,
        'minimal_m2_nonassociativity_statement': False,
        'positivity_or_unitarity_statement': False,
    }


def shapovalov_channel_norm_squared(
    coefficients,
    gram=None,
):
    r"""Quadratic Shapovalov channel norm for finite-level coefficients.

    The bootstrap-positive object is the Gram norm of a whole finite-level
    channel vector, not a coordinate-wise bound on individual OPE
    coefficients.  This helper uses real/symmetric coordinates; it is the
    algebraic guard needed by the manuscript statement.
    """
    coeffs = tuple(S(x) for x in coefficients)
    size = len(coeffs)
    if gram is None:
        matrix = [
            [S.One if i == j else S.Zero for j in range(size)]
            for i in range(size)
        ]
    else:
        matrix = [[S(entry) for entry in row] for row in gram]
        if len(matrix) != size or any(len(row) != size for row in matrix):
            raise ValueError("Gram matrix must be square with vector size")

    return simplify(
        sum(coeffs[i] * matrix[i][j] * coeffs[j]
            for i in range(size) for j in range(size))
    )


def shapovalov_projected_channel_norm_squared(
    coefficients,
    kept_indices,
    gram=None,
):
    r"""Quadratic norm after projection to a coordinate-orthogonal channel.

    The caller is responsible for using a Shapovalov-orthogonal finite-level
    basis.  In that basis this is the finite-dimensional model for the
    quotient contraction \(\|P_N v\|^2\le \|v\|^2\).
    """
    coeffs = tuple(S(x) for x in coefficients)
    kept = tuple(kept_indices)
    if any(index < 0 or index >= len(coeffs) for index in kept):
        raise ValueError("kept_indices must refer to entries of coefficients")
    projected = [S.Zero for _ in coeffs]
    for index in kept:
        projected[index] = coeffs[index]
    return shapovalov_channel_norm_squared(projected, gram=gram)


def shapovalov_bootstrap_scope() -> Dict[str, Any]:
    r"""Scope guard for Shapovalov positivity in the bootstrap paragraph."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'delta'),
        'hypotheses': (
            'finite-level Shapovalov normalisation',
            'hypAmbientWtCpl',
            'unitary Hilbert quotient',
            'chiral OPE-channel comparison',
        ),
        'object': 'finite-level channel Gram norm',
        'channel_vector': 'v_ij^(N) in M_(c,h)[N]',
        'propagator_source': 'inverse finite-level Shapovalov Gram matrix',
        'kac_denominators': 'det G_M(c,h) for finite levels M<=k',
        'unitary_statement': 'orthogonal quotient projection is contractive',
        'coordinatewise_ope_bound_asserted': False,
        'channel_norm_contraction_asserted': True,
        'unitary_minimal_models_only': 'c=1-6/(m(m+1)), m>=2',
        'generic_minimal_model_unitarity_asserted': False,
        'raw_verma_bar_implies_bpz_constraints': False,
        'raw_hpl_summands_may_have_kac_poles': True,
        'scalar_shadow_projection_is_separate': True,
    }


def virasoro_large_c_shadow_asymptotics(max_r: int = 9) -> Dict[int, Dict[str, Any]]:
    r"""Large-\(c\) asymptotics of the scalar Virasoro shadow branch.

    For \(r\ge4\), the scalar-projected branch satisfies

        S_r = 8(-6)^(r-4)/(r c^(r-2)) + O(c^(1-r)).

    This is the scalar identity/contact lane after scalar shadow projection,
    not the full large-\(c\) conformal bootstrap.
    """
    if max_r < 4:
        raise ValueError("large-c shadow asymptotics are recorded for r >= 4")
    c_sym = Symbol('c')
    coeffs = virasoro_shadow_metric_coefficients(max_r=max_r, c=c_sym)
    out: Dict[int, Dict[str, Any]] = {}
    for r in range(4, max_r + 1):
        leading_constant = simplify(limit(coeffs[r] * c_sym**(r - 2), c_sym, oo))
        out[r] = {
            'coefficient': coeffs[r],
            'leading_constant': leading_constant,
            'expected_leading_constant': simplify(S(8) * (-6)**(r - 4) / r),
            'decay_power': 2 - r,
            'leading_term': f'8*(-6)^({r}-4)/({r} c^({r}-2))',
            'error_order': f'O(c^({1-r}))',
        }
    return out


def large_c_bootstrap_scope() -> Dict[str, Any]:
    r"""Scope guard for the scalar large-c shadow/contact proposition."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('gamma', 'epsilon', 'delta'),
        'hypotheses': (
            'hypAmbientWtCpl',
            'effScalarShadowProj',
            'identity-block contact comparison',
            'Brown-Henneaux normalisation c=3ell/(2G_N)',
        ),
        'object': 'scalar Virasoro T-line shadow branch',
        'Q_Vir': '(c+6t)^2 + 80 t^2/(5c+22)',
        'H': 't^2 sqrt(Q_Vir(t))',
        'coefficient_rule': 'S_r=[t^r]H(t)/r',
        'large_c_asymptotic': 'S_r=8*(-6)^(r-4)*c^(2-r)/r+O(c^(1-r))',
        'S4_leading': '2/c^2',
        'S5_leading': '-48/(5c^3)',
        'scalar_radius': 'c*sqrt((5c+22)/(180c+872)) ~ c/6',
        'full_large_c_bootstrap_statement': False,
        'non_vacuum_blocks_controlled': False,
        'single_valued_crossing_asserted': False,
        'positivity_or_spectral_density_asserted': False,
        'identity_contact_lane_only_after_comparison': True,
        'old_coefficient_exponent': 'not [t^r]H=O(c^(3-r))',
    }


def otoc_braiding_phase(h_p, h_v, h_w, orientation: str = 'positive'):
    r"""Diagonal Virasoro braiding phase for an OTO block channel.

    This is only the multiplicity-free diagonal specialization of the
    conformal-block monodromy matrix.  The full OTOC requires thermal
    coefficients, antiholomorphic blocks, and analytic continuation data.
    """
    sign_by_orientation = {
        'positive': S.One,
        'inverse': S.NegativeOne,
    }
    if orientation not in sign_by_orientation:
        raise ValueError("orientation must be 'positive' or 'inverse'")
    exponent = sign_by_orientation[orientation] * 2 * pi * I * (
        S(h_p) - S(h_v) - S(h_w)
    )
    return simplify(exp(exponent))


def otoc_r_matrix_scope() -> Dict[str, Any]:
    r"""Scope guard for the OTOC/R-matrix monodromy theorem."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': (
            'hypAmbientWtCpl',
            'hypVirTorusBlock',
            'conformal-block local system',
            'chosen OTO continuation path',
            'thermal trace coefficients',
            'spectral-line R-matrix comparison',
        ),
        'object': 'chiral conformal-block monodromy operator',
        'formula': 'F_OTOC^(chi,mon)=sum_pq a_p^beta rho(gamma_OTO)_p^q F_q(z_OTO)',
        'diagonal_phase': 'exp(2*pi*i*(h_p-h_V-h_W)) for positive orientation',
        'inverse_orientation': 'inverse phase',
        'full_normalized_otoc_asserted': False,
        'thermal_coefficients_computed_by_r_matrix': False,
        'anti_holomorphic_sector_included': False,
        'lyapunov_exponent_determined': False,
        'scalar_phase_without_diagonalization': False,
        'boltzmann_weight_formula_asserted': False,
        'monodromy_matrix_required_generically': True,
    }


def mss_bound_value(beta):
    r"""Maldacena--Shenker--Stanford Lyapunov bound \(2\pi/\beta\).

    The bound belongs to the normalized physical thermal OTOC satisfying
    the MSS strip hypotheses.  The annular bar complex supplies the
    periodicity/curvature datum used in the comparison, not the analytic
    boundedness hypothesis.
    """
    beta_val = S(beta)
    if beta_val == 0:
        raise ValueError("beta must be nonzero")
    return simplify(2 * pi / beta_val)


def mss_annular_bar_scope() -> Dict[str, Any]:
    r"""Scope guard for the MSS/annular-bar theorem."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': (
            'hypAmbientWtCpl',
            'hypVirTorusBlock',
            'hypModularCardy',
            'MSS normalized thermal OTOC strip hypotheses',
            'unitarity',
            'factorisation',
            'boundedness on the half-strip',
            'thermal analyticity of width beta',
            'HHLL identity-block dominance for saturation',
        ),
        'object': 'normalized physical thermal OTOC strip bound',
        'analytic_theorem': 'Maldacena-Shenker-Stanford',
        'bound': 'lambda_L <= 2*pi/beta',
        'saturation': 'lambda_L=2*pi/beta only under HHLL identity-block dominance',
        'identity_block_scale': 'exp(2*pi*t/beta)',
        'bar_complex_role': 'genus-1 annular curvature and wrap-around periodicity datum',
        'bar_curvature_alone_proves_bound': False,
        'boundedness_or_positivity_from_bar_complex': False,
        'thermal_coefficients_computed_by_bar_complex': False,
        'wraparound_mode_alone_equals_lyapunov': False,
        'full_normalized_otoc_from_annular_bar_alone': False,
        'identity_block_saturation_requires_HHLL': True,
        'identity_block_saturation_requires_vacuum_dominance': True,
    }


def scrambling_time_from_amplitude(beta, c, amplitude=1):
    r"""Scrambling scale from \(A c^{-1}\exp(2\pi t/\beta)\asymp 1\).

    The returned value is the exact solution of the normalized threshold
    \(A c^{-1}\exp(2\pi t/\beta)=1\).  The manuscript only uses this up
    to the \(O_\beta(1)\) ambiguity in the physical threshold.
    """
    beta_val = S(beta)
    c_val = S(c)
    amplitude_val = S(amplitude)
    if beta_val == 0:
        raise ValueError("beta must be nonzero")
    if amplitude_val == 0:
        raise ValueError("amplitude must be nonzero")
    return simplify(beta_val * log(c_val / amplitude_val) / (2 * pi))


def scrambling_time_scope() -> Dict[str, Any]:
    r"""Scope guard for the conditional scrambling-time proposition."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('gamma', 'epsilon', 'delta'),
        'hypotheses': (
            'hypAmbientWtCpl',
            'effScalarShadowProj',
            'hypVirTorusBlock',
            'hypModularCardy',
            'HHLL identity-block dominance',
            'nonzero O(1) probe normalization',
            'pre-scrambling expansion F=1-u+O(u^2)',
        ),
        'object': 'HHLL identity-block scrambling scale',
        'expansion_parameter': 'u(t)=A_id*c^(-1)*exp(2*pi*t/beta)',
        'threshold_equation': 'u(t_*) asymp 1',
        'scale_with_amplitude': 't_*=(beta/(2*pi))*log(c/A_id)+O_beta(1)',
        'leading_scale': 't_*=(beta/(2*pi))*log(c)+O_beta(1)',
        'exact_log_c_without_amplitude_asserted': False,
        'shadow_tower_alone_proves_scrambling': False,
        'full_thermal_otoc_from_scalar_shadow_asserted': False,
        'nth_physical_otoc_correction_from_Sr_asserted': False,
        'all_fm4_strata_equal_at_tstar_asserted': False,
        'post_tstar_interior_dominance_asserted': False,
        'fully_scrambled_from_bar_complex_alone': False,
        'finite_degree_failure_scale_on_selected_branch': True,
    }


def w3_w_line_shadow_coefficient(arity: int, c=None):
    r"""Normalized scalar \(W\)-line shadow coefficient for \(\mathcal W_3\).

    The scalar branch is

        Q_W(t) = (4 c^2 / 9) (1 + delta_3 t^2),
        delta_3 = 122880/(c^2(5c+22)^3),
        H_W(t) = t^2 sqrt(Q_W(t)),  sqrt(Q_W(0)) = 2c/3.

    Odd arities vanish by \(W\mapsto -W\).  This is the scalar \(W\)-line
    projection, not the full two-variable \(\mathcal W_3\) shadow tensor.
    """
    if arity < 2:
        raise ValueError("shadow coefficients start at arity 2")
    c_val = Symbol('c') if c is None else S(c)
    if arity % 2 == 1:
        return S.Zero
    half_arity = arity // 2
    if half_arity == 1:
        return simplify(c_val / 3)
    cat = factorial(2 * half_arity - 2) / (
        factorial(half_arity - 1) * factorial(half_arity)
    )
    return simplify(
        (-1)**half_arity
        * cat
        * S(30720)**(half_arity - 1)
        / (
            3
            * (2 * half_arity - 3)
            * c_val**(2 * half_arity - 3)
            * (5 * c_val + 22)**(3 * (half_arity - 1))
        )
    )


def w3_w_line_closed_form_scope() -> Dict[str, Any]:
    r"""Scope guard for the \(\mathcal W_3\) scalar W-line closed form."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'gamma', 'epsilon'),
        'hypotheses': ('hypAmbientWtCpl', 'effScalarShadowProj'),
        'object': 'normalized scalar W-line shadow coefficient',
        'line_branch': 'x_T=0, x_W != 0',
        'branch': 'sqrt(Q_W(0))=2c/3',
        'excluded_central_charges': (S.Zero, -S(22) / 5),
        'Q_W': '(4c^2/9)(1 + 122880 t^2/(c^2(5c+22)^3))',
        'H_W': 't^2 sqrt(Q_W(t))',
        'coefficient_rule': 'S_n^W=[t^n]H_W(t)/n',
        'S2': 'c/3',
        'odd_coefficients_vanish': True,
        'closed_form_even': (
            'S_(2r)^W=(-1)^r C_(r-1) 30720^(r-1) / '
            '(3(2r-3)c^(2r-3)(5c+22)^(3(r-1)))'
        ),
        'raw_two_variable_shadow_tensor_formula': False,
        'raw_complementarity_for_S': False,
        'pole_cleared_normalization_is_constant': True,
        'verified_degrees': tuple(range(2, 15)),
    }


def gravity_infinite_soft_shadow_hierarchy_scope() -> Dict[str, Any]:
    r"""Scope guard for the Virasoro infinite soft-shadow hierarchy.

    The infinite class-M statement is an algebraic support-shadow
    non-termination statement in the completed ambient.  Physical higher
    soft theorems require the celestial comparison datum, and each
    degree-r channel contains the normalized shadow projection together
    with lower Maurer--Cartan bracket terms.  It is not the assertion
    that the raw operation m_r alone is a physical soft operator.
    """
    coeffs = virasoro_shadow_metric_coefficients(max_r=5)
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta'),
        'hypotheses': ('hypAmbientWtCpl', 'hypCelSoft'),
        'ambient_object': 'completed class-M Virasoro support-shadow packet',
        'soft_order_to_arity': 'p maps to r=p+2 after CelSoft comparison',
        'degree_r_channel': 'normalized Sh(Theta_r) plus lower MC brackets',
        'raw_operation_alone_controls_soft_order': False,
        'higher_physical_soft_theorem_asserted_without_celsoft': False,
        'independent_primitive_classes_asserted': False,
        'nontermination_statement': (
            'no finite support cutoff for the scalar metric branch'
        ),
        'shadow_metric_branch': (
            'H(t)=t^2 c sqrt(1+12t/c+(36+80/(5c+22))t^2/c^2)'
        ),
        'S2': coeffs[2],
        'S3': coeffs[3],
        'S4': coeffs[4],
        'S5': coeffs[5],
        'S5_formula': '-48/(c^2(5c+22))',
        'quintic_bracket_subchannel': '20/(c^2(5c+22))',
        'large_c_S5': '-48/(5c^3)',
        'bh_scaling_requires': 'c=3 ell/(2 G_N)',
        'table_source_language': (
            'Theta_(p+2) projection with lower MC bracket terms, not m_(p+2) alone'
        ),
    }


def virasoro_shadow_metric_pole_profile(max_r: int = 12) -> Dict[int, Dict[str, Any]]:
    r"""Pole divisors of normalized Virasoro scalar shadow coefficients.

    This records the nonconstant denominator support in Q[c].  Rational
    scalar denominators such as 3, 7, 11 are kept separately; they are
    not central-charge poles.
    """
    if max_r < 4:
        raise ValueError("the pole-divisor formula starts at arity r=4")
    c_sym = Symbol('c')
    coeffs = virasoro_shadow_metric_coefficients(max_r=max_r, c=c_sym)
    out: Dict[int, Dict[str, Any]] = {}
    for r in range(4, max_r + 1):
        numerator, denominator = fraction(simplify(coeffs[r]))
        expected = c_sym**(r - 3) * (5 * c_sym + 22)**((r - 2) // 2)
        ratio = simplify(denominator / expected)
        out[r] = {
            'coefficient': coeffs[r],
            'numerator': numerator,
            'denominator': denominator,
            'pole_divisor': expected,
            'rational_scalar_denominator': ratio,
            'pole_divisor_matches': ratio.is_number,
        }
    return out


def stasheff_scalar_pole_cancellation_scope() -> Dict[str, Any]:
    r"""Scope guard for Stasheff cancellation on the scalar shadow branch.

    The cancellation theorem is a statement about the normalized scalar
    shadow projection.  It says the central-charge pole divisor has only
    c and 5c+22 factors.  It does not claim that raw HPL summands or the
    full transferred operation m_r have no Kac-table denominators.
    """
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('gamma', 'epsilon'),
        'hypotheses': ('hypAmbientWtCpl', 'effScalarShadowProj'),
        'object': 'normalized scalar shadow pole divisor',
        'branch': 'sqrt(Q_Vir(0))=c',
        'Q_Vir': '(c+6t)^2 + 80 t^2/(5c+22)',
        'coefficient_rule': 'S_r=[t^r]H(t)/r',
        'nonconstant_pole_divisor': 'c^(r-3)(5c+22)^floor((r-2)/2)',
        'rational_scalar_denominators_ignored': True,
        'raw_hpl_summands_may_have_kac_poles': True,
        'full_transferred_mr_pole_cancellation_asserted': False,
        'kac_table_poles_cancel_after_scalar_projection': True,
        'planarity_literal_nonplanar_tree_claim': False,
        'operator_to_scalar_projection_constructed_all_arities': False,
        'verified_degrees': tuple(range(4, 13)),
    }


def maloney_witten_scope_profile() -> Dict[str, Any]:
    r"""Scope guard for the Maloney--Witten and Page/Stokes readings.

    The universal holography functor produces the chain-level scalar
    trace and its genus expansion.  It does not produce the
    SL_2(Z)/Gamma_infty orbit sum over thermal saddles from
    Z_der^ch(Vir_c) alone, and the convergent FP scalar tower does not
    contain a Page or BTZ Stokes jump.
    """

    return {
        'bar_trace': 'Tr_Bord(Vir_c) exp(Theta)',
        'bar_trace_role': 'perturbative thermal-AdS/BTZ seed',
        'maloney_witten_sum': (
            'Z_MW(tau,taubar) = sum_{gamma in SL2Z/Gamma_inf} '
            'Z_thermal(gamma tau, gamma taubar)'
        ),
        'phi_hol_outputs_maloney_witten_sum': False,
        'bar_trace_equals_maloney_witten_sum': False,
        'requires_for_maloney_witten_sum': (
            'Brown-Henneaux dictionary',
            'saddle labelling',
            'saddle set',
            'ensemble prescription',
            'modular invariance',
        ),
        'scalar_series_variable': 'hbar^2',
        'scalar_series_radius': '4*pi^2',
        'ordinary_gevrey1_borel_transform': 'entire',
        'page_or_btz_stokes_from_scalar_tower': False,
        'page_stokes_requires': (
            'raw two-sector transseries',
            'sectorial Borel summability',
            'Stokes automorphism',
            'modular invariance',
            'vacuum dominance',
            'saddle extraction',
        ),
    }


def page_curve_profile(c, s_bh) -> Dict[str, Any]:
    r"""Two-sector Page-profile algebra in the real comparison window.

    The model is conditional on the raw gravitational Page transseries,
    the modular-Cardy entropy functional, and the same-family window
    0 < c < 26.  It records only the scalar branch crossing
    c t / 6 = S_BH - (26-c)t/6; it does not derive the transseries or
    the entropy functional from the FP scalar tower.
    """

    c_val = S(c)
    s_val = S(s_bh)

    if c_val.is_number:
        if not (bool(c_val > 0) and bool(c_val < 26)):
            raise ValueError("Page profile requires the real window 0 < c < 26")
    if s_val.is_number:
        if not bool(s_val > 0):
            raise ValueError("Page profile requires S_BH > 0")

    hawking_rate = simplify(c_val / 6)
    island_rate = simplify((26 - c_val) / 6)
    t_page = simplify(3 * s_val / 13)
    s_page = simplify(c_val * s_val / 26)
    t_evap = simplify(6 * s_val / (26 - c_val))
    crossing_balance = simplify(
        hawking_rate * t_page - (s_val - island_rate * t_page)
    )

    return {
        'claim_status': 'Conditional',
        'window': '0 < c < 26',
        'interval': '0 <= t <= t_evap',
        'hawking_rate': hawking_rate,
        'island_rate': island_rate,
        't_page': t_page,
        's_page': s_page,
        't_evap': t_evap,
        'crossing_balance': crossing_balance,
        'branch_equation': 'c*t/6 = S_BH - (26-c)*t/6',
        'profile': 'min(c*t/6, S_BH - (26-c)*t/6)',
        'scalar_tower_alone_derives_page_curve': False,
        'borel_singularity_gives_exp_page_time_formula': False,
    }


def page_curve_scope() -> Dict[str, Any]:
    r"""Scope guard for the conditional Page-profile proposition."""

    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta', 'epsilon'),
        'hypotheses': (
            'hypAmbientWtCpl',
            'effScalarShadowProj',
            'hypModularCardy',
            'raw gravitational Page transseries',
            'real same-family window 0<c<26',
            'two-sector saddle completeness',
            'entropy functional',
        ),
        'object': 'two-sector model radiation entropy profile',
        'time_interval': '0 <= t <= t_evap = 6*S_BH/(26-c)',
        'profile': 'min(c*t/6, S_BH-(26-c)*t/6)',
        't_page': '3*S_BH/13',
        's_page': 'c*S_BH/26',
        'scalar_tower_alone_derives_page_curve': False,
        'valid_for_c_ge_26_without_extra_hypothesis': False,
        'post_evaporation_profile_asserted': False,
        'borel_singularity_gives_exp_page_time_formula': False,
        'line_verdier_sector_not_naive_koszul_slogan': True,
    }


def desitter_central_charge(ell_dS, g_newton) -> Any:
    r"""Real-section de Sitter Brown--Henneaux scalar normalization."""

    ell_val = S(ell_dS)
    g_val = S(g_newton)
    if ell_val.is_number and not bool(ell_val > 0):
        raise ValueError("de Sitter radius must be positive")
    if g_val.is_number and not bool(g_val > 0):
        raise ValueError("Newton constant must be positive")
    return simplify(3 * ell_val / (2 * g_val))


def desitter_horizon_entropy_from_radius(ell_dS, g_newton) -> Any:
    r"""Gibbons--Hawking entropy of the dS3 cosmological horizon."""

    ell_val = S(ell_dS)
    g_val = S(g_newton)
    if ell_val.is_number and not bool(ell_val > 0):
        raise ValueError("de Sitter radius must be positive")
    if g_val.is_number and not bool(g_val > 0):
        raise ValueError("Newton constant must be positive")
    return simplify(pi * ell_val / (2 * g_val))


def desitter_horizon_entropy_from_c(c_dS) -> Any:
    r"""dS3 entropy written in the real scalar central-charge parameter."""

    c_val = S(c_dS)
    if c_val.is_number and not bool(c_val > 0):
        raise ValueError("de Sitter scalar central charge must be positive")
    return simplify(pi * c_val / 3)


def desitter_shadow_profile(ell_dS, g_newton) -> Dict[str, Any]:
    r"""Scope and formula data for the de Sitter scalar shadow lane."""

    c_val = desitter_central_charge(ell_dS, g_newton)
    entropy_radius = desitter_horizon_entropy_from_radius(ell_dS, g_newton)
    entropy_c = desitter_horizon_entropy_from_c(c_val)
    kappa = simplify(c_val / 2)

    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'gamma', 'epsilon'),
        'hypotheses': (
            'de Sitter real-section metric normalization',
            'hypAmbientWtCpl',
            'effScalarShadowProj',
        ),
        'c_dS': c_val,
        'kappa_dS': kappa,
        'entropy_from_radius': entropy_radius,
        'entropy_from_c': entropy_c,
        'entropy_balance': simplify(entropy_radius - entropy_c),
        'horizon_length': simplify(2 * pi * S(ell_dS)),
        'genus_coefficient': (
            'F_g^{dS,sc}(Vir_{c_dS}) = '
            'kappaChHodge(Vir_{c_dS}) lambda_g^FP'
        ),
        'fixed_point_c': S(13),
        'fixed_point_kappa': Rational(13, 2),
        'fixed_point_entropy': simplify(13 * pi / 3),
        'nariai_geometry_constructed': False,
        'desitter_hilbert_space_constructed': False,
        'dscft_correlator_functor_constructed': False,
        'banks_dimension_theorem': False,
        'literal_complex_wick_rotation_gives_real_coefficients': False,
    }


def jt_wp_spectral_curve_y(z) -> Any:
    r"""Compact Weil--Petersson normalization of the JT spectral curve."""

    z_val = S(z)
    return simplify(sin(2 * pi * z_val) / (4 * pi))


def jt_disk_density(E) -> Any:
    r"""Physical JT disk density after the energy-contour continuation."""

    e_val = S(E)
    if e_val.is_number and not bool(e_val >= 0):
        raise ValueError("JT disk energy must be nonnegative")
    return simplify(sinh(2 * pi * sqrt(e_val)) / (4 * pi**2))


def jt_wp_to_density_balance(E) -> Any:
    r"""Check y_WP(i sqrt(E))/(i pi) against the physical disk density."""

    e_val = S(E)
    continued = jt_wp_spectral_curve_y(I * sqrt(e_val)).rewrite(sinh)
    density_from_curve = simplify(continued / (I * pi))
    return simplify(density_from_curve - jt_disk_density(e_val))


def jt_schwarzian_scope() -> Dict[str, Any]:
    r"""Scope guard for the conditional Schwarzian/JT scalar comparison."""

    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma', 'delta', 'epsilon'),
        'hypotheses': (
            'Schwarzian comparison datum',
            'hypAmbientWtCpl',
            'effScalarShadowProj',
            'JT spectral-curve normalization',
            'Eynard-Orantin recursion kernel',
            'Weil-Petersson volume comparison',
            'energy-contour convention',
            'matrix-integral or Stokes completion datum',
        ),
        'compact_curve': 'x=z^2, y_WP(z)=sin(2*pi*z)/(4*pi)',
        'physical_density': 'rho_0(E)=sinh(2*pi*sqrt(E))/(4*pi^2)',
        'branch_of_sqrt_x_required': True,
        'physical_density_requires_contour': True,
        'scalar_shadow_alone_proves_wp_volumes': False,
        'scalar_shadow_alone_completes_jt_series': False,
        'physical_jt_density_from_sine_curve_without_contour': False,
        'bernoulli_decay_is_jt_borel_completion': False,
        'schur_q_series_supplies_uv_completion': False,
        'miki_stress_sheaf_derives_jt_measure': False,
    }


def colored_partition_number(colors: int, n: int) -> int:
    r"""Coefficient of p^n in product_{m>=1}(1-p^m)^(-colors)."""

    if colors < 0:
        raise ValueError("number of colors must be nonnegative")
    if n < 0:
        raise ValueError("partition degree must be nonnegative")

    coeffs = [0] * (n + 1)
    coeffs[0] = 1
    for part in range(1, n + 1):
        for _ in range(colors):
            for degree in range(part, n + 1):
                coeffs[degree] += coeffs[degree - part]
    return coeffs[n]


def k3_borcherds_scalar_bridge_scope() -> Dict[str, Any]:
    r"""Scope and normalization data for the K3 x E Borcherds scalar bridge."""

    d5_factor = Rational(1, 64)
    phi10_op_factor = d5_factor**2
    return {
        'claim_status': 'ProvedElsewhere',
        'licensing_tags': ('beta', 'gamma', 'epsilon'),
        'hypotheses': (
            'Gritsenko-Nikulin half-K3 Jacobi normalization',
            'Borcherds singular-theta lift',
            'DMVV second-quantized K3 elliptic-genus scalar',
            'Gottsche Hilbert-scheme specialization',
            'Bruinier Heegner projection in the half-K3 convention',
            'hypAmbientWtCpl',
            'effScalarShadowProj',
        ),
        'phi10_un': 'Delta_5^2',
        'd5_factor': d5_factor,
        'phi10_op_factor': phi10_op_factor,
        'bps_scalar_prefactor_delta5_inverse_square': S(1),
        'op_dt_scalar_prefactor_delta5_inverse_square': S(-4096),
        'p24_5': colored_partition_number(24, 5),
        'hilb5_k3_euler': S(176256),
        'bruinier_reduced_c2_triangle': S(0),
        'bruinier_reduced_c3_half_k3': S(-64),
        'phi_minus2_1_reduced_c3_other_convention': S(-8),
        'uses_phi_minus2_1_reduced_convention': False,
        'p24_5_is_bruinier_obstruction_coefficient': False,
        'scalar_identity_promotes_to_gravity_line_trace': False,
        'scalar_identity_proves_maloney_witten_equality': False,
        'one_variable_hilbert_series_is_full_siegel_form': False,
        'op_normalization_changes_bkm_denominator_algebra': False,
    }


def k3_scalar_no_promotion_scope() -> Dict[str, Any]:
    r"""No-promotion data for the K3 x E scalar shadow.

    The scalar coordinates determine the Borcherds shadow, not a filtered
    SC^{ch,top} morphism to the gravity-line boundary algebra.
    """

    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('beta', 'gamma', 'epsilon'),
        'ambient': 'hypAmbientWtCpl',
        'effectiveness': 'effScalarShadowProj',
        'uses_scalar_non_faithfulness_lemma': True,
        'scalar_character_functor_faithful': False,
        'acyclic_filtered_summand_preserves_scalar_character': True,
        'acyclic_filtered_summand_preserves_chain_object': False,
        'scalar_data': {
            'phi10_un': 'Delta_5^2',
            'z_bps': '(Phi_10^{un})^{-1}',
            'kappa_bkm': S(5),
            'kappa_tuple_fourfold': (S(0), S(3), S(5), S(24)),
        },
        'missing_beta_comparison_data': (
            'positive-half Hall-Borcherds E1-chiral bialgebra morphism',
            'completed Drinfeld-double extension through the current envelope on E',
            'filtered SC^{ch,top} morphism to the Virasoro gravity-line boundary algebra',
            'derived-centre trace compatibility with the BPS scalar character',
        ),
        'determines_filtered_sc_chtop_morphism': False,
        'determines_ordered_virasoro_bar_trace': False,
        'determines_maloney_witten_sum': False,
        'scalar_shadow_is_object_equivalence': False,
    }


def heptagon_growth_bound_constant(generator_rank: int, ope_norm=1):
    r"""Constant in the finite-type heptagon Gevrey-1 bound.

    The bound depends on the finite strong-generator rank R.  W fixes the
    weight window; R counts labels inside that window.
    """

    if generator_rank < 1:
        raise ValueError("generator rank must be positive")
    norm = S(ope_norm)
    if norm < 0:
        raise ValueError("OPE norm must be nonnegative")
    norm_factor = max(S.One, norm)
    return 4 * S(generator_rank) * norm_factor * exp(pi * sqrt(S(2) / 3))


def heptagon_growth_bound_scope() -> Dict[str, Any]:
    r"""Scope for the heptagon ordered-bar Gevrey-1 coefficient bound."""

    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('gamma', 'epsilon'),
        'hypotheses': (
            'hypAmbientWtCpl',
            'effKoszul',
            'finite strong-generator profile (W,R)',
            'Lyndon-PBW basis',
            'finite-type OPE norm M',
        ),
        'constant_parameters': ('W', 'R', 'M'),
        'constant_formula': '4*R*max(1,M)*exp(pi*sqrt(2/3))',
        'depends_only_on_W_M': False,
        'requires_generator_rank': True,
        'tree_shape_bound': 'Catalan_n <= 4^n',
        'weight_decomposition_bound': 'p(n) <= exp(pi*sqrt(2*n/3))',
        'partition_absorbed_exponentially': True,
        'pbw_symmetrization_bound': 'n!',
        'borel_transform_local_radius': '>= 1/C(W,R,M)',
        'proves_sectorial_continuation': False,
        'proves_borel_singularity_location': False,
        'proves_maloney_witten_interpretation': False,
    }


def factorization_l0_gevrey_constant(
    bar_generator_rank: int,
    bv_generator_rank: int,
    ope_norm=1,
    propagator_norm=1,
):
    r"""Constant for the finite-profile factorisation-BV Gevrey-1 bound."""

    if bar_generator_rank < 1:
        raise ValueError("bar generator rank must be positive")
    if bv_generator_rank < 1:
        raise ValueError("BV generator rank must be positive")
    m_val = S(ope_norm)
    p_val = S(propagator_norm)
    if m_val < 0:
        raise ValueError("OPE norm must be nonnegative")
    if p_val < 0:
        raise ValueError("propagator norm must be nonnegative")
    norm_factor = max(S.One, m_val, p_val)
    return (
        8
        * S(bar_generator_rank)
        * S(bv_generator_rank)
        * norm_factor
        * exp(pi * sqrt(S(2) / 3))
    )


def factorization_l0_gevrey_scope() -> Dict[str, Any]:
    r"""Scope for the factorisation-BV L0 Gevrey-1 coefficient bound."""

    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('gamma', 'epsilon'),
        'hypotheses': (
            'hypAmbientWtCpl',
            'effKoszul',
            'finite strong-generator profile (W,R)',
            'finite factorisation-BV graph profile (U,D)',
            'finite OPE norm M',
            'finite regularised propagator/vertex norm P',
        ),
        'constant_parameters': ('W', 'U', 'R', 'D', 'M', 'P'),
        'constant_formula': '8*R*D*max(1,M,P)*exp(pi*sqrt(2/3))',
        'same_as_heptagon_constant': False,
        'splitting_bound': 'n+1 <= 2^n',
        'bar_tree_shape_bound': '4^j',
        'bv_graph_skeleton_bound': '4^i',
        'factorial_absorption': 'i!*j! <= n!',
        'borel_transform_local_radius': (
            '>= 1/C_fact(W,U,R,D,M,P)'
        ),
        'proves_sectorial_continuation': False,
        'proves_borel_singularity_location': False,
        'proves_bcov_instanton': False,
        'proves_maloney_witten_interpretation': False,
    }


def factorization_bcov_candidate_location():
    r"""Candidate Borel location for the conditional BCOV saddle datum."""

    return S.One / (2 * pi)


def factorization_bcov_candidate_residue(c=None):
    r"""Candidate one-loop residue c(c-25)/24 for the BCOV saddle datum."""

    c_val = Symbol('c') if c is None else S(c)
    return simplify(c_val * (c_val - 25) / 24)


def factorization_bcov_instanton_scope() -> Dict[str, Any]:
    r"""Scope guard for the conditional factorisation-BV BCOV saddle datum.

    The Borel location is the renormalised action difference in the chosen
    Borel coordinate.  The Virasoro beta coefficient c(c-25)/24 is the
    candidate one-loop residue, not the action.
    """

    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'gamma', 'delta'),
        'hypotheses': (
            'Borel-coordinate normalisation',
            'hypAmbientWtCpl',
            'hypStokes',
            'nondegenerate factorisation-BV saddle',
            'saddle-to-ordered-bar comparison',
            'oriented one-loop determinant normalisation',
        ),
        'action_difference': 'A_BV(Phi_*)=S_eff(Phi_*)-S_eff(Phi_0)',
        'candidate_location': '1/(2*pi)',
        'candidate_residue': 'c*(c-25)/24',
        'action_location_and_residue_are_separate': True,
        'action_equals_beta_residue': False,
        'old_action_formula_asserted': False,
        'local_qme_proves_singularity': False,
        'gevrey_bound_proves_singularity': False,
        'requires_sectorial_comparison': True,
        'requires_one_loop_determinant': True,
        'proves_maloney_witten_interpretation': False,
    }


def heisenberg_zwegers_shadow_scope() -> Dict[str, Any]:
    r"""Scope for the rank-one Heisenberg inverse-eta shadow theorem.

    The proved statement concerns the standard Fock vacuum trace of the
    rank-one Heisenberg algebra at nonzero level.  A determinant-line power
    eta^(-r) is a separate scalar shadow and needs its own branch and
    multiplier datum before a modular-weight assertion is meaningful.
    """

    return {
        'claim_status': 'ProvedHere',
        'licensing_tags': ('gamma', 'epsilon'),
        'hypotheses': (
            'hypAmbientWtCpl',
            'effScalarShadowProj',
            'rank-one standard Fock vacuum trace',
            'nonzero Heisenberg level',
            'inverse Dedekind eta multiplier',
        ),
        'raw_character': 'eta(tau)^(-1)',
        'coefficient_sequence': 'p(n)',
        'modular_weight': Rational(-1, 2),
        'zwegers_shadow_vanishes': True,
        'mock_completion_required': False,
        'borel_transform': 'entire of exponential type zero',
        'finite_borel_singularities': False,
        'level_changes_oscillator_count': False,
        'unbranched_real_eta_power_modular_asserted': False,
        'determinant_line_power_requires_branch_multiplier': True,
        'determinant_power_weight_formula': '-r/2 after branch/multiplier datum',
        'proves_maloney_witten_path_integral': False,
    }


def partition_number(n: int) -> int:
    r"""Ordinary partition number p(n), with p(0)=1."""

    if n < 0:
        raise ValueError("partition degree must be nonnegative")
    return colored_partition_number(1, n)


def virasoro_vacuum_verma_coefficient(n: int) -> int:
    r"""Weight-n coefficient of the generic Virasoro vacuum Verma module.

    The vacuum relation L_{-1}1=0 removes partitions with a part 1.
    Hence the coefficient of prod_{m>=2}(1-q^m)^(-1) is p(n)-p(n-1).
    """

    if n < 0:
        raise ValueError("Virasoro vacuum degree must be nonnegative")
    previous = 0 if n == 0 else partition_number(n - 1)
    return partition_number(n) - previous


def virasoro_vacuum_verma_asymptotic(n: int):
    r"""Hardy--Ramanujan leading term for p(n)-p(n-1)."""

    if n <= 0:
        raise ValueError("asymptotic degree must be positive")
    n_val = S(n)
    return (
        pi / (12 * sqrt(2))
        * n_val ** Rational(-3, 2)
        * exp(pi * sqrt(2 * n_val / 3))
    )


def virasoro_hardy_ramanujan_cardy_scope() -> Dict[str, Any]:
    r"""Scope for the Virasoro vacuum coefficient and Cardy comparison.

    The PBW/Hardy--Ramanujan statement is a coefficient theorem for the
    universal vacuum module.  Cardy growth is a separate theorem for a
    modular-invariant physical CFT under the modular/vacuum/saddle package.
    """

    return {
        'claim_status_vacuum_coefficients': 'ProvedHere',
        'claim_status_cardy_density': 'ProvedElsewhere',
        'licensing_tags': ('gamma', 'delta'),
        'hypotheses': (
            'hypAmbientWtCpl',
            'generic c>1 vacuum Virasoro module',
            'only vacuum null L_{-1}1=0',
            'hypModularCardy for physical CFT comparison',
        ),
        'vacuum_character': 'prod_{m>=2}(1-q^m)^(-1)',
        'coefficient_formula': 'p(n)-p(n-1)',
        'asymptotic_prefactor': 'pi/(12*sqrt(2))',
        'asymptotic_exponential': 'exp(pi*sqrt(2*n/3))',
        'independent_of_c': True,
        'cardy_density_formula': (
            'log rho_phys(Delta) ~ 2*pi*sqrt(c_eff*Delta/6)'
        ),
        'cardy_is_verma_coefficient_identity': False,
        'bare_verma_module_has_modular_invariance': False,
        'vacuum_pqw_equals_physical_density': False,
        'zwegers_shadow_from_verma_coefficients_alone': False,
    }


# =========================================================================
# 5. QUARTIC CONTACT INVARIANT
# =========================================================================

def quartic_contact_virasoro(c=None):
    r"""Quartic contact invariant Q^contact_Vir = 10/(c(5c+22)).

    This is the first nonlinear modular shadow coefficient beyond kappa.
    It measures the quartic resonance obstruction in the shadow obstruction
    tower for the Virasoro algebra.

    Poles: c = 0 and c = -22/5.
    Positive for c > 0.
    The quartic contact invariant controls the genus-2 shadow correction.

    Ground truth: Vol I, thm:nms-virasoro-quartic-explicit.
    """
    c_val = Symbol('c') if c is None else S(c)
    return 10 / (c_val * (5 * c_val + 22))


def quartic_contact_virasoro_exact(c_num, c_den=1):
    """Exact rational quartic contact invariant.

    Uses fractions.Fraction for exact arithmetic.
    """
    c_val = Fraction(c_num, c_den)
    return Fraction(10, 1) / (c_val * (5 * c_val + 22))


# =========================================================================
# 5. GRAVITATIONAL KAPPA
# =========================================================================

def gravity_kappa(c=None):
    r"""Modular characteristic kappa(Vir_c) = c/2.

    This is the genus-1 curvature: d_B^2 = kappa(A) * omega_1.
    For the Virasoro algebra, kappa = c/2.

    Ground truth: Vol I Theorem D, landscape_census.tex.
    """
    c_val = Symbol('c') if c is None else S(c)
    return c_val / 2


def gravity_kappa_exact(c_num, c_den=1):
    """Exact rational kappa using fractions.Fraction."""
    return Fraction(c_num, c_den) / 2


# =========================================================================
# 6. KOSZUL DUAL CENTRAL CHARGE
# =========================================================================

def koszul_dual_central_charge(c=None):
    r"""Same-family Virasoro comparison central charge: c_dual = 26 - c.

    The line-side homotopy dual of Vir_c is represented by Vir_{26-c}
    only on the strict same-family comparison surface; this function
    computes that representative's central charge.

    CRITICAL: Self-dual at c = 13, NOT c = 26.
    This is the Feigin-Frenkel involution for Virasoro.
    """
    c_val = Symbol('c') if c is None else S(c)
    return 26 - c_val


# =========================================================================
# 7. R-MATRIX POLE STRUCTURE
# =========================================================================

def gravity_laplace_kernel_poles(c=None):
    r"""Pre-dlog Laplace/OPE kernel for the Virasoro algebra.

    This is the Laplace transform of the lambda-bracket, equivalently
    the OPE generating kernel:

      r^L(z) = (c/2)/z^4 + 2T/z^2 + dT/z

    It is not the bar collision r-matrix.

    Returns a dict {pole_order: residue_description}.
    """
    c_val = Symbol('c') if c is None else S(c)

    return {
        4: {'residue': c_val / 2, 'description': 'c/2 (quartic OPE coefficient)'},
        2: {'residue': '2T', 'description': '2T (energy-momentum)'},
        1: {'residue': 'dT', 'description': 'dT (derivative)'},
    }


def gravity_r_matrix_poles(c=None):
    r"""Bar collision r-matrix pole structure for the Virasoro algebra.

    The r-matrix r(z) = Res^coll_{0,2}(Theta_A) from the binary
    genus-0 shadow of the universal MC element. The d-log kernel in the
    bar construction absorbs one pole of the pre-dlog Laplace/OPE
    kernel:

      r^L(z) = (c/2)/z^4 + 2T/z^2 + dT/z
      r^{coll}(z) = (c/2)/z^3 + 2T/z

    The derivative term dT becomes the regular z^0 part and is not a
    singular collision pole.

    Returns a dict {pole_order: residue_description}.
    """
    c_val = Symbol('c') if c is None else S(c)

    return {
        3: {'residue': c_val / 2, 'description': 'c/2 (central collision residue)'},
        1: {'residue': '2T', 'description': '2T (energy-momentum collision residue)'},
    }


def gravity_r_matrix_leading_residue(c=None):
    """Leading residue of the collision r-matrix: c/2 at the order-3 pole."""
    c_val = Symbol('c') if c is None else S(c)
    return c_val / 2


def virasoro_collision_cybe_scope():
    r"""Scope guard for the Virasoro collision-residue CYBE.

    The gravitational CYBE is the arity-three MC/Arnold cancellation after
    collision-residue extraction.  It is not the ordinary finite-dimensional
    Casimir identity for k*Omega/z, and it is not the strict quantum YBE for
    the Ponsot-Teschner fusion kernel before the associator comparison is fixed.
    """
    return {
        'kernel': 'r_coll_Vir(z) = (c/2)/z^3 + 2T/z',
        'relation': 'classical Yang-Baxter relation in completed line-side endomorphisms',
        'pre_comparison_form': 'arity-three MC equation for Theta_Vir_c',
        'proof_mechanism': 'Arnold relation on FM_3(C) after collision-residue comparison',
        'licensing_tags': ('alpha', 'beta', 'gamma'),
        'requires': ('BHdict', 'completed/pro ambient', 'collision-residue comparison'),
        'is_finite_dimensional_casimir_cybe': False,
        'is_strict_quantum_ybe_for_fusion_kernel': False,
    }


def virasoro_ds_hpl_transfer_scope() -> Dict[str, Any]:
    r"""Scope of the DS-HPL Virasoro spectral transfer theorem.

    The HPL theorem transfers products, morphisms, and the collision
    target twist to a homotopy-coherent Virasoro spectral package.  It
    does not by itself prove the strict or dg-shifted Yangian axioms.
    """
    return {
        'claim_status': 'Conditional',
        'constructed_object': 'homotopy-coherent Virasoro spectral package',
        'package_name': 'Y_Vir^HPL',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'delta'),
        'requires': (
            'principal DS chart at k != -2',
            'source affine dg-shifted Yangian comparison',
            'hypAmbientWtCpl',
            'hypKZSDR',
            'linear DS SDR',
            'complete finite-weight filtration',
        ),
        'proved_core': (
            'HPL-transferred A_infinity products',
            'HPL-transferred A_infinity morphism Delta_z',
            'finite-weight convergence of the transfer series',
            'homotopy-coherent A_infinity Yang-Baxter relation',
            'collision-residue target twist r_coll_Vir(z)',
        ),
        'collision_kernel': 'r_coll_Vir(z) = (c/2)/z^3 + 2T/z',
        'not_proved_here': (
            'strict dg-shifted Yangian presentation of H_Vir',
            'closed sign normalisation for transferred products in degree >= 3',
            'translation compatibility of the transferred spectral parameter',
            'rational strictification of the collision-residue r-matrix',
        ),
        'is_strict_or_dg_shifted_yangian_assertion': False,
        'promotes_after_open_axioms': True,
    }


def ds_ordered_bar_intertwining_scope() -> Dict[str, Any]:
    r"""Scope of the filtered principal DS comparison for the ordered bar."""
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('alpha', 'beta', 'gamma', 'delta'),
        'requires': (
            'principal DS chart at non-critical level',
            'hypAmbientWtCpl',
            'hypKZSDR',
            'finite-weight BRST concentration',
            'ordered BRST-bar bicomplex',
            'HPL-transferred R-descent datum',
            'strongly convergent PBW/ghost spectral sequence',
        ),
        'object': 'filtered ordered-bar comparison in the homotopy category',
        'comparison': (
            'B^ord(W_k(g)) -> H^0_DS(B^ord(V_k(g))) filtered quasi-isomorphism'
        ),
        'base_space': 'Conf^<_bullet(R) x FM_bullet(C)',
        'descent_square': 'commutes up to transferred R-descent homotopy',
        'not_asserted': (
            'unconditional ordered DS functor',
            'strict equality of ordered bar dg-coalgebras',
            'R-descent unchanged by DS reduction',
            'strict dg-shifted Yangian presentation',
            'physical gravitational Yangian identification',
            'ordered classification before R-descent datum',
        ),
        'h0_requires_concentration': True,
        'r_matrix_transferred_not_copied': True,
        'yangian_promotion_separate': True,
    }


def principal_ds_coproduct_primitivity_scope() -> Dict[str, Any]:
    r"""Scope of the principal-DS coproduct primitivity theorem.

    The degree-two Virasoro computation proves generator primitivity for
    T on the ghost-zero projection surface.  The all-degree principal
    W statement is conditional on the signed ghost-defect lemma for
    every nonlinear HPL morphism tree.
    """
    return {
        'claim_status': 'Conditional',
        'base_case': 'Delta_{z,2}^Vir(T,T) = 0',
        'linear_generator_coproduct': (
            'Delta_z^Vir(T) = tau_z(T) tensor 1 + 1 tensor T'
        ),
        'all_degree_statement': 'generator-level principal DS primitivity',
        'licensing_tags': ('alpha', 'gamma', 'delta'),
        'requires': (
            'principal DS chart at non-critical level',
            'ghost-zero BRST projection',
            'hypAmbientWtCpl',
            'hypKZSDR',
            'signed ghost-defect lemma for HPL morphism trees',
        ),
        'proved_here': (
            'degree-two Virasoro source-tree ghost defect',
            'degree-two Virasoro target-tree ghost defect',
            'projection kills tensor ghost bidegrees not equal to (0,0)',
        ),
        'not_asserted': (
            'ordinary Hopf primitivity for arbitrary composites',
            'unconditional all-principal-W primitivity without the signed tree count',
            'the slogan gh(h)=-1 as a complete proof',
        ),
        'degree_two_source_tree_ghost_shift': -1,
        'degree_two_target_tree_ghost_shift': -1,
        'projection_accepts_only_bidegree': (0, 0),
    }


def virasoro_primary_ward_connection(c=None, h=None, z=None):
    r"""Scalar highest-weight projection of the Virasoro collision connection.

    The singular collision residue

      r^{coll}(z) = (c/2)/z^3 + 2T/z

    acts on a primary highest-weight line by T -> h. The resulting
    scalar Ward connection one-form has coefficient

      2h/z + c/(2z^3).
    """
    c_val = Symbol('c') if c is None else S(c)
    h_val = Symbol('h') if h is None else S(h)
    z_val = Symbol('z') if z is None else S(z)
    return 2 * h_val / z_val + c_val / (2 * z_val**3)


def virasoro_primary_ward_log(c=None, h=None, z=None):
    r"""Logarithm of the scalar primary-state Ward factor.

    This is a primitive of ``virasoro_primary_ward_connection``:

      2h log(z) - c/(4z^2).
    """
    c_val = Symbol('c') if c is None else S(c)
    h_val = Symbol('h') if h is None else S(h)
    z_val = Symbol('z') if z is None else S(z)
    return 2 * h_val * log(z_val) - c_val / (4 * z_val**2)


def virasoro_primary_ward_even_coefficients(c=None, max_k=3):
    r"""Coefficients of exp(-c/(4z^2)) as a series in z^{-2}.

    Returns {k: (-c/4)^k/k!} for

      exp(-c/(4z^2)) = sum_{k>=0} coeff[k] z^{-2k}.
    """
    if max_k < 0:
        raise ValueError(f"max_k must be nonnegative, got {max_k}")

    c_val = Symbol('c') if c is None else S(c)
    return {
        k: expand((-c_val / 4)**k / factorial(k))
        for k in range(max_k + 1)
    }


def virasoro_primary_ward_even_exponents(max_k=3):
    """Exponents appearing after the z^(2h) factor is removed."""
    if max_k < 0:
        raise ValueError(f"max_k must be nonnegative, got {max_k}")
    return [-2 * k for k in range(max_k + 1)]


# =========================================================================
# 8. COMPLEMENTARITY CONSTANT
# =========================================================================

def complementarity_constant_virasoro():
    r"""Complementarity constant for Virasoro: kappa + kappa' = 13.

    kappa(Vir_c) + kappa(Vir_{26-c}) = c/2 + (26-c)/2 = 13.

    This is the Virasoro instance of Theorem C (complementarity).
    The constant 13 is twice kappa(Vir_13), the comparison-fixed value.

    For comparison:
      - Heisenberg: kappa + kappa' = 0
      - Affine KM: kappa + kappa' = 0
      - W_3: kappa + kappa' = 250/3
    """
    return S(13)


def verify_complementarity(c=None):
    """Verify kappa(Vir_c) + kappa(Vir_{26-c}) = 13."""
    c_val = Symbol('c') if c is None else S(c)
    kappa = gravity_kappa(c_val)
    kappa_dual = gravity_kappa(26 - c_val)
    total = simplify(kappa + kappa_dual)
    return {
        'kappa': kappa,
        'kappa_dual': kappa_dual,
        'sum': total,
        'equals_13': simplify(total - 13) == 0,
    }


def verdier_modular_s_shadow_scope() -> Dict[str, Any]:
    r"""Scope of the formal Verdier--S square at genus one.

    The square is a statement in the completed algebraic shadow
    complex. Verdier changes the Virasoro parameter c -> 26-c; the
    shadow S-pullback changes tau -> -1/tau and acts on the Hodge
    coefficient line. The physical torus partition function requires
    the Dedekind multiplier, E2 anomaly/completion, modular invariance,
    and vacuum dominance separately.
    """
    return {
        'claim_status': 'Conditional',
        'licensing_tags': ('beta', 'gamma'),
        'ambient': 'completed genus-one shadow complex (hypAmbientWtCpl)',
        'verdier_action': 'c -> 26 - c',
        'shadow_s_action': 'tau -> -1/tau on Hodge/multiplier coefficient line',
        'degree_zero_input': 'Theta^(1)_0(Vir_c) = (c/2) omega_1',
        'degree_zero_output': '((26-c)/2) S_sh(omega_1)',
        'commuting_square': True,
        'physical_torus_partition_function_statement': False,
        'dedekind_multiplier_required_for_partition_function': True,
        'e2_anomaly_required_for_partition_function': True,
        'cardy_expression_scope': (
            'requires ModInv, VacDom, Cardy/BTZ comparison, 0 <= c <= 26, Delta >= 0'
        ),
        'paired_cardy_expression': (
            '2*pi*sqrt(c*Delta/6) + 2*pi*sqrt((26-c)*Delta/6)'
        ),
        'not_asserted': (
            'physical torus modular invariance',
            'vacuum dominance',
            'BTZ entropy theorem',
            'real Cardy expression for c outside [0,26]',
            'Dedekind eta multiplier is trivial',
            'E2 is a modular form without completion',
        ),
    }


# =========================================================================
# 9. GENUS GENERATING FUNCTION (A-HAT GENUS)
# =========================================================================

def _lambda_fp(g):
    r"""Faber-Pandharipande intersection number lambda_g^FP.

    lambda_g^FP = (2^{2g-1} - 1) / 2^{2g-1} * |B_{2g}| / (2g)!

    where B_{2g} is the 2g-th Bernoulli number.

    These are the coefficients in:
      A-hat(x) = (x/2)/sinh(x/2) = sum_{g>=0} (-1)^g lambda_g^FP x^{2g}
    with lambda_0^FP := 1.
    """
    if g < 1:
        raise ValueError(f"Genus must be >= 1, got {g}")
    B_2g = bernoulli(2 * g)
    numerator = (2**(2*g - 1) - 1) * Abs(B_2g)
    denominator = 2**(2*g - 1) * factorial(2 * g)
    return Rational(numerator, denominator)


def _ahat_coefficient(g):
    r"""Coefficient of x^{2g} in A-hat(x) = (x/2)/sinh(x/2).

    Returns (-1)^g * lambda_g^FP for g >= 1, and 1 for g = 0.
    """
    if g == 0:
        return Rational(1)
    return (-1)**g * _lambda_fp(g)


def genus_generating_function_coefficients(c=None, max_genus=5):
    r"""Genus expansion free energy coefficients F_g for Virasoro.

    The genus-g scalar coefficient from the Wick-rotated A-hat formula:
      F_g = kappa_eff * lambda_g^FP

    where kappa_eff = (c - 26)/2 is the EFFECTIVE curvature that
    accounts for the Koszul dual shift. The physical free energy
    uses kappa_eff rather than kappa = c/2 because the partition
    function involves the relative curvature.

    The scalar free energy at genus g is:
      F_g = kappa * lambda_g^FP
    where kappa = c/2 and lambda_g^FP > 0.

    So F_g > 0 iff kappa > 0 iff c > 0.

    The EFFECTIVE (shifted) version uses kappa_eff = (c-26)/2.

    We return BOTH: the raw F_g = kappa * lambda_g^FP and the
    shifted F_g^eff = kappa_eff * lambda_g^FP.

    Explicit values:
      lambda_1 = 1/24
      lambda_2 = 7/5760
      lambda_3 = 31/967680
      lambda_4 = 127/154828800

    Returns
    -------
    dict with keys 'raw' (using kappa=c/2) and 'effective' (using
    kappa_eff=(c-26)/2), each mapping genus g to F_g.
    """
    c_val = Symbol('c') if c is None else S(c)
    kappa = c_val / 2
    kappa_eff = (c_val - 26) / 2

    raw = {}
    effective = {}
    lambda_fp_values = {}

    for g in range(1, max_genus + 1):
        lfp = _lambda_fp(g)
        lambda_fp_values[g] = lfp
        raw[g] = expand(kappa * lfp)
        effective[g] = expand(kappa_eff * lfp)

    return {
        'raw': raw,
        'effective': effective,
        'lambda_fp': lambda_fp_values,
        'kappa': kappa,
        'kappa_eff': kappa_eff,
    }


def ahat_series_coefficients(max_genus=10):
    r"""Return coefficients of x^{2g} in A-hat(x) = (x/2)/sinh(x/2).

    Standard expansion:
      A-hat(x) = 1 - x^2/24 + 7x^4/5760 - 31x^6/967680 + ...

    The sign pattern is (-1)^g for the genus-g coefficient.
    """
    result = {0: Rational(1)}
    for g in range(1, max_genus + 1):
        result[g] = _ahat_coefficient(g)
    return result


def verify_ahat_series(max_genus=6):
    r"""Verify A-hat coefficients against direct series expansion of (x/2)/sinh(x/2)."""
    x = Symbol('x')
    s = series(x / 2 / sinh(x / 2), x, 0, 2 * max_genus + 2)

    results = {}
    all_match = True

    for g in range(max_genus + 1):
        coeff_from_series = Rational(s.coeff(x, 2 * g))
        coeff_from_formula = _ahat_coefficient(g)
        match = simplify(coeff_from_series - coeff_from_formula) == 0
        all_match = all_match and match
        results[g] = {
            'series': coeff_from_series,
            'formula': coeff_from_formula,
            'match': match,
        }

    return {
        'all_match': all_match,
        'by_genus': results,
    }


# =========================================================================
# 10. SHADOW DEPTH AND CLASSIFICATION
# =========================================================================

def virasoro_shadow_class():
    """Shadow depth class for Virasoro: M (mixed, r_max = infinity).

    The Virasoro algebra has infinite shadow depth because the quintic
    obstruction o^(5)_Vir != 0 (the shadow obstruction tower never terminates).

    Contrast with:
      Heisenberg: G (Gaussian, r_max = 2)
      Affine KM:  L (Lie/tree, r_max = 3)
      Beta-gamma: C (contact/quartic, r_max = 4)
    """
    return 'M'


def virasoro_shadow_depth():
    """Shadow depth r_max for Virasoro: infinity."""
    return float('inf')


# =========================================================================
# 11. GENUS-1 HESSIAN CORRECTION
# =========================================================================

def genus1_hessian_correction_virasoro(c=None):
    r"""Genus-1 Hessian correction delta_H^(1)_Vir.

    delta_H^(1)_Vir = 120 / [c^2 (5c + 22)] * x^2

    This is the genus-1 correction to the shadow obstruction tower from the
    quartic contact invariant. It lives in the weight-4 sector.

    Returns the coefficient (without the x^2 factor).
    """
    c_val = Symbol('c') if c is None else S(c)
    return 120 / (c_val**2 * (5 * c_val + 22))


# =========================================================================
# 12. SUMMARY / DIAGNOSTIC
# =========================================================================

def gravity_diagnostic(c_val):
    """Comprehensive diagnostic for 3d gravity at central charge c.

    Returns a dict summarizing all gravitational invariants.
    """
    c = S(c_val)
    c_dual = koszul_dual_central_charge(c)
    kappa = gravity_kappa(c)
    kappa_dual = gravity_kappa(c_dual)
    q_contact = quartic_contact_virasoro(c)
    genus_data = genus_generating_function_coefficients(c, max_genus=4)

    return {
        'c': c,
        'c_dual': c_dual,
        'line_comparison_fixed': simplify(c - 13) == 0,
        'kappa': kappa,
        'kappa_dual': kappa_dual,
        'kappa_sum': simplify(kappa + kappa_dual),
        'Q_contact': simplify(q_contact),
        'shadow_class': virasoro_shadow_class(),
        'shadow_depth': virasoro_shadow_depth(),
        'F_g_raw': {g: simplify(v) for g, v in genus_data['raw'].items()},
        'delta_H1': simplify(genus1_hessian_correction_virasoro(c)),
    }
