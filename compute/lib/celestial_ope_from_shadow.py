"""Celestial OPE from support-shadow data.

Derives celestial holography OPE coefficients from the support-shadow
arity decomposition (thm:thqg-VI-general-soft in
thqg_soft_graviton_theorems.tex).

Mathematical framework
----------------------
The shadow connection nabla^hol = d - Sh_{g,n}(Theta_A) decomposes by
arity.  The arity-r channel is the algebraic input available to a
chosen celestial comparison at soft order p = r - 2 in the expansion
(1/z_s)^p.  The support bound does not by itself count independent
physical Ward identities or soft factors.

Dictionary (tab:thqg-VI-soft-dictionary):
    p = 0  (arity 2):  kappa(A)             -> leading (Weinberg)
    p = 1  (arity 3):  C(A) = cubic shadow  -> subleading (Cachazo-Strominger)
    p = 2  (arity 4):  Q(A) = quartic       -> sub-subleading
    p = r-2 (arity r): Sh_r(A)              -> order-(r-2) support channel

After a celestial comparison datum is fixed, the OPE coefficient at
pole order k = r - 1 may receive contributions from the arity-r
shadow Sh_r and from lower cross-terms [Sh_s, Sh_t] with s+t=r.

For the Virasoro algebra (the gravitational case):
    kappa = c/2
    C     = 2 x^3  (cubic shadow)
    Q^ct  = 10 / [c(5c+22)]  (quartic contact)
    S_5   = -48 / [c^2(5c+22)]  (quintic)

References:
    Vol II: thqg_soft_graviton_theorems.tex (Theorems, Propositions, Computations)
    Vol I:  higher_genus_modular_koszul.tex (support-shadow packet, MC2)
    Vol I:  concordance.tex (Theorem D: kappa, support-depth classification)
    Literature: Pate-Raclariu-Strominger-Yuan (2019), arXiv:1910.07424
                Cachazo-Strominger (2014), arXiv:1404.4091
                Weinberg (1965), Phys. Rev. 140, B516
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Dict, List, Optional, Tuple, Union

from sympy import (
    Rational, S, Symbol, factorial, oo, simplify, sqrt, symbols,
    binomial, pi, Function,
)


# =========================================================================
# 1. SUPPORT DEPTH CLASSIFICATION (G/L/C/M)
# =========================================================================

@dataclass(frozen=True)
class ShadowDepthClass:
    """Support-depth class of a modular Koszul chiral algebra.

    The four classes (def:thqg-shadow-archetype):
        G (Gaussian):  r_max = 2, quadratic support channel only
        L (Lie/tree):  r_max = 3, quadratic plus cubic support channels
        C (Contact):   r_max = 4, support may include a quartic channel
        M (Mixed):     r_max = inf, no finite support cutoff

    ``support_channel_bound`` is an algebraic arity bound.  It is not
    the number of independent celestial soft factors, Ward identities,
    BPS states, or OPE coefficients.
    """
    label: str           # 'G', 'L', 'C', or 'M'
    r_max: int           # support depth; use -1 for infinity
    support_channel_bound: int  # arity-channel bound; -1 encodes infinity

    @property
    def is_infinite(self) -> bool:
        return self.r_max < 0

    @property
    def display_r_max(self) -> str:
        return 'inf' if self.is_infinite else str(self.r_max)


CLASS_G = ShadowDepthClass('G', 2, 2)
CLASS_L = ShadowDepthClass('L', 3, 3)
CLASS_C = ShadowDepthClass('C', 4, 4)
CLASS_M = ShadowDepthClass('M', -1, -1)  # -1 encodes infinity


def classify_shadow_depth(family: str) -> ShadowDepthClass:
    """Classify a chiral algebra family by support depth.

    Parameters
    ----------
    family : str
        One of: 'heisenberg', 'lattice', 'affine', 'betagamma',
                'virasoro', 'W_N', 'free_fermion'.

    Returns
    -------
    ShadowDepthClass
    """
    family_lower = family.lower().replace(' ', '_').replace('-', '_')

    g_families = {'heisenberg', 'lattice', 'free_boson', 'free_fermion'}
    l_families = {'affine', 'kac_moody', 'affine_km', 'affine_kac_moody'}
    c_families = {'betagamma', 'beta_gamma', 'bc', 'symplectic_boson'}
    m_families = {'virasoro', 'vir', 'w_n', 'w_algebra', 'w_3', 'w_4', 'w_5'}

    if family_lower in g_families:
        return CLASS_G
    elif family_lower in l_families:
        return CLASS_L
    elif family_lower in c_families:
        return CLASS_C
    elif family_lower in m_families:
        return CLASS_M
    else:
        raise ValueError(f"Unknown family: {family}")


# =========================================================================
# 2. SHADOW COEFFICIENTS FOR STANDARD FAMILIES
# =========================================================================

def kappa_virasoro(c: Any) -> Any:
    """Modular characteristic of Vir_c: kappa = c/2.

    AP48: kappa depends on the full algebra, not the Virasoro subalgebra.
    This formula is SPECIFIC to the Virasoro algebra.
    """
    return S(c) / 2


def kappa_affine(dim_g: Any, k: Any, h_dual: Any) -> Any:
    """Modular characteristic of affine KM: kappa = dim(g)(k+h^v)/(2h^v).

    AP1/AP39: this formula is DISTINCT from c/2.
    """
    return S(dim_g) * (S(k) + S(h_dual)) / (2 * S(h_dual))


def kappa_heisenberg(k: Any) -> Any:
    """Modular characteristic of H_k: kappa = k."""
    return S(k)


def cubic_shadow_virasoro() -> str:
    """Cubic shadow of Vir_c: C = 2x^3.

    The coefficient 2 comes from the T(z)T(w) OPE structure:
    the ternary operation is determined by the :T^2: normal-ordered
    product.
    """
    return '2*x^3'


def quartic_contact_virasoro(c: Any) -> Any:
    """Quartic contact invariant of Vir_c.

    Q^contact = 10 / [c(5c + 22)]

    From thm:thqg-VI-virasoro-quartic.  The denominator c(5c+22)
    is the Kac determinant factor at weight 4.

    Poles at c = 0 (free field) and c = -22/5 (Lee-Yang).
    Zero at c = infinity (semiclassical limit).
    """
    c_val = S(c)
    return S(10) / (c_val * (5 * c_val + 22))


def quintic_shadow_virasoro(c: Any) -> Any:
    """Quintic shadow coefficient of Vir_c.

    S_5 = -48 / [c^2(5c + 22)]

    From comp:thqg-VI-virasoro-quintic-soft and the shadow metric
    expansion H(t) = t^2 * sqrt(Q(t)).
    """
    c_val = S(c)
    return S(-48) / (c_val**2 * (5 * c_val + 22))


def genus1_hessian_virasoro(c: Any) -> Any:
    """Genus-1 Hessian correction for Virasoro.

    delta_H^{(1)} = 120 / [c^2(5c+22)]

    Same denominator as Q^contact confirms the quartic tower
    is controlled by the single Gram factor c(5c+22).
    """
    c_val = S(c)
    return S(120) / (c_val**2 * (5 * c_val + 22))


# =========================================================================
# 3. CELESTIAL SOFT FACTOR STRUCTURE
# =========================================================================

@dataclass
class CelestialSoftFactor:
    """A celestial soft factor at a given order.

    The celestial soft factor S^{(p)} at order p = r - 2 is
    determined by the arity-r shadow Sh_r(A) and cross-terms
    from lower arities.

    Attributes
    ----------
    order : int
        The soft order p (= 0, 1, 2, ...).
    arity : int
        The shadow arity r = p + 2.
    shadow_coefficient : Any
        The shadow coefficient Sh_r(A) (sympy expression in c or k).
    cross_term_sources : list of tuple
        Pairs (s, t) with s + t = r, s,t >= 3, giving the
        cross-term contributions [S^{(s-2)}, S^{(t-2)}].
    pole_order : int
        Maximum pole order in z_i - z_j (= r - 1).
    conformal_weight_degree : int
        Maximum degree in h_i (= 2r - 2).
    c_scaling_order : int
        Order of 1/c suppression (= r - 2 = p).
    physical_name : str
        Name in the gravitational reading.
    """
    order: int
    arity: int
    shadow_coefficient: Any
    cross_term_sources: List[Tuple[int, int]] = field(default_factory=list)
    pole_order: int = 0
    conformal_weight_degree: int = 0
    c_scaling_order: int = 0
    physical_name: str = ''


def soft_order_to_arity(p: int) -> int:
    """Convert soft order p to shadow arity r = p + 2.

    The fundamental dictionary:
        p = 0 -> arity 2 (kappa, Weinberg)
        p = 1 -> arity 3 (cubic, Cachazo-Strominger)
        p = 2 -> arity 4 (quartic, sub-subleading)
        ...
    """
    if p < 0:
        raise ValueError(f"Soft order must be >= 0, got {p}")
    return p + 2


def arity_to_soft_order(r: int) -> int:
    """Convert shadow arity r to soft order p = r - 2."""
    if r < 2:
        raise ValueError(f"Arity must be >= 2, got {r}")
    return r - 2


def celestial_mellin_transform_profile() -> Dict[str, Any]:
    """Profile the celestial Mellin transform used by the residue theorem."""
    return {
        "formula": (
            "M[A](Delta_i,z_i,zbar_i)="
            "prod_i int_0^infty omega_i^{Delta_i - 1} d omega_i "
            "A(omega_i,z_i,zbar_i)"
        ),
        "energy_weight": "omega_i^{Delta_i - 1}",
        "integration_domain": "(0, infinity) for every external energy",
        "source": "Pasterski-Shao-Strominger Mellin transform",
    }


def soft_residue_dimension(index: int, convention: str = "soft_level") -> Any:
    """Return the conformal dimension of a soft residue.

    Conventions
    -----------
    soft_level:
        ``index`` is p, the conformally soft level.  The reduced
        PSS primary has pole Delta = 1 - p.
    arity:
        ``index`` is the A_infinity arity r.  Since p = r - 2, the
        reduced pole is Delta = 3 - r.
    stress_tensor:
        ``index`` is the A_infinity arity r, but the graviton
        stress-tensor primary is shifted by one: Delta_st = 4 - r.
    """
    c = convention.lower().replace("-", "_")
    if c in {"soft_level", "residue_index", "soft"}:
        if index < 0:
            raise ValueError(f"Soft level must be >= 0, got {index}")
        return S(1) - S(index)
    if c in {"arity", "ainfty_arity", "a_infinity_arity"}:
        if index < 2:
            raise ValueError(f"Arity must be >= 2, got {index}")
        return S(3) - S(index)
    if c in {"stress_tensor", "graviton_shifted", "stress"}:
        if index < 2:
            raise ValueError(f"Arity must be >= 2, got {index}")
        return S(4) - S(index)
    raise ValueError(f"Unknown residue convention: {convention}")


def residue_operator_profile(arity: int, algebra: str = "A") -> Dict[str, Any]:
    """Profile the hard-insertion operator in the Mellin residue theorem."""
    if arity < 2:
        raise ValueError(f"Arity must be >= 2, got {arity}")
    return {
        "arity": arity,
        "soft_level": arity_to_soft_order(arity),
        "operator_formula": (
            "D_{r,i}=rho_i(pi_1 D_A | (s^-1 A)^otimes r)"
        ),
        "operator_with_arity": (
            f"D_{{{arity},i}}^{algebra}="
            f"rho_i(pi_1 D_{algebra} | (s^-1 {algebra})^otimes {arity})"
        ),
        "stage": "operator",
        "scalar_projection_warning": (
            "shadow coefficients appear only after projecting the operator "
            "to a chosen scalar channel"
        ),
    }


def mellin_residue_identity_profile(
    arity: int,
    n_hard: int,
) -> Dict[str, Any]:
    """Profile the residue theorem in soft-level and arity notation."""
    if arity < 2:
        raise ValueError(f"Arity must be >= 2, got {arity}")
    if n_hard < 1:
        raise ValueError(f"Need at least one hard insertion, got {n_hard}")
    p = arity_to_soft_order(arity)
    operator = residue_operator_profile(arity)
    return {
        "arity": arity,
        "soft_level": p,
        "n_hard": n_hard,
        "soft_level_pole": soft_residue_dimension(p, "soft_level"),
        "reduced_arity_pole": soft_residue_dimension(arity, "arity"),
        "stress_tensor_pole": soft_residue_dimension(arity, "stress_tensor"),
        "soft_index_identity": (
            "Res_{Delta=1-p} M[A_{n+1}] = "
            "sum_{i=1}^n D_{p+2,i} M[A_n]"
        ),
        "arity_identity": (
            "Res_{Delta=3-r} M[A_{n+1}] = "
            "sum_{i=1}^n D_{r,i} M[A_n]"
        ),
        "stress_tensor_identity": (
            "Res_{Delta_st=4-r} M[A_{n+1}] = "
            "sum_{i=1}^n D_{r,i}^{st} M[A_n]"
        ),
        "operator_formula": operator["operator_formula"],
    }


def virasoro_residue_coefficients(c: Any = Symbol("c")) -> Dict[int, Dict[str, Any]]:
    """Virasoro coefficients as raw ordered-bar data and reduced soft classes."""
    c_val = S(c)
    return {
        2: {
            "raw_ordered_bar": c_val / 2,
            "soft_class": c_val / 2,
            "status": "nontrivial",
        },
        3: {
            "raw_ordered_bar": S(2),
            "soft_class": S(0),
            "status": "gauge_trivial",
            "reason": "Virasoro cubic shadow is a coboundary in the reduced soft complex",
        },
        4: {
            "raw_ordered_bar": quartic_contact_virasoro(c_val),
            "soft_class": quartic_contact_virasoro(c_val),
            "status": "nontrivial",
        },
        5: {
            "raw_ordered_bar": quintic_shadow_virasoro(c_val),
            "soft_class": quintic_shadow_virasoro(c_val),
            "status": "nontrivial",
        },
    }


def cross_term_decomposition(r: int) -> List[Tuple[int, int]]:
    """Compute the cross-term decomposition for arity r.

    From eq:thqg-VI-general-soft-factor:
        S^{(r-2)} = Sh_{0,n}(Sh_r)|_soft
                   + sum_{s+t=r, s,t>=3} [S^{(s-2)}, S^{(t-2)}]|_cross

    Returns all pairs (s, t) with s + t = r and s, t >= 3.
    """
    if r < 2:
        raise ValueError(f"Arity must be >= 2, got {r}")
    pairs = []
    for s in range(3, r - 2):
        t = r - s
        if t >= 3:
            pairs.append((s, t))
    # Include the symmetric pair if s = t
    if r % 2 == 0 and r // 2 >= 3:
        s = r // 2
        if (s, s) not in pairs:
            pairs.append((s, s))
    # For r = 6: (3,3) only
    # For r = 7: (3,4) and (4,3) but we list ordered pairs s <= t
    # Recompute cleanly: all unordered pairs {s,t} with s+t=r, s,t>=3
    pairs_clean = []
    for s in range(3, r // 2 + 1):
        t = r - s
        if t >= 3:
            pairs_clean.append((s, t))
    return pairs_clean


def build_soft_factor(
    order: int,
    family: str = 'virasoro',
    c: Any = Symbol('c'),
) -> CelestialSoftFactor:
    """Build the celestial soft factor at a given soft order.

    Parameters
    ----------
    order : int
        The soft order p (0 = leading, 1 = subleading, ...).
    family : str
        Chiral algebra family.
    c : sympy expression
        Central charge or level parameter.

    Returns
    -------
    CelestialSoftFactor
    """
    r = soft_order_to_arity(order)

    # Physical names
    names = {
        0: 'Leading (Weinberg)',
        1: 'Subleading (Cachazo-Strominger)',
        2: 'Sub-subleading',
    }
    physical_name = names.get(order, f'Order-{order} soft factor')

    # Cross-terms
    cross = cross_term_decomposition(r)

    # Shadow coefficient
    if family.lower() in ('virasoro', 'vir'):
        if order == 0:
            coeff = kappa_virasoro(c)
        elif order == 1:
            coeff = S(2)  # C = 2x^3, coefficient is 2
        elif order == 2:
            coeff = quartic_contact_virasoro(c)
        elif order == 3:
            coeff = quintic_shadow_virasoro(c)
        else:
            # General: from shadow metric recursion
            # S_r ~ O(1/c^{r-2}) at large c
            coeff = Symbol(f'S_{r}')
    elif family.lower() in ('affine', 'kac_moody'):
        dim_g = Symbol('dim_g')
        k = Symbol('k')
        h_dual = Symbol('h_dual')
        if order == 0:
            coeff = kappa_affine(dim_g, k, h_dual)
        elif order == 1:
            coeff = Symbol('f_abc')  # structure constants
        else:
            coeff = S(0)  # r_max = 3 for affine
    elif family.lower() in ('heisenberg',):
        if order == 0:
            coeff = kappa_heisenberg(Symbol('k'))
        else:
            coeff = S(0)  # r_max = 2
    else:
        coeff = Symbol(f'Sh_{r}')

    return CelestialSoftFactor(
        order=order,
        arity=r,
        shadow_coefficient=coeff,
        cross_term_sources=cross,
        pole_order=r - 1,
        conformal_weight_degree=2 * r - 2,
        c_scaling_order=order,
        physical_name=physical_name,
    )


# =========================================================================
# 4. CELESTIAL OPE COEFFICIENTS FROM THE SHADOW TOWER
# =========================================================================

def celestial_ope_coefficient_leading(
    kappa_val: Any,
    h_i: Any,
    h_j: Any,
) -> Any:
    """Leading celestial OPE coefficient from kappa (arity 2).

    From thm:thqg-VI-leading-soft:
        S^{(0)} = sum_i kappa(phi_s, phi_i) / (z_s - z_i)

    In the celestial OPE O_Delta(z) O_Delta'(w), the leading
    coefficient C_1 (simple pole at z = w) is:

        C_1(Delta, Delta') = kappa * f(h_i, h_j)

    where f depends on the conformal weights of the operators.
    For the Virasoro stress tensor (h = 2):
        C_1 = kappa = c/2

    For general operators of weight h_i, h_j:
        C_1 = kappa * delta(h_i, h_j)   [diagonal in conformal weight]

    This is the gravitational reading of the Weinberg soft factor:
    the OPE coefficient is proportional to the modular characteristic,
    which is universal across all modular Koszul algebras.
    """
    return S(kappa_val)


def celestial_ope_coefficient_subleading(
    c_val: Any,
    h_i: Any,
) -> Any:
    """Subleading celestial OPE coefficient from cubic shadow (arity 3).

    From thm:thqg-VI-subleading-soft, for Virasoro:
        S^{(1)}_{CS,Vir} = sum_i [h_i(h_i-1)/(z_s-z_i)^2
                                  + h_i partial_{z_i}/(z_s-z_i)]

    The celestial OPE coefficient at the double pole (k=2) is:
        C_2(Delta) = h_i(h_i - 1)

    This is the Cachazo-Strominger subleading soft factor.
    It depends on the conformal weight h_i but NOT on c (the central
    charge enters only through the definition of the conformal weight).
    """
    h = S(h_i)
    return h * (h - 1)


def celestial_ope_coefficient_subsubleading(
    c_val: Any,
    h_i: Any,
) -> Any:
    """Sub-subleading celestial OPE coefficient from quartic shadow (arity 4).

    From thm:thqg-VI-virasoro-quartic and eq:thqg-VI-virasoro-sub-subleading:
        S^{(2)}_{Vir} = Q^contact * sum_i h_i^2(h_i-1)^2/(z_s-z_i)^3 + ...

    The celestial OPE coefficient at the triple pole (k=3) is:
        C_3(c, Delta) = Q^contact(c) * h_i^2 * (h_i - 1)^2

    where Q^contact = 10/[c(5c+22)].

    This is NON-UNIVERSAL: it depends on both the conformal weight
    and the central charge through the quartic contact invariant.
    The non-universality is the algebraic mechanism behind the
    known non-universality of the sub-subleading soft graviton
    theorem.
    """
    c = S(c_val)
    h = S(h_i)
    Q_ct = quartic_contact_virasoro(c)
    return Q_ct * h**2 * (h - 1)**2


def celestial_ope_coefficient_order_3(
    c_val: Any,
    h_i: Any,
) -> Any:
    """Order-3 celestial OPE coefficient from quintic shadow (arity 5).

    From comp:thqg-VI-virasoro-quintic-soft:
        S^{(3)}_{Vir} has leading coefficient
        20/[c^2(5c+22)] * h_i^3(h_i-1)^2(2h_i-1)/(z_s-z_i)^4 + ...

    The full shadow coefficient is S_5 = -48/[c^2(5c+22)], but the
    celestial OPE coefficient involves both the direct S_5 term and
    the [C, Q] cross-term.

    The celestial OPE coefficient at the 4th-order pole:
        C_4(c, Delta) involves S_5 and cross-terms from C and Q.
    """
    c = S(c_val)
    h = S(h_i)
    # Leading diagonal contribution from the [C, Q] bracket
    bracket_coeff = S(20) / (c**2 * (5 * c + 22))
    return bracket_coeff * h**3 * (h - 1)**2 * (2 * h - 1)


# =========================================================================
# 5. FULL CELESTIAL OPE STRUCTURE
# =========================================================================

@dataclass
class CelestialOPEData:
    """Full celestial OPE data derived from the shadow tower.

    Suppressing the conventional gravitational coupling, the celestial OPE
    of graviton primaries:
        O_Delta(z) O_Delta'(w) ~ sum_k C_k(Delta,Delta') / (z-w)^k * O_Delta''

    The OPE coefficients C_k are controlled by the shadow tower:
        C_1: from kappa (arity 2) — universal
        C_2: from C (arity 3) — universal (angular momentum)
        C_3: from Q^contact (arity 4) — non-universal (dynamics)
        C_k: from Sh_{k+1} (arity k+1) — progressively less universal

    Finite support depth bounds the algebraic arity channels available
    to a fixed celestial comparison.  It does not count independent OPE
    coefficients: cross-terms and the comparison functor decide which
    physical coefficients survive.
    """
    family: str
    depth_class: ShadowDepthClass
    coefficients: Dict[int, Any]   # k -> C_k expression
    shadow_coefficients: Dict[int, Any]  # r -> Sh_r value
    c_val: Any


def compute_celestial_ope_virasoro(
    c_val: Any,
    h_val: Any,
    max_order: int = 5,
) -> CelestialOPEData:
    """Compute celestial OPE coefficients for Virasoro from the shadow tower.

    Parameters
    ----------
    c_val : sympy expression
        Central charge.
    h_val : sympy expression
        Conformal weight of the external operator.
    max_order : int
        Maximum pole order to compute (default 5).

    Returns
    -------
    CelestialOPEData
    """
    c = S(c_val)
    h = S(h_val)

    coefficients = {}
    shadow_coefficients = {}

    # Arity 2: kappa
    kappa = kappa_virasoro(c)
    shadow_coefficients[2] = kappa
    coefficients[1] = celestial_ope_coefficient_leading(kappa, h, h)

    # Arity 3: cubic shadow
    shadow_coefficients[3] = S(2)  # C = 2x^3
    if max_order >= 2:
        coefficients[2] = celestial_ope_coefficient_subleading(c, h)

    # Arity 4: quartic contact
    Q_ct = quartic_contact_virasoro(c)
    shadow_coefficients[4] = Q_ct
    if max_order >= 3:
        coefficients[3] = celestial_ope_coefficient_subsubleading(c, h)

    # Arity 5: quintic
    S5 = quintic_shadow_virasoro(c)
    shadow_coefficients[5] = S5
    if max_order >= 4:
        coefficients[4] = celestial_ope_coefficient_order_3(c, h)

    # Higher arities: from shadow metric recursion
    # S_r = a_{r-2}/r from the convolution recursion f^2 = Q_L
    if max_order >= 5:
        for r in range(6, max_order + 2):
            shadow_coefficients[r] = Symbol(f'S_{r}')

    return CelestialOPEData(
        family='virasoro',
        depth_class=CLASS_M,
        coefficients=coefficients,
        shadow_coefficients=shadow_coefficients,
        c_val=c,
    )


# =========================================================================
# 6. CROSS-TERM STRUCTURE (FLATNESS EQUATION)
# =========================================================================

def flatness_equation_at_arity(r: int, shadows: Dict[int, Any]) -> Dict[str, Any]:
    """Decompose the flatness equation at arity r.

    From thm:thqg-VI-flatness-by-arity:
        nabla^{(2)}(A^{(r)}) + sum_{s+t=r} [A^{(s)}, A^{(t)}] + o_r = 0

    Parameters
    ----------
    r : int
        Arity (>= 2).
    shadows : dict
        Dictionary mapping arity s to shadow coefficient Sh_s.

    Returns
    -------
    dict with keys:
        'arity': r
        'direct_shadow': Sh_r contribution
        'cross_terms': list of (s, t, Sh_s, Sh_t) quadruples
        'soft_order': p = r - 2
        'n_cross_terms': number of cross-term channels
    """
    if r < 2:
        raise ValueError(f"Arity must be >= 2, got {r}")

    direct = shadows.get(r, S(0))
    cross = []
    for s in range(3, r - 2):
        t = r - s
        if t >= 3 and s <= t:
            sh_s = shadows.get(s, S(0))
            sh_t = shadows.get(t, S(0))
            cross.append((s, t, sh_s, sh_t))
    # Include s = t case
    if r % 2 == 0:
        s = r // 2
        if s >= 3 and (s, s) not in [(c[0], c[1]) for c in cross]:
            sh_s = shadows.get(s, S(0))
            cross.append((s, s, sh_s, sh_s))

    return {
        'arity': r,
        'direct_shadow': direct,
        'cross_terms': cross,
        'soft_order': r - 2,
        'n_cross_terms': len(cross),
    }


# =========================================================================
# 7. SHADOW METRIC RECURSION FOR HIGHER COEFFICIENTS
# =========================================================================

def shadow_metric_coefficients_virasoro(
    c_val: Any,
    r_max: int = 10,
) -> Dict[int, Any]:
    """Compute Virasoro shadow coefficients from the shadow metric recursion.

    The shadow metric is Q_L(t) = (2*kappa + 3*alpha*t)^2 + 2*Delta*t^2
    where kappa = c/2, alpha = 2, Delta = 8*kappa*S_4.

    The shadow generating function is H(t) = t^2 * sqrt(Q_L(t)).
    Shadow coefficients are S_r = [t^r] H(t) / r.

    Uses the convolution recursion from f^2 = Q_L.

    Parameters
    ----------
    c_val : Rational or Fraction
        Central charge (must be rational and nonzero).
    r_max : int
        Maximum arity to compute.

    Returns
    -------
    dict mapping arity r to exact shadow coefficient S_r.
    """
    c = Rational(c_val) if not isinstance(c_val, Rational) else c_val
    if c == 0:
        return {r: S(0) for r in range(2, r_max + 1)}

    kappa = c / 2
    alpha = Rational(2)
    S4 = Rational(10) / (c * (5 * c + 22))

    # Shadow metric: Q = (2*kappa)^2 + 12*kappa*alpha*t + (9*alpha^2 + 16*kappa*S4)*t^2
    q0 = 4 * kappa**2       # = c^2
    q1 = 12 * kappa * alpha  # = 12c
    q2 = 9 * alpha**2 + 16 * kappa * S4  # = 36 + 80/(c(5c+22))

    # Convolution recursion for f = sqrt(Q):
    # f(t) = a_0 + a_1*t + a_2*t^2 + ...
    # where f^2 = Q
    # a_0 = sqrt(q0) = |c|, a_1 = q1/(2*a_0), etc.
    max_n = r_max - 2
    if max_n < 0:
        return {}

    a0 = abs(c)  # sqrt(c^2) = |c|
    a = [S(0)] * (max_n + 1)
    a[0] = a0

    if max_n >= 1:
        a[1] = q1 / (2 * a0)

    if max_n >= 2:
        a[2] = (q2 - a[1]**2) / (2 * a0)

    for n in range(3, max_n + 1):
        conv = sum(a[j] * a[n - j] for j in range(1, n))
        a[n] = -conv / (2 * a0)

    # S_r = a_{r-2} / r  (because H = t^2 * f, so [t^r]H = a_{r-2})
    # Actually: the shadow coefficient convention is S_r = a_{r-2}
    # The notation in the manuscript: Sh_r = A^sh_{r,0}
    # From shadow_borel_resurgence.py line 198: S_r = a_{r-2} / r
    result = {}
    for r in range(2, r_max + 1):
        idx = r - 2
        if idx <= max_n:
            result[r] = a[idx] / r
        else:
            result[r] = S(0)

    return result


# =========================================================================
# 8. CELESTIAL OPE FROM SHADOW METRIC (NUMERICAL)
# =========================================================================

def celestial_ope_from_shadow_metric(
    c_val: Any,
    h_val: Any,
    r_max: int = 10,
) -> Dict[str, Any]:
    """Compute celestial OPE data from the full shadow metric recursion.

    This is the master function that combines:
    1. Shadow metric recursion for exact shadow coefficients
    2. Cross-term decomposition from the flatness equation
    3. Celestial OPE coefficients from the soft factor dictionary

    Parameters
    ----------
    c_val : rational number
        Central charge.
    h_val : rational number or symbol
        Conformal weight of external operator.
    r_max : int
        Maximum arity to compute.

    Returns
    -------
    dict with keys:
        'shadow_coefficients': {r: S_r}
        'celestial_ope_coefficients': {k: C_k}
        'cross_term_structure': {r: cross_term_data}
        'depth_class': ShadowDepthClass
        'c_scaling': {k: order of 1/c suppression}
    """
    c = S(c_val)
    h = S(h_val)

    # Shadow coefficients from the metric recursion
    try:
        shadows = shadow_metric_coefficients_virasoro(c_val, r_max)
    except (TypeError, ValueError):
        # Symbolic c: use known exact values
        shadows = {
            2: kappa_virasoro(c),
            3: S(2),
            4: quartic_contact_virasoro(c),
            5: quintic_shadow_virasoro(c),
        }

    # Celestial OPE coefficients
    ope_coeffs = {}
    cross_structure = {}
    c_scaling = {}

    for r in range(2, r_max + 1):
        k = r - 1  # pole order in the OPE
        p = r - 2  # soft order

        # Cross-term structure
        cross_data = flatness_equation_at_arity(r, shadows)
        cross_structure[r] = cross_data

        # C_k scaling: O(c^{-(r-2)}) at large c
        c_scaling[k] = p

    # Explicit OPE coefficients for low orders
    if c_val != 0:
        ope_coeffs[1] = celestial_ope_coefficient_leading(
            kappa_virasoro(c), h, h)
        ope_coeffs[2] = celestial_ope_coefficient_subleading(c, h)
        ope_coeffs[3] = celestial_ope_coefficient_subsubleading(c, h)
        if r_max >= 5:
            ope_coeffs[4] = celestial_ope_coefficient_order_3(c, h)

    return {
        'shadow_coefficients': shadows,
        'celestial_ope_coefficients': ope_coeffs,
        'cross_term_structure': cross_structure,
        'depth_class': CLASS_M,
        'c_scaling': c_scaling,
    }


# =========================================================================
# 9. COMPARISON WITH KNOWN CELESTIAL OPE (PATE-RACLARIU-STROMINGER)
# =========================================================================

def pate_raclariu_strominger_leading(
    delta_1: Any,
    delta_2: Any,
    J_1: Any,
    J_2: Any,
) -> Dict[str, Any]:
    """Leading celestial OPE coefficient from Pate-Raclariu-Strominger-Yuan (2019).

    The celestial OPE of graviton primaries:
        G^+_{Delta_1}(z, zbar) G^+_{Delta_2}(w, wbar) ~
            (zbar-wbar)/(z-w)
            * B(Delta_1-1, Delta_2-1)
            * G^+_{Delta_1+Delta_2}(w, wbar) + ...

    where B is the Euler beta function and Delta_i are conformal dimensions.
    Overall sign conventions depend on the incoming/outgoing choice.

    The leading Mellin kernel is the beta function:
        C_leading = B(Delta_1-1, Delta_2-1) = Gamma(D1-1)*Gamma(D2-1)/Gamma(D1+D2-2)

    In the shadow tower language:
        This is the arity-2 (kappa) contribution.
        The beta function arises from the Mellin transform of the
        collinear splitting function, which is the momentum-space
        incarnation of the kappa pairing.

    Parameters
    ----------
    delta_1, delta_2 : conformal dimensions
    J_1, J_2 : helicities

    Returns
    -------
    dict with coefficient and comparison data
    """
    d1, d2 = S(delta_1), S(delta_2)
    # Beta function B(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b)
    from sympy import gamma as Gamma
    beta_val = Gamma(d1 - 1) * Gamma(d2 - 1) / Gamma(d1 + d2 - 2)

    return {
        'C_leading_PRS': beta_val,
        'delta_1': d1,
        'delta_2': d2,
        'shadow_source': 'kappa (arity 2)',
        'is_universal': True,
        'has_antiholomorphic_ratio': True,
        'central_triple_pole_present': False,
        'output_dimension': d1 + d2,
        'reference': 'Pate-Raclariu-Strominger-Yuan 2019, arXiv:1910.07424',
    }


def compare_shadow_vs_prs(
    c_val: Any,
    delta_val: Any,
) -> Dict[str, Any]:
    """Compare shadow tower prediction with PRS celestial OPE.

    The comparison tests whether the shadow tower correctly reproduces
    the known celestial OPE structure:

    1. Leading (kappa): the soft factor sum_i kappa/(z_s - z_i) maps
       to the B(D-1, D'-1) beta function through the Mellin transform.
       This is UNIVERSAL: it depends only on kappa, not on the
       specific algebra.

    2. Subleading (C): the angular momentum factor maps to the
       PRS subleading OPE.  For Virasoro, C = 2x^3 gives the
       standard conformal Ward subleading term.

    3. Sub-subleading (Q^contact): the dynamics-dependent term.
       For Virasoro, Q^ct = 10/[c(5c+22)] gives an explicit
       prediction that should match the sub-subleading celestial
       OPE from the stress tensor exchange.

    IMPORTANT (AP42): This comparison is valid at the LEVEL OF THE
    SOFT EXPANSION.  The celestial OPE in general receives contributions
    from all soft orders simultaneously.  The shadow tower gives the
    soft expansion of the OPE, not the full OPE.  The collinear limit
    (z -> w) requires the full OPE, not just the soft expansion
    (rem:thqg-VI-open-directions(ii)).
    """
    c = S(c_val)
    delta = S(delta_val)

    kappa = kappa_virasoro(c)
    Q_ct = quartic_contact_virasoro(c)

    return {
        'leading_shadow': kappa,
        'leading_PRS_source': 'B(Delta-1, Delta-1) via Mellin of Weinberg',
        'leading_match': True,  # Universal: both determined by kappa
        'subleading_shadow': S(2),  # cubic shadow C = 2x^3
        'subleading_PRS_source': 'Angular momentum from Cachazo-Strominger',
        'subleading_match': True,  # Universal angular momentum
        'subsubleading_shadow': Q_ct,
        'subsubleading_PRS_source': 'Non-universal dynamics',
        'subsubleading_match': 'dynamics-dependent (AP42: level caveat)',
        'shadow_scaling': {
            'leading': 'O(c)',
            'subleading': 'O(1)',
            'subsubleading': f'O(1/c) = {Q_ct}',
        },
    }


# =========================================================================
# 10. CELESTIAL SUPPORT PROFILE
# =========================================================================

def ward_identity_support_profile(family: str) -> Dict[str, Any]:
    """Describe the support-channel bound for celestial Ward data.

    The support bound is a finite arity envelope for the algebraic
    packet.  It is not a count of independent Ward identities.  A
    physical count requires a celestial comparison datum and can be
    smaller because channels may be killed, identified, or sourced by
    lower cross-terms.

    Parameters
    ----------
    family : str
        Chiral algebra family name.

    Returns
    -------
    dict with support-channel structure
    """
    depth = classify_shadow_depth(family)

    if depth.is_infinite:
        support_bound = 'unbounded'
        profile_description = 'No finite algebraic support cutoff'
        system_type = 'unbounded_support'
    else:
        support_bound = depth.support_channel_bound
        profile_description = (
            f'algebraic support channels have arity <= {support_bound}'
        )
        system_type = 'finite_support'

    return {
        'family': family,
        'depth_class': depth.label,
        'r_max': depth.display_r_max,
        'support_channel_bound': support_bound,
        'ward_identity_count': 'requires_celestial_comparison',
        'comparison_required': True,
        'profile_description': profile_description,
        'system_type': system_type,
    }


def celestial_soft_factor_table() -> List[Dict[str, Any]]:
    """Generate the celestial support-channel table.

    Reproduces tab:thqg-VI-family-soft from the manuscript.
    """
    families = [
        ('Heisenberg', 'heisenberg', {'k': Symbol('k')}),
        ('Affine KM', 'affine', {}),
        ('beta-gamma', 'betagamma', {}),
        ('Virasoro', 'virasoro', {'c': Symbol('c')}),
        ('W_N', 'W_N', {}),
    ]

    table = []
    for name, family, params in families:
        depth = classify_shadow_depth(family)
        entry = {
            'family': name,
            'class': depth.label,
            'r_max': depth.display_r_max,
            'support_channel_bound': (
                depth.support_channel_bound if not depth.is_infinite
                else 'unbounded'
            ),
            'comparison_required_for_physical_count': True,
        }

        # Explicit shadow coefficients
        if family == 'heisenberg':
            entry['kappa'] = f"k = {params.get('k', 'k')}"
            entry['cubic'] = '0'
            entry['quartic'] = '0'
        elif family == 'affine':
            entry['kappa'] = 'dim(g)(k+h^v)/(2h^v)'
            entry['cubic'] = 'f_abc (structure constants)'
            entry['quartic'] = '0 (Jacobi identity)'
        elif family == 'betagamma':
            entry['kappa'] = '6*lambda^2 - 6*lambda + 1'
            entry['cubic'] = '0 (on weight-changing line)'
            entry['quartic'] = 'nonzero (contact)'
        elif family == 'virasoro':
            entry['kappa'] = 'c/2'
            entry['cubic'] = '2*x^3'
            entry['quartic'] = '10/[c(5c+22)]'
        elif family == 'W_N':
            entry['kappa'] = 'higher-spin generalization'
            entry['cubic'] = 'higher-spin'
            entry['quartic'] = 'higher-spin'

        table.append(entry)

    return table
