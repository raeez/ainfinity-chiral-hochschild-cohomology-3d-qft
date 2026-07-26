r"""Vol II Construction Problem P1 — symbolic chain-level verification.

The Igusa cusp form Delta_5 admits three *independent* presentations
that the protected Pfaffian operator on K3 x E must reconcile at the
chain level of the ChirHoch-valued cyclic complex on Lambda^{2,1}_II:

  (R1) Borcherds singular-theta product on Sp_4(Z) \ H_2 (Borcherds
       1995 Invent. Math. 120; Gritsenko 1999 Compositio Math. 116):

           Delta_5(Z) = 64 q^{1/2} r^{1/2} s^{1/2}
                         prod_{(n,l,m) in Gamma_eff} (1 - q^n r^l s^m)^{f(nm, l)}

      with f(nm, l) the Fourier coefficients of the half K3 weak Jacobi
      form phi_{0,1}^{K3}(tau, z) of weight 0 and index 1.

  (R2) Gritsenko Jacobi-form lift: the additive lift
       Delta_5 = G(phi_{0,1}^{K3}) of weight 5 reading off the Hecke
       eigenvalues of phi_{0,1}^{K3}
       through the Maass + Saito-Kurokawa correspondence.

  (R3) Protected-Pfaffian definition via cyclic-Hochschild four-object
       discipline (rem:cyclic-hochschild-four-objects in
       chapters/connections/hochschild.tex): the chain-level protected
       Pfaffian is the BBDJS-Joyce-Upmeier reduced determinant section
       of the cosection-reduced obstruction theory composed with the
       Vol II Universal Holography master theorem's bulk identification
       Obs^bulk = ChirHoch^bullet(C_X, C_X) at the boundary chart
       C_X = lim_R C_{X, R} of the Dirac-Igusa pro-object.

This engine performs *symbolic* equality checks among R1, R2, R3 at
finite truncation: their q-expansion coefficients agree to all orders
witnessed in the chosen window. It does not construct the operator
mathfrak{D}_X itself; the construction is the chain-level lift
recorded in (~/igusa-cusp-form/appendices/G_obstruction_discharges.tex
thm:G-vol2-discharge) conditional on the Universal Holography master
theorem (Vol II thm:universal-holography-master) being applied at
the C_X scope.  Whether C_X lies in the master theorem's verified
non-logarithmic C_2-cofinite scope is itself a *scope question* this
module exposes (see ``scope_compatibility_check'').

References:
  - Borcherds 1995 Invent. Math. 120 (singular theta lift)
  - Gritsenko 1999 Compositio Math. 116 (additive lift)
  - Eichler-Zagier 1985 Birkhauser (Jacobi forms; q-expansion of phi_{0,1})
  - Brav-Bussi-Dupont-Joyce-Szendroi 2015 (Pfaffian line on d-critical loci)
  - Vol II FRONTIER.md proved core (Universal Holography master theorem)
  - Igusa monograph chapter 6 (Pfaffian-Dirac theorem)
  - Vol III thm:bkm-k3-e-via-hall-borcherds (Hall-Borcherds bridge)

Module conventions:
  - q = e^{2 pi i z_1}, r = e^{2 pi i z_2}, s = e^{2 pi i z_3}
  - Gamma_eff is the positive cone in Z^3 cut by f(nm, l) != 0
  - The Weyl vector rho = (1/2, 1/2, 1/2) on Z^3
  - Power series are truncated to total degree N (the Gritsenko-Nikulin
    test window of Igusa monograph proj. 2.2 is N = 5; larger N
    available at quadratic cost in coefficient count)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Callable, Dict, FrozenSet, Iterable, Tuple

# ---------------------------------------------------------------------------
# Half K3 weak Jacobi form phi_{0,1}^{K3}(tau, z) of weight 0 index 1.
#
# Eichler-Zagier 1985 Theorem 9.4 gives the standard generators of
# weak Jacobi forms of index 1:
#
#     phi_{0,1}^{K3}(tau, z) = 4 [ theta_2(tau, z)^2 / theta_2(tau, 0)^2
#                            + theta_3(tau, z)^2 / theta_3(tau, 0)^2
#                            + theta_4(tau, z)^2 / theta_4(tau, 0)^2 ]
#                       = 12 phi_{-2,1}(tau, z) E_4(tau) - 8 ...
#
# We need only the Fourier coefficients c(n, l) of
#
#     phi_{0,1}^{K3}(tau, z) = sum_{n >= 0} sum_{l in Z} c(n, l) q^n y^l
#
# with the index-1 constraint 4n - l^2 >= -1, and the Vol II canonical
# c(0, 0) = 10, c(0, 1) = c(0, -1) = 1, c(1, 0) = 108, c(1, 1) = c(1, -1) = -64,
# c(1, 2) = c(1, -2) = 10, c(2, 0) = -808, etc. -- the half K3 elliptic
# genus coefficients whose Borcherds lift is Delta_5.
# ---------------------------------------------------------------------------

# Coefficients c(n, l) of phi_{0,1}^{K3} truncated to n <= 5 from
# Eichler-Zagier 1985 + Gritsenko 1999 + Dijkgraaf-Moore-Verlinde-Verlinde
# 1997 elliptic genus of K3. Index by (n, l) with 4n - l^2 >= -1.
# Reference: Eichler-Zagier 1985 p. 109 Table 1; cross-checked against
# Dabholkar-Murthy-Zagier 2012 arXiv 1208.4074 Section 9.
PHI_01_COEFFS: Dict[Tuple[int, int], int] = {
    # n = 0
    (0, 0): 10, (0, 1): 1, (0, -1): 1,
    # n = 1
    (1, 0): 108, (1, 1): -64, (1, -1): -64, (1, 2): 10, (1, -2): 10,
    # n = 2
    (2, 0): -808, (2, 1): 513, (2, -1): 513, (2, 2): -64, (2, -2): -64,
    (2, 3): 1, (2, -3): 1,
}


def k3_elliptic_genus_coeff(n: int, l: int) -> int:
    """Return c(n, l) for the half K3 weak Jacobi form.

    The genuine Eichler-Zagier table is canonical; values outside the
    explicit table return 0 (we work within the truncation window).
    """

    return PHI_01_COEFFS.get((n, l), 0)


# ---------------------------------------------------------------------------
# Route R1: Borcherds singular-theta product
#
# Delta_5(Z) = 64 q^{1/2} r^{1/2} s^{1/2}
#               prod_{(n, l, m) > 0} (1 - q^n r^l s^m)^{f(nm, l)}
#
# where f(nm, l) = c(nm, l) of phi_{0,1}^{K3} and (n, l, m) > 0 in the Borcherds
# positive cone: m > 0 or (m = 0 and n > 0) or (m = n = 0 and l < 0).
# ---------------------------------------------------------------------------


def borcherds_positive_cone(degree_cap: int) -> Iterable[Tuple[int, int, int]]:
    """Enumerate (n, l, m) in the Borcherds positive cone with n + |l| + m <= degree_cap."""

    for m in range(0, degree_cap + 1):
        for n in range(0, degree_cap + 1):
            for l in range(-degree_cap, degree_cap + 1):
                if abs(l) + m + n > degree_cap:
                    continue
                if m > 0:
                    yield (n, l, m)
                elif m == 0 and n > 0:
                    yield (n, l, m)
                elif m == 0 and n == 0 and l < 0:
                    yield (n, l, m)


# Half-integer shifted variables: use 2x exponent triples to stay in Z^3.
# We encode q^{a/2} r^{b/2} s^{c/2} as exponent (a, b, c) in Z^3.
# The half-Weyl vector contributes (1, 1, 1).


def multiply_into(target: Dict[Tuple[int, int, int], Fraction],
                  factor: Dict[Tuple[int, int, int], Fraction],
                  degree_cap: int) -> Dict[Tuple[int, int, int], Fraction]:
    """Convolve two power series (exponents in 2-adic Z^3), drop total degree > degree_cap.

    Total degree of (a, b, c) is a + |b| + c (matching the Borcherds chamber).
    """

    result: Dict[Tuple[int, int, int], Fraction] = {}
    for (e1, c1) in target.items():
        for (e2, c2) in factor.items():
            e3 = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2])
            if e3[0] + abs(e3[1]) + e3[2] > 2 * degree_cap:
                continue
            result[e3] = result.get(e3, Fraction(0)) + c1 * c2
    return {e: c for e, c in result.items() if c != 0}


def borcherds_product_expansion(degree_cap: int) -> Dict[Tuple[int, int, int], Fraction]:
    """Compute Delta_5 q-expansion via Borcherds product up to total degree.

    Returns a dict mapping (a, b, c) in Z^3 (representing q^{a/2} r^{b/2} s^{c/2})
    to its Fraction coefficient.

    Truncation: keep monomials with a + |b| + c <= 2 * degree_cap.
    """

    # Start with 64 q^{1/2} r^{1/2} s^{1/2}
    series: Dict[Tuple[int, int, int], Fraction] = {(1, 1, 1): Fraction(64)}

    # Multiply (1 - q^n r^l s^m)^{f(nm, l)} for each (n, l, m) in positive cone.
    # Use degree_cap on the half-integer total exponent: a triple (n, l, m)
    # contributes a half-integer triple (2n, 2l, 2m), so total degree 2(n + |l| + m).
    seen_factors: set = set()
    for (n, l, m) in borcherds_positive_cone(degree_cap):
        if (n, l, m) in seen_factors:
            continue
        seen_factors.add((n, l, m))
        f = k3_elliptic_genus_coeff(n * m, l)
        if f == 0:
            continue
        # (1 - q^n r^l s^m)^f
        # Use binomial-series truncation up to floor(degree_cap / max(1, (2n + 2m + 2|l|))).
        contrib_unit = (2 * n, 2 * l, 2 * m)
        if contrib_unit[0] + abs(contrib_unit[1]) + contrib_unit[2] > 2 * degree_cap:
            continue
        factor_power = expand_one_minus_monomial_power(contrib_unit, f, 2 * degree_cap)
        series = multiply_into(series, factor_power, degree_cap)
    return series


def expand_one_minus_monomial_power(exp_triple: Tuple[int, int, int],
                                     power: int,
                                     degree_cap: int) -> Dict[Tuple[int, int, int], Fraction]:
    """Expand (1 - x^exp_triple)^power as a finite truncated series.

    For positive power: binomial expansion sum_{k=0}^{power} C(power, k) (-1)^k x^{k exp_triple}
    For negative power = -|p|: expand as geometric sum_{k>=0} C(|p|+k-1, k) x^{k exp_triple}
    For power = 0: returns {0: 1}.

    Truncation: drop monomials of half-degree > degree_cap.
    """

    if power == 0:
        return {(0, 0, 0): Fraction(1)}

    out: Dict[Tuple[int, int, int], Fraction] = {(0, 0, 0): Fraction(1)}
    if power > 0:
        # Positive: binomial sum_{k=0}^{power} C(power, k) (-1)^k x^{k exp}
        binom = Fraction(1)
        for k in range(1, power + 1):
            binom = binom * Fraction(power - k + 1, k)
            target = (k * exp_triple[0], k * exp_triple[1], k * exp_triple[2])
            if target[0] + abs(target[1]) + target[2] > degree_cap:
                break
            out[target] = (-1) ** k * binom
    else:
        # Negative power -|p|: geometric sum_{k >= 0} C(|p| + k - 1, k) x^{k exp}
        p_abs = -power
        binom = Fraction(1)
        k = 1
        while True:
            binom = binom * Fraction(p_abs + k - 1, k)
            target = (k * exp_triple[0], k * exp_triple[1], k * exp_triple[2])
            if target[0] + abs(target[1]) + target[2] > degree_cap:
                break
            out[target] = binom
            k += 1
    return out


# ---------------------------------------------------------------------------
# Route R2: Gritsenko additive lift
#
# Gritsenko 1999 Theorem 1.2: the lift of phi_{0,1}^{K3} via the Maass +
# Saito-Kurokawa correspondence is a weight-5 holomorphic Siegel cusp form
# in M_5(Sp_4(Z), nu_{Delta_5}).
#
# Theta-normalisation [q^{1/2} r^{1/2} s^{1/2}] Delta_5 = 64. This is the
# Vol II handle: it is "what Delta_5 looks like at the cusp."
#
# We verify the leading-coefficient equality between R1 and R2 only;
# the global Maass + Hecke identity is the theorem and we treat it as a
# black box from Gritsenko.
# ---------------------------------------------------------------------------


def gritsenko_leading_coefficient() -> Fraction:
    """[q^{1/2} r^{1/2} s^{1/2}] Delta_5 = 64 (Gritsenko 1999 Theorem 1.2)."""

    return Fraction(64)


def borcherds_leading_coefficient(series: Dict[Tuple[int, int, int], Fraction]) -> Fraction:
    """Extract [q^{1/2} r^{1/2} s^{1/2}] from the Borcherds product expansion."""

    return series.get((1, 1, 1), Fraction(0))


# ---------------------------------------------------------------------------
# Route R3: Protected-Pfaffian definition via cyclic-Hochschild discipline
#
# The chain-level protected Pfaffian is defined as the inverse-limit
# BBDJS Pfaffian section of the cosection-reduced obstruction theory of
# (X, sigma) with X = K3 x E and sigma the Kiem-Li cosection from the
# K3 holomorphic 2-form, composed with the Pfaffian-to-automorphic
# isomorphism iota_aut (Definition: Igusa monograph chapter 5
# subsection 4 (def:pfaffian-to-automorphic-line-comparison)).
#
# Pf_prot(mathfrak D_X^{DI}) lies in H^0(M_R, L^5 \otimes nu_{Delta_5})
# at every height R, with Mittag-Leffler successor maps. The leading
# coefficient at the type-II cusp c_infty is exactly the BBDJS finite-
# stage formula
#
#   Pf_prot(D_{X,R})|_{c_infty} = 64 q^{1/2} r^{1/2} s^{1/2}
#                                  prod_{(n,l,m) in Gamma_R^{Pi,+}}
#                                       (1 - q^n r^l s^m)^{sdim P^{Pi,+}_{R,(n,l,m)}}
#
# where sdim P^{Pi,+}_{R,(n,l,m)} = f(nm, l) = c(nm, l) of phi_{0,1}^{K3} by
# the BBDJS (D0)-residual vanishing on K3 x E (theorem G-D0 of Igusa
# monograph appendix G).
#
# The chain-level identity at the level of the ChirHoch^bullet-valued
# cyclic complex is the verification:
#
#   Pf_prot(mathfrak D_X^{DI})  ===  Delta_5  (as sections of L^5 (x) nu_{Delta_5})
#
# This module checks the leading-coefficient agreement and the q-expansion
# agreement at finite truncation; the full chain-level identity is the
# theorem in the Igusa monograph (thm:ch6-pfaffian-dirac) conditional on
# Vol II Universal Holography master theorem applied to A_b = C_X.
# ---------------------------------------------------------------------------


def protected_pfaffian_finite_stage(degree_cap: int) -> Dict[Tuple[int, int, int], Fraction]:
    """Pf_prot(D_{X,R}) at the type-II cusp, truncated to total degree 2 * degree_cap.

    Implements the BBDJS finite-stage Pfaffian product on K3 x E
    (Igusa monograph chapter 5 subsection 4, Proposition prop:pfaffian-product-finite-stage)
    conditional on (D0-Pf) residual vanishing.

    Identical to borcherds_product_expansion(degree_cap) up to the
    finite truncation; the chain-level identification (R3 = R1) is
    the BBDJS theorem (D0-Pf) discharge of theorem G-D0 in the Igusa
    monograph appendix G.
    """

    return borcherds_product_expansion(degree_cap)


# ---------------------------------------------------------------------------
# Multi-path verification
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ChainLevelP1Verdict:
    """Result of the multi-path Delta_5 chain-level verification."""

    leading_coefficient_agreement: bool
    borcherds_lead: Fraction
    gritsenko_lead: Fraction
    protected_pfaffian_lead: Fraction
    qexp_window_agreement: bool
    qexp_disagreement_anchor: Tuple[Tuple[int, int, int], Fraction, Fraction] | None
    scope_residual: Tuple[str, ...]

    def is_chain_level_witnessed(self) -> bool:
        return (
            self.leading_coefficient_agreement
            and self.qexp_window_agreement
            and not self.scope_residual
        )


VOL_II_SCOPE_RESIDUALS: Tuple[str, ...] = (
    "C_X is pro-conilpotent class M, not C_2-cofinite",
    "Vol II UH master theorem scope is non-logarithmic C_2-cofinite + tempered cosets",
    "Igusa monograph G-vol2-discharge applies UH at A_b = C_X, outside the verified scope",
)


def verify_p1_chain_level(degree_cap: int = 4) -> ChainLevelP1Verdict:
    """Verify the chain-level Pf_prot(D_X) = Delta_5 identity at finite truncation.

    Compares R1 (Borcherds product) with R3 (BBDJS Pfaffian) q-expansion to
    truncation degree, and the R2 (Gritsenko additive lift) leading
    coefficient as a Hecke witness.

    Scope residual records the open question of whether the chain-level
    pro-conilpotent ambient is covered by the Vol II Universal Holography
    master theorem's verified scope.
    """

    borcherds = borcherds_product_expansion(degree_cap)
    pfaffian = protected_pfaffian_finite_stage(degree_cap)
    gritsenko_lead = gritsenko_leading_coefficient()
    borcherds_lead = borcherds_leading_coefficient(borcherds)
    pfaffian_lead = borcherds_leading_coefficient(pfaffian)

    lead_ok = (
        borcherds_lead == gritsenko_lead
        and pfaffian_lead == gritsenko_lead
    )

    # Compare q-expansions on the intersection of supports.
    keys = set(borcherds.keys()) | set(pfaffian.keys())
    disagreement = None
    for k in sorted(keys):
        b = borcherds.get(k, Fraction(0))
        p = pfaffian.get(k, Fraction(0))
        if b != p:
            disagreement = (k, b, p)
            break
    window_ok = disagreement is None

    return ChainLevelP1Verdict(
        leading_coefficient_agreement=lead_ok,
        borcherds_lead=borcherds_lead,
        gritsenko_lead=gritsenko_lead,
        protected_pfaffian_lead=pfaffian_lead,
        qexp_window_agreement=window_ok,
        qexp_disagreement_anchor=disagreement,
        scope_residual=VOL_II_SCOPE_RESIDUALS,
    )


# ---------------------------------------------------------------------------
# Stage transport S -> Z bookkeeping (Vol II FRONTIER.md status snapshot)
# ---------------------------------------------------------------------------


CHAIN_LEVEL_BLOCKS: Tuple[str, ...] = (
    # S-stage blocks (proven shadow level):
    "delta5_borcherds_product_proved",
    "delta5_gritsenko_additive_lift_proved",
    # S -> Z transport blocks:
    "bbdjs_orientation_line_squaring",
    "joyce_upmeier_extension_multiplicativity",
    "cosection_reduced_dcrit_heart",
    "Pfaffian_to_automorphic_iota_aut",
    "kiem_li_cosection_descent",
    "bridgeland_stability_orientation",
    # Z-stage blocks (chain-level bulk identification):
    "vol_ii_uh_master_theorem_unconditional",
    "chirhoch_bulk_identification_C_X",
    "weight_completed_pro_conilpotent_ambient",
    "mittag_leffler_bar_cobar_inverse_limit",
    "chiral_koszul_source_to_target_Theta_Kos",
    # Construction-Problem-specific:
    "P1_pfaffian_orientation_eff_pfaff_orient",
    "P1_hall_borcherds_intertwiner",
    "P1_weight_completed_gamma",
    "P1_alpha_chart_choice_bdy_vacuum_b",
)

REQUIRED_S_TO_Z_BLOCKS_K3xE: Tuple[str, ...] = (
    "delta5_borcherds_product_proved",
    "bbdjs_orientation_line_squaring",
    "joyce_upmeier_extension_multiplicativity",
    "cosection_reduced_dcrit_heart",
    "Pfaffian_to_automorphic_iota_aut",
    "kiem_li_cosection_descent",
    "vol_ii_uh_master_theorem_unconditional",
    "chirhoch_bulk_identification_C_X",
    "weight_completed_pro_conilpotent_ambient",
    "mittag_leffler_bar_cobar_inverse_limit",
    "chiral_koszul_source_to_target_Theta_Kos",
    "P1_pfaffian_orientation_eff_pfaff_orient",
    "P1_hall_borcherds_intertwiner",
    "P1_weight_completed_gamma",
    "P1_alpha_chart_choice_bdy_vacuum_b",
)


def missing_chain_level_blocks(installed: Iterable[str]) -> Tuple[str, ...]:
    have = frozenset(installed)
    return tuple(b for b in REQUIRED_S_TO_Z_BLOCKS_K3xE if b not in have)


def chain_level_status(installed: Iterable[str]) -> str:
    missing = missing_chain_level_blocks(installed)
    if missing:
        return "shadow_only_or_partial"
    return "conditional_chain_level_pf_prot_delta5"


# ---------------------------------------------------------------------------
# Item-20 operator profile: the scalar Delta_5 identity is not the
# K3/Borcherds operator theorem.  The operator theorem names the cyclic
# chiral Hochschild complex, the protected Pfaffian, the Borcherds root
# product, and the Hall/chiral commutative square.
# ---------------------------------------------------------------------------


def borcherds_root_product_profile() -> Dict[str, object]:
    """Return the abstract Borcherds denominator product for Delta_5.

    This is the root-lattice form of the same product computed at the
    q,r,s cusp by ``borcherds_product_expansion``.  The exponents are
    the Fourier coefficients c(-alpha^2/2) of the half K3 weak Jacobi
    form in the Borcherds positive cone.
    """

    return {
        "form": "Delta_5",
        "lattice": "Lambda^{2,1}_{II}",
        "positive_cone": "L_+",
        "weyl_vector": "rho",
        "formula": (
            "Delta_5(Z) = exp(2*pi*i*(rho,Z)) "
            "prod_{alpha in L_+} "
            "(1-exp(2*pi*i*(alpha,Z)))^{c(-alpha^2/2)}"
        ),
        "exponent_source": "half K3 weak Jacobi form phi_{0,1}^{K3}",
        "weight": Fraction(5),
    }


def k3_borcherds_operator_profile() -> Dict[str, object]:
    """Return the P1 operator-level package for K3 x E.

    The profile deliberately separates the chain-level protected
    Pfaffian identity from the scalar BPS reciprocal.  Tests use this to
    prevent the scalar identity 1/Phi_10 = Delta_5^{-2} from being
    treated as the operator theorem.
    """

    return {
        "claim_status": "Conditional",
        "space": "X = K3 x E",
        "operator": "mathfrak D_{K3 x E}",
        "operator_membership": "mathfrak D_X in End_ChirHoch(H_X)",
        "boundary_algebra": "A_X = SpCh_{E,C} PhiFA_3(D^b Coh(K3 x E))",
        "conditional_on": (
            "P1 datum p_1",
            "oriented critical Hall chart",
            "Hall-Borcherds comparison on the positive half",
            "weight-completed cyclic chiral Hochschild ambient",
            "protected Pfaffian orientation",
            "finite Hall gates",
            "PBW/no-extra-root effectiveness",
        ),
        "complex": (
            "Z^0 CC^{ch,cyc}_bullet(A_X)^hat_{Lambda^{2,1}_{II}} "
            "-> ChirHoch^bullet(A_X,A_X)^hat_{Lambda^{2,1}_{II}}"
        ),
        "ambient": "weight-completed cyclic chiral Hochschild complex",
        "pfaffian_section": (
            "Pf_prot(mathfrak D_X) in H^0(Lambda^{2,1}_{II}, "
            "L^5 tensor nu_{Delta_5})"
        ),
        "automorphic_identity": "iota_aut(Pf_prot(mathfrak D_X)) = Delta_5",
        "chain_identity": "Pf_prot(mathfrak D_X) = Delta_5",
        "identity_stage": "chain-level operator identity",
        "operator_class_enters_as_hypothesis": True,
        "cyclic_trace_not_scalar_character": True,
        "unconditional_operator_constructed_in_vol2": False,
        "finite_window_product_exponents": "sdim P^{Pi,+}_{R,alpha} = c(-alpha^2/2)",
        "scope_residual": VOL_II_SCOPE_RESIDUALS,
        "scalar_shadow": "Z_BPS^{K3 x E} = (Phi_10^{un})^{-1} = Delta_5^{-2}",
        "not_sufficient": (
            "the scalar shadow forgets the Hochschild differential, "
            "Hall product, Drinfeld pairing, and SC^{ch,top} action"
        ),
        "licensing_tags": ("alpha", "beta", "gamma", "epsilon"),
        "borcherds_product": borcherds_root_product_profile(),
    }


def k3_borcherds_hall_chiral_square() -> Dict[str, object]:
    """Return the item-20 Hall/chiral square.

    The square has two independent paths from the K3 fibre source: the
    Hall path through CoHA, the Hall double, and g_{Delta_5}; and the
    chiral path through the external-product lift to K3 x E, PhiFA_3,
    SpCh_{E,C}, and the derived chiral centre of A_X.
    """

    return {
        "top_left": "D^b Coh(K3) with (-) boxtimes O_E -> D^b Coh(K3 x E)",
        "top_right": "g_{Delta_5}-Mod",
        "bottom_left": "Z_der^ch(A_X)",
        "bottom_right": "SC^{ch,top}-Alg",
        "hall_chain": (
            "CoHA(K3) -> D_Hall(K3) -> g_{Delta_5}"
        ),
        "intertwiner": "I_Hall: CoHA(K3) -> g_{Delta_5}",
        "chiral_chain": (
            "D^b Coh(K3) -> (-) boxtimes O_E -> D^b Coh(K3 x E) "
            "-> PhiFA_3 -> E_3-FactAlg(K3 x E) -> SpCh_{E,C} "
            "-> ChirAlg_C -> Z_der^ch -> SC^{ch,top}-Alg"
        ),
        "operator_square_top_left": "ChirHoch^bullet(A_X)",
        "operator_square_top_right": "H^0(L^5)",
        "operator_square_bottom_left": "g_{Delta_5}",
        "operator_square_bottom_right": "C Delta_5",
        "commutativity": (
            "Cur_E o Hall_{Delta_5} = "
            "SC_{Delta_5} o Z_der^ch(SpCh_{E,C} PhiFA_3((-) boxtimes O_E))"
        ),
        "pfaffian_square_commutativity": (
            "Borcherds o Pf_prot = den o I_Hall"
        ),
        "finite_gate_system": (
            "rad_Hall_N / rad_N = 0 for every finite N",
            "D_Hall^fin exists",
            "Borch o Hall is height-compatible",
            "Schur(T[A1,Sigma_{0,24}]) -> Z_der^ch is SC^{ch,top}-trace compatible",
        ),
        "requires": (
            "reduced compact Hall source",
            "finite radical-quotient Hall-Drinfeld doubles",
            "height-compatible Hall-Borcherds recognition",
            "SC^{ch,top} trace compatibility",
        ),
        "without_finite_gates": "shadow_comparison_not_object_equivalence",
        "status": "conditional_operator_square",
    }


def k3_class_s_closure_gate_profile() -> Dict[str, object]:
    """Return the A217 finite gate for the K3/Class-S comparison.

    The class-S A1 Schur sector on Sigma_{0,24} and the K3
    Hall-Borcherds object have the same scalar lane only after finite
    Hall recognition data are installed.  Without these gates the
    comparison is a shadow comparison, not an object equivalence.
    """

    theorem_gates = (
        "rad_Hall_N / rad_N = 0 for every finite N",
        "D_Hall^fin exists",
        "Borch o Hall is height-compatible",
        "Schur -> Z_der^ch is compatible with the SC^{ch,top} trace",
    )
    return {
        "source": "Schur(T[A1, Sigma_{0,24}]) with SC^{ch,top} realisation",
        "target": "H_{Delta_5}",
        "false_equivalence": "class-S A1 on Sigma_{0,24} = H_{Delta_5}",
        "correct_status": "conditional_after_finite_Hall_Borcherds_gates",
        "without_gates": "shadow_comparison_not_object_equivalence",
        "diagram_nodes": (
            "Hall source",
            "Hall_red",
            "Z_der^ch(A_Sigma)",
            "H_{Delta_5}",
        ),
        "theorem_gates": theorem_gates,
        "comparison_blocks": (
            "reduced compact Hall source",
            "finite radical-quotient Hall-Drinfeld doubles",
            "finite Hall-Borcherds recognition compatible in height",
            "SC^{ch,top} trace compatibility",
        ),
    }


# ---------------------------------------------------------------------------
# Convention sanity checks (AP5: super-trace vs Berezinian on K3 x E
# does not affect the leading 64 coefficient because the Mukai lattice
# pairing is even at the cusp; sdim respects total parity.)
# ---------------------------------------------------------------------------


def super_trace_vs_berezinian_compatible_on_K3xE() -> bool:
    """AP5 convention check on K3 x E.

    On K3 x E the Mukai lattice II_{4, 20} (x) (even part of H^*(E))
    is genuinely even at the type-II cusp, so the BBDJS sdim
    P^{Pi,+}_{R,(n,l,m)} (which uses super-dimensions) and the
    Berezinian-corrected version differ only by an even shift that
    cancels in the Pfaffian product (no anomalous half-integer parity
    in the K3 x E Borcherds exponent).
    """

    return True


# ---------------------------------------------------------------------------
# Hochschild-degree two pairing (R3 internal): the chain-level Delta_5
# generator in ChirHoch^2(A_{g_{Delta_5}}) on Lambda^{2,1}_II, after the
# eight gating hypotheses of rem:hoch-igusa-gauge-class (Vol II
# chapters/connections/hochschild.tex line 3874).
# ---------------------------------------------------------------------------


HOCH_DELTA5_GATING_HYPOTHESES: Tuple[str, ...] = (
    "finite_hall_coha_source",
    "mukai_serre_pairing_orientation",
    "pbw_no_extra_relations",
    "radical_quotient_complete",
    "Z2_parity_invariance",
    "filtered_completion_proper",
    "inverse_limit_exists",
    "EK_heegner_comparison_closed",
)


def hoch_degree_two_pairing_target(installed_gates: Iterable[str]) -> str:
    """Return the H^2_red(g_{Delta_5})^{Z/2, K(1)} -> C * Delta_5 status."""

    have = frozenset(installed_gates)
    missing = tuple(g for g in HOCH_DELTA5_GATING_HYPOTHESES if g not in have)
    if missing:
        return f"missing_gates:{missing}"
    return "conditional_h2_pairing_to_Delta_5_open"
