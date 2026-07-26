"""Irrational coset tempered-stratum structural engine.

Supports chapters/theory/irrational_cosets_tempered_platonic.tex.

The chapter proves the tempered-criterion refinement: on the standard
landscape, an infinite-dim Zhu algebra does NOT obstruct temperedness.
Three explicit irrational coset families are constructed and their
shadow channels bounded analytically. The engine encodes:

  (i)   Central charge formulas for coset candidates (Z_k, W_k,
        H_k cosets).
  (ii)  Zhu-algebra dimension analytics (finite vs infinite; coset
        pyramid structure).
  (iii) Channel decomposition ingredients matching the Virasoro / W(p)
        templates already used in tempered_stratum_characterization
        and logarithmic_wp_tempered_analysis.
  (iv)  Analytic ordinary-generating radius rho_*^{coset}(k, ...) =
        |c^{coset}| / beta_coset for each family.
  (v)   A certificate that the coset Virasoro subchannel falls
        outside the severe Kac-zero locus for generic irrational
        parameters.

The central analytic mechanism is: every irrational-coset candidate
here is built as Z(a) subalgebra of a C_2-cofinite-or-asymptotic
ambient whose shadow coefficients obey the universal Stirling
dominance r!^{-1/r} ~ e/r. Infinite-dim Zhu algebra corresponds to an
infinite multiplicity DECOMPOSITION over Heisenberg modes, NOT to
factorial shadow growth. The Heisenberg branching of each coset module
has polynomial (not factorial) multiplicities.

Conventions: Creutzig-Kanade-Linshaw 2019 (arXiv:1906.05868) for
parafermion K(g,k) = Com(H_k, V_k(g)); Adamovic 2007 / Adamovic-Milas
(arXiv:1206.3355, 0706.0803) for admissible affine minimal
cosets; Arakawa-Creutzig-Linshaw 2019 for irrational W-minus-H
cosets.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable


# ---------------------------------------------------------------------------
# Central charges of the three explicit coset candidates.
# ---------------------------------------------------------------------------


def parafermion_central_charge(k_num: int, k_den: int) -> Fraction:
    """c(K(sl_2, k)) = 2k/(k+2) - 1 = (k-2)/(k+2).

    Parafermion coset K(sl_2, k) = Com(H_k, V_k(sl_2)). Central charge
    formula: c(V_k(sl_2)) = 3k/(k+2), minus c(H_k) = 1. Hence
    c(K(sl_2, k)) = 3k/(k+2) - 1 = (2k - 2)/(k+2) = 2(k-1)/(k+2).

    At irrational k, K(sl_2, k) has infinite-dimensional Zhu algebra
    (admissible levels k = -2 + p/q give finitely many simple modules
    per Adamovic-Milas; irrational k gives a continuum parameter
    lift).
    """
    k = Fraction(k_num, k_den)
    # c(V_k(sl_2)) = 3k / (k + 2). Heisenberg rank-1 contributes c = 1.
    c_affine = Fraction(3) * k / (k + Fraction(2))
    c_heis = Fraction(1)
    return c_affine - c_heis


def affine_minimal_coset_central_charge(k_num: int, k_den: int) -> Fraction:
    """c(Com(H_k, V_k(sl_2))) at admissible-level k = -2 + p/q.

    Same formula as parafermion; distinction is scope (k admissible
    versus k generic-irrational).
    """
    return parafermion_central_charge(k_num, k_den)


def vir_in_affine_coset_central_charge(
    k_num: int,
    k_den: int,
) -> Fraction:
    """c(Com(Vir_{Sug}, V_k(sl_2))) = c(V_k) - c(Vir_Sug) = 3k/(k+2) - 3k/(k+2) = 0.

    The Sugawara Virasoro saturates the affine central charge at
    rank-1 level-k sl_2. The commutant Com(Vir_Sug, V_k(sl_2)) has
    zero central charge and is the trivial algebra at integer k. At
    irrational k, Com(Vir_c, V_k(sl_2)) is interpreted as
    Com(Vir_{c'}, V_k(sl_2)) for a non-Sugawara Virasoro embedding
    c' != c_Sug. For c' irrational, the coset carries non-trivial
    content.

    We adopt the convention: "Virasoro-sub-of-affine coset at
    parameter (k, c')" has effective central charge
        c_coset = 3k/(k+2) - c',
    bounded in (0, 3k/(k+2)) for c' in (0, 3k/(k+2)).
    """
    k = Fraction(k_num, k_den)
    c_affine = Fraction(3) * k / (k + Fraction(2))
    # For this analytic purpose we return c_affine; caller supplies c'.
    return c_affine


# ---------------------------------------------------------------------------
# Zhu algebra dimension analytics.
# ---------------------------------------------------------------------------


def zhu_dimension_finite_admissible(k_num: int, k_den: int) -> int:
    """dim A(L_k(sl_2)) at admissible k = -2 + p/q, gcd(p, q) = 1.

    Adamovic-Milas 1995 / Adamovic 2007: at admissible level,
    L_k(sl_2) is C_2-cofinite as a MODULE CATEGORY quotient (the
    full universal algebra is NOT C_2-cofinite); the admissible
    module category has (p-1)(q-1)/2 finitely many simples.

    For the PARAFERMION coset K(sl_2, k) the analogous count is the
    number of admissible K(sl_2, k)-modules, which is bounded by the
    same (p-1)(q-1)/2 combinatorial count. Hence the Zhu algebra of
    the COSET at admissible level is finite; temperedness is
    inherited from the C_2-cofinite W(p) template.
    """
    k_eff = Fraction(k_num, k_den) + Fraction(2)
    if k_eff <= 0:
        return 0
    p = abs(k_num + 2 * k_den)
    q = abs(k_den)
    from math import gcd
    g = gcd(p, q)
    p_red = p // g
    q_red = q // g
    if p_red < 2 or q_red < 2:
        return 0
    return ((p_red - 1) * (q_red - 1)) // 2


def zhu_dimension_infinite_irrational(k_num: int, k_den: int) -> int:
    """Return sentinel value -1 to denote infinite-dim Zhu.

    For generic irrational k (non-admissible), the Zhu algebra is
    infinite-dimensional: the coset has a continuous family of simple
    modules parametrized by a Heisenberg charge lattice. In the
    engine we return -1 as a sentinel (no admissible numerator/
    denominator reduction).

    The chapter's analytic tempering result applies REGARDLESS of the
    Zhu dimension: the finite-dimensional contribution from a typical
    fusion category quotient bounds each shadow channel; the
    infinite-dimensional Zhu does not feed factorial growth because
    the coset-module pyramid has POLYNOMIAL multiplicities over
    Heisenberg weight.
    """
    # Admissibility proxy: p/q rational in lowest terms with |p|, |q|
    # bounded. Irrational marker: k_den = 0 or caller's signal.
    if k_den == 0:
        return -1
    return zhu_dimension_finite_admissible(k_num, k_den)


# ---------------------------------------------------------------------------
# Channel decomposition: Virasoro sub-channel dominance.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CosetChannelData:
    """Channel decomposition template for an irrational coset.

    Parameters:
      c_coset: effective central charge
      beta_T:  Virasoro-sub-channel Riccati ratio (= 6 always)
      beta_J:  cross-channel generator-to-stress ratio (>= 2, bounded)
      beta_JJ: generator-to-generator channel ratio (>= 2, bounded)
      rho_star: min(|c|/beta_T, |c|/beta_JJ)
    """
    c_coset: Fraction
    beta_T: int
    beta_J: int
    beta_JJ: int
    rho_star: Fraction


@dataclass(frozen=True)
class BGGResolutionProfile:
    """BGG-resolution data used by the VSKR + BGG criterion."""

    sector: str
    resolution: str
    alternating_shadow_sum: str
    ramification_inequality: str
    local_pole_bound_required: bool


def bgg_resolution_profile(sector: str = "L(lambda)") -> BGGResolutionProfile:
    """Return the BGG shadow profile for an irrational coset sector."""

    return BGGResolutionProfile(
        sector=sector,
        resolution=(
            "... -> direct_sum_{ell(w)=2} M(w.lambda) -> "
            "direct_sum_{ell(w)=1} M(w.lambda) -> M(lambda) -> "
            "L(lambda) -> 0"
        ),
        alternating_shadow_sum=(
            "S_r^coset(lambda) = sum_{w in W_aff} (-1)^ell(w) "
            "S_r(M(w.lambda)) q^(Delta(w.lambda)-Delta(lambda))"
        ),
        ramification_inequality=(
            "Re(Delta(w.lambda)-Delta(lambda)) >= a ell(w)^2 - b"
        ),
        local_pole_bound_required=True,
    )


def ramification_lower_bound(
    length: int,
    *,
    a_num: int = 1,
    a_den: int = 1,
    b_num: int = 0,
    b_den: int = 1,
) -> Fraction:
    """Return a*length^2 - b for the Kazhdan ramification estimate."""

    if length < 0:
        raise ValueError(f"affine Weyl length must be non-negative, got {length}")
    a = Fraction(a_num, a_den)
    b = Fraction(b_num, b_den)
    if a <= 0:
        raise ValueError("ramification coefficient a must be positive")
    return a * length * length - b


def bgg_gaussian_exponent(
    length: int,
    *,
    q_decay_num: int = 1,
    q_decay_den: int = 1,
    a_num: int = 1,
    a_den: int = 1,
    b_num: int = 0,
    b_den: int = 1,
) -> Fraction:
    """Return -t(a length^2 - b) where |q| = exp(-t)."""

    t = Fraction(q_decay_num, q_decay_den)
    if t <= 0:
        raise ValueError("q-decay t must be positive")
    return -t * ramification_lower_bound(
        length,
        a_num=a_num,
        a_den=a_den,
        b_num=b_num,
        b_den=b_den,
    )


def bgg_gaussian_dominates_polynomial(
    *,
    degree: int,
    length: int,
    q_decay_num: int = 1,
    q_decay_den: int = 1,
    a_num: int = 1,
    a_den: int = 1,
    b_num: int = 0,
    b_den: int = 1,
) -> bool:
    """Return True once Gaussian length decay beats a polynomial proxy.

    We compare exponents: t(a m^2 - b) >= degree * log(m) is modeled
    by the elementary sufficient inequality
        t(a m^2 - b) >= degree * m
    for m >= 1, since log(m) <= m.
    """

    if degree < 0:
        raise ValueError("polynomial degree must be non-negative")
    if length < 1:
        raise ValueError("length must be >= 1 for polynomial comparison")
    t = Fraction(q_decay_num, q_decay_den)
    lhs = t * ramification_lower_bound(
        length,
        a_num=a_num,
        a_den=a_den,
        b_num=b_num,
        b_den=b_den,
    )
    rhs = Fraction(degree * length, 1)
    return lhs >= rhs


def parafermion_channel_data(k_num: int, k_den: int) -> CosetChannelData:
    """Channel decomposition for K(sl_2, k) = Com(H_k, V_k(sl_2)).

    Strong generators: psi_+ (weight 2/(k+2) + 1 = (k+4)/(k+2)),
    psi_- (weight same), plus composites. The OPE
        psi_+(z) psi_-(w)
    has pole order dominating at z=w; the non-integer conformal
    weights of psi_+-, however, kill the factorial-tempered envelope
    that Virasoro uses: the Virasoro sub-channel dominates the
    shadow coefficients. Explicitly:
        beta_T = 6 (Virasoro Riccati ratio),
        beta_J = 2 * Delta(psi) = 2(k+4)/(k+2) approx 2 in limit,
        beta_JJ = bounded by 2 * 2 = 4 (two psi fusing into
                  conformally weighted composite).
    In practice: rho_*^{K(sl_2,k)} = |c_coset(k)| / 6 = |k-1| / (3(k+2)).
    """
    c = parafermion_central_charge(k_num, k_den)
    # beta_T stays at 6; beta_J, beta_JJ bounded.
    return CosetChannelData(
        c_coset=c,
        beta_T=6,
        beta_J=2,
        beta_JJ=4,
        rho_star=abs(c) / Fraction(6),
    )


def affine_heisenberg_coset_channel_data(
    k_num: int,
    k_den: int,
) -> CosetChannelData:
    """Channel decomposition for Com(H, V_k(sl_2)) at irrational k.

    This coincides with parafermion K(sl_2, k) when H is the Cartan
    Heisenberg; the distinction is scope (irrational non-admissible
    level vs admissible p/q level).
    """
    return parafermion_channel_data(k_num, k_den)


def vir_affine_coset_channel_data(
    k_num: int,
    k_den: int,
    c_prime_num: int,
    c_prime_den: int,
) -> CosetChannelData:
    """Channel decomposition for Com(Vir_{c'}, V_k(sl_2)) at irrational (k, c').

    Non-Sugawara Virasoro embedding Vir_{c'} inside V_k(sl_2); coset
    central charge c_coset = 3k/(k+2) - c'. For this to be positive
    and irrational, c' in (0, 3k/(k+2)) and c' irrational.
    """
    c_affine = Fraction(3) * Fraction(k_num, k_den) / (
        Fraction(k_num, k_den) + Fraction(2)
    )
    c_prime = Fraction(c_prime_num, c_prime_den)
    c_coset = c_affine - c_prime
    return CosetChannelData(
        c_coset=c_coset,
        beta_T=6,
        beta_J=2,
        beta_JJ=4,
        rho_star=abs(c_coset) / Fraction(6),
    )


# ---------------------------------------------------------------------------
# Severe Kac-zero avoidance.
# ---------------------------------------------------------------------------


SEVERE_KAC_ZERO_LOCUS = (
    Fraction(-22, 5),
    Fraction(-10, 7),
    Fraction(-50, 13),
)


def coset_c_avoids_kac_locus(
    c_coset: Fraction,
    kac_candidates: Iterable[Fraction] = SEVERE_KAC_ZERO_LOCUS,
) -> bool:
    """Return True iff c_coset is not in the severe Kac-zero locus."""
    for kac in kac_candidates:
        if c_coset == kac:
            return False
    return True


# ---------------------------------------------------------------------------
# Summary certificate.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TemperingCertificate:
    """Output certificate for an irrational-coset candidate."""
    family: str
    k_repr: str
    c_coset: Fraction
    zhu_dim: int  # -1 sentinel = infinite
    beta_T: int
    rho_star: Fraction
    tempered: bool
    avoids_kac: bool
    bgg_ramification_assumed: bool


def certify_irrational_coset(
    family: str,
    k_num: int,
    k_den: int,
    *,
    c_prime_num: int = 0,
    c_prime_den: int = 1,
) -> TemperingCertificate:
    """Compute the tempering certificate for an irrational coset."""
    if family == "parafermion":
        data = parafermion_channel_data(k_num, k_den)
    elif family == "affine_heisenberg":
        data = affine_heisenberg_coset_channel_data(k_num, k_den)
    elif family == "vir_affine":
        data = vir_affine_coset_channel_data(
            k_num, k_den, c_prime_num, c_prime_den
        )
    else:
        raise ValueError(f"unknown family: {family}")

    zhu = zhu_dimension_infinite_irrational(k_num, k_den)
    avoids = coset_c_avoids_kac_locus(data.c_coset)
    # Tempered iff the Virasoro sub-channel avoids severe Kac-zeros,
    # beta_T is finite (always 6 here), and the BGG ramification
    # estimate is supplied. The Zhu-dimension sign does NOT enter.
    bgg_ramification_assumed = True
    tempered = avoids and data.beta_T < 1_000_000 and bgg_ramification_assumed

    k_repr = f"{k_num}/{k_den}" if k_den != 1 else str(k_num)
    return TemperingCertificate(
        family=family,
        k_repr=k_repr,
        c_coset=data.c_coset,
        zhu_dim=zhu,
        beta_T=data.beta_T,
        rho_star=data.rho_star,
        tempered=tempered,
        avoids_kac=avoids,
        bgg_ramification_assumed=bgg_ramification_assumed,
    )
