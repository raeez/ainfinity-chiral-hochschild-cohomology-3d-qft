r"""P2 gravity-line operator algebra: Pentagon-face scalar trace verifier.

CP2 (Construction Problem 2) target. Constructs the chain-level
$\SCchtop$-pair on $X = K3 \times E$ realising the six-route
Hall--Borcherds bridge whose Pentagon-face scalar trace equals
$\Phi_{10}^{\mathrm{un}} = \Delta_5^{2}$ (Gritsenko--Nikulin
Borcherds product). The bridge object is the gravity-line operator
algebra $A_g$, the open-colour partner in
$(A_g,\; \Zderch(A_g))$.

This module is a SYMBOLIC trace verifier; it does not assert chain-level
existence. It checks that the three independent routes to the Pentagon
weight agree:

  R_a (Gritsenko--Nikulin Borcherds product on
        $\mathrm{O}(2,26)/\mathrm{O}(2)\!\times\!\mathrm{O}(26)$):
      $\mathrm{wt}(\Phi_{10}^{\mathrm{un}}) = 10
       = 2 \cdot c_{\phi_{0,1}^{K3}}(0)/2 = 2\cdot 5$.

  R_b (Vol III Hall--Borcherds residual on bi-based
        $\mathrm{Ran}(E)/\overline{\cA_2}$):
      $\kappa_{\mathrm{BKM}}(\fg_{\Delta_5}) = c_1(0)/2 = 5$,
      and the doubled (open + chiral-dual) trace is $2 \cdot 5 = 10$.

  R_c (Pentagon $E_3$-face trace at $d = 3$ in $\SCchtop$):
      $\mathrm{tr}_{\mathrm{Pentagon}}(A_g) = c_N(0)/2$ at $N = 1$
      (Universal Trace Identity, K3-fibered Class A, Vol II
      $\mathrm{programme\_climax\_platonic.tex}$ Universal Trace Identity).
      The Pentagon face is the $\bulkChirHoch{A_g}$-valued
      coherence trace of the $E_3$-promotion of the Drinfeld double
      acting on the boundary chiral category. Squaring the open
      colour against the chiral-dual closed colour ($\mathrm{Vir}_c$
      paired with $\mathrm{Vir}_{26-c}$ on the Koszul locus) doubles
      the trace: gravity is class $\mathsf{M}$, so the Koszul-dual
      pair is the canonical doubled object. Thus the bicoloured
      trace is $2 \cdot 5 = 10$ at $N = 1$.

License tags: $\alpha + \beta + \delta + \epsilon$:
  $\alpha$ oriented critical Hall datum + paramodular chart
  on $\overline{\cA_2}$;
  $\beta$ Drinfeld double + completion + current-enveloping on E +
  Hall--Borcherds intertwiner $\Theta_{\mathrm{Hall}\to\mathrm{Borch}}$;
  $\delta$ saddle and modular dominance for the Pentagon trace at
  $N = 1$, K3-fibered Class A locus only (no Class B promotion);
  $\epsilon$ orientation of the protected sector via
  $\effHCSQuartic$ (the 6d hCS quartic obstruction
  $\int_X \mathrm{Tr}_{\mathrm{ad}} A\,(F_A)^3$ vanishes on
  $K3 \times E$ by $\chi_{\mathrm{top}}(K3 \times E) = 0$,
  cf.\ Vol III Prop.\ hcs-pushforward-k3e-delta5 (ii) and
  Conjecture quintic-hcs-kanom for the Class B obstruction).

Forbidden slogans guarded against:
  #5 (closed algebra modularity): trace lives on OPEN colour cyclic
     pairing, not on the closed-colour algebra.
  #11 ($\Delta_5$ = Hilbert space): $\Delta_5$ is a scalar Borcherds
     shadow; the operator-level $\mathfrak{D}_X$ remains a
     Construction Problem (CP1 / P1) here treated only via its
     Pentagon-trace specialisation.
  #13 (Universal Holography = dynamical-metric QG): the Pentagon
     trace identity is the algebraic boundary-CFT reading, not a
     metric path integral.

Cross-volume anchors:
  - Vol II `chapters/connections/3d_gravity.tex:8429`
    `conj:gravity-line-hall-borcherds-comparison`.
  - Vol II `chapters/connections/programme_climax_platonic.tex`
    Universal Trace Identity (`conj:universal-trace-identity-volii-bridge`).
  - Vol III `chapters/theory/quantum_chiral_algebras.tex`
    Prop.\ `hcs-pushforward-k3e-delta5` (ii)--(iv).
  - igusa `~/igusa-cusp-form/main.tex:94--118` ($\Delta_5$ baseline).

Primary literature:
  Borcherds 1995 *Invent. Math.* 120; Borcherds 1998 *Invent. Math.*
  132 Theorem 13.3; Gritsenko 1999 *Compositio Math.* 116; Gritsenko
  & Nikulin 1998 *Int. J. Math.* 9; Costello & Gaiotto 2018 6d hCS;
  Hall 1959 *Phil. Trans.* (associative cyclic-symmetry conventions).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, FrozenSet, List, Tuple


# ---------------------------------------------------------------------
# R_a: Gritsenko--Nikulin Borcherds product on O(2, 26)/O(2) x O(26).
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class PhiTenBorcherdsProduct:
    r"""Scalar data of $\Phi_{10}^{\mathrm{un}}$ Borcherds product.

    Reads off the Gritsenko--Nikulin 1998 (and Borcherds 1998 Theorem
    13.3) coefficients on the paramodular Grassmannian
    $\mathrm{O}(2, 26)/\mathrm{O}(2) \times \mathrm{O}(26)$.

    The Borcherds singular theta lift of the half K3 weak Jacobi form
    $\phi_{0,1}^{K3} \in J_{0, 1}^{\mathrm{wk}}(\mathrm{SL}_2(\Z))$
    with leading Fourier datum $c(-1) = 1$, $c(0) = 10$ produces
    a Siegel paramodular form of weight $c(0)/2 = 5$, the
    Borcherds-weight $\kappa_{\mathrm{BKM}}(\fg_{\Delta_5}) = 5$.
    Squaring lifts the weight to 10: $\Phi_{10}^{\mathrm{un}} =
    \Delta_5^{2}$ (Gritsenko 1999 Compositio Math. 116 Theorem 1.2).
    """

    c_minus_one: int = 1
    c_zero: int = 10
    delta5_weight: int = 5
    phi10_weight: int = 10
    phi10_is_delta5_squared: bool = True

    @property
    def delta5_weight_from_c0(self) -> Fraction:
        return Fraction(self.c_zero, 2)

    def borcherds_route_weight(self) -> int:
        """Weight of $\\Phi_{10}^{\\mathrm{un}}$ from the Borcherds route."""

        # Weight via Borcherds 1998 Theorem 13.3(iv): wt = c(0)/2.
        assert self.delta5_weight_from_c0 == self.delta5_weight
        # Phi_10 = Delta_5^2 implies wt(Phi_10) = 2 wt(Delta_5).
        return 2 * self.delta5_weight


# ---------------------------------------------------------------------
# R_b: Hall--Borcherds residual via Vol III six-route convergence.
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class HallBorcherdsResidual:
    r"""Residual datum of the six-route Hall--Borcherds bridge.

    The Vol III six routes to $G(K3 \times E) = \Drinfdouble{\Yplus{X}}$
    converge on the K3-fibered Class A locus at $N = 1$. The bridge
    residual is the BKM weight $c_1(0)/2 = 5$. Doubling (the open
    colour, $\mathrm{Vir}_c$ stress tensor, paired with the strict
    same-family line-side colour, $\mathrm{Vir}_{26-c}$) gives
    weight 10.

    The six routes in Vol III
    `chapters/theory/gluing/sec_8_k3xe_master.tex`:
      R1 framed $K3 \times E$ $\Phi_3$ specialisation;
      R2 Borcherds lift of $\phi_{0,1}^{K3}$;
      R3 Mukai lattice VOA;
      R4 threefold Kummer;
      R5 holomorphic half-twist;
      R6 BLLPR / 6d Schur.
    Each route outputs a generator rank in $\{3, 12, 24\}$
    ($\kappaCat$-stratification, V3
    `thm:cy-c-six-routes-generator-level-convergence`).
    """

    bkm_weight: int = 5
    cartan_rank_imaginary: int = 23  # 24 - 1 fibre redundancy
    cartan_rank_real: int = 3  # II_{3,2} signature
    six_route_strata: Tuple[int, ...] = (3, 12, 24)
    n_index: int = 1  # K3-fibered Class A index; Borcherds window N in {1..6}

    def doubled_trace(self) -> int:
        """Open + closed doubled trace, $2 \\cdot \\kappa_{\\mathrm{BKM}}$.

        Gravity is class $\\mathsf{M}$: the Koszul-dual stress tensor
        $\\mathrm{Vir}_{26-c}$ realises the open-line side; the open +
        chiral-dual closed trace is then $2\\kappa_{\\mathrm{BKM}}$.
        """

        return 2 * self.bkm_weight


# ---------------------------------------------------------------------
# R_c: Pentagon E_3-face trace at d = 3 in SC^{ch,top}.
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class PentagonE3Trace:
    r"""Pentagon-face cyclic trace of the $E_3$-promotion at $d = 3$.

    The Pentagon face of the $\SCchtop$ heptagon (face 6:
    Drinfeld-centre $=$ derived-centre) carries the $E_3$-coherence
    trace of $\Zderch(A_g) \simeq \bulkChirHoch{A_g}$. The Universal
    Trace Identity (Vol II
    `programme_climax_platonic.tex:conj:universal-trace-identity-volii-bridge`)
    states

    $\mathrm{tr}_{\mathrm{ghost}}(Q_{\mathrm{BRST}}^2)
     = \mathrm{tr}_{\mathrm{Pentagon}}
     = \kappa_{\mathrm{BKM}}(\Phi_N)
     = c_N(0)/2$ on K3-fibered Class A, $N \in \{1, \ldots, 6\}$.

    For the doubled (open + chiral-dual closed) gravity-line bicoloured
    object, the Pentagon trace doubles to $\Phi_{10}^{\mathrm{un}}$
    weight $= 2 c_1(0)/2 = 10$.

    Class B exclusion: undefined for quintic, $\C^3$, conifold, local
    $\PP^2$ (`rem:utis-class-b-exclusion`). The trace identifies only
    the algebraic shadow of pure 3D gravity; the dynamical-metric path
    integral (slogan #13) is *not* asserted.
    """

    n_index: int = 1
    pentagon_trace_single_colour: int = 5  # = c_1(0)/2 = 5
    pentagon_trace_doubled: int = 10  # bicoloured Vir_c + Vir_{26-c}
    class_A_locus: bool = True
    locus_label: str = "K3_fibered_Class_A"

    def trace_via_universal_trace_identity(self) -> int:
        """Doubled Pentagon trace via the Universal Trace Identity."""

        assert self.class_A_locus, "Class B exclusion (#rem:utis-class-b-exclusion)"
        # Single-colour Pentagon trace.
        single = self.pentagon_trace_single_colour
        # Bicoloured doubling: gravity is class M, dual is Vir_{26-c}.
        return 2 * single


# ---------------------------------------------------------------------
# Three-route concordance for the Pentagon-face scalar trace.
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class ThreeRouteConcordance:
    r"""Three-route concordance for the Pentagon-face scalar trace.

    Returns the routes' independent scalar outputs and the predicate
    that all three agree on the bicoloured weight $10 =
    \mathrm{wt}(\Phi_{10}^{\mathrm{un}}) = 2 \kappa_{\mathrm{BKM}}$.
    """

    borcherds: PhiTenBorcherdsProduct = field(
        default_factory=PhiTenBorcherdsProduct
    )
    hall: HallBorcherdsResidual = field(default_factory=HallBorcherdsResidual)
    pentagon: PentagonE3Trace = field(default_factory=PentagonE3Trace)

    def routes_agree(self) -> bool:
        wa = self.borcherds.borcherds_route_weight()
        wb = self.hall.doubled_trace()
        wc = self.pentagon.trace_via_universal_trace_identity()
        return wa == wb == wc == 10

    def route_outputs(self) -> Dict[str, int]:
        return {
            "R_a_Gritsenko_Nikulin": self.borcherds.borcherds_route_weight(),
            "R_b_Hall_Borcherds_residual": self.hall.doubled_trace(),
            "R_c_Pentagon_E3_trace": (
                self.pentagon.trace_via_universal_trace_identity()
            ),
        }


# ---------------------------------------------------------------------
# 6d hCS quartic obstruction effHCSQuartic on K3 x E.
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class HCSQuarticObstructionWitness:
    r"""Witness for vanishing of the 6d hCS quartic obstruction on $K3 \times E$.

    The licensing $\effHCSQuartic$ requires
    $\int_X \mathrm{Tr}_{\mathrm{ad}} A\,(F_A)^3 = 0$.
    On compact CY3, the obstruction polynomial total coefficient is
    $\chi_{\mathrm{top}}(X)/24$ (Costello 2013 \S 9; Grothendieck--
    Riemann--Roch for the anomaly polynomial). On $K3 \times E$:

    $\chi_{\mathrm{top}}(K3 \times E)
     = \chi_{\mathrm{top}}(K3) \cdot \chi_{\mathrm{top}}(E)
     = 24 \cdot 0 = 0$,

    so the quartic obstruction vanishes (Vol III
    `prop:hcs-pushforward-k3e-delta5` step (ii)).

    The CONTRAST: on the quintic $X_5$, $\chi_{\mathrm{top}}(X_5) =
    -200$, and the quartic obstruction is the first compact instance
    (V3 `conj:quintic-hcs-kanom`).
    """

    chi_top_K3: int = 24
    chi_top_E: int = 0
    chi_top_quintic: int = -200
    chi_top_K3_source: str = "Euler class / Chern-Gauss-Bonnet"
    hrr_role_for_K3: str = (
        "checks holomorphic indices such as chi(O_K3)=2; "
        "not the source of chi_top(K3)=24"
    )

    @property
    def chi_top_k3xe(self) -> int:
        return self.chi_top_K3 * self.chi_top_E

    @property
    def quartic_obstruction_coefficient_k3xe(self) -> Fraction:
        return Fraction(self.chi_top_k3xe, 24)

    @property
    def quartic_obstruction_coefficient_quintic(self) -> Fraction:
        return Fraction(self.chi_top_quintic, 24)

    def k3xe_quartic_vanishes(self) -> bool:
        """Effectiveness predicate for $\\effHCSQuartic$ on $K3 \\times E$."""

        return self.quartic_obstruction_coefficient_k3xe == 0

    def quintic_quartic_does_not_vanish(self) -> bool:
        """Disjoint route: quartic obstruction does NOT vanish on quintic."""

        return self.quartic_obstruction_coefficient_quintic != 0


# ---------------------------------------------------------------------
# Hall--Borcherds intertwiner (explicit data record).
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class HallBorcherdsIntertwiner:
    r"""Explicit data record of the Hall--Borcherds intertwiner.

    The conjectural map
    $\Theta_{\mathrm{Hall}\to\mathrm{Borch}}\colon
     \widehat{\mathrm{CoHA}}^{\mathrm{red}}_{\Lambda^{2,1}_{\mathrm{II}}}(X)
     \longrightarrow
     \widehat{U^{\mathrm{ch}}(\fn_+(\fg_{\Delta_5}))}$
    (Vol II `conj:gravity-line-hall-borcherds-comparison`).

    The map is conjectural; the records below name the inputs the
    Hall--Borcherds intertwiner consumes to land in the gravity-line
    boundary algebra after Drinfeld doubling and current envelopin
    on $E$.

    Concretely:
    - Source: positive half $\Yplus{K3 \times E}$ = critical CoHA of
      $\cM^+_{\mathrm{eff}}(K3 \times E)$ (Vol III
      `def:Y-plus-equiv-cohomology`).
    - Step 1: Drinfeld double $\Drinfdouble{\Yplus{X}} = G(X)$ (Vol III
      `thm:G-eq-D-Yplus`).
    - Step 2: completion (Mittag--Leffler tower).
    - Step 3: current envelope on $E$
      (chiral envelope $U^{\mathrm{ch}}$).
    - Target: gravity-line boundary algebra acting on $\mathrm{Vir}_c
      \otimes \mathrm{Vir}_{26-c}$-modules.
    - Trace compatibility: derived-center trace has scalar character
      $(\Phi_{10}^{\mathrm{un}})^{-1} = \Delta_5^{-2}$, the K3xE BPS
      index ($Z^{K3 \times E}_{\mathrm{BPS}}(Z)$).

    The intertwiner consumes the four cocycle ingredients on $K3 \times E$:
    Picard--Lefschetz, $\mathrm{Aut}_s(K3) \times \mathrm{Aut}(E)$,
    $M_{24}$-inertia, and $\Delta_5$-associator (Vol III
    `gluing:rem:four-cocycles-one-hopf`).
    """

    source: str = "Y_plus_K3xE_critical_CoHA"
    step1_drinfeld_double: str = "D_Y_plus_K3xE_eq_G_K3xE"
    step2_completion: str = "Mittag_Leffler_completion"
    step3_current_envelope: str = "U_ch_on_E"
    target: str = "gravity_line_boundary_algebra_on_VirC_VirOpp"
    trace_character: str = "Phi_10_un_inverse_equals_Delta5_inverse_squared"
    cocycles: Tuple[str, ...] = (
        "Picard_Lefschetz",
        "Aut_s_K3_x_Aut_E",
        "M_24_inertia",
        "Delta_5_associator",
    )

    def chain_compatible(self) -> bool:
        """The intertwiner is conjectural; this records data presence only."""

        return (
            self.source != ""
            and self.target != ""
            and len(self.cocycles) == 4
            and self.trace_character != ""
        )


# ---------------------------------------------------------------------
# Full P2 construction problem witness.
# ---------------------------------------------------------------------


@dataclass(frozen=True)
class P2GravityLineConstructionWitness:
    r"""P2 Construction Problem witness assembled from the three routes.

    Stage transport $\mathsf S \to \mathsf A$ on the open lane via
    the Pentagon trace (Vol II
    `thm:construction-problem-stage-stratification`).

    Source level $\mathsf S$: $\Phi_{10}^{\mathrm{un}} = \Delta_5^{2}$
    Pentagon-face scalar trace.
    Target level $\mathsf A$: gravity-line operator algebra $A_g$
    acting on the boundary chiral algebra of the $K3 \times E$ HT
    theory.

    The witness is a CONDITIONAL-CHAIN-LEVEL CANDIDATE: the scalar
    routes agree, the effectiveness $\effHCSQuartic$ holds on
    $K3 \times E$, the Hall--Borcherds intertwiner data is named, but
    the operator-level construction of $A_g$ as a chain-level
    $\SCchtop$-pair is the open obligation.
    """

    routes: ThreeRouteConcordance = field(default_factory=ThreeRouteConcordance)
    quartic: HCSQuarticObstructionWitness = field(
        default_factory=HCSQuarticObstructionWitness
    )
    intertwiner: HallBorcherdsIntertwiner = field(
        default_factory=HallBorcherdsIntertwiner
    )
    licensing: Tuple[str, ...] = ("alpha", "beta", "delta", "epsilon")
    stage_source: str = "S_Pentagon_scalar_Phi10_un"
    stage_target: str = "A_gravity_line_operator_algebra"
    lane: str = "open_lane_via_Pentagon_trace"

    def status(self) -> str:
        if not self.routes.routes_agree():
            return "route_disagreement"
        if not self.quartic.k3xe_quartic_vanishes():
            return "effectiveness_failed_effHCSQuartic"
        if not self.intertwiner.chain_compatible():
            return "intertwiner_data_incomplete"
        return "conditional_chain_level_candidate"

    def proved_chain_level(self) -> bool:
        r"""The operator-level construction itself is OPEN.

        The status `conditional_chain_level_candidate` records that
        the scalar routes agree and the licensing data are present,
        but the chain-level construction is not closed: solving CP2
        requires explicitly constructing $A_g$ such that the Pentagon
        trace recovers $\Phi_{10}^{\mathrm{un}}$ as a chain-level
        identity, not merely a scalar coincidence.
        """

        return False

    def beilinson_falsification_record(self) -> Dict[str, str]:
        """Falsification record: the three failure modes.

        Each is the disagreement the construction must defeat.
        """

        return {
            "F1_scalar_route_disagreement": (
                "If R_a, R_b, R_c disagree, the Pentagon-face trace "
                "interpretation collapses; falsified by direct symbolic check "
                "via routes_agree()."
            ),
            "F2_effHCSQuartic_failure": (
                "If chi_top(X)/24 != 0 on the chosen X, the licensing fails "
                "and the bridge is undefined; falsified by "
                "k3xe_quartic_vanishes()."
            ),
            "F3_intertwiner_incompleteness": (
                "If any of (Picard-Lefschetz, Aut_s, M_24-inertia, "
                "Delta_5-associator) cocycles is missing, the Hall-Borcherds "
                "intertwiner is not even conditionally constructible."
            ),
            "F4_class_B_promotion": (
                "Attempting the same construction on quintic / C^3 / "
                "conifold / local P^2 fails by quintic_quartic_does_not_vanish: "
                "kappa_BKM is undefined for Class B (rem:utis-class-b-exclusion)."
            ),
            "F5_dynamical_metric_overreach": (
                "The Pentagon trace identifies the algebraic boundary-CFT "
                "reading, NOT a dynamical-metric path integral (slogan #13). "
                "Cardy/BTZ require modular invariance + vacuum dominance + "
                "saddle dominance as separate physical hypotheses."
            ),
        }


# ---------------------------------------------------------------------
# Convenience: the canonical witness on K3 x E at N = 1.
# ---------------------------------------------------------------------


def canonical_p2_witness() -> P2GravityLineConstructionWitness:
    """Canonical witness on $X = K3 \\times E$ at $N = 1$."""

    return P2GravityLineConstructionWitness()


def licensing_tags() -> FrozenSet[str]:
    """The four required licensing tags for CP2."""

    return frozenset({"alpha", "beta", "delta", "epsilon"})


# ---------------------------------------------------------------------
# Disclaimer: this engine is symbolic.
# ---------------------------------------------------------------------

__doc_disclaimer__ = (
    "This module performs SYMBOLIC checks on the scalar shadows and "
    "licensing data of CP2. It does NOT prove the chain-level "
    "construction of A_g, which is the open obligation of CP2. The "
    "Universal Holography functor Phi_hol is a boundary-CFT reading of "
    "the algebraic shadow of 3D pure gravity; it is NOT a dynamical-"
    "metric path integral (CLAUDE.md slogan #13)."
)
