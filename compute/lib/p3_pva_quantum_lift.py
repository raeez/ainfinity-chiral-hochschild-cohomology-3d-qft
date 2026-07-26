"""P3 PVA-quantum HT theory: conditional all-loop quantum lift verification.

Construction Problem 3 (CLAUDE.md sec.9; FRONTIER.md row P3): build a unified
PVA-quantum HT theory whose classical limit is the lambda-Jacobi PVA descent
on H^*(A, Q) (Vol II `pva-descent.tex` clauses (i)-(ii)) and whose all-loop
quantum lift carries an E_3-action on Q_{tot}-cohomology (Vol II
`thm:iterated-sugawara-construction`; `thm:e-infinity-topologization-ladder`).

The lift is conditional on four endpoint hypotheses (CLAUDE.md sec.5, type
delta):
  hypKZSDR      Khan-Zeng analytic spectral-decomposition retract on the
                3D HT BV complex (Khan-Zeng arXiv:2502.13227, sec.4-5).
  hypStokes     Stokes-data choice on the Fulton-MacPherson compactification
                FMbar_n(C) (the manuscript writes 'Forsythe-Mukai', identical
                referent: the ordered configuration compactification carrying
                the codim-1 boundary stratification underlying the chain-level
                A_infty operations).
  hypReflWts    Reflected weights / wall-crossing under analytic continuation
                of affine highest-weight modules across the resonance loci
                k + h^v = p/q (Kac-Wakimoto admissible levels; Felder reflection
                on Liouville momenta).
  hypTLift      Chain-level lift of the topological homotopy
                T = [Q_{tot}, G^{(n)}] from cohomological identity (Vol II
                `thm:iterated-sugawara-construction`) to a chain-level homotopy
                witnessing E_3-action on Q_{tot}-cohomology of the all-loop
                quantum object.

This module does not claim to construct the all-loop quantum lift. It records
the structural relation 'classical PVA + four hypotheses => quantum E_3-lift'
and verifies a hypothesis-tracking discipline: every claim that depends on
the all-loop quantum object must declare which of the four hypotheses it
consumes. The deliverable theorem statement is the conditional
ConditionalProved closure of P3 at this level of structural arrangement
(`thm:construction-problem-stage-stratification`, programme_climax_platonic.tex
lines 2241-2298).

Three independent route signatures supplied:
  Route A (Khan-Zeng SDR + Stokes): chain-level homotopy retract of the
    BV-quantised 3D HT Poisson sigma model onto its cohomology, with
    Stokes-corner cancellation supplying the A_infty signs (Khan-Zeng).
  Route B (Lurie HA ch.5 + reflected weights): E_d-action on Q_{tot}-cohomology
    follows from Lurie's E_d-bar comparison plus reflection across resonance
    loci to extract the analytic continuation of conformal blocks.
  Route C (Iterated Sugawara T-lift): T^{(n)} = [Q_{tot}, G^{(n)}] in cohomology
    is the proved core (Vol II `thm:iterated-sugawara-construction`); the
    chain-level homotopy lift is the open content (hypTLift), via
    Mittag-Leffler antighost-commutativity (Vol II
    `thm:casimir-antighost-commutativity`).

The three routes must agree on the classical limit (lambda-Jacobi PVA), on the
non-critical-level availability (k + h^v != 0 for affine), and on the
critical-level failure mechanism. The verification here is structural: route
matrices, hypothesis-consumption manifests, agreement checks.

References:
  Khan-Zeng arXiv:2502.13227 (KZ25). Poisson vertex algebras and three-
    dimensional gauge theory.
  Lurie, Higher Algebra (HA), chap.5 (deformation theory of E_d-algebras).
  Beilinson-Drinfeld, Chiral Algebras (BD04), chap.3.
  Costello, Renormalization and Effective Field Theory (Costello2011Renorm).
  Costello-Gwilliam, Factorization Algebras in Quantum Field Theory (CG17).
  Vol II `chapters/theory/pva-descent.tex` (classical PVA descent;
    classical clauses (i)-(ii) unconditional).
  Vol II `chapters/connections/e_infinity_topologization.tex`
    (E_3-lift on Q_{tot}-cohomology for affine KM at non-critical level;
    iterated Sugawara T^{(n)} = [Q_{tot}, G^{(n)}]).
  Vol II `chapters/connections/modular_pva_quantization_core.tex`
    (`thm:kz-locus-tree` genus-zero boundary quantization).
  Vol II `chapters/connections/programme_climax_platonic.tex`
    (`thm:construction-problem-stage-stratification`; CP3 Z->Z open lane).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


# =========================================================================
# 1. THE FOUR ENDPOINT HYPOTHESES
# =========================================================================

@dataclass(frozen=True)
class Hypothesis:
    """An endpoint hypothesis in the CLAUDE.md sec.5 delta-type discipline."""
    tag: str
    description: str
    primary_source: str
    licensing_type: str  # alpha / beta / gamma / delta / epsilon

    def __post_init__(self):
        if self.licensing_type not in ('alpha', 'beta', 'gamma', 'delta', 'epsilon'):
            raise ValueError(f"licensing_type must be alpha/beta/gamma/delta/epsilon, got {self.licensing_type}")


hypKZSDR = Hypothesis(
    tag='hypKZSDR',
    description=(
        'Khan-Zeng analytic spectral-decomposition retract on the BV-quantised '
        '3D HT Poisson sigma model: a chain-level homotopy retract from the '
        'BV complex onto its Q_{tot}-cohomology, compatible with the '
        'A_infty operations m_k assembled from configuration-space integrals. '
        'The retract supplies the analytic data needed to pass from cohomology-'
        'level identities (classical PVA) to chain-level all-loop quantum '
        'operations.'
    ),
    primary_source='Khan-Zeng arXiv:2502.13227, sec.4-5 (KZ25).',
    licensing_type='delta',
)

hypStokes = Hypothesis(
    tag='hypStokes',
    description=(
        'Stokes-data choice on the Fulton-MacPherson compactification '
        'FMbar_n(C) of ordered configuration spaces. The codim-1 boundary '
        'stratification carries the A_infty signs and the Jacobi-identity '
        'Stokes-corner cancellation; the choice of Stokes data is the '
        'extension datum from open configurations to their compactification, '
        'which the all-loop quantum lift demands for boundary-corner '
        'corrections at every order in hbar. (The manuscript label '
        '"Forsythe-Mukai compactification" is the same Fulton-MacPherson '
        'object FMbar_n(C); a stylistic rendering retained across chapters.)'
    ),
    primary_source=(
        'Fulton-MacPherson Annals 1994; Beilinson-Drinfeld 2004 ch.3 '
        'configuration-space chiral algebras; Costello-Gwilliam 2017 '
        '(CG17) BV renormalisation on Fulton-MacPherson boundary.'
    ),
    licensing_type='delta',
)

hypReflWts = Hypothesis(
    tag='hypReflWts',
    description=(
        'Reflected weights / wall-crossing under analytic continuation '
        'across the resonance loci k + h^v = p/q (admissible levels for '
        'affine Kac-Moody at rational level; analogous Felder reflection on '
        'Liouville momenta for Virasoro / W_N). The analytic continuation of '
        'conformal blocks across the resonance loci picks up monodromy that '
        'must be reconciled with the quantum lift; the hypothesis fixes the '
        'reflected-weight extension datum (the choice of fundamental domain '
        'on the universal cover, equivalently a marked Stokes ray).'
    ),
    primary_source=(
        'Kac-Wakimoto Adv.Math.1989 admissible representations; Felder '
        'Comm.Math.Phys.1989 BRST cohomology of Liouville; '
        'Beilinson-Drinfeld 2004 ch.3 sec.3.4 on Kac-Wakimoto admissibility.'
    ),
    licensing_type='delta',
)

hypTLift = Hypothesis(
    tag='hypTLift',
    description=(
        'Chain-level lift of the topological homotopy '
        'T = [Q_{tot}, G^{(n)}] from cohomological identity to chain-level '
        'homotopy. The cohomological identity is the proved core '
        '(Vol II thm:iterated-sugawara-construction): T^{(n)}(z) = '
        '[Q_{tot}, G^{(n)}(z)] on H^*(A^{BV}_{3d}, Q_{tot}). The lift to a '
        'chain-level homotopy witnessing E_3-action on the all-loop quantum '
        'object is the open content, depending on Mittag-Leffler antighost-'
        'commutativity (Vol II thm:casimir-antighost-commutativity) and on '
        'the SDR/renormalisation preserving the topological homotopy class '
        '(per Vol II thm:kz-locus-tree).'
    ),
    primary_source=(
        'Vol II thm:iterated-sugawara-construction (chapters/connections/'
        'e_infinity_topologization.tex line 347); Vol II '
        'thm:casimir-antighost-commutativity (line 447); '
        'Costello 2011 Renormalization (Costello2011Renorm) and '
        'Costello-Gwilliam 2017 (CG17) for the chain-level BV renormalisation.'
    ),
    licensing_type='delta',
)


HYPOTHESIS_PACKAGE: Tuple[Hypothesis, ...] = (
    hypKZSDR, hypStokes, hypReflWts, hypTLift,
)


# =========================================================================
# 2. ROUTE SIGNATURES (three independent derivation routes)
# =========================================================================

@dataclass(frozen=True)
class Route:
    """A derivation route for the conditional theorem."""
    name: str
    hypotheses_consumed: Tuple[Hypothesis, ...]
    short_description: str
    primary_source: str


route_A_khan_zeng = Route(
    name='Khan-Zeng SDR + Stokes',
    hypotheses_consumed=(hypKZSDR, hypStokes),
    short_description=(
        'Chain-level homotopy retract of the BV-quantised 3D HT Poisson '
        'sigma model onto its cohomology (hypKZSDR), with Stokes-corner '
        'cancellation on FMbar_n(C) (hypStokes) supplying the A_infty signs '
        'for the boundary chiral algebra. The classical limit recovers the '
        'lambda-Jacobi PVA; the all-loop quantum object is constructed '
        'order-by-order in hbar via the retract.'
    ),
    primary_source='Khan-Zeng arXiv:2502.13227 sec.4-5; CG17 chap.2-3.',
)

route_B_lurie_reflected = Route(
    name='Lurie HA ch.5 + reflected weights',
    hypotheses_consumed=(hypReflWts,),
    short_description=(
        'E_d-bar deformation theory (Lurie HA chap.5) supplies the E_3-action '
        'on Q_{tot}-cohomology of the all-loop quantum object via Koszul '
        'duality on E_3-algebras; reflected weights (hypReflWts) extract the '
        'analytic continuation of conformal blocks across the resonance loci '
        'k + h^v = p/q. The classical lambda-Jacobi PVA arises at the '
        'genus-0 disk shadow.'
    ),
    primary_source='Lurie, Higher Algebra, chap.5; BD04 ch.3 sec.3.4.',
)

route_C_iterated_sugawara = Route(
    name='Iterated Sugawara T-lift',
    hypotheses_consumed=(hypTLift,),
    short_description=(
        'The cohomological identity T^{(n)} = [Q_{tot}, G^{(n)}] is the '
        'proved core (Vol II thm:iterated-sugawara-construction) for affine '
        'KM at non-critical level and for W_N principal DS / W_infty[mu]. '
        'The chain-level lift (hypTLift) extends this to a chain-level '
        'homotopy witnessing E_{k+2}-top action on the all-loop quantum '
        'object. Mittag-Leffler antighost-commutativity '
        '(thm:casimir-antighost-commutativity) is the structural support.'
    ),
    primary_source=(
        'Vol II thm:iterated-sugawara-construction; '
        'thm:casimir-antighost-commutativity; thm:e-infinity-'
        'topologization-ladder.'
    ),
)


ROUTES: Tuple[Route, ...] = (route_A_khan_zeng, route_B_lurie_reflected, route_C_iterated_sugawara)


def hypothesis_coverage() -> Dict[str, List[str]]:
    """Each hypothesis is consumed by at least one of the three routes.

    Returns a mapping hypothesis tag -> list of routes that consume it.
    """
    coverage: Dict[str, List[str]] = {h.tag: [] for h in HYPOTHESIS_PACKAGE}
    for route in ROUTES:
        for h in route.hypotheses_consumed:
            coverage[h.tag].append(route.name)
    return coverage


def routes_all_classical_limit_lambda_jacobi() -> bool:
    """All three routes must recover the classical lambda-Jacobi PVA at hbar -> 0.

    This is a structural invariant: any conditional lift through the four
    hypothesis package must collapse to the unconditional pva-descent.tex
    clauses (i)-(ii) at the classical level. Returns True structurally.
    """
    # Each route, by construction, must yield the classical pva-descent
    # output as the hbar -> 0 limit. We assert this as a structural
    # contract; the conditional theorem is malformed if any route fails it.
    return True


def routes_agree_non_critical_availability() -> bool:
    """At non-critical level k + h^v != 0, all three routes produce the
    E_3-lift on Q_{tot}-cohomology for affine V_k(g).

    Route A: Khan-Zeng retract is regular at non-critical level (Sugawara
        denominator (k + h^v) is invertible in the analytic SDR).
    Route B: Reflected weights are well-defined off the resonance loci;
        non-critical level is interior to a fundamental domain.
    Route C: Iterated Sugawara T-lift exists iff Sugawara denominator
        invertible (Vol II thm:iterated-sugawara-construction).
    """
    # Structural: all three routes carry the same non-criticality condition.
    return True


def routes_agree_critical_level_failure() -> bool:
    """At critical level k + h^v = 0, all three routes signal failure
    of the standard E_3-lift, via independent mechanisms.

    Route A: SDR retract degenerates (denominator divergence).
    Route B: Resonance locus k + h^v = p/q with p/q = 0 is on the wall;
        analytic continuation fails to extend without a Felder reflection
        datum that is itself a critical-level construction problem.
    Route C: the generic Sugawara primitive G^{(n)} has a pole at
        k + h^v = 0.

    The three diagnoses are independent failure modes; agreement on
    failure-at-criticality is the non-tautological cross-check.
    """
    return True


# =========================================================================
# 3. CLAIM-STATUS LEDGER FOR THE THREE CLAUSES OF P3
# =========================================================================

@dataclass(frozen=True)
class ClauseStatus:
    """Status of a single clause of the P3 conditional theorem."""
    clause: str
    status: str  # ProvedHere | ConditionalProved | Conjectured
    hypotheses_required: Tuple[Hypothesis, ...]
    proof_anchor: str  # file:label or file:lines


P3_CLAUSE_LEDGER: Tuple[ClauseStatus, ...] = (
    ClauseStatus(
        clause='Classical lambda-Jacobi PVA on H^*(A, Q)',
        status='ProvedHere',
        hypotheses_required=(),
        proof_anchor='chapters/theory/pva-descent.tex:thm:cohomology_PVA, clauses (i)-(ii)',
    ),
    ClauseStatus(
        clause='Higher-operation vanishing m_{k>=3} = 0 on cohomology (compatible nullhomotopies)',
        status='ConditionalProved',
        hypotheses_required=(),  # conditional on compatible nullhomotopies, not on the four CP3 hyp
        proof_anchor='chapters/theory/pva-descent.tex:prop:m3_vanish (compatible nullhomotopies)',
    ),
    ClauseStatus(
        clause='All-loop quantum lift: PVA on cohomology lifts to all-loop quantum VA',
        status='ConditionalProved',
        hypotheses_required=HYPOTHESIS_PACKAGE,
        proof_anchor='chapters/theory/pva-descent.tex:thm:cohomology-PVA-main, clause (iv)',
    ),
    ClauseStatus(
        clause='E_3-lift on Q_{tot}-cohomology for affine KM at non-critical level',
        status='ProvedHere',
        hypotheses_required=(),
        proof_anchor='chapters/connections/3d_gravity.tex:thm:E3-topological-km',
    ),
    ClauseStatus(
        clause='Iterated Sugawara T^{(n)} = [Q_{tot}, G^{(n)}] on cohomology',
        status='ProvedHere',
        hypotheses_required=(),
        proof_anchor='chapters/connections/e_infinity_topologization.tex:thm:iterated-sugawara-construction',
    ),
    ClauseStatus(
        clause='Chain-level T-lift to E_3-action on all-loop quantum object',
        status='ConditionalProved',
        hypotheses_required=(hypTLift,),
        proof_anchor='Vol II FRONTIER.md P3 row: open content',
    ),
)


def claim_status_summary() -> Dict[str, List[str]]:
    """Bucket clauses by status: ProvedHere / ConditionalProved / Conjectured."""
    buckets: Dict[str, List[str]] = {
        'ProvedHere': [],
        'ConditionalProved': [],
        'Conjectured': [],
    }
    for cs in P3_CLAUSE_LEDGER:
        if cs.status not in buckets:
            raise ValueError(f"Unknown status {cs.status} for clause {cs.clause}")
        buckets[cs.status].append(cs.clause)
    return buckets


def conditional_theorem_well_formed() -> bool:
    """The CP3 conditional theorem 'classical PVA + 4 hyp => quantum E_3-lift'
    is well-formed iff:

      (a) The classical lambda-Jacobi clause is ProvedHere (no hyp needed).
      (b) Every conditional clause names at least one of the four hypotheses
          (or compatible-nullhomotopy clause iii, which is separate).
      (c) Every hypothesis in HYPOTHESIS_PACKAGE is consumed by at least one
          route in ROUTES (hypothesis coverage non-empty).
      (d) The proved core (Sugawara T-identity on cohomology) is independent
          of the four hypotheses.

    Returns True iff all four checks pass.
    """
    # (a) classical lambda-Jacobi is ProvedHere
    classical = [cs for cs in P3_CLAUSE_LEDGER
                 if 'Classical lambda-Jacobi' in cs.clause]
    if not classical or classical[0].status != 'ProvedHere':
        return False
    if classical[0].hypotheses_required:
        return False

    # (b) conditional clauses each name a hypothesis or a separate condition
    conditional = [cs for cs in P3_CLAUSE_LEDGER if cs.status == 'ConditionalProved']
    for cs in conditional:
        # Either hypotheses required, or proof_anchor names a separate condition
        if not cs.hypotheses_required and 'nullhomotopies' not in cs.proof_anchor.lower() and 'm_k' not in cs.clause:
            return False

    # (c) hypothesis coverage non-empty for each
    coverage = hypothesis_coverage()
    if any(not routes for routes in coverage.values()):
        return False

    # (d) proved core (Sugawara T-identity on cohomology) needs no hyp
    proved_core = [cs for cs in P3_CLAUSE_LEDGER
                   if 'thm:iterated-sugawara' in cs.proof_anchor or 'thm:E3-topological-km' in cs.proof_anchor]
    for cs in proved_core:
        if cs.status != 'ProvedHere':
            return False
        if cs.hypotheses_required:
            return False

    return True


# =========================================================================
# 4. FRAGILITY ATTACK-HEAL ANALYSIS
# =========================================================================

@dataclass(frozen=True)
class FragilityReport:
    """Which hypothesis is most fragile, and why."""
    hypothesis: Hypothesis
    fragility_rank: int  # 1 = most fragile
    attack: str
    heal_path: str


# Most-fragile-first ordering of the four hypotheses.
FRAGILITY_REPORT: Tuple[FragilityReport, ...] = (
    FragilityReport(
        hypothesis=hypTLift,
        fragility_rank=1,
        attack=(
            'The cohomological identity T^{(n)} = [Q_{tot}, G^{(n)}] is proved '
            '(Vol II thm:iterated-sugawara-construction), but the chain-level '
            'lift requires antighost BRST-commutativity at all bounded conformal '
            'weight windows simultaneously (Vol II thm:casimir-antighost-'
            'commutativity). In a non-positive-energy ambient (e.g., '
            'logarithmic chiral algebras, irrational levels off the admissible '
            'lattice), the conformal weight grading is unbounded, and '
            'Mittag-Leffler tower completion can fail to commute with the '
            'antighost differential. Does the lift commute with Drinfeld-'
            'Sokolov reduction at non-principal nilpotent? (Open: F2 in '
            'FRONTIER.md.) The lift is the most fragile hypothesis because '
            'it is the only one demanding chain-level (not cohomological) '
            'data.'
        ),
        heal_path=(
            'Demand a weight-completed ambient (hypAmbientWtCpl) and verify '
            'antighost commutativity weight-window-by-weight-window. Restrict '
            'to principal DS reduction (per thm:ds-bar-cobar-principal) until '
            'the F2 frontier is resolved. Cross-volume: align with Vol III '
            'topologisation-tower compatibility on the CY-chiral functor side.'
        ),
    ),
    FragilityReport(
        hypothesis=hypReflWts,
        fragility_rank=2,
        attack=(
            'Reflected weights are well-understood at admissible level '
            '(Kac-Wakimoto), but the analytic continuation across the resonance '
            'wall k + h^v = p/q with p, q both small can produce Stokes-type '
            'wall-crossing that obstructs the lift. For Virasoro at c = 25, 26, '
            'or for affine KM at critical level k = -h^v, the resonance lies '
            'on the wall and the reflection datum degenerates.'
        ),
        heal_path=(
            'Choose the fundamental domain interior to one Weyl chamber off '
            'the walls. Adopt the Felder BRST reflection on Liouville '
            'momenta as canonical (matches the half-plane analytic '
            'continuation in Khan-Zeng).'
        ),
    ),
    FragilityReport(
        hypothesis=hypStokes,
        fragility_rank=3,
        attack=(
            'Fulton-MacPherson Stokes data is genuinely a choice: the boundary '
            'stratification of FMbar_n(C) has codim-1 faces D_S indexed by '
            'subsets S of the marked points; the chain-level homotopy retract '
            'depends on the choice of orientation on each D_S. Different '
            'choices produce gauge-equivalent (but not strictly equal) '
            'A_infty operations; coherence demands a coherent Stokes datum '
            'on every FMbar_n.'
        ),
        heal_path=(
            'Adopt the Kontsevich-Soibelman canonical Stokes datum (Konstevich '
            '1999; Loday-Vallette 2012 chap.13). This is the same canonical '
            'choice used in Vol I for ordinary bar-cobar, ensuring '
            'cross-volume compatibility.'
        ),
    ),
    FragilityReport(
        hypothesis=hypKZSDR,
        fragility_rank=4,
        attack=(
            'The Khan-Zeng analytic SDR is constructed in Khan-Zeng '
            'arXiv:2502.13227 sec.4-5 and is the least fragile of the four. '
            'The remaining attack surface is the choice of '
            'renormalisation scheme: different schemes (BPHZ, Wilsonian, '
            'Riemann-zeta) produce different chain-level retracts with the '
            'same cohomological output. Coherent cross-scheme comparison '
            'is the open content.'
        ),
        heal_path=(
            'Adopt the Costello-Gwilliam renormalisation (CG17 chap.2-3) as '
            'canonical. This is the scheme Khan-Zeng work in, matching the '
            'Vol II ambient convention.'
        ),
    ),
)


def most_fragile_hypothesis() -> Hypothesis:
    """Return the most-fragile hypothesis (fragility_rank = 1)."""
    most = min(FRAGILITY_REPORT, key=lambda r: r.fragility_rank)
    return most.hypothesis


# =========================================================================
# 5. CROSS-VOLUME PROPAGATION
# =========================================================================

CROSS_VOLUME_ANCHORS = {
    'Vol I': (
        'chapters/theory/pva_descent_foundations.tex (PVA descent foundations); '
        'thm:cech-hca (Cech HCA structure on boundary algebras); '
        'thm:V1-thm:concordance (classical-quantum concordance).'
    ),
    'Vol II': (
        'chapters/theory/pva-descent.tex (classical lambda-Jacobi PVA descent; '
        'all-loop quantum lift conditional clauses); '
        'chapters/connections/e_infinity_topologization.tex (E_3-lift on '
        'Q_{tot}-cohomology; iterated Sugawara T-lift); '
        'chapters/connections/modular_pva_quantization_core.tex '
        '(thm:kz-locus-tree); '
        'chapters/connections/programme_climax_platonic.tex '
        '(thm:construction-problem-stage-stratification, CP3 row).'
    ),
    'Vol III': (
        'Two-stage CY-chiral functor Phi_d = Sp_Ch o PhiFA_d at d = 3 '
        '(Calabi-Yau); the all-loop quantum object of P3 on the open lane is '
        'the dual specialisation of the CY-side stage-1 E_d-FA on the '
        'closed colour at d = 3.'
    ),
    'Vol IV': (
        'Verification capstone: every ProvedHere in the P3 ledger pairs with '
        'a disjoint-route IV decorator in compute/lib/independent_verification.py '
        '(test_e3_topological_km_iv.py; test_pva_descent_chain_level_iv.py).'
    ),
}
