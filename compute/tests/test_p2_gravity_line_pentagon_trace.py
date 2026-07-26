"""Tests for P2 gravity-line Pentagon-trace verifier (CP2)."""

from fractions import Fraction
import os
import sys


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.p2_gravity_line_pentagon_trace import (  # noqa: E402
    HCSQuarticObstructionWitness,
    HallBorcherdsIntertwiner,
    HallBorcherdsResidual,
    P2GravityLineConstructionWitness,
    PentagonE3Trace,
    PhiTenBorcherdsProduct,
    ThreeRouteConcordance,
    canonical_p2_witness,
    licensing_tags,
)


# ---------------------------------------------------------------------
# R_a: Gritsenko--Nikulin Borcherds product.
# ---------------------------------------------------------------------


def test_borcherds_route_recovers_phi10_weight_10():
    b = PhiTenBorcherdsProduct()
    # Borcherds 1998 Theorem 13.3(iv): wt(Phi) = c(0)/2.
    assert b.delta5_weight_from_c0 == Fraction(5)
    assert b.delta5_weight == 5
    # Phi_10 = Delta_5^2 doubles the weight.
    assert b.borcherds_route_weight() == 10
    assert b.phi10_weight == 10
    assert b.phi10_is_delta5_squared


def test_borcherds_route_pinned_by_K3_elliptic_genus_polar_term():
    b = PhiTenBorcherdsProduct()
    # Polar term c(-1) = 1 in the Eichler-Zagier normalisation;
    # constant term c(0) = 10.
    assert b.c_minus_one == 1
    assert b.c_zero == 10


# ---------------------------------------------------------------------
# R_b: Hall--Borcherds residual via Vol III six-route convergence.
# ---------------------------------------------------------------------


def test_hall_borcherds_residual_yields_kappa_BKM_5():
    h = HallBorcherdsResidual()
    # Vol III: kappa_BKM(g_{Delta_5}) = c_1(0)/2 = 5.
    assert h.bkm_weight == 5
    assert h.doubled_trace() == 10


def test_imaginary_cartan_rank_23_real_3():
    h = HallBorcherdsResidual()
    # 24 - 1 fibre redundancy = 23 (V3
    # `gluing:prop:cartan-count`).
    assert h.cartan_rank_imaginary == 23
    # II_{3, 2} signature (3, 2): the real-root Cartan.
    assert h.cartan_rank_real == 3


def test_six_route_strata_match_V3_kappa_cat_stratification():
    h = HallBorcherdsResidual()
    # V3 `thm:cy-c-six-routes-generator-level-convergence`:
    # generator-level stratification {3, 12, 24}.
    assert set(h.six_route_strata) == {3, 12, 24}


def test_n_index_is_one_for_class_A_K3xE_canonical():
    h = HallBorcherdsResidual()
    # The N = 1 K3-fibered Class A row carries kappa_BKM = 5.
    assert h.n_index == 1


# ---------------------------------------------------------------------
# R_c: Pentagon E_3-face trace at d = 3 in SC^{ch,top}.
# ---------------------------------------------------------------------


def test_pentagon_trace_via_universal_trace_identity_doubles_to_10():
    p = PentagonE3Trace()
    # Single-colour Pentagon trace at N = 1, K3-fibered Class A.
    assert p.pentagon_trace_single_colour == 5
    # Bicoloured trace doubles: Vir_c open paired with Vir_{26-c} dual.
    assert p.trace_via_universal_trace_identity() == 10
    assert p.pentagon_trace_doubled == 10


def test_pentagon_trace_undefined_on_class_B():
    # Class B exclusion (rem:utis-class-b-exclusion) is essential:
    # quintic/C^3/conifold/local P^2 have undefined kappa_BKM.
    p = PentagonE3Trace(class_A_locus=False, locus_label="Class_B_quintic")
    try:
        p.trace_via_universal_trace_identity()
    except AssertionError:
        return
    raise AssertionError(
        "Pentagon trace should be UNDEFINED on Class B; assertion missing."
    )


# ---------------------------------------------------------------------
# Three-route concordance (the three independent verifications).
# ---------------------------------------------------------------------


def test_three_routes_agree_on_bicoloured_weight_10():
    c = ThreeRouteConcordance()
    outputs = c.route_outputs()
    assert outputs["R_a_Gritsenko_Nikulin"] == 10
    assert outputs["R_b_Hall_Borcherds_residual"] == 10
    assert outputs["R_c_Pentagon_E3_trace"] == 10
    assert c.routes_agree()


def test_route_outputs_disagree_under_route_perturbation():
    # Perturbation: take the wrong N-index; the routes desynchronise.
    bad_pentagon = PentagonE3Trace(
        pentagon_trace_single_colour=7,  # not c_1(0)/2 = 5
        pentagon_trace_doubled=14,
        class_A_locus=True,
    )
    c = ThreeRouteConcordance(pentagon=bad_pentagon)
    assert not c.routes_agree()


def test_borcherds_perturbation_breaks_concordance():
    # Borcherds with c(0) = 12 (e.g. Fake Monster Phi_12 dual-indexing
    # AP5 confusion); this is NOT the K3xE row.
    bad_borcherds = PhiTenBorcherdsProduct(
        c_zero=12,  # Fake Monster, not paramodular K3 x E
        delta5_weight=6,
        phi10_weight=12,
    )
    c = ThreeRouteConcordance(borcherds=bad_borcherds)
    # If we don't fix the cross-checks, routes disagree.
    assert not c.routes_agree()


# ---------------------------------------------------------------------
# Effectiveness: effHCSQuartic on K3 x E.
# ---------------------------------------------------------------------


def test_effHCSQuartic_vanishes_on_K3xE():
    q = HCSQuarticObstructionWitness()
    # chi_top(K3 x E) = chi_top(K3) * chi_top(E) = 24 * 0 = 0.
    assert q.chi_top_k3xe == 0
    # Quartic obstruction coefficient = chi_top(X)/24 = 0.
    assert q.quartic_obstruction_coefficient_k3xe == Fraction(0)
    assert q.k3xe_quartic_vanishes()


def test_effHCSQuartic_does_NOT_vanish_on_quintic_class_B():
    q = HCSQuarticObstructionWitness()
    # chi_top(X_5) = -200 (V3 quintic-hcs-kanom).
    assert q.chi_top_quintic == -200
    # Quartic obstruction coefficient = -200/24 != 0.
    assert q.quartic_obstruction_coefficient_quintic == Fraction(-200, 24)
    assert q.quintic_quartic_does_not_vanish()
    # This is the Beilinson falsification: P2 promotion off K3 x E fails.


def test_effHCSQuartic_chi_top_K3_is_24():
    q = HCSQuarticObstructionWitness()
    assert q.chi_top_K3 == 24
    assert q.chi_top_K3_source == "Euler class / Chern-Gauss-Bonnet"
    assert "not the source of chi_top(K3)=24" in q.hrr_role_for_K3


def test_effHCSQuartic_chi_top_E_is_0():
    # Compact CY1: chi_top(E) = 0.
    q = HCSQuarticObstructionWitness()
    assert q.chi_top_E == 0


# ---------------------------------------------------------------------
# Hall--Borcherds intertwiner data record.
# ---------------------------------------------------------------------


def test_intertwiner_has_four_cocycles():
    iv = HallBorcherdsIntertwiner()
    # V3 gluing:rem:four-cocycles-one-hopf: PL, Aut_s, M_24-inertia,
    # Delta_5-associator.
    assert len(iv.cocycles) == 4
    assert "Picard_Lefschetz" in iv.cocycles
    assert "M_24_inertia" in iv.cocycles
    assert "Delta_5_associator" in iv.cocycles
    assert iv.chain_compatible()


def test_intertwiner_source_is_Yplus_critical_CoHA():
    iv = HallBorcherdsIntertwiner()
    # Vol III `def:Y-plus-equiv-cohomology`.
    assert "Y_plus" in iv.source
    assert "critical" in iv.source.lower() or "CoHA" in iv.source


def test_intertwiner_trace_character_is_Delta5_inverse_squared():
    iv = HallBorcherdsIntertwiner()
    # Conj. gravity-line-hall-borcherds-comparison: trace character =
    # (Phi_10_un)^{-1} = Delta_5^{-2} = K3 x E BPS index.
    assert (
        "Delta5"
        in iv.trace_character.replace(" ", "").replace("-", "_")
        .replace("__", "_")
    )


def test_intertwiner_missing_cocycle_breaks_data_record():
    iv = HallBorcherdsIntertwiner(cocycles=("Picard_Lefschetz",))
    assert not iv.chain_compatible()


# ---------------------------------------------------------------------
# Full P2 witness.
# ---------------------------------------------------------------------


def test_canonical_witness_is_conditional_chain_level_candidate():
    w = canonical_p2_witness()
    assert w.status() == "conditional_chain_level_candidate"
    # The chain-level construction itself is OPEN.
    assert not w.proved_chain_level()


def test_canonical_witness_licensing_tags():
    w = canonical_p2_witness()
    # CP2 licensing: alpha + beta + delta + epsilon.
    assert set(w.licensing) == licensing_tags()
    assert "gamma" not in w.licensing  # CP2 does not need gamma


def test_stage_transport_is_S_to_A_open_lane():
    w = canonical_p2_witness()
    # CP2 stratification: S (Pentagon scalar) -> A (gravity-line operator).
    assert "S" in w.stage_source
    assert "A" in w.stage_target
    assert "open" in w.lane.lower()
    assert "Pentagon" in w.lane


def test_beilinson_falsification_record_lists_five_failure_modes():
    w = canonical_p2_witness()
    rec = w.beilinson_falsification_record()
    assert "F1_scalar_route_disagreement" in rec
    assert "F2_effHCSQuartic_failure" in rec
    assert "F3_intertwiner_incompleteness" in rec
    assert "F4_class_B_promotion" in rec
    assert "F5_dynamical_metric_overreach" in rec


def test_F4_class_B_promotion_falsified_by_quintic_quartic():
    # On Class B, the quartic obstruction is nonzero.
    q = HCSQuarticObstructionWitness()
    assert q.quintic_quartic_does_not_vanish()


def test_F5_dynamical_metric_disclaimer_text_present():
    w = canonical_p2_witness()
    rec = w.beilinson_falsification_record()
    f5 = rec["F5_dynamical_metric_overreach"]
    # Slogan #13 (Universal Holography != dynamical-metric QG).
    assert "boundary-CFT" in f5
    assert "saddle" in f5.lower()
    assert "modular" in f5.lower()
    assert "Cardy" in f5 or "BTZ" in f5


# ---------------------------------------------------------------------
# Cross-volume AP5 + AP49 checks.
# ---------------------------------------------------------------------


def test_AP5_dual_indexing_fake_monster_vs_paramodular():
    # AP5 dual-indexing: Vol I Fake Monster Phi_12 has c(0) = 24,
    # weight 12; Vol III paramodular Delta_5 has c(0) = 10, weight 5.
    # The CP2 witness uses the PARAMODULAR row, not the Fake Monster.
    b_paramod = PhiTenBorcherdsProduct(c_zero=10, delta5_weight=5)
    assert b_paramod.borcherds_route_weight() == 10  # Delta_5^2
    # If someone confuses indexing: Fake Monster Phi_12 has weight 12,
    # NOT compatible with CP2 (the paramodular K3xE row).
    b_fakemonster = PhiTenBorcherdsProduct(
        c_zero=24, delta5_weight=12, phi10_weight=24
    )
    assert b_fakemonster.borcherds_route_weight() == 24  # NOT 10
    # The two rows are distinct kappa_BKM strata; conflation is AP5.


def test_class_B_undefined_kappa_BKM():
    # Class B (quintic / C^3 / conifold / local P^2) has undefined
    # kappa_BKM (rem:utis-class-b-exclusion).
    p = PentagonE3Trace(class_A_locus=False)
    try:
        p.trace_via_universal_trace_identity()
    except AssertionError:
        return
    raise AssertionError("Class B should not admit a Pentagon trace.")
