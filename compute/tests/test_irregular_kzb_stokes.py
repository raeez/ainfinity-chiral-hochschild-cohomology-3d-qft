"""Tests for irregular KZB formal type and Stokes cocycle profiles."""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.irregular_kzb_stokes import (
    boundary_connection_profile,
    clutching_profile,
    covered_level_locus,
    formal_gauge_profile,
    formal_type_profile,
    kzb_vs_curved_dunn_profile,
    stokes_sector_profile,
    uncovered_level_locus,
)


def test_boundary_connection_profile_records_q_normal_form():
    profile = boundary_connection_profile(3)

    assert profile["chart"] == "alpha=(q,t_alpha)"
    assert profile["polar_terms"] == ("A_3/q^3", "A_2/q^2", "A_1/q")
    assert profile["regular_q_term"] == "A_0 dq"
    assert profile["transverse_terms"] == "sum_alpha B_alpha(q) dt_alpha"
    assert profile["irregular_requires_higher_pole"] is True
    assert boundary_connection_profile(1)["irregular_requires_higher_pole"] is False


def test_formal_type_and_gauge_profile_record_jmu_shape():
    formal_type = formal_type_profile(4)
    gauge = formal_gauge_profile(5)

    assert formal_type["terms"] == (
        "A_2/(1 q^1)",
        "A_3/(2 q^2)",
        "A_4/(3 q^3)",
        "A_1 log q",
    )
    assert "Theta_partial(q)=sum_{r=1}^{m-1}" in formal_type["theta"]
    assert gauge["transformation"] == "G(q) in Aut(V)[[q^(1/5)]]"
    assert gauge["normal_form"] == "nabla_q ~= d_q-d_q Theta_partial-Lambda dq/q"


def test_stokes_sector_and_clutching_profiles_record_cocycle():
    stokes = stokes_sector_profile()
    clutching = clutching_profile()

    assert stokes["sectorial_solution"] == (
        "Y_ell(q)=H_ell(q) exp(Theta_partial(q)) q^Lambda"
    )
    assert stokes["stokes_matrix"] == "S_{ell,ell'}=Y_ell^{-1}Y_ell'"
    assert stokes["cocycle"] == "S_{ell1,ell2} S_{ell2,ell3} S_{ell3,ell1}=1"
    assert "Theta_{partial,Gamma1#Gamma2}" in clutching["formal_type"]
    assert "S^{Gamma1#Gamma2}_{ell,ell'}" in clutching["stokes_product"]
    assert clutching["order"] == "tangential basepoint convention"


def test_level_loci_are_exactly_the_covered_rows():
    assert covered_level_locus("integrable") == {
        "level": "k in Z_{>=0}",
        "mechanism": "KZ pentagon + Kazhdan-Lusztig regularity",
    }
    assert covered_level_locus("generic_nonrational") == {
        "level": "k in C\\Q",
        "mechanism": "irregular KZB formal type + wild groupoid",
    }
    assert covered_level_locus("critical") == {
        "level": "k=-h^vee",
        "mechanism": "Feigin-Frenkel centre",
    }
    assert uncovered_level_locus()["required_input"] == (
        "admissible-level comparison theorem"
    )


def test_kzb_stokes_is_not_curved_dunn_h2_vanishing():
    profile = kzb_vs_curved_dunn_profile()

    assert profile["kzb_stokes_supplies"] == (
        "boundary formal type",
        "mixed Stokes matrices",
        "wild-groupoid sector-crossing cocycle",
    )
    assert profile["curved_dunn_h2_requires"] == (
        "modular-bootstrap H^2 acyclicity",
        "degree-2 comparison map",
        "transport of zero cohomology class",
    )
    assert not set(profile["kzb_stokes_supplies"]) & set(
        profile["curved_dunn_h2_requires"]
    )


def test_profiles_reject_impossible_orders():
    with pytest.raises(ValueError):
        boundary_connection_profile(0)
    with pytest.raises(ValueError):
        formal_type_profile(0)
    with pytest.raises(ValueError):
        formal_gauge_profile(0)


def test_manuscript_contains_explicit_kzb_formal_type_labels():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    source = os.path.join(
        repo_root, 'chapters', 'theory', 'modular_swiss_cheese_operad.tex'
    )
    curved_source = os.path.join(
        repo_root, 'chapters', 'theory', 'curved_dunn_higher_genus.tex'
    )

    with open(source, encoding='utf-8') as handle:
        text = handle.read()
    with open(curved_source, encoding='utf-8') as handle:
        curved_text = handle.read()

    assert 'prop:irregular-kzb-q-formal-type-stokes-cocycle' in text
    assert 'eq:irregular-kzb-boundary-connection-form' in text
    assert 'eq:irregular-kzb-boundary-formal-type' in text
    assert 'eq:irregular-kzb-formal-gauge-form' in text
    assert 'eq:irregular-kzb-sectorial-solution' in text
    assert 'eq:irregular-kzb-stokes-matrix' in text
    assert 'eq:irregular-kzb-clutching-formal-type' in text
    assert 'eq:irregular-kzb-clutching-stokes-product' in text
    assert 'eq:irregular-kzb-stokes-cocycle' in text
    assert 'rational nonintegral' in text
    assert 'admissible-level comparison' in text
    assert 'it is not the curved-Dunn \\(H^{2}\\)-acyclicity mechanism' in curved_text
