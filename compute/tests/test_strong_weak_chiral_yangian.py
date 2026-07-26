"""Tests for the strong/weak chiral Yangian comparison obstruction."""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.strong_weak_chiral_yangian import (
    comparison_map_profile,
    monodromy_rmatrix_statement_status,
    strong_chiral_yangian_profile,
    weak_chiral_yangian_profile,
)


def test_weak_chiral_yangian_profile_is_rtt_factorization():
    profile = weak_chiral_yangian_profile()

    assert profile["datum"] == "(A,R(z),tensor_z)"
    assert "R12(z-w)R13(z)R23(w)" in profile["relation"]
    assert profile["deformation_problem"] == "RTTDef(A,R)"


def test_strong_chiral_yangian_profile_is_modular_mc():
    profile = strong_chiral_yangian_profile()

    assert "{m_k}_{k>=2}" in profile["datum"]
    assert "Theta_mod" in profile["datum"]
    assert profile["mc_equation"] == "d Theta_mod + 1/2[Theta_mod,Theta_mod]=0"
    assert profile["deformation_problem"] == "MC_mod(A)"


def test_comparison_map_tangent_and_obstruction():
    profile = comparison_map_profile()

    assert profile["map"] == "Psi_RTT_to_mod: RTTDef(A,R) -> MC_mod(A)"
    assert profile["tangent"] == "dPsi_RTT_to_mod(dotR(z))=Theta_1"
    assert "d_{Theta_0}Theta_2" in profile["second_obstruction"]
    assert "1/2[Theta_1,Theta_1]" in profile["second_obstruction"]
    assert "H_mod^2(A)" in profile["second_obstruction"]


def test_monodromy_equals_rmatrix_needs_all_strong_inputs():
    assert monodromy_rmatrix_statement_status(True, True, True) == "strong modular-MC theorem"
    assert (
        monodromy_rmatrix_statement_status(False, True, True)
        == "weak RTT/factorization statement only"
    )
    assert (
        monodromy_rmatrix_statement_status(True, False, True)
        == "weak RTT/factorization statement only"
    )
    assert (
        monodromy_rmatrix_statement_status(True, True, False)
        == "weak RTT/factorization statement only"
    )


def test_line_operator_source_contains_strong_weak_obstruction_theorem():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    source = os.path.join(repo_root, 'chapters', 'connections', 'line-operators.tex')
    with open(source, encoding='utf-8') as handle:
        text = handle.read()

    assert 'thm:strong-weak-chiral-yangian-comparison-obstruction' in text
    assert 'eq:strong-chiral-yangian-modular-mc' in text
    assert 'eq:rttdef-to-modmc-comparison-map' in text
    assert 'eq:rttdef-tangent-to-theta-one' in text
    assert 'eq:strong-weak-chiral-yangian-ob2' in text
    assert r'\operatorname{RTTDef}(A,R)' in text
    assert r'\operatorname{MC}_{\mathrm{mod}}(A)' in text
    assert r'd\Psi_{\mathrm{RTT}\to\mathrm{mod}}(\dot R(z))=\Theta_1' in text


def test_cross_reference_propagation_to_existing_two_definition_remarks():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    paths = [
        os.path.join(repo_root, 'chapters', 'connections', 'ordered_associative_chiral_kd_core.tex'),
        os.path.join(repo_root, 'chapters', 'connections', 'dg_shifted_factorization_bridge.tex'),
        os.path.join(repo_root, 'chapters', 'connections', 'ht_physical_origins.tex'),
    ]

    for path in paths:
        with open(path, encoding='utf-8') as handle:
            text = handle.read()
        assert 'thm:strong-weak-chiral-yangian-comparison-obstruction' in text
