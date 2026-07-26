"""Tests for the shifted RTT rank-one Kleinian test."""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.shifted_rtt_kleinian_test import (
    kleinian_boundary_relations,
    pairing_annihilator_profile,
    quantum_determinant_profile,
    shifted_generator_profile,
    shifted_rtt_candidate_passes_rank_one_test,
)


def test_shifted_generator_uses_chamber_exponent():
    profile = shifted_generator_profile(mu_i=4, mu_j=1)

    assert profile["u_exponent"] == -3
    assert "u^{-mu_i+mu_j}" in profile["formula"]
    assert "delta_ij+sum_{r>=1}T_ij^(r)u^{-r}" in profile["formula"]


def test_pairing_annihilator_profile_records_coideal_condition():
    profile = pairing_annihilator_profile()

    assert "I_mu={x in Y_hbar^RTT" in profile["ideal"]
    assert "Y_hbar^{RTT,vee}_{>=mu}" in profile["ideal"]
    assert profile["coideal_condition"] == "Delta(I_mu) subset I_mu tensor Y_hat + Y_hat tensor I_mu"
    assert profile["interpretation"] == "without the coideal condition the quotient is not a line algebra"


def test_quantum_determinant_has_factorial_terms_and_boundary_polynomial():
    profile = quantum_determinant_profile(4)

    assert profile["term_count"] == 24
    assert profile["shifts"] == ("u-0hbar", "u-1hbar", "u-2hbar", "u-3hbar")
    assert "prod_{j=1}^N" in profile["formula"]
    assert profile["boundary_value"] == "qdet T(u)=P_mu(u)"


def test_rank_one_kleinian_relations_and_associated_graded():
    profile = kleinian_boundary_relations(5)

    assert profile["commutators"] == ("[z,x]-hbar x", "[z,y]+hbar y")
    assert profile["yx_factors"] == (
        "z-0hbar",
        "z-1hbar",
        "z-2hbar",
        "z-3hbar",
        "z-4hbar",
    )
    assert profile["xy_factors"] == (
        "z+1hbar",
        "z+2hbar",
        "z+3hbar",
        "z+4hbar",
        "z+5hbar",
    )
    assert profile["associated_graded"] == "C[x,y,z]/(xy-z^5)"
    assert profile["kleinian_type"] == "A_4"


def test_candidate_requires_coideal_and_kleinian_deformation():
    assert shifted_rtt_candidate_passes_rank_one_test(
        has_coideal=True,
        gr_relation="C[x,y,z]/(xy-z^3)",
        m=3,
    )
    assert not shifted_rtt_candidate_passes_rank_one_test(
        has_coideal=False,
        gr_relation="C[x,y,z]/(xy-z^3)",
        m=3,
    )
    assert not shifted_rtt_candidate_passes_rank_one_test(
        has_coideal=True,
        gr_relation="C[x,y,z]/(xy-z^2)",
        m=3,
    )


def test_invalid_rank_and_kleinian_m_rejected():
    with pytest.raises(ValueError):
        quantum_determinant_profile(0)
    with pytest.raises(ValueError):
        kleinian_boundary_relations(1)


def test_shifted_rtt_source_contains_pairing_annihilator_kleinian_test():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    source = os.path.join(
        repo_root,
        'chapters',
        'connections',
        'shifted_rtt_duality_orthogonal_coideals.tex',
    )
    with open(source, encoding='utf-8') as handle:
        text = handle.read()

    assert 'thm:shifted-rtt-pairing-annihilator-kleinian-test' in text
    assert 'eq:shifted-rtt-pairing-annihilator' in text
    assert 'eq:shifted-rtt-coideal-condition' in text
    assert 'eq:shifted-rtt-local-generators' in text
    assert 'eq:shifted-rtt-qdet-polynomial' in text
    assert 'eq:rank-one-quantized-kleinian-boundary' in text
    assert 'eq:rank-one-kleinian-associated-graded-general' in text
    assert r'I_\mu' in text
    assert r'\operatorname{qdet}T(u)' in text
    assert r'\mathbb C[x,y,z]/(xy-z^m)' in text
