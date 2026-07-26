"""Tests for the Type-A Baxter-Rees RTT obstruction tower."""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.typeA_baxter_rees_obstruction import (
    baxter_rees_family_profile,
    beta_polynomial_terms,
    rtt_component_relation,
    rtt_obstruction_tower_profile,
    weightwise_continuation_profile,
)


def test_rtt_component_relation_records_matrix_and_components():
    relation = rtt_component_relation()

    assert relation["r_matrix"] == "R(u)=1-hbar P/u"
    assert relation["matrix_relation"] == "R(u-v)T_1(u)T_2(v)=T_2(v)T_1(u)R(u-v)"
    assert "(u-v)[T_ij(u),T_kl(v)]" in relation["component_relation"]
    assert "T_kj(u)T_il(v)-T_kj(v)T_il(u)" in relation["component_relation"]


def test_baxter_rees_family_profile_has_two_fibres():
    profile = baxter_rees_family_profile()

    assert profile["weight"] == "wt(T_ij^(r))=r"
    assert "beta^d F_d" in profile["rees_algebra"]
    assert profile["formal_family"] == "Spf completed R_beta Y_hbar(gl_N) -> Spf C[[beta]]"
    assert profile["generic_fiber"] == "RTT Yangian Y_hbar(gl_N)"
    assert profile["special_fiber"] == "associated graded RTT boundary algebra"
    assert "negative-prefundamental" in profile["boundary_packet"]


def test_weightwise_continuation_is_finite_in_each_window():
    profile = weightwise_continuation_profile(3)

    assert profile["finite_R_degree_bound"] == 3
    assert profile["finite_Theta_degree_bound"] == 3
    assert profile["r_terms"] == (
        "beta^0 R_0^(<= 3)(u)",
        "beta^1 R_1^(<= 3)(u)",
        "beta^2 R_2^(<= 3)(u)",
        "beta^3 R_3^(<= 3)(u)",
    )
    assert profile["rtt_scope"] == "satisfies RTT modulo F_{>w}"
    assert profile["reason"] == "positive RTT weights make each fixed window finite"


def test_obstruction_tower_records_tangent_and_second_class():
    tower = rtt_obstruction_tower_profile()

    assert tower["differential"] == "d_{Theta_0}=[Theta_0,-]"
    assert tower["tangent_equation"] == "d_{Theta_0} dotTheta = 0"
    assert tower["second_obstruction"] == "1/2[dotTheta,dotTheta]+d_{Theta_0}Theta_2"
    assert tower["cohomology"] == "H_RTT^2"
    assert tower["tower_meaning"] == "boundary tensor geometry is the RTT obstruction tower"


def test_beta_polynomial_terms_rejects_negative_degree():
    assert beta_polynomial_terms("Theta", 2) == (
        "beta^0 Theta_0",
        "beta^1 Theta_1",
        "beta^2 Theta_2",
    )
    with pytest.raises(ValueError):
        beta_polynomial_terms("Theta", -1)
    with pytest.raises(ValueError):
        weightwise_continuation_profile(-1)


def test_typeA_appendix_contains_explicit_rtt_rees_obstruction_tower():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    source = os.path.join(repo_root, 'chapters', 'connections', 'typeA_baxter_rees_theta.tex')
    with open(source, encoding='utf-8') as handle:
        text = handle.read()

    assert 'eq:typeA-rtt-component-relation' in text
    assert 'thm:typeA-rtt-baxter-rees-obstruction-tower' in text
    assert 'eq:typeA-rtt-rees-algebra' in text
    assert 'eq:typeA-rtt-baxter-rees-family' in text
    assert 'eq:typeA-weightwise-R-continuation' in text
    assert 'eq:typeA-weightwise-theta-continuation' in text
    assert 'eq:typeA-boundary-dot-theta' in text
    assert 'eq:typeA-boundary-tangent-cocycle' in text
    assert 'eq:typeA-boundary-second-obstruction' in text
    assert r'\frac12[\dot\Theta,\dot\Theta]+d_{\Theta_0}\Theta_2' in text
    assert 'After the polynomial $R$-matrix input' in text
    assert 'formal braided boundary germ relative to those' in text
    assert 'principles after the polynomial tensor data have been supplied' in text
    assert 'Formal braided boundary germ from polynomial continuation' in text
    assert 'Assume the hypotheses of' in text
    assert 'The full deformation complex is left to subsequent development' not in text
    assert 'yielding a canonical formal braided boundary germ' not in text
    assert 'Unconditional internal theorems whose proofs are given here' not in text
