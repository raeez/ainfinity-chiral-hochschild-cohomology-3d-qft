"""Tests for the punctured-disk versus circle E_2/E_1 gap."""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.punctured_disk_circle_gap import (
    angular_projection_effect,
    punctured_disk_circle_gap_profile,
    rotation_bv_operator,
)


class TestPuncturedDiskCircleGap:
    def test_rotation_bv_operator_profile(self):
        profile = rotation_bv_operator(3)

        assert profile["operator"] == "Delta_rot = sum_i iota_{partial_arg_z_i}"
        assert profile["contractions"] == (
            "iota_partial_arg_z_1",
            "iota_partial_arg_z_2",
            "iota_partial_arg_z_3",
        )
        assert profile["survives_on"] == "D^* punctured-disk chain model"
        assert profile["circle_image"] == "ordinary cyclic rotation"

    def test_gap_profile_separates_qiso_scope_from_full_operadic_module(self):
        profile = punctured_disk_circle_gap_profile(2)

        assert profile["quasi_isomorphism_scope"] == "angular E_1 quotient only"
        assert profile["full_operadic_module_equivalence"] is False
        assert "theta_i = d arg z_i" in profile["gap_generators"]
        assert "omega_ij = d log(z_i - z_j)" in profile["gap_generators"]
        assert "Delta_rot = sum_i iota_{partial_arg_z_i}" in profile["gap_generators"]

    def test_circle_and_punctured_disk_compute_different_chain_models(self):
        profile = punctured_disk_circle_gap_profile(4)

        assert profile["circle_computes"] == "annular trace"
        assert profile["punctured_disk_computes"] == "chiral centre with rotation"
        assert "Conf_n(D^*)" in profile["symmetric_model"]
        assert "Conf_n(S^1)" in profile["circle_model"]

    def test_angular_projection_effect(self):
        effect = angular_projection_effect()

        assert "cyclic order" in effect["preserves"]
        assert "wraparound monodromy" in effect["preserves"]
        assert "radial collision data" in effect["forgets"]
        assert "rotation BV operator as D^* chain operation" in effect["forgets"]

    def test_nonpositive_arity_rejected(self):
        with pytest.raises(ValueError):
            rotation_bv_operator(0)
        with pytest.raises(ValueError):
            punctured_disk_circle_gap_profile(-1)

    def test_hochschild_source_contains_rotation_gap(self):
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        source = os.path.join(repo_root, 'chapters', 'connections', 'hochschild.tex')
        with open(source, encoding='utf-8') as handle:
            text = handle.read()

        assert 'prop:punctured-disk-circle-rotation-gap' in text
        assert 'eq:punctured-disk-symmetric-bar' in text
        assert 'eq:punctured-disk-rotation-bv' in text
        assert r'\Delta_{\mathrm{rot}}' in text
        assert r'\theta_i=d\arg z_i' in text
        assert r'\omega_{ij}=d\log(z_i-z_j)' in text
        assert 'angular' in text
        assert '\\(E_1\\)-quotient' in text

    def test_raviolo_source_uses_angular_quotient(self):
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        source = os.path.join(repo_root, 'chapters', 'theory', 'raviolo-restriction.tex')
        with open(source, encoding='utf-8') as handle:
            text = handle.read()

        assert r'C^{\mathrm{ch,ang}}(D^\times,\cA)' in text
        assert r'\Delta_{\mathrm{rot}}=\sum_i\iota_{\partial_{\arg z_i}}' in text
        assert 'It does not declare the full \\(E_2/E_1\\) relative' in text
