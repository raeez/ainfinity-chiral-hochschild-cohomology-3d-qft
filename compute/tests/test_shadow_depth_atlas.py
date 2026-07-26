from __future__ import annotations

import os
import sys
from math import inf
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.shadow_depth_atlas import (
    WILD_BOUNDARY,
    class_profile,
    koszul_shadow_depth_values,
    shadow_depth_defined,
    wild_boundary_statement,
)


def test_koszul_atlas_has_exact_shadow_depth_values():
    assert koszul_shadow_depth_values() == (2, 3, 4, inf)
    assert class_profile("G").chirally_koszul is True
    assert class_profile("L").shadow_depth == 3
    assert class_profile("C").algebraic_depth == 2
    assert class_profile("M").shadow_depth == inf


def test_wild_boundary_is_not_class_m_refinement():
    wild = class_profile("W")
    mixed = class_profile("M")

    assert mixed.chirally_koszul is True
    assert mixed.shadow_depth == inf
    assert wild == WILD_BOUNDARY
    assert wild.chirally_koszul is False
    assert wild.shadow_depth is None
    assert wild.algebraic_depth is None
    assert shadow_depth_defined("W") is False
    assert "Kronecker K_m" in wild.example


def test_wild_boundary_statement_names_the_real_open_problem():
    statement = wild_boundary_statement()

    assert statement["koszul_sequent"] == (
        "A in Kosz_ch => r_sh(A) in {2,3,4,infty}"
    )
    assert statement["wild_depth"] is None
    assert statement["wild_depth_text"] == "r_sh(A) undefined"
    assert statement["not_open_problem"] == "whether infinite depth exists"
    assert statement["class_m_examples"] == ("Virasoro", "W_N")
    assert "refuses the Koszul shadow" in statement["open_problem"]


def test_shadow_depth_atlas_source_guard():
    root = Path(__file__).resolve().parents[2]
    intro = (root / "chapters/theory/introduction.tex").read_text()
    rft = (root / "chapters/connections/relative_feynman_transform.tex").read_text()

    combined = " ".join((intro + "\n" + rft).split())
    assert "A\\in\\operatorname{Kosz}_{\\mathrm{ch}}" in combined
    assert "r_{\\mathrm{sh}}(A)\\in\\{2,3,4,\\infty\\}" in combined
    assert "\\mathbf{W}\\ \\text{(Wild)}" in combined
    assert "outside the chirally Koszul locus" in combined
    assert "shadow depth is undefined" in combined
    assert "The open problem is not the existence of infinite depth" in combined
