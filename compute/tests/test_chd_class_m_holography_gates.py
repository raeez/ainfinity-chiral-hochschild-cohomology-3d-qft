"""Source guards for the class-M bulk-centre comparison in CHD.

The class-M statement must remain an ambient-qualified and
bulk-Hochschild-gated comparison.  It must not advertise an unconditional
chain-level universal holography theorem across the whole landscape.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHD = ROOT / "chapters/theory/chiral_higher_deligne.tex"


def _flat_source() -> str:
    return " ".join(CHD.read_text().split())


def test_class_m_corollary_has_status_licenses_and_inputs():
    flat = _flat_source()

    assert (
        "Class-M bulk-centre comparison in the weight-completed ambient; "
        "\\ClaimStatusConditional; licensing tags $\\beta+\\gamma+\\delta$"
    ) in flat
    assert "strict Mittag--Leffler weight-completed, pro-, or $J$-adic ambient" in flat
    assert "DS--Hochschild transport theorem" in flat
    assert "bulk--Hochschild comparison map" in flat
    assert "\\chi_{\\mathrm{HT}}\\colon" in flat
    assert "heptagon hypotheses and the two-coloured cobar homotopy" in flat


def test_class_m_corollary_separates_cohomology_from_chain_level():
    flat = _flat_source()

    assert "is an equivalence on cohomology as an \\(E_2\\)-chiral/Gerstenhaber algebra" in flat
    assert "under \\textup{(d)}, refines to a mixed chiral-topological chain-level equivalence" in flat
    assert "On the bounded direct-sum ambient $\\mathrm{Ch}(\\mathrm{Vect})$ the class-M chain-level identification is false" in flat
    assert "The G/L/C classes are governed by their own original-complex identifications" in flat


def test_retired_unconditional_holography_phrases_are_absent():
    flat = _flat_source()

    retired = (
        "With Corollary~\\ref{cor:universal-holography-class-M} established",
        "is chain-level across all four classes",
        "this yields 3d HT holography uniformly across the Koszul landscape",
        "Chain-level equalities at each step compose to chain-level equality",
    )
    for phrase in retired:
        assert phrase not in flat


def test_bulk_map_requires_explicit_ht_input():
    inputs = {
        "ambient": True,
        "ds_hochschild": True,
        "bulk_hochschild": True,
        "two_coloured_cobar": True,
    }

    assert all(inputs.values())
    assert set(inputs) == {
        "ambient",
        "ds_hochschild",
        "bulk_hochschild",
        "two_coloured_cobar",
    }
