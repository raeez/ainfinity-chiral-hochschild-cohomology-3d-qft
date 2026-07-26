"""Source guards for the semisimple Pixton/MC tower theorem.

The theorem uses the FSZ/PPZ all-r-spin route.  A fixed scalar value of
kappa cannot generate the whole Pixton ideal by that route; the manuscript
must keep the statement on a scalar family containing every r-spin point.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HT_BULK_BOUNDARY_LINE = ROOT / "chapters/connections/ht_bulk_boundary_line.tex"


def _flat_source() -> str:
    return " ".join(HT_BULK_BOUNDARY_LINE.read_text().split())


def test_pixton_theorem_has_scalar_family_status_and_licenses():
    flat = _flat_source()

    assert (
        "Pixton ideal from the scalar MC tower on the semisimple "
        "$r$-spin locus; \\ClaimStatusProvedHere{} under the displayed "
        "scalar-family hypotheses; "
        "licensing tags $\\alpha+\\beta+\\gamma+\\delta$"
    ) in flat
    assert "ClaimStatusProvedHereConditional" not in flat
    assert "\\kappa_r=\\frac{r-2}{2r}" in flat
    assert "U\\supset \\left\\{\\kappa_r" in flat
    assert "r\\ge 3" in flat


def test_pixton_theorem_keeps_fsz_ppz_hypotheses_explicit():
    flat = _flat_source()

    assert "semisimple rank-one CohFT with flat unit" in flat
    assert "Givental--Teleman classification" in flat
    assert "Faber--Shadrin--Zvonkine $r$-spin $R$-matrix" in flat
    assert "Pandharipande--Pixton--Zvonkine" in flat
    assert "$3,4,5,\\ldots$" in flat
    assert "the Pixton ideal" in flat


def test_pixton_theorem_does_not_advertise_all_standard_families():
    flat = _flat_source()

    assert "This covers all standard families at generic level" not in flat
    assert "For every chirally Koszul algebra with semisimple genus-$0$" not in flat
    assert "This is a scalar-lane statement" in flat
    assert "does not assert a full tensor-channel theorem" in flat
    assert (
        "does not cover admissible-level, resonant, or logarithmic "
        "non-semisimple specialisations"
    ) in flat


def test_fixed_kappa_cannot_supply_the_all_r_ppz_route():
    def kappa_r(r: int) -> Fraction:
        return Fraction(r - 2, 2 * r)

    fixed_family_values = {kappa_r(3)}
    ppz_required_values = {kappa_r(r) for r in range(3, 10)}

    assert len(ppz_required_values) == 7
    assert fixed_family_values < ppz_required_values
    assert kappa_r(3) != kappa_r(4)
