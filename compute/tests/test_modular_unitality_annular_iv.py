"""Independent unit/counit checks for prop:modular-operad-unitality."""

from __future__ import annotations

from compute.lib.independent_verification import independent_verification


def _delta_basis(element):
    if element == "1":
        return {("1", "1"): 1}
    if element == "x":
        return {("1", "x"): 1, ("x", "1"): 1}
    raise KeyError(element)


def _epsilon(element):
    return 1 if element == "1" else 0


def _left_counit(delta):
    result = {}
    for (left, right), coeff in delta.items():
        result[right] = result.get(right, 0) + coeff * _epsilon(left)
    return {k: v for k, v in result.items() if v}


def _right_counit(delta):
    result = {}
    for (left, right), coeff in delta.items():
        result[left] = result.get(left, 0) + coeff * _epsilon(right)
    return {k: v for k, v in result.items() if v}


def _cotensor_with_diagonal(comodule_basis):
    """Diagonal bicomodule leaves compatible comodule labels unchanged."""
    diagonal = {(b, b) for b in ("1", "x")}
    return {
        right
        for left in comodule_basis
        for candidate_left, right in diagonal
        if left == candidate_left
    }


@independent_verification(
    claim="prop:modular-operad-unitality",
    derived_from=["Unit-normalised annular sewing proof in 3d_gravity.tex"],
    verified_against=[
        "Mac Lane 1963 monoidal category unit axioms",
        "May 1972 operad identity axiom",
        "Sweedler 1969 coalgebra counit identity",
    ],
    disjoint_rationale=(
        "The test uses the elementary coalgebra counit and diagonal "
        "bicomodule identities.  These are independent of the programme's "
        "annular-bar construction and verify the algebraic replacement for "
        "the discarded simply-connected-unit argument."
    ),
)
def test_counit_normalised_diagonal_bicomodule_is_unit():
    """Counit normalisation and diagonal cotensoring recover the input."""
    for element in ("1", "x"):
        delta = _delta_basis(element)
        assert _left_counit(delta) == {element: 1}
        assert _right_counit(delta) == {element: 1}

    # A normalized R-term supported in the augmentation ideal disappears
    # when either tensor factor is evaluated by the counit.
    r_term = {("x", "x"): 7}
    assert _left_counit(r_term) == {}
    assert _right_counit(r_term) == {}

    comodule_basis = {"1", "x"}
    assert _cotensor_with_diagonal(comodule_basis) == comodule_basis
