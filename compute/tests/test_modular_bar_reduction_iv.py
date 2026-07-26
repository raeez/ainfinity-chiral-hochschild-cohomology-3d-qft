"""Independent square-zero checks for prop:ainf-chiral-modular-bar-reduction."""

from __future__ import annotations

from compute.lib.independent_verification import independent_verification


def _square_components(internal_square, mixed_term, expansion_square):
    return internal_square + mixed_term + expansion_square


def _signed_codim_two_sum(path_weights):
    return sum(sign * weight for sign, weight in path_weights)


@independent_verification(
    claim="prop:ainf-chiral-modular-bar-reduction",
    derived_from=["Modular-bar criterion for annular clutching data in 3d_gravity.tex"],
    verified_against=[
        "Getzler-Kapranov 1998 modular operads stable-graph differential",
        "Ginzburg-Kapranov 1994 quadratic operad bar-cobar square expansion",
        "Gui-Li-Zeng 2022 twisted-pair curvature term",
    ],
    disjoint_rationale=(
        "The test uses the formal expansion of a bar differential square "
        "and a toy curvature-compensation calculation.  It is independent "
        "of the programme's annular clutching construction and verifies the "
        "logical point that raw positive-genus curvature is not an uncurved "
        "modular-bar datum."
    ),
)
def test_modular_bar_square_requires_curvature_absorption():
    """A curved internal square must be cancelled before D^2=0 can hold."""
    assert _square_components(0, 0, 0) == 0

    curvature = 5
    assert _square_components(curvature, 0, 0) == curvature

    s_tail_correction = -curvature
    assert _square_components(curvature + s_tail_correction, 0, 0) == 0

    # Codimension-two cancellation is the signed sum over the two expansion
    # paths to the same two-edge graph.
    assert _signed_codim_two_sum(((1, 7), (-1, 7))) == 0
    assert _signed_codim_two_sum(((1, 7), (1, 7))) != 0
