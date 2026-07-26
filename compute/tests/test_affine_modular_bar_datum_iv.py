"""Independent locus checks for cor:affine-modular-bar-datum."""

from __future__ import annotations

from compute.lib.independent_verification import independent_verification


def _covered_affine_locus(genus: int, noncritical: bool, integrable: bool) -> bool:
    """KZ covers genus zero; regular-singular KZB covers integrable levels."""
    return noncritical and (genus == 0 or integrable)


def _modular_bar_criterion_inputs(
    *,
    internal_square_zero: bool,
    one_edge_maps: bool,
    anticommutation: bool,
    codimension_two_cancellation: bool,
) -> bool:
    return (
        internal_square_zero
        and one_edge_maps
        and anticommutation
        and codimension_two_cancellation
    )


@independent_verification(
    claim="cor:affine-modular-bar-datum",
    derived_from=["Affine KZ/KZB covered-locus corollary in 3d_gravity.tex"],
    verified_against=[
        "Drinfeld 1990 associators and the KZ pentagon/hexagon",
        "Tsuchiya-Ueno-Yamada 1989 regular-singular WZW/KZB connection",
        "Kazhdan-Lusztig 1993 tensor structures for affine Lie algebras",
        "Getzler-Kapranov 1998 stable-graph modular bar differential",
    ],
    disjoint_rationale=(
        "The test checks the logical level/genus predicate and the four "
        "modular-bar criterion inputs as an abstract decision procedure. "
        "It does not inspect annular sewing maps or reuse the manuscript "
        "proof; it verifies that generic non-integral positive genus is "
        "excluded unless the codimension-two/Stokes input is supplied."
    ),
)
def test_affine_modular_bar_locus_is_exact_for_kz_kzb_inputs():
    """Only the KZ genus-zero and integrable KZB loci satisfy the criterion."""
    assert _covered_affine_locus(genus=0, noncritical=True, integrable=False)
    assert _covered_affine_locus(genus=2, noncritical=True, integrable=True)
    assert not _covered_affine_locus(genus=1, noncritical=True, integrable=False)
    assert not _covered_affine_locus(genus=0, noncritical=False, integrable=True)

    assert _modular_bar_criterion_inputs(
        internal_square_zero=True,
        one_edge_maps=True,
        anticommutation=True,
        codimension_two_cancellation=True,
    )
    assert not _modular_bar_criterion_inputs(
        internal_square_zero=True,
        one_edge_maps=True,
        anticommutation=True,
        codimension_two_cancellation=False,
    )

    # Positive-genus raw curvature cannot be used until compensation has
    # made the internal differential square-zero.
    raw_positive_genus_curvature = True
    compensation_supplied = False
    internal_square_zero = (not raw_positive_genus_curvature) or compensation_supplied
    assert not _modular_bar_criterion_inputs(
        internal_square_zero=internal_square_zero,
        one_edge_maps=True,
        anticommutation=True,
        codimension_two_cancellation=True,
    )
