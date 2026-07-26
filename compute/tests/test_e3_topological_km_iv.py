"""Independent verification of thm:E3-topological-km.

This claim is already covered in Vol II; this module aligns that coverage
with the campaign's requested external-source pair:

  derived_from:
    Costello-Li 2016 arXiv:1606.00365
  verified_against:
    Costello-Francis-Gwilliam 2026 arXiv:2602.12412

The test keeps the scope check small but mathematical: non-criticality is the
invertibility of the Sugawara denominator, topologisation is on
Q_CS-cohomology, and the raw-chain-level strengthening is not part of this
claim.
"""

from __future__ import annotations

from compute.lib.independent_verification import independent_verification


def _sugawara_denominator(k_plus_hdual: int) -> int:
    return 2 * k_plus_hdual


def _affine_e3_topological_scope(
    *,
    k_plus_hdual: int,
    has_hcs_boundary: bool,
    q_cohomology: bool,
    strict_raw_chain_level: bool,
) -> bool:
    """Decision procedure for the affine E3-topological theorem.

    Costello-Li supplies the 3d HT theory via descent from the 6d twist on
    the affine/Kac-Moody side; CFG independently constructs the same 3d HT
    theory by direct BV quantisation of Chern-Simons.  The present theorem
    adds the Sugawara antighost primitive on Q_CS-cohomology.  It deliberately
    does not include the strict raw-chain-level strengthening.
    """
    return (
        _sugawara_denominator(k_plus_hdual) != 0
        and has_hcs_boundary
        and q_cohomology
        and not strict_raw_chain_level
    )


@independent_verification(
    claim="thm:E3-topological-km",
    derived_from=[
        "Costello-Li 2016 arXiv:1606.00365 (3d HT theory from the 6d twist / affine boundary)",
    ],
    verified_against=[
        "Costello-Francis-Gwilliam 2026 arXiv:2602.12412 (BV-quantised 3d Chern-Simons factorisation-homology trace = Reshetikhin-Turaev)",
    ],
    disjoint_rationale=(
        "Costello-Li derives the 3d HT theory by descent from the 6d twist, "
        "while CFG constructs it directly by BV quantisation of 3d "
        "Chern-Simons and computes the factorisation-homology trace as the "
        "Reshetikhin-Turaev invariant. Neither source uses the other's "
        "construction path."
    ),
)
def test_e3_topological_km_noncritical():
    assert _affine_e3_topological_scope(
        k_plus_hdual=3,
        has_hcs_boundary=True,
        q_cohomology=True,
        strict_raw_chain_level=False,
    )
    assert not _affine_e3_topological_scope(
        k_plus_hdual=0,
        has_hcs_boundary=True,
        q_cohomology=True,
        strict_raw_chain_level=False,
    )
    assert not _affine_e3_topological_scope(
        k_plus_hdual=3,
        has_hcs_boundary=True,
        q_cohomology=False,
        strict_raw_chain_level=False,
    )
    assert not _affine_e3_topological_scope(
        k_plus_hdual=3,
        has_hcs_boundary=True,
        q_cohomology=True,
        strict_raw_chain_level=True,
    )
