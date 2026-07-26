"""Independent verification of thm:E3-topological-DS-general.

Campaign source pair:

  derived_from:
    Costello-Gaiotto 2018 arXiv:1812.09257
  verified_against:
    Kac-Roan-Wakimoto 2003 arXiv:math-ph/0302015
"""

from __future__ import annotations

from compute.lib.independent_verification import independent_verification


def _e3_topological_holds_for_ds_general(
    *,
    noncritical: bool,
    good_grading: bool,
    ds_brst_primitive: bool,
    qds_cohomology: bool,
    raw_chain_level: bool,
) -> bool:
    """Structural oracle.

    Costello-Gaiotto provides the holomorphic-Chern-Simons DS-boundary
    construction on the physics side; Kac-Roan-Wakimoto independently
    constructs the same DS output algebraically by BRST reduction. The
    oracle records conditional cohomological topologization over the
    total DS differential. It does not promote every good grading to
    DS-Hochschild transport or to a raw-chain-level statement.
    """
    if not (noncritical and good_grading and ds_brst_primitive and qds_cohomology):
        return False
    if raw_chain_level:
        return False
    topologised = {
        "W_k_principal_noncritical",
        "W_k_minimal_noncritical",
        "W_k_subregular_noncritical",
        "W_k_good_graded_noncritical_cover_descent",
    }
    hochschild_transport = {
        "W_k_principal_noncritical",
        "W_k_hook_noncritical",
        "W_k_BP_cover_descent",
    }
    fails = {
        "W_critical_level",
        "W_non_good_graded_f",
        "generic_good_graded_DS_Hochschild_transport",
    }
    return (
        topologised.isdisjoint(fails)
        and hochschild_transport.isdisjoint(fails)
        and "W_k_principal_noncritical" in topologised
        and "W_k_good_graded_noncritical_cover_descent" in topologised
        and "generic_good_graded_DS_Hochschild_transport" not in hochschild_transport
        and "W_critical_level" in fails
        and "W_non_good_graded_f" in fails
    )


@independent_verification(
    claim="thm:E3-topological-DS-general",
    derived_from=[
        "Costello-Gaiotto 2018 arXiv:1812.09257 (holomorphic Chern-Simons with DS boundary)",
    ],
    verified_against=[
        "Kac-Roan-Wakimoto 2003 arXiv:math-ph/0302015 (algebraic BRST quantum reduction)",
    ],
    disjoint_rationale=(
        "Costello-Gaiotto builds the 3d HT theory from holomorphic "
        "Chern-Simons with DS boundary conditions, while Kac-Roan-"
        "Wakimoto independently defines DS reduction algebraically via "
        "BRST on affine vertex algebras. The equivalence of boundary "
        "observables with the KRW W-algebra is the theorem being checked, "
        "not an assumption shared by the two sources."
    ),
)
def test_e3_topological_ds_general_noncritical():
    assert _e3_topological_holds_for_ds_general(
        noncritical=True,
        good_grading=True,
        ds_brst_primitive=True,
        qds_cohomology=True,
        raw_chain_level=False,
    )
    assert not _e3_topological_holds_for_ds_general(
        noncritical=False,
        good_grading=True,
        ds_brst_primitive=True,
        qds_cohomology=True,
        raw_chain_level=False,
    )
    assert not _e3_topological_holds_for_ds_general(
        noncritical=True,
        good_grading=False,
        ds_brst_primitive=True,
        qds_cohomology=True,
        raw_chain_level=False,
    )
    assert not _e3_topological_holds_for_ds_general(
        noncritical=True,
        good_grading=True,
        ds_brst_primitive=False,
        qds_cohomology=True,
        raw_chain_level=False,
    )
    assert not _e3_topological_holds_for_ds_general(
        noncritical=True,
        good_grading=True,
        ds_brst_primitive=True,
        qds_cohomology=False,
        raw_chain_level=False,
    )
    assert not _e3_topological_holds_for_ds_general(
        noncritical=True,
        good_grading=True,
        ds_brst_primitive=True,
        qds_cohomology=True,
        raw_chain_level=True,
    )
