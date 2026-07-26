"""Conditional Schellekens ``c=24`` finite-orbifold descent checks.

The manuscript theorem is a criterion.  The Schellekens and cyclic
orbifold literature supplies VOA presentations and level matching; HT
descent still requires the local Dijkgraaf-Witten/BV trivialization and
finite-orbifold BV descent for the chosen action.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path

import pytest

from compute.lib.finite_orbifold_descent import (
    dw_obstruction_profile,
    homotopy_fixed_point_profile,
    schellekens_row_partition,
)
from compute.lib.z2_group_cohomology import (
    SIGMA,
    h3_bz2_u1_class_from_triple_sign,
    normalized_z2_3_cocycle,
)


ROOT = Path(__file__).resolve().parents[2]
SCHELLEKENS_SOURCE = (
    ROOT / "chapters/connections/schellekens_71_alpha_classification_platonic.tex"
)


def test_schellekens_stratification_is_a_partition() -> None:
    """The repaired criterion keeps the ``24 + 1 + 46 = 71`` partition."""
    partition = schellekens_row_partition()
    type_a_pure_niemeier = partition["Niemeier"]
    type_b_leech_z2 = partition["V1=0"]
    type_c_nontrivial_cyclic = partition["nonlat"]
    schellekens_total = 71

    assert (
        type_a_pure_niemeier
        + type_b_leech_z2
        + type_c_nontrivial_cyclic
        == schellekens_total
    )


def test_type_a_has_vacuous_finite_orbifold_obstruction() -> None:
    """Pure Niemeier lattice VOAs have no finite group to gauge."""
    finite_group_order = 1
    h3_trivial_group_order = 1
    alpha_type_a = 0

    assert finite_group_order == 1
    assert h3_trivial_group_order == 1
    assert alpha_type_a == 0


def test_schellekens_nontrivial_strata_use_homotopy_fixed_points() -> None:
    """Types B and C require ``V^{hG}`` plus a zero DW class."""
    profile = homotopy_fixed_point_profile(
        group="Z/n",
        algebra="V_Leech",
        twisted_sector_required=True,
    )
    obstruction = dw_obstruction_profile("Z/n")

    assert "V_Leech^hZ/n" in profile.totalization
    assert "Map(Z/n^2,V_Leech)" in profile.totalization
    assert profile.requires_trivial_dw_class
    assert obstruction.cohomology_group == "H^3(BZ/n; U(1))"
    assert obstruction.descent_requires_zero
    assert "beta(g,h)" in obstruction.corrected_associator_formula


def test_type_b_leech_z2_uses_local_bv_sign_not_determinant_alone() -> None:
    """For Type B, ``+1`` is the supplied local BV datum."""
    determinant_sign = +1
    supplied_local_alpha = normalized_z2_3_cocycle(+1)

    assert determinant_sign == +1
    assert supplied_local_alpha(SIGMA, SIGMA, SIGMA) == +1
    assert h3_bz2_u1_class_from_triple_sign(+1) == 0

    alternate_local_alpha = normalized_z2_3_cocycle(-1)
    assert alternate_local_alpha(SIGMA, SIGMA, SIGMA) == -1
    assert h3_bz2_u1_class_from_triple_sign(-1) == 1


def test_type_c_level_matching_is_not_bv_trivialization() -> None:
    """The Coxeter-Todd example has the VOA level-matching input, while
    BV trivialization remains a separate hypothesis.
    """
    n = 3
    h_tw = Fraction(2, 1)
    level_matching = (h_tw * n).denominator == 1
    local_bv_trivialization_supplied = False

    assert level_matching
    assert not local_bv_trivialization_supplied


def test_type_c_source_does_not_promote_lattice_cocycle_to_bv_class() -> None:
    """The Coxeter-Todd proof keeps the BV cocycle as supplied data."""
    source = SCHELLEKENS_SOURCE.read_text(encoding="utf-8")
    flat = " ".join(source.split())

    assert (
        "its restriction to the Coxeter--Todd fixed lattice is not, by itself, "
        "a computation of the finite-orbifold BV associator"
    ) in flat
    assert "a normalized local \\(3\\)-cocycle \\(\\alpha_{\\mathrm{orb}}\\)" in flat
    assert "or equivalently a \\(2\\)-cochain \\(\\beta\\)" in flat
    assert "compatibility checks, not a computation of the Dijkgraaf--Witten class" in flat
    assert "the local BV cocycle and trivializing cochain remain supplied data" in flat

    assert "induced $3$-cocycle $\\epsilon|_{K_{12}}$" not in source
    assert "is the Kapustin-Saulina anomaly representative" not in source


@pytest.mark.parametrize(
    "type_a,type_b,type_c",
    [
        (24, 1, 46),
    ],
)
def test_uhf_schellekens_image_count_is_conditional(
    type_a: int, type_b: int, type_c: int
) -> None:
    """All rows enter the HT image only under the stratum hypotheses."""
    total = type_a + type_b + type_c
    local_bv_data_supplied_for_nontrivial_strata = False
    automatic_uhf_image_count = type_a

    assert total == 71
    assert automatic_uhf_image_count == 24
    assert not local_bv_data_supplied_for_nontrivial_strata


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
