"""Source guards for W_N self-dual halving scope.

The W_N self-dual point halves the kappa curvature.  That immediately
halves the scalar lane.  It halves the full multi-weight genus tower only
after a separate cross-channel equality, which the manuscript must not
silently assert.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FRONTIER = ROOT / "chapters/examples/w-algebras-frontier.tex"
STABLE = ROOT / "chapters/examples/w-algebras-stable.tex"


def _flat(path: Path) -> str:
    return "".join(path.read_text(encoding="utf-8").split())


def test_live_frontier_halving_is_scalar_lane_not_full_multi_weight():
    source = _flat(FRONTIER)

    assert "inheritascalar-lanepropertyfromthetable" in source
    assert r"F_1^{\mathrm{sc}}(\cW_{N,c})=\kappaChHodge(\cW_{N,c})\lambda_1^{\mathrm{FP}}" in source
    assert r"F_1^{\mathrm{sc}}(c^*)=\tfrac{1}{2}F_1^{\mathrm{sc}}(\alpha_N)" in source
    assert r"ForVirasoro($N=2$,uniform-weight),thescalarlaneisthefullgenustower" in source
    assert r"F_g(\cW_{N,c})=F_g^{\mathrm{sc}}(\cW_{N,c})+\deltaF_g^{\mathrm{cross}}(\cW_{N,c})" in source
    assert r"\deltaF_g^{\mathrm{cross}}(\cW_{N,c^*})=\tfrac12\,\deltaF_g^{\mathrm{cross}}(\cW_{N,\alpha_N})" in source
    assert "whichisnotassertedhere" in source
    assert r"F_1(c^*)=\tfrac{1}{2}\,F_1(\alpha_N)" not in source


def test_stable_copy_quarantines_same_halving_scope():
    source = _flat(STABLE)

    assert r"\begin{remark}[Scalar-laneself-dualhalving]" in source
    assert r"F_1^{\mathrm{sc}}(c^*)=\tfrac{1}{2}\,F_1^{\mathrm{sc}}(\alpha_N)" in source
    assert r"Fullmulti-weighthalvingwouldrequiretheseparateequality" in source
    assert r"F_1(c^*)=\tfrac{1}{2}\,F_1(\alpha_N)" not in source
    assert r"F_1(c^*)=(\alpha_N/4)\lambda_1^{\mathrm{FP}}" not in source
    assert r"F_1^{\mathrm{sc}}(c^*)=(\alpha_N/4)\lambda_1^{\mathrm{FP}}" in source


def test_cross_channel_equality_is_independent_of_scalar_halving():
    def scalar_halves(kappa_at_critical: Fraction, kappa_at_ghost: Fraction) -> bool:
        return kappa_at_critical == kappa_at_ghost / 2

    def full_halves(
        kappa_at_critical: Fraction,
        kappa_at_ghost: Fraction,
        cross_at_critical: Fraction,
        cross_at_ghost: Fraction,
    ) -> bool:
        return scalar_halves(kappa_at_critical, kappa_at_ghost) and (
            cross_at_critical == cross_at_ghost / 2
        )

    assert scalar_halves(Fraction(50), Fraction(100))
    assert full_halves(Fraction(50), Fraction(100), Fraction(7), Fraction(14))
    assert not full_halves(Fraction(50), Fraction(100), Fraction(7), Fraction(13))
