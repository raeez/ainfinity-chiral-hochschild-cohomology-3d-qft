"""Logarithmic boundary-amplitude obligations for triplet W(p).

The W(p) tempering problem has two independent pieces:

* the regular semisimple TW/WW shadow bound, recorded by
  ``conj:wp-regular-sector-amplitude-bound``;
* the logarithmic boundary-changing estimate for the phi_{0,1}
  sector, recorded by
  ``conj:wp-logarithmic-boundary-amplitude-bound``.

Finite-dimensional Zhu data and C2-cofiniteness do not imply the
second bound.  This module only encodes the exact proof obligation and
the elementary Stirling consequence of a subfactorial logarithmic
bound.  It does not compute the actual triplet shadow tower.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import exp, lgamma, log


@dataclass(frozen=True)
class LogarithmicShadowComponent:
    """One summand in S_r(W(p)) = S_reg + S_nil + S_log."""

    name: str
    source: str
    status: str
    description: str


def reg_nil_log_decomposition_profile() -> dict[str, object]:
    """Return the semisimple/nilpotent/logarithmic shadow split."""

    return {
        "formula": "S_r(W(p)) = S_r^reg + S_r^nil + S_r^log",
        "regular_refinement": (
            "S_r^reg = S_r^TT + S_r^TW + S_r^WW on the semisimple OPE sector"
        ),
        "components": (
            LogarithmicShadowComponent(
                name="S_r^reg",
                source="conj:wp-regular-sector-amplitude-bound",
                status="conjectural outside TT",
                description="ordinary semisimple propagation",
            ),
            LogarithmicShadowComponent(
                name="S_r^nil",
                source="Jordan nilpotent L0 insertions",
                status="requires fixed nilpotency-length control",
                description="nilpotent Jordan-block insertions",
            ),
            LogarithmicShadowComponent(
                name="S_r^log",
                source="conj:wp-logarithmic-boundary-amplitude-bound",
                status="conjectural",
                description="explicit powers of log(z_i-z_j)",
            ),
        ),
        "finite_zhu_implies_log_bound": False,
        "full_tempering_status": "frontier until regular and logarithmic bounds are proved",
    }


def boundary_changing_correlator_profile() -> dict[str, object]:
    """Return the phi_{0,1} logarithmic amplitude estimate."""

    return {
        "field": "phi_{0,1}",
        "correlator": "<phi_{0,1}(z_1) ... phi_{0,1}(z_r)>",
        "cutoffs": "configuration-space cutoffs and branch sectors fixed",
        "estimate": (
            "abs(correlator) <= C^r (r!)^(1-epsilon) "
            "prod_{i<j}|z_i-z_j|^{-N_ij} "
            "(1 + sum_{i<j}|log|z_i-z_j||)^M"
        ),
        "shadow_equivalent": (
            "limsup_r (|S_r(W(p))|/r!)^(1/r) = 0 and "
            "|S_r^log| <= C^r (r!)^(1-epsilon)"
        ),
        "zhu_dimension_sufficient": False,
    }


@dataclass(frozen=True)
class LogShadowSubfactorialBound:
    """Bound |S_r^log| <= C^r (r!)^(1-epsilon)."""

    base: int
    epsilon: Fraction

    def __post_init__(self) -> None:
        if self.base < 1:
            raise ValueError("base C must be positive")
        if self.epsilon <= 0 or self.epsilon > 1:
            raise ValueError("epsilon must satisfy 0 < epsilon <= 1")

    def log_normalized_root_rate(self, r: int) -> float:
        """Return log(((C^r (r!)^(1-eps))/r!)^(1/r))."""

        if r < 2:
            raise ValueError(f"arity must be >= 2, got {r}")
        return log(self.base) - float(self.epsilon) * lgamma(r + 1) / r

    def normalized_root_rate(self, r: int) -> float:
        """Return ((C^r (r!)^(1-eps))/r!)^(1/r)."""

        return exp(self.log_normalized_root_rate(r))

    def tends_to_zero_certificate(self, r_small: int, r_large: int) -> bool:
        """Check the monotone Stirling trend on two large arities."""

        if r_large <= r_small:
            raise ValueError("r_large must be greater than r_small")
        return (
            self.log_normalized_root_rate(r_large)
            < self.log_normalized_root_rate(r_small)
            and self.log_normalized_root_rate(r_large) < 0
        )


def standard_log_shadow_bound() -> LogShadowSubfactorialBound:
    """Return a representative subfactorial bound profile."""

    return LogShadowSubfactorialBound(base=2, epsilon=Fraction(1, 4))
