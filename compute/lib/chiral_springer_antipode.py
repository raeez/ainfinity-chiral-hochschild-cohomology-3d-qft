"""Finite witnesses for the chiral Springer antipode package.

The functions here do not construct the derived Steinberg stack or the
conjectural chiral Drinfeld double.  They record the exact formal data that
the manuscript uses as a checkable contract: the complementary chiral
Steinberg carrier, convolution correspondence, orientation-reversal/Verdier
antipode, Hopf identities, and the higher Maurer-Cartan homotopy obstruction.
"""

from __future__ import annotations

from dataclasses import dataclass


HOPF_ANTIPODE_IDENTITIES: tuple[str, ...] = (
    "m(S tensor 1) Delta = eta epsilon",
    "m(1 tensor S) Delta = eta epsilon",
    "Delta S = (S tensor S) Delta^op",
    "epsilon S = epsilon",
)

MC_ORIENTATION_REVERSAL_EQUATION = "S(Theta) = -Theta + d_Theta(Xi)"


@dataclass(frozen=True)
class ChiralSpringerAntipodeData:
    """Formal contract for the complementary chiral Steinberg antipode."""

    steinberg_object: str
    lagrangian_left: str
    lagrangian_right: str
    convolution_space: str
    convolution_formula: str
    point_involution: str
    sheaf_antipode: str
    hopf_identities: tuple[str, ...]
    mc_obstruction_equation: str
    fallback_without_homotopy: str


def chiral_springer_antipode_profile() -> ChiralSpringerAntipodeData:
    """Return the formula package required by the A214 correction."""

    return ChiralSpringerAntipodeData(
        steinberg_object="St^ch(A) = L_A x^h_{M_vac} L_A^vee",
        lagrangian_left="L_A",
        lagrangian_right="L_A^vee (Verdier/Koszul-complementary)",
        convolution_space="L_A x^h_M L_A^vee x^h_M L_A",
        convolution_formula=(
            "F star G = p_13,*(p_12^*F tensor p_23^*G)"
        ),
        point_involution="iota(x, y) = (y, x)",
        sheaf_antipode="S(F) = D_Verdier(iota^*F)",
        hopf_identities=HOPF_ANTIPODE_IDENTITIES,
        mc_obstruction_equation=MC_ORIENTATION_REVERSAL_EQUATION,
        fallback_without_homotopy=(
            "associated graded anti-equivalence, not full chiral Drinfeld antipode"
        ),
    )


def mc_orientation_reversal_status(
    *,
    class_label: str,
    has_homotopy: bool,
) -> dict[str, str | bool]:
    """Classify whether orientation reversal gives a full antipode.

    Class G has no higher MC tower in this finite diagnostic, so the homotopy
    obstruction is vacuous.  For classes L/C/M, the condition
    S(Theta) = -Theta + d_Theta(Xi) is the gate between an associated-graded
    anti-equivalence and a full antipode of the chiral Drinfeld double.
    """

    normalized = class_label.strip().upper()
    if normalized not in {"G", "L", "C", "M"}:
        raise ValueError("class_label must be one of G, L, C, M")

    if normalized == "G":
        return {
            "class": normalized,
            "homotopy_required": False,
            "condition": MC_ORIENTATION_REVERSAL_EQUATION,
            "status": "verified_gaussian_antipode",
            "full_antipode": True,
        }

    if has_homotopy:
        return {
            "class": normalized,
            "homotopy_required": True,
            "condition": MC_ORIENTATION_REVERSAL_EQUATION,
            "status": "conditional_full_antipode",
            "full_antipode": True,
        }

    return {
        "class": normalized,
        "homotopy_required": True,
        "condition": MC_ORIENTATION_REVERSAL_EQUATION,
        "status": "associated_graded_anti_equivalence_only",
        "full_antipode": False,
    }


def antipode_identity_set_complete(identities: tuple[str, ...]) -> bool:
    """Check that the four Hopf antipode identities are present exactly."""

    return tuple(identities) == HOPF_ANTIPODE_IDENTITIES
