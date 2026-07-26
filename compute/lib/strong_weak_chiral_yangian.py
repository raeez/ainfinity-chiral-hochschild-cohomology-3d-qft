"""Strong versus weak chiral Yangian comparison data."""
from __future__ import annotations

from typing import Dict


def weak_chiral_yangian_profile() -> Dict[str, str]:
    """Return the E1-factorization plus RTT datum."""
    return {
        "datum": "(A,R(z),tensor_z)",
        "relation": "R12(z-w)R13(z)R23(w)=R23(w)R13(z)R12(z-w)",
        "deformation_problem": "RTTDef(A,R)",
    }


def strong_chiral_yangian_profile() -> Dict[str, str]:
    """Return the modular-MC chiral Yangian datum."""
    return {
        "datum": "(A,Delta_z,R(z),S(z),epsilon,{m_k}_{k>=2},Theta_mod)",
        "mc_equation": "d Theta_mod + 1/2[Theta_mod,Theta_mod]=0",
        "deformation_problem": "MC_mod(A)",
    }


def comparison_map_profile() -> Dict[str, str]:
    """Return the required RTT-to-modular comparison map."""
    return {
        "map": "Psi_RTT_to_mod: RTTDef(A,R) -> MC_mod(A)",
        "tangent": "dPsi_RTT_to_mod(dotR(z))=Theta_1",
        "second_obstruction": "ob_2(dotR)=[d_{Theta_0}Theta_2+1/2[Theta_1,Theta_1]] in H_mod^2(A)",
    }


def monodromy_rmatrix_statement_status(
    has_comparison_map: bool,
    obstruction_vanishes: bool,
    monodromy_lands_in_strong_object: bool,
) -> str:
    """Classify whether 'monodromy equals R-matrix' is strong or weak."""
    if has_comparison_map and obstruction_vanishes and monodromy_lands_in_strong_object:
        return "strong modular-MC theorem"
    return "weak RTT/factorization statement only"
