"""Irregular KZB boundary formal type and Stokes cocycle profiles."""
from __future__ import annotations

from typing import Dict, Literal, Tuple


LevelLocus = Literal["integrable", "generic_nonrational", "critical"]


def boundary_connection_profile(pole_order: int) -> Dict[str, object]:
    """Return the ramified q-coordinate KZB boundary connection profile."""
    if pole_order < 1:
        raise ValueError("pole_order must be positive")

    polar_terms = tuple(
        f"A_{j}/q^{j}" if j > 1 else "A_1/q"
        for j in range(pole_order, 0, -1)
    )
    return {
        "pole_order": pole_order,
        "chart": "alpha=(q,t_alpha)",
        "connection": (
            "nabla=d-(A_m/q^m+...+A_1/q+A_0)dq-"
            "sum_alpha B_alpha(q) dt_alpha"
        ),
        "polar_terms": polar_terms,
        "regular_q_term": "A_0 dq",
        "transverse_terms": "sum_alpha B_alpha(q) dt_alpha",
        "irregular_requires_higher_pole": pole_order > 1,
    }


def formal_type_profile(pole_order: int) -> Dict[str, object]:
    """Return the formal type Theta_partial(q) for the given pole order."""
    if pole_order < 1:
        raise ValueError("pole_order must be positive")

    terms = tuple(
        f"A_{r + 1}/({r} q^{r})" for r in range(1, pole_order)
    ) + ("A_1 log q",)
    return {
        "pole_order": pole_order,
        "theta": "Theta_partial(q)=sum_{r=1}^{m-1} A_{r+1}/(r q^r)+A_1 log q",
        "terms": terms,
        "regular_singular_when": "pole_order=1",
    }


def formal_gauge_profile(ramification_index: int) -> Dict[str, str]:
    """Return the Levelt-Turrittin/JMU formal gauge shape."""
    if ramification_index <= 0:
        raise ValueError("ramification_index must be positive")

    return {
        "transformation": f"G(q) in Aut(V)[[q^(1/{ramification_index})]]",
        "normal_form": "nabla_q ~= d_q-d_q Theta_partial-Lambda dq/q",
        "category": "formal/sectorial wild category",
    }


def stokes_sector_profile() -> Dict[str, str]:
    """Return the sectorial solution, transition matrix, and cocycle."""
    return {
        "sectorial_solution": "Y_ell(q)=H_ell(q) exp(Theta_partial(q)) q^Lambda",
        "stokes_matrix": "S_{ell,ell'}=Y_ell^{-1}Y_ell'",
        "cocycle": "S_{ell1,ell2} S_{ell2,ell3} S_{ell3,ell1}=1",
        "meaning": "generic nonrational associativity is Stokes cocycle coherence",
    }


def clutching_profile() -> Dict[str, str]:
    """Return the mixed clutching identities for formal types and Stokes factors."""
    return {
        "formal_type": (
            "Theta_{partial,Gamma1#Gamma2}=Theta_{partial,Gamma1}"
            " direct_sum Theta_{partial,Gamma2} direct_sum Theta_node"
        ),
        "stokes_product": (
            "S^{Gamma1#Gamma2}_{ell,ell'}=S^{Gamma1}_{ell,ell'} "
            "S^{Gamma2}_{ell,ell'} S^{node}_{ell,ell'}"
        ),
        "order": "tangential basepoint convention",
    }


def covered_level_locus(kind: LevelLocus) -> Dict[str, str]:
    """Return the exact covered level-locus mechanism."""
    loci = {
        "integrable": {
            "level": "k in Z_{>=0}",
            "mechanism": "KZ pentagon + Kazhdan-Lusztig regularity",
        },
        "generic_nonrational": {
            "level": "k in C\\Q",
            "mechanism": "irregular KZB formal type + wild groupoid",
        },
        "critical": {
            "level": "k=-h^vee",
            "mechanism": "Feigin-Frenkel centre",
        },
    }
    return loci[kind]


def uncovered_level_locus() -> Dict[str, str]:
    """Return the rational nonintegral noncritical proof obligation."""
    return {
        "level": "rational nonintegral noncritical",
        "status": "outside covered loci",
        "required_input": "admissible-level comparison theorem",
    }


def kzb_vs_curved_dunn_profile() -> Dict[str, Tuple[str, ...]]:
    """Separate KZB Stokes composition from curved-Dunn H^2 vanishing."""
    return {
        "kzb_stokes_supplies": (
            "boundary formal type",
            "mixed Stokes matrices",
            "wild-groupoid sector-crossing cocycle",
        ),
        "curved_dunn_h2_requires": (
            "modular-bootstrap H^2 acyclicity",
            "degree-2 comparison map",
            "transport of zero cohomology class",
        ),
    }
