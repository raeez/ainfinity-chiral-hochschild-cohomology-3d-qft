"""Deletion ledger for maximal Universal Holography overclaims."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DeletionLedgerEntry:
    """One corrected claim in the A218 deletion ledger."""

    claim: str
    correct_status: str
    replacement: str


DELETION_LEDGER: tuple[DeletionLedgerEntry, ...] = (
    DeletionLedgerEntry(
        claim="SC^{ch,top} recognition is global",
        correct_status="Conditional",
        replacement="H4 product-formal local-shadow recognition",
    ),
    DeletionLedgerEntry(
        claim="Obs^bulk ~= CH^bullet(A,A) for arbitrary logarithmic SC",
        correct_status="NotProved",
        replacement=(
            "requires a chosen HT realisation satisfying H1-H2 and "
            "exact-sector hypotheses"
        ),
    ),
    DeletionLedgerEntry(
        claim="class M raw direct-sum chain theorem",
        correct_status="False",
        replacement="weight-completed/pro or finite-propagation ambient",
    ),
    DeletionLedgerEntry(
        claim="E3-PBW proves concentration",
        correct_status="Conjectured",
        replacement=(
            "conj:H-concentration-via-E3-rigidity is the precise route; "
            "PBW is a filtration input, not the established proof"
        ),
    ),
    DeletionLedgerEntry(
        claim="DS-Hochschild all nilpotents",
        correct_status="False",
        replacement=(
            "principal, hook, and named cover-descent cases; other "
            "nilpotents need case checks"
        ),
    ),
    DeletionLedgerEntry(
        claim="chain-level associator-free mixed structure",
        correct_status="False",
        replacement="Phi-dependent chain representative",
    ),
    DeletionLedgerEntry(
        claim="bar scalar trace = Maloney-Witten sum",
        correct_status="False",
        replacement="bar trace is perturbative seed; modular saddles are extra",
    ),
    DeletionLedgerEntry(
        claim="scalar genus tower determines full tensor channel",
        correct_status="False",
        replacement=(
            "requires channel splitting, scalar diagonalisation, and "
            "non-scalar component data"
        ),
    ),
    DeletionLedgerEntry(
        claim="K3/Class-S comparison automatic",
        correct_status="Conditional",
        replacement="finite Hall-Borcherds gate",
    ),
)


def deletion_ledger_claims() -> tuple[str, ...]:
    """Return the claims corrected by the deletion ledger."""

    return tuple(entry.claim for entry in DELETION_LEDGER)


def deletion_ledger_status_map() -> dict[str, str]:
    """Map each deleted claim to its correct status."""

    return {entry.claim: entry.correct_status for entry in DELETION_LEDGER}


def deletion_ledger_replacement_map() -> dict[str, str]:
    """Map each deleted claim to its replacement theorem form."""

    return {entry.claim: entry.replacement for entry in DELETION_LEDGER}


def maximal_theorem_status_alphabet() -> tuple[str, ...]:
    """The status alphabet allowed in the corrected maximal theorem."""

    return ("ProvedHere", "ProvedElsewhere", "Conditional", "Conjectured")


def false_claims() -> tuple[str, ...]:
    """Return the entries whose corrected status is False."""

    return tuple(
        entry.claim for entry in DELETION_LEDGER
        if entry.correct_status == "False"
    )


def corrected_maximal_theorem_form() -> dict[str, object]:
    """Return the A219 sequent calculus for the maximal theorem."""

    return {
        "sequent": (
            "Xi(A) |-_{A(A)} Phi_hol(A)=T_A; "
            "Obs_boundary(T_A) ~= A; Obs_bulk(T_A) ~= Z_der_ch(A)"
        ),
        "status_alphabet": maximal_theorem_status_alphabet(),
        "ambients": {
            "G/L/C": "Ch(Vect) on stated non-critical loci",
            "M": "Ch_hat_wt_rho or pro-Ch(Vect)",
            "W_infty": "bounded-weight pro-window tower",
        },
        "arbitrary_logarithmic": (
            "algebraic package only unless H1, H2, and exact-sector "
            "comparison are supplied"
        ),
        "critical_affine": "k = -h^vee is outside Phi_hol non-critical source",
        "k3_borcherds": (
            "scalar shadow only unless Hall-Borcherds operator gates are constructed"
        ),
        "raw_direct_sum": "not a theorem in class M",
    }


def ds_nilpotent_scope_profile() -> dict[str, object]:
    """Separate DS-BRST existence from DS-Hochschild transport."""

    return {
        "good_grading_role": (
            "necessary input for forming the Kazhdan DS-BRST complex"
        ),
        "good_grading_sufficient_for_ds_hochschild": False,
        "blanket_all_nilpotents_transport": False,
        "cohomological_topologization": (
            "good-graded DS data use branched-cover/Galois descent where "
            "the Costello-Gaiotto boundary theorem applies"
        ),
        "hochschild_transport_verified": (
            "principal",
            "hook",
            "Bershadsky-Polyakov cover descent",
        ),
        "case_by_case": (
            "subregular outside checked hooks",
            "minimal outside named cover-descent fibres",
            "exceptional good-graded nilpotents",
            "exotic nilpotents",
        ),
        "required_transport_gates": (
            "non-critical level",
            "DS HPL special deformation retract",
            "mu_q-invariant transferred braces",
            "mu_q-invariant transported antighost",
            "Hochschild comparison",
        ),
    }


def associator_chain_scope_profile() -> dict[str, object]:
    """Separate cohomological associator-independence from chain data."""

    return {
        "chain_level_associator_free": False,
        "cohomological_associator_free": True,
        "bar_side_associator_free_invariants": (
            "kappa",
            "shadow tower",
            "Koszul depth",
        ),
        "chain_level_object": "Phi-dependent representative",
        "torsor": "GRT_1(Q) Drinfeld-associator torsor",
        "deleted_form": "associator-independent chain-level collapse",
    }


def maloney_witten_bridge_scope_profile() -> dict[str, object]:
    """Separate the raw bar trace from the Maloney-Witten orbit sum."""

    return {
        "raw_bar_trace_equals_maloney_witten_sum": False,
        "raw_bar_trace_role": "chain-level perturbative seed",
        "maloney_witten_object": "completed modular orbit over saddles",
        "requires": (
            "Borel summability",
            "Zwegers completion",
            "Brown-Henneaux dictionary",
            "saddle labelling",
            "ensemble prescription",
        ),
        "heisenberg_orbit_scope": (
            "Rademacher orbit object, not a pure-gravity path integral "
            "without extra saddle interpretation"
        ),
    }


def arbitrary_logarithmic_bulk_scope_profile() -> dict[str, object]:
    """Scope of bulk/Hochschild comparison for logarithmic SC data."""

    return {
        "abstract_log_sc_has_physical_bulk": False,
        "abstract_log_sc_package": "algebraic local-shadow package",
        "bulk_observable_identification_requires": (
            "chosen HT prefactorization realization",
            "product-formal local shadow",
            "H1-H2 physics bridge",
            "boundary-linear exact-sector comparison",
        ),
        "bulk_hochschild_scope": (
            "HT prefactorization realization with exact-sector hypotheses"
        ),
    }


def recognition_scope_profile() -> dict[str, object]:
    """Scope of the SC^{ch,top} recognition theorem."""

    return {
        "global_recognition_theorem": False,
        "recognized_surface": "product-formal local-shadow rectangle data",
        "not_recognized": (
            "arbitrary global Ran-space factorization data",
            "physical HT bulk without chosen prefactorization realization",
        ),
        "output": "local C_*(W(SC^{ch,top})) algebra shadow",
    }


def pbw_concentration_scope_profile() -> dict[str, object]:
    """Separate proved concentration from the conjectural E3-PBW route."""

    return {
        "e3_pbw_proves_concentration": False,
        "established_mechanism": "Arnold-Orlik-Solomon/FM local proof",
        "pbw_role": "filtration and associated-graded input",
        "proved_route": (
            "ordered bar",
            "chiral PBW/Koszul collapse",
            "Arnold-Orlik-Solomon/FM concentration",
        ),
        "conjectural_route_requires": (
            "filtered E3 envelope",
            "Free_E3 associated graded",
            "Rees-flat no-hidden-extension lift",
            "convergent PBW spectral sequence",
            "polynomial-growth/amplitude bounds",
            "E1-page support in total degrees <= 2",
        ),
    }


def scalar_tensor_channel_scope_profile() -> dict[str, object]:
    """Separate scalar genus traces from full tensor-channel data."""

    return {
        "scalar_genus_tower_determines_full_tensor_channel": False,
        "scalar_tower_scope": (
            "uniform-weight closed-sector trace of the modular bar curvature"
        ),
        "rank_one_exception": (
            "Heisenberg/Gaussian rank-one abelian channel, where the tensor "
            "channel has no mixed entries"
        ),
        "full_tensor_channel_requires": (
            "chosen channel splitting",
            "scalar diagonalisation",
            "basis of conformal blocks",
            "off-diagonal stable-graph component integrals",
            "non-scalar ordered or field-valued coefficients",
        ),
        "trace_forgets": (
            "mixed tensor entries",
            "ordered depth filtration",
            "field-valued coefficients",
            "off-diagonal cross-channel corrections",
        ),
    }
