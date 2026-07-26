from __future__ import annotations

import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.deletion_ledger import (
    DELETION_LEDGER,
    associator_chain_scope_profile,
    arbitrary_logarithmic_bulk_scope_profile,
    corrected_maximal_theorem_form,
    deletion_ledger_claims,
    deletion_ledger_replacement_map,
    deletion_ledger_status_map,
    ds_nilpotent_scope_profile,
    false_claims,
    maloney_witten_bridge_scope_profile,
    maximal_theorem_status_alphabet,
    pbw_concentration_scope_profile,
    recognition_scope_profile,
    scalar_tensor_channel_scope_profile,
)


def test_deletion_ledger_has_pdf_rows():
    assert deletion_ledger_claims() == (
        "SC^{ch,top} recognition is global",
        "Obs^bulk ~= CH^bullet(A,A) for arbitrary logarithmic SC",
        "class M raw direct-sum chain theorem",
        "E3-PBW proves concentration",
        "DS-Hochschild all nilpotents",
        "chain-level associator-free mixed structure",
        "bar scalar trace = Maloney-Witten sum",
        "scalar genus tower determines full tensor channel",
        "K3/Class-S comparison automatic",
    )
    assert len(DELETION_LEDGER) == 9


def test_deletion_ledger_statuses_and_replacements():
    statuses = deletion_ledger_status_map()
    replacements = deletion_ledger_replacement_map()

    assert statuses["SC^{ch,top} recognition is global"] == "Conditional"
    assert statuses["class M raw direct-sum chain theorem"] == "False"
    assert statuses["E3-PBW proves concentration"] == "Conjectured"
    assert statuses["K3/Class-S comparison automatic"] == "Conditional"
    assert "product-formal local-shadow" in replacements[
        "SC^{ch,top} recognition is global"
    ]
    assert "H1-H2" in replacements[
        "Obs^bulk ~= CH^bullet(A,A) for arbitrary logarithmic SC"
    ]
    assert "weight-completed/pro" in replacements[
        "class M raw direct-sum chain theorem"
    ]
    assert "conj:H-concentration-via-E3-rigidity" in replacements[
        "E3-PBW proves concentration"
    ]
    assert "named cover-descent" in replacements[
        "DS-Hochschild all nilpotents"
    ]
    assert "modular saddles" in replacements[
        "bar scalar trace = Maloney-Witten sum"
    ]
    assert "non-scalar component data" in replacements[
        "scalar genus tower determines full tensor channel"
    ]


def test_false_claims_are_exactly_the_deleted_overclaims():
    assert false_claims() == (
        "class M raw direct-sum chain theorem",
        "DS-Hochschild all nilpotents",
        "chain-level associator-free mixed structure",
        "bar scalar trace = Maloney-Witten sum",
        "scalar genus tower determines full tensor channel",
    )
    assert maximal_theorem_status_alphabet() == (
        "ProvedHere",
        "ProvedElsewhere",
        "Conditional",
        "Conjectured",
    )


def test_corrected_maximal_theorem_form_records_ambient_sequent():
    profile = corrected_maximal_theorem_form()

    assert profile["sequent"] == (
        "Xi(A) |-_{A(A)} Phi_hol(A)=T_A; "
        "Obs_boundary(T_A) ~= A; Obs_bulk(T_A) ~= Z_der_ch(A)"
    )
    assert profile["status_alphabet"] == maximal_theorem_status_alphabet()
    assert profile["ambients"] == {
        "G/L/C": "Ch(Vect) on stated non-critical loci",
        "M": "Ch_hat_wt_rho or pro-Ch(Vect)",
        "W_infty": "bounded-weight pro-window tower",
    }
    assert "H1, H2, and exact-sector" in profile["arbitrary_logarithmic"]
    assert profile["critical_affine"] == (
        "k = -h^vee is outside Phi_hol non-critical source"
    )
    assert profile["raw_direct_sum"] == "not a theorem in class M"
    assert "Hall-Borcherds operator gates" in profile["k3_borcherds"]


def test_concentration_route_status_is_conjectural_in_source():
    root = Path(__file__).resolve().parents[2]
    climax = (root / "chapters/connections/programme_climax_platonic.tex").read_text()
    chd = (root / "chapters/theory/chiral_higher_deligne.tex").read_text()
    flat_climax = " ".join(climax.split())

    assert r"\ClaimStatusConjectured &" in climax
    assert r"Conjecture~\ref{conj:H-concentration-via-E3-rigidity}" in climax
    assert r"\ClaimStatusNeedsVerification" not in climax
    assert r"\{\ClaimStatusProvedHere,\ \ClaimStatusProvedElsewhere, \ClaimStatusConditional,\ \ClaimStatusConjectured\}" in flat_climax

    assert r"\label{conj:H-concentration-via-E3-rigidity}" in chd
    assert r"\label{thm:H-concentration-via-E3-rigidity}" not in chd


def test_ds_nilpotent_scope_separates_good_grading_from_transport():
    profile = ds_nilpotent_scope_profile()

    assert "Kazhdan DS-BRST complex" in profile["good_grading_role"]
    assert profile["good_grading_sufficient_for_ds_hochschild"] is False
    assert profile["blanket_all_nilpotents_transport"] is False
    assert profile["hochschild_transport_verified"] == (
        "principal",
        "hook",
        "Bershadsky-Polyakov cover descent",
    )
    assert "subregular outside checked hooks" in profile["case_by_case"]
    assert "mu_q-invariant transferred braces" in profile[
        "required_transport_gates"
    ]


def test_ds_nilpotent_scope_manuscript_source_guard():
    root = Path(__file__).resolve().parents[2]
    axioms = (root / "chapters/theory/axioms.tex").read_text()
    chd = (root / "chapters/theory/chiral_higher_deligne.tex").read_text()
    gravity = (root / "chapters/connections/3d_gravity.tex").read_text()
    flat_axioms = " ".join(axioms.split())
    flat_chd = " ".join(chd.split())
    flat_gravity = " ".join(gravity.split())

    assert "This is the scope for forming the DS-BRST complex" in flat_axioms
    assert "not a DS--Hochschild transport theorem" in flat_axioms
    assert "not a chain-level $\\SCchtop$ lift" in flat_axioms
    assert "not a Koszul-preservation statement" in flat_axioms
    assert "principal and hook-type data" in flat_axioms
    assert "Bershadsky--Polyakov checked separately" in flat_axioms
    assert "case-by-case BRST admissibility" in flat_axioms
    assert "principal or hook-type" in flat_chd
    assert "Exotic nilpotents remain open" in flat_chd
    assert "without the cover descent datum" in flat_chd
    assert "This is a cohomological topologization statement" in flat_gravity
    assert "stronger DS--Hochschild chain-level transport" in flat_gravity


def test_associator_scope_separates_cohomology_from_chain_level():
    profile = associator_chain_scope_profile()

    assert profile["chain_level_associator_free"] is False
    assert profile["cohomological_associator_free"] is True
    assert profile["bar_side_associator_free_invariants"] == (
        "kappa",
        "shadow tower",
        "Koszul depth",
    )
    assert profile["chain_level_object"] == "Phi-dependent representative"
    assert "GRT_1" in profile["torsor"]
    assert profile["deleted_form"] == (
        "associator-independent chain-level collapse"
    )


def test_associator_scope_manuscript_source_guard():
    root = Path(__file__).resolve().parents[2]
    synthesis = (root / "chapters/frame/part_viii_synthesis.tex").read_text()
    chd = (root / "chapters/theory/chiral_higher_deligne.tex").read_text()
    flat_synthesis = " ".join(synthesis.split())
    flat_chd = " ".join(chd.split())

    assert (
        "associator-independent chain-level formulation is not a theorem"
        in flat_synthesis
    )
    assert "non-trivial $\\mathrm{GRT}_{1}(\\Q)$-torsor" in flat_synthesis
    assert "associator-free statement lives on cohomology" in flat_synthesis
    assert "chain-level mixed object is a \\(\\Phi\\)-dependent representative" in flat_synthesis
    assert "\\(\\Phi\\)-dependent chain-level Deligne--Tamarkin torsor" in flat_synthesis
    assert "associator-independent chain-level collapse is false" in flat_synthesis
    assert "The cohomological derived centre is associator-free" in flat_chd
    assert "The chain-level derived centre depends on $\\Phi$" in flat_chd


def test_pbw_concentration_scope_separates_routes():
    profile = pbw_concentration_scope_profile()

    assert profile["e3_pbw_proves_concentration"] is False
    assert profile["established_mechanism"] == (
        "Arnold-Orlik-Solomon/FM local proof"
    )
    assert profile["pbw_role"] == "filtration and associated-graded input"
    assert profile["proved_route"] == (
        "ordered bar",
        "chiral PBW/Koszul collapse",
        "Arnold-Orlik-Solomon/FM concentration",
    )
    assert profile["conjectural_route_requires"] == (
        "filtered E3 envelope",
        "Free_E3 associated graded",
        "Rees-flat no-hidden-extension lift",
        "convergent PBW spectral sequence",
        "polynomial-growth/amplitude bounds",
        "E1-page support in total degrees <= 2",
    )


def test_pbw_concentration_manuscript_source_guard():
    root = Path(__file__).resolve().parents[2]
    concordance = (root / "chapters/connections/concordance.tex").read_text()
    line_ext = (
        root / "chapters/connections/thqg_line_operators_extensions.tex"
    ).read_text()
    line_ops = (root / "chapters/connections/line-operators.tex").read_text()
    hochschild = (root / "chapters/connections/hochschild.tex").read_text()
    anomaly = (
        root / "chapters/connections/thqg_anomaly_extensions.tex"
    ).read_text()
    flat_concordance = " ".join(concordance.split())
    flat_line_ext = " ".join(line_ext.split())
    flat_line_ops = " ".join(line_ops.split())
    flat_hochschild = " ".join(hochschild.split())
    flat_anomaly = " ".join(anomaly.split())

    assert (
        "PBW/Koszul collapse is the input; the vanishing mechanism is "
        "ordered-bar Arnold--Orlik--Solomon/FM concentration"
        in flat_concordance
    )
    assert "MC1 (PBW/Koszul input + ordered-bar concentration)" in flat_line_ext
    assert "DK-2 & MC1 (PBW/Koszul input) + DK-1" in flat_line_ext
    assert "PBW/Koszul concentration input alone does not" in flat_line_ops
    assert "not a PBW-only proof of concentration" in flat_hochschild
    assert "PBW/Koszul filtration input plus ordered-bar" in flat_anomaly
    assert "PBW concentration for all standard families" not in flat_concordance
    assert "MC1 (PBW concentration)" not in flat_line_ext


def test_scalar_tensor_channel_scope_separates_trace_from_tensor():
    profile = scalar_tensor_channel_scope_profile()

    assert profile["scalar_genus_tower_determines_full_tensor_channel"] is False
    assert "closed-sector trace" in profile["scalar_tower_scope"]
    assert "rank-one abelian channel" in profile["rank_one_exception"]
    assert profile["full_tensor_channel_requires"] == (
        "chosen channel splitting",
        "scalar diagonalisation",
        "basis of conformal blocks",
        "off-diagonal stable-graph component integrals",
        "non-scalar ordered or field-valued coefficients",
    )
    assert "mixed tensor entries" in profile["trace_forgets"]
    assert "off-diagonal cross-channel corrections" in profile["trace_forgets"]


def test_scalar_tensor_channel_manuscript_source_guard():
    root = Path(__file__).resolve().parents[2]
    rosetta = (root / "chapters/examples/rosetta_stone.tex").read_text()
    gravity = (root / "chapters/connections/3d_gravity.tex").read_text()
    movements = (
        root / "chapters/connections/thqg_3d_gravity_movements_vi_x.tex"
    ).read_text()
    reconstruction = (
        root / "chapters/connections/thqg_holographic_reconstruction.tex"
    ).read_text()
    climax = (
        root / "chapters/connections/programme_climax_platonic.tex"
    ).read_text()
    tensor_d = (
        root / "chapters/theory/theorems_C_D_native_vol2_platonic.tex"
    ).read_text()
    flat_rosetta = " ".join(rosetta.split())
    flat_gravity = " ".join(gravity.split())
    flat_movements = " ".join(movements.split())
    flat_reconstruction = " ".join(reconstruction.split())
    flat_climax = " ".join(climax.split())
    flat_tensor_d = " ".join(tensor_d.split())

    assert "rank-one abelian channel" in flat_rosetta
    assert "does not extend to a multi-channel tensor claim" in flat_rosetta
    assert (
        "closed scalar observable is a "
        "\\(\\kappaChHodge(\\mathrm{Vir}_c)\\)-polynomial"
        in flat_gravity
    )
    assert "not a reconstruction of the full tensor channel" in flat_gravity
    assert "Read as the complete scalar genus tower" in flat_movements
    assert "rank-one Gaussian shadow" in flat_reconstruction
    assert "not the full tensor-channel package" in flat_reconstruction
    assert "not the full tensor-channel object" in flat_climax
    assert "mixed entries $K_{TW}$ and $K_{WT}$ are invisible" in flat_tensor_d


def test_maloney_witten_scope_separates_trace_from_orbit_sum():
    profile = maloney_witten_bridge_scope_profile()

    assert profile["raw_bar_trace_equals_maloney_witten_sum"] is False
    assert profile["raw_bar_trace_role"] == "chain-level perturbative seed"
    assert profile["maloney_witten_object"] == "completed modular orbit over saddles"
    assert profile["requires"] == (
        "Borel summability",
        "Zwegers completion",
        "Brown-Henneaux dictionary",
        "saddle labelling",
        "ensemble prescription",
    )
    assert "not a pure-gravity path integral" in profile["heisenberg_orbit_scope"]


def test_maloney_witten_bridge_manuscript_source_guard():
    root = Path(__file__).resolve().parents[2]
    gravity = (root / "chapters/connections/3d_gravity.tex").read_text()
    climax = (root / "chapters/connections/programme_climax_platonic.tex").read_text()
    flat_gravity = " ".join(gravity.split())
    flat_climax = " ".join(climax.split())

    assert (
        "Conditional orbit-sum comparison with the Maloney--Witten partition function"
        in flat_gravity
    )
    assert "not an assertion that the raw trace" in flat_gravity
    assert "the Maloney--Witten object is the completed modular orbit" in flat_gravity
    assert "Brown--Henneaux, saddle-labelling, and ensemble-prescription" in flat_gravity
    assert "not a pure-gravity path integral" in flat_gravity
    assert "The bar trace is a perturbative seed" in flat_climax
    assert "separate orbit sum, ensemble prescription, and dominance data" in flat_climax


def test_arbitrary_logarithmic_bulk_scope_requires_realization():
    profile = arbitrary_logarithmic_bulk_scope_profile()

    assert profile["abstract_log_sc_has_physical_bulk"] is False
    assert profile["abstract_log_sc_package"] == "algebraic local-shadow package"
    assert profile["bulk_observable_identification_requires"] == (
        "chosen HT prefactorization realization",
        "product-formal local shadow",
        "H1-H2 physics bridge",
        "boundary-linear exact-sector comparison",
    )
    assert "exact-sector" in profile["bulk_hochschild_scope"]


def test_arbitrary_logarithmic_bulk_scope_manuscript_source_guard():
    root = Path(__file__).resolve().parents[2]
    hochschild = (root / "chapters/connections/hochschild.tex").read_text()
    brace = (root / "chapters/connections/brace.tex").read_text()
    flat_hochschild = " ".join(hochschild.split())
    flat_brace = " ".join(brace.split())

    assert "Package content and scope" in flat_hochschild
    assert "chosen HT prefactorization realization" in flat_hochschild
    assert "boundary-linear exact-sector hypotheses" in flat_hochschild
    assert "abstract logarithmic $\\SCchtop$-algebra without such a realization" in flat_hochschild
    assert "does not determine a physical bulk observable complex" in flat_hochschild
    assert "only after a chosen HT prefactorization realization" in flat_brace
    assert "product-formal local shadow" in flat_brace


def test_standard_family_bulk_computation_setup_is_scoped():
    root = Path(__file__).resolve().parents[2]
    source = (root / "chapters/connections/hochschild.tex").read_text()
    flat = " ".join(source.split())

    assert (
        "The algebraic object computed below is the chiral Hochschild "
        "cochain complex"
        in flat
    )
    assert (
        "When a standard family is realized by a chosen 3d HT "
        "prefactorization model"
        in flat
    )
    assert (
        "Theorem~\\ref{thm:bulk_hochschild} identifies that complex with "
        "the corresponding physical bulk observables"
        in flat
    )
    assert (
        "For each standard family, $C^*(\\cA,\\cA)$ is computed and "
        "identified with the bulk observables."
        not in source
    )


def test_bulk_hochschild_summary_echoes_are_scoped():
    root = Path(__file__).resolve().parents[2]
    files = {
        "conclusion": root / "chapters/connections/conclusion.tex",
        "examples_computing": root / "chapters/examples/examples-computing.tex",
        "examples_worked": root / "chapters/examples/examples-worked.tex",
        "introduction": root / "chapters/theory/introduction.tex",
        "raviolo": root / "chapters/theory/raviolo.tex",
    }
    text = {name: path.read_text() for name, path in files.items()}
    flat = {name: " ".join(value.split()) for name, value in text.items()}

    assert (
        "the algebraic closed-sector object is the chiral Hochschild complex"
        in flat["conclusion"]
    )
    assert (
        "physical bulk observables identify with it only under the chosen "
        "HT prefactorisation and boundary-linear exact-sector comparison"
        in flat["conclusion"]
    )
    assert (
        "Under the same field-theoretic realization and exact-sector scope"
        in flat["examples_computing"]
    )
    assert (
        "This is the universal closed-sector complex of the boundary algebra"
        in flat["examples_worked"]
    )
    assert (
        "only after the HT prefactorisation and exact-sector comparison"
        in flat["examples_worked"]
    )
    assert (
        "computes the algebraic closed-sector Hochschild complex"
        in flat["introduction"]
    )
    assert (
        "the physical bulk identification is the scoped HT "
        "prefactorisation/exact-sector comparison"
        in flat["introduction"]
    )
    assert (
        "Under the HT prefactorisation and boundary-linear exact-sector "
        "hypotheses of Theorem~\\ref{thm:bulk_hochschild}"
        in flat["raviolo"]
    )
    assert (
        "relates physical bulk observables, under the comparison hypotheses"
        in flat["raviolo"]
    )
    stale_forms = (
        "bulk observables identify with chiral Hochschild cochains",
        "identifies bulk observables with chiral Hochschild cochains",
        "identifies bulk observables with Hochschild cochains",
        "The bulk observables are chiral Hochschild cochains",
        "Section~\\ref{sec:chiral_hochschild_expanded} relates bulk "
        "observables to chiral Hochschild cochains.",
    )
    for value in text.values():
        for stale in stale_forms:
            assert stale not in value


def test_bulk_algebra_computation_headings_are_closed_sector():
    root = Path(__file__).resolve().parents[2]
    files = {
        "main": root / "main.tex",
        "concordance": root / "chapters/connections/concordance.tex",
        "hochschild": root / "chapters/connections/hochschild.tex",
        "introduction": root / "chapters/theory/introduction.tex",
        "w3": root / "chapters/examples/w-algebras-w3.tex",
        "stable": root / "chapters/examples/w-algebras-stable.tex",
    }
    text = {name: path.read_text() for name, path in files.items()}
    flat = {name: " ".join(value.split()) for name, value in text.items()}
    theorem_start = text["hochschild"].index(
        "\\subsection{The bulk--boundary factorization theorem}"
    )
    hochschild_computation_band = text["hochschild"][
        text["hochschild"].index(
            "\\subsection{Closed sector from boundary: chiral Hochschild computations}"
        ):
        theorem_start
    ]
    hochschild_assembly_band = text["hochschild"][
        theorem_start:
        text["hochschild"].index("\\end{itemize}", theorem_start)
    ]
    repaired_surface = (
        hochschild_computation_band
        + "\n"
        + hochschild_assembly_band
        + "\n"
        + text["main"]
        + "\n"
        + text["concordance"]
        + "\n"
        + text["introduction"]
        + "\n"
        + text["w3"]
        + "\n"
        + text["stable"]
    )

    assert "Closed sector from boundary: chiral Hochschild computations" in text[
        "hochschild"
    ]
    assert (
        "Chiral Hochschild closed sector for Heisenberg" in text["hochschild"]
    )
    assert (
        "Chiral Hochschild closed sector for affine "
        "$\\widehat{\\mathfrak{sl}}_2$"
        in text["hochschild"]
    )
    assert (
        "Chiral Hochschild closed sector for Virasoro" in text["hochschild"]
    )
    assert (
        "Chiral Hochschild closed sector for $\\mathcal W_3$"
        in text["hochschild"]
    )
    assert (
        "Chiral Hochschild closed sector for $\\mathcal W_4$"
        in text["hochschild"]
    )
    assert (
        "fixed-level Hochschild closed-sector state algebra"
        in flat["hochschild"]
    )
    assert (
        "closed-sector-to-boundary comparison map" in flat["hochschild"]
    )
    assert (
        "physical bulk $\\simeq$ chiral Hochschild cochains only in the "
        "scoped HT prefactorization and boundary-linear exact-sector "
        "comparison"
        in flat["hochschild"]
    )
    assert (
        "algebraic closed sector $\\simeq$ chiral Hochschild; physical bulk "
        "only after the scoped comparison"
        in flat["hochschild"]
    )
    assert (
        "Algebraic closed sector $\\simeq$ chiral Hochschild cochains; "
        "physical bulk only after the scoped HT comparison"
        in flat["introduction"]
    )
    assert (
        "Physical bulk $\\simeq$ chiral Hochschild only under the comparison "
        "hypotheses"
        in flat["concordance"]
    )
    assert (
        "chiral Hochschild closed sector; physical bulk scoped"
        in text["main"]
    )
    assert (
        "Chiral Hochschild closed sector for $\\mathcal{W}_N$"
        in text["w3"]
    )
    assert (
        "Chiral Hochschild closed sector for $\\mathcal{W}_N$"
        in text["stable"]
    )

    stale_forms = (
        "Bulk algebra for Heisenberg",
        "Bulk algebra for affine",
        "Bulk algebra for Virasoro",
        "Bulk algebra for $\\mathcal W_3$",
        "Bulk algebra for $\\mathcal W_4$",
        "Bulk algebra for $\\mathcal{W}_N$",
        "fixed-level bulk state algebra",
        "The bulk algebra at the Hochschild level",
        "The bulk algebra is",
        "The bulk algebra jumps",
        "bulk deformation directions",
        "\\item Bulk $\\simeq$ chiral Hochschild cochains",
        "Bulk $\\simeq$ chiral Hochschild",
        "chiral Hochschild = bulk",
    )
    for stale in stale_forms:
        assert stale not in repaired_surface


def test_lagrangian_hkr_bulk_is_not_ambient_restriction():
    root = Path(__file__).resolve().parents[2]
    hochschild = (root / "chapters/connections/hochschild.tex").read_text()
    remark = hochschild[
        hochschild.index("\\begin{remark}[Lagrangian self-intersection interpretation]"):
        hochschild.index("\\end{remark}", hochschild.index(
            "\\begin{remark}[Lagrangian self-intersection interpretation]"
        ))
    ]
    flat = " ".join(remark.split())

    assert (
        "\\cO\\bigl(\\widehat{\\cM_{\\mathrm{vac}}}_{\\cL_b}\\bigr)_{\\mathrm{Darboux}}"
        in remark
    )
    assert "It is not ordinary restriction" in flat
    assert "\\cO(\\cM_{\\mathrm{vac}})|_{\\cL_b}=\\cO(\\cL_b)" in flat
    assert "boundary chart" in flat
    assert (
        "the bulk cochain algebra keeps the shifted cotangent fibre"
        in flat
    )
    assert (
        "functions on the formal Lagrangian self-intersection"
        in flat
    )
    stale_forms = (
        "\\cO(\\cM_{\\mathrm{vac}})\\big|_{\\cL_b}",
        "The bulk algebra is the ambient",
        "ambient $(-2)$-shifted symplectic stack seen from the Lagrangian",
        "this self-intersection is the ambient stack restricted to $\\cL_b$",
    )
    for stale in stale_forms:
        assert stale not in remark


def test_affine_cs_current_algebra_is_boundary_not_bulk():
    root = Path(__file__).resolve().parents[2]
    files = {
        "spectral_core": root / "chapters/connections/spectral-braiding-core.tex",
        "spectral_split": root / "chapters/connections/spectral-braiding.tex",
        "examples_proved": root / "chapters/examples/examples-complete-proved.tex",
        "examples_full": root / "chapters/examples/examples-complete.tex",
        "examples_core": root / "chapters/examples/examples-complete-core.tex",
    }
    text = {name: path.read_text() for name, path in files.items()}
    flat = {name: " ".join(value.split()) for name, value in text.items()}

    assert (
        "the boundary current chart is the Heisenberg/affine "
        "Kac--Moody algebra $\\widehat{\\mathfrak{u}(1)}_k$"
        in flat["spectral_core"]
    )
    assert (
        "the associated physical bulk object is the chiral derived centre"
        in flat["spectral_split"]
    )
    assert (
        "this kernel comes from the boundary current OPE, not from an "
        "identification of the boundary chart with bulk"
        in flat["spectral_core"]
    )
    for key in ("examples_proved", "examples_full", "examples_core"):
        assert (
            "\\mathcal{A}_{\\partial} = \\text{(affine Kac--Moody PVA "
            "for $\\mathfrak{g}$ at level $k$)}"
            in text[key]
        )
        assert "as the boundary current chart" in flat[key]
        assert (
            "the associated bulk object is the chiral derived centre"
            in flat[key]
        )

    repaired_surface = "\n".join(text.values())
    stale_forms = (
        "the bulk algebra is the affine Kac--Moody",
        "the bulk algebra is the affine Kac--Moody "
        "$\\widehat{\\mathfrak{u}(1)}_k$",
        "\\mathcal{A}_{\\text{bulk}} = \\text{(affine Kac--Moody PVA "
        "for $\\mathfrak{g}$ at level $k$)}",
    )
    for stale in stale_forms:
        assert stale not in repaired_surface


def test_n4_symmetry_matching_is_obstruction_gated():
    root = Path(__file__).resolve().parents[2]
    files = {
        "frontier": root / "chapters/connections/celestial_holography_frontier.tex",
        "split": root / "chapters/connections/celestial_holography.tex",
    }
    text = {name: path.read_text() for name, path in files.items()}
    flat = {name: " ".join(value.split()) for name, value in text.items()}

    assert (
        "This is a perturbative coalgebra-level matching" in flat["frontier"]
    )
    assert (
        "It becomes a bulk-as-derived-centre statement only after the mixed "
        "bulk--boundary comparison class"
        in flat["frontier"]
    )
    assert (
        "\\gamma_{\\cN=4}\\colon \\cA_{\\mathrm{bulk}} \\longrightarrow "
        "Z_{\\mathrm{der}}(\\cA_{\\partial})"
        in flat["frontier"]
    )
    assert (
        "before that vanishing it is comparison data, not an identification"
        in flat["frontier"]
    )
    assert (
        "\\gamma_{\\cN=4}\\colon \\A_{\\mathrm{bulk}} \\longrightarrow "
        "Z_{\\mathrm{der}}(\\A_{\\partial})"
        in flat["split"]
    )
    assert (
        "It becomes a quasi-isomorphism only after the mixed obstruction class"
        in flat["split"]
    )

    repaired_surface = "\n".join(text.values())
    stale_forms = (
        "The bulk algebra itself is the chiral derived centre; the displayed "
        "coalgebra is its perturbative BV resolution.",
        "The global symmetry matching is the bulk--boundary--line triangle",
        "\\cA_{\\mathrm{bulk}} \\;\\simeq\\; Z_{\\mathrm{der}}(\\cA_{\\partial}).",
        "\\A_{\\mathrm{bulk}} \\;\\simeq\\; Z_{\\mathrm{der}}(\\A_{\\partial}).",
        "the global symmetry matching is the statement",
    )
    for stale in stale_forms:
        assert stale not in repaired_surface


def test_boundary_linear_bulk_line_is_morita_scoped():
    root = Path(__file__).resolve().parents[2]
    files = {
        "core": root / "chapters/connections/ht_bulk_boundary_line_core.tex",
        "split": root / "chapters/connections/ht_bulk_boundary_line.tex",
        "frontier": root / "chapters/connections/ht_bulk_boundary_line_frontier.tex",
    }
    text = {name: path.read_text() for name, path in files.items()}
    flat = {name: " ".join(value.split()) for name, value in text.items()}

    for key in ("core", "split"):
        assert (
            "the bulk dg algebra \\(\\cO(\\dCrit(W))\\) is canonically "
            "quasi-isomorphic to \\(\\Zder(B_{L,W})\\)"
            in flat[key]
        )
        assert "not a literal equality of presentations" in flat[key]
        assert (
            "A pointed line algebra such as \\(K_\\kappa\\) is a Morita "
            "presentation of the local line category"
            in flat[key]
        )
        assert (
            "the model bulk algebra is compared with the boundary "
            "derived centre in the derived Morita category"
            in flat[key]
        )
        assert (
            "not reconstruction of an arbitrary physical bulk from a "
            "boundary algebra alone"
            in flat[key]
        )
        assert (
            "the completed model bulk algebra is identified, by Morita "
            "comparison, with \\(\\HH^\\bullet(K_\\kappa)\\)"
            in flat[key]
        )

    for key in files:
        if key in {"core", "split"}:
            assert (
                "model bulk compared with the shifted cotangent of the "
                "derived boundary zero locus"
                in flat[key]
            )
        assert "derived-centre/Morita comparison" in flat[key]

    repaired_surface = "\n".join(text.values())
    stale_forms = (
        "In particular the derived center of the boundary algebra equals "
        "the bulk algebra in the boundary-linear sector.",
        "The line algebra does not equal the bulk algebra; the bulk is "
        "the derived center of the line algebra",
        "the bulk is the derived center of the line algebra",
        "Boundary-linear LG theorem: bulk $=$ shifted cotangent",
        "the bulk is recovered from the boundary algebra",
        "bulk is recovered from lines by Hochschild cochains",
        "the local bulk algebra is recovered from $K_\\kappa$ by Hochschild cochains",
        "the completed local bulk algebra is computed by "
        "\\(\\HH^\\bullet(K_\\kappa)\\)",
    )
    for stale in stale_forms:
        assert stale not in repaired_surface


def test_examples_worked_derived_centre_bulk_is_scoped():
    root = Path(__file__).resolve().parents[2]
    examples = (root / "chapters/examples/examples-worked.tex").read_text()
    flat = " ".join(examples.split())

    assert (
        "is the closed-sector derived-centre vertex. Its identification "
        "with the universal physical bulk algebra uses the "
        "Kodaira--Spencer HT realisation and the scoped bulk--Hochschild "
        "comparison"
        in flat
    )
    assert (
        "In this scoped 5d HT realisation, bulk observables are represented "
        "by the chiral derived center"
        in flat
    )
    assert (
        "On the bulk--Hochschild comparison surface it is represented by "
        "the chiral derived center"
        in flat
    )

    stale_forms = (
        "computes the universal bulk algebra",
        "The bulk algebra is the chiral derived center",
        "The bulk algebra is \\emph{not} the bar complex",
        "It is the chiral derived center",
    )
    for stale in stale_forms:
        assert stale not in examples


def test_gravity_derived_centre_bulk_is_brown_henneaux_scoped():
    root = Path(__file__).resolve().parents[2]
    gravity = (root / "chapters/connections/3d_gravity.tex").read_text()
    flat = " ".join(gravity.split())

    assert (
        "is the algebraic closed-colour vertex attached to the Virasoro "
        "boundary chart"
        in flat
    )
    assert (
        "Its physical gravitational reading uses the Virasoro HT "
        "prefactorization realisation, the Brown--Henneaux boundary chart "
        "\\(c=3\\ell/(2G_N)\\), and the \\(\\Sigma_n\\)-closed-sector descent"
        in flat
    )
    assert (
        "In pure 3d gravity this reading further passes to the "
        "central-extension quotient \\(\\C[\\![c]\\!]\\)"
        in flat
    )
    assert "it is not the whole derived centre" in flat

    stale_forms = (
        "universal bulk algebra of the gravitational theory",
        "the derived center classifies bulk observables",
        "Three-dimensional gravity is the closed-sector "
        "($\\Sigma_n$-averaged) projection",
    )
    for stale in stale_forms:
        assert stale not in gravity


def test_concordance_derived_centre_bulk_summaries_are_scoped():
    root = Path(__file__).resolve().parents[2]
    files = {
        "concordance": root / "chapters/connections/concordance.tex",
        "hochschild": root / "chapters/connections/hochschild.tex",
        "extension": root / "chapters/connections/thqg_ht_bbl_extensions.tex",
        "preface_trimmed": root / "chapters/frame/preface_trimmed.tex",
    }
    text = {name: path.read_text() for name, path in files.items()}
    flat = {name: " ".join(value.split()) for name, value in text.items()}

    assert (
        "This is the algebraic closed-sector vertex of the boundary chart. "
        "It becomes a physical bulk observable algebra only after a chosen "
        "HT prefactorization realisation satisfies the bulk--Hochschild "
        "comparison hypotheses"
        in flat["concordance"]
    )
    assert (
        "the boundary-linear exact-sector comparison identifies physical "
        "bulk with the Morita-invariant derived centre of the boundary "
        "category"
        in flat["hochschild"]
    )
    assert (
        "the derived centre is the closed-colour comparison vertex"
        in flat["extension"]
    )
    assert (
        "A physical bulk reading requires the same HT prefactorization and "
        "exact-sector comparison data"
        in flat["extension"]
    )
    assert (
        "the derived centre is the algebraic closed-sector vertex; the "
        "physical bulk reading requires the HT comparison data"
        in flat["preface_trimmed"]
    )

    repaired_surface = "\n".join(text.values())
    stale_forms = (
        "serves as the\n universal bulk algebra",
        "serves as the universal bulk algebra",
        "boundary-linear bulk is the derived center of boundary",
        "the bulk is the derived center",
        "the bulk is the derived centre",
        "bulk is the derived center",
        "bulk is the derived centre",
    )
    for stale in stale_forms:
        assert stale not in repaired_surface


def test_active_bulk_shorthand_is_closed_sector_scoped():
    root = Path(__file__).resolve().parents[2]
    files = {
        "preface": root / "chapters/frame/preface.tex",
        "introduction": root / "chapters/theory/introduction.tex",
        "super_yangian": root / "chapters/theory/super_chiral_yangian.tex",
    }
    text = {name: path.read_text() for name, path in files.items()}
    flat = {name: " ".join(value.split()) for name, value in text.items()}

    assert (
        "The algebraic closed-sector vertex is "
        "\\(\\Zderch(\\cH_k)=\\bulkChirHoch{\\cH_k}\\)"
        in flat["preface"]
    )
    assert (
        "it is not the assertion that the boundary chart \\(\\cH_k\\) "
        "literally equals the bulk"
        in flat["preface"]
    )
    assert (
        "The level-$\\mathsf{Z}$ object is the algebraic derived chiral centre"
        in flat["introduction"]
    )
    assert (
        "its physical bulk reading requires the scoped HT open/closed "
        "comparison"
        in flat["introduction"]
    )
    assert (
        "The Virasoro algebraic closed-sector vertex is the chiral "
        "Hochschild object"
        in flat["introduction"]
    )
    assert "passes through the Brown--Henneaux comparison" in flat[
        "introduction"
    ]
    assert (
        "whose algebraic closed sector is the super derived chiral centre"
        in flat["super_yangian"]
    )
    assert (
        "the physical bulk reading is part of the same HT open/closed "
        "comparison package"
        in flat["super_yangian"]
    )

    repaired_surface = "\n".join(text.values())
    stale_forms = (
        "$\\Abulk = \\cH_k$: the bulk is the abelian chiral algebra",
        "The bulk is the level-$\\mathsf{Z}$ derived chiral centre",
        "The Virasoro bulk is the chiral Hochschild object",
        "whose bulk is the super derived chiral centre",
    )
    for stale in stale_forms:
        assert stale not in repaired_surface


def test_ht_physical_origin_requires_constructed_realisation():
    root = Path(__file__).resolve().parents[2]
    ht_origins = (
        root / "chapters/connections/ht_physical_origins.tex"
    ).read_text()
    flat = " ".join(ht_origins.split())

    assert (
        "The physical-origin lane of the E$_1$ core contains those chiral "
        "algebras for which a holomorphic-topological field theory in the "
        "sense of Costello--Li, or an explicit open/closed comparison "
        "datum, has been constructed"
        in flat
    )
    assert (
        "an abstract chiral algebra in the E$_1$ core does not thereby "
        "acquire a physical HT origin"
        in flat
    )
    assert (
        "the algebraic closed-sector vertex is the chiral derived center"
        in flat
    )
    assert (
        "Its physical bulk reading requires the chosen HT open/closed "
        "comparison"
        in flat
    )

    stale_forms = (
        "Every chiral algebra in the E$_1$ core arises as the boundary "
        "theory of a holomorphic-topological field theory",
        "the bulk algebra is a separate object, the chiral derived center",
    )
    for stale in stale_forms:
        assert stale not in ht_origins


def test_ht_path_integral_claims_are_perturbative_bv_scoped():
    root = Path(__file__).resolve().parents[2]
    ht_origins = (
        root / "chapters/connections/ht_physical_origins.tex"
    ).read_text()
    bv_ht = (root / "chapters/connections/bv_ht_physics.tex").read_text()
    flat = " ".join(ht_origins.split())
    flat_bv = " ".join(bv_ht.split())

    assert "Perturbative localization to holomorphic data" in ht_origins
    assert (
        "Fix a Costello--Li holomorphic twist, a BV gauge fixing, and a "
        "renormalisation scheme"
        in flat
    )
    assert (
        "it is not an ordinary convergent measure-theoretic path integral"
        in flat
    )
    assert "Its underived truncation is the moduli space" in flat
    assert (
        "Determined by the renormalised boundary factorization product of "
        "the gauge-fixed bulk theory"
        in flat
    )
    assert (
        "an unrenormalised analytic HCS path integral is not being invoked"
        in flat
    )
    assert (
        "This is computed from the gauge-fixed HCS BV Feynman rules"
        in flat
    )
    assert (
        "it is not a bare measure-theoretic path-integral identity"
        in flat
    )
    assert "\\subsubsection{Perturbative BV functional}" in ht_origins
    assert "\\subsection{Perturbative BV functional}" in bv_ht
    assert "\\begin{definition}[Gauge-fixed perturbative BV functional]" in ht_origins
    assert "\\begin{definition}[Gauge-fixed perturbative BV functional]" in bv_ht
    assert (
        "the perturbative BV functional is modeled by the bar-cobar "
        "pairing"
        in flat
    )
    assert "The gauge-fixed perturbative BV functional is formally denoted by" in flat
    assert "The gauge-fixed perturbative BV functional is formally denoted by" in flat_bv
    assert (
        "[D\\phi]$ denotes the induced formal BV Berezinian density "
        "after the renormalisation prescription has been fixed"
        in flat
    )
    assert (
        "It is not an ordinary measure on a space of fields"
        in flat
    )
    assert (
        "After the configuration-space density, gauge-fixing "
        "Lagrangian, and BV counterterm scheme are fixed, the "
        "perturbative BV functional is modeled by the bar--cobar pairing"
        in flat_bv
    )
    assert "\\emph{Step~2: Formal BV Berezinian density.}" in ht_origins
    assert "\\emph{Step~2: Formal BV Berezinian density.}" in bv_ht
    assert (
        "The induced BV Berezinian density on $\\mathcal{L}$ is formally "
        "written"
        in flat
    )
    assert (
        "Geometrically, this is represented by a configuration-space "
        "density"
        in flat
    )
    assert "Perturbative BV functional & Bar-cobar pairing" in ht_origins
    assert "Perturbative BV functional & Bar-cobar pairing" in bv_ht
    assert (
        "the algebraic configuration-space model for the perturbative BV "
        "integral"
        in flat
    )
    assert (
        "the algebraic configuration-space model for the perturbative BV "
        "integral"
        in flat_bv
    )
    assert "\\begin{remark}[Configuration-space density input]" in ht_origins
    assert "\\begin{remark}[Configuration-space density input]" in bv_ht
    assert (
        "The comparison between the formal perturbative BV functional and "
        "the bar-cobar pairing requires a configuration-space density "
        "model for BV gauge-fixing data"
        in flat
    )
    assert (
        "The comparison between the formal perturbative BV functional and "
        "the bar--cobar pairing requires a configuration-space density "
        "model for BV gauge-fixing data"
        in flat_bv
    )

    stale_forms = (
        "After the holomorphic twist \\cite{costello-renormalization}, "
        "the path integral localizes to:",
        "Since $\\{Q, V\\} \\geq 0$, the path integral localizes to",
        "Determined by bulk path integral with boundary insertions",
        "The OPE arises from short-distance singularities of the HCS "
        "path integral with boundary insertions",
        "This is derivable from the HCS path integral.",
        "\\subsubsection{BV path integral}",
        "\\subsection{BV path integral}",
        "\\begin{definition}[BV partition function]",
        "The BV path integral is realized by the bar-cobar pairing",
        "is exactly the BV path integral with gauge fixing",
        "The BV partition function is formally written as",
        "The BV partition function is:",
        "The BV path integral is represented by the bar--cobar pairing",
        "The identification of the BV path integral with the bar-cobar "
        "pairing",
        "In that model the BV partition function $Z_{\\mathrm{BV}}$ is "
        "the bar--cobar pairing",
        "\\emph{Step~2: BV Measure.}",
        "where $\\mathcal{L}$ is a Lagrangian submanifold (gauge fixing)",
        "and the integration uses the BV measure (Berezinian)",
        "The BV measure on $\\mathcal{L}$ is:",
        "Geometrically, this is the measure on configuration space",
        "Path integral & Bar-cobar pairing & Configuration integrals",
    )
    for stale in stale_forms:
        assert stale not in ht_origins
        assert stale not in bv_ht


def test_thqg_bv_extensions_use_graph_shadow_not_bv_partition_function():
    root = Path(__file__).resolve().parents[2]
    target = (
        root / "chapters/connections/thqg_bv_ht_extensions.tex"
    ).read_text()
    flat = " ".join(target.split())

    assert "\\section{Perturbative BV graph functionals}" in target
    assert "The construction below is narrower and algebraic" in flat
    assert "\\subsection{Graph-weight density on the bar complex}" in target
    assert "\\begin{definition}[Graph-weight BV density]" in target
    assert "graph-weight BV density on $\\barBch(\\cA)$" in flat
    assert (
        "\\begin{definition}[Genus-$g$ perturbative BV graph shadow]"
        in target
    )
    assert "Z_g^{\\mathrm{BV,gr}}(\\cA)" in target
    assert (
        "\\begin{proposition}[Perturbative graph shadow as "
        "bar-complex integral]"
        in target
    )
    assert "not an independent non-perturbative BV path integral" in flat
    assert "\\subsection{Heisenberg: the Gaussian determinant shadow}" in target
    assert "ordinary field-space measure on maps" in flat
    assert "\\begin{construction}[Graph density and Mumford isomorphism]" in target
    assert "\\subsection{Anomaly-cancelled scalar shadows}" in target
    assert (
        "Perturbative BV graph shadow $Z_g^{\\mathrm{BV,gr}}$ &"
        in target
    )
    assert "\\subsection{Virasoro: gravity scalar shadow}" in target
    assert "\\begin{construction}[Virasoro graph shadow from bar data]" in target
    assert "Z_g^{\\mathrm{BV,gr}}(\\mathrm{Vir}_c)" in target
    assert (
        "not a construction of the gravitational path integral over "
        "hyperbolic $3$-manifolds"
        in flat
    )
    assert (
        "perturbative shadow corrections in the bar complex"
        in flat
    )

    stale_forms = (
        "\\section{BV integration and partition functions}",
        "The conjecture that the BV path integral equals the "
        "bar-cobar pairing",
        "The key construction is the BV measure on the bar complex",
        "The partition function at genus $g$ is then a specific "
        "integral over",
        "\\subsection{BV measure on the bar complex}",
        "\\begin{definition}[Graph-weight BV measure]",
        "The \\emph{BV measure on $\\barBch(\\cA)$ at genus $g$}",
        "\\begin{definition}[BV partition function]",
        "The \\emph{BV partition function at genus $g$} is",
        "\\begin{proposition}[BV partition function as bar-complex "
        "integral]",
        "Then the BV partition function $Z_g(\\cA)$ equals the "
        "bar-cobar pairing",
        "The graph-weight BV measure",
        "\\subsection{Heisenberg: the Gaussian partition function}",
        "\\begin{theorem}[Heisenberg partition function]",
        "The graph-weight sum reduces to a Gaussian integral",
        "replacing the path integral over maps",
        "\\begin{construction}[BV measure and Mumford isomorphism]",
        "the BV partition function at genus $g$ is related to the "
        "Mumford isomorphism",
        "The BV measure $\\mu^{\\mathrm{BV}}_g$ is a section",
        "\\subsection{Anomaly-cancelled partition functions}",
        "The anomaly-cancelled partition function is:",
        "BV partition function $Z_g$ &",
        "The BV partition function at genus $g \\geq 2$",
        "\\subsection{Virasoro: gravity partition function}",
        "\\begin{construction}[Virasoro partition function from bar data]",
        "The Virasoro partition function at genus $g$ presents",
        "At genus $g$, the partition function is:",
        "Z_g(\\mathrm{Vir}_c)",
        "gives a precise version of their modular sum",
        "the combinatorial analogue of the sum over hyperbolic "
        "3-manifolds in the gravitational path integral",
        "``gravitational instantons'' of the Maloney--Witten approach",
    )
    for stale in stale_forms:
        assert stale not in target


def test_bv_brst_heisenberg_uses_determinant_line_scalar_not_bv_partition():
    root = Path(__file__).resolve().parents[2]
    vol1_root = root.parent / "chiral-bar-cobar"
    vol2 = (root / "chapters/connections/bv_brst.tex").read_text()
    vol1 = (vol1_root / "chapters/connections/bv_brst.tex").read_text()
    vol1_compute = (
        vol1_root / "compute/lib/heisenberg_bv_bar_proof.py"
    ).read_text()
    flat2 = " ".join(vol2.split())
    flat1 = " ".join(vol1.split())

    required_manuscript = (
        "The zero-mode-reduced gauge-fixed free-boson determinant-line "
        "scalar",
        "Z_g^{\\mathrm{det}}(\\cH_\\kappa;\\,\\Sigma_g)",
        "no ordinary measure on the field space is being asserted",
        "The Heisenberg determinant-line scalar factorizes as",
        "Z_g^{\\mathrm{det}}(\\cH_\\kappa;\\Sigma_g)",
        "whose determinant-line scalar is determined by a "
        "zeta-regularized determinant",
    )
    for required in required_manuscript:
        assert required in flat2
        assert required in flat1

    assert (
        "the moduli-space density, gauge-fixing data, and Costello "
        "renormalization framework"
        in flat2
    )
    assert "not an unconstructed field-space measure" in flat2
    assert (
        "field-theoretic comparison with the renormalized BV-BRST "
        "complex of the higher-genus WZW model"
        in flat2
    )
    assert "it is not supplied by an unconstructed field-space measure" in flat1
    assert "The BV side is the determinant-line scalar" in vol1_compute
    assert "Z_g^{det}(H_k) = (det'_zeta dbar_{Sigma_g})^{-k}" in vol1_compute
    assert "The determinant-line scalar uses zeta-regularized determinants" in vol1_compute

    stale_forms = (
        "The BV partition function of $\\kappa$ free bosons",
        "Gaussian functional integral",
        "the path integral measure on $\\overline{\\mathcal{M}}_g$",
        "Z_g^{\\mathrm{BV}}(\\cH_\\kappa;\\,\\Sigma_g)",
        "The Heisenberg partition function factorizes as",
        "whose partition function is determined by a zeta-regularized",
        "the WZW path integral at higher genus",
    )
    for stale in stale_forms:
        assert stale not in vol2
        assert stale not in vol1
        assert stale not in vol1_compute


def test_factorization_swiss_cheese_k3_delta5_is_scalar_anomaly_characteristic():
    root = Path(__file__).resolve().parents[2]
    target = (
        root / "chapters/theory/factorization_swiss_cheese.tex"
    ).read_text()
    flat = " ".join(target.split())

    assert (
        "The perturbative Swiss-cheese anomaly line at $K3 \\times T^2$ "
        "has a scalar Deligne characteristic"
        in flat
    )
    assert (
        "After the Pfaffian orientation, Hall--Borcherds comparison, "
        "and K3 elliptic-genus normalisation are fixed"
        in flat
    )
    assert "the primitive Borcherds trivialising section" in flat
    assert "the BKM denominator $\\Delta_5$" in flat
    assert (
        "\\kappa_{\\mathrm{BKM}}(\\mathfrak g_{\\Delta_5})"
        in target
    )
    assert "c_{\\phi_{0,1}^{K3}}(0)/2" in target

    stale_forms = (
        "The Swiss-cheese BV partition function at $K3 \\times T^2$ "
        "receives a",
        "paramodular anomaly contribution in "
        "$H^2(\\overline{\\mathcal{A}_2}, \\mathcal{O}^*)$; "
        "the unique trivialiser",
        "the unique trivialiser in Deligne cohomology",
    )
    for stale in stale_forms:
        assert stale not in target


def test_ordered_drinfeld_boundary_pairing_is_gated_scalar_shadow():
    root = Path(__file__).resolve().parents[2]
    target = (
        root / "chapters/connections/ordered_associative_chiral_kd_frontier.tex"
    ).read_text()
    flat = " ".join(target.split())

    assert "\\textbf{Boundary pairing comparison.}" in target
    assert (
        "after the boundary polarization, localization density, and "
        "hemisphere comparison map are fixed"
        in flat
    )
    assert (
        "the Dimofte hemisphere scalar amplitude maps to the cyclic "
        "form on $B^{\\mathrm{ord}}(\\cA)$"
        in flat
    )
    assert (
        "This is not an identification of a physical partition function "
        "with the raw ordered bar complex"
        in flat
    )
    assert "the scalar boundary-pairing shadow required" in flat
    assert "When all four data are supplied" in flat

    stale_forms = (
        "\\textbf{Physics equals algebra at the boundary pairing.}",
        "the hemisphere partition function equals the cyclic pairing "
        "on the ordered bar complex",
        "The physical pairing (Dimofte hemisphere) equals the algebraic "
        "pairing",
        "This is the final consistency check: the reconstructor "
        "constructed algebraically agrees with the reconstructor "
        "constructed physically",
        "Taken together, (a) existence, (b) uniqueness up to natural "
        "transformation, (c) correctness, and (d) physics matches "
        "algebra, establish",
    )
    for stale in stale_forms:
        assert stale not in target


def test_modular_pva_heisenberg_hemisphere_uses_scalar_amplitude():
    root = Path(__file__).resolve().parents[2]
    target = (
        root / "chapters/connections/modular_pva_quantization_frontier.tex"
    ).read_text()
    flat = " ".join(target.split())

    assert "\\index{hemisphere scalar amplitude|textbf}" in target
    assert "\\label{subsec:hemisphere-scalar-amplitude}" in target
    assert "\\index{hemisphere pairing!cyclic comparison datum}" in target
    assert (
        "\\index{abelian Chern--Simons!hemisphere scalar amplitude}"
        in target
    )
    assert (
        "localised hemisphere scalar amplitude of a $3$d "
        "holomorphic-topological theory"
        in flat
    )
    assert (
        "after a boundary polarisation, an Omega-background "
        "localisation density, and a one-loop determinant "
        "regularisation have been chosen"
        in flat
    )
    assert (
        "After a boundary polarization, Omega-background localization "
        "data, and one-loop determinant regularization are fixed"
        in flat
    )
    assert (
        "Z_{D^3}^{\\mathrm{scal}}\\bigl[\\mathrm{U}(1)_k\\bigr]"
        in target
    )
    assert (
        "formal Gaussian determinant-line prediction after the "
        "Omega-background localisation and determinant regularisation "
        "are fixed"
        in flat
    )

    stale_forms = (
        "\\index{hemisphere partition function!equals cyclic pairing}",
        "\\index{abelian Chern--Simons!hemisphere partition function}",
        "\\label{subsec:hemisphere-partition-function}",
        "the hemisphere partition function is the perturbative path "
        "integral around the Gaussian saddle",
        "the hemisphere partition function is",
        "expression above is obtained as a direct perturbative Gaussian "
        "integral",
        "computation of the $D^3$ hemisphere partition function for "
        "abelian CS",
        "Z_{D^3}\\bigl[\\mathrm{U}(1)_k\\bigr]",
    )
    for stale in stale_forms:
        assert stale not in target


def test_ht_agt_and_class_s_localization_are_scoped():
    root = Path(__file__).resolve().parents[2]
    ht_origins = (
        root / "chapters/connections/ht_physical_origins.tex"
    ).read_text()
    bv_ht = (root / "chapters/connections/bv_ht_physics.tex").read_text()
    flat = " ".join(ht_origins.split())
    flat_bv = " ".join(bv_ht.split())

    assert (
        "Equivariant localization computes the Nekrasov partition "
        "function by fixed-point data on the instanton/Hitchin moduli "
        "problem"
        in flat
    )
    assert "Under the AGT comparison hypotheses" in flat
    assert (
        "it does not turn the unregularised 4d path integral into a "
        "literal integral over \\(\\mathcal{M}_{\\mathrm{Hit}}\\)"
        in flat
    )
    assert (
        "followed by BV gauge fixing and dimensional reduction to a "
        "Riemann surface"
        in flat
    )
    assert (
        "produces the perturbative holomorphic BF theory whose derived "
        "classical moduli problem is the \\(\\bar\\partial\\)-connection "
        "stack"
        in flat
    )
    assert "Step~2: Perturbative BV localization" in ht_origins
    assert (
        "the perturbative factorization algebra is computed in the formal "
        "neighbourhood of the \\(Q\\)-fixed locus"
        in flat
    )
    assert (
        "This OPE is computed from the renormalised holomorphic "
        "Chern--Simons/BF BV Feynman rules"
        in flat
    )
    assert (
        "produces the perturbative holomorphic BF theory whose derived "
        "classical moduli problem is the \\(\\bar\\partial\\)-connection "
        "stack"
        in flat_bv
    )
    assert "Step~2: Perturbative BV localization" in bv_ht
    assert (
        "This OPE is computed from the renormalised holomorphic "
        "Chern--Simons/BF BV Feynman rules"
        in flat_bv
    )

    stale_forms = (
        "The $\\Omega$-background localizes the path integral to "
        "$\\mathcal{M}_{\\mathrm{Hit}}$",
        "localizes to $\\bar{\\partial}$-connections on a Riemann surface",
        "By localization, path integral reduces to:",
        "This OPE is computed via the holomorphic Chern--Simons path integral",
    )
    for stale in stale_forms:
        assert stale not in ht_origins
        assert stale not in bv_ht


def test_chiral_ce_one_loop_uses_renormalised_feynman_residue():
    root = Path(__file__).resolve().parents[2]
    chiral_ce = (
        root / "chapters/connections/chiral_ce_factalg_gen_rel.tex"
    ).read_text()
    flat = " ".join(chiral_ce.split())

    assert "Fix the Costello renormalisation scheme and the Bergman propagator" in flat
    assert "The coefficient is the renormalised BV Feynman-rule residue" in flat
    assert (
        "configuration-space integral on $\\P^1\\times\\P^1$, not by an "
        "unqualified measure-theoretic path integral over fields"
        in flat
    )
    assert (
        "$\\zeta$-regularisation~\\cite[\\S10]{CostelloRenormalization}"
        in chiral_ce
    )

    stale_forms = (
        "Computing the coefficient via the path integral",
        "via the path integral (Costello renormalisation",
    )
    for stale in stale_forms:
        assert stale not in chiral_ce


def test_feynman_connection_one_loop_uses_zeta_regularised_determinant():
    root = Path(__file__).resolve().parents[2]
    feynman = (root / "chapters/connections/feynman_connection.tex").read_text()
    flat = " ".join(feynman.split())

    assert (
        "The zero-mode-normalised, zeta-regularised Gaussian determinant "
        "of the free boson on $E_\\tau$ is"
        in flat
    )
    assert "Z_{E_\\tau}=(\\operatorname{Im}\\tau)^{-1/2}|\\eta(\\tau)|^{-2}" in flat
    assert "analytic determinant-line amplitude of the free theory" in flat
    assert (
        "or construct it as a field-space measure integral"
        in flat
    )

    stale_forms = (
        "The Gaussian path integral on $E_\\tau$ gives",
        "Gaussian path integral on $E_\\tau$",
    )
    for stale in stale_forms:
        assert stale not in feynman


def test_spectral_braiding_cfg_uses_perturbative_bv_feynman_expansion():
    root = Path(__file__).resolve().parents[2]
    spectral = (
        root / "chapters/connections/spectral-braiding-core.tex"
    ).read_text()
    flat = " ".join(spectral.split())

    assert "Comparison with Costello--Francis--Gwilliam" in spectral
    assert (
        "the filtered structure records the loop expansion of the "
        "renormalised perturbative BV Feynman rules"
        in flat
    )
    assert "It is not a non-perturbative BV path integral" in flat

    stale_forms = (
        "the filtered structure reflects the loop expansion of the\nBV path integral",
        "the filtered structure reflects the loop expansion of the BV path integral",
    )
    for stale in stale_forms:
        assert stale not in spectral
        assert stale not in flat


def test_affine_generic_center_not_promoted_to_cs_path_integral():
    root = Path(__file__).resolve().parents[2]
    proved = (
        root / "chapters/examples/examples-complete-proved.tex"
    ).read_text()
    flat = " ".join(proved.split())

    assert (
        "At generic non-critical level, $Z(\\widehat{\\fg}_k)=\\C$ is a "
        "degree-zero algebraic statement"
        in flat
    )
    assert (
        "the chiral Hochschild closed-sector complex has no central local "
        "class beyond the unit"
        in flat
    )
    assert (
        "Under a chosen perturbative Chern--Simons/BV boundary comparison "
        "this unit is the vacuum summand of local closed observables"
        in flat
    )
    assert (
        "It is not a theorem about the full non-perturbative "
        "Chern--Simons functional integral"
        in flat
    )
    assert (
        "global flat-connection sectors, framing dependence, and "
        "finite-level modular data require separate $3$d input"
        in flat
    )
    assert "support shadow depth" in flat
    assert "r^\\Theta_{\\max}=3" in proved
    assert (
        "This is a support-packet statement, not an assertion that "
        "obstruction-class vanishing determines the packet"
        in flat
    )
    assert (
        "\\(r^\\Theta_{\\max}\\) fixes the support envelope of the "
        "canonical bar-intrinsic MC packet"
        in flat
    )
    assert (
        "It is not an obstruction-depth invariant and it does not by "
        "itself determine the homological length of \\(B(A)\\)"
        in flat
    )
    assert (
        "The bar weight profile is fixed only after the OPE-residue data"
        in flat
    )

    stale_forms = (
        "unique bulk vacuum of nonabelian Chern--Simons",
        "the CS path integral $\\int \\mathcal{D}A\\, "
        "e^{ikS_{\\mathrm{CS}}(A)}$ admits no adjustable parameter "
        "beyond the level",
        "$\\rmax",
        "\\rmax",
        "r_{\\max}",
        "r_max",
        "rmax",
        "truncation of the shadow obstruction tower",
        "homological length of $B(A)$ on the Koszul locus",
        "The shadow obstruction tower terminates at degree~$3$",
        "the MC element terminates at finite degree",
        "infinite shadow depth",
    )
    for stale in stale_forms:
        assert stale not in proved
    assert not re.search(r"\\rmax|r_\{\\max\}|r_max|rmax", proved)


def test_bv_brst_hcs_finiteness_is_perturbative_regulator_statement():
    root = Path(__file__).resolve().parents[2]
    bv_brst = (root / "chapters/connections/bv_brst.tex").read_text()
    flat = " ".join(bv_brst.split())

    assert (
        "asserts that the renormalised perturbative 6d hCS/BV amplitudes "
        "on $\\bR\\times K3\\times E$ are finite at every genus in the "
        "Fulton--MacPherson compactified configuration-space model"
        in flat
    )
    assert (
        "This is a statement about gauge-fixed counterterms and "
        "factorisation algebra amplitudes after the regulator has been "
        "installed"
        in flat
    )
    assert (
        "not an ordinary non-perturbative analytic hCS measure integral"
        in flat
    )

    stale_forms = (
        "asserts that the 6D hCS path integral on $\\bR\\times K3\\times E$ "
        "is finite at every genus with Fulton--MacPherson as geometric "
        "regulator",
        "6D hCS path integral on $\\bR\\times K3\\times E$ is finite",
    )
    for stale in stale_forms:
        assert stale not in bv_brst


def test_gravity_cs_torsion_is_semiclassical_not_full_path_integral():
    root = Path(__file__).resolve().parents[2]
    gravity = (
        root / "chapters/connections/thqg_3d_gravity_movements_vi_x.tex"
    ).read_text()
    flat = " ".join(gravity.split())

    assert (
        "formulation gives the same semiclassical one-loop datum after a "
        "flat connection and boundary polarization have been fixed"
        in flat
    )
    assert (
        "the quadratic BV determinant is the Ray--Singer/Reidemeister "
        "torsion"
        in flat
    )
    assert (
        "For the thermal AdS$_3$ vacuum this torsion is the "
        "eta/vacuum-character determinant factor after the chosen "
        "boundary-mode normalization"
        in flat
    )
    assert (
        "It is not the full Chern--Simons path integral over all flat "
        "connections or the complete modular sum"
        in flat
    )

    stale_forms = (
        "the CS path integral is computed by the Reidemeister torsion",
        "$Z_{\\mathrm{CS}} = |\\eta(\\tau)|^{-2}$ at leading order",
    )
    for stale in stale_forms:
        assert stale not in gravity


def test_gravity_effective_scalar_shadow_series_not_partition_function():
    root = Path(__file__).resolve().parents[2]
    gravity = (
        root / "chapters/connections/thqg_3d_gravity_movements_vi_x.tex"
    ).read_text()
    flat = " ".join(gravity.split())

    assert (
        "Spectral decomposition of the genus-$g$ effective scalar "
        "shadow series"
        in flat
    )
    assert "\\index{spectral decomposition!genus scalar shadow series}" in gravity
    assert "The genus-$g$ effective scalar shadow series decomposes as" in flat
    assert (
        "the effective scalar shadow series is computed by the curved "
        "scalar sector alone"
        in flat
    )
    assert (
        "The scalar generating series is "
        "$Z_{\\mathrm{grav}}^{\\mathrm{scal}} = "
        "\\exp(\\cF_{\\mathrm{grav}}^{\\mathrm{scal}})$"
        in flat
    )
    assert (
        "the effective scalar shadow series is the exponential series "
        "generated by the normalized $\\hat A$-class weighted by $-13/2$"
        in flat
    )
    assert (
        "Self-duality means invariance under the reflection "
        "$\\cS\\colon c\\mapsto 26-c$, not inversion of the scalar "
        "series"
        in flat
    )

    stale_forms = (
        "Spectral decomposition of the genus-$g$ effective scalar "
        "partition function",
        "\\index{spectral decomposition!genus partition function}",
        "The genus-$g$ effective scalar partition function decomposes as",
        "the effective scalar partition function is computed by the "
        "curved scalar sector alone",
        "The exponentiation: $Z_{\\mathrm{grav}}^{\\mathrm{scal}}",
        "the effective scalar partition function is the exponential of "
        "the $\\hat A$-genus",
        "The self-duality means\n$Z_{13} = Z_{13}^{-1}$",
        "The self-duality means $Z_{13} = Z_{13}^{-1}$",
    )
    for stale in stale_forms:
        assert stale not in gravity


def test_symplectic_polarization_holography_uses_comparison_datum():
    root = Path(__file__).resolve().parents[2]
    symplectic = (
        root / "chapters/connections/thqg_symplectic_polarization.tex"
    ).read_text()
    flat = " ".join(symplectic.split())

    assert (
        "\\subsection{Holographic comparison datum: genus-$1$ capacity "
        "and BTZ}"
        in flat
    )
    assert (
        "\\index{Lagrangian capacity!holographic comparison|textbf}"
        in symplectic
    )
    assert (
        "physical AdS$_3$/CFT$_2$ comparison datum is supplied"
        in flat
    )
    assert (
        "completed boundary-CFT genus-$g$ partition function together "
        "with a chosen bulk saddle expansion, or with a gravitational "
        "functional integral after the bulk measure and contour "
        "prescription have been specified"
        in flat
    )
    assert (
        "Complementarity supplies the algebraic bulk--boundary "
        "polarization entering such a comparison"
        in flat
    )
    assert (
        "it does not construct the bulk measure, prove modular "
        "invariance, select saddles, or imply a physical entropy formula"
        in flat
    )
    assert (
        "\\subsubsection{Genus-$1$ complementarity and the BTZ "
        "comparison package}"
        in flat
    )
    assert (
        "\\subsubsection{BTZ comparison package and the complementarity "
        "potential}"
        in flat
    )
    assert "\\index{BTZ black hole!comparison package}" in symplectic

    stale_forms = (
        "the physical AdS$_3$/CFT$_2$ dictionary is supplied",
        "the genus-$g$ partition function of a boundary CFT is "
        "compared with a gravitational path integral over bulk "
        "geometries",
        "Complementarity provides the algebraic bulk--boundary "
        "polarization; it does not by itself prove saddle dominance, "
        "modular invariance, or a physical entropy formula",
        "Genus-$1$ complementarity and the BTZ partition function",
        "BTZ partition function and the complementarity potential",
        "\\index{bulk--boundary entanglement|textbf}",
        "\\index{BTZ black hole!partition function}",
    )
    for stale in stale_forms:
        assert stale not in symplectic


def test_symplectic_polarization_entropy_is_capacity_bound_not_entropy():
    root = Path(__file__).resolve().parents[2]
    symplectic = (
        root / "chapters/connections/thqg_symplectic_polarization.tex"
    ).read_text()
    flat = " ".join(symplectic.split())

    assert (
        "the construction reduces to a one-line Lagrangian-capacity "
        "computation"
        in flat
    )
    assert (
        "A physical BTZ or entanglement interpretation requires the "
        "AdS$_3$/CFT$_2$ comparison datum, Hilbert completion, modular "
        "invariance, vacuum dominance, and saddle analysis"
        in flat
    )
    assert "single scalar degree-$2$ channel and no nonlinear shadow jets" in flat
    assert "\\rank H_\\cA^{(1)} = \\dim Q_1(\\cA)=1" in symplectic
    assert (
        "the genus-$1$ cubic and quartic shadow projections vanish"
        in flat
    )
    assert (
        "bulk--boundary pairing & comparison channel"
        in flat
    )
    assert (
        "\\subsubsection{Lagrangian capacity bound from the "
        "complementarity potential}"
        in flat
    )
    assert (
        "\\begin{proposition}[Genus-$1$ Lagrangian capacity bound; "
        "\\ClaimStatusConditional]"
        in flat
    )
    assert (
        "physical comparison datum realizes \\(Q_1(\\cA)\\) and "
        "\\(Q_1(\\cA^!)\\) as trace-class Hilbert factors"
        in flat
    )
    assert (
        "Complementarity itself proves only the algebraic capacity and "
        "pairing"
        in flat
    )
    assert (
        "it does not construct a density matrix, a partial trace, or a "
        "Hilbert tensor product"
        in flat
    )
    assert (
        "The Lagrangian polarization gives a direct-sum decomposition"
        in flat
    )
    assert "not a tensor product" in flat
    assert (
        "\\begin{remark}[Ryu--Takayanagi comparison datum]"
        in flat
    )
    assert "the rank of the Hessian is a Lagrangian capacity, not an area" in flat
    assert (
        "requires a boundary subregion, a bulk metric, a homology "
        "constraint, and a variational minimal-surface problem"
        in flat
    )
    assert (
        "it is not itself \\(\\mathrm{Area}(\\gamma_A)/(4G_N)\\)"
        in flat
    )
    assert "BTZ comparison capacity" in flat

    stale_forms = (
        "bulk--boundary entanglement computation",
        "bulk--boundary pairing & entanglement",
        "The complementarity potential at genus $1$ is quadratic",
        "linear in the coordinate on $Q_1(\\cA)$",
        "cotangent normal form gives a linear function",
        "Entanglement entropy from the complementarity potential",
        "Entanglement entropy from the Hessian",
        "the entanglement entropy of the bulk--boundary system at "
        "genus~$1$ is determined by the Hessian",
        "-\\Tr(H_\\cA^{(1)} \\log H_\\cA^{(1)})",
        "The entanglement entropy measures the correlation between",
        "the reduced density matrix on\n$Q_1(\\cA)$",
        "the entanglement entropy is nonzero and related to the",
        "the ``area'' is the rank of the Hessian",
        "The complementarity potential $S_\\cA$ provides the\n"
        "``gravitational action'' whose critical locus determines",
        "BTZ bulk--boundary entanglement",
    )
    for stale in stale_forms:
        assert stale not in symplectic


def test_symplectic_polarization_potential_is_comparison_data_not_action():
    root = Path(__file__).resolve().parents[2]
    symplectic = (
        root / "chapters/connections/thqg_symplectic_polarization.tex"
    ).read_text()
    gravity = (root / "chapters/connections/3d_gravity.tex").read_text()
    fm_calculus = (
        root / "chapters/connections/thqg_fm_calculus_extensions.tex"
    ).read_text()
    flat = " ".join(symplectic.split())
    gravity_flat = " ".join(gravity.split())
    fm_flat = " ".join(fm_calculus.split())

    assert (
        "\\begin{remark}[The complementarity potential under a "
        "gravitational comparison datum]"
        in flat
    )
    assert "\\index{complementarity potential!comparison datum}" in symplectic
    assert (
        "The complementarity potential $S_\\cA$ is the formal function "
        "whose graph is the dual Lagrangian"
        in flat
    )
    assert (
        "A physical interpretation of $S_\\cA$ as an action requires "
        "a comparison datum"
        in flat
    )
    assert (
        "a bulk field space, BV/gauge-fixing data, boundary "
        "polarization, a normalization of observables, and a "
        "renormalization scheme"
        in flat
    )
    assert (
        "Only after such data have been fixed can the algebraic shadow "
        "channels be read as metric or Chern--Simons vertices"
        in flat
    )
    assert "\\textbf{Algebraic role}" in symplectic
    assert "\\textbf{After comparison}" in symplectic
    assert "Newton-coupling channel" in flat
    assert "possible cubic vertex channel" in flat
    assert "possible quartic/contact channel" in flat
    assert "higher vertex or composite channels" in flat
    assert "G class: only the quadratic shadow channel survives" in flat
    assert "L class: a cubic shadow channel is present" in flat
    assert "C class: a quartic/contact shadow channel is present" in flat
    assert "M class: the shadow obstruction tower is unbounded" in flat
    assert (
        "It is not a theorem that three-dimensional gravity is "
        "perturbatively non-renormalizable"
        in flat
    )
    assert (
        "does not produce gravitational counterterms unless a physical "
        "renormalization problem and comparison map have first been "
        "specified"
        in flat
    )
    assert "The four terms have algebraic collision-channel roles" in gravity_flat
    assert (
        "After a Brown--Henneaux/metric comparison datum is fixed, they "
        "may be read as"
        in gravity_flat
    )
    assert "central cubic vertex comparison channel" in gravity_flat
    assert (
        "to \\(1/G_N\\) only after the Brown--Henneaux normalization has "
        "been installed"
        in gravity_flat
    )
    assert "the central cubic comparison channel" in fm_flat
    assert (
        "After a metric comparison datum has been fixed, this channel "
        "can be compared with a cubic gravitational vertex"
        in fm_flat
    )

    stale_forms = (
        "The complementarity potential as gravitational\naction",
        "plays the role of the gravitational action",
        "Its Taylor jets are gravitational coupling constants",
        "\\textbf{Gravitational}",
        "three-graviton coupling",
        "four-graviton coupling",
        "higher-point graviton couplings",
        "determines the perturbative\ncomplexity of the gravitational dual",
        "gravity is free (Gaussian, no graviton\n interactions)",
        "gravity has cubic coupling only (tree-level\n GR)",
        "gravity has quartic coupling (one-loop\n GR)",
        "gravity has all-order couplings (full quantum\n gravity)",
        "algebraic counterpart of the perturbative non-renormalizability "
        "of\nthree-dimensional quantum gravity",
        "gravitational counterterm at $r$-point order",
    )
    for stale in stale_forms:
        assert stale not in symplectic

    propagation_stale_forms = (
        "The four terms have direct gravitational content",
        "is the three-graviton coupling ($\\propto 1/G$)",
        "the \\emph{three-graviton coupling}",
    )
    for stale in propagation_stale_forms:
        assert stale not in gravity
        assert stale not in fm_calculus


def test_modular_pva_gv_one_loop_is_gauge_fixed_seed():
    root = Path(__file__).resolve().parents[2]
    modular = (
        root / "chapters/connections/modular_pva_quantization_core.tex"
    ).read_text()
    flat = " ".join(modular.split())

    assert (
        "matches the gauge-fixed semiclassical one-loop determinant of "
        "the Chern--Simons BV complex at the trivial flat connection"
        in flat
    )
    assert "after the same boundary-mode normalization" in flat
    assert "It is not the full matrix-model functional integral" in flat
    assert (
        "The full Gopakumar--Vafa result requires summing over all flat "
        "connections"
        in flat
    )

    stale_forms = (
        "matches the one-loop determinant of the CS path integral around "
        "the trivial flat connection",
        "CS path integral around the trivial flat connection",
    )
    for stale in stale_forms:
        assert stale not in modular


def test_heisenberg_recursion_uses_gaussian_determinant_not_path_integral():
    root = Path(__file__).resolve().parents[2]
    bootstrap = (
        root / "chapters/connections/thqg_modular_bootstrap.tex"
    ).read_text()
    flat = " ".join(bootstrap.split())

    assert (
        "The genus-$g$ Gaussian determinant-line amplitude of the "
        "gauge-fixed free-boson factorization algebra is the $k$-th "
        "determinant-line power"
        in flat
    )
    assert (
        "This is the perturbative determinant model for the Heisenberg "
        "sector, not an unqualified measure-theoretic path integral"
        in flat
    )
    assert (
        "Products of lower-genus terms are disconnected "
        "determinant-line amplitudes"
        in flat
    )
    assert "Thus the connected class is $k\\omega_g$" in flat
    assert (
        "The zero-mode Gaussian determinant-line scalar carried by a "
        "period matrix $\\Omega_g$ is"
        in flat
    )
    assert "\\chi_{g,\\mathrm{det}}^{(d)}(\\Omega_g)" in bootstrap
    assert (
        "This analytic determinant-line scalar is a representative of "
        "the Hodge determinant line; it is not the exponential of the "
        "Faber--Pandharipande scalar free energy"
        in flat
    )
    assert (
        "the integrated connected Hodge free energy and the zero-mode "
        "determinant-line scalar, not a physical partition function"
        in flat
    )

    stale_forms = (
        "The genus-$g$ Heisenberg path integral is the $k$-th "
        "determinant-line power",
        "genus-$g$ Heisenberg path integral",
        "contributions to the partition function and disappear after "
        "the",
        "The genus-$g$ partition function is:",
        "The partition function is the exponential of the free energy,",
        "Z_g^{(d)}",
        "eq:thqg-VII-rank-d-Zg",
        "gives $Z_g = (\\det(\\mathrm{Im}\\,\\Omega))^{-dk/2}$",
    )
    for stale in stale_forms:
        assert stale not in bootstrap


def test_modular_bootstrap_outputs_scalar_shadow_not_gravity_partition():
    root = Path(__file__).resolve().parents[2]
    bootstrap = (
        root / "chapters/connections/thqg_modular_bootstrap.tex"
    ).read_text()
    flat = " ".join(bootstrap.split())

    assert (
        "complete for the fixed bar--modular shadow: it determines "
        "the connected scalar shadow class at each genus"
        in flat
    )
    assert "hence the genus-$g$ scalar shadow class is determined" in flat
    assert "One-loop scalar shadow from genus-$1$ MC" in flat
    assert "\\index{one-loop gravitational scalar shadow}" in bootstrap
    assert "Z_{1}^{\\mathrm{scal}}(\\cA)" in bootstrap
    assert (
        "The genus-$1$ scalar shadow contribution is "
        "$Z_1^{\\mathrm{scal}} = \\exp(F_1)$"
        in flat
    )
    assert (
        "\\mathcal{A}_g^{\\mathrm{scal,conn}}(\\cA)"
        in bootstrap
    )
    assert (
        "\\kappaChHodge(\\cA)\\cdot \\lambda_g^{\\mathrm{FP}}"
        in bootstrap
    )
    assert "The disconnected scalar shadow generating series is:" in flat
    assert "Z_{\\mathrm{sh}}^{\\mathrm{scal}}(\\cA;\\hbar)" in bootstrap
    assert (
        "No claim about a physical gravitational partition function "
        "is made here"
        in flat
    )
    assert (
        "The output is the completed scalar shadow generating series"
        in flat
    )
    assert (
        "conditional on a separate comparison between \\(\\Theta_\\cA\\), "
        "a bulk saddle expansion, and the completed physical partition "
        "function"
        in flat
    )
    assert "No finite support cutoff for mixed-type algebras" in bootstrap
    assert "bar-intrinsic support packet can force degeneration" in flat
    assert (
        "Actual nonvanishing of a differential \\(d_r\\) remains a "
        "separate cohomological statement"
        in flat
    )
    assert "Spectral sequence differentials and support projections" in flat
    assert (
        "non-scalar sources from the bar-intrinsic support packet"
        in flat
    )
    assert (
        "This support cutoff is stronger than eventual vanishing of "
        "obstruction classes"
        in flat
    )
    assert (
        "not equivalent to a statement that every corresponding "
        "spectral-sequence differential is nonzero"
        in flat
    )
    assert "It does not say \\(o_{N+1}(\\cA)\\neq0\\)" in flat
    assert (
        "\\(r^\\Theta_{\\max}(\\cA)\\) (support shadow depth)"
        in bootstrap
    )
    assert "Support-degree bounds on gravitational amplitudes" in bootstrap
    assert "not a graph-degree theorem" in flat
    assert "Support depth bounds modular-bootstrap complexity" in bootstrap
    assert "no depth argument can force finite-page degeneration" in flat

    stale_forms = (
        "it determines the partition\nfunction at each genus",
        "it determines the partition function at each genus",
        "hence the genus-$g$ partition function is determined",
        "partition function is uniquely determined",
        "One-loop gravity from genus-$1$ MC",
        "\\index{one-loop gravitational correction}",
        "gravitational correction to the partition function",
        "The genus-$1$ contribution to the partition function is",
        "\\mathcal{A}_g(\\cA)\n\\;=\\;\n\\kappaChHodge(\\cA)^g",
        "\\kappaChHodge(\\cA)^g \\cdot F_g",
        "The partition function at all genera is:",
        "Z_g(\\cA)",
        "The output is the complete gravitational partition function",
        "determined by the boundary data without free parameters.\n\n"
        "Within the fixed HT construction",
        "shadow depth $r_{\\max}",
        "Shadow depth $r_{\\max}$",
        "r_{\\max}",
        "r_max",
        "rmax",
        "the shadow obstruction tower terminates at weight~$2$",
        "the shadow obstruction tower is infinite",
        "shadow obstruction tower $\\Theta^{\\leq r}$",
        "the $d_r$ differential is nontrivial for each~$r$",
        "if $\\cA$ has shadow depth $r_{\\max}$",
        "$o_{r+1}(\\cA) = 0$ for $r \\geq r_{\\max}$",
        "the genus spectral sequence collapses at $E_2$",
        "the genus spectral sequence collapses at $E_3$",
        "the genus spectral sequence does not degenerate at any finite page",
        "The shadow depth $r_{\\max}(\\cA)$ determines the collapse",
        "The shadow obstruction tower $\\Theta^{\\leq r_{\\max}} = "
        "\\Theta_\\cA$ stabilizes",
    )
    for stale in stale_forms:
        assert stale not in bootstrap
    assert not re.search(r"r_\{\\max\}|r_max|rmax", bootstrap)


def test_thqg_finiteness_shadow_free_energy_not_path_integral():
    root = Path(__file__).resolve().parents[2]
    finiteness = (
        root / "chapters/connections/thqg_perturbative_finiteness.tex"
    ).read_text()
    flat = " ".join(finiteness.split())

    assert (
        "The shadow free energy $F_g(\\cA)$ is the genus-$g$ scalar "
        "vacuum-energy shadow of the twisted holographic theory"
        in flat
    )
    assert (
        "connected, no-insertion coefficient of the algebraic shadow "
        "determinant expansion"
        in flat
    )
    assert (
        "Its physical reading as a genus-$g$ cosmological constant "
        "requires a chosen bulk comparison datum"
        in flat
    )
    assert (
        "not a contribution to an independently constructed path integral"
        in flat
    )
    assert (
        "complete, explicit, genus-by-genus calculation of the "
        "scalar-shadow partition series"
        in flat
    )
    assert "\\index{BTZ black hole!scalar shadow comparison}" in finiteness
    assert (
        "not the full gravitational partition function over bulk "
        "geometries"
        in flat
    )
    assert (
        "$F_g(\\cA)$ (shadow free energy) & genus-$g$ scalar vacuum "
        "amplitude"
        in flat
    )
    assert "metric one-loop determinant" in flat

    stale_forms = (
        "is the genus-$g$ vacuum energy of the twisted holographic theory: "
        "it is the contribution to the path integral",
        "contribution to the path integral from genus-$g$ worldsheets",
        "calculation of the gravitational partition function",
        "$F_g(\\cA)$ (shadow free energy) & genus-$g$ path integral",
        "gravity path integral (Giombi--Maloney--Yin",
        "\\index{BTZ black hole!partition function}",
    )
    for stale in stale_forms:
        assert stale not in flat


def test_thqg_degree_summed_series_not_full_gravity_partition_function():
    root = Path(__file__).resolve().parents[2]
    finiteness = (
        root / "chapters/connections/thqg_perturbative_finiteness.tex"
    ).read_text()
    fredholm = (
        root / "chapters/connections/thqg_fredholm_partition_functions.tex"
    ).read_text()
    flat = " ".join(finiteness.split())
    fredholm_flat = " ".join(fredholm.split())

    assert "Scalar shadow genus series" in finiteness
    assert "Z_{\\mathrm{sh}}^{\\mathrm{scal}}(T;\\,\\hbar)" in finiteness
    assert (
        "The definition is the degree-$0$ trace of the shadow "
        "Maurer--Cartan element"
        in flat
    )
    assert (
        "not a full partition function over bulk metrics or all shadow "
        "degrees"
        in flat
    )
    assert (
        "Completed degree-summed shadow partition series and Gaussian "
        "scalar convergence"
        in flat
    )
    assert "the completed degree-summed shadow partition series satisfies" in flat
    assert "Z_g^{\\mathrm{sh,deg}}(T)" in finiteness
    assert (
        "Completed degree-summed shadow two-parameter convergence" in flat
    )
    assert "\\index{shadow partition series!two-parameter convergence}" in finiteness
    assert "fixed-genus degree-summed shadow expression" in flat
    assert "Z_g^{\\mathrm{sh,deg}}(T;\\,t)" in finiteness
    assert "degree-summed shadow two-parameter series" in flat
    assert (
        "Z_{\\mathrm{sh}}^{\\mathrm{deg}}(T;\\,\\hbar,t)"
        in finiteness
    )
    assert (
        "bidisk for the completed degree-summed shadow genus expansion"
        in flat
    )
    assert "finite support shadow depth" in flat
    assert "r^\\Theta_{\\max}(\\cA)<\\infty" in finiteness
    assert "r^\\Theta_{\\max}(\\cA)=\\infty" in finiteness
    assert "bar-intrinsic support packet" in flat
    assert "This is stronger than eventual vanishing of the obstruction classes" in flat
    assert "does not kill the corresponding \\(H^1\\)-lift coordinate" in flat
    assert "support-depth class" in flat
    assert "\\(r^\\Theta_{\\max}(\\cA)\\) (support shadow depth)" in finiteness
    assert (
        "the degree-weighted shadow trace into an analytic function "
        "of~$t$"
        in flat
    )
    assert (
        "Z_{\\mathrm{sh}}^{\\mathrm{deg}}(T;\\,\\hbar) \\;=\\; "
        "Z_{\\mathrm{sh}}^{\\mathrm{scal}}(T;\\,\\hbar)"
        in flat
    )
    assert "degree-summed shadow partition series" in fredholm_flat
    assert "Z_g^{\\mathrm{sh,deg}} = Z_g^{\\mathrm{G}} + \\delta Z_g" in fredholm
    assert (
        "A comparison with a physical gravitational partition function "
        "requires additional bulk data"
        in fredholm_flat
    )

    stale_forms = (
        "\\begin{definition}[Gravitational partition function]",
        "the \\emph{gravitational partition function} is",
        "The \\emph{scalar gravitational partition function} is the "
        "restriction to the degree-$0$ scalar channel",
        "The full gravitational partition function includes the "
        "higher-degree contributions",
        "Completed full partition function and Gaussian scalar convergence",
        "full gravitational partition function (including all degrees) "
        "satisfies",
        "Z_g^{\\mathrm{full}}(T)",
        "Z_{\\mathrm{grav}}(T;\\,\\hbar) \\;=\\; "
        "Z_{\\mathrm{grav}}^{\\mathrm{scal}}(T;\\,\\hbar)",
        "the full partition function equals the scalar partition function",
        "Completed degree series and Gaussian two-parameter convergence",
        "\\index{gravitational partition function!2D convergence}",
        "Z_{\\mathrm{grav}}(T;\\,\\hbar,t)",
        "bidisk for the full degree-summed genus expansion",
        "the scalar trace into an analytic function of~$t$",
        "finite shadow depth",
        "infinite shadow depth",
        "The shadow depth $r_{\\max}(\\cA)$ bounds the degree of contributing graphs",
        "obstruction class $o_{r+1}(\\cA) = 0 and the shadow obstruction tower stabilizes",
        "Part~(iii) is the definition of shadow depth",
        "$r_{\\max}$",
        "r_{\\max}",
        "r_max",
    )
    for stale in stale_forms:
        assert stale not in finiteness
    assert not re.search(r"r_\{\\max\}|r_max|rmax|d_\\infty|f_\\infty", finiteness)

    fredholm_stale_forms = (
        "The full partition function",
        "Z_g^{\\mathrm{full}} = Z_g^{\\mathrm{G}} + \\delta Z_g",
        "is the analytic realization of the complete MC element",
    )
    for stale in fredholm_stale_forms:
        assert stale not in fredholm


def test_thqg_scalar_closed_forms_are_shadow_data_before_comparison():
    root = Path(__file__).resolve().parents[2]
    finiteness = (
        root / "chapters/connections/thqg_perturbative_finiteness.tex"
    ).read_text()
    claims = (root / "metadata/claims.jsonl").read_text()
    graph = (root / "metadata/dependency_graph.dot").read_text()
    flat = " ".join(finiteness.split())
    metadata = claims + "\n" + graph

    assert "Convergence radius of the scalar shadow closed form" in flat
    assert "\\index{convergence radius!scalar shadow closed form}" in finiteness
    assert "Z_{\\mathrm{sh}}^{\\mathrm{scal}}(\\cA;\\,\\hbar)" in finiteness
    assert (
        "Completed scalar shadow series on the proved scalar lane" in flat
    )
    assert (
        "Brown--Henneaux normalization of the scalar shadow closed form"
        in flat
    )
    assert (
        "Substituting $c = 3\\ell/(2G_N)$ and "
        "$\\hbar = 2G_N/(3\\ell)$ in the scalar shadow closed form gives"
        in flat
    )
    assert "Non-perturbative comparison datum" in flat
    assert (
        "only after a bulk completion, contour, measure, "
        "unitarity/positivity structure, and saddle sector have been "
        "supplied"
        in flat
    )
    assert (
        "Positivity datum: positivity for real $\\hbar > 0$ belongs "
        "to the completed physical partition function"
        in flat
    )
    assert "positivity structure, and saddle sector" in flat
    assert "shadow entropy" in flat
    assert "choose a branch of the logarithm" in flat
    assert "Z_{\\mathrm{sh}}^{\\mathcal{W}_N,\\mathrm{scal}}" in finiteness
    assert "R_{\\mathrm{sh}}^{\\mathrm{scal}}" in finiteness
    assert (
        "quartic/contact comparison channel"
        in flat
    )
    assert (
        "Convergence radius of the scalar shadow closed form" in metadata
    )
    assert (
        "Leading Brown--Henneaux expansion of the reduced scalar "
        "shadow series"
        in metadata
    )
    assert (
        "Brown--Henneaux normalization of the scalar shadow closed form"
        in metadata
    )
    assert "Non-perturbative comparison datum" in metadata

    stale_forms = (
        "Convergence radii of the scalar partition function",
        "The scalar gravitational partition function "
        "$Z^{\\mathrm{scal}}_{\\mathrm{grav}}",
        "\\index{convergence radius!scalar partition function}",
        "Leading gravitational coupling expansion of the reduced "
        "scalar series",
        "\\index{gravitational coupling!large-$c$ expansion}",
        "Substituting $c = 3\\ell/(2G)$ and "
        "$\\hbar = 2G/(3\\ell)$, the gravitational partition "
        "function becomes",
        "Newton's constant expansion",
        "\\index{Newton's constant!expansion}",
        "\\index{non-perturbative!gravitational partition function}",
        "The meromorphic function $Z^{\\mathrm{scal}}_{\\mathrm{grav}}",
        "admits a non-perturbative completion satisfying",
        "Positivity: for real $\\hbar > 0$, the full "
        "non-perturbative partition function $Z > 0$",
        "Z_{\\mathrm{grav,scal}}",
        "R_{\\mathrm{grav}}",
        "two-loop graviton self-energy",
    )
    for stale in stale_forms:
        assert stale not in finiteness
        assert stale not in metadata


def test_thqg_faber_pandharipande_scalar_lane_uses_kappa_ch_hodge():
    root = Path(__file__).resolve().parents[2]
    finiteness = (
        root / "chapters/connections/thqg_perturbative_finiteness.tex"
    ).read_text()
    start = finiteness.index(
        "\\begin{proposition}[Shadow scalar coefficient in closed form"
    )
    end = finiteness.index(
        "% SUBSECTION 3: CONVERGENCE OF THE SCALAR SHADOW SERIES",
        start,
    )
    scalar_lane = finiteness[start:end]
    flat = " ".join(scalar_lane.split())

    assert (
        "F_g^{\\mathrm{sc}}(\\cA)=\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}"
        in flat
    )
    assert (
        "\\sum_{g=1}^{\\infty} F_g^{\\mathrm{sc}}(\\cA)\\,x^{2g} &\\;=\\; "
        "\\kappaChHodge(\\cA) \\sum_{g=1}^{\\infty}"
        in flat
    )
    assert (
        "\\kappaChHodge(\\cA)\\bigl(x/(2\\sin(x/2)) - 1\\bigr)"
        in scalar_lane
    )
    assert (
        "\\frac{2|\\kappaChHodge(\\cA)|}{(2\\pi)^{2g}}"
        in scalar_lane
    )
    assert (
        "8(-1)^n \\kappaChHodge(\\cA)n^2\\pi^2"
        in scalar_lane
    )
    assert (
        "F_1^{\\mathrm{sc}}(\\cA) &= \\frac{\\kappaChHodge(\\cA)}{24}"
        in scalar_lane
    )
    assert (
        "\\caption{Scalar shadow free energies $F_g^{\\mathrm{sc}}(\\cA)="
        "\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}$"
        in scalar_lane
    )

    stale_forms = (
        "F_g = \\kappa \\cdot \\lambda_g^{\\mathrm{FP}}",
        "value of $\\kappa$",
        "\\sum_{g=1}^{\\infty} F_g\\,x^{2g} &\\;=\\; \\kappa \\sum",
        "&\\;=\\; \\kappa \\left( \\frac{x/2}{\\sin(x/2)} - 1 \\right)",
        "\\kappaChHodge(x/(2\\sin(x/2)) - 1)",
        "$\\kappa \\cdot (\\frac{\\sqrt{\\hbar}/2}{\\sin(\\sqrt{\\hbar}/2)} - 1)$",
        "2|\\kappa|",
        "4|\\kappa|",
        "F_1(\\cA) &= \\frac{\\kappa}{24}",
        "\\frac{7\\kappa}{5760}",
        "F_g = \\kappa \\cdot \\frac{2^{2g-1}-1}",
        "\\caption{Shadow free energies $F_g(\\cA) = "
        "\\kappa \\cdot \\lambda_g^{\\mathrm{FP}}",
        "converging to $\\kappa \\cdot",
    )
    for stale in stale_forms:
        assert stale not in scalar_lane

    whole_file_stale_forms = (
        "\\kappa \\cdot \\lambda_g^{\\mathrm{FP}}",
        "\\kappa/24",
        "|\\kappa|",
    )
    for stale in whole_file_stale_forms:
        assert stale not in finiteness


def test_live_vol2_fp_lane_echoes_use_kappa_ch_hodge():
    root = Path(__file__).resolve().parents[2]
    paths = []
    paths.extend(sorted((root / "chapters").rglob("*.tex")))
    appendices = root / "appendices"
    if appendices.exists():
        paths.extend(sorted(appendices.rglob("*.tex")))
    paths.extend(sorted((root / "compute" / "lib").rglob("*.py")))

    live_surface = "\n".join(
        f"% {path.relative_to(root)}\n{path.read_text()}" for path in paths
    )

    stale_exact_forms = (
        "\\kappa \\cdot \\lambda_g^{\\mathrm{FP}}",
        "\\kappa\\cdot\\lambda_g^{\\mathrm{FP}}",
        "\\kappa/24",
        "\\kappa \\cdot c_g",
    )
    for stale in stale_exact_forms:
        assert stale not in live_surface

    stale_patterns = (
        r"F_g(?:\([^)]*\))?\s*=\s*\\kappa(?![A-Za-z_])"
        r"(?:\s*\\cdot|\s*\\,)?\s*\\lambda_g\^\{\\mathrm\{FP\}\}",
        r"F_1(?:\([^)]*\))?\s*=\s*\\kappa(?![A-Za-z_])\s*/24",
        r"\\kappa(?![A-Za-z_])\s*\\lambda_g\^\{\\mathrm\{FP\}\}",
        r"\\kappa(?![A-Za-z_])\s*\\cdot\s*c_g",
        r"\\kappa(?![A-Za-z_])\s*\\cdot\s*B_2",
    )
    for pattern in stale_patterns:
        assert re.search(pattern, live_surface) is None


def test_gravity_effective_fp_lane_names_hodge_scalar_source():
    root = Path(__file__).resolve().parents[2]
    gravity = root / "chapters/connections/3d_gravity.tex"
    movements = root / "chapters/connections/thqg_3d_gravity_movements_vi_x.tex"

    gravity_text = gravity.read_text()
    movements_text = movements.read_text()
    combined = gravity_text + "\n" + movements_text

    required_gravity_forms = (
        "\\kappa_{\\mathrm{eff}}\n"
        "\\;:=\\;\n"
        "\\kappaChHodge(\\mathrm{Vir}_c \\otimes bc_2)",
        "F_g^{\\mathrm{eff}}(\\mathrm{Vir}_c \\otimes bc_2)\n"
        "= \\kappaChHodge(\\mathrm{Vir}_c \\otimes bc_2)\n"
        "  \\lambda_g^{\\mathrm{FP}}",
        "F_g^{\\mathrm{dS},\\mathrm{sc}}(\\mathrm{Vir}_{c_{\\mathrm{dS}}})\n"
        " = \\kappaChHodge(\\mathrm{Vir}_{c_{\\mathrm{dS}}})\n"
        " \\lambda_g^{\\mathrm{FP}}",
    )
    for required in required_gravity_forms:
        assert required in gravity_text

    required_movement_forms = (
        "F_g^{\\mathrm{eff}}(\\mathrm{Vir}_{26}\\otimes bc_2)\n"
        "= \\kappaChHodge(\\mathrm{Vir}_{26}\\otimes bc_2)"
        "\\lambda_g^{\\mathrm{FP}}",
        "F_g^{\\mathrm{eff}}(\\mathrm{Vir}_{26}\\otimes bc_2)\n"
        "= \\kappaChHodge(\\mathrm{Vir}_{26}\\otimes bc_2)"
        "\\lambda_g^{\\mathrm{FP}}\n"
        "= \\kappa_{\\mathrm{eff}}\\lambda_g^{\\mathrm{FP}}$ vanishes",
        "\\sum_{g \\geq 1} F_g^{\\mathrm{eff}}\\hbar^{2g}\n"
        "= \\kappaChHodge(\\mathrm{Vir}_c\\otimes bc_2)"
        "(\\hat A(i\\hbar)-1)$",
    )
    for required in required_movement_forms:
        assert required in movements_text

    stale_forms = (
        "F_g = \\kappa_{\\mathrm{eff}}\\lambda_g^{\\mathrm{FP}}",
        "F_g = \\kappa_{\\mathrm{eff}} \\cdot \\lambda_g^{\\mathrm{FP}}",
        "F_g^{\\mathrm{eff}} = \\kappa_{\\mathrm{eff}}"
        "\\lambda_g^{\\mathrm{FP}}",
        "F_g^{\\mathrm{eff}} = \\kappa_{\\mathrm{eff}} \\cdot",
        "\\sum_{g \\geq 1} F_g \\hbar^{2g} = "
        "\\kappa_{\\mathrm{eff}}(\\hat A(i\\hbar)-1)",
        "F_g = \\kappa_{\\mathrm{dS}}\\cdot\\lambda_g^{\\mathrm{FP}}",
    )
    for stale in stale_forms:
        assert stale not in combined


def test_thqg_handlebody_example_is_comparison_not_partition_function():
    root = Path(__file__).resolve().parents[2]
    finiteness = (
        root / "chapters/connections/thqg_perturbative_finiteness.tex"
    ).read_text()
    flat = " ".join(finiteness.split())

    assert "Scalar-shadow comparison with selected 3-manifold saddles" in flat
    assert "the thermal AdS$_3$ saddle in the usual comparison" in flat
    assert "Z_{\\mathrm{char}}(\\cA;\\tau)" in finiteness
    assert (
        "agrees with the holomorphic one-loop graviton determinant of "
        "\\cite{GMY08} after the thermal-AdS boundary-mode normalization "
        "has been fixed"
        in flat
    )
    assert (
        "The scalar shadow value $F_1=c/48="
        "\\kappaChHodge(\\mathrm{Vir}_c)/24$ is the Hodge trace"
        in flat
    )
    assert (
        "is a candidate bulk filling after a handlebody saddle "
        "prescription has been chosen"
        in flat
    )
    assert (
        "A physical handlebody partition function requires a "
        "conformal-block pairing, a measure or sum over intermediate "
        "labels, and a metric normalization of the handlebody action"
        in flat
    )
    assert "Z^{\\mathrm{block}}_{H_2}(\\Omega)" in finiteness
    assert "S_{H_2}^{\\mathrm{met}}(\\Omega)" in finiteness
    assert (
        "Vacuum dominance and the metric action "
        "\\(S_{H_2}^{\\mathrm{met}}\\) are part of the physical "
        "comparison datum"
        in flat
    )
    assert "F_2^{\\mathrm{scal}}(\\mathrm{Vir}_c)" in finiteness
    assert (
        "is the connected genus-$2$ Hodge-trace shadow coefficient on "
        "the uniform-weight scalar lane"
        in flat
    )
    assert (
        "It is not a computed two-loop handlebody determinant and not a "
        "physical handlebody partition function"
        in flat
    )

    stale_forms = (
        "The partition function on specific 3-manifolds",
        "The partition function is $Z_{H_2}",
        "the dominant contribution is the vacuum block",
        "Z_{H_2} \\approx |\\det(\\operatorname{Im}\\Omega)|^{c/2}",
        "S_{\\mathrm{grav}}(\\Omega)",
        "where $S_{\\mathrm{grav}}$ is the gravitational action of the "
        "handlebody filling",
        "The shadow free energy $F_2 = 7c/11520$ is the two-loop "
        "correction from the genus-$2$ moduli integral",
    )
    for stale in stale_forms:
        assert stale not in finiteness


def test_recognition_scope_blocks_global_reading():
    profile = recognition_scope_profile()

    assert profile["global_recognition_theorem"] is False
    assert profile["recognized_surface"] == (
        "product-formal local-shadow rectangle data"
    )
    assert "arbitrary global Ran-space factorization data" in profile["not_recognized"]
    assert "chosen prefactorization realization" in profile["not_recognized"][1]
    assert profile["output"] == "local C_*(W(SC^{ch,top})) algebra shadow"


def test_recognition_scope_manuscript_source_guard():
    root = Path(__file__).resolve().parents[2]
    examples = (root / "chapters/examples/examples-complete.tex").read_text()
    examples_core = (root / "chapters/examples/examples-complete-core.tex").read_text()
    bar_cobar = (root / "chapters/connections/bar-cobar-review.tex").read_text()
    wave17 = (root / "compute/tests/test_climax_theorems_wave17_iv.py").read_text()
    flat_examples = " ".join(examples.split())
    flat_examples_core = " ".join(examples_core.split())
    flat_bar_cobar = " ".join(bar_cobar.split())

    assert "product-formal local-shadow recognition theorem" in flat_examples
    assert "product-formal local-shadow recognition theorem" in flat_examples_core
    assert "product-formal local-shadow recognition conditions" in flat_bar_cobar
    assert "product-formal local-shadow recognition theorem" in flat_bar_cobar
    assert "Programme product-formal local-shadow recognition theorem" in wave17
    assert "Programme recognition theorem for log SC" not in wave17


def test_deletion_ledger_manuscript_source_guard():
    root = Path(__file__).resolve().parents[2]
    source = (root / "chapters/connections/programme_climax_platonic.tex").read_text()
    flat = " ".join(source.split())

    assert "\\label{tab:universal-holography-deletion-ledger}" in source
    assert "Deletion ledger for the master theorem" in source
    assert "\\(\\SCchtop{}\\) recognition is global" in source
    assert "class \\(\\mathbf M\\) raw direct-sum chain theorem" in source
    assert "\\(E_3\\)-PBW proves concentration" in source
    assert "DS--Hochschild all nilpotents" in source
    assert "bar scalar trace \\(=\\) Maloney--Witten sum" in source
    assert "K3/Class-\\(\\mathcal S\\) comparison automatic" in source
    assert "product-formal local-shadow" in flat
    assert "weight-completed/pro or finite-propagation ambient" in flat
    assert "Finite Hall--Borcherds gate" in flat


def test_corrected_theorem_form_manuscript_source_guard():
    root = Path(__file__).resolve().parents[2]
    source = (root / "chapters/connections/programme_climax_platonic.tex").read_text()
    flat = " ".join(source.split())

    assert "\\label{rem:universal-holography-sequent-form}" in source
    assert "\\Xi(A)\\ \\vdash_{\\mathcal A(A)}" in source
    assert "\\Phihol(A,\\Xi_A)=T_{A,\\Xi}" in source
    assert "\\mathrm{Obs}^{\\partial}(T_{A,\\Xi})\\simeq A" in source
    assert "\\mathrm{Obs}^{\\mathrm{bulk}}(T_{A,\\Xi})\\simeq" in source
    assert "\\sigma(A)\\in" in source
    assert "\\ClaimStatusProvedElsewhere" in source
    assert "A\\in\\{\\mathbf G,\\mathbf L,\\mathbf C\\}" in source
    assert "\\widehat{\\mathrm{Ch}}_{\\mathrm{wt},\\rho}" in source
    assert "\\mathrm{pro}\\text{-}\\mathrm{Ch}(\\mathrm{Vect})" in source
    assert "\\varprojlim_N\\mathcal A_{\\le N}^{\\mathrm{wt}}" in source
    assert "raw direct-sum statement is not a theorem" in flat
    assert "For the non-\\(W_{\\infty}\\) class-\\(\\mathbf M\\) families" in flat
    assert "replace the uniform Banach algebra" in flat
    assert "arbitrary logarithmic \\(\\SCchtop\\)-algebras" in flat
    assert "\\(k=-h^{\\vee}\\), the affine object is outside" in flat
    assert "bounded-weight pro-window tower is the object" in flat
    assert "Hall--Borcherds operator gates" in flat


def test_holographic_reconstruction_is_xi_realized_comparison():
    root = Path(__file__).resolve().parents[2]
    reconstruction = (
        root / "chapters/connections/thqg_holographic_reconstruction.tex"
    ).read_text()
    part_vi = (
        root / "chapters/connections/part_vi_platonic_introduction.tex"
    ).read_text()
    flat_reconstruction = " ".join(reconstruction.split())
    flat_part_vi = " ".join(part_vi.split())

    assert "Universal Holography, admissible shadow projection with realization datum" in reconstruction
    assert "\\mathrm{ChirAlg}^{\\omega,\\mathrm{adm},\\Xi}_X" in reconstruction
    assert "\\mathrm{HT\\text{-}QFT}^{\\Xi}_{X\\times\\R}" in reconstruction
    assert "\\Xi_\\cA=(\\mathcal F_{\\cA,\\Xi},\\eta^\\partial_{\\cA}," in reconstruction
    assert "\\mathcal{T}_{\\cA,\\Xi}:=\\mathcal F_{\\cA,\\Xi}" in reconstruction
    assert "\\chi_{\\mathrm{HT},\\cA}\\colon" in reconstruction
    assert "not a construction of a physical bulk from \\(\\cA\\) alone" in flat_reconstruction
    assert "Class~$\\mathsf{M}$ cohomological and completed bulk--centre comparison" in reconstruction
    assert "a triangle of comparison maps in the declared ambient" in flat_reconstruction
    assert "\\cT_{\\cA,\\Xi}" in part_vi
    assert "\\eta^\\partial_\\cA" in part_vi
    assert "\\chi_{\\mathrm{HT},\\cA}\\colon" in part_vi
    assert "the compared bulk is $\\Einfty$-topological" in part_vi
    assert "\\cT_{\\mathrm{Vir}_c,\\Xi^{\\mathrm{BH}}}" in part_vi

    stale_forms = (
        "Then the admissible datum determines a canonical relative",
        "there is a canonical relative $3$d holomorphic-topological gauge theory",
        "the bulk is the chiral derived centre",
        "bulk as derived chiral centre",
        "The universal bulk identification",
        "the bulk is $Z^{\\mathrm{der}}_{\\mathrm{ch}}(\\cW^k(\\fg))$",
        "the bulk is $\\Einfty$-topological",
        "whose bulk is $\\Zderch(\\cA)",
        "$\\Phi_{\\mathrm{hol}}^{-1}$ on the derived chiral centre image",
    )
    repaired_surface = "\n".join((reconstruction, part_vi))
    for stale in stale_forms:
        assert stale not in repaired_surface


def test_shadow_tower_reconstruction_requires_lift_coordinates():
    root = Path(__file__).resolve().parents[2]
    source = (
        root / "chapters/connections/thqg_holographic_reconstruction.tex"
    ).read_text()
    celestial = (
        root / "chapters/connections/thqg_celestial_holography_extensions.tex"
    ).read_text()
    flat = " ".join(source.split())
    flat_celestial = " ".join(celestial.split())

    assert "Filtered shadow-tower reconstruction of the bar-intrinsic MC gauge class" in source
    assert "Obstruction classes alone decide extendability" in flat
    assert "they do not select a point of the \\(H^1\\)-torsor of lifts" in flat
    assert "\\mathfrak S(\\cA)=" in source
    assert "\\{(o_n(\\cA),\\lambda_n)\\}_{n\\geq 5}" in source
    assert "\\lambda_n\\in H^1(J^n(\\cA),d_2)" in source
    assert "bar-intrinsic gauge slice" in flat
    assert "Obstruction depth and support depth of the shadow tower" in source
    assert "d_{\\mathrm{obs}}(\\cA)" in source
    assert "r^\\Theta_{\\max}(\\cA)" in source
    assert "d_{\\mathrm{obs}}(\\cA)\\le r^\\Theta_{\\max}(\\cA)" in source
    assert "finite reconstruction is governed by \\(r^\\Theta_{\\max}\\)" in flat
    assert "The lifts above that degree may still vary by the \\(H^1(J^n(\\cA),d_2)\\)-torsors" in flat
    assert "Support shadow depth" in source
    assert "The obstruction depth \\(d_{\\mathrm{obs}}\\) records only eventual \\(H^2\\)-vanishing" in flat
    assert "The $k$-invariants $\\{k_r(\\cA)\\}_{r \\geq 2}$ determine the" in source
    assert "obstruction-depth profile \\(d_{\\mathrm{obs}}(\\cA)\\)" in flat
    assert "Mittag--Leffler for the bar-intrinsic shadow branch" in source
    assert "No surjectivity statement is made for arbitrary components" in flat
    assert "\\varprojlim^1_r E^\\Theta_\\cA(r) = 0" in source
    assert "\\MC(\\widehat{\\gAmod})_{\\Theta}" in source
    assert "a triangle of comparison maps" in flat
    assert "\\mathcal T_{\\mathcal A,\\Xi}" in source
    assert "\\chi_{\\mathrm{HT},\\mathcal A}" in source
    assert "comparison equivalence of $E_3$-topological" in flat
    assert "Finite support shadow depth controls the $q=0$ reconstruction" in celestial
    assert "support depth \\(r_\\Theta:=r^\\Theta_{\\max}(\\cA)\\)" in flat_celestial
    assert "all \\(H^1\\)-lift coordinates above degree \\(r_\\Theta\\) vanish" in flat_celestial
    assert "Finite-jet characterization of support shadow depth" in celestial
    assert "the obstruction classes \\(o_n(\\cA)\\) and the lift coordinates" in flat_celestial
    assert "\\(r^\\Theta_{\\max}(\\mathcal{H}_k) = 2\\)" in celestial
    assert "\\(r^\\Theta_{\\max}(\\widehat{\\mathfrak{g}}_k) = 3\\)" in celestial
    assert "\\(r^\\Theta_{\\max}(V_{\\beta\\gamma}) = 4\\)" in celestial
    assert "\\(r^\\Theta_{\\max}(\\mathrm{Vir}_c) = \\infty\\)" in celestial
    assert "scalar modular shadow series" in celestial
    assert "the dimension \\(\\dim(\\mathfrak{g})\\)" in flat_celestial
    assert "\\kappaChHodge(\\widehat{\\mathfrak{g}}_k) = \\dim(\\mathfrak{g})(k + h^\\vee)/(2h^\\vee)" in celestial
    assert "completed shadow reconstruction packet, including the lift coordinates" in flat_celestial
    assert "finite contact packet determines the scalar modular shadow series" in flat_celestial
    assert "$r^\\Theta_{\\max}$ & $2$ & $3$ & $4$ & $\\infty$" in source

    stale_forms = (
        "The full datum $\\Theta_\\cA$ is determined by the shadow-tower data:",
        "The reconstruction is:",
        "Universal Holography, strict $E_3$ projection]\\label{thm:holographic-recon-e3}",
        "\\mathrm{Obs}^{\\mathrm{bulk}}(\\mathcal T_{\\mathcal A}) \\simeq",
        "Obstruction classes alone determine",
        "Infinite data } \\{S_r\\}_{r \\geq 2} \\text{ required.}",
        "r_{\\max}(\\cA) < \\infty\n\\quad\\Longleftrightarrow\\quad\n\\Theta_\\cA^{\\leq r_{\\max}} = \\Theta_\\cA",
        "If $r_{\\max} < \\infty$, then\n$\\Theta_\\cA^{\\leq r} = \\Theta_\\cA^{\\leq r_{\\max}}$",
        "the component $\\Theta_{r+1} = 0 at every level above",
        "The $k$-invariants $\\{k_r(\\cA)\\}_{r \\geq 2}$ determine\n$\\Theta_\\cA$ up to gauge equivalence.",
        "$r_{\\max}$ equals the $L_\\infty$ obstruction depth",
        "\\varprojlim^1_r E_\\cA(r) = 0",
        "\\MC(\\gAmod) \\to \\varprojlim_r \\MC(\\gAmod / F^{r+1})",
        "Surjectivity for general",
    )
    for stale in stale_forms:
        assert stale not in source

    stale_celestial = (
        "Finite shadow depth controls the $q=0$ reconstruction",
        "the shadow obstruction tower stabilizes at degree $r_{\\max}$, then $\\Theta_\\cA^{\\leq r_{\\max}} = \\Theta_\\cA$",
        "\\item $r_{\\max}(\\cA) < \\infty$;",
        "\\item $J^{r_{\\max}}(\\cA) = J^\\infty(\\cA)$;",
        "the obstruction classes $o_{r+1}(\\cA) = 0$ for all $r \\geq r_{\\max}$",
        "r_{\\max}(\\mathcal{H}_k)",
        "r_{\\max}(\\widehat{\\mathfrak{g}}_k)",
        "r_{\\max}(V_{\\beta\\gamma})",
        "r_{\\max}(\\mathrm{Vir}_c)",
        "scalar modular partition function",
        "full scalar modular partition function",
        "rank $\\dim(\\mathfrak{g})$",
        "Four numbers ($\\kappa",
    )
    for stale in stale_celestial:
        assert stale not in celestial


def test_support_depth_is_not_obstruction_depth_in_gravity_and_soft_hierarchy():
    root = Path(__file__).resolve().parents[2]
    gravity = (
        root / "chapters/connections/thqg_gravitational_complexity.tex"
    ).read_text()
    soft = (
        root / "chapters/connections/thqg_soft_graviton_theorems.tex"
    ).read_text()
    flat_gravity = " ".join(gravity.split())
    flat_soft = " ".join(soft.split())

    assert "\\rho_{\\mathrm{grav}}(\\cA)" in gravity
    assert "r^\\Theta_{\\max}(\\cA)" in gravity
    assert "d_{\\mathrm{obs}}(\\cA):=" in gravity
    assert "d_{\\mathrm{obs}}(\\cA)\\le r^\\Theta_{\\max}(\\cA)" in flat_gravity
    assert "Obstruction depth only records eventual \\(H^2\\)-vanishing" in flat_gravity
    assert "d^\\Theta_\\infty(\\cA)" in gravity
    assert "f^\\Theta_\\infty(\\cA)" in gravity
    assert (
        "The \\(H^2\\)-obstruction depth \\(d_{\\mathrm{obs}}\\) is a "
        "separate invariant"
        in flat_gravity
    )
    assert (
        "If \\(r^\\Theta_{\\max}(\\cA)<\\infty\\), the bar-intrinsic MC packet"
        in flat_gravity
    )
    assert "Standard-family support-depth classification" in gravity
    assert "Support depth as a gravitational invariant" in gravity
    assert (
        "canonical bar-intrinsic branches of the standard landscape"
        in flat_gravity
    )
    assert (
        "not a classification of arbitrary components of the full filtered "
        "MC moduli problem"
        in flat_gravity
    )
    assert "obstruction depth and support depth may differ" in flat_gravity
    assert "canonical standard-family gauge slice" in flat_gravity
    assert "no monotonicity statement is made without controlling the \\(H^1\\)-lift coordinates" in flat_gravity
    assert "No classification of arbitrary components of the filtered MC moduli problem is asserted here" in flat_gravity
    assert "Support depth as stratification depth of the Steinberg" in gravity

    assert "support soft depth" in flat_soft
    assert "r_\\Theta(\\cA):=r^\\Theta_{\\max}(\\cA)" in soft
    assert "d_{\\mathrm{obs}}(\\cA):=" in soft
    assert (
        "records only eventual \\(H^2\\)-vanishing and may be strictly smaller"
        in flat_soft
    )
    assert "Let \\(r_\\Theta:=r^\\Theta_{\\max}(\\cA)\\)" in flat_soft
    assert (
        "The system of shadow Ward identities at degrees \\(2\\le r\\le "
        "r_\\Theta\\)"
        in flat_soft
    )
    assert (
        "support shadow depth equals \\(A_\\infty\\)-support depth equals "
        "\\(L_\\infty\\)-support depth"
        in flat_soft
    )
    assert "not by obstruction vanishing alone" in flat_soft
    assert "\\mathfrak{C}(\\cA)=\\Theta_{\\cA,3}\\neq0" in soft
    assert "\\Theta_{\\cA,4}=\\mathfrak Q(\\cA)" in soft
    assert "\\bigoplus_{2\\le r\\le r^\\Theta_{\\max}(\\cA)}" in soft
    assert "formal infinite OPE expansion" in flat_soft
    assert "Support shadow depth controls the truncation of the celestial OPE" in flat_soft

    retired_forms = (
        "r_{\\max}(\\cA) := \\sup",
        "The shadow obstruction tower stabilises at degree $r_{\\max}(\\cA)$",
        "r_{\\max}(\\cA) = d_\\infty(\\cA) = f_\\infty(\\cA)",
        "Shadow termination degree equals $A_\\infty$-depth equals",
        "The obstruction tower forces four mutually exclusive",
        "For an algebra of shadow depth $r_{\\max}(\\cA)$",
        "the system of shadow Ward identities\nat degrees $r = 2, 3, \\dotsc, r_{\\max}$",
        "The classification into exactly four classes is forced by the\nobstruction tower",
        "At critical values, $r_{\\max}$ can only decrease",
        "Let $\\cA$ have shadow depth",
        "The shadow depth r_max(A) controls",
    )
    repaired_surface = "\n".join((gravity, soft))
    for stale in retired_forms:
        assert stale not in repaired_surface
    assert not re.search(r"r_\{\\max\}|r_max|rmax|d_\\infty|f_\\infty", repaired_surface)
