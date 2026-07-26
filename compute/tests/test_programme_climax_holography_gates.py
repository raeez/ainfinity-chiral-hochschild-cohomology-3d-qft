"""Regression guards for the Universal Holography climax theorem.

The live theorem is a sequent from explicit HT realization data.  It must
not advertise a bare canonical 3d HT theory for the whole standard landscape.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "chapters/connections/programme_climax_platonic.tex"
INTRO_SOURCE = ROOT / "chapters/theory/introduction.tex"
PREFACE_SOURCE = ROOT / "chapters/frame/preface.tex"
LOG_WP_SOURCE = ROOT / "chapters/theory/logarithmic_wp_tempered_analysis_platonic.tex"


def _source() -> str:
    return SOURCE.read_text()


def _flat() -> str:
    return " ".join(_source().split())


def test_master_functor_is_xi_relative_with_bulk_hochschild_gate():
    source = _source()
    flat = _flat()

    assert "\\ChirAlgclimax^{\\omega, \\mathrm{BL,adm},\\Xi}_{X}" in source
    assert "\\HTQFTclimax^{\\Xi}_{X \\times \\mathbb{R}}" in source
    assert "\\Xi(A,\\omega,\\mathrm{BL})" in source
    assert "Costello--Gwilliam HT BV/factorization model" in flat
    assert "bulk--Hochschild comparison map in the declared ambient" in flat
    assert "\\chi_{\\mathrm{HT},A}\\colon" in source
    assert "\\hypAmbientWtCpl" in source
    assert "\\hypProchazka,\\hypCKL,\\hypPRSh,\\hypYamada" in source
    assert "Without \\(\\Xi\\) the theorem asserts" in flat
    assert "not a bare existence theorem" in flat


def test_pointwise_climax_theorem_is_not_ungated_canonical_existence():
    source = _source()
    flat = _flat()

    assert "There exists a canonical three-" not in source
    assert "There exists a canonical three-dimensional" not in flat
    assert "Assume that \\(A\\) is equipped with the pointwise" in flat
    assert "universal-holography datum \\(\\Xi(A)\\)" in flat
    assert "\\eta^\\partial_A\\colon" in source
    assert "\\chi_{\\mathrm{HT},A}\\) identifies" in flat
    assert "\\mathcal A(A)" in source
    assert "not a bare existence theorem for a canonical physical theory" in flat


def test_functoriality_preserves_reductions_only_when_xi_commutes():
    flat = _flat()

    assert "morphisms preserve the boundary comparison, bulk comparison, and ambient data" in flat
    assert "whenever the corresponding \\(\\Xi\\)-data commute with those operations" in flat
    assert "\\SCchtop\\text{-}\\mathsf{Alg}^{\\omega, \\mathrm{BL},\\Xi}_{X \\times \\R}" in _source()
    assert "\\mathsf{BicNT}^{\\,\\mathsf{SC}\\text{-coh},\\Xi}_{(C, \\R)}" in _source()


def test_sequent_display_and_rungs_are_xi_typed():
    source = _source()
    flat = _flat()

    assert "\\Phihol(A,\\Xi_A)=T_{A,\\Xi}" in source
    assert "\\mathrm{Obs}^{\\partial}(T_{A,\\Xi})\\simeq A" in source
    assert "\\mathrm{Obs}^{\\mathrm{bulk}}(T_{A,\\Xi})\\simeq" in source
    assert "with its realization datum \\(\\Xi_A\\)" in flat
    assert "\\Phihol(H_{k},\\omega_H,\\mathrm{BL}_{H_k},\\Xi_{H_k})" in source
    assert "\\Phihol(V_{k}(\\mathfrak{g}),\\omega_{\\mathrm{Sug}},\\mathrm{BL}_{V_k},\\Xi_{V_k(\\mathfrak g)})" in source
    assert "\\Phihol(\\Vir_{c},\\omega_{\\Vir}," in source
    assert "\\mathrm{BL}_{\\mathrm{BH}},\\Xi_{\\Vir_c}^{\\mathrm{BH}})" in source
    assert "\\Phihol(W_{N,c},\\omega_{W_N}," in source
    assert "\\Phihol(W_{\\infty}[\\lambda]," in source
    assert "Specialise \\cref{thm:universal-holography-master} at" not in source
    assert "\\Phihol(H_{k}))" not in source
    assert "\\Phihol(\\Vir_{c}))" not in source
    assert "\\Phihol(A))" not in source


def test_finite_zhu_is_not_used_as_massey_bound_or_heisenberg_input():
    programme = _flat()
    intro = " ".join(INTRO_SOURCE.read_text(encoding="utf-8").split())
    preface = " ".join(PREFACE_SOURCE.read_text(encoding="utf-8").split())
    log_wp = " ".join(LOG_WP_SOURCE.read_text(encoding="utf-8").split())
    combined = "\n".join([programme, intro, preface, log_wp])

    retired = (
        "finite-dimensionality forces bounded Massey",
        "finite $\\dim\\mathrm{Zhu}(\\cA)$ forces polynomial-in-$r$ Massey",
        "Every $C_2$-cofinite $\\cA$ is tempered",
        "follows from finite-dim Zhu and the absence of higher Massey products",
        "Rigidity pillar (bounded Zhu + refinement)",
        "Tempering criterion: bounded Zhu implies tempered",
        "Equivalently (contrapositive form): every non-tempered VOA",
        "the non-logarithmic standard landscape FAILS for $\\mathcal{W}(p)$",
    )
    for phrase in retired:
        assert phrase not in combined

    assert "Rigidity pillar (finite amplitude data, not Zhu alone)" in programme
    assert "finite-dimensionality of $\\mathrm{Zhu}(A)$ does not bound" in programme
    assert "no finite-Zhu input is used" in programme
    assert "Rigidity (finite amplitude data + refinement)" in preface
    assert "Finite-Zhu amplitude criterion for tempering" in log_wp
    assert "Finite-dimensionality of $A(\\cA)$ is not itself the bound" in log_wp
    assert "finite-envelope non-logarithmic rows does not apply to $\\mathcal{W}(p)$" in programme


def test_high_genus_gravity_trace_is_scalar_weight_completed_only():
    source = _source()
    flat = _flat()
    start = source.index(r"\label{rem:gravity-climax-genus-g-through-F5}")
    end = source.index(r"\section{What universal holography does not say}", start)
    block = source[start:end]
    block_flat = " ".join(block.split())

    assert "High-genus scalar trace closure routes through F5 for the original complex" in block
    assert "only as a scalar-shadow trace comparison in the weight-completed category" in block_flat
    assert "\\(\\hypAmbientWtCpl+\\effScalarShadowProj\\)" in block
    assert "This is not a statement on the original bar complex" in block_flat
    assert "Open Frontier~F5" in block

    retired = (
        "High-genus partition function closure routes through F5",
        "genus $g \\ge 2$ in the weight-completed category unconditionally",
    )
    for phrase in retired:
        assert phrase not in flat
