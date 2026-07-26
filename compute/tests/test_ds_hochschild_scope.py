"""Source guards for the DS--Hochschild bridge.

The bridge is a conditional coderivation-transfer theorem.  It is not a
consequence of generic C2-cofiniteness, and it is not obtained by
applying Hochschild cochains to a DS algebra retract.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHD = ROOT / "chapters" / "theory" / "chiral_higher_deligne.tex"
MAIN = ROOT / "main.tex"
HOCHSCHILD = ROOT / "chapters" / "connections" / "hochschild.tex"
THQG_RECON = ROOT / "chapters" / "connections" / "thqg_holographic_reconstruction.tex"
UCH = ROOT / "chapters" / "connections" / "universal_celestial_holography.tex"
UHF = ROOT / "chapters" / "connections" / "universal_holography_functor.tex"
CLIMAX = ROOT / "chapters" / "connections" / "programme_climax_platonic.tex"


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return " ".join(text.split())


def _ds_block() -> str:
    text = _source(CHD)
    start = text.index(r"\label{thm:chd-ds-hochschild}")
    end = text.index(r"\begin{remark}[Subregular and minimal nilpotents]", start)
    return text[start:end]


def _summary_block() -> str:
    text = _source(CHD)
    start = text.index(r"\item \textbf{thm:chd-ds-hochschild}")
    end = text.index(r"\section{Proof obligations not discharged here}", start)
    return text[start:end]


def _hochschild_duplicate_block() -> str:
    text = _source(HOCHSCHILD)
    start = text.index(r"\begin{proposition}[DS--Hochschild class dichotomy")
    end = text.index(r"\begin{remark}[HKR scope and chiral analogue]", start)
    return text[start:end]


def _uhf_ds_block() -> str:
    text = _source(UHF)
    start = text.index(r"\begin{proposition}[DS-Hochschild compatibility under completed")
    end = text.index(r"\section{Monster finite-orbifold BV datum}", start)
    return text[start:end]


def _uch_gravity_block() -> str:
    text = _source(UCH)
    start = text.index(r"\label{thm:uch-gravity-chain-level}")
    end = text.index(r"\subsubsection*{(d) Yang--Mills", start)
    return text[start:end]


def _thqg_class_m_block() -> str:
    text = _source(THQG_RECON)
    start = text.index(r"\label{thm:class-M-chain-bulk}")
    end = text.index(r"\begin{remark}[Admissible construction]", start)
    return text[start:end]


def _thqg_functoriality_block() -> str:
    text = _source(THQG_RECON)
    start = text.index(r"\emph{Step} (iii) \emph{-- functoriality and DS intertwining.}")
    end = text.index(r"\begin{corollary}[Class~$\mathsf{M}$", start)
    return text[start:end]


def _climax_deletion_block() -> str:
    text = _source(CLIMAX)
    start = text.index("DS--Hochschild all nilpotents")
    end = text.index("chain-level associator-free mixed structure", start)
    return text[start:end]


def _climax_rung3_block() -> str:
    text = _source(CLIMAX)
    start = text.index(r"\label{cor:rung-virasoro}")
    end = text.index(r"\begin{corollary}[\label{cor:rung-w-N}", start)
    return text[start:end]


def _generic_level_is_not_lisse(level: str) -> bool:
    """Return whether the level should be treated as a lisse/C2-finite input."""
    return level in {"nondegenerate-admissible", "positive-integral-rational"}


def test_ds_hochschild_uses_coderivation_transfer_package():
    block = _flat(_ds_block())

    assert "completed DS bar-coalgebra SDR on the reduced chiral bar coalgebras" in block
    assert r"\Coder_0(B^{\mathrm{ch}}(-))" in block
    assert "bar-coalgebra SDR in \\textup{(b)} induces the displayed coderivation SDR" in block
    assert "Hochschild cochains are not functorial in an algebra map" in block
    assert "is not obtained by applying \\(\\ChirHoch^\\bullet(-,-)\\) to the DS algebra retract" in block
    assert "completed DS homotopy equivalence of reduced chiral bar coalgebras" in block
    assert "coderivation HPL give the displayed Hochschild SDR" in block


def test_ds_hochschild_convergence_does_not_use_generic_c2_cofiniteness():
    block = _flat(_ds_block())
    summary = _flat(_summary_block())

    assert "Finiteness of the HPL sum does not come from generic \\(C_2\\)-cofiniteness" in block
    assert "bounded-shift estimate in hypothesis \\textup{(c)}" in block
    assert "direct-sum complex one still needs the corresponding Banach or bounded-shift estimate" in block
    assert "completed HPL with the bounded-shift/weightwise-convergence estimate" in summary
    assert "not the source of functoriality" in summary
    assert not _generic_level_is_not_lisse("generic-irrational")
    assert _generic_level_is_not_lisse("nondegenerate-admissible")


def test_retired_ds_hochschild_false_routes_are_absent_from_live_statement():
    block = _ds_block()
    summary = _summary_block()
    duplicate = _hochschild_duplicate_block()
    main = _source(MAIN)
    thqg = _thqg_class_m_block()
    thqg_functoriality = _thqg_functoriality_block()
    uch = _uch_gravity_block()
    uhf = _uhf_ds_block()
    climax = "\n".join([_climax_deletion_block(), _climax_rung3_block()])

    retired = (
        "Applying\n\\(\\ChirHoch^\\bullet(-,-)\\) and tensoring with the ghost Fock space gives",
        "Arakawa's \\(C_2\\)-cofiniteness\nensures that each conformal-weight stratum is finite-dimensional",
        "Arakawa 2015 $C_2$-cofiniteness",
        "presentation on the $C_2$-cofinite locus",
        "unconditional cohomological core via the chiral Hochschild--DS",
        "retract plus Arakawa $C_2$-cofiniteness",
        "cohomologically by Arakawa $C_2$-cofiniteness",
        "by Arakawa $C_2$-cofiniteness lifted through HPL transfer",
        "HPL transfer through the de Boer-Tjin DS strong deformation",
        "leveraging Arakawa $C_2$-cofiniteness and chiral HKR",
        "fibres such as\nBershadsky--Polyakov are verified",
        r"Hochschild functoriality $C^\bullet_{\mathrm{ch}}(\varphi,\varphi)$",
        "chiral Hochschild functoriality",
    )
    live = "\n".join([block, summary, duplicate, main, thqg, thqg_functoriality, uch, uhf, climax])
    for phrase in retired:
        assert phrase not in live

    assert "conditional cohomological core via the chiral Hochschild--DS" in main
    assert "completed Drinfeld--Sokolov bar-coalgebra SDR" in main


def test_hochschild_duplicate_bridge_has_same_scope_discipline():
    duplicate = _flat(_hochschild_duplicate_block())

    assert "ClaimStatusProvedHere as a conditional implication" not in duplicate
    assert "ClaimStatusProvedHere{} as a class-stratified criterion" in duplicate
    assert "smooth chiral-HKR/derived-vacuum identification" in duplicate
    assert "completed DS coderivation transfer package" in duplicate
    assert "Generic non-critical level is not a $C_2$-cofinite input" in duplicate
    assert "finite-weight Kazhdan/PBW profile" in duplicate
    assert "no generic \\(C_2\\)-cofiniteness assertion is used" in duplicate


def test_holography_surfaces_inherit_completed_bar_coderivation_scope():
    thqg = _flat(_thqg_class_m_block())
    thqg_functoriality = _flat(_thqg_functoriality_block())
    uch = _flat(_uch_gravity_block())
    uhf = _flat(_uhf_ds_block())
    climax = _flat("\n".join([_climax_deletion_block(), _climax_rung3_block()]))

    assert "smooth chiral-HKR/derived-vacuum package" in thqg
    assert "completed DS coderivation transfer" in thqg
    assert "Kazhdan/PBW finite-weight profile" in thqg
    assert "generic non-critical level is not being used as a \\(C_2\\)-cofinite hypothesis" in thqg
    assert "not automatic functoriality of cochains along the algebra map \\(\\varphi\\)" in thqg_functoriality

    assert "completed bar-coalgebra/coderivation-transfer package" in uch
    assert r"\Coder_0(B^{\mathrm{ch}}(-))" in uch
    assert "not obtained by applying Hochschild cochains to the DS algebra retract" in uch
    assert "bounded-shift/weight-completed estimate" in uch
    assert "Vol~I PBW/Koszul class-$\\mathbf L$ finite-propagation theorem" in uch
    assert "not by a lisse \\(C_2\\)-finite hypothesis" in uch

    assert "ClaimStatusProvedHereConditional" not in uhf
    assert (
        "\\ClaimStatusProvedHere{} under the completed SDR/HPL hypotheses"
        in uhf
    )
    assert "completed DS bar-coalgebra SDR on reduced chiral bar coalgebras" in uhf
    assert r"\Coder_0(B^{\mathrm{ch}}(-))" in uhf
    assert "smooth chiral-HKR/derived-vacuum comparison" in uhf
    assert "bounded-shift HPL convergence estimates" in uhf

    assert "completed bar-coalgebra/coderivation-transfer package" in climax
    assert "bounded-shift HPL package is supplied" in climax
