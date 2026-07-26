"""Source guards for the negative-cyclic/K-theory distinction.

Negative cyclic homology is the target of a Chern character from
K-theory.  It is not equivalent to K-theory under a Bridgeland
orientation.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HOCHSCHILD = ROOT / "chapters" / "connections" / "hochschild.tex"


def _source() -> str:
    return HOCHSCHILD.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return " ".join(text.split())


def _four_object_block() -> str:
    text = _source()
    start = text.index(r"\begin{remark}[Four-object cyclic-Hochschild discipline")
    end = text.index(r"\begin{remark}[Bar $H^\bullet$ versus cyclic", start)
    return text[start:end]


def _stage_stratification_block() -> str:
    text = _source()
    start = text.index(r"\label{thm:cyclic-hochschild-stage-stratification}")
    end = text.index(r"\begin{remark}[Construction Problem 4", start)
    return text[start:end]


def _base_field_profiles() -> dict[str, str]:
    """Algebraic K-theory and negative cyclic profiles for A = C."""
    return {
        "K0(C)": "abelian-group:Z",
        "HCminus0(C)": "C-vector-space",
        "K1(C)": "multiplicative-group:Cx",
        "HCminus1(C)": "zero",
    }


def test_negative_cyclic_receives_chern_character_not_equivalence():
    four = _flat(_four_object_block())
    stage = _flat(_stage_stratification_block())

    assert "negative cyclic target of the Chern character from $K$-theory" in four
    assert "it does not make $\\mathrm{HC}^{-}(A)$ equal to $K$-theory" in four
    assert r"\operatorname{ch}^{-}_{\sigma}" in stage
    assert (
        r"K^{\mathrm{ch}}\bigl(\mathrm{Mod}(A)\bigr) "
        r"\;\xrightarrow{\;\operatorname{ch}^{-}_{\sigma}\;}\; "
        r"\mathrm{HC}^{-}(A)[\sigma]"
        in stage
    )
    assert "not a level-$\\mathsf{A}$ object" in stage
    assert "not the $K$-theory object itself" in stage


def test_base_field_obstruction_to_absolute_equivalence_is_recorded():
    stage = _flat(_stage_stratification_block())
    profiles = _base_field_profiles()

    assert "Already for $A=\\C$ in algebraic $K$-theory" in stage
    assert r"K_0(\C)=\Z" in stage
    assert r"\mathrm{HC}^{-}_0(\C)" in stage
    assert r"K_1(\C)=\C^\times" in stage
    assert r"\mathrm{HC}^{-}_1(\C)=0" in stage
    assert profiles["K0(C)"] != profiles["HCminus0(C)"]
    assert profiles["K1(C)"] != profiles["HCminus1(C)"]


def test_retired_negative_cyclic_k_theory_equivalence_phrases_are_absent():
    text = _source()
    stage = _stage_stratification_block()

    retired = [
        "computes $K$-theory",
        "computes K-theory",
        r"\mathrm{HC}^{-}(A)[\sigma] " r"\;\xrightarrow{\;\sim\;}\;",
        r"identifies $\mathrm{HC}^{-}$ with chiral $K$-theory",
        "K-theoretic Hochschild-Bridgeland equivalence",
        r"\cite{BZN12} \S~6",
        r"\cite{HKKP11}",
    ]
    for phrase in retired:
        assert phrase not in text

    assert "BZN12" not in stage
    assert "HKKP11" not in stage
    assert "Bridgeland-oriented module category and its\n$K$-theory object" in stage
