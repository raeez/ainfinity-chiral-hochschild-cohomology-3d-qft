import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FM_PROOFS = ROOT / "chapters/theory/fm-proofs.tex"
ROSETTA = ROOT / "chapters/examples/rosetta_stone.tex"


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


def test_fm4_full_identity_term_count_and_minimal_scope_are_correct():
    text = FM_PROOFS.read_text()
    body = compact(text)

    required = (
        "The ten displayed full-chain terms are supported on the six",
        r"D_{\{1,2,3,4\}}$contributions",
        r"D_{\{1,2,3\}}$and$D_{\{2,3,4\}}$contributions",
        r"D_{\{1,2\}}$,$D_{\{2,3\}}$,$D_{\{3,4\}}$contributions",
        r"m_2\circm_3+m_3\circm_2=0",
        "The nine-term count belongs instead",
        r"m_2\circm_4+m_3\circm_3+m_4\circm_2=0",
    )
    for needle in required:
        haystack = body if "\\" in needle else text
        assert needle in haystack

    retired = (
        "The nine terms correspond to the nine non-vanishing boundary strata",
        r"$D_{\{1,2\}}\cong\FM_3(\C)\times\FM_2(\C)$",
        r"$D_{\{1,2,3\}}\cong\FM_2(\C)\times\FM_3(\C)$",
        r"\epsilon(1,3)=(3-1)|a_1|",
    )
    for phrase in retired:
        assert phrase not in text


def test_arnold_relation_is_printed_as_pva_jacobi_residue():
    body = compact(FM_PROOFS.read_text())

    required = (
        r"\label{cor:arnold-pva-jacobi-residue}",
        r"Ataritythree,theformidentity\eqref{eq:AOS}istheJacobiidentityforthesingularpartof\(m_2\)",
        r"\Omega_3^{\mathrm{sing}\text{-}\mathrm{sing}}(a,b,c)",
        r"D_{\{1,2\}},\qquadD_{\{2,3\}},\qquadD_{\{1,3\}}",
        r"withincidencesigns\(+,+,-\)",
        r"J_{\lambda,\mu}(a,b,c)",
        r"orientedboundarycontributions",
        r"\{\{a_\lambdab\}_{\lambda+\mu}c\}",
        r"\{a_\lambda\{b_\muc\}\}",
        r"-\,(-1)^{(\degree{a}+1)(\degree{b}+1)}",
        r"-\{\{a_\lambdab\}_{\lambda+\mu}c\}",
        r"\{b_\mu\{a_\lambdac\}\}",
        r"Theorem~\ref{thm:Jacobi}",
    )
    for needle in required:
        assert needle in body


def test_rosetta_class_l_does_not_claim_m4_is_killed_directly_by_jacobi():
    text = ROSETTA.read_text()
    body = compact(text)

    required = (
        r"Cohomologicaloperations:$m_3^H\ne0$,$m_4^H=0$",
        r"Thearity-foursource\(m_2\circm_3+m_3\circm_2\)istheLie-Jacobiator",
        r"thisisnotachain-leveldeterminationof\(m_4\)bythearity-fouridentity",
        r"thearity-fourLie-Jacobisourcestillvanishesoncohomology",
        r"nocohomological\(m_4^H\)appears",
    )
    for needle in required:
        assert needle in body

    retired = (
        "The quartic operation $m_4 = 0$ by the Jacobi identity",
        "The quartic operation $m_4 = 0$\nby the Jacobi identity",
        "the Jacobi identity still kills $m_4$",
    )
    for phrase in retired:
        assert phrase not in text
