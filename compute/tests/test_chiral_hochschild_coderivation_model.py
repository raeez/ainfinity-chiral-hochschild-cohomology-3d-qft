import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HOCHSCHILD = ROOT / "chapters" / "connections" / "hochschild.tex"
FOUNDATIONS = ROOT / "chapters" / "theory" / "foundations.tex"
AXIOMS = ROOT / "chapters" / "theory" / "axioms.tex"


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


def test_hochschild_chapter_uses_coderivation_model_and_commutator_differential():
    text = HOCHSCHILD.read_text()
    body = compact(text)

    required = (
        r"B^{\mathrm{ch}}(A)=k\oplus\barBch(A)",
        r"D_A=\sum_{r\geq1}D_{m_r}",
        r"C^\bullet_{\mathrm{ch}}(A,A):=\Coder_0\bigl(B^{\mathrm{ch}}(A)\bigr)",
        r"\simeq\Coder\bigl(\barBch(A)\bigr)",
        r"d_{\mathrm{Hoch}}f=[D_A,f]",
        r"=D_A\circf-(-1)^{|f|}f\circD_A",
        r"C^\bullet_{\mathrm{ch}}(A,A)=\Coder_0(B^{\mathrm{ch}}(A))",
        r"Foran\(\Ainf\)-chiralalgebrathehigher\(m_r\)termsarepartofthedifferential",
        r"\Coder_0\bigl(B^{\mathrm{ch}}(A_\partial)\bigr)",
        r"differential\(d_{\mathrm{Hoch}}=[D_{A_\partial},-]\)",
    )
    for needle in required:
        assert needle in body

    retired = (
        r"C^\bullet_{\text{ch}}(A, A) = \prod_{n \geq 0} \Hom_{\text{VA}}",
        r"where $\Hom_{\text{VA}}$ denotes vertex algebra homomorphisms",
        r"d_{\text{ch}} = \delta_Q + \delta_{\text{Hoch}}",
        "the Hochschild differential inserts $m_2$ at adjacent positions",
    )
    for phrase in retired:
        assert phrase not in text


def test_bulk_hochschild_theorem_is_comparison_not_definition():
    text = HOCHSCHILD.read_text()
    flat = " ".join(text.split())

    required = (
        "Bulk--Hochschild comparison for a chosen physical prefactorization model",
        "\\Psi_{\\mathrm{Ran}}\\colon",
        "\\chi_{\\mathrm{HT},A_\\partial}\\colon",
        "it is not a definition of a physical bulk from \\(A_\\partial\\) alone",
        "The inverse homotopy class is \\(\\chi_{\\mathrm{HT},A_\\partial}\\)",
        "The first equivalence is the \\(\\Psi_{\\mathrm{Ran}}\\) comparison",
    )
    for needle in required:
        assert needle in flat

    retired = (
        "Bulk $=$ chiral Hochschild cochains",
        "bulk local operators form the chiral derived centre",
        "Bulk $=$ derived center",
        "the bulk is identified from the boundary chiral Hochschild cochains",
    )
    for phrase in retired:
        assert phrase not in text


def test_foundations_and_axioms_use_same_bch_coderivation_bridge():
    foundations = FOUNDATIONS.read_text()
    axioms = AXIOMS.read_text()
    f_body = compact(foundations)
    a_body = compact(axioms)

    foundations_required = (
        r"\operatorname{Coder}_0\bigl(\mathrmB^{\mathrm{ch}}(A_b)\bigr)",
        r"d_{\mathrm{Hoch}}=[D_{A_b},-]",
        r"D_{A_b}=\sum_{r\ge1}D_{m_r}",
        r"thefull\(\Ainf\)-chiraldifferentialusesall\(D_{m_r}\)",
    )
    for needle in foundations_required:
        assert needle in f_body

    axioms_required = (
        r"\mathrm{C}^{\bullet}_{\mathrm{ch}}(A,A)\;\simeq\;\mathrm{Coder}_0\bigl(\mathrmB^{\mathrm{ch}}(A)\bigr)",
        r"\simeq\;\mathrm{Coder}\bigl(\barBch(A)\bigr)",
        r"d_{\mathrm{Hoch}}\;=\;[D_A,\,-\,]",
        r"\Coder_0(\mathrmB^{\mathrm{ch}}(A))\simeq\Zderch{A}",
    )
    for needle in axioms_required:
        assert needle in a_body

    retired = (
        r"\mathrm{Coder}\bigl(\BarTwc{A}\bigr),",
        r"d_{\mathrm{Hoch}} \;=\; [d_{B},\, -\,],",
        r"d_{\mathrm{ch}} = \delta_Q + \delta_{\mathrm{Hoch}}^{\mathrm{ch}}",
    )
    for phrase in retired:
        assert phrase not in foundations
        assert phrase not in axioms
