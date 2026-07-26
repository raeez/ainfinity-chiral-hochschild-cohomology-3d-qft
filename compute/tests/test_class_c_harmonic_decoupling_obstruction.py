from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BV_BRST = ROOT / "chapters/connections/bv_brst.tex"
BV_CONSTRUCTION = ROOT / "chapters/theory/bv-construction.tex"
PVA_DESCENT = ROOT / "chapters/theory/pva-descent-repaired.tex"


def _flat(path: Path) -> str:
    return " ".join(path.read_text(encoding="utf-8").split())


def test_class_c_harmonic_decoupling_is_a_named_obstruction_class():
    source = _flat(BV_BRST)

    assert r"\label{def:class-c-harmonic-decoupling-obstruction}" in source
    assert r"\label{eq:class-c-harmonic-decoupling-obstruction}" in source
    assert r"\mathfrak h_{g,r}^{\mathsf C}(\cA)" in source
    assert r"\pi_{\mathrm{ncomp}}\delta_{g,r}^{\mathrm{harm}}" in source
    assert r"\mathcal F_{\mathrm{comp}}^{g,r}(\cA)" in source
    assert r"J_j=:\!\beta\,\partial^j\gamma\!:" in source
    assert r"d_{\operatorname{Hom}}(\phi)" in source
    assert "vanishing and compatible lifting of these classes" in source


def test_class_c_theorem_uses_obstruction_vanishing_not_vague_decoupling():
    source = _flat(BV_BRST)

    theorem_window = source[
        source.index(r"\label{thm:bv-bar-coderived}") :
        source.index(r"\begin{proof}", source.index(r"\label{thm:bv-bar-coderived}"))
    ]
    assert r"\mathfrak h_{g,r}^{\mathsf C}(\cA)" in theorem_window
    assert "vanish compatibly for all \\(r\\ge4\\)" in theorem_window
    assert "and harmonic decoupling holds" not in theorem_window

    proof_window = source[
        source.index("The same argument applies to class~$\\mathsf{C}$") :
        source.index("Now assume $\\cA$ lies in class~$\\mathsf{M}$")
    ]
    assert r"\mathfrak h_{g,r}^{\mathsf C}(\cA)" in proof_window
    assert "non-composite quotient contains no residual harmonic component" in proof_window
    assert "composite-channel representatives are null-homotopic" in proof_window
    assert "under harmonic decoupling" not in proof_window


def test_genus_one_betagamma_result_is_only_the_first_obstruction():
    source = _flat(BV_BRST)
    remark = source[
        source.index(r"\label{rem:bv-bar-class-c-proof}") :
        source.index(r"\label{def:coacyclic-fact}")
    ]

    assert r"\mathfrak h_{1,4}^{\mathsf C}(\beta\gamma)=0" in remark
    assert "genus~$1$" in remark
    assert "all-genera harmonic-decoupling hypothesis" in remark
    assert "for every \\(g\\ge1\\) and every \\(r\\ge4\\)" in source


def test_downstream_surfaces_reference_the_obstruction_class():
    bv = _flat(BV_CONSTRUCTION)
    pva = _flat(PVA_DESCENT)

    assert r"\mathfrak h_{g,r}^{\mathsf C}" in bv
    assert r"\ref{def:class-c-harmonic-decoupling-obstruction}" in bv
    assert "class C/M closes in the weight-completed category" not in bv
    assert "class-C harmonic-decoupling datum" not in bv
    assert "QME after harmonic-decoupling" not in bv
    assert r"QME after vanishing/lifting of \(\mathfrak h_{g,r}^{\mathsf C}\)" in bv
    assert "raw direct-sum chain-level statement is false" in bv

    assert r"\mathfrak h_{g,r}^{\mathsf C}" in pva
    assert r"\ref{def:class-c-harmonic-decoupling-obstruction}" in pva
    assert "classes \\(\\mathfrak h_{g,r}^{\\mathsf C}\\) vanish and lift" in pva
