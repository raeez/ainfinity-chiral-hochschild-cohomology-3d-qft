import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AXIOMS = ROOT / "chapters/theory/axioms.tex"
FM_CALCULUS = ROOT / "chapters/theory/fm-calculus.tex"
FM_PROOFS = ROOT / "chapters/theory/fm-proofs.tex"
ORIENTATIONS = ROOT / "chapters/theory/orientations.tex"


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


def test_cluster_factorization_uses_residue_not_verbal_sign_appeal():
    text = AXIOMS.read_text()
    body = compact(text)

    required = (
        r"\label{eq:cluster-fm-residue}",
        r"\Res_{D_S}(\omega_k)=\beta_S\big|_{\varepsilon_S=0}=\lim_{\varepsilon_S\to0}\bigl(\varepsilon_S\iota_{\partial_{\varepsilon_S}}\omega_k\bigr)\big|_{D_S}",
        r"\label{eq:cluster-fm-orientation}",
        r"\operatorname{or}\bigl(\FM_k(\C)\bigr)=(-d\varepsilon_S)\wedge\operatorname{or}(D_S)",
        r"\label{eq:cluster-collision-differential}",
        r"D_{\mathrm{coll}}=\sum_{\substack{S\subset\{1,\ldots,k\}\\|S|\ge2}}\epsilon_{D_S}\,\Res_{D_S}\otimesm_S",
        r"D_{\mathrm{coll}}^2=0",
        r"\partial^2[\FM_k(\C)]=0",
        r"d\varepsilon_T\wedged\varepsilon_S=-d\varepsilon_S\wedged\varepsilon_T",
    )
    for needle in required:
        assert needle in body

    assert "The naive contraction" in text
    assert "is not used" in text

    retired = (
        "The orientation and sign conventions are those of",
        "up to shuffle signs",
    )
    for phrase in retired:
        assert phrase not in text


def test_cluster_factorization_uses_reduced_inner_fm_strata():
    body = compact(AXIOMS.read_text())

    assert (
        r"D_\pi\cong\FM_r(\C)\times\prod_{j=1}^r\FM_{k_j}^{\mathrm{red}}(\C)"
        in body
    )
    assert (
        r"D_\pi\cong\FM_r(\C)\times\prod_{j=1}^r\FM_{k_j}(\C)"
        not in body
    )


def test_fm_calculus_and_appendix_share_residue_orientation_conventions():
    calculus = compact(FM_CALCULUS.read_text())
    proofs = compact(FM_PROOFS.read_text())

    shared_required = (
        r"\Res_{D_S}(\omega):=\beta\big|_{\varepsilon_S=0}",
        r"\Res_{D_S}(\omega)=\lim_{\varepsilon_S\to0}\bigl(\varepsilon_S\cdot\iota_{\partial_{\varepsilon_S}}\omega\bigr)\big|_{D_S}",
        r"nottheordinarynormalcontraction",
    )
    for needle in shared_required:
        assert needle in calculus

    appendix_required = (
        r"\text{or}(\FM_k(\C))=(-d\varepsilon_S)\wedge\text{or}(D_S)",
        r"\Res_{D_S}(\omega)=\beta\big|_{\varepsilon_S=0}",
        r"\Res_{D_S}(\omega)=\lim_{\varepsilon_S\to0}\bigl(\varepsilon_S\cdot\iota_{\partial_{\varepsilon_S}}\omega\bigr)\big|_{D_S}",
        "naivecontraction",
        "mustnotbeused",
    )
    for needle in appendix_required:
        assert needle in proofs


def test_fm_reduction_notation_and_planar_projection_scope_are_explicit():
    calculus_text = FM_CALCULUS.read_text()
    proofs_text = FM_PROOFS.read_text()
    orientations_text = ORIENTATIONS.read_text()

    calculus = compact(calculus_text)
    proofs = compact(proofs_text)
    orientations = compact(orientations_text)

    required_calculus = (
        r"\FM_k^{\mathrm{tr}}(\C):=\FM_{\C}(k)/\C",
        r"\FM_k^{\mathrm{red}}(\C):=\FM_{\C}(k)/(\C\rtimes\R_{>0})",
        r"\FM_k^{\mathrm{aff}}(\C):=\FM_{\C}(k)/(\C\rtimes\C^\times)",
        r"\dim_\R\FM_k(\C)=2(k-1)",
        r"withrealdimension\(2k-4\)",
        r"\FM_k^{\mathrm{aff}}(\C)",
        "whichwoulderasetheangularresidue",
        r"D_S\cong\FM_r^{\mathrm{tr}}(\C)\times\FM_{|S|}^{\mathrm{red}}(\C)",
        r"codimension-\(2\)stratumof\(\FM_k^{\mathrm{tr}}(\C)\)",
        r"vanishonlyafterapplyingtheopenplanar\(E_1\)-orderedprojection",
        r"TheclosedchiralcollisioncalculusandthedressedPVAdescentretainthenon-consecutivedivisors",
    )
    for needle in required_calculus:
        assert needle in calculus

    required_proofs = (
        r"\FM_3(\C)=\FM_3^{\mathrm{tr}}(\C)",
        r"\FM_4(\C)=\FM_4^{\mathrm{tr}}(\C)",
        r"\FM_2^{\mathrm{tr}}(\C)\times\FM_2^{\mathrm{red}}(\C)",
        r"\FM_1^{\mathrm{tr}}(\C)\times\FM_3^{\mathrm{red}}(\C)",
        r"\FM_1^{\mathrm{tr}}(\C)\times\FM_4^{\mathrm{red}}(\C)",
        r"Thedivisor\(D_{\{1,3\}}\)ispresentintheclosedchiralFMspace",
        r"theplanarresiduevanishes",
    )
    for needle in required_proofs:
        assert needle in proofs

    assert (
        r"\FM_k^{\mathrm{tr}}(\C)=\FM_{\C}(k)/\C" in orientations
    )
    assert "thefullcomplex-affinequotientisnotusedintheStokessignsbelow" in orientations

    retired = (
        "After quotienting $\\C^3$ by translations and dilations",
        r"$D_{\{1,2\}} \cong \FM_2(\C) \times \FM_2^{\mathrm{red}}(\C)$",
        r"$D_{\{2,3\}} \cong \FM_2(\C) \times \FM_2^{\mathrm{red}}(\C)$",
        r"$D_{\{1,2,3\}} \cong \FM_1(\C) \times \FM_3(\C)",
        "This is \\emph{not} a codimension-2 corner",
        "zero by planarity",
        "contribute zero by the ordering constraint",
    )
    combined = "\n".join([calculus_text, proofs_text, orientations_text])
    for phrase in retired:
        assert phrase not in combined
