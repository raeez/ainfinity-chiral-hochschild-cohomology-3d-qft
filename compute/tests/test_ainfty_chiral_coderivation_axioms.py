import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "chapters/theory/axioms.tex"
MAIN = ROOT / "main.tex"
FOUNDATIONS = ROOT / "chapters/theory/foundations.tex"
ROSETTA = ROOT / "chapters/examples/rosetta_stone.tex"


def read() -> str:
    return TARGET.read_text()


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text)


def test_ainfty_definition_starts_from_chiral_bar_coalgebra():
    text = read()
    body = compact(text)

    assert text.index(r"\mathrm B^{\mathrm{ch}}(A):=T^c") < text.index(
        r"\label{eq:mk-type}"
    )

    required = (
        r"\mathrmB^{\mathrm{ch}}(A):=T^c(s^{-1}\barA)=\bigoplus_{n\ge0}(s^{-1}\barA)^{\otimesn}",
        r"\barBch(A):=\overline{\mathrmB}^{\mathrm{ch}}(A)=\bigoplus_{n\ge1}(s^{-1}\barA)^{\otimesn}",
        r"The$n=0$summandis$k$",
        r"allcoderivationsbelowpreserve$\barBch(A)$",
    )
    for needle in required:
        assert needle in body


def test_ainfty_coderivation_records_corestrictions_and_suspended_sign():
    body = compact(read())

    required = (
        r"m_k\colonA^{\otimesk}\longrightarrowA((\lambda_1))\cdots((\lambda_{k-1}))",
        r"\degm_k=2-k",
        r"Theunaryoperationis$m_1=Q$",
        r"\bara_i:=|a_i|-1=|s^{-1}a_i|",
        r"\piD_A\bigl[s^{-1}a_1|\cdots|s^{-1}a_k\bigr]=s^{-1}m_k(a_1,\ldots,a_k)",
        r"\pi\colon\barBch(A)\tos^{-1}\barA",
        r"D_A\bigl[s^{-1}a_1|\cdots|s^{-1}a_n\bigr]",
        r"\sum_{\substack{k\ge1\\1\lei\len-k+1}}",
        r"(-1)^{\epsilon(i,k)}",
        r"\epsilon(i,k)=\sum_{r<i}\bara_r",
    )
    for needle in required:
        assert needle in body


def test_ainfty_structure_is_square_zero_coalgebra_coderivation():
    body = compact(read())

    required = (
        r"\DeltaD_A=(D_A\otimes\id+\id\otimesD_A)\Delta",
        r"Thesingleequation$D_A^2=0$",
        r"Equation~\eqref{eq:axioms-ainfty-coderivation}isthesuspendedbarconvention",
        r"thesuspendedcoderivationequation$D_A^2=0$",
    )
    for needle in required:
        assert needle in body


def test_rectification_uses_reduced_and_completed_reduced_bar_coalgebra():
    body = compact(read())

    required = (
        r"completed}reducedcofreecoalgebra$\widehat{\barBch}(A)=\prod_{n\ge1}(s^{-1}\barA)^{\otimesn}$",
        r"\mathcal{C}_A\;:=\;\barBch(A)\;=\;\bigoplus_{n\ge1}(s^{-1}\bar{A})^{\otimesn}",
        r"thecoaugmentedcoalgebrais$\mathrmB^{\mathrm{ch}}(A)=k\oplus\mathcalC_A$",
    )
    for needle in required:
        assert needle in body

    retired = (
        r"completed}cofreecoalgebra$\widehat{T}^c(s^{-1}\bar{A})",
        r"\mathcal{C}_A\;:=\;T^c\!\bigl(s^{-1}\bar{A}\bigr)\;=\;\bigoplus_{n\ge1}",
    )
    for needle in retired:
        assert needle not in body


def test_spectral_substitution_remains_part_of_main_stasheff_identity():
    body = compact(read())

    required = (
        r"\label{eq:ainfty-relation-raw}",
        r"m_j(a_{s+1},\dots,a_{s+j}\,;\,\lambda_{s+1},\dots,\lambda_{s+j-1})",
        r"\Lambda_{[s+1,\dots,s+j]}",
        r"\bm\lambda^{\mathrm{out}}=\bigl(\lambda_1,\dots,\lambda_s,\\Lambda_{[s+1,\dots,s+j]},\\lambda_{s+j+1},\dots,\lambda_{n-1}\bigr)",
    )
    for needle in required:
        assert needle in body


def test_sesquilinearity_uses_translation_modes_not_formal_lambda_derivatives():
    text = read()
    body = compact(text)

    required = (
        r"$m_k(\ldots,\partiala_k)=(\partial+\sum_{j=1}^{k-1}j\,\lambda_j)\,m_k(\ldots)$",
        r"m_k(a_1,\ldots,a_{k-1},\partiala_k)&=\Bigl(\partial+\sum_{j=1}^{k-1}j\,\lambda_j\Bigr)\,m_k(a_1,\ldots,a_k)",
        r"m_{k,\mathbfn}(a_1,\ldots,\partiala_i,\ldots,a_k)=-\,n_i\,m_{k,\mathbfn-\mathbfe_i}(a_1,\ldots,a_k)",
        r"m_{k,\mathbfn}(\ldots,\partiala_i,\ldots)=-n_im_{k,\mathbfn-\mathbfe_i}(\ldots)",
        r"(Ta)_{(n)}=-n\,a_{(n-1)}",
        r"Multiplyingbythedividedpowersandsumminggives",
        r"\left(\partial_A+\sum_{r=1}^{k-1}\mu_r\right)m_k(a_1,\ldots,a_k;\mu)",
        r"\sum_{i=1}^{k-1}\mu_i=\sum_{i=1}^{k-1}\sum_{j=i}^{k-1}\lambda_j=\sum_{j=1}^{k-1}j\,\lambda_j",
        r"Noformalderivative\(\partial_{\lambda_i}\)appearsinthe\(\lambda\)-bracketnormalisation",
    )
    for needle in required:
        assert needle in body

    retired = (
        "Via the spectral parameters:",
        "F(\\varepsilon)",
        "partial_{\\lambda_1} m_k",
        "m_k(\\partial a_1, a_2, \\ldots, a_k)\n = \\partial_{\\lambda_1}",
        "the proof tracks how a shift of one insertion point",
    )
    for needle in retired:
        assert needle not in text


def test_live_summary_surfaces_do_not_collapse_reduced_and_coaugmented_bar():
    main = compact(MAIN.read_text())
    foundations = compact(FOUNDATIONS.read_text())
    rosetta = compact(ROSETTA.read_text())

    main_required = (
        r"Thereducedbarcomplex\[\barBch(\cA)=\bigoplus_{n\ge1}(s^{-1}\bar\cA)^{\otimesn}",
        r"\mathrmB^{\mathrm{ch}}(\cA)=k\oplus\barBch(\cA)",
        r"Proposition~\ref{prop:sc-koszul-duality-nonselfdual}",
        r"\phantomsection\label{V1-cor:def-obs-exchange-genus0}",
        r"\phantomsection\label{V1-rem:critical-level-lie-vs-chirhoch}",
        r"\phantomsection\label{V1-rem:conv-strict-vs-homotopy}",
        r"\phantomsection\label{V1-prop:koszul-dual-tensor-product}",
    )
    for needle in main_required:
        assert needle in main

    foundations_required = (
        r"\mathrmB^{\mathrm{ch}}(A):=\bigoplus_{n\ge0}(s^{-1}\barA)^{\otimesn}",
        r"\barBch(A):=\bigoplus_{n\ge1}(s^{-1}\barA)^{\otimesn}",
        r"\pi_1D_A\big|_{(s^{-1}\barA)^{\otimesk}}",
    )
    for needle in foundations_required:
        assert needle in foundations

    rosetta_required = (
        r"\mathrmB^{\mathrm{ch,ord}}(A)=\bigoplus_{n\ge0}A^{\boxtimesn}|_{\mathrm{FM}_n(C)^{\mathrm{ord}}}",
        r"\barB^{\mathrm{ch,ord}}(A)=\bigoplus_{n\ge1}\barA^{\boxtimesn}|_{\mathrm{FM}_n(C)^{\mathrm{ord}}}",
        r"Thereducedcoidealisfunctorialin$E_1$-chiralmaps",
    )
    for needle in rosetta_required:
        assert needle in rosetta

    retired = (
        r"\barBch(\cA)=T^c(s^{-1}\bar\cA)",
        "thm:SC-" + "self-duality",
        r"\mathrmB^{\mathrm{ch}}(A):=T^c(s^{-1}\barA)",
        r"\barBch(A)=\bigoplus_{n\ge0}A^{\boxtimesn}|_{\mathrm{FM}_n(C)^{\mathrm{ord}}}",
    )
    for needle in retired:
        assert needle not in main
        assert needle not in foundations
        assert needle not in rosetta
