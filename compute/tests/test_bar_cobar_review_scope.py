"""Guards for the bar-cobar review theorem-package scope.

The review chapter must track the current Vol I theorem package: MC2
uses the completed coinvariant carrier, Theorems A--D keep their
Koszul/coderived/genus/spectral-scalar scope, and the genus-zero bar
differential is not itself Main Theorem A.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BAR_COBAR = ROOT / "chapters" / "connections" / "bar-cobar-review.tex"
HT_BBL_EXTENSIONS = ROOT / "chapters" / "connections" / "thqg_ht_bbl_extensions.tex"


def _source() -> str:
    return BAR_COBAR.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return " ".join(text.split())


def _between(start: str, end: str) -> str:
    text = _source()
    lo = text.index(start)
    hi = text.index(end, lo)
    return text[lo:hi]


def _vol1_projection_block() -> str:
    return _between(
        r"\label{constr:vol2-five-theorems-projections}",
        r"\begin{remark}[Transition to line operators]",
    )


def _bar_differential_block() -> str:
    return _between(
        r"\label{thm:bar_differential_squared}",
        r"\subsubsection{Structural properties of the $\Ainf$ tower}",
    )


def test_bar_cobar_adjunction_label_is_live_and_not_stale():
    text = _source()
    ht_bbl = HT_BBL_EXTENSIONS.read_text(encoding="utf-8")

    assert r"\label{thm:bar_cobar_adjunction}" in text
    assert "thm:bar-cobar-adjunction" not in text
    assert "V1-thm:bar-cobar-adjunction" not in ht_bbl
    assert r"Theorem~\ref{thm:bar_cobar_adjunction}" in ht_bbl


def test_bar_differential_is_genus_zero_input_not_main_theorem_a():
    block = _flat(_bar_differential_block())

    assert "the one-input-$m_1$ part of the full arity-$k$ $A_\\infty$ identity" in block
    assert "not only the binary relation" in block
    assert "genus-zero geometric presentation of the bar coderivation" in block
    assert "genus-zero bar-construction identity $d_B^2=0$" in block
    assert "Main Theorem~A uses this identity inside the bar--cobar adjunction and Verdier lane" in block
    assert "homotopy Koszul dual factorization algebra" in block

    retired = (
        "which follows from the $k=2$ case of the $A_\\infty$ identity",
        "independent geometric verification of Vol~I's algebraic $d_B^2 = 0$",
        "the two results are equivalent, not merely analogous",
    )
    for phrase in retired:
        assert phrase not in block


def test_vol1_projection_block_uses_completed_mc2_carrier():
    block = _flat(_vol1_projection_block())

    assert r"\MC\bigl(G^1(\Defcyc(\cA)\widehat{\otimes}\Gmod)\bigr)" in block
    assert r"G^1(\Defcyc(\cA)\widehat{\otimes}\Gmod)\cong\gAmod" in block
    assert "completed coinvariant modular convolution algebra" in block
    assert "not the uncompleted direct sum" in block
    assert "not the ordered $E_1$ carrier" in block
    assert "ordered lift is a separate theorem before averaging" in block

    retired = (
        r"\Theta_\cA \in \mathrm{MC}(\gAmod)",
        "projections of a single Maurer--Cartan element",
    )
    for phrase in retired:
        assert phrase not in block


def test_vol1_theorem_rows_carry_current_scope_qualifiers():
    block = _flat(_vol1_projection_block())

    assert r"\Barch_X\dashv\Omegach_X" in block
    assert "homotopy Koszul dual factorization algebra" in block
    assert "small classical dual only on the Koszul-effective formal locus" in block
    assert r"\Omegach_X\Barch_X(\cA)\xrightarrow{\sim}\cA" in block
    assert "weight-completed coderived/contraderived Positselski ambient" in block
    assert "not in the ordinary derived category" in block
    assert r"H^*(\overline{\mathcal M}_g,\mathcal Z(\cA))=Q_g(\cA)\oplus Q_g(\cA^!)" in block
    assert "perfect-duality/Lagrangian clause is asserted for $g\\ge1$" in block
    assert "genus zero is the pointed curved base case" in block
    assert r"D_\cA^2=0\Rightarrow\Theta_\cA" in block
    assert "scalar-shadow/determinant-line extraction" in block
    assert "uniform-weight lane, with genus one unconditional" in block

    retired = (
        r"$B_\kappa \dashv \Omega_\kappa$ on $\mathrm{Ran}(X)$",
        r"$\Omega(B(\cA)) \xrightarrow{\sim} \cA$ on Koszul locus",
        r"$Q_g(\cA) + Q_g(\cA^!) = H^*(\overline{\mathcal{M}}_g, Z(\cA))$",
        r"$\kappaChHodge(\cA)$ universal, additive, $\hat{A}$-genus",
        r"\Bbarch",
    )
    for phrase in retired:
        assert phrase not in block


def test_vol1_rectification_remark_records_current_firewall():
    block = _flat(_vol1_projection_block())

    assert r"\label{rem:vol1-theorem-package-rectification}" in block
    assert "bar coalgebra is a twisting-morphism classifier, not the bulk and not the small dual" in block
    assert "Verdier object" in block
    assert "finite-dualizable Koszul-effective/formality hypotheses" in block
    assert "ordinary bar--cobar inversion on the Koszul locus" in block
    assert "completed coderived/contraderived Positselski ambient off it" in block
    assert "$g\\ge1$ Verdier pairing" in block
    assert "pointed curved base" in block
    assert "bar-intrinsic MC element through the scalar-shadow determinant line" in block
    assert "not from a circular appeal to the family index formula" in block
