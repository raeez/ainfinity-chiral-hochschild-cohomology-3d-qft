"""Source guards for the two-coloured Boardman--Vogt homotopy.

The Research Paper Refinement pass requires the missing two-coloured
coherence in Chiral Higher Deligne to be an explicit signed tree
homotopy, not an unnamed appeal to higher coherence.
"""
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHD = ROOT / "chapters" / "theory" / "chiral_higher_deligne.tex"
FOUNDATIONS = ROOT / "chapters" / "theory" / "foundations.tex"
FM_BOUNDARY = ROOT / "compute" / "lib" / "fm_boundary.py"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_main_theorem_advertises_conditional_two_coloured_scope():
    text = _read(CHD)

    assert "conditional two-coloured lift" in text
    assert r"\ClaimStatusConditional on clauses" in text
    assert r"licensing tags" in text
    assert r"$\alpha+\beta+\gamma+\delta$" in text
    assert r"\textup{(2)} is conditional on the signed two-coloured cobar" in text
    assert (
        r"contraction of Proposition~\ref{prop:chd-two-coloured-cobar-obstruction}"
        in text
    )


def test_signed_boardman_vogt_tree_formula_is_explicit():
    text = _read(CHD)

    assert "def:chd-two-coloured-cobar-tree-homotopy" in text
    assert r"\mathcal C_{\mathrm{oc}}:=\Omega((\SCchtop)^!)" in text
    assert r"\mathfrak{o}(T)=\det(\Bbbk^{E_{\mathrm{int}}(T)})" in text
    assert r"degree \(-1\) two-coloured homotopy" in text
    assert r"\label{eq:chd-hoc-tree-formula}" in text
    assert r"h_{\mathrm{oc}}^{\mathrm{tree}}(T)" in text
    assert r"\sum_{e\in E_{\mathrm{int}}(T)}" in text
    assert r"(-1)^{\sigma(e,T)}T_e" in text
    assert r"\sigma(e_j,T)" in text
    assert r"(j-1)+\sum_{v\prec e_j}|x_v|\pmod 2" in text


def test_open_to_closed_void_and_orientation_ordering_are_guarded():
    text = _read(CHD)

    assert r"\mathcal C_{\mathrm{oc}}(\mathsf{cl}^{p},\mathsf{op}^{q};\mathsf{cl})" in text
    assert r"\bigr)=0\qquad(q>0)." in text
    assert "Changing the total order changes both" in text
    assert "orientation-line expression rather than an ordering-dependent" in text


def test_retract_equations_and_obstruction_class_are_precise():
    text = _read(CHD)

    assert "prop:chd-two-coloured-cobar-obstruction" in text
    assert r"p\circ i=\mathrm{id}_{\mathcal C_{\mathrm{min}}}" in text
    assert r"\label{eq:chd-hoc-contraction}" in text
    assert r"\mathrm{id}-i\circ p" in text
    assert r"\label{eq:chd-hoc-side-conditions}" in text
    assert r"p\,h_{\mathrm{oc}}^{\mathrm{tree}}=0" in text
    assert r"h_{\mathrm{oc}}^{\mathrm{tree}}\,i=0" in text
    assert r"\bigl(h_{\mathrm{oc}}^{\mathrm{tree}}\bigr)^2=0" in text
    assert r"\label{eq:chd-hoc-delta-mix}" in text
    assert r"\Delta_{\mathrm{mix}}h_{\mathrm{oc}}^{\mathrm{tree}}(x)" in text
    assert r"\label{eq:chd-oc-obstruction-class}" in text
    assert r"H^{1}\operatorname{Cone}" in text
    assert r"\End(Z^{\mathrm{der}}_{\mathrm{ch}}(\cA),\cA)" in text


def test_low_arity_tree_sign_calculation_is_printed():
    text = _read(CHD)

    assert "The low-arity calculation fixes the sign convention." in text
    assert r"E_{\mathrm{int}}(T)=\{e_1\}" in text
    assert r"a=\sum_{v\prec e_1}|x_v|" in text
    assert r"h_{\mathrm{oc}}^{\mathrm{tree}}(T)=(-1)^aT_{e_1}" in text
    assert r"d_{\mathrm{cube}}T_{e_1}=(-1)^a\bigl(T-i p(T)\bigr)" in text
    assert r"E_{\mathrm{int}}(T)=\{e_1<e_2\}" in text
    assert r"\sigma(e_1,T)+\sigma(e_2,T_{e_1})" in text
    assert r"\sigma(e_2,T)+\sigma(e_1,T_{e_2})+1" in text
    assert "no new sign occurs in higher arity" in text
    assert "the proof obligation here is" not in text
    assert "extra datum is exactly the two-coloured homotopy" in text


def test_foundations_and_compute_summary_share_the_same_sdr_formula():
    foundations = _read(FOUNDATIONS)
    fm_boundary = _read(FM_BOUNDARY)

    assert r"dh_{\mathrm{oc}}+h_{\mathrm{oc}}d=\mathrm{id}-i\circ p" in foundations
    assert r"p h_{\mathrm{oc}}=h_{\mathrm{oc}}i=h_{\mathrm{oc}}^2=0" in foundations
    assert "d h_oc + h_oc d = id - i circ p" in fm_boundary
    assert "p h_oc = h_oc i = h_oc^2 = 0" in fm_boundary
    assert ("1 - i" " p") not in fm_boundary


def test_ambiguous_retract_notation_is_not_used_on_this_surface():
    text = "\n".join([_read(CHD), _read(FOUNDATIONS), _read(FM_BOUNDARY)])

    assert ("1" "-ip") not in text
    assert ("1 - i" " p") not in text
    assert r"\mathrm{id}-ip" not in text
