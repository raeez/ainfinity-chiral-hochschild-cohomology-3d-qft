"""Guards for the three-Hochschild comparison theorem.

The three comparison maps do not identify the three Hochschild-type
invariants in all regimes.  The critical affine Feigin-Frenkel centre is
larger than ordinary mode Hochschild degree zero, and Gel'fand-Fuchs
cohomology is not isomorphic to Weyl Hochschild cohomology.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HOCHSCHILD = ROOT / "chapters" / "connections" / "hochschild.tex"


def _source() -> str:
    return HOCHSCHILD.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return " ".join(text.split())


def _three_way_scope_block() -> str:
    text = _source()
    start = text.index(r"\label{rem:thm-H-three-way-scope}")
    end = text.index(r"\subsection{The Zhu centre shadow", start)
    return text[start:end]


def _three_hochschild_theorem_block() -> str:
    text = _source()
    start = text.index(r"\label{thm:three-hochschild-unification}")
    end = text.index(
        r"\begin{proposition}[Algebraic--geometric chiral Hochschild models are equivalent",
        start,
    )
    return text[start:end]


def _three_invariants_block() -> str:
    text = _source()
    start = text.index(r"\label{rem:three-hochschild-invariants}")
    end = text.index(r"\begin{remark}[Chiral Hochschild comparison at the K3", start)
    return text[start:end]


def test_eta_mode_is_not_critical_degree_zero_isomorphism():
    scope = _flat(_three_way_scope_block())
    theorem = _flat(_three_hochschild_theorem_block())
    combined = scope + " " + theorem

    assert "degree-$0$ isomorphism only on the non-critical scalar-centre locus" in scope
    assert "Feigin--Frenkel zero-mode quotient with infinite-dimensional kernel" in scope
    assert "isomorphism in degree $0$ on the non-critical scalar-centre locus" in theorem
    assert "Feigin--Frenkel zero-mode/one-point specialization" in theorem
    assert "with infinite-dimensional kernel" in theorem
    assert "not an isomorphism in this uncompleted mode-algebra ambient" in theorem
    assert "Feigin--Frenkel--Reshetikhin zero-mode/one-point specialization" in theorem

    retired = (
        "is an isomorphism in degree $0$, always",
        "Degree-$0$ iso and degree-$1$ injection are direct computations",
        "Agreement occurs only on the cohomology of the comparison map",
    )
    for phrase in retired:
        assert phrase not in combined


def test_zeta_is_character_not_weyl_isomorphism():
    theorem = _flat(_three_hochschild_theorem_block())

    assert "WRONG" not in theorem
    assert "Theorem~H has no classical analogue" not in theorem
    assert "opposite shape from a chiral-only uniqueness claim" in theorem
    assert "spectral and higher-configuration data" in theorem
    assert "CE/HKR character" in theorem
    assert r"\HH^\bullet(\cA_{\mathrm{mode}}, \cA_{\mathrm{mode}})=k" in theorem
    assert "identifies only the degree-$0$ quotient" in theorem
    assert "all positive-degree Gel'fand--Fuchs classes lie in its kernel" in theorem
    assert "deformation-quantisation trace character" in theorem
    assert "kills the positive-degree GF generators" in theorem

    retired = (
        r"is an isomorphism when \(\cA_{\mathrm{mode}}\) is a Weyl algebra",
        r"is an isomorphism when $\cA_{\mathrm{mode}}$ is a Weyl algebra",
        "HKR for Weyl algebras is classical; continuous variant",
    )
    for phrase in retired:
        assert phrase not in theorem


def test_three_invariants_remark_uses_positive_scope_not_false_marker():
    remark = _flat(_three_invariants_block())

    assert "has the wrong target: Weyl Hochschild cohomology is concentrated in degree $0$" in remark
    assert "The unbounded object is (iii)" in remark
    assert "Theorem~\\ref{thm:three-hochschild-unification} makes the comparison precise" in remark

    retired = (
        "is FALSE",
        "The genuine ``fails-to-concentrate'' object",
        "while THH is infinite'' is FALSE",
    )
    for phrase in retired:
        assert phrase not in _three_invariants_block()
