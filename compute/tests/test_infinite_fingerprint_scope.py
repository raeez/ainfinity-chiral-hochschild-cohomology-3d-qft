"""Source guards for the six-slot fingerprint reconstruction scope."""
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _read(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def _flat(rel: str) -> str:
    return " ".join(_read(rel).split())


def test_six_slot_fingerprint_requires_branchwise_reconstruction_datum():
    source = _read("chapters/theory/infinite_fingerprint_classification.tex")
    flat = " ".join(source.split())

    assert r"\label{def:fingerprint-branch-reconstruction-datum}" in source
    assert r"\mathfrak R_c\colon" in source
    assert "The datum is not recovered from the six scalar slots alone" in flat
    assert "Six-slot branchwise reconstruction criterion" in source
    assert r"\ClaimStatusConditional" in source
    assert "Without the datum \\(\\mathfrak R_c\\)" in flat
    assert "not an isomorphism theorem" in flat

    assert "unique carrier of a complete fingerprint" not in source
    assert "If $\\varphi'(A_b) = \\varphi'(A_{b'})$, then $A_b \\simeq A_{b'}$" not in source
    assert r"\begin{proof}[Sketch]" not in source


def test_mkosz_uses_conditional_branchwise_coordinates():
    source = _read("chapters/theory/koszulness_moduli_M_kosz.tex")
    flat = " ".join(source.split())

    assert "Shadow-class stratification with branchwise coordinates" in source
    assert r"\ClaimStatusConditional" in source
    assert "branchwise coordinates only after the reconstruction datum" in flat
    assert r"\mathfrak R_c" in source
    assert "six slots become reconstructive only after" in flat

    assert "complete coordinates on each" not in source
    assert "is a complete coordinate on each component" not in source


def test_fingerprint_summaries_do_not_export_global_completeness():
    paths = (
        "chapters/frame/preface.tex",
        "chapters/frame/part_viii_synthesis.tex",
        "chapters/theory/introduction.tex",
        "main.tex",
    )
    combined = "\n".join(_flat(path) for path in paths)

    assert "branchwise reconstruction invariant" in combined
    assert "fingerprints becomes reconstructive only after the branchwise datum" in combined
    assert "conditional on the class-\\(\\mathbf M\\) reconstruction datum" in combined
    assert "branchwise datum \\(\\mathfrak R_c\\) is supplied" in combined

    stale_phrases = (
        "complete invariant of the Koszul-bar-complex structure",
        "equality of fingerprints forces isomorphism of chiral algebras",
        "sixth slot restores completeness",
        "Six-slot fingerprint completeness",
    )
    for phrase in stale_phrases:
        assert phrase not in combined


def test_standard_five_slot_bar_quillen_theorem_remains_separate():
    source = _read("chapters/examples/examples-complete-proved.tex")
    flat = " ".join(source.split())

    assert r"\label{thm:fingerprint-completeness-standard}" in source
    assert "bar complexes are Quillen equivalent" in flat
    assert "no bar-Quillen-completeness statement is asserted outside" in flat
