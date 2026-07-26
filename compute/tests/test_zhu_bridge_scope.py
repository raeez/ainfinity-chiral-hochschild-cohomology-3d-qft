"""Source guards for the Zhu centre shadow in the Hochschild chapter."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HOCHSCHILD = ROOT / "chapters" / "connections" / "hochschild.tex"


def _source() -> str:
    return HOCHSCHILD.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return " ".join(text.split())


def _zhu_block() -> str:
    text = _source()
    start = text.index(r"\label{subsec:zhu-bridge}")
    end = text.index(r"\subsection{\texorpdfstring{Chiral Hochschild cochains", start)
    return text[start:end]


def _semisimple_hh_profile(number_of_simple_blocks: int) -> dict[int, int]:
    """Hochschild cohomology dimensions for a finite semisimple algebra."""
    return {0: number_of_simple_blocks, 1: 0, 2: 0}


def test_zhu_bridge_is_degree_zero_centre_shadow_not_full_mode_quotient():
    text = _source()
    block = _zhu_block()
    flat = _flat(block)

    assert r"\mathrm{Zhu}^0\colon Z_{\mathrm{ch}}(V)\to Z(A(V))" in text
    assert "not a quotient of the full mode algebra" in block
    assert "not an algebra homomorphism on the full mode algebra" in flat
    assert "requires extra bimodule or derived Morita data" in flat
    assert r"\cA_{\mathrm{mode}}\twoheadrightarrow A(V)" not in block
    assert "quasi-isomorphism in degrees" not in block
    assert "chiral Higher-Deligne concentration" not in block
    assert r"\mathrm{Zhu}^*\colon" not in block
    assert "Zhu centre shadow" in flat


def test_rational_zhu_case_records_semisimple_block_obstruction():
    block = _zhu_block()

    assert r"\HH^0(A(V))=Z(A(V))\simeq \C^{I(V)}" in block
    assert r"\HH^j(A(V))=0\quad(j>0)" in block
    assert r"\C\to\C^3" in block
    assert "It is not a\n quasi-isomorphism" in block

    ising_hh = _semisimple_hh_profile(3)
    assert ising_hh[0] == 3
    assert ising_hh[1] == 0
    assert ising_hh[2] == 0
    assert 1 != ising_hh[0]


def test_critical_affine_zhu_map_is_surjective_with_kernel():
    block = _zhu_block()
    flat = _flat(block)

    assert "Frenkel--Zhu gives \\(A(V)\\simeq U(\\fg)\\)" in flat
    assert "Feigin--Frenkel--Reshetikhin" in block
    assert "It is surjective and has infinite-dimensional kernel" in block
    assert "it is not an inclusion" in block
    assert "has infinite codimension" not in block
    assert r"\hookrightarrow Z(U(\fg))" not in block
    assert "cokernel of" not in block


def test_zhu_folklore_remark_excludes_hochschild_quasi_isomorphism():
    block = _zhu_block()
    flat = _flat(block)

    assert "not for a quasi-isomorphism from chiral Hochschild cochains" in flat
    assert "finite semisimple block centre" in flat
    assert "maps onto \\(Z(U(\\fg))\\) with infinite kernel" in flat
    assert "degree-zero centre shadow" in flat
    assert "semisimple block-trace target" in flat
    assert "explicitly supplied derived Morita comparison" in flat
