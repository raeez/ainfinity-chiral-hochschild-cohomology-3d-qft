from pathlib import Path
import sys

from sympy import S, Symbol, expand, factor, symbols


ROOT = Path(__file__).resolve().parents[2]
GRAVITY = ROOT / "chapters" / "connections" / "3d_gravity.tex"

sys.path.insert(0, str(ROOT / "compute" / "lib"))


def between(source: str, start: str, end: str) -> str:
    assert start in source, f"missing start marker: {start}"
    tail = source.split(start, 1)[1]
    assert end in tail, f"missing end marker: {end}"
    return start + tail.split(end, 1)[0]


def test_gravity_m4_statement_is_gauge_fixed_and_licensed():
    source = GRAVITY.read_text()
    section = between(
        source,
        r"\begin{proposition}[Quaternary Virasoro $\Ainf$ operation in the",
        r"\subsubsection*{The quinary operation $m_5$",
    )
    flat = " ".join(section.split())

    required = (
        "Stasheff gauge; \\ClaimStatusProvedHere; licensing",
        r"$\alpha+\gamma+\varepsilon$ via consecutive collision coordinates",
        r"$\hypAmbientWtCpl$",
        r"BRST contraction \(h_{\mathrm{Vir}}\)",
        r"m_4^{h_{\mathrm{Vir}}}(T,T,T,T)",
        r"=-h_{\mathrm{Vir}}\mathrm{Obs}_4(T,T,T,T)",
        "chain-level gauge-fixed formula",
        "not a claim that the arity-four minimal identity canonically determines a cohomology operation",
        "source-tree witness is `compute/lib/symbolic_stasheff.py`",
        "`stasheff_rhs_arity4` and `m4_virasoro_symbolic`",
        "`compute/field_sector_generating_function.py::m4_factorization_verify`",
    )
    for needle in required:
        assert needle in flat

    retired = (
        r"\begin{proposition}[Quaternary Virasoro $\Ainf$ operation]",
        r"The quaternary $\Ainf$ operation for the Virasoro algebra is",
        r"Solve the degree-$4$ Stasheff relation~\eqref{eq:gravity-associator}",
        "verified by direct symbolic computation.",
    )
    for needle in retired:
        assert needle not in section


def test_symbolic_stasheff_m4_matches_displayed_formula():
    from symbolic_stasheff import m4_virasoro_symbolic

    l1, l2, l3, c = symbols("l1 l2 l3 c")
    result = m4_virasoro_symbolic(l1, l2, l3, c)
    m4 = result["m4"]
    common = (l1 - l3) * (l1 - l2 + l3)
    sigma = l1 + l2 + l3

    assert expand(m4.get("dT", S.Zero) - 2 * common) == 0
    assert expand(m4.get("T", S.Zero) - 4 * common * sigma) == 0
    assert expand(m4.get("1", S.Zero) - c * common * sigma**3 / 6) == 0
    assert m4.get("d2T", S.Zero) == 0
    assert m4.get("d3T", S.Zero) == 0
    assert expand(factor(m4["T"]) - 4 * (l1 - l3) * (l1 - l2 + l3) * sigma) == 0
    assert result["consistency"]["stasheff_gauge_representative"] is True
    assert "no homotopy needed" not in result["consistency"]["note"]
