from pathlib import Path
import sys

from sympy import S, expand, symbols


ROOT = Path(__file__).resolve().parents[2]
GRAVITY = ROOT / "chapters" / "connections" / "3d_gravity.tex"

sys.path.insert(0, str(ROOT / "compute" / "lib"))


def between(source: str, start: str, end: str) -> str:
    assert start in source, f"missing start marker: {start}"
    tail = source.split(start, 1)[1]
    assert end in tail, f"missing end marker: {end}"
    return start + tail.split(end, 1)[0]


def test_gravity_m5_statement_is_gauge_fixed_licensed_and_complete():
    source = GRAVITY.read_text()
    section = between(
        source,
        r"\begin{proposition}[Quinary Virasoro $\Ainf$ operation in the",
        r"\subsubsection*{The degree pattern",
    )
    flat = " ".join(section.split())

    required = (
        r"Stasheff gauge; \ClaimStatusProvedHere; licensing",
        r"$\alpha+\gamma+\varepsilon$ via consecutive collision coordinates",
        r"$\hypAmbientWtCpl$",
        r"BRST contraction",
        r"m_5^{h_{\mathrm{Vir}}}(T^5;",
        r"where the non-leading polynomials are",
        r"P_1(\ell)=",
        r"P_0(\ell)=",
        r"Q_5(\ell)=",
        r"Q_5(\lambda,\ldots,\lambda)=-60\lambda^6",
        r"\mathrm{Obs}_5",
        r"m_5^{h_{\mathrm{Vir}}}=-h_{\mathrm{Vir}}\mathrm{Obs}_5",
        "source-tree witness is `compute/lib/symbolic_stasheff.py`",
        "`stasheff_rhs_arity5`",
        "`m5_virasoro_symbolic`",
        "`compute/tests/test_gravity_m5_stasheff_gauge.py`",
    )
    for needle in required:
        assert needle in flat

    coefficient_markers = (
        r"2\ell_1^3 - 15\ell_1^2\ell_2 + 5\ell_1^2\ell_4",
        r"- 8\ell_3^2\ell_4 - 3\ell_4^3",
        r"-4\ell_1^4 - 8\ell_1^3\ell_2",
        r"- 2\ell_4^4",
        r"-2\ell_1^6 - 6\ell_1^5\ell_2",
        r"- 3\ell_3\ell_4^5 - \ell_4^6",
    )
    for needle in coefficient_markers:
        assert needle in section

    retired = (
        r"\begin{proposition}[Quinary Virasoro $\Ainf$ operation]",
        r"The quinary $\Ainf$ operation for the Virasoro algebra is",
        r"where $P_1$, $P_0$ are explicit polynomials",
        r"Verified by symbolic computation.",
        r"Same method as Proposition~\ref{prop:gravity-m4}",
    )
    for needle in retired:
        assert needle not in section


def test_symbolic_stasheff_m5_matches_displayed_formula():
    from symbolic_stasheff import m5_virasoro_symbolic

    l1, l2, l3, l4, c = symbols("l1 l2 l3 l4 c")
    result = m5_virasoro_symbolic(l1, l2, l3, l4, c)
    m5 = result["m5"]

    P1 = (
        2 * l1**3
        - 15 * l1**2 * l2
        + 5 * l1**2 * l4
        - 22 * l1 * l2**2
        - 20 * l1 * l2 * l3
        - 38 * l1 * l2 * l4
        - 4 * l1 * l3**2
        - 6 * l1 * l3 * l4
        - 6 * l1 * l4**2
        - 4 * l2**3
        - 6 * l2**2 * l3
        - 12 * l2**2 * l4
        - 3 * l2 * l3**2
        - l2 * l3 * l4
        - 9 * l2 * l4**2
        - 4 * l3**3
        - 8 * l3**2 * l4
        - 3 * l4**3
    )
    P0 = (
        -4 * l1**4
        - 8 * l1**3 * l2
        - 8 * l1**3 * l3
        + 12 * l1**3 * l4
        - 6 * l1**2 * l2**2
        - 12 * l1**2 * l2 * l3
        - 12 * l1**2 * l2 * l4
        + 12 * l1**2 * l3 * l4
        + 10 * l1**2 * l4**2
        - 12 * l1 * l2**2 * l3
        - 26 * l1 * l2**2 * l4
        - 16 * l1 * l2 * l3 * l4
        - 22 * l1 * l2 * l4**2
        + 8 * l1 * l3**3
        - 8 * l1 * l3**2 * l4
        - 4 * l1 * l4**3
        - 6 * l2**3 * l3
        - 2 * l2**3 * l4
        - 6 * l2**2 * l3**2
        + 6 * l2**2 * l3 * l4
        - 4 * l2**2 * l4**2
        + 2 * l2 * l3**3
        + 8 * l2 * l3 * l4**2
        - 4 * l2 * l4**3
        + 4 * l3**4
        - 10 * l3**3 * l4
        - 8 * l3**2 * l4**2
        - 2 * l3 * l4**3
        - 2 * l4**4
    )
    Q5 = (
        -2 * l1**6
        - 6 * l1**5 * l2
        - 8 * l1**5 * l3
        - 4 * l1**5 * l4
        - 4 * l1**4 * l2**2
        - 18 * l1**4 * l2 * l3
        - 8 * l1**4 * l2 * l4
        - 10 * l1**4 * l3**2
        - 12 * l1**4 * l3 * l4
        + 2 * l1**3 * l2**3
        - 8 * l1**3 * l2**2 * l3
        - 12 * l1**3 * l2 * l3**2
        - 16 * l1**3 * l2 * l3 * l4
        - 8 * l1**3 * l3**2 * l4
        + 10 * l1**3 * l4**3
        + l1**2 * l2**4
        + 2 * l1**2 * l2**3 * l4
        + 12 * l1**2 * l2 * l3**3
        + 10 * l1**2 * l3**4
        + 8 * l1**2 * l3**3 * l4
        + 10 * l1**2 * l3 * l4**3
        + 5 * l1**2 * l4**4
        - 2 * l1 * l2**4 * l3
        + l1 * l2**4 * l4
        - 4 * l1 * l2**3 * l3**2
        + 4 * l1 * l2**2 * l3**3
        - 14 * l1 * l2**2 * l4**3
        + 16 * l1 * l2 * l3**4
        + 16 * l1 * l2 * l3**3 * l4
        - 8 * l1 * l2 * l3 * l4**3
        - 11 * l1 * l2 * l4**4
        + 8 * l1 * l3**5
        + 14 * l1 * l3**4 * l4
        + 4 * l1 * l3**3 * l4**2
        - 10 * l1 * l3**2 * l4**3
        - 4 * l1 * l3 * l4**4
        - 2 * l1 * l4**5
        - l2**5 * l3
        + l2**5 * l4
        - 4 * l2**4 * l3**2
        + l2**4 * l3 * l4
        + 3 * l2**4 * l4**2
        - 6 * l2**3 * l3**3
        - 4 * l2**3 * l3**2 * l4
        + 6 * l2**3 * l3 * l4**2
        + l2**2 * l3**4
        - 2 * l2**2 * l3**3 * l4
        + 4 * l2**2 * l3 * l4**3
        - 5 * l2**2 * l4**4
        + 5 * l2 * l3**5
        + 8 * l2 * l3**4 * l4
        + 2 * l2 * l3**3 * l4**2
        - 4 * l2 * l3**2 * l4**3
        - 4 * l2 * l3 * l4**4
        - 4 * l2 * l4**5
        + 2 * l3**6
        + 5 * l3**5 * l4
        + 3 * l3**4 * l4**2
        - 8 * l3**3 * l4**3
        - 7 * l3**2 * l4**4
        - 3 * l3 * l4**5
        - l4**6
    )

    assert expand(m5.get("d4T", S.Zero) + 1) == 0
    assert expand(m5.get("d3T", S.Zero) + 4 * l1 + 5 * l2 + 2 * l3 + 3 * l4) == 0
    assert expand(
        m5.get("d2T", S.Zero)
        + (
            2 * l1**2
            + 19 * l1 * l2
            + 6 * l1 * l3
            + 9 * l1 * l4
            + 8 * l2**2
            + 6 * l2 * l3
            + 12 * l2 * l4
            + 3 * l3**2
            + 3 * l3 * l4
            + 3 * l4**2
        )
    ) == 0
    assert expand(m5.get("dT", S.Zero) - P1) == 0
    assert expand(m5.get("T", S.Zero) - P0) == 0
    assert expand(m5.get("1", S.Zero) - c * Q5 / 12) == 0
    assert result["consistency"]["stasheff_gauge_representative"] is True


def test_symbolic_stasheff_m5_symmetric_point_and_numerical_engine_agree():
    from compute.m7_m10_depth_frontier import StasheffEngine
    from symbolic_stasheff import m5_virasoro_symbolic, mk_exact_numerical

    l1, l2, l3, l4, c, lam = symbols("l1 l2 l3 l4 c lam")
    m5 = m5_virasoro_symbolic(l1, l2, l3, l4, c)["m5"]
    subs = {l1: lam, l2: lam, l3: lam, l4: lam}

    assert expand(m5["d4T"].subs(subs) + 1) == 0
    assert expand(m5["d3T"].subs(subs) + 14 * lam) == 0
    assert expand(m5["d2T"].subs(subs) + 71 * lam**2) == 0
    assert expand(m5["dT"].subs(subs) + 154 * lam**3) == 0
    assert expand(m5["T"].subs(subs) + 120 * lam**4) == 0
    assert expand(m5["1"].subs(subs) + 5 * c * lam**6) == 0

    exact = mk_exact_numerical(5, [1, 1, 1, 1], 12)
    assert exact["d4T_coeff"] == -1.0
    assert exact["d3T_coeff"] == -14.0
    assert exact["d2T_coeff"] == -71.0
    assert exact["dT_coeff"] == -154.0
    assert exact["T_coeff"] == -120.0
    assert exact["scalar_coeff"] == -60.0

    engine = StasheffEngine(c_val=12.0)
    numerical = engine.mk((1.0, 1.0, 1.0, 1.0))
    expected_by_order = {4: -1.0, 3: -14.0, 2: -71.0, 1: -154.0, 0: -120.0, -1: -60.0}
    for order, expected in expected_by_order.items():
        assert abs(numerical[order] - expected) < 1e-9


def test_quintic_shadow_paragraph_separates_raw_m5_from_shadow_coefficient():
    source = GRAVITY.read_text()
    section = between(
        source,
        r"\subsubsection*{The quintic shadow:",
        r"\subsubsection*{The shadow table through $S_9$}",
    )
    flat = " ".join(section.split())

    required = (
        "two scalar-shadow routes",
        "operator-level convention check",
        r"$-5c\lambda^6$",
        "it is not the shadow coefficient",
        r"\Pi^{(5)}_{\mathrm{sh}}",
        "that map is not constructed in this computation",
        "The two scalar routes agree",
    )
    for needle in required:
        assert needle in flat

    retired = (
        "three independent routes",
        "three independent verifications",
        r"-(17c/6)\lambda^6",
        "All three routes agree.",
        "division by the Shapovalov norm at levels $\\le 5$",
    )
    for needle in retired:
        assert needle not in section
