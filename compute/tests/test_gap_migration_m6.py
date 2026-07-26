from pathlib import Path

from compute.m7_m10_depth_frontier import StasheffEngine


ROOT = Path(__file__).resolve().parents[2]
GRAVITY = ROOT / "chapters" / "connections" / "3d_gravity.tex"
INTRO = ROOT / "chapters" / "theory" / "introduction.tex"
W3 = ROOT / "chapters" / "examples" / "w-algebras-w3.tex"
W_STABLE = ROOT / "chapters" / "examples" / "w-algebras-stable.tex"
FRONTIER = ROOT / "chapters" / "examples" / "w-algebras-frontier.tex"


def between(source: str, start: str, end: str) -> str:
    assert start in source, f"missing start marker: {start}"
    tail = source.split(start, 1)[1]
    assert end in tail, f"missing end marker: {end}"
    return start + tail.split(end, 1)[0]


def field_depths(result: dict[int, float], arity: int, threshold: float = 1e-9) -> set[int]:
    """Convert derivative-order keys to spectral depths in the field sector."""
    depths = set()
    for order, coeff in result.items():
        if order < 0 or abs(coeff) <= threshold:
            continue
        depths.add(arity - 1 - order)
    return depths


def test_m6_numerical_stasheff_depth_spectrum_has_no_depth_0_or_1():
    samples = (
        (1.0, 2.0, -1.0, 0.5, 3.0),
        (2.0, -3.0, 5.0, 7.0, -11.0),
        (0.25, 1.5, -2.5, 3.5, -4.5),
    )

    for lams in samples:
        result = StasheffEngine(c_val=12.0).mk(lams)
        assert field_depths(result, arity=6) == {2, 3, 4, 5}
        assert 5 not in result or abs(result.get(5, 0.0)) < 1e-9
        assert 4 not in result or abs(result.get(4, 0.0)) < 1e-9
        assert abs(result[-1]) > 1e-9


def test_gap_migration_theorem_and_m6_table_are_scoped_and_correct():
    source = GRAVITY.read_text()
    section = between(
        source,
        r"\begin{theorem}[Virasoro weight-depth gap",
        r"\begin{theorem}[Symmetric-point period-$2$ parity",
    )
    flat = " ".join(section.split())

    required = (
        r"\ClaimStatusConditional",
        r"$\hypAmbientWtCpl+\effKoszul$",
        r"\operatorname{Spec}(m_6|_T)=\{2,3,4,5\}",
        r"\operatorname{Spec}(m_6)",
        r"\{2,\, 3,\, 4,\, 5,\, 7\}",
        "secondary cancellation at depths $d=0,1$",
        r"$0$ & $5$ & $\partial^5 T$ & $0$ \textup{(secondary cancellation)}",
        r"$1$ & $4$ & $\partial^4 T$ & $0$ \textup{(secondary cancellation)}",
        r"$6$ & $-1$ & \textup{--} & $0$ \textup{(structural gap)}",
        "compute/m7_m10_depth_frontier.py",
        "compute/tests/test_gap_migration_m6.py",
        r"T-propagated binding sector",
        r"d_{\mathrm{gap}}^{\mathrm{bind}}(\mathcal{W}_N,n)",
        "not a statement about the union over all sectors",
    )
    for needle in required:
        assert needle in flat

    retired = (
        r"\begin{theorem}[Gap migration]",
        r"\{1,\, 2,\, 3,\, 4,\, 5,\, 7\}",
        r"$1$ & $4$ & $\partial^4 T$ & $5$",
        r"$6$ & $\{1, 2, 3, 4, 5\}$",
        "The monomial counts in the scalar sector",
        r"d_{\mathrm{gap}}(\mathcal{W}_N, n) = 2N + n - 4",
        r"binding sector \((W_N,\ldots,W_N)\)",
    )
    for needle in retired:
        assert needle not in section


def test_gap_migration_propagation_surfaces_use_binding_sector_scope():
    surfaces = {
        INTRO: ("T-propagated binding-sector gap migrates",),
        W3: (
            "T-propagated binding-sector form",
            "conditional",
            "principal T-propagated binding-sector statement",
        ),
        W_STABLE: (
            "T-propagated binding-sector gap migration",
            r"d_{\mathrm{gap}}^{\mathrm{bind}}(\mathcal{W}_N,\, n)",
        ),
        FRONTIER: (
            "T-propagated binding-sector formula",
            r"not for the raw",
        ),
        ROOT / "chapters" / "examples" / "rosetta_stone.tex": (
            "T-propagated",
            r"d_{\mathrm{gap}}^{\mathrm{bind}}=2N+n-4",
        ),
    }
    for path, needles in surfaces.items():
        text = path.read_text()
        for needle in needles:
            assert needle in text, f"{needle!r} missing from {path}"

    stale_phrases = (
        "was proved in general in",
        r"formula $d_{\mathrm{gap}} = 2N{-}2$ is proved for all principal",
        r"The gap migration formula remains valid as a\nsector-level",
        r"The Virasoro depth column $\{0,\ldots,n{-}1,n{+}1\}$",
    )
    for path in surfaces:
        text = path.read_text()
        for phrase in stale_phrases:
            assert phrase not in text, f"stale phrase survived in {path}"
