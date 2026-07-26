"""Source guards for the higher-genus Hochschild bridge.

The genus-g bridge is a filtered curved comparison.  It is not an
ordinary quasi-isomorphism of curved complexes unless the scalar curvature
vanishes.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HOCHSCHILD = ROOT / "chapters" / "connections" / "hochschild.tex"
THOLOG = (
    ROOT
    / "chapters"
    / "connections"
    / "anomaly_completed_topological_holography.tex"
)


def _source(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return " ".join(text.split())


def _higher_genus_block() -> str:
    text = _source(HOCHSCHILD)
    start = text.index(r"\label{thm:hochschild-bridge-higher-genus}")
    end = text.index(r"\begin{proposition}[Moduli action", start)
    return text[start:end]


def _full_bridge_remark() -> str:
    text = _source(HOCHSCHILD)
    start = text.index(r"\begin{remark}[The full Hochschild bridge]")
    end = text.index(r"\begin{corollary}[Bridge to Volume~I Theorem~H", start)
    return text[start:end]


def test_higher_genus_bridge_defines_curved_comparison():
    text = _source(HOCHSCHILD)
    block = _higher_genus_block()

    assert r"\label{def:filtered-curved-comparison}" in text
    assert "filtered curved comparison" in block
    assert r"\Psi_g\colon" in block
    assert r"\ChirHoch^\bullet_g" in block
    assert r"\Theta_g\in H^2_{\mathrm{cyc}}(\cA,\cA)" in block
    assert (
        "not an ordinary quasi-isomorphism unless the curvature vanishes"
        in _flat(block)
    )


def test_higher_genus_bridge_rejects_ordinary_curved_quasi_isomorphism():
    text = _source(HOCHSCHILD)
    block = _higher_genus_block()
    remark = _full_bridge_remark()

    assert "canonical filtered quasi-isomorphism of complexes" not in block
    assert "quasi-isomorphism of complexes" not in block
    assert r"H^2_{\mathrm{cyc}}(\ChirHoch^\bullet(\cA))" not in text
    assert r"biholomorphic to $\C$ (uniformization)" not in text
    assert "ARE the chiral Hochschild cochains" not in text
    assert "With this bridge closed" not in text
    flat_remark = _flat(remark)
    assert "filtered curved \\(R_g\\)-sense" in flat_remark
    assert "only when the curvature vanishes" in flat_remark


def test_downstream_anomaly_heuristic_uses_curved_comparison():
    text = _source(THOLOG)
    start = text.index(r"\begin{remark}[Heuristic justification]")
    end = text.index("Centrality:", start)
    heuristic = text[start:end]

    assert "filtered curved comparison" in heuristic
    assert "not an ordinary quasi-isomorphism at nonzero curvature" in heuristic
    assert "provides a\nfiltered quasi-isomorphism" not in heuristic
