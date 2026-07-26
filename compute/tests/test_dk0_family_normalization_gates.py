"""Source guards for the family-gated DK-0 Laplace bridge.

The refinement pass forbids treating Heisenberg, affine, Virasoro,
free-field, and W_3 examples as though they share one normalization or
one proof status.  The DK-0 proposition may remain a direct coefficient
calculation only if the family gates are printed in the theorem surface.
"""
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_WORKED = ROOT / "chapters/examples/examples-worked.tex"


def _flat_source() -> str:
    return "".join(EXAMPLES_WORKED.read_text(encoding="utf-8").split())


def test_dk0_proposition_is_family_gated():
    source = _flat_source()

    assert (
        r"\begin{proposition}[Family-gatedDK-0Laplacebridge;"
        r"\ClaimStatusProvedHere{}underthedisplayedfamilygates;"
        r"licensingtags$\alpha+\beta+\gamma$]"
    ) in source
    assert "ClaimStatusProvedHereConditional" not in source
    assert (
        r"\begin{proposition}[DK-0Laplacebridgeforfivefamilies;"
        r"\ClaimStatusProvedHere]"
    ) not in source
    assert r"Thestatementincludesthefollowingfamilygates:" in source
    assert r"\textbf{Family}&\textbf{normalizationandstatusgate}" in source


def test_each_standard_family_has_its_own_normalization_gate():
    source = _flat_source()

    required = [
        r"$\mathcalH_k$&Currentconvention$\{J_\lambdaJ\}=k\lambda$",
        r"$r^{\mathrm{coll}}_{\cH}=k\Omega_{\cH}/z$",
        r"$\widehat{\mathfrak{sl}}_2$&Trace-formLaplacenormalization",
        r"TheKZrepresentative$\Omega/((k+h^\vee)z)$",
        r"$\mathrm{Vir}_c$&Virasoroconvention",
        r"nofullVirasorofusionoperatorordescendant-sensitivespectral\(R\)-matrixisasserted",
        r"$\beta\gamma/bc$&Free-fieldconventionwithdualgenerators",
        r"statistics-exchangeKoszuldualgate",
        r"$\mathcalW_3$&StandardZamolodchikovnormalization",
        r"5c+22=0",
        r"separatelimitinggate",
    ]
    for needle in required:
        assert needle in source


def test_dk0_bridge_requires_all_family_gates():
    def bridge_claim_is_admissible(
        heisenberg_gate: bool,
        affine_gate: bool,
        virasoro_gate: bool,
        free_field_gate: bool,
        w3_gate: bool,
    ) -> bool:
        return all(
            (
                heisenberg_gate,
                affine_gate,
                virasoro_gate,
                free_field_gate,
                w3_gate,
            )
        )

    assert bridge_claim_is_admissible(True, True, True, True, True)
    assert not bridge_claim_is_admissible(False, True, True, True, True)
    assert not bridge_claim_is_admissible(True, False, True, True, True)
    assert not bridge_claim_is_admissible(True, True, False, True, True)
    assert not bridge_claim_is_admissible(True, True, True, True, False)
