"""Source guard for the CFG ordinary-CS lane versus Vol II hCS lane."""
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CFG_SIDE_BY_SIDE = ROOT / "appendices/cfg_side_by_side.tex"
SPECTRAL_BRAIDING_CORE = ROOT / "chapters/connections/spectral-braiding-core.tex"


def test_cfg_ordinary_cs_is_not_identified_with_hcs_chiral_avatar():
    source = " ".join(CFG_SIDE_BY_SIDE.read_text(encoding="utf-8").split())

    assert "ordinary topological Chern--Simons factorisation lane" in source
    assert "Vol~II hCS and chiral-avatar lanes supply Dolbeault BV fields" in source
    assert "CFG2026 remains the ordinary Chern--Simons source lane" in source
    assert "not a CFG2026 theorem" in source
    assert "comparison functors have been applied" in source
    assert "comparison surface, not a dictionary theorem" in source


def test_cfg_hcs_global_comparison_is_componentwise_datum():
    source = CFG_SIDE_BY_SIDE.read_text(encoding="utf-8")
    flat = " ".join(source.split())

    assert r"\label{def:cfg-hcs-comparison-datum}" in source
    assert r"\ClaimStatusDefinitional" in source
    assert r"\mathfrak C_{\mathrm{CFG}\to\hCS}" in source
    for component in (
        r"\mathcal R_{\mathrm{ell}}",
        r"\rho_{\mathrm{BV}}",
        r"\rho_{\mathrm{fac}}",
        r"\tau_{\mathrm{tr}}",
        r"\{h_n,c_n\}_{n\ge 1}",
        r"\nu_{\mathrm{Mell}}",
        r"\eta_{\mathrm{K3}}",
    ):
        assert component in source

    assert r"\mathfrak o_{\mathrm{CFG}\to\hCS}" in source
    for component in (
        "o_{\\mathrm{ell}}",
        "o_{\\mathrm{BV}}",
        "o_{\\mathrm{fac}}",
        "o_{\\mathrm{tr}}",
        "o_{\\mathrm{reg}}",
        "o_{\\mathrm{Mell}}",
        "o_{\\mathrm{K3}}",
    ):
        assert component in source

    assert "A global CFG-to-hCS comparison theorem is exactly the assertion" in flat
    assert "the rows below do not prove that theorem" in flat.lower()
    assert "row-wise shadows and the components of the datum" in flat
    assert "there is no canonical map between these complexes before this replacement is supplied" in flat


def test_spectral_cfg_comparison_is_trace_shadow_not_global_functor():
    source = " ".join(SPECTRAL_BRAIDING_CORE.read_text(encoding="utf-8").split())

    assert "not a theorem identifying the CFG $E_3$-algebra with the hCS" in source
    assert "After choosing the trace component $\\tau_{\\mathrm{tr}}$" in source
    assert "Definition~\\ref{def:cfg-hcs-comparison-datum}" in source
    assert "common perturbative trace shadow" in source
    assert "requiring the trace component $\\tau_{\\mathrm{tr}}$" in source
    assert "between traces, not between global $E_3$-algebras" in source
    assert "does not produce a global comparison functor" in source

    assert "structure arises non-perturbatively from a different route" not in source
    assert "is the perturbative shadow of the $E_3$-topological structure" not in source
    assert "both of which model the same underlying Chern--Simons factorisation algebra" not in source
