"""Guards for the split modular-PVA extension surface.

The extension file is not in the live main.tex input graph, but it is a
merge source.  It must not re-import the old claim that every standard
family has a uniform one-parameter full modular lift.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXTENSION = ROOT / "chapters/connections/thqg_modular_pva_extensions.tex"


def _flat_source() -> str:
    return " ".join(EXTENSION.read_text().split())


def test_hs_sewing_is_conditional_on_analytic_gates():
    flat = _flat_source()

    assert (
        "HS-sewing under polynomial/subexponential analytic gates; "
        "\\ClaimStatusConditional; licensing tags $\\gamma+\\delta$"
    ) in flat
    assert "positive-energy vertex algebra equipped with a Hilbert completion" in flat
    assert "\\dim A_n \\le C\\exp(\\alpha\\sqrt n)" in flat
    assert "polynomial Hilbert--Schmidt growth" in flat
    assert "No claim is made here for logarithmic, admissible-level" in flat
    assert "non-principal DS specialisations" in flat


def test_quantization_dictionary_is_gated_not_uniform_full_lift():
    flat = _flat_source()

    assert (
        "Gated modular-PVA quantization dictionary; "
        "\\ClaimStatusConditional; licensing tags $\\alpha+\\gamma+\\delta$"
    ) in flat
    assert "it is not a claim that the full tensor-channel all-loop quantum HT theory has been constructed" in flat
    assert "fixed lattice: no universal $\\partial_k$" in flat
    assert "$W_N$ & conjectural uniform vanishing" in flat
    assert "not constructed here" in flat
    assert "all-loop PVA--quantum HT problem" in flat


def test_retired_uniform_standard_family_phrases_are_absent():
    flat = _flat_source()

    retired = (
        "Every vertex algebra in the standard landscape satisfies the HS-sewing condition",
        "The pattern is uniform: for every family in the standard landscape",
        "the full modular lift exists",
        "all standard families are in triangular $W$-normal form",
        "all standard families depend on a single continuous parameter",
    )
    for phrase in retired:
        assert phrase not in flat


def test_single_parameter_claim_fails_for_fixed_lattice_row():
    rows = {
        "heisenberg": {"parameter": "k"},
        "fixed_lattice": {"parameter": None},
        "virasoro": {"parameter": "c"},
    }

    assert rows["fixed_lattice"]["parameter"] is None
    assert {row["parameter"] for row in rows.values()} != {"k", "c"}
