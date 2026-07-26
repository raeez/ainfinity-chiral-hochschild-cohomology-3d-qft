from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SURFACES = (
    ROOT / "chapters" / "connections",
    ROOT / "chapters" / "theory",
    ROOT / "chapters" / "examples",
)


def iter_tex_sources():
    for surface in SURFACES:
        yield from surface.rglob("*.tex")


def squashed(path: Path) -> str:
    return " ".join(path.read_text().split())


def test_active_surfaces_do_not_identify_full_free_energy_with_scalar_trace():
    stale_forms = (
        "F_g(\\cA)=\\kappaChHodge",
        "F_g(\\cA) = \\kappaChHodge",
        "F_g(\\cA_T)=\\kappaChHodge",
        "F_g(\\cA_T) = \\kappaChHodge",
        "F_g=\\kappaChHodge",
        "F_g = \\kappaChHodge",
        "F_g(\\widehat{\\mathfrak{g}}_k) = \\kappaChHodge",
        "scalar formula receives",
        "genus-$g$ free energy acquires a cross-channel correction",
    )

    offenders: list[str] = []
    for path in iter_tex_sources():
        text = path.read_text()
        for stale in stale_forms:
            if stale in text:
                rel = path.relative_to(ROOT)
                offenders.append(f"{rel}: contains {stale!r}")

    assert offenders == []


def test_repaired_surfaces_state_scalar_trace_and_full_decomposition():
    anchors = {
        "chapters/connections/holomorphic_topological.tex": (
            (
                "F_g^{\\mathrm{sc}}(\\cA_T) "
                "=\\kappaChHodge(\\cA_T)\\lambda_g^{\\mathrm{FP}}"
            ),
            (
                "F_g(\\cA_T)=F_g^{\\mathrm{sc}}(\\cA_T) "
                "+\\delta F_g^{\\mathrm{cross}}(\\cA_T)"
            ),
        ),
        "chapters/connections/thqg_perturbative_finiteness.tex": (
            (
                "F_g^{\\mathrm{sc}}(\\cA) \\;=\\; "
                "\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}"
            ),
            (
                "F_g(\\cA)=F_g^{\\mathrm{sc}}(\\cA) "
                "+\\delta F_g^{\\mathrm{cross}}(\\cA)"
            ),
        ),
        "chapters/connections/modular_pva_quantization_core.tex": (
            (
                "F_g^{\\mathrm{sc}}(\\cA)=\\kappaChHodge(\\cA)"
                "\\lambda_g^{\\mathrm{FP}}"
            ),
            (
                "F_g(\\cA)=F_g^{\\mathrm{sc}}(\\cA)"
                "+\\delta F_g^{\\mathrm{cross}}(\\cA)"
            ),
        ),
    }

    for rel_path, required_forms in anchors.items():
        text = squashed(ROOT / rel_path)
        for required in required_forms:
            assert required in text
