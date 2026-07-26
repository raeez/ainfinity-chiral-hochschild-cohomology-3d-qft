from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TARGET = "chapters/connections/bv_brst.tex"


def read_squashed() -> str:
    return " ".join((ROOT / TARGET).read_text().split())


def test_bv_brst_theorem_d_summary_uses_scalar_free_energy():
    squashed = read_squashed()

    required_forms = (
        (
            "F_g^{\\mathrm{sc}}(\\cA) "
            "=\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}"
        ),
        "multi-weight correction $\\delta F_g^{\\mathrm{cross}}$",
        (
            "class-$\\mathsf{C}$ scalar pattern "
            "\\(F_g^{\\mathrm{sc}}(\\cA) "
            "=\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}\\)"
        ),
    )
    for required in required_forms:
        assert required in squashed


def test_bv_brst_theorem_d_summary_has_no_bare_fp_free_energy_slogan():
    squashed = read_squashed()

    stale_forms = (
        "F_g = \\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}",
        "F_g = \\kappaChHodge(\\cA)\\cdot \\lambda_g^{\\mathrm{FP}}",
        "F_g(\\cA) = \\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}",
        "F_g(\\cA)=\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}",
        "F_g(\\cA)=\\kappaChHodge(\\cA)\\cdot\\lambda_g^{\\mathrm{FP}}",
        (
            "pattern \\(F_g = \\kappaChHodge(\\cA)\\cdot "
            "\\lambda_g^{\\mathrm{FP}}\\)"
        ),
    )
    for stale in stale_forms:
        assert stale not in squashed
