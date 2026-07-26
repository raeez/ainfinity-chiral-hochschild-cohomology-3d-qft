from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TARGET = "chapters/theory/foundations.tex"


def read_squashed() -> str:
    return " ".join((ROOT / TARGET).read_text().split())


def test_foundations_theorem_d_summary_splits_scalar_and_all_weight_free_energy():
    squashed = read_squashed()

    required_forms = (
        (
            "F_g^{\\mathrm{sc}}(\\cA) "
            "=\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}"
        ),
        (
            "F_g(\\cA)=F_g^{\\mathrm{sc}}(\\cA) "
            "+\\delta F_g^{\\mathrm{cross}}(\\cA)"
        ),
        "F_1^{\\mathrm{sc}}(\\cA)=\\kappaChHodge(\\cA)/24",
        "with $F_g(\\cA)=F_g^{\\mathrm{sc}}(\\cA) +\\delta F_g^{\\mathrm{cross}}(\\cA)$ in all weights",
    )
    for required in required_forms:
        assert required in squashed


def test_foundations_theorem_d_summary_has_no_bare_fp_free_energy_slogan():
    squashed = read_squashed()

    stale_forms = (
        "F_g = \\kappaChHodge(\\cA) \\cdot \\lambda_g^{\\mathrm{FP}}",
        "F_g=\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}",
        "F_g(\\cA)=\\kappaChHodge(\\cA)\\lambda_g^{\\mathrm{FP}}",
        "F_g(\\cA)=\\kappaChHodge(\\cA) \\cdot \\lambda_g^{\\mathrm{FP}}",
        "F_1 = \\kappaChHodge(\\cA)/24",
    )
    for stale in stale_forms:
        assert stale not in squashed
