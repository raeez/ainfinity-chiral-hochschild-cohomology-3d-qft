from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TARGET = "chapters/connections/holomorphic_topological.tex"


def read_squashed() -> str:
    return " ".join((ROOT / TARGET).read_text().split())


def test_ht_mc_theorem_splits_scalar_trace_from_all_weight_free_energy():
    squashed = read_squashed()

    required_forms = (
        (
            "F_g^{\\mathrm{sc}}(\\cA_T) "
            "=\\kappaChHodge(\\cA_T)\\lambda_g^{\\mathrm{FP}}"
        ),
        (
            "F_g(\\cA_T)=F_g^{\\mathrm{sc}}(\\cA_T) "
            "+\\delta F_g^{\\mathrm{cross}}(\\cA_T)"
        ),
        "multi-weight cross-channel trace",
    )
    for required in required_forms:
        assert required in squashed


def test_ht_costello_gaiotto_comparisons_use_scalar_uniform_weight_formula():
    squashed = read_squashed()

    required_forms = (
        (
            "F_g^{\\mathrm{sc}}(\\cA)=\\kappaChHodge(\\cA)"
            "\\lambda_g^{\\mathrm{FP}}"
        ),
        (
            "F_g^{\\mathrm{sc}}(\\cA_T)=\\kappaChHodge(\\cA_T)"
            "\\lambda_g^{\\mathrm{FP}}"
        ),
        (
            "F_g^{\\mathrm{sc}}(\\widehat{\\mathfrak{g}}_k) = "
            "\\kappaChHodge(\\widehat{\\mathfrak{g}}_k) "
            "\\lambda_g^{\\mathrm{FP}}"
        ),
    )
    for required in required_forms:
        assert required in squashed


def test_ht_burns_example_names_uniform_weight_equality_explicitly():
    squashed = read_squashed()

    required_forms = (
        (
            "F_g^{\\mathrm{sc}}(\\cA_{\\mathrm{Burns}}) "
            "=\\kappaChHodge(\\cA_{\\mathrm{Burns}}) "
            "\\lambda_g^{\\mathrm{FP}}"
        ),
        (
            "F_g(\\cA_{\\mathrm{Burns}}) "
            "=F_g^{\\mathrm{sc}}(\\cA_{\\mathrm{Burns}})"
        ),
        (
            "F_2^{\\mathrm{sc}}(\\cA_{\\mathrm{Burns}}) \\;=\\; "
            "F_2(\\cA_{\\mathrm{Burns}})"
        ),
    )
    for required in required_forms:
        assert required in squashed


def test_ht_surface_has_no_bare_fp_free_energy_slogan():
    squashed = read_squashed()

    stale_forms = (
        "F_g=\\kappaChHodge(\\cA_T) \\cdot \\lambda_g^{\\mathrm{FP}}",
        (
            "F_g(\\cA) = \\kappaChHodge(\\cA) "
            "\\cdot \\lambda_g^{\\mathrm{FP}}"
        ),
        (
            "F_g(\\cA_T) = \\kappaChHodge(\\cA_T) "
            "\\cdot \\lambda_g^{\\mathrm{FP}}"
        ),
        (
            "F_g(\\widehat{\\mathfrak{g}}_k) = "
            "\\kappaChHodge(\\widehat{\\mathfrak{g}}_k)"
            "\\cdot \\lambda_g^{\\mathrm{FP}}"
        ),
        (
            "gives $F_g(\\cA_{\\mathrm{Burns}}) = "
            "\\kappaChHodge(\\cA_{\\mathrm{Burns}})"
        ),
    )
    for stale in stale_forms:
        assert stale not in squashed


def test_ht_support_depth_is_not_a_bps_multiplet_counter():
    squashed = read_squashed()

    required_forms = (
        "Support-depth archetype classification",
        (
            "support depth \\(r^\\Theta_{\\max}(\\cA_T)\\) "
            "of the canonical bar-intrinsic MC packet"
        ),
        "only after a BPS comparison datum has been fixed",
        "stability condition, a central-charge map",
        "It does not count protected BPS multiplet types by itself.",
        "not on support depth alone",
        "Support class~$\\mathsf{L}$, \\(r^\\Theta_{\\max}=3\\)",
        "Support class~$\\mathsf{C}$, \\(r^\\Theta_{\\max}=4\\)",
        "Support class~$\\mathsf{M}$, \\(r^\\Theta_{\\max}=\\infty\\)",
        "Boundary VOA classification by support depth",
        "Formality and support depth",
        "support packet detects the \\emph{residual} quartic anomaly",
        "genus-$0$ support-packet projections",
    )
    for required in required_forms:
        assert required in squashed

    stale_forms = (
        "r_{\\max}",
        "r_max",
        "rmax",
        "\\rmax",
        "Shadow archetype",
        "shadow archetype",
        "shadow class",
        "shadow classes",
        "shadow depth",
        "shadow obstruction tower",
        "tower terminates",
        "counts the number of protected BPS multiplet types",
        "BPS interpretation depends",
        "while the BPS interpretation depends on~$r_{\\max}$ alone",
        "protected BPS multiplet types whose contributions survive",
        "Shadow class~$\\mathsf{L}$",
        "Shadow class~$\\mathsf{C}$",
        "Shadow class~$\\mathsf{M}$",
    )
    for stale in stale_forms:
        assert stale not in squashed
