"""Source guards for the 3d HT bridge BV datum.

These tests keep the physics bridge from drifting back to an unnamed
field/action package.  The manuscript theorem and construction must display
the BV tuple, the product HT field complex, the kernel homotopy identity for
the propagator, and the scale quantum master equation.
"""
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAVIOLO = ROOT / "chapters/theory/raviolo.tex"
BV_CONSTRUCTION = ROOT / "chapters/theory/bv-construction.tex"
BV_HT_PHYSICS = ROOT / "chapters/connections/bv_ht_physics.tex"
HT_PHYSICAL_ORIGINS = ROOT / "chapters/connections/ht_physical_origins.tex"
RAVIOLO_RESTRICTION = ROOT / "chapters/theory/raviolo-restriction.tex"
SIX_D_HCS_AVATAR = (
    ROOT / "chapters/connections/six_d_hcs_e3_chiral_avatar_platonic.tex"
)


def _squash(path: Path) -> str:
    return "".join(path.read_text(encoding="utf-8").split())


def _flat(path: Path) -> str:
    return " ".join(path.read_text(encoding="utf-8").split())


def _window(source: str, needle: str, radius: int = 2500) -> str:
    index = source.index(needle)
    return source[max(0, index - radius): index + radius]


def test_physics_bridge_theorem_starts_from_bv_datum():
    source = _squash(RAVIOLO)

    assert (
        r"\begin{theorem}[Bridgefromphysicstoalgebra;"
        r"\ClaimStatusConditional;licensingtags$\alpha+\gamma+\delta$]"
    ) in source
    assert r"(\mathcalE,Q,\omega_{\mathrm{BV}},I,B)" in source
    assert (
        r"\mathcalE\simeq\Omega^\bullet(\R_t)\widehat\otimes"
        r"\Omega^{0,\bullet}(\C_z)\otimes\mathfraka[1]"
    ) in source
    assert r"Q=d_t+\dbar_z+d_{\mathfraka}" in source
    assert r"S_{\mathrm{cl}}=\frac12\omega_{\mathrm{BV}}(A,QA)+I(A)" in source
    assert r"\{S_{\mathrm{cl}},S_{\mathrm{cl}}\}_{\mathrm{BV}}=0" in source


def test_physics_bridge_theorem_requires_kernel_homotopy():
    source = _squash(RAVIOLO)

    assert r"P_{\epsilon<L}=\int_\epsilon^LQ^*K_s\,ds" in source
    assert r"Q_{\boxtimes}P_{\epsilon<L}=K_L-K_\epsilon" in source
    assert r"where$Q_{\boxtimes}$isthedifferentialinducedby$Q$onkernels" in source
    assert r"P_{\mathrm{sing}}(z,t)=\sum_{\alpha=1}^{N}" in source
    assert r"moduloasmooth$Q_{\boxtimes}$-exactkernel" in source
    assert r"QG=\delta_{\C}\boxtimes\delta_{\R}" not in source


def test_bv_construction_displays_scale_qme_and_kernel_identity():
    source = _squash(BV_CONSTRUCTION)

    assert r"(\mathcalE,Q,\omega_{\mathrm{BV}},I,B)" in source
    assert r"\label{eq:3dht-fields}" in source
    assert r"\label{eq:3dht-bv-pairing}" in source
    assert r"\label{eq:3dht-interaction}" in source
    assert r"\label{eq:3dht-heat-propagator}" in source
    assert r"Q_{\boxtimes}P_{\epsilon<L}=K_L-K_\epsilon" in source
    assert r"\label{eq:3dht-propagator-kernel-homotopy}" in source
    assert r"[Q,P_{\epsilon<L}]=K_L-K_\epsilon" in source
    assert r"QI[L]+\hbar\Delta_LI[L]+\frac12\{I[L],I[L]\}_L=0" in source
    assert r"Q_{\mathrm{BRST}}^L=Q+\{I[L],-\}_L+\hbar\Delta_L" in source


def test_bv_bar_boundary_comparison_is_transferred_not_bare_equality():
    """BV graph residues match the bar coderivation only through transfer data."""
    flat = _flat(BV_CONSTRUCTION)

    assert r"\label{eq:3dht-bv-bar-transfer-comparison}" in flat
    assert "boundary strong deformation retract" in flat
    assert r"i_\partial" in flat
    assert r"p_\partial" in flat
    assert r"h_\partial Q_{\mathrm{BRST}}^L+Q_{\mathrm{BRST}}^Lh_\partial" in flat
    assert r"D_A^{\mathrm{BV}}" in flat
    assert "strict equality is asserted only in the transferred boundary model" in flat
    assert "raw graph residues are not identified with the bar coderivation" in flat
    assert r"m_k^{\mathrm{BV}}=m_k^{\mathrm{bar}}\) under the bridge hypotheses" not in flat
    assert "equality of the operations" not in flat


def test_bv_construction_records_standard_family_realisation_rows():
    """The standard-family BV datum must include fields, boundary, QME, maps."""
    source = _flat(BV_CONSTRUCTION)

    assert r"\label{prop:standard-family-bv-realisation-data}" in source
    assert "field complex" in source
    assert "BRST operator" in source
    assert "boundary Lagrangian" in source
    assert r"\eta_A^\partial" in source
    assert r"\chi_A" in source
    assert "Heisenberg/lattice" in source
    assert r"affine \(V_k(\fg)\)" in source
    assert r"\(\beta\gamma\)/symplectic boson" in source
    assert r"\(\mathrm{Vir}_c,\mathcal W_N\)" in source
    assert "principal DS BRST reduction of the affine row" in source
    assert "the algebra \\(A\\) alone never constructs its BV theory" in source


def test_bridge_applies_only_with_full_bv_kernel_qme_package():
    def bridge_applies(
        bv_tuple: bool,
        kernel_homotopy: bool,
        factorized_singular_part: bool,
        scale_qme: bool,
        one_loop_finite: bool,
    ) -> bool:
        return all(
            (
                bv_tuple,
                kernel_homotopy,
                factorized_singular_part,
                scale_qme,
                one_loop_finite,
            )
        )

    assert bridge_applies(True, True, True, True, True)
    assert not bridge_applies(False, True, True, True, True)
    assert not bridge_applies(True, False, True, True, True)
    assert not bridge_applies(True, True, True, False, True)


def test_bv_ht_physics_theorem_heads_are_licensed_in_split_and_active_copy():
    """The BV/HT theorem surface carries the hypothesis package explicitly."""
    for path in (BV_HT_PHYSICS, HT_PHYSICAL_ORIGINS):
        source = _flat(path)
        assert "HT theory from 4d" in source
        assert r"\ClaimStatusConditional; licensing $\alpha+\gamma$" in source
        assert r"\hypKMHTBV" in source
        assert (
            "Costello holomorphic-twist datum, BV gauge fixing, "
            "dimensional reduction"
        ) in source

        boundary = _window(
            source,
            r"\ClaimStatusConditional; licensing $\alpha+\beta+\gamma$ "
            "via a chosen HCS boundary condition",
        )
        assert r"\ClaimStatusConditional; licensing $\alpha+\beta+\gamma$" in boundary
        assert "bulk-to-boundary OPE comparison" in boundary
        assert "level normalisation" in boundary
        assert "whose boundary factorisation algebra" in boundary

        central = _window(
            source,
            r"\ClaimStatusProvedElsewhere; licensing $\alpha+\beta$",
        )
        assert r"\ClaimStatusProvedElsewhere; licensing $\alpha+\beta$" in central
        assert r"\hypDSBRST" in central
        assert "Drinfeld--Sokolov BRST reduction" in central

        observables = _window(
            source,
            r"\ClaimStatusConditional; licensing $\beta+\gamma$ via",
        )
        assert r"\ClaimStatusConditional; licensing $\beta+\gamma$" in observables
        assert "BRST complex" in observables
        assert "BV/bar comparison" in observables


def test_early_ht_localization_is_conditional_bv_statement():
    """The first localization theorem must not read as a bare path integral."""
    source = _flat(HT_PHYSICAL_ORIGINS)
    theorem = _window(source, r"\label{thm:ht-localization}", radius=1800)

    assert "Perturbative localization to holomorphic data" in theorem
    assert r"\ClaimStatusConditional; licensing $\alpha+\gamma$" in theorem
    assert "Costello--Li holomorphic-twist datum" in theorem
    assert "BV gauge fixing" in theorem
    assert "renormalisation scheme" in theorem
    assert "it is not an ordinary convergent measure-theoretic path integral" in theorem
    assert r"\begin{proof}[Sketch]" not in theorem
    assert "Since $\\{Q, V\\} \\geq 0$, the path integral localizes to" not in theorem


def test_early_ht_boundary_chiral_theorem_uses_supplied_boundary_datum():
    """The early CDG boundary theorem must match the later BV hypothesis package."""
    source = _flat(HT_PHYSICAL_ORIGINS)
    theorem = _window(source, r"\label{thm:ht-boundary-chiral-algebra}", radius=2800)

    assert "Boundary chiral algebra from supplied HT boundary data" in theorem
    assert r"\ClaimStatusConditional; licensing $\alpha+\beta+\gamma$" in theorem
    assert "chosen boundary chart" in theorem
    assert "renormalised boundary factorization algebra" in theorem
    assert "bulk-to-boundary OPE comparison" in theorem
    assert "level normalisation" in theorem
    assert "Costello--Dimofte--Gaiotto comparison" in theorem
    assert r"\Obs^\partial_B|_{C_D}" in theorem
    assert r"z \in C_D" in theorem
    assert "unrenormalised analytic HCS path integral is not being invoked" in theorem

    assert "An HT boundary condition supports a chiral algebra" not in theorem
    assert "Boundary operators $\\mathcal{O}(z) = \\lim" not in theorem
    assert "Determined by bulk path integral with boundary insertions" not in theorem
    assert r"\begin{proof}[Sketch]" not in theorem


def test_raviolo_higgs_coulomb_claim_requires_comparison_datum():
    """Higgs/Coulomb quantization is not automatic from the HT bridge alone."""
    source = _flat(RAVIOLO)
    datum = _window(
        source,
        r"\label{def:higgs-coulomb-raviolo-comparison-datum}",
        radius=3200,
    )
    theorem = _window(source, r"\label{prop:Higgs-Coulomb-compat}", radius=3800)

    assert "Higgs--Coulomb raviolo comparison datum" in datum
    assert r"\mathfrak h_{\mathcal T}^{\mathrm{rav}}" in datum
    assert r"F_H,F_C,\Psi_H,\Psi_C,\psi_{\mathrm{cur}},k_{\mathrm{cur}}" in datum
    assert "For a general \\(3\\)d~\\(\\mathcal N=4\\) theory its existence is not asserted here" in datum

    assert "Higgs/Coulomb compatibility under supplied comparison datum" in theorem
    assert r"\ClaimStatusConditional; licensing tags $\alpha+\beta+\gamma$" in theorem
    assert "associated gradeds are identified by \\(\\Psi_H\\) and \\(\\Psi_C\\)" in theorem
    assert "current component \\((\\psi_{\\mathrm{cur}},k_{\\mathrm{cur}})\\)" in theorem
    assert "proves no existence theorem for" in theorem
    assert r"\begin{proof}[Sketch]" not in theorem

    stale = (
        "is simultaneously a quantization of the Higgs and Coulomb branch Poisson algebras",
        "expected $\\widehat{\\mathfrak g}$",
        "Direct calculation in three families verifies this:",
    )
    for needle in stale:
        assert needle not in theorem


def test_raviolo_restriction_echo_keeps_higgs_coulomb_comparison_conditional():
    """The expanded raviolo corollary must not re-export the old global claim."""
    source = _flat(RAVIOLO_RESTRICTION)
    corollary = _window(source, r"\label{cor:agreement-HT}", radius=2600)

    assert "Agreement with HT constructions under comparison data" in corollary
    assert r"\ClaimStatusConditional; licensing tags $\alpha+\beta+\gamma$" in corollary
    assert "secondary-product comparison datum" in corollary
    assert r"\mathfrak h_{\mathcal T}^{\mathrm{rav}}" in corollary
    assert "factorization-homology comparison along \\(\\R\\)" in corollary
    assert (
        "For HT twists of 3d $\\mathcal N=2,4$ theories whose observables "
        "lie in the scope"
        not in corollary
    )


def test_bv_ht_split_e3_surface_requires_dolbeault_koszul_reduction():
    """The split echo must not advertise raw C^3 hCS as E3 before reduction."""
    source = _flat(BV_HT_PHYSICS)

    e3 = _window(source, r"\label{thm:E3-action-6d-hCS}")
    assert r"\ClaimStatusConditional; licensing $\gamma$" in e3
    assert "Dolbeault--Koszul reduction" in e3
    assert r"\hypKMHTBV" in e3
    assert "not asserted on the unreduced six-real-dimensional Dolbeault complex" in e3

    dunn = _window(source, r"\label{thm:dunn-E3-6d-hCS}")
    assert r"\ClaimStatusProvedElsewhere; licensing $\gamma$" in dunn
    assert "After the Dolbeault--Koszul reduction" in dunn
    assert "This theorem counts real locally-constant directions" in dunn

    p3 = _window(source, r"\label{thm:P3-obs-6d-hCS}")
    assert r"\ClaimStatusConditional; licensing $\gamma$" in p3
    assert "Dolbeault-reduced classical observables" in p3
    assert "degree \\(1-n\\)" in p3


def test_six_d_hcs_avatar_e3_surface_is_dolbeault_reduced():
    """The live 6d hCS chapter states E3 only after Dolbeault reduction."""
    source = _flat(SIX_D_HCS_AVATAR)

    prologue = _window(source, "On $X = \\C^3$")
    assert "Dolbeault--Koszul cohomological reduction" in prologue
    assert "three residual locally constant $E_1$ directions" in prologue
    assert "unreduced six-real-dimensional Dolbeault complex" in prologue

    theorem = _window(source, r"\label{thm:hcs-e3-from-translation}")
    assert r"\ClaimStatusConditional; licensing $\gamma$" in theorem
    assert "Dolbeault--Koszul reduction" in theorem
    assert r"\hypKMHTBV" in theorem
    assert "translation-equivariant sector" in theorem
    assert "not an $E_3$ assertion on the unreduced six-real-dimensional" in theorem

    assert "E_3$ after three Dolbeault--Koszul contractions" in source
