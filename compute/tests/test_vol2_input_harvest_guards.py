"""Source guards for rigorous corrections harvested from the June 2026 PDFs."""
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _read(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def _flat(rel: str) -> str:
    return " ".join(_read(rel).split())


def _assert_all(source: str, needles: tuple[str, ...]) -> None:
    for needle in needles:
        assert needle in source


def test_spectral_discriminants_are_cyclic_modules_and_divisor_cores():
    text = _read("chapters/connections/casimir_divisor_core_transport.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\label{thm:divisor-core-calculus}",
            r"\label{thm:minimal-intrinsic-realization}",
            r"\label{constr:casimir-divisor-shadow}",
            r"H^{1,\alg}_{\red}(A):=\mathcal R_A^\vee",
            r"M(p_A),\qquad p_A(t)=t^{d_A}\Delta_A(t^{-1}).",
        ),
    )
    assert "cyclic transport object" in flat
    assert (
        r"A divisor of the discriminant does \emph{not} in general "
        r"define a projector; it defines a canonical subquotient."
    ) in flat


def test_relative_hochschild_duality_is_cech_duality_with_bridge_data():
    text = _read("chapters/connections/ym_higher_body_couplings.tex")

    _assert_all(
        text,
        (
            r"\label{def:relative-hochschild-cochains}",
            r"\label{thm:relative-duality-theorem}",
            r"\label{def:w-derived-center-cech}",
            r"\label{thm:w-boundary-cech-duality}",
            r"\label{def:higher-body-ym-bridge-datum}",
            "Novikov-completed deformation complex",
            "screening positivity",
            "domination estimates",
        ),
    )


def test_instanton_completion_and_screenings_are_novikov_hodge_gap_data():
    text = _read("chapters/connections/ym_instanton_screening.tex")

    _assert_all(
        text,
        (
            r"\label{thm:novikov-completion-theorem}",
            r"\label{thm:instanton-completed-tangent-center}",
            r"\label{def:central-koszul-screening-complex}",
            r"\label{thm:screening-hodge-theorem}",
            r"\label{cor:screening-spectral-gap-criterion}",
            r"\label{def:internal-screening-gap-domination}",
            r"\label{thm:mass-gap-reduction-internal-screening}",
        ),
    )


def test_celestial_transfer_is_obstruction_tower_and_residue_identity():
    text = _read("chapters/connections/celestial_boundary_transfer_core.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\label{def:cbt-filtered-celestial}",
            r"\label{def:cbt-obstruction-classes}",
            r"\label{thm:cbt-recursive-linearization}",
            r"\label{thm:cbt-spectral-obstruction}",
            r"\label{thm:cbt-lowest-nonlinearity}",
            r"\label{thm:cbt-bf-ladder}",
            r"\label{thm:cbt-airy-witt}",
            r"\label{def:cbt-pt-residue-datum}",
            r"\label{thm:cbt-mhv-bar-integral}",
        ),
    )
    assert "The theorem is a residue identity on the cyclic FM boundary" in flat
    assert "not construct the Lorentzian scattering theory" in flat


def test_chiral_ce_generators_relations_and_6d_hcs_reduction_are_printed():
    text = _read("chapters/connections/chiral_ce_factalg_gen_rel.tex")

    _assert_all(
        text,
        (
            r"\label{def:ch-ce-cochains}",
            r"\label{def:ch-ce-generators}",
            r"\label{def:ch-ce-factalg}",
            r"\label{thm:ch-ce-weiss-descent}",
            r"\label{thm:ch-ce-disk3-algebra}",
            r"\label{thm:ch-ce-fubini-disk3-to-disk1}",
            r"\label{def:sunrise-integrand}",
            r"\label{thm:ch-ce-sunrise-sl2-minus-eighth}",
            r"\label{thm:ch-ce-6d-hcs-reduction}",
        ),
    )


def test_six_dimensional_hcs_coefficients_are_graph_complex_coefficients():
    text = _read("chapters/theory/six_d_hcs_feynman_coefficients.tex")

    _assert_all(
        text,
        (
            r"\label{thm:6d-hcs-one-loop-coefficient}",
            r"\label{thm:6d-hcs-two-loop-coefficient}",
            r"\label{thm:6d-hcs-three-loop-coefficient}",
            r"\label{prop:cfg-6d-comparison}",
            r"\label{thm:four-channel-match}",
            r"\begin{theorem}[One-loop scalar graph weight; \ClaimStatusProvedHere; licensing tags $\alpha+\beta+\gamma$]",
            r"\begin{theorem}[Two-loop scalar graph weight; \ClaimStatusProvedHere; licensing tags $\alpha+\beta+\gamma$]",
            r"\begin{theorem}[Third-obstruction scalar graph weight; \ClaimStatusProvedHere; licensing tags $\alpha+\beta+\gamma$]",
            r"\section{All-loop BV anomaly conjecture}",
            r"\begin{conjecture}[All-loop 6d hCS BV anomaly; \ClaimStatusConjectured; licensing tags $\beta+\gamma+\epsilon$]",
            "No all-loop theorem is claimed here.",
            "is conjectured to compare with the",
            "graph complex",
            "sunrise",
            "Kontsevich",
        ),
    )
    assert r"\section{The theorem to be proved at loop $n$ in general}" not in text
    assert "is expected to compare with the" not in text


def test_drinfeld_center_bulk_comparison_remains_conjectural():
    text = _read("chapters/connections/hochschild.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\label{conj:drinfeld-center-equals-bulk}",
            r"\ClaimStatusConjectured",
            "Drinfeld centre is conjectured in",
            "below conjecturally identifies the center of the chiral double",
            "does not assert that bar--cobar ``produces the bulk.''",
        ),
    )
    assert "Drinfeld centre is expected to compare with the bulk derived centre" not in flat
    assert "The analogous chiral statement identifies the center" not in flat


def test_kontsevich_integral_is_restricted_operation_not_global_feynman_slogan():
    text = _read("chapters/connections/kontsevich_integral.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\label{thm:kontsevich-invariant}",
            r"\label{thm:bar-weight-systems}",
            r"\label{prop:propagator-restriction}",
            r"\label{prop:feynman-graph-complex}",
            r"\label{def:real-loop-obstruction}",
            r"\label{constr:kontsevich-rectification}",
            "standard framing normalization",
            "4T relation",
        ),
    )
    assert "This restriction-to-$S^1$ step is part of the theorem." in flat
    assert "not a consequence of the Feynman-transform isomorphism alone" in flat


def test_non_principal_w_algebras_have_orbit_sensitive_bar_data():
    text = _read("chapters/connections/fractional_ghost_chain_level_platonic.tex")

    _assert_all(
        text,
        (
            r"\label{def:good-1-over-d-grading}",
            r"\label{def:branched-cover-integralization}",
            r"\label{thm:branched-cover-integralization}",
            r"\label{thm:sugawara-antighost-branched-cover}",
            r"\label{thm:E3-topological-DS-general-explicit-BP}",
            r"\label{cor:fractional-ghost-healed-non-principal}",
            r"\label{prop:three-lane-fractional-ghost}",
            "good grading",
            "fractional conformal weight",
            "Galois-descent",
        ),
    )


def test_w3_lambda_channel_and_gaudin_scope_are_harvested():
    text = _read("chapters/examples/w-algebras-w3.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\label{eq:w3-lambda-WW}",
            r"\label{eq:w3-lambda-projector}",
            r"\label{eq:w3-unit-lambda-channel}",
            r"\label{eq:w3-lambda-boundary-kernel}",
            r"\label{eq:w3-m4-lambda-channel}",
            r"\label{eq:w3-m4-literal-lambda-channel}",
            r"\label{eq:w3-primary-lambda-eigenvalue}",
            r"\label{thm:w3-CYBE}",
            r"K^{W_3}(u, v) = \begin{pmatrix}",
            "No non-dynamical CYBE for the",
            "field-dependent matrix",
        ),
    )
    assert (
        r"P_\Lambda\,m_2(W,W;\lambda) = \beta\,\lambda\,\Lambda"
        in flat
    )
    assert (
        r"\Lambda_0|h\rangle = \left(h^2-\frac{3h}{5}\right)|h\rangle"
        in flat
    )
    assert (
        "The commutativity $[H_i^{\\mathrm{Gaud}}, H_j^{\\mathrm{Gaud}}] = 0$ "
        "is a theorem only after a $W_3$ Gaudin algebra with the corresponding "
        "quantum $RLL$ relation has been constructed."
    ) in flat
    assert "does not by itself construct the commuting Hamiltonians" in flat
    assert (
        "The $\\mathcal{W}_3$ Bethe-equation comparison is conjectural here"
        in flat
    )
    assert "requires a supplied $\\mathcal W_3$ Gaudin/RLL algebra" in flat
    assert "The target obligation is that the holonomic $D$-module" in flat
    assert "open recognition diagram is" in flat
    assert "target-recognition data, not a construction claim" in flat
    assert "$\\mathcal{W}_3$ Bethe equations are expected to be" not in flat
    assert "is expected to carry this order-$8$ class" not in flat
    assert "diagram to be proved is" not in flat


def test_w3_m4_odd_sector_is_parity_vanishing_not_a_sketch():
    text = _read("chapters/examples/w-algebras-w3.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\label{thm:m4-WWWW-tree}",
            r"\ClaimStatusProvedHere; licensing tags $\alpha+\gamma$",
            r"\sigma(T)=T",
            r"\sigma(W)=-W",
            r"\label{eq:m4-WWWW-odd-vanishing}",
            r"\label{eq:m4-WWWW-equivariant-odd-zero}",
        ),
    )
    assert (
        r"\Pi_{\cW_3^{\mathrm{odd}}}\Pi_{9}A_4^{WW}(\ell,\mu,\nu)=0"
        in text
    )
    assert (
        "the actual tree representative also depends on the chosen contraction"
        in flat
    )
    assert "The proof is a direct but lengthy computation" not in text
    assert "sketched below" not in text


def test_soft_graviton_phantom_citation_cluster_is_removed():
    active_sources = {
        "main.tex": _read("main.tex"),
        "soft": _read("chapters/connections/soft_graviton_mellin_shadow_bridge_platonic.tex"),
        "universal": _read("chapters/connections/universal_celestial_holography.tex"),
        "core": _read("chapters/connections/celestial_holography_core.tex"),
        "thqg": _read("chapters/connections/thqg_soft_graviton_theorems.tex"),
        "gravity": _read("chapters/connections/3d_gravity.tex"),
        "winfty": _read("chapters/connections/w_infty_e_infty_endpoint_platonic.tex"),
    }
    joined = "\n".join(active_sources.values())

    _assert_all(
        active_sources["main.tex"],
        (
            r"\bibitem{HamadaShiuSubsubSoft}",
            r"\bibitem{GuevaraHimwichPateStrominger21}",
            "arXiv:1801.05528",
            "arXiv:2103.03961",
        ),
    )
    _assert_all(
        active_sources["soft"],
        (
            r"not an external amplitude derivation",
            r"not used as coefficient oracles",
            r"external comparison surface",
            r"physical celestial central charge is a separate comparison datum",
        ),
    )
    _assert_all(
        active_sources["universal"],
            (
                "The leading soft-graviton theorem gives the supertranslation current.",
                r"The celestial \(\wtc\)-algebra is a separate soft-current algebra",
                "coefficient oracle for this rational function",
            ),
        )
    _assert_all(
        active_sources["gravity"],
        (
            "Cachazo--Strominger / Hamada--Shiu sub-subleading",
            r"\cite{CachazoStrominger14,HamadaShiuSubsubSoft}",
        ),
    )

    retired = (
        "LiStromingerHigherSoft",
        "Li--Strominger",
        "Li-Strominger",
        "Hamada--Shiu--Li--Strominger",
        "Cachazo--He--Yuan",
        "He--Mitra--Pasterski--Strominger",
        "arXiv:XXXX.XXXXX",
        r"\cite[(3.11)--(3.12)]{HamadaShiuSubsubSoft}",
        r"\cite[(4.23)]{HamadaShiuSubsubSoft}",
        "leading one-loop order",
        "genuine cross-derivation",
    )
    for needle in retired:
        assert needle not in joined


def test_celestial_ope_transfer_separates_bar_residue_from_physical_ope():
    core = _read("chapters/connections/celestial_holography_core.tex")
    main = _read("main.tex")
    flat = " ".join(core.split())

    _assert_all(
        main,
        (
            r"\bibitem{PateRaclariuStromingerYuan19}",
            "arXiv:1910.07424",
        ),
    )
    _assert_all(
        core,
        (
            "internal collision-residue data",
            r"r^{\mathrm{KM}}_{\mathrm{lev}}(z)=\frac{k\,\kappa^{ab}\,t_a\otimes t_b}{z}",
            "This tensor is symmetric in the two colour labels.",
            "antisymmetric structure constant",
            r"\frac{B(\Delta_1-1,\Delta_2-1)f^{ab}{}_c}{z_1 - z_2}",
            r"\cite{FanPasterskiShaoCelOPE,PateRaclariuStromingerYuan19}",
            r"\frac{\bar z_{12}}{z_{12}}",
            r"B(\Delta_1-1,\Delta_2-1)",
            r"\(\Delta_1+\Delta_2-3\)",
        ),
    )
    assert "These residues are not, by themselves, the physical celestial OPE coefficients." in flat
    assert "the bar level residue and the physical celestial OPE coefficient are distinct pieces" in flat
    assert "No tree-level physical celestial OPE contains" in flat
    assert (
        "the bar complex alone does not extract \\(f^{ab}{}_{c}\\) from \\(k\\Omega\\)"
        in flat
    )

    retired = (
        "colour projection of $k\\Omega$ is $f^{ab}{}_{c}$",
        "bar-complex form of the celestial gluon OPE obtained by direct Mellin transform",
        "The two-pole structure matches the Pasterski--Shao--Strominger graviton OPE.",
        r"\frac{c/2}{(z_1 - z_2)^3}",
        r"\mathcal{O}^+_{\Delta_1 + \Delta_2 - 3}",
    )
    for needle in retired:
        assert needle not in core


def test_celestial_moonshine_bridge_is_virtual_conditional_and_not_global_m24():
    cmb = _read("chapters/connections/celestial_moonshine_bridge.tex")
    main = _read("main.tex")
    endpoint = _read("chapters/connections/w_infty_e_infty_endpoint_platonic.tex")
    gravity = _read("chapters/connections/3d_gravity.tex")
    active_celestial_sources = "\n".join((main, cmb, endpoint, gravity))
    flat = " ".join(cmb.split())

    _assert_all(
        cmb,
        (
            "Gannon's virtual graded module",
            "not a symmetry group of a single K3 sigma model",
            r"\begin{theorem}[Celestial--moonshine bridge; \ClaimStatusConditional]",
            r"\begin{theorem}[Mathieu--celestial correspondence; \ClaimStatusConditional]",
            r"7A/7B",
            r"23A/23B",
            r"\begin{theorem}[Shadow-tower moonshine comparison; \ClaimStatusConditional]",
            r"\(-2\)",
            r"\mathcal R_n^{\mathrm{K3}}",
            "Duncan--Griffin--Ono",
            r"\begin{proposition}[Umbral celestial correspondence; \ClaimStatusConditional]",
            r"\mathfrak U_N",
        ),
    )
    assert "The Frame shape is a useful label for eta-products, but it is not an injective" in flat
    assert (
        r"\begin{corollary}[Celestial placement of Mathieu twining characters; \ClaimStatusConditional]"
        in flat
    )
    assert "not \\(-\\Tr_{\\widetilde H(K3,\\ZZ)}(g)/12\\)" in flat
    assert "there is no theorem identifying bar depth \\(r\\) with a Rademacher" in flat

    retired_main = (
        "DonnayFotopoulosPasterskiTaylor21",
        "DonnayFotopoulosPasterskiTaylor22",
        "Celestial Glauber phase from asymptotic symmetries",
    )
    for needle in retired_main:
        assert needle not in active_celestial_sources

    retired_cmb = (
        "Donnay--Fotopoulos--Pasterski--Taylor",
        "DFPT",
        r"$\MtwFour$-action on the combined",
        "the $24$ K3-realized conjugacy classes",
        "indexed by the $24$ K3-realized Frame shapes",
        r"S_r^{(g)}(\cA^{\mathcal{N}{=}4}_{K3})",
        r"\mathrm{Rad}_r\bigl[h_g(\tau)\bigr]",
        r"$-\Tr_{\widetilde H(K3,\ZZ)}(g) / 12$ for general $g$",
        "specialization of Gannon's theorem to Niemeier lattices",
        "where $G_N$ embeds into the",
        "Compareerence",
        "Moonshine is an algebraic consequence of the celestial framework",
        "Mathieu moonshine is a theorem of the celestial framework",
        "outputs of the celestial dictionary",
        r"\ClaimStatusProvedHere\ \text{in celestial sector",
        "Mathieu moonshine is an output of celestial holography",
        "with moonshine Rademacher coefficients",
    )
    for needle in retired_cmb:
        assert needle not in cmb


def test_class_m_banach_and_ds_hochschild_hpl_are_harvested():
    class_m = _read("chapters/theory/topologization_class_m_original_complex_platonic.tex")
    class_m_flat = " ".join(class_m.split())
    chd = _read("chapters/theory/chiral_higher_deligne.tex")

    _assert_all(
        class_m,
        (
            r"\label{def:class-m-banach-algebra-body}",
            r"\label{def:class-m-banach-operation-constants}",
            r"C_k(\rho)",
            r"\sum_{k\ge2}C_k(\rho)<\infty",
            r"\label{prop:class-m-banach-mc-convergence}",
            r"\label{thm:chain-level-e3-original-complex-conditional}",
            r"\mathfrak{o}_{\mathrm{Ban}}(\cA,\rho)=0",
            r"\mathfrak{o}_{\oplus}(\cA,\rho)=0",
        ),
    )
    assert (
        "a Banach statement; it does not imply that the same operations preserve "
        "the raw direct sum."
    ) in class_m_flat
    assert "descent from the completed model to the raw direct sum requires" in class_m_flat

    _assert_all(
        chd,
        (
            r"\label{eq:ds-brst-w-algebra}",
            r"\label{eq:ds-brst-charge}",
            r"\label{eq:ds-hochschild-sdr}",
            r"\label{thm:chd-ds-hochschild}",
            r"\ClaimStatusConditional",
            "principal or hook-type",
            "or hook type",
            "bounded-shift/weightwise-convergence",
            r"\label{eq:chd-ds-hh}",
            r"\label{eq:ds-hochschild-hpl-brace}",
            r"Q_{\mathrm{DS}}h+hQ_{\mathrm{DS}}=\mathrm{id}-i\circ p",
            r"\sum_{T\in\mathrm{PRT}_n}",
        ),
    )


def test_modular_graph_arakelov_and_tensor_curvature_are_harvested():
    modular = _read("chapters/connections/modular_pva_quantization_core.tex")
    arakelov = _read("chapters/theory/factorization_swiss_cheese.tex")
    tensor = _read("chapters/theory/theorems_C_D_native_vol2_platonic.tex")
    tensor_flat = " ".join(tensor.split())

    _assert_all(
        modular,
        (
            r"\label{thm:explicit-relative-modular-graph-complex}",
            r"\mathfrak{o}_{\mathrm{GK}}(\Gamma)",
            r"\mathfrak{o}_{\mathrm{Hdg}}(\Gamma)",
            r"\FT_{\mathrm{rel}}(C)",
            r"D=D_{\mathrm{int}}+D_{\mathrm{sep}}+D_{\mathrm{nsep}}",
            r"D_0^2=0,\qquad",
            r"D_1^2=0,\qquad",
            r"D_0D_1+D_1D_0=0",
            r"[o_g]",
            r"H^2\!\left(\gr_F^g L,\ D_0+[\Theta^{(0)},-]\right)",
        ),
    )
    _assert_all(
        arakelov,
        (
            r"\label{def:elliptic-arakelov-green-kernel}",
            r"\partial\bar\partial G_\tau",
            "The residue-one chiral propagator used in this chapter is",
            r"P_\tau(z)",
            r"\bar\partial P_\tau=-2\pi i\,\mu_\tau",
            r"\label{thm:scalar-arakelov-genus-tower}",
            r"\lambda_g^{\mathrm{FP}}",
            r"\delta F_2^{\mathrm{cross}}(\cW_3)",
        ),
    )
    _assert_all(
        tensor,
        (
            r"\label{def:theoremC-total-gluing-datum}",
            r"\mathfrak s_{\mathcal A}^{C}",
            r"\begin{theorem}[Total shifted-symplectic consequence of Theorem~C;",
            r"\ClaimStatusConditional; licensing tags $\alpha+\beta+\gamma+\epsilon$]",
            r"\label{def:tensor-arakelov-graph-datum}",
            r"\mathfrak k_{\mathcal A}^{D}",
            r"\begin{theorem}[Tensor-Arakelov consequence of Theorem~D;",
            r"\ClaimStatusConditional; licensing tags $\alpha+\beta+\gamma+\epsilon$]",
            r"\label{thm:theoremD-tensor-arakelov}",
            r"\label{eq:tensor-arakelov-component-graph-integral}",
            r"P_{\alpha_e\beta_e}",
            r"C_v",
            r"\operatorname{tr}_{\cF}K(\cA)",
            r"\label{ex:theoremD-W3-four-channel-tensor}",
            "K_{TT} & K_{TW}",
            r"\delta F_2^{\mathrm{cross}}(\cW_3)=\frac{c+204}{16c}",
        ),
    )
    assert "The scalar identity \\(\\operatorname{tr}_{\\mathcal F}K=\\kappaChHodge(\\mathcal A)\\omega_g\\) does not determine it." in tensor_flat
    assert (
        "The scalar Vol~I Theorem~C identity holds on each "
        "\\textsc{uniform-weight} fibre"
        in tensor_flat
    )
    assert (
        "The tensor statement of "
        "Theorem~\\ref{thm:theoremD-tensor-arakelov} additionally "
        "requires the graph datum"
        in tensor_flat
    )
    assert r"\begin{proof}[Sketch]" not in tensor
    assert "Tensor-Arakelov upgrade of Theorem~D;\n\\ClaimStatusProvedHere" not in tensor
    assert "Total shifted-symplectic upgrade of Theorem~C;\n\\ClaimStatusProvedHere" not in tensor


def test_anomaly_completion_and_k3_borcherds_are_scoped_not_promoted():
    anomaly = _read("chapters/connections/anomaly_completed_core.tex")
    gravity = _read("chapters/connections/3d_gravity.tex")
    gravity_flat = " ".join(gravity.split())

    _assert_all(
        anomaly,
        (
            r"\label{def:tholog-curved-mc-anomaly-completion}",
            r"\label{eq:curved-mc-anomaly-completion}",
            r"d\eta+\frac12[\eta,\eta]+\Theta=0",
            r"d\eta+\eta^{2}+\Theta=0",
            "not the strict central-shadow transgression algebra",
            r"\label{ex:virasoro-curved-mc-anomaly-completion}",
            r"d\eta+\eta^{2}+\Theta_{\mathrm{Vir}}=0",
        ),
    )
    _assert_all(
        gravity,
        (
            r"\label{thm:k3-borcherds-operator-square}",
            r"\ClaimStatusConditional",
            r"\hypAmbientWtCpl",
            r"\effPfaffOrient+\effPBWnoExtra",
            "Assume a P1 datum",
            "cyclic-trace condition",
            "not a scalar-character statement",
            r"\label{eq:k3-borcherds-operator-pfaffian}",
            r"\label{eq:k3-borcherds-hall-chain}",
            r"\label{eq:k3-borcherds-chiral-chain}",
            "The finite Hall gate system is",
            "shadow comparison, not an object equivalence",
        ),
    )
    assert "not a consequence of the scalar identity" in gravity_flat


def test_supersymmetric_indices_are_annular_bar_supertraces():
    text = _read("chapters/connections/ordered_associative_chiral_kd_frontier.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\label{sec:susy-index-annular-bar}",
            r"\label{def:refined-annular-bar}",
            r"\label{def:annular-index}",
            r"\label{eq:annular-index-trace}",
            r"\label{thm:elliptic-genus-annular}",
            r"\label{thm:weak-jacobi-annular}",
            r"\label{def:n4-superconformal-index}",
            r"\label{thm:n4-index-annular}",
            r"\operatorname{str}_{B^{\mathrm{ann}}_\bullet(\cA)}",
        ),
    )
    assert "the physical index is the Euler characteristic of a chain complex" in flat


def test_chiral_springer_resolves_singular_steinberg_before_yangian_claim():
    text = _read("chapters/connections/spectral-braiding-frontier.tex")

    _assert_all(
        text,
        (
            r"\label{subsec:chiral-springer-resolution}",
            r"\label{constr:springer-sl2}",
            r"\label{constr:springer-type-A}",
            r"\label{constr:convolution-resolved}",
            r"\label{conj:yangian-springer}",
            r"\widetilde{\Steinberg}_b \to \Steinberg_b",
            r"\ClaimStatusConjectured",
        ),
    )


def test_yang_mills_boundary_lane_has_central_formality_and_mixed_couplings():
    text = _read("chapters/connections/ym_boundary_theory.tex")

    _assert_all(
        text,
        (
            r"\label{def:twisted-ym-boundary-datum}",
            r"\label{thm:twisted-ym-boundary-brst}",
            r"\label{thm:twisted-ym-tangent-center}",
            r"\label{thm:twisted-ym-central-formality}",
            r"\label{thm:twisted-ym-binary-mixed-couplings}",
            r"\label{cor:twisted-ym-multibody-couplings}",
            "central formality",
            "mixed couplings controlled by tensor products of reduced dual centers",
        ),
    )


def test_tail_false_patterns_are_scoped_by_product_formal_and_completed_ambient():
    foundations = _read("chapters/theory/foundations.tex")
    foundations_flat = " ".join(foundations.split())
    concordance = _read("chapters/connections/concordance.tex")
    concordance_flat = " ".join(concordance.split())
    e_infinity = _read("chapters/connections/e_infinity_topologization.tex")

    _assert_all(
        foundations,
        (
            r"\label{thm:recognition-foundations}",
            "product-formal rectangle prefactorization algebras",
            "on the raw direct sum",
            r"Theorem~\ref{thm:weight-completed-topologization-class-m}",
        ),
    )
    assert (
        "does not assert that the local operad recovers arbitrary Ran-space factorization data"
        in foundations_flat
    )
    assert "Its globalization is exactly" in foundations_flat
    assert "no theorem of global boundary Ran descent is claimed here" in foundations_flat
    assert "formal model suggests the conjectural geometric lift" in foundations_flat
    assert "not asserted as a global derived-stack theorem here" in foundations_flat
    assert "conjectural \\'etale shadow" in foundations_flat
    assert "The conjectural Koszul--\\'etale surface asserts" in foundations_flat
    assert "If this conjectural \\'etale surface is present" in foundations_flat
    assert "is expected to globalize by general categorical" not in foundations_flat
    assert "is expected to record the failure" not in foundations_flat
    assert "is expected to be \\'etale" not in foundations_flat
    assert "both source and target maps are expected to be quasi-isomorphisms" not in foundations_flat
    assert "Koszul inversion is then expected" not in foundations_flat
    _assert_all(
        concordance,
        (
            "BV datum",
            "Factorized logarithmic HT parametrix",
            "FM compactification, logarithmic forms, AOS relations, Stokes exactness",
            "Product-formal rectangle compatibility",
            "identifies the product-formal rectangle shadow, not arbitrary global Ran-space factorization data",
        ),
    )
    assert "BV--BRST landing on the chiral Koszul locus is the effectiveness clause" in concordance_flat
    assert r"\(\effKoszul\)" in concordance
    assert "must be carried as an \\(\\epsilon\\)-type datum wherever it is used" in concordance_flat
    assert "BV-BRST is expected to land on the chiral Koszul locus" not in concordance_flat
    _assert_all(
        e_infinity,
        (
            r"\label{thm:e-infinity-weightwise-inverse-limit}",
            "bounded conformal-weight window",
            "weight-completed inverse limit is degreewise eventually constant",
        ),
    )


def test_standard_fingerprint_does_not_claim_nonstandard_completeness():
    text = _read("chapters/examples/examples-complete-proved.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\label{thm:fingerprint-completeness-standard}",
            r"\label{rem:fingerprint-nonstandard-open}",
            "the five-slot fingerprint is not a complete invariant",
            "requires either verified free strong generation",
            "no bar-Quillen-completeness statement is asserted outside",
        ),
    )
    assert "a sixth coordinate recording twisted-sector data" in flat
    assert "The proof obligation reduces to" not in flat
    assert "proof uses two inputs not supplied by the five-slot fingerprint" not in flat


def test_tail_virasoro_rmatrix_and_pt_comparison_remain_conditional():
    core = _read("chapters/connections/spectral-braiding-core.tex")
    frontier = _read("chapters/connections/spectral-braiding-frontier.tex")

    _assert_all(
        core,
        (
            r"\label{eq:virasoro-unitarity-failure}",
            r"\label{prop:virasoro-first-quantum}",
            "operator-valued term",
            r"\Delta r^{(2)}_{\Ainf}(z)",
            r"\ClaimStatusConditional",
        ),
    )
    _assert_all(
        frontier,
        (
            r"\label{conj:R-equals-PT}",
            r"\label{thm:R-equals-PT-conditional}",
            "same four-point line channel",
            "channel, contour, normalization, and external momenta",
            r"\ClaimStatusConditional",
        ),
    )


def test_tail_page_desitter_and_borel_claims_require_raw_physical_data():
    gravity = _read("chapters/connections/3d_gravity.tex")
    flat = " ".join(gravity.split())

    _assert_all(
        gravity,
        (
            r"\label{conj:gravity-page-raw-transseries}",
            r"\label{prop:gravity-page-curve}",
            r"\label{conj:page-borel-singularity}",
            r"\label{rem:stokes-is-page}",
        ),
    )
    assert "It is not a consequence of scalar complementarity alone." in flat
    assert "has no factorial growth and supplies no Page Borel singularity" in flat
    assert "requires the physical entropy functional" in flat


def test_tail_harmonic_beta_law_is_conjectural_and_finite_beta_tempering_is_separate():
    beta = _read("chapters/theory/beta_N_closed_form_all_platonic.tex")
    closure = _read("chapters/theory/wn_tempered_closure_platonic.tex")

    _assert_all(
        beta,
        (
            "scaling-law frontier",
            r"\label{conj:beta-N-harmonic-closed-form}",
            r"\ClaimStatusConjectured",
            "all-rank proof obligation",
            "not polynomial in $N$",
        ),
    )
    _assert_all(
        closure,
        (
            r"\label{thm:wn-tempered-all-N}",
            "requires only existence of a finite",
            "closed form of $\\beta_N$ is conditional",
        ),
    )


def test_tail_m4_cone_one_wheel_and_superconnection_replacements_are_printed():
    rosetta = _read("chapters/examples/rosetta_stone.tex")
    celestial = _read("chapters/connections/celestial_holography_core.tex")
    log_ht = _read("chapters/connections/log_ht_monodromy_core.tex")

    assert "this is not a chain-level determination of" in " ".join(rosetta.split())
    _assert_all(
        celestial,
        (
            r"\label{thm:relative-bulk-boundary-reconstruction}",
            r"H^1\bigl(\Cone(F_\pi)\bigr)=H^2\bigl(\Cone(F_\pi)\bigr)=0",
            r"\label{def:one-wheel-sum-and-class}",
            r"\label{thm:one-wheel-reduction}",
            "completed one-wheel sum",
        ),
    )
    _assert_all(
        log_ht,
        (
            r"\label{def:superconnection}",
            r"\label{thm:master-curvature}",
            r"\label{thm:support-hierarchy}",
            r"\label{cor:recursive-support}",
            "bar-valued superconnection",
            "curvature is the Maurer-Cartan defect",
        ),
    )


def test_retired_input_pdf_slogans_do_not_reappear_on_guarded_surfaces():
    active_sources = (
        "chapters/connections/casimir_divisor_core_transport.tex",
        "chapters/connections/celestial_boundary_transfer_core.tex",
        "chapters/connections/kontsevich_integral.tex",
        "chapters/connections/fractional_ghost_chain_level_platonic.tex",
        "chapters/connections/spectral-braiding-frontier.tex",
        "chapters/connections/ym_boundary_theory.tex",
    )
    combined = "\n".join(_read(path) for path in active_sources)

    for retired in (
        "bar complex resembles scattering",
        "Chiral Springer gives Yangian",
        "the discriminant controls transport",
        "DS preserves spectral data",
        "Kontsevich integral follows from the Feynman transform",
    ):
        assert retired not in combined


def test_source_truth_bibliography_residue_is_closed():
    main = _read("main.tex")
    axioms = _read("chapters/theory/axioms.tex")
    einfty = _read("chapters/connections/e_infinity_topologization.tex")
    thqg_bv = _read("chapters/connections/thqg_bv_ht_extensions.tex")
    thqg_lines = _read("chapters/connections/thqg_line_operators_extensions.tex")
    active = "\n".join((main, axioms, einfty, thqg_bv, thqg_lines))
    einfty_flat = " ".join(einfty.split())

    _assert_all(
        main,
        (
            r"\bibitem{deBoerTjinRelation93}",
            "arXiv:hep-th/9302006",
            r"\bibitem{CLNS25}",
            "Duality via convolution of $\\mathcal{W}$-algebras",
            "arXiv:2203.01843",
            "unpublished draft, January 2009",
            r"\bibitem{FF92}",
        ),
    )
    _assert_all(
        axioms,
        (
            r"\cite{deBoerTjin93,deBoerTjinRelation93}",
            "Bouwknegt--Schoutens",
        ),
    )
    _assert_all(
        einfty_flat,
        (
            "Bouwknegt--Schoutens",
            "Linshaw/Bouwknegt--Schoutens OPE",
            "primary-field convention reviewed by Bouwknegt--Schoutens",
        ),
    )
    assert r"\cite{CLNS25}" in thqg_bv
    assert r"\cite{CLNS25}" in thqg_lines

    retired = (
        "TODO: librarian verification",
        "Bouwknegt--McCarthy--Pilch",
        "Bouwknegt-McCarthy-Pilch",
        "Linshaw/BMP",
        "BMP review",
        "BMP lower-Casimir",
        "BMP's Miura",
        r"\bibitem{CLNS24}",
        r"\cite{CLNS24}",
        "arXiv:2404.02410",
        r"\bibitem{FF90}",
        r"\cite{FF90}",
        "based on preprints circa 1990",
        "preprint identifier / published version",
        "final arXiv identifier",
        "exact preprint date/venue",
        "alternative candidate Bouwknegt",
    )
    for needle in retired:
        assert needle not in active


def test_log_ht_reduction_datum_is_hodge_relative_not_canonical_retract():
    core = _read("chapters/connections/log_ht_monodromy_core.tex")
    duplicate = _read("chapters/connections/log_ht_monodromy.tex")
    spectral = _read("chapters/connections/spectral-braiding-core.tex")
    combined = "\n".join((core, duplicate, spectral))

    assert r"\label{thm:hodge-hbar-retract}" in core
    for text in (core, duplicate):
        flat = " ".join(text.split())
        _assert_all(
            text,
            (
                r"Hodge $\h$-adic strong deformation retract",
                "After choosing a Hermitian metric, equivalently a tree-level Hodge contraction",
                "canonical only relative to the chosen tree-level Hodge contraction",
                "different contractions give homotopic retract data",
                r"an $\h$-adic retract is induced by any chosen tree-level Hodge contraction",
                "without any further quantum concentration hypothesis",
                "functorially in its chosen reduction data",
            ),
        )
        _assert_all(
            flat,
            (
                r"Q_{0,n}h_{0,n}+h_{0,n}Q_{0,n}=\id-i_{0,n}p_{0,n}",
                r"tree-level degree-zero resolution condition on $(E_n,Q_{0,n})$, together with a compatible tree-level Hodge contraction",
            ),
        )
    spectral_flat = " ".join(spectral.split())
    _assert_all(
        spectral_flat,
        (
            "fix a compatible tree-level Hodge contraction defining the reduced line category",
            "applies after choosing this tree-level Hodge contraction",
            r"admits an $\h$-adic strong deformation retract $(p,i,h)$ relative to that choice",
        ),
    )

    retired = (
        "Canonical $\\h$-adic strong deformation retract",
        "canonical $\\h$-adic strong deformation retract",
        "canonical $\\h$-adic retract",
        "the canonical $\\h$-adic retract exists",
        "or equivalently a canonical reduction datum",
        "determines canonically",
        "without any further concentration or retract hypothesis",
        "Transfer this superconnection along the canonical retract",
        "identifies canonically with $H_n",
        "construct a canonical $\\h$-adic strong deformation retract",
        "reduces canonically to an ordinary logarithmic connection",
        "thm:canonical-retract",
    )
    for needle in retired:
        assert needle not in combined


def test_affine_monodromy_is_braided_comparison_not_pointwise_r_matrix_equality():
    core = _read("chapters/connections/log_ht_monodromy_core.tex")
    iv = _read("compute/tests/test_affine_monodromy_identification_iv.py")
    frontier = _read("chapters/connections/spectral-braiding-frontier.tex")
    active = "\n".join((core, iv, frontier))
    core_flat = " ".join(core.split())
    frontier_flat = " ".join(frontier.split())

    _assert_all(
        core,
        (
            "After choosing the KZ/Drinfeld associator and the corresponding",
            "braided-monoidally equivalent to the quantum-group braiding",
            "through this monodromy functor",
        ),
    )
    assert "not by a bare equality of functions" in core_flat
    assert "proved on the fixed reduced evaluation comparison surface" in frontier_flat
    assert "a compatible Hodge reduction transports the reduced HT monodromy" in frontier_flat
    assert "It is not a pointwise equality of spectral and quantum-group" in frontier_flat
    _assert_all(
        iv,
        (
            "programme spectral half-monodromy on evaluation modules is compared",
            "only up to the Drinfeld associator",
            "comparison is braided-monoidal, not pointwise equality of",
            "HT monodromy should compare with q-braiding",
        ),
    )

    retired = (
        "spectral R-matrix restricted to evaluation modules equals",
        "equals the quantum-group R-matrix",
        "HT monodromy = KZ monodromy = q-R-matrix",
        "should equal spectral R",
        "bare equality of meromorphic R-functions",
        "For the affine lineage, this is proved unconditionally",
    )
    for needle in retired:
        assert needle not in active


def test_abelian_cs_example_uses_fixed_reduction_data_for_affine_bridge():
    text = _read("chapters/examples/examples-worked.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            "under the displayed Gaussian and reduction data",
            r"licensing tags $\alpha+\beta+\gamma$",
            "fix the Gaussian HT gauge data, a compatible Hodge reduction",
            "KZ/Drinfeld associator and",
            "The four ingredients are available on the fixed abelian comparison",
            r"r^{\mathrm{coll}}_{\cH}(z)=\frac{k\,\Omega_{\cH}}{z}",
        ),
    )
    assert (
        "The reduced line category on the evaluation core is compared with "
        "$\\mathrm{Rep}_q(\\mathfrak{g})$ by Theorem~\\ref{thm:affine-monodromy-identification}"
        in flat
    )
    assert (
        "The comparison is braided-monoidal; it is not a pointwise equality of "
        "spectral and quantum-group $R$-matrices."
    ) in flat

    retired = (
        "cross-volume bridge applies unconditionally",
        "The four ingredients are each unconditionally available",
        r"\mathcal{C}_{\mathrm{line}}^{\mathrm{red}}\simeq\mathrm{Rep}_q(\mathfrak{g})$ on the evaluation core is proved unconditionally",
        "HT monodromy = KZ monodromy = q-R-matrix",
    )
    for needle in retired:
        assert needle not in text


def test_bicoloured_nine_tuple_reconstructs_sc_algebra_not_physical_ht_qft():
    text = _read("chapters/theory/foundations.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\label{cor:bicoloured-primitive-universality-3dHT}",
            r"canonical $\SCchtop$-algebra",
            "local algebraic",
            "bulk--boundary shadow of a $3$d HT theory",
            "separate BV/HT",
            "realisation, descent, and anomaly data",
            "does not by",
            "itself supply a BV action, QME solution, analytic integration cycle, or",
            "global descent package",
        ),
    )
    assert (
        "every coherent bicoloured nine-tuple extends to a canonical "
        "$\\SCchtop$-algebra"
    ) in flat

    retired = (
        "every coherent bicoloured nine-tuple extends to a\ncanonical $3$d HT theory",
        "every coherent bicoloured nine-tuple extends to a canonical $3$d HT theory",
        "canonical $3$d HT theory on $C \\times \\R$",
        "by itself supply a physical HT QFT",
    )
    for needle in retired:
        assert needle not in text


def test_m2_zeta_scalar_is_regularized_benchmark_not_unconditional_kappa():
    text = _read("chapters/connections/thqg_ht_bbl_extensions.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            r"\begin{theorem}[Components of $\Theta_{M2}$; \ClaimStatusProvedHere]",
            r"\label{def:m2-zeta-regularization-datum}",
            r"\label{prop:m2-zeta-regularized-kappa}",
            r"\label{conj:m2-zeta-comparison}",
            r"\ClaimStatusConjectured",
            r"\kappa_{M2}^{\zeta}=K^2\zeta(-1)=-\frac{K^2}{12}",
            r"F_1^{\zeta}(A_{M2,\infty})",
            r"F_1^{\zeta}=\kappa_{M2}^{\zeta}/24",
            r"\kappa_{M2}^{\zeta}=-K^2/12$; comparison conjectural",
        ),
    )
    _assert_all(
        flat,
        (
            "The theorem does not identify this characteristic with a number.",
            r"identifying $\kappa_{M2}$ with $\kappa_{M2}^{\zeta}$ is the comparison conjecture",
            "Proposition~\\ref{prop:m2-zeta-regularized-kappa} proves the value of the regularized scalar once the regularization datum has been installed.",
            "replacing it by the actual \\(\\kappa_{M2}\\) requires Conjecture~\\ref{conj:m2-zeta-comparison}",
        ),
    )

    retired = (
        r"\ClaimStatusProvedHere{} except genus-$1$ value",
        r"\ClaimStatusConditional modulo $\kappa_{M2}$ value",
        r"$\kappa_{M2} = -K^2/12$ heuristic",
        r"\kappa_{M2} = -K^2/12$ is heuristic",
        r"F_1=\kappa/24",
        "well-motivated heuristic, not a theorem",
        "unique prescription compatible with modular invariance",
    )
    for needle in retired:
        assert needle not in text


def test_gravitational_s_duality_exactness_is_algebraic_not_physical_duality():
    text = _read("chapters/connections/thqg_gravitational_s_duality.tex")
    flat = " ".join(text.split())

    _assert_all(
        text,
        (
            "algebraic gravitational $S$-duality of the",
            r"\begin{theorem}[Algebraic gravitational $S$-duality; \ClaimStatusProvedHere]",
            "Algebraic $S$-duality as exact symmetry",
            "The algebraic gravitational $S$-duality constructed in this",
            "chapter is exact on the chiral Koszul--Verdier surface",
            "not a Montonen--Olive-style physical duality conjecture",
            "definition of algebraic gravitational $S$-duality in this framework",
        ),
    )
    assert (
        "real forms, analytic sewing, and dynamical-metric interpretations are separate inputs"
        in flat
    )

    retired = (
        r"\begin{theorem}[Gravitational $S$-duality; \ClaimStatusProvedHere]",
        r"\subsection{$S$-duality as exact symmetry: the all-genera theorem}",
        r"The gravitational $S$-duality is \emph{exact}, not conjectural",
        "it is a theorem ($\\sigma(\\Theta_\\cA) = \\Theta_{\\cA^!}$), not a duality conjecture",
        "definition of gravitational $S$-duality in this framework",
    )
    for needle in retired:
        assert needle not in text


def test_summary_layer_keeps_harvested_comparison_boundaries():
    main = _read("main.tex")
    preface = _read("chapters/frame/preface.tex")
    intro = _read("chapters/theory/introduction.tex")
    combined = "\n".join((main, preface, intro))
    flat = " ".join(combined.split())
    intro_flat = " ".join(intro.split())

    _assert_all(
        main,
        (
            r"the candidate $\Phi_4^{\mathrm{HK}}$",
            "asks for the $3$d HT-QFT climax through twistor-$S^1$ reduction",
            "are comparison targets until the twistor",
            "descent and $11$d lift are constructed",
        ),
    )
    _assert_all(
        preface,
        (
            r"CY-$d$-enlarged sparse support",
            r"$\{0,1,2,d\}$",
            r"$\{0,1,2,4\}$",
            r"\Phi_4^{\mathrm{HK}}$ candidate becomes",
        ),
    )
    assert "no degree-$3$ occupancy asserted in the CY-$4$ case" in flat
    _assert_all(
        intro_flat,
        (
            "remaining datum is a moduli-stratification construction",
            "conditional on the Verdier-helicity comparison",
            "tree-level Mellin/Parke--Taylor propagator comparison",
            r"not an unconditional consequence of \(\SCchtop\) alone",
            "The required comparison datum is the map between",
            "the displayed equality is a comparison problem, not a theorem",
            r"candidate $\Phi_4^{\mathrm{HK}}$",
            "with three independent comparison targets",
            "These routes are evidence and verification targets, not a proof",
        ),
    )
    assert r"\paragraph{The CY-$4$ HK-restricted $\Phi_4$ comparison problem.}" in intro
    assert r"$\{0,1,2,3,4\}$" not in combined
    assert r"$\{0,1,\ldots,d\}$" not in combined

    retired = (
        "realises the $3$d HT-QFT climax via twistor-$S^1$ reduction",
        "reached through DT-$4$ sheaf counting",
        "The Vol~III functor $\\Phi$ extends from CY-$3$",
        "HK-restricted $\\Phi_4$ proof obligation",
        "pairing formula to be proved",
        "verified along three disjoint verification paths",
        "Its proof obligation is",
        "the proof obligation is the form-type comparison",
        "proof obligation is to construct the corresponding moduli",
        "The remaining proof obligation is the cross-volume comparison map",
    )
    for needle in retired:
        assert needle not in flat
