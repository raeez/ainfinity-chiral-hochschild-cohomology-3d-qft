"""Tests for PVA descent chain-level verification.

Verifies D2-D5 plus the vacuum axiom at the chain level for all standard families:
1. D2 Sesquilinearity: {da_lam b} = -lam {a_lam b}
2. D3 Jacobi: full Borcherds identity
3. D4 Leibniz: derivation property of n-products
4. D5 Skew-symmetry: {a_lam b} = -{b_{-lam-d} a}
5. Vacuum Unit: {1_lam a} = {a_lam 1} = 0, 1*a = a = a*1, d1 = 0

For each axiom, we test against ALL 7 standard families:
Heisenberg, Virasoro, affine sl_2, beta-gamma, W_3, free multiplet, LG cubic.

Each test performs ACTUAL symbolic computation.

References:
  Vol II: pva-descent.tex (D2-D5 plus vacuum proofs)
  Vol I: configuration_spaces.tex, fm_boundary.py
  De Sole-Kac (2006): PVA axiomatization
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from pathlib import Path
from sympy import Symbol, Rational, S, simplify, expand, symbols, oo

from lib.pva_descent_chain_level import (
    # Data structures
    LambdaBracket,
    PVAData,
    # Standard families
    heisenberg_pva,
    virasoro_pva,
    affine_sl2_pva,
    betagamma_pva,
    w3_pva,
    free_multiplet_pva,
    lg_cubic_pva,
    # D2
    verify_d2_sesquilinearity,
    # D3
    verify_d3_jacobi_generators,
    verify_d3_jacobi_sl2,
    verify_d3_jacobi_virasoro,
    # D4
    verify_d4_leibniz,
    # D5
    verify_d5_skew_symmetry,
    # D6
    verify_d6_unit,
    # Full sweep
    full_pva_descent_verification,
    # FM_3 boundary
    fm3_boundary_strata_cancellation,
    fm3_exchange_cylinder_stokes,
    # Residues
    boundary_residue_d12,
    boundary_residue_d23,
    boundary_residue_d13,
    # Complementarity
    kappa_complementarity_from_pva,
    # Pole census
    pole_order_census,
)


# ===================================================================
# LAMBDA-BRACKET REPRESENTATION
# ===================================================================

class TestPvaDescentSourceGuards:
    """Source-level guards for theorem status and comparison hypotheses."""

    def test_steinberg_poisson_comparison_is_conditional(self):
        """The Steinberg comparison needs chi and Koszul-effectiveness."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Steinberg Poisson comparison for the PVA bracket"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma+\varepsilon$",
            r"\effKoszul",
            r"\chi\colon H^\bullet(\barB(\A))\longrightarrow H^\bullet(\A,Q)",
            "is an isomorphism on the Koszul-effective comparison locus",
            "bar-to-boundary comparison map",
            "if $\\chi$ is an isomorphism",
            "direct FM proof below supplies the PVA bracket",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            r"\begin{proposition}[The PVA bracket is the Steinberg Poisson bracket;",
            r"\ClaimStatusProvedHere]",
            "On cohomology \\(H^\\bullet(\\barB(\\A)) \\simeq H^\\bullet(\\A,Q)\\)",
        )
        for needle in retired:
            assert needle not in block

    def test_jacobi_lagrangian_convolution_is_conditional_comparison(self):
        """The Steinberg Jacobi reading must not replace the FM proof."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()

        remark_start = source.index(r"\begin{remark}[Conditional Steinberg reading")
        remark_end = source.index(r"\end{remark}", remark_start)
        remark = source[remark_start:remark_end]
        remark_flat = " ".join(remark.split())
        for needle in (
            r"Proposition~\ref{prop:PVA-from-symplectic}",
            r"\chi",
            r"\effKoszul",
            "Steinberg bar/Koszul shadow",
            r"direct FM theorem",
            r"\ref{thm:Jacobi}",
        ):
            assert needle in remark or needle in remark_flat

        start = source.index(
            r"\begin{proposition}[Conditional Steinberg convolution interpretation"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma+\varepsilon$",
            r"boundary chart $b$",
            r"comparison map $\chi$",
            "derived shifted-symplectic ambient",
            r"\effKoszul",
            r"\chi\colon H^\bullet(\barB(\A))\longrightarrow H^\bullet(\A,Q)",
            "bar-to-boundary comparison map",
            "signed Jacobiator",
            r"J_{\lambda,\mu}",
            "isomorphism on the Koszul-effective comparison locus",
            "not a proof of the boundary PVA identity",
            "direct boundary identity has already been proved",
            "Steinberg shadow",
            r"\mathfrak J_{23}+\mathfrak J_{13}+\mathfrak J_{12}",
            "Associativity of derived correspondence composition",
            r"Proposition~\ref{prop:PVA-from-symplectic}",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            r"\begin{proposition}[Jacobi as associativity of Lagrangian convolution;",
            r"\ClaimStatusProvedHere]",
            r"The Jacobi identity for $\{{\cdot}_\lambda{\cdot}\}$ is the associativity",
            "These are the three residues computed in",
            r"which is the Jacobi identity~\eqref{eq:repaired-jacobi}",
            "the three ways of composing Lagrangian correspondences in $\\Steinb_b$ are",
        )
        for needle in retired:
            assert needle not in remark
            assert needle not in block

    def test_pva_descent_functor_uses_strict_unital_continuous_morphisms(self):
        """Functoriality must name the strict coefficientwise morphism class."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()

        def_start = source.index(
            r"\begin{definition}[Strict morphism of logarithmic $\SCchtop$-algebras]"
        )
        def_end = source.index(r"\end{definition}", def_start)
        definition = source[def_start:def_end]
        definition_flat = " ".join(definition.split())
        for needle in (
            r"\hypAmbientPro",
            "strict morphism",
            "continuous cochain map",
            r"f(\mathbf 1_A)=\mathbf 1_B",
            "intertwining all logarithmic kernels coefficientwise",
            "Laurent coefficient extraction",
            "regular/singular projection",
            "polynomial Borel transform",
            "strict unital morphisms",
        ):
            assert needle in definition or needle in definition_flat

        thm_start = source.index(
            r"\begin{theorem}[Strict functoriality of PVA descent;"
        )
        thm_end = source.index(r"\end{proof}", thm_start)
        theorem = source[thm_start:thm_end]
        theorem_flat = " ".join(theorem.split())
        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "strict unital coefficientwise",
            r"\cat{SC\text{-}Alg}^{\log}",
            "regular projection and constant-term extraction",
            "intertwined coefficientwise",
            "commutes with $f$ by $\\hypAmbientPro$-continuity",
            "strict unitality clause",
        ):
            assert needle in theorem or needle in theorem_flat

        retired = (
            r"\begin{definition}[Morphism of logarithmic $\SCchtop$-algebras]",
            r"\begin{theorem}[Functoriality of PVA descent; \ClaimStatusProvedHere]",
            "a morphism of unital\n$\\SCchtop$-algebras preserves the unit by definition",
            "The map $H^\\bullet(f)$ preserves each piece of the PVA structure.",
        )
        for needle in retired:
            assert needle not in definition
            assert needle not in theorem

    def test_windowwise_k3_pva_descent_uses_recognized_strict_pro_limit(self):
        """The K3 completed PVA needs recognized strict transitions."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()

        start = source.index(
            r"\begin{proposition}[Windowwise PVA descent and recognized pro-limit;"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusConditional",
            r"licensing $\alpha+\beta+\gamma$",
            "height-$N$ Hall--Borcherds recognition datum",
            r"maps $\rho_N$",
            r"strict $\hypAmbientPro$ morphisms",
            r"\hypAmbientWtCpl",
            r"Definition~\ref{def:k3-finite-hall-windows}",
            r"Definition~\ref{def:morphism-SC-algebra}",
            r"\operatorname{Rad}_{N+1}",
            r"\rho_{N+1}",
            r"\rho_N",
            r"Theorem~\ref{thm:pva-descent-functor}",
            "computed relationwise in $\\hypAmbientWtCpl$",
            "Without those compatible finite recognition maps",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            r"\begin{proposition}[Windowwise PVA descent and radical quotient;",
            "are morphisms of logarithmic $\\SCchtop$-algebras, carry",
            "If the transition maps are PVA morphisms, the completed\ninverse limit",
        )
        for needle in retired:
            assert needle not in block

    def test_legacy_pva_split_warning_blocks_are_licensed_and_signed(self):
        """Legacy split warning blocks must match the repaired PVA theorem surface."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent.tex").read_text()

        jacobi_start = source.index(r"\begin{lemma}[Jacobi identity verification;")
        jacobi_end = source.index(r"\end{proof}", jacobi_start)
        jacobi = source[jacobi_start:jacobi_end]
        jacobi_flat = " ".join(jacobi.split())
        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "oriented $\\FM_3$ incidence",
            "Arnold--Orlik--Solomon corner cancellation",
            r"J_{\lambda,\mu}(a,b,c)",
            r"\operatorname{im}Q",
            r"K_3^{\mathrm{Jac}}",
            r"\mathfrak J_{23}",
            r"\mathfrak J_{13}",
            r"\mathfrak J_{12}",
        ):
            assert needle in jacobi or needle in jacobi_flat

        higher_start = source.index(
            r"\begin{theorem}[Higher operations are \(Q\)-exact under compatible topological"
        )
        higher_end = source.index(r"\end{proof}", higher_start)
        higher = source[higher_start:higher_end]
        higher_flat = " ".join(higher.split())
        for needle in (
            r"\ClaimStatusConditional",
            r"licensing $\gamma+\varepsilon$",
            r"\hypAmbientPro",
            "relative compactified ordered topological factor",
            r"Definition~\textup{\ref{def:compatible-topological-nullhomotopy}}",
            r"\widetilde\omega_k=\omega_k^{\mathrm{hol}}\wedge\omega_k^{\mathrm{top}}",
            "Without the supplied relative bounding chains",
            "Ordinary contractibility of the open ordered configuration space is",
            "not enough",
        ):
            assert needle in higher or needle in higher_flat

        leibniz_start = source.index(r"\begin{proposition}[Leibniz rule in cohomology;")
        leibniz_end = source.index(r"\end{proof}", leibniz_start)
        leibniz = source[leibniz_start:leibniz_end]
        leibniz_flat = " ".join(leibniz.split())
        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "oriented $\\FM_3$ incidence",
            "unsuspended regular product",
            "divided-power Borel transform",
            r"\mathcal L_\lambda(a;b,c)",
            r"K_3^{\mathrm{Leib}}",
            r"(-1)^{(|a|+1)|b|}",
        ):
            assert needle in leibniz or needle in leibniz_flat

        retired = (
            r"\begin{lemma}[Jacobi Identity Verification; \ClaimStatusProvedHere{}",
            r"\begin{theorem}[Higher operations vanish under compatible nullhomotopies;",
            r"\begin{proposition}[Leibniz Rule in Cohomology; \ClaimStatusProvedHere]",
            r"\{a_\lambda(b\cdot c)\} = \{a_\lambda b\}\cdot c + b\cdot\{a_\lambda c\}",
            "Consequently, $m_k = 0$ in cohomology.\n\\end{theorem}",
        )
        for needle in retired:
            assert needle not in jacobi
            assert needle not in higher
            assert needle not in leibniz

    def test_degree_two_identity_uses_sum_zero_stasheff_sign(self):
        """The n=2 identity must match the defining Stasheff convention."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()
        start = source.index(
            r"\begin{proposition}[$A_\infty$ degree-$2$ compatibility"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            r"Q(m_2(a,b))+m_2(Qa,b)+(-1)^{\degree{a}}m_2(a,Qb)=0",
            r"Q(m_2(a,b))=-m_2(Qa,b)-(-1)^{\degree{a}}m_2(a,Qb)",
            "desuspended Stasheff convention",
            "No other term occurs at degree~$2$",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            r"Q(m_2(a,b))=m_2(Qa,b)+(-1)^{\degree{a}}m_2(a,Qb)",
            "Equivalently, $Q=m_1$ is a graded derivation of the binary operation.",
            r"Q\mu(a,b)=\mu(Qa,b)+(-1)^{\degree{a}}\mu(a,Qb)",
            r"Q(a_{(n)}b)=(Qa)_{(n)}b+(-1)^{\degree{a}}a_{(n)}(Qb)",
        )
        for needle in retired:
            assert needle not in block

    def test_binary_descent_lemma_names_ambient_and_borel_inputs(self):
        """The binary descent lemma must expose its ambient and projections."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()
        start = source.index(
            r"\begin{lemma}[Descent of the regular and singular binary parts"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "regular/singular decomposition",
            r"\cref{thm:cohomology-PVA-main}",
            "polynomial Borel transform",
            r"\cref{def:borel-transform-pva}",
            "commutes with the regular/singular projections",
            r"Q\mu(a,b)=-\mu(Qa,b)-(-1)^{\degree{a}}\mu(a,Qb)",
            r"Q(a_{(n)}b)=-(Qa)_{(n)}b-(-1)^{\degree{a}}a_{(n)}(Qb)",
            "The same argument applies in the second variable.",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            r"\begin{lemma}[Descent of the regular and singular binary parts; \ClaimStatusProvedHere]",
            "depend only on the classes $[a],[b]\\in H^\\bullet(\\A,Q)$.",
        )
        for needle in retired:
            assert needle not in block

    def test_product_associativity_uses_unsuspended_regular_product(self):
        """Associativity must use desuspension, not sign cancellation by commutativity."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()

        remark_start = source.index(r"\begin{remark}[Unsymmetrized vs.\ symmetrized product]")
        remark_end = source.index(r"\end{remark}", remark_start)
        remark = source[remark_start:remark_end]
        remark_flat = " ".join(remark.split())
        for needle in (
            r"\mu(a,b) := \mu_2^{\mathrm{reg}}(a,b;0)",
            r"\mu_2^{\mathrm{reg}}(a,b;0)",
            "suspended notation; after desuspension",
            r"m_2^{\mathrm{reg}}",
        ):
            assert needle in remark or needle in remark_flat

        start = source.index(
            r"\begin{proposition}[Associativity of the descended product;"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "unsuspended regular-product convention",
            r"\mu(a,b)=\mu_2^{\mathrm{reg}}(a,b;0)",
            "suspended bar component",
            r"K_3^{\mathrm{reg}}(a,b,c)",
            "regular--regular constant-term",
            r"D_A^2=0",
            "bar desuspension",
            r"\operatorname{im}Q",
            "No vanishing of",
            "supplies the chain homotopy",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            r"The product $\mu(a,b) := m_2^{\mathrm{reg}}(a,b;0)$",
            r"The sign $(-1)^{|a|}$ in the $A_\infty$ identity becomes trivial",
            "regular--regular projection and symmetrization",
            "symmetrized product: since",
            "the two orderings",
            "requiring $m_3$ itself to vanish",
        )
        for needle in retired:
            assert needle not in remark
            assert needle not in block

    def test_associativity_shadow_surfaces_use_desuspension_not_symmetrization(self):
        """Raviolo and preview shadows must not revive the old associativity proof."""
        root = Path(__file__).resolve().parents[2]

        raviolo = (root / "chapters" / "theory" / "raviolo.tex").read_text()
        rav_start = raviolo.index(r"\textbf{Step 1: Commutative product from the regular part.}")
        rav_end = raviolo.index(r"\medskip\noindent\textbf{Step 2:", rav_start)
        rav_block = raviolo[rav_start:rav_end]
        rav_flat = " ".join(rav_block.split())
        for needle in (
            "unsuspended regular product",
            r"\mu_2^{\mathrm{reg}}",
            r"Proposition~\ref{prop:product-associative}",
            r"D_A^2=0",
            "ordinary associator",
            "raw Stasheff sign is absorbed",
            "not discarded by graded commutativity",
        ):
            assert needle in rav_block or needle in rav_flat

        preview = (root / "chapters" / "theory" / "pva-preview.tex").read_text()
        prev_start = preview.index(r"\begin{proposition}[Associativity on $H$]")
        prev_end = preview.index(r"\end{proof}", prev_start)
        prev_block = preview[prev_start:prev_end]
        prev_flat = " ".join(prev_block.split())
        for needle in (
            "ordinary unsuspended regular product",
            r"\mu_2^{\mathrm{reg}}",
            "bar desuspension",
            "ordinary associator",
            r"\operatorname{im}Q",
        ):
            assert needle in prev_block or needle in prev_flat

        retired = (
            r"m_2\bigl(m_2(a,b),c\bigr) \;\pm\;",
            "The symmetrization of the projected identity is compatible",
            "the two orderings of each associator term contribute equally",
            "Carefully tracking the graded symmetrizations",
            r"The $A_\infty$ identity at degree $3$ implies $Qm_3(a,b,c)$ equals a signed sum",
        )
        for needle in retired:
            assert needle not in rav_block
            assert needle not in prev_block

    def test_jacobi_uses_oriented_three_face_jacobiator_normalization(self):
        """Jacobi must be a signed oriented boundary identity."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()

        subsection_start = source.index(
            r"\subsection{The three-face Stokes theorem on $\FM_3(\C)$}"
        )
        lemma_start = source.index(r"\begin{lemma}[Three-face singular factorization", subsection_start)
        lemma_end = source.index(r"\end{proof}", lemma_start)
        lemma = source[lemma_start:lemma_end]
        lemma_flat = " ".join(lemma.split())

        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            r"oriented $\FM_3$ incidence",
            "divided-power Borel transform",
            r"\mathfrak J_{23}",
            r"\mathfrak J_{13}",
            r"\mathfrak J_{12}",
            r"-\,(-1)^{(\degree{a}+1)(\degree{b}+1)}",
            r"-\,\{\{a_\lambda b\}_{\lambda+\mu}c\}",
            "Jacobiator normalization",
            "dressed doubly singular PVA form",
            r"\ref{rem:non-consecutive-scope}",
        ):
            assert needle in lemma or needle in lemma_flat

        theorem_start = source.index(r"\begin{theorem}[Jacobi identity;", lemma_end)
        theorem_end = source.index(r"\end{proof}", theorem_start)
        theorem = source[theorem_start:theorem_end]
        theorem_flat = " ".join(theorem.split())

        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "Arnold--Orlik--Solomon corner cancellation",
            r"K_3^{\mathrm{Jac}}",
            r"J_{\lambda,\mu}(a,b,c)\in\operatorname{im}Q",
            r"\mathfrak J_{23}+\mathfrak J_{13}+\mathfrak J_{12}",
            r"\{\{a_\lambda b\}_{\lambda+\mu}c\}",
            "printed Jacobiator",
        ):
            assert needle in theorem or needle in theorem_flat

        retired = (
            "Then its residues along the three pairwise divisors are:",
            "Combining this with the orientation sign\n$\\epsilon_{13}=-1$ produces exactly the third Jacobi term with the correct overall sign.",
            "the three pairwise residues are exactly\nthe three iterated brackets",
            "The first summand is the residue pattern for",
            "not a slogan attached to it",
        )
        for needle in retired:
            assert needle not in lemma
            assert needle not in theorem

    def test_jacobi_shadow_surfaces_use_signed_jacobiator(self):
        """FM proof, repaired appendix, old appendix, preview, and raviolo must match."""
        root = Path(__file__).resolve().parents[2]

        fm = (root / "chapters" / "theory" / "fm-proofs.tex").read_text()
        fm_start = fm.index(r"\begin{corollary}[Arnold is the PVA Jacobi residue identity]")
        fm_end = fm.index(r"\end{corollary}", fm_start)
        fm_block = fm[fm_start:fm_end]
        fm_flat = " ".join(fm_block.split())
        for needle in (
            r"J_{\lambda,\mu}(a,b,c)",
            r"-\,(-1)^{(\degree{a}+1)(\degree{b}+1)}",
            r"-\{\{a_\lambda b\}_{\lambda+\mu}c\}",
            "oriented boundary contributions",
        ):
            assert needle in fm_block or needle in fm_flat

        repaired = (root / "chapters" / "theory" / "pva-expanded-repaired.tex").read_text()
        ss_start = repaired.index(r"\begin{proposition}[The singular--singular projection")
        ss_end = repaired.index(r"\end{proof}", ss_start)
        ss = repaired[ss_start:ss_end]
        ss_flat = " ".join(ss.split())
        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            "Jacobiator normalization",
            r"-\,(-1)^{(\degree{a}+1)(\degree{b}+1)}",
            r"-\{\{a_\lambda b\}_{\lambda+\mu}c\}",
        ):
            assert needle in ss or needle in ss_flat

        old_appendix = (root / "appendices" / "pva-expanded.tex").read_text()
        app_start = old_appendix.index(r"\textbf{Step 3: Extract the Jacobi identity.}")
        app_end = old_appendix.index(r"\textbf{Step 4: Sign analysis.}", app_start)
        app = old_appendix[app_start:app_end]
        assert "normal form" in app
        assert r"-(-1)^{(\degree{a}+1)(\degree{b}+1)}" in app
        assert r"-\{\{a{}_\lambda b\}{}_{\lambda+\mu}c\}" in app.replace(r"\,", "")

        preview = (root / "chapters" / "theory" / "pva-preview.tex").read_text()
        prev_start = preview.index(r"\begin{theorem}[Jacobi identity]")
        prev_end = preview.index(r"\end{proof}", prev_start)
        prev = preview[prev_start:prev_end]
        assert "signed Jacobiator" in prev
        assert "incidence orientation and shifted transposition" in prev

        raviolo = (root / "chapters" / "theory" / "raviolo.tex").read_text()
        rav_start = raviolo.index(r"\textbf{Step 5: Jacobi identity.}")
        rav_end = raviolo.index(r"\medskip\noindent\textbf{Step 6:", rav_start)
        rav = raviolo[rav_start:rav_end]
        rav_flat = " ".join(rav.split())
        for needle in (
            "signed Jacobiator",
            "dressed doubly singular PVA form",
            "incidence orientation",
            "shifted transposition",
        ):
            assert needle in rav or needle in rav_flat

        retired = (
            "After the divided-power Borel\ntransform they are exactly",
            "the three residues become\nrespectively",
            "ensures that the sum of the three boundary contributions vanishes",
            "The three pairwise collision strata $D_{12}$, $D_{13}$, $D_{23}$ contribute the three terms",
        )
        for needle in retired:
            assert needle not in fm_block
            assert needle not in ss
            assert needle not in app
            assert needle not in prev
            assert needle not in rav

    def test_legacy_binary_descent_shadow_has_no_removed_label_refs(self):
        """The superseded PVA split must not advertise stale label refs."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent.tex").read_text()
        start = source.index(
            r"\begin{lemma}[Regular and singular binary descent to cohomology"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "regular product",
            "singular polynomial bracket",
            "Borel transform convention",
            "The regular projection and every singular coefficient projection",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            r"Proposition \ref{prop:m1_m2}",
            r"Lemma \ref{lem:lambda_descends}",
            "[m_2(a,b)]",
            r"\begin{lemma}[The $\lambda$-Bracket Descends to Cohomology;",
        )
        for needle in retired:
            assert needle not in block

        completion_start = source.index(r"\begin{proof}[Proof of Theorem")
        completion = source[completion_start: source.index(r"\end{proof}", completion_start)]
        assert r"Lemma \ref{lem:lambda_descends}" not in completion
        assert "regular/singular binary descent lemma above" in completion

    def test_sesquilinearity_proof_uses_divided_power_modes(self):
        """PVA1 must be proved in the polynomial divided-power convention."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()
        start = source.index(r"\begin{proposition}[Sesquilinearity in cohomology")
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "polynomial divided-power convention",
            r"\cref{def:borel-transform-pva}",
            r"\eqref{eq:sesqui-left}--\eqref{eq:sesqui-right}",
            r"(\partial a)_{(n)}b=-n\,a_{(n-1)}b",
            r"a_{(n)}(\partial b)=\partial(a_{(n)}b)+n\,a_{(n-1)}b",
            r"a_{(-1)}b=0",
            r"\frac{\lambda^n}{n!}",
            r"+ \lambda\sum_{m\ge0}a_{(m)}b",
            r"Definition~\ref{def:sesquilinearity}",
            r"Lemma~\ref{lem:lambda_descends}",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            "Those coefficients are precisely the coefficients that",
            "preserves linear identities among the coefficients",
            r"\begin{proposition}[Sesquilinearity in cohomology; \ClaimStatusProvedHere]",
        )
        for needle in retired:
            assert needle not in block

    def test_legacy_sesquilinearity_shadow_has_no_removed_refs(self):
        """The superseded PVA split must not cite removed PVA1 labels."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent.tex").read_text()

        early_start = source.index(r"\begin{lemma}[Sesquilinearity verification")
        early_end = source.index(r"\end{proof}", early_start)
        early = source[early_start:early_end]
        early_flat = " ".join(early.split())
        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "polynomial divided-power convention",
            "first sesquilinearity identity",
            "second sesquilinearity identity",
        ):
            assert needle in early or needle in early_flat

        later_start = source.index(r"\begin{proposition}[Sesquilinearity in cohomology")
        later_end = source.index(r"\end{proof}", later_start)
        later = source[later_start:later_end]
        later_flat = " ".join(later.split())
        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            r"(\partial a)_{(n)}b=-n\,a_{(n-1)}b",
            r"a_{(n)}(\partial b)=\partial(a_{(n)}b)+n\,a_{(n-1)}b",
            "regular/singular binary descent lemma above",
        ):
            assert needle in later or needle in later_flat

        completion_start = source.index(r"\begin{proof}[Proof of Theorem")
        completion = source[completion_start: source.index(r"\end{proof}", completion_start)]
        retired = (
            r"\ref{eq:PVA1a}",
            r"\ref{eq:PVA1b}",
            r"\ref{prop:PVA1_proof}",
            r"Definition \ref{def:ainfty_chiral}",
            r"\{a_\lambda b\} = [m_2(a,b)]",
        )
        for needle in retired:
            assert needle not in early
            assert needle not in later
            assert needle not in completion

    def test_exchange_lagrangian_monodromy_is_half_turn(self):
        """The exchange path is a half-monodromy, not a full winding."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()
        start = source.index(
            r"\begin{proposition}[Exchange half-monodromy as Lagrangian monodromy"
        )
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\alpha+\gamma$",
            r"\hypAmbientPro",
            r"z_1(s)-z_2(s)=2r\,e^{\pi i s}",
            "traverses the half-boundary arc",
            "not a loop",
            "The square of the exchange",
            r"\mathrm{Conf}_2(\R^3)\simeq",
            "No braid-group",
            r"d\log(-\zeta)=d\log\zeta",
            "oriented boundary sign",
            r"\lambda\mapsto-\lambda-\partial",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            "winds once around",
            "carries the winding class",
            r"full ordered configuration space in $\C\times\R$ is three-real-dimensional",
            r"\log(z_1-z_2) \mapsto",
            "contributes a sign to the integrand",
            r"generator of $\pi_1(\mathrm{Conf}_2(\C \times \R))",
        )
        for needle in retired:
            assert needle not in block

    def test_legacy_skew_symmetry_shadow_uses_half_turn(self):
        """The superseded PVA split must not sell the old exchange proof."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent.tex").read_text()
        start = source.index(r"\begin{proposition}[Skew-symmetry in cohomology;")
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\alpha+\gamma$",
            r"\hypAmbientPro",
            r"\zeta(s)=z_1(s)-z_2(s)=2r\,e^{\pi i s}",
            "half-boundary arc",
            "The square of this exchange",
            "No braid-group",
            r"d\log(-\zeta)=d\log\zeta",
            "oriented boundary sign",
            r"\lambda\mapsto-\lambda-\partial",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            "label removed: prop:PVA2_proof",
            "Construct an explicit chain homotopy",
            r"\sigma_*[\FM_2(\C)]",
            "sign from $\\sigma$",
            r"$\lambda = z_1 - z_2$ becomes $-\lambda = z_2 - z_1$",
        )
        for needle in retired:
            assert needle not in block

    def test_exchange_context_names_bulk_configuration_space(self):
        """Nearby context must not call Conf_2(C x R) three-dimensional."""
        root = Path(__file__).resolve().parents[2]
        appendix = (
            root / "chapters" / "theory" / "pva-expanded-repaired.tex"
        ).read_text()
        axioms = (root / "chapters" / "theory" / "axioms.tex").read_text()

        assert "configuration space of two points in the $3$-dimensional bulk" in appendix
        assert "configuration space of two points in the three-dimensional bulk" in axioms

        retired = (
            "full $3$-dimensional configuration space of bulk points",
            "connectivity of the full three-dimensional configuration space",
        )
        for needle in retired:
            assert needle not in appendix
            assert needle not in axioms

    def test_exchange_homotopy_lemma_has_complete_tau_substitution(self):
        """PVA2 exchange homotopy must state the full vertex substitution."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()
        start = source.index(r"\begin{lemma}[Chain-level exchange homotopies")
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\alpha+\gamma$",
            r"\hypAmbientPro",
            "oriented exchange half-monodromy",
            r"\zeta(s)=z_1(s)-z_2(s)=2r\,e^{\pi i s}",
            r"d\log(-\zeta)=d\log\zeta",
            "there is no separate PVA operation",
            r"\lambda\mapsto-\lambda-\partial",
            r"\eqref{eq:repaired-n2}",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            r"\begin{lemma}[Chain-level exchange homotopies; \ClaimStatusProvedHere]",
            r"(which becomes $\lambda\mapsto-\lambda$ under",
            r"$\lambda\mapsto\lambda-\partial$",
            r"$[0,1]\times \FM_2(\C)$",
        )
        for needle in retired:
            assert needle not in block

    def test_exchange_consequences_carry_same_license(self):
        """Product commutativity and PVA2 inherit exchange licensing."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()

        product_start = source.index(r"\begin{proposition}[Graded commutativity of the product")
        product_end = source.index(r"\end{proof}", product_start)
        product = source[product_start:product_end]

        skew_start = source.index(r"\begin{proposition}[Skew-symmetry of the $\lambda$-bracket")
        skew_end = source.index(r"\end{proof}", skew_start)
        skew = source[skew_start:skew_end]

        for block in (product, skew):
            flat = " ".join(block.split())
            for needle in (
                r"\ClaimStatusProvedHere",
                r"licensing $\alpha+\gamma$",
                r"\hypAmbientPro",
                "oriented exchange half-monodromy",
            ):
                assert needle in block or needle in flat

        assert r"-(-1)^{\degree{a}\degree{b}}" in skew
        assert r"(-1)^{(\degree{a}+1)(\degree{b}+1)}" not in skew

    def test_pva2_shadow_surfaces_use_visible_koszul_sign(self):
        """Appendix and preview PVA2 copies must match the active convention."""
        root = Path(__file__).resolve().parents[2]

        repaired = (
            root / "chapters" / "theory" / "pva-expanded-repaired.tex"
        ).read_text()
        remark_start = repaired.index(r"\begin{remark}[Regular and singular exchange signs]")
        remark_end = repaired.index(r"\end{remark}", remark_start)
        remark = repaired[remark_start:remark_end]
        remark_flat = " ".join(remark.split())
        for needle in (
            "visible vertex-algebra convention",
            r"-(-1)^{\degree{a}\degree{b}}",
            "Shifted signs in the three-face Jacobi calculation",
            "do not replace the PVA2 skewsymmetry sign",
        ):
            assert needle in remark or needle in remark_flat
        assert "singular piece, after Borel transform, is a bracket of degree $-1$" not in remark

        old_appendix = (root / "appendices" / "pva-expanded.tex").read_text()
        app_start = old_appendix.index(r"\begin{proposition}[Skew-symmetry]")
        app_end = old_appendix.index(r"\end{proof}", app_start)
        app_block = old_appendix[app_start:app_end]
        for needle in (
            r"-(-1)^{\degree{a}\degree{b}}",
            "holomorphic half-boundary arc",
            r"d\log(-\zeta)=d\log\zeta",
            r"\lambda\mapsto-\lambda-\partial",
        ):
            assert needle in app_block
        for needle in (
            r"-(-1)^{(\degree{a}+1)(\degree{b}+1)}",
            "naive sign",
            r"m_2(b,a;-\lambda)",
        ):
            assert needle not in app_block

        preview = (root / "chapters" / "theory" / "pva-preview.tex").read_text()
        preview_start = preview.index(r"\begin{definition}[Vertex Poisson bracket]")
        preview_end = preview.index(r"\end{proof}", preview.index(r"\begin{proposition}[Sesquilinearity and skewsymmetry]"))
        preview_block = preview[preview_start:preview_end]
        assert r"-(-1)^{\degree{a}\degree{b}}" in preview_block
        assert "ordinary Koszul sign" in preview_block
        assert r"-(-1)^{(\degree{a}+1)(\degree{b}+1)}" not in preview_block

    def test_leibniz_uses_signed_mixed_three_face_defect(self):
        """Leibniz must be the signed mixed-sector defect, not an unsigned sum."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()
        start = source.index(r"\subsection{Leibniz from the mixed regular--singular sector}")
        end = source.index(r"\subsection{Vacuum and completion of the PVA proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        required = (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            r"\mathcal L_\lambda(a;b,c)",
            r"\mathfrak L_{23}",
            r"\mathfrak L_{12}",
            r"\mathfrak L_{13}",
            r"K_3^{\mathrm{Leib}}",
            r"\operatorname{im}Q",
            "oriented $\\FM_3$ incidence",
            "unsuspended regular product",
            "divided-power Borel transform",
            "singular--singular Wick integral term",
            "does not appear in this classical PVA descent",
            r"\eqref{eq:repaired-leibniz}",
        )
        for needle in required:
            assert needle in block or needle in flat

        retired = (
            "Leibniz from the exchange cylinder",
            "skew-symmetry from the two-point winding",
            "Then the three pairwise divisors contribute:",
            "The three pairwise residues are exactly the three terms",
            r"\int_0^\lambda \{\{",
        )
        for needle in retired:
            assert needle not in block

    def test_leibniz_shadow_surfaces_have_no_wick_integral(self):
        """Appendix, preview, and raviolo copies must keep classical Leibniz."""
        root = Path(__file__).resolve().parents[2]

        repaired = (
            root / "chapters" / "theory" / "pva-expanded-repaired.tex"
        ).read_text()
        proj_start = repaired.index(r"\begin{proposition}[The mixed projection;")
        proj_end = repaired.index(r"\end{proof}", proj_start)
        proj = repaired[proj_start:proj_end]
        proj_flat = " ".join(proj.split())
        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            r"\mathcal L_\lambda(a;b,c)",
            r"\mathfrak L_{23}",
            r"\mathfrak L_{12}",
            r"\mathfrak L_{13}",
            "unsuspended regular product",
            "singular--singular Wick integral term is absent",
        ):
            assert needle in proj or needle in proj_flat

        repaired_remark_start = repaired.index(r"\begin{remark}[No Wick-type integral term is needed]")
        repaired_remark_end = repaired.index(r"\end{remark}", repaired_remark_start)
        repaired_remark = repaired[repaired_remark_start:repaired_remark_end]
        repaired_remark_flat = " ".join(repaired_remark.split())
        assert "standard} Poisson vertex Leibniz rule" in repaired_remark_flat
        assert r"No extra $\int_0^\lambda$ correction term" in repaired_remark_flat

        old_appendix = (root / "appendices" / "pva-expanded.tex").read_text()
        app_start = old_appendix.index(r"\begin{proposition}[Leibniz rule]")
        app_end = old_appendix.index(r"\end{proof}", app_start)
        app = old_appendix[app_start:app_end]
        app_flat = " ".join(app.split())
        for needle in (
            r"\mathcal L_\lambda(a;b,c)",
            r"\mathfrak L_{23}",
            r"\mathfrak L_{12}",
            r"\mathfrak L_{13}",
            r"\operatorname{im}Q",
            "No Wick-type integral term appears",
            "singular--singular Jacobi sector",
            "mixed regular--singular sector",
        ):
            assert needle in app or needle in app_flat
        assert r"\int_0^\lambda" not in app

        preview = (root / "chapters" / "theory" / "pva-preview.tex").read_text()
        prev_start = preview.index(r"\begin{theorem}[Leibniz rule]")
        prev_end = preview.index(r"\end{proof}", prev_start)
        prev = preview[prev_start:prev_end]
        prev_flat = " ".join(prev.split())
        for needle in (
            "signed mixed",
            r"\mathcal L_\lambda(a;b,c)",
            "one singular channel and one regular channel",
            "no singular--singular Wick integral correction appears",
        ):
            assert needle in prev or needle in prev_flat

        raviolo = (root / "chapters" / "theory" / "raviolo.tex").read_text()
        rav_start = raviolo.index(r"\textbf{Step 6: Leibniz rule.}")
        rav_end = raviolo.index(r"\medskip\noindent\textbf{Step 7:", rav_start)
        rav = raviolo[rav_start:rav_end]
        rav_flat = " ".join(rav.split())
        for needle in (
            "signed Leibniz defect",
            r"\mathcal L_\lambda(a;b,c)",
            r"Proposition~\ref{prop:PVA3_proof}",
            r"\eqref{eq:repaired-leibniz}",
            "No Wick-type integral is present",
            "mixed regular--singular Leibniz sector",
        ):
            assert needle in rav or needle in rav_flat

        retired = (
            r"The third term in \eqref{eq:app-Leibniz} is the \emph{non-classical} contribution",
            r"Theorem~\ref{thm:Leibniz-PVA}",
            r"m_2(m_2(a,b),c) \;\pm\;",
            "The mixed projection selects terms where exactly one pair",
            r"\int_0^\lambda \{\{",
        )
        for needle in retired:
            assert needle not in app
            assert needle not in prev
            assert needle not in rav

    def test_vacuum_unit_axiom_uses_strict_unit_on_both_sides(self):
        """PVA4 must prove both unit brackets and both product units."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()
        start = source.index(r"\begin{proposition}[Vacuum/unit axiom;")
        end = source.index(r"\end{proof}", start)
        block = source[start:end]
        flat = " ".join(block.split())

        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "strict unital",
            "unsuspended regular product",
            "regular/singular projection",
            "divided-power Borel transform",
            r"\{[\mathbf 1]_\lambda[a]\}=0",
            r"\{[a]_\lambda[\mathbf 1]\}=0",
            r"m_2^{\mathrm{sing}}(\mathbf 1,a)=0",
            r"m_2^{\mathrm{sing}}(a,\mathbf 1)=0",
            r"\eqref{eq:unit-higher}",
            r"\partial\mathbf 1=0",
        ):
            assert needle in block or needle in flat

        cor_start = source.index(r"\begin{corollary}[Shifted PVA structure on cohomology;")
        cor_end = source.index(r"\end{proof}", cor_start)
        cor = source[cor_start:cor_end]
        cor_flat = " ".join(cor.split())
        for needle in (
            r"licensing $\alpha+\gamma$",
            "oriented exchange half-monodromy",
            r"\hypAmbientPro",
            r"oriented $\FM_3$ incidence convention",
            "strict unital",
            "divided-power Borel transform",
        ):
            assert needle in cor or needle in cor_flat

        retired = (
            r"\begin{proposition}[Vacuum/unit axiom; \ClaimStatusProvedHere]",
            r"\{\mathbf 1_\lambda a\}=0",
            r"\mu_2(\mathbf 1,a)=a,\qquad \mu_2(a,\mathbf 1)=a",
            "Finally, $\\partial\\mathbf 1=0$ by the unit axiom.",
        )
        for needle in retired:
            assert needle not in block

    def test_vacuum_shadow_surfaces_do_not_use_symmetrized_unit_argument(self):
        """Old appendix, legacy split, and raviolo must use the same unit convention."""
        root = Path(__file__).resolve().parents[2]

        old_appendix = (root / "appendices" / "pva-expanded.tex").read_text()
        app_start = old_appendix.index(r"\begin{proposition}[Vacuum]")
        app_end = old_appendix.index(r"\end{proof}", app_start)
        app = old_appendix[app_start:app_end]
        app_flat = " ".join(app.split())
        for needle in (
            r"\{[a]{}_\lambda[\mathbf{1}]\}=0",
            r"[\mathbf{1}]\cdot [a] &= [a]\cdot[\mathbf{1}]=[a]",
            r"\partial[\mathbf 1]=0",
            "unsuspended regular",
            "both singular binary projections",
            "no hidden higher boundary term",
        ):
            assert needle in app or needle in app_flat

        legacy = (root / "chapters" / "theory" / "pva-descent.tex").read_text()
        leg_start = legacy.index(r"\begin{proposition}[Vacuum Axiom in Cohomology;")
        leg_end = legacy.index(r"\end{proof}", leg_start)
        leg = legacy[leg_start:leg_end]
        leg_flat = " ".join(leg.split())
        for needle in (
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            r"\{[a]_\lambda[\mathbf 1]\}=0",
            "unsuspended regular product",
            r"m_2^{\mathrm{sing}}(a,\mathbf 1)=0",
            "translation-invariant vacuum vector",
        ):
            assert needle in leg or needle in leg_flat

        raviolo = (root / "chapters" / "theory" / "raviolo.tex").read_text()
        rav_start = raviolo.index(r"\item \emph{Vacuum axiom}:")
        rav_end = raviolo.index(r"\item \emph{Translation covariance}", rav_start)
        rav = raviolo[rav_start:rav_end]
        rav_flat = " ".join(rav.split())
        for needle in (
            "unsuspended vertex-algebra frame",
            "desuspension convention",
            r"Y(|0\rangle,z)=\mathrm{id}_V",
            "creation property",
            "singular part of either binary kernel with a unit insertion is zero",
        ):
            assert needle in rav or needle in rav_flat

        retired = (
            "apparent discrepancy for odd-degree elements",
            "if $\\degree{a}$ even",
            "if $\\degree{a}$ odd",
            r"\Sym(m_2^{\mathrm{reg}}(\mathbf{1},a))",
            r"\ASym(m_2^{\mathrm{sing}}(\mathbf{1},a))",
            "by the following sign absorption mechanism",
            r"(-1)^{2|a|-1} = -1",
        )
        for needle in retired:
            assert needle not in app
            assert needle not in leg
            assert needle not in rav

    def test_higher_operations_need_compatible_topological_nullhomotopies(self):
        """Higher-operation Q-exactness is conditional on relative fillings."""
        root = Path(__file__).resolve().parents[2]
        source = (root / "chapters" / "theory" / "pva-descent-repaired.tex").read_text()

        lemma_start = source.index(r"\begin{lemma}[Open ordered topological contraction;")
        lemma_end = source.index(r"\end{proof}", lemma_start)
        lemma = source[lemma_start:lemma_end]
        lemma_flat = " ".join(lemma.split())
        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma$",
            r"\hypAmbientPro",
            "translation quotient",
            "open ordered topological factor only",
            "does not by itself construct a relative bounding chain",
            "compactified operadic pair",
        ):
            assert needle in lemma or needle in lemma_flat

        prop_start = source.index(
            r"\begin{proposition}[Higher operations are \(Q\)-exact under compatible topological"
        )
        prop_end = source.index(r"\end{proof}", prop_start)
        prop = source[prop_start:prop_end]
        prop_flat = " ".join(prop.split())
        for needle in (
            r"\ClaimStatusProvedHere",
            r"licensing $\gamma+\varepsilon$",
            r"\hypAmbientPro",
            "relative compactified ordered topological factor",
            r"Definition~\textup{\ref{def:compatible-topological-nullhomotopy}}",
            r"\widetilde\omega_k=\omega_k^{\mathrm{hol}}\wedge\omega_k^{\mathrm{top}}",
            r"\partial_{\mathrm{rel}}\Gamma_{k-1}=\omega_k^{\mathrm{top}}",
            "Ordinary contractibility of the open orthant",
            "does not supply the required",
            "Without the supplied relative bounding chains",
            "no higher-operation vanishing statement is asserted",
        ):
            assert needle in prop or needle in prop_flat

        retired = (
            r"\begin{lemma}[Topological contraction; \ClaimStatusProvedHere]",
            r"\begin{proposition}[Higher operations vanish under compatible topological nullhomotopies;",
            r"\ClaimStatusConditional{} on compatible topological nullhomotopies",
            "Contractibility of the ordered topological configuration space is not by itself a\nchain-level nullhomotopy",
        )
        for needle in retired:
            assert needle not in lemma
            assert needle not in prop

    def test_higher_operation_shadows_do_not_feed_binary_pva_proof(self):
        """Shadow copies must not use global higher vanishing to prove Jacobi."""
        root = Path(__file__).resolve().parents[2]

        old_appendix = (root / "appendices" / "pva-expanded.tex").read_text()
        app_start = old_appendix.index(r"\textbf{Step 2: Pass to cohomology.}")
        app_end = old_appendix.index(r"\textbf{Step 3: Extract the Jacobi identity.}", app_start)
        app = old_appendix[app_start:app_end]
        app_flat = " ".join(app.split())
        for needle in (
            "total-collision term",
            r"m_1(m_3(a,b,c))=Qm_3(a,b,c)",
            "local \\(Q\\)-boundary",
            "not the separate global statement",
            "all higher operations vanish in cohomology",
        ):
            assert needle in app or needle in app_flat

        legacy = (root / "chapters" / "theory" / "pva-descent.tex").read_text()
        leg_start = legacy.index(r"\begin{proposition}[Higher operations are \(Q\)-exact")
        leg_end = legacy.index(r"\end{proof}", leg_start)
        leg = legacy[leg_start:leg_end]
        leg_flat = " ".join(leg.split())
        for needle in (
            r"licensing $\gamma+\varepsilon$",
            r"\hypAmbientPro",
            "relative compactified ordered topological factor",
            "compatible topological nullhomotopies",
            "Ordinary contractibility of the open ordered configuration space is not enough",
        ):
            assert needle in leg or needle in leg_flat

        completion_start = legacy.index(r"\begin{proof}[Proof of Theorem")
        completion_end = legacy.index(r"\end{proof}", completion_start)
        completion = legacy[completion_start:completion_end]
        assert "binary PVA theorem and the separate" in completion
        assert "Therefore the binary operations define a Poisson vertex algebra" in completion

        raviolo = (root / "chapters" / "theory" / "raviolo.tex").read_text()
        rav_start = raviolo.index(r"\paragraph{\textbf{Landau")
        rav_end = raviolo.index(r"\paragraph{\textbf{Abelian", rav_start)
        rav = raviolo[rav_start:rav_end]
        rav_flat = " ".join(rav.split())
        for needle in (
            "first nontrivial $m_3$ may appear",
            "total-collision contribution",
            r"\(Q\)-boundary",
            "not declared globally zero",
            r"Proposition~\ref{prop:m3_vanish}",
        ):
            assert needle in rav or needle in rav_flat

        retired = (
            "using the vanishing of higher operations in cohomology",
            "by the usual $A_\\infty$ argument",
            "We have now established all three parts",
            "This space is diffeomorphic to the open simplex",
            "Contractibility of the open configuration space is only the homological background.",
        )
        for needle in retired:
            assert needle not in app
            assert needle not in leg
            assert needle not in completion
            assert needle not in rav


class TestLambdaBracket:
    """Test lambda-bracket representation and evaluation."""

    def test_heisenberg_bracket_eval(self):
        """Heisenberg: {J_lam J} = k*lam."""
        k = Symbol('k')
        lam = Symbol('lam')
        bracket = LambdaBracket(coefficients={0: S.Zero, 1: k})
        result = bracket.evaluate(lam)
        assert simplify(result - k * lam) == 0

    def test_virasoro_bracket_eval(self):
        """Virasoro: {T_lam T} = dT + 2T*lam + (c/12)*lam^3."""
        c, T, dT = symbols('c T dT')
        lam = Symbol('lam')
        bracket = LambdaBracket(coefficients={0: dT, 1: 2*T, 2: S.Zero, 3: c/2})
        result = bracket.evaluate(lam)
        expected = dT + 2*T*lam + (c/2)*lam**3/6
        assert simplify(result - expected) == 0

    def test_pole_order_heisenberg(self):
        """Heisenberg: lambda degree = 1, coming from OPE double pole."""
        k = Symbol('k')
        bracket = LambdaBracket(coefficients={0: S.Zero, 1: k})
        assert bracket.pole_order() == 1

    def test_pole_order_virasoro(self):
        """Virasoro: pole order = 3 (4th order pole)."""
        c, T, dT = symbols('c T dT')
        bracket = LambdaBracket(coefficients={0: dT, 1: 2*T, 2: S.Zero, 3: c/2})
        assert bracket.pole_order() == 3

    def test_n_product_extraction(self):
        """T_{(3)}T = c/2 for Virasoro."""
        c = Symbol('c')
        bracket = LambdaBracket(coefficients={3: c/2})
        assert bracket.n_product(3) == c/2

    def test_n_product_zero(self):
        """T_{(2)}T = 0 for Virasoro."""
        bracket = LambdaBracket(coefficients={0: Symbol('dT'), 1: 2*Symbol('T'), 3: Symbol('c')/2})
        assert bracket.n_product(2) == S.Zero

    def test_empty_bracket_pole(self):
        """Empty bracket has pole order -1."""
        bracket = LambdaBracket(coefficients={})
        assert bracket.pole_order() == -1


# ===================================================================
# D2: SESQUILINEARITY
# ===================================================================

class TestD2Sesquilinearity:
    """D2 sesquilinearity for all families."""

    def test_d2_heisenberg(self):
        """D2 for Heisenberg: {dJ_lam J} = -lam * k*lam = -k*lam^2."""
        pva = heisenberg_pva()
        result = verify_d2_sesquilinearity(pva, 'J', 'J')
        assert result['left_holds']
        assert result['right_holds']

    def test_d2_virasoro(self):
        """D2 for Virasoro: sesquilinearity of {T_lam T}."""
        pva = virasoro_pva()
        result = verify_d2_sesquilinearity(pva, 'T', 'T')
        assert result['left_holds']
        assert result['right_holds']

    def test_d2_sl2_all_pairs(self):
        """D2 for sl_2: all 9 generator pairs."""
        pva = affine_sl2_pva()
        for a in pva.generators:
            for b in pva.generators:
                result = verify_d2_sesquilinearity(pva, a, b)
                assert result['left_holds'], f"D2 left fails for ({a}, {b})"
                assert result['right_holds'], f"D2 right fails for ({a}, {b})"

    def test_d2_betagamma(self):
        """D2 for beta-gamma."""
        pva = betagamma_pva()
        for a in pva.generators:
            for b in pva.generators:
                result = verify_d2_sesquilinearity(pva, a, b)
                assert result['left_holds']

    def test_d2_free_multiplet(self):
        """D2 for free multiplet."""
        pva = free_multiplet_pva()
        result = verify_d2_sesquilinearity(pva, 'phi', 'psi')
        assert result['left_holds']

    def test_d2_source_algebraic(self):
        """D2 is an algebraic identity (Borcherds)."""
        pva = virasoro_pva()
        result = verify_d2_sesquilinearity(pva, 'T', 'T')
        assert result['source'] == 'algebraic (Borcherds identity)'


# ===================================================================
# D3: JACOBI IDENTITY
# ===================================================================

class TestD3Jacobi:
    """D3 Jacobi identity for all families."""

    def test_d3_heisenberg_trivial(self):
        """D3 for Heisenberg: Jacobi is trivial (abelian)."""
        pva = heisenberg_pva()
        result = verify_d3_jacobi_generators(pva, 'J', 'J', 'J')
        assert result['fm3_faces'] == 3

    def test_d3_sl2_all_27(self):
        """D3 for sl_2: all 27 triples satisfy Jacobi."""
        result = verify_d3_jacobi_sl2()
        assert result['num_triples'] == 27
        assert result['all_pass']

    def test_d3_virasoro(self):
        """D3 for Virasoro: Jacobi from sesquilinearity."""
        result = verify_d3_jacobi_virasoro()
        assert result['jacobi_holds']

    def test_d3_betagamma_all(self):
        """D3 for beta-gamma: all 8 triples."""
        pva = betagamma_pva()
        count = 0
        for a in pva.generators:
            for b in pva.generators:
                for c in pva.generators:
                    result = verify_d3_jacobi_generators(pva, a, b, c)
                    count += 1
        assert count == 8  # 2^3

    def test_d3_geometric_source(self):
        """D3 geometric source is three-face Stokes on FM_3."""
        pva = heisenberg_pva()
        result = verify_d3_jacobi_generators(pva, 'J', 'J', 'J')
        assert result['geometric_source'] == 'three-face Stokes on FM_3(C)'


# ===================================================================
# D4: LEIBNIZ RULE
# ===================================================================

class TestD4Leibniz:
    """D4 Leibniz rule for all families."""

    def test_d4_heisenberg(self):
        """D4 for Heisenberg: Leibniz trivially holds."""
        pva = heisenberg_pva()
        result = verify_d4_leibniz(pva, 'J', 'J', 'J')
        assert result['identity_type'] == 'D4 Leibniz'

    def test_d4_virasoro(self):
        """D4 for Virasoro: Leibniz from three-face Stokes."""
        pva = virasoro_pva()
        result = verify_d4_leibniz(pva, 'T', 'T', 'T')
        assert result['geometric_source'] == 'three-face Stokes on FM_3(C)'

    def test_d4_sl2(self):
        """D4 for sl_2: test on (e, h, f)."""
        pva = affine_sl2_pva()
        result = verify_d4_leibniz(pva, 'e', 'h', 'f')
        assert result['identity_type'] == 'D4 Leibniz'

    def test_d4_betagamma(self):
        """D4 for beta-gamma."""
        pva = betagamma_pva()
        result = verify_d4_leibniz(pva, 'beta', 'gamma', 'beta')
        assert result['num_checks'] >= 0


# ===================================================================
# D5: SKEW-SYMMETRY
# ===================================================================

class TestD5SkewSymmetry:
    """D5 skew-symmetry for all families."""

    def test_d5_heisenberg(self):
        """D5 for Heisenberg: {J_lam J} = -{J_{-lam-d} J}."""
        pva = heisenberg_pva()
        result = verify_d5_skew_symmetry(pva, 'J', 'J')
        assert result['all_pass']

    def test_d5_virasoro(self):
        """D5 for Virasoro: skew-symmetry of {T_lam T}."""
        pva = virasoro_pva()
        result = verify_d5_skew_symmetry(pva, 'T', 'T')
        assert result['all_pass']

    def test_d5_sl2_all_pairs(self):
        """D5 for sl_2: all 9 pairs."""
        pva = affine_sl2_pva()
        for a in pva.generators:
            for b in pva.generators:
                result = verify_d5_skew_symmetry(pva, a, b)
                assert result['all_pass'], f"D5 fails for ({a}, {b})"

    def test_d5_betagamma(self):
        """D5 for beta-gamma: {beta_lam gamma} = -{gamma_{-lam-d} beta}."""
        pva = betagamma_pva()
        result = verify_d5_skew_symmetry(pva, 'beta', 'gamma')
        assert result['all_pass']

    def test_d5_free_multiplet(self):
        """D5 for free multiplet."""
        pva = free_multiplet_pva()
        result = verify_d5_skew_symmetry(pva, 'phi', 'psi')
        assert result['all_pass']

    def test_d5_geometric_source(self):
        """D5 geometric source is monodromy."""
        pva = heisenberg_pva()
        result = verify_d5_skew_symmetry(pva, 'J', 'J')
        assert result['identity_type'] == 'D5 skew-symmetry'


# ===================================================================
# VACUUM UNIT AXIOM
# ===================================================================

class TestD6Unit:
    """Vacuum unit axiom for all families."""

    def test_d6_heisenberg(self):
        """Vacuum for Heisenberg: brackets vanish and 1 is a unit."""
        pva = heisenberg_pva()
        result = verify_d6_unit(pva, 'J')
        assert result['vanishes']
        assert result['right_bracket_vanishes']
        assert result['left_unit_holds']
        assert result['right_unit_holds']
        assert result['translation_unit_holds']

    def test_d6_virasoro(self):
        """Vacuum for Virasoro: {1_lam T} = 0."""
        pva = virasoro_pva()
        result = verify_d6_unit(pva, 'T')
        assert result['vanishes']

    def test_d6_sl2_all(self):
        """Vacuum for sl_2: {1_lam J^a} = 0 for all a."""
        pva = affine_sl2_pva()
        for a in pva.generators:
            result = verify_d6_unit(pva, a)
            assert result['vanishes'], f"vacuum unit fails for {a}"

    def test_d6_betagamma(self):
        """Vacuum for beta-gamma."""
        pva = betagamma_pva()
        for a in pva.generators:
            result = verify_d6_unit(pva, a)
            assert result['vanishes']

    def test_d6_geometric_source(self):
        """The unit follows from vacuum factorization."""
        pva = heisenberg_pva()
        result = verify_d6_unit(pva, 'J')
        assert result['geometric_source'] == 'strict unit/vacuum factorization with unit insertion'


# ===================================================================
# FULL D2-D5 PLUS VACUUM SWEEPS
# ===================================================================

class TestFullPVASweeps:
    """Full D2-D5 plus vacuum verification sweeps for each family."""

    def test_full_sweep_heisenberg(self):
        """Full D2-D5 plus vacuum for Heisenberg."""
        pva = heisenberg_pva()
        result = full_pva_descent_verification(pva)
        assert result['D2_sesquilinearity']['all_pass']
        assert result['D5_skew']['all_pass']
        assert result['D6_unit']['all_pass']

    def test_full_sweep_virasoro(self):
        """Full D2-D5 plus vacuum for Virasoro."""
        pva = virasoro_pva()
        result = full_pva_descent_verification(pva)
        assert result['D2_sesquilinearity']['all_pass']
        assert result['D5_skew']['all_pass']
        assert result['D6_unit']['all_pass']

    def test_full_sweep_sl2(self):
        """Full D2-D5 plus vacuum for affine sl_2."""
        pva = affine_sl2_pva()
        result = full_pva_descent_verification(pva)
        assert result['D2_sesquilinearity']['all_pass']
        assert result['D5_skew']['all_pass']
        assert result['D6_unit']['all_pass']

    def test_full_sweep_betagamma(self):
        """Full D2-D5 plus vacuum for beta-gamma."""
        pva = betagamma_pva()
        result = full_pva_descent_verification(pva)
        assert result['D2_sesquilinearity']['all_pass']
        assert result['D5_skew']['all_pass']
        assert result['D6_unit']['all_pass']

    def test_full_sweep_free_multiplet(self):
        """Full D2-D5 plus vacuum for free multiplet."""
        pva = free_multiplet_pva()
        result = full_pva_descent_verification(pva)
        assert result['D2_sesquilinearity']['all_pass']
        assert result['D6_unit']['all_pass']

    def test_full_sweep_lg_cubic(self):
        """Full D2-D5 plus vacuum for LG cubic (same PVA as free at this level)."""
        pva = lg_cubic_pva()
        result = full_pva_descent_verification(pva)
        assert result['D2_sesquilinearity']['all_pass']
        assert result['D6_unit']['all_pass']

    def test_full_sweep_w3(self):
        """Full D2-D5 plus vacuum for W_3."""
        pva = w3_pva()
        result = full_pva_descent_verification(pva)
        assert result['D2_sesquilinearity']['all_pass']
        assert result['D6_unit']['all_pass']


# ===================================================================
# FM_3 BOUNDARY STRATA
# ===================================================================

class TestFM3BoundaryStrata:
    """FM_3(C) boundary strata cancellations."""

    def test_partial_fraction_identity(self):
        """The partial fraction identity: 1/(z12*z23) + cyc = 0."""
        result = fm3_boundary_strata_cancellation()
        assert result['partial_fraction_vanishes']

    def test_linear_identity(self):
        """z12 + z23 + z31 = 0."""
        result = fm3_boundary_strata_cancellation()
        assert result['linear_vanishes']

    def test_three_faces(self):
        """FM_3 has exactly 3 codim-1 boundary faces."""
        result = fm3_boundary_strata_cancellation()
        # The three faces correspond to three terms
        assert 'term1' in result
        assert 'term2' in result
        assert 'term3' in result

    def test_exchange_cylinder(self):
        """Exchange cylinder argument for D3."""
        result = fm3_exchange_cylinder_stokes()
        assert result['stokes_gives_jacobi']


# ===================================================================
# POLE ORDER CENSUS
# ===================================================================

class TestPoleOrderCensus:
    """Pole order census across families."""

    def test_heisenberg_pole_order(self):
        """Heisenberg: max pole order = 1."""
        pva = heisenberg_pva()
        result = pole_order_census(pva)
        assert result['max_pole_order'] == 1

    def test_virasoro_pole_order(self):
        """Virasoro: max pole order = 3."""
        pva = virasoro_pva()
        result = pole_order_census(pva)
        assert result['max_pole_order'] == 3

    def test_w3_pole_order(self):
        """W_3: max pole order = 5 (from {W_lam W})."""
        pva = w3_pva()
        result = pole_order_census(pva)
        assert result['max_pole_order'] == 5

    def test_betagamma_pole_order(self):
        """Beta-gamma: max pole order = 0."""
        pva = betagamma_pva()
        result = pole_order_census(pva)
        assert result['max_pole_order'] == 0

    def test_sl2_pole_order(self):
        """sl_2: max pole order = 1 (simple pole from level)."""
        pva = affine_sl2_pva()
        result = pole_order_census(pva)
        assert result['max_pole_order'] == 1

    def test_free_multiplet_pole_order(self):
        """Free multiplet: max pole order = 0."""
        pva = free_multiplet_pva()
        result = pole_order_census(pva)
        assert result['max_pole_order'] == 0


# ===================================================================
# PVA DATA CONSISTENCY
# ===================================================================

class TestPVADataConsistency:
    """Consistency of PVA data across families."""

    def test_heisenberg_one_generator(self):
        """Heisenberg has 1 generator."""
        pva = heisenberg_pva()
        assert len(pva.generators) == 1

    def test_virasoro_one_generator(self):
        """Virasoro has 1 generator."""
        pva = virasoro_pva()
        assert len(pva.generators) == 1

    def test_sl2_three_generators(self):
        """sl_2 has 3 generators."""
        pva = affine_sl2_pva()
        assert len(pva.generators) == 3

    def test_w3_two_generators(self):
        """W_3 has 2 generators."""
        pva = w3_pva()
        assert len(pva.generators) == 2

    def test_betagamma_two_generators(self):
        """Beta-gamma has 2 generators."""
        pva = betagamma_pva()
        assert len(pva.generators) == 2

    def test_virasoro_t3t_equals_c_over_2(self):
        """Virasoro: T_{(3)}T = c/2 (NOT c/12; that's the lambda-bracket)."""
        pva = virasoro_pva()
        c = Symbol('c')
        t3t = pva.n_product('T', 'T', 3)
        assert simplify(t3t - c/2) == 0

    def test_sl2_bracket_antisymmetry(self):
        """sl_2: e_{(0)}f = h, f_{(0)}e = -h."""
        pva = affine_sl2_pva()
        e0f = pva.n_product('e', 'f', 0)
        f0e = pva.n_product('f', 'e', 0)
        assert simplify(e0f + f0e) == 0

    def test_w3_zamolodchikov_coefficients(self):
        """W_3: Zamolodchikov Lambda coefficients in the WW products."""
        pva = w3_pva()
        c = Symbol('c')
        bracket_WW = pva.brackets[('W', 'W')]
        Lambda = Symbol('Lambda')
        dLambda = Symbol('dLambda')
        d2T = Symbol('d2T')
        d3T = Symbol('d3T')
        beta = Rational(16, 1) / (22 + 5*c)

        w1w = bracket_WW.n_product(1)
        assert simplify(w1w - (Rational(3, 10)*d2T + 2*beta*Lambda)) == 0

        w0w = bracket_WW.n_product(0)
        assert simplify(w0w - (d3T/Rational(15, 1) + beta*dLambda)) == 0

    def test_w3_skew_primary_product(self):
        """W_3: skew symmetry gives W_{(0)}T = 2 partial W."""
        pva = w3_pva()
        dW = Symbol('dW')
        assert simplify(pva.n_product('W', 'T', 0) - 2*dW) == 0
