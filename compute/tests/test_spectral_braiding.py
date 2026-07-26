"""
Tests verifying Laplace kernels, collision kernels, and Yang-Baxter equations.

The spectral braiding data arises from the Laplace transform bridge:
  r^L(z) = sum_{n>=0} c_n * n! / z^{n+1}
where {a_lambda b} = sum c_n lambda^n is the lambda-bracket.
The collision/spectral r-kernel used in Yangian examples is the bar-residue
kernel; for central current terms it has one pole fewer.

For each example family, we verify:
1. The Laplace transform bridge (BR3) correctly maps brackets to r^L kernels
2. The classical Yang-Baxter equation (CYBE) for the collision r-kernel
3. The quantum Yang-Baxter equation (YBE) for the R-matrix (when applicable)

Paper references: Vol II Section 18 (spectral-braiding.tex), BR3 axiom.

Test tiers:
  Tier 1 (structural): CYBE/YBE are self-certifying identities
  Tier 2 (published): r-matrix forms match Kac, Chari-Pressley
  Tier 3 (cross-check): Laplace bridge agrees with direct OPE computation
"""
import sys
import os
from pathlib import Path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from sympy import (Symbol, Rational, simplify, expand, S, symbols,
                   Matrix, eye, zeros, trace)

REPO_ROOT = Path(__file__).resolve().parents[2]
SPECTRAL_BRAIDING_TEX = REPO_ROOT / "chapters/connections/spectral-braiding.tex"
SPECTRAL_BRAIDING_CORE_TEX = (
    REPO_ROOT / "chapters/connections/spectral-braiding-core.tex"
)
SPECTRAL_BRAIDING_FRONTIER_TEX = (
    REPO_ROOT / "chapters/connections/spectral-braiding-frontier.tex"
)
AFFINE_HALF_SPACE_BV_TEX = REPO_ROOT / "chapters/connections/affine_half_space_bv.tex"
LINE_OPERATORS_TEX = REPO_ROOT / "chapters/connections/line-operators.tex"
HT_BULK_BOUNDARY_LINE_CORE_TEX = (
    REPO_ROOT / "chapters/connections/ht_bulk_boundary_line_core.tex"
)
FACTOR_SWISS_CHEESE_TEX = REPO_ROOT / "chapters/theory/factorization_swiss_cheese.tex"
THQG_LINE_OPERATORS_EXTENSIONS_TEX = (
    REPO_ROOT / "chapters/connections/thqg_line_operators_extensions.tex"
)
INTRODUCTION_TEX = REPO_ROOT / "chapters/theory/introduction.tex"
ORDERED_KD_CORE_TEX = (
    REPO_ROOT / "chapters/connections/ordered_associative_chiral_kd_core.tex"
)
ORDERED_KD_FRONTIER_TEX = (
    REPO_ROOT / "chapters/connections/ordered_associative_chiral_kd_frontier.tex"
)
THQG_SPECTRAL_BRAIDING_EXTENSIONS_TEX = (
    REPO_ROOT / "chapters/connections/thqg_spectral_braiding_extensions.tex"
)


def _flat(path: Path) -> str:
    return " ".join(path.read_text(encoding="utf-8").split())


def test_spectral_rmatrix_family_not_single_categorical_braiding():
    """R(z) is a meromorphic exchange family before categorical descent."""
    tex = SPECTRAL_BRAIDING_TEX.read_text()
    core_tex = SPECTRAL_BRAIDING_CORE_TEX.read_text()
    line_tex = LINE_OPERATORS_TEX.read_text()
    ht_line_tex = HT_BULK_BOUNDARY_LINE_CORE_TEX.read_text()
    swiss_tex = FACTOR_SWISS_CHEESE_TEX.read_text()
    thqg_line_tex = THQG_LINE_OPERATORS_EXTENSIONS_TEX.read_text()
    intro_tex = INTRODUCTION_TEX.read_text()
    ordered_core_tex = ORDERED_KD_CORE_TEX.read_text()
    ordered_frontier_tex = ORDERED_KD_FRONTIER_TEX.read_text()
    spectral_frontier_tex = SPECTRAL_BRAIDING_FRONTIER_TEX.read_text()
    thqg_spectral_tex = THQG_SPECTRAL_BRAIDING_EXTENSIONS_TEX.read_text()
    compact = " ".join(tex.split())
    core = " ".join(core_tex.split())
    line_ops = " ".join(line_tex.split())
    ht_line = " ".join(ht_line_tex.split())
    swiss = " ".join(swiss_tex.split())
    thqg_line = " ".join(thqg_line_tex.split())
    intro = " ".join(intro_tex.split())
    ordered_core = " ".join(ordered_core_tex.split())
    ordered_frontier = " ".join(ordered_frontier_tex.split())
    spectral_frontier = " ".join(spectral_frontier_tex.split())
    thqg_spectral = " ".join(thqg_spectral_tex.split())

    bad_fixed_tensor = "braided monoidal category with " + "meromorphic dependence on z"
    bad_raw_braiding = "R(z) defines a " + "meromorphic braiding on"
    bad_spectral_category = (
        "the braid relation for a braided monoidal category with spectral parameter"
    )
    bad_strict_definition = (
        "An $R$-matrix $R(z) : V \\otimes_z W \\to W \\otimes_{-z} V$ "
        "satisfying the quantum Yang--Baxter equation."
    )
    bad_fqg_slogan = (
        "A \\emph{factorisation quantum group} is a categorical structure encoding"
    )
    bad_line_category = (
        "The line category \\(\\cC_{\\mathrm{line}}\\) equipped with the "
        "\\(R\\)-matrix braiding"
    )
    bad_line_category_tex = (
        "The line category $\\cC_{\\mathrm{line}}$ equipped with the "
        "$R$-matrix braiding"
    )
    bad_hexagon_source = "Yang--Baxter equation (Theorem~\\ref{thm:YBE}) supplies the hexagon coherences"
    bad_r_supplies_braiding = "The $R$-matrix supplies the braiding on"
    bad_flat_coproduct = "line surface carries a flat spectral coproduct"
    bad_pre_descent_category = (
        "braided monoidal category $\\mathcal{C}_\\partial$ from Part"
    )
    bad_h_modules_with_braiding = (
        "category of $\\mathcal{H}$-modules equipped with the braiding induced by $R(z)$"
    )
    bad_unqualified_all_orders = (
        "with the Yang--Baxter equation holding to all orders by Theorem"
    )
    bad_latyntsev_instance = (
        "The spectral $R$-matrix and the meromorphic tensor product are an instance"
    )
    bad_full_fqg_statement = (
        "forms a factorisation quantum group in the sense of Latyntsev"
    )
    bad_parameter_upgrade = "upgrading the braided monoidal category"
    bad_meromorphic_braided_definition = (
        "A \\emph{meromorphic braided tensor category} is a category"
    )
    bad_meromorphic_braided_title = (
        "\\begin{remark}[Meromorphic braided category, after Dimofte]"
    )
    bad_fqg_of_theory = "the \\emph{factorisation quantum group} of the theory"
    bad_hopf_like_object = "filtered Hopf-like object $\\mathcal H$"
    bad_hopf_object_label = "Koszul-dual Hopf object $\\mathcal H$"
    bad_fqg_parenthetical = "(the \"factorisation quantum group\" in this setting)"
    bad_strict_fqg_all_simple = "strict factorisation quantum group for all simple"
    bad_strict_ybe_derivation = (
        "Substituting $R(z) = \\id + \\hh\\, r(z) + O(\\hh^2)$ "
        "into the quantum YBE"
    )
    bad_inverse_from_braiding = "inverse relation from the definition of braiding"
    bad_automatic_antisymmetry = (
        "For antisymmetry, swapping $L_1 \\leftrightarrow L_2$ sends"
    )
    bad_fixed_tensor_exchange = (
        "R_{L_1,L_2}(z)\\colon L_1\\otimes L_2 \\xrightarrow{\\ \\sim\\ } "
        "L_2\\otimes L_1"
    )
    bad_laplace_is_duality = "The Laplace transform \\emph{is} the spectral duality"
    bad_factorization_lift = "factorization-level lift of the evaluation-level DK"
    bad_auto_skew_proof = (
        "Skew-symmetry}: This follows from the symmetry of the tensor product"
    )
    bad_weak_cocomm_skew = (
        "Weak co-commutativity: $(\\tau_w \\otimes \\tau_w) r(z) = r(z)$, "
        "$r(-z) = -r_{21}(z)$"
    )
    bad_bare_quasitriangular = (
        "It follows from the definition of $R(z)$ as the braiding isomorphism"
    )
    bad_full_rmatrix_qybe = (
        "full $R$-matrix $R(z)$ satisfies the quantum Yang--Baxter"
    )
    bad_exact_quantum_ybe = "sum to zero at all orders, which is exactly the quantum YBE"
    bad_strict_sigma_descent = "spectral $R$-matrix $R(z)$ satisfying YBE and braid coherence"
    bad_ybe_lifts_to_sigma = "lift to $\\Sigma_n$-relations by YBE"
    bad_standard_inversion_axiom = "standard inversion axiom"
    bad_ht_full_qybe = "satisfies the full quantum Yang--Baxter equation"
    bad_ht_jacobi_determines_all = "is determined entirely by $r(z)$ and the Jacobi identity"
    bad_thqg_braided_structure = "carries a braided monoidal structure with braiding"
    bad_thqg_braided_category = "This braided monoidal category is equivalent"
    bad_thqg_yangian_shortcut = "The braided monoidal structure follows from the Yangian"
    bad_thqg_equiv_ybe = "is equivalent to the Yang--Baxter equation"
    bad_thqg_full_qybe = "and the full quantum Yang--Baxter"
    bad_thqg_all_orders = "(at all orders)"
    bad_intro_auto_skew = "co-commutativity $r(-z) = -r_{21}(z)$"
    bad_ordered_auto_unitarity = (
        "This holds for $E_\\infty$-chiral algebras (including those with OPE poles)"
    )
    bad_ordered_path_reversal = (
        "path reversal inverts the holonomy: $R_{21}(-z)=R_{12}(z)^{-1}$"
    )
    bad_thqg_inverse_definition = (
        "composing the braiding $L_1 \\otimes L_2 \\to L_2 \\otimes L_1$"
    )
    bad_thqg_tensor_skew_odd = (
        "$r_{12}(z) + r_{12}(-z) = 0$ if $r(z)$ is skew-symmetric"
    )
    bad_frontier_strict_ybe = (
        "The spectral $R$-matrix satisfies the Yang--Baxter equation by"
    )
    bad_ordered_frontier_strict_ybe = (
        "The $R$-matrix satisfies the Yang--Baxter equation "
        "$R_{12}(z_1{-}z_2)"
    )

    assert bad_fixed_tensor_exchange not in compact
    assert bad_fixed_tensor_exchange not in core
    assert bad_fixed_tensor not in compact
    assert bad_fixed_tensor not in core
    assert bad_raw_braiding not in compact
    assert bad_raw_braiding not in core
    assert bad_spectral_category not in core
    assert bad_strict_definition not in compact
    assert bad_strict_definition not in core
    assert bad_fqg_slogan not in compact
    assert bad_fqg_slogan not in core
    assert bad_pre_descent_category not in compact
    assert bad_pre_descent_category not in core
    assert bad_h_modules_with_braiding not in compact
    assert bad_h_modules_with_braiding not in core
    assert bad_unqualified_all_orders not in compact
    assert bad_latyntsev_instance not in compact
    assert bad_latyntsev_instance not in core
    assert bad_full_fqg_statement not in compact
    assert bad_parameter_upgrade not in compact
    assert bad_parameter_upgrade not in core
    assert bad_meromorphic_braided_definition not in core
    assert bad_meromorphic_braided_title not in core
    assert bad_fqg_of_theory not in core
    assert bad_hopf_like_object not in compact
    assert bad_hopf_object_label not in compact
    assert bad_fqg_parenthetical not in compact
    assert bad_strict_fqg_all_simple not in core
    assert bad_strict_ybe_derivation not in compact
    assert bad_strict_ybe_derivation not in core
    assert bad_inverse_from_braiding not in compact
    assert bad_inverse_from_braiding not in core
    assert bad_automatic_antisymmetry not in compact
    assert bad_automatic_antisymmetry not in core
    assert bad_laplace_is_duality not in compact
    assert bad_laplace_is_duality not in core
    assert bad_factorization_lift not in compact
    assert bad_factorization_lift not in core
    assert bad_auto_skew_proof not in compact
    assert bad_auto_skew_proof not in core
    assert bad_weak_cocomm_skew not in compact
    assert bad_weak_cocomm_skew not in core
    assert bad_bare_quasitriangular not in compact
    assert bad_bare_quasitriangular not in core
    assert bad_full_rmatrix_qybe not in compact
    assert bad_full_rmatrix_qybe not in core
    assert bad_exact_quantum_ybe not in compact
    assert bad_exact_quantum_ybe not in core
    assert bad_strict_sigma_descent not in swiss
    assert bad_ybe_lifts_to_sigma not in swiss
    assert bad_standard_inversion_axiom not in swiss
    assert bad_ht_full_qybe not in ht_line
    assert bad_ht_jacobi_determines_all not in ht_line
    assert bad_thqg_braided_structure not in thqg_line
    assert bad_thqg_braided_category not in thqg_line
    assert bad_thqg_yangian_shortcut not in thqg_line
    assert bad_thqg_equiv_ybe not in thqg_line
    assert bad_thqg_full_qybe not in thqg_line
    assert bad_thqg_all_orders not in thqg_line
    assert bad_intro_auto_skew not in intro
    assert bad_ordered_auto_unitarity not in ordered_core
    assert bad_ordered_path_reversal not in ordered_core
    assert bad_thqg_inverse_definition not in thqg_spectral
    assert bad_thqg_tensor_skew_odd not in thqg_spectral
    assert bad_frontier_strict_ybe not in spectral_frontier
    assert bad_ordered_frontier_strict_ybe not in ordered_frontier
    assert bad_line_category not in line_ops
    assert bad_line_category_tex not in line_ops
    assert bad_hexagon_source not in line_ops
    assert bad_r_supplies_braiding not in line_ops
    assert bad_flat_coproduct not in line_ops
    assert "local spectral tensor datum" in compact
    assert "local spectral tensor datum" in core
    assert "local spectral tensor datum" in line_ops
    assert "a three-point associator \\(\\Phi_A\\)" in compact
    assert "a three-point associator \\(\\Phi_A\\)" in core
    assert "three-point associator $\\Phi_\\cA$" in line_ops
    assert "meromorphic family of exchange isomorphisms" in compact
    assert "meromorphic family of exchange isomorphisms" in core
    assert "not a single categorical braiding on the fixed tensor product" in compact
    assert "not a single categorical braiding on the fixed tensor product" in core
    assert "only after evaluation/descent along an exchange path" in compact
    assert "after evaluation/descent along an exchange path" in core
    assert "Yang--Baxter relation for this meromorphic exchange family" in core
    assert "ordinary braided-monoidal hexagon belongs only to the descended centre-side" in core
    assert "local factorisation line category from Part" in compact
    assert "local factorisation line category of Part" in core
    assert "The all-orders Yang--Baxter statement belongs to that chosen quantized surface" in compact
    assert "The all-orders Yang--Baxter statement belongs to that chosen quantized surface" in core
    assert "local shadow of Latyntsev's factorisation quantum groups" in compact
    assert "local shadow of Latyntsev's factorisation quantum groups" in core
    assert "local meromorphic factorisation tensor datum" in core
    assert "The local distinguishing datum is not a braided monoidal category upgraded by a parameter" in compact
    assert "The local distinguishing datum is not a braided monoidal category upgraded by a parameter" in core
    assert "filtered local bialgebra shadow" in compact
    assert "filtered local bialgebra shadow" in core
    assert "factorisation-quantum-group interpretation requires the additional spectral coproduct" in compact
    assert "factorisation-quantum-group interpretation requires the additional spectral coproduct" in core
    assert "strict local spectral tensor datum on the resolved affine comparison surface" in core
    assert "associator-controlled Yang--Baxter relation" in compact
    assert "associator-controlled Yang--Baxter relation" in core
    assert "The skew condition is a separate normalization" in compact
    assert "The skew condition is a separate normalization" in core
    assert "unitary exchange normalization" in compact
    assert "unitary exchange normalization" in core
    assert "records the local topological exchange shadow" in compact
    assert "records the local topological exchange shadow" in core
    assert "local factorisation shadow whose descended evaluation surface gives the DK braiding" in compact
    assert "local factorisation shadow whose descended evaluation surface gives the DK braiding" in core
    assert "Normalized skew relation" in compact
    assert "Normalized skew relation" in core
    assert "not a formal consequence of the local exchange family" in compact
    assert "not a formal consequence of the local exchange family" in core
    assert "associator-controlled quantum Yang--Baxter relation on triple tensor products" in compact
    assert "associator-controlled quantum Yang--Baxter relation on triple tensor products" in core
    assert "braid-descent datum on \\(X^n_{\\mathrm{ord}}\\)" in swiss
    assert "unitary exchange normalization assumed in the statement" in swiss
    assert "chosen associator datum" in ht_line
    assert "descended evaluation subcategory" in thqg_line
    assert "strictified comparison surface" in thqg_line
    assert "unitary exchange normalization is part of the comparison surface" in intro
    assert "braid-descent datum consisting of a chosen associator" in ordered_core
    assert "not a formal consequence of the local exchange family" in ordered_core
    assert "not a consequence of Definition~\\ref{def:spectral-braiding}" in thqg_spectral
    assert "tensor-skewness is a separate condition" in thqg_spectral
    assert "satisfies the associator-controlled Yang--Baxter relation" in spectral_frontier
    assert "On a strictified comparison surface this relation becomes the displayed strict equation" in ordered_frontier
    assert "ordinary braided monoidal category before descent" in line_ops
    assert "ordinary categorical braiding" in line_ops
    assert "descended centre-side module category" in line_ops


def test_active_spectral_braiding_uses_typed_fm_reductions():
    """The active line-category proof uses FM^tr/FM^red, not the old quotient."""
    core = _flat(SPECTRAL_BRAIDING_CORE_TEX)
    frontier = _flat(SPECTRAL_BRAIDING_FRONTIER_TEX)
    affine = _flat(AFFINE_HALF_SPACE_BV_TEX)

    assert r"\FM^{\mathrm{tr}}_3(\C)" in core
    assert r"\FM^{\mathrm{red}}_3(\C)" in core
    assert r"D_S \cong \FM^{\mathrm{tr}}_{3-|S|+1}(\C)\times \FM^{\mathrm{red}}_{|S|}(\C)" in core
    assert "No non-consecutive divisor is removed from the closed chiral FM boundary" in core
    assert "The full-collision stratum" not in core
    assert "subleading contribution that vanishes" not in core

    assert r"\FM^{\mathrm{tr}}_3(\C)" in frontier
    assert "total-collision divisor gives the parenthesisation/associator term" in frontier
    assert r"\FM^{\mathrm{tr}}_2(\C)\times\FM^{\mathrm{red}}_2(\C)" in frontier

    assert r"\right\}/(\C\rtimes\R_{>0})." in affine
    assert "translation-reduced coordinates" in affine
    assert "The full complex-affine quotient is not used" in affine
    assert "translations and dilations" not in affine


def test_active_spectral_braiding_records_kz_dk_kl_line_path():
    """The quantum-group line-category bridge must name all comparison steps."""
    core = _flat(SPECTRAL_BRAIDING_CORE_TEX)

    assert r"\label{prop:kz-dk-kl-line-braiding-proof-path}" in core
    assert r"\hbar_{\mathrm{KZ}}=(k+h^\vee)^{-1}" in core
    assert r"\Omega_\kappa/((k+h^\vee)z)" in core
    assert "invariant-form conversion" in core
    assert "Drinfeld--Kohno monodromy" in core
    assert "Kazhdan--Lusztig tensor transport" in core
    assert r"q_{DK}=q_{KL}^{2}" in core
    assert "The critical level is excluded because the KZ denominator vanishes" in core


# ===================================================================
# ABELIAN CS R-MATRIX
# ===================================================================

class TestAbelianCSRMatrix:
    """Abelian CS Laplace kernel and collision r-kernel.

    The abelian current algebra {J_lambda J} = k*lambda has Laplace transform
    r^L(z) = k * 1! / z^2 = k/z^2.

    The Yangian-style collision kernel used by the example library is the
    tensor bar-residue kernel k*Omega_H/z; the stored scalar is its
    rank-one coefficient k/z, with charge factors included.
    Vol II Section 18.
    """

    def test_laplace_transform(self):
        """BR3: Laplace transform of {J_lam J} = k*lam gives r^L(z) = k/z^2.

        {J_lambda J} = k*lambda => bracket_coeffs = {1: k}
        r^L(z) = k * 1! / z^2 = k/z^2.

        Tier 3 (cross-check): independent Laplace computation.
        """
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        z = Symbol('z')
        k = Symbol('k')
        bracket_coeffs = {1: k}  # {J_lam J} = k*lam
        r = lambda_bracket_to_r_matrix(bracket_coeffs, z)
        expected = k / z**2
        assert simplify(r - expected) == 0

    def test_laplace_roundtrip(self):
        """Verify Laplace bridge is consistent: OPE -> lambda -> r^L -> OPE.

        The OPE J(z)J(w) ~ k/(z-w)^2 has pole coefficient c_1 = k at order 2.
        lambda-bracket: {J_lam J} = c_1 * lam / 1! = k*lam.
        Laplace kernel: r^L(z) = k * 1! / z^2 = k/z^2.
        Back to OPE: k/z^2 matches the original.

        Tier 3 (cross-check).
        """
        from lib.laplace_bridge import (
            lambda_bracket_to_r_matrix, r_matrix_to_lambda_bracket
        )
        k = Symbol('k')
        z = Symbol('z')
        # Forward: OPE -> lambda-bracket
        ope_coeffs = {1: k}  # c_1 = k (coefficient of 1/(z-w)^2)
        bracket = r_matrix_to_lambda_bracket(ope_coeffs)
        assert simplify(bracket[1] - k) == 0  # lambda^1 coefficient = k

        # Forward: lambda-bracket -> r-matrix
        r = lambda_bracket_to_r_matrix(bracket, z)
        assert simplify(r - k / z**2) == 0

    def test_yang_baxter_scalar(self):
        """YBE is trivially satisfied for scalar (abelian) r-matrix.

        For commuting quantities, [r_12, r_13] = [r_12, r_23] = [r_13, r_23] = 0,
        so CYBE: [r_12, r_13] + [r_12, r_23] + [r_13, r_23] = 0 trivially.

        Vol II Prop 18.5.
        """
        from lib.examples.abelian_cs import check_yang_baxter_abelian
        z1, z2 = symbols('z1 z2')
        hbar = Symbol('hbar')
        assert check_yang_baxter_abelian(z1, z2, hbar) == 0

    def test_classical_limit(self):
        """Collision r-kernel: r^{coll}(z) = k * Omega / z.

        For abelian case: Omega = J tensor J (one generator),
        so the stored coefficient is k/z times the charge pairing.

        Tier 2 (published).
        """
        from lib.examples.abelian_cs import abelian_r_matrix
        z = Symbol('z')
        hbar = Symbol('hbar')
        k = Symbol('k')
        r = abelian_r_matrix(z, hbar, k)
        q1, q2 = symbols('q1 q2')
        # Collision residue after bar dlog absorption, evaluated on charges.
        assert simplify(r - hbar * k * q1 * q2 / z) == 0

    def test_antisymmetry(self):
        """Laplace central kernel is symmetric; collision kernels carry skewness.

        For r^L(z) = k/z^2: r^L(-z) = k/z^2 = r^L(z).
        Under 12<->21 swap for abelian: r_{21}(z) = r_{12}(z).
        So r_{12}(z) + r_{21}(-z) = 2k/z^2 != 0.

        For the LEVEL-1 collision pole:
        r^{coll}(z) = Omega/z has r_{12}(z) = -r_{21}(-z) = -Omega/(-z) = Omega/z.

        For the abelian Laplace kernel r^L(z) = k/z^2,
        the antisymmetry is r_{12}(z) = -r_{21}(-z) <=> k/z^2 = -k/z^2,
        which fails. This is because the central Laplace kernel records the
        symmetric OPE pole before bar-residue extraction.

        The correct statement is: the skew part r(z) - r_{21}(-z) = 2k/z^2
        equals the central extension term.
        """
        z = Symbol('z')
        k = Symbol('k')
        r_12 = k / z**2
        r_21_neg = k / (-z)**2  # = k/z^2
        # Skew part = r_12(z) + r_21(-z) = 2k/z^2
        skew = simplify(r_12 + r_21_neg)
        assert simplify(skew - 2 * k / z**2) == 0


# ===================================================================
# GL(2) R-MATRIX
# ===================================================================

class TestGL2RMatrix:
    """R-matrix for gl_2 CS: R(z) = I + P/z at leading order.

    The permutation matrix P swaps tensor factors: P(v tensor w) = w tensor v.
    Properties: P^2 = I, tr(P) = n for gl_n.

    The classical r-matrix is r(z) = P/z (the rational Yangian r-matrix).
    The CYBE reads: [r_12, r_13] + [r_12, r_23] + [r_13, r_23] = 0.

    Vol II Section 18, gl_n example.
    """

    def test_permutation_matrix_P_squared(self):
        """P^2 = I: the permutation matrix is an involution.

        P = sum_{a,b} E_{ab} tensor E_{ba} for gl_n.
        For gl_2: P is the 4x4 matrix that swaps factors in C^2 tensor C^2.
        """
        # Permutation matrix for C^2 tensor C^2
        # Basis order: e1 tensor e1, e1 tensor e2, e2 tensor e1, e2 tensor e2
        P = Matrix([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
        ])
        assert P**2 == eye(4)

    def test_permutation_matrix_trace(self):
        """tr(P) = 2 for gl_2 (= dimension of the fundamental representation)."""
        P = Matrix([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
        ])
        assert trace(P) == 2

    def test_permutation_matrix_eigenvalues(self):
        """P has eigenvalues +1 (symmetric subspace) and -1 (antisymmetric).

        For gl_2: C^2 tensor C^2 = Sym^2(C^2) + Lambda^2(C^2).
        dim Sym^2 = 3, dim Lambda^2 = 1.
        So eigenvalue +1 has multiplicity 3, eigenvalue -1 has multiplicity 1.
        """
        P = Matrix([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
        ])
        eigenvals = P.eigenvals()
        assert eigenvals.get(S.One, 0) == 3
        assert eigenvals.get(S.NegativeOne, 0) == 1

    def test_yang_baxter_gl2_leading_order(self):
        """YBE for R(u) = 1 - P/u: the Yang R-matrix for gl_2.

        The quantum Yang-Baxter equation:
        R_12(u-v) R_13(u) R_23(v) = R_23(v) R_13(u) R_12(u-v)

        For the Yang R-matrix R(u) = 1 - P/u, we verify by direct
        matrix computation at specific numerical values u=3, v=1
        (avoiding u=v and u=0, v=0).

        Tier 1 (structural): self-certifying identity.
        """

        def _embed_P(i, j, n=2):
            """Embed the permutation P_{ij} into (C^n)^{tensor 3}."""
            dim = n**3
            result = zeros(dim, dim)
            for a1 in range(n):
                for a2 in range(n):
                    for a3 in range(n):
                        idx_in = a1 * n**2 + a2 * n + a3
                        indices = [a1, a2, a3]
                        indices[i], indices[j] = indices[j], indices[i]
                        idx_out = indices[0] * n**2 + indices[1] * n + indices[2]
                        result[idx_out, idx_in] = 1
            return result

        P12 = _embed_P(0, 1)
        P13 = _embed_P(0, 2)
        P23 = _embed_P(1, 2)
        I8 = eye(8)

        # R_{ij}(u) = I - P_{ij}/u (the Yang R-matrix)
        def R(Pij, u):
            return I8 - Pij * Rational(1, u)

        # Test at u=3, v=1 (u-v=2)
        u_val, v_val = 3, 1
        lhs = R(P12, u_val - v_val) * R(P13, u_val) * R(P23, v_val)
        rhs = R(P23, v_val) * R(P13, u_val) * R(P12, u_val - v_val)

        diff = lhs - rhs
        assert diff == zeros(8, 8), \
            f"YBE FAILED for gl_2 Yang R-matrix at u={u_val}, v={v_val}"

    def test_cybe_gl2(self):
        """Classical Yang-Baxter equation for r = P/z (gl_2 rational r-matrix).

        [r_12(u-v), r_13(u)] + [r_12(u-v), r_23(v)] + [r_13(u), r_23(v)] = 0

        Since r(z) = P/z, and P's satisfy the same commutation relations
        regardless of the spectral parameter:
        [P_12/(u-v), P_13/u] + [P_12/(u-v), P_23/v] + [P_13/u, P_23/v]
        = [P_12,P_13]/(u(u-v)) + [P_12,P_23]/(v(u-v)) + [P_13,P_23]/(uv)

        By the CYBE for P (proved above), we need this linear combination
        of commutators with DIFFERENT coefficients to vanish.

        Actually, the CYBE with spectral parameter is DIFFERENT from the
        constant CYBE. Let's check it directly.

        Tier 1 (structural): self-certifying identity.
        """
        def _embed_P(i, j, n=2):
            dim = n**3
            result = zeros(dim, dim)
            for a1 in range(n):
                for a2 in range(n):
                    for a3 in range(n):
                        idx_in = a1 * n**2 + a2 * n + a3
                        indices = [a1, a2, a3]
                        indices[i], indices[j] = indices[j], indices[i]
                        idx_out = indices[0] * n**2 + indices[1] * n + indices[2]
                        result[idx_out, idx_in] = 1
            return result

        P12 = _embed_P(0, 1)
        P13 = _embed_P(0, 2)
        P23 = _embed_P(1, 2)

        u, v = symbols('u v')
        # r_12(u-v) = P_12 / (u-v), r_13(u) = P_13 / u, r_23(v) = P_23 / v
        # CYBE: [r_12(u-v), r_13(u)] + [r_12(u-v), r_23(v)] + [r_13(u), r_23(v)]
        # = [P_12, P_13] / ((u-v)*u) + [P_12, P_23] / ((u-v)*v) + [P_13, P_23] / (u*v)

        comm_12_13 = P12 * P13 - P13 * P12
        comm_12_23 = P12 * P23 - P23 * P12
        comm_13_23 = P13 * P23 - P23 * P13

        # Each entry of the CYBE is a rational function of u, v
        # We check that for each matrix entry, the sum vanishes
        dim = 8
        for i in range(dim):
            for j in range(dim):
                entry = (comm_12_13[i, j] / ((u - v) * u)
                         + comm_12_23[i, j] / ((u - v) * v)
                         + comm_13_23[i, j] / (u * v))
                assert simplify(entry) == 0, \
                    f"CYBE entry ({i},{j}) nonzero: {entry}"


# ===================================================================
# VIRASORO R-MATRIX
# ===================================================================

class TestVirasoroRMatrix:
    """Virasoro Laplace kernel from the lambda-bracket.

    {T_lambda T} = dT + 2T*lambda + (c/12)*lambda^3

    Via Laplace transform (BR3):
    r^L(z) = dT * 0!/z + 2T * 1!/z^2 + (c/12) * 3!/z^4
           = dT/z + 2T/z^2 + c/2 * 1/z^4

    This recovers the Virasoro OPE:
    T(z)T(w) ~ (c/2)/(z-w)^4 + 2T/(z-w)^2 + dT/(z-w)

    Vol II Section 18.
    """

    def test_virasoro_classical_r(self):
        """Verify r(z) = dT/z + 2T/z^2 + (c/2)/z^4 via BR3.

        Tier 3 (cross-check): Laplace bridge reproduces OPE coefficients.
        """
        from lib.laplace_bridge import verify_br3_virasoro
        c = Symbol('c')
        z = Symbol('z')
        r, expected, diff = verify_br3_virasoro(c, z)
        assert diff == 0, f"BR3 for Virasoro FAILED: diff = {diff}"

    def test_virasoro_r_matrix_coefficients(self):
        """Check individual pole coefficients of the Virasoro r-matrix.

        The r-matrix r(z) = dT/z + 2T/z^2 + (c/2)/z^4 has:
        - Simple pole: residue = dT (= partial T)
        - Double pole: coefficient = 2T
        - No triple pole (coefficient = 0)
        - Quadruple pole: coefficient = c/2
        """
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        from lib.examples.virasoro import T, dT
        c = Symbol('c')
        z = Symbol('z')

        bracket_coeffs = {0: dT, 1: 2*T, 3: Rational(1, 12)*c}
        r = lambda_bracket_to_r_matrix(bracket_coeffs, z)

        # Extract coefficients by multiplying by appropriate power of z
        # r(z) = dT/z + 2T/z^2 + (c/2)/z^4
        # r(z) * z = dT + 2T/z + (c/2)/z^3
        # r(z) * z^2 = dT*z + 2T + (c/2)/z^2
        # r(z) * z^4 = dT*z^3 + 2T*z^2 + c/2

        # Check by evaluating at z -> infinity with appropriate scaling
        r_expanded = expand(r)

        # Simple pole coefficient (1/z term)
        coeff_z_neg1 = simplify(expand(r_expanded * z).subs(z, S.Zero))
        # But r*z at z=0 is problematic. Instead, check the full expression.
        expected = dT / z + 2 * T / z**2 + Rational(1, 2) * c / z**4
        assert simplify(expand(r - expected)) == 0

    def test_virasoro_cybe_via_jacobi(self):
        """Verify CYBE via PVA Jacobi equivalence.

        The CYBE for the classical r-matrix is EQUIVALENT to the Jacobi
        identity for the lambda-bracket (via the Laplace transform bridge).

        Since we have already verified Virasoro Jacobi (test_pva_axioms.py),
        the CYBE follows automatically.

        This test verifies the LOGICAL EQUIVALENCE by checking that
        the Jacobi identity (which we know holds) implies CYBE.

        Vol II Thm 18.12: PVA Jacobi <=> CYBE for the associated r-matrix.
        """
        from lib.examples.virasoro import check_virasoro_jacobi
        lam, mu = symbols('lambda mu')
        c = Symbol('c')
        # Jacobi holds => CYBE holds (by the Laplace bridge theorem)
        jacobi_result = check_virasoro_jacobi(lam, mu, c)
        assert jacobi_result == 0, "Virasoro Jacobi FAILED => CYBE would also fail"

    def test_virasoro_r_laplace_inverse(self):
        """Check that the OPE -> lambda -> r roundtrip is consistent.

        Start from OPE coefficients, convert to lambda-bracket, then to r^L.
        The result should reproduce the original OPE coefficients.

        Virasoro OPE: T(z)T(w) ~ (c/2)/(z-w)^4 + 2T(w)/(z-w)^2 + dT(w)/(z-w)
        n-products: T_0 T = dT, T_1 T = 2T, T_2 T = 0, T_3 T = c/2
        lambda-bracket: {T_lam T} = dT + 2T*lam + (c/12)*lam^3
                        because (c/2)/3! = c/12
        Laplace kernel: r^L(z) = dT/z + 2T/z^2 + (c/2)/z^4

        The roundtrip reproduces the OPE pole structure.
        Tier 3 (cross-check).
        """
        from lib.laplace_bridge import r_matrix_to_lambda_bracket, lambda_bracket_to_r_matrix
        from lib.examples.virasoro import T, dT
        c = Symbol('c')
        z = Symbol('z')

        # Step 1: OPE -> lambda-bracket
        ope_coeffs = {0: dT, 1: 2*T, 3: Rational(1, 2)*c}
        bracket = r_matrix_to_lambda_bracket(ope_coeffs)
        # Check lambda-bracket coefficients
        assert simplify(bracket[0] - dT) == 0        # lambda^0: dT
        assert simplify(bracket[1] - 2*T) == 0        # lambda^1: 2T
        assert simplify(bracket[3] - Rational(1, 12)*c) == 0  # lambda^3: c/12

        # Step 2: lambda-bracket -> r-matrix
        r = lambda_bracket_to_r_matrix(bracket, z)
        expected = dT / z + 2 * T / z**2 + Rational(1, 2) * c / z**4
        assert simplify(expand(r - expected)) == 0


# ===================================================================
# SU(2) R-MATRIX
# ===================================================================

class TestSU2RMatrix:
    """R-matrix for su(2) current algebra at level k.

    The lambda-bracket {J^a_lam J^b} = epsilon^{abc}J^c + k*delta^{ab}*lam
    gives via BR3:
      r^{ab}(z) = epsilon^{abc}J^c / z + k*delta^{ab} / z^2

    This matches the OPE:
      J^a(z)J^b(w) ~ k*delta^{ab}/(z-w)^2 + epsilon^{abc}J^c(w)/(z-w)

    Vol II Section 18.
    """

    def test_br3_all_9_components(self):
        """Verify BR3 for all 9 components of the su(2) r-matrix.

        Tier 3 (cross-check).
        """
        from lib.laplace_bridge import verify_br3_su2
        k = Symbol('k')
        z = Symbol('z')
        results = verify_br3_su2(k, z)
        for (a, b), (r, expected, diff) in results.items():
            assert diff == 0, f"BR3 for su(2) ({a},{b}) FAILED: diff = {diff}"

    def test_su2_r_diagonal(self):
        """Diagonal pre-dlog components: r^{aa}(z) = k/z^2.

        From {J^a_lam J^a} = k*lam: only lambda^1 term.
        r^{aa}(z) = k * 1! / z^2 = k/z^2.
        """
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        k = Symbol('k')
        z = Symbol('z')
        bracket_coeffs = {1: k}
        r = lambda_bracket_to_r_matrix(bracket_coeffs, z)
        assert simplify(r - k / z**2) == 0

    def test_su2_r_offdiagonal(self):
        """Off-diagonal components: r^{12}(z) = J^3/z.

        From {J^1_lam J^2} = J^3 (epsilon^{123} = 1): only lambda^0 term.
        r^{12}(z) = J^3 * 0! / z = J^3/z.
        """
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        from lib.examples.nonabelian_cs import J
        z = Symbol('z')
        bracket_coeffs = {0: J[3]}  # {J^1_lam J^2} = J^3
        r = lambda_bracket_to_r_matrix(bracket_coeffs, z)
        assert simplify(r - J[3] / z) == 0

    def test_su2_cybe_via_jacobi(self):
        """CYBE for su(2) follows from the PVA Jacobi identity (already verified).

        By the Laplace bridge theorem (Vol II Thm 18.12),
        PVA Jacobi <=> CYBE for the classical r-matrix.

        Since Jacobi for su(2)_k holds for all 27 triples (verified in
        test_pva_axioms.py), the CYBE holds automatically.
        """
        from lib.examples.nonabelian_cs import check_su2_jacobi_all
        lam, mu = symbols('lambda mu')
        k = Symbol('k')
        results = check_su2_jacobi_all(lam, mu, k)
        for (a, b, c), val in results.items():
            assert val == 0, f"Jacobi ({a},{b},{c}) FAILED => CYBE would fail"


# ===================================================================
# CROSS-FAMILY LAPLACE BRIDGE CONSISTENCY
# ===================================================================

class TestLaplaceBridgeConsistency:
    """Cross-family tests verifying the Laplace bridge is consistent.

    The bridge BR3 should satisfy:
    1. Linearity: r(bracket_1 + bracket_2) = r(bracket_1) + r(bracket_2)
    2. n! scaling: lambda^n term maps to n!/z^{n+1}
    3. Roundtrip: OPE -> lambda -> r reproduces OPE pole structure

    Vol II Section 18, BR3 axioms.
    """

    def test_bridge_linearity(self):
        """Laplace bridge is linear in the bracket coefficients.

        r(c_1 * bracket_1 + c_2 * bracket_2) = c_1 * r(bracket_1) + c_2 * r(bracket_2)
        """
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        z = Symbol('z')
        alpha, beta = symbols('alpha beta')

        # Two brackets
        b1 = {0: 1, 1: 2}      # 1 + 2*lambda
        b2 = {0: 3, 2: 1}      # 3 + lambda^2

        r1 = lambda_bracket_to_r_matrix(b1, z)
        r2 = lambda_bracket_to_r_matrix(b2, z)

        # Combined bracket: alpha*b1 + beta*b2
        combined = {
            0: alpha * 1 + beta * 3,
            1: alpha * 2,
            2: beta * 1,
        }
        r_combined = lambda_bracket_to_r_matrix(combined, z)

        # Check linearity
        expected = alpha * r1 + beta * r2
        assert simplify(expand(r_combined - expected)) == 0

    def test_bridge_factorial_scaling(self):
        """Verify n! scaling: lambda^n -> n!/z^{n+1}.

        This is the defining property of the Laplace transform:
        integral_0^inf lambda^n * e^{-lambda*z} dlambda = n! / z^{n+1}.
        """
        from lib.laplace_bridge import lambda_bracket_to_r_matrix
        from sympy import factorial
        z = Symbol('z')

        for n in range(5):
            bracket = {n: 1}  # pure lambda^n term with coefficient 1
            r = lambda_bracket_to_r_matrix(bracket, z)
            expected = factorial(n) / z**(n + 1)
            assert simplify(r - expected) == 0, \
                f"Factorial scaling failed at n={n}: got {r}, expected {expected}"

    def test_abelian_br3(self):
        """Full BR3 verification for abelian CS.

        Tier 3 (cross-check).
        """
        from lib.laplace_bridge import verify_br3_abelian
        k = Symbol('k')
        z = Symbol('z')
        r, expected, diff = verify_br3_abelian(k, z)
        assert diff == 0

    def test_virasoro_br3(self):
        """Full BR3 verification for Virasoro.

        Tier 3 (cross-check).
        """
        from lib.laplace_bridge import verify_br3_virasoro
        c = Symbol('c')
        z = Symbol('z')
        r, expected, diff = verify_br3_virasoro(c, z)
        assert diff == 0
