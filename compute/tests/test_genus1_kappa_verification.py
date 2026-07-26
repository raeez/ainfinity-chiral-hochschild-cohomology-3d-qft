r"""
test_genus1_kappa_verification.py -- Tests that kappa(A) is DERIVED from
bar-complex self-sewing, not hardcoded.

Each test constructs OPE data from first principles, computes kappa via
the genus-1 self-loop graph amplitude, and verifies against the independently
known formula.

The four test families:
  1. Heisenberg H_k:     kappa = k                          (from J_{(1)} J = k)
  2. Virasoro Vir_c:     kappa = c/2                        (from T_{(3)} T = c/2)
  3. Affine sl_2 at k:   kappa = 3(k+2)/4                   (from Casimir trace + Sugawara)
  4. AP19 pole absorption: bar residue poles one order below OPE poles

Plus:
  5. W_3:                kappa = 5c/6                        (from T + W contributions)
  6. W_N for N=2..6:     kappa = c * (H_N - 1)              (multi-generator sum)
  7. Affine g_k for all simple types: kappa = dim(g)(k+h^v)/(2h^v)
  8. Cross-family additivity: kappa(A tensor B) = kappa(A) + kappa(B)
  9. Koszul dual complementarity: kappa + kappa' checks for specific families
  10. Beta-gamma cross-sewing: kappa = 1
  11. Free fermion: kappa = 1/2
"""

import pytest
from fractions import Fraction

from compute.lib.genus1_kappa_verification import (
    # OPE data constructors
    heisenberg_ope,
    virasoro_ope,
    affine_sl2_ope,
    affine_ope,
    w3_ope,
    w_n_ope,
    betagamma_ope,
    free_fermion_ope,
    # Computation functions
    kappa_from_self_sewing,
    kappa_betagamma_from_sewing,
    self_sewing_residue,
    extract_leading_scalar_n_product,
    verify_pole_absorption,
    bar_residue_pole_orders,
    verify_kappa_from_first_principles,
    # Expected values (independent computation)
    expected_kappa_heisenberg,
    expected_kappa_virasoro,
    expected_kappa_affine_sl2,
    expected_kappa_affine,
    expected_kappa_w_n,
)


# ================================================================
# 1. Heisenberg: kappa(H_k) = k
# ================================================================

class TestHeisenbergKappa:
    """kappa(Heisenberg) from J_{(1)}J = k (the leading scalar n-product)."""

    def test_heisenberg_k1(self):
        """kappa(H_1) = 1."""
        ope = heisenberg_ope(Fraction(1))
        result = verify_kappa_from_first_principles(
            ope, expected_kappa_heisenberg(Fraction(1)),
        )
        assert result['match'], result['note']

    def test_heisenberg_k_generic(self):
        """kappa(H_k) = k for several values of k."""
        for k_val in [1, 2, 3, 5, 10, -1, Fraction(1, 2), Fraction(7, 3)]:
            k = Fraction(k_val)
            ope = heisenberg_ope(k)
            kappa = kappa_from_self_sewing(ope)
            assert kappa == k, f"kappa(H_{k}) = {kappa}, expected {k}"

    def test_heisenberg_self_sewing_residue(self):
        """The self-sewing residue R_{JJ} = J_{(1)}J = k."""
        ope = heisenberg_ope(Fraction(7))
        R = self_sewing_residue(ope, 'J')
        assert R == Fraction(7), f"R_JJ = {R}, expected 7"

    def test_heisenberg_leading_n_product(self):
        """The leading scalar n-product of J is J_{(1)}J."""
        ope = heisenberg_ope(Fraction(3))
        val = extract_leading_scalar_n_product(ope, 'J', 'J')
        assert val == Fraction(3), f"J_{{(1)}}J = {val}, expected 3"


# ================================================================
# 2. Virasoro: kappa(Vir_c) = c/2
# ================================================================

class TestVirasoroKappa:
    """kappa(Virasoro) from T_{(3)}T = c/2 (the leading scalar n-product)."""

    def test_virasoro_c1(self):
        """kappa(Vir_1) = 1/2."""
        ope = virasoro_ope(Fraction(1))
        result = verify_kappa_from_first_principles(
            ope, expected_kappa_virasoro(Fraction(1)),
        )
        assert result['match'], result['note']

    def test_virasoro_c2(self):
        """kappa(Vir_2) = 1."""
        ope = virasoro_ope(Fraction(2))
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(1), f"kappa(Vir_2) = {kappa}, expected 1"

    def test_virasoro_c25(self):
        """kappa(Vir_25) = 25/2."""
        ope = virasoro_ope(Fraction(25))
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(25, 2), f"kappa(Vir_25) = {kappa}, expected 25/2"

    def test_virasoro_c26(self):
        """kappa(Vir_26) = 13.  The critical string theory central charge."""
        ope = virasoro_ope(Fraction(26))
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(13), f"kappa(Vir_26) = {kappa}, expected 13"

    def test_virasoro_comparison_fixed_c13(self):
        """kappa(Vir_13) = 13/2 at the comparison fixed point."""
        ope = virasoro_ope(Fraction(13))
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(13, 2), f"kappa(Vir_13) = {kappa}, expected 13/2"

    def test_virasoro_c0(self):
        """kappa(Vir_0) = 0.  The bar complex is uncurved at c = 0."""
        ope = virasoro_ope(Fraction(0))
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(0), f"kappa(Vir_0) = {kappa}, expected 0"

    def test_virasoro_negative_c(self):
        """kappa(Vir_{-2}) = -1.  Negative central charge (ghost system)."""
        ope = virasoro_ope(Fraction(-2))
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(-1), f"kappa(Vir_-2) = {kappa}, expected -1"

    def test_virasoro_generic(self):
        """kappa(Vir_c) = c/2 for many values of c."""
        for c_val in [1, 2, 13, 25, 26, 0, -2, -26, Fraction(1, 2),
                      Fraction(7, 5), Fraction(22, 5)]:
            c = Fraction(c_val)
            ope = virasoro_ope(c)
            kappa = kappa_from_self_sewing(ope)
            expected = c / 2
            assert kappa == expected, f"kappa(Vir_{c}) = {kappa}, expected {expected}"

    def test_virasoro_self_sewing_residue(self):
        """R_{TT} = T_{(3)}T = c/2."""
        ope = virasoro_ope(Fraction(26))
        R = self_sewing_residue(ope, 'T')
        assert R == Fraction(13), f"R_TT = {R}, expected 13"

    def test_virasoro_leading_n_product(self):
        """Leading scalar n-product of T is T_{(3)}T = c/2.

        The n=3 n-product is the coefficient of (z-w)^{-4} in T(z)T(w),
        which IS the two-point function coefficient c/2.
        """
        for c_val in [1, 26, Fraction(1, 2)]:
            c = Fraction(c_val)
            ope = virasoro_ope(c)
            val = extract_leading_scalar_n_product(ope, 'T', 'T')
            assert val == c / 2, f"T_{{(3)}}T at c={c}: {val}, expected {c/2}"

    def test_virasoro_field_valued_n_products_not_scalar(self):
        """T_{(1)}T = 2T and T_{(0)}T = dT are field-valued, not scalar.

        The extraction function should skip these and return the highest
        SCALAR n-product, which is T_{(3)}T = c/2.
        """
        ope = virasoro_ope(Fraction(26))
        # Check that lower n-products are not returned
        val = extract_leading_scalar_n_product(ope, 'T', 'T')
        # Should be c/2 = 13, NOT the field-valued T_{(1)}T or T_{(0)}T
        assert val == Fraction(13), (
            f"Should return T_{{(3)}}T = 13, got {val}"
        )


# ================================================================
# 3. Affine sl_2: kappa = 3(k+2)/4
# ================================================================

class TestAffineSl2Kappa:
    """kappa(sl_2_k) = 3(k+2)/4 from Casimir trace + Sugawara shift."""

    def test_sl2_k1(self):
        """kappa(sl_2, k=1) = 3*3/4 = 9/4."""
        ope = affine_sl2_ope(Fraction(1))
        kappa = kappa_from_self_sewing(ope)
        expected = expected_kappa_affine_sl2(Fraction(1))
        assert kappa == expected == Fraction(9, 4), (
            f"kappa(sl_2, k=1) = {kappa}, expected {expected}"
        )

    def test_sl2_k2(self):
        """kappa(sl_2, k=2) = 3*4/4 = 3."""
        ope = affine_sl2_ope(Fraction(2))
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(3), f"kappa = {kappa}, expected 3"

    def test_sl2_k_generic(self):
        """kappa(sl_2, k) = 3(k+2)/4 for several values."""
        for k_val in [1, 2, 3, 4, 10, Fraction(1, 2), Fraction(5, 3)]:
            k = Fraction(k_val)
            ope = affine_sl2_ope(k)
            kappa = kappa_from_self_sewing(ope)
            expected = Fraction(3) * (k + 2) / 4
            assert kappa == expected, (
                f"kappa(sl_2, k={k}) = {kappa}, expected {expected}"
            )

    def test_sl2_critical_level(self):
        """At k = -h^v = -2 (critical level), kappa = 0."""
        ope = affine_sl2_ope(Fraction(-2))
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(0), f"kappa at critical level = {kappa}"

    def test_sl2_decomposition(self):
        """Verify the decomposition: kappa = raw_trace + sugawara_shift.

        raw_trace = dim(g) * k / (2h^v) = 3k/4
        sugawara_shift = dim(g) / 2 = 3/2
        total = 3k/4 + 3/2 = 3(k+2)/4
        """
        k = Fraction(5)
        ope = affine_sl2_ope(k)

        # Raw trace (without Sugawara shift)
        raw = Fraction(0)
        for gen in ope.generators:
            R = self_sewing_residue(ope, gen.name)
            eta_inv = Fraction(1) / ope.cyclic_pairing[(gen.name, gen.name)]
            raw += eta_inv * R
        # raw = 3 * (1/4) * 5 = 15/4
        assert raw == Fraction(15, 4), f"raw trace = {raw}"

        # Sugawara shift = 3/2
        assert ope.sugawara_shift == Fraction(3, 2)

        # Total
        kappa = kappa_from_self_sewing(ope)
        assert kappa == raw + ope.sugawara_shift == Fraction(21, 4)


# ================================================================
# 4. AP19 pole absorption: r-matrix poles one order below OPE
# ================================================================

class TestPoleAbsorption:
    """Verify AP19: bar residue has poles one order lower than OPE."""

    def test_virasoro_pole_absorption(self):
        """Virasoro OPE: max pole 4.  Bar residue: max pole 3."""
        ope = virasoro_ope(Fraction(26))
        result = verify_pole_absorption(ope, 'T', 'T')
        assert result['absorption_verified'], result['note']
        assert result['ope_max_pole'] == 4
        assert result['rmatrix_max_pole'] == 3

    def test_heisenberg_pole_absorption(self):
        """Heisenberg OPE: max pole 2.  Bar residue: max pole 1."""
        ope = heisenberg_ope(Fraction(1))
        result = verify_pole_absorption(ope, 'J', 'J')
        assert result['absorption_verified'], result['note']
        assert result['ope_max_pole'] == 2
        assert result['rmatrix_max_pole'] == 1

    def test_w3_T_pole_absorption(self):
        """W_3 T-channel: same as Virasoro."""
        ope = w3_ope(Fraction(10))
        result = verify_pole_absorption(ope, 'T', 'T')
        assert result['absorption_verified']
        assert result['ope_max_pole'] == 4

    def test_w3_W_pole_absorption(self):
        """W_3 W-channel OPE: max pole 6.  Bar residue: max pole 5."""
        ope = w3_ope(Fraction(10))
        result = verify_pole_absorption(ope, 'W', 'W')
        assert result['absorption_verified']
        assert result['ope_max_pole'] == 6
        assert result['rmatrix_max_pole'] == 5

    def test_sl2_pole_absorption(self):
        """Affine sl_2 OPE: max pole 2.  Bar residue: max pole 1."""
        ope = affine_sl2_ope(Fraction(1))
        result = verify_pole_absorption(ope, 'J1', 'J1')
        assert result['absorption_verified']
        assert result['ope_max_pole'] == 2
        assert result['rmatrix_max_pole'] == 1

    def test_pole_absorption_rule(self):
        """General rule: for weight h, OPE max pole is 2h, r-matrix is 2h-1."""
        for c_val in [1, 26]:
            ope = virasoro_ope(Fraction(c_val))
            h = 2  # weight of T
            result = verify_pole_absorption(ope, 'T', 'T')
            assert result['ope_max_pole'] == 2 * h
            assert result['rmatrix_max_pole'] == 2 * h - 1

    def test_bar_residue_pole_orders_virasoro(self):
        """The bar residue r(z) for Virasoro has poles at z^{-3} and z^{-1}.

        OPE: c/2 * z^{-4} + 2T * z^{-2} + dT * z^{-1}
        After d-log: r(z) has poles at z^{-3} (from c/2), z^{-1} (from 2T)

        AP19: the r-matrix has NO z^{-4} pole; the d-log absorbs it.
        For bosonic algebras, even-order poles in the OPE become
        odd-order poles in the r-matrix.
        """
        ope = virasoro_ope(Fraction(26))
        residues = bar_residue_pole_orders(ope, 'T', 'T')
        # n=3 gives r-matrix pole at z^{-3} with coefficient c/2 = 13
        assert 3 in residues
        assert residues[3] == Fraction(13)
        # n=1 gives r-matrix pole at z^{-1} with field-valued coefficient 'T'
        assert 1 in residues
        assert residues[1] == 'T'


# ================================================================
# 5. W_3: kappa = 5c/6
# ================================================================

class TestW3Kappa:
    """kappa(W_3) = c/2 + c/3 = 5c/6 from T and W contributions."""

    def test_w3_generic(self):
        """kappa(W_3, c) = 5c/6 for several values of c."""
        for c_val in [1, 6, 10, 26, Fraction(22, 5)]:
            c = Fraction(c_val)
            ope = w3_ope(c)
            kappa = kappa_from_self_sewing(ope)
            expected = Fraction(5) * c / 6
            assert kappa == expected, (
                f"kappa(W_3, c={c}) = {kappa}, expected {expected}"
            )

    def test_w3_decomposition(self):
        """T contributes c/2, W contributes c/3.  Total = 5c/6."""
        c = Fraction(12)
        ope = w3_ope(c)

        R_T = self_sewing_residue(ope, 'T')
        R_W = self_sewing_residue(ope, 'W')

        assert R_T == Fraction(6), f"R_TT = {R_T}, expected 6 (= c/2)"
        assert R_W == Fraction(4), f"R_WW = {R_W}, expected 4 (= c/3)"

        kappa = R_T + R_W  # eta^{ii} = 1 for both
        assert kappa == Fraction(10), f"kappa = {kappa}, expected 10 (= 5*12/6)"


# ================================================================
# 6. W_N: kappa = c * (H_N - 1)
# ================================================================

class TestWNKappa:
    """kappa(W_N) = c * sum_{s=2}^N 1/s = c * (H_N - 1)."""

    def test_w2_is_virasoro(self):
        """W_2 = Virasoro.  kappa = c/2."""
        c = Fraction(26)
        ope = w_n_ope(c, N=2)
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(13)

    def test_w3_from_w_n(self):
        """W_3 from the general W_N constructor.  kappa = 5c/6."""
        c = Fraction(6)
        ope = w_n_ope(c, N=3)
        kappa = kappa_from_self_sewing(ope)
        assert kappa == Fraction(5)  # 5 * 6 / 6

    def test_w4(self):
        """kappa(W_4) = c * (1/2 + 1/3 + 1/4) = 13c/12."""
        c = Fraction(12)
        ope = w_n_ope(c, N=4)
        kappa = kappa_from_self_sewing(ope)
        expected = expected_kappa_w_n(c, 4)
        assert expected == Fraction(13)
        assert kappa == expected

    def test_w5(self):
        """kappa(W_5) = c * (1/2 + 1/3 + 1/4 + 1/5) = 77c/60."""
        c = Fraction(60)
        ope = w_n_ope(c, N=5)
        kappa = kappa_from_self_sewing(ope)
        expected = expected_kappa_w_n(c, 5)
        assert expected == Fraction(77)
        assert kappa == expected

    def test_w6(self):
        """kappa(W_6) = c * (1/2 + 1/3 + 1/4 + 1/5 + 1/6) = 29c/20."""
        c = Fraction(20)
        ope = w_n_ope(c, N=6)
        kappa = kappa_from_self_sewing(ope)
        expected = expected_kappa_w_n(c, 6)
        # H_6 - 1 = 1/2 + 1/3 + 1/4 + 1/5 + 1/6 = 29/20
        assert expected == Fraction(29)
        assert kappa == expected

    @pytest.mark.parametrize("N", [2, 3, 4, 5, 6, 7, 8])
    def test_w_n_formula(self, N):
        """kappa(W_N) = c * (H_N - 1) for N = 2, ..., 8."""
        c = Fraction(1)
        ope = w_n_ope(c, N=N)
        kappa = kappa_from_self_sewing(ope)
        expected = expected_kappa_w_n(c, N)
        assert kappa == expected, (
            f"kappa(W_{N}, c=1) = {kappa}, expected {expected}"
        )


# ================================================================
# 7. Affine g_k for all simple types
# ================================================================

class TestAffineAllTypes:
    """kappa(g_k) = dim(g)(k+h^v)/(2h^v) for all simple Lie algebras."""

    @pytest.mark.parametrize("g_name,dim_g,h_vee", [
        ('sl2', 3, 2),
        ('sl3', 8, 3),
        ('sl4', 15, 4),
        ('sl5', 24, 5),
        ('so5', 10, 3),
        ('g2', 14, 4),
        ('so8', 28, 6),
        ('f4', 52, 9),
        ('e6', 78, 12),
        ('e7', 133, 18),
        ('e8', 248, 30),
    ])
    def test_affine_kappa(self, g_name, dim_g, h_vee):
        """kappa(g_k) = dim(g)(k+h^v)/(2h^v) at k=1."""
        k = Fraction(1)
        ope = affine_ope(g_name, k)
        kappa = kappa_from_self_sewing(ope)
        expected = expected_kappa_affine(g_name, k)
        assert kappa == expected, (
            f"kappa({g_name}, k=1) = {kappa}, expected {expected}"
        )

    def test_e8_k1(self):
        """kappa(E_8, k=1) = 248 * 31 / 60."""
        k = Fraction(1)
        kappa = kappa_from_self_sewing(affine_ope('e8', k))
        expected = Fraction(248) * (1 + 30) / 60
        assert kappa == expected == Fraction(248 * 31, 60)

    def test_critical_level_all_types(self):
        """At k = -h^v (critical level), kappa = 0 for all types."""
        for g_name, (dim_g, h_vee) in [
            ('sl2', (3, 2)), ('sl3', (8, 3)), ('g2', (14, 4)),
            ('e8', (248, 30)),
        ]:
            k = Fraction(-h_vee)
            ope = affine_ope(g_name, k)
            kappa = kappa_from_self_sewing(ope)
            assert kappa == Fraction(0), (
                f"kappa({g_name}, k=-h^v={-h_vee}) = {kappa}, should be 0"
            )


# ================================================================
# 8. Cross-family additivity
# ================================================================

class TestAdditivity:
    """kappa is additive: kappa(A tensor B) = kappa(A) + kappa(B)."""

    def test_virasoro_plus_heisenberg(self):
        """kappa(Vir_c tensor H_k) = c/2 + k."""
        c, k = Fraction(26), Fraction(1)
        kappa_vir = kappa_from_self_sewing(virasoro_ope(c))
        kappa_heis = kappa_from_self_sewing(heisenberg_ope(k))
        assert kappa_vir + kappa_heis == c / 2 + k == Fraction(14)

    def test_w3_decomposition_additivity(self):
        """W_3 kappa decomposes as T-contribution + W-contribution.

        This tests additivity at the generator level: each generator
        contributes independently to kappa.
        """
        c = Fraction(6)
        ope = w3_ope(c)
        kappa_T = self_sewing_residue(ope, 'T')
        kappa_W = self_sewing_residue(ope, 'W')
        kappa_total = kappa_from_self_sewing(ope)
        assert kappa_total == kappa_T + kappa_W


# ================================================================
# 9. Koszul dual complementarity (AP24)
# ================================================================

class TestComplementarity:
    """kappa + kappa' complementarity checks.

    WARNING (AP24): kappa + kappa' = 0 for KM/free fields/lattice,
    but kappa(Vir_c) + kappa(Vir_{26-c}) = 13, NOT 0.
    """

    def test_virasoro_complementarity(self):
        """kappa(Vir_c) + kappa(Vir_{26-c}) = 13 for all c.

        The same-family line-side representative is Vir_{26-c}.
        kappa(Vir_c) = c/2, kappa(Vir_{26-c}) = (26-c)/2.
        Sum = c/2 + (26-c)/2 = 13.
        """
        for c_val in [0, 1, 13, 25, 26, Fraction(1, 2)]:
            c = Fraction(c_val)
            kappa_A = kappa_from_self_sewing(virasoro_ope(c))
            kappa_dual = kappa_from_self_sewing(virasoro_ope(26 - c))
            total = kappa_A + kappa_dual
            assert total == Fraction(13), (
                f"kappa(Vir_{c}) + kappa(Vir_{26-c}) = {total}, expected 13"
            )

    def test_heisenberg_anti_symmetry(self):
        """kappa(H_k) + kappa(H_{-k}) = 0.

        This is the scalar sign-flip row for the open-colour model.
        The chiral Koszul dual is Sym^ch(V*), not H_{-k}.
        """
        for k_val in [1, 2, Fraction(1, 2)]:
            k = Fraction(k_val)
            kappa_A = kappa_from_self_sewing(heisenberg_ope(k))
            kappa_dual = kappa_from_self_sewing(heisenberg_ope(-k))
            assert kappa_A + kappa_dual == Fraction(0)

    def test_sl2_anti_symmetry(self):
        """kappa(sl_2_k) + kappa(sl_2_{-k-4}) = 0.

        Feigin-Frenkel: the dual level is k' = -k - 2h^v = -k - 4.
        kappa(sl_2, k) = 3(k+2)/4.
        kappa(sl_2, -k-4) = 3(-k-4+2)/4 = 3(-k-2)/4 = -3(k+2)/4.
        Sum = 0.
        """
        for k_val in [1, 2, 5, Fraction(1, 3)]:
            k = Fraction(k_val)
            k_dual = -k - 4
            kappa_A = kappa_from_self_sewing(affine_sl2_ope(k))
            kappa_dual = kappa_from_self_sewing(affine_sl2_ope(k_dual))
            assert kappa_A + kappa_dual == Fraction(0), (
                f"kappa(sl2, {k}) + kappa(sl2, {k_dual}) = "
                f"{kappa_A + kappa_dual}"
            )

    def test_virasoro_line_comparison_fixed_point(self):
        """At c = 13, the same-family Virasoro comparison is fixed.

        kappa(Vir_13) = 13/2.
        The complementarity sum is 13/2 + 13/2 = 13, consistent with AP24.
        """
        kappa = kappa_from_self_sewing(virasoro_ope(Fraction(13)))
        assert kappa == Fraction(13, 2)


# ================================================================
# 10. Beta-gamma cross-sewing
# ================================================================

class TestBetaGamma:
    """kappa(beta-gamma) = 1 from cross-sewing."""

    def test_betagamma_kappa(self):
        """kappa(betagamma) = 1 via the cross-sewing formula."""
        ope = betagamma_ope()
        kappa = kappa_betagamma_from_sewing(ope)
        assert kappa == Fraction(1), f"kappa(betagamma) = {kappa}, expected 1"

    def test_betagamma_verify(self):
        """Full verification for beta-gamma."""
        ope = betagamma_ope()
        result = verify_kappa_from_first_principles(
            ope, Fraction(1), use_cross_sewing=True,
        )
        assert result['match'], result['note']


# ================================================================
# 11. Free fermion
# ================================================================

class TestFreeFermion:
    """kappa(free fermion) = 1/2."""

    def test_fermion_kappa(self):
        """kappa(psi) = 1/2, with c = 1/2."""
        ope = free_fermion_ope()
        kappa = kappa_from_self_sewing(ope)
        # psi_{(0)} psi = 1, cyclic pairing = 2, so contribution = 1/2
        assert kappa == Fraction(1, 2), (
            f"kappa(fermion) = {kappa}, expected 1/2"
        )


# ================================================================
# 12. First-principles consistency: kappa determines F_1
# ================================================================

class TestKappaDeterminesF1:
    """F_1(A) = kappa(A) / 24 (the first free energy).

    This is the Faber-Pandharipande integral: int_{M_{1,1}} lambda_1 = 1/24.
    So F_1 = kappa * lambda_1^FP = kappa / 24.
    """

    def test_virasoro_F1(self):
        """F_1(Vir_c) = c/48."""
        for c_val in [1, 2, 26]:
            c = Fraction(c_val)
            kappa = kappa_from_self_sewing(virasoro_ope(c))
            F1 = kappa * Fraction(1, 24)
            assert F1 == c / 48, f"F_1(Vir_{c}) = {F1}, expected {c}/48"

    def test_heisenberg_F1(self):
        """F_1(H_k) = k/24."""
        k = Fraction(1)
        kappa = kappa_from_self_sewing(heisenberg_ope(k))
        F1 = kappa * Fraction(1, 24)
        assert F1 == Fraction(1, 24)

    def test_sl2_F1(self):
        """F_1(sl_2, k=1) = 9/96 = 3/32."""
        k = Fraction(1)
        kappa = kappa_from_self_sewing(affine_sl2_ope(k))
        F1 = kappa * Fraction(1, 24)
        assert F1 == Fraction(9, 96) == Fraction(3, 32)


# ================================================================
# 13. Self-consistency: OPE determines all of kappa
# ================================================================

class TestSelfConsistency:
    """The OPE data alone determines kappa: no external input needed."""

    def test_no_hardcoded_kappa(self):
        """Verify that kappa_from_self_sewing does NOT access any
        hardcoded kappa lookup table.  It uses ONLY:
          - Generator weights (from OPE data)
          - n-product values (from OPE data)
          - Cyclic pairing (from OPE data)
          - Sugawara shift (from OPE data)

        This is verified by constructing a CUSTOM algebra with made-up
        OPE data and checking that the formula applies correctly.
        """
        # Construct a custom algebra: single generator of weight 3,
        # with self-OPE leading coefficient 7/11
        custom_ope = ChiralAlgebraOPE(
            name='custom',
            generators=[GeneratorData('X', Fraction(3))],
            n_products={
                ('X', 'X'): {
                    5: Fraction(7, 11),  # X_{(5)} X = 7/11 (scalar)
                    3: 'Y',              # field-valued
                },
            },
            cyclic_pairing={('X', 'X'): Fraction(1)},
        )
        kappa = kappa_from_self_sewing(custom_ope)
        assert kappa == Fraction(7, 11), (
            f"Custom algebra: kappa = {kappa}, expected 7/11"
        )


# Import for the custom algebra test
from compute.lib.genus1_kappa_verification import ChiralAlgebraOPE, GeneratorData
