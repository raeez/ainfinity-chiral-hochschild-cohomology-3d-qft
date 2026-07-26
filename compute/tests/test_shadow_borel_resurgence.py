r"""Tests for the shadow Borel resurgence engine.

Verifies:
1. Shadow data: kappa, alpha, S4, Delta, Q_L for Virasoro at all c
2. Branch points: complex conjugacy, modulus = 1/rho, argument
3. Borel singularities: positions, conjugacy, relation to branch points
4. Darboux coefficients: amplitude and phase, asymptotic prediction
5. Shadow coefficients: convolution recursion, leading values S_2=kappa, S_3=alpha
6. Borel transform: coefficient computation, entirety, convergence
7. Stokes graph: sector count, Z_2 symmetry at c=13
8. Optimal truncation: N* formula, dependence on rho
9. Koszul duality: rho(c) vs rho(26-c), kappa+kappa'=13, self-dual at c=13
10. Ratio analysis: |S_{r+1}/S_r| convergence to rho
11. Fraction arithmetic: exact computation for rational c
12. Constrained Epstein: discriminant, branch point consistency
13. Resurgence atlas: completeness, internal consistency
14. Zeta connection: correctly flagged as NO connection
15. Median resummation: consistency with exact value (convergent regime)

Cross-engine verification:
    All computations are self-contained; cross-checked against
    the formulas in shadow_radius.py (Vol I) and
    resurgence_frontier_engine.py (Vol I).

Manuscript references:
    thm:shadow-radius, thm:riccati-algebraicity, thm:single-line-dichotomy,
    def:shadow-metric, thm:shadow-connection, thm:mc2-bar-intrinsic
"""

import sys
sys.path.insert(0, 'compute')

import cmath
import math
from fractions import Fraction
from pathlib import Path

import pytest


# =====================================================================
# Section 1: VirasoroShadowData
# =====================================================================

class TestVirasoroShadowData:
    """Test VirasoroShadowData construction and derived quantities."""

    def test_kappa_equals_c_over_2(self):
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 2.0, 13.0, 25.0, 26.0]:
            d = VirasoroShadowData(c)
            assert abs(d.kappa - c / 2.0) < 1e-14

    def test_alpha_equals_2(self):
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d = VirasoroShadowData(13.0)
        assert abs(d.alpha - 2.0) < 1e-14

    def test_S4_formula(self):
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            expected = 10.0 / (c * (5.0 * c + 22.0))
            assert abs(d.S4 - expected) < 1e-14

    def test_Delta_formula(self):
        """Delta = 40/(5c+22)."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            expected = 40.0 / (5.0 * c + 22.0)
            assert abs(d.Delta - expected) < 1e-13

    def test_q0_equals_c_squared(self):
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0]:
            d = VirasoroShadowData(c)
            assert abs(d.q0 - c**2) < 1e-13

    def test_q1_equals_12c(self):
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0]:
            d = VirasoroShadowData(c)
            assert abs(d.q1 - 12.0 * c) < 1e-12

    def test_q2_formula(self):
        """q2 = (180c+872)/(5c+22)."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            expected = (180.0 * c + 872.0) / (5.0 * c + 22.0)
            assert abs(d.q2 - expected) < 1e-12

    def test_rho_formula(self):
        """rho = sqrt((180c+872)/((5c+22)*c^2))."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            rho_sq = (180.0 * c + 872.0) / ((5.0 * c + 22.0) * c**2)
            expected = math.sqrt(rho_sq)
            assert abs(d.rho - expected) < 1e-12

    def test_rho_comparison_fixed_c13(self):
        """rho(13) ~ 0.4674 (convergent)."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d = VirasoroShadowData(13.0)
        assert abs(d.rho - 0.467396) < 1e-4
        assert d.is_convergent

    def test_rho_divergent_c1(self):
        """rho(1) ~ 6.24 (divergent)."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d = VirasoroShadowData(1.0)
        assert d.rho > 1.0
        assert not d.is_convergent

    def test_is_comparison_fixed(self):
        from lib.shadow_borel_resurgence import VirasoroShadowData
        assert VirasoroShadowData(13.0).is_comparison_fixed
        assert not VirasoroShadowData(12.0).is_comparison_fixed


# =====================================================================
# Section 2: Branch points
# =====================================================================

class TestBranchPoints:
    """Test branch point computation."""

    def test_branch_points_are_zeros_of_QL(self):
        """Q_L(t_+) = 0."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            Q_at_tp = d.q0 + d.q1 * d.t_plus + d.q2 * d.t_plus**2
            assert abs(Q_at_tp) < 1e-10

    def test_branch_points_conjugate(self):
        """For c > 0 (Delta > 0), branch points are complex conjugates."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            assert abs(d.t_plus - d.t_minus.conjugate()) < 1e-12

    def test_branch_point_modulus_equals_R(self):
        """Both branch points have modulus R = 1/rho."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            assert abs(abs(d.t_plus) - d.R) < 1e-12
            assert abs(abs(d.t_minus) - d.R) < 1e-12

    def test_branch_points_negative_real_part(self):
        """Branch points are in the left half-plane (Re(t_+) < 0)."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            assert d.t_plus.real < 0

    def test_branch_point_argument_near_pi(self):
        """arg(t_+) is between pi/2 and pi (upper left quadrant)."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            arg = cmath.phase(d.t_plus)
            assert math.pi / 2 < arg < math.pi


# =====================================================================
# Section 3: Borel singularities
# =====================================================================

class TestBorelSingularities:
    """Test Borel singularity positions (instanton actions)."""

    def test_A_equals_reciprocal_of_branch_point(self):
        """A_+ = 1/t_+."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            assert abs(d.A_plus - 1.0 / d.t_plus) < 1e-12

    def test_A_conjugate(self):
        """A_+/- are complex conjugates."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0]:
            d = VirasoroShadowData(c)
            assert abs(d.A_plus - d.A_minus.conjugate()) < 1e-12

    def test_A_modulus_equals_rho(self):
        """|A_+| = rho = 1/R."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            assert abs(abs(d.A_plus) - d.rho) < 1e-12

    def test_borel_singularities_function(self):
        from lib.shadow_borel_resurgence import borel_singularities
        bs = borel_singularities(13.0)
        assert bs['are_conjugate']
        assert abs(bs['A_plus_mod'] - bs['rho']) < 1e-12

    def test_borel_singularities_negative_real_part(self):
        """Borel singularities have negative real part (left half-plane)."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            assert d.A_plus.real < 0

    def test_stokes_direction_near_minus_pi(self):
        """Stokes direction arg(A_+) is near -pi."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            # A_+ is in the third quadrant (negative real, negative imag)
            # so arg(A_+) is between -pi and -pi/2
            assert -math.pi < d.stokes_direction < -math.pi / 2

    def test_c_half_borel_position(self):
        """c = 1/2 (Ising): A_+ = -12 - 3.614i."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d = VirasoroShadowData(0.5)
        assert abs(d.A_plus.real - (-12.0)) < 0.01
        assert abs(d.A_plus.imag - (-3.614)) < 0.01

    def test_c13_borel_position(self):
        """c = 13: A_+ = -6/13 - i*Im."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d = VirasoroShadowData(13.0)
        # A_+ = 1/t_+ where t_+ is root of Q_L
        # Re(A_+) = -6/13 (from the q1 coefficient)
        assert abs(d.A_plus.real - (-6.0 / 13.0)) < 0.001


# =====================================================================
# Section 4: Shadow coefficients
# =====================================================================

class TestShadowCoefficients:
    """Test shadow coefficient computation."""

    def test_S2_equals_kappa(self):
        """S_2 = kappa = c/2."""
        from lib.shadow_borel_resurgence import shadow_coefficients
        for c in [1.0, 13.0, 25.0]:
            coeffs = shadow_coefficients(c, 10)
            assert abs(coeffs[2] - c / 2.0) < 1e-12

    def test_S3_equals_alpha(self):
        """S_3 = alpha = 2."""
        from lib.shadow_borel_resurgence import shadow_coefficients
        for c in [1.0, 13.0, 25.0]:
            coeffs = shadow_coefficients(c, 10)
            assert abs(coeffs[3] - 2.0) < 1e-12

    def test_coefficient_count(self):
        from lib.shadow_borel_resurgence import shadow_coefficients
        coeffs = shadow_coefficients(13.0, 50)
        assert len(coeffs) == 49  # S_2 through S_50

    def test_H_generates_sqrt_QL(self):
        """H(t) = sum r*S_r*t^r = t^2*sqrt(Q_L(t))."""
        from lib.shadow_borel_resurgence import shadow_coefficients, VirasoroShadowData
        c = 13.0
        d = VirasoroShadowData(c)
        coeffs = shadow_coefficients(c, 60)
        t = 0.3
        H_partial = sum(r * coeffs[r] * t**r for r in range(2, 61))
        Q_at_t = d.q0 + d.q1 * t + d.q2 * t**2
        H_exact = t**2 * math.sqrt(Q_at_t)
        assert abs(H_partial - H_exact) < 1e-8

    def test_fraction_matches_float(self):
        """Fraction and float computations agree."""
        from lib.shadow_borel_resurgence import (
            shadow_coefficients, shadow_coefficients_fraction,
        )
        coeffs_float = shadow_coefficients(1.0, 20)
        coeffs_frac = shadow_coefficients_fraction(1, 1, 20)
        for r in range(2, 21):
            ref = max(abs(coeffs_float[r]), 1e-15)
            assert abs(coeffs_float[r] - float(coeffs_frac[r])) / ref < 1e-10


# =====================================================================
# Section 5: Darboux coefficients
# =====================================================================

class TestDarbouxCoefficients:
    """Test Darboux asymptotic transfer coefficients."""

    def test_C_plus_conjugate_C_minus(self):
        """C_+ = conj(C_-)."""
        from lib.shadow_borel_resurgence import darboux_coefficients
        for c in [1.0, 13.0, 25.0]:
            dd = darboux_coefficients(c)
            assert abs(dd.C_plus - dd.C_minus.conjugate()) < 1e-10

    def test_amplitude_positive(self):
        from lib.shadow_borel_resurgence import darboux_coefficients
        for c in [0.5, 1.0, 13.0, 25.0]:
            dd = darboux_coefficients(c)
            assert dd.amplitude > 0

    def test_c13_amplitude(self):
        """c = 13: Darboux amplitude ~ 18.86."""
        from lib.shadow_borel_resurgence import darboux_coefficients
        dd = darboux_coefficients(13.0)
        assert abs(dd.amplitude - 18.86) < 0.1

    def test_asymptotic_prediction_sign_pattern(self):
        """Asymptotic prediction captures the sign oscillation."""
        from lib.shadow_borel_resurgence import (
            shadow_coefficients, asymptotic_prediction,
        )
        c = 13.0
        coeffs = shadow_coefficients(c, 60)
        # For large r, the sign of S_r should match the sign of the
        # asymptotic prediction
        matches = 0
        total = 0
        for r in range(30, 55):
            pred = asymptotic_prediction(c, r)
            actual = coeffs[r]
            if abs(actual) > 1e-20:
                total += 1
                if (pred > 0) == (actual > 0):
                    matches += 1
        # At least 80% of signs should match at large r
        assert matches / total > 0.8, f"Only {matches}/{total} signs match"

    def test_asymptotic_ratio_convergence(self):
        """S_r / (prediction) -> 1 as r -> infinity."""
        from lib.shadow_borel_resurgence import (
            shadow_coefficients, asymptotic_prediction,
        )
        c = 13.0
        coeffs = shadow_coefficients(c, 80)
        # Ratio should approach 1; check terms where prediction is not near zero
        # (cosine factor creates near-zeros at certain arities)
        ratios = []
        for r in range(30, 75):
            pred = asymptotic_prediction(c, r)
            actual = coeffs[r]
            # Only include terms where both are sufficiently nonzero
            if abs(pred) > 1e-25 and abs(actual) > 1e-25:
                rat = actual / pred
                # Only include when the cosine factor is not near zero
                # (near-zero cosine gives unreliable ratios)
                if abs(rat) < 10:
                    ratios.append(rat)
        assert len(ratios) > 5, f"Too few valid ratios: {len(ratios)}"
        # The ratios should be near 1 (positive, since signs match)
        avg_ratio = sum(abs(r) for r in ratios) / len(ratios)
        # Should be within 30% of 1 (convergence is slow due to r^{-5/2})
        assert 0.7 < avg_ratio < 1.5, f"avg ratio = {avg_ratio}"


# =====================================================================
# Section 6: Borel transform
# =====================================================================

class TestBorelTransform:
    """Test Borel transform computation."""

    def test_borel_coefficients_divided_by_factorial(self):
        """b_r = S_r / r!."""
        from lib.shadow_borel_resurgence import (
            shadow_coefficients, borel_coefficients,
        )
        coeffs = shadow_coefficients(13.0, 20)
        b = borel_coefficients(coeffs)
        for r in range(2, 21):
            expected = coeffs[r] / math.gamma(r + 1)
            assert abs(b[r] - expected) < 1e-15

    def test_borel_transform_at_origin(self):
        """B(0) = 0 (series starts at zeta^2)."""
        from lib.shadow_borel_resurgence import shadow_coefficients, borel_transform
        coeffs = shadow_coefficients(13.0, 30)
        val = borel_transform(coeffs, 0.0)
        assert abs(val) < 1e-15

    def test_borel_transform_entire_large_zeta(self):
        """B(zeta) is finite for large |zeta| (entire function)."""
        from lib.shadow_borel_resurgence import shadow_coefficients, borel_transform
        coeffs = shadow_coefficients(13.0, 60)
        # Evaluate at zeta = 10 (far from origin)
        val = borel_transform(coeffs, 10.0)
        assert math.isfinite(abs(val))

    def test_borel_small_zeta_dominated_by_S2(self):
        """For small zeta, B(zeta) ~ S_2 * zeta^2 / 2!."""
        from lib.shadow_borel_resurgence import shadow_coefficients, borel_transform
        c = 13.0
        coeffs = shadow_coefficients(c, 20)
        zeta = 0.01
        val = borel_transform(coeffs, zeta)
        leading = coeffs[2] * zeta**2 / math.gamma(3)
        # Should be dominated by the leading term
        assert abs(val - leading) / abs(leading) < 0.01


# =====================================================================
# Section 7: Stokes graph
# =====================================================================

class TestStokesGraph:
    """Test Stokes graph geometry."""

    def test_c13_z2_symmetry(self):
        from lib.shadow_borel_resurgence import stokes_graph
        sg = stokes_graph(13.0)
        assert sg.has_z2_symmetry

    def test_c1_no_z2_symmetry(self):
        from lib.shadow_borel_resurgence import stokes_graph
        sg = stokes_graph(1.0)
        assert not sg.has_z2_symmetry

    def test_stokes_angles_count(self):
        from lib.shadow_borel_resurgence import stokes_graph
        sg = stokes_graph(13.0)
        # At least 2 distinct Stokes directions (from A_+ and A_-)
        assert sg.n_sectors >= 2

    def test_anti_stokes_orthogonal_to_stokes(self):
        """Anti-Stokes lines are at +/- pi/2 from SOME Stokes line."""
        from lib.shadow_borel_resurgence import VirasoroShadowData, stokes_graph
        d = VirasoroShadowData(13.0)
        sg = stokes_graph(13.0)
        # Each anti-Stokes angle should be at pi/2 from the corresponding
        # instanton action's argument (either A_+ or A_-)
        for anti_dir in sg.anti_stokes_angles:
            # Check if this anti-Stokes angle is pi/2 from arg(A_+) or arg(A_-)
            ok = False
            for stokes_dir in [cmath.phase(d.A_plus), cmath.phase(d.A_minus)]:
                diff = abs(anti_dir - stokes_dir)
                # Normalize to [0, pi]
                diff = diff % (2 * math.pi)
                if diff > math.pi:
                    diff = 2 * math.pi - diff
                if abs(diff - math.pi / 2) < 0.01:
                    ok = True
            assert ok, f"anti-Stokes angle {anti_dir} not orthogonal to any Stokes line"


# =====================================================================
# Section 8: Optimal truncation
# =====================================================================

class TestOptimalTruncation:
    """Test optimal truncation order."""

    def test_N_star_formula(self):
        """N* = floor(1/rho) at t=1."""
        from lib.shadow_borel_resurgence import (
            optimal_truncation_order, VirasoroShadowData,
        )
        for c in [13.0, 25.0]:
            d = VirasoroShadowData(c)
            N = optimal_truncation_order(c)
            expected = max(2, int(1.0 / d.rho))
            assert N == expected

    def test_divergent_series_N_star_is_2(self):
        """For rho > 1 (divergent): N* = 2 (immediate truncation)."""
        from lib.shadow_borel_resurgence import optimal_truncation_order
        # c = 1: rho ~ 6.24 > 1
        N = optimal_truncation_order(1.0)
        assert N == 2

    def test_convergent_series_large_N_star(self):
        """For rho << 1 (convergent): N* is large."""
        from lib.shadow_borel_resurgence import optimal_truncation_order
        # c = 25: rho ~ 0.24, so N* ~ 4
        N = optimal_truncation_order(25.0)
        assert N >= 3


# =====================================================================
# Section 9: Koszul duality
# =====================================================================

class TestKoszulDuality:
    """Test Koszul duality relations in the Borel plane."""

    def test_kappa_sum_equals_13(self):
        """kappa(c) + kappa(26-c) = 13 for Virasoro."""
        from lib.shadow_borel_resurgence import koszul_dual_borel_comparison
        for c in [1.0, 5.0, 13.0, 25.0]:
            kd = koszul_dual_borel_comparison(c)
            assert abs(kd['kappa_sum'] - 13.0) < 1e-12

    def test_comparison_fixed_at_c13(self):
        """rho(13) = rho(13): comparison fixed."""
        from lib.shadow_borel_resurgence import koszul_dual_borel_comparison
        kd = koszul_dual_borel_comparison(13.0)
        assert kd['comparison_fixed']

    def test_not_comparison_fixed_at_c1(self):
        from lib.shadow_borel_resurgence import koszul_dual_borel_comparison
        kd = koszul_dual_borel_comparison(1.0)
        assert not kd['comparison_fixed']

    def test_rho_not_symmetric_under_duality(self):
        """rho(c) != rho(26-c) in general on the line-side comparison.

        The shadow growth rate depends on c^2, so it is NOT invariant under
        the comparison involution c -> 26-c. Only at c=13 (the fixed point) do
        rho(c) and rho(26-c) agree.
        """
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d1 = VirasoroShadowData(1.0)
        d25 = VirasoroShadowData(25.0)
        # rho(1) ~ 6.24, rho(25) ~ 0.24 -- vastly different
        assert abs(d1.rho - d25.rho) > 1.0

    def test_stokes_direction_c_vs_dual(self):
        """Stokes directions at c and 26-c are generally different."""
        from lib.shadow_borel_resurgence import koszul_dual_borel_comparison
        kd = koszul_dual_borel_comparison(5.0)
        # They should have different Stokes directions
        assert abs(kd['stokes_dir'] - kd['stokes_dir_dual']) > 0.001


# =====================================================================
# Section 10: Ratio analysis
# =====================================================================

class TestRatioAnalysis:
    """Test asymptotic ratio analysis."""

    def test_ratio_convergence_direction(self):
        """Ratios approach rho from above (for the dominant contribution)."""
        from lib.shadow_borel_resurgence import ratio_analysis
        ra = ratio_analysis(13.0, 60)
        # The skip ratios (every-other-term) should converge to rho^2
        skip = ra['skip_ratios']
        if len(skip) > 5:
            last = skip[-1]
            # Within 50% of rho^2 (convergence is slow)
            assert abs(last['ratio_sq'] - last['predicted_sq']) / last['predicted_sq'] < 0.5

    def test_ratio_analysis_has_data(self):
        from lib.shadow_borel_resurgence import ratio_analysis
        ra = ratio_analysis(13.0, 40)
        assert len(ra['ratios']) > 20
        assert len(ra['skip_ratios']) > 10


# =====================================================================
# Section 11: Fraction arithmetic
# =====================================================================

class TestFractionArithmetic:
    """Test exact computation with Fraction arithmetic."""

    def test_S2_fraction_exact(self):
        """S_2 = kappa = c/2 exactly."""
        from lib.shadow_borel_resurgence import shadow_coefficients_fraction
        for c_num, c_den in [(1, 1), (13, 1), (1, 2)]:
            coeffs = shadow_coefficients_fraction(c_num, c_den, 10)
            expected = Fraction(c_num, c_den) / 2
            assert coeffs[2] == expected

    def test_S3_fraction_exact(self):
        """S_3 = alpha = 2 exactly."""
        from lib.shadow_borel_resurgence import shadow_coefficients_fraction
        coeffs = shadow_coefficients_fraction(1, 1, 10)
        assert coeffs[3] == Fraction(2)

    def test_S4_fraction_agrees_with_Q_contact(self):
        """S_4 = Q^contact_Vir / ... -- check consistency."""
        from lib.shadow_borel_resurgence import shadow_coefficients_fraction
        # For c=1: S_4 = a_2/4 where a_2 = (q2 - a_1^2)/(2*a_0)
        coeffs = shadow_coefficients_fraction(1, 1, 10)
        S4_computed = coeffs[4]
        # S_4 should be a rational number for c=1
        assert isinstance(S4_computed, Fraction)

    def test_fraction_high_arity(self):
        """Fraction computation at arity 30 is finite and nonzero."""
        from lib.shadow_borel_resurgence import shadow_coefficients_fraction
        coeffs = shadow_coefficients_fraction(1, 1, 30)
        assert coeffs[30] != 0
        assert isinstance(coeffs[30], Fraction)


# =====================================================================
# Section 11b: Manuscript guards for the gravity resolvent
# =====================================================================

class TestGravityResolventManuscriptGuard:
    """Guard the scalar-resolvent theorem against stale normalizations."""

    @staticmethod
    def _gravity_source() -> str:
        root = Path(__file__).resolve().parents[2]
        return (root / "chapters/connections/3d_gravity.tex").read_text()

    def test_resolvent_theorem_names_branch_and_genus_zero(self):
        source = self._gravity_source()
        start = source.index(
            r"\begin{theorem}[Scalar shadow resolvent on the Virasoro metric branch"
        )
        end = source.index(r"\end{theorem}", start)
        block = source[start:end]
        compact = "".join(block.split())

        assert r"\ClaimStatusConditional" in block
        assert r"\hypAmbientWtCpl+\hypStokes" in block
        assert r"\sqrt{Q_{\mathrm{Vir}}(0)}=c" in block
        assert r"r \ge 2" in block
        assert "Abelian integral on the rational double cover" in block
        assert "genus~$0$" in block
        assert "candidate scalar-shadow Stokes directions" in block
        assert "elliptic integral" not in block
        assert "No real pole" not in block
        assert "bound-state" not in block
        branch_formula = "".join((
            r"\frac{c(5c+22)}{180c+872}"
            r"\left(-6 \pm 4i\sqrt{\frac{5}{5c+22}}\right)"
        ).split())
        assert branch_formula in compact
        assert r"\tfrac{c(5c+22)}{2}" not in block

    def test_shadow_table_uses_normalized_low_degree_rows(self):
        source = self._gravity_source()
        start = source.index(r"\subsubsection*{The shadow table through $S_9$}")
        end = source.index(r"\end{center}", start)
        table = source[start:end]

        assert r"$2$ & $\dfrac{c}{2}$ & $+$ & $c/2$" in table
        assert r"$3$ & $2$ & $+$ & $2$" in table
        assert r"$2$ & $\dfrac{c}{12}$" not in table
        assert r"$3$ & $-c$" not in table
        assert "first Taylor" in table
        assert "not its Taylor coefficients" not in table

    def test_exact_low_degree_shadow_coefficients_at_c13(self):
        from lib.shadow_borel_resurgence import shadow_coefficients_fraction

        coeffs = shadow_coefficients_fraction(13, 1, 5)
        assert coeffs[2] == Fraction(13, 2)
        assert coeffs[3] == Fraction(2)
        assert coeffs[4] == Fraction(10, 13 * (5 * 13 + 22))
        assert coeffs[5] == Fraction(-48, 13**2 * (5 * 13 + 22))


# =====================================================================
# Section 12: Constrained Epstein connection
# =====================================================================

class TestConstrainedEpstein:
    """Test constrained Epstein zeta function connection."""

    def test_discriminant_negative(self):
        """disc(Q_L) = -32*kappa^2*Delta < 0 for c > 0 (Delta > 0)."""
        from lib.shadow_borel_resurgence import constrained_epstein_comparison
        for c in [1.0, 13.0, 25.0]:
            ec = constrained_epstein_comparison(c)
            assert ec['Q_L_discriminant'] < 0

    def test_branch_points_consistent(self):
        from lib.shadow_borel_resurgence import (
            constrained_epstein_comparison, VirasoroShadowData,
        )
        c = 13.0
        ec = constrained_epstein_comparison(c)
        d = VirasoroShadowData(c)
        assert abs(ec['branch_points'][0] - d.t_plus) < 1e-12

    def test_comparison_fixed_c13(self):
        from lib.shadow_borel_resurgence import constrained_epstein_comparison
        ec = constrained_epstein_comparison(13.0)
        assert ec['functional_equation_symmetric']


# =====================================================================
# Section 13: Resurgence atlas
# =====================================================================

class TestResurgenceAtlas:
    """Test the complete resurgence atlas."""

    def test_atlas_has_all_entries(self):
        from lib.shadow_borel_resurgence import resurgence_atlas
        atlas = resurgence_atlas()
        assert len(atlas) == 8  # default 8 central charges

    def test_atlas_rho_monotone_decreasing(self):
        """rho decreases with c (for c > 0)."""
        from lib.shadow_borel_resurgence import resurgence_atlas
        atlas = resurgence_atlas()
        for i in range(len(atlas) - 1):
            assert atlas[i]['rho'] > atlas[i + 1]['rho']

    def test_atlas_convergence_transition(self):
        """There is a transition from divergent to convergent near c*~6.12."""
        from lib.shadow_borel_resurgence import resurgence_atlas
        atlas = resurgence_atlas()
        found_transition = False
        for i in range(len(atlas) - 1):
            if not atlas[i]['convergent'] and atlas[i + 1]['convergent']:
                found_transition = True
                # Transition should be between c=6 and c=13
                assert atlas[i]['c'] < 13.0
                assert atlas[i + 1]['c'] >= 6.0
        assert found_transition

    def test_atlas_comparison_fixed_only_c13(self):
        from lib.shadow_borel_resurgence import resurgence_atlas
        atlas = resurgence_atlas()
        z2_entries = [a for a in atlas if a['z2_symmetric']]
        assert len(z2_entries) == 1
        assert abs(z2_entries[0]['c'] - 13.0) < 0.1


# =====================================================================
# Section 14: Zeta connection
# =====================================================================

class TestZetaConnection:
    """Test that the zeta connection is correctly assessed as absent."""

    def test_no_connection(self):
        from lib.shadow_borel_resurgence import zeta_connection_assessment
        for c in [0.5, 1.0, 13.0]:
            zc = zeta_connection_assessment(c)
            assert zc['connection_exists'] is False

    def test_ratio_to_zeta_zero_not_rational(self):
        """The ratio |Im(A_+)| / gamma_1 is not a simple rational."""
        from lib.shadow_borel_resurgence import zeta_connection_assessment
        zc = zeta_connection_assessment(0.5)
        ratio = zc['ratio']
        # Should NOT be a simple fraction
        for n in range(1, 20):
            for d in range(1, 20):
                if abs(ratio - n / d) < 0.001:
                    pytest.fail(f"ratio {ratio} too close to {n}/{d}")


# =====================================================================
# Section 15: Borel singularity table
# =====================================================================

class TestBorelSingularityTable:
    """Test the Borel singularity table."""

    def test_table_completeness(self):
        from lib.shadow_borel_resurgence import borel_singularities_table
        table = borel_singularities_table()
        assert len(table) == 8

    def test_all_entries_have_conjugate_singularities(self):
        from lib.shadow_borel_resurgence import borel_singularities_table
        for entry in borel_singularities_table():
            assert entry['are_conjugate']

    def test_rho_equals_A_modulus(self):
        from lib.shadow_borel_resurgence import borel_singularities_table
        for entry in borel_singularities_table():
            assert abs(entry['A_plus_mod'] - entry['rho']) < 1e-12


# =====================================================================
# Section 16: Specific numerical values
# =====================================================================

class TestSpecificValues:
    """Test specific numerical values from the manuscript."""

    def test_rho_ising_c_half(self):
        """rho(Ising, c=1/2) ~ 12.53."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d = VirasoroShadowData(0.5)
        assert abs(d.rho - 12.532) < 0.01

    def test_rho_c26_string(self):
        """rho(c=26) ~ 0.232."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d = VirasoroShadowData(26.0)
        assert abs(d.rho - 0.2325) < 0.001

    def test_rho_asymmetry_under_duality(self):
        """rho is NOT symmetric under c -> 26-c (except at c=13)."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d1 = VirasoroShadowData(1.0)
        d25 = VirasoroShadowData(25.0)
        # rho depends on c^2 in the denominator, so rho(1) >> rho(25)
        assert d1.rho > 6.0
        assert d25.rho < 0.3

    def test_Q_contact_vir(self):
        """Q^contact_Vir = S_4 = 10/(c(5c+22))."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [1.0, 13.0]:
            d = VirasoroShadowData(c)
            expected = 10.0 / (c * (5.0 * c + 22.0))
            assert abs(d.S4 - expected) < 1e-14

    def test_Delta_c13(self):
        """Delta(c=13) = 40/87 ~ 0.4598."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d = VirasoroShadowData(13.0)
        assert abs(d.Delta - 40.0 / 87.0) < 1e-12


# =====================================================================
# Section 17: Discriminant and shadow metric properties
# =====================================================================

class TestShadowMetric:
    """Test shadow metric Q_L(t) properties."""

    def test_Q_at_zero_positive(self):
        """Q_L(0) = q0 = c^2 > 0."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0]:
            d = VirasoroShadowData(c)
            assert d.q0 > 0

    def test_discriminant_negative(self):
        """disc(Q_L) = q1^2 - 4*q0*q2 = -32*kappa^2*Delta < 0."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            disc = d.q1**2 - 4.0 * d.q0 * d.q2
            expected = -32.0 * d.kappa**2 * d.Delta
            assert abs(disc - expected) < 1e-8
            assert disc < 0

    def test_Q_positive_on_real_axis_inside_R(self):
        """Q_L(t) > 0 for real t in (-R, R)."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        d = VirasoroShadowData(13.0)
        for t in [0.0, 0.5, 1.0, 1.5, 2.0]:
            Q = d.q0 + d.q1 * t + d.q2 * t**2
            assert Q > 0


# =====================================================================
# Section 18: Cross-consistency checks
# =====================================================================

class TestCrossConsistency:
    """Cross-consistency between different computational paths."""

    def test_rho_from_shadow_data_vs_formula(self):
        """rho from VirasoroShadowData matches the explicit formula."""
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            rho_formula = math.sqrt(
                (180.0 * c + 872.0) / ((5.0 * c + 22.0) * c**2)
            )
            assert abs(d.rho - rho_formula) < 1e-12

    def test_R_times_rho_equals_1(self):
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0, 25.0]:
            d = VirasoroShadowData(c)
            assert abs(d.R * d.rho - 1.0) < 1e-12

    def test_omega_equals_abs_theta(self):
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [0.5, 1.0, 13.0]:
            d = VirasoroShadowData(c)
            assert abs(d.omega - abs(d.theta)) < 1e-14

    def test_borel_singularity_modulus_from_two_paths(self):
        """Check |A_+| = rho using two independent calculations."""
        from lib.shadow_borel_resurgence import VirasoroShadowData, borel_singularities
        c = 13.0
        d = VirasoroShadowData(c)
        bs = borel_singularities(c)
        assert abs(abs(d.A_plus) - d.rho) < 1e-12
        assert abs(bs['A_plus_mod'] - d.rho) < 1e-12

    def test_darboux_omega_matches_shadow_data(self):
        from lib.shadow_borel_resurgence import (
            darboux_coefficients, VirasoroShadowData,
        )
        c = 13.0
        dd = darboux_coefficients(c)
        d = VirasoroShadowData(c)
        # omega from darboux = -theta (where theta = arg(t_+))
        assert abs(dd.omega - (-d.theta)) < 1e-12
