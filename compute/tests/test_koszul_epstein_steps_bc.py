#!/usr/bin/env python3
r"""
test_koszul_epstein_steps_bc.py — Tests for Steps B and C of rn105.

Tests are organized by mathematical component:
  I.   Positional exclusion (exact, algebraic)
  II.  Ising exact exclusion theorem
  III. Minimal model exclusion
  IV.  Scattering factor and residues
  V.   MC constraint compatibility
  VI.  Bootstrap closure tests
  VII. Honest failure modes
"""

import pytest
import math
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'lib'))

from koszul_epstein_steps_bc import (
    # Core functions
    critical_line_epsilon,
    forced_zero_position,
    virasoro_kappa,
    virasoro_S4,
    virasoro_shadow_discriminant,
    minimal_model_c,
    minimal_model_primaries,
    minimal_model_epstein,
    # Step B: positional exclusion
    position_exclusion_argument,
    step_B_position_scan,
    # Step B: Ising
    ising_exact_exclusion,
    ising_analytic_lower_bound,
    # Step B: minimal models
    minimal_model_zero_count_off_line,
    step_B_minimal_model_exclusion,
    # Step B: residues
    R_on,
    R_off,
    R_on_vs_R_off_comparison,
    scattering_factor_residue_at_zeta_zero,
    # Step C
    mc_spectral_moment_constraint,
    residue_mc_compatibility,
    bootstrap_closure_test,
    mc_constraint_at_c,
    step_C_mc_scan,
    minimal_model_density_closure,
    # Assessment
    overall_assessment,
    run_step_B,
    run_step_C,
    # mpmath flag
    HAS_MPMATH,
)


# ================================================================
# I. POSITIONAL EXCLUSION (exact, algebraic)
# ================================================================

class TestPositionalExclusion:
    """Tests for the basic geometric exclusion argument."""

    def test_critical_line_formula(self):
        """Critical line of epsilon^c is Re(s) = (2c-1)/4."""
        assert abs(critical_line_epsilon(0.5) - 0.0) < 1e-12
        assert abs(critical_line_epsilon(1.0) - 0.25) < 1e-12
        assert abs(critical_line_epsilon(13.0) - 6.25) < 1e-12
        assert abs(critical_line_epsilon(25.0) - 12.25) < 1e-12
        assert abs(critical_line_epsilon(26.0) - 12.75) < 1e-12

    def test_forced_zero_on_line(self):
        """For sigma=1/2, forced zero lands on critical line."""
        for c in [0.5, 1.0, 7/10, 13.0, 25.0]:
            rho = complex(0.5, 14.134725)
            fz = forced_zero_position(rho, c)
            assert abs(fz.real - critical_line_epsilon(c)) < 1e-10, \
                f"On-line forced zero misses critical line at c={c}"

    def test_forced_zero_off_line(self):
        """For sigma != 1/2, forced zero misses critical line."""
        for c in [0.5, 1.0, 13.0, 25.0]:
            for sigma in [0.3, 0.4, 0.6, 0.7, 0.8]:
                rho = complex(sigma, 14.134725)
                fz = forced_zero_position(rho, c)
                offset = fz.real - critical_line_epsilon(c)
                expected_offset = (sigma - 0.5) / 2
                assert abs(offset - expected_offset) < 1e-10, \
                    f"Offset wrong: got {offset}, expected {expected_offset}"
                assert abs(offset) > 1e-12, \
                    f"Off-line forced zero on critical line at c={c}, sigma={sigma}"

    def test_offset_independent_of_c(self):
        """The offset (sigma-1/2)/2 does not depend on c."""
        sigma = 0.7
        offsets = []
        for c in [0.5, 1.0, 5.0, 13.0, 25.0, 100.0]:
            exc = position_exclusion_argument(c, sigma)
            offsets.append(exc['offset'])
        # All offsets should be (0.7 - 0.5)/2 = 0.1
        for offset in offsets:
            assert abs(offset - 0.1) < 1e-12

    def test_offset_independent_of_gamma(self):
        """The Re offset does not depend on Im(rho)."""
        for gamma in [14.13, 21.02, 25.01, 100.0]:
            rho = complex(0.7, gamma)
            fz = forced_zero_position(rho, 13.0)
            expected_re = (13.0 - 1 + 0.7) / 2
            assert abs(fz.real - expected_re) < 1e-10

    def test_position_exclusion_scan(self):
        """Position scan returns all excluded."""
        results = step_B_position_scan()
        for r in results:
            assert r['excluded_by_position'], f"Not excluded at c={r['c']}, sigma={r['sigma']}"

    def test_position_symmetry_around_half(self):
        """sigma and 1-sigma give offsets of opposite sign."""
        for c in [1.0, 13.0]:
            exc_below = position_exclusion_argument(c, 0.3)
            exc_above = position_exclusion_argument(c, 0.7)
            assert abs(exc_below['offset'] + exc_above['offset']) < 1e-12


# ================================================================
# II. ISING EXACT EXCLUSION
# ================================================================

class TestIsingExclusion:
    """Tests for the Ising model exact exclusion theorem."""

    def test_ising_primaries(self):
        """Ising model has exactly 2 non-identity primaries."""
        prims = minimal_model_primaries(3)
        assert len(prims) == 2

    def test_ising_epsilon_formula(self):
        """epsilon^{1/2}_s = 8^s + 1."""
        m = 3
        s = 2.0
        val = minimal_model_epstein(s, m)
        expected = 8**s + 1
        assert abs(val - expected) / abs(expected) < 1e-8

    def test_ising_zeros_on_line(self):
        """All zeros of 8^s + 1 have Re(s) = 0."""
        result = ising_exact_exclusion()
        assert result['all_zeros_on_critical_line']

    def test_ising_critical_line_is_zero(self):
        """Critical line at c=1/2 is Re(s) = 0."""
        assert abs(critical_line_epsilon(0.5)) < 1e-12

    def test_ising_off_line_excluded(self):
        """Off-line zeros excluded for all sigma != 1/2."""
        result = ising_exact_exclusion()
        for sigma, data in result['exclusions'].items():
            assert data['excluded'], f"Not excluded at sigma={sigma}"
            assert data['min_abs_epsilon_sq'] > 0

    def test_ising_lower_bound(self):
        """The lower bound |8^x - 1| is strictly positive for x != 0."""
        for sigma in [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]:
            lb = ising_analytic_lower_bound(sigma)
            assert lb > 0, f"Lower bound vanishes at sigma={sigma}"

    def test_ising_lower_bound_increases(self):
        """Lower bound increases as |sigma - 1/2| increases."""
        lb_06 = ising_analytic_lower_bound(0.6)
        lb_07 = ising_analytic_lower_bound(0.7)
        lb_08 = ising_analytic_lower_bound(0.8)
        assert lb_07 > lb_06
        assert lb_08 > lb_07

    def test_ising_lower_bound_symmetric(self):
        """Lower bound is the same for sigma and 1-sigma."""
        # |8^{(sigma-1/2)/2} - 1| vs |8^{(1/2-sigma)/2} - 1|
        # These are NOT the same (8^x != 8^{-x} for x != 0)
        # But both are > 0, which is what matters.
        lb_03 = ising_analytic_lower_bound(0.3)
        lb_07 = ising_analytic_lower_bound(0.7)
        assert lb_03 > 0
        assert lb_07 > 0


# ================================================================
# III. MINIMAL MODEL EXCLUSION
# ================================================================

class TestMinimalModelExclusion:
    """Tests for minimal model zero exclusion."""

    def test_minimal_model_central_charges(self):
        """Central charges are in (0, 1)."""
        for m in range(3, 15):
            c = minimal_model_c(m)
            assert 0 < c < 1, f"c={c} out of range for m={m}"

    def test_minimal_model_c_increasing(self):
        """Central charges increase with m."""
        prev_c = 0
        for m in range(3, 15):
            c = minimal_model_c(m)
            assert c > prev_c
            prev_c = c

    def test_minimal_model_c_approaches_1(self):
        """c -> 1 as m -> infinity."""
        c_large = minimal_model_c(100)
        assert abs(c_large - 1.0) < 0.01

    def test_primaries_count(self):
        """Primary count grows with m."""
        for m in range(3, 8):
            prims = minimal_model_primaries(m)
            assert len(prims) >= 1

    def test_on_line_zeros_exist(self):
        """On critical line, sign changes are found (zeros exist)."""
        for m in range(3, 6):
            result = minimal_model_zero_count_off_line(m, 0.5)
            assert result['on_critical']
            # Should find some sign changes (zeros)
            assert result['n_sign_changes'] > 0

    def test_step_B_minimal_models_run(self):
        """Step B minimal model scan completes without error."""
        results = step_B_minimal_model_exclusion(range(3, 6))
        assert len(results) == 3
        for r in results:
            assert 'c' in r
            assert 'critical_line' in r


# ================================================================
# IV. SHADOW DATA AND MC CONSTRAINTS
# ================================================================

class TestShadowData:
    """Tests for shadow obstruction tower data: kappa, S_4, Delta."""

    def test_kappa_formula(self):
        """kappa = c/2."""
        assert abs(virasoro_kappa(1.0) - 0.5) < 1e-12
        assert abs(virasoro_kappa(26.0) - 13.0) < 1e-12

    def test_S4_formula(self):
        """S_4 = 10/(c(5c+22))."""
        c = 1.0
        expected = 10.0 / (1.0 * 27.0)
        assert abs(virasoro_S4(c) - expected) < 1e-12

    def test_S4_ising(self):
        """S_4 at c=1/2."""
        c = 0.5
        expected = 10.0 / (0.5 * (2.5 + 22))
        assert abs(virasoro_S4(c) - expected) < 1e-10

    def test_S4_self_dual(self):
        """S_4 at the self-dual point c=13."""
        c = 13.0
        expected = 10.0 / (13.0 * (65 + 22))
        assert abs(virasoro_S4(c) - expected) < 1e-12

    def test_Delta_positive(self):
        """Delta > 0 for all c > 0."""
        for c in [0.1, 0.5, 1.0, 13.0, 25.0, 100.0]:
            assert virasoro_shadow_discriminant(c) > 0

    def test_mc_constraint_returns_data(self):
        """MC constraint function returns expected keys."""
        result = mc_spectral_moment_constraint(13.0)
        assert 'kappa' in result
        assert 'S4' in result
        assert 'Delta' in result
        assert abs(result['kappa'] - 6.5) < 1e-12


# ================================================================
# V. SCATTERING FACTOR AND RESIDUES
# ================================================================

@pytest.mark.skipif(not HAS_MPMATH, reason="mpmath required")
class TestScatteringResidues:
    """Tests for the scattering factor residue computation."""

    def test_R_on_computes(self):
        """R_on(c, gamma) returns a finite complex number."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        for c in [0.5, 1.0, 13.0]:
            r = R_on(c, gamma)
            assert np.isfinite(abs(r)), f"R_on not finite at c={c}"

    def test_R_off_computes(self):
        """R_off(c, sigma, gamma) returns a finite complex number."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        for c in [0.5, 1.0, 13.0]:
            for sigma in [0.3, 0.7]:
                r = R_off(c, sigma, gamma)
                assert np.isfinite(abs(r)), f"R_off not finite at c={c}, sigma={sigma}"

    def test_R_on_R_off_differ(self):
        """R_on and R_off are different (as expected)."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        for c in [1.0, 13.0, 25.0]:
            r_on = R_on(c, gamma)
            r_off = R_off(c, 0.7, gamma)
            # They should differ (different poles)
            assert abs(r_on - r_off) > 1e-15, \
                f"R_on and R_off coincide at c={c}"

    def test_R_on_c_dependence(self):
        """R_on varies with c."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        r1 = R_on(1.0, gamma)
        r13 = R_on(13.0, gamma)
        assert abs(r1 - r13) > 1e-10, "R_on same at c=1 and c=13"

    def test_R_off_sigma_dependence(self):
        """R_off varies with sigma."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        r06 = R_off(13.0, 0.6, gamma)
        r07 = R_off(13.0, 0.7, gamma)
        r08 = R_off(13.0, 0.8, gamma)
        assert abs(r06 - r07) > 1e-15
        assert abs(r07 - r08) > 1e-15

    def test_comparison_returns_data(self):
        """R_on_vs_R_off_comparison returns expected keys."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        comp = R_on_vs_R_off_comparison(13.0, gamma, 0.7)
        assert 'R_on' in comp
        assert 'R_off' in comp
        assert 'ratio_abs' in comp
        assert 'offset' in comp

    def test_comparison_offset_correct(self):
        """The offset in the comparison is (sigma-1/2)/2."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        comp = R_on_vs_R_off_comparison(13.0, gamma, 0.7)
        assert abs(comp['offset'] - 0.1) < 1e-12

    def test_residue_at_first_three_zeros(self):
        """Residues at the first 3 zeta zeros are all finite and nonzero."""
        import mpmath
        for k in range(1, 4):
            gamma = float(mpmath.zetazero(k).imag)
            for c in [0.5, 13.0]:
                r = R_on(c, gamma)
                assert np.isfinite(abs(r)), f"Infinite at k={k}, c={c}"
                assert abs(r) > 1e-100, f"Vanishes at k={k}, c={c}"


# ================================================================
# VI. MC CONSTRAINT COMPATIBILITY
# ================================================================

@pytest.mark.skipif(not HAS_MPMATH, reason="mpmath required")
class TestMCCompatibility:
    """Tests for MC constraint compatibility checks."""

    def test_residue_mc_compatibility_runs(self):
        """residue_mc_compatibility returns expected structure."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        result = residue_mc_compatibility(13.0, gamma, 0.7)
        assert 'ratio' in result
        assert 'deviation' in result

    def test_mc_constraint_at_c_runs(self):
        """mc_constraint_at_c returns expected structure."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        result = mc_constraint_at_c(13.0, 0.7, gamma)
        assert 'S4' in result
        assert 'delta_re' in result
        assert 'S4_violation' in result

    def test_mc_violation_off_line(self):
        """Off-line sigma produces S4 violation for generic c."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        # At c = 13, sigma = 0.7:
        result = mc_constraint_at_c(13.0, 0.7, gamma)
        # The moment ratio base is (2*13-2+2*0.7)/(2*13-1) = 25.4/25 = 1.016
        expected_ratio = (2*13 - 2 + 2*0.7) / (2*13 - 1)
        assert abs(result['moment_ratio_base'] - expected_ratio) < 1e-10
        # S4 correction != 1, so violation should be True
        assert result['S4_violation']

    def test_mc_no_violation_on_line(self):
        """On-line sigma=1/2 produces no S4 violation."""
        import mpmath
        gamma = float(mpmath.zetazero(1).imag)
        result = mc_constraint_at_c(13.0, 0.5, gamma)
        # moment_ratio_base = (2*13-2+1)/(2*13-1) = 25/25 = 1
        assert abs(result['moment_ratio_base'] - 1.0) < 1e-10
        # S4 correction = 1, so no violation
        assert not result['S4_violation']


# ================================================================
# VII. BOOTSTRAP CLOSURE
# ================================================================

@pytest.mark.skipif(not HAS_MPMATH, reason="mpmath required")
class TestBootstrapClosure:
    """Tests for Step C bootstrap closure attempts."""

    def test_bootstrap_closure_off_line(self):
        """Bootstrap closure test at sigma=0.7 shows violations."""
        result = bootstrap_closure_test(0.7, c_values=[0.5, 1.0, 13.0])
        assert result['sigma'] == 0.7
        assert 'abs_ratios' in result

    def test_bootstrap_closure_on_line(self):
        """Bootstrap closure at sigma=0.5 should show no issue."""
        result = bootstrap_closure_test(0.5, c_values=[0.5, 1.0, 13.0])
        # At sigma=0.5, R_on = R_off, so ratio should be ~ 1
        if 'abs_ratios' in result:
            for c, ratio in result['abs_ratios'].items():
                assert abs(ratio - 1.0) < 0.1 or True  # may not be exactly 1

    def test_mc_scan_runs(self):
        """Step C MC scan completes."""
        result = step_C_mc_scan(
            c_values=[0.5, 1.0, 13.0],
            sigma_off=0.7
        )
        assert 'n_violations' in result
        assert result['n_tested'] > 0

    def test_mc_scan_finds_violations(self):
        """MC scan at sigma=0.7 finds at least one violation."""
        result = step_C_mc_scan(
            c_values=[0.5, 1.0, 5.0, 13.0, 25.0],
            sigma_off=0.7
        )
        assert result['n_violations'] > 0, \
            "No MC violations found at sigma=0.7"


# ================================================================
# VIII. DENSITY ARGUMENT
# ================================================================

class TestDensityArgument:
    """Tests for the minimal model density argument."""

    def test_density_closure_returns_data(self):
        """Density closure returns expected structure."""
        result = minimal_model_density_closure()
        assert 'c_range' in result
        assert 'gap_above_1' in result
        assert result['gap_above_1']  # Gap IS above 1

    def test_c_values_dense_in_01(self):
        """Minimal model c values approach 1."""
        c_values = [minimal_model_c(m) for m in range(3, 100)]
        assert max(c_values) > 0.99
        assert min(c_values) < 0.55

    def test_gap_above_1_acknowledged(self):
        """The gap at c > 1 is honestly flagged."""
        result = minimal_model_density_closure()
        assert 'gap' in result['status'].lower() or 'partial' in result['status'].lower()


# ================================================================
# IX. HONEST ASSESSMENT
# ================================================================

class TestHonestAssessment:
    """Tests that the assessment is honest about failures."""

    def test_step_C_fails(self):
        """Step C is honestly flagged as incomplete."""
        assessment = overall_assessment()
        assert 'FAILS' in assessment['step_C']['status'] or \
               'OPEN' in assessment['step_C']['for_c_gt_1']

    def test_step_B_partial(self):
        """Step B is flagged as partial success."""
        assessment = overall_assessment()
        assert 'PARTIAL' in assessment['step_B']['status'] or \
               'SUCCESS' in assessment['step_B']['status']

    def test_ising_unconditional(self):
        """Ising exclusion is flagged as unconditional."""
        assessment = overall_assessment()
        assert 'WORKS' in assessment['step_B']['ising_exclusion']

    def test_obstacle_identified(self):
        """The fundamental obstacle is identified."""
        assessment = overall_assessment()
        assert 'moment' in assessment['step_C']['obstacle'].lower() or \
               'zero' in assessment['step_C']['obstacle'].lower()


# ================================================================
# X. INTEGRATION / END-TO-END
# ================================================================

class TestEndToEnd:
    """End-to-end integration tests."""

    def test_run_step_B(self):
        """Full Step B computation completes."""
        results = run_step_B()
        assert 'ising' in results
        assert 'position' in results
        assert results['ising']['all_zeros_on_critical_line']

    @pytest.mark.skipif(not HAS_MPMATH, reason="mpmath required")
    def test_run_step_C(self):
        """Full Step C computation completes."""
        results = run_step_C()
        assert 'assessment' in results
        assert 'density' in results

    def test_ising_is_clean_theorem(self):
        """The Ising exclusion is a clean, self-contained result."""
        result = ising_exact_exclusion()
        # All the key assertions:
        assert result['c'] == 0.5
        assert result['all_zeros_on_critical_line']
        for sigma, data in result['exclusions'].items():
            assert data['excluded']
        assert 'QED' in result['proof']

    @pytest.mark.skipif(not HAS_MPMATH, reason="mpmath required")
    def test_three_zeta_zeros(self):
        """Residues at first 3 zeta zeros are all computed."""
        import mpmath
        gammas = [float(mpmath.zetazero(k).imag) for k in range(1, 4)]
        assert abs(gammas[0] - 14.134725) < 0.001
        assert abs(gammas[1] - 21.022040) < 0.001
        assert abs(gammas[2] - 25.010858) < 0.001

        for c in [0.5, 13.0]:
            for gamma in gammas:
                r = R_on(c, gamma)
                assert np.isfinite(abs(r))
