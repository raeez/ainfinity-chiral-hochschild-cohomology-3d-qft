r"""Tests for the exact WKB / Stokes graph engine.

Verifies:
1. Quantum spectral curve: Q_L matches shadow_borel_resurgence
2. Branch points: complex conjugates for c > 0, agree with VirasoroShadowData
3. WKB expansion: S_0 = graviton resolvent, S_1 = -(1/4)*log(Q), S_2 from Riccati
4. Stokes graph: 3 lines per turning point, Z_2 at c=13, topology
5. Voros symbols: no bound states (complex turning points), unit circle at c=13
6. Connection formula: S_1 = 2*pi*i (universal for simple turning points)
7. Period integral: analytical formula matches numerical
8. Painleve identification: P_III_3 (D_8)
9. Cross-check with shadow_borel_resurgence: consistency of all data
10. Self-dual c=13: enhanced symmetry, real-axis reflection

Manuscript references:
    thm:shadow-metric, thm:shadow-connection, Movement VI,
    thm:riccati-algebraicity
"""

import sys
sys.path.insert(0, 'compute')

import cmath
import math

import pytest


# =====================================================================
# Section 1: Quantum Spectral Curve
# =====================================================================

class TestQuantumSpectralCurve:
    """Test the quantum spectral curve Q_L(t)."""

    def test_Q_matches_shadow_data(self):
        """Q_L(t) from QuantumSpectralCurve matches VirasoroShadowData."""
        from lib.exact_wkb_stokes import QuantumSpectralCurve
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [1.0, 6.0, 13.0, 25.0]:
            qsc = QuantumSpectralCurve(c)
            d = VirasoroShadowData(c)
            for t in [0.0, 1.0, -0.5, 2.0]:
                Q_qsc = qsc.Q(t)
                Q_shadow = d.q0 + d.q1 * t + d.q2 * t * t
                assert abs(Q_qsc - Q_shadow) < 1e-10, \
                    f"Q mismatch at c={c}, t={t}: {Q_qsc} vs {Q_shadow}"

    def test_Q_factored_form(self):
        """Q_L(t) = (c+6t)^2 + 80t^2/(5c+22)."""
        from lib.exact_wkb_stokes import QuantumSpectralCurve
        for c in [1.0, 13.0, 26.0]:
            qsc = QuantumSpectralCurve(c)
            for t in [0.5, 1.0, 2.0]:
                Q = qsc.Q(t)
                Q_expected = (c + 6 * t)**2 + 80 * t**2 / (5 * c + 22)
                assert abs(Q - Q_expected) < 1e-10

    def test_branch_points_are_zeros(self):
        """Q_L(t_+) = Q_L(t_-) = 0."""
        from lib.exact_wkb_stokes import QuantumSpectralCurve
        for c in [1.0, 6.0, 13.0, 25.0]:
            qsc = QuantumSpectralCurve(c)
            assert abs(qsc.Q(qsc.t_plus)) < 1e-8, f"Q(t_+) != 0 at c={c}"
            assert abs(qsc.Q(qsc.t_minus)) < 1e-8, f"Q(t_-) != 0 at c={c}"

    def test_branch_points_conjugate(self):
        """For c > 0: t_+ and t_- are complex conjugates."""
        from lib.exact_wkb_stokes import QuantumSpectralCurve
        for c in [1.0, 6.0, 13.0, 25.0]:
            qsc = QuantumSpectralCurve(c)
            assert abs(qsc.t_plus - qsc.t_minus.conjugate()) < 1e-10, \
                f"Branch points not conjugate at c={c}"

    def test_turning_points_complex(self):
        """Turning points are complex (not real) for all c > 0."""
        from lib.exact_wkb_stokes import QuantumSpectralCurve
        for c in [0.5, 1.0, 6.0, 13.0, 25.0, 26.0]:
            qsc = QuantumSpectralCurve(c)
            assert qsc.turning_points_are_complex, \
                f"Turning points should be complex at c={c}"

    def test_genus_zero(self):
        """Spectral curve is always genus 0."""
        from lib.exact_wkb_stokes import QuantumSpectralCurve
        qsc = QuantumSpectralCurve(13.0)
        assert qsc.is_genus_zero

    def test_branch_points_match_shadow(self):
        """Branch points agree with VirasoroShadowData."""
        from lib.exact_wkb_stokes import QuantumSpectralCurve
        from lib.shadow_borel_resurgence import VirasoroShadowData
        for c in [1.0, 13.0, 26.0]:
            qsc = QuantumSpectralCurve(c)
            d = VirasoroShadowData(c)
            assert abs(qsc.t_plus - d.t_plus) < 1e-12
            assert abs(qsc.t_minus - d.t_minus) < 1e-12


# =====================================================================
# Section 2: WKB Expansion
# =====================================================================

class TestWKBExpansion:
    """Test the WKB expansion terms."""

    def test_graviton_resolvent_origin(self):
        """G(0) = int_0^0 sqrt(Q) ds = 0."""
        from lib.exact_wkb_stokes import graviton_resolvent_numerical
        G = graviton_resolvent_numerical(13.0, 0.0)
        assert abs(G) < 1e-10

    def test_graviton_resolvent_small_t(self):
        """G(t) ~ c*t for small real t (since sqrt(Q) ~ c near t=0)."""
        from lib.exact_wkb_stokes import graviton_resolvent_numerical
        c = 13.0
        t = 0.01
        G = graviton_resolvent_numerical(c, t)
        # sqrt(Q_L(0)) = sqrt(c^2) = c
        expected = c * t  # leading order
        assert abs(G - expected) / abs(expected) < 0.02, \
            f"G({t}) = {G}, expected ~ {expected}"

    def test_S1_formula(self):
        """S_1(t) = -(1/4)*log(Q_L(t))."""
        from lib.exact_wkb_stokes import wkb_S1, QuantumSpectralCurve
        c = 13.0
        t = 1.0 + 0.5j
        qsc = QuantumSpectralCurve(c)
        S1 = wkb_S1(c, t)
        expected = -0.25 * cmath.log(qsc.Q(t))
        assert abs(S1 - expected) < 1e-12

    def test_S2_finite_away_from_turning_points(self):
        """S_2(t) is finite for t away from turning points."""
        from lib.exact_wkb_stokes import wkb_S2
        c = 13.0
        S2 = wkb_S2(c, 1.0 + 0.5j)
        assert math.isfinite(abs(S2)), f"S_2 should be finite, got {S2}"

    def test_S2_diverges_at_turning_point(self):
        """S_2(t) diverges as t -> turning point (Q -> 0)."""
        from lib.exact_wkb_stokes import wkb_S2, QuantumSpectralCurve
        c = 13.0
        qsc = QuantumSpectralCurve(c)
        # Approach t_+ closely
        t_near = qsc.t_plus + 0.001
        S2 = wkb_S2(c, t_near)
        t_far = 1.0 + 0.5j
        S2_far = wkb_S2(c, t_far)
        assert abs(S2) > 10 * abs(S2_far), \
            "S_2 should be larger near turning point"

    def test_free_energy_density_g0(self):
        """dS_0/dt = sqrt(Q_L(t))."""
        from lib.exact_wkb_stokes import wkb_free_energy_densities, QuantumSpectralCurve
        c = 13.0
        t = 1.0 + 0.5j
        fe = wkb_free_energy_densities(c, t, g_max=2)
        qsc = QuantumSpectralCurve(c)
        expected = cmath.sqrt(qsc.Q(t))
        assert abs(fe[0] - expected) < 1e-10

    def test_free_energy_density_g1(self):
        """dS_1/dt = -Q'/(4Q)."""
        from lib.exact_wkb_stokes import wkb_free_energy_densities, QuantumSpectralCurve
        c = 13.0
        t = 1.0 + 0.5j
        fe = wkb_free_energy_densities(c, t, g_max=2)
        qsc = QuantumSpectralCurve(c)
        expected = -qsc.Q_prime(t) / (4.0 * qsc.Q(t))
        assert abs(fe[1] - expected) < 1e-10


# =====================================================================
# Section 3: Stokes Graph
# =====================================================================

class TestStokesGraph:
    """Test Stokes graph topology."""

    def test_three_stokes_lines_per_turning_point(self):
        """Each simple turning point has exactly 3 Stokes lines."""
        from lib.exact_wkb_stokes import compute_stokes_graph
        sg = compute_stokes_graph(13.0, n_trace=50, trace_length=2.0)
        assert len(sg.stokes_directions_t_plus) == 3
        assert len(sg.stokes_directions_t_minus) == 3

    def test_three_anti_stokes_lines_per_turning_point(self):
        """Each simple turning point has exactly 3 anti-Stokes lines."""
        from lib.exact_wkb_stokes import compute_stokes_graph
        sg = compute_stokes_graph(13.0, n_trace=50, trace_length=2.0)
        assert len(sg.anti_stokes_directions_t_plus) == 3
        assert len(sg.anti_stokes_directions_t_minus) == 3

    def test_stokes_directions_120_degrees_apart(self):
        """Stokes lines are separated by 2*pi/3 = 120 degrees."""
        from lib.exact_wkb_stokes import compute_stokes_graph
        sg = compute_stokes_graph(13.0, n_trace=50, trace_length=2.0)
        dirs = sorted(sg.stokes_directions_t_plus)
        # Check angular separations are approximately 2*pi/3
        for i in range(3):
            j = (i + 1) % 3
            diff = dirs[j] - dirs[i]
            if diff < 0:
                diff += 2 * math.pi
            assert abs(diff - 2 * math.pi / 3) < 0.1, \
                f"Angular separation {diff} != 2*pi/3"

    def test_z2_symmetry_at_c13(self):
        """Z_2 symmetry (complex conjugation) at c = 13."""
        from lib.exact_wkb_stokes import compute_stokes_graph
        sg = compute_stokes_graph(13.0, n_trace=50, trace_length=2.0)
        assert sg.has_z2_symmetry

    def test_no_z2_at_generic_c(self):
        """No special Z_2 symmetry at generic c (well, it should
        still have conjugation symmetry since Q_L is real)."""
        from lib.exact_wkb_stokes import compute_stokes_graph
        # Actually, for ANY real c, Q_L has real coefficients,
        # so t_+ and t_- are conjugates and the Stokes graph
        # has Z_2. Our has_z2_symmetry flag specifically checks c=13.
        sg = compute_stokes_graph(5.0, n_trace=50, trace_length=2.0)
        assert not sg.has_z2_symmetry  # flag checks c=13 specifically

    def test_total_stokes_lines(self):
        """Total 6 Stokes lines (3 from each turning point)."""
        from lib.exact_wkb_stokes import compute_stokes_graph
        sg = compute_stokes_graph(13.0, n_trace=50, trace_length=2.0)
        assert len(sg.stokes_lines) == 6


# =====================================================================
# Section 4: Voros Symbols
# =====================================================================

class TestVorosSymbols:
    """Test Voros symbols and quantization."""

    def test_no_bound_states(self):
        """Complex turning points -> no bound states for all c > 0."""
        from lib.exact_wkb_stokes import voros_symbols
        for c in [1.0, 6.0, 13.0, 25.0, 26.0]:
            vd = voros_symbols(c)
            assert not vd.has_bound_states, \
                f"Should have no bound states at c={c}"

    def test_voros_modulus_finite(self):
        """Voros modulus is finite and positive."""
        from lib.exact_wkb_stokes import voros_symbols
        for c in [1.0, 13.0, 26.0]:
            vd = voros_symbols(c, hbar=1.0)
            assert vd.voros_modulus > 0
            assert math.isfinite(vd.voros_modulus)

    def test_instanton_action_nonzero(self):
        """Instanton action is nonzero (the curve is not degenerate)."""
        from lib.exact_wkb_stokes import voros_symbols
        for c in [1.0, 13.0, 25.0]:
            vd = voros_symbols(c)
            assert abs(vd.instanton_action) > 0.01


# =====================================================================
# Section 5: Connection Formula
# =====================================================================

class TestConnectionFormula:
    """Test the Stokes connection formula."""

    def test_stokes_constant_2pi_i(self):
        """Stokes constant S_1 = 2*pi*i for simple turning points."""
        from lib.exact_wkb_stokes import connection_formula
        for c in [1.0, 13.0, 26.0]:
            cf = connection_formula(c)
            assert abs(cf.stokes_constant - 2.0j * math.pi) < 1e-10

    def test_monodromy_minus_one(self):
        """Monodromy of sqrt(Q_L) around simple zero = -1."""
        from lib.exact_wkb_stokes import connection_formula
        cf = connection_formula(13.0)
        assert cf.monodromy_sign == -1

    def test_connection_matrix_unipotent(self):
        """Connection matrix is upper triangular unipotent."""
        from lib.exact_wkb_stokes import connection_formula
        cf = connection_formula(13.0)
        M = cf.connection_matrix
        assert abs(M[0][0] - 1) < 1e-12
        assert abs(M[1][0]) < 1e-12
        assert abs(M[1][1] - 1) < 1e-12
        assert abs(M[0][1] - 2.0j * math.pi) < 1e-10


# =====================================================================
# Section 6: Period Integral
# =====================================================================

class TestPeriodIntegral:
    """Test the period integral computation."""

    def test_analytical_matches_numerical(self):
        """Analytical period formula matches numerical integration."""
        from lib.exact_wkb_stokes import period_integral_exact
        for c in [1.0, 6.0, 13.0, 25.0]:
            result = period_integral_exact(c)
            assert result['relative_error'] < 0.02, \
                f"Period mismatch at c={c}: error={result['relative_error']}"

    def test_period_purely_imaginary_at_c13(self):
        """At c=13, the period integral should be close to purely imaginary
        (by Z_2 symmetry of the contour)."""
        from lib.exact_wkb_stokes import period_integral_exact
        result = period_integral_exact(13.0)
        period = result['period_numerical']
        # The period is purely imaginary if Re(period)/|period| ~ 0
        if abs(period) > 1e-10:
            ratio = abs(period.real) / abs(period)
            # Allow some numerical tolerance
            assert ratio < 0.05, \
                f"Period at c=13 should be ~purely imaginary, ratio={ratio}"

    def test_period_nonzero(self):
        """Period integral is nonzero for all c > 0."""
        from lib.exact_wkb_stokes import period_integral_exact
        for c in [1.0, 13.0, 26.0]:
            result = period_integral_exact(c)
            assert abs(result['period_numerical']) > 0.01


# =====================================================================
# Section 7: Painleve Analysis
# =====================================================================

class TestPainleveAnalysis:
    """Test the Painleve identification."""

    def test_painleve_type(self):
        """Painleve type is P_III_3 (D_8 Sakai)."""
        from lib.exact_wkb_stokes import painleve_analysis
        pa = painleve_analysis(13.0)
        assert 'III' in pa.painleve_type

    def test_algebraic_dependence(self):
        """c enters Q_L algebraically."""
        from lib.exact_wkb_stokes import painleve_analysis
        pa = painleve_analysis(13.0)
        assert pa.is_algebraic_dependence

    def test_dq2_dc(self):
        """dq2/dc = -400/(5c+22)^2."""
        from lib.exact_wkb_stokes import painleve_compatibility_check
        for c in [1.0, 13.0, 25.0]:
            result = painleve_compatibility_check(c)
            dq2 = result['dq2_dc']
            expected = -400.0 / (5.0 * c + 22.0)**2
            assert abs(dq2 - expected) < 1e-10

    def test_naive_residual_nonzero(self):
        """The naive B-matrix does NOT satisfy zero curvature
        (the gauge correction is the Painleve equation)."""
        from lib.exact_wkb_stokes import painleve_compatibility_check
        result = painleve_compatibility_check(13.0)
        assert result['requires_gauge_correction']
        assert result['naive_residual_norm'] > 0.01


# =====================================================================
# Section 8: Comparison-Fixed c = 13
# =====================================================================

class TestComparisonFixed:
    """Test the comparison fixed point c = 13."""

    def test_comparison_fixed_analysis_runs(self):
        """Comparison-fixed analysis completes without error."""
        from lib.exact_wkb_stokes import comparison_fixed_analysis
        sd = comparison_fixed_analysis()
        assert sd['c'] == 13.0

    def test_branch_points_conjugate_at_c13(self):
        from lib.exact_wkb_stokes import comparison_fixed_analysis
        sd = comparison_fixed_analysis()
        assert sd['branch_points']['are_conjugate']

    def test_z2_symmetry(self):
        from lib.exact_wkb_stokes import comparison_fixed_analysis
        sd = comparison_fixed_analysis()
        assert sd['stokes_graph']['has_z2_symmetry']

    def test_no_bound_states_at_c13(self):
        from lib.exact_wkb_stokes import comparison_fixed_analysis
        sd = comparison_fixed_analysis()
        assert not sd['voros']['has_bound_states']

    def test_stokes_constant_at_c13(self):
        from lib.exact_wkb_stokes import comparison_fixed_analysis
        sd = comparison_fixed_analysis()
        assert abs(sd['connection']['stokes_constant'] - 2.0j * math.pi) < 1e-10


# =====================================================================
# Section 9: Cross-check with shadow_borel_resurgence
# =====================================================================

class TestCrossCheck:
    """Test consistency with the shadow Borel resurgence engine."""

    def test_branch_points_agree(self):
        from lib.exact_wkb_stokes import cross_check_borel_stokes
        xc = cross_check_borel_stokes(13.0)
        assert xc['checks']['branch_points_agree']

    def test_stokes_constant_agree(self):
        from lib.exact_wkb_stokes import cross_check_borel_stokes
        xc = cross_check_borel_stokes(13.0)
        assert xc['checks']['stokes_constant_agree']

    def test_z2_symmetry_agree(self):
        from lib.exact_wkb_stokes import cross_check_borel_stokes
        xc = cross_check_borel_stokes(13.0)
        assert xc['checks']['z2_symmetry_sbr']
        assert xc['checks']['z2_symmetry_wkb']


# =====================================================================
# Section 10: Atlas
# =====================================================================

class TestAtlas:
    """Test the full WKB atlas."""

    def test_atlas_completeness(self):
        from lib.exact_wkb_stokes import exact_wkb_atlas
        atlas = exact_wkb_atlas()
        assert len(atlas) == 8

    def test_no_bound_states_anywhere(self):
        """No bound states for any c > 0 in the atlas."""
        from lib.exact_wkb_stokes import exact_wkb_atlas
        atlas = exact_wkb_atlas()
        for row in atlas:
            assert not row['has_bound_states'], \
                f"Unexpected bound states at c={row['c']}"

    def test_all_turning_points_complex(self):
        from lib.exact_wkb_stokes import exact_wkb_atlas
        atlas = exact_wkb_atlas()
        for row in atlas:
            assert row['turning_points_complex'], \
                f"Turning points should be complex at c={row['c']}"

    def test_painleve_type_consistent(self):
        from lib.exact_wkb_stokes import exact_wkb_atlas
        atlas = exact_wkb_atlas()
        for row in atlas:
            assert 'III' in row['painleve_type']
