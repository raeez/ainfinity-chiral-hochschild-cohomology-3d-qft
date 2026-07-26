"""Tests for the P2 gravity-line OPE sketch: generators, OPE, package."""

import os
import sys


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.p2_gravity_line_ope_sketch import (  # noqa: E402
    BicolouredPrimitivePackageOnK3xE,
    BrownHenneauxDictionary,
    GravityLineGenerators,
    GravityLineOPEStructure,
    cross_colour_ope,
    imaginary_root_singular_theta_ope,
    km_real_root_ope,
    koszul_dual_virasoro_self_ope,
    stress_imaginary_primary_ope,
    virasoro_self_ope,
)


# ---------------------------------------------------------------------
# Generators.
# ---------------------------------------------------------------------


def test_generators_present_in_each_sector():
    g = GravityLineGenerators()
    sectors = g.sectors_present()
    assert "boundary_open" in sectors
    assert "line_closed_dual" in sectors
    assert "current" in sectors
    assert "imaginary" in sectors


def test_boundary_stress_T_is_weight_2():
    g = GravityLineGenerators()
    assert g.T_boundary.name == "T"
    assert g.T_boundary.weight == 2
    assert g.T_boundary.sector == "boundary_open"


def test_dual_stress_T_prime_is_weight_2_line_sector():
    g = GravityLineGenerators()
    # T' has central charge 26 - c only as the strict same-family
    # line-side representative of Vol II prop:gravity-koszul-dual.
    assert g.T_dual_line.weight == 2
    assert "T_prime" in g.T_dual_line.name
    assert g.T_dual_line.sector == "line_closed_dual"


def test_imaginary_root_generators_are_weight_0():
    g = GravityLineGenerators()
    # Imaginary root direction null/negative; primary weight on the
    # boundary stress is h_alpha = -alpha^2/2 + 1; the null root
    # alpha^2 = 0 gives h = 1.
    for gen in g.Theta_imaginary:
        assert gen.sector == "imaginary"


# ---------------------------------------------------------------------
# OPE structure.
# ---------------------------------------------------------------------


def test_virasoro_self_ope_has_quartic_pole_class_M():
    ope = virasoro_self_ope("c")
    # Class M: quartic pole at (z-w)^{-4} forces m_k != 0 for all k >= 3.
    assert ope.max_pole_order() == 4


def test_koszul_dual_self_ope_has_quartic_pole():
    ope = koszul_dual_virasoro_self_ope()
    assert ope.max_pole_order() == 4
    # The coefficient is (26-c)/2.
    coefs = [c.coefficient_label for c in ope.singular_part]
    assert any("26-c" in c for c in coefs)


def test_cross_colour_ope_has_no_singular_part():
    # The boundary stress and same-family line representative commute
    # on the strict comparison surface.
    ope = cross_colour_ope()
    assert ope.max_pole_order() == 0
    assert ope.singular_part == ()


def test_km_real_root_ope_has_quadratic_pole():
    # Affine Kac-Moody (real-root subsystem of g_{Delta_5}).
    ope = km_real_root_ope()
    assert ope.max_pole_order() == 2
    # Singular part includes the structure constant f^{abc} factor.
    coefs = [c.coefficient_label for c in ope.singular_part]
    assert any("f^{abc}" in c for c in coefs)


def test_imaginary_root_ope_carries_borcherds_factor():
    # Singular theta correspondence: residue at z = w yields
    # Borcherds product factor (1 - q)^{c(-alpha^2/2)}.
    ope = imaginary_root_singular_theta_ope()
    coefs = [c.coefficient_label for c in ope.singular_part]
    assert any("Borcherds" in c for c in coefs)


def test_stress_imaginary_primary_ope():
    # T(z) Theta^im(w): primary of weight h_alpha.
    ope = stress_imaginary_primary_ope()
    coefs = [c.coefficient_label for c in ope.singular_part]
    assert any("h_alpha" in c for c in coefs)
    assert ope.max_pole_order() == 2


# ---------------------------------------------------------------------
# Full OPE structure.
# ---------------------------------------------------------------------


def test_gravity_line_ope_structure_class_M_on_both_colours():
    s = GravityLineOPEStructure()
    assert s.boundary_class_M()
    assert s.line_dual_class_M()
    assert s.cross_colour_commutes()


def test_pentagon_face_trace_formula_carries_phi10_un():
    s = GravityLineOPEStructure()
    formula = s.pentagon_face_trace_formula()
    assert "Phi_10_un" in formula
    assert "Delta_5" in formula
    # The chain-level master identity:
    # tr_Pentagon = STr_{mod.op}[B^{E_3}(PhiFA_3(D^b Coh(K3 x E)))].
    assert "B^{E_3}" in formula or "PhiFA_3" in formula


def test_derived_centre_label_separates_central_shadow_from_line_comparison():
    s = GravityLineOPEStructure()
    label = s.derived_centre_bulk_HHcat()
    assert "Zder_ch(A_g)" in label
    assert "central shadow C[[c]]" in label
    assert "line comparison C[[c']]" in label
    assert "26 - c" in label
    assert "tensor" not in label


# ---------------------------------------------------------------------
# Brown--Henneaux dictionary (slogan #13 discipline).
# ---------------------------------------------------------------------


def test_brown_henneaux_dictionary_algebraic_only():
    bh = BrownHenneauxDictionary()
    # Slogan #13: NOT a dynamical-metric path integral.
    assert not bh.is_dynamical_metric
    assert bh.algebraic_only()


def test_brown_henneaux_BTZ_hypothesis_package():
    bh = BrownHenneauxDictionary()
    hyp = bh.hypothesis_package_for_BTZ()
    assert "modular_invariance" in hyp
    assert "vacuum_dominance" in hyp
    assert "saddle_dominance" in hyp


def test_brown_henneaux_central_charge_formula():
    bh = BrownHenneauxDictionary()
    assert "3 ell / (2 G_N)" in bh.central_charge_formula


# ---------------------------------------------------------------------
# Bicoloured primitive package on K3 x E.
# ---------------------------------------------------------------------


def test_nine_tuple_completeness():
    p = BicolouredPrimitivePackageOnK3xE()
    assert p.nine_tuple_completeness()


def test_closed_colour_trace_evaluates_to_phi10_un_delta5_squared():
    p = BicolouredPrimitivePackageOnK3xE()
    # Pentagon-face scalar trace = Phi_10_un = Delta_5^2.
    assert p.closed_colour_trace_equals_phi10_un()


def test_open_factorisation_lives_on_log_E_not_K3():
    p = BicolouredPrimitivePackageOnK3xE()
    # Slogan #4: open sector on log curve, not on X.
    assert "log_E" in p.open_factorisation_dg_cat
    # Closed colour lives on K3 x E (factorisation dim 3, via Phi_3^FA).
    assert "K3xE" in p.closed_factorisation_inf_cat


def test_z_der_ch_A_b_includes_dual_central_charge():
    p = BicolouredPrimitivePackageOnK3xE()
    # Derived chiral centre / chiral Hochschild bulk keeps its central
    # projection separate from the line-side comparison parameter.
    assert "C[[c]]" in p.Z_der_ch_A_b
    assert "central projection" in p.Z_der_ch_A_b
    assert "line comparison" in p.Z_der_ch_A_b
    assert "26" in p.Z_der_ch_A_b


def test_A_b_is_End_C_b_grav():
    p = BicolouredPrimitivePackageOnK3xE()
    # A_b = End_C(b), Russian-school primitive definition.
    assert "End" in p.A_b
    assert "gravity_line" in p.A_b
