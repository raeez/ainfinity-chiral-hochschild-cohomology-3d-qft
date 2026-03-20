"""Tests for planted-forest obstruction theory.

Verifies the genuinely modular corrections that distinguish the
modular bar from the genus-0 bar: FM boundary types, planted-forest
coefficient algebra, MC dictionary, residue splitting, cubic gauge
triviality, and D^2 = 0.

Ground truth:
  - higher_genus_modular_koszul.tex (Vol I): G_pf, d_pf, D^2=0
  - nonlinear_modular_shadows.tex (Vol I): shadow tower
  - fm3_planted_forest_synthesis.tex (Vol II): FM_3 strata
  - Mok25: log FM tropicalization, thm:ambient-d-squared-zero
  - modular_bar.py (Vol I): PlantedForestType, fm3_planted_forest_types

Tier 1 (structural): all tests are self-certifying against known counts.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from fractions import Fraction

from lib.planted_forest_obstruction import (
    RootedTree,
    PlantedForest,
    make_leaf,
    make_binary_tree,
    make_corolla,
    FMBoundaryType,
    fm3_boundary_types,
    fm3_planted_forest_strata,
    fm4_boundary_types,
    planted_forest_coefficient,
    planted_forest_coefficient_from_tree,
    d_pf_arity3,
    d_pf_arity4,
    mc_dictionary_algebraic,
    mc_dictionary_geometric,
    mc_dictionary_verify,
    residue_splitting,
    cubic_gauge_triviality_check,
    genuinely_planted_forest_cubic,
    d_squared_zero_check,
    codimension_table,
    mok_tropicalization_verify,
    boundary_operator_matrix,
)


# ===================================================================
# PLANTED FOREST / ROOTED TREE REPRESENTATION
# ===================================================================

class TestRootedTree:
    """Basic rooted tree construction and properties."""

    def test_leaf_construction(self):
        """A leaf has 1 leaf, 0 internal vertices, 0 edges."""
        leaf = make_leaf(1)
        assert leaf.is_leaf
        assert leaf.num_leaves == 1
        assert leaf.num_internal_vertices == 0
        assert leaf.num_edges == 0

    def test_binary_tree_two_leaves(self):
        """Binary tree with 2 leaves has 1 internal vertex, 2 edges."""
        tree = make_binary_tree(make_leaf(1), make_leaf(2))
        assert not tree.is_leaf
        assert tree.num_leaves == 2
        assert tree.num_internal_vertices == 1
        assert tree.num_edges == 2

    def test_binary_tree_three_leaves(self):
        """Binary tree with 3 leaves: ((1,2),3) has 2 internal nodes, 4 edges."""
        left = make_binary_tree(make_leaf(1), make_leaf(2))
        tree = make_binary_tree(left, make_leaf(3))
        assert tree.num_leaves == 3
        assert tree.num_internal_vertices == 2
        # 4 edges: root->left_subtree, root->leaf_3, left->leaf_1, left->leaf_2
        assert tree.num_edges == 4

    def test_corolla_construction(self):
        """A corolla on n leaves has 1 internal vertex, n edges."""
        for n in [2, 3, 4, 5]:
            tree = make_corolla(list(range(1, n + 1)))
            assert tree.num_leaves == n
            assert tree.num_internal_vertices == 1
            assert tree.num_edges == n

    def test_leaf_labels(self):
        """all_leaves returns the correct label set."""
        tree = make_binary_tree(make_leaf(3), make_leaf(7))
        assert tree.all_leaves == frozenset({3, 7})

    def test_three_leaf_labels(self):
        """Nested binary tree preserves leaf labels."""
        left = make_binary_tree(make_leaf(1), make_leaf(2))
        tree = make_binary_tree(left, make_leaf(3))
        assert tree.all_leaves == frozenset({1, 2, 3})

    def test_total_genus_zero(self):
        """Default genus is 0 at all vertices."""
        tree = make_binary_tree(make_leaf(1), make_leaf(2))
        assert tree.total_genus == 0

    def test_total_genus_nonzero(self):
        """Genus accumulates over internal vertices."""
        left = RootedTree(children=(make_leaf(1), make_leaf(2)), genus=1)
        tree = RootedTree(children=(left, make_leaf(3)), genus=2)
        assert tree.total_genus == 3  # genus 1 + genus 2

    def test_automorphism_leaf(self):
        """Aut of a single leaf is trivial."""
        assert make_leaf(1).automorphism_order() == 1

    def test_automorphism_binary_asymmetric(self):
        """Aut of ((1,2),3): root children are non-isomorphic, but
        the inner node (1,2) has Aut = 2 (swap leaves). Total |Aut| = 2."""
        left = make_binary_tree(make_leaf(1), make_leaf(2))
        tree = make_binary_tree(left, make_leaf(3))
        # Root: children are (0:.,.) and . -- different types, so no root swap.
        # Inner binary node: two leaf children are isomorphic, giving factor 2.
        assert tree.automorphism_order() == 2

    def test_automorphism_binary_symmetric(self):
        """Aut of a symmetric binary tree with two leaves is 2."""
        # (1, 2): the two leaves are at the same level but labeled,
        # so the tree automorphism swaps them -> |Aut| = 2
        # (when ignoring labels, which _tree_canonical_form does)
        tree = make_binary_tree(make_leaf(1), make_leaf(2))
        assert tree.automorphism_order() == 2


# ===================================================================
# PLANTED FOREST
# ===================================================================

class TestPlantedForest:
    """Planted forest: collection of rooted trees."""

    def test_single_tree_forest(self):
        """A forest with one tree."""
        tree = make_binary_tree(make_leaf(1), make_leaf(2))
        forest = PlantedForest(trees=(tree,))
        assert forest.num_trees == 1
        assert forest.total_leaves == 2

    def test_two_tree_forest(self):
        """A forest with two disjoint trees."""
        t1 = make_binary_tree(make_leaf(1), make_leaf(2))
        t2 = make_binary_tree(make_leaf(3), make_leaf(4))
        forest = PlantedForest(trees=(t1, t2))
        assert forest.num_trees == 2
        assert forest.total_leaves == 4
        assert forest.all_leaves == frozenset({1, 2, 3, 4})

    def test_forest_coefficient_single_corolla(self):
        """Coefficient of a single corolla tree: 1/|Aut|."""
        # Corolla on 3 leaves: all 3 children are leaves (isomorphic)
        # |Aut| = 3! = 6
        tree = make_corolla([1, 2, 3])
        coeff = planted_forest_coefficient_from_tree(tree)
        assert coeff == Fraction(1, 6)

    def test_forest_coefficient_binary(self):
        """Coefficient of a binary tree with 2 leaves."""
        tree = make_binary_tree(make_leaf(1), make_leaf(2))
        # |Aut| = 2 (swap the two leaves)
        coeff = planted_forest_coefficient_from_tree(tree)
        assert coeff == Fraction(1, 2)

    def test_forest_coefficient_three_leaf_binary(self):
        """((1,2),3) has |Aut| = 2 (inner swap), so coefficient = 1/2."""
        left = make_binary_tree(make_leaf(1), make_leaf(2))
        tree = make_binary_tree(left, make_leaf(3))
        coeff = planted_forest_coefficient_from_tree(tree)
        assert coeff == Fraction(1, 2)


# ===================================================================
# FM_3 BOUNDARY TYPES
# ===================================================================

class TestFM3BoundaryTypes:
    """The 7 boundary types of FM_3."""

    def test_total_count(self):
        """FM_3 has exactly 7 boundary types."""
        types = fm3_boundary_types()
        assert len(types) == 7

    def test_codim1_count(self):
        """FM_3 has 4 codimension-1 boundary strata."""
        types = fm3_boundary_types()
        codim1 = [t for t in types if t.codimension == 1]
        assert len(codim1) == 4

    def test_codim2_count(self):
        """FM_3 has 3 codimension-2 boundary strata (nested collisions)."""
        types = fm3_boundary_types()
        codim2 = [t for t in types if t.codimension == 2]
        assert len(codim2) == 3

    def test_pair_collisions(self):
        """3 pair collisions: (12), (23), (13)."""
        types = fm3_boundary_types()
        pair_names = sorted(t.name for t in types if len(t.name) == 2)
        assert pair_names == ["12", "13", "23"]

    def test_triple_collision(self):
        """1 triple collision: (123)."""
        types = fm3_boundary_types()
        triples = [t for t in types if t.name == "123"]
        assert len(triples) == 1
        assert triples[0].codimension == 1

    def test_nested_collisions(self):
        """3 nested collisions are genuinely planted-forest."""
        types = fm3_boundary_types()
        nested = [t for t in types if t.is_planted_forest]
        assert len(nested) == 3
        nested_names = sorted(t.name for t in nested)
        assert nested_names == ["12<123", "13<123", "23<123"]

    def test_planted_forest_strata_count(self):
        """fm3_planted_forest_strata returns exactly the 3 nested strata."""
        strata = fm3_planted_forest_strata()
        assert len(strata) == 3
        assert all(s.is_planted_forest for s in strata)
        assert all(s.codimension == 2 for s in strata)

    def test_codim1_are_not_planted_forest(self):
        """Codimension-1 strata are NOT planted-forest corrections."""
        types = fm3_boundary_types()
        codim1 = [t for t in types if t.codimension == 1]
        assert all(not t.is_planted_forest for t in codim1)


# ===================================================================
# FM_4 BOUNDARY TYPES
# ===================================================================

class TestFM4BoundaryTypes:
    """Boundary types of FM_4."""

    def test_codim1_count(self):
        """FM_4 has 11 codimension-1 strata: 6 pairs + 4 triples + 1 quadruple."""
        types = fm4_boundary_types()
        codim1 = [t for t in types if t.codimension == 1]
        assert len(codim1) == 11

    def test_codim1_breakdown(self):
        """Verify 6+4+1 = 11 breakdown of codim-1 strata."""
        types = fm4_boundary_types()
        codim1 = [t for t in types if t.codimension == 1]
        pairs = [t for t in codim1 if len(t.subset) == 2]
        triples = [t for t in codim1 if len(t.subset) == 3]
        quads = [t for t in codim1 if len(t.subset) == 4]
        assert len(pairs) == 6
        assert len(triples) == 4
        assert len(quads) == 1

    def test_codim2_has_nested_and_disjoint(self):
        """FM_4 codim-2 strata include both nested and disjoint types."""
        types = fm4_boundary_types()
        codim2 = [t for t in types if t.codimension == 2]
        nested = [t for t in codim2 if "<" in t.name]
        disjoint = [t for t in codim2 if "|" in t.name]
        assert len(nested) > 0
        assert len(disjoint) > 0

    def test_disjoint_pair_count(self):
        """3 disjoint-pair strata in FM_4: {12|34, 13|24, 14|23}."""
        types = fm4_boundary_types()
        codim2 = [t for t in types if t.codimension == 2]
        disjoint = [t for t in codim2 if "|" in t.name]
        assert len(disjoint) == 3

    def test_codim3_exists(self):
        """FM_4 has codimension-3 strata (doubly nested chains)."""
        types = fm4_boundary_types()
        codim3 = [t for t in types if t.codimension == 3]
        assert len(codim3) > 0

    def test_all_codim2_are_planted_forest(self):
        """All codimension >= 2 strata are planted-forest corrections."""
        types = fm4_boundary_types()
        higher_codim = [t for t in types if t.codimension >= 2]
        assert all(t.is_planted_forest for t in higher_codim)


# ===================================================================
# d_pf CORRECTIONS
# ===================================================================

class TestDPFCorrections:
    """The planted-forest corrections d_pf at arities 3 and 4."""

    def test_d_pf_arity3_structure(self):
        """d_pf at arity 3 = 3 * S2 * S3."""
        result = d_pf_arity3(Fraction(1), Fraction(1))
        assert result == Fraction(3)

    def test_d_pf_arity3_zero_cubic(self):
        """d_pf vanishes when cubic shadow is zero (Heisenberg)."""
        result = d_pf_arity3(Fraction(1, 2), Fraction(0))
        assert result == Fraction(0)

    def test_d_pf_arity3_scaling(self):
        """d_pf at arity 3 scales linearly in both S2 and S3."""
        for s2 in [Fraction(1, 2), Fraction(1), Fraction(3, 2)]:
            for s3 in [Fraction(0), Fraction(1), Fraction(2)]:
                result = d_pf_arity3(s2, s3)
                assert result == 3 * s2 * s3

    def test_d_pf_arity4_channels(self):
        """d_pf at arity 4 has 5 distinct channels."""
        channels = d_pf_arity4(Fraction(1), Fraction(1), Fraction(1))
        assert "disjoint_pair" in channels
        assert "binary_ternary" in channels
        assert "ternary_binary" in channels
        assert "doubly_nested" in channels
        assert "quartic" in channels
        assert "total" in channels

    def test_d_pf_arity4_disjoint_pair(self):
        """Disjoint-pair channel has coefficient 3 (three pairings)."""
        channels = d_pf_arity4(Fraction(1), Fraction(0), Fraction(0))
        assert channels["disjoint_pair"] == Fraction(3)

    def test_d_pf_arity4_all_zero(self):
        """When all shadows vanish, d_pf at arity 4 is zero."""
        channels = d_pf_arity4(Fraction(0), Fraction(0), Fraction(0))
        assert channels["total"] == Fraction(0)


# ===================================================================
# MC DICTIONARY
# ===================================================================

class TestMCDictionary:
    """MC dictionary: algebraic <-> geometric."""

    def test_arity2_obstruction_zero(self):
        """At arity 2, the MC obstruction is always zero (kappa closed)."""
        obs = mc_dictionary_algebraic(Fraction(1), Fraction(1), Fraction(1))
        assert obs[2] == Fraction(0)

    def test_arity3_obstruction(self):
        """At arity 3, obstruction = K_2^2."""
        obs = mc_dictionary_algebraic(Fraction(3), Fraction(0), Fraction(0))
        assert obs[3] == Fraction(9)  # 3^2

    def test_arity4_obstruction(self):
        """At arity 4, obstruction = 2*K_2*K_3."""
        obs = mc_dictionary_algebraic(Fraction(2), Fraction(5), Fraction(0))
        assert obs[4] == Fraction(20)  # 2*2*5

    def test_geometric_arity3(self):
        """Geometric side at arity 3: 7 total, 3 planted-forest."""
        geo = mc_dictionary_geometric(3)
        assert geo["total"] == 7
        assert geo["planted_forest_count"] == 3
        assert geo["codim1_count"] == 4
        assert geo["codim2_count"] == 3

    def test_dictionary_verification_holds(self):
        """MC dictionary verification passes at arities 3 and 4."""
        results = mc_dictionary_verify(max_arity=4)
        assert results[3]["dictionary_holds"]
        assert results[4]["dictionary_holds"]

    def test_dictionary_arity3_mechanism(self):
        """Arity 3: 3 nested strata match 3 terms in K_2*K_2."""
        results = mc_dictionary_verify(max_arity=3)
        assert results[3]["geometric_planted_forest_count"] == 3


# ===================================================================
# RESIDUE SPLITTING
# ===================================================================

class TestResidueSplitting:
    """Residue splitting theorem: algebraic x geometric."""

    def test_splitting_product(self):
        """Total residue = algebraic * geometric."""
        tree = make_binary_tree(make_leaf(1), make_leaf(2))
        forest = PlantedForest(trees=(tree,))
        result = residue_splitting(forest, Fraction(3))
        # geometric weight = 1/|Aut| = 1/2
        assert result["geometric_weight"] == Fraction(1, 2)
        assert result["total_residue"] == Fraction(3, 2)

    def test_splitting_with_override(self):
        """Override geometric weight."""
        tree = make_leaf(1)
        forest = PlantedForest(trees=(tree,))
        result = residue_splitting(forest, Fraction(5), geometric_weight=Fraction(1, 3))
        assert result["total_residue"] == Fraction(5, 3)

    def test_splitting_trivial_forest(self):
        """Single leaf forest has |Aut| = 1, so geometric weight = 1."""
        forest = PlantedForest(trees=(make_leaf(1),))
        result = residue_splitting(forest, Fraction(7))
        assert result["geometric_weight"] == Fraction(1)
        assert result["total_residue"] == Fraction(7)


# ===================================================================
# CUBIC GAUGE TRIVIALITY
# ===================================================================

class TestCubicGaugeTriviality:
    """Cubic planted-forest correction o_3^{pf}."""

    def test_heisenberg_trivially_zero(self):
        """Heisenberg: o_3^{pf} = 0 (no cubic shadow at all)."""
        result = cubic_gauge_triviality_check("heisenberg")
        assert result["gauge_trivial"]
        assert result["o3_pf"] == Fraction(0)

    def test_affine_zero_by_jacobi(self):
        """Affine: o_3^{pf} = 0 by Jacobi identity."""
        result = cubic_gauge_triviality_check("affine")
        assert result["gauge_trivial"]
        assert result["o3_pf"] == Fraction(0)
        assert "Jacobi" in result["mechanism"]

    def test_betagamma_zero(self):
        """Beta-gamma: o_3^{pf} = 0 (no cubic, only quartic)."""
        result = cubic_gauge_triviality_check("betagamma")
        assert result["gauge_trivial"]

    def test_virasoro_principal_ds(self):
        """Virasoro: principal DS => o_3^{pf} = 0 at cubic level."""
        result = cubic_gauge_triviality_check("virasoro")
        assert result["gauge_trivial"]
        assert "principal" in result["mechanism"]

    def test_w3_principal_ds(self):
        """W_3 (principal): o_3^{pf} = 0."""
        result = cubic_gauge_triviality_check("w3")
        assert result["gauge_trivial"]

    def test_non_principal_nonzero(self):
        """Non-principal W_3: o_3^{pf} != 0."""
        result = cubic_gauge_triviality_check("w3_non_principal")
        assert not result["gauge_trivial"]
        assert result["o3_pf"] != 0

    def test_shadow_depth_heisenberg(self):
        """Heisenberg shadow depth = 2 (Gaussian)."""
        result = cubic_gauge_triviality_check("heisenberg")
        assert result["shadow_depth"] == 2
        assert result["archetype"] == "Gaussian"

    def test_shadow_depth_affine(self):
        """Affine shadow depth = 3 (Lie/tree)."""
        result = cubic_gauge_triviality_check("affine")
        assert result["shadow_depth"] == 3
        assert result["archetype"] == "Lie/tree"

    def test_shadow_depth_betagamma(self):
        """Beta-gamma shadow depth = 4 (Contact/quartic)."""
        result = cubic_gauge_triviality_check("betagamma")
        assert result["shadow_depth"] == 4
        assert result["archetype"] == "Contact/quartic"

    def test_shadow_depth_virasoro_infinite(self):
        """Virasoro shadow depth = infinite (Mixed modular)."""
        result = cubic_gauge_triviality_check("virasoro")
        assert result["shadow_depth"] == -1  # -1 encodes infinite
        assert result["archetype"] == "Mixed modular"


# ===================================================================
# GENUINELY PLANTED-FOREST CUBIC
# ===================================================================

class TestGenuinelyPlantedForestCubic:
    """The first non-gauge-trivial cubic correction."""

    def test_heisenberg_zero(self):
        """Heisenberg: no genuine cubic correction."""
        result = genuinely_planted_forest_cubic("heisenberg")
        assert not result["nonzero"]
        assert result["o3_pf"] == Fraction(0)

    def test_non_principal_nonzero(self):
        """Non-principal W_3: genuinely nonzero cubic correction."""
        result = genuinely_planted_forest_cubic("w3_non_principal")
        assert result["nonzero"]
        assert result["num_nested_strata"] == 3

    def test_strata_names_correct(self):
        """All families see the same 3 nested strata."""
        for fam in ["heisenberg", "affine", "virasoro", "w3_non_principal"]:
            result = genuinely_planted_forest_cubic(fam)
            assert sorted(result["strata_names"]) == ["12<123", "13<123", "23<123"]


# ===================================================================
# D^2 = 0
# ===================================================================

class TestDSquaredZero:
    """Verification of D^2 = 0 at low arity."""

    def test_arity2(self):
        """D^2 = 0 at arity 2 is trivial."""
        results = d_squared_zero_check(max_arity=2)
        assert results[2]["d_squared_zero"]

    def test_arity3(self):
        """D^2 = 0 at arity 3 via codim-2 face pairing."""
        results = d_squared_zero_check(max_arity=3)
        assert results[3]["d_squared_zero"]
        assert results[3]["num_cancelling_pairs"] == 3

    def test_arity4(self):
        """D^2 = 0 at arity 4 via Mok25."""
        results = d_squared_zero_check(max_arity=4)
        assert results[4]["d_squared_zero"]

    def test_arity3_codim_counts(self):
        """At arity 3: 4 codim-1 strata, 3 codim-2 strata."""
        results = d_squared_zero_check(max_arity=3)
        assert results[3]["codim1_count"] == 4
        assert results[3]["codim2_count"] == 3


# ===================================================================
# CODIMENSION TABLE
# ===================================================================

class TestCodimensionTable:
    """Codimension distribution of FM_n boundary strata."""

    def test_fm2(self):
        """FM_2 has 1 codim-1 stratum."""
        table = codimension_table(2)
        assert table.get(1, 0) == 1

    def test_fm3(self):
        """FM_3: 4 codim-1, 3 codim-2."""
        table = codimension_table(3)
        assert table[1] == 4
        assert table[2] == 3

    def test_fm4_codim1(self):
        """FM_4 has 11 codimension-1 strata."""
        table = codimension_table(4)
        assert table[1] == 11

    def test_fm5_codim1_formula(self):
        """FM_5 has 2^5 - 5 - 1 = 26 codim-1 strata."""
        table = codimension_table(5)
        assert table[1] == 26


# ===================================================================
# MOK TROPICALIZATION
# ===================================================================

class TestMokTropicalization:
    """G_pf = Trop(FM_n(C|D)) verification."""

    def test_n2(self):
        """Tropicalization at n=2: trivially holds."""
        result = mok_tropicalization_verify(2)
        assert result["tropicalization_holds"]

    def test_n3_count(self):
        """Tropicalization at n=3: 7 boundary types = 7 planted-forest cells."""
        result = mok_tropicalization_verify(3)
        assert result["tropicalization_holds"]
        assert result["fm_boundary_types"] == 7
        assert result["planted_forest_cells"] == 7

    def test_n4_holds(self):
        """Tropicalization at n=4 holds."""
        result = mok_tropicalization_verify(4)
        assert result["tropicalization_holds"]

    def test_n3_codim_split(self):
        """At n=3: 4 codim-1, 3 codim-2."""
        result = mok_tropicalization_verify(3)
        assert result["codim1_strata"] == 4
        assert result["codim2_strata"] == 3


# ===================================================================
# BOUNDARY OPERATOR MATRIX
# ===================================================================

class TestBoundaryOperatorMatrix:
    """Incidence structure of FM_n boundary."""

    def test_fm2_trivial(self):
        """FM_2 boundary operator matrix is trivial."""
        result = boundary_operator_matrix(2)
        assert result["d_squared_zero"]

    def test_fm3_incidence(self):
        """FM_3: each codim-2 face has exactly 2 cofaces."""
        result = boundary_operator_matrix(3)
        assert result["d_squared_zero"]
        assert result["num_cancelling_pairs"] == 3

    def test_fm3_incidence_detail(self):
        """FM_3: (12)<(123) is a face of D_{12} and D_{123}."""
        result = boundary_operator_matrix(3)
        inc = result["incidence"]
        assert "12<123" in inc
        coface_names = [name for name, _ in inc["12<123"]]
        assert "12" in coface_names
        assert "123" in coface_names

    def test_fm4_d_squared(self):
        """FM_4: D^2 = 0 holds."""
        result = boundary_operator_matrix(4)
        assert result["d_squared_zero"]


# ===================================================================
# CROSS-CHECKS WITH VOL I modular_bar.py
# ===================================================================

class TestCrossVolumeConsistency:
    """Cross-checks with Vol I modular_bar.py conventions."""

    def test_fm3_count_matches_vol1(self):
        """Our 7 boundary types match Vol I fm3_planted_forest_types count."""
        # Vol I has 7 types too (modular_bar.py)
        types = fm3_boundary_types()
        assert len(types) == 7

    def test_codim1_pair_count_matches(self):
        """3 pair collisions match Vol I convention."""
        types = fm3_boundary_types()
        pairs = [t for t in types if t.codimension == 1 and len(t.subset) == 2]
        assert len(pairs) == 3

    def test_nested_count_matches(self):
        """3 nested collisions match Vol I convention."""
        types = fm3_boundary_types()
        nested = [t for t in types if t.is_planted_forest]
        assert len(nested) == 3

    def test_shadow_archetype_consistency(self):
        """Shadow archetypes match Vol I classification (G/L/C/M)."""
        # Vol I: shadow_archetype(cubic_nonzero, quartic_nonzero)
        # Heisenberg: (False, False) -> Gaussian
        # Affine: (True, False) -> Lie/tree
        # Beta-gamma: (False, True) -> Contact/quartic
        # Virasoro: (True, True) -> Mixed modular
        assert cubic_gauge_triviality_check("heisenberg")["archetype"] == "Gaussian"
        assert cubic_gauge_triviality_check("affine")["archetype"] == "Lie/tree"
        assert cubic_gauge_triviality_check("betagamma")["archetype"] == "Contact/quartic"
        assert cubic_gauge_triviality_check("virasoro")["archetype"] == "Mixed modular"
