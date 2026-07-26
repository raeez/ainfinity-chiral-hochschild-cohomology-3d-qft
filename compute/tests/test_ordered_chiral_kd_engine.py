"""Tests for Ordered Associative Chiral KD Engine.

Verifies the ordered (E_1) sector of chiral Koszul duality:
1. Deconcatenation coproduct: all ordered splits
2. Bar differential: consecutive-collapse formula
3. d^2 = 0: associativity implies d^2 = 0
4. Shuffle product: count, associativity, symmetry
5. Opposite involution: involution property

Each test performs ACTUAL computation, not lookup.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from fractions import Fraction
from math import comb

from lib.ordered_chiral_kd_engine import (
    deconcatenation_coproduct,
    bar_differential_ordered,
    d_squared_check,
    shuffle_product,
    opposite_involution,
    bordered_diagonal_bicomodule_profile,
    annular_cyclic_rotation_sign,
    annular_diagonal_closure_profile,
    commutative_m2,
    free_associative_m2,
)
from lib.fm_boundary import (
    local_residue_convention,
    ordered_chiral_bar_residue_skeleton,
)


# ===================================================================
# 1. DECONCATENATION COPRODUCT
# ===================================================================

class TestDeconcatenation:
    """Test the deconcatenation coproduct on words."""

    def test_length_1_empty(self):
        """A single-letter word has no nontrivial splits."""
        result = deconcatenation_coproduct(['a'])
        assert result == []

    def test_length_0_empty(self):
        """Empty word has no splits."""
        result = deconcatenation_coproduct([])
        assert result == []

    def test_length_2(self):
        """[a, b] splits only as ([a], [b])."""
        result = deconcatenation_coproduct(['a', 'b'])
        assert len(result) == 1
        assert result[0] == (['a'], ['b'])

    def test_length_3(self):
        """[a, b, c] has exactly 2 splits."""
        result = deconcatenation_coproduct(['a', 'b', 'c'])
        assert len(result) == 2
        assert (['a'], ['b', 'c']) in result
        assert (['a', 'b'], ['c']) in result

    def test_length_4(self):
        """[a, b, c, d] has exactly 3 splits."""
        result = deconcatenation_coproduct(['a', 'b', 'c', 'd'])
        assert len(result) == 3
        expected = [
            (['a'], ['b', 'c', 'd']),
            (['a', 'b'], ['c', 'd']),
            (['a', 'b', 'c'], ['d']),
        ]
        for e in expected:
            assert e in result

    def test_split_count_general(self):
        """A word of length n has exactly n-1 splits."""
        for n in range(1, 8):
            word = [chr(ord('a') + i) for i in range(n)]
            result = deconcatenation_coproduct(word)
            assert len(result) == max(0, n - 1)

    def test_splits_reconstruct(self):
        """Each split (left, right) concatenates back to the original word."""
        word = ['x', 'y', 'z', 'w']
        for left, right in deconcatenation_coproduct(word):
            assert left + right == word


# ===================================================================
# 2. BAR DIFFERENTIAL
# ===================================================================

class TestBarDifferential:
    """Test the ordered bar differential d_C."""

    def test_d_of_length_1(self):
        """d of a single element is 0."""
        result = bar_differential_ordered(['a'], commutative_m2)
        assert result == {}

    def test_d_of_pair(self):
        """d[a|b] = m2(a, b), one term, positive sign."""
        result = bar_differential_ordered(['a', 'b'], commutative_m2)
        # m2(a, b) = 'ab' with coefficient 1
        # Sign: (-1)^0 = +1
        assert len(result) == 1
        # commutative_m2 sorts: m2('a','b') = {'ab': 1}
        assert result[('ab',)] == Fraction(1)

    def test_d_of_triple_has_two_terms(self):
        """d[a|b|c] has exactly 2 terms (consecutive collapse at positions 0 and 1)."""
        result = bar_differential_ordered(['a', 'b', 'c'], commutative_m2)
        # i=0: +1 * [m2(a,b) | c] = +[ab | c]
        # i=1: -1 * [a | m2(b,c)] = -[a | bc]
        assert len(result) == 2
        assert result[('ab', 'c')] == Fraction(1)
        assert result[('a', 'bc')] == Fraction(-1)

    def test_d_of_triple_noncommutative(self):
        """d[a|b|c] with free associative m2."""
        result = bar_differential_ordered(['a', 'b', 'c'], free_associative_m2)
        # i=0: +1 * [ab | c]
        # i=1: -1 * [a | bc]
        assert result[('ab', 'c')] == Fraction(1)
        assert result[('a', 'bc')] == Fraction(-1)

    def test_d_of_quadruple_has_three_terms(self):
        """d[a|b|c|d] has 3 terms."""
        result = bar_differential_ordered(['a', 'b', 'c', 'd'], free_associative_m2)
        # i=0: +1 * [ab | c | d]
        # i=1: -1 * [a | bc | d]
        # i=2: +1 * [a | b | cd]
        assert len(result) == 3
        assert result[('ab', 'c', 'd')] == Fraction(1)
        assert result[('a', 'bc', 'd')] == Fraction(-1)
        assert result[('a', 'b', 'cd')] == Fraction(1)

    def test_sl2_low_degree_source_has_fixed_desuspension_signs(self):
        """The low-degree sl2 ordered bar computation must not use free ± signs."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        for rel_path in (
            'chapters/connections/ordered_associative_chiral_kd_core.tex',
            'chapters/connections/ordered_associative_chiral_kd.tex',
        ):
            source = os.path.join(repo_root, rel_path)
            with open(source, encoding='utf-8') as handle:
                text = handle.read()

            start = text.index(r'For $\widehat{\mathfrak{sl}}_2$ at level~$k$')
            end = text.index(r'Bar cohomology at degree~$1$:', start)
            block = text[start:end]
            flat = ' '.join(block.split())

            required = (
                r'(-1)^{i-1}',
                r'[e,f]=h,\ [h,e]=2e,\ [h,f]=-2f',
                r's^{-1}(e_{(0)}f)\;=\;s^{-1}h',
                r's^{-1}(e_{(0)}h)\;=\;-2\,s^{-1}e',
                r's^{-1}(h_{(0)}e)\;=\;2\,s^{-1}e',
                r's^{-1}e\otimes s^{-1}h',
                r'-\,s^{-1}h\otimes s^{-1}e',
                'because the two desuspended weight-one currents are odd',
            )
            for needle in required:
                assert ' '.join(needle.split()) in flat

            retired = (
                r'\pm s^{-1}(e_{(0)}f)',
                r'\pm s^{-1}(e_{(0)}h)',
                r'\pm s^{-1}(h_{(0)}e)',
                r'h\otimes e$ up to sign',
            )
            for needle in retired:
                assert needle not in block

    def test_ordered_low_degree_source_has_fixed_face_signs(self):
        """The low-degree construction uses face signs for mu_2, not the Lie bracket."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        for rel_path in (
            'chapters/connections/ordered_associative_chiral_kd_core.tex',
            'chapters/connections/ordered_associative_chiral_kd.tex',
        ):
            source = os.path.join(repo_root, rel_path)
            with open(source, encoding='utf-8') as handle:
                text = handle.read()

            start = text.index(r'\mu_2(e_I,e_J)')
            end = text.index('This vanishes if and only if', start)
            block = text[start:end]
            flat = ' '.join(block.split())

            required = (
                r'\mu_2(e_I,e_J)=\sum_K\mu^K_{IJ}e_K',
                'This symbol is not the vertex-algebra',
                r'then \(n=0\) is the Lie bracket',
                r'not to the associative \(E_1\) proof below',
                r'd_{\mathrm{res}}[s^{-1}a_1|\cdots|s^{-1}a_m]',
                r'\sum_{i=1}^{m-1}(-1)^{i-1}',
                r's^{-1}\mu_2(a_i,a_{i+1})',
                r'\sum_K\mu^K_{IJ}\;s^{-1}e_K',
                r'\sum_L \mu^L_{IJ}\;',
                r'\;-\;',
                r'\sum_M \mu^M_{JK}\;',
                r'\sum_{L,M,N}\bigl(',
                r'\mu^N_{LK}\,\mu^L_{IJ}',
                r'-\mu^N_{IM}\,\mu^M_{JK}',
                r's^{-1}e_N',
                r'not an \(E_1\)-associative differential by itself',
            )
            for needle in required:
                assert ' '.join(needle.split()) in flat

            retired = (
                r'\pm\,c^K_{IJ;n}',
                r'\pm\,c^L_{IJ;0}',
                r'\pm\,c^M_{JK;0}',
                r'(-1)^{|e_I|}',
                r's^{-1}e_?',
                r'\sum_{i=1}^{m-1}\sum_{n\ge0}(-1)^{i-1}',
                r'\sum_{K,n}c^K_{IJ;n}',
                r'c^N_{LK;0}\,c^L_{IJ;0}',
                r'c^N_{IM;0}\,c^M_{JK;0}',
                'For the simple-pole residue',
                r'associativity of the $0$-th product',
                'In the ordered complex, no Arnold relation is needed',
            )
            for needle in retired:
                assert needle not in block

            heisenberg_retired = (
                r'\pm k\,s^{-1}(1)',
                r'\pm s^{-1}(\alpha_{(1)}\alpha)',
                r'\pm s^{-1}(k)',
            )
            for needle in heisenberg_retired:
                assert needle not in text

            assert r'=k\,s^{-1}(1)' in text

            if rel_path.endswith('_core.tex'):
                assert r'associativity of the open-colour product \(\mu_2\)' in text
                assert r'The Lie bracket is not the open \(E_1\) product' in text

    def test_ybe_boundary_face_proof_uses_mu2_associativity(self):
        """The YBE-from-d^2 proof must not use the vertex 0-mode as the E_1 product."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        for rel_path in (
            'chapters/connections/ordered_associative_chiral_kd_core.tex',
            'chapters/connections/ordered_associative_chiral_kd.tex',
        ):
            source = os.path.join(repo_root, rel_path)
            with open(source, encoding='utf-8') as handle:
                text = handle.read()

            start = text.index('There are exactly two pathways')
            end = text.index('equality is $d^2=0$.', start)
            block = text[start:end]
            flat = ' '.join(block.split())

            required = (
                r'applying \(\mu_2(a,b)\)',
                r'applying \(\mu_2(\mu_2(a,b),c)\)',
                r'applying \(\mu_2(b,c)\)',
                r'applying \(\mu_2(a,\mu_2(b,c))\)',
                r'associativity of',
                r'open-colour product \(\mu_2\)',
                'monodromy part',
            )
            for needle in required:
                assert ' '.join(needle.split()) in flat

            retired = (
                r'extracting $a_{(0)}b$',
                r'extracting $(a_{(0)}b)_{(0)}c$',
                r'extracting $b_{(0)}c$',
                r'extracting $a_{(0)}(b_{(0)}c)$',
                r'the zeroth product ($E_1$ axiom)',
            )
            for needle in retired:
                assert needle not in block

    def test_concentration_theorem_separates_mu2_from_vertex_depth(self):
        """The concentration theorem must not make Lie homology the open E_1 page."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        core_path = os.path.join(
            repo_root,
            'chapters/connections/ordered_associative_chiral_kd_core.tex',
        )
        with open(core_path, encoding='utf-8') as handle:
            text = handle.read()

        start = text.index(r'\begin{theorem}[Bar cohomology concentration;')
        end = text.index(r'\begin{computation}[Concentration data', start)
        block = text[start:end]
        flat = ' '.join(block.split())

        required = (
            r'first differential \(d_\mu\)',
            r'associative product \(\mu_2\)',
            'does not have Lie algebra homology as its first page',
            r'a local vertex-algebra shadow',
            r'[a,b]=a_{(0)}b',
            r'converges to the symmetric chiral\slash Francis--Gaitsgory shadow',
            r'not to the open-colour \(E_1\) proof',
            r'extracts the open product \(\mu_2(a,b)\)',
            r'V=Q_{\mu_2}(\bar\cA)',
            r'\Barchord(\bar\cA,\mu_2)',
            r'Remark~\ref{rem:depth-decomposition-d-squared}',
        )
        for needle in required:
            assert ' '.join(needle.split()) in flat

        retired = (
            r'zeroth products) converges to $H^*(\Barchord)$',
            r'extracts the zeroth product',
            r'a_{(0)}b$; its image is the space of decomposable',
            r'The $E_0$ page uses only the zeroth product',
        )
        for needle in retired:
            assert needle not in block

        table_start = text.index(r'\begin{computation}[Concentration data')
        table_end = text.index(r'\begin{computation}[Dimension table', table_start)
        table_block = text[table_start:table_end]
        table_flat = ' '.join(table_block.split())
        assert 'vertex-shadow depth page' in table_flat
        assert 'The last column records the separate vertex-shadow spectral sequence' in table_flat
        assert r'not the open-colour \(E_1\) concentration proof' in table_flat

        for rel_path in (
            'chapters/connections/ordered_associative_chiral_kd_core.tex',
            'chapters/connections/ordered_associative_chiral_kd.tex',
        ):
            source = os.path.join(repo_root, rel_path)
            with open(source, encoding='utf-8') as handle:
                source_text = handle.read()
            r_start = source_text.index(r'The \emph{$R$-matrix} is the monodromy')
            r_end = source_text.index(r'The $R$-matrix intertwines', r_start)
            r_block = source_text[r_start:r_end]
            r_flat = ' '.join(r_block.split())

            assert 'vertex-shadow OPE matrix' in r_flat
            assert 'classical $r$-matrix} of the local vertex-algebra shadow' in r_flat
            assert 'It is not the open-colour associative product' in r_flat
            assert r'\mu_2' in r_block
            assert 'of zeroth products (Lie brackets)' not in r_block

    def test_yangian_double_bar_is_associated_graded_shadow(self):
        """The Yangian double-bar claim is PBW/current shadow, not full equality."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        core_path = os.path.join(
            repo_root,
            'chapters/connections/ordered_associative_chiral_kd_core.tex',
        )
        with open(core_path, encoding='utf-8') as handle:
            text = handle.read()

        start = text.index(r'\begin{proposition}[Associated-graded current shadow')
        end = text.index(r'\end{proposition}', start)
        prop = text[start:end]
        prop_flat = ' '.join(prop.split())

        required = (
            'Drinfeld',
            'PBW filtration',
            r'\operatorname{gr}_F Y_\hbar(\fg)',
            r'U(\fg[t])',
            'primitive collision-residue projection',
            r'[-,-]_{\fg[t]}',
            r'not an isomorphism',
            r'\Barchord(Y_\hbar(\fg))\simeq U(\fg[t])',
            'full ordered bar objects',
            'central extension',
            'associated-graded current shadow',
        )
        for needle in required:
            assert ' '.join(needle.split()) in prop_flat

        retired_in_prop = (
            r'Then the open-colour double bar',
            r'recovers the \emph{current algebra} $U(\fg[t])$',
            r'\Barchord\bigl(Y_\hbar(\fg)\bigr) \;\xrightarrow{\;\sim\;} U(\fg[t])',
        )
        for needle in retired_in_prop:
            assert ' '.join(needle.split()) not in prop_flat

        shadow_start = text.index(r'The open-colour Koszul dual of the Yangian itself')
        shadow_end = text.index(r'\noindent\emph{(5)', shadow_start)
        shadow = text[shadow_start:shadow_end]
        shadow_flat = ' '.join(shadow.split())

        assert 'not computed by discarding the deformation parameter' in shadow_flat
        assert 'associated-graded current-shadow calculation' in shadow_flat
        assert 'not an isomorphism of full' in shadow_flat
        assert 'This chain records a shadow, not a second full Koszul duality' in shadow_flat
        assert 'the full ordered bar' in shadow_flat
        assert 'The current shadow is not the Yangian' in shadow_flat

        summary_start = text.index(r'\textbf{Summary at degree $2$:}')
        summary_end = text.index(r'\medskip', summary_start)
        summary = text[summary_start:summary_end]
        summary_flat = ' '.join(summary.split())
        assert r'\operatorname{gr}_F d_{\mathrm{bar}}^{\mathrm{prim,coll}}' in summary
        assert 'associated-graded current-shadow clause' in summary_flat
        assert r'\Barchord(Y_\hbar(\mathfrak{sl}_2)) \simeq U(\mathfrak{sl}_2[t])' not in summary

        sl3_guard = os.path.join(repo_root, 'compute/tests/test_double_bar_sl3.py')
        with open(sl3_guard, encoding='utf-8') as handle:
            sl3_text = handle.read()
        assert 'Associated-graded collision shadow' in sl3_text
        assert 'This is not an isomorphism' in sl3_text
        assert 'B^{ord,ch}(Y_hbar(sl_N)) = U(sl_N[t])' not in sl3_text

        quantum_note = os.path.join(repo_root, 'compute/quantum_lattice_voa_bar.py')
        with open(quantum_note, encoding='utf-8') as handle:
            quantum_text = handle.read()
        assert 'deformation-parameter loss is an associated-graded' in quantum_text
        assert 'full_yangian_hbar_visible' in quantum_text
        assert 'yangian_primitive_shadow_loses_hbar' in quantum_text
        assert 'quantum_lattice_collision_shadow_hbar_visible' in quantum_text
        assert 'The full ordered Yangian bar and the RTT product retain hbar' in quantum_text
        assert 'not an isomorphism' in quantum_text
        assert 'Previous session established: B^ord(Y_hbar(sl_2)) = U(sl_2[t])' not in quantum_text
        assert 'yangian_hbar_invisible' not in quantum_text
        assert 'quantum_lattice_hbar_invisible' not in quantum_text
        assert 'B^ord sees only Omega' not in quantum_text

    def test_quantum_lattice_oracle_is_finite_window_not_firstness_claim(self):
        """The quantum-lattice oracle is a finite A1 benchmark, not a global theorem."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        quantum_note = os.path.join(repo_root, 'compute/quantum_lattice_voa_bar.py')
        with open(quantum_note, encoding='utf-8') as handle:
            quantum_text = handle.read()

        required = (
            'Finite-window A_1 check',
            'A1 finite-window oracle',
            'not a theorem that every higher operation or every quantum lattice',
            'Run the focused quantum-lattice VOA bar oracle',
            'explicit',
            'A_1 benchmark here for genuinely E_1-chiral ordered-bar behaviour',
            'FINITE-WINDOW q-STABILITY',
            'the tested m4 window',
            'extra quasi-triangular input beyond',
        )
        for needle in required:
            assert needle in quantum_text

        retired = (
            'The quantum lattice VOA V_q (Etingof-Kazhdan type) is the FIRST',
            'first genuinely E_1-chiral algebra',
            'the complete q-stability dichotomy',
            'Run the complete test suite for quantum lattice VOA bar complex',
            'This is the deepest structural fact',
            'topological invariant',
            'm_4 = 0 -> m_4 = 0 (tower truncation)',
            'Poincare series, tower truncation) are q-STABLE',
            'The trigonometric R-matrix is NOT derivable from a local OPE',
            'it is independent input (tier (iii) of the three-tier picture)',
        )
        for needle in retired:
            assert needle not in quantum_text

    def test_quantum_affine_to_yangian_requires_additive_scaling(self):
        """Bare q->1 gives the affine/loop limit, not the Yangian."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        quantum_note = os.path.join(repo_root, 'compute', 'quantum_lattice_voa_bar.py')
        with open(quantum_note, encoding='utf-8') as handle:
            quantum_text = handle.read()

        required = (
            'not by the bare specialization q -> 1',
            'bare limit',
            'q -> 1 at fixed z gives the undeformed multiplicative/affine',
            'additive scaling z = exp(eps*u), q = exp(eps*hbar/2), eps -> 0',
            'bare_q_to_1_limit',
            'rational_scaling_limit',
            'Bare q->1: undeformed affine/loop object',
            'Yangian: additive rational scaling degeneration',
            'Bare q->1 gives the undeformed',
            'not the Yangian',
        )
        for needle in required:
            assert needle in quantum_text

        retired = (
            'U_q(sl_2-hat) -> Y_hbar(sl_2) as q -> 1',
            'Classical limit: U_q -> Y_hbar as q -> 1',
            'At q -> 1: R(z) -> P (the permutation), recovering the rational limit',
            'leading-order approximation to the full trigonometric R-matrix',
        )
        for needle in retired:
            assert needle not in quantum_text

    def test_yangian_oracles_use_mu2_not_vertex_zero_mode_product(self):
        """Executable Yangian notes must not call the vertex 0-mode the E_1 product."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

        sl2_path = os.path.join(repo_root, 'compute', 'ordered_e1_shadow_sl2.py')
        with open(sl2_path, encoding='utf-8') as handle:
            sl2_text = handle.read()
        sl2_flat = ' '.join(sl2_text.split())

        required_sl2 = (
            'For E_1-chiral algebras: d^2 = 0 <=> associativity of mu_2.',
            'The product mu_2 is not the vertex 0-mode a_{(0)}b.',
            'Its RTT associative multiplication supplies the open-colour product mu_2.',
            'The primitive collision-residue shadow is the current Lie bracket, not the product.',
            'd^2 = 0 on B^ord_3(Y_hbar) holds by mu_2-associativity.',
            'Yangian RTT multiplication supplies mu_2;',
            'the V_k(sl_2) vertex 0-mode Lie bracket does not.',
        )
        for needle in required_sl2:
            assert needle in sl2_flat

        retired_sl2 = (
            'For E_1-chiral algebras: d^2 = 0 <=> associativity of zeroth product.',
            'The Yangian Y_hbar(sl_2) IS E_1-chiral. Its zeroth product IS associative.',
            'The zeroth product IS associative (it is the Yangian product).',
            r'$a_{(0)}b$ is the Yangian multiplication',
            'd^2 = 0 requires associative zeroth product',
        )
        for needle in retired_sl2:
            assert needle not in sl2_text

        sl3_path = os.path.join(repo_root, 'compute', 'tests', 'test_double_bar_sl3.py')
        with open(sl3_path, encoding='utf-8') as handle:
            sl3_text = handle.read()
        sl3_flat = ' '.join(sl3_text.split())

        required_sl3 = (
            'bar differential extracts the primitive collision-residue bracket',
            'For the associated-graded primitive shadow computation',
            'The primitive collision residue is the sl_3 Lie bracket.',
            'The primitive collision-residue bracket for Y_hbar(sl_3).',
            'At the level of the collision-residue shadow, the bracket is not the open E_1 product',
            'The associated-graded primitive shadow extends uniformly to sl_N.',
        )
        for needle in required_sl3:
            assert needle in sl3_flat

        retired_sl3 = (
            'bar differential extracts zeroth product (Lie bracket)',
            'The zeroth product (residue at simple pole)',
            'The COLLISION RESIDUE (zeroth product)',
            'The zeroth product (collision residue) for Y_hbar(sl_3).',
            'zeroth product = Lie bracket',
            'The degree-3 consistency of the double bar',
            'Lie bracket recovered by the double bar',
        )
        for needle in retired_sl3:
            assert needle not in sl3_text


# ===================================================================
# 3. D-SQUARED = 0
# ===================================================================

class TestDSquared:
    """Test d^2 = 0 (equivalent to associativity of m2)."""

    def test_d_squared_arity_2_commutative(self):
        """d^2 = 0 on [a|b] with commutative m2."""
        result = d_squared_check(['a', 'b'], commutative_m2)
        assert result == {}

    def test_d_squared_arity_3_commutative(self):
        """d^2 = 0 on [a|b|c] with commutative m2."""
        result = d_squared_check(['a', 'b', 'c'], commutative_m2)
        assert result == {}

    def test_d_squared_arity_4_commutative(self):
        """d^2 = 0 on [a|b|c|d] with commutative m2."""
        result = d_squared_check(['a', 'b', 'c', 'd'], commutative_m2)
        assert result == {}

    def test_d_squared_arity_3_free_assoc(self):
        """d^2 = 0 on [a|b|c] with free associative m2."""
        result = d_squared_check(['a', 'b', 'c'], free_associative_m2)
        assert result == {}

    def test_d_squared_arity_4_free_assoc(self):
        """d^2 = 0 on [a|b|c|d] with free associative m2."""
        result = d_squared_check(['a', 'b', 'c', 'd'], free_associative_m2)
        assert result == {}

    def test_d_squared_arity_5_commutative(self):
        """d^2 = 0 at arity 5."""
        word = ['a', 'b', 'c', 'd', 'e']
        result = d_squared_check(word, commutative_m2)
        assert result == {}

    def test_d_squared_repeated_generators(self):
        """d^2 = 0 even with repeated generators."""
        result = d_squared_check(['a', 'a', 'a'], commutative_m2)
        assert result == {}

    def test_d_squared_fails_nonassociative(self):
        """d^2 != 0 for a non-associative product.

        Define a product where m2(m2(a,b),c) != m2(a,m2(b,c)) by
        making left-association and right-association produce distinct labels.
        """
        def nonassoc_m2(a: str, b: str):
            # m2 that is NOT associative:
            # Always wraps with explicit parentheses, so
            # m2(m2(a,b), c) = ((ab)c) but m2(a, m2(b,c)) = (a(bc))
            return {'(' + a + b + ')': Fraction(1)}

        # For [a|b|c]:
        # d[a|b|c] = [(ab)|c] - [a|(bc)]
        # d(d[a|b|c]) = m2((ab),c) - m2(a,(bc))
        #             = ((ab)c) - (a(bc))  -- different strings!
        result = d_squared_check(['a', 'b', 'c'], nonassoc_m2)
        # Should NOT be zero for this non-associative product
        assert result != {}


# ===================================================================
# 4. SHUFFLE PRODUCT
# ===================================================================

class TestShuffle:
    """Test the shuffle product on words."""

    def test_shuffle_count_1_1(self):
        """Sh([a], [b]) has C(2,1) = 2 shuffles."""
        result = shuffle_product(['a'], ['b'])
        assert len(result) == comb(2, 1)
        assert ('a', 'b') in result
        assert ('b', 'a') in result

    def test_shuffle_count_2_1(self):
        """Sh([a,b], [c]) has C(3,2) = 3 shuffles."""
        result = shuffle_product(['a', 'b'], ['c'])
        assert len(result) == comb(3, 2)

    def test_shuffle_count_2_2(self):
        """Sh([a,b], [c,d]) has C(4,2) = 6 shuffles."""
        result = shuffle_product(['a', 'b'], ['c', 'd'])
        assert len(result) == comb(4, 2)

    def test_shuffle_count_3_2(self):
        """Sh of length 3 and 2 has C(5,3) = 10 shuffles."""
        result = shuffle_product(['a', 'b', 'c'], ['d', 'e'])
        assert len(result) == comb(5, 3)

    def test_shuffle_count_general(self):
        """Sh(w1, w2) has exactly C(|w1|+|w2|, |w1|) elements."""
        for p in range(1, 5):
            for q in range(1, 5):
                w1 = [f'a{i}' for i in range(p)]
                w2 = [f'b{j}' for j in range(q)]
                result = shuffle_product(w1, w2)
                assert len(result) == comb(p + q, p)

    def test_shuffle_preserves_order_w1(self):
        """Each shuffle preserves the internal order of w1."""
        w1 = ['a', 'b', 'c']
        w2 = ['x', 'y']
        for shuf in shuffle_product(w1, w2):
            # Extract positions of w1 elements
            positions = [i for i, s in enumerate(shuf) if s in ['a', 'b', 'c']]
            extracted = [shuf[i] for i in positions]
            assert extracted == ['a', 'b', 'c']

    def test_shuffle_preserves_order_w2(self):
        """Each shuffle preserves the internal order of w2."""
        w1 = ['a', 'b']
        w2 = ['x', 'y', 'z']
        for shuf in shuffle_product(w1, w2):
            positions = [i for i, s in enumerate(shuf) if s in ['x', 'y', 'z']]
            extracted = [shuf[i] for i in positions]
            assert extracted == ['x', 'y', 'z']

    def test_shuffle_empty_w1(self):
        """Shuffle with empty w1 returns just w2."""
        result = shuffle_product([], ['a', 'b'])
        assert len(result) == 1
        assert result[0] == ('a', 'b')

    def test_shuffle_empty_w2(self):
        """Shuffle with empty w2 returns just w1."""
        result = shuffle_product(['a', 'b'], [])
        assert len(result) == 1
        assert result[0] == ('a', 'b')


# ===================================================================
# 5. OPPOSITE INVOLUTION
# ===================================================================

class TestOpposite:
    """Test the opposite (reversal) involution."""

    def test_involution_squared_is_identity(self):
        """(w^op)^op = w for several words."""
        words = [
            ['a'],
            ['a', 'b'],
            ['a', 'b', 'c'],
            ['a', 'b', 'c', 'd'],
            ['x', 'y', 'z', 'w', 'v'],
        ]
        for w in words:
            assert opposite_involution(opposite_involution(w)) == w

    def test_reversal_of_pair(self):
        """[a, b]^op = [b, a]."""
        assert opposite_involution(['a', 'b']) == ['b', 'a']

    def test_reversal_of_triple(self):
        """[a, b, c]^op = [c, b, a]."""
        assert opposite_involution(['a', 'b', 'c']) == ['c', 'b', 'a']

    def test_palindrome_fixed(self):
        """A palindrome is a fixed point of the involution."""
        w = ['a', 'b', 'a']
        assert opposite_involution(w) == w

    def test_empty_word(self):
        """Empty word is a fixed point."""
        assert opposite_involution([]) == []


# ===================================================================
# 6. ORDERED CHIRAL BAR RESIDUE SKELETON
# ===================================================================

class TestOrderedChiralBarResidueSkeleton:
    """Guard the local collision-residue skeleton in foundations.tex."""

    def test_local_residue_and_orientation_convention(self):
        """The FM convention extracts the dlog coefficient with outward normal."""
        convention = local_residue_convention()

        assert (
            convention['residue_formula']
            == 'Res_D_S(dlog(epsilon_S) wedge beta + alpha) = beta|epsilon_S=0'
        )
        assert (
            convention['orientation_formula']
            == 'or(FM_k(C)) = (-d epsilon_S) wedge or(D_S)'
        )
        assert convention['outward_normal'] == '-partial_epsilon_S'

    def test_coderivation_projection_formula(self):
        """The cogenerator projection is pi_1 D_A = s^-1 m_k s^otimes k."""
        profile = ordered_chiral_bar_residue_skeleton(4)

        assert profile['coalgebra'] == 'B^ch(A) = T^c(s^-1 bar A)'
        assert profile['coderivation'] == 'D_A = sum_{k>=1} D_{m_k}'
        assert (
            profile['projection_formula']
            == 'pi_1 D_A | (s^-1 bar A)^otimes 4 = s^-1 m_4 s^otimes 4'
        )

    def test_collision_residue_and_d_squared_sources(self):
        """D_A^2 is guarded by FM boundary-of-boundary and Arnold relation."""
        profile = ordered_chiral_bar_residue_skeleton(3)

        assert 'Res_D_S prod_{e in E(S)}' in profile['collision_residue']
        assert 'dlog(z_{s(e)}-z_{t(e)}) tensor OPE_S' in profile['collision_residue']
        assert 'partial^2[FM_k(C)] = 0' in profile['d_squared_sources']
        assert any(
            'omega_ij wedge omega_jk' in source
            for source in profile['d_squared_sources']
        )

    def test_derived_centre_brace_and_mixed_operations(self):
        """The centre is Coder(B^ch(A)) with ordered-insertion braces."""
        profile = ordered_chiral_bar_residue_skeleton(2)

        assert profile['derived_centre'] == 'Z_der_ch(A) = Coder(B^ch(A))'
        assert profile['centre_differential'] == 'd = [D_A,-]'
        assert 'sum_ordered_insertions' in profile['brace_formula']
        assert profile['mixed_operations'].startswith('nu_k^(m): Z_der_ch(A)^otimes k')

    def test_strictification_topologisation_and_class_m_completion(self):
        """Strict W and E_3 clauses are conditional on named hypotheses."""
        profile = ordered_chiral_bar_residue_skeleton(5)

        assert profile['strictification_requires'] == (
            'Drinfeld associator',
            'two-coloured contracting homotopy',
        )
        assert 'd h_oc + h_oc d = id - i circ p' in profile['contracting_homotopy']
        assert 'p h_oc = h_oc i = h_oc^2 = 0' in profile['contracting_homotopy']
        assert 'T=[Q,G]' in profile['topologisation_condition']
        assert 'translations are Q-exact' in profile['topologisation_condition']
        assert 'sum_k C_k(rho) < infinity' in profile['class_m_completion']

    def test_nonpositive_arity_rejected(self):
        """The Taylor arity must be positive."""
        with pytest.raises(ValueError):
            ordered_chiral_bar_residue_skeleton(0)

    def test_foundations_source_contains_ordered_bar_skeleton(self):
        """Source guard: the active foundations chapter carries the theorem."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        source = os.path.join(repo_root, 'chapters', 'theory', 'foundations.tex')
        with open(source, encoding='utf-8') as handle:
            text = handle.read()

        assert 'thm:ordered-chiral-bar-residue-skeleton' in text
        assert r'\pi_1D_A\big|_{(s^{-1}\bar A)^{\otimes k}}' in text
        assert r'\operatorname{Res}_{D_S}' in text
        assert r'\operatorname{or}\bigl(\FM_k(C)\bigr)' in text
        assert r'\Zderch{A}' in text
        assert r'h_{\mathrm{oc}}' in text


# ===================================================================
# 7. BORDERED DIAGONAL BICOMODULE AND ANNULAR CLOSURE
# ===================================================================

class TestBorderedDiagonalBicomodule:
    """Guard the bordered FM chain model for the diagonal bicomodule."""

    def test_three_bordered_stratum_types_present(self):
        """Conf_{p,q}(H) has interior, boundary, and mixed divisors."""
        profile = bordered_diagonal_bicomodule_profile(2, 2)
        types = {entry['type'] for entry in profile['stratum_types']}

        assert types == {
            'interior-interior',
            'boundary-boundary',
            'interior-boundary',
        }
        equations = {entry['local_equation'] for entry in profile['stratum_types']}
        assert 'z_i -> z_j' in equations
        assert 'x_a -> x_b' in equations
        assert 'Im z_i -> 0 and Re z_i -> x_a' in equations

    def test_diagonal_bimodule_approach_convention(self):
        """The two A-actions come from the two boundary approach directions."""
        profile = bordered_diagonal_bicomodule_profile(1, 1)

        assert profile['diagonal_bimodule'] == 'Delta_A = _A A_A'
        assert profile['left_action'] == 'negative boundary approach'
        assert profile['right_action'] == 'positive boundary approach'
        assert profile['forest_roots'] == ('rho_bulk', 'rho_boundary')

    def test_residue_differential_decomposition(self):
        """The diagonal model splits into interior, boundary, and mixed residues."""
        profile = bordered_diagonal_bicomodule_profile(3, 2)
        differential = profile['differential']

        assert differential['d_int'] == 'sum_S Res_{D_S^int} tensor m_S'
        assert differential['d_boundary'] == 'sum_I Res_{D_I^boundary} tensor mu_I'
        assert differential['d_mix'] == 'sum_{S,I} Res_{D_{S,I}^mix} tensor nu_{S,I}'
        assert 'C_log(Confbar_{p,q}(H))' in profile['chain_model']

    def test_missing_stratum_types_are_not_fabricated(self):
        """The profile only lists geometrically available codimension-one types."""
        profile = bordered_diagonal_bicomodule_profile(0, 3)
        types = {entry['type'] for entry in profile['stratum_types']}

        assert types == {'boundary-boundary'}

    def test_negative_point_counts_rejected(self):
        """Point counts are arities and must be nonnegative."""
        with pytest.raises(ValueError):
            bordered_diagonal_bicomodule_profile(-1, 0)

    def test_annular_cyclic_sign_even_and_odd_cases(self):
        """The cyclic sign uses desuspended degrees |a|-1."""
        assert annular_cyclic_rotation_sign([0, 0, 0]) == 1
        assert annular_cyclic_rotation_sign([1, 0, 0]) == -1
        assert annular_cyclic_rotation_sign([2, 3, 3]) == 1

    def test_empty_cyclic_sign_rejected(self):
        """There is no generator to rotate in the empty word."""
        with pytest.raises(ValueError):
            annular_cyclic_rotation_sign([])

    def test_annular_closure_profile(self):
        """The annular closure takes cyclic coinvariants and carries genus-one data."""
        profile = annular_diagonal_closure_profile(4)

        assert 'Conf_n(S^1 x R_{>=0})' in profile['annular_closure']
        assert profile['coinvariants'] == 'Z/n cyclic coinvariants'
        assert profile['cyclic_sign_formula'] == (
            '(-1)^((|a_n|-1) * sum_{i<n} (|a_i|-1))'
        )
        assert profile['first_genus_one_data'] == (
            'R-matrix monodromy',
            'trace pairing',
        )

    def test_active_source_contains_bordered_diagonal_model(self):
        """Source guard: active ordered KD core contains the bordered model."""
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
        source = os.path.join(
            repo_root,
            'chapters',
            'connections',
            'ordered_associative_chiral_kd_core.tex',
        )
        with open(source, encoding='utf-8') as handle:
            text = handle.read()

        assert 'thm:bordered-diagonal-bicomodule-chain-model' in text
        assert 'eq:bordered-diagonal-chain' in text
        assert 'eq:bordered-diagonal-differential' in text
        assert 'eq:bordered-annular-closure' in text
        assert 'eq:annular-cyclic-sign' in text
        assert r'\operatorname{Res}_{D_{S,I}^{\mathrm{mix}}}' in text
