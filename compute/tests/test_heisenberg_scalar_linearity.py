"""Heisenberg connected scalar genus amplitudes are linear in the level."""

from __future__ import annotations

import unittest
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ROSETTA = ROOT / "chapters" / "examples" / "rosetta_stone.tex"
PREFACE = ROOT / "chapters" / "frame" / "preface.tex"
MODULAR_BOOTSTRAP = ROOT / "chapters" / "connections" / "thqg_modular_bootstrap.tex"


def between(source: str, start: str, end: str) -> str:
    begin = source.index(start)
    finish = source.index(end, begin)
    return source[begin:finish]


@lru_cache(maxsize=None)
def bernoulli(n: int) -> Fraction:
    if n == 0:
        return Fraction(1)
    if n == 1:
        return Fraction(-1, 2)
    if n > 1 and n % 2 == 1:
        return Fraction(0)

    total = Fraction(0)
    for k in range(n):
        total += Fraction(comb(n + 1, k)) * bernoulli(k)
    return -total / Fraction(n + 1)


def lambda_fp(g: int) -> Fraction:
    if g < 1:
        raise ValueError("lambda_g^FP requires g >= 1")
    return (
        Fraction(2 ** (2 * g - 1) - 1, 2 ** (2 * g - 1))
        * abs(bernoulli(2 * g))
        / Fraction(factorial(2 * g))
    )


def heisenberg_connected_free_energy(level: Fraction, rank: int, g: int) -> Fraction:
    return rank * level * lambda_fp(g)


class TestHeisenbergScalarLinearity(unittest.TestCase):
    def test_genus_two_is_linear_in_level(self):
        self.assertEqual(lambda_fp(2), Fraction(7, 5760))
        self.assertEqual(
            heisenberg_connected_free_energy(Fraction(3), 1, 2),
            Fraction(21, 5760),
        )

    def test_rank_d_is_additive(self):
        self.assertEqual(
            heisenberg_connected_free_energy(Fraction(5), 4, 3),
            Fraction(20) * Fraction(31, 967680),
        )

    def test_power_law_is_disconnected_partition_bookkeeping(self):
        level = Fraction(3)
        connected = heisenberg_connected_free_energy(level, 1, 2)
        wrong_power = level**2 * lambda_fp(2)
        self.assertNotEqual(connected, wrong_power)

    def test_rosetta_opening_scopes_scalar_line_and_bulk_claims(self):
        source = ROSETTA.read_text()
        preface = PREFACE.read_text()
        flat = " ".join(source.split())
        preface_flat = " ".join(preface.split())

        required = (
            r"The Heisenberg algebra~$\cH_k$ is the Gaussian support-depth test case",
            r"$r^\Theta_{\max}(\cH_k)=2$",
            r"support packet of $\Theta^{\mathrm{oc}}$ has no primitive component above degree~$2$",
            r"Every scalar projection on this rank-one abelian channel",
            r"$F_g^{\mathrm{sc}}(\cH_k)=k\,\lambda_g^{\mathrm{FP}}$ for $g\geq1$",
            r"completed abelian line comparison models the relevant open-colour line sector",
            r"this is not the chiral Koszul dual",
            r"physical free-boson bulk requires the HT bulk--Hochschild comparison map",
            r"The scalar $k$ determines only the rank-one abelian scalar/OPE data",
        )
        for needle in required:
            self.assertIn(needle, flat)

        retired = (
            r"has shadow depth $r_{\max} = 2$ (class~G)",
            r"so $\Theta^{\mathrm{oc}}$ terminates at degree~2",
            r"genus tower $F_g = k\,\lambda_g^{\mathrm{FP}}",
            r"The line category is $\cC_{\mathrm{line}} \simeq \cH_{-k}\text{-mod}$",
            "the derived center is the free boson bulk",
            r"The single scalar~$k$ determines the Heisenberg package",
            r"Koszulness ($m_{k \ge 3} = 0$) is cleanness",
        )
        for needle in retired:
            self.assertNotIn(needle, source)
        self.assertIn("the derived centre is the algebraic abelian bulk candidate", preface_flat)
        self.assertIn("the physical bulk reading uses the HT bulk--Hochschild comparison", preface_flat)
        self.assertIn("higher support operations vanish", preface_flat)
        self.assertNotIn("bulk $=$ boundary", preface)

    def test_rosetta_ordered_shadow_uses_support_depth_and_log_transport(self):
        source = ROSETTA.read_text()
        section = between(
            source,
            r"\label{comp:heisenberg-e1-ordered-shadow}",
            r"\begin{computation}[Heisenberg classifying space]",
        )
        flat = " ".join(section.split())

        required = (
            "support calculation: the canonical bar-intrinsic support packet",
            r"$r^\Theta_{\max}(\cH_k)=2$",
            "the only primitive support operation is the binary residue",
            r"The support packet terminates: $m^H_j = 0$ for all $j \ge 3$",
            r"Support depth $r^\Theta_{\max}$",
            r"\mathcal R_\hbar(z;z_0)",
            r"\exp(k\hbar\,\xi)",
            r"not a Taylor series in \(1/z\)",
            r"$F_g^{\mathrm{sc}}(\cH_k)=k \cdot \lambda_g^{\mathrm{FP}}$",
            "The determinant-line shadow of the Heisenberg bar complex",
            "not a multi-channel tensor classification",
        )
        for needle in required:
            self.assertIn(needle, flat)

        retired = (
            r"all data reduce to a single scalar~$k$",
            r"the shadow tower terminates at degree~$2$",
            "the entire ordered structure is computed in closed form at every genus",
            r"The Taylor expansion to order~$50$ in $u = 1/z$",
            r"e^{1/z}",
            r"entire function of $1/z$",
            r"Shadow depth $r_{\max}$",
            r"The single scalar~$k$ determines the entire $E_1$ ordered shadow",
            r"The genus tower is $F_g = k \cdot \lambda_g^{\mathrm{FP}}$",
        )
        for needle in retired:
            self.assertNotIn(needle, section)

    def test_rosetta_hydrogen_atom_separates_support_degree_from_cogenerators(self):
        source = ROSETTA.read_text()
        section = between(
            source,
            r"\label{comp:heisenberg-hydrogen-atom}",
            r"\label{comp:heisenberg-e1-ordered-shadow}",
        )
        flat = " ".join(section.split())

        required = (
            r"rank-one scalar/OPE data are controlled by the level~$k$",
            r"canonical support packet has no primitive operation above arity~$2$",
            r"All higher support operations vanish: $m^H_j = 0$ for $j \ge 3$",
            r"The support degree is \(r^\Theta_{\max}(\cH_k)=2\)",
            r"$\barBch(\cH_k)$ is cogenerated in degree~$1$",
            r"higher support packet contributes no correction because $m^H_j=0$",
            r"$F_g^{\mathrm{sc}}(\cH_k)=k\lambda_g^{\mathrm{FP}}$",
            "The determinant-line shadow of the Heisenberg bar complex records",
            r"$(w,d,n)=(0,1,2)$",
        )
        for needle in required:
            self.assertIn(needle, flat)

        retired = (
            r"all data reduce to a single scalar~$k$",
            "all higher complexity vanishes",
            "The complete package follows",
            r"All $m_k = 0$ for $k \ge 3$",
            r"bar cogenerators live in degree~$2$ only",
            r"genus tower contains no shadow corrections because $m_n = 0$",
            r"scalar coefficient $F_g=k\lambda_g^{\mathrm{FP}}$",
            "The Heisenberg bar complex, through the genus tower, knows",
        )
        for needle in retired:
            self.assertNotIn(needle, section)

    def test_heisenberg_fp_formula_is_scalar_lane_across_exports(self):
        sources = {
            "rosetta": ROSETTA.read_text(),
            "preface": PREFACE.read_text(),
            "modular_bootstrap": MODULAR_BOOTSTRAP.read_text(),
        }
        required = r"F_g^{\mathrm{sc}}(\cH_k)=k\lambda_g^{\mathrm{FP}}"
        retired = (
            r"F_g(\cH_k)=k\lambda_g^{\mathrm{FP}}",
            r"F_g(\cH_k) = k\lambda_g^{\mathrm{FP}}",
            r"F_g=k\lambda_g^{\mathrm{FP}}",
        )
        for name, source in sources.items():
            flat = " ".join(source.split())
            self.assertIn(required, flat, name)
            for needle in retired:
                self.assertNotIn(needle, flat, name)

    def test_rosetta_simple_pole_examples_scope_generator_span_vanishing(self):
        source = ROSETTA.read_text()
        section = between(
            source,
            r"\label{subsec:betagamma-ordered-bar}",
            r"\label{sec:shadow-tower-atlas}",
        )
        flat = " ".join(section.split())

        required = (
            r"support depth $r^\Theta_{\max}=4$ by rank-$1$ abelian rigidity",
            r"preserves $r^\Theta_{\max}=4$ on the composite sector",
            r"m^{\mathrm{gen}}_j=0$ for $j\ge3$",
            r"The full support packet reaches $r^\Theta_{\max}=4$ for $\beta\gamma$",
            "the generator-span vanishing statement does not classify the closed composite algebra",
            r"On the generator span it is class~$\mathbf{G}$",
            "It is not a statement about the algebra after closure under normal-ordered composites",
            r"support depth $r^\Theta_{\max}=4$ on composites",
            "This generator-span statement does not classify the closed composite algebra",
            r"support depth $r^\Theta_{\max}=4$ (the ghost current",
            "not the support depth of the closed composite algebra",
            r"full support depth $r^\Theta_{\max}=4$",
            "Generator collision depth and full support depth",
        )
        for needle in required:
            self.assertIn(needle, flat)

        retired = (
            r"All $m_k = 0$ for $k \ge 3$ on generators: the OPE",
            r"All $m_k = 0$ for $k \ge 3$ on generators. The OPE",
            r"On generators, all $m_k = 0$ for",
            r"terminates at depth~$4$",
            r"preserves $r_{\max} = 4$",
            r"The full depth spectrum including composites is $\{0,1,4\}$",
            r"with depth~$4$ on composites",
            r"with depth~$4$ (the ghost current",
            r"quartic exchange at depth~$4$",
            r"Depth spectrum: $\{0\}$ on generators",
            r"Depth spectrum: $\{0\}$ on generators; depth~$4$ full",
        )
        for needle in retired:
            self.assertNotIn(needle, section)

    def test_rosetta_pole_structure_theorem_uses_composite_closed_support(self):
        source = ROSETTA.read_text()
        section = between(
            source,
            r"\begin{theorem}[Composite-closed support-depth stratification]",
            r"\subsection{The genus-1 atlas theorem}",
        )
        flat = " ".join(section.split())

        required = (
            "Composite-closed support-depth stratification",
            "composite-closed chiral generating package",
            r"canonical bar-intrinsic support packet~$\Theta_\cA$",
            r"r^\Theta_{\max}(\cA)",
            "not by the raw pole order of a chosen generator span",
            r"Class~$\mathbf{G}$: the composite-closed support packet is binary",
            r"Class~$\mathbf{L}$: the support packet carries the finite Lie--Jacobi cubic transfer",
            r"Class~$\mathbf{C}$: after composite closure, a quartic contact support component survives",
            r"$\beta\gamma$, $bc$, free fermion after adjoining its",
            r"Class~$\mathbf{M}$: the support packet has no finite cutoff",
            r"\(r^\Theta_{\max}\) is undefined",
            "a simple-pole generator calculation may give",
            "while the composite-closed algebra still has a quartic support component",
        )
        for needle in required:
            self.assertIn(needle, flat)

        retired = (
            "Pole-structure dichotomy",
            r"whose OPE has maximal pole order\/~$d$",
            "The pole structure determines the Swiss-cheese formality",
            r"If\/ $d = 2$",
            r"If\/ $d = 3$",
            r"If\/ $d \ge 4$",
            r"shadow depth\/~$r_{\max}",
            "The pole order $d$ determines the maximal degree",
            "bar cogenerators appear",
            r"the bar complex has non-vanishing $m_k$ for all $k \ge 3$",
            "Shadow depth classifies",
            "shadow depth is undefined",
        )
        for needle in retired:
            self.assertNotIn(needle, section)

    def test_rosetta_support_packet_atlas_defines_composite_closed_depth(self):
        source = ROSETTA.read_text()
        section = between(
            source,
            r"\section{The support-packet atlas}",
            r"\subsection{Virasoro: the infinite tower}",
        )
        flat = " ".join(section.split())

        required = (
            r"\Theta_\cA=\sum_{r\ge2}\Theta_{\cA,r}",
            r"\(\operatorname{Sh}_r(\cA)=\pi_{\mathrm{sc}}(\Theta_{\cA,r})\)",
            "are projections of this composite-closed packet, not raw generator pole orders",
            r"\begin{definition}[Support depth as excess intersection",
            r"r^\Theta_{\max}(\cA)",
            r"\max\{r \ge 2 : \Theta_{\cA,r} \ne 0\}",
            r"\begin{computation}[Support-packet atlas]",
            r"S_r=\pi_{\mathrm{sc}}(\Theta_{\cA,r})/x^r",
            "of the composite-closed support packet",
        )
        for needle in required:
            self.assertIn(needle, flat)

        retired = (
            r"\section{The shadow obstruction tower atlas}",
            r"\begin{definition}[Shadow depth as excess intersection",
            r"The \emph{shadow depth} of a chiral algebra~$\cA$",
            r"\max\{r \ge 2 : \operatorname{Sh}_r(\cA) \ne 0\}",
            r"\begin{computation}[Shadow obstruction tower atlas]",
            "are raw generator pole orders",
        )
        for needle in retired:
            self.assertNotIn(needle, section)


if __name__ == "__main__":
    unittest.main()
