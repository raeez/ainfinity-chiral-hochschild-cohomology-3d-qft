# The integrated monograph

`make` → `out/platonic.pdf` — **495 pages, 0 LaTeX errors, 0 undefined
references, 0 undefined citations.** Default target: a plain build produces the
current, merged monograph.

## This repository

**Volume II — Mixed Holomorphic–Topological Theory**: Hamiltonian jets, shifted cotangent fields, universal branes, Duflo quantization, and the all-arity open–closed map. Placed first.

## The programme, and what stands outside it

| | subject | repository |
|---|---|---|
| **I** | Ordered Chiral Geometry | `chiral-bar-cobar` |
| **II** | Mixed Holomorphic–Topological Theory | `chiral-bar-cobar-vol2` |
| **III** | Calabi–Yau Quantum Groups | `calabi-yau-quantum-groups` |
| — | Universal Chiral BV and Einstein Completion | no repository yet |
| — | *The Igusa Square Root* | `igusa-cusp-form`, standalone |
| — | *Remarks on Mixed HT Strings* | `mixed-holomorphic-topological-strings`, standalone |

The upstream spine ships five files named `Volume_I` … `Volume_V`, with
Igusa–Borcherds at `III` and Calabi–Yau at `IV`. **Those filenames are
upstream identifiers, not the programme numbering.** They are left untouched so
the files stay byte-identical across repositories; do not renumber them.

## One spine

The five volume files, `platonic.sty`, `integrated_macros.tex`, and
`references.bib` are byte-identical in every repository. Repositories differ
only in their own frame chapter and in which volume they read first. Editing a
volume file here forks it from every other copy.

## Merged in — `Volume_I_audit_and_obstructions.tex`

The spine as shipped had no audit of the legacy planar tables and no
obstruction constitution: `Motzkin`, `Riordan`, `Goncharova`,
`pentagonal`, and `no-go` all occurred zero times across all five volumes.
Now merged in, in the spine's own notation: the bar–Chevalley comparison with
$\mathfrak{sl}_2\to1+t^3$, $\mathfrak{sl}_3\to1+t^3+t^5+t^8$,
$\mathfrak g_2\to1+t^3+t^{11}+t^{14}$; pentagonal rigidity
($\dim H^n(\mathrm{Bar}^\perp U(L_1))=2$ in the weights $(3n^2\mp n)/2$,
Euler character $\prod(1-q^n)$ saturating it, verified $n=1,2,3,4$); the
identity $M(n)=R(n)+R(n+1)$ dissolving the Drinfeld–Sokolov attribution; the
carrier theorem; the location of a mixed distributive law on the aligned
idempotent's image; and seven obstruction theorems, the last of which retracts
the $5\times5$ $\kappa$-matrix and the Universal Trace Identity rather than
merely omitting them.

## What is discarded, and what is kept for the record

The retracted architecture — Open Beilinson tower, Theorems A/B/C/D/H, the five
archetypes, the $5\times5$ $\kappa$-matrix, the Universal Trace Identity — is
**not** in the default build. It remains tracked, still builds under its own
target, and its disposition is recorded file by file. Nothing was deleted; it
is out of the monograph and available as history.

## Verification

`chiral-bar-cobar/compute/` carries two independent harnesses, cross-checked:
exact `Fraction` sparse elimination with $d^2=0$ asserted per space, and
sympy dense matrices. They share one claim, $H^*(\mathfrak{sl}_2)=1,0,0,1$,
and agree with no shared code, basis, sign convention, or backend. 23 tests.

## Still open

No $\lambda$ is constructed for any example, so the doubly noncommutative
object remains a definition without a model; the anticommutation of the two
edge-contraction differentials rests on determinant-line independence that no
witness proves; the genus-two-and-above analytic material is unported; and the
408-page revision still lacks 46 of its 53 sources.
