# CLAUDE.md — Volume II: A∞ Chiral Algebras and 3D Holomorphic-Topological QFT

## What the Engine Computes

Volume I built the categorical logarithm — the bar construction B(A) for chiral algebras on curves, with five theorems proving its existence, invertibility, branch structure, leading coefficient, and coefficient ring. Volume II reads the output in three dimensions.

The bar complex carries two structures: a **differential** d_B from OPE residues on FM_k(ℂ), encoding the holomorphic chiral product, and a **coproduct** Δ from ordered deconcatenation on Conf_k(ℝ), encoding the topological interval-cutting. The differential lives in the ℂ-direction; the coproduct lives in the ℝ-direction. Together, a bar element of degree k is parametrized by FM_k(ℂ) × Conf_k(ℝ) — the product of holomorphic and topological configuration spaces.

This product is the operadic fingerprint of a 3d holomorphic-topological QFT on ℂ_z × ℝ_t, where observables factorize holomorphically in z and associatively in t. The two-colored Swiss-cheese operad SC^{ch,top} has operation spaces FM_k(ℂ) × E₁(m). The bar differential is the closed (holomorphic) color. The bar coproduct is the open (topological) color. The no-open-to-closed rule reflects that bulk interactions restrict to boundaries but not conversely. **The bar complex presents the Swiss-cheese algebra, as the Steinberg variety presents the Hecke algebra.**

At genus g ≥ 1: curved Swiss-cheese with curvature κ(A)·ω_g from the Hodge bundle. The non-vanishing of higher A∞ operations IS the curved bar structure d² = κ(A)·ω₁ — formality fails precisely because the logarithm acquires monodromy.

## The Monograph

Two volumes by Raeez Lorgat. Vol I (~1,960pp, ~/chiral-bar-cobar) proves the machine. Vol II (~450pp, this repo) shows what it computes and how to read the output.

Every theorem proved, every physical identification precise, every construction functorial. When claims outrun proofs, strengthen the proof first. Target: Annals/Astérisque grade.

## Vol I Theorems Used Here

Every chapter depends on Vol I's five theorems. Cross-references to Vol I labels resolve as "undefined" — expected for a multi-volume work.

| Vol I Theorem | What it supplies to Vol II |
|---------------|---------------------------|
| **(A)** Bar-cobar adjunction | The bar complex exists as a factorization coalgebra; reinterpreted as SC^{ch,top}-algebra structure — the bridge theorem |
| **(B)** Koszul inversion | Bar-cobar equivalence on the Koszul locus; lifted to raviolo VA setting and completed towers |
| **(C)** Complementarity | Genus-g obstructions decompose as complementary Lagrangians; the bulk-boundary-line triangle inherits this (−1)-shifted symplectic structure |
| **(D)** Leading coefficient | Curvature κ(A)·ω_g governs the genus tower; curved Swiss-cheese = Swiss-cheese + Hodge deformation |
| **(H)** Hochschild ring | BV-BRST origin of the deformation ring; bulk ≃ chiral Hochschild (Theorem H gets its physical explanation) |

## Five Parts

**I. From the Bar Complex to the Swiss-Cheese Operad.** The bridge from Vol I to three dimensions. SC^{ch,top} is constructed: closed color = FM_k(ℂ), open color = FM_k(ℂ) × E₁(m). The closed color recovers BD chiral algebras (ProvedHere). Recognition theorem: an A∞ chiral algebra satisfying (H1)–(H4) is an SC^{ch,top}-algebra (ProvedHere). Homotopy-Koszulity of SC^{ch,top} proved via Kontsevich formality + transfer from classical Swiss-cheese (ProvedHere). Raviolo descent: SC-algebra → raviolo vertex algebra → Poisson vertex algebra on cohomology (ProvedHere). The raviolo is the algebraic avatar of the ℂ × ℝ geometry. Axiomatics: sesquilinearity, A∞ relations with spectral substitution, cluster factorization from Stokes on FM compactifications.

**II. The Descent Calculus.** Two descent mechanisms: cohomological (bar → PVA via Arnold/Stokes) and genus (curved bar over M̄_g). PVA on cohomology H•(A,Q) is a (−1)-shifted Poisson vertex algebra — all axioms verified (ProvedHere). FM calculus, chiral Hochschild cohomology, and the bar-cobar review (Quillen equivalence, filtered Koszul duality).

**III. Dualities and the Bulk-Boundary-Line Triangle.** The **bulk-boundary-line Koszul triangle**: bulk algebra A, boundary algebra A! (Koszul dual), line operators as A!-modules — three vertices of a single triangle, with Koszul duality mediating each edge. Corrected form: bulk ≃ derived center of boundary, NOT bulk = boundary. Spectral braiding: R(z) from bulk-boundary composition solves Yang-Baxter via Stokes on FM₃; classical limit r(z) = Laplace of λ-bracket. Celestial transfer: the boundary PVA descends to a celestial OPE algebra.

**III. The Standard Landscape.** Three worked examples verify (H1)–(H4) and ground the abstract machinery:
- *Free multiplet*: all m_{k≥3} = 0 (no interaction vertices). The trivial case that calibrates signs.
- *Landau-Ginzburg cubic*: truncation at m₃ by degree counting. The simplest interacting example.
- *Abelian Chern-Simons*: boundary = û(1)_k, OPE computed directly from propagator.
- Plus: *Virasoro* (genuinely infinite A∞ -- all m_k != 0 via wheel diagrams, Q² = 0, recursive m_k from BV master equation), *W₃* (classical YBE from λ-bracket Jacobi, central charge shift from ghost counting).

**IV. The Standard Landscape.** Three worked examples verify (H1)–(H4) and ground the abstract machinery: free multiplet, LG cubic, abelian CS, plus Virasoro and W₃.

**V. Quantization and Holography.** PVA quantization via the modular bar complex. Affine half-space BV (solved for the affine case). Planted-forest L∞ obstructions. YM boundary packages. Celestial holography. Logarithmic HT monodromy. The **anomaly-completed Koszul triangle**: transgression algebra B_Θ, secondary anomaly u = η², genus-Clifford dichotomy — the holographic dictionary is presented by the transgression algebra, as the Swiss-cheese algebra is presented by the bar complex.

## Standing Hypotheses

**No conjectural algebraic inputs remain.** Both Recognition and Homotopy-Koszulity are proved. All algebraic results are unconditional. The only conditional inputs are four analytic axioms, verified in the three worked examples but not proved in general:

| | Content | Status |
|---|---------|--------|
| (H1) | BV data, HT gauge fixing, one-loop finiteness | Physical axiom |
| (H2) | Propagator: meromorphic in ℂ, exponential decay in ℝ | Physical axiom |
| (H3) | FM compactification, logarithmic forms, AOS relations, Stokes exactness | Physical axiom |
| (H4) | Factorization compatibility with C_*(W(SC^{ch,top})) | Physical axiom |

## Critical Pitfalls

**Inherited from Vol I (NEVER VIOLATE):**
- Four objects: A (algebra), B(A) (bar coalgebra), A^i (dual coalgebra), A^! (dual algebra). NEVER conflate.
- Ω(B(A)) = A is inversion. A^! via Verdier duality. Cobar does NOT produce A^!.
- COHOMOLOGICAL grading (|d| = +1). Bar uses DESUSPENSION. m₁²(a) = [m₀, a] (commutator, MINUS sign). Bar d² = 0 always.
- Sugawara UNDEFINED at critical level k = −h∨ (not "c diverges"). Feigin-Frenkel: k ↔ −k−2h∨.
- Virasoro self-dual at c = 13, NOT c = 26. Vir_c^! = Vir_{26−c}.

**Specific to Vol II:**
- Swiss-cheese **directionality is strict**: open-to-closed is EMPTY. No open inputs produce closed outputs. This is the mathematical expression of bulk→boundary directionality.
- PVA on cohomology H•(A,Q) is **(−1)-shifted**: the λ-bracket has shifted parity relative to classical PVA conventions.
- The R-matrix R(z) comes from **bulk-boundary composition**, NOT from the universal R-matrix of a quantum group (though they agree on the evaluation locus — this is DK-0).
- Formality fails at d' = 1: this is NOT a defect. The non-vanishing of higher A∞ operations IS the curved bar structure d² = κ(A)·ω₁ from Vol I.
- The corrected bulk/boundary/line triangle: **bulk ≃ derived CENTER of boundary**, NOT bulk = boundary.
- Chiral Koszulness from physics is **OPEN**: (H1)–(H4) + BV-BRST should imply it, but this is not proved.
- **Homotopy-Koszulity of SC^{ch,top} is PROVED** (Theorem thm:homotopy-Koszul): via Kontsevich formality + transfer from classical Swiss-cheese (Livernet). ALL formerly conditional results (bar-cobar Quillen equivalence, filtered Koszul duality, C_line ≃ A!-mod, dg-shifted Yangian) are now unconditional.
- The **Koszul dual is the boundary**, not the bulk: A! lives on the boundary ℝ, not in the bulk ℂ × ℝ.

## Cross-Volume Bridges

| Bridge | Vol II claim | Vol I anchor | Status |
|--------|-------------|--------------|--------|
| Bar-cobar | SC^{ch,top} bar-cobar specializes Vol I Thm A when curve = ℂ, topological = ℝ | Theorem A | Proved |
| DS-bar | Bar-cobar commutes with DS reduction | Theorem ds-koszul-intertwine | Proved (Vol I) |
| Hochschild | BV-BRST origin of Vol I's Theorem H complex | Theorem H | Proved (all genera) |
| DK/YBE | r(z) = ∫₀^∞ e^{-λz}{·_λ·}dλ provides DK-0 shadow | MC3 (DK extension) | Proved (Laplace) |
| PVA-Coisson | PVA descent at X = pt recovers Coisson structure | Deformation theory | Proved |
| W-algebras | Feynman-diagrammatic m_k matches bar differential at all genera | MC5 (BRST = bar) | **Proved** (all genera) |

## Build

```
pkill -9 -f pdflatex 2>/dev/null || true; sleep 2; make    # Full build (single pass usually suffices)
```

Same engine as Vol I: memoir, EB Garamond, newtxmath, thmtools, microtype.

## LaTeX Rules

- Use `\providecommand` (not `\newcommand`) for macros — many pre-defined in the compatibility block
- Do NOT add `\newtheorem` in chapter files — all theorem environments defined in main.tex preamble
- Claim tags match Vol I: \ClaimStatusProvedHere, \ClaimStatusProvedElsewhere, \ClaimStatusConjectured, \ClaimStatusHeuristic, \ClaimStatusOpen
- Key macros: \cA, \Ainf, \Linf, \barB, \Omegach, \hh, \HH, \Sym, \End
- Chapters moved from Vol I retain mathematical content verbatim; cite keys mapped to Vol II bibliography
- Do not add packages without checking preamble
- Do not create new files when content belongs in existing chapters
- Do not duplicate definitions from Vol I — reference them textually

## Git — HARD RULE

All commits authored by Raeez Lorgat. **Never credit an LLM.** No "co-authored-by", no "generated by", no AI attribution anywhere in the repo.

## File Map

**Theory** (chapters/theory/, 10 files): The operadic foundations and equivalence theorems.
- foundations: **the bridge** — bar complex as Swiss-cheese algebra, SC^{ch,top} operad construction
- raviolo, raviolo-restriction: algebraic avatar of ℂ × ℝ geometry, raviolo VA from SC-algebra
- pva-descent, pva-preview: PVA on cohomology, (−1)-shifted Poisson structure, all axioms verified
- axioms: sesquilinearity, A∞ relations with spectral substitution, cluster factorization
- equivalence: operad ⟹ axioms and axioms ⟹ operad (rectification)
- bv-construction: chain-level A∞ from BV-BRST, conditional on (H1)–(H4)
- fm-calculus: A∞ from Stokes on FM compactifications, AOS corner cancellations
- locality: HT prefactorization structure

**Examples** (chapters/examples/, 4 files): Free multiplet, LG, CS, Virasoro, W₃ — complete computations.

**Connections** (chapters/connections/, 19 files): The programme and frontier.
- concordance: **status ledger** and cross-volume bridges (constitutional for Vol II)
- bar-cobar-review: bar-cobar adjunction and filtered Koszul duality in the SC setting
- hochschild: brace algebra on HC cochains, Gerstenhaber bracket, bulk ≃ chiral Hochschild
- line-operators: C_line ≃ A!-mod (unconditional; homotopy-Koszulity proved)
- spectral-braiding: R(z) from bulk-boundary, Yang-Baxter from Stokes on FM₃
- ht_bulk_boundary_line: the Koszul triangle — bulk, boundary, line
- celestial_holography, celestial_boundary_transfer: celestial OPE and boundary transfer
- anomaly_completed_topological_holography: B_Θ, secondary anomaly u = η², genus-Clifford
- log_ht_monodromy: logarithmic monodromy in the HT direction
- modular_pva_quantization, affine_half_space_bv: quantization programme
- fm3_planted_forest_synthesis: planted-forest L∞ obstructions
- bv_ht_physics, physical_origins, holomorphic_topological: physics bridges
- ym_synthesis, brace: YM boundary, brace structures
- conclusion: synthesis and open problems

## The Aesthetic

Show, don't tell. Every construction should feel inevitable — a single phenomenon viewed from different angles, not parts assembled. The Steinberg variety presents the Hecke algebra; the bar complex presents the Swiss-cheese algebra; the transgression algebra presents the holographic dictionary; the complementarity potential presents the nonlinear modular shadows. Synthesize disparate mathematical domains to bring out their inner wonder and music.
