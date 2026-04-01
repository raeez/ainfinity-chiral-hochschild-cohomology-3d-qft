# CLAUDE.md — Volume II: A∞ Chiral Algebras and 3D Holomorphic-Topological QFT

## What the Engine Computes

Volume I built the categorical logarithm — the bar construction B(A) for chiral algebras on curves, with five theorems proving its existence, invertibility, branch structure, leading coefficient, and coefficient ring. Volume II reads the output in three dimensions.

The bar complex carries two structures: a **differential** d_B from OPE residues on FM_k(ℂ), encoding the holomorphic chiral product, and a **coproduct** Δ from ordered deconcatenation on Conf_k(ℝ), encoding the topological interval-cutting. The differential lives in the ℂ-direction; the coproduct lives in the ℝ-direction. Together, a bar element of degree k is parametrized by FM_k(ℂ) × Conf_k(ℝ) — the product of holomorphic and topological configuration spaces.

This product is the operadic fingerprint of a 3d holomorphic-topological QFT on ℂ_z × ℝ_t, where observables factorize holomorphically in z and associatively in t. The two-colored Swiss-cheese operad SC^{ch,top} has operation spaces FM_k(ℂ) × E₁(m). The bar differential is the closed (holomorphic) color. The bar coproduct is the open (topological) color. The no-open-to-closed rule reflects that bulk interactions restrict to boundaries but not conversely. **The bar complex presents the Swiss-cheese algebra, as the Steinberg variety presents the Hecke algebra.**

**CRITICAL DISTINCTION (AP-OC):** The bar complex classifies *twisting morphisms* (universal couplings between A and A!). The bulk algebra — the observables of the 3d theory — is a DIFFERENT object: the chiral derived center Z^der_ch(A) = C^•_ch(A_b, A_b), realized as chiral Hochschild cochains of any boundary chart A_b (Vol I, thm:thqg-swiss-cheese). Bar classifies couplings; derived center classifies bulk observables. The primitive object is the open-sector factorization dg-category C_op; the boundary algebra A_b = End(b) is a chart, not an essence (Morita invariance: thm:thqg-local-global-bridge(iii)).

At genus g ≥ 1: curved Swiss-cheese with curvature κ(A)·ω_g from the Hodge bundle. The non-vanishing of higher A∞ operations IS the curved bar structure d² = κ(A)·ω₁ — formality fails precisely because the logarithm acquires monodromy.

## The Monograph

Two volumes by Raeez Lorgat. Vol I (~2,200pp, ~/chiral-bar-cobar) proves the machine. Vol II (~900pp, this repo) shows what it computes and how to read the output.

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

## Ten Parts (Treatise Architecture)

**I. The Open Primitive** (part:swiss-cheese). The primitive datum is a category, not an algebra: the open/closed factorization dg-category C on a tangential log curve. SC^{ch,top} constructed, recognition theorem, homotopy-Koszulity proved. Chapters: foundations, locality, axioms, equivalence, BV construction, factorization Swiss-cheese.

**II. The Universal Bulk** (part:bbl-core). The bulk is the chiral derived center Z^der_ch(A), not the bar complex. Chiral Hochschild cochains, brace dg algebra, bar-cobar review, line operators, spectral braiding (proved core), Koszul triangle (proved core), celestial boundary transfer (proved core).

**III. Modularity as Trace and Clutching** (part:modularity). Modularity arises from trace and clutching on the open sector, not as an axiom on the closed algebra. Modular Swiss-cheese operad, relative Feynman transform, modular PVA quantization (core), physical origins.

**IV. Descent and the Classical Shadow** (part:descent). Extracts classical PVA from the quantum SC structure. Raviolo VA, raviolo restriction, FM calculus, PVA descent (D2-D6 all proved).

**V. The Standard HT Landscape** (part:examples). Worked examples: Rosetta stone, free multiplet, LG, CS (proved), Virasoro, W_3.

**VI. Three-Dimensional Quantum Gravity** (part:quantization). Modular Koszul duality of Virasoro IS 3d quantum gravity. Affine half-space BV, planted-forest synthesis, gravitational complexity, 3d gravity movements, perturbative finiteness, soft graviton theorems, critical string dichotomy, symplectic polarization.

**VII. The Ordered Sector and Factorization Transport** (part:ordered). The E_1 wing: ordered associative chiral KD (core), dg-shifted factorization bridge, gravitational Yangian.

**VIII. Holographic and Celestial Frontier** (part:holography). YM synthesis (core), celestial holography (core), logarithmic HT monodromy (core), anomaly-completed holography (core), modular bootstrap.

**IX. Extensions, Conditional Results, and Frontier** (part:frontier). All frontier/conjectural material from chapter splits. No earlier part depends on this part.

**Conclusion and Aftermatter** (unnumbered). Conclusion, appendices (brace signs, orientations, FM proofs, PVA expanded).

## Standing Hypotheses — MADE EXPLICIT

**The algebraic framework is unconditional.** The former standing hypotheses (H1)–(H4) are no longer background axioms. They have been made explicit: (H1)–(H2) are now conditions of the physics bridge theorem (Theorem thm:physics-bridge), (H3) is a theorem, and (H4) is the recognition theorem. The pure-algebraic theory requires only a single definition:

**Definition (Logarithmic SC^{ch,top}-algebra, Definition def:log-SC-algebra):** A C_*(W(SC^{ch,top}))-algebra whose closed-colour A∞ operations are defined by logarithmic weight forms factoring as ω_k = ω_k^hol ⊗ ω_k^top on FM_k(ℂ) × Conf_k(ℝ).

The former axioms are now derived:

| | Content | New Status |
|---|---------|--------|
| (H1) | BV data, one-loop finiteness | Condition of the bridge theorem (Theorem thm:physics-bridge) — applies only to physical realisations |
| (H2) | Propagator: meromorphic in ℂ, exponential decay in ℝ | Consequence of Q = ∂̄ + d_t splitting — Green's function is Cauchy kernel × Heaviside |
| (H3) | FM compactification, logarithmic forms, AOS relations, Stokes exactness | Theorem of configuration space geometry (Theorem thm:FM-calculus) |
| (H4) | Factorization compatibility with C_*(W(SC^{ch,top})) | Recognition theorem (Theorem thm:recognition-SC) — already proved |

All results in Parts I–VII hold unconditionally for any logarithmic SC^{ch,top}-algebra. Physical theories (gauge theories satisfying Theorem thm:physics-bridge) provide the standard class of examples.

## Beilinson Rectification Loop

The Vol I `CLAUDE.md` supplies the source protocol. For Vol II, the Codex-native enactment now lives locally in `AGENTS.md`: use `update_plan` for loop state, local RED/BLUE/GREEN passes by default, `multi_tool_use.parallel` for shell fan-out, `exec_command` for build/test/worktree steps, `apply_patch` for fixes, and `compute/audit/linear_read_notes.md` as the findings register. Run the Beilinson loop on the live Vol II surface or on a named target chapter, and treat convergence as "no actionable findings at severity MODERATE or above after re-audit plus the relevant verification passes."

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
- The corrected bulk/boundary/line triangle: **bulk ≃ derived CENTER of boundary**, NOT bulk = boundary. **Scope:** Proved rigorously in the boundary-linear exact sector (Theorem thm:boundary-linear-bulk-boundary). The global triangle for all HT theories is conjectural — it reduces formally to hypotheses (compact generator, derived center quasi-isomorphism) that are verified only in the boundary-linear case and for inputs from CDG20/GKW24.
- Chiral Koszulness from physics: **RESOLVED for the affine lineage** by the loop-order criterion (Theorem thm:one-loop-koszul), which shows one-loop exactness of the BV-BRST differential implies Koszulness. DS reduction (Theorem thm:ds-koszul-obstruction) explains the failure mechanism: Koszulness descends through DS iff the BRST cohomology inherits one-loop exactness. General case beyond the affine lineage remains open.
- **Homotopy-Koszulity of SC^{ch,top} is PROVED** (Theorem thm:homotopy-Koszul): via Kontsevich formality + transfer from classical Swiss-cheese (Livernet). ALL formerly conditional results (bar-cobar Quillen equivalence, filtered Koszul duality, C_line ≃ A!-mod, dg-shifted Yangian) are now unconditional.
- **Spectral Drinfeld strictification is PROVED** (Theorem thm:complete-strictification in dg_shifted_factorization_bridge.tex) for all filtrations and all simple Lie algebras: root multiplicity one forces the spectral Drinfeld class to vanish at every filtration, via root-space one-dimensionality (Theorem thm:root-space-one-dim) and the Jacobi collapse lemma (Lemma lem:jacobi-collapse). The true remaining frontier is Kac-Moody algebras with root multiplicities > 1, where root-space one-dimensionality fails and the strictification mechanism requires new input.
- The **Koszul dual is the boundary**, not the bulk: A! lives on the boundary ℝ, not in the bulk ℂ × ℝ.

## Anti-Patterns (inherited from Vol I + Vol II-specific)

All Vol I anti-patterns AP1-AP32 apply here. The following encode deep mathematical errors found propagated across both volumes:

**AP-OC — Bar = bulk conflation.** The bar complex B(A) classifies TWISTING MORPHISMS (universal couplings between A and A!). The BULK OBSERVABLES are a DIFFERENT object: the chiral derived center Z^der_ch(A) = C^•_ch(A_b, A_b). The PRIMITIVE OBJECT is the open-sector factorization dg-category C_op; the boundary algebra A_b = End(b) is a CHART (Morita-invariant). Modularity belongs to TRACE + CLUTCHING on the open sector, not to the closed algebra alone. Found at: ht_physical_origins.tex lines 359-361 (corrected), celestial_holography_frontier.tex line 846 (corrected), thqg_bv_ht_extensions.tex lines 218/302/704 (corrected), plus 4 superseded files.

**AP19 — The bar kernel absorbs a pole.** The collision residue r(z) has pole orders ONE LESS than the OPE — the d log(z-w) kernel absorbs one power. Virasoro r-matrix: (c/2)/z³ + 2T/z, NOT (c/2)/z⁴ + 2T/z² + ∂T/z.

**AP20 — An invariant of one algebra is not an invariant of a system.** κ(A) is intrinsic; κ_eff = κ(matter) + κ(ghost) is composite; κ(B) where B = A! is the dual's. F_g always uses κ(A). State WHICH algebra's κ.

**AP21 — A class is not a scalar; Clifford ≠ exterior.** u = η² = λ = κ(B)·ω_g is LINEAR in κ (a class), NOT κ² (a scalar). The gravity dichotomy (c≠26 vs c=26) is whether this class vanishes. Squaring κ destroys the bifurcation.

**AP22 — Generating function index mismatch.** Â(iℏ)-1 starts at ℏ². If F_1 ≠ 0, the pairing Σ F_g ℏ^{power} must use ℏ^{2g} (not ℏ^{2g-2}) unless an explicit 1/ℏ² appears.

**AP23 — Flat section vs weighted transport.** √(Q_L) is the flat section of ∇^sh. The shadow generating function H(t) = t²√(Q_L) is NOT horizontal — the t² is the arity offset.

**AP24 — The complementarity sum is not universally zero.** κ + κ! = 0 for KM/free fields/lattice/principal W. For Virasoro: κ + κ! = 13 ≠ 0. The anti-symmetry was overclaimed in 20+ locations. The Feigin-Frenkel involution ensures anti-symmetry for KM; the Virasoro involution c ↦ 26−c is anti-symmetric around 13, not 0.

**AP25 — Three functors, three outputs: bar ≠ Verdier dual ≠ cobar.** B(A) = coalgebra. D_Ran(B(A)) ≃ B(A!) = Verdier dual (an ALGEBRA). Ω(B(A)) ≃ A = cobar (recovers the ORIGINAL). The Koszul dual A! is obtained by Verdier duality, NOT by cobar. Found conflated in 16 files.

**AP26 — Fock inner product ≠ BPZ inner product.** At weight ≥ 4 for rank ≥ 3 W-algebras, dim(Fock) > dim(W-algebra) and ⟨W₄|Λ⟩_Fock ≠ 0 even though ⟨W₄|Λ⟩_BPZ = 0. W-algebra decompositions MUST use BPZ (Wick contractions), not the free-field dot product.

**Meta-principle:** Every error traces to confusing two objects that share a name, formula, or special-case coincidence. The meta-rule: never trust a coincidence. Verify at the most general case, highest weight, most general level, most general family.

## Cross-Volume Bridges

| Bridge | Vol II claim | Vol I anchor | Status |
|--------|-------------|--------------|--------|
| Bar-cobar | SC^{ch,top} bar-cobar specializes Vol I Thm A when curve = ℂ, topological = ℝ | Theorem A | Proved |
| DS-bar | Bar-cobar commutes with DS reduction | Theorem ds-koszul-intertwine | Proved (Vol I) |
| Hochschild | BV-BRST origin of Vol I's Theorem H complex | Theorem H | Proved (all genera) |
| DK/YBE | r(z) = ∫₀^∞ e^{-λz}{·_λ·}dλ provides DK-0 shadow | MC3 (DK extension) | Proved (Laplace) |
| PVA-Coisson | PVA descent at X = pt recovers Coisson structure | Deformation theory | Proved |
| W-algebras | Feynman-diagrammatic m_k matches bar differential at all genera | MC5 (BRST = bar) | **Proved** (all genera) |
| Affine monodromy | Reduced HT monodromy = quantum group R-matrix; C_line^red ≃ Rep_q(𝔤) on eval modules; Jones polynomial from bar complex | Thm A + affine half-space BV + Drinfeld-Kohno | **Proved** (affine lineage) |

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

**Theory** (chapters/theory/, 10 files): Parts I–II.
- foundations, locality, axioms, equivalence, bv-construction: Part I (Swiss-Cheese)
- raviolo, raviolo-restriction: Part I (promoted from appendix)
- fm-calculus, pva-descent-repaired: Part II (Descent Calculus)
- introduction: global introduction

**Examples** (chapters/examples/): Part IV.
- rosetta_stone, examples-computing, examples-worked: proved core
- examples-complete-proved: proved computations (split from examples-complete)
- examples-complete-conditional: conditional computations (split, Part VIII)
- w-algebras-stable: general framework (split from w-algebras)
- w-algebras-conditional: Virasoro/W₃ conditional (split, Part VIII)

**Connections — Core** (used in Parts III, V, VI, VII):
- hochschild, brace, bar-cobar-review, line-operators: Part III core
- spectral-braiding-core, ht_bulk_boundary_line_core, celestial_boundary_transfer_core: Part III (split)
- ht_physical_origins: Part III (merged from physical_origins + holomorphic_topological + bv_ht_physics)
- modular_pva_quantization_core, affine_half_space_bv, fm3_planted_forest_synthesis, 3d_gravity: Part V
- ordered_associative_chiral_kd_core, dg_shifted_factorization_bridge: Part VI
- ym_synthesis_core, celestial_holography_core, log_ht_monodromy_core, anomaly_completed_core: Part VII

**Connections — Frontier** (Part VIII, all splits):
- spectral-braiding-frontier, ht_bulk_boundary_line_frontier, celestial_boundary_transfer_frontier
- modular_pva_quantization_frontier, ordered_associative_chiral_kd_frontier
- ym_synthesis_frontier, celestial_holography_frontier, log_ht_monodromy_frontier, anomaly_completed_frontier

**Aftermatter**: conclusion (Part IX). Concordance removed (external compile route, matching Vol I).

**Appendices**: brace-signs, orientations, fm-proofs, pva-expanded-repaired

**Superseded files** (still in repo, no longer \input'd):
- spectral-braiding.tex, ht_bulk_boundary_line.tex, celestial_boundary_transfer.tex, modular_pva_quantization.tex
- ordered_associative_chiral_kd.tex, ym_synthesis.tex, celestial_holography.tex, log_ht_monodromy.tex
- anomaly_completed_topological_holography.tex, examples-complete.tex, w-algebras.tex
- physical_origins.tex, holomorphic_topological.tex, bv_ht_physics.tex, concordance.tex

## The Aesthetic

Show, don't tell. Every construction should feel inevitable — a single phenomenon viewed from different angles, not parts assembled. The Steinberg variety presents the Hecke algebra; the bar complex presents the Swiss-cheese algebra; the transgression algebra presents the holographic dictionary; the complementarity potential presents the nonlinear modular shadows. Synthesize disparate mathematical domains to bring out their inner wonder and music.
