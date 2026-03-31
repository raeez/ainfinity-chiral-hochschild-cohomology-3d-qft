# AGENTS.md — Volume II: A∞ Chiral Algebras and 3D Holomorphic-Topological QFT

## What the Engine Computes

Volume I built the categorical logarithm — the bar construction B(A) for chiral algebras on curves, with five theorems proving its existence, invertibility, branch structure, leading coefficient, and coefficient ring. Volume II reads the output in three dimensions.

The bar complex carries two structures: a **differential** d_B from OPE residues on FM_k(ℂ), encoding the holomorphic chiral product, and a **coproduct** Δ from ordered deconcatenation on Conf_k(ℝ), encoding the topological interval-cutting. The differential lives in the ℂ-direction; the coproduct lives in the ℝ-direction. Together, a bar element of degree k is parametrized by FM_k(ℂ) × Conf_k(ℝ) — the product of holomorphic and topological configuration spaces.

This product is the operadic fingerprint of a 3d holomorphic-topological QFT on ℂ_z × ℝ_t, where observables factorize holomorphically in z and associatively in t. The two-colored Swiss-cheese operad SC^{ch,top} has operation spaces FM_k(ℂ) × E₁(m). The bar differential is the closed (holomorphic) color. The bar coproduct is the open (topological) color. The no-open-to-closed rule reflects that bulk interactions restrict to boundaries but not conversely. **The bar complex presents the Swiss-cheese algebra, as the Steinberg variety presents the Hecke algebra.**

At genus g ≥ 1: curved Swiss-cheese with curvature κ(A)·ω_g from the Hodge bundle. The non-vanishing of higher A∞ operations IS the curved bar structure d² = κ(A)·ω₁ — formality fails precisely because the logarithm acquires monodromy.

## The Monograph

Two volumes by Raeez Lorgat. Vol I (~1,960pp, ~/chiral-bar-cobar) proves the machine. Vol II (~450pp, this repo) shows what it computes and how to read the output.

Every theorem proved, every physical identification precise, every construction functorial. When claims outrun proofs, strengthen the proof first. Target: Annals/Astérisque grade.

## Research Audit Doctrine — The Beilinson Principle

This repository is in a long correction-and-rectification phase. The default stance for every agent is frontier mathematical skepticism, not trust in inherited prose.

**Core doctrine:** the main obstacle to progress is not lack of ideas but failure to discard false ones quickly enough. In this repo, confidence, elegance, repetition, and prior status labels are never evidence. A claim is true only after it has been traced to ground.

What counts as tracing to ground:
- a complete local proof whose nontrivial steps can be named and checked;
- an exact citation to a result that really states the needed claim;
- a computation or test that genuinely verifies the mathematical statement being used.

Repository-wide rules:
- Every claim is provisional until verified, even if the text states it confidently.
- Downgrading a statement, narrowing hypotheses, or isolating a missing lemma is progress.
- Always distinguish sharply between: proved here, proved elsewhere, consequence of Definition~`def:log-SC-algebra`, conditional on Theorem~`thm:physics-bridge`, conjectural, heuristic, and open.
- Treat contradictions between theorem statements, surrounding prose, examples, tests, appendices, build logs, and split/superseded files as mathematical bugs.
- Prefer the live manuscript surface (`main.tex`, files currently `\input`'d, current diffs, actual test/build output) over stale notes or historical summaries.
- Never let tone hide a gap. If a proof does not close, either strengthen it, weaken the claim, mark the missing input explicitly, or move the statement to the frontier.

### Default Audit Protocol

For any nontrivial task, follow this order:
1. Identify the exact statement, labels, hypotheses, and dependencies.
2. Check conventions and invariants first: grading, signs, bar/cobar objects, Swiss-cheese directionality, shifted PVA parity, and Vol I/Vol II terminology.
3. Trace the argument to ground line by line. Do not compress "standard" steps unless they are genuinely standard and correctly cited.
4. Search for collisions across the repo: introduction, chapter-local restatements, examples, compute code, tests, appendices, and superseded split files.
5. Propagate corrections everywhere the statement is used or advertised.
6. Run the narrowest relevant verification before declaring success: targeted `pytest`, targeted grep for labels/refs, `make fast`, or log inspection.

### Specific Failure Modes To Resist

- Do not mistake a status tag for a proof.
- Do not infer truth from repeated appearance across multiple files; duplication can amplify error.
- Do not trust stale programme notes over current source, git history, and executable verification.
- Do not preserve a beautiful but false formulation when a narrower true statement is available.
- Do not "repair" by rhetoric alone; if the mathematical content changes, update theorem status, hypotheses, and downstream uses.

## Deep Beilinson Audit — Codex Procedure

When instructed to perform a deep Beilinson audit on Vol II, execute the loop on the **live manuscript surface**: `main.tex`, the files currently `\input`'d, the dirty git diff, the active build logs, and the narrowest relevant compute/tests slice.

### Core stance

- Falsification comes first. The task is to break claims before trusting them.
- Treat `AGENTS.md` and `CLAUDE.md` as operational guides, not mathematical ground truth.
- Prefer a smaller true statement to a larger false one.
- If a correction is not independently verified, mark the gap explicitly rather than silently "repairing" it.

### Initialization and loop state

- Start with a short `commentary` update naming the target and first verification step.
- Use `update_plan` to track the target, iteration number, current stage, and open findings count. Update the plan at each stage boundary.
- Record findings in `compute/audit/linear_read_notes.md`. Each entry should include date, target, severity, class, exact file/line, issue, fix, and status.
- Unless the user explicitly asks for sub-agents or parallel agent work, keep the audit local. In Codex, parallel shell inspection happens through `multi_tool_use.parallel`; agent fan-out via `spawn_agent` is reserved for explicit user-authorized delegation.

### Per-claim verification protocol

For every nontrivial claim on the target surface:

1. Identify the exact statement, labels, hypotheses, and dependencies.
2. Check conventions first: grading, signs, bar/cobar objects, Swiss-cheese directionality, shifted PVA parity, and Vol I/Vol II terminology.
3. Recompute formulas independently. Do not correct by pattern or majority vote.
4. Trace proof steps line by line. Verify cited results really state what is being used and that their hypotheses are satisfied.
5. Check scope honestly: genus, arity, filtration, algebra family, open/closed colour, ambient versus convolution level.
6. Propagate across the repo: introduction, current chapter, examples, appendices, compute layer, superseded splits, and both volumes when the claim is cross-volume.
7. Run the narrowest relevant verification before declaring success: targeted `pytest`, targeted grep, log inspection, or `make fast`.

### The Six Hostile Examiners, Codex edition

Apply these lenses explicitly on dense passages:

- **Beilinson**: is the claim really proved, or merely asserted?
- **Witten**: does the mathematics match the physical claim, and is the scope qualified?
- **Costello**: are the factorization and Ran-space constructions actually well-formed?
- **Gaiotto**: do the VOA and W-algebra identifications survive edge cases and reductions?
- **Drinfeld**: are the quantum-group and transport statements structurally correct?
- **Kontsevich**: are the operadic/formality statements precise, functorial, and non-circular?

### Severity and classification

- Severity: `CRITICAL`, `SERIOUS`, `MODERATE`, `MINOR`
- Class:
  - `A` logical/circular
  - `B` formula/computation
  - `C` structural/label/dependency
  - `D` status/scope qualification
  - `E` editorial wording that affects truth conditions or mathematical readability

## Beilinson Rectification Loop — Codex Protocol

When instructed to run the Beilinson loop on a target chapter or on the live Vol II surface, iterate until convergence. Convergence means: no actionable findings at severity `MODERATE` or above remain on the modified surface, and the narrowest relevant verification passes.

### Stage 0 — Initialize

- Register the target in `update_plan`.
- Read the target file, its neighboring context, the active `\input` map from `main.tex`, and the current git diff before editing.
- Prefer the active split/core files over superseded files, but audit superseded files when they advertise the same claim.

### Stage 1 — Audit

Run three passes, locally by default:

- **RED pass**: logic, formulas, signs, scope, status tags, dependency honesty.
- **BLUE pass**: consistency, stale status prose, duplicate or conflicting formulations, AP5 propagation failures, build-log collisions, compute/test mismatches.
- **GREEN pass**: missing definitions, dangling forward references, unproved lemmas, over-conditional statements that can now be strengthened, and load-bearing structural gaps.

Codex execution rules:

- Use `multi_tool_use.parallel` for independent reads, greps, log checks, and targeted test invocations that can run concurrently.
- Use `exec_command` for sequential shell work, `git diff`, `git log`, `make fast`, and `pytest`.
- Only if the user explicitly authorizes delegation may RED/BLUE/GREEN be split across `spawn_agent`; otherwise they are three local passes in one session.

Merge the findings into a single register, deduplicate them, and classify them by severity and class.

### Stage 2 — Rectify

Fix findings in dependency order, starting with severity `CRITICAL` and `SERIOUS`, then `MODERATE`.

For each fix:

1. Re-read the local context.
2. Compute the correction independently.
3. Announce the intended edit in `commentary`.
4. Edit with `apply_patch`.
5. Grep for all variant forms across the active Vol II surface and, when relevant, `~/chiral-bar-cobar`.
6. Run the narrowest verification that can actually falsify the change.

Build discipline:

- After every three fixes, or after any load-bearing theorem/proof rewrite, run `make fast` unless a narrower check is clearly sufficient.
- If a structural rewrite needs isolation, create a temporary `git worktree` with `exec_command`, inspect the diff, and only then port the successful change back. Never overwrite unrelated user edits.

### Stage 3 — Re-audit

- Re-run RED/BLUE/GREEN on every modified location and on any nearby statements whose truth conditions changed.
- Treat Stage 3 as hostile verification of your own rewrite, not a rubber stamp.
- If any actionable finding at severity `MODERATE` or higher remains, start the next iteration with those findings as input.

### Finalize

- Run the narrowest defensible closing verification. For broad manuscript changes, that usually means `make fast`; for compute-facing claims, the relevant `pytest` slice as well.
- Update `compute/audit/linear_read_notes.md` with final statuses.
- Mark the `update_plan` entry complete only after verification, not after edits.

### Execution rules

- Never edit without reading first.
- Never trust repetition as evidence.
- Never upgrade a claim to `\ClaimStatusProvedHere` in the same pass as its first unchecked proof draft.
- After every substantive correction, search for downstream advertisements of the old claim and update them in the same session.
- Treat contradictions between manuscript, compute, tests, and build logs as mathematical bugs, not as "documentation drift."

## Vol I Theorems Used Here

Every chapter depends on Vol I's five theorems. Cross-references to Vol I labels resolve as "undefined" — expected for a multi-volume work.

| Vol I Theorem | What it supplies to Vol II |
|---------------|---------------------------|
| **(A)** Bar-cobar adjunction | The bar complex exists as a factorization coalgebra; reinterpreted as SC^{ch,top}-algebra structure — the bridge theorem |
| **(B)** Koszul inversion | Bar-cobar equivalence on the Koszul locus; lifted to raviolo VA setting and completed towers |
| **(C)** Complementarity | Genus-g obstructions decompose as complementary Lagrangians; the bulk-boundary-line triangle inherits this (−1)-shifted symplectic structure |
| **(D)** Leading coefficient | Curvature κ(A)·ω_g governs the genus tower; curved Swiss-cheese = Swiss-cheese + Hodge deformation |
| **(H)** Hochschild ring | BV-BRST origin of the deformation ring; bulk ≃ chiral Hochschild (Theorem H gets its physical explanation) |

## Nine Parts (Treatise Architecture)

**I. From the Bar Complex to the Swiss-Cheese Operad.** The bridge from Vol I to three dimensions. SC^{ch,top} constructed, recognition theorem, homotopy-Koszulity proved. Raviolo descent: SC-algebra → raviolo VA → PVA on cohomology. Axiomatics: sesquilinearity, A∞ relations, cluster factorization.

**II. The Descent Calculus.** Two descent mechanisms: cohomological (bar → PVA via Arnold/Stokes) and genus (curved bar over M̄_g). FM calculus and PVA descent.

**III. The Bulk-Boundary-Line Core.** Chiral Hochschild, brace algebra, bar-cobar review, line operators, spectral braiding (proved core), Koszul triangle (proved core), celestial boundary transfer (proved core), physical origins (merged chapter).

**IV. The Standard HT Landscape.** Worked examples: Rosetta stone, free multiplet, LG, CS (proved), W-algebras (stable framework).

**V. Quantization and Obstruction Theory.** Modular PVA quantization (core), affine half-space BV, planted-forest L∞, 3d gravity.

**VI. The Ordered/Open Sector and Transport.** Ordered associative chiral KD (core), dg-shifted factorization bridge. The E₁ wing as equal partner.

**VII. Holographic and Celestial Frontier.** YM synthesis (core), celestial holography (core), logarithmic HT monodromy (core), anomaly-completed holography (core).

**VIII. Extensions, Conditional Results, and Frontier.** All frontier/conjectural material from splits. No earlier part depends on this part.

**IX. Conclusion and Aftermatter.** Conclusion, appendices (brace signs, orientations, FM proofs, PVA expanded).

## Standing Hypotheses — REPLACED BY DEFINITION

The former standing hypotheses (H1)–(H4) have been replaced by a single definition, which axiomatises the class of algebras for which the full theory holds:

**Definition (Logarithmic SC^{ch,top}-algebra, Definition def:log-SC-algebra):** A C_*(W(SC^{ch,top}))-algebra whose closed-colour A∞ operations are defined by logarithmic weight forms factoring as ω_k = ω_k^hol ⊗ ω_k^top on FM_k(ℂ) × Conf_k(ℝ).

The former hypotheses now have precise logical status:

| | Content | Status |
|---|---------|--------|
| (H1) | BV data, one-loop finiteness | Condition of the bridge theorem (Theorem thm:physics-bridge) — applies only to physical realisations |
| (H2) | Propagator: meromorphic in ℂ, exponential decay in ℝ | Consequence of the bridge theorem for physical realisations |
| (H3) | FM compactification, logarithmic forms, AOS relations, Stokes exactness | Core content axiomatised in Definition def:log-SC-algebra; analytic consequences are Theorem thm:FM-calculus |
| (H4) | Factorization compatibility with C_*(W(SC^{ch,top})) | Recognition theorem (Theorem thm:recognition-SC) — already proved |

All results in Parts I–VII hold for any logarithmic SC^{ch,top}-algebra. Physical theories (gauge theories satisfying Theorem thm:physics-bridge) provide the standard class of examples. The word "derived" is not used for (H3): the definition axiomatises the logarithmic form data; Theorem thm:FM-calculus shows this implies the needed analytic properties.

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
- Chiral Koszulness from physics: **RESOLVED for the affine lineage** by the loop-order criterion (Theorem thm:one-loop-koszul), which shows one-loop exactness of the BV-BRST differential implies Koszulness. DS reduction (Theorem thm:ds-koszul-obstruction) explains the failure mechanism: Koszulness descends through DS iff the BRST cohomology inherits one-loop exactness. General case beyond the affine lineage remains open.
- **Homotopy-Koszulity of SC^{ch,top} is PROVED** (Theorem thm:homotopy-Koszul): via Kontsevich formality + transfer from classical Swiss-cheese (Livernet). ALL formerly conditional results (bar-cobar Quillen equivalence, filtered Koszul duality, C_line ≃ A!-mod, dg-shifted Yangian) are now unconditional.
- **Spectral Drinfeld strictification is COMPLETE** for all filtrations and all simple Lie algebras, via the complete strictification theorem. The former "filtration >= 4" frontier is fully resolved. The true remaining frontier is Kac-Moody algebras with root multiplicities > 1, where the strictification mechanism requires new input.
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
