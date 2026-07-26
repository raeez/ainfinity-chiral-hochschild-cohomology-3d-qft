# Beilinson Audit -- Vol II Preface (Wave 11)

Target: /Users/raeez/chiral-bar-cobar-vol2/chapters/frame/preface.tex
Lines: 700 (merged Wave 4-1 version)
Protocol: six-hostile-examiner deep Beilinson audit (read-only)
Scope: structural, mathematical, typographic, AP compliance

Severity tags: CRITICAL / HIGH / MODERATE / LOW / NIT.

## Executive summary

The preface is structurally strong. The Steinberg analogy (one coalgebra, three projections: cobar = inversion, Verdier = Koszul, Hochschild = bulk) is precise and survives examiner pressure. The three-leap architecture is well motivated. The AP-OC distinction is stated correctly. The slab-is-bimodule (RS-9) language is clean. All five Dimofte remark labels exist in ht_bulk_boundary_line_core.tex. The six integrated AP fixes (AP24, AP25/AP34, AP126, AP117, AP32, V2-AP26) are all present and correctly applied.

Two findings of MODERATE-to-HIGH severity survived the pass: a Virasoro "uncurved at c=13" mis-characterization that conflates uncurvedness with self-duality, and an AP125 prefix/environment mismatch on two Cardy/BTZ references written as "Theorems" but pointing at `conj:` labels. Several LOW/NIT items are logged below.

---

## Findings by examiner

### Formalist (proof structure, label hygiene, environment/tag discipline)

- **HIGH -- AP125 prefix mismatch (line 610--613).** Prose: "the two-sided Cardy formula for these MC states is conjectural (Theorems~\ref{conj:cardy-mc} and \ref{conj:gravity-btz-entropy})." Both labels resolve to `\begin{conjecture}` environments (3d_gravity.tex lines 3805 and 3855), tagged `\ClaimStatusConjectured`. Per AP125, the reference word must match the environment: these should be written as "Conjectures~\ref{conj:cardy-mc} and \ref{conj:gravity-btz-entropy}" -- writing "Theorems" here invites a future agent grepping for `conj:` (to find conjectures) to be misdirected. Self-consistent with the immediately adjacent "is conjectural" phrase, so factually honest, but the word "Theorems" violates AP125. Note: the same file has `conj:gravity-cardy` aliased to `conj:cardy-mc` two lines apart in 3d_gravity.tex (lines 3807--3808), which is also a separate AP124 risk but outside the preface scope.

- **HIGH -- unchecked universal quantifier on "Framework-level results are unconditional" (line 660).** The XIV signposts section groups thirteen-plus theorems and then writes "Framework-level results are unconditional; Vol~I Theorem~D and Theorem~C2 inputs carry uniform-weight scope at g >= 2." This is correct as stated, but the line above lists "complete spectral Drinfeld strictification (Theorem~\ref{thm:complete-strictification})" which per Vol II CLAUDE.md is PROVED for all simple Lie (frontier: Kac-Moody root mult > 1). The preface does not flag the Kac-Moody frontier caveat. MODERATE/HIGH depending on whether the body chapter carries the caveat at the theorem statement.

- **MODERATE -- "Theorem C2" is cited as a Vol I input (line 268--277).** The preface says "Theorem~D and the scalar part~C2 of Theorem~C are tagged uniform-weight." This is AP32-compliant tagging: every occurrence of obs_g, F_g, lambda_g downstream "inherits the explicit tag from its Vol~I source." Good. No violation here -- logging as a positive finding that AP32 was properly absorbed.

- **MODERATE -- "inversion (not Koszul duality, and the left adjoint of bar)" parenthetical (line 37).** Correct content, but the parenthetical has two clauses stacked: "(X, and Y)". Either drop the "and" or split into two parentheticals. Cosmetic; content is right.

### Topologist (operadic structure, geometric claims, fiber-base levels)

- **CRITICAL -- "The chiral dual Vir_{26-c} is uncurved at c = 13" (line 613--614).** This is a mathematical mis-characterization. Uncurvedness means kappa = 0. For Virasoro kappa(Vir_c) = c/2, so Vir_{26-c} is uncurved only at c = 26 (the matter-ghost critical dimension), NOT at c = 13. At c = 13, the strict same-family line-side comparison representative is fixed: Vir_{26-13} = Vir_{13}. Both Vir_{13} and its comparison representative have kappa = 13/2, both are curved. The preface conflates "uncurved" with "fixed point of the line-side comparison."
    - The very next clause ("...the structural antipode to the matter-ghost critical dimension c = 26") is the correct phrasing but for the WRONG structural fact: c = 13 is the self-dual point (half of 26), not the uncurved point.
    - AP8 explicitly says: NEVER "self-dual" unqualified; specify which duality. Here the author knows the duality (Koszul) and knows c = 13 is special, but wrote "uncurved" instead of "self-dual." This is CRITICAL because the body chapter (3d_gravity.tex) depends on the self-dual interpretation for the Page-time derivation.
    - The Page time argument t_P = 3 S_BH / 13 depending on kappa(Vir_c) + kappa(Vir_{26-c}) = 13 is CORRECT (verified algebraically: c/2 + (26-c)/2 = 13), and the body chapter (3d_gravity.tex lines 7488--7529) derives it self-consistently. The CRITICAL is strictly local to the word "uncurved."
    - **Recommendation (do not modify in this audit):** replace "uncurved at c = 13" with "self-dual at c = 13" in the next write-pass.

- **HIGH -- "Vir_{6k}" notation (line 600).** The preface writes "\Bbound = \mathrm{Vir}_{6k} at c = 6k." This is the central-charge-as-subscript convention used throughout Vol I and Vol II, and 3d_gravity.tex uses it consistently. The convention itself is fine, but at first appearance in the preface (where no prior convention is set) this reads as "Virasoro at level 6k" to a reader coming from Kac-Moody notation. Logging as HIGH because this is a preface and the convention should be either stated or use Vir_c with "c = 6k" explicit.

- **HIGH -- "the holomorphic direction is the Koszul differential of the Lagrangian embedding, the topological direction is the groupoid diagonal" (line 191--194).** Section IV calls the bar complex the "Koszul model of the Lagrangian self-intersection in a (-2)-shifted symplectic stack." The (-2)-shifted claim is precise, but the groupoid-diagonal description of the topological direction is thin. Is this proved in Theorem~\ref{thm:bar-is-self-intersection}? The preface cites the theorem but does not explain why groupoid diagonal ↔ deconcatenation coproduct. For a preface this is acceptable (body does the work), but a hostile reader will note the lack of definitional precision.

- **MODERATE -- AP130 fiber/base check on omega_g.** The preface writes d_fib^2 = kappa(A) * omega_g (lines 82--83, 158, 536, 688). omega_g is correctly the Hodge class on M̄_g (the base) while d_fib is the FIBERWISE differential. The equation d_fib^2 = kappa * omega_g lives on the TOTAL space of the bar complex fibered over M̄_g; kappa is a scalar (arithmetic) and omega_g pulls back to the total space. This is fine but the preface never says so explicitly. NIT/MODERATE: a one-line "the square lives on the total space, with omega_g pulled back from M̄_g" would close the gap.

- **MODERATE -- "the (-2)-shifted symplectic stack" (line 190).** Is this genuinely (-2)-shifted, or is the Lagrangian in a (-1)-shifted ambient with the self-intersection carrying the shift? PTVV convention: a Lagrangian in an n-shifted symplectic stack has an (n-1)-shifted symplectic structure on the derived self-intersection. So a Lagrangian in an (n=-1)-shifted ambient gives (-2)-shifted self-intersection -- consistent. But the preface says the stack itself is (-2)-shifted, which would put the Lagrangian self-intersection at (-3). Check thm:bar-is-self-intersection for the precise shift. Logging as MODERATE because a single-letter shift error is the kind of thing this audit catches.

### Physicist (physical claims, dimensional analysis, boundary/bulk structure)

- **LOW -- Brown-Henneaux formula.** c = 3ℓ/(2G) is the standard Brown-Henneaux central charge for AdS_3/CFT_2. The preface also identifies it with c = 6k where k = ℓ/(4G) is the CS level. Verification: 6 · ℓ/(4G) = 3ℓ/(2G). CORRECT.

- **LOW -- Heisenberg R-matrix R(z) = exp(k*hbar/z).** Matches V2-AP7. Correct with level prefix; vanishes at k = 0. AP126-clean.

- **MODERATE -- "Heisenberg class G, formal" (line 349).** The "formal" tag is correct (G is formal by definition), but a hostile reader will ask: "formal" in which sense -- SC-formal (m_k^SC = 0 for k >= 3), or dg-algebra formal (H*(A) quasi-iso to A)? AP14 explicitly distinguishes: "Koszulness != SC formality." In the Heisenberg case both hold, but the preface does not say which it means. For a preface this is LOW.

- **MODERATE -- "F_1 = -k log eta(tau)" for Heisenberg (line 361).** Verification: -k log eta = -k log(q^{1/24} prod(1 - q^n)) = -k/24 * log(q) + analytic corrections. As a numerical scalar at leading order, this gives F_1 = k/24 = kappa(H_k)/24, consistent with AP120's sanity-check identity F_1 = kappa/24. But "F_1" conventionally denotes a NUMBER (the first free energy), while "-k log eta(tau)" is a FUNCTION of tau. The preface is using "F_1" loosely to denote the whole genus-1 free energy. Either use a different symbol (F^(1)(tau)) for the function, or make explicit that F_1 here means the one-loop partition function. LOW/MODERATE.

- **MODERATE -- "the secondary anomaly u := eta^2" (line 543).** The letter eta is overloaded in the preface: (i) Dedekind eta(tau) in Heisenberg F_1 (line 361); (ii) the Arnold form eta_12 = d log(z_1 - z_2) (lines 72, 86); (iii) the transgression-algebra odd generator eta with d eta = Theta (line 541); and (iv) u = eta^2 as secondary anomaly (line 543). Four distinct uses in one preface. A careful reader can keep them apart by context, but a tired reader cannot. MODERATE.

- **LOW -- "infinite-depth class M" for Virasoro 3d gravity (line 601).** Correct: Virasoro is class M (infinite A-inf tower, d_alg = infinity). The preface appropriately does not conflate d_gen with d_alg here (no AP131 violation).

### Number Theorist (kappa identities, level arithmetic, modular data)

- **POSITIVE -- AP24 Virasoro Feigin-Frenkel (lines 382--386).** "kappa(Vir_c) + kappa(Vir_{26-c}) = 13, not zero: the two families share the duality name but not the involution data." This is AP24 stated correctly and FORCEFULLY. The comparison with affine Kac-Moody "kappa(g_k) + kappa(g_{-k-2h^v}) = 0" is the right pedagogical move. The preface even calls out the shared name / different involution data. PASS.

- **POSITIVE -- AP126 level-prefix on all r-matrices.** Heisenberg r_H(z) = k/z (line 355), Kac-Moody r(z) = k*Omega/z (line 373), "at k = 0 the r-matrix vanishes" sanity check stated explicitly (lines 346, 355, 376, 531). AP126/AP141 fully absorbed. PASS.

- **POSITIVE -- AP117 connection-vs-Arnold (lines 86--92, 373--377).** "The form eta_12 = d log(z_1 - z_2) is a bar-construction coefficient... it is not the connection form of a flat bundle. The KZ connection on conformal blocks... is the one-form Omega dz/(k + h^v), with a rational one-form in z, not a logarithmic one-form in z_1 - z_2." This is AP117 stated perfectly. PASS.

- **LOW -- kappa(Kac-Moody) = dim(g)(k + h^v)/(2 h^v) (line 379).** Matches CLAUDE.md constants. CORRECT.

- **LOW -- W_3 coefficient 32/(22 + 5c) (line 394--396).** Cross-verified against w-algebras-w3.tex line 51, w-algebras-conditional.tex line 185, w3_shadow_coefficients.py, and w3_shadow_closed_form.py. All write 32/(5c+22), which is the same polynomial. PASS.

- **LOW -- completion entropy h_K ~ 0.772 for W_3 (line 397--398).** Cross-verified against thqg_modular_pva_extensions.tex tables (lines 638, 1581) and standalone preface_full_survey.tex line 596. PASS.

- **LOW -- "quartic contact Q^{contact}_Vir = 10/[c(5c+22)]" (line 606--607).** Per AP129 this is exactly the reciprocal trap (previously written as -(5c+22)/(10c)). Verified here: at c = 1, denominator = 1*(27) = 27, so Q = 10/27 > 0. Preface writes the CORRECT form. PASS.

- **MODERATE -- kappa_eff = (c - 26)/2 "at one loop" (line 604).** The Vol I / Vol II convention: kappa_eff = kappa(matter) + kappa(ghost) vanishes at c = 26 (ghost cancellation). Matter kappa = c/2, ghost kappa = -13 (bc system). Sum = c/2 - 13 = (c - 26)/2. CORRECT. But the preface does not identify kappa(ghost) = -13 anywhere, and a hostile reader will want to see the origin. LOW/MODERATE.

### Adversarial Chef (rhetorical traps, AI slop, prose laws)

- **POSITIVE -- zero em dashes, zero AI slop.** Grepped the preface for "notably", "crucially", "remarkably", "interestingly", "furthermore", "moreover", "delve", "leverage", "tapestry", "cornerstone", "it is worth noting", "We now". ZERO hits. Em dashes: ZERO. Passive hedging: minimal. V2-AP29 compliant. PASS.

- **LOW -- "Unconditionally, Abulk identifies with chiral Hochschild cochains of the boundary algebra" (line 457--459, 508--509).** Stated TWICE verbatim (Section IX "Six direct faces" and Section XI "Bulk-boundary-line triangle"). Prose law "state once, prove once" gently violated. The repetition is defensible (each occurrence anchors a different construction), but a hostile Editor will flag it.

- **LOW -- "every chapter of this volume is a projection of the single coalgebra barB_ch(A) through one of the three" (line 45--47).** Repeated almost verbatim at line 664--666: "every chapter is a different projection of the single coalgebra barB_ch(A)." This is the explicit preface thesis, and the repetition IS intentional (it frames Section I and closes Section XIV). Consider it rhetorical bookending, not violation. NIT.

- **LOW -- "The climax of the volume" (line 597).** Meta-commentary. The Vol II CLAUDE.md identity says gravity IS the climax (Part VI, per RS-13). The prose laws warn against "This chapter..." meta-openings (AP106, AP108, AP111). "The climax" is a mild version of this. The rest of Part VI description is concrete, so this is borderline. LOW.

- **MODERATE -- Section V mixes AP-OC, slab, and interface Dimofte content in one section.** The section title is "The open primitive and the AP-OC distinction" and contains three distinct sub-topics: (a) the categorical open sector (lines 201--212); (b) the three independent objects under AP-OC (lines 214--230); (c) the slab fiber functor and five Dimofte cross-refs (lines 231--250); (d) Morita invariance and modularity (lines 253--259). Each is a paragraph or less. A hostile Editor would split this into two sections (AP-OC proper, and Slab) or at least bold-subtitle each paragraph. But all content is genuine and tightly worded. MODERATE.

- **POSITIVE -- Six direct faces table (lines 423--448) is clean and dense.** Exactly the CG "decomposition table" structural move. PASS.

- **LOW -- "The two-structure thesis of this volume" (line 120).** Introduces a name ("two-structure thesis") that does not appear to be used anywhere else in the preface or body. Logging as minor terminological orphan. NIT.

### Editor (cross-references, consistency, Part numbers, label sync)

- **POSITIVE -- V2-AP26 compliance.** All Part references use \ref{part:...}: part:swiss-cheese (line 565), part:e1-core (line 569), part:bbl-core (line 575), part:examples (lines 584, 697), part:holography (line 588), part:gravity (line 596), part:frontier (line 620). ZERO hardcoded Part numbers. V2-AP26 fully absorbed. PASS.

- **POSITIVE -- all five Dimofte cross-references have targets.** Grepped for rem:slab-fiber-functor, rem:dimofte-interface-generalization, rem:dimofte-k-matrix, rem:dimofte-double-bosonization, rem:dimofte-meromorphic-braided across vol2. All five appear in ht_bulk_boundary_line_core.tex (plus the frontier + spectral-braiding-core partners). PASS.

- **POSITIVE -- key theorem labels resolve.** Spot-checked thm:rosetta-swiss-cheese, thm:homotopy-Koszul, thm:cohomology_PVA, thm:YBE, thm:boundary-linear-bulk-boundary, thm:lines_as_modules, thm:3d-universal-mc, thm:affine-monodromy-identification, thm:bar-is-self-intersection, thm:holographic-reconstruction, thm:complete-strictification, thm:general-half-space-bv, thm:doubling-rwi, thm:bulk-boundary-line-factorization, thm:modular-bar, thm:formal-genus-expansion, thm:synthesis, thm:resolvent-principle, prop:sc-koszul-duality-nonselfdual. All resolve via grep across chapters/. PASS.

- **HIGH -- V1-thm:shadow-formality-identification (line 285).** Prefixed "V1-" suggests a cross-volume reference to Vol I. Does Vol II's build phantomsection-stub this label, or does it pull from Vol I's label space? Per AP127, cross-refs to migrated chapters must have \phantomsection\label{} stubs. Out-of-scope to verify build-side here, but flagging for the next rectification pass: does \ref{V1-thm:shadow-formality-identification} resolve in the Vol II build?

- **MODERATE -- "scoped identifications" phrasing (lines 449--456).** "The scoped identifications C_line ≃ A^!_line-mod and Abulk ≃ Z_der(Bbound) ≃ HH^bullet_ch(A^!_line) assemble these faces into the bulk-boundary-line triangle (Theorems~\ref{thm:boundary-linear-bulk-boundary} and \ref{thm:lines_as_modules}) on the boundary-linear exact sector and on the chirally Koszul locus respectively." The phrase "boundary-linear exact sector... chirally Koszul locus respectively" pairs the two conditions with the two identifications -- but which scope goes with which? The reader has to match by order. MODERATE: rewrite as "... on the boundary-linear exact sector (for Abulk) and on the chirally Koszul locus (for C_line)" would disambiguate. The content is right, the pairing is unclear.

- **LOW -- "Definition~\ref{def:btz-deformation}" (line 610).** Verified the label exists in 3d_gravity.tex. PASS.

- **LOW -- Steinberg section numbering.** The preface uses I -- XV section headers, all \section* (unnumbered). Consistent. PASS.

---

## AP-fix integration audit (the six claimed Wave 4-1 fixes)

| AP | Where | Status |
|----|-------|--------|
| AP24 (Vir FF anti-symmetry = 13, not 0) | lines 382--386 | INTEGRATED, forceful |
| AP25 / AP34 (dbar sign; divided-power lambda-bracket c/12 lambda^3) | line 394--395 (T-T lambda bracket (c/12) lambda^3) | INTEGRATED |
| AP126 (level-prefixed r-matrix; k=0 vanishing) | lines 346, 355, 373, 376, 530--531 | INTEGRATED, five explicit vanishing checks |
| AP117 (Arnold eta is bar coefficient, not connection; KZ = Omega dz/(k+h^v)) | lines 86--92, 373--377 | INTEGRATED, correctly separated |
| AP32 (obs_g / F_g / lambda_g carry explicit uniform-weight tag) | lines 268--277, 537--539, 660--663 | INTEGRATED, all three occurrences tagged |
| V2-AP26 (no hardcoded Part numbers, use \ref{part:...}) | lines 565, 569, 575, 584, 588, 596, 620, 697 | INTEGRATED, zero hardcoded Part numbers |

All six claimed fixes are present and correct. PASS.

---

## Steinberg analogy -- depth check

The Steinberg analogy in Section I claims: the Steinberg variety is one geometric object with THREE distinct extractions (BM-homology = Hecke algebra, K-theory = Springer correspondence, equivariant derived category = KL basis), and the chiral bar complex barB_ch(A) is ONE coalgebra with THREE distinct projections (cobar = inversion = left adjoint; Verdier on Ran = Koszul duality; Hochschild cochains = derived centre = bulk).

Is the analogy mathematically precise?

- **The three Steinberg extractions ARE genuinely distinct** (BM != K != D^b_eq) and produce genuinely different algebras. PASS.
- **The three bar projections ARE genuinely distinct** (cobar is a left adjoint; Verdier is a duality functor on Ran; Hochschild is a homotopical center computation). Omega B(A) ≃ A is inversion (Vol I architecture). D_Ran(barB(A)) ≃ barB(A^!) is Verdier. C^bullet_ch(A,A) ≃ Z^der_ch(A) is Hochschild. These are three genuinely distinct functors from coalgebras (or from the Ran-equivariant category) and the preface correctly says "no two of them agree." PASS.
- **The analogy is SUGGESTIVE but NOT a precise parallel.** In the Steinberg case, the three extractions are cohomology theories applied to the SAME geometric object. In the bar case, the three projections are three different CONSTRUCTIONS on the same homotopical object, and only Hochschild is naturally described as "cohomology." Cobar is a functor in a different direction; Verdier is a duality involution. The analogy is a rhetorical frame, not a functorial equivalence.
- The preface does not overclaim. It says "the Steinberg analogy is not decoration: the three projections play three genuinely different roles, and never coincide." This is accurate. The analogy is suggestive-plus, not suggestive-only. PASS with a note that a hostile reader will remark on the category-theoretic asymmetry.

---

## Three conceptual leaps -- precision check

- **(1) Com -> E_1/A_inf.** "Classical Koszul duality works for commutative/associative algebras. The chiral setting replaces Com by the Swiss-cheese operad SC^{ch,top}." Precise. The Stasheff hierarchy on FM_k(C) is correctly identified. The "deconcatenation coproduct on Conf_k^<(R)" is correctly identified. PASS.

- **(2) Genus 0 -> modular.** "Classical Koszul duality is a genus-0 computation: one differential, one cohomology. The modular setting lifts the bar complex over M̄_g,n and tracks the curved differential d_fib^2 = kappa(A) * omega_g." Precise. AP32 tag inherited from Section VI. PASS.

- **(3) Local -> nonlocal.** "The E_inf layer (locality, Sigma_n-equivariance, OPE with arbitrary pole order) sits inside a genuinely E_1 layer. All standard vertex algebras, including those with poles, are E_inf: locality is compatible with OPE poles. The genuinely E_1 layer is where the spectral R-matrix is an INDEPENDENT input (not derived from the local OPE) and the line-side open-colour dual A^!_line is a dg-shifted Yangian." This is V2-AP1/V2-AP12 stated perfectly. "The distinction is provenance, not pole order" is V2-AP2 stated perfectly. PASS.

All three leaps are stated precisely and correctly. The leaps act independently and together produce the seven-part architecture. PASS.

---

## AP-OC distinction -- check

Section V states: "the bar complex classifies twisting morphisms, the bulk is the derived centre, and the two live on different levels of the categorical hierarchy. The curved genus-g extension deforms the bar differential by kappa(A) * omega_g without disturbing the coproduct: the deformation is purely closed-colour."

- AP-OC stated PRECISELY: bar = twisting morphism classifier, bulk = derived centre. The hierarchy (categorical levels) is correctly identified.
- "Curved genus-g extension does not disturb the coproduct" is the claim that the deformation is closed-colour only. Verify: d_fib^2 = kappa * omega_g is a curvature of the DIFFERENTIAL; the deconcatenation coproduct Delta is unchanged. This is consistent with eq preface2-coderivation (Delta circ D = (D x id + id x D) circ Delta): adding a closed-colour curvature to D leaves this identity intact modulo the standard curved-coderivation signs. PASS.

---

## Slab = bimodule (RS-9) -- check

Line 233--236: "The SC^{ch,top} open disk admits a concrete geometric realisation as a slab with transverse boundary conditions: a bimodule geometry with two boundary components, NOT a Swiss-cheese disk."

- RS-9 ("The slab is a bimodule, NOT a Swiss-cheese disk. Two boundary components") stated EXPLICITLY with the capitalized NOT. PASS.
- The subsequent sentence identifies the slab universal algebra as the Drinfeld double of the two boundary algebras. This matches the Dimofte integration story in MEMORY.md.
- Vol I Theorem C (complementarity) is called "the algebraic manifestation of transversality" -- this connection is a new framing and is stated crisply.

PASS.

---

## Summary of findings

| Severity | Count | Summary |
|----------|-------|---------|
| CRITICAL | 1 | "Vir_{26-c} uncurved at c=13" (should be "self-dual at c=13") |
| HIGH | 4 | AP125 "Theorems" vs `conj:` labels; "Vir_{6k}" convention at first use; Lagrangian-diagonal thinness; V1-thm prefix label build check |
| MODERATE | 9 | parenthetical stacking; shift-grading check on (-2); fiber/base statement on d_fib; heisenberg "formal" ambiguity; F_1 scalar vs function; eta overload (four senses); "scoped identifications" pairing clarity; Section V density; kappa_eff = (c-26)/2 origin |
| LOW | 8 | Brown-Henneaux, Heisenberg R, kappa-KM formula, W_3 coefficient, h_K 0.772, quartic contact reciprocal, "state once" repetition, "climax" meta, BTZ label |
| NIT | 2 | "two-structure thesis" orphan; Steinberg thesis bookending |
| POSITIVE | 12 | V2-AP26 Part refs; zero AI slop; AP126 level prefix; AP117 Arnold vs KZ; AP24 forceful; AP32 tags; all six Wave 4-1 fixes absorbed; Steinberg analogy precise; three leaps precise; AP-OC precise; RS-9 slab-is-bimodule explicit; decomposition table |

## Recommendations for next rectification pass

1. **CRITICAL fix**: line 613--614, replace "is uncurved at c = 13" with "the same-family line-side comparison is fixed at c = 13". Do NOT touch the Page time argument which is correct.
2. **HIGH fix**: line 610--613, replace "Theorems" with "Conjectures" to match `conj:` labels (AP125).
3. **HIGH fix**: line 600, either introduce the Vir_c convention at first use or write "Vir_c at c = 6k" instead of "Vir_{6k}".
4. **HIGH check**: verify \ref{V1-thm:shadow-formality-identification} resolves in Vol II build (AP127 phantom stub).
5. **MODERATE fixes**: disambiguate "scoped identifications" pairing; clarify F_1 symbol (scalar vs function); consider splitting Section V into AP-OC proper and Slab subsections.
6. **LOW fixes**: eta overloading footnote; "state once" dedup between Section IX and Section XI for the Hochschild = Abulk identification.

The preface is close to publication-ready. The CRITICAL uncurved/self-dual conflation is the one finding that MUST be fixed before any further build; every other finding is a polish pass.

## Files referenced

- /Users/raeez/chiral-bar-cobar-vol2/chapters/frame/preface.tex (target)
- /Users/raeez/chiral-bar-cobar-vol2/CLAUDE.md (V2-AP constraints)
- /Users/raeez/chiral-bar-cobar/CLAUDE.md (canonical AP index)
- /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/3d_gravity.tex (Cardy conjecture labels, Page time derivation, BTZ def)
- /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/ht_bulk_boundary_line_core.tex (all five Dimofte rem: labels)
- /Users/raeez/chiral-bar-cobar-vol2/chapters/examples/w-algebras-w3.tex (32/(5c+22) cross-check)
- /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/thqg_modular_pva_extensions.tex (h_K 0.772 cross-check)
- /Users/raeez/chiral-bar-cobar-vol2/compute/w3_shadow_coefficients.py (W_3 beta cross-check)
- /Users/raeez/chiral-bar-cobar/chapters/examples/w_algebras_deep.tex (KM Feigin-Frenkel cross-check)

End of audit.
