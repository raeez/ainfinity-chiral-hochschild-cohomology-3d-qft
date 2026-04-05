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
- The **pole-order dichotomy**: double poles → class L (formal SC structure), quartic poles → class M (non-formal, genuinely infinite A∞). DS reduction transports L→M via Sugawara: the affine double-pole OPE becomes the Virasoro quartic-pole OPE, escalating shadow depth from finite to infinite.
- The **self-dual point c* ≠ critical string point c_crit**. For W_N: c* = α_N/2 (Koszul self-duality), c_crit = α_N (matter-ghost cancellation: c = |c_ghost| = α_N). For Virasoro: c*=13, c_crit=26. For W₃: c*=50, c_crit=100. NEVER conflate these. The formula c_crit = α_N/(2(H_N−1)) is WRONG for N≥3; the correct ghost curvature is κ_ghost = −α_N·(H_N−1), not −α_N/2.

## Anti-Patterns (inherited from Vol I + Vol II-specific)

All Vol I anti-patterns AP1-AP34 apply here. In particular, AP34 (bar-cobar inversion ≠ open-to-closed passage) is especially load-bearing in Vol II: the three functors on B(A) are (1) Ω(B(A)) ≅ A (reconstruction), (2) Ω(D_Ran(B(A))) ≅ A! (Koszul duality), (3) RHom(Ω(B(A)), A) = C^•_ch(A,A) (derived centre = universal bulk). NEVER write "bar-cobar inversion produces the bulk." The following encode deep mathematical errors found propagated across both volumes:

**AP-OC — Bar = bulk conflation.** The bar complex B(A) classifies TWISTING MORPHISMS (universal couplings between A and A!). The BULK OBSERVABLES are a DIFFERENT object: the chiral derived center Z^der_ch(A) = C^•_ch(A_b, A_b). The PRIMITIVE OBJECT is the open-sector factorization dg-category C_op; the boundary algebra A_b = End(b) is a CHART (Morita-invariant). Modularity belongs to TRACE + CLUTCHING on the open sector, not to the closed algebra alone. Found at: ht_physical_origins.tex lines 359-361 (corrected), celestial_holography_frontier.tex line 846 (corrected), thqg_bv_ht_extensions.tex lines 218/302/704 (corrected), plus 4 superseded files.

**AP19 — The bar kernel absorbs a pole.** The collision residue r(z) has pole orders ONE LESS than the OPE — the d log(z-w) kernel absorbs one power. Virasoro r-matrix: (c/2)/z³ + 2T/z, NOT (c/2)/z⁴ + 2T/z² + ∂T/z.

**AP20 — An invariant of one algebra is not an invariant of a system.** κ(A) is intrinsic; κ_eff = κ(matter) + κ(ghost) is composite; κ(B) where B = A! is the dual's. F_g always uses κ(A). State WHICH algebra's κ.

**AP21 — A class is not a scalar; Clifford ≠ exterior.** u = η² = λ = κ(B)·ω_g is LINEAR in κ (a class), NOT κ² (a scalar). The gravity dichotomy (c≠26 vs c=26) is whether this class vanishes. Squaring κ destroys the bifurcation.

**AP22 — Generating function index mismatch.** Â(iℏ)-1 starts at ℏ². If F_1 ≠ 0, the pairing Σ F_g ℏ^{power} must use ℏ^{2g} (not ℏ^{2g-2}) unless an explicit 1/ℏ² appears.

**AP23 — Flat section vs weighted transport.** √(Q_L) is the flat section of ∇^sh. The shadow generating function H(t) = t²√(Q_L) is NOT horizontal — the t² is the arity offset.

**AP24 — The complementarity sum is not universally zero.** κ + κ! = 0 for KM/free fields/lattice/principal W. For Virasoro: κ + κ! = 13 ≠ 0. The anti-symmetry was overclaimed in 20+ locations. The Feigin-Frenkel involution ensures anti-symmetry for KM; the Virasoro involution c ↦ 26−c is anti-symmetric around 13, not 0.

**AP25 — Three functors, three outputs: bar ≠ Verdier dual ≠ cobar.** B(A) = coalgebra. D_Ran(B(A)) ≃ B(A!) = Verdier dual (an ALGEBRA). Ω(B(A)) ≃ A = cobar (recovers the ORIGINAL). The Koszul dual A! is obtained by Verdier duality, NOT by cobar. Found conflated in 16 files.

**AP26 — Fock inner product ≠ BPZ inner product.** At weight ≥ 4 for rank ≥ 3 W-algebras, dim(Fock) > dim(W-algebra) and ⟨W₄|Λ⟩_Fock ≠ 0 even though ⟨W₄|Λ⟩_BPZ = 0. W-algebra decompositions MUST use BPZ (Wick contractions), not the free-field dot product.

**AP33 — Koszul duality ≠ Feigin-Frenkel duality ≠ negative-level substitution.** Three operations share surface similarities but differ. Koszul: A ↦ A^! = (H*(B(A)))^v. FF involution: k ↦ -k-2h^v within same family. Negative-level: H_k ↦ H_{-k}. For Heisenberg: κ(H_k^!) = -k = κ(H_{-k}), but H_k^! ≠ H_{-k} as chiral algebras. H_k^! = Sym^ch(V*). NEVER write H_k^! = H_{-k}.

**AP-CHR — Pole order ≠ chromatic height; the classical theory is height 0.** The holomorphic weight filtration F^p on SC^{ch,top} classifies algebras by algebraic depth (shadow depth r_max: classes G/L/C/M). This is NOT chromatic height. The classical bar complex B(A) is a chain complex over ℂ, hence an Hℂ-module, hence rational: L_{K(n)}(B(A)) = 0 for all n ≥ 1. The ENTIRE classical bar-cobar theory is chromatic height 0. No OPE pole produces positive-height chromatic data in the linear theory. The spectral-parameter composition law is the additive formal group Ĝ_a (height ∞ over F_p, not height 1 — do NOT confuse with Ĝ_m). Chromatic height ≥ 1 enters only after spectralization (chain complexes → spectra); height ≥ 2 enters at genus ≥ 1 through the formal group of E_τ. Conf_k(ℝ) is discrete (k! contractible components) and chromatically inert. A quartic pole does NOT mean "height 2."

**AP-RED — Koszul duality does not raise chromatic height; clutching raises p-div height by 2.** Koszul duality A ↦ A^! reflects the curvature (κ ↦ −κ for KM, κ ↦ 13−κ for Virasoro) but does NOT raise chromatic height. The genus-raising clutching D₁ adds one handle, increasing dim(J) by 1 and the p-divisible group height ht(J[p^∞]) by 2 (since ht = 2·dim for abelian varieties). NEVER write "one unit of p-divisible group height per handle" — it is two. The FF family {V_k(g)} interpolates between chromatic regimes as k varies, but formal-group height is a p-local invariant — it requires reduction mod p, not a complex parameter. NEVER write "Koszul duality is chromatic redshift." NEVER assign formal-group heights to chiral algebras over ℂ without specifying a prime.

**AP-TOWER — The chromatic tower uses E(n), not K(n).** The chromatic convergence theorem (Hopkins-Ravenel) uses the tower ... → L_{E(n)}X → L_{E(n-1)}X → ... where L_{E(n)} = L_{E(n)} is Johnson-Wilson localization. The monochromatic K(n)-localizations L_{K(n)} do NOT form a tower with natural transition maps (they live in different categories). NEVER write a "chromatic tower" indexed by L_{K(n)} → L_{K(n-1)}.

**AP-WICK — The modular S-transform is NOT the Wick rotation of the topological direction.** The S-transform τ → −1/τ exchanges the A-cycle and B-cycle of E_τ: it acts on the CLOSED colour (complex structure of the worldsheet). The Wick rotation of ℝ_t acts on the OPEN colour (the E₁ ordering). They are sourced by the same B-cycle quasi-periodicity defect but couple to DIFFERENT algebraic data: the S-transform couples to κ(A) (highest-pole coefficient), while the braiding monodromy couples to c₀ = {a_{(0)}b} (lowest-pole coefficient). Heisenberg proves the distinction: κ = k ≠ 0 (S nontrivial) while c₀ = 0 (braiding trivial). The two colours can entangle at genus ≥ 1 when c₀ ≠ 0 (Theorem thm:elliptic-spectral-dichotomy). The identification of S with "Wick rotation" is a 2d worldsheet perspective that does NOT survive the 3d HT lift.

**AP35 — E_∞-chiral includes ALL vertex algebras, not just pole-free ones.** Every BD chiral algebra (vertex algebra) — including Kac-Moody V_k(𝔤), Virasoro Vir_c, Heisenberg H_k — is E_∞-chiral. The E_∞ condition is LOCALITY: the chiral product is defined on the UNORDERED configuration space Conf_n(X) with Σ_n-equivariant factorization structure. OPE poles are compatible with E_∞ — they are singularities of a LOCAL product, not a failure of commutativity. "E₁-chiral algebra" means NONLOCAL: the factorization structure is defined on ORDERED configurations as primitive data, with no underlying symmetric structure. The standard class of genuinely E₁ objects is Etingof-Kazhdan quantum vertex algebras, Yangians, and quantum groups. NEVER say "E_∞ means no OPE poles." NEVER say "Heisenberg/Virasoro/Kac-Moody is not E_∞." NEVER contrast E₁ and E_∞ as "poles vs no poles" — the contrast is NONLOCAL vs LOCAL.

**AP36 — R(z) ≠ τ does NOT imply genuinely E₁; the discriminant is provenance, not value.** The R-matrix R(z) on the ordered bar complex can be nontrivial for BOTH E_∞ and E₁ algebras. For E_∞-chiral algebras (all vertex algebras), R(z) is DERIVED from the local OPE via analytic continuation — it is the monodromy of the flat connection on configurations, determined by the OPE pole residues. For genuinely E₁-chiral algebras (quantum vertex algebras), R(z) is INDEPENDENT INPUT — part of the defining data, not derivable from a local symmetric structure. The Heisenberg H_k has R(z) ≠ τ (nontrivial spectral dependence from OPE poles), but it IS E_∞ because R(z) is determined by the local OPE J(z)J(w) ~ k/(z-w)². The condition R(z) = τ (Koszul-signed flip) characterizes the pole-free subclass (commutative chiral algebras with no OPE singularities), which is a STRICT SUBCLASS of E_∞-chiral, not all of it. NEVER write "E_∞ implies R(z) = τ." The correct three-tier picture: (i) pole-free commutative: R(z) = τ, (ii) vertex algebras with poles: R(z) ≠ τ but derived from locality, (iii) genuinely E₁: R(z) ≠ τ and independent input. Both (i) and (ii) are E_∞.

**AP37 — Three bar complexes, three objects: FG bar ≠ full symmetric bar ≠ ordered bar.** (a) The FG bar B^{FG}(A) (Francis-Gaitsgory) uses ONLY the zeroth product a_{(0)}b (the chiral Lie bracket), ignoring all higher OPE poles. It is the bar complex of A viewed as a chiral Lie algebra. (b) The full symmetric bar B^{Σ}(A) uses ALL OPE products {a_{(n)}b}_{n≥0} and takes Σ_n-coinvariants — this is the bar complex of Vol I Theorem A. (c) The ordered bar B^{ord}(A) uses ALL OPE products but retains the linear ordering — no Σ_n quotient. These are THREE DIFFERENT chain complexes producing three different Koszul duals: B^{FG} gives the chiral Lie Koszul dual, B^{Σ} gives the full chiral Koszul dual A^! of Vol I, B^{ord} gives the ordered (associative) chiral Koszul dual of Part VII. There is a natural map B^{ord} → B^{Σ} (take coinvariants) and a filtration on B^{Σ} whose associated graded recovers B^{FG} (retain only zeroth poles). NEVER conflate these. In particular, B^{FG}(A) does NOT see higher poles and therefore misses the curvature κ(A).

**AP38 — Ordered-to-unordered descent is R-matrix twisted, not naive, even for E_∞.** For E_∞-chiral algebras with OPE poles (all interesting vertex algebras), the ordered bar B^{ord}(A) and the symmetric bar B^{Σ}(A) are DIFFERENT chain complexes: one lives on ordered configurations, the other on Σ_n-coinvariants. The descent B^{ord} → B^{Σ} requires the R-matrix as twisting datum: B^{Σ}_n ≃ (B^{ord}_n)^{R-Σ_n} (Proposition prop:r-matrix-descent). For pole-free commutative algebras, R(z) = τ and the descent is the naive Σ_n-quotient. For vertex algebras with poles (Heisenberg, Kac-Moody, Virasoro), R(z) ≠ τ and the descent is genuinely twisted — this is nontrivial mathematics even though both sides are E_∞. The ordered bar complex of an E_∞-chiral algebra with poles carries STRICTLY MORE INFORMATION than the symmetric bar: the R-matrix monodromy is the surplus datum. NEVER write "ordered = unordered for E_∞" without specifying whether you mean pole-free E_∞ (true, R = τ) or E_∞ with poles (false without R-twist). The manuscript's three-tier picture at ordered_associative_chiral_kd_core.tex lines 2204-2214 is correct and should be preserved.

**AP39 — NEVER equate "E_∞-chiral" with "no OPE poles" or "commutative chiral algebra in the BD sense."** BD's "commutative chiral algebra" (= no poles, = D-scheme) is a STRICT SUBCLASS of E_∞-chiral, not a synonym. E_∞-chiral = LOCAL = Σ_n-equivariant factorization = ALL vertex algebras. The BD commutative subclass is the pole-free stratum within E_∞. The false identification "E_∞ = no poles" was introduced into the manuscript during the 2 April 2026 editing session and caused cascading errors. Specific false glosses that must NEVER appear: "E_∞-chiral algebras (commutative chiral algebras in the sense of BD — those whose chiral product has no OPE singularities)" — this parenthetical restricts E_∞ to the pole-free subclass. The correct parenthetical is the original: "E_∞-chiral algebras (commutative chiral algebras)" — where "commutative" means LOCAL/Σ_n-equivariant, not pole-free.

**AP40 — BD do not study E₁-chiral algebras. E₁-chiral is a NEW concept.** All of Beilinson-Drinfeld's chiral algebras are E_∞-chiral (local, Σ_n-equivariant). The concept of "E₁-chiral algebra" (nonlocal, ordered configurations, no Σ_n symmetry) is introduced in THIS manuscript (Vol II, Part VII) as a genuinely new algebraic structure. It is NOT a BD concept. NEVER say "BD chiral algebras include E₁" or "some BD chiral algebras are E₁." NEVER say "E₁-chiral algebra is a standard notion." It is the manuscript's contribution, motivated by the Etingof-Kazhdan theory of quantum vertex algebras.

**AP41 — The Heisenberg R-matrix is R(z) = exp(kℏ/z), NOT trivial.** The Heisenberg OPE J(z)J(w) ~ k/(z-w)² has a double pole. By AP19, the d log kernel absorbs one power: the COLLISION RESIDUE is r^{coll}(z) = k/z (a simple pole in the collision residue, from a double pole in the OPE). The connection on Conf₂^ord(ℂ) is ∇ = d - k·d log(z), and the monodromy is exp(-2πik) — a nontrivial scalar phase. The spectral R-matrix on modules is R(z) = exp(kℏ/z) (confirmed in rosetta_stone.tex line 1102). A previous version of this anti-pattern (and multiple agent computations on 2 April 2026) INCORRECTLY claimed R = 1 for Heisenberg by integrating the LAPLACE kernel k/z² instead of the COLLISION RESIDUE k/z. This error arose from forgetting AP19: the bar complex uses d log propagators, not dz propagators. NEVER compute R-matrix monodromy from the raw OPE kernel — always apply the d log absorption first.

**AP42 — NEVER introduce parenthetical glosses that RESTRICT a defined term to a subclass.** The specific error pattern: the manuscript defines "E_∞-chiral algebras" correctly (= all vertex algebras = local). An LLM assistant then adds a parenthetical like "(commutative chiral algebras in the sense of BD — those whose chiral product has no OPE singularities)" which RESTRICTS the meaning to the pole-free subclass. This is a metamathematical error: the parenthetical narrows a broader term to a special case, causing all subsequent reasoning to apply only to the special case while thinking it applies generally. When adding parenthetical clarifications, ALWAYS check: does the parenthetical PRESERVE the full scope of the term, or does it RESTRICT it? If it restricts, it is an error.

**AP43 — NEVER say "ĝ_k is not E_∞-chiral" or "Virasoro is not E_∞-chiral" or "Heisenberg is not E_∞-chiral."** All three ARE E_∞-chiral. They are local vertex algebras with Σ_n-equivariant operations. Having OPE poles does not break E_∞. This error was committed multiple times on 2 April 2026 and caused cascading damage.

**AP44 — NEVER say "E_∞ implies R(z) = τ" without the pole-free qualifier.** The correct statement: "For pole-free E_∞-chiral algebras, R(z) = τ." For E_∞-chiral algebras WITH OPE poles (all interesting vertex algebras), R(z) ≠ τ. The E_∞ condition constrains R(z) to be DERIVABLE from the local OPE, not to be trivial.

**AP45 — NEVER conflate "E_∞-chiral" with "BD commutative chiral algebra."** BD's commutative chiral algebras (Chapter 4 of Chiral Algebras) are the pole-free subclass: chiral algebras whose bracket extends across the diagonal without poles. This is STRICTER than E_∞-chiral. All BD commutative chiral algebras are E_∞, but not all E_∞-chiral algebras are BD-commutative. The gap is precisely the vertex algebras with poles (ĝ_k, H_k, Vir_c, etc.).

**AP46 — NEVER say "the distinction between E₁ and E_∞ is about poles."** The distinction is about LOCALITY. E_∞ = local (the chiral operations are defined on unordered configuration spaces). E₁ = nonlocal (the operations are defined on ordered configurations only, with no underlying symmetric equivariance). A vertex algebra with poles (like ĝ_k) is E_∞ because it is LOCAL — the poles do not break locality. A quantum vertex algebra (Etingof-Kazhdan) is E₁ because it is NONLOCAL — the braiding is independent input, not derived from a local OPE.

**AP47 — NEVER treat the literature agent's claim "vertex algebras are NOT E_∞" as authoritative.** On 2 April 2026, an LLM agent consulted the BD literature and concluded "A vertex algebra is NOT E_∞-chiral — it has OPE poles." This conclusion was WRONG. The agent confused BD's "commutative chiral algebra" (= no poles, a subclass of E_∞) with E_∞-chiral itself (= all local chiral algebras). The BD book studies E_∞-chiral algebras throughout; their "commutative" subclass is the special case where the Lie*-bracket vanishes. ALL BD chiral algebras are E_∞. BD do not study E₁-chiral algebras at all.

**AP48 — NEVER "fix" pre-existing manuscript text by adding restrictive parentheticals.** The specific destructive pattern: the manuscript says "E_∞-chiral algebras (commutative chiral algebras)" and an LLM adds "in the sense of BD — those whose chiral product has no OPE singularities." This NARROWS the parenthetical from the correct meaning (commutative = local = Σ_n-equivariant) to the incorrect meaning (commutative = no poles). The pre-existing parenthetical was correct. The "fix" broke it. Before adding any parenthetical gloss to a mathematical term, verify that the gloss PRESERVES the original scope.

**AP49 — NEVER oscillate between conventions within a single session.** On 2 April 2026, an LLM assistant changed its interpretation of "E_∞-chiral" FOUR TIMES in a single session: (1) initially treated vertex algebras as E_∞ (correct), (2) then concluded vertex algebras are NOT E_∞ because they have poles (wrong), (3) then reverted to "E_∞ = pole-free" based on a literature agent (wrong), (4) then accepted the author's correction that vertex algebras ARE E_∞ (correct). Each oscillation introduced new errors into the manuscript. The rule: once a convention is established by the AUTHOR, lock it in and do not revisit based on agent outputs or literature searches. The author's word is final on conventions.

**AP50 — NEVER edit E₁/E_∞ language without explicit author confirmation.** This anti-pattern exists because every autonomous attempt to "clarify" or "correct" E₁/E_∞ terminology on 2 April 2026 made things worse. The error surface is too large and the conventions too subtle for autonomous editing. Before touching ANY sentence containing E₁, E_∞, "commutative chiral," "ordered," "unordered," "locality," or "R(z) = τ," the assistant MUST state the proposed edit and wait for explicit author approval. No exceptions.

**AP51 — The three-tier picture is WITHIN E_∞, not a division between E_∞ and E₁.** The manuscript correctly identifies three tiers of R-matrix behavior: (i) pole-free commutative (R = τ), (ii) vertex algebras with poles (R ≠ τ, derived from local OPE), (iii) genuinely E₁/nonlocal (R ≠ τ, independent input). Tiers (i) and (ii) are BOTH E_∞-chiral. Only tier (iii) is E₁. The three-tier picture is a refinement WITHIN the E_∞ class plus one class outside it. NEVER present it as "E_∞ vs E₁" — that is only the (i+ii) vs (iii) division.

**AP52 — CONSEQUENCE: NEVER revert correct manuscript text based on a false premise.** On 2 April 2026, the assistant reverted the entire file ordered_associative_chiral_kd_core.tex to its committed state, destroying ~200 lines of correct insertions (sl₃ computation, two-colour architecture Remark, small vs homotopy distinction, gravitational Yangian Remark, annular bar differential), because it FALSELY believed the pre-existing E₁/E_∞ text was wrong. The pre-existing text was correct. The assistant's own parenthetical additions were the problem. The correct action was surgical removal of the contaminated parentheticals, not a full file revert. NEVER use `git checkout HEAD -- file` on a file with mixed correct/incorrect edits. Use targeted Edit operations instead.

**AP53 — CONSEQUENCE: NEVER trust an agent's literature search over the author's explicit statement.** On 2 April 2026, a literature agent concluded "vertex algebras are NOT E_∞-chiral" based on reading BD Chapter 4 (which defines "commutative chiral algebras" as pole-free). The author explicitly corrected this: "ALL vertex algebras are E_∞-chiral." The assistant then went back and forth between the agent's conclusion and the author's correction, causing multiple rounds of damage. The author's explicit mathematical statements about their OWN conventions ALWAYS override agent literature searches. An agent can report what a reference says; it cannot override the author's definition of terms in the author's own manuscript.

**AP54 — CONSEQUENCE: NEVER propagate a correction to multiple files before verifying it is correct.** On 2 April 2026, the assistant made an incorrect "correction" (replacing "E_∞" with "abelian Lie bracket" as the R-matrix triviality condition) and then propagated it to 6 locations across 3 files before verifying it was right. Each propagation amplified the damage. Then it "corrected" the correction (replacing "abelian Lie bracket" with "no OPE singularities") and propagated THAT to the same 6 locations — creating a second wave of damage. The rule: make ONE edit, verify with the author, THEN propagate. Never batch-propagate an unverified change.

**AP55 — CONSEQUENCE: NEVER add "in the sense of BD" or "in the sense of [reference]" as a parenthetical without verifying the reference actually supports the gloss.** The assistant added "in the sense of BD — those whose chiral product has no OPE singularities" as a gloss on "commutative chiral algebras." BD's Chapter 4 does define "commutative chiral algebra" as pole-free, but "commutative chiral algebra" in the MANUSCRIPT's usage means something broader (= E_∞ = local = all vertex algebras). The BD reference was cited accurately but applied to the WRONG term — the manuscript's "commutative" is not BD's "commutative." Citing a reference correctly does not mean the citation is appropriate in context.

**AP56 — CONSEQUENCE: The working notes contain ~800 lines of content added on 2 April 2026. Some of this content contains contaminated E₁/E_∞ language that has been partially but not fully corrected.** Before using any content from working_notes.tex sections dated 2 April 2026 (sections on "dg-shifted Yangians," "two-colour Koszul duality architecture," "curvature of the Yangian," etc.), verify that every E₁/E_∞ statement matches the author's convention: E_∞ = local = all vertex algebras, E₁ = nonlocal = quantum vertex algebras.

**AP57 — PVA ≠ P_∞-chiral algebra. These are COMPLETELY DIFFERENT objects at different levels of the hierarchy.** A Poisson vertex algebra (PVA) sits between a commutative associative algebra and an E_∞-chiral algebra: it is the CLASSICAL SHADOW obtained by passing to cohomology of an E_∞-chiral algebra, retaining the commutative product (from the regular OPE) and the Poisson bracket (from the singular OPE), with the PVA axioms encoding their compatibility. A P_∞-chiral algebra sits between an E_∞-chiral algebra and an E₁-chiral algebra: it is a HOMOTOPY-COHERENT intermediate structure that remembers some but not all of the Σ_n-equivariance. The hierarchy is: commutative associative ⊂ PVA ⊂ E_∞-chiral ⊂ P_∞-chiral ⊂ E₁-chiral. NEVER conflate PVA (a classical/cohomological object on the LEFT side of this chain) with P_∞-chiral (a homotopy-coherent object on the RIGHT side). NEVER say "the PVA is the P_∞-chiral algebra" or "the P_∞ structure gives a PVA." The PVA is obtained by DESCENDING from E_∞-chiral to cohomology. The P_∞-chiral algebra is obtained by RELAXING E_∞-chiral toward E₁. These are opposite directions.

**AP58 — The full hierarchy of chiral operadic structures.** From most to least symmetric: commutative associative (no poles, no bracket) → PVA (commutative product + Poisson bracket, classical shadow) → E_∞-chiral = vertex algebra (local, Σ_n-equivariant, OPE with poles) → P_∞-chiral (partial symmetry, homotopy-coherent intermediate) → E₁-chiral = nonlocal/quantum vertex algebra (ordered only, no Σ_n). The bar complex and Koszul duality operate at the E_∞ and E₁ levels. The PVA is the associated graded / classical limit / cohomological shadow. NEVER place PVA at the same level as E_∞-chiral or P_∞-chiral — it lives one categorical level down (on cohomology, not on chain complexes).

**AP35 — Accidentally correct theorem (false proof, true conclusion).** A correct answer is NOT evidence of a correct proof. When two errors cancel (orientation double-error, factorial/derivative double-error), fix BOTH — the cancellation is accidental and will break under generalization. In Vol II: orientation convention +∂_{ε_S} vs correct -∂_{ε_S} in fm-proofs appendix — double error cancelled giving correct final signs. Elliptic r-matrix: factorial AND derivative order both off by one — errors partially cancelled (correct pole structure, wrong regular part). **Rule: verify proof steps independently of the conclusion.**

**AP36 — Biconditional overclaim (⟹ proved, ⟺ claimed).** "Convolution formality = one-channel" retracted to one-directional. DS-KD intertwining claimed for "all Koszul at generic level" but only proved on abelian-n+ lane. "Filtration formality" conflated with "dg formality." **Rule: before writing "iff," verify BOTH directions independently.**

**AP37 — Spectral sequence page from pole order alone.** W_N collapse claimed at E₄; correct is E_{2N} for N ≥ 3 (pole order 2N drives d_{2N-1} differentials). Ordered bar E₁ page identified as Lie homology; correct is Hochschild homology (tensor ≠ exterior construction). **Rule: compute spectral sequence pages from full differential structure, not pole-order heuristics.**

**AP38 — Literature normalization convention in hardcoded values.** Faber-Pandharipande test expectations: λ₂ = 1/1152 was wrong; correct is 7/5760. Engine correct but test hardcoded wrong expected value. F₁ coefficient: -(k/24)log η(τ) → -k log η(τ). **Rule: record source paper and normalization convention in comments. Derive expected values independently (AP10).**

**AP39 — κ ≠ S₂ for non-Virasoro families.** S₂ = c/2 ≠ κ = dim(g)·(k+h∨)/(2h∨) for affine KM. They coincide only for Virasoro and Heisenberg. 8+ table values corrected. **Rule: table headers must distinguish S₂ (arity-2 shadow) from κ (modular characteristic). NEVER copy between families.**

**AP40 — LaTeX environment contradicts claim status tag.** Three frontier conjectures in \\begin{theorem} despite \\ClaimStatusConjectured. thm:Koszul_dual_Yangian tagged ProvedHere but cited "(Dimofte-Niu-Py, Thm 5.5)" in the statement. **Rule: environment MUST match claim status. Systematic check: grep for mismatches.**

**AP41 — Prose mechanism ≠ mathematical mechanism.** "The residue extracts the simple-pole coefficient" — wrong, the bar kernel extracts ALL modes via d log(z-w). Formulas correct; English description wrong. **Rule: verify prose descriptions match formulas. One-sentence summaries that skip steps are lies of omission.**

**AP42 — Correct at sophisticated level, false at naive level.** Slogans like "CY-A works for all d" or "scattering = shadow tower" capture deep truths but fail when instantiated naively. **Rule: state the level of validity explicitly from the first occurrence.**

**AP43 — Central object defined by aspiration, not by axioms.** Objects used in theorem statements without formal \\begin{definition}. A property list is a conjecture about an object, not a definition. **Rule: the central object MUST be formally defined before use.**

**Meta-principle:** Every error traces to confusing two objects that share a name, formula, or special-case coincidence. The meta-rule: never trust a coincidence. Verify at the most general case, highest weight, most general level, most general family. **The meta-meta-meta-rule (from AP35-AP43): the same error can recur at different EPISTEMIC levels.** AP35-AP39 catch errors in the VERIFICATION LAYER (false proofs, biconditional overclaims, convention mismatches). AP40-AP43 catch errors in the COMMUNICATION LAYER (environment/tag mismatches, prose lies, slogans without scope, undefined objects).

## Cross-Volume Bridges

| Bridge | Vol II claim | Vol I anchor | Status |
|--------|-------------|--------------|--------|
| Bar-cobar | SC^{ch,top} bar-cobar specializes Vol I Thm A when curve = ℂ, topological = ℝ | Theorem A | Proved |
| DS-bar | Bar-cobar commutes with DS reduction | Theorem ds-koszul-intertwine | Proved (Vol I) |
| Hochschild | BV-BRST origin of Vol I's Theorem H complex | Theorem H | Proved (all genera) |
| DK/YBE | r(z) = ∫₀^∞ e^{-λz}{·_λ·}dλ provides DK-0 shadow | MC3 (DK extension) | Proved (Laplace) |
| PVA-Coisson | PVA descent at X = pt recovers Coisson structure | Deformation theory | Proved |
| W-algebras | Feynman-diagrammatic m_k matches bar differential at genus 0 | MC5 (BRST = bar) | Proved (genus 0); conjectural at g≥1 (conj:master-bv-brst) |
| Affine monodromy | Reduced HT monodromy = quantum group R-matrix; C_line^red ≃ Rep_q(𝔤) on eval modules; Jones polynomial from bar complex | Thm A + affine half-space BV + Drinfeld-Kohno | **Proved** (affine lineage) |
| Soft theorems | Shadow tower controls soft graviton hierarchy via soft order p ↔ arity r=p+2 | Theorem H + Movement I | Proved (genus 0) |
| W_N Koszul | α_N = 2(N−1)(2N²+2N+1) generalises Virasoro c→26−c to all W-algebras | Theorem B + DS-bar | Proved |
| Wick anomaly | Genus tower measures Wick rotation breaking: F_g ≠ 0 ⟺ κ_eff ≠ 0 ⟺ S-transform anomalous; quartic pole → κ → d² = κ_eff·ω_g → F_g | Theorem D + Movement IV | Proved (genus tower); Conjectured (Cardy extraction) |
| Two-colour architecture | ordered bar → A^!_line, symmetric bar → A^!_ch; R-matrix descent | rem:two-colour-architecture, thm:two-color-master | Proved |
| Annular bar-Hochschild | B^{ann}(A) computes HH^{ch}_•(A); genus-1 ordered sector | thm:annular-bar-differential | Proved |
| FG-shadow-stratification | Commutator filtration spectral sequence; E_1-page = FG bar | conj:FG-shadow (now theorem) | Proved |
| Gauge-gravity dichotomy | m_k=0 (gauge) vs m_k≠0 (gravity); DS transports L→M | rem:gauge-gravity-yangian-dichotomy | Proved |

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
