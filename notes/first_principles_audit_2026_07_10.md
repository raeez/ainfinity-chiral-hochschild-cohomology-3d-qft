# First-Principles Audit of the Mathematics and Physics — Vol II

**Date:** 2026-07-10.
**State audited:** HEAD `525accd` ("release pdf", 2026-06-18) plus the
working tree (383 files modified vs HEAD at audit time; no mid-audit
drift observed in the audited files). Line numbers are as of this read.
**Method:** two independent adversarial deep reads — (i) the
SC^{ch,top} operad and the topologisation tower; (ii) the Part VI
3d-gravity climax, the Theorem C conductor set, Brown–Henneaux, the
chiral quantum group Q_g^{k,f,μ}, and celestial holography — each
instructed to ignore claim-status macros and referee the actual
statements and proof bodies. Load-bearing findings were independently
confirmed at the main line (marked ✓M below) by direct reading or
symbolic recomputation. Companion audits: Vol I
(`~/chiral-bar-cobar/notes/first_principles_audit_2026_07_10.md`),
Vol III, MHT, igusa-cusp-form; ecosystem copies in
`~/ecosystem/swarm-reports/`.

---

## 1. Executive verdict

**Overall grade: C. The epistemic scaffolding is unusually good; the
two headline structures underneath it are damaged** — one by a
foundational type defect (the operad's composition is not closed),
one by a fabricated constant in the flagship invariant set (98/3).

- Conditional claims mostly carry their hypothesis packages
  *in-statement*; the gravity climax explicitly denies producing a
  dynamical-metric QFT, a microscopic Hilbert space, or a proof that
  Z^der_ch is the QG Hilbert space; the celestial soft brackets are
  exactly right; obstruction classes are named.
- But the defining composition of SC^{ch,top} is undefined on
  mixed-into-open insertions (✓M); the volume maintains two
  inequivalent definitions of its central operad and proves theorems
  about whichever is convenient; the Koszulity theorem is a citation
  transfer to literature about the *other* model; and the Theorem C
  value 98/3 comes from an invented "Arakawa-normalized" BP central
  charge whose Feigin–Frenkel sum is 196, while the true BP conductor
  — proved in Vol I — is 50 (✓M).

## 2. Independently verified at the main line (✓M)

1. **The BP falsification.** True Bershadsky–Polyakov central charge
   c(k) = 25 − 24/(k+3) − 6(k+3) (= −(2k+3)(3k+1)/(k+3)); c(0) = −1;
   c(k) + c(−k−6) ≡ 50 identically (the 1/(k+3) and (k+3) terms are
   odd under k ↦ −k−6). Ghost check: spins (2, 3/2, 3/2, 1) give
   −26 − 11 − 11 − 2 = −50, confirming criticality at 50. The
   inscribed c^Ar(k) = 2 − 24(k+1)²/(k+3)
   (`chapters/theory/bp_chain_level_strict_platonic.tex:794–830`,
   read verbatim ✓M) has c^Ar(0) = −6 ≠ −1 — it is not the BP central
   charge in any standard convention — and its FF sum is
   4 + 24[(k+5)² − (k+1)²]/(k+3) = 4 + 192 = 196, whence
   98/3 = 196/6.
2. **The missing composition case.** `chapters/theory/foundations.tex:2930–2955`
   (read verbatim ✓M): composition is defined in exactly three
   bullets — closed-into-closed, closed-into-open, and open-into-open
   *where every inserted open-output operation is a pure E₁(n_j)
   element*. Insertion of a mixed operation (FM_{k'}(ℂ) × E₁(n),
   k' ≥ 1) into an open slot — required for any bulk–boundary
   calculus and for every tree of the bar/cobar of the operad — is
   absent; componentwise it would need a merge
   FM_k(ℂ) × FM_{k'}(ℂ) → FM_{k+k'}(ℂ), which does not exist (the
   relative position of the two clusters is lost by the
   translation/dilation quotient).
3. W_N conductor formula: complementary-pair sum
   c_N(b) + c_N(ib) = 2(N−1)(2N²+2N+1) = 4N³−2N−2, giving 26 (Vir),
   100 (W₃), 246 (W₄) — so the set entries 13 = 26/2 and
   250/3 = (5/6)·100 are correct arithmetic (✓M, symbolic).

## 3. Pillar-by-pillar findings

### 3.1 SC^{ch,top} and the topologisation tower — grade **C−**

Definition sites: product model
`foundations.tex:2857–2986` (+ duplicate `locality.tex:199–230`);
Voronov half-plane model `sc_chtop_heptagon.tex:232–315`;
"product = only depth-zero associated graded"
`factorization_swiss_cheese.tex:773–798`.

- **F-II.1 (composition not closed; two inequivalent definitions).**
  Mixed-into-open insertion is undefined in the product model (✓M,
  §2.2). The heptagon's face (1) uses Voronov's coupled FM_k(ℍ) model
  and `factorization_swiss_cheese.tex:792–796` states the separated
  product is *not* the global parameter space — directly
  contradicting `foundations.tex:2891–2920`, where the product with
  componentwise composition *is* the operad "forced by locality".
  The heptagon therefore identifies seven presentations of an object
  whose two definitional models are inequivalent. (Also the fibre
  product at `sc_chtop_heptagon.tex:242` is literally FM_k(ℍ) —
  broken notation.)
- **F-II.2 ("dioperad" is a category error).** Gan dioperads are
  multi-output structures; a one-output dioperad is an operad. The
  justification "Σ_n-equivariance fails across colours"
  (`foundations.tex:3017–3024`) misreads coloured-operad
  equivariance, which permutes inputs together with their colour
  signature. Voronov's Swiss-cheese is a 2-coloured symmetric operad
  with exactly the claimed directional restriction. The restriction
  itself IS closed under composition (`locality.tex:229`) — that part
  is sound (trivially).
- **F-II.3 (Koszulity by citation transfer; Livernet uncited).**
  `line-operators.tex:84–201` asserts Koszulity of "the classical
  Swiss-cheese operad" via Voronov/Hoefel/Hoefel–Livernet — literature
  about the honest *coupled* operad — while the comparison map
  Φ = (formality) × id exists only for the decoupled product model.
  Livernet's non-formality of Swiss-cheese (2015) is cited nowhere in
  the volume and obstructs any formality-based identification for the
  coupled model. Either reading leaves a hole. `prop:gr-chiral`
  (`line-operators.tex:283–303`, "gr ≅ P_ch ⊔ E₁, no mixed
  operations") is false as stated: gr of their own definition has
  components gr(FM_k) ⊗ E₁(m); "decouple" ≠ "vanish". This
  proposition feeds the bar–cobar adjunction, filtered-Koszul, the
  (Lie, Ass) dual, and the Yangian recognition chain.
- **F-II.4 (tower assumes its hard part).** G^{(n)} is genuinely
  constructed (`e_infinity_topologization.tex:349–357`, Koszul–Tate
  contraction :463–491) and the pairwise-commutativity computation is
  real (the referee verified the Kirillov–Kostant step and the
  1/(p+q) normalisation). But: affine KM E₃ is **cohomology-only** by
  its own statement (`3d_gravity.tex:7739–7742`), with the chain-level
  [Q, G̃] = T_Sug "a genuine computation not carried out in the
  present volume" (`e_infinity_topologization.tex:706–707`) — MA-15's
  T1 open at chain level in the flagship case, T4 nowhere addressed.
  The ladder rung E_{k+2} for k ≥ 2 assumes hypTLift ("k pairwise
  commuting topological translation operators with coherent higher
  homotopies", :754–759), and the proof concedes the new direction
  exists "only because that independence is part of the hypotheses"
  (:810–812). A translation operator plus homotopy is not an
  E₁-multiplication; no factorization structure along the W-flow
  directions is constructed for any example. Even the ProvedHere W_N
  specialisation assumes hypTLift.
- **F-II.5 (chiral higher Deligne: two coherences proved, all
  asserted).** `chiral_higher_deligne.tex:326–349` proves the first
  two Stasheff coherences; :351–370 self-authorises "through all
  degrees" by a "uniform Stokes template" remark. Classically the
  brace identities are strict; here they hold up to homotopy, so full
  B_∞ coherence is genuinely infinite data — asserted, not built.
- **F-II.6 (curved-Dunn H² = 0 at g ≥ 2 is a conditional transport
  with a broken step).** `curved_dunn_higher_genus.tex:241–291` is
  honestly Conditional (if H²(bootstrap) = 0 and a bridge Φ exists,
  then H²(TCo_g) = 0), with neither hypothesis discharged in any
  example (:340–345). In the bridge proof (:153–166) the step "adjusts
  ξ̃ by a coboundary to make it d-closed" is wrong — adding a
  coboundary does not change dξ̃; salvageable from hypothesis (3),
  but the written argument is broken. The CLAUDE.md headline
  overstates a conditional implication.
- **F-II.7 (Dunn non-additivity is a remark, not a result).**
  True as a type observation (`foundations.tex:3176–3223`), but no
  theorem, no invariant, no counterexample; its supporting
  "dioperad/equivariance" rationale is F-II.2.
- **Sound and worth keeping:** the Casimir–antighost calculus; the
  W_∞[λ] → E_∞ endpoint carrying hypProchazka/CKL/PRSh/Yamada
  *in-statement* (`e_infinity_topologization.tex:1090–1205`, threshold
  N₀(w) = 2w−1 table :1399–1414) — with one slip:
  `sc_chtop_heptagon.tex:3211–3213` calls the Conditional endpoint
  "proved".

### 3.2 Part VI climax, conductors, quantum group, celestial — grade **C+**

- **F-II.8 (the fabricated 98/3; CONTRADICTS VOL I).**
  `bp_chain_level_strict_platonic.tex:794–859`
  (`prop:bp-kappa-conductor`, tagged ProvedElsewhere): defines
  c^Ar_BP(k) := 2 − 24(k+1)²/(k+3) ("Arakawa-normalized"), derives
  K^Ar ≡ 196, multiplies by an underived "generator-profile
  coefficient ϱ_BP = 1/6" to get 98/3, and cites Vol I's BP theorem
  as "matching" — but Vol I
  (`chapters/examples/bershadsky_polyakov.tex:365–411`) proves the
  conductor ≡ **50** with the correct central charge, and Vol I's
  census (`landscape_census.tex:2853–2875`) declares the BP κ-sum
  **Open**. My recomputation (§2.1): c^Ar is not the BP central
  charge in any convention (c^Ar(0) = −6 vs true −1); the true FF sum
  is 50; the ghost anomaly independently gives 50; ϱ_BP = 1/6
  contradicts the volume's own Σ_{h}1/h generator-profile convention
  (spins 1, 3/2, 3/2, 2 give 17/6, as Vol I records). The proposition
  itself concedes in-body that the Fateev–Lukyanov convention gives
  K^FL ≡ 50. The 98/3 value is propagated through ≥10 files including
  `CLAUDE.md:31`, `preface.tex:1429–1433` (misattributed to Vol I),
  `conclusion.tex:1358, 1376–1384` (a second "witness"
  W₂^{su(2|1)} at c = 49/3 that has no derivation anywhere in the
  corpus), `sc_chtop_heptagon.tex:846–859`,
  `bar-cobar-review.tex:2184–2213`, and numerology at
  `introduction.tex:2308–2311`. **Repair (forced):** restore the BP
  slot to Open with conductor 50; derive or delete the su(2|1) row;
  restate the set as {0, 8, 13, 250/3} ∪ {open BP slot}.
- **F-II.9 (the set is not a bijection with the families; B-row
  double-booked).** Sound entries: 0 for G, L, C (FF-duality
  κ(V_{−k−2h∨}) = −κ; H_k: k + (−k); βγ: −1+1 —
  `conclusion.tex:1350–1353`); 13 for M-Virasoro (c + c^! = 26);
  250/3 for M-W₃ (three paths, `landscape_census.tex:2835–2851` Vol I).
  So G/L/C share 0 and M carries {13, 250/3, (open BP)}. The B-row 8
  (`conclusion.tex:1372–1374`, K^κ(H_{Δ₅}) = 2c₊(II₄,₂₀) = 8 via
  Bruinier) is never derived as κ^Hodge(A) + κ^Hodge(A^!) for any
  algebra, and `theorems_C_D_native_vol2_platonic.tex:347–361`
  computes ρ_K(K3×E) = **0** on archetype B, calling the Mukai 8 "a
  separate sixth invariant, outside the κ-tuple rows". The honest
  passage exists (`introduction.tex:2250–2254`: a "combined value
  set"); the CLAUDE.md phrasing "∈ {0,8,13,250/3,98/3} on G/L/C/M/B"
  is unfaithful shorthand.
- **F-II.10 (Φ_hol and the gravity reading — sound as labelled,
  borderline tautological).** Source/target ∞-categories defined
  (`universal_holography_functor.tex:339–436`; master theorem
  `programme_climax_platonic.tex:127–229`, Conditional), but the HT
  model, boundary comparison, and bulk–Hochschild comparison are all
  part of the input datum Ξ — the functor largely repackages supplied
  realization data. The scope remark `rem:climax-qg-scope`
  (:1367–1423) explicitly denies (a) a dynamical-metric QFT, (b) a
  microscopic Hilbert space, (c) Z^der_ch = QG Hilbert space; the
  Maloney–Witten convergence problem is correctly flagged
  (`thqg_perturbative_finiteness.tex:4655–4659`). Weakness: Part VI
  leans on Witten 2007; the modern Virasoro-TQFT literature
  (Collier–Eberhardt–Mühlmann–Zhang) is cited only in Part VII and
  never confronted in the gravity part; genus-2 partition-function
  subtleties unengaged.
- **F-II.11 (Brown–Henneaux hygiene — sound).** c = 3ℓ/2G_N is an
  assumed dictionary in a setup environment
  (`thqg_perturbative_finiteness.tex:1414–1420`); F₁ = c/48 = ℓ/32G
  checks; BTZ/Cardy inputs (modular invariance, vacuum dominance,
  saddle) are named as separate hypotheses. Nowhere is semiclassical
  gravity claimed to be derived from Koszul duality.
- **F-II.12 (Q_g^{k,f,μ} — programmatic; "unified" overclaims).**
  `thqg_gravitational_yangian.tex:127–202`: coproduct and R-matrix
  are "the separate Drinfeld/RTT input"; antipode excluded;
  uniqueness only on the verified fibre; the existence of Δ_z as an
  E₁-chiral coproduct on W_k(g,f) is asserted via an unproved
  "descends unambiguously" through KRW BRST. The eight-fibre table is
  known objects. The fine print is honest; the "unified chiral
  quantum group" headline is a template, not a construction.
- **F-II.13 (celestial — sound, one discipline gap).** The soft
  bracket [w^p_m, w^q_n] = (m(q−1) − n(p−1))w^{p+q−2}_{m+n}
  (`w_infty_e_infty_endpoint_platonic.tex:695, 755`) matches
  Strominger exactly; classical-w vs quantum-W_{1+∞} (Moyal) handled
  correctly. The universal theorem
  (`universal_celestial_holography.tex:238–293`) is genuinely scoped
  (QME anomaly + gluing class as hypotheses; iff on Mellin cone
  obstructions). Gap: the file has zero claim-status tags, violating
  the repo's own convention; "category of 4d QFTs" in clause (i) is
  not pinned down.

## 4. Systemic diagnosis

Vol II's failure mode differs from Vol I's. Vol I promotes
definitions to theorems; Vol II (i) **maintains two inequivalent
models of its central object** and lets each theorem cite the
convenient one, and (ii) **let one manufactured number into the
flagship invariant set** — the exact Beilinson-cascade failure the
constitution warns about, complete with a false cross-volume
citation. Meanwhile the *conditional* architecture (hypothesis
packages in-statement, named obstructions, honest gravity scoping) is
the best in the constellation after MHT. The volume needs structural
repair of SC^{ch,top} (adopt the Voronov coupled model as THE
definition; re-prove or re-scope Koszulity against it; cite and
confront Livernet) more than it needs more theorems.

## 5. Triage (ordered)

1. **Excise 98/3**: restore the BP conductor to 50 with the slot
   Open (as Vol I has it); delete or derive the c = 49/3 su(2|1)
   witness; fix `CLAUDE.md:31`, `preface.tex:1429–1433`,
   `conclusion.tex:1350–1384`, `sc_chtop_heptagon.tex:846–859`,
   `bar-cobar-review.tex:2184–2213`, `introduction.tex:2308–2311`.
2. **Fix the B-row double-booking**: 8 is the Mukai/Bruinier
   conductor, not a ρ_K value; state the set as
   {0, 8(Mukai), 13, 250/3} ∪ {open BP}.
3. **One definition for SC^{ch,top}** (the coupled FM_k(ℍ) model);
   define mixed-into-open composition; delete the "dioperad"
   framing; re-scope `thm:homotopy-Koszul` and fix
   `prop:gr-chiral`; cite Livernet 2015 and address non-formality.
4. **Downgrade the ladder**: E₃-on-cohomology for KM stays (with its
   package); rungs k ≥ 2 restated as conditional on hypTLift with
   the E₁-structure gap named; inscribe the chain-level [Q,G̃] = T_Sug
   computation or keep T1 open explicitly.
5. Fix the curved-Dunn bridge-proof step (:153–166) and the
   "proved" slip at `sc_chtop_heptagon.tex:3211–3213`; add status
   tags to `universal_celestial_holography.tex`.

**Provenance.** ✓M items verified by direct main-line reading or
symbolic recomputation; all other findings are referee reports whose
decisive quotes were spot-checked where marked. Line numbers drift
under the active working tree; committed state `525accd` is the
stable reference.

---

# Part II — Mathematical yield (fresh-eyes pass, same date)

Stricter second pass: referee forbidden from reading CLAUDE.md,
FRONTIER.md, notes/, status appendices, or Part I; graded only
**true + proved + new**, with hypothesis-contains-conclusion counting
as zero. Note: the working tree drifted mid-audit
(`hochschild.tex` `comp:bulk-virasoro` flipped from ProvedHere to
Conditional between two reads — the concurrent editing session is
active in this repo too).

**Yield grade: C−.**

1. **Fatal gap in the load-bearing theorem.**
   `thm:casimir-antighost-commutativity`
   (`e_infinity_topologization.tex:429–456`): to apply the contraction
   to [G⁽ⁿ⁾, G⁽ᵐ⁾] the proof must show it is Q-*closed*; instead it
   concludes "[Q_tot,[G⁽ⁿ⁾,G⁽ᵐ⁾]] = 0 in cohomology" (:507–532) —
   vacuous, since any [Q,X] is exact. Strict closure is false:
   Q[G⁽ⁿ⁾,G⁽ᵐ⁾] = [T⁽ⁿ⁾,G⁽ᵐ⁾] ± [G⁽ⁿ⁾,T⁽ᵐ⁾] ≠ 0. The bridging
   identity gr[Q,[G,G]] = h({Pₙ,Pₘ}_KK) is asserted, never derived.
   This theorem discharges axiom (T5) for the E∞-topologisation
   ladder — the volume's single thesis — so the gap propagates to the
   climax. (The Koszul contraction, the 1/(p+q) normalization, and
   the Kirillov–Kostant vanishing were re-verified and are correct;
   the defect is precisely the closure step.)
2. **The advertised Virasoro derived-centre profile is unproved and
   mismatched.** main.tex (~2586–2592) advertises (1,0,1) in degrees
   {0,1,2}; the only theorem behind it (`hochschild.tex:370–405`)
   *cites* Bakalov–De Sole–Kac — whose computation has support
   {0,2,3} with dims (1,0,1,1) — plus a comparison hypothesis that
   contains the conclusion. BDSK has no bibliography entry.
3. **Outright false ProvedHere:** `comp:bulk-heisenberg`
   (`hochschild.tex:2684–2712`) claims H¹ = 0 ("all derivations
   inner") — refuted by the explicit outer derivation α ↦ |0⟩ (no
   weight-0 zero mode realizes it) and by the BDSK profile (2,1,0,…)
   recorded in the *same file* (:407–430, :5832). Its opening
   sentence ("concentrated in degree 0") contradicts its own H² ≠ 0
   three lines later.
4. **Curved-Dunn:** no unconditional H² computation of any named
   complex for any algebra at any genus exists in the chapter; the
   single unconditional item unwinds the definition ρ_K := κ + κ^!
   (c/2 + (26−c)/2 = 13). The citation "Markl–Shnider–Stasheff 2002
   Theorem 3.42 sharpens to n ≤ g+1" (:1686–1688) is decorative —
   that theorem concerns nothing of the sort.
5. **F₁ = c/48 is definitional** (F₁ := κ∫λ₁ = (c/2)(1/24); both
   factors classical); the surrounding torus identities (Zhu
   recursion, P₂ = ℘ + E₂/12, q∂_q log η^{−c/2} = −(c/48)E₂) are
   correct and standard. Nothing approaches Maloney–Witten.

**The genuinely new proved mathematics of the volume** (complete list
found): the symmetric-point combinatorics of the recursively defined
Virasoro Stasheff tower (`3d_gravity.tex:700–1985`) — even-arity
vanishing by telescoping, the odd-arity **Catalan factorization**
φ₂ₙ₊₃(x) = (−1)ⁿ Cₙ (x+2)⋯(x+2n+3), the explicit gauge-fixed m₄/m₅
formulas, c-linearity of the scalar sector, and the weight-depth
grading. The referee reproduced m₃ by hand from the Virasoro
λ-bracket, the Obs₄ factorization by an independent sympy
implementation, and the full Catalan induction (checked k = 4–7):
correct, real, new — and lemma-grade, gauge-dependent (the
identification with an A∞ minimal model of Vir_c is not established;
with m₁ = 0 the tower does not satisfy the standard A∞ relations, as
the proposition itself concedes at :808–812, and the recursion is
verified computationally only through arity 13).

**Consequences for Part I:** Part I's praise of the antighost
calculus as "real mathematics" stands for the contraction and KK
steps but is **overturned at the theorem level** — the
commutativity theorem has a fatal, load-bearing gap. Triage
additions (top priority): repair or weaken
`thm:casimir-antighost-commutativity` (zero-mode charges, or a
statement modulo the lower-Casimir ideal); retract or re-prove
`comp:bulk-heisenberg` H¹ = 0; reconcile the (1,0,1) advertisement
with BDSK's {0,2,3} and add the missing bibliography entry; delete
the false MSS citation.

---

# Part III — Healing ledger (2026-07-10, same date)

All seven triage items executed against the working tree; every edit
re-read its file immediately before writing (the tree is under
concurrent edit — `comp:bulk-virasoro` was found already healed to
Conditional with the BDSK {0,2,3} content and was preserved, and
`factorization_swiss_cheese.tex` drifted mid-session and was
re-read before its two edits). Statuses moved only by proof;
hypothesis-containing conclusions split into proved part + named
Conjecture. Line numbers below are post-edit.

## 1. 98/3 excised (F-II.8) — the BP slot is OPEN, the conductor is 50

- `chapters/theory/bp_chain_level_strict_platonic.tex:105–147`
  (`rem:anchor-bp-self-duality`): the fabricated "Arakawa convention
  c^Ar = 2 − 24(k+1)²/(k+3), K^Ar ≡ 196" deleted; the remark now
  states the true FL central charge c_BP(k) = 25 − 6(k+3) − 24/(k+3)
  = −(2k+3)(3k+1)/(k+3) (c_BP(0) = −1), the odd-under-t↦−t proof of
  K^c ≡ 50, and the ghost cross-check −26−11−11−2 = −50; the false
  claim that Vol I `standalone/bp_self_duality.tex` proves 196 (it
  proves 50, verbatim checked) removed.
- `…:794–884` (`prop:bp-kappa-conductor` + proof): full replacement.
  Now states K^c_BP ≡ 50 (ProvedElsewhere, Vol I
  `thm:bp-koszul-conductor-polynomial`), ghost-anomaly confirmation,
  and the κ-complementarity slot as OPEN: no derivation of the BP
  κ-sum exists; ϱ_BP = 1/6 deleted (the volume's own Σ1/h convention
  on weights (1,3/2,3/2,2) gives 17/6 and is itself not a
  derivation — stated in-proof). Label `eq:bp-kappa-complementarity`
  removed (no dangling refs; verified).
- `rem:bp-role-of-K-196` → `rem:bp-role-of-K-50` (no external refs);
  end-of-chapter summary (`sec:bp-summary-forward`) updated to 50 +
  open slot.
- Propagation, every site verified clean by final grep (0 hits for
  98/3, 49/3, BP-context 196, ϱ_BP):
  `curved_dunn_higher_genus.tex:1537–1560` (BP bullet of
  `prop:chiral-deligne-conductor-identity` + proof: slot open,
  conductor 50 named as distinct); `modular_swiss_cheese_operad.tex`
  (`cor:kappa-zero-crit` hypothesis set); `introduction.tex` (M-row
  paragraph, family-wise-ceiling paragraph, Schur-sector clause (C),
  H_{Δ5}/Lee–Yang 245/33 numerology clause deleted);
  `sc_chtop_heptagon.tex` (chain-homotopy normalisation passage —
  also repaired the deeper conflation: eq now uses the Vol I
  Verdier–Ran-on-bar norm N(A) = κ + κ^!_LV with witness values
  N(G)=1, N(L)=2(k+h∨), N(M)=c(5c+22)/10, N(B)=8, stated distinct
  from ρ_K with computed spectrum {0,13,250/3});
  `theorems_C_D_native_vol2_platonic.tex:288–297`; `preface.tex`
  (W₃ paragraph — mis-attribution "Vol I five-archetype bucket"
  removed — and ceiling paragraph);
  `thqg_perturbative_finiteness.tex` (`rem:thqg-pf-Kkappa-8`);
  `conclusion.tex` (complementarity table: W₂^{su(2|1)} c = 49/3 row
  DELETED entire — underived everywhere in the corpus — and the
  two-invariant B-row phrasing installed);
  `bar-cobar-review.tex` (`rem:bar-cobar-2-theoremC-bucket`
  retitled and rewritten: Mukai 8 = sixth invariant, su(2|1) witness
  deleted); `ordered_associative_chiral_kd_core.tex`
  (`rem:ordered-kd-core-conductor-bridge`: K(BP) = 50, ϱ_BP removed,
  bridge scoped to principal W_N);
  `programme_climax_platonic.tex:584–592`;
  `anomaly_completed_frontier.tex` (four sites: c_sum(3,1) = 50 with
  true formula; c−25 odd, midpoint at k = −3±2i; κ-sum open; false
  "computed in Volume I" for κ(B^k) → open bar-complex computation;
  `rem:ac-comparison-mukai-bruinier` + section title);
  `anomaly_completed_topological_holography.tex` (two sites; false
  "BP duality conjecture … = 196" → Verdier-companion assumption +
  proved 50; decorative `V1-comp:sl3-ds-hierarchy` citation removed
  after verifying it computes no κ);
  `thqg_gravitational_s_duality.tex` (Trinity value list);
  `thqg_anomaly_extensions.tex` (family table row: 50, 25, κ open);
  `ordered_associative_chiral_kd_frontier.tex` (true c_BP; the
  cross-family SCA–BP sum RECOMPUTED from the true formula:
  (6k²+67k+156)/(k+3) = 6k+49+9/(k+3), still not level-independent,
  conclusion unchanged; same-family identity 50);
  `thqg_symplectic_polarization.tex` (bridge scoped to principal
  W_N; K_BP = 50; no derived BP coefficient);
  `chapters/examples/w-algebras.tex` (κ_W(BP) = 50 ≠ 100);
  `line-operators.tex:2468` (K = 50 line-operator bridge);
  `CLAUDE.md` §3 Theorem C line; `main.tex:2427` Part-IV overview.

## 2. B-row double-booking (F-II.9)

Every site presenting {0,8,13,250/3,…} as "ρ_K on G/L/C/M/B" now
carries the two-invariant phrasing: ρ_K(K3×E) = 0 (Hodge-row sum;
`theorems_C_D_native_vol2_platonic.tex:347–361` was already honest
and is now matched by its own :288–297) and
K^{κ} = 2c₊(II₄,₂₀) = 8 (Mukai/Bruinier, a distinct sixth
invariant). Sites: introduction (table B-row now shows both
invariants; ceiling paragraphs), preface (two ceiling paragraphs),
conclusion (`rem:complementarity-table`), bar-cobar-review,
thqg_perturbative_finiteness, anomaly_completed_frontier,
programme_climax_platonic, main.tex, CLAUDE.md. The computed
ρ_K spectrum reads {0, 13, 250/3} + open BP slot everywhere.

## 3. SC^{ch,top} unification (F-II.1/2/3, F-II.7)

- (a) Coupled model is THE definition:
  `foundations.tex` `def:SC-operations` now defines mixed operations
  as FM_{k,m}(ℍ) (Voronov half-plane compactification of k interior
  + m boundary points mod boundary-affine group);
  `def:SC-composition` now has the full coloured composition
  including **mixed-into-open boundary insertion carrying the
  transverse position**, with the explicit statement that the
  decoupled model cannot express it (the would-be merge
  FM_k×FM_{k'}→FM_{k+k'} forgets relative position);
  `rem:product-forced-by-locality` rewritten: the product
  FM_k(ℂ)×E₁(m) is the depth-zero associated graded of the boundary
  stratification (as `factorization_swiss_cheese.tex`
  `def:mixed-ht-ran-geometry` already proves), and the three-step
  locality argument now proves the gr splitting, not the operad
  product. `locality.tex` `def:SC`/`prop:SC-operad` rebuilt on the
  coupled model with gr⁰ as the computational model;
  `def:SC-alg` distinguishes coupled strict algebras from
  gr⁰-algebras. Duplicate-definition inconsistency resolved:
  heptagon Face (1) now cross-references `def:SC-operations`.
  Downstream sites asserting the product as THE operad re-scoped to
  gr⁰: `bv-construction.tex:353`, `3d_gravity.tex` three-colour
  extension (three sites), heptagon `C^log` dg object,
  `introduction.tex` construction section (full rewrite of
  mixed-to-open + composition paragraphs).
- (b) "Dioperad" framing deleted volume-wide (foundations
  `prop:SC-well-defined` + `rem:sc-chtop-not-operad`, locality
  intro, and 20+ loose occurrences renamed operad; Gan bib entry
  kept). The false "Σ_n-equivariance fails across colours" rationale
  replaced by: coloured-operad equivariance permutes inputs together
  with their colour signature; the directional restriction is
  emptiness of components, closed under composition by the
  empty-target argument (kept). Dunn-additivity non-applicability
  restated on the true ground (coupling, not equivariance failure).
- (c) `thm:homotopy-Koszul` (`line-operators.tex`) restated as
  **gr-level**: the decoupled operad gr SC^{ch,top} is
  homotopy-Koszul (distributive-law Koszulity of the decoupled
  homology operad + Kontsevich formality componentwise + bar-cobar
  transfer — the map Φ = φ×id exists only for the decoupled model,
  now said explicitly). NEW `conj:coupled-sc-koszul` (coupled
  homotopy-Koszulity) and NEW `rem:livernet-obstruction` citing
  Livernet 2015 non-formality as the obstruction; `\bibitem{Liv15}`
  added to main.tex; the coupled/homology-level Hoefel–Livernet
  Koszulity quarantined to the homology level in the proof and the
  three-inputs remark.
- (d) `prop:gr-chiral` corrected: gr splits componentwise
  (gr(FM_k)⊗E₁(m)); mixed components do NOT vanish; the free product
  P_ch ⊔ E₁ appears only after further forgetting the mixed
  components — proof Steps 2–3 rewritten accordingly. Dependents
  conditioned on the corrected statement:
  `thm:scchtop-bar-cobar-adjunction` (equivalence proved for
  gr-algebras; coupled conditional on the conjecture),
  `thm:filtered-koszul` (graded equivalence proved; filtered coupled
  statement conditional on the conjecture + complete exhaustive
  filtration), `rmk:hKoszul-dependence` (Yangian recognition chain
  `thm:lines_as_modules` / `thm:Koszul_dual_Yangian` inherit the
  conditioning), plus every external citation site reconditioned:
  axioms.tex (3), equivalence.tex, modular_swiss_cheese_operad.tex
  (2: incl. `thm:modular-hkoszul-SC`), introduction.tex (main
  theorem restatement + summary), preface.tex (3),
  factorization_swiss_cheese.tex (2), foundations.tex
  (`constr:vol2-en-formality`, `rem:two-color-koszul-duality-operadic`).

## 4. thm:casimir-antighost-commutativity (Part II finding 1 — the fatal gap)

- `e_infinity_topologization.tex:429–…`: theorem REPLACED by the two
  provable statements. (i) Zero-mode statement:
  [Q_tot, G⁽ⁿ⁾₀] = T⁽ⁿ⁾₀ = L⁽ⁿ⁾_{−n+1} (zero-mode functor commutes
  with the BRST contour); (ii) {P_n,P_m}_KK = 0 with the corrected
  invariant-theory computation (type error "[dF_ξ,ξ] = 0" fixed to
  ad*_{dF_ξ}ξ = 0); (iii) gr-level exactness of the zero-mode
  brackets modulo lower filtration, with the bridging identity now
  DERIVED at gr (single J–J structure-constant contractions; Q₀
  reassembles h({P_n,P_m}_KK)₀; level terms lower filtration) rather
  than asserted. Status ProvedHere for exactly (i)–(iii); licensing
  β+γ (hypTLift removed from the theorem, it now lives in the
  conjecture).
- NEW `conj:casimir-antighost-commutativity`: chain-level pairwise
  commutativity [G⁽ⁿ⁾(z),G⁽ᵐ⁾(w)] = [Q_tot,H⁽ⁿᵐ⁾] with coherent
  higher homotopies ≡ the pairwise-commutativity clause of hypTLift
  for Casimir towers. NEW
  `rem:casimir-commutativity-obstruction`: names the precise defect
  (Q[G,G] = [W⁽ⁿ⁾,G⁽ᵐ⁾] ± [G⁽ⁿ⁾,W⁽ᵐ⁾] ≠ 0; "[Q,[G,G]] = 0 in
  cohomology" carries no information; the descending recursion is
  the conjecture's content).
- Honesty propagated: `prop:closed-normal-ordered-antighost-homotopies`
  re-scoped (contraction identity ProvedHere; the commutativity
  reading conditional on the conjecture — proof's closure step now
  says "If … the conjectural input");
  `rem:axiom-T5-scope` + `rem:frontier-antighost-…` now state (T5)
  is NOT discharged for any tower (gr-level zero-mode shadow proved,
  chain level = named Conjecture);
  `thm:e-infinity-topologization-ladder` hypothesis list names the
  conjecture explicitly for the Casimir specialisations;
  `thm:e-infinity-specialisation-WN` proof's (T5) step rewritten
  (enters as hypothesis, not theorem); external claim sites fixed:
  `sc_chtop_heptagon.tex` (class-M edge status: proved →
  conditional-on-conjecture), `preface.tex` ((T5) "discharged" claim
  removed, Vir depth-2 vacuity noted, W_N conditional),
  `part_vi_platonic_introduction.tex`,
  `w_infty_e_infty_endpoint_platonic.tex` (endpoint step (ii)).

## 5. comp:bulk-heisenberg (Part II finding 3) + BDSK reconciliation (finding 2)

- `hochschild.tex` `comp:bulk-heisenberg`: the false H¹ = 0 ("all
  derivations inner") and the self-contradictory "concentrated in
  degree 0 … H² = k⟨ν_k⟩" REPLACED by the BDSK-consistent profile
  (2,1,0,…), support {0,1}, Hilbert polynomial 2+t, with the outer
  derivation D: α ↦ |0⟩ exhibited and proved outer
  ((:α^m:)₍₀₎α = mk ∂:α^{m−1}: is never the vacuum); H² = 0 with the
  level-rescaling exactness (now consistent with the same file's
  own :420 and :5852 BDSK records); status ProvedElsewhere for the
  bounded profile, chart reading conditional on χ^bd; the
  fixed-fibre reduced presentation (`comp:tamarkin-e2-heisenberg`)
  cross-referenced as a different complex with the comparison as
  the grading obligation of `comp:drinfeld-center-heisenberg`.
- `main.tex:2588–2599` advertisement reconciled: "(1,0,1) in degrees
  {0,1,2}" → dims (1,0,1,1) in degrees {0,2,3} (BDSK), transported
  to the chart exactly when χ_c^bd is a quasi-isomorphism —
  matching the current `thm:chirhoch-virasoro-concentration` and
  `comp:bulk-virasoro` (found already Conditional; preserved,
  semantic merge).
- `\bibitem{BDSK21}` (Bakalov–De Sole–Kac, Computation of cohomology
  of vertex algebras, arXiv:2106.00688, + Jpn. J. Math. operadic
  companion) added to main.tex; all five plain-text "Bakalov–De
  Sole–Kac (2021, Thm 7.2/7.4)" sites in hochschild.tex wired to
  \cite{BDSK21}.

## 6. Curved-Dunn (F-II.6, Part II finding 4)

- `curved_dunn_higher_genus.tex` bridge proof
  (`prop:modular-bootstrap-to-curved-dunn-bridge`): the broken
  "adjust ξ̃ by a coboundary to make it d-closed" step REPLACED by
  the correct route through hypothesis (3): lift [ξ] through the
  iso H²(Φ), kill it upstairs by H² = 0, conclude [ξ] = 0.
- Decorative citation "Markl–Shnider–Stasheff 2002 Theorem 3.42
  sharpens to n ≤ g+1" DELETED (both occurrences: tower-bound
  remark and the stray \cite in the ML remarks); the n ≤ g+1 bound
  now honestly attributed to the comparison-theorem hypothesis
  package, with the proved bound n ≤ 3g retained.
- g = 9, 10 remarks retitled "Formal … expressions: unverified
  comparison data" with an explicit caveat sentence; the Borcherds
  legs (Φ₁₀^un)⁵/η¹²⁰ and (Φ₁₀^un)^{11/2}/η¹³² labelled numerology
  (exponents match weights, nothing more claimed); section retitled
  "Formal obs_g expressions…".

## 7. Small fixes

- `sc_chtop_heptagon.tex` `conj:heptagon-collapse-general`: the
  W_∞[μ] endpoint "proved in Theorem…" → "established conditionally
  in …" with the four endpoint hypotheses named.
- `sc_chtop_heptagon.tex:242`: broken fibre product
  FM_k(ℍ)×_{Conf_m(ℝ)}Conf_m(ℝ) → FM_{k,m}(ℍ) with the honest
  description and the composition cases corrected (closed-into-open
  at interior points; mixed-into-open at boundary points with
  transverse position).
- `universal_celestial_holography.tex`: claim-status tags added to
  all 12 theorem-environments (uch-main Conditional; bulk-recovery
  Conditional; sdgauge ProvedElsewhere; ds-matter Conditional;
  gravity-leading ProvedElsewhere; higher-partial ProvedHere;
  chain-level class M Conditional; uch-ym Conditional;
  mellin-shadow ProvedHere; soft-hierarchy Conditional; celestial
  duals of Heisenberg and Virasoro ProvedHere). Source category
  pinned in the theorem head and clause (i): objects are the BV–HT
  ansatz pairs of Definition `def:uch-ht-twist` (perturbative local
  BV theories on ℝ⁴ with HT-twist data), not arbitrary 4d QFTs.

## Not completed / left open (with reasons)

- The chain-level proof of `conj:casimir-antighost-commutativity`
  and of `conj:coupled-sc-koszul`: these are genuine open
  mathematics (the audit's point); they are now named conjectures
  with their obstructions stated, not silently assumed.
- A first-principles computation of κ^Hodge(BP_k) (would close the
  open BP slot): not attempted; recorded as open in
  `prop:bp-kappa-conductor`.
- `preface_trimmed.tex` still contains four pre-healing Koszulity
  sentences; it is not \input by main.tex (verified) and is treated
  as an inactive draft.
- F-II.4's chain-level [Q,G̃] = T_Sug for affine KM (T1 at chain
  level) remains open and explicitly so in
  `rem:frontier-class-L-strict-chain-level` (unchanged; already
  honest).
- Environment balance verified on all nine heavily edited files
  (begin/end counts match); no full `make` run (session-end builds
  are user-opt-in per CLAUDE.md).
