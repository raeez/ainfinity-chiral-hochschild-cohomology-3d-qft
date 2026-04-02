# ADVERSARIAL AUDIT: FG-Shadow Conjecture Infrastructure
## ordered_associative_chiral_kd_frontier.tex, lines 43–448

### CRITICAL FINDING: The conjecture statement and proof strategy are **MISALIGNED**.

---

## Issue 1: "Exhaustive separated commutator filtration" is **NEVER FORMALLY DEFINED**

**Location:** Conjecture 1 (line 46–48)

**The problem:**
- Line 47: "Let $A$ be a filtered strongly admissible associative chiral algebra with exhaustive separated **commutator filtration** $F_{\mathrm{com}}^\bullet A$."
- The phrase "exhaustive separated commutator filtration" appears **only** in the conjecture hypothesis. 
- No formal definition of what this filtration **is** appears in the text before line 369.
- The conjecture is stated as a blank check: "assume there exists a filtration with property X, then Y follows."

**Where it IS defined (eventually):**
- Lines 369–374: "Specifically, define the \emph{commutator depth} of a field $a$ as the minimal $d$ such that..."
- Lines 390–403: Reformulation of "the commutator filtration on the bar complex" in terms of tensor factors and their depths.

**Status:** The conjecture statement (Conjecture 1, lines 45–63) **assumes the existence of a commutator filtration without defining what it is**. The definition comes **after** the conjecture, in Step 6 of the proof of Theorem 1 (line 195). This is a **logical inversion**: the conjecture cannot be stated precisely without a prior definition.

**Recommendation:** Define commutator depth and the commutator filtration on $A$ **before** Conjecture 1. The current path is: **Conjecture (undefined) → Remark → Construction → Proposition → Corollary → Proof (with definition embedded in Step 6)**.

---

## Issue 2: The "pole-order stratification" is **NOT identical** to the commutator filtration

**Location:** Construction (lines 87–109)

**The problem:**

The Construction claims (lines 92–93):
> "The commutator filtration $F_{\mathrm{com}}^\bullet A$ **induces** a stratification of the configuration factor $\FM_k(\mathbb{C})\times\Conf_k(\mathbb{R})$ by **pole order of the OPE**."

But then lines 369–403 **reveal a conflict**:
- The commutator filtration is **not** the pole-order filtration.
- Descendants $\partial^k e_I$ belong to the **same** commutator-filtration level as $e_I$ (line 364).
- The commutator filtration measures **iterated 0-th products (Lie brackets)**, not OPE pole order.
- Translation is an **inner derivation**, not a commutator operation (lines 365–377).

**What actually happens:**

The pole-order stratification in Construction (Definition of $S_p$, lines 95–103) reads:
$$S_p := \left\{ (z_1, \ldots, z_k; t_1 < \cdots < t_k) \,\big|\, \text{maximal pole order at any collision} \le p \right\}$$

This stratifies the **configuration space** (the domain of the bar complex) by **OPE pole orders**.

**But the commutator filtration** (lines 390–393) stratifies the bar complex **itself** by the sum of commutator depths of tensor factors:
$$F^p_{\mathrm{com}}\overline{B}^{\mathrm{ch,ord}}(\cA) := \bigoplus_k \left(\bigotimes_{i=1}^k F^{p_i}_{\mathrm{com}}\bar\cA\right) \otimes C_\bullet(\Conf_k^{\mathrm{ord}})$$
where $\sum_i p_i \le p$.

**Are these the same?**
- The Construction defines $S_p$ as a **subset of the configuration space** based on which fields can collide with which others.
- The bar complex filtration (lines 390–393) defines levels based on **how many times the 0-th product (Lie bracket) has been applied**.
- These **are not the same object**. One is a geometric stratification of the configuration space. The other is an algebraic filtration of a chain complex.

**How the proof claims they connect (line 84–85):**
> "using the pole-order stratification **induced by** the commutator filtration of Conjecture~\ref{conj:FG-shadow}"

But **what does "induced by" mean?** The text never explains the precise map. It jumps from:
1. The commutator filtration on $A$ (depth of a field)
2. To a stratification of $\FM_k(\mathbb{C})\times\Conf_k(\mathbb{R})$ (by maximal collision pole order)

**The missing piece:** How does the commutator depth on $A$ **constrain** which fields can have high pole orders?

**Hypothesis:** If $a \in F^d_{\mathrm{com}}A$ (depth $d$), then $a$ is a sum of $d$-fold iterated 0-th products of generators. Since $T_{(0)}\,a$ does increase pole order (line 355), maybe the claim is that iterating the Lie bracket $d$ times bounds the maximum pole order. But this is **not stated formally**.

**Status:** The bridge between the commutator filtration (on $A$) and the pole-order stratification (on $\FM_k \times \Conf_k$) is **logically undefined**. The construction assumes it exists but does not prove it.

---

## Issue 3: Hypothesis (a) is **trivially true and is NOT automatic**

**Location:** Conjecture 1, hypothesis (a) (line 49)

**The problem:**

Hypothesis (a) states:
> "$\gr_{\mathrm{com}}A$ is $E_\infty$-chiral"

The text claims (lines 106–108):
> "For the associated graded $\gr_{\mathrm{com}}A$ (which is $E_\infty$-chiral, i.e.\ has only simple poles), only the stratum $S_1$ contributes to the bar differential..."

**Is this automatic from pole-order killing?**

The text seems to suggest: "If we kill all poles ≥ 2 in the associated graded, we're left with only simple poles, which is $E_\infty$."

But:
- If $A$ has fields with pole orders up to $p$, then $\gr_{\mathrm{com}}A$ is the associated graded of the commutator filtration, which **kills commutators** (the Lie bracket), not poles.
- Hypothesis (a) **is not automatically true** from the pole-order interpretation.
- Rather, hypothesis (a) is an **assumption on the structure of** $\gr_{\mathrm{com}}A$: it says that the associated graded happens to be $E_\infty$-chiral.

**Example:** For Virasoro:
- The generator $T$ has depth 0 in the commutator filtration ($F^0_{\mathrm{com}}$).
- Its pole order is 4.
- The derivative $\partial T$ is also in $F^0_{\mathrm{com}}$ (translation is inner).
- But $\partial T$ has pole order 5 with $T$.
- In $\gr_{\mathrm{com}}$, we kill commutators, but $\partial T$ and $T$ still collide with pole order ≥ 2.
- So $\gr_{\mathrm{com}}\text{Vir}_c$ is **not** $E_\infty$-chiral unless we restrict further (e.g., to the commutative subalgebra).

**Status:** Hypothesis (a) is **not automatic**. It is a genuine constraint on the input algebra $A$. The text conflates "associated graded of commutator filtration" with "killing high poles" — these are different operations.

---

## Issue 4: Stratification-dependence on bar element is **NOT addressed**

**Location:** Construction (lines 95–103) and Proposition (lines 112–124)

**The problem:**

The stratum $S_p$ is defined as:
$$S_p := \left\{ (z_1, \ldots, z_k; t_1 < \cdots < t_k) \,\big|\, \text{maximal pole order at any collision} \le p \right\}$$

But "maximal pole order at any collision" **depends on the bar element**:
- A bar element is a tensor $e_{I_1} \otimes \cdots \otimes e_{I_k}$ of fields.
- The collision between $e_{I_i}$ and $e_{I_j}$ has pole order $N(e_{I_i}, e_{I_j})$.
- Different bar elements have different field types, so **different pole orders**.

**Example:** In Virasoro:
- The bar element $T \otimes T$ (two Virasoro generators) has collision pole order 4 between them.
- The bar element $\partial T \otimes T$ has collision pole order 5 between them.
- Are these in the same stratum $S_4$? The definition says "maximal pole order ≤ 4," so $T \otimes T \in S_4$, but $\partial T \otimes T \notin S_4$ (since $N(\partial T, T) = 5$).

**What the text actually means:**

Lines 136–137 clarify:
> "More precisely, the quotient $S_p/S_{p-1}$ carries chains supported on configurations whose maximal pole order is **exactly** $p$."

So $S_p$ is understood as a **subset of the geometric space** $\FM_k(\mathbb{C})\times\Conf_k(\mathbb{R})$, and the bar complex decomposes as:
$$\overline{B}^{\mathrm{ch,ord}}(A) = \bigoplus_k \left( \text{chains on } \FM_k(\mathbb{C})\times\Conf_k(\mathbb{R}) \right)$$

The stratification is **geometric** and **does not depend on which fields are placed at which points**. Rather, it depends on the configuration of points $(z_1, \ldots, z_k; t_1, \ldots, t_k)$ themselves and what pole orders **could** appear for any choice of fields.

**Status:** The stratification is on the **configuration space**, not on individual bar elements. This is correct but **could be much clearer**. The definition should say: "For a fixed bar element $e_{I_1} \otimes \cdots \otimes e_{I_k}$, the stratum $S_p$ consists of those configurations $(z_1, \ldots, z_k; t_1, \ldots, t_k)$ where maximal pole order across all pairwise collisions of fields in the element is ≤ $p$." Currently it reads like a definition of a fixed geometric object independent of the element.

---

## Issue 5: The "$E_1$ page identification" (Proposition 1) is **conceptually unclear**

**Location:** Proposition (lines 112–124)

**The claim:**
$$E_1^{p,q} \cong H^{p+q}_{\BM}(S_p/S_{p-1}) \cong \overline{B}^{FG}_{p+q}(\gr_{\mathrm{com}}A)$$

**The problem:**

1. **Left side:** $E_1^{p,q} = H^{p+q}_{\BM}(S_p/S_{p-1})$ is the Borel–Moore homology of the quotient stratum.
2. **Right side:** $\overline{B}^{FG}_{p+q}(\gr_{\mathrm{com}}A)$ is the Francis–Gaitsgory bar complex of the associated graded.

**Is this isomorphism **literally** the FG differential?**

The proof (lines 127–133) says:
> "On the associated graded, all OPE poles are simple, so the bar differential $d_{\mathrm{bar}}$ **restricts** to the simple-pole collision stratum $S_1$. The Borel–Moore chains $C_*^{\BM}(S_1)$ form precisely the $E_\infty$ configuration space underlying the Francis–Gaitsgory bar construction: the FG differential is the **residue map along** $\partial S_1$, i.e.\ the Stokes extraction of the simple-pole coefficient at each pairwise collision boundary."

**Interpretation:** 
- The $E_1$ page has only $p=1$ contributing (for $\gr_{\mathrm{com}}A$).
- The differential on $E_1$ is the boundary map of the Borel–Moore complex, which extracts residues at the boundary $\partial S_1$.
- This is the **definition** of the FG differential (residue at pairwise collisions).
- So they are the **same object**, not just isomorphic.

**But the text uses the word "residue map along $\partial S_1$"** — is this the same as the FG differential from the literature?

The reference is Francis–Gaitsgory, § 4.2. The current text does not cite what the FG differential actually is; it assumes the reader knows.

**Status:** The identification is **correct in spirit** but **relies on unstated equivalences** between:
- Borel–Moore boundary maps on $\partial S_1$
- Residue extraction at pairwise collisions
- The Francis–Gaitsgory differential from their paper

The text should be explicit: "The FG differential extracts the residue of the OPE at each first-order pole (pairwise collision boundary). This is exactly what the Borel–Moore boundary map $\delta: C_*^{\BM}(S_1) \to C_{*-1}^{\BM}(\partial S_1)$ computes."

---

## Issue 6: Convergence claim uses a **non-standard** criterion

**Location:** Corollary (lines 146–170)

**The claim:**
> "Strong convergence of the filtration (hypothesis (c)) gives convergence of the spectral sequence."

**The issue:**

A spectral sequence associated to a filtered complex $F^\bullet C_\bullet$ converges when the filtration is:
1. **Exhaustive:** $\bigcup_p F^p = C_\bullet$
2. **Separated (Hausdorff):** $\bigcap_p F^p = 0$
3. **Either:** (a) the complex is bounded below and the filtration is bounded above, **or** (b) the filtration is "strongly convergent" in the sense that the spectral sequence stabilizes.

**Hypothesis (c) of Conjecture 1 (line 51):**
> "the induced filtration on $\Barch(A)$ is complete and strongly convergent"

**What does "complete and strongly convergent" mean?**

In the context of pro-nilpotent filtered algebras:
- **Complete** means the inverse limit of finite filtration quotients recovers the original object.
- **Strongly convergent** for a spectral sequence means $E_\infty^{p,q} \Rightarrow \text{gr}^p(C_{\bullet+q})$ (the abutment is the associated graded).

**Does this guarantee convergence?**

Standard homological algebra says: **Yes**, if the filtration is exhaustive, separated, and the object is pro-nilpotent (or more generally, if the associated graded is bounded or the filtration has finite width).

**Status:** The corollary **asserts** convergence from hypothesis (c) but **does not prove that hypothesis (c) implies convergence**. It is an appeal to general principle, which is fine, but it should cite the theorem (e.g., "by completeness and Hausdorff property, standard spectral sequence convergence applies"). Currently the corollary reads like a tautology: "strong convergence gives convergence."

---

## Issue 7: The **duality line** in Corollary uses **unstated assumptions**

**Location:** Corollary (lines 156–169)

**The claim:**
> "On the Koszul locus, $\gr_{\mathrm{com}}(A^!)\simeq(\gr_{\mathrm{com}}A)^!_{FG}$."

**The proof (lines 164–169):**
> "the bar-cobar equivalence (Vol I, Theorem A) gives $A^!\simeq\overline{B}(A)^\vee$, and the associated graded of the dual is the dual of the associated graded (the filtration is exhaustive and separated by hypothesis (c)), so... $\gr_{\mathrm{com}}(A^!)\simeq ... =(\gr_{\mathrm{com}}A)^!_{FG}$."

**The issue:**

The claim "associated graded of the dual = dual of the associated graded" requires **finite-dimensionality of each graded piece** (Definition of strongly admissible, line 168 of core chapter).

If $\dim F^p_{\mathrm{com}}A/F^{p-1}_{\mathrm{com}}A = \infty$, then the dual might not preserve grading.

**Status:** The argument is **correct** given hypothesis (b) of the Definition of strongly admissible (line 168 of core: "each graded/conformal piece... is finite-dimensional"). But the Corollary proof **does not restate this assumption**, so a reader checking the logic would need to trace back to the definition of strongly admissible.

**Recommendation:** Quote hypothesis (b) explicitly in the proof.

---

## Issue 8: **Critical gap in Step 6** — the "derivative correction" is **introduced but not proven to work**

**Location:** Step 6 (lines 351–403)

**What happens:**

Step 5 shows that the **naive pole-order filtration fails** for Virasoro: the derivative $\partial T$ has pole order 5 with $T$, but $T$ has pole order 4 with itself, so $d_{\mathrm{bar}}(S_4) \not\subset S_4$.

Step 6 **redefines** the filtration to use "commutator depth" instead of "pole order":
- Line 369–374: Define commutator depth as iterated 0-th products.
- Lines 375–377: Claim that translation acts on each depth level.
- Lines 379–403: Reformulate the bar complex filtration by tensor factor depths.

**The problem:**

The reformulation (lines 390–403) claims:
> "the commutator filtration on the bar complex... respects this filtration because: (i) The $n=0$ term $a_{(0)}\,b$ has commutator depth $\le\max(d(a),d(b))+1$ but goes into the next filtration level... (ii) The $n\ge 1$ terms $a_{(n)}\,b$ preserve commutator depth: $d(a_{(n)}\,b)\le d(a)+d(b)$ for $n\ge 1$."

**But this is the definition of "the bar differential respects the filtration," not a proof of it.**

For this to hold, the text must prove:
1. The 0-th product $a_{(0)}\,b$ (Lie bracket) is a commutator, so it **increases depth by 1**.
2. The $n\ge 1$ residues $a_{(n)}\,b$ (for $n\ge 1$) **do not introduce new Lie brackets**, so they preserve depth.

**Where is this proof?**

- For (1): The text says "it is the Lie bracket, the defining operation that the commutator filtration measures" (line 399). This is true but needs justification: **why is the Lie bracket in the definition the same as the $n=0$ product?**
- For (2): The text says $a_{(n)}\,b$ with $n\ge 1$ is "a 'higher symmetry' extraction, not a commutator" (lines 383–384). But **why** is extracting higher OPE residues not a commutator? The higher residues can still be sums of Lie brackets applied to lower modes.

**Status:** Step 6 **asserts the key commutator-depth properties without proof**. It is heuristically correct — the idea that Lie brackets form a graded space while residues are "other operations" is sound — but it is **not rigorously justified in the text**.

---

## Issue 9: Step 7 **does not close the gap** for general Virasoro weights

**Location:** Step 7 (lines 414–432)

**What it does:**

For affine Kac–Moody $\hat\mathfrak{g}_k$ with weight-1 generators:
- The bar differential produces $J^a{}_{(0)}\,J^b = f^{ab}_c J^c$ (a Lie bracket) and $J^a{}_{(1)}\,J^b = k\kappa^{ab}$ (a scalar).
- The text claims: "The $n=0$ output is the Lie bracket with $N(J^c,J^d)\le 2$ for all $c,d$, so $d_{\mathrm{bar}}(S_2)\subset S_2$" (lines 420–421).

**The issue:**

This works because:
- Affine generators have weight 1.
- Their pole orders are at most 2 (since $N(a,b) \le h_a + h_b$).
- The Lie bracket doesn't increase pole order (since it's a generator, not a derivative).

But for **Virasoro** (weight 2):
- The bar differential produces $T_{(0)}\,T = \partial T$, which is a **derivative**, not a new generator.
- The derivative **increases pole order by 1** (line 355).

**The text then claims (lines 435–448):**
> "The pole non-increase estimate... holds in its stated form (with $N$ the OPE pole order) for the affine lineage and for any vertex algebra whose generators all have conformal weight 1. For generators of weight ≥ 2, the 0-th product can produce descendants with higher pole orders. However, the **commutator filtration** (which is the filtration relevant to Conjecture~\ref{conj:FG-shadow}) is always respected."

**So the claim is:**
- For weight ≥ 2, the **bare pole-order filtration fails**.
- But the **commutator-depth filtration** works.

**The problem:**

The text **asserts** this but **does not prove** that the commutator-depth filtration works for Virasoro. It proves:
- For affine (weight 1): the pole-order filtration works.
- For all weights: the commutator-depth filtration is defined and its $n=0$ term increases depth while $n\ge 1$ terms preserve depth (lines 396–402).

But **where is the proof that Steps 1–5 of Theorem 1 still hold** when we switch to commutator depth?

The Theorem statement (lines 195–215) is titled "Pole non-increase for the bar differential" and its conclusion (lines 206–211) is stated in terms of pole order $N$, not commutator depth. So **the theorem proves pole non-increase, not commutator-depth preservation**.

**Status:** **Critical gap:** Step 8 (Conclusion) claims "the commutator filtration is always respected" (line 440) and thus "completes the reduction of Conjecture~\ref{conj:FG-shadow}" (lines 444–445), but **the commutator-depth version is never formally proved**. The main Theorem 1 proves pole order is non-increasing, which **fails for Virasoro**. The claimed fix (use commutator depth instead) is **conceptually correct but mathematically unjustified**.

---

## Summary of Issues

| Issue | Severity | Type | Location |
|-------|----------|------|----------|
| 1. "Exhaustive separated commutator filtration" undefined before use | **CRITICAL** | Logical inversion | Conjecture 1 (line 45), definition at line 369 |
| 2. "Pole-order stratification induced by commutator filtration" has no bridge | **CRITICAL** | Gap in logic | Construction (lines 84–93) |
| 3. Hypothesis (a) claimed automatic but is nontrivial | **MODERATE** | Imprecision | Conjecture 1 (line 49), lines 106–108 |
| 4. Stratification definition appears independent of bar element | **MODERATE** | Clarity issue | Construction (lines 95–103) |
| 5. $E_1$ identification relies on FG differential literature | **MODERATE** | Incomplete reference | Proposition 1 (lines 112–133) |
| 6. Convergence asserted from "strong convergence" (tautology) | **LOW** | Expository | Corollary (lines 160–162) |
| 7. Duality proof omits finite-dimensionality hypothesis | **LOW** | Technical omission | Corollary (lines 164–169) |
| 8. Commutator-depth properties asserted without proof | **CRITICAL** | Gap in logic | Step 6 (lines 379–403) |
| 9. Main theorem doesn't prove commutator-depth version | **CRITICAL** | Mismatch between claim and proof | Theorem 1 + Step 8 |

---

## Recommendations for Repair

### Priority 1: Reorder the exposition
1. Define commutator depth (lines 369–374) **before** Conjecture 1.
2. Add a Definition section: "Definition (Exhaustive separated commutator filtration): A filtration $F^\bullet A$ on an associative chiral algebra $A$ by commutator depth, meaning..."
3. Then state Conjecture 1 with a well-defined object.

### Priority 2: Justify the bridge
1. After Conjecture 1, prove (or cite as Lemma): "If $A$ is equipped with a commutator filtration, then every field $a \in F^d_{\mathrm{com}}A$ satisfies [some bound on its pole orders with all other fields]."
2. Use this Lemma to justify Construction: "Thus the commutator filtration on $A$ induces a stratification on $\FM_k(\mathbb{C})\times\Conf_k(\mathbb{R})$..."

### Priority 3: Separate pole-order from commutator-depth
1. Rename Theorem 1: "Theorem 1 (Pole non-increase for weight-1 generators):" or "Theorem 1 (Pole non-increase for affine algebras)."
2. After Step 5 (Virasoro failure), add a **new Theorem 2**: "Theorem 2 (Commutator-depth preservation): For any vertex algebra, the commutator-filtered bar complex satisfies $d_{\mathrm{bar}}(F^p_{\mathrm{com}}\Barch) \subset F^{p+1}_{\mathrm{com}}\Barch$."
3. Then use Theorem 2 in the conclusion.

### Priority 4: Clarify the FG differential
1. Add a remark after Proposition 1: "The Francis–Gaitsgory differential is defined in [FG] as the residue extraction at each pairwise collision in the $E_\infty$ configuration space. The Borel–Moore boundary map $\delta: C_*^{\mathrm{BM}}(S_1) \to C_{*-1}^{\mathrm{BM}}(\partial S_1)$ computes exactly this, via the stokes formula on the boundary stratum."

### Priority 5: Explicit hypothesis tracing
1. At the start of Step 7, add: "Hypothesis (b) of Definition (strongly admissible) ensures each graded piece is finite-dimensional. Thus... [continue]."
2. After the Virasoro section, add: "We have now verified the pole-order version for weight-1 (Step 7) and shown the failure at weight 2 (Step 5). To complete the proof for all weights, we adopt the commutator-depth filtration (Step 6) and prove the following:"

---

## Conclusion

The FG-shadow conjecture infrastructure is **conceptually sound but pedagogically inverted and logically gapped**:

1. The main hypothesis is undefined at first use.
2. The bridge between the hypothesis and the construction is not explicit.
3. The proof strategy (Theorem 1) proves the wrong thing (pole non-increase, which fails).
4. The fix (commutator-depth version) is asserted without a separate proof.
5. Hypothesis (a) of the conjecture is presented as automatic when it is nontrivial.

**Fixes are substantial but straightforward:** reorder definitions, separate the affine and general cases, and prove the commutator-depth version as a standalone theorem.

The **mathematical content is likely correct**. The issue is **presentation order and logical gaps**, not mathematical error.
