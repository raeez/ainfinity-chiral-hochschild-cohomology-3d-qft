# Route B: Factorization Swiss-Cheese — THE PRIMARY TREATMENT

## Scope
This is the foundational chapter AND the rewrite of foundations.tex and introduction.tex to make factorization primary throughout. This route creates the new chapter AND corrects the hierarchy in existing chapters.

## Read First
1. /Users/raeez/chiral-bar-cobar-vol2/CLAUDE.md
2. /Users/raeez/chiral-bar-cobar-vol2/.claude/specs/master.md (THE SIX-LAYER FRAMEWORK — read this completely)
3. /Users/raeez/chiral-bar-cobar-vol2/chapters/theory/foundations.tex (FULL FILE — understand current operadic framing)
4. /Users/raeez/chiral-bar-cobar-vol2/chapters/theory/introduction.tex (FULL FILE — understand current framing)
5. /Users/raeez/chiral-bar-cobar-vol2/chapters/theory/axioms.tex lines 1-120 (operadic formulation remark)
6. /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/dg_shifted_factorization_bridge.tex (existing factorization chapter)
7. /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/line-operators.tex lines 1-225 (current homotopy-Koszulity)
8. /Users/raeez/chiral-bar-cobar/chapters/theory/configuration_spaces.tex lines 2393-2443 (Vol I D^oc_A)
9. /Users/raeez/chiral-bar-cobar/appendices/coderived_models.tex lines 1-100 (Vol I Positselski)

## Part 1: New Chapter

### Output File
/Users/raeez/chiral-bar-cobar-vol2/chapters/theory/factorization_swiss_cheese.tex

### Mathematical Content

#### Section 1: The factorization substrate

**Opening paragraph.** State the ground truth directly: a chiral algebra is not an algebra over an operad. It is a factorizable D-module on the Ran space. The operad is extracted from the factorization structure by formal completion at collision points — it is the local shadow, not the governing structure. At genus 0 on P^1, local ~ global. At genus g >= 1, local != global, and the divergence IS the curvature.

**Definition (Ran space and factorization pattern).** For a smooth curve X:
- Ran(X) = colim_{S in fSet^{surj}} X^S — the space of finite non-empty subsets
- The factorization structure: for disjoint S_1, S_2 in Ran(X), the inclusion U_{S_1} x U_{S_2} -> Ran(X) induces the factorization map
- A factorizable D-module F on Ran(X): a D-module with factorization isomorphism F|_{U_I} ~ box_{i in I} F on complements of diagonals

**Definition (Chiral operations from factorization).** The chiral operations are NOT maps between fibers at points. They are:
- The restriction F^{box n}|_{Delta_{ij}} -> F^{box (n-1)} to the diagonal Delta_{ij} = {z_i = z_j} in X^n
- These are maps between D-modules supported on different diagonals
- The composition law comes from iterating diagonal restrictions — this is the factorization structure, not an operad composition

**Remark (The local shadow).** The topological operad SC^{ch,top} arises by:
1. Restricting to a formal neighborhood of a point x in X (where X ~ C)
2. Replacing the D-module structure by its formal completion (which is just a vector space with a differential)
3. Identifying the formal neighborhood of the small diagonal in X^n with FM_n(C)
This extraction loses global information: monodromy around cycles, D-module connections, period integrals. At genus 0 on P^1, the loss is trivial (P^1 - pt = C). At genus g >= 1, the loss is the curvature kappa * omega_g.

#### Section 2: BD factorization Swiss-cheese

**Definition (BD factorization Swiss-cheese algebra on Sigma_g x R).** Full definition with:
- Closed sector: factorizable D-module F_cl on Ran(Sigma_g)
- Open sector: factorization algebra F_op on Ran(R) ~ E_1-algebra
- Mixed sector: factorization module structure for configurations with closed points on Sigma_g and open points on boundary component of Sigma_g x R_{>= 0}
- Directionality: bulk -> boundary only (no open-to-closed maps)

**The genus tower from families of factorizable D-modules.** For varying g, the collection {F_cl^{(g)}} defines a family over M-bar_g. The D-module structure on Ran(Sigma_g) varies with the complex structure of Sigma_g. The curvature d_fib^2 = kappa * omega_g arises from the monodromy of the D-module connection around B-cycles of Sigma_g. This is a theorem about D-modules, not about topological operads.

**Bar complex as factorization coalgebra.** B(A) is the factorization coalgebra on Ran(X) Koszul-dual to the chiral algebra A:
- Bar differential = chiral operations dualized via factorization structure
- Bar coproduct = factorization isomorphism dualized
- At genus 0: d^2 = 0 (honest dg coalgebra)
- At genus g >= 1: d_fib^2 = kappa * omega_g (curved factorization coalgebra)

#### Section 3: CG factorization Swiss-cheese

**Definition (CG prefactorization algebra on Sigma_g x R).** Full definition with:
- Holomorphic factorization in Sigma_g direction
- Locally constant factorization in R direction
- BV quantization origin: action functional S = S_hol + S_top + S_mix
- Cosheaf descent making it a genuine factorization algebra

**Factorization homology.** integral_{Sigma_g x R} F = global sections. For varying g, produces the genus-g bar complex. The genus tower from dependence on complex structure.

**CG Koszul duality.** Ayala-Francis framework: Koszul duality for locally constant factorization algebras on R^n exchanges E_n-algebras and E_n-coalgebras. For HT product Sigma_g x R:
- R-direction: E_1 Koszul duality
- Sigma_g-direction: chiral Koszul duality
- Combined: Swiss-cheese Koszul duality

#### Section 4: Equivalence of BD and CG

**Theorem (BD-CG equivalence; ProvedHere conditional on standard comparison results).** State with precise adjectives:
- Characteristic 0
- X smooth projective
- Holomorphically translation-invariant on X, locally constant on R
- Completeness: ind-coherent on BD side, cosheaf topology on CG side

Proof through three layers:
(i) Closed (chiral): BD Thm 4.5.4 / Francis-Gaitsgory
(ii) Open (associative): Lurie HA 5.4.5.9
(iii) Mixed (Swiss-cheese): comparison of D-module and cosheaf descriptions of mixed Ran space

#### Section 5: Factorization Koszul duality and derived-coderived equivalence

**Theorem (Factorization Swiss-cheese Koszul duality; ProvedHere).** For a factorization Swiss-cheese algebra F on Sigma_g x R:
(i) Bar construction is factorization coalgebra
(ii) Cobar recovers F at genus 0 (by local homotopy-Koszulity)
(iii) Feynman transform recovers F at all genera (by Getzler-Kapranov involutivity on closed color)
(iv) Derived-coderived equivalence at all genera

**Corollary (Three models from factorization).** The three chain-level models are three gauges of the factorization structure:
- Flat = algebraic gauge (formal completion at collision points)
- Corrected holomorphic = flat gauge (universal cover)
- Curved = single-valued gauge (Sigma_g itself)
Reference Remark rem:three-models, which should be updated to cite this corollary.

#### Section 6: Aside on Latyntsev factorization quantum groups

**Remark (Connection to Latyntsev).** Develop:
- Latyntsev's Phi_q(g) as factorization algebra on Ran(P^1) -> quantum group representations via factorization homology
- Our factorization Swiss-cheese algebra -> dg-shifted Yangian (open-color projection) and chiral Hochschild (closed-color projection)
- The bridge: Latyntsev works on closed color, producing quantum groups via Ran(X). Our Yangian works on open color, producing shifted Yangians via E_1. The Swiss-cheese mixed structure mediates: R-matrix from bulk-boundary composition.
- The spectral R-matrix connection: Latyntsev's quantum R-matrix and our r(z) agree on the evaluation locus (DK-0 shadow).
- What's new: the two-colored factorization framework places both as projections of a single object.

## Part 2: Rewrites of Existing Chapters

### foundations.tex Rewrites

**Line 27: "operadic fingerprint"**
REPLACE: "This product is the operadic fingerprint of a three-dimensional holomorphic--topological (HT) field theory"
WITH: "This product is the local algebraic shadow of the factorization structure of a three-dimensional holomorphic--topological (HT) field theory"

**Lines 30-31: "governing structure is the two-colored operad"**
REPLACE: "The governing structure is the two-colored Swiss-cheese operad $\SCchtop$ with operation spaces $\FM_k(\C) \times E_1(m)$."
WITH: "The local model for this factorization structure is the two-colored Swiss-cheese operad $\SCchtop$ with operation spaces $\FM_k(\C) \times E_1(m)$, extracted from the factorization pattern by formal completion at collision points."

**Lines 364-366: "homotopy-Koszulity...operadic input"**
REPLACE: "The homotopy-Koszulity of $\SCchtop$ (Theorem~\ref{thm:homotopy-Koszul}) is the operadic input to"
WITH: "The homotopy-Koszulity of the local operadic model $\SCchtop$ (Theorem~\ref{thm:homotopy-Koszul}), reflecting the Koszulity of the underlying factorization structure (Chapter~\ref{ch:factorization-swiss-cheese}), provides the local input to"

**Lines 394-395: "SC operad is homotopy-Koszul"**
ADD qualifier after sentence: "(This is a property of the local model; the global factorization structure is the primary datum from which it is derived.)"

### introduction.tex Rewrites

**Lines 87-91: "bar complex is SCchtop-algebra"**
REPLACE "algebra" with "coalgebra" on every instance where the bar complex is called an "SC^{ch,top}-algebra." The bar complex is an SC^{ch,top}-COALGEBRA (it has a coproduct, not a product). The Koszul dual A^! carries the algebra structure.

**Line 158-159: "bar complex is a Swiss-cheese algebra"**
REPLACE: "the bar complex is a Swiss-cheese algebra, and a formal one"
WITH: "the bar complex is a Swiss-cheese coalgebra (the local model of a factorization coalgebra on $\operatorname{Ran}(X)$), and a formal one"

**Line 350: "governed by...operad"**
REPLACE: "The structures above are governed by a single two-coloured operad $\SCchtop$"
WITH: "The structures above are locally modelled by a single two-coloured operad $\SCchtop$, extracted from the factorization structure on $\operatorname{Ran}(X) \times \operatorname{Ran}(\R)$"

**Lines 281-285: "MC element lives at genus 0"**
ADD at end of paragraph: "At genus~$0$, the local operadic model suffices; the global family structure over $\overline{\cM}_g$ requires the factorization-algebraic framework of Chapter~\ref{ch:factorization-swiss-cheese}."

### axioms.tex Rewrite

**Lines 59-101: Operadic formulation remark**
ADD at end of remark (before \end{remark}): "This operadic description is a local shadow of the factorization-categorical formulation (Chapter~\ref{ch:factorization-swiss-cheese}), obtained by evaluating the factorization structure at formal completions of collision points. At genus~$0$ on $\mathbb{P}^1$, the two descriptions coincide; at genus~$g \geq 1$, the factorization formulation carries the D-module monodromy data that produces the curvature $\dfib^{\,2} = \kappa \cdot \omega_g$, which is invisible to the local operadic description."

## main.tex Edit
Add after line 807 (Part I, after raviolo-restriction):
```latex
\input{chapters/theory/factorization_swiss_cheese}
```

## CLAUDE.md Edits
1. File map: Add "- factorization_swiss_cheese: Part I (BD and CG factorization, equivalence, Koszul duality, Latyntsev)"
2. Add to Critical Pitfalls: "- **Factorization is primary, operad is local shadow**: A chiral algebra is a factorizable D-module on Ran(X). The SC^{ch,top} operad is the formal completion at collision points — the local model. At genus 0 on P^1, local = global. At genus g >= 1, local != global; the curvature is the global excess. NEVER present the operad as the 'governing structure' — it is the local approximation."

## Definition of Done
- [ ] New chapter chapters/theory/factorization_swiss_cheese.tex exists and compiles
- [ ] Contains Sections 1-6 as described above
- [ ] All internal \ref targets resolve
- [ ] foundations.tex lines 27, 30-31, 364-366, 394-395 rewritten
- [ ] introduction.tex lines 87-91, 158-159, 281-285, 350 rewritten (algebra->coalgebra + local shadow qualifier)
- [ ] axioms.tex lines 59-101 qualified
- [ ] main.tex \input added
- [ ] CLAUDE.md updated (file map + pitfall)
- [ ] Full build passes: pkill -9 -f pdflatex 2>/dev/null || true; sleep 2; make
