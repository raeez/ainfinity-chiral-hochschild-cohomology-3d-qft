# Route A: Modular Swiss-Cheese Operad — THE LOCAL APPROXIMATION

## Scope
This route creates the local-model chapter AND rewrites passages in line-operators.tex and bar-cobar-review.tex to be honest about the operadic scope. Route A is honest: operadic structure captures the formal-disc data and misses the global D-module structure that produces curvature.

## Read First
1. /Users/raeez/chiral-bar-cobar-vol2/CLAUDE.md
2. /Users/raeez/chiral-bar-cobar-vol2/.claude/specs/master.md (THE SIX-LAYER FRAMEWORK)
3. /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/line-operators.tex lines 1-275 (current homotopy-Koszulity proof + consequences)
4. /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/bar-cobar-review.tex lines 1100-1900 (Swiss-cheese framing in bar-cobar review)
5. /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/modular_pva_quantization_core.tex lines 425-850 (modular bar, D0/D1)
6. /Users/raeez/chiral-bar-cobar-vol2/chapters/examples/rosetta_stone.tex lines 90-150 (Heisenberg Swiss-cheese)
7. /Users/raeez/chiral-bar-cobar/chapters/theory/higher_genus_modular_koszul.tex lines 331-438 (Vol I modular Koszul)
8. /Users/raeez/chiral-bar-cobar/chapters/theory/bar_cobar_adjunction_curved.tex lines 5910-5958 (Vol I FCom)
9. /Users/raeez/chiral-bar-cobar/chapters/theory/higher_genus_foundations.tex lines 5918-5967 (Vol I Feynman involution)

## Part 1: New Chapter

### Output File
/Users/raeez/chiral-bar-cobar-vol2/chapters/theory/modular_swiss_cheese_operad.tex

### Mathematical Content

#### Opening: Honest framing

The chapter opens by stating its scope explicitly: this is the LOCAL STORY. The topological operad SC^{ch,top} models what happens in a formal disc around a collision point, where X ~ C. It does not see global information: D-module monodromy, periods of Sigma_g, the factorization pattern on Ran(Sigma_g). The global story requires the factorization framework of Chapter ref{ch:factorization-swiss-cheese}.

The local story is still essential: it explains WHY the corrected holomorphic differential D^(g)^2 = 0 is flat (the Fay trisecant identity is a LOCAL identity on configuration spaces), while the curvature d_fib^2 = kappa * omega_g is a GLOBAL phenomenon. The corrected holomorphic model lives in the operadic world; the curved geometric model lives in the factorization world; the modular bar complex mediates between them.

**Remark (What the local model sees and what it misses).**
- SEES: collision data (OPE coefficients), the Arnold relation and its Fay trisecant generalization, the Koszul dual A^!, the Swiss-cheese directionality, the formal-disc propagator
- MISSES: D-module monodromy around cycles, period corrections from H^1(Sigma_g), the Arakelov propagator's non-holomorphic correction, the curvature kappa * omega_g
- The corrected holomorphic differential D^(g) uses ONLY local data (the holomorphic propagator, which is the prime form) and is therefore visible to the local operad. The curved differential d_fib uses GLOBAL data (the Arakelov propagator, which depends on the period matrix of Sigma_g) and is therefore invisible to the local operad.

#### Section 1: The partially modular Swiss-cheese operad

**Definition (Partially modular two-colored operad).** As in the original Route A spec, but with explicit acknowledgment that this is the local model extracted from the factorization structure by formal completion.

**Definition (SC^{ch,top}_mod).** The modular chiral Swiss-cheese operad:
- Closed: C_*(FM_k(Sigma_g)) for all g, k — these are LOCAL models (FM compactifications of configuration spaces)
- Open: E_1(m) — genus 0 only (no genus for the interval)
- Mixed: C_*(FM_{k|m}(Sigma_g, partial))
- Directionality: no open-to-closed

**Remark (Extraction from factorization).** SC^{ch,top}_mod is extracted from the BD factorization Swiss-cheese algebra (Definition in factorization_swiss_cheese.tex) by:
1. Restricting the factorizable D-module to formal neighborhoods of the small diagonal
2. Replacing D-module data by formal completion data (vector spaces + differentials)
3. Identifying the resulting local data with chains on FM compactifications
This extraction is a functor from factorization algebras to operadic algebras, and it is faithful at genus 0 on P^1 but not at genus g >= 1.

#### Section 2: Modular homotopy-Koszulity

**Definition (Modular homotopy-Koszulity for partially modular operads).** As before.

**Theorem (Modular homotopy-Koszulity of SC^{ch,top}_mod; ProvedHere).**

Four-step proof:
1. Genus-0 homotopy-Koszulity (Livernet + Kontsevich + Berger-Moerdijk — existing proof)
2. Closed-color modular Koszulity (formality at each genus via Fay trisecant + Feynman involutivity)
3. Mixed-operation compatibility (Swiss-cheese directionality ensures D_1 acts only on closed color)
4. Assembly via spectral sequence convergence

**Remark (What modular homotopy-Koszulity achieves and what it doesn't).**
- ACHIEVES: bar-cobar equivalence for the flat models (D_0^2 = 0 and D^(g)^2 = 0)
- DOES NOT ACHIEVE: the curved bar-cobar equivalence (d_fib^2 = kappa * omega_g) — this requires the factorization input (global D-module data) from Chapter ref{ch:factorization-swiss-cheese}
- The derived-coderived equivalence at genus g >= 1 is NOT a consequence of modular homotopy-Koszulity alone. It requires: (1) modular homotopy-Koszulity (this chapter) for the flat model, (2) factorization Koszulity (Chapter ref{ch:factorization-swiss-cheese}) for the curved model, (3) the chiral Riemann-Hilbert correspondence to connect them.

#### Section 3: Consequences and scope

**Corollary (Flat-model equivalence; ProvedHere).** For any logarithmic SC^{ch,top}_mod-algebra A, the flat modular bar complex B_mod(A) with differential D = D_0 + D_1 has D^2 = 0, and the modular bar-cobar adjunction is a Quillen equivalence. This is purely operadic — it lives in the derived category.

**Corollary (Spectral sequence convergence; ProvedHere).** The genus-filtration spectral sequence converges for any SC^{ch,top}_mod-algebra. The curvature kappa * omega_g is the d_1 obstruction.

**Non-corollary (Curved-model equivalence).** The curved bar-cobar equivalence (in the coderived category) does NOT follow from operadic modular homotopy-Koszulity alone. It requires factorization input. See Chapter ref{ch:factorization-swiss-cheese}, Theorem ref{thm:fact-sc-koszul-duality}.

**Remark (Relation to the three models).** Of the three models in Remark rem:three-models:
- Model 1 (flat associated graded, D_0^2 = 0): fully captured by the local operad
- Model 2 (corrected holomorphic, D^(g)^2 = 0): fully captured by the local operad (Fay is local)
- Model 3 (curved geometric, d_fib^2 = kappa * omega_g): NOT captured by the local operad (curvature is global)
The derived-coderived equivalence connecting Models 1-2 (derived) with Model 3 (coderived) requires the factorization framework of Chapter ref{ch:factorization-swiss-cheese}.

## Part 2: Rewrites of Existing Chapters

### bar-cobar-review.tex Rewrites

**Line 1159: "form an algebra over the Swiss-cheese operad"**
REPLACE "algebra" with "coalgebra"

**Line 1178: "is an algebra over the Swiss-cheese operad"**
REPLACE "algebra" with "coalgebra"

**Lines 1860-1872: "SC is homotopy-Koszul...property of the operad"**
ADD after "not of individual algebras": "It is a property of the \emph{local model} — the formal completion of the factorization structure at collision points. The global factorization Koszulity (Chapter~\ref{ch:factorization-swiss-cheese}) is the primary datum from which operadic homotopy-Koszulity is extracted."

### rosetta_stone.tex Rewrites

**Lines 99, 126, 138: "Swiss-cheese algebra" (Heisenberg)**
REPLACE all instances of "Swiss-cheese algebra" referring to the bar complex with "Swiss-cheese coalgebra". Specifically:
- Line 99 subsection header: "Swiss-cheese coalgebra"
- Line 126 theorem title: "Swiss-cheese coalgebra"
- Line 138: "is a coalgebra over $\SCchtop$"

### line-operators.tex Rewrite

**After Theorem thm:homotopy-Koszul (after line 54), add a remark:**

```
\begin{remark}[Scope: local model]
\label{rem:hkoszul-scope}
Theorem~\ref{thm:homotopy-Koszul} is a theorem about the local
model $\SCchtop$ --- the formal completion of the factorization
structure at collision points.  It captures the Koszulity of the
collision data (the OPE coefficients and their Arnold relations)
but does not see global information: D-module monodromy around
cycles, period corrections from $H^1(\Sigma_g)$, or the
Arakelov propagator's non-holomorphic correction that produces
the curvature $\dfib^{\,2} = \kappa \cdot \omega_g$.  The
global factorization Koszulity, from which this local property
is derived, is the content of
Chapter~\textup{\ref{ch:factorization-swiss-cheese}}.
\end{remark}
```

## main.tex Edit
Add after factorization_swiss_cheese (Part I):
```latex
\input{chapters/theory/modular_swiss_cheese_operad}
```

## CLAUDE.md Edit
Add to Theory file map:
```
- modular_swiss_cheese_operad: Part I (local operadic approximation, modular SC^{ch,top}_mod, honest about scope)
```

## Definition of Done
- [ ] chapters/theory/modular_swiss_cheese_operad.tex exists and compiles
- [ ] Opens with honest framing: this is the LOCAL story
- [ ] Contains: extraction from factorization, partially modular definition, modular homotopy-Koszulity theorem, scope limitations
- [ ] bar-cobar-review.tex algebra->coalgebra fixes (lines 1159, 1178)
- [ ] bar-cobar-review.tex homotopy-Koszulity qualified as local (lines 1860-1872)
- [ ] rosetta_stone.tex algebra->coalgebra fixes (lines 99, 126, 138)
- [ ] line-operators.tex scope remark added after thm:homotopy-Koszul
- [ ] main.tex \input added
- [ ] CLAUDE.md updated
- [ ] Full build passes
