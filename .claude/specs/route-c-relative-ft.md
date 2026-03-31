# Route C: Relative Feynman Transform — THE ALGEBRAIC SKELETON

## Scope
This route creates the algebraic-skeleton chapter AND updates the modular PVA quantization chapter, the preface, main.tex CLAUDE.md, and Remark rem:three-models to reflect the full picture. Route C is the connective tissue: it names the algebraic structure that the factorization approach (Route B) and the operadic approach (Route A) both produce.

## Read First
1. /Users/raeez/chiral-bar-cobar-vol2/CLAUDE.md
2. /Users/raeez/chiral-bar-cobar-vol2/.claude/specs/master.md (THE SIX-LAYER FRAMEWORK)
3. /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/modular_pva_quantization_core.tex lines 425-900 (modular bar coalgebra, D0/D1, genus completion, three-models remark — READ CAREFULLY, this is the chapter you're extending)
4. /Users/raeez/chiral-bar-cobar-vol2/chapters/connections/line-operators.tex lines 1-275
5. /Users/raeez/chiral-bar-cobar-vol2/chapters/theory/foundations.tex lines 2195-2260 (genus tower as coderived)
6. /Users/raeez/chiral-bar-cobar-vol2/chapters/frame/preface.tex lines 1355-1400
7. /Users/raeez/chiral-bar-cobar/chapters/theory/configuration_spaces.tex lines 2393-2443 (Vol I D^oc_A)
8. /Users/raeez/chiral-bar-cobar/chapters/theory/bar_cobar_adjunction_curved.tex lines 5910-5958
9. /Users/raeez/chiral-bar-cobar/chapters/theory/higher_genus_foundations.tex lines 5918-5967
10. /Users/raeez/chiral-bar-cobar/chapters/theory/higher_genus_modular_koszul.tex lines 10490-10540

## Part 1: New Chapter

### Output File
/Users/raeez/chiral-bar-cobar-vol2/chapters/connections/relative_feynman_transform.tex

### Mathematical Content

#### Opening

The modular bar complex B_mod is already an algebra over a relative Feynman transform. This chapter names the structure, proves its involutivity, and connects it to the factorization (Chapter ref{ch:factorization-swiss-cheese}) and operadic (Chapter ref{ch:modular-sc-operad}) approaches.

The relative Feynman transform stands intermediate: it is the algebraic skeleton shared by the global (factorization) and local (operadic) viewpoints. The factorization approach produces it via families of factorizable D-modules over M-bar_g. The operadic approach produces it via the modular bar construction on SC^{ch,top}_mod-coalgebras. The relative Feynman transform is the abstract structure that both flesh out.

#### Section 1: The relative Feynman transform

**Definition (Relative modular extension of a two-colored operad).** As before: modular operad M on closed color, operad P on open color, mixed operations, compatibility condition at genus 0.

**Definition (Relative Feynman transform FT_{M/P}).** The Feynman transform of M taken relative to P:
- FT_{M/P}(A) is a coalgebra over P^! with modular coalgebra structure on closed color
- Total differential D = D_P + D_M where D_P = D_0 (P-coalgebra structure, genus-preserving) and D_M = D_1 (M-modular structure, genus-raising)
- Bicomplex: D_P^2 = 0, D_M^2 = 0, D_P D_M + D_M D_P = 0

**Key observation.** FT_{Com_mod / SC^{ch,top}} = B_mod (the existing modular bar coalgebra):
- D_P = D_0 = D_int + D_sep
- D_M = D_1 = D_nsep
This is Proposition prop:D0D1 reinterpreted.

#### Section 2: Recognition theorem

**Theorem (Recognition; ProvedHere).** B_mod(A) is naturally an algebra over FT_{Com_mod / SC^{ch,top}}.

Proof: reinterpretation of Theorem thm:modular-bar, Proposition prop:D0D1, Theorem thm:genus-completion.

**Remark (Local vs global incarnations).** The relative Feynman transform has two incarnations:
- LOCAL (operadic): extracted from SC^{ch,top}_mod by formal completion — this produces the flat models (D_0^2 = 0 and D^(g)^2 = 0)
- GLOBAL (factorization): extracted from the factorization Swiss-cheese algebra by the Ran space structure — this produces both flat and curved models
The relative Feynman transform is the abstract structure common to both. It sees the bicomplex (D_0, D_1) and the genus spectral sequence, but does not by itself determine whether the fiberwise model is flat (operadic) or curved (factorization).

#### Section 3: Homotopy-involutivity

**Definition (Relative homotopy-involutivity).** FT_{M/P}(FT_{M/P}(A)) ~ A.

**Theorem (Homotopy-involutivity of FT_{Com_mod / SC^{ch,top}}; ProvedHere).**

Proof by three steps:
1. P-involutivity at genus 0 (from operadic homotopy-Koszulity, Theorem thm:homotopy-Koszul)
2. M-involutivity on closed color (from Feynman involutivity, Vol I Theorem thm:feynman-involution)
3. Assembly via spectral sequence convergence (completeness + pronilpotence from Theorem thm:modular-bar(iii))

**Remark (What involutivity achieves).** Relative homotopy-involutivity gives:
- Bar-cobar equivalence for the flat models: Omega B_mod(A) ~ A in the derived category
- Spectral sequence convergence for the genus filtration
- But NOT the derived-coderived equivalence at genus g >= 1 — that requires the factorization input (the global D-module data that produces the curved model)

#### Section 4: The full picture

**Theorem (Derived-coderived equivalence: full statement; ProvedHere).** Stated here as the conjunction of three results:
1. Relative homotopy-involutivity of FT_{Com_mod / SC^{ch,top}} (this chapter)
2. Factorization Koszul duality for the curved model (Chapter ref{ch:factorization-swiss-cheese})
3. The chiral Riemann-Hilbert correspondence connecting flat and curved models (Vol I)

This is the single structural principle: the three chain-level models (Remark rem:three-models) are three views of the relative Feynman transform, seen through the local (operadic), global (factorization), and mediating (algebraic) lenses.

**Remark (The three routes unified).**
- Route B (factorization) provides the global truth: factorization on Ran(Sigma_g) x Ran(R) over M-bar_g
- Route A (operadic) provides the local approximation: SC^{ch,top}_mod at formal completions
- Route C (this chapter) provides the algebraic skeleton: FT_{Com_mod / SC^{ch,top}} as the common abstraction
The factorization structure projects to the relative Feynman transform by forgetting D-module data; the operadic structure injects into the relative Feynman transform by formal completion. The relative Feynman transform is the categorical image.

## Part 2: Rewrites of Existing Files

### Remark rem:three-models (modular_pva_quantization_core.tex)

This remark (which we wrote earlier today) needs updating to reflect the full picture. Find the remark labeled rem:three-models and ADD at the end, before \end{remark}:

```
\smallskip\noindent
\textbf{Structural origin.}\enspace
The three models are three views of the relative Feynman transform
$\mathrm{FT}_{\Com_{\mathrm{mod}} / \SCchtop}$
(Chapter~\textup{\ref{ch:relative-feynman}}):
\begin{enumerate}[label=\textup{(\roman*)}]
\item The flat associated graded is the \emph{local} (operadic) view,
  extracted from the factorization structure by formal completion at
  collision points
  (Chapter~\textup{\ref{ch:modular-sc-operad}}).
\item The corrected holomorphic model uses only local data
  (the holomorphic propagator) and is visible to the local operad.
\item The curved geometric model uses global data (the Arakelov
  propagator, D-module monodromy, period corrections) and is visible
  only to the factorization framework
  (Chapter~\textup{\ref{ch:factorization-swiss-cheese}}).
\end{enumerate}
The derived--coderived equivalence connecting models~(i)--(ii)
(derived) with model~(iii) (coderived) is the conjunction of
relative homotopy-involutivity
(Theorem~\textup{\ref{thm:relative-ft-involutive}}) and
factorization Koszul duality
(Theorem~\textup{\ref{thm:fact-sc-koszul-duality}}).
```

### Remark rem:curvature-spectral-sequence (modular_pva_quantization_core.tex)

Find the remark we added today and ADD at the end:

```
The curvature is a \emph{global} datum: it arises from
the D-module monodromy of the factorizable sheaf on
$\operatorname{Ran}(\Sigma_g)$, not from the local
operadic structure. The spectral sequence trades this
global curvature for the genus-mixing operator~$D_1$,
which IS visible to the local operad.
```

### preface.tex Update

Find the paragraph we added today ("The modular bar complex is flat..."). REPLACE with:

```
The modular bar complex is flat ($D^2 = 0$) because it couples all
genera simultaneously via the nonseparating one-edge operator~$D_1$;
the genus-$g$ bar complex on a fixed Riemann surface is curved
($\dfib^{\,2} = \kappa \cdot \omega_g$) because restricting to a
single genus decouples this connection.  The flatness is visible
to the local operadic model; the curvature is a global phenomenon
requiring the factorization framework on $\operatorname{Ran}(\Sigma_g)$.
The equivalence of these two perspectives is the conjunction of
relative homotopy-involutivity
(Chapter~\ref{ch:relative-feynman}) and factorization Koszul
duality (Chapter~\ref{ch:factorization-swiss-cheese}).
```

### main.tex

**Line 873** (part epigraph): If it says "Swiss-cheese algebra is formal", change "algebra" to "coalgebra".

Add after line 909 (Part V, after 3d_gravity):
```latex
\input{chapters/connections/relative_feynman_transform}
```

### CLAUDE.md Updates

1. File map: "- relative_feynman_transform: Part V (relative FT, recognition, involutivity, three routes unified)"

2. Update the "What the Engine Computes" section (very top of CLAUDE.md): the second paragraph currently says "This product is the operadic fingerprint..." — this MUST be changed to match the factorization-primary framing. Change to:
"This product is the local algebraic shadow of the factorization structure of a 3d holomorphic-topological QFT on C_z x R_t..."

3. Update the Critical Pitfall about three models to reference the relative Feynman transform.

## Definition of Done
- [ ] chapters/connections/relative_feynman_transform.tex exists and compiles
- [ ] Contains: relative FT definition, recognition theorem, homotopy-involutivity, full derived-coderived equivalence, three-routes-unified remark
- [ ] Remark rem:three-models updated with structural origin
- [ ] Remark rem:curvature-spectral-sequence updated with global-datum note
- [ ] preface.tex updated
- [ ] main.tex updated (coalgebra fix + \input)
- [ ] CLAUDE.md updated (file map + What the Engine Computes + three-models pitfall)
- [ ] Full build passes
