# Master Spec: The Factorization Foundation of Modular Swiss-Cheese Koszulity

## The Ground Truth

A chiral algebra is not an algebra over an operad. It is a factorizable D-module on the Ran space. The operad is a local shadow — what you see when you pass to formal completions at collision points. At genus 0 on P^1, the local shadow determines the global object (because P^1 minus a point is C). At genus g >= 1, it does not: the periods of Sigma_g, the monodromy of the D-module around cycles, and the factorization pattern on Ran(Sigma_g) carry genuinely more information than the local data. The curvature d_fib^2 = kappa * omega_g is the measure of this excess — it is global information invisible to the local operad.

The manuscript currently papers over this. It proves homotopy-Koszulity of the topological operad SC^{ch,top} and invokes it for the all-genus story. This is a category error: the operadic theorem is about formal discs; the genus tower lives on Ran(Sigma_g) over M-bar_g. The two are connected but not identical, and the connection is itself a theorem requiring factorization-algebraic input.

The correct hierarchy has six layers:

### Layer 0: The geometric substrate
Smooth (or stable nodal) curves X. The Ran space Ran(X) — finite non-empty subsets of X. The moduli stack M-bar_{g,n} of stable curves. The universal Ran space: Ran(X) varying in families over M-bar_g.

### Layer 1: Factorizable D-modules (the BD structure)
A chiral algebra A on X is a factorizable D-module F_A on Ran(X). The factorization isomorphism F_A|_{U_I} ~ box_{i in I} F_A on the complement of diagonals is the structural datum. Chiral operations A^{box n}|_Delta -> A are restrictions to diagonals — maps between D-modules, not maps between fibers at points. The composition of operations goes through the factorization structure on Ran(X). This is the Beilinson-Drinfeld framework [BD04].

### Layer 2: The HT product geometry
The holomorphic-topological manifold Sigma_g x R. A factorizable D-module on Ran(Sigma_g) (holomorphic, closed color) tensored with a locally constant factorization algebra on R (topological, open color). The Swiss-cheese structure: holomorphic in z, associative in t, with bulk -> boundary directionality. The half-space Sigma_g x R_{>= 0} with boundary Sigma_g x {0}.

### Layer 3: The bar complex as factorization coalgebra
The bar complex B(A) is the factorization coalgebra on Ran(X) Koszul-dual to the chiral algebra A. The bar differential = chiral operations dualized via the factorization structure. The bar coproduct = factorization isomorphism dualized. At genus 0: honest dg coalgebra (d^2 = 0). At genus g >= 1: curved factorization coalgebra (d_fib^2 = kappa * omega_g), the curvature arising from the monodromy of the D-module around B-cycles of Sigma_g — global information that the local formal-disc model does not see.

### Layer 4: Koszul duality as factorization duality
The bar-cobar adjunction is factorization algebra <-> factorization coalgebra duality. At genus 0, it reduces to the operadic bar-cobar (the local shadow). At genus g >= 1, it requires factorization on Ran(Sigma_g) over M-bar_g. The Feynman involution FT^2 ~ id is the involutivity of factorization duality for the modular operad {M-bar_{g,n}}.

### Layer 5: Three models as three gauges
The three chain-level models of the genus-g bar complex (Remark rem:three-models) are three gauges of the factorization structure:
- **Flat associated graded** (D_0^2 = 0): the factorization coalgebra in algebraic gauge — formal completion at collision points, where the local operad acts
- **Corrected holomorphic** (D^(g)^2 = 0): the factorization coalgebra in flat gauge — on the universal cover of Sigma_g, where the holomorphic propagator (prime form) defines a flat connection
- **Curved geometric** (d_fib^2 = kappa * omega_g): the factorization coalgebra in single-valued gauge — on Sigma_g itself, where the Arakelov propagator is single-valued but curved
The passage between gauges is the chiral Riemann-Hilbert correspondence. The curvature is the monodromy of the gauge transformation.

### Layer 6: The local shadow
The topological operad SC^{ch,top} = the formal completion of the factorization structure at collision points. Its operation spaces C_*(FM_k(C)) x E_1(m) are the local models of the factorization pattern near the small diagonal of Ran(X). Homotopy-Koszulity of SC^{ch,top} (Theorem thm:homotopy-Koszul) is the theorem that the formal completion faithfully encodes the collision data — the local Koszul duality works. This is a theorem about the local model, valid because FM_k(C) is the local model of Ran(X) near the thin diagonal.

## The Three Chapters

The six layers organize into three chapters, ordered by depth:

### Chapter 1 (Route B): Factorization Swiss-Cheese [PRIMARY — Layers 0-5]
**File:** chapters/theory/factorization_swiss_cheese.tex (Part I, after line 807)

The foundational development. Constructs the Swiss-cheese structure directly from factorization data, without passing through operads. Two parallel treatments:
- **BD factorization** (Layers 0-1): factorizable D-modules on Ran(X), chiral operations as diagonal restrictions, factorization on families over M-bar_g. The genus tower as a family of factorizable D-modules. Curvature as D-module monodromy.
- **CG factorization** (Layer 2): prefactorization algebras on Sigma_g x R, BV quantization, factorization homology. The genus tower from the dependence of the prefactorization algebra on the complex structure of Sigma_g.
- **Equivalence**: precise theorem with full adjectives (char 0, smooth projective, hol TI, locally constant). The equivalence passes through three layers: chiral/closed (BD Thm 4.5.4 / Francis-Gaitsgory), associative/open (Lurie HA 5.4.5.9), mixed/Swiss-cheese (comparison of D-module and cosheaf descriptions).
- **Factorization Koszul duality** (Layers 3-4): bar complex as factorization coalgebra, Koszul duality from Ayala-Francis / Lurie, derived-coderived equivalence at all genera as a consequence of factorization duality.
- **Three models from factorization** (Layer 5): the three gauges derived from the factorization structure.
- **Latyntsev aside**: factorization quantum groups as the closed-color projection; dg-shifted Yangians as the open-color projection; R-matrix as the mixed projection.

### Chapter 2 (Route A): Modular Swiss-Cheese Operad [LOCAL APPROXIMATION — Layer 6]
**File:** chapters/theory/modular_swiss_cheese_operad.tex (Part I, after Route B)

The local/operadic formalization. Honest about its scope: this is what you see at collision points. Defines SC^{ch,top}_mod as a partially modular operad (closed color modular, open color genus-0). Proves modular homotopy-Koszulity. But explicitly acknowledges:
- This captures formal-disc data and misses global D-module structure
- The curvature d_fib^2 = kappa * omega_g is INVISIBLE to the local operad — it requires the factorization input of Chapter 1
- The corrected differential D^(g)^2 = 0 (flat, by Fay trisecant) IS visible to the local operad — the Fay identity is a local identity on configuration spaces
- The passage from local to global (operadic to factorization) is the localization functor, whose properties are a theorem (not a tautology)

The chapter should open by extracting the local operad from the factorization data of Chapter 1, making the derivation explicit. Then prove modular homotopy-Koszulity as a property of this extracted local object.

### Chapter 3 (Route C): Relative Feynman Transform [ALGEBRAIC SKELETON — all layers]
**File:** chapters/connections/relative_feynman_transform.tex (Part V, after line 909)

Formalizes the common algebraic structure. The modular bar complex B_mod is an algebra over FT_{Com_mod / SC^{ch,top}} (the relative Feynman transform). This is the skeleton that both the factorization approach (Chapter 1) and the operadic approach (Chapter 2) flesh out. Recognition theorem: B_mod IS FT_{Com_mod / SC^{ch,top}}. Homotopy-involutivity theorem. The chapter stands intermediate, connecting the global (factorization) and local (operadic) viewpoints through a shared algebraic framework.

## What Each Chapter Must Address

Every chapter must explicitly state its relationship to the local/global tension:
- Chapter 1 resolves it (factorization is the global truth)
- Chapter 2 localizes it (operadic = formal-disc approximation, honest about what it misses)
- Chapter 3 algebraicizes it (relative Feynman transform = the formal skeleton shared by both)

## Shared Conventions

### What Vol I Provides (agents MUST read)
- **Convention 3.1** (higher_genus_foundations.tex:194-258): Three differentials
- **Theorem thm:bar-modular-operad** (bar_cobar_adjunction_curved.tex:5910-5958): Bar as FCom-algebra
- **Theorem thm:feynman-involution** (higher_genus_foundations.tex:5918-5967): FT^2 ~ id
- **Definition def:modular-koszul-chiral** (higher_genus_modular_koszul.tex:331-438): MK1-MK3
- **D^oc_A** (configuration_spaces.tex:2393-2443): Modular open-closed bar differential
- **Coderived models** (appendices/coderived_models.tex:1-100): Positselski CDG framework

### What Vol II Provides (agents MUST read)
- **CLAUDE.md**: All conventions, pitfalls, file map
- **Theorem thm:homotopy-Koszul** (line-operators.tex:39-173): Operadic homotopy-Koszulity
- **Theorem thm:modular-bar** (modular_pva_quantization_core.tex:463-502): B_mod with D^2=0
- **Proposition prop:D0D1** (modular_pva_quantization_core.tex:524-537): Bicomplex D_0, D_1
- **Remark rem:three-models** (modular_pva_quantization_core.tex): Three chain-level models
- **Theorem thm:genus-completion** (modular_pva_quantization_core.tex:655-701): Genus spectral sequence
- **foundations.tex:466-552**: Current HT prefactorization discussion

### LaTeX
- \providecommand for new macros (never \newcommand in chapter files)
- No \newtheorem — all environments in main.tex preamble
- Existing macros: \SCchtop, \Bmod, \barB, \Omegach, \dfib, \dzero, \Dg{g}, \FM, \Conf, \Fact, \Ainf, \Linf, \Eone, \Etwo, \Ass, \Com, \Lie
- Claim tags: \ClaimStatusProvedHere, \ClaimStatusProvedElsewhere, \ClaimStatusConjectured
- Build: pkill -9 -f pdflatex 2>/dev/null || true; sleep 2; make

### Git
All commits by Raeez Lorgat. NEVER credit LLM. Do NOT commit — leave uncommitted for review.
