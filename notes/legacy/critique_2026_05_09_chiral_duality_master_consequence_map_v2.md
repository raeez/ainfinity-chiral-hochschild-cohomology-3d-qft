# Critique of 2026-05-09 — Architectural Reconstitution Map (v2, Strongest Form)

**Supersedes**: `critique_2026_05_09_chiral_duality_master_consequence_map.md` (v1).
**Why a v2**: v1 enumerated consequences; v2 is the architectural theorem that *generates* the enumeration. The five Beilinson cuts of v1 are not five independent cuts — they are five *types of licensing data* for **one** cut. The two structural chains the critique writes on its last page (open lane and CY lane) are not two separate architectures — they are two instances of **one universal stage chain** whose level-3 stages are unified by SC$^{\mathrm{ch,top}}$ heptagon face (6). v2 takes those identifications seriously and rewrites the consequences accordingly.

**Voice note**: this is workshop floor (`notes/`), so AP language, wave language, taxonomy, and bookkeeping vocabulary are allowed. The prose itself is meant to read as mathematics, not as a configuration manual.

---

## 0. The reflexive note (what kind of object this is)

This map is **a recipe, not a deliverable**. The deliverable is the rectified manuscripts. The map is a workshop scaffold whose job is to expire — once the rectification lands, the discipline is in the prose, the macros, the APs, the IV tests, and this note becomes redundant. If a future session is reading this map *instead of* the rectified prose, the rectification has not landed.

The same reflexive rule applies to the critique itself. The critique is a *triage*, not a theorem. Its job is to make space; the theorems are in the chapters.

---

## 1. The Master Architectural Theorem (the universal chain)

> **Every working object in the five programmes (Vol I, Vol II, Vol III, mixed-HT-strings, igusa-cusp-form) occupies one of five stages of the universal chain**
>
> $$
>   \underbrace{\mathsf{P}}_{\text{primitive}}
>     \;\xrightarrow{\;\alpha\;}\;
>   \underbrace{\mathsf{C}}_{\text{chart}}
>     \;\xrightarrow{\;\beta\;}\;
>   \underbrace{\mathsf{S}}_{\text{shadow}}
>     \;\xrightarrow{\;\gamma\;}\;
>   \underbrace{\mathsf{Z}}_{\text{centre / bulk}}
>     \;\xrightarrow{\;\delta\;}\;
>   \underbrace{\mathsf{A}}_{\text{acting object}}
> $$
>
> with arrows $\alpha$ (chart-projection), $\beta$ (shadow / forgetful), $\gamma$ (centre / Hochschild), $\delta$ (extension to maximal action). Each stage has its own legitimate operations. Cross-stage collapse — identifying objects at different stages without naming the licensing arrow — **is the master pattern `shadow = object`** the critique names.

The 17 dismissals are 17 cross-stage collapses. The Beilinson cut is the *licensing declaration* that the chain has stages. The five "cuts" of v1 are five *types* of licensing data, indexed by which arrow ($\alpha, \beta, \gamma, \delta$) the collapsing slogan was eliding, plus a fifth arrow $\epsilon$ for *ambient declarations* that sit outside the chain (completion / topology / pro-system) but are required to make any of the chain stages well-defined.

The universal chain is not metaphorical: each manuscript-level slogan can be located on one stage, and the licensing arrow back to the previous stage can be named.

---

## 2. The two lanes are parallel instances of the universal chain

The critique closes (page 9) with two structural chains:

$$
\textbf{Open lane}\colon\quad
  \openFactCat
    \;\rightsquigarrow\;
  A_b
    \;\rightsquigarrow\;
  B(A_b)
    \;\rightsquigarrow\;
  \Zderch{A_b}
    \;\rightsquigarrow\;
  \mathrm{Line}(A_b).
$$

$$
\textbf{CY lane}\colon\quad
  \mathrm{CY}_d\text{-cat}
    \;\rightsquigarrow\;
  \PhiFA_d(A)
    \;\rightsquigarrow\;
  \SpCh\bigl(\PhiFA_d(A)\bigr)
    \;\rightsquigarrow\;
  Y^+(X)
    \;\rightsquigarrow\;
  D\bigl(Y^+(X)\bigr).
$$

These are **two instances** of the universal chain $\mathsf{P} \to \mathsf{C} \to \mathsf{S} \to \mathsf{Z} \to \mathsf{A}$ — not two different architectures.

| Stage         | Open lane                        | CY lane                                  | Igusa lane                            |
| ------------- | -------------------------------- | ---------------------------------------- | ------------------------------------- |
| $\mathsf{P}$ primitive  | open factorisation dg-cat $\openFactCat$ on $\logCurve$ | CY$_d$-category                          | K3-Heisenberg datum $(L, \theta_L)$    |
| $\mathsf{C}$ chart      | $A_b = \mathrm{End}_\cC(b)$, after choosing vacuum $b$  | $\PhiFA_d(A) = $ stage-1 $E_d$-FA on CY target | $K_0$-determinant package              |
| $\mathsf{S}$ shadow     | $\BarTwc{A_b} = $ twisting / coupling coalgebra         | $\SpCh = $ specialisation along $(\Sigma_{d-1}, C)$ to $E_1$-chiral | $\Delta_5$ Borcherds denominator       |
| $\mathsf{Z}$ centre / bulk | $\Zderch{A_b} \simeq \bulkChirHoch{A_b}$               | chiral shadow on $C$ + its derived centre        | (open: $\mathfrak{D}_X$ centre)        |
| $\mathsf{A}$ acting object | line operators $\mathrm{Line}(A_b)$ acting on bulk     | $D(Y^+(X)) = $ Drinfeld double = full quantum group | $\mathfrak{D}_X$ operator-level (open) |

**The level-3 fact** that unifies the two lanes is the SC$^{\mathrm{ch,top}}$ heptagon face (6) Drinfeld-centre identification:

$$
Z\bigl(\mathrm{Rep}_{\mathrm{fact}}(A)\bigr) \;\simeq\; \mathrm{Rep}_{\mathrm{fact}}\bigl(\Zderch{A}\bigr)^{E_2}.
$$

Reading: the categorified Drinfeld double on the CY-lane Hopf-side and the derived chiral centre on the open-lane Hochschild-side are the *same* level-$\mathsf{Z}$ object viewed through two functorial windows. This is *already proved* (Vol II `chapters/connections/...` heptagon face (6); see FRONTIER.md DEFINITIVE STATUS) — but it is not currently exhibited as the structural unification of the two lanes. v2's first move is to make this unification visible.

A **third lane** is implicit in the 2026-04-22 sharpenings: the Pentagon-trace / ghost-trace / Borcherds-weight three-factor lane.

$$
\textbf{Trace lane}\colon\quad
  \fg_{\mathrm{ghost}}
    \;\rightsquigarrow\;
  Q_{\mathrm{BRST}}^{2}
    \;\rightsquigarrow\;
  \mathrm{tr}_{\mathrm{ghost}}\bigl(Q_{\mathrm{BRST}}^{2}\bigr) = \mathrm{tr}_{\mathrm{Pentagon}} = \omega_{\mathrm{Borcherds}} = c_N(0)/2.
$$

Stages: $\mathsf{P}$ = ghost Lie algebra; $\mathsf{C}$ = chosen BRST nilpotent $Q_{\mathrm{BRST}}$; $\mathsf{S}$ = scalar trace; $\mathsf{Z}$ = Pentagon-trace face of the SC$^{\mathrm{ch,top}}$ closed colour (single-colour MC5 Pentagon); $\mathsf{A}$ = the Borcherds product / lattice-vertex-operator algebra. The eight-form spread ($w(N) \in \{5, 2, 1, 1, 1/2, 1, 1/4, 0\}$) is a parameter on the $\mathsf{S} \to \mathsf{Z}$ arrow indexed by $N \in \{1, 2, 3, 4, 6\}$ + half/quarter-integer continuations.

**Three lanes, one architecture.** The unification: each stage of each lane carries the same five-licensing-type discipline. The 17 dismissals scatter across the three lanes:

| Lane     | Dismissals                                              |
| -------- | ------------------------------------------------------- |
| Open     | 1, 2, 3, 4, 5, 13, 15, 17                               |
| CY       | 6, 7, 8, 9, 10                                          |
| Trace    | 11, 12 (and the Pentagon-trace anchor of #6, #14)        |
| Cross-lane | 14, 16 (endpoint hypotheses bleed into all three)      |

This is the hidden symmetry v1 missed: every dismissal is a cross-stage collapse on a known lane, not a unique error.

---

## 3. Five licensing types (the indexing of the one cut)

Every false slogan is *one* slogan with *one* fix: name the licensing arrow. The arrows partition into five types, *indexed by what kind of data the arrow is*:

### Type α — Chart / scope / log decoration

A chart-projection $\alpha : \mathsf{P} \to \mathsf{C}$ requires choosing a structure on $\mathsf{P}$:

- *Boundary vacuum* $b$ for the open factorisation category $\openFactCat$ (gives $A_b = \mathrm{End}_\cC(b)$).
- *Scope label* for $\kappa$ (which row of the κ-tuple is meant).
- *Tangential log data* $(D, \tau)$ for the curve (gives $\logCurve$ as opposed to bare $X$).
- *BRST nilpotent* $Q_{\mathrm{BRST}}$ for the ghost Lie algebra (gives $\mathrm{tr}_{\mathrm{ghost}}(Q^2)$).
- *Branch / Stokes datum* on the FM compactification (gives Stokes choice in PVA quantisation).

**Dismissals**: 1, 4, 5, 6, partially 16. **Macro skeleton**: `\primPkg`, `\logCurve`, `\kappaCat / \kappaChHodge / \kappaChHeis / \kappaBKM / \kappaFiber`, plus `\bdyVac{b}`.

### Type β — Comparison / functor / natural transformation

A shadow / forgetful / Hochschild / Pfaffian arrow $\beta : \mathsf{C} \to \mathsf{S}$ or $\gamma : \mathsf{S} \to \mathsf{Z}$ requires a comparison map:

- $\mathrm{ChirHoch}^\bullet(A, A)$ functor (Lurie HA.5.5 + chiral Deligne-Tamarkin).
- Drinfeld-double construction $D : Y^+ \mapsto D(Y^+)$ (Hall pairing + completion + integral form + stable-envelope transport + descent).
- Protected Pfaffian $\mathrm{Pf}_{\mathrm{prot}} : \mathfrak{D}_X \mapsto \Delta_5$.
- Hall--Borcherds residual: scalar trace $\to$ gravity-line partition function.
- Stage-2 specialisation $\SpCh = \int_{\Sigma_{d-1}}$ : $E_d$-FA on CY target $\to$ $E_1$-chiral on $C$.
- MC bijection $\mathrm{Hom}(A, B) \hookrightarrow \mathrm{MC}(A^! \otimes B)$ (Gui-Li-Zeng injection; bijectivity Koszul-effective-special).

**Dismissals**: 2, 3, 7, 8, 11, 12, 17. **Macro skeleton**: `\bulkOf`, `\bulkChirHoch`, `\BarTwc`, `\PhiFA`, `\SpCh`, `\Yplus`, `\Drinfdouble`, `\protectedPfaff`, `\operatorPrim`.

### Type γ — Ambient declaration

An ambient $\gamma$-declaration is **outside** the universal chain but conditions every stage of it. It is the choice of:

- *Chain-level ambient* (ordinary $\mathrm{Ch}(\mathrm{Vect})$ vs weight-completed vs pro-object vs $J$-adic).
- *Categorical lane* (chain-level vs $(\infty, 1)$-categorical via Lurie HA.5.5 — both are load-bearing per CLAUDE.md, never substitutable).
- *Topological completion* (analytic / Schwartz / formal) on the Stokes / FM side.

**Dismissals**: 15. **Macro skeleton**: `\hypAmbientWtCpl`, `\hypAmbientPro`, `\hypAmbientJadic`, plus a parametric `\hypothesis{tag}{description}` macro for indefinite expansion.

The chain-level / $(\infty,1)$-categorical lane discipline is *a Type γ declaration* (CLAUDE.md "Chain-level and $(\infty,1)$-categorical: equal status"). v1 missed this: ambient is not just *a* type of licensing, it is the *prerequisite* for every other licensing type. Without ambient, the centre $\Zderch{A}$ is undefined.

### Type δ — Endpoint hypothesis / convergence

An $\mathsf{S} \to \mathsf{Z}$ or $\mathsf{Z} \to \mathsf{A}$ arrow that survives a *limit* operation requires hypotheses:

- $W_\infty[\lambda] \Rightarrow E_\infty$ : Prochazka triangular truncation + Creutzig-Kanade-Linshaw parafermion compatibility + Pope-Romans-Shen / Bakas + Yamada weight-window.
- All-loop quantum HT from PVA λ-Jacobi : KZ analytic SDR + Stokes choices + reflected weights + lift of $T = [Q_{\mathrm{tot}}, G]$.
- Universal Holography $\Rightarrow$ dynamical metric : Brown-Henneaux dictionary + modular invariance + vacuum dominance + saddle dominance.
- Mittag-Leffler / antighost-commutativity for the iterated Sugawara ladder.
- Class M chain-level via inverse-system convergence of $S_4$.

**Dismissals**: 13, 14, 16, partially 10. **Macro skeleton**: `\hypProchazka`, `\hypCKL`, `\hypPRSh`, `\hypYamada`, `\hypKZSDR`, `\hypStokes`, `\hypReflWts`, `\hypTLift`, plus parametric `\hypothesis`.

### Type ε — Effectiveness / orientation / non-degeneracy

A $\mathsf{Z} \to \mathsf{A}$ arrow that requires a special-case effectiveness condition:

- *Koszul effectiveness* : MC injection $\to$ MC bijection (Gui-Li-Zeng vs Vol I Theorem B Positselski).
- *Holomorphic de Rham vanishing* : formal-local Hamiltonian model $\to$ global compact background ($H^1(X, \C_X)$ obstruction).
- *Orientation* / Pfaffian sign datum : $\Delta_5^2$ (squared scalar) $\to$ $\Delta_5$ (signed Pfaffian).
- *6d hCS quartic obstruction vanishing* : $\PhiFA_3$ on verified loci ($\chi_{\mathrm{top}}(K3 \times E) = 0$ kills the quartic anomaly polynomial integration).
- *PBW / no-extra-relations* / parity / centre / associator / Mittag-Leffler defect vanishing for the K3$\times E$ $\Delta_5$ recognition (Wave 13 R1).
- *MO Yang-Baxter unitarity* / cocycle condition for Maulik-Okounkov R-matrix as gluing residue.

**Dismissals**: 9, 10, 11, 17. **Macro skeleton**: `\effKoszul`, `\effHTGlobalDR`, `\effPfaffOrient`, `\effHCSQuartic`, `\effPBWnoExtra`, `\effMOcocycle`.

### How the five types index the dismissals

| Dismissal | Type | Stage collapsed                                                             |
| --------- | :--: | --------------------------------------------------------------------------- |
| 1         |  α   | $\mathsf{P} \to \mathsf{C}$ (primitive treated as chart)                    |
| 2         |  β   | $\mathsf{S} \to \mathsf{Z}$ (shadow treated as centre)                      |
| 3         |  β   | $\mathsf{S} \to \mathsf{Z}$ (bar mechanism stands in for Deligne-Tamarkin)  |
| 4         |  α   | log $\to$ bare (curve loses tangential decoration)                          |
| 5         |  β   | shadow trace $\to$ adjective (open trace projected onto closed shadow)      |
| 6         |  α   | tuple $\to$ scalar (κ-row scope label dropped)                              |
| 7         |  β   | composite-functor $\to$ direct-functor (two-stage Φ collapsed)              |
| 8         |  β / δ | $\mathsf{C} \to \mathsf{A}$ (positive half identified with full quantum group) |
| 9         |  ε   | analogy treated as identification (3d-CS knot $\to$ 6d hCS without quartic) |
| 10        |  ε   | local $\to$ global (formal Darboux without de Rham obstruction)             |
| 11        |  β   | shadow $\to$ operator ($\Delta_5$ as Hilbert space)                         |
| 12        |  β   | scalar trace $\to$ operator algebra ($Z_{\mathrm{BPS}}$ as path integral)   |
| 13        |  δ   | functor $\to$ dynamical theory (Universal Holography $\to$ metric integral) |
| 14        |  δ   | endpoint $\to$ limit theorem ($W_\infty \to E_\infty$ unconditional)        |
| 15        |  γ   | weight-completed $\to$ ordinary (class M ambient dropped)                   |
| 16        |  δ   | classical $\to$ quantum (PVA Jacobi $\to$ all-loop)                         |
| 17        |  ε   | injection $\to$ bijection (quadratic dual $\to$ Koszul duality)             |

**The fact**: every dismissal is exactly one cross-stage collapse, with the licensing arrow correctly typed. The classification has no leftover and no overlap (8 dismissals are β-type, 4 are α, 3 are δ, 2 are ε, 1 is γ; #8 is dual α/β at the chart-functor boundary).

This typology *replaces* the v1 "five Beilinson cuts." It is denser, more falsifiable, and directly enforceable at the macro level: every theorem statement carries an explicit licensing-type tag.

---

## 4. Heptagon face ↔ licensing-type correspondence

The SC$^{\mathrm{ch,top}}$ heptagon has seven faces (`chapters/theory/sc_chtop_heptagon.tex:197`). The licensing types map onto faces:

| Heptagon face                              | Licensing type | Cut in v1 | Dismissals it gates                  |
| ------------------------------------------ | :------------: | :-------: | ------------------------------------ |
| (1) two-coloured operadic primitive         |       α        |     A     | 1, 4, 5                              |
| (2) quantum complementarity (κ + κ⁺)        |       α        |     C     | 6                                    |
| (3) bar-cobar adjoint inversion             |       β        |     B     | 2, 3, 17                             |
| (4) brace / Deligne-Tamarkin action         |       β        |     B     | 13                                   |
| (5) topologisation ladder ($E_n$ promotion) |       δ        |     E     | 14                                   |
| (6) Drinfeld centre = derived centre        |       β        |     B + D | 8                                    |
| (7) PTVV derived-AG witness                 |       γ + ε    |     E     | 9, 10, 11, 12, 15                    |

This mapping is **already available** in the manuscript: face (6) is the unification we made explicit in §2. The architectural value is that the seven heptagon faces and the five licensing types are **two presentations of one structure**:

- The heptagon is the *operadic / categorical* presentation.
- The licensing-type taxonomy is the *epistemic / proof-discipline* presentation.

Each is the dual of the other under: face = stage-pair to be licensed, licensing-type = data needed to license that stage-pair.

This is the **deep symmetry** v1 missed.

---

## 5. The cyclic Hochschild discipline

Cut B identifies $\mathrm{Bar}(A) \neq \mathsf{bulk}$ and $\Zderch{A} \simeq \mathrm{ChirHoch}^\bullet(A, A)$. v1 stopped there. v2 sharpens: there are **four parallel chiral Hochschild objects**, and all-but-one is a shadow of the centre:

$$
\mathrm{ChirHoch}^\bullet(A, A) \quad
\mathrm{ChirHoch}_\bullet(A, A) \quad
\mathrm{HC}^{\mathrm{cyc}}(A) \quad
\mathrm{HC}^{-}(A).
$$

The bulk is $\mathrm{ChirHoch}^\bullet(A, A)$ (cohomology, *not* homology). The other three are shadows. Their identification with the bulk requires further licensing data:

- $\mathrm{ChirHoch}_\bullet(A, A) \simeq \mathrm{ChirHoch}^\bullet(A, A)$: Calabi-Yau / dualising-line compatibility. License = orientation datum.
- $\mathrm{HC}^{\mathrm{cyc}}(A) \to \mathrm{ChirHoch}^\bullet(A, A)$: $S^1$-equivariant Hochschild spectral sequence. License = $S^1$-action on bar complex (not free in chain-level ordinary; weight-completed required).
- $\mathrm{HC}^{-}(A)$: negative cyclic is the target of the Chern character from the $K$-theory of the CY module category; it is not equal to $K$-theory. License = compact CY trace datum + Bridgeland orientation/chamber.

The current Vol II `chapters/connections/hochschild.tex` is correct on $\mathrm{ChirHoch}^\bullet(A, A) = $ bulk and now articulates the four-object discipline. Reconstitution work: **keep the four objects and their forward arrows separated, especially the Chern-character direction $K \to \mathrm{HC}^{-}$**.

This is a **β-type sharpening** within Cut B that v1 missed.

---

## 6. The closed-colour primitive package

v1 Cut A = the seven-tuple $\primPkg{X}{D}{\tau}{b}$ for the **open** colour. The closed colour has its own primitive package, which v1 silently elided.

The closed-colour primitive on a smooth complex curve $X$ (no log decoration: closed sector is global, not boundary):

$$
\bigl(X;\;\; \cD^{\mathrm{cl}},\;\; F^{\mathrm{cl}},\;\; \mathrm{tr}^{\mathrm{closed}}_{X}\bigr)
$$

with $\cD^{\mathrm{cl}}$ a closed factorisation $\infty$-category (chain-level: factorisation D-modules on $\mathrm{Ran}(X)$ at the BD-side), $F^{\mathrm{cl}}$ a closed colour vacuum factorisation algebra ($E_2$-on-curve factorisation algebra), $\mathrm{tr}^{\mathrm{closed}}_X$ the *closed-colour trace* (which is precisely what evaluates to the Pentagon-trace value $c_N(0)/2$ at the modular endpoint).

The **directional restriction** of $\mathsf{SC}^{\mathrm{ch,top}}$ — closed colour cannot map to open colour, but open can clutch through closed — is the operadic shadow of the **stage-1 vs stage-2 distinction**. Stage-1 is a closed-colour primitive on the CY target; stage-2 is the open-colour primitive on the curve. The two-stage factorisation $\Phi_d = \SpCh \circ \PhiFA_d$ is a **forgetful** from closed-colour stage-1 to open-colour stage-2.

**Reconstitution work** v1 missed: write the closed-colour primitive package as an explicit Architectural Definition next to the open-colour one. Together they form the *bicoloured primitive package* of the heptagon. This is a **Cut A sharpening** that resolves dismissal 5 (modularity) more cleanly: modularity is the closed-colour trace evaluating to $c_N(0)/2$, with open-colour trace + clutching as the *operadic origin* of the closed-colour modular structure.

The bicoloured primitive package macro:

```latex
\providecommand{\primPkgBicolour}[2]{
  \bigl(
    \openFactCat\bigl|_{\logCurve}\,;\;
    \cD^{\mathrm{cl}}\bigl|_{X}\,;\;
    \bdyVac,\;\,
    \bdyAlg,\;\,
    F^{\mathrm{cl}};\;\,
    \HalfBraid,\;\,
    \mathrm{tr}^{\mathrm{cl}}_{X};\;\,
    \TraceC
  \bigr)
}
```

This nine-tuple replaces the v1 seven-tuple. It is the *correct primitive package* for the Vol II $\mathsf{SC}^{\mathrm{ch,top}}$ + 3D HT QFT setting.

---

## 7. The four new construction problems the cut creates space for

The critique closes by saying *primitive objects first, shadows second, scalar modular forms last.* The implied generative question: once shadows are correctly demoted, what new operator-level primitives are being asked for?

Four explicit construction problems crystallise:

### Problem 1 — $\mathfrak{D}_X$ for K3$\times E$ such that $\mathrm{Pf}_{\mathrm{prot}}(\mathfrak{D}_X) = \Delta_5$

Source: Igusa `~/igusa-cusp-form/main.tex:110-118`. The missing primitive is a first-order protected operator algebra on K3$\times E$ whose protected Pfaffian, after a Pfaffian-to-automorphic line comparison, yields $\Delta_5$. The construction problem is:

> Construct $\mathfrak{D}_X$ as a sheaf of $\Z/2$-graded ChirAlg objects on $X = $ K3$\times E$ whose protected Pfaffian module $\mathrm{Pf}_{\mathrm{prot}}(\mathfrak{D}_X)$ is canonically isomorphic to the Borcherds product $\Delta_5(\Omega)$ on the Siegel upper half-space, with orientation datum producing the sign distinguishing $\Delta_5$ from $-\Delta_5$.

This is *the* operator-level construction problem of the Igusa programme. Vol II Wave 13 retraction R1 says it is *recognition-conditional* on five compact-provenance source-matrix gates. The construction itself is not yet attempted.

**Cross-volume locus**: igusa $\to$ Vol III $\to$ Vol II climax.

### Problem 2 — Gravity-line operator algebra whose scalar trace is $Z_{\mathrm{BPS}}^{K3 \times E} = \Delta_5^{-2}$

Source: Vol II `3d_gravity.tex:8429` `conj:gravity-line-hall-borcherds-comparison`. The construction problem is:

> Construct a chain-level $\mathsf{SC}^{\mathrm{ch,top}}$-pair $(\mathrm{boundary}_{\mathrm{grav}}, \Zderch{\mathrm{boundary}_{\mathrm{grav}}})$ on K3$\times E$ such that the scalar trace of the Pentagon face on the closed colour evaluates to $\Phi_{10}^{\mathrm{un}} = \Delta_5^{2}$, with explicit Hall-Borcherds intertwiner from the K3$\times E$ Hall-Drinfeld double constructed in Wave 13.

This is the **gravity-line construction problem**. It absorbs Problem 1 (since $\mathrm{Pf}_{\mathrm{prot}}(\mathfrak{D}_X)^{-2}$ is a candidate scalar trace).

### Problem 3 — Unified PVA-quantum HT theory

Source: dismissal 16. The construction problem is:

> Construct a *single* HT theory $T_{HT}$ such that its classical limit recovers Khan-Zeng's PVA $\lambda$-Jacobi structure and its all-loop quantum extension carries the $E_3$-lift on $Q_{\mathrm{tot}}$-cohomology. The bridge is the lift $T = [Q_{\mathrm{tot}}, G]$ and the analytic data $(\hypKZSDR, \hypStokes, \hypReflWts)$.

The PVA-classical-vs-quantum gap (current `pva-descent.tex`) is precisely the asymmetry where the licensing data $\hypKZSDR, \hypStokes, \hypReflWts, \hypTLift$ is the construction obligation.

### Problem 4 — Koszul theorem at chiral generality

Source: dismissal 17. Gui-Li-Zeng give MC injection; Vol I Theorem B (Positselski coderived) gives full Koszul duality at *associative* generality. The construction problem is:

> Construct chiral-categorical Koszulness extending Vol I Theorem B: state and prove a chiral-Positselski theorem whose specialisation to ordinary algebras recovers Vol I, and whose specialisation to quadratic chiral algebras recovers Gui-Li-Zeng MC bijection.

This is the **chiral Positselski problem**. It is the natural completion of Vol I Theorem B at chiral generality; the *gap* is precisely what the critique's #17 names.

The four problems are coherent: they form a single research line on operator-level constructions of objects whose scalar / shadow / classical / injection-only forms are currently named. **They are the true ideas the cut creates space for.** Reconstitution work should explicitly inscribe them as Construction Problems in `FRONTIER.md` Wave 14+ queue.

---

## 8. Bidirectional cross-volume propagation

v1 mapped Vol II corrections forward to Vol I and Vol III. v2 acknowledges: when Vol II tightens its discipline, the corresponding Vol I and Vol III prose may need *back-edits*. Examples:

- **Cut A back-edit to Vol I**: Vol I `bar-cobar` chapter currently begins "Let $A$ be a chiral algebra." Per Cut A, Vol I should *also* state $A = A_b$ for some boundary chart $b$ in $\openFactCat$. The Vol I bar-cobar duality is then an inversion *on the level of charts*, not on primitives.
- **Cut C back-edit to Vol I**: Vol I five-archetype ceiling table $\{0, 8, 13, 250/3, 98/3\}$ uses $\kappa_{\mathrm{ch}}$ in the Hodge-supertrace sense. Vol I should explicitly tag every $\kappa$ usage as $\kappaChHodge$ to disambiguate from $\kappaChHeis$, $\kappaBKM$.
- **Cut D back-edit to Vol III**: Vol III `quantum_chiral_algebras.tex` should state $G(X) = D(Y^+(X))$ at theorem level (currently main-text prose at lines 550-564), and elevate `Y^+(X) = H^\bullet_{\mathrm{eq}}(\cM^+_{\mathrm{eff}}(X), \phi_W)` to a load-bearing theorem.
- **Cut E back-edit to Vol I**: Vol I Theorem B (Positselski coderived) is the candidate Koszul theorem. Vol I should mark Theorem B as *the chiral Koszul candidate* for Problem 4 (chiral Positselski).
- **Cut E back-edit to mixed-HT-strings**: discipline at `~/mixed-holomorphic-topological-strings/main.tex:3200-3210, 3232-3233` is correct. Vol II `bv_ht_physics.tex` and `holomorphic_topological.tex` should *cross-cite* this discipline at every HT import. Bidirectional: mixed-HT-strings should *cross-cite* Vol II's HT-bulk-boundary chapters where the import sites are.
- **Cut B back-edit to igusa**: igusa `main.tex:110-118` names the construction problem; igusa should also state the four-object Hochschild discipline (§5 of v2) as the level-$\mathsf{Z}$ stage for $\Delta_5$.

The bidirectional propagation is a single audit: for each primitive theorem in each volume, check whether the licensing-type taxonomy exposes a stronger / sharper / more stage-disciplined statement. Add a `propagate` skill invocation across all five programmes in Phase 7.

---

## 9. Build-system layer: `make verify-licensing`

A new Makefile target:

```makefile
verify-licensing:
	@bash scripts/verify_licensing.sh
```

The script runs a mechanical grep audit:

1. **Cross-stage collapse detection**: for each forbidden phrasing in v2 §10 below, grep Vol II `chapters/`. Any match outside a `\begin{remark}{historical-frame}` or `\begin{verbatim}` block is a violation. List with file:line.
2. **Macro discipline**: every `\kappa` not within `\kappaCat`/`\kappaChHodge`/`\kappaChHeis`/`\kappaBKM`/`\kappaFiber`/`\Kkappa`/`\kappaTuple` is a violation.
3. **Hypothesis tag presence**: every theorem statement in `e_infinity_topologization.tex`, `pva-descent.tex`, `bv_ht_physics.tex`, `wn_tempered_closure_platonic.tex`, `programme_climax_platonic.tex` must carry at least one `\hyp...` tag in its statement (not its proof).
4. **Heptagon face cross-reference**: every theorem citing `face (6)` must also cite `\ref{rem:heptagon-face-6-drinfeld-derived}` (and parallel for face 7).
5. **Two-stage Φ discipline**: every `\Phi_{d}` *not* in a `\providecommand` line must appear within the same paragraph as `\PhiFA` or `\SpCh`.

The script is a $O(\text{file count} \times \text{phrase count})$ grep; runtime ~30 seconds on Vol II (199 chapters as of 2026-05-09).

This is the *enforcement layer*. v1 had APs (caught by `beilinson-gate.sh` PostToolUse hook, file-by-file); v2 adds the *whole-volume* enforcement target, runnable on demand and as a pre-commit gate.

**Cross-volume**: the same `verify_licensing.sh` runs against Vol I, Vol III, mixed-HT-strings, igusa with an `--volume vol1` flag. Phase 1 of reconstitution includes scaffolding this script.

---

## 10. The 17-line voice table (refined from v1)

This is the prose-level grep target for `verify_licensing.sh`. Each row: forbidden $\to$ allowed.

| #  | Forbidden                                                                  | Allowed                                                                                                                                                           |
| -- | -------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1  | "Let $A$ be a chiral algebra"                                             | "Let $\primPkgBicolour{X}{(D, \tau)}$ be the bicoloured primitive package on the tangentially-decorated curve $\logCurve$, and $A_b = \mathrm{End}_{\openFactCat}(b)$" |
| 2  | "$\mathrm{Bar}(A)$ is the bulk"                                            | "$\bulkOf{A} = \bulkChirHoch{A}$ is the bulk; $\BarTwc{A}$ is the twisting/coupling coalgebra"                                                                     |
| 3  | "the $E_1$-bar direction explains $2d \to 3d$"                             | "boundary $A_\infty$-chiral $\rightsquigarrow \bulkOf{A}$ via chiral Deligne-Tamarkin / Swiss-cheese"                                                              |
| 4  | "open sector on $X$" / "open colour on the curve $X$"                      | "open sector on $\logCurve = (X, D, \tau)$"                                                                                                                       |
| 5  | "the closed chiral algebra is modular"                                     | "the open category carries cyclic trace + clutching ($\TraceC$); the closed-colour shadow has modular consequences via $\mathrm{tr}^{\mathrm{cl}}_X$"               |
| 6  | "$\kappa = $ \dots" without row tag                                        | "$\kappaTuple{A}$ takes values $(\kappaCat, \kappaChHodge, \kappaChHeis, \kappaBKM, \kappaFiber)(A) = \dots$"                                                       |
| 7  | "$\Phi_d : \mathrm{CY}_d\text{-Cat} \to \mathrm{ChirAlg}$"                  | "$\Phi_d^{(\Sigma_{d-1}, C)} = \SpCh \circ \PhiFA_d$"                                                                                                              |
| 8  | "$\mathrm{CoHA}(\C^3) = W_{1+\infty}$"                                      | "$\CoHA{\C^3} = \Yplus{\fgl_1}$; $\Wonepinf$ is the Fock representation of $\Drinfdouble{\Yplus{\fgl_1}}$"                                                          |
| 9  | "6d hCS = 3d Chern-Simons in disguise"                                     | "6d hCS realises $\PhiFA_3$ on verified loci; one-loop obstruction quartic $\int_X \mathrm{Tr}_{\mathrm{ad}} A(F_A)^3$; 3d-CS = analogy"                          |
| 10 | "formal local HT $\Rightarrow$ compact global theory"                       | "formal local HT + $\effHTGlobalDR$ + descent + QME + anomaly + locality $\Rightarrow$ candidate compact background"                                              |
| 11 | "$\Delta_5$ = compact BPS Hilbert space"                                  | "$\Deltafive$ = Borcherds denominator shadow of $\operatorPrim{X}$; $\protectedPfaff{\operatorPrim{X}} = \Deltafive$ is the construction problem"                  |
| 12 | "$Z_{\mathrm{BPS}}$ = gravitational path integral"                          | "$\ZBPS = \Phitenun^{-1} = \Deltafive^{-2}$ = protected scalar shadow; gravity-line via Hall-Borcherds residual"                                                  |
| 13 | "Universal Holography constructs 3d quantum gravity"                       | "Universal Holography functor identifies $(\bdyAlg, \bulkOf{A}, \SCchtop\text{-brace})$; for $A = \Vir_c$ it is the holographic boundary-CFT reading"             |
| 14 | "$W_\infty[\lambda] \Rightarrow E_\infty$"                                  | "$W_\infty[\lambda] \Rightarrow E_\infty$ conditional on $\hypProchazka, \hypCKL, \hypPRSh, \hypYamada$"                                                            |
| 15 | "class M chain-level holds"                                                | "class M chain-level holds in $\hypAmbientWtCpl$; fails in ordinary $\mathrm{Ch}(\mathrm{Vect})$"                                                                  |
| 16 | "PVA $\lambda$-Jacobi $\Rightarrow$ quantum theory"                         | "PVA $\lambda$-Jacobi $\Rightarrow$ classical gauge invariance; quantum conditional on $\hypKZSDR, \hypStokes, \hypReflWts, \hypTLift$"                            |
| 17 | "quadratic chiral duality = Koszul duality"                                | "quadratic dual gives candidate dual + MC injection; Koszulness conditional on $\effKoszul$ (Vol I Theorem B / Positselski)"                                       |

The voice table replaces the v1 §7 table with refined conditional-licensing tags, in the §3 typology.

---

## 11. Licensing decision tree (cheat sheet)

For any new statement / theorem / remark, ask the questions in order:

```
1. WHAT STAGE is the statement at?
   P (primitive) ........... open-fact-cat / CY-cat / Lie-algebra-of-ghosts / CY-datum
   C (chart) ............... A_b / Φ^FA_d(A) / Q_BRST / K_0-determinant
   S (shadow) .............. Bar(A) / Sp^ch(stage-1) / scalar trace / κ-row / MC-injection
   Z (centre / bulk) ....... Z^der_ch(A) = ChirHoch^•(A,A) / chiral shadow on C / Pentagon face / D(Y^+) at Hopf-Hochschild level
   A (acting object) ....... Line(A_b) / D(Y^+) full / boundary CFT at Brown-Henneaux

2. WHICH LANE does the object belong to?
   Open / CY / Trace / Cross-lane

3. WHAT LICENSING TYPE bridges the previous stage?
   α (chart / scope / log)       — choose b, name κ-row, name (X, D, τ)
   β (comparison / functor)      — name ChirHoch / D(·) / Pf_prot / specialisation
   γ (ambient declaration)       — declare weight-completed / pro / J-adic / chain-level vs (∞,1)
   δ (endpoint hypothesis)       — name Prochazka / CKL / PRSh / Yamada / KZ SDR / Stokes / etc.
   ε (effectiveness condition)   — name Koszul-effective / vanishing H¹(X, C_X) / orientation / quartic vanishing / PBW

4. IS the statement licensed?
   - All applicable types named at theorem statement (not introduction): yes / no
   - Cross-volume macros used: yes / no
   - Heptagon face cited (if relevant): yes / no
   - APs flagged or refuted: yes / no
   - IV test exists or scheduled: yes / no

5. IF licensed: write the statement with all tags inline.
   IF not licensed: do not write the statement; identify the missing licensing data and either supply it or state the construction problem.
```

This tree is the **operating manual** for the Vol II writer (human or agent) post-reconstitution. It absorbs the §3 typology, §6 bicoloured primitive, §7 construction problems, §10 voice table, and the v1 antipattern catalogue.

---

## 12. Honest phasing (25-35 sessions)

v1 estimated 14-20 sessions. v2 corrects: the realistic budget is **25-35 sessions** including:

- Phase 1 macro + script scaffolding: 2-3 sessions (was 1-2). The bicoloured primitive package adds a non-trivial macro layer; `verify_licensing.sh` requires test coverage.
- Phase 2 architectural Definition + canonical statements: 3-4 sessions. Adding the closed-colour primitive adds a chapter-level inscription.
- Phase 3 AP + cache append: 1 session.
- Phase 4 GAP-closing (4 GAPs, parallelisable): 4-6 sessions (4 agents × 1-1.5 sessions each).
- Phase 5 voice sweep via chriss-ginzburg-rectify: 4-6 sessions (one CG-rectify run per chapter touched in Phase 4 plus the bicoloured primitive package chapters; each rectify is 30-60 min agent runtime).
- Phase 6 IV scaffolding: 4-6 sessions (8 test files; each requires disjoint-route witness design).
- Phase 7 cross-volume integrity audit + bidirectional back-edits: 3-4 sessions.
- Phase 8 final Beilinson audit + audit-skill invocation: 2-3 sessions.
- Phase 9 (new in v2) Construction-Problem inscription in FRONTIER.md Wave 14+: 1 session.
- Phase 10 (new in v2) Cyclic-Hochschild discipline articulation: 1 session.
- Phase 11 (new in v2) Heptagon-face / licensing-type duality remark in `sc_chtop_heptagon.tex`: 1 session.

**Parallelism**: Phase 4 (4 agents), Phase 5 (per-chapter parallel), Phase 6 (per-test-file parallel), Phase 7 back-edits (per-volume parallel) — the wallclock can compress to 12-15 sessions if maximum-parallel.

**Strict sequentiality**: Phases 1 → 2 → 3, then Phase 4-6 in parallel, then Phase 7, then Phases 8-11 in any order, with Phase 8 last.

Time honesty matters because v1 phasing under-budgeted Phase 5 (CG-rectify is expensive) and Phase 6 (disjoint-route IV witnesses are non-trivial), and entirely omitted Phases 9-11.

---

## 13. Adversarial audit of v2 itself

A reflexive cut: is v2 itself an instance of `shadow = object`?

**Audit 1**: Is the Master Architectural Theorem (§1) a primitive or a shadow? It is a *taxonomic statement* about the manuscripts. If treated as an *intrinsic theorem* (a structural fact about $E_n$-operadic semantics, say), it would be a shadow. v2 treats it as a *meta-organisational claim*: every working object the manuscripts study can be located on the chain. The chain itself is not a Vol I/II/III theorem; it is a *reading discipline*. This passes the cut.

**Audit 2**: Are the five licensing types (§3) primitive or a shadow? They are a *typology*, not an intrinsic structure. The actual primitive is the seven-faced heptagon (Vol II) plus the two-stage CY-chiral functor (Vol III). The typology projects onto these via the §4 correspondence. v2 declares the typology to be a *shadow of the heptagon + two-stage functor*. This is honest.

**Audit 3**: Are the four construction problems (§7) primitive or shadows? They are *primitive construction obligations* — each names an operator-level object to construct, with explicit licensing data. v2 inscribes them as Wave 14+ FRONTIER queue entries; they are operator-level primitives (Type-β arrows from existing shadows to constructed centres). This passes.

**Audit 4**: Is the bicoloured primitive package (§6) a true sharpening of Cut A or a notational complication? Test: does it absorb dismissal 5 (modularity) more cleanly? Yes — modularity becomes the closed-colour trace $\mathrm{tr}^{\mathrm{cl}}_X$ evaluating to $c_N(0)/2$, with the open-colour trace + clutching as the operadic origin. This is a *substantive* sharpening, not a notational one.

**Audit 5**: Is the cyclic-Hochschild discipline (§5) a real low-hanging fruit or made-up? Test: do the four objects $\mathrm{ChirHoch}^\bullet, \mathrm{ChirHoch}_\bullet, \mathrm{HC}^{\mathrm{cyc}}, \mathrm{HC}^-$ all appear in Vol II / Vol III prose? Vol II `hochschild.tex` references $\mathrm{ChirHoch}^\bullet$ extensively; $\mathrm{HC}^-$ appears in CY contexts (Bridgeland orientation). Yes, real. The discipline is a sharpening.

**Audit 6**: Is the bidirectional propagation (§8) genuine or a make-work fan-out? Test: are the proposed Vol I and Vol III back-edits *distinct* corrections, or just restatements of forward propagation? The Vol I bar-cobar back-edit (insert $A = A_b$) is a real change, not a restatement. The Vol III $G(X) = D(Y^+(X))$ elevation is a real change (theorem-grade vs prose). The Vol I κ-tag is a real change (every $\kappa$ in Vol I currently bare). Genuine.

**Audit 7**: Does v2 unify what should be unified, or merge what should stay distinct? The heptagon face / licensing-type duality (§4) is a *bijection*, not a conflation. The open-lane / CY-lane parallelism (§2) is a *parallelism via heptagon face (6)*, not an identification. The trace lane is *parallel* to but *not a forgetful image of* the open or CY lanes. v2 does not over-merge.

**Audit 8**: Does v2 invent any new mathematical objects without licensing? The bicoloured primitive package nine-tuple is a *combination* of Vol II's existing open and closed colour primitives (both in `sc_chtop_heptagon.tex`). The four construction problems are *named at* the existing manuscript loci (igusa main.tex:110-118; Vol II 3d_gravity.tex:8429; pva-descent.tex; bar-cobar-review.tex). The licensing-type taxonomy is *new* but is a meta-structure, not a mathematical object. No unlicensed inventions.

**Audit 9**: Where does v2 still fall short? Two known gaps:
- The mathematics of the *passing-along-stages* arrows ($\alpha, \beta, \gamma, \delta, \epsilon$) is not stated as a Vol II theorem. It should be. Specifically, the meta-theorem "the universal chain $\mathsf{P} \to \mathsf{C} \to \mathsf{S} \to \mathsf{Z} \to \mathsf{A}$ is functorial in the data" is *not proved* in any of the five programmes. v2 punts on this.
- The v2 typology is a *static* classification. The dynamics of how an object *flows* along the chain — e.g. how $\PhiFA_d$ specialises through $\SpCh$ to $E_1$-chiral — is implicit, not articulated as a *deformation* or *RG-flow*-style theorem. There may be a deeper RG-style architecture beneath v2.

These are explicit further research directions, named for honesty.

**Audit 10**: Is v2 clean enough to retire v1? Yes. v2 supersedes v1 in every section: §1 (Master theorem) replaces v1's §0 plus §1; §3 (five licensing types) replaces v1's "five Beilinson cuts" cleanly; §4 (heptagon ↔ licensing duality) is new; §5-7 are new sharpenings; §8-12 refine v1's §3-9. The v2 phasing (§12) is a strict superset of v1's. Mark v1 as *first-pass; superseded by v2*.

---

## 14. What the cut actually purchases (executive summary)

After Phases 1-11 land, the five-volume programme has the following invariants — stated operationally (each invariant is mechanically checkable):

1. **Stage discipline**: every primitive object names its stage on the universal chain (§1) and its lane (§2). No statement is allowed to span stages without naming the licensing arrow.
2. **Licensing-type tags**: every theorem statement carries explicit α/β/γ/δ/ε licensing-type tags at the *statement* level. `verify_licensing.sh` enforces.
3. **Bicoloured primitive package**: every "Let $A$ be \dots" replaced by `\primPkgBicolour{X}{(D, \tau)}`.
4. **κ-tuple discipline**: bare `\kappa` forbidden; tuple components used.
5. **Two-stage Φ + heptagon-face citation**: every $\Phi_d$ accompanied by $\PhiFA_d, \SpCh$, and a heptagon face reference.
6. **Shadow-vs-operator manifest**: eight shadows named with operator-level primitives. Four open Construction Problems inscribed.
7. **Conditional-locus discipline**: four endpoint chapters carry $\hyp...$ tags at theorem statements.
8. **Cyclic-Hochschild discipline**: four-object remark in `hochschild.tex`.
9. **Bidirectional cross-volume back-edits**: Vol I, Vol III, mixed-HT-strings, igusa updated to mirror Vol II discipline.
10. **Build-system enforcement**: `make verify-licensing` passes; pre-commit hook runs the script.
11. **CG voice sweep applied**: 17-line voice table enforced in every rectified chapter.
12. **IV layer**: 8 disjoint-route IV tests pass; `make verify-independence` green.
13. **Master theorem inscribed**: a single Architectural Definition states the universal chain in `chapters/frame/architectural_definition.tex`.
14. **Heptagon face / licensing-type duality remark**: stated in `sc_chtop_heptagon.tex` as the unification.
15. **Reflexive note carried**: this map is a recipe, expires when the discipline is in the prose.

The cut purchases: a programme whose architectural commitments are visible at every layer (preamble macro, theorem statement, chapter prose, AP catalogue, IV test, build-system gate, cross-volume manifest), whose 17 historical false slogans are mechanically caught, whose 4 construction problems are named, and whose mathematical content is *unchanged but more cleanly licensed*.

The cut does **not** purchase: new theorems. It is a compression, not an expansion. The new theorems remain the work of the manuscripts; the cut creates the architectural surface on which they will be written.

---

## Appendix A — Diff against v1

| v1 section                                | v2 disposition                                                                                                                                                                                  |
| ----------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| §0 Architectural fact                      | Promoted to §1 Master Architectural Theorem with explicit universal stage chain                                                                                                                |
| §1 Five Beilinson cuts                     | Replaced by §3 Five licensing types (chart / comparison / ambient / endpoint / effectiveness); the "cuts" are now *one cut* indexed by type                                                    |
| §2 17-dismissal cascade table              | Refined as §3 closing table with stage-collapse column + licensing-type column                                                                                                                  |
| §3 Cross-volume propagation matrix         | Retained, refined as §8 with bidirectionality                                                                                                                                                   |
| §4 Macro layer                             | Refined; bicoloured primitive package macro added (§6); `\hypothesis` parametric                                                                                                                |
| §5 AP + cache rows                         | Retained; V2-AP-SHADOW-EQUALS-OBJECT demoted from AP to *licensing rule*                                                                                                                        |
| §6 Compute / IV                            | Retained, refined                                                                                                                                                                               |
| §7 Voice / CG-rectification                | Refined as §10 with conditional-licensing tags                                                                                                                                                  |
| §8 Phasing                                 | Retained, honest re-budgeting (25-35 sessions); Phases 9-11 added                                                                                                                               |
| §9 Already-inscribed                       | Retained as part of §13 audit                                                                                                                                                                   |
| §10 Reconstitution invariant               | Refined as §14 What the cut purchases                                                                                                                                                            |
| Appendix A one-line recipes                | Retained as part of §10 voice table                                                                                                                                                              |
| Appendix B status flag                     | Refined as §14 invariants                                                                                                                                                                        |
| (new) §2 Two lanes parallel                | New: open / CY / trace lanes via heptagon face (6)                                                                                                                                              |
| (new) §4 Heptagon ↔ licensing duality       | New: deep symmetry articulating heptagon faces and licensing types                                                                                                                              |
| (new) §5 Cyclic Hochschild discipline       | New: four-object discipline at level $\mathsf{Z}$                                                                                                                                                |
| (new) §6 Closed-colour primitive package    | New: bicoloured nine-tuple replacing open-only seven-tuple                                                                                                                                      |
| (new) §7 Four construction problems         | New: $\mathfrak{D}_X$, gravity-line, PVA-quantum, chiral Positselski                                                                                                                            |
| (new) §9 Build-system layer                 | New: `make verify-licensing`                                                                                                                                                                    |
| (new) §11 Licensing decision tree           | New: cheat sheet for writers                                                                                                                                                                    |
| (new) §13 Adversarial audit of v2           | New: reflexive cut                                                                                                                                                                              |

## Appendix B — Hidden symmetries v2 makes explicit (the "yearned-for" inner form)

1. The two structural chains the critique closes with are **two presentations of one universal chain**, unified at level $\mathsf{Z}$ by SC$^{\mathrm{ch,top}}$ heptagon face (6). (§2)
2. The heptagon's seven faces and the licensing-type taxonomy's five types are **two presentations of one architectural commitment** (operadic / categorical vs epistemic / proof-discipline). (§4)
3. The 17 dismissals are **17 instances of one cross-stage collapse** along the universal chain, parametrised by lane × stage-pair × licensing-type. (§3 table)
4. The two-stage CY-chiral factorisation $\Phi_d = \SpCh \circ \PhiFA_d$ and the open-lane chart-projection $\openFactCat \to A_b = \mathrm{End}(b)$ are **two instances of the same chart functor** (§2 table) at different lanes.
5. The Pentagon trace, ghost trace, and Borcherds weight three-factor identity is **the level-$\mathsf{Z}$ closed-colour trace** of the trace lane, parallel to the open-lane derived centre and the CY-lane Drinfeld double. (§2)
6. The four construction problems (§7) are **the four open level-$\mathsf{C}\to\mathsf{Z}$ arrows** that current shadows lack their operator-level domain. They are coherent — collectively, they are the operator-level reconstitution of the protected-shadow content of the five programmes.
7. The cyclic-Hochschild four-object discipline at level $\mathsf{Z}$ (§5) is the **finer presentation of the bulk** that current Vol II prose treats as a single object.
8. The bicoloured primitive package (§6) is the **finer presentation of the primitive** that v1's Cut A treated only as an open-colour seven-tuple. Modularity (dismissal 5) is now naturally a closed-colour trace evaluation, not a "property of the closed algebra."

These eight inner symmetries are not new mathematics. They are the *form that the existing mathematics yearns to be presented in*. The Beilinson cut is what makes them visible. v2's job is to ensure the rectification work brings them to the surface; v1's job was to enumerate the consequences. v2 is the deeper move.

---

## End-state

The map is now strongest in the following operational sense: every consequence in v1 is preserved, every hidden symmetry the subject was carrying is named, every adversarial audit of the map itself is run with explicit gaps named, the construction problems are inscribed, the build-system layer is scoped, the phasing is honest, and the reflexive note is present.

The next move is to start Phase 1 (macro + script scaffolding). The map is ready to expire.
