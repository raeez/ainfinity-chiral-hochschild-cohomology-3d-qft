# FRONTIER — Vol II Open Frontiers

The platonic-ideal frontier of Vol II: open mathematical problems whose resolution would advance the $\mathsf{SC}^{\mathrm{ch,top}}$ programme, the topologisation ladder, the bulk–boundary pair $(A, Z^{\mathrm{der}}_{\mathrm{ch}}(A))$, the Universal Holography functor, or the chiral lane of the five theorems A, B, C, D, H. Each frontier carries a hypothesis package, a reconstruction theorem, and a heal path.

**Architectural anchors.** `notes/legacy/critique_2026_05_09_chiral_duality_master_consequence_map_v2.md` (universal stage chain $\mathsf{P}\to\mathsf{C}\to\mathsf{S}\to\mathsf{Z}\to\mathsf{A}$, five licensing types α/β/γ/δ/ε, four Construction Problems, 17-line voice table) · `notes/legacy/vol2_platonic_architecture.md` (the seven-part form Vol II yearns to be) · `notes/_METACOGNITIVE_CANON.md` (live vs archival index). Drafting-history cascades, session-dated snapshots, and historical wave-N stratification live in `archive/2026-04/FRONTIER_layered.md`; the manuscript carries the current state.

The four open Construction Problems named by the architectural cut (operator-level constructions whose scalar / shadow / classical / injection-only forms are currently the point of attack) recur across the frontiers below: $\mathfrak{D}_X$ for K3$\times E$ with $\mathrm{Pf}_{\mathrm{prot}}(\mathfrak{D}_X) = \Delta_5$ (igusa programme); the gravity-line operator algebra with Pentagon-face scalar trace $= \Phi_{10}^{\mathrm{un}}$ (F1, `3d_gravity.tex:8429`); the unified PVA-quantum HT theory with classical $\lambda$-Jacobi limit and $E_3$-lift on $Q$-cohomology (F8 layers + `pva-descent.tex`); and the chiral Positselski theorem extending Vol I Theorem B at chiral generality (F1, F6, F7).

---

## The proved core

| Theorem | Statement | Stage | Status |
|:-------:|-----------|:----:|--------|
| **A** | bar–cobar inversion $\Omega(\mathrm{Bar}(A)) \xrightarrow{\sim} A$ | $\mathsf{C} \leftrightarrow \mathsf{S}$ | ProvedHere on the chiral lane (Quillen / Lefèvre-Hasegawa / Loday–Vallette in $\mathrm{ChirAlg}^{E_1}$) |
| **B** | $\Omega(\mathrm{Bar}(C)) \xrightarrow{\sim} C$ in $D^{\mathrm{co}}_{\mathrm{ch}}(X)$ | $\mathsf{C} \to \mathsf{Z}$ | Conditional on Vol I's Theorem-B package: Koszul/collapse locus plus strict Mittag-Leffler completed chiral Positselski CP1-CP3; class M chain-level original-complex is false in $\mathrm{Ch}(\mathrm{Vect})$ and only lives in completed / pro / weight-completed / $J$-adic ambients |
| **C** | $\rho_K = \kappa_{\mathrm{ch}}^{\mathrm{Hodge}}(A) + \kappa_{\mathrm{ch}}^{\mathrm{Hodge}}(A^!) \in \{0, 8, 13, 250/3, 25/3\}$ on $\mathsf{G}/\mathsf{L}/\mathsf{C}/\mathsf{M}/\mathsf{B}$ five-archetype ceiling | $\mathsf{Z}$ (per stratum) | ProvedHere on Koszul locus + $\mathsf{B}$-row K3 Mukai-doubling witness; the Bershadsky–Polyakov value $25/3$ is conditional, pending the genus-one curvature computation of $\kappa_{\mathrm{BP}}$ (its central-charge conductor $c_{\mathrm{BP}}(k)+c_{\mathrm{BP}}(-k-6) \equiv 50$ is proved) |
| **D** | obstruction-tower universality: Hodge shadow $\mathrm{obs}^{\mathrm{Hdg}}_g = (-1)^g\,\kappa_{\mathrm{ch}}^{\mathrm{Hodge}} \cdot \lambda_g \in H^{2g}$ of the K-theoretic class $\mathfrak{O}_g^K = \kappa\,\lambda_{-1}(\mathbb{E}_g)$, uniform-weight, all $g \ge 1$; the deformation-theoretic object $\mathrm{Obs}^{\mathrm{def}}_g \in H^2(\mathrm{Def}_g)$ is distinct, bridged by clause $H_D^{\mathrm{tr}}$ | $\mathsf{Z}$ on $\overline{M}_{g,n}$ | Holds under $H_D^1$; K-theoretic identification under $H_D^K$; multi-weight extension carries cross-channel correction $\delta F_g^{\mathrm{cross}}$ |
| **H** | $\mathrm{ChirHoch}^\bullet \subset \{0, 1, 2\}$ on the PBW/Koszul generic finite-window locus | $\mathsf{Z}$ | Conditional on Vol I's Theorem-H package: PBW chiral Koszulness, $E_\infty$-chiral completion, finite-type/perfect diagonal Hochschild complexes, strict Mittag-Leffler passage, and ordered twisted-tensor acyclicity/averaging; Feigin-Frenkel critical level is outside the finite lane |

**Vol II-specific load-bearing closures.**

- $\mathsf{SC}^{\mathrm{ch,top}}$ heptagon, seven faces inscribed: bicoloured operadic primitive · quantum complementarity · bar–cobar adjoint inversion · brace / chiral Deligne–Tamarkin · topologisation ladder · Drinfeld-centre $=$ derived-centre $\bigl(Z(\mathrm{Rep}_{\mathrm{fact}}(A)) \simeq \mathrm{Rep}_{\mathrm{fact}}(Z^{\mathrm{der}}_{\mathrm{ch}}(A))^{E_2}\bigr)$ · PTVV derived-AG.
- Chiral Higher Deligne (`thm:chiral-higher-deligne`): $Z^{\mathrm{der}}_{\mathrm{ch}}(A)$ is the universal $E_3$-chiral acted on by $\mathsf{SC}^{\mathrm{ch,top}}$ only under the stated two-coloured cobar and topologisation hypotheses. The route from $E_3$-rigidity + PBW collapse to Theorem H is conjectural and requires the derived-centre $E_3$-PBW package; Vol I's Theorem H route is the conditional ordered-bar / twisted-tensor acyclicity / averaging package. `thm:chd-ds-hochschild` is conditional on the completed DS bar-coalgebra SDR, the induced coderivation transfer, and bounded-shift HPL convergence; with that package it gives $\mathrm{ChirHoch}^\bullet(W_k(\mathfrak{g})) \simeq H^\bullet_{\mathrm{DS}}(\mathrm{ChirHoch}^\bullet(V_k(\mathfrak{g})))$ as an $E_2$-chiral Gerstenhaber comparison.
- Curved-Dunn $H^2 = 0$ at $g \ge 2$: `thm:curved-dunn-H2-vanishing-all-genera`. Combined with `thm:irregular-singular-kzb-regularity` (Jimbo–Miwa) replacing KZ Stokes, closes modular operad composition at generic non-integral level for all $g \ge 1$.
- $E_\infty^{\mathrm{top}}$ topologisation ladder (`thm:e-infinity-topologization-ladder`): $k$ inner stress tensors $T^{(n)} = [Q_{\mathrm{tot}}, G^{(n)}]$ promote $E_1$ to $E_{k+2}^{\mathrm{top}}$ on $Q$-cohomology. Affine Kac–Moody at non-critical level $\to E_3^{\mathrm{top}}$. $\mathcal{W}_N \to E_{N+1}^{\mathrm{top}}$. $\mathcal{W}_\infty[\lambda] \to E_\infty^{\mathrm{top}}$ Platonic endpoint conditional on $\hypProchazka$, $\hypCKL$, $\hypPRSh$, $\hypYamada$ + antighost-commutativity axiom.
- Universal celestial holography (`thm:uch-main`): $\mathsf{SC}^{\mathrm{ch,top}}$ on $(A^{\mathrm{cel}}, Z^{\mathrm{der}}_{\mathrm{ch}}(A^{\mathrm{cel}}))$ + celestial OPE $=$ chiral factorisation homology on $\mathbb{P}^1_{\mathrm{cel}}$. Coverage: self-dual gauge (KM), gauge + matter (DS), gravity (Virasoro $+$ $w_{1+\infty}$), YM (Beem–Rastelli $\chi$-functor). Class M $g \ge 1$ chain-level retained as `conj:uch-gravity-chain-level`.
- Unified $Q^{k, f, \mu}_g$ chiral quantum group (`chapters/theory/unified_chiral_quantum_group.tex`): nine specialisation fibres cover the non-gauge-theoretic landscape.
- Principal $\mathcal{W}_N$ tempering (`thm:wn-tempered-all-N`):
  finite-\(\beta_N\) Stirling dominance is proved for every finite
  principal rank.  The harmonic identification
  `conj:beta-N-harmonic-closed-form`,
  \(\beta_N=12(H_N-1)=\sum_{s=2}^{N}12/s\), is not a closure; it is the
  open \(\kappa\)-ratio scaling-law input.  Conditional on that input,
  \(\beta_2,\beta_3,\beta_4\) are integral and \(\beta_N\) is
  non-integral for every \(N\ge5\).
- MC5 chain-level class M weight-completed (`prop:bv-bar-class-m-weight-completed`): direct-sum class M chain-level genuinely false in $\mathrm{Ch}(\mathrm{Vect})$ (Vol I `thm:mc5-class-m-chain-level-pro-ambient` shows the same on pro-object and $J$-adic ambients).

**Universal Holography master theorem.** $\Phi_{\mathrm{hol}}: \mathsf{SC}^{\mathrm{ch,top}}\text{-Alg}_A \to \text{HT-pair}(A, Z^{\mathrm{der}}_{\mathrm{ch}}(A))$. For $A = \mathrm{Vir}_c$ in the Brown–Henneaux dictionary $c = 3\ell / (2 G_N)$, the master theorem identifies the algebraic holographic HT sector of pure 3D quantum gravity. The theorem does not construct a dynamical-metric path integral; the BTZ / Cardy promotion requires modular-invariance, vacuum-dominance, and saddle-dominance hypotheses, declared at every site that crosses from algebraic identification to physical interpretation.

**Cross-volume bridges.**

- The Vol III chain fusion conjecture $A_X = A_{b(X, \Sigma, C)}$ identifies Vol III's Stage-2 chiral output on a curve $C$ with Vol II's chart algebra $A_b$ in an open factorisation dg-category on $(C, D_C, \tau_C)$. Verified at $\mathbb{C}^3$, local $\mathbb{P}^2$, conifold, $K3 \times E$. The conjecture is the bridge from Vol III's CY input to Vol II's centre / universal-holography output.
- The Universal Trace Identity $\mathrm{tr}_{\mathrm{ghost}}(Q_{\mathrm{BRST}}^2) = \mathrm{tr}_{\mathrm{Pentagon}} = \omega_{\mathrm{Borcherds}} = c_N(0)/2$ at $N \in \{1, 2, 3, 4, 6\}$: ghost trace on the BRST anomaly $=$ Pentagon trace (Vol II's single-colour coherence trace of the $E_3$-algebra at $d = 3$) $=$ Borcherds weight (Vol III's $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$). At $N = 1$ the three readings coincide on the K3-Heisenberg witness: $\{5, 5, 5\}$.
- Vol II's BKM-side K3 quantum-group lane targets the Hall–Drinfeld double $\mathbf{H}_{\Delta_5} = \mathcal{D}_\hbar(\mathrm{CoHA}_{K3 \times E})$ on the bi-based $\mathrm{Ran}(E) / \overline{\mathcal{A}_2}$ datum with CY-2 $[2]$ shift, class-$\mathcal{S}$ parent $\mathcal{T}[A_1, \Sigma_{0, 24}]$, and $\Delta_5$ as 1-loop-forced output. *K3 Yangian* without the Hall–Drinfeld qualifier denotes the separate Mukai self-mirror branch $Y_\hbar(\mathfrak{so}(4 \mid 20))$, not the BKM lane.

---

## Open frontiers

### F1 — BV/BRST $=$ bar at chain level for class $\mathsf{M}$, $g \ge 2$

Class $\mathsf{M}$ (Virasoro at $c \ne 25, 26$; principal $\mathcal{W}_N$; Bershadsky–Polyakov; non-trivial-cocycle algebras) carries chain-level coderived equivalence $\mathrm{BV} \simeq \mathrm{Bar}$ in $D^{\mathrm{co}}(A)$ at $g \ge 2$ via MC5 weight-completed pattern + chiral Higher Deligne; the original-complex direct-sum form remains chain-level genuinely false in $\mathrm{Ch}(\mathrm{Vect})$. The remaining frontier is explicit chain-level $A_\infty$ coherence at $g \ge 2$ for non-formal class $\mathsf{M}$, and the conceptual identification of $D^{\mathrm{co}}$ as the home where curvature is absorbed by the coderived differential.

- **Hypothesis package**: chiral Higher Deligne + curved-Dunn $H^2 = 0$ + MC5 weight-completed; Mittag-Leffler antighost commutativity; class M ambient declaration.
- **Reconstruction theorem**: explicit $A_\infty$-quasi-isomorphism $\mathrm{BV}(A) \xrightarrow{\sim} \mathrm{Bar}(A)$ on $\overline{M}_{g, n}$ at $g \ge 2$ in the weight-completed / pro / $J$-adic ambient, with curvature $m_0 = \kappa_{\mathrm{ch}}^{\mathrm{Hodge}} \cdot \omega_g$ absorbing the harmonic discrepancy via Fay trisecant.
- **Heal path**: explicit coacyclicity at genus 1 for $\mathrm{Vir}_c$ at $c \in \{1, 25, 26\}$; extend the Felder BRST + Zamolodchikov norm computation to $g = 2$.

### F2 — The $(3, 2)$ nilpotent in $\mathfrak{sl}_5$: gateway to non-principal DS–KD

Drinfeld–Sokolov reduction commutes with bar–cobar in the principal case (`thm:ds-bar-cobar-principal`); for non-principal nilpotents the proof breaks where $\mathfrak{n}_+$ is non-abelian. The $(3, 2)$ partition of 5 is the first non-abelian case in type A (dim$\,\mathfrak{n}_+ = 8$, 2-step nilpotent, four nonzero commutators). The ghost-ghost BRST corrections $Q_{\mathrm{gh}} \ne 0$ defeat the Kazhdan-filtration argument; $E_1$-degeneration is the obstruction.

- **Hypothesis package**: 2-step nilpotent in type A; Kazhdan-filtration grading; ghost-number decomposition.
- **Reconstruction theorem**: $E_1$-degeneration at $(3, 2)$ implies DS–KD for all 2-step nilpotents in type A; failure mode identifies the precise obstruction to the non-principal direction.
- **Heal path**: the BRST complex has matrix sizes $\le 3000 \times 3000$ at the hardest weight; computationally accessible by ghost-number decomposition (17 sectors). Build `brst_sl5_subregular_engine.py`.

### F3 — Genus-5 cross-channel: the Borel-determining computation

The cross-channel generating function $\delta F_g^{\mathrm{cross}}$ on multi-weight algebras (e.g., $\mathcal{W}_3$) grows factorially in $g$ while the scalar tower $F_g^{\mathrm{scal}} = \kappa_{\mathrm{ch}}^{\mathrm{Hodge}} \cdot \lambda_g^{\mathrm{FP}}$ is Gevrey-0. Three data points $(g = 2, 3, 4)$ yield $A_{\mathrm{cross}}/A_{\mathrm{scalar}} \in [1.7, 3.1]$; the genus-5 datum determines the Gevrey shift $b$, hence $A_{\mathrm{cross}}$ uniquely. The genus-5 verification is also a closed-form test of CohFT-weighted topological recursion on the $A_2$ Frobenius manifold.

- **Hypothesis package**: stable graphs at $g = 5$ ($\sim 4000$–$5000$); rational arithmetic Newton interpolation at $\sim 12$ integer $c$ values; CohFT-weighted topological recursion structure on the $A_2$ Frobenius.
- **Reconstruction theorem**: the symbolic Givental–Stokes extraction on the $A_2$ Frobenius manifold yields $A_{\mathrm{cross}}(c)$ in closed form, by Givental–Teleman (DBSS 2019 / DYZ 2019 / Teleman / FSZ10 / PPZ19).
- **Heal path**: the analytic route suffices; genus-5 graph enumeration is verification, not determination. Carry out the symbolic Givental–Stokes extraction.

### F4 — Admissible $\mathfrak{sl}_3$ Koszulness: low-hanging fruit via periodic CDG

For admissible levels of $\widehat{\mathfrak{sl}}_2$ the simple quotient $L_k(\mathfrak{sl}_2)$ is Koszul (single-weight null-vector argument + Kac–Wakimoto). For $\mathfrak{sl}_3$ the null-vector ideal carries multi-weight contributions (highest root $\theta$ at grade $(p - 2)q$, simple roots $\alpha_1, \alpha_2$ at grade $(p - 1)q$); the single-generator argument fails. Universal $V_k(\mathfrak{sl}_3)$ is Koszul at all levels (`prop:pbw-universality`); $L_k(\mathfrak{sl}_3)$ Koszulness at admissible levels is open at rank $\ge 2$.

- **Hypothesis package**: admissible level $k = -3 + p/q$, $\gcd(p, q) = 1$, $p \ge 2$; periodic CDG filtration tooling (Vol I `chapters/theory/periodic_cdg_admissible.tex`); Arakawa 2015 *J. AMS* 22:1 Theorem 5.4.4 + Corollary 5.4.5: $L_k(\mathfrak{sl}_n)$ is $C_2$-cofinite iff $q \ge h^\vee = n$, so the lisse branch for $\mathfrak{sl}_3$ is $q \ge 3$.
- **Reconstruction theorem**: explicit Li-bar $E_2$ degeneration at $k = -5/3$ ($p = 4$, $q = 3$), the first admissible **lisse** level for $\widehat{\mathfrak{sl}}_3$ where nulls enter the bar range ($h_\theta = 6$, $h_{\alpha_i} = 9$, $c = -10$); $R_{L_k}$ finite-dimensional Artinian, dimension bounded by the structural decomposition $d_C^2 \cdot d_R^6 = 1024$. Three independent routes: Li-bar via ghost-number filtration + Lie–Poisson surjectivity; Adamović–Milas screening with quantum Serre at $q_{qg} = e^{3\pi i/4}$; GKO coset $L_{-5/3}(\mathfrak{sl}_3) = \mathrm{Com}(V_{\mathrm{lat}}, V_c(\mathrm{Vir}) \otimes \mathcal{W}_3)$.
- **Beilinson cut.** The level $k = -3/2$ ($p = 3$, $q = 2$) is **quasi-lisse, not lisse**: $X_{L_{-3/2}(\mathfrak{sl}_3)} = \overline{\mathcal{O}_{\min}} \subset \mathfrak{sl}_3^*$ is 4-complex-dimensional (Adamović–Perse 2014, Beem–Rastelli 2018), and $R_{L_{-3/2}}$ is infinite-dimensional. The five quasi-lisse admissible levels $k \in \{-3/2, -1/2, 1/2, 3/2, 5/2\}$ (all at $q = 2$) fall outside the structural Kunneth argument; Koszulness analysis there requires the Beem–Rastelli sheaf-of-vertex-algebras framework on $\overline{\mathcal{O}_{\min}}$.
- **Heal path**: engines `admissible_koszul_rank2_engine.py` and `theorem_admissible_sl3_libar_engine.py` (Vol I); `compute/lib/f4_admissible_sl3_libar_vol2.py` (Vol II Frontier sweep, 34 admissible $(p,q)$ verified: 29 lisse-Koszul, 5 quasi-lisse open).

### F5 — Restricted DK-4 on the evaluation-generated core

Drinfeld–Kohno extends from finite-dim representations to the formal moduli of line operators in 3D HT theory. MC3 is proved on the evaluation-generated core for all simple types (`thm:categorical-cg-all-types`); the type-A reduction (`prop:yangian-dk4-typea-frontier`) reduces DK-4 to a single mixed-tensor coefficient identity, satisfied on the factorisation side. The pointwise data is confirmed for $\mathfrak{sl}_2$ through $\mathfrak{sl}_8$. The missing step is the pointwise-to-global passage: proving $\mathfrak{g}_A = Y^{\mathrm{dg}}_A$ as a filtered complete dg Lie algebra.

- **Hypothesis package**: evaluation-generated core; pointwise Ext data; type A; Mittag-Leffler tower completion.
- **Reconstruction theorem**: $\mathrm{Fact}_{E_1}(X; A) \xrightarrow{\sim} \mathrm{Mod}^{\mathrm{comp}}(Y^{\mathrm{dg}}_A)$ on the evaluation-generated core.
- **Heal path**: extend existing engines to compute $\mathrm{Ext}^*(V_\omega(a), V_\omega(b))$ for $\mathfrak{sl}_3$ (first rank-2 case); compute the degree-2 seed comparison.

### F6 — DK-5: categorical $E_1$ primacy

The full triple bridge $\mathrm{Fact}_{E_1}(X; A) \simeq \mathrm{Mod}^{\mathrm{comp}}(Y^{\mathrm{dg}}_A) \simeq \mathrm{Rep}^{\mathrm{spec}}(\mathrm{QG}^{\mathrm{spec}}(R_A))^{\mathrm{op}}$ unifies factorisation algebras of local operators, dg-shifted Yangian module categories, and spectral quantum group representations. Bridge Criterion Theorem (`thm:bridge-criterion`): B1 + B2 + B4 $\Rightarrow$ full bridge. B1 is full $\mathcal{O}$-Koszulness beyond the eval core; B2 is tower completion (Mittag-Leffler proved, algebraic identification open); B4 is the spectral quantum group comparison with Latyntsev.

- **Hypothesis package**: B1 full $\mathcal{O}$-Koszulness; B2 Mittag-Leffler + algebraic identification; B4 Latyntsev spectral-QG.
- **Reconstruction theorem**: full triple bridge, equivalence of model categories.
- **Heal path**: B2 algebraic identification via reduced-resolution Yangian + filtered comparison; B4 via Latyntsev's spectral construction. B1 via beyond-eval Vermas + dual modules.

### F7 — The Grand Completion

The modular cumulant transform packages the full bar–cobar machine — modular MC element, genus tower, shadow obstruction tower, $R$-matrix — into the completed pronilpotent modular cumulant coalgebra, equivalent to the chiral algebra up to homotopy. Two sub-conjectures: (a) cumulant recognition: the resonance-graded associated graded of the completed bar is the cofree coalgebra on primitive cumulants. (b) Drinfeld spectral filtration: reduced-weight-$q$ bar windows determine the Yangian $r$-matrix through jet order $z^{-(q-1)}$ modulo the regular and lower-order jets already determined by $K_{<q}$. The AP19 d-log shift identifies OPE pole order $n$ with $r$-matrix jet order $z^{-(n-1)}$; the naive reading "$K_q \to z^{-q}$" is off by one. The corrected principle holds vacuously on K3-Heisenberg (Gaussian, $r_{\max} = 2$, $K_2 \to z^{-1}$ matches the Mukai Casimir $\Omega_{\mathrm{Muk}} = G^{II}_{ij} h^i \otimes h^j$) but fails in the naive form on Virasoro (Class M, $r_{\max} = \infty$, OPE poles at orders $4, 2, 1$ with $r$-matrix jets $z^{-3}, z^{-1}, z^0$ — the cubic shadow at $K_3$ does not lift to a $z^{-2}$ jet).

- **Hypothesis package**: completed pronilpotent modular cumulant coalgebra; resonance-graded filtration; primitive cumulant generators; jet-order weight filtration.
- **Reconstruction theorem**: equivalence of model categories extending the proved genus-0 Quillen equivalence.
- **Heal path**: very hard; the principal open structural problem of Vol II. Sub-conjecture (b) is potentially accessible via explicit computation on the K3-Heisenberg.

### F8 — Analytic realisation: three-layer gap

A vertex algebra is an algebraic skeleton; the physical theory requires convergent correlation functions, partition functions, sewing amplitudes. HS-sewing for the standard landscape (`thm:general-hs-sewing`), Heisenberg sewing (`thm:heisenberg-sewing`), lattice sewing (`thm:lattice-sewing`) are proved; three layers remain.

- **Layer 1 (sewing envelope for interacting algebras)**: sl$_2$ at $k = 1$ closed only algebraically (FKS); vertex operators are not Hilbert–Schmidt on a single Fock; requires $\sim 3$ functional-analytic theorems.
- **Layer 2 (metric independence on the conformally flat 2-disk)**: closed only for Heisenberg via Moriwaki 2026b. Anomaly cancellation at chain level open.
- **Layer 3 (higher-genus coderived shadow)**: downstream of layers 1–2.

Heal path: layer 1 via explicit kernel computation on the Segal–Banach completion of the Fock module at sl$_2$, $k = 1$; layer 2 via heat-kernel regularisation on the conformally flat disk.

### F9 — $E_1$ Verdier on ordered configurations

Verdier duality on $\mathrm{Ran}(X)$ intertwines $\mathrm{Bar}(A)$ and $\mathrm{Bar}(A^!)$ — the algebraic incarnation of $S$-duality in the HT theory. The ordered bar $B^{\mathrm{ord}}$ lives on $\mathrm{Conf}^<(X)$, not $\mathrm{Ran}(X)$; pushforward to $\mathrm{Ran}$ loses the ordering. The two-colour double Koszul duality theorem (`thm:two-color-master`) confirms: closed colour uses Verdier / Ran; open colour uses linear duality. The intertwining $B^{\mathrm{ord}}(A^{\mathrm{op}}) = B^{\mathrm{ord}}(A)^{\mathrm{cop}}$ proves the open-colour Quillen equivalence at the chiral-algebraic level.

- **Hypothesis package**: $D_{\mathrm{Conf}^<}$ Verdier duality on ordered configuration spaces; or a ribbon Ran space.
- **Reconstruction theorem**: open-colour $E_1$ Verdier with full six-functor package.
- **Chiral-algebraic content (already proved).** `thm:opposite` in `chapters/connections/ordered_associative_chiral_kd_core.tex:451` supplies the algebraic shadow of the literal Verdier identity at the level of bar coalgebras: $R_A: B^{\mathrm{ord}}(A^{\mathrm{op}}) \xrightarrow{\sim} B^{\mathrm{ord}}(A)^{\mathrm{cop}}$, induced by the order-reversal involution on $\mathrm{Conf}^<_n(\R)$, with the Koszul sign $(-1)^{\sum_{i < j}(|a_i| + 1)(|a_j| + 1)}$ absorbing the parity. Verified at finite degree by `compute/lib/f9_verdier_ordered_engine.py` (face-bijection, coproduct intertwining $\Delta^{\mathrm{cop}} \circ R_A = (R_A \otimes R_A) \circ \Delta$, involutivity $R_A \circ R_A = \mathrm{id}$ on words up to length $6$--$7$ over a three-letter alphabet).
- **Three obstructions to the six-functor upgrade**: (G1) general curves admit no global $\mathrm{Conf}^<$, only ordered configurations per open chart glued along monodromy; (G2) the six-functor package on $\mathrm{Conf}^<$ has no published construction (Costello--Gwilliam 2021 Vol II §7.3 builds ordered prefactorization algebras, Lurie HA §5.5 builds factorization on stratified spaces with $(f^*, f_*, f_!, f^!)$, but neither extracts a Verdier $D$ on $\mathrm{Conf}^<$); (G3) ribbon-Ran exists for $\R, S^1, \mathbb{H}$ via BD04 §3.3 / Kontsevich--Soibelman bordered Swiss-cheese / Horel, not for general algebraic curves.
- **Heal path**: two candidate routes for the six-functor upgrade. (a) Mok 2025 logarithmic Fulton--MacPherson compactification on a smooth curve combined with a cyclic-operad orientation gives a ribbon-FM, then push to ribbon Ran by a colimit over symmetric subgroups inducing cyclic orderings. (b) Choose a real structure / antiholomorphic involution on $C$ producing a real Ran subspace where ordering is canonical near real points, allowing a partial $D_{\mathrm{Conf}^<}$ restricted to real-Ran loci.

### F10 — Resurgence: pin down $A_{\mathrm{cross}}$ from genus-5

The cross-channel instanton action $A_{\mathrm{cross}}$ controls the large-order behaviour of the multi-weight genus expansion; determines Borel summability and non-perturbative effects on the exact partition function. Scalar instanton action $A_{\mathrm{scalar}} = (2\pi)^2$ from the $\hat{A}$ genus; $A_{\mathrm{cross}}$ from the multi-weight OPE. Current bounds $A_{\mathrm{cross}} / A_{\mathrm{scalar}} \in [1.7, 3.1]$; F3 above resolves the Gevrey shift $b$, hence $A_{\mathrm{cross}}$ uniquely.

### F11 — Cross-channel generating function

No closed-form $\hat{A}$-type generating function exists for $\delta F_g^{\mathrm{cross}}$ (`prop:cross-channel-no-closed-form`, ProvedHere, five independent obstruction paths). Three obstructions: inhomogeneous $c$-scaling ($\mathcal{O}(1)$ at $g = 2$ vs $\mathcal{O}(c)$ for $g \ge 3$); super-linear ratio growth; irreducible numerators. If a generating function exists it must be bivariate in $(c, \hbar)$ and non-separable.

### F12 — Scalar saturation beyond algebraic families

Layer 1 ($\dim H^2_{\mathrm{cyc}} = 1$) for all algebraic families with rational OPE coefficients. Layer 2 ($\Gamma_A = \kappa_{\mathrm{ch}}^{\mathrm{Hodge}} \cdot \Lambda$) on the uniform-weight lane; fails for multi-weight at $g \ge 2$. Residual content: layer 1 for non-algebraic-family modular Koszul algebras. Three candidate families need checking — (1) non-GKO cosets, (2) 4D $\mathcal{N} = 2$ quiver VOAs, (3) admissible-level simple quotients at rank $\ge 2$. No counterexample known.

### F13 — Super-Yangian complementarity

$\kappa_{\mathrm{ch}}^{\mathrm{Hodge}}(Y(\mathfrak{sl}(m \mid n))) + \kappa_{\mathrm{ch}}^{\mathrm{Hodge}}(Y(\mathfrak{sl}(n \mid m))^!) = \max(m, n)$ at Sugawara-shifted dual level (symbolic small-rank verification); the Verdier pairing inscription depends on the convention reconciliation between super-trace and Berezinian, currently in flight on the cross-volume side.

- **Hypothesis package**: super-trace vs Berezinian convention; Nazarov quantum-Berezinian centrality.
- **Reconstruction theorem**: a single canonical Verdier pairing for $Y(\mathfrak{sl}(m \mid n))$ with super-trace and Berezinian readings related by a named comparison map.
- **Heal path**: Vol I `chapters/frame/programme_overview_platonic.tex` overview-vs-foundations reconciliation; bridge via Nazarov.

### F14 — $W(p)$ logarithmic triplet placement

The logarithmic triplet $W(p)$ ($p \ge 2$) sits outside the $\mathsf{G} / \mathsf{L} / \mathsf{C} / \mathsf{M}$ four-class partition. $p = 2$ closed trivially (symplectic-fermion $r_{\max} = 4$); the $TT$ sub-channel proved tempered for all $p \ge 2$ via Virasoro generic-$c$ at $c(p) = 1 - 6(p - 1)^2 / p$. The $TW$ / $WW$ sub-channels at $p \ge 3$ are genuinely open.

- **Hypothesis package**: extend the partition to a five-class $\mathsf{G} / \mathsf{L} / \mathsf{C} / \mathsf{M} / \mathsf{log}$ including unbounded-Massey LCFTs.
- **Reconstruction theorem**: replace Zhu-bounded-Massey with direct character-amplitude growth bound.
- **Heal path**: inscribe the Adamović–Milas amplitude bound; verify $|S_r(W(3))|$ for $r \le 8$ via Adamović–Milas $\phi_{0, 1}$.

### F15 — Arithmetic typed separation: $H^{\bullet}_{\mathrm{Den}}$ vs $H^{\bullet}_{\mathrm{ord}}(\operatorname{Spec}\mathbb{Z}, \mathcal{H}_{\mathrm{Ar}})$ vs $H^{\bullet}_{\partial}(\mathrm{ACS}_{\mathrm{qu}})$

The arithmetic branch of the programme — the ordered-incidence chiral candidate over $\operatorname{Spec}\mathbb{Z}$, the conjectural Deninger cohomology of the Riemann zeta, and the boundary $\partial$-homology of arithmetic Chern–Simons — has been written informally as a single object. The platonic-ideal form distinguishes three carriers that must not be conflated without naming the comparison morphism:

\begin{enumerate}[label=\textup{(\arabic*)}]
\item $H^{\bullet}_{\mathrm{ord}}(\operatorname{Spec}\mathbb{Z}, \mathcal{H}_{\mathrm{Ar}}) = H^{\bullet}(B^{\mathrm{ord}}_{\bullet}(V_{\mathrm{prim}}))$ — the \emph{ordered-incidence bar homology} on $\operatorname{Spec}\mathbb{Z}$ with arithmetic primitive coefficient bundle $V_{\mathrm{prim}}$. Construction: the ordered Fulton–MacPherson tower over $\operatorname{Spec}\mathbb{Z}$ with Bost–Connes finite Euler terms. \emph{Status:} this is the only carrier that admits a Vol II-style chain-level construction.
\item $H^{\bullet}_{\Delta, \chi}(\operatorname{Spec}\mathbb{Z}) = H^{\bullet}_{\mathrm{Den}, \chi}(\operatorname{Spec}\mathbb{Z})$ — the Deninger candidate cohomology obtained from (1) by adjoining completed continuous duals and formal Tate polar lines in degrees $0$ and $2$, then imposing the spectral-flow operator $\Theta_{\xi} = \tfrac{1}{2}\mathrm{id} + i H_{\mathrm{HP}}$ with determinant $\xi_{\mathbb{R}}$. \emph{Status:} target package; conditional on $\Xi_{\mathrm{BD/CG}}$ and $\Theta_{\xi}$.
\item $H^{\bullet}_{\partial}(\mathrm{ACS}_{\mathrm{qu}})$ — the boundary $\partial$-homology of the arithmetic Chern–Simons quantum theory of Kim–Chung–Kim–Park–Yoo. \emph{Status:} the natural physical target; requires a Vol IV / arithmetic-CS comparison.
\end{enumerate}

The forbidden slogan ``Deninger cohomology is constructed'' replaces these three carriers with a single object. The allowed form is the conditional packaging below.

- **Hypothesis package**: $\Xi_{\mathrm{BD/CG}}\colon H^{\bullet}_{\mathrm{ord}} \to \int^{\mathrm{arith\text{-}ch}}_{\operatorname{Spec}\mathbb{Z}} \mathcal{H}_{\mathrm{Ar}}$, the Beilinson–Drinfeld / Costello–Gwilliam factorization-homology comparison over $\operatorname{Spec}\mathbb{Z}$; $\Theta_{\xi} = \tfrac{1}{2}\mathrm{id} + i H_{\mathrm{HP}}$, the spectral-flow operator with determinant identity $\det \Theta_{\xi} = \xi_{\mathbb{R}}$; the Deninger trace formula $\sum_{\rho} h(\rho) = \mathrm{polar} + \mathrm{arch} - \sum_{p, k}(\log p)\, p^{-k/2}\, \hat{h}(k\log p)$.
- **Reconstruction theorem** (\ClaimStatusConditional): under $\Xi_{\mathrm{BD/CG}}$ + $\Theta_{\xi}$ + the trace formula, the three carriers fit into the diagram
\[
 H^{\bullet}_{\mathrm{ord}}(\operatorname{Spec}\mathbb{Z}, \mathcal{H}_{\mathrm{Ar}})
 \xrightarrow{\;\Xi_{\mathrm{BD/CG}}\;}
 H^{\bullet}_{\Delta, \chi}(\operatorname{Spec}\mathbb{Z})
 \xrightarrow{\;\partial\text{-comp}\;}
 H^{\bullet}_{\partial}(\mathrm{ACS}_{\mathrm{qu}}),
\]
the Deninger package; without $\Xi_{\mathrm{BD/CG}}$ + $\Theta_{\xi}$ + trace, $H^{\bullet}_{\Delta, \chi}$ is an ordered-incidence candidate plus polar completion, not constructed Deninger cohomology.
- **Heal path**: prove $\Xi_{\mathrm{BD/CG}}$ as a BD/CG quasi-isomorphism over $\operatorname{Spec}\mathbb{Z}$ using the chiral Hochschild local structure $\Zderch{A_{\mathrm{Ar}}}$ adapted to arithmetic factorization homology (typed master functor of `chapters/theory/foundations.tex` `thm:typed-boundary-holographic-realisation` applied to $A_{\mathrm{Ar}}$); construct $\Theta_{\xi}$ via the Connes–Deninger zero-side input on the boundary curve $\operatorname{Spec}\mathbb{Z}$. Cross-volume locus: Vol IV verification capstone or the dedicated arithmetic-chiral branch.

The Vol II contribution is the typed-separation discipline itself, propagated from `def:six-typed-carriers` in `chapters/theory/foundations.tex`: arithmetic objects are mapped through the same six carriers (C1)–(C6), with $H^{\bullet}_{\mathrm{ord}}$ at (C3) (chiral Hochschild candidate over $\operatorname{Spec}\mathbb{Z}$), $H^{\bullet}_{\Delta, \chi}$ at (C5) (the L-function shadow), and $H^{\bullet}_{\partial}(\mathrm{ACS}_{\mathrm{qu}})$ at (C6) (the gravity-line analogue for the arithmetic curve).

---

## Open Construction Problems

The architectural cut (`notes/legacy/critique_2026_05_09_chiral_duality_master_consequence_map_v2.md` §7) names four operator-level Construction Problems whose currently available statements inhabit only the scalar / shadow / classical / injection face of the universal stage chain $\mathsf{P}\to\mathsf{C}\to\mathsf{S}\to\mathsf{Z}\to\mathsf{A}$. Lifting each to the acting-object level $\mathsf{A}$ is the stated research target. The chapter-level inscription is `chapters/connections/programme_climax_platonic.tex` `prop:four-construction-problems` (also `rem:four-construction-problems`); the four problems and their cross-volume anchors:

- **(P1)** $\mathfrak{D}_X$ for $X = K3 \times E$ with $\mathrm{Pf}_{\mathrm{prot}}(\mathfrak{D}_X) = \Delta_5$ as a chain-level identity in the $\mathrm{ChirHoch}^\bullet$-valued cyclic complex on $\Lambda^{2,1}_{\mathrm{II}}$. Licensing: $\alpha$ + $\beta$ + $\gamma$ + $\epsilon$. Anchor: `~/igusa-cusp-form/main.tex:94--118` (Igusa 1962); Vol III K3$\times E$ programme; primary literature Borcherds 1995 *Invent. Math.* 120 and Gritsenko 1999 *Compositio Math.* 116.
- **(P2)** Gravity-line operator algebra: chain-level $\mathsf{SC}^{\mathrm{ch,top}}$-pair on K3$\times E$ with Pentagon-face scalar trace $\Phi_{10}^{\mathrm{un}} = \Delta_5^{2}$ and an explicit Hall–Borcherds intertwiner. Licensing: $\alpha$ + $\beta$ + $\delta$ + $\epsilon$. Anchor: `chapters/connections/3d_gravity.tex:8429` `conj:gravity-line-hall-borcherds-comparison`; Vol III six-route Hall–Borcherds bridge.
- **(P3)** Unified PVA-quantum HT theory with classical $\lambda$-Jacobi limit and $E_3$-lift on $Q_{\mathrm{tot}}$-cohomology, conditional on $\hypKZSDR$, $\hypStokes$, $\hypReflWts$, $\hypTLift$. Licensing: $\alpha$ + $\gamma$ + $\delta$ + $\epsilon$. Anchor: F8 layers above; `chapters/theory/pva-descent.tex` and the iterated-Sugawara ladder `thm:iterated-sugawara-construction`; primary literature Khan–Zeng (`KZ25`), Lurie HA ch. 5.
- **(P4)** Chiral Positselski theorem: $\Omega(\mathrm{Bar}^{\mathrm{ch}}(A)) \xrightarrow{\sim} A$ in $D^{\mathrm{co}}_{\mathrm{ch}}(C)$ for arbitrary chiral algebras on the open factorisation dg-category, specialising to ordinary $\to$ Vol I Theorem B and to quadratic chiral $\to$ Gui–Li–Zeng MC bijection. Licensing: $\alpha$ + $\beta$ + $\gamma$ + $\epsilon$ ($\effKoszul$). Anchor: `thm:hochschild-bridge-higher-genus`; primary literature Positselski 2011 SMF, Beilinson–Drinfeld 2004 ch. 3.

Solving any one of (P1)–(P4) collapses one face of the universal chain into a chain-level identity at the acting-object stage and is recorded as a ProvedHere upgrade in the proved-core table above.

The companion structural discipline is `rem:cyclic-hochschild-four-objects` in `chapters/connections/hochschild.tex` (cyclic-Hochschild four-object discipline: bulk $\mathrm{ChirHoch}^\bullet$ vs. shadow $\mathrm{ChirHoch}_\bullet$ vs. cyclic $\mathrm{HC}^{\mathrm{cyc}}$ vs. negative cyclic $\mathrm{HC}^-$; licensing $\gamma + \epsilon$).

**Status snapshot.** Per CP, the gap between the current proved face and the target acting-object construction is:

| CP | Current face (proved) | Target face (open) | Gap |
|:---:|----------------------|--------------------|-----|
| **P1** | $\Delta_5(\Omega) \in M_5(\mathrm{Sp}_4(\mathbb Z))$ scalar Borcherds product (Borcherds 1995, Gritsenko 1999); chain-level cyclic complex $\mathrm{ChirHoch}^\bullet$ on $\Lambda^{2,1}_{\mathrm{II}}$ identified (Vol II `thm:chirhoch-virasoro-concentration`) | $\mathfrak{D}_X$ on $X = K3 \times E$ with $\mathrm{Pf}_{\mathrm{prot}}(\mathfrak{D}_X) = \Delta_5$ as chain-level identity | construct the protected-Pfaffian generator and exhibit the $\mathsf{S}\to\mathsf{Z}$ stage transport |
| **P2** | $\Phi_{10}^{\mathrm{un}} = \Delta_5^{2}$ Borcherds product on $\mathrm{O}(2, 26)/\mathrm{O}(2)\times\mathrm{O}(26)$ (Gritsenko–Nikulin) carrying Pentagon-face scalar trace; six-route Hall–Borcherds bridge in Vol III | gravity-line operator algebra realising the bridge at chain level | construct the $\mathsf{SC}^{\mathrm{ch,top}}$-pair on K3$\times E$ and the explicit Hall–Borcherds intertwiner; verify $\effHCSQuartic$ |
| **P3** | classical $\lambda$-Jacobi PVA limit (`pva-descent.tex`); $E_3$-lift on $Q_{\mathrm{tot}}$-cohomology proved for affine KM at non-critical level (`thm:e-infinity-topologization-ladder`); iterated-Sugawara $E_{k+2}^{\mathrm{top}}$ ladder | unified PVA-quantum HT theory with both limits realised on a single chain-level object | install $\hypKZSDR + \hypStokes + \hypReflWts + \hypTLift$; close the all-loop quantum lift |
| **P4** | Vol I Theorem B on positive-energy / Koszul class (`thm:chiral-positselski-cross-volume`, listed above); MC bijection of Gui–Li–Zeng on quadratic specialisation | chiral Positselski on arbitrary chiral algebras with the open factorisation dg-category | refined $\hypAmbientWtCpl$ ambient covering non-positive-energy, log-decorated, non-rational chiral algebras |

Each row inhabits a single stage transport in $\mathsf{P}\to\mathsf{C}\to\mathsf{S}\to\mathsf{Z}\to\mathsf{A}$: P1 is $\mathsf{S} \to \mathsf{Z}$ (trace + CY lane); P2 is $\mathsf{S} \to \mathsf{A}$ (open lane via Pentagon trace); P3 is $\mathsf{Z} \to \mathsf{Z}$ inside the open lane (classical to quantum); P4 is $\mathsf{S} \to \mathsf{Z}$ on the open lane outside the Vol~I positive-energy ambient. The stratification is `thm:construction-problem-stage-stratification` in `chapters/connections/programme_climax_platonic.tex`.

---

## Cross-volume cross-references

- **Vol I** (`~/chiral-bar-cobar`): five-theorem core; $E_1$–$E_1$ operadic Koszul duality; tangential log curves; modular open-closed convolution; canonical $\kappa$ / $r(z)$ / $S_r$ in `chapters/examples/landscape_census.tex`.
- **Vol III** (`~/calabi-yau-quantum-groups`): CY-to-chiral functor $\Phi_d$, Hall–Drinfeld double $\mathbf{H}_{\Delta_5}$, BKM $\mathfrak{g}_{\Delta_5}$, $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$, 6d hCS quartic obstruction, three-tier $r_{\mathrm{CY}}$ structure.
- **Vol IV** (`~/chiral-bar-cobar-vol4`): verification capstone; pairs every Vol II ProvedHere with an independent-verification decorator and disjoint compute engine.
- **igusa-cusp-form** (`~/igusa-cusp-form`): $\Delta_5$ Borcherds denominator; terminal scalar shadow at level 4 of the universal stage chain. Disclaimer at `main.tex:96` carried through Vol II's BKM-side references.
- **mixed-holomorphic-topological-strings** (`~/mixed-holomorphic-topological-strings`): formal-Darboux stalk theorem; holomorphic de Rham obstruction. The level-1 physical face of the CY-side stage-1 chain.

---

## Reading guide

This frontier is the platonic-ideal mathematical state: open problems with hypothesis packages and heal paths. Closures and proved theorems live in their home chapters; this file refers to them, it does not re-prove. Drafting-history cascades, session-dated snapshots, wave-N stratification, and adversarial-audit logs live in `archive/2026-04/FRONTIER_layered.md` (preserved lossless per `~/ecosystem/INVARIANTS.md §IX`).

When a frontier closes — by proof, restriction, or honest downgrade — its entry moves to the proved-core table or to the discarded-conjecture appendix in the closing chapter. Wave numbers, drafting-history language, and session-dated annotations do not enter this document.
