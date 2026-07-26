# Vol II Anti-Patterns Catalogue

Author: Raeez Lorgat.

## Canonical mathematical values (locked audit anchors, 2026-04-21)

When any entry in this catalogue or the first-principles cache appears
to assert one value while another entry asserts a different value, the
registry below pins the canonical mathematical value. Older values are
retained only as historical audit context, explicitly flagged as
"wrong claim", "old draft", or "retracted".

| Quantity | Canonical mathematical value | Retracted earlier values | Lock |
|---|---|---|---|
| $(c_{4d}, c_{2d})$ of $\mathcal T[A_1, \Sigma_{0, 24}]$ | $(107/6, -214)$ (Wave 15 Gaiotto) | $(26, -312)$ (Wave 14 error); $(23/4, -69)$ (Wave 25 Gaiotto agent error) | WOV-2 |
| Trinion/tube count at $n=24$ | $(n_v, n_h) = (63, 88)$ | $(21, 27)$ (Wave 25 confusion of trinion/tube count with multiplicities) | Chacaltana--Distler 2010 §5.14 |
| Universal $c_{4d}$ formula, $A_1$ genus 0 | $c_{4d} = (5n - 13)/6$ | $(12(g-1) + 7n)/6$ (Wave 14 error — fails $\mathrm{SU}(2)$ $N_f = 4$) | Shapere--Tachikawa 2008 |
| Monster BKM Cartan rank | $2$ (on $\mathrm{II}_{1,1}$) | $26$ (Wave 16 confusion with Fake-Monster) | Borcherds 1992 *Invent Math* 109 |
| Fake-Monster Cartan rank | $26$ (on $\mathrm{II}_{25,1}$) | --- | Borcherds 1992 |
| K3-BKM Cartan rank | $3$ (on $\Lambda^{2,1}_{II}$) | --- | Gritsenko--Nikulin 1998 §3 |
| $c_3$ (Borcherds coefficient) | $-8$ in Bruinier reduced-class convention | $176256$ (= $p_{24}(5)$, Wave 16 error; also the Gritsenko--Nikulin Cartan-matrix value with conversion factor $-22032$) | WOV-5 |
| $\zeta(3,3,3,3)$ numerical value | $0.000295999\ldots$ (depth-4, weight 12) | $0.0028565$ (Wave 17 draft, 10× error) | Brown 2012 motivic basis |
| $c(1, 2, \pm 2)$ heterotic lift coefficient | $-2$ (Wave 22 Witten) | $+1$ (Wave 21 error) | $\eta^{18} \theta_1^2$ direct expansion + DMZ 2012 + holomorphic-anomaly (three-path) |
| Humbert $H_4$ arithmetic description | $(2, 2)$-isogeny quotient of $E_1 \times E_2$ with $\mathrm{End} \supset \mathbb{Z}[2i]$; monodromy order 2 | $\mathbb{Q}(\sqrt 2)$-RM locus (Wave 15 imprecise; that is $H_8$, not $H_4$) | van der Geer 1988 Ch IX |
| Theorem B scope (Koszul locus) | $\overline{\mathcal A_2} \setminus \bigcup_{n \,\mathrm{admissible}} H_n$ (all admissible Heegner divisors) | $\overline{\mathcal A_2} \setminus (H_1 \cup H_4)$ (Wave 15 narrow) | Wave 18 Beilinson tightening |
| K3-BKM Weyl denominator | $\Delta_5 = \mathrm{Grit}(\eta^9 \vartheta_1) \in S_5(K(1))$ (additive) | $\Phi_{12}$ attempted (Fake-Monster, $\mathrm{II}_{25,1}$, multiplicative — different object) | Gritsenko 1999 Thm 6.1 |
| Umbral $A_{N-1}$ labelling rule | $(N-1) \mid 24$; at $N = 6$ reanchor to Niemeier $6 D_4$ | $N \mid 24$ (Wave 18 error) | Wave 19 Gaiotto |
| Leech lattice simple root norm | $2$; $r_\lambda = (\lambda; 1, 1 - \lambda^2/2) = (\lambda; 1, -1)$ at $\lambda^2 = 4$ | $6$ (Wave 24 Nekrasov formula error, using $r_\lambda = (\lambda; e + (\lambda^2 - 2)/2\,f)$) | Conway 1983 *Proc R Soc Lond A* 384; SPLAG Ch 27 |
| Fake-Monster Weyl vector | $\rho = e = (0_\Lambda; 1, 0)$ lightlike, $\rho_\Lambda = 0$ | Sign-inconsistent Wave-23 draft | Borcherds 1992 Thm 10.1 |
| $\Phi_{12}$ automorphic home | Orthogonal Shimura variety $\mathcal{D}_{\mathrm{II}_{26,2}} = O(26, 2)^+ / (O(26) \times O(2))$, complex dim 26, singular weight $12 = (26 - 2)/2$ | Siegel $\mathrm{Sp}_{26}(\mathbb Z)$ (impossible, $\mathbb H_{26}$ has dim 351); Jacobi-in-26-variables (no torus direction) | Borcherds 1992 |
| Restriction sublattice for $\Phi_{12} \to \Phi_{10}$ | Primitive $\mathrm{II}_{2,2} \hookrightarrow \mathrm{II}_{25,1}$, giving $O^+(\mathrm{II}_{2,2}) \cong \mathrm{Sp}_4(\mathbb Z)$ | $\mathrm{II}_{2,1}$ (Wave 23 error — sig $(2,1)$ carries no holomorphic modular form) | Gritsenko--Nikulin 1998 |
| Master $L$-value identification | $\log Z^{(1)}_{\mathbf H_{\Delta_5}} = -\log \Delta_5 - \kappa_{\mathrm{BGS}} \cdot L'(0, \Delta_5, \mathrm{std}) + \log C$, $\kappa_{\mathrm{BGS}} = 24$ | $L'(0, \Delta_{10}, \mathrm{ad}^0)$ (Wave 24 Costello — CAP-reducible, conflates three $L$-functions); $L'(0, \mathrm{ad}^0 \rho_{\Delta_{12}})$ (Wave 25 Gaiotto agent alternative, incompatible regulator) | Bruinier--Kühn 2003 Thm 4.11; Yoshikawa 2004 Thm 5.7 |
| Bloch--Kato Selmer on representation of $\mathbf H_{\Delta_5}$ deformation | $\dim H^1_f(\mathbb Q, \mathrm{std}\,\rho_{\Delta_5}) = 1$ (paramodular cyclotomic Hida family tangent) | $\dim H^1_f(\mathrm{ad}^0 \rho_{\Delta_{10}}) = 1$ (Wave 24 — right number, wrong representation; the adjoint spinor is rigid, dim 0) | Pilloni 2011 + Urban 2011 + Poor--Yuen 2015 + Thorne 2020 |
| Bridgeland $\dim \mathrm{Stab}(K3 \times E)$ | $48 = \mathrm{rk}\,\mathcal N(K3 \times E) = 24 \cdot 2$ | $26$ (Wave 23 claimed "codim-0 slice" — actual Künneth image has codim $22 = \mathrm{rk}\,T_{K3}$) | Bridgeland 2007 Thm 1.2 |
| Codim of $\mathrm{Stab}^\Phi$ inside ambient | $22$ | $0$ (Wave 23 error) | via four independent paths |
| Rank of $u_{\zeta_8}^{\mathrm{tilt}}$-mod | $162 = 27 \cdot 6 = (\ell' - 1)^3 \cdot \lvert S_3 \rvert$ | --- | Andersen--Polo--Wen 1994 |
| Structure of rank-$162$ MTC | $(A_1)_{k=2}^{\otimes 3} \rtimes S_3$ ($S_3$-crossed braided fusion) | $(A_1)_{k=2}^{\otimes 3} \boxtimes \mathbb Z[S_3]$ (Wave 23 tensor claim — $\mathrm{Vec}_{S_3}$ non-modular) | Turaev 2000; ENO 2010 |
| Conway $V^{s\natural}$ central charge | $c = 12$ | $c = 24$ (Wave 23 inherited from Monster error) | Duncan 2007 *Duke Math J* 139 |
| Conway $V^{s\natural}$ $\Psi$-placement | $\Psi^{\mathrm{metap}}$-image on metaplectic $\overline{\mathcal A_2^{(2)}}$ branch | Bosonic $\Psi$-image with $(K, \hbar^2) = (2, -1/2)$ (Wave 19/23 — Leech has no hyperbolic plane, universal identity out of scope) | Scheithauer 2008 |
| "Four is all" citation chain | Scheithauer 2017 + Dittmann--Ma--Scheithauer 2021 + Scheithauer 2006 (three papers) | Scheithauer 2017 alone (Wave 23/24 incomplete) | Scope: GN-reflective signature-$(2, n)$ with $n \ge 3$ |
| Fifth Borcherds product outside GN-scope | $24 A_1$ Niemeier product (Borcherds 1995 *Invent Math* 120 §13), reflective automorphic on sig-$(2, 24)$ singular weight 12 but fails GN-reflectivity (divisor non-rational-quadratic-hyperplane components) | --- | Borcherds 1995 |
| Pentagon admissibility congruence variable | $D_n = (n - 3)/2 \pmod 4 \in \{0, 1\}$, equivalently $n \equiv 3, 5 \pmod 8$ | $n \not\equiv 0, 3 \pmod 4$ (loose restatement) | Eichler--Zagier 1985 Thm 9.1 |
| Humbert--Heegner admissible $n \in [3, 36]$ | $\{3, 5, 11, 13, 19, 21, 27, 29, 35\}$ (nine values) | All Padovan-positive $n$ (loose) | Eichler--Zagier 1985 Thm 9.3 + Gritsenko--Nikulin 1998 |
| First admissible non-vanishing $\phi^{(n)}$ | $\phi^{(5)} = -2 \cdot [\mathrm{gen}]^{\otimes 5}$ | --- | Gritsenko--Nikulin 1998 Table 2 |
| Padovan $d_n$ reference table ($n \le 12$) | $(d_3, \dots, d_{12}) = (1, 0, 1, 1, 1, 2, 2, 2, 3, 4)$ with $d_n = d_{n-2} + d_{n-3}$, seeds $(1, 0, 1)$ | Fibonacci recurrence | Brown 2012 *Ann Math* 175 Thm 1 |
| $\phi^{(n)} \ne 0$ on K3-Humbert iff | (i) $n \equiv 3, 5 \pmod 8$ AND (ii) $d_n > 0$ AND (iii) $D_n \le 1$ (polar cutoff) | Padovan-only (AP-V2-24 / V2-AP127) | Eichler--Zagier Thm 9.3; Gritsenko--Nikulin 1998 |
| $\Psi$-functor completeness | **Four** sibling functors $\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$ (disjoint union surjective) | $\Psi$ alone surjective (Wave 23 claim falsified) | Baily--Borel--Freitag stratification of $\overline{\mathcal A_2}$ |
| $\mathrm{grt}_1^{(1/2)}$ extension status | Non-split, non-abelian, non-central (three distinct properties) | "Split because central" (Wave 24 conflation) | Lie $H^2_{\mathrm{Lie}}(\mathfrak q; \mathrm{grt}_1)$ via van Est transgression |
| Obstruction cohomology type | Lie $H^2_{\mathrm{Lie}}$ via van Est transgression from group $H^1_{\mathrm{grp}}$ (Saito--Kurokawa Eichler cocycle) | Group $H^1$ only (Wave 24 misattribution) | Costello--Gwilliam 2017 BV realisation |
| $[\widetilde\sigma_2, \widetilde\sigma_2]$ | $288 \widetilde\sigma_4$ via $\tau(2)^2 / 2 = 288$ | --- | $\tau(2) = -24$ |
| $e_5$ vs $W_5$ | $e_5 = W_5$ identically at every $c$ | --- | Wang 1998 *Prog Theor Phys* Prop 4.2 three-leg uniqueness |
| $e_4$ at $c = -214$ | $:T \partial^2 T: - (3/2):(\partial T)^2: + (321/10) \partial^4 T + \hbar\mathrm{qt}(J^{(4)})$; pairing $(65193/10) \mathrm{Vol}(E) (2\pi i)^4$ | --- | Pope--Romans--Shen 1990; Wang 1998 |
| Chiral-Hochschild $e_k$ motivic home | Single-valued $\mathrm{zv}^{\mathrm{sv}}_{3k}$ (Brown 2014; Schnetz 2013) | Full motivic $\mathrm{Per}^{\mathrm{mot}}_{3k}$ (Wave 23 naive Costello cross-cut) | BGS analytic torsion on Shimura varieties forces SV |
| Pseudo-character target | Chenevier 2014 determinant $D_{\Delta_{10}}$ of dim 4 | Pseudo-representation (loses mod-$\ell^n$ Cayley--Hamilton for reducible $\rho_{\Delta_{10}}$) | Chenevier 2014 *Ann Inst Fourier* 64 |
| Hecke field of $\Delta_{10}$ | $\mathbb Q$ (since $\dim S_{26}(\mathrm{SL}_2(\mathbb Z)) = 1$); minimal coefficient ring $\mathbb Z$ | --- | Standard |
| Humbert divisor / AD correspondence | Humbert divisor $= $ Argyres--Douglas point at $E_{\tau_1} \times E_{\tau_2}$ (SW curve degenerates to pair of elliptic curves → $(A_1, A_{2N-1})$) | --- | GMN 2009 *Adv Theor Math Phys* 13 Ex 8.3 |
| K3 $\times$ E spectrum and $\kappa_\bullet$ indexing | $\{0,3,5,24\}$ from four distinct constructions: $\kappa_{\mathrm{cat}}(K3\times E)=0$ on the total space, $\kappa_{\mathrm{ch}}^K=3$, $\kappa_{\mathrm{BKM}}(\Delta_5)=5$, and the Mukai/Heisenberg--Fock $24$-trace. $\chi(\mathcal O_{K3})=2$ is only the K3-fibre witness, not a total-space value. | Naive "$\kappa_{\mathrm{BKM}} = \kappa_{\mathrm{ch}} + \chi(\mathcal O_{\mathrm{fibre}})$"; treating $2$ as $\kappa_{\mathrm{cat}}(K3\times E)$ or as a spectrum entry. | Always name the construction and the $\kappa_\bullet$ index |
| $\Delta_5/\Phi_{10}/\Phi_{12}$ denominator split | $\Delta_5$ is the primitive K3 paramodular denominator with $\kappa_{\mathrm{BKM}}=5$; $\Phi_{10}^{\mathrm{un}}=\Delta_5^2$ is the unnormalised Igusa square; OP/DT uses \(D_5=64^{-1}\Delta_5\) and \(\chi_{10}^{\mathrm{OP}}=D_5^2\); $\Phi_{12}$ is the Fake-Monster denominator with weight $12$. | Bare cross-volume values $5$ or $12$ without the input denominator; bare OP scalar \(-\Delta_5^{-2}\). | Gritsenko 1999 + OP normalisation cache |
| CoHA vs chiral-algebra type | $\mathrm{CoHA}_{K3 \times E}$ is $E_1$-associative (Hall product); chiralisation via $\Phi_3$-arrow gives $\mathbf H_{\Delta_5} = \Phi_3(\mathcal D_\hbar(\mathcal Y^{\mathrm{Hall}}(\mathrm{CoHA}_{K3 \times E})))$ | "CoHA is a chiral algebra" (type error) | Schiffmann--Vasserot + $\Phi_d$-functor framework |
| $\mathrm{CoHA}(\mathbb C^3)$ identification | $Y^+(\widehat{\mathfrak{gl}}_1)$ — the **positive half** of the affine Yangian. A $W_{1+\infty}$ vertex algebra appears only after the Hall--Drinfeld double, derived/chiral centre, and Fock/evaluation representation layers. | Directly identifying $\mathrm{CoHA}(\mathbb C^3)$ with $W_{1+\infty}$ or with the full Yangian $Y$. | Schiffmann--Vasserot |
| $W_{1+\infty}$ vs $W_\infty[c]$ | $W_{1+\infty} = W_\infty[c] \otimes \mathcal H$ (Heisenberg); different objects | Conflating the two | Pope--Romans--Shen 1990 |
| Two-$\hbar$ discipline | $\hbar^{\mathrm{Drinfeld}} = 2\pi i / \ell$ (root-of-unity) vs $\hbar^{\mathrm{BV}}$ (loop-counting); semantically distinct but agree numerically at $\hbar^2 = -1/8$ for $\ell = 8$ | Using bare $\hbar$ without subscript | AP151 bridge |
| Bar cohomology class-$\mathbf M$ at $E_3$-level | $6^g$ (cohomological dim at $g \in \{1, 2, 3\}$; chain level is infinite) | "Infinite at cohomological level" | Wave 12 Vol III |
| $\kappa_{\mathrm{cat}}(K3 \times E)$ | $0$ (total space, Künneth multiplicative $\chi(\mathcal O_{K3}) \chi(\mathcal O_E) = 2 \cdot 0$) | $2$ (fiber only) | AP-CY68/AP234 discipline |

**Every cache/AP entry citing an older value must appear in a "Wrong
Claim" / "retracted" / "earlier draft" context. Any standalone assertion
of an older value is itself a latent AP violation; report it via
Gate 0 on sight.**


## Compatible dual readings (non-contradictions)

Some statements appear in different waves with apparently conflicting
framings yet describe the **same mathematical object** through two
compatible frames. An agent reading only one framing may mistake the
pair for a contradiction and attempt to "fix" one side, destroying
content. Every entry below is a certified **non-contradiction**: both
framings are correct simultaneously and reconcile through the named
mechanism. Do not flatten either reading. Do not treat this as an AP
retraction list; treat it as a frame-compatibility registry.

| Statement A (frame 1) | Statement B (frame 2) | Compatibility mechanism | Primary reference |
|---|---|---|---|
| Conway $V^{s\natural}$ is the $\Psi^{\mathrm{metap}}$-image on the metaplectic $\overline{\mathcal A_2^{(2)}}$ branch (canonical preamble, row 41). | Conway $V^{s\natural}$ is the $\mathbb Z/2$-super-twin of $V^\natural$ via Duncan's commutative orbifolding diamond (W20.2 healing). | Same VOA-super-algebra read from the $\Psi$-functor / $\overline{\mathcal A_2}$-stratification side and from the commutative-orbifolding-diamond side; the metaplectic $\Psi$-image equals the super-twin of the Leech CFT via the Duncan 2007 Frame-change. | Duncan 2007 *Duke Math J* 139; Scheithauer 2008. |
| $\kappa_{\mathrm{BKM}}(\mathbf H_{\Delta_5}) = 12$ (Vol I three-faces identity, Fake-Monster $\Phi_{12}$ input). | $\kappa_{\mathrm{BKM}}(\mathbf H_{\Delta_5}) = 5$ (Vol III Borcherds-weight count on paramodular $\Delta_5$). | Denominator split: the two values index different input denominators ($\Phi_{12}$ vs $\Phi_{10}^{\mathrm{un}}=\Delta_5^2$) under the same $c_N(0)/2$ rule. Each site must name the input denominator, and OP/DT sites must name \(D_5=64^{-1}\Delta_5\); the two BKM weights are $(c_{\Phi_{12}}(0)/2, c_{\Delta_5}(0)/2) = (12, 5)$. | Borcherds 1992; Gritsenko 1999 Thm 6.1; canonical preamble row 57. |
| $e_4 = W_4 - (107/11)\Lambda_Z$ in the $W$-basis (Vol III abstract, $c = -214$). | $e_4 = {:}T\partial^2 T{:} - (3/2){:}(\partial T)^2{:} + (321/10)\partial^4 T + \hbar\,\mathrm{qt}(J^{(4)})$ in the Virasoro-composite basis (canonical preamble row 51). | Same cocycle class, two bases. Pope--Romans--Shen projection $W_k = e_k + (-c/22)\Lambda_Z$ converts between $W$-basis and Virasoro-composite basis; the scalar $-107/11 = c/22$ at $c = -214$ is the frame-change constant. | Pope--Romans--Shen 1990; Wang 1998 Prop 4.2. |
| Archetype $\mathbf B$ row: $\kappa + \kappa^! = 8$ (Vol I five-archetype table, $r(z)$-families). | K3 Mukai pairing $K = 2c_+ = 8$ (Vol III Mukai form on $\Lambda^{3,19}$). | Three faces of the same canonical $K = 8$: (i) derived-centre complementarity sum $\kappa + \kappa^!$; (ii) Mukai pairing $K = 2c_+$ with $c_+ = 4$; (iii) Lusztig reflection length $\ell_{\mathrm{Lusztig}} = 8$ at $\zeta_8$. All three equal $8$ because $K3 \times \mathrm{elliptic}$ pairs with $\mathbf H_{\Delta_5}$ through a single self-dual lattice indexing. | Borcherds--Mukai pairing; Lusztig 1989 at $\zeta_8$; Gritsenko--Nikulin 1998. |
| Padovan counts $d_n = d_{n-2} + d_{n-3}$ give MZV basis dimensions in weight $n$ (Brown 2012 conjecture). | Humbert-Heegner divisor filter selects arithmetically admissible $n \equiv 3, 5 \pmod 8$ on $\phi^{(n)}$ (Vol I / Vol III pentagon admissibility row). | Orthogonal filters. Padovan counts basis vectors; Humbert-Heegner selects admissible embedding points. The two do not compete: Padovan $d_n$ governs the motivic dimension; Humbert-Heegner governs the arithmetic non-vanishing locus. Compatible because $\phi^{(n)}$ lives in a weight-$n$ MZV-period space filtered further by admissibility. | Brown 2012; Eichler--Zagier 1985 Thm 9.1 (canonical row 44). |
| $(\infty, 1)$-obstruction tower (W21.2) gives the Maurer--Cartan tower of universal obstructions in an $(\infty, 1)$-stable category. | $(\infty, 2)$-adjunction AR structure (W22.1) exhibits bar--cobar as an adjunction in an $(\infty, 2)$-category with an Auslander--Reiten sequence. | The $(\infty, 2)$ structure is the $2$-categorical lift of the $(\infty, 1)$ obstruction tower: AR sequence at object level descends to the MC obstruction tower at morphism level. Both proved, both load-bearing. | Lurie HA 5.5; Riehl--Verity 2020. |
| Six-path $\chi_3$ (W22.6): six independent verification paths for the third Chern character. | Seven-path $\chi_3$ (W25.4): seven paths including a derived-deformation check. | Six independent paths plus one additional derived-deformation path, consistent with prior comparison. The seventh path refines the six without contradicting them. Every path reports the same $\chi_3$ value. | Beilinson multi-path verification discipline; W25.4 cache entry. |
| Four sibling functors $\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$ parametrise the Baily--Borel--Freitag stratification of $\overline{\mathcal A_2}$ (canonical preamble row 46). | Five $\Psi$-image rows including the Conway $V^{s\natural}$ row on the metaplectic branch (W20.2 table). | The four siblings parametrise the four strata; each sibling has one image row, and the $\Psi^{\mathrm{metap}}$ sibling hosts two rows (Conway $V^{s\natural}$ and a second metaplectic image). Row count $5 = 4$ siblings with one extra Conway row on the metaplectic branch. | Baily--Borel--Freitag; Duncan 2007. |
| Absolute Kuznetsov HPD blocked on $K3 \times E$ (obstruction from the non-trivial Brauer class on $K3 \times E$). | Relative HPD over $E$ works (the relative $E$-base kills the obstruction by descent). | Obstruction at the absolute level; healing at the relative level. Both statements are proved theorems; relative HPD is the correct home for the $K3 \times E$ Kuznetsov lift. | Kuznetsov 2007; Perry 2019 relative HPD. |
| Pseudo-character framework (Taylor 1991) is sufficient on reduced rings for Galois deformations of $\rho_{\Delta_{10}}$. | Chenevier determinant framework (2014) is necessary on non-reduced rings (where mod-$\ell^n$ Cayley--Hamilton fails for reducible $\rho$). | Compatible hierarchy, not competing frameworks. Taylor's pseudo-character suffices on reduced; Chenevier's determinant is the correct lift to non-reduced. Each dominates in its natural domain; neither refutes the other. Canonical preamble row 53 pins the non-reduced case to Chenevier; reduced-ring calculations may still cite Taylor. | Taylor 1991 *Invent Math* 116; Chenevier 2014 *Ann Inst Fourier* 64. |

**Operating rule.** If you encounter two inscriptions that seem to
contradict, before editing either side **check this table**. If the
pair matches an entry, both are correct. If not, escalate through
Gate 0: the confusion is a candidate for a new compatibility mechanism,
not a candidate for retraction.


---

Scope: Vol II ($A_\infty$ Chiral Algebras and 3D HT QFT). Patterns discovered during Waves 14-24 attack-heal programme.

*Relocated to `notes/` per Manuscript Metadata Hygiene: anti-pattern
catalogues are working-notes scaffolding and do not belong in the
typeset manuscript.*

## Severity levels

| Level    | Meaning                                                    | Action                                 |
| -------- | ---------------------------------------------------------- | -------------------------------------- |
| Critical | Theorem status wrong (conjecture $\to$ theorem) or scope inflation at climax | Immediate fix; audit all instances     |
| High     | Numerical or structural error propagates                   | Fix before next build                  |
| Medium   | Convention clash or ambiguity                              | Fix in current session                 |
| Low      | Cosmetic or cross-reference staleness                      | Fix in batch                           |

## Section 1: Bar-cobar structure and the Swiss-cheese octachotomy

### AP-V2-1: Swiss-cheese octachotomy as colouring
**Severity**: Critical
**Trigger**: Assertion that the eight bar-cobar ambients
(fibrewise / single-monodromy / bi-unipotent / tri-unipotent /
weight-completed / $A_\infty$-corrected / $(\infty,1)$-$\mathrm{Perf}$ /
Kontsevich-formal) are interchangeable or that the octachotomy is a
mere indexing convention.
**Ghost**: Each of the 8 ambients corresponds to a specific
$\mathsf{SC}^{\mathrm{ch,top}}$ colouring stratum on the Swiss-cheese
operadic cover of the boundary-bulk pair; the octachotomy is a
stratification of the coloured operadic cover, not a choice of
formalism.
**Error**: Category error: moving a bar-cobar statement across
ambients (e.g., from bi-unipotent to Kontsevich-formal) requires
crossing a coloured stratum and may incur $\mathrm{grt}_1$-torsor
obstructions that do not appear in any single stratum.
**Correction**: Construction: label each bar-cobar statement with its
ambient $\mathcal A^{(i)}$, $i \in \{1,\dots,8\}$; $\mathsf{SC}^{\mathrm{ch,top}}$
carries a 2-colouring $\mathfrak c: \{e, o\} \to \{1,\dots,8\}$ whose
strata refine Dunn's additive dichotomy. Remark: Pattern 277 is the
local-to-global restatement; each colouring stratum carries its own
adjunction $\Omega^{\mathrm{ch}}_{(i)} \dashv B^{\mathrm{ch}}_{(i)}$
and the crossing maps are non-trivial natural transformations
$\Omega^{\mathrm{ch}}_{(i)} \Rightarrow \Omega^{\mathrm{ch}}_{(j)}$.
**Primary**: Francis--Gaitsgory 2012 *Selecta Math* 18 §2.4 (Ran-space
factorisation stratification); Voronov 1998 *Math Ann* 312
(Swiss-cheese operads); Tamarkin 2003 *Lett Math Phys* 66 (GRT torsor
structure).
**Inscription**: `chapters/theory/sc_chtop_octachotomy.tex`, Thm 2.1
(\texttt{thm:sc-octachotomy-colouring}).
**Cross-ref**: Vol I Pattern 269 (adjunction-strictness conflation);
Vol III AP-CY33 (chain vs rational).

### AP-V2-2: Universal $k$-tower genus formula
**Severity**: High
**Trigger**: Omitting the genus-dependent upper bound
$k_{\max}(g) = g(g+1)/2$ on the Swiss-cheese coloured coherence tower,
or asserting a genus-independent bound.
**Ghost**: The coloured coherence tower on $\overline{\mathcal M_g}$
truncates at height $k_{\max}(g) = g(g+1)/2$, controlled by the rank
of the intersection form on $H^1(\Sigma_g, \mathbb Z)$ via
Mumford--Morita--Miller classes.
**Error**: Claiming chain-level $E_n$-coherences up to arbitrary $n$
at all genera ignores the symplectic rank bottleneck; the higher
coherences simply do not exist above $k_{\max}(g)$.
**Correction**: Construction: for each genus $g$, the Swiss-cheese
coloured coherence tower stabilises at $k_{\max}(g) = g(g+1)/2$;
at $g=0$ trivial, $g=1$ height $1$, $g=2$ height $3$, $g=3$ height
$6$. Remark: matches triangular numbers $T_g$ and equals
$\dim_{\mathbb Q} H^1(\Sigma_g, \mathbb Q)^{\wedge 2} / (\text{image of } \cup)$.
**Primary**: Morita 1989 *Invent Math* 90 (Mumford--Morita--Miller);
Costello--Gwilliam 2021 Vol 2 §7.3 (genus-scaling of BV anomaly);
Hain 1997 *Invent Math* 128 (mapping-class group cohomology).
**Inscription**: `chapters/theory/higher_genus_modular_koszul.tex`,
Prop 4.2 (\texttt{prop:coloured-tower-genus-bound}).
**Cross-ref**: Vol I AP27 (propagator weight); Vol III AP-CY12
(shadow class from full tower).

### AP-V2-3: Codim-4 admissible Heegner quadruple emptiness
**Severity**: High
**Trigger**: Asserting emptiness of admissible Heegner quadruples in
codimension 4 on $\overline{\mathcal A_2}$ by dimensional argument
alone ($\dim = 3$, so codim-4 strata are empty).
**Ghost**: Dimension count gives the correct answer (codim-4 strata
are indeed empty on $\overline{\mathcal A_2}$), but the genuine proof
is arithmetic: Bruinier--Humbert rank obstruction via regularised
theta lift rules out admissible quadruples by $L$-function
non-vanishing at $s=0$.
**Error**: Dimensional argument is insufficient because the stratum
closure locus can exceed ambient dimension under Baily--Borel
compactification; only the arithmetic obstruction rules out
admissible quadruples unconditionally.
**Correction**: Construction: Bruinier 2002 *Lecture Notes 1780*
Thm 5.11 regularised theta-lift non-vanishing + Humbert divisor
arithmetic gives admissibility rank bound; a codim-4 admissible
quadruple would produce an $L$-value of the wrong sign at $s=0$.
Remark: dimensional argument is a heuristic shadow of the arithmetic
content.
**Primary**: Bruinier 2002 *Lecture Notes Math* 1780 §5.3;
Humbert 1899 *Liouville J* 5 Ser 5 (Humbert surfaces);
Borcherds--Katzarkov--Pantev--Shepherd-Barron 1998 *Compositio* 122
(Picard-modular singular weight).
**Inscription**: `chapters/theory/codim4_heegner_emptiness.tex`,
Thm 3.2 (\texttt{thm:codim4-bruinier-humbert}).
**Cross-ref**: Vol III AP-CY42 (pentagon admissibility mod 4).

## Section 2: $(\infty,1)$ and $(\infty,2)$-categorical lifts

### AP-V2-4: Compatibility data canonicality
**Severity**: Critical
**Trigger**: Writing the bar-cobar compatibility homotopy as "chosen"
data or as a choice made by the author, rather than as the image of a
canonical torsor.
**Ghost**: The space of compatibility data is the canonical torsor
$\mathrm{Data}(A) \simeq B\mathrm{Aut}^h_{L_\infty}(\mathfrak g^{\infty,1}_A)$
(the classifying space of $L_\infty$-homotopy automorphisms of the
underlying $(\infty,1)$-algebra), not an arbitrary choice.
**Error**: Headline-vs-full-structure confusion: stating the
compatibility homotopy as "a choice" truncates the higher structure to
the 0-truncation and loses the $L_\infty$-automorphism $(\infty,1)$-group
action that governs deformation theory.
**Correction**: Construction: $\mathrm{Data}(A)$ is an
$(\infty,0)$-groupoid, canonically equivalent to
$B\mathrm{Aut}^h_{L_\infty}(\mathfrak g^{\infty,1}_A)$ where
$\mathfrak g^{\infty,1}_A = \mathrm{Der}^\bullet(A)[1]$ is the tangent
$L_\infty$-algebra. Any two compatibility data are canonically
equivalent through a unique-up-to-contractible $L_\infty$-automorphism.
Remark: the full $(\infty,1)$-structure is essential for the $\mathrm{grt}_1$
action on the bar-cobar GRT-torsor.
**Primary**: Lurie 2017 *Higher Algebra* §5.5.4 (canonical equivalence
in presentable $\infty$-categories); Pridham 2010 *Adv Math* 224
(deformation theory via $L_\infty$); Hinich 1997 *Comm Algebra* 25
(homotopy automorphisms).
**Inscription**: `chapters/theory/compatibility_torsor.tex`, Thm 2.4
(\texttt{thm:compatibility-torsor-canonical}).
**Cross-ref**: Vol I AP8 (canonical-vs-chosen);
Vol III AP-CY30 (factored $\neq$ solved).

### AP-V2-5: $(\infty,2)$-adjunction headline vs 2-morphism content
**Severity**: Critical
**Trigger**: Presenting the bar-cobar $(\infty,2)$-adjunction only via
its unit/counit (1-morphism data) without the Auslander--Reiten
triangular structure on 2-morphisms.
**Ghost**: The $(\infty,2)$-adjunction $\Omega^{\mathrm{ch}} \dashv B^{\mathrm{ch}}$
lives in an $(\infty,2)$-category whose 2-morphisms carry
Auslander--Reiten triangles; the 1-truncation captures only the
headline, losing the triangular Hom-functor structure that encodes
the derived centre.
**Error**: 1-truncation throws away the AR-triangle content; the unit
and counit are 1-morphisms, but the content lies in the 2-morphism
$\mathrm{Hom}$-spaces and their triangular composition.
**Correction**: Construction: state the adjunction in
$\mathrm{Cat}^{(\infty,2)}_{\mathrm{chiral}}$; 2-morphisms form
Auslander--Reiten triangles
$\tau X \to \mathrm{Hom}(X, Y) \to X \oplus Y$
where $\tau$ is the AR translation induced by the chiral Serre
functor. Remark: the 1-truncation is $\mathrm{Cat}^{(\infty,1)}$
and suffices for the adjunction headline; the 2-morphism content is
essential for the derived-centre complementarity theorem (Theorem C).
**Primary**: Gaitsgory--Rozenblyum 2017 *A Study in Derived Algebraic
Geometry I* §I.1.1 ($(\infty,2)$-categories);
Auslander--Reiten 1975 *Comm Algebra* 3 (AR triangles);
Lurie 2009 *$(\infty,2)$-Categories and the Goodwillie Calculus I*.
**Inscription**: `chapters/theory/infty2_adjunction_ar.tex`, Thm 3.1
(\texttt{thm:adjunction-infty2-ar-triangles}).
**Cross-ref**: Vol I AP24 (three-functor conflation); Vol III AP-CY3
($E_2 \neq$ commutative).

## Section 3: Chiral Kontsevich formality

### AP-V2-6: Kontsevich smooth-manifold transport
**Severity**: Critical
**Trigger**: Direct transport of Kontsevich 1997 smooth-manifold
formality across $\overline{\mathcal A_2}$ without restricting to the
Koszul admissible locus $U^{\mathrm{adm}}$.
**Ghost**: Kontsevich formality requires a smooth ambient manifold;
$\overline{\mathcal A_2}$ is stratified with Humbert divisor
singularities, so the direct transport fails. The valid transport
goes through $U^{\mathrm{adm}} \subset \overline{\mathcal A_2}$ via
Francis--Gaitsgory $\infty$-operadic and Tamarkin operadic formality.
**Error**: Smoothness hypothesis violation: $\overline{\mathcal A_2}$
carries codimension-1 Humbert divisors and the symmetric-square
stratum along the hyperelliptic locus; Kontsevich's construction does
not extend to stratified ambients without per-stratum correction.
**Correction**: Construction: restrict chiral formality to
$U^{\mathrm{adm}} = \overline{\mathcal A_2} \setminus
\bigcup_n H_n^{\mathrm{non-adm}}$; on $U^{\mathrm{adm}}$, combine
Francis--Gaitsgory $\infty$-operadic formality with Tamarkin operadic
formality to obtain a quasi-isomorphism
$B^{\mathrm{ch}}(A)|_{U^{\mathrm{adm}}} \simeq \mathrm{Tam}(A)|_{U^{\mathrm{adm}}}$.
Remark: on Humbert divisors, the formality fails with explicit
obstruction class in $H^2_{\mathrm{Kont}}|_{H_n}$.
**Primary**: Kontsevich 1997 *Lett Math Phys* 66 §7 (smooth manifolds
only); Francis--Gaitsgory 2012 *Selecta Math* 18;
Tamarkin 2003 *Lett Math Phys* 66; Willwacher 2015 *Invent Math* 200
(stratified-ambient obstructions).
**Inscription**: `chapters/theory/chiral_kontsevich_formality.tex`,
Thm 4.3 (\texttt{thm:chiral-formality-admissible-locus}).
**Cross-ref**: Vol I AP34 (stratified-ambient scope).

### AP-V2-7: Per-wall vs global Kontsevich obstruction sheaf
**Severity**: High
**Trigger**: Conflating the per-Humbert-divisor obstruction
$\mathrm{ob}^{\mathrm{Kont}}_n$ (local on $H_n$) with the global
pushforward sheaf $\mathrm{ob}^{\mathrm{Kont}}$ on
$\overline{\mathcal A_2}$.
**Ghost**: $\mathrm{ob}^{\mathrm{Kont}}_n$ is a sheaf on the Humbert
divisor $H_n$, cohomology class in $H^2(H_n)$; the global
$\mathrm{ob}^{\mathrm{Kont}} = \bigoplus_n (\iota_n)_* \mathrm{ob}^{\mathrm{Kont}}_n$
on $\overline{\mathcal A_2}$ is the pushforward, cohomology class in
$H^2_{\{H_n\}}(\overline{\mathcal A_2})$.
**Error**: Scope error: writing "the Kontsevich obstruction vanishes"
as a global statement when the vanishing holds only per-wall, or vice
versa, confuses local and global contributions.
**Correction**: Construction: $\mathrm{ob}^{\mathrm{Kont}}_n
\in H^2(H_n, \mathcal N^\vee_{H_n/\overline{\mathcal A_2}})$ (local);
$\mathrm{ob}^{\mathrm{Kont}} = \sum_n (\iota_n)_* \mathrm{ob}^{\mathrm{Kont}}_n
\in H^2_{\{H_n\}}(\overline{\mathcal A_2}, \mathcal N^\vee)$ (global).
Remark: per-wall vanishing does not imply global vanishing; a
residual class in $H^3(\overline{\mathcal A_2})$ can appear via
$E_2^{\mathrm{Leray}}$ spectral sequence.
**Primary**: Hartshorne 1977 *Algebraic Geometry* III.8 (local
cohomology); Willwacher 2015 *Invent Math* 200 §4 (per-wall
obstructions); Bruinier--Funke 2004 *J Reine Angew Math* 570
(regularised theta-lift on Humbert divisors).
**Inscription**: `chapters/theory/chiral_kontsevich_formality.tex`,
Prop 4.7 (\texttt{prop:kont-obstruction-local-vs-global}).
**Cross-ref**: Vol III AP-CY13 (stale Part references, scope error
as cousin).

## Section 4: BD factorisation compatibility

### AP-V2-8: Bi-unipotent pro-ambient vs BD $\mathcal D$-module
**Severity**: Critical
**Trigger**: Asserting the bi-unipotent bar-cobar ambient is merely a
bounded-derived pro-object (e.g., $\mathrm{Pro}(D^b_{\mathrm{coh}})$)
without the genuine BD factorisation-$\mathcal D$-module structure.
**Ghost**: The bi-unipotent ambient carries a genuine
Beilinson--Drinfeld factorisation-$\mathcal D$-module structure on
$\mathrm{Ran}\,X$; the pro-object description is a shadow that loses
the factorisation compatibility and the 1-loop counterterm content.
**Error**: Type error: $\mathrm{Pro}(D^b)$ is a category, not a
$\mathcal D$-module; the factorisation structure
$\Delta_* \mathcal A \simeq \mathcal A \boxtimes \mathcal A$ along
the diagonal is essential and is not captured by the pro-object
description.
**Correction**: Construction: the bi-unipotent ambient is a BD
factorisation $\mathcal D$-module $\mathcal A^{\mathrm{BD}}_{\mathrm{bi-uni}}$
on $\mathrm{Ran}\,X$ with factorisation isomorphism
$\mathcal A^{\mathrm{BD}}_{\mathrm{bi-uni}}|_{\mathrm{Ran}\,X_1 \sqcup \mathrm{Ran}\,X_2}
\simeq \mathcal A^{\mathrm{BD}}_{\mathrm{bi-uni}}|_{\mathrm{Ran}\,X_1}
\boxtimes \mathcal A^{\mathrm{BD}}_{\mathrm{bi-uni}}|_{\mathrm{Ran}\,X_2}$
along disjointness opens. Remark: Costello--Gwilliam 2021 Vol 2
Thm 7.1.1 gives the 1-loop counterterm $= 8$ which verifies the BD
structure (not the pro-object shadow); the specific value $8$ is the
Euler characteristic contribution of the bi-unipotent stratum.
**Primary**: Beilinson--Drinfeld 2004 *Chiral Algebras* Ch 3;
Francis--Gaitsgory 2012 *Selecta Math* 18; Costello--Gwilliam 2021
Vol 2 *Factorization Algebras in Quantum Field Theory* Thm 7.1.1;
Gwilliam 2012 Northwestern thesis.
**Inscription**: `chapters/theory/bi_unipotent_bd_structure.tex`,
Thm 5.2 (\texttt{thm:bi-uni-bd-factorisation}).
**Cross-ref**: Vol I AP165 (B(A) is E_1 coalgebra, not SC);
Vol III AP-CY30 (factored $\neq$ solved).

## Section 5: Theorem H Vol II $\leftrightarrow$ Vol III hinge

### AP-V2-9: Theorem H concentration $d$-independence
**Severity**: Critical
**Trigger**: Asserting chiral Hochschild concentration
$\mathrm{ChirHoch}^\bullet(A) \in \{0,1,2\}$ holds at all dimensions
$d$ of the CY input, without scope qualifier.
**Ghost**: Theorem H concentration $\{0,1,2\}$ is the **ambient-curve
scope** statement (Vol II, $X$ an algebraic curve); on the
$\Phi_d$-enlarged ambient with a CY-$d$ transverse direction, the
concentration enlarges to $\{0,1,2,d\}$ (Vol III, scope for
$d \geq 2$).
**Error**: Scope inflation: the Vol II statement is tight on curves
($d=1$); transporting to CY-$d$ without enlarging to $\{0,1,2,d\}$
misses the transverse direction contribution and gives wrong
Hochschild spectral-sequence page.
**Correction**: Construction: Vol II Theorem H (on curves):
$\mathrm{ChirHoch}^\bullet(A) \in \{0,1,2\}$ for $A$ a chiral algebra
on curve $X$. Vol III enlargement via $\Phi_d$: for $d \geq 2$,
$\mathrm{ChirHoch}^\bullet(\Phi_d(\mathcal C)) \in \{0,1,2,d\}$ where
the extra $d$-entry is the transverse-direction contribution.
Remark: the hinge is the Koszul-locus equivalence
$\mathrm{ChirHoch}^\bullet_{\mathrm{Vol II}}|_{U^{\mathrm{Koszul}}}
\simeq \mathrm{ChirHoch}^\bullet_{\mathrm{Vol III}}|_{U^{\mathrm{Koszul}}}$
which holds on the Koszul locus only.
**Primary**: Francis 2013 *Compositio* 149 §3 (Hochschild
concentration for $E_n$-algebras); Gaitsgory--Lurie 2019 *Weil
Conjecture for Function Fields* Ch 4 (chiral Hochschild); Costello
2011 *Renormalization and Effective Field Theory* §5.
**Inscription**: `chapters/theory/theorem_h_hinge_vol2_vol3.tex`,
Thm 6.1 (\texttt{thm:theorem-h-scope-hinge}).
**Cross-ref**: Vol I Theorem H statement; Vol III AP-CY17 (MF CY
dimension, transverse-direction analogue).

## Section 6: Cyclic chiral homology (Feigin--Tsygan 5th verification path)

### AP-V2-10: Classical Feigin--Tsygan verbatim transport
**Severity**: Critical
**Trigger**: Applying classical Feigin--Tsygan cyclic homology
directly to chiral bialgebras via verbatim transport of the classical
constructions.
**Ghost**: Classical Feigin--Tsygan is built for associative algebras;
chiral bialgebras are non-associative (the chiral product is defined
only on disjoint support), so verbatim transport fails. Five chiral
steps are required for the genuine chiral-cyclic homology.
**Error**: Associativity mismatch: Connes $B$ operator on a chiral
bialgebra requires an $S^1$-action on the bar complex that respects
the factorisation structure; the classical $B$ does not respect
factorisation.
**Correction**: Construction: five-step chiral Feigin--Tsygan: (i)
Francis--Gaitsgory $E_1$-operadic upgrade of the chiral bar complex;
(ii) $S^1$-action providing Connes $B^{\mathrm{ch}}$ on $B^{\mathrm{ch}}(A)$;
(iii) SBI invertibility on $U^{\mathrm{adm}}$ via Kontsevich formality
spectral-sequence collapse at $E_2$-page; (iv) CoHA-to-cyclic matching
through admissible-graph integral
$\int_{\mathrm{Graphs}^{\mathrm{adm}}_n} \omega_{\mathrm{Kont}}$;
(v) explicit Serre-dual numerical evaluation giving
$\chi_3^{\mathrm{ch}} = c(5c+22)/10$-normalised cyclic
characteristic. Remark: the five steps are independently checkable;
each provides a verification path.
**Primary**: Feigin--Tsygan 1987 *K-Theory* 1 (classical cyclic);
Francis--Gaitsgory 2012 *Selecta Math* 18 §5 ($E_1$-operadic upgrade);
Kassel 1987 *J Pure Appl Algebra* 48 (Connes $B$); Kontsevich--Soibelman
2009 arXiv:0606241 (CoHA cyclic structure).
**Inscription**: `chapters/theory/cyclic_chiral_homology.tex`, Thm 7.2
(\texttt{thm:chiral-feigin-tsygan-five-step}).
**Cross-ref**: Vol I AP35--43 (epistemic anti-patterns);
Vol III AP-CY23 ($E_1$ not $E_\infty$ bialgebra).

### AP-V2-11: Three cyclic theories conflation
**Severity**: High
**Trigger**: Conflating any two of: ordinary algebraic cyclic
$\mathrm{HC}^\bullet_{\mathrm{alg}}$, topological
$\mathrm{TC}/\mathrm{TP}$, chiral cyclic $\mathrm{HC}^\bullet_{\mathrm{ch}}$.
**Ghost**: The three theories are distinct and their distinctness is
load-bearing: conflating any two collapses the fifth $\chi_3$
verification path (the one that uses independence of the three).
**Error**: Structural collapse: ordinary $\mathrm{HC}^\bullet_{\mathrm{alg}}$
is $\pi$-exact for DG-algebras; $\mathrm{TC}/\mathrm{TP}$ is defined
on spectra via Bokstedt; chiral $\mathrm{HC}^\bullet_{\mathrm{ch}}$ is
defined on $\mathrm{Ran}\,X$ via factorisation. They agree on
commutative rational input but differ for non-commutative or
non-rational input.
**Correction**: Construction: state for each use which cyclic theory
is invoked; triangulation via all three (direct algebraic / topological
via Bokstedt / chiral via factorisation) gives a fifth $\chi_3$ path
independent of Serre duality. Remark: the three theories fit into a
diagram of natural transformations
$\mathrm{HC}^\bullet_{\mathrm{alg}} \to \mathrm{TC} \to
\mathrm{HC}^\bullet_{\mathrm{ch}}$; each arrow is an equivalence on
the respective semiprime locus but not elsewhere.
**Primary**: Loday 1998 *Cyclic Homology* §4 (algebraic);
Bokstedt--Hsiang--Madsen 1993 *Invent Math* 111 (topological);
Francis 2013 *Compositio* 149 §5 (chiral).
**Inscription**: `chapters/theory/three_cyclic_theories.tex`, Prop 7.5
(\texttt{prop:three-cyclic-distinct}).
**Cross-ref**: Vol I AP152 ("ordered" ambiguity, cousin-collapse).

## Section 7: NC Hodge Mukai--Hochschild

### AP-V2-12: Kaledin degeneration scope
**Severity**: Critical
**Trigger**: Transporting Kaledin 2008 noncommutative Hodge
degeneration from smooth-proper DG-scheme scope to
factorisation-$\infty$-categorical scope on $\overline{\mathcal A_2}$
without checking the Humbert-stratum obstructions.
**Ghost**: Kaledin degeneration holds in smooth-proper DG-scope; on
$\overline{\mathcal A_2}$ the Humbert strata carry non-zero
Heegner--Bruinier obstruction classes that feed into higher-page
$d_r$ of the Hodge-to-de-Rham spectral sequence, so global
degeneration fails.
**Error**: Scope inflation: Kaledin's Cor 5.5 requires the ambient to
be smooth-proper; $\overline{\mathcal A_2}$ is not, and the
Hodge-to-de-Rham $E_2$-page carries non-zero $d_2$ on Humbert divisors.
**Correction**: Construction: Hodge-to-de-Rham degenerates on the
complement $\overline{\mathcal A_2} \setminus \bigcup H_n$;
on Humbert divisors $H_n$, the $d_2$ differential is the
Heegner--Bruinier class
$d_2([\alpha]) = [\alpha] \smile \mathrm{ob}^{\mathrm{HB}}_n$.
Remark: higher-page differentials $d_r$ for $r \geq 3$ can be
non-zero even on the complement, via arithmetic of iterated Humbert
intersections.
**Primary**: Kaledin 2008 *J Pure Appl Algebra* 212 Thm 1.1 (smooth
proper DG); Kaledin 2010 *Ann Inst Fourier* 60 Cor 5.5 (Serre
duality); Bruinier 2002 *Lecture Notes Math* 1780 §5 (Heegner
obstruction).
**Inscription**: `chapters/theory/nc_hodge_mukai_hochschild.tex`,
Thm 8.3 (\texttt{thm:kaledin-scope-humbert}).
**Cross-ref**: Vol III AP-CY1 (CY dim $\neq$ cpx dim).

### AP-V2-13: Mukai--Hochschild pairing triangulation
**Severity**: High
**Trigger**: Computing the Mukai--Hochschild pairing
$\langle[\chi_3], [\chi_3^\vee]\rangle = 2\,\mathrm{Vol}(E)(2\pi i)^3$
by a single-path argument, without independent triangulation.
**Ghost**: The pairing value $2\,\mathrm{Vol}(E)(2\pi i)^3$ is the
correct numerical answer, but its derivation requires triangulation
via three independent paths: (i) Kaledin 2010 Cor 5.5 Serre duality;
(ii) BV Serre-dual cup product; (iii) Siegel--Eisenstein period
cross-check.
**Error**: Single-path insufficiency: the factor
$2\,\mathrm{Vol}(E)$ arises from $\chi_{\mathrm{top}}(E) \cdot \mathrm{Vol}$,
the $(2\pi i)^3$ from degree-$3$ NC Hodge structure; a single path
can produce the correct numeric by chance but fails multi-path
triangulation.
**Correction**: Construction: triangulate via (i) Kaledin 2010
Cor 5.5 applied to $\mathrm{Perf}(E) \times \mathrm{Perf}(K3)$
Serre-dual trace; (ii) BV Serre-dual cup product in
$\mathrm{HH}^\bullet_{\mathrm{BV}}(\mathcal A)$ using
Costello--Gwilliam 2021 Vol 2 Thm 5.2.3; (iii) Siegel--Eisenstein
period cross-check through Pitale--Saha--Schmidt 2014 $L$-value at
$s=0$. All three paths converge on $2\,\mathrm{Vol}(E)(2\pi i)^3$.
Remark: absence of any one path indicates a computation error.
**Primary**: Kaledin 2010 *Ann Inst Fourier* 60 Cor 5.5;
Costello--Gwilliam 2021 Vol 2 *Factorization Algebras* Thm 5.2.3;
Pitale--Saha--Schmidt 2014 *Memoirs AMS* 232 §7;
Shimura 1976 *Comm Pure Appl Math* 29 (Eisenstein periods).
**Inscription**: `chapters/theory/mukai_hochschild_pairing.tex`,
Thm 8.5 (\texttt{thm:mh-pairing-triangulation}).
**Cross-ref**: Vol III AP-CY40 (master $L$-value triangulation).

## Section 8: CY-4 HK-restricted $\Phi_4$

### AP-V2-14: Generic CY-4 $\Phi_4$ existence
**Severity**: Critical
**Trigger**: Asserting that the CY-to-chiral functor $\Phi_4$ exists
on generic CY-4 input without scope restriction.
**Ghost**: $\Phi_4$ existence is BLOCKED on generic CY-4 by
Kapustin--Rozansky--Saulina 3d/4d dichotomy: the Rozansky--Witten
twist producing 3d HT content requires hyperkähler structure on the
CY-4 target (not merely holomorphic symplectic). Scope-restricted
existence on hyperkähler CY-4 via twistor-$S^1$ reduction holds.
**Error**: Scope inflation: Rozansky--Witten produces a 3d B-twist on
hyperkähler targets; generic CY-4 (holomorphic symplectic but not
hyperkähler) lacks the $S^1$-twistor family needed for the reduction.
**Correction**: Construction: $\Phi_4^{\mathrm{HK}}$ defined on
hyperkähler CY-4 via Rozansky--Witten 1997 + twistor-$S^1$ reduction
+ Kapustin--Rozansky--Saulina 2009 3d/4d dichotomy. Generic CY-4
$\Phi_4$ conjectural; hyperkähler scope unconditional.
Remark: hyperkähler CY-4 examples: $K3 \times K3$, $K3^{[2]}$,
generalised Kummer $K_2(A)$.
**Primary**: Rozansky--Witten 1997 *Selecta Math* 3 (3d B-twist);
Kapustin--Rozansky--Saulina 2009 *Nucl Phys B* 823 (3d/4d dichotomy);
Huybrechts 1999 *Invent Math* 135 (hyperkähler topology).
**Inscription**: `chapters/theory/phi4_hk_restricted.tex`, Thm 9.1
(\texttt{thm:phi4-hk-scope}).
**Cross-ref**: Vol III AP-CY6 ($A_X$ at $d=3$, parallel scope).

### AP-V2-15: K3$\times$K3 pairing triangulation
**Severity**: High
**Trigger**: Computing the CY-4 HK-restricted pairing
$\langle \Phi_4(K3_1), \Phi_4(K3_2)\rangle = 4\,\mathrm{Vol}(K3_1)\mathrm{Vol}(K3_2)(2\pi i)^4$
by single-path argument.
**Ghost**: The pairing value is correct but requires triangulation via
three independent paths: (i) Cao--Kool--Monavari DT-4 MacMahon vertex;
(ii) Nekrasov 8d prepotential equivariant localisation;
(iii) Dijkgraaf--Moore--Verlinde--Verlinde 1997 second-quantised
elliptic genus.
**Error**: Single-path insufficiency analogous to AP-V2-13; the
coefficient $4 = 2 \cdot 2$ combines $\chi_{\mathrm{top}}(K3) = 24$
with the hyperkähler double-cover factor.
**Correction**: Construction: triangulate via (i) Cao--Kool--Monavari
2021 *J Differ Geom* 117 DT-4 MacMahon vertex on $K3_1 \times K3_2$;
(ii) Nekrasov 2003 8d prepotential equivariant localisation (genus
$0$ contribution only, by HK-restriction);
(iii) DMVV 1997 *Commun Math Phys* 185 second-quantised elliptic
genus at $K3 \times K3$. All three paths converge on
$4\,\mathrm{Vol}(K3_1)\mathrm{Vol}(K3_2)(2\pi i)^4$. Remark: the
$(2\pi i)^4$ encodes degree-$4$ NC Hodge structure; the $4$
encodes chi-top halving.
**Primary**: Cao--Kool--Monavari 2021 *J Differ Geom* 117 DT-4
MacMahon; Nekrasov 2003 *Adv Theor Math Phys* 7 8d prepotential;
Dijkgraaf--Moore--Verlinde--Verlinde 1997 *Commun Math Phys* 185.
**Inscription**: `chapters/theory/k3k3_pairing_triangulation.tex`,
Thm 9.4 (\texttt{thm:k3k3-pairing-three-path}).
**Cross-ref**: Vol III AP-CY15 (README inflation), AP-CY40 ($L$-value
triangulation).

## Section 9: Plancherel super-quasi-Hopf

### AP-V2-16: Plancherel measure on simples vs projectives
**Severity**: High
**Trigger**: Defining the Plancherel measure on simples $L_\lambda$
of a non-semisimple chiral quantum group, without replacing by
projective covers $P_\lambda$.
**Ghost**: Non-semisimple Kerler--Lyubashenko MTCs force the
Plancherel measure to be supported on projective covers $P_\lambda$,
not simples $L_\lambda$; simples have zero quantum dimension in the
non-semisimple case, so the measure degenerates.
**Error**: Semisimple-only error: in a semisimple MTC, simples equal
their own projective covers; in non-semisimple MTCs (e.g., small
quantum group at root of unity), $P_\lambda \neq L_\lambda$ and
$\dim_q(L_\lambda) = 0$ for non-projective simples.
**Correction**: Construction: $\mu^{\mathrm{Pl}}_{\mathrm{KL}}(\lambda)
= \dim_q(P_\lambda) / \sum_\mu \dim_q(P_\mu)$; for small quantum
group at root of unity, this gives a non-degenerate measure on
$\{P_\lambda\}$. Remark: Kerler--Lyubashenko modified trace
$\mathrm{tr}^{\mathrm{mod}}_{P_\lambda}$ encodes the projective
measure.
**Primary**: Kerler--Lyubashenko 2001 *Non-Semisimple Topological
Quantum Field Theories for 3-Manifolds with Corners* (projective covers);
Geer--Kujawa--Patureau-Mirand 2018 *Selecta Math* 24 (modified trace);
Costantino--Geer--Patureau-Mirand 2014 *Quantum Topol* 5.
**Inscription**: `chapters/theory/plancherel_superquasihopf.tex`,
Thm 10.1 (\texttt{thm:plancherel-projective}).
**Cross-ref**: Vol III AP-CY37 (non-semisimple tilting MTC).

### AP-V2-17: Convergence mode
**Severity**: High
**Trigger**: Claiming the Plancherel decomposition
$\mathbbm 1 = \sum_\lambda \mu^{\mathrm{Pl}}(\lambda)\,\Pi_\lambda$
converges in norm or weak-$*$ sense.
**Ghost**: The convergence is Mittag--Leffler pro-stabilisation in
$\mathrm{Pro}(\mathrm{Mod}_H)$ (category of pro-objects), not norm or
weak-$*$; the pro-structure is essential for infinite-dimensional
chiral representations.
**Error**: Topological-convergence error: for infinite-dimensional
chiral representations, norm convergence is not available
(chiral-module norms do not exist canonically) and weak-$*$ requires
a dual pairing that the chiral category does not carry.
**Correction**: Construction: the Plancherel decomposition converges
as a pro-object limit
$\mathbbm 1 = \lim_{N \to \infty} \sum_{\lambda : |\lambda| \leq N}
\mu^{\mathrm{Pl}}(\lambda)\,\Pi_\lambda$ in $\mathrm{Pro}(\mathrm{Mod}_H)$
with Mittag--Leffler stability. Remark: Mittag--Leffler stabilisation
holds because the truncated sums are eventually constant in
$\mathrm{Pro}$ on any compactly-supported object.
**Primary**: Artin--Mazur 1969 *Etale Homotopy* (pro-objects);
Beilinson 1986 *K-Theory* 2 (pro-$\mathcal D$-modules);
Lurie 2017 *Higher Algebra* §6 (pro-categories).
**Inscription**: `chapters/theory/plancherel_convergence.tex`,
Prop 10.3 (\texttt{prop:plancherel-pro-ml}).
**Cross-ref**: Vol I AP142 (pro-object discipline); Vol III AP-CY33
(chain-level pro-structure).

### AP-V2-18: MO-extension super-quasi-Hopf via dynamical twist
**Severity**: Medium
**Trigger**: Attempting to extend the Maulik--Okounkov construction
to the super-quasi-Hopf setting without the Etingof--Kazhdan
dynamical twist.
**Ghost**: The MO-extension goes through Etingof--Kazhdan dynamical
twist $\partial_Z \log \Phi_{10}$ where $\Phi_{10}$ is the Igusa cusp
form; this specific twist is the unique $H^2$-rigidity class (up to
gauge).
**Error**: Non-canonical choice error: a generic dynamical twist
would break super-quasi-Hopf associativity; the EK-twist
$\partial_Z \log \Phi_{10}$ is uniquely determined (up to gauge) by
$H^2$-rigidity of the super-quasi-Hopf deformation.
**Correction**: Construction: EK dynamical twist $J(Z) = \exp(\hbar
\partial_Z \log \Phi_{10}(Z))$ where $Z \in \mathbb H_2$ Siegel upper
half-space; the twist satisfies the dynamical pentagon
equation and realises the unique $H^2$-rigidity class on super-quasi-Hopf
deformations. Remark: $H^2$-rigidity follows from
$\dim H^2(\mathfrak g; \mathfrak g^\vee) = 1$ with generator
$\partial_Z \log \Phi_{10}$.
**Primary**: Etingof--Kazhdan 1996 *Selecta Math* 2 (dynamical
twist); Maulik--Okounkov 2012 arXiv:1211.1287 (MO geometric
R-matrix); Gritsenko--Nikulin 1998 *Invent Math* 129 ($\Phi_{10}$
Igusa form).
**Inscription**: `chapters/theory/mo_superquasihopf.tex`, Thm 10.6
(\texttt{thm:mo-ek-dynamical}).
**Cross-ref**: Vol III AP-CY39 ($\Phi_{10}$ automorphic home).

## Section 10: Universal ratio-of-levels and three-faces identity

### AP-V2-19: Universal ratio-of-levels scope
**Severity**: High
**Trigger**: Asserting the universal ratio-of-levels identity
$\ell_X/\ell_Y = c_+(L_X)/c_+(L_Y)$ holds uniformly across all Mukai
lattice data, including the Leech--Conway row.
**Ghost**: The universal ratio ties Lusztig-level quantisation $\ell$
to Mukai-signature data $c_+$; the identity holds for the four
Gritsenko--Nikulin reflective rows plus the $24A_1$ Niemeier row, but
**breaks for Leech--Conway** (no Fricke involution on the positive-definite
Leech lattice).
**Error**: Scope inflation: the identity requires a Fricke involution
to compare $c_+(L_X)$ across rows; Leech (positive-definite) and
Conway $V^{s\natural}$ (super-lattice) do not admit Fricke involutions,
so the ratio is undefined.
**Correction**: Construction: universal ratio
$\ell_X/\ell_Y = c_+(L_X)/c_+(L_Y)$ holds for rows
$\{\mathrm{II}_{2,2}, \mathrm{II}_{10,2}, \mathrm{II}_{18,2}, \mathrm{II}_{26,2}\}$
(reflective Lorentzian, Fricke-involution available) plus $24A_1$
Niemeier. Fails on Leech row (no Fricke). Remark: the failure is
structural, not an error; it distinguishes reflective from
positive-definite data.
**Primary**: Lusztig 1990 *Geom Dedicata* 35 (levels); Mukai 1988
*Nagoya Math J* 116 (Mukai vector); Borcherds 1992 *Invent Math* 109
(Fricke involution on reflective lattices); Scheithauer 2017
*Commun Math Phys* 350.
**Inscription**: `chapters/theory/universal_ratio_levels.tex`,
Thm 11.1 (\texttt{thm:universal-ratio-levels-scope}).
**Cross-ref**: Vol III AP-CY34 (Conway $V^{s\natural}$ scope).

### AP-V2-20: Three-faces universal identity undercounting
**Severity**: Critical
**Trigger**: Treating any of the three faces (Vol I conductor $K$,
Vol II one-loop Quillen exponent, Vol III Borcherds weight) of the
universal identity $\hbar^2 \cdot K^{\kappa_{\mathrm{ch}}} = -1$ in
isolation.
**Ghost**: The three faces are three independent cross-volume
realisations of the **same** universal identity; at the distinguished
point $(K, \hbar^2) = (8, -1/8)$ all three reconcile on
$\hbar^2 K^{\kappa_{\mathrm{ch}}} = -1$, and any isolated-face
treatment is structural undercounting.
**Error**: Triangulation failure: treating only the Vol II face
loses the Vol I conductor arithmetic and the Vol III Borcherds weight
content; all three are needed for multi-path verification.
**Correction**: Construction: state the identity with all three faces
named: (a) Vol I conductor face $K = \mathrm{cond}(A)$, chiral algebra
$A$; (b) Vol II one-loop Quillen exponent
$\hbar^2 = \exp(2\pi i \cdot Q_{\mathrm{Quillen}})$;
(c) Vol III Borcherds weight $\kappa_{\mathrm{ch}} = \mathrm{wt}_{\mathrm{Borch}}(\Phi)$.
Distinguished point $(K, \hbar^2) = (8, -1/8)$:
$\hbar^2 \cdot K^{\kappa_{\mathrm{ch}}} = (-1/8) \cdot 8^{\kappa_{\mathrm{ch}}} = -1$
forces $\kappa_{\mathrm{ch}} = 1$ on the distinguished stratum.
Remark: the three faces are pairwise independent; their mutual
constraint determines $(K, \hbar^2, \kappa_{\mathrm{ch}})$ uniquely
on each reflective row.
**Primary**: Borcherds 1992 *Invent Math* 109 (Borcherds weight);
Quillen 1985 *Invent Math* 82 (one-loop determinant);
Scheithauer 2006 *Commun Math Phys* 266 (conductor).
**Inscription**: `chapters/theory/three_faces_universal.tex`, Thm 11.4
(\texttt{thm:three-faces-identity}).
**Cross-ref**: Vol I AP151 ($\hbar$ convention); Vol III AP-Vol-III-prop-1
(bare $\kappa$).

## Section 11: Chain-level vs cohomological

### AP-V2-21: Class M $E_3$-bar dimensionality
**Severity**: High
**Trigger**: Asserting the class-M $E_3$-bar cohomology is
infinite-dimensional at all lanes, or equivalently, claiming it is
$6^g$ at chain level.
**Ghost**: Class M $E_3$-bar cohomology is $6^g$ at the
**cohomological** level (finite, controlled by Euler-characteristic data
of $E_3$ to the power $g$), and **infinite-dimensional** at the
**chain** level (because the differential is not rationally exact);
mixing the lanes gives wrong genus-scaling.
**Error**: Lane conflation: Vol III AP-CY21 notes that $E_3$-bar for
class M is infinite-dimensional at chain level because $d_4$ survives;
at cohomological level, $d_4$ kills everything except a finite
$6^g$-dimensional quotient by the universal class-M coherence
spectral sequence collapse.
**Correction**: Construction: class M $E_3$-bar at chain level:
infinite-dimensional (Pattern 236 chain-level scope); at cohomological
level: $6^g$-dimensional. Remark: $6$ is the Euler-characteristic-weighted
$E_3$ contribution on a class-M generator with sign discipline;
$6^g$ matches the genus-$g$ tensor power structure.
**Primary**: Francis 2013 *Compositio* 149 §4 (chain vs cohomological
for $E_n$); Costello--Gwilliam 2021 Vol 2 §9 (BV chain vs
cohomological).
**Inscription**: `chapters/theory/class_m_e3_bar.tex`, Thm 12.1
(\texttt{thm:class-m-e3-bar-lane}).
**Cross-ref**: Vol III AP-CY21 ($E_3$ bar class M); Vol I AP33 (chain
vs rational).

## Section 12: Five-archetype landscape cross-volume

### AP-V2-22: Four-archetype vs five-archetype landscape
**Severity**: Critical
**Trigger**: Treating the G/L/C/M four-archetype landscape as
exhaustive of the $\kappa + \kappa^! \in \{0, 13, 250/3, 98/3\}$
complementarity table.
**Ghost**: The four-archetype classification extends to **five
archetypes** G/L/C/M/**B** via a BKM-crown row B at $\kappa + \kappa^! = 8$;
the B-row is populated by Borcherds--Kac--Moody algebras ($\Phi_{12}$
fake-Monster; $\Phi_{10}$ Igusa) and extends the Virasoro-ended
classification by one reflective row.
**Error**: Cross-volume scope inflation: the four-archetype table is
scoped to Vol I's standard landscape; the BKM crown row B appears
**only** when Vol III's CY-to-chiral functor $\Phi$ imports
Borcherds-product data, enlarging the complementarity table. Treating
four as exhaustive misses the crown content and undercounts by one
row.
**Correction**: Construction: extended complementarity table
$\kappa + \kappa^! \in \{0, 8, 13, 250/3, 98/3\}$ over five
archetypes G ($= 0$, Heisenberg), B ($= 8$, BKM crown), L
($= 13$, affine Kac--Moody), C ($= 98/3$, $\beta\gamma$), M
($= 250/3$, Virasoro) **[RETRACTED AP5 per W25.C6 cross-volume audit
2026-04-20: the per-row assignment L = 13 / C = 98/3 / M = 250/3 is
WRONG relative to Vol I `chapters/examples/landscape_census.tex` L1576-
1579 canonical: $\kappa(\widehat{\fg}_k) + \kappa(\widehat{\fg}_k^!) = 0$,
$\kappa(\beta\gamma_\lambda) + \kappa((\beta\gamma_\lambda)^!) = 0$,
$\kappa(\mathrm{Vir}_c) + \kappa(\mathrm{Vir}_c^!) = 13$; the values
$250/3$ and $98/3$ are W-extensions on the M-row
($\mathsf{M}$-ext $\mathcal{W}_3^k$ and $\mathrm{BP}_k$ per
landscape\_census.tex L1756, L1758), NOT the native L/C/M assignments.
Canonical per-row: G = 0 ($r_{\max} = 2$), L = 0 ($r_{\max} = 3$),
C = 0 ($r_{\max} = 4$), M = 13 ($r_{\max} = \infty$), B = 8
($r_{\max} = \infty$); $\mathsf{M}$-ext values $250/3, 98/3$ sit on
the M-row via principal-W and minimal-W specialisations. See Vol I
AP357 L1316-1327 + Vol II `first_principles_cache_comprehensive.md`
L3207-3211 for the correct per-row table.]** Remark: the B-row is
Vol II's contribution to the cross-volume synthesis; its value $8$
matches the Costello--Gwilliam 1-loop counterterm of AP-V2-8 and the
$(K, \hbar^2) = (8, -1/8)$ distinguished point of AP-V2-20.
**Primary**: Borcherds 1992 *Invent Math* 109 (BKM algebras);
Gritsenko--Nikulin 1998 *Invent Math* 129 ($\Phi_{10}$ Igusa);
Scheithauer 2017 *Commun Math Phys* 350 (reflective crown);
Costello--Gwilliam 2021 Vol 2 §7 (1-loop value).
**Inscription**: `chapters/theory/five_archetype_bkm_crown.tex`,
Thm 13.2 (\texttt{thm:five-archetype-landscape}).
**Cross-ref**: Vol I Theorem C (derived-centre complementarity);
Vol III AP-CY8 (Borcherds $\neq$ bar Euler), AP-Vol-III-prop-2
($N=1$ coincidence inflation).

## Statistics

| Category                                              | Count |
| ----------------------------------------------------- | ----: |
| Bar-cobar structure (AP-V2-1--3)                      |     3 |
| $(\infty,1)/(\infty,2)$-categorical (AP-V2-4--5)      |     2 |
| Chiral Kontsevich formality (AP-V2-6--7)              |     2 |
| BD factorisation compatibility (AP-V2-8)              |     1 |
| Theorem H Vol II $\leftrightarrow$ Vol III (AP-V2-9)  |     1 |
| Cyclic chiral homology (AP-V2-10--11)                 |     2 |
| NC Hodge Mukai--Hochschild (AP-V2-12--13)             |     2 |
| CY-4 HK-restricted $\Phi_4$ (AP-V2-14--15)            |     2 |
| Plancherel super-quasi-Hopf (AP-V2-16--18)            |     3 |
| Universal ratio and three-faces (AP-V2-19--20)        |     2 |
| Chain-level vs cohomological (AP-V2-21)               |     1 |
| Five-archetype landscape (AP-V2-22)                   |     1 |
| **Total catalogued (Wave 20-24)**                     |  **22** |
| Critical severity                                     |    10 |
| High severity                                         |    11 |
| Medium severity                                       |     1 |
| Low severity                                          |     0 |

Three Wave-20-24 anti-patterns have already driven 5+ independent
fixes across the Vol II chapter corpus: AP-V2-1 (Swiss-cheese
octachotomy, 7 fixes across `sc_chtop_octachotomy.tex` +
`heptagon_edge_34_*.tex` + downstream cross-references); AP-V2-9
(Theorem H scope hinge, 6 fixes reconciling Vol II/Vol III ambient
statements); AP-V2-22 (five-archetype B-row, 5 fixes enlarging the
complementarity table across `chiral_center_theorem.tex` and
cross-volume references).

---

## Historical K3-BKM / 3d HT-QFT anti-patterns (Waves 14-22, Vol II frame)

- **V2-AP56 -- $(c_{4d}, c_{2d})$ central-charge reversal (Critical).**
  Wave-14 retraction to $(26, -312)$ fails $\mathrm{SU}(2)$ $N_f = 4$
  cross-check. Canonical $(107/6, -214)$ at $n = 24$ via $(5n-13)/6$;
  WOV-2 locks.
- **V2-AP57 -- Monster vs Fake-Monster rank (Critical).**
  Monster rank 2 ($\mathrm{II}_{1,1}$), Fake-Monster rank 26
  ($\mathrm{II}_{25,1}$), K3-BKM rank 3 ($\Lambda^{2,1}_{II}$).
- **V2-AP58 -- $c_3$ normalisation (High).**
  Bruinier reduced-class: $c_3 = -8$. Conversion factor $-22032$ to
  Gritsenko--Nikulin Cartan-matrix convention.
- **V2-AP59 -- Umbral Niemeier labelling (High).**
  $(N-1) \mid 24$ for $A_{N-1}$ umbral; $N = 6$ reanchored to $6D_4$.
- **V2-AP60 -- $\zeta(3,3,3,3)$ value (High).**
  $0.000295999\ldots$, not $0.0028565$. Depth-4 MZV at weight 12.
- **V2-AP61 -- $c(1,2,\pm 2) = -2$ (High).**
  Integer heterotic lift; sign-flip corrected via three-path
  cross-check.
- **V2-AP62 -- $H_4$ vs $H_8$ (High).**
  $H_4$ is $(2,2)$-isogeny with $\mathrm{End} \supset \mathbb{Z}[2i]$.
- **V2-AP63 -- Theorem B scope (Medium).**
  Exclude all admissible $H_n$, not just $H_1 \cup H_4$.
- **V2-AP64 -- Borcherds vs Gritsenko weight (Medium).**
  K3-BKM denominator $\Delta_5$ (Grit additive), not $\Phi_{12}$
  (Borch multiplicative).
- **V2-AP65 -- Two-$\hbar$ (High).**
  $\hbar^{\mathrm{Drinfeld}} \ne \hbar^{\mathrm{BV}}$ semantically.
  Name the $\hbar$ at every site.
- **V2-AP66 -- Sylvester vs Feingold--Frenkel (Medium).**
  Use eigenvalue signature $\{+4, +4, -2\}$, not Sylvester
  $(2, 0, -32)$ trap.
- **V2-AP67 -- Schmidt parameter pair (Medium).**
  $(17/2, 15/2)$ for $\Delta_{10}$; $(7/2, 5/2) \otimes \mathrm{sgn}_\mathbb{R}$
  for $\Delta_5$ Maass-spin cover.
- **V2-AP68 -- CoHA $\ne$ chiral (Critical).**
  $\mathbf{H}_{\Delta_5} = \Phi_3(\mathcal{D}_\hbar(\mathcal{Y}^{\mathrm{Hall}}(\mathrm{CoHA})))$.
- **V2-AP69 -- $\kappa_\bullet$ Künneth (High).**
  $\kappa_{\mathrm{cat}}(K3 \times E) = 0$ (multiplicative), not 2.


---

## V2-AP1 through V2-AP20: E_1/E_inf Locality Hierarchy (re-extracted 2026-04-21)

The canonical Vol II register of locality and pole-structure anti-patterns,
re-extracted from the E_1/E_inf Locality Hierarchy section of the legacy
CLAUDE.md. These were the scrub-cascade casualty identified in V2-AP51.

**Canonical statement.** $E_1$ vs $E_\infty$ is about LOCALITY, not poles.
$E_\infty$ = LOCAL = $\Sigma_n$-equivariant; $E_1$ = NONLOCAL. OPE poles
are compatible with $E_\infty$. Three-tier picture WITHIN $E_\infty$:
(i) pole-free ($R = \tau$, e.g. BD "commutative chiral algebra" — STRICT
SUBCLASS of $E_\infty$); (ii) VA with poles ($R \ne \tau$ but DERIVED
from local OPE); (iii) genuinely $E_1$ ($R \ne \tau$, independent input).
Tiers (i)+(ii) BOTH $E_\infty$. **Kac--Moody, Virasoro, Heisenberg, $W$-algebras
are ALL $E_\infty$** — poles do not break $E_\infty$. Discriminant
$E_1$ vs $E_\infty$ is PROVENANCE of $R(z)$, not its value. BD do NOT study $E_1$.
$E_1$-chiral algebra is a NEW concept of THIS manuscript.

- **V2-AP1**: $E_\infty$-chiral INCLUDES ALL vertex algebras. KM, Virasoro,
  Heisenberg are ALL $E_\infty$-chiral. NEVER write "VAs are not
  $E_\infty$-chiral."
- **V2-AP2**: $R(z) \ne \tau$ does NOT imply genuinely $E_1$-chiral. For
  $E_\infty$ with poles, $R(z)$ is derived from local OPE; for $E_1$-chiral,
  $R(z)$ is independent input. Discriminant is PROVENANCE, not value.
- **V2-AP3**: Three bars: $B^{FG}$ (zeroth pole only) $\ne B^\Sigma$ (all poles
  $+$ coinvariants) $\ne B^{\mathrm{ord}}$ (all poles $+$ ordering).
  Maps: $B^{\mathrm{ord}} \to B^\Sigma$ (coinvariants);
  $\mathrm{gr}(B^\Sigma) \to B^{FG}$ (filtration).
- **V2-AP4**: Ordered-to-unordered descent is $R$-matrix twisted:
  $B^\Sigma_n = (B^{\mathrm{ord}}_n)^{R\text{-}\Sigma_n}$.
  Naive quotient only for pole-free.
- **V2-AP5**: NEVER equate $E_\infty$-chiral with "no OPE poles." BD
  "commutative chiral algebra" (no poles) is a STRICT SUBCLASS of
  $E_\infty$-chiral.
- **V2-AP6**: BD do NOT study $E_1$. $E_1$-chiral algebra = NEW concept
  of THIS manuscript.
- **V2-AP7**: Heisenberg $R$-matrix $= \exp(k \hbar / z)$, NOT trivial.
  Collision residue $k/z$. Monodromy $\exp(-2\pi i k)$.
- **V2-AP8**: NEVER add restrictive parenthetical glosses. "$E_\infty$-chiral
  ($=$ BD commutative $=$ no poles)" NARROWS the term.
- **V2-AP9**: NEVER say VA "is not $E_\infty$-chiral." KM, Virasoro,
  Heisenberg, $W$-algebras are ALL $E_\infty$-chiral. Poles do not break
  $E_\infty$-chirality.
- **V2-AP10**: NEVER "$E_\infty$-chiral implies $R(z) = \tau$" without
  pole-free qualifier. Correct: "For POLE-FREE $E_\infty$-chiral, $R(z) = \tau$."
- **V2-AP11**: NEVER conflate $E_\infty$-chiral with BD "commutative."
  BD Chapter 4 "commutative" = pole-free = strict subclass.
- **V2-AP12**: $E_1$- vs $E_\infty$-chiral is about LOCALITY, not poles.
- **V2-AP13**: NEVER trust agent claim "VAs are not $E_\infty$-chiral."
  This exact error caused cascading damage.
- **V2-AP14**: NEVER oscillate between conventions in single session.
- **V2-AP15**: NEVER edit $E_1$/$E_\infty$ language without author confirmation.
- **V2-AP16**: Three-tier picture is WITHIN $E_\infty$-chiral, not a
  division between $E_\infty$-chiral and $E_1$-chiral. (i) + (ii) both
  $E_\infty$-chiral. Only (iii) is $E_1$-chiral.
- **V2-AP17**: NEVER revert file based on false premise. Surgical removal only.
- **V2-AP18**: Author's explicit statements override agent literature searches.
- **V2-AP19**: NEVER batch-propagate unverified corrections. ONE edit,
  verify, THEN propagate.
- **V2-AP20**: NEVER add "in the sense of [reference]" without verification.

---

## V2-AP21 through V2-AP48: legacy extension (re-extracted 2026-04-21)

Continuation of the legacy V2-AP register re-extracted from
`notes/claude_md_legacy_20260418.md` and
`~/chiral-bar-cobar/notes/cross_volume_aps.md`.

### Section L.3: Hierarchy and chromatic (V2-AP21--V2-AP24)

- **V2-AP21**: PVA $\ne$ $P_\infty$-chiral. PVA is the CLASSICAL
  shadow (descend to cohomology); $P_\infty$-chiral is the HOMOTOPY
  intermediate. Opposite directions.
- **V2-AP22**: Full hierarchy: commutative associative $\subset$ PVA
  $\subset$ $E_\infty$-chiral $\subset$ $P_\infty$-chiral $\subset$
  $E_1$-chiral. Bar/Koszul at the $E_\infty$ and $E_1$ levels only.
- **V2-AP23**: Chromatic: the classical theory is height $0$.
  $L_{K(n)}(B(A)) = 0$ for $n \ge 1$. Pole order $\ne$ chromatic
  height.
- **V2-AP24**: $S$-transform (closed, complex structure) $\ne$
  Wick rotation of $R$ (open, $E_1$ ordering). Different algebraic
  data.

### Section L.4: Empirical anti-patterns from error archaeology (V2-AP25--V2-AP39)

- **V2-AP25**: Complex-analytic sign. $\mathrm{Im}(f) = (f - \bar f)
  / (2 i)$; $\bar\partial \,\mathrm{Im}(f) = (i/2) \bar\partial \bar f$,
  NOT $1 / (2 i)$. Identity $-1 / (2 i) = i / 2$ is a common
  confusion. Verify at EACH propagation site.
- **V2-AP26**: NEVER hardcode Part/chapter numbers. Always
  `\ref{part:...}`. 24+ stale refs after $10 \to 7$ Part
  restructuring.
- **V2-AP27**: Duplicated content across files FORBIDDEN. Use
  `\input{}` or `\ref{}`, never copy-paste theorem environments.
- **V2-AP28**: Test expected values from 2+ independent sources with
  documented derivation. Engine and test from the same mental model
  share the same error. $\lambda_3 = 1/82944$ was WRONG (correct:
  $31/967680$) because both used the same faulty computation.
- **V2-AP29**: AI-slop cleanup mandatory post-generation. Grep for:
  moreover, additionally, notably, crucially, remarkably, "it is
  worth noting", em dashes, "We now", passive "can be shown."
- **V2-AP30**: After architecture restructuring: grep
  `Part~[IVXL]` in `chapters/`; also `\ref{part:` to verify targets.
- **V2-AP31**: AP4 at write time. Before `\begin{proof}`, verify the
  preceding environment is theorem/proposition/lemma with
  `\ClaimStatusProvedHere`. If conjecture: use
  `\begin{remark}[Evidence]` instead. 25-instance fix commit.
- **V2-AP32**: Standalone-document artifact leak. Chapter files
  `\input`'d into `main.tex` MUST NOT contain `\title{}`,
  `\begin{abstract}`, `\tableofcontents`, `\date{}`, `\author{}`.
- **V2-AP33**: `RECTIFICATION-FLAG` must NOT become permanent debt.
  Zero tolerance for unresolved flags at session end.
- **V2-AP34**: Divided-power convention in $\lambda$-brackets. Vol II
  uses $\{T_\lambda T\} = (c / 12) \lambda^3$ (divided power). OPE
  mode $T_{(3)} T = c/2 \mapsto (c/2)/3! = c/12$. Grep
  `c/2.*lambda^3`; if found, almost certainly wrong (should be
  $c/12$). $\mathcal W_3$: $c/3 \cdot \lambda^5$ wrong; correct
  $c/360$.
- **V2-AP35**: Unresolved logical connectives after correction.
  Audit all "therefore"/"hence"/"it follows" within 5 lines.
- **V2-AP36**: Terminology rename atomicity. Enumerate all variant
  forms; grep all three volumes; complete replacements in a SINGLE
  commit. "shadow Postnikov tower" $\to$ "shadow obstruction tower"
  needed 5 commits.
- **V2-AP37**: Arakelov form canonical normalisation.
  $\omega_1 = (i / (2 \,\mathrm{Im}\,\tau)) dz \wedge d\bar z$
  (integral $= +1$);
  $\omega_{\mathrm{Ar}} = -(\pi / \mathrm{Im}\,\tau) dz \wedge d\bar z$
  (integral $= -1$); $\omega_1 = -\omega_{\mathrm{Ar}} / (2\pi)$.
  Same error fixed THREE times.
- **V2-AP38**: Phantom label retirement. After chapter migration,
  366 phantoms across 2 commits after Vol I $\to$ Vol II migration.
- **V2-AP39**: Macro portability after migration. After migrating
  Vol I $\to$ Vol II: compile; grep "Undefined control sequence";
  add `\providecommand` for each. 7 macros across 2 commits.

### Section L.5: Wave-17-19 synthesis legacy (V2-AP40--V2-AP48)

- **V2-AP40**: NO anti-pattern tags or metadata leakage into the
  typeset manuscript. AP identifiers, FM identifiers, session codes,
  commit hashes, meta-stamps, visible-label leaks, session dates
  in body prose MUST NOT appear in any `.tex` file compiled into
  the monograph PDF or standalone paper PDF. Six subclauses:
  V2-AP40a (bibkey), V2-AP40b (label), V2-AP40c (index/hyperref),
  V2-AP40d (environment/section title), V2-AP40e (monospace
  filename), V2-AP40f (math-mode label leak).
- **V2-AP41**: Internal theorem-statement-vs-proof formula
  inconsistency. A main theorem states one formula; the immediate
  proof cites a different formula diverging at boundary cases.
  Example: $\beta_{W_N} = (N+1)(N+2)/2$ stated, Stirling-dominance
  proof cites $\beta_N = 12(H_N - 1)$; agree at $N = 2, 3$;
  diverge $N \ge 4$.
- **V2-AP42**: $(\infty, 2)$-adjunction bar-cobar. 1-truncation does
  NOT capture Auslander--Reiten triangular structure on
  compatibility 2-morphisms. Specify $(\infty, 1)$ or
  $(\infty, 2)$; never elide 2-morphism coherence.
- **V2-AP43**: Chiral Kontsevich formality fails on stratified
  $\overline{\mathcal A_2}$. Francis--Gaitsgory $+$ Tamarkin
  transport holds ONLY on Koszul locus $U^{\mathrm{adm}}$; Humbert
  strata are OBSTRUCTED. Name the ambient and the locus at every
  invocation.
- **V2-AP44**: Swiss-cheese colouring is an OCTACHOTOMY. Eight
  bar-cobar ambients correspond to eight specific
  $\mathsf{SC}^{\mathrm{ch,top}}$ colourings; treating them as
  interchangeable is a category error.
- **V2-AP45**: CY-4 HK-restricted $\Phi_4$. Generic CY-4 $\Phi_4$
  BLOCKED by Kapustin--Rozansky--Saulina 3d/4d dichotomy;
  scope-restricted to hyperkähler CY-4 via twistor-$S^1$ reduction.
- **V2-AP46**: Francis--Gaitsgory cyclic chiral homology needs
  $E_1$-operadic ambient upgrade; classical Feigin--Tsygan verbatim
  transport is a FIVE-STEP FAILURE.
- **V2-AP47**: NC Hodge Kaledin degeneration scope. Smooth-proper DG
  scope $\ne$ factorisation $\infty$-category scope; Humbert strata
  carry non-zero Heegner--Bruinier obstruction classes feeding
  higher-page $d_r$.
- **V2-AP48**: Theorem H hinge enlarges under $\Phi_d$. Ambient-curve
  $\mathrm{ChirHoch}$ concentration $\{0, 1, 2\}$ ENLARGES to
  $\{0, 1, 2, d\}$ under $\Phi_d$.

---

## Waves 1-13 pre-synthesis anti-patterns (V2-AP100--V2-AP125)

Distilled from the Vol II pre-Wave-14 programme record: SC foundations,
Dimofte integration, Swiss-cheese / seven faces, and the 3D QG climax.
Numbering is disjoint from `V2-AP1--V2-AP48` and `AP-V2-1--AP-V2-22`.

### Section W.1: $\mathsf{SC}^{\mathrm{ch,top}}$ foundations (V2-AP100--V2-AP105)

- **V2-AP100**: $\mathsf{SC}^{\mathrm{ch,top}}$ is a TWO-COLOURED
  dioperad with directional restriction (closed output cannot
  receive open input). Dunn additivity does NOT apply; the correct
  promotion is through Voronov's degeneration plus the
  product-forced locality hypothesis, NOT Dunn.
- **V2-AP101**: $B(A)$ is an $E_1$ COASSOCIATIVE coalgebra, NOT a
  $\mathsf{SC}^{\mathrm{ch,top}}$-coalgebra (cf. Vol I AP165). The
  $\mathsf{SC}$ structure emerges in the DERIVED CENTRE PAIR
  $(C^\bullet_{\mathrm{ch}}(A, A), A)$, not on $B(A)$.
- **V2-AP102**: Bar differential on the curve is HOLOMORPHIC
  FACTORISATION; coproduct on $\mathbb R$ is TOPOLOGICAL
  FACTORISATION. The $\mathsf{SC}^{\mathrm{ch,top}}$ datum is the
  pairing over $C \times \mathbb R$. Conflation collapses
  bulk-boundary.
- **V2-AP103**: Homotopy retract is DATA. Chain homotopy,
  Mittag--Leffler witness, or explicit MC element must be named;
  disc $\ne$ point, $\mathbb P^1 \ne \mathbb A^1$ (cf. Vol I
  AP142 / B53 / FM23).
- **V2-AP104**: Costello--Gwilliam bar differential $\ne$ chiral
  bar differential verbatim. CG's holomorphic factorisation
  framework provides the AMBIENT; the chiral bar differential is
  constructed via explicit OPE poles on the curve.
- **V2-AP105**: $E_{k+2}^{\mathrm{top}}$ topologisation ladder
  requires conformal vector at NON-CRITICAL level. Affine KM
  reaches $E_3^{\mathrm{top}}$ (PROVED; cf. Vol I
  AP167 / B59 / FM28); general chiral-with-conformal is
  CONJECTURAL. $\mathcal W_N$ reaches $E_{N+1}^{\mathrm{top}}$;
  $\mathcal W_\infty$ is the Platonic endpoint.

### Section W.2: Dimofte integration and Drinfeld double (V2-AP106--V2-AP112)

- **V2-AP106**: Slab $=$ bimodule. Dimofte slab geometry on
  $C \times [0, 1]$ produces an $(A, B)$-bimodule structure; the
  two ends carry distinct chiral algebras.
- **V2-AP107**: Drinfeld double identification is a CONJECTURE.
  $Z^{\mathrm{der}}_{\mathrm{ch}}(A) \simeq \mathcal D(A)$ only on
  Koszul locus; beyond the locus, open.
- **V2-AP108**: Line operators live on $\mathbb R$, not on the
  curve. Topological-factorisation coproduct along $\mathbb R$
  generates line operators; holomorphic-factorisation bar
  differential on $C$ does not.
- **V2-AP109**: $\hat Z$-invariants are CONJECTURAL cross-volume
  objects. Gukov--Pei--Putrov--Vafa $\hat Z$ sit at the 3d-boundary
  of Vol II's HT setup; identification with Vol III's
  Borcherds-side quantities is open.
- **V2-AP110**: Holomorphic blocks $\ne$ conformal blocks. A
  holomorphic block on $C \times \mathbb R$ is a chiral bundle
  section; a conformal block is a BD-factorisation-sheaf section
  on $\mathrm{Ran}\,X$. Coincide on genus-$0$ degeneration only.
- **V2-AP111**: 3D quantum gravity CLIMAX (Part VI) requires ALL
  prior parts. Asserting the climax from a subset is a
  dependency-inflation error.
- **V2-AP112**: Genus-1 Bernard--Felder propagator involves theta
  functions and quasi-periodicity under $z \to z + \tau$. The
  genus-$0$ limit is a DEGENERATION, not a specialisation.

### Section W.3: Early Swiss-cheese and seven faces (V2-AP113--V2-AP119)

- **V2-AP113**: Voronov 1998 Swiss-cheese $\ne$
  $\mathsf{SC}^{\mathrm{ch,top}}$ of Vol II. Voronov uses
  $\mathrm{FM}_{k,m}(\mathbb H, \partial \mathbb H)$ (boundary
  fixed to $\partial \mathbb H$); $\mathsf{SC}^{\mathrm{ch,top}}$
  uses $\mathrm{FM}_k(\mathbb C) \times E_1(m)$ (product). The
  latter is a DEGENERATION of the former.
- **V2-AP114**: Seven faces of $r(z)$ chain F1--F7 is a STAR
  GRAPH: F1 hub, six spokes; plus external F5 $\leftrightarrow$ F6
  (Drinfeld--STS 1983--85). NEVER "seven equivalent."
- **V2-AP115**: F5 (Yangian, level-independent) vs F7 (Gaudin,
  level-shifted by $k + h^\vee$) CONFLATED at multiple sites. At
  $k = 0$: Yangian $r(z) \ne 0$; Gaudin $r(z) = \Omega / h^\vee
  \ne 0$; level-stripped $r = 0$.
- **V2-AP116**: F1 $\Leftrightarrow$ F4 is INJECTION, not
  bijection. $r(z)$ determines the KZ 1-form; the Drinfeld
  associator $\Phi_{\mathrm{KZ}}$ has MZV coefficients invisible
  to $r(z)$.
- **V2-AP117**: Heptagon edge 3-4 is associator-INDEPENDENT.
  Associator enters in edge 2-3 through $\mathrm{GRT}_1$-equivariant
  lift.
- **V2-AP118**: Heptagon edge 3-4 $\mathrm{GRT}$-equivariant heal.
  Associator-dependence in edge 2-3 is healed by passage to
  $\mathrm{GRT}_1$-equivariant category.
- **V2-AP119**: Genus-1 seven faces require KZB connection as the
  elliptic analogue of KZ. Felder propagator provides F6 at
  genus 1; genus $\ge 2$ degenerates via pair-of-pants
  factorisation.

### Section W.4: Drinfeld-double gap and frontier (V2-AP120--V2-AP125)

- **V2-AP120**: Drinfeld double gap.
  $Z^{\mathrm{der}}_{\mathrm{ch}}(A) \simeq \mathcal D(A)$ OPEN in
  general; proved on Koszul locus only. Gap is the OBSTRUCTION
  CLASS in $H^2(\mathrm{grt}_1; \widehat{\mathrm{Imag}})$.
- **V2-AP121**: Conway $V^{s\natural}$ is NOT an independent
  $\Psi$-image row. Duncan 2007 \emph{Duke Math.\ J.}\ 139
  constructs $V^{s\natural} = A(\Lambda_{24})^+ \oplus
  A(\Lambda_{24})^{\mathrm{tw}, +}$ as $\mathbb Z / 2$-orbifold of
  the 24-generator fermionic VOA on the Leech lattice;
  super-twin of $V^\natural$ on the Monster row. Three earlier
  defects: VENUE (Duke 139, not MRL 14), CONSTRUCTION (Leech
  fermionic VOA, not $E_8$ super-lattice), SIGNATURE
  ($c_+(\Lambda_{24}) = 24$, not $0$).
- **V2-AP122**: Six routes to $G(K3 \times E)$ are six DIFFERENT
  constructions ($\Phi$, Borcherds lift, lattice VOA, Kummer,
  sigma model, BLLPR), NOT six applications of the same functor.
- **V2-AP123**: Class M $E_3$-bar cohomology is $6^g$
  (COHOMOLOGICAL level), infinite-dimensional (CHAIN level).
  Mixing lanes gives WRONG genus-scaling. $g = 1$ bar profile
  $[0, 3, 3, 0]$.
- **V2-AP124**: $\zeta(3, 3, 3, 3) = 0.000295999\ldots$, NOT
  $0.0028565$. Depth-4 MZV at weight 12. Brown 2011 motivic MZV
  basis confirmation.
- **V2-AP125**: Fricke LDP variance $\sigma_k^2 = \pi / (2 \sin^2
  \theta_k^*)$ (inverse-density), NOT $(\pi^2 / 2) \cos^{-2}
  \theta_k^*$ (density-curvature, DIVERGES at finite-density
  nodes).

### Section W.5: Single-valued MZV scope on chiral-Hochschild periods (V2-AP126)

- **V2-AP126 -- Single-valued MZV scope of chiral-Hochschild
  periods (Critical).**
  (a) **Ghost.** The Deligne-Goncharov 2005 *Ann Sci ENS* 38
  mixed-Tate motivic framework gives the full motivic Galois
  $\mathrm{grt}_1^{\mathrm{mot}}$ acting on the motivic MZV ring
  $\mathrm{MZV}^{\mathrm{mot}}$; the chiral-Hochschild period
  $2\mathrm{Vol}(E)(2\pi i)^3$ is a genuine period with a natural
  motivic home. The Vol II one-loop Quillen exponent (Path E) and
  the Costello-Gwilliam factorisation-cyclic-homology $\chi_3$
  pairing on $\mathsf{SC}^{\mathrm{ch,top}}$-algebras $a\text{-}priori$
  expand in $\mathrm{MZV}^{\mathrm{mot}}$.
  (b) **Precise error.** Asserting that the chiral-Hochschild
  period identity $\chi_3 = 2\mathrm{Vol}(E)(2\pi i)^3$ lies in the
  full motivic ring $\mathrm{MZV}^{\mathrm{mot}}$ is a scope
  inflation. The Arnold forms $\eta_{ij} = d\log|z_{ij}|^2$ that
  witness the chain-level $\chi_3$ cocycle are **single-valued
  real**; the period pairing factors through Brown 2013 *Ann Math*
  175 projection $\mathrm{proj}: \mathrm{MZV}^{\mathrm{mot}} \to
  \mathrm{MZV}^{\mathrm{sv}}$. At weight 2, $\zeta^{\mathrm{sv}}(2)
  = 0$; conflating the two rings overcounts admissible periods and
  would predict $\zeta(2)$-weighted Hochschild contributions that
  do not survive the single-valued projection.
  (c) **Correct.** Chiral-Hochschild periods live in
  $\zeta^{\mathrm{sv}}$ (Brown 2013 single-valued MZVs), **not**
  in $\mathrm{grt}_1^{\mathrm{mot}}$-stable full motivic MZVs.
  Three sites must be distinguished:
  **chain-level** (explicit $\eta_{ij}$-integrals on
  $\mathrm{Conf}_n(X)$, rational-coefficient),
  **motivic** ($\mathrm{MZV}^{\mathrm{mot}}$ target of the period
  map),
  **single-valued** ($\zeta^{\mathrm{sv}}$ image under Brown's
  projection). Canonical identifications:
  $\zeta^{\mathrm{sv}}(2) = 0$,
  $\zeta^{\mathrm{sv}}(2k+1) = 2\zeta(2k+1)$ at odd weight, and
  $\zeta^{\mathrm{sv}}$ is a **proper** subring at depth $\ge 2$
  (Schnetz 2014 *Commun Num Theor Phys* 8; Panzer 2015
  *Commun Num Theor Phys* 9).
  **Vol II reading.** The single-valued scope constrains Vol II's
  one-loop Quillen exponent / Costello-Gwilliam
  factorisation-cyclic-homology reading of the $\chi_3$ pairing.
  The Vol II cyclic chiral homology (Path E) lands in
  $\zeta^{\mathrm{sv}}$ via the Francis-Gaitsgory
  cyclic-factorisation trace composed with the Brown single-valued
  projection. The Theorem H amplitude bound
  $\mathrm{ChirHoch}^\bullet \in \{0, 1, 2\}$ is recovered as a
  **single-valued consequence** of $\zeta^{\mathrm{sv}}(2) = 0$,
  not imposed as a separate concentration axiom.
  **Counter**: never equate chiral-Hochschild periods with
  $\mathrm{MZV}^{\mathrm{mot}}$; always name the Brown 2013
  projection and land in $\zeta^{\mathrm{sv}}$.
  Cross-reference Vol I AP901 and Theorem
  `thm:sv-scope-restriction-chiralhoch` in
  `/Users/raeez/chiral-bar-cobar/chapters/theory/motivic_shadow_tower.tex`.
  Related: AP888 (shadow-ChirHoch bridge) and the seven-path
  $\chi_3$ comparison theorem.
  **Primary citations**: Brown 2013 "Mixed Tate motives over
  $\mathbb Z$" *Ann Math* 175; Brown 2013 *Ann Sci ENS* 46
  single-valued multiple polylogarithms; Schnetz 2014 single-valued
  zeta; Deligne-Goncharov 2005; Panzer 2015 single-valued
  algorithms.

---

## Cross-reference index

- Legacy `V2-AP1--V2-AP48` (this file: sections directly above
  plus re-extracted V2-AP1-20 block) is verbatim from
  `notes/claude_md_legacy_20260418.md` and
  `~/chiral-bar-cobar/notes/cross_volume_aps.md`. Number preservation:
  `V2-APN` is the historical trigger; `AP-V2-N` (top of this file)
  is the current Platonic-synthesis scheme. Independent numberings.
- Waves 1-13 (sections W.1--W.4, `V2-AP100--V2-AP125`) distilled
  from session archives in `~/chiral-bar-cobar/notes/session_*`
  $+$ `project_dimofte_full_integration.md`
  $+$ `project_sc_chtop_concrete_programme.md`
  $+$ `project_preface_geometric_escalation.md`.
- Waves 14-19 entries appear as `V2-AP56--V2-AP69` distilled from
  `~/calabi-yau-quantum-groups/notes/ADJUDICATION_LEDGER_WAVES_14_TO_19.md`.
- Waves 20-24 entries appear as `AP-V2-1--AP-V2-22` distilled from
  `~/chiral-bar-cobar/notes/GRAND_SYNTHESIS_WAVES_20_22.md`.

## Per-wave coverage table

| Wave(s)     | Count | Scheme                | Focus                                                      |
| ----------- | ----: | --------------------- | ---------------------------------------------------------- |
| Waves 1-9   |    16 | V2-AP1--V2-AP16       | $E_1 / E_\infty$ locality hierarchy                        |
| Waves 1-9   |     4 | V2-AP17--V2-AP20      | Process and workflow guards                                |
| Waves 10-13 |     4 | V2-AP21--V2-AP24      | Hierarchy and chromatic                                    |
| Waves 10-13 |    15 | V2-AP25--V2-AP39      | Empirical error archaeology                                |
| Waves 17-19 |     9 | V2-AP40--V2-AP48      | Wave-17-19 synthesis legacy                                |
| Waves 1-13  |    26 | V2-AP100--V2-AP125    | SC foundations, Dimofte, Swiss-cheese, seven faces, climax |
| Waves 14-19 |    14 | V2-AP56--V2-AP69      | K3-BKM / 3d HT-QFT adjudication                            |
| Waves 20-24 |    22 | AP-V2-1--AP-V2-22     | Octachotomy, Theorem H hinge, BKM crown, five-archetype    |
| Wave 27     |     1 | AP-V2-23              | Chenevier determinant (not pseudo-character)               |
| Wave 28     |     1 | V2-AP126              | Single-valued MZV scope on chiral-Hochschild periods       |
| Wave 29     |     1 | AP-V2-24 / V2-AP127   | Humbert--Heegner admissibility filter $n\equiv 3,5\pmod 8$ |
| Wave 14 cross | 5 | AP-V2-25--AP-V2-29 / V2-AP128--V2-AP132 | hCS / KT formality torsor, MNOP centre trace, Swiss-cheese compatibility, $\mathbb E_n$ ambient vs category, non-isolated LG critical locus (Vol III AP-CY160--172 cross-append) |
| **Total**   | **118** | legacy + synthesis  | Vol II anti-patterns                                       |

---

### AP-V2-23: Chenevier determinant, not Taylor-Wiles pseudo-character (arithmetic anchor for $\mathbf H_{\Delta_5}$)

**Trigger**: "Taylor 1991 pseudo-character $S^{\mathrm{ps}} :
\mathbb T^{\mathrm{par}}_1 \to \mathcal O_E$ with axioms (symmetry /
multiplicativity / dimension)" attached to the Galois-side invariants
of the BKM crown algebra $\mathbf H_{\Delta_5}$.

**Ghost**: the Taylor-Wiles pseudo-character $S^{\mathrm{ps}}$ is a
real object (Taylor 1991 Duke 63 Thm 2.1; Rouquier 1996); the
Hecke-algebra 4-tuple $(S_1, S_2, S_3, S_4)$ computed from the
Saito-Kurokawa lift Satake parameters of $\Delta_{10}$ is correct
data; the bridge to a Galois representation $\rho_{\Delta_{10}} :
\mathrm{Gal}(\overline{\mathbb Q} / \mathbb Q) \to \mathrm{GSp}_4(\mathcal O_E)$
works on reduced rings (Chenevier 2014 Thm 2.12 equivalence
pseudo-characters $\leftrightarrow$ determinants on reduced rings).

**Precise error**: conflates the Taylor-Wiles pseudo-character
(older, weaker, multilinear symmetric trace functions) with the
Chenevier 2014 determinant (newer, stronger, single homogeneous
polynomial law with multiplicativity, unitality, Cayley-Hamilton).
On non-reduced rings — exactly the deformation rings that Vol II's
one-loop Quillen norm and BV 3d HT-QFT partition function feel — the
two objects differ; the determinant captures strictly more structure
(nilpotent Cayley-Hamilton witnesses that the pseudo-character
silently drops).

**Correct**: $D^{\mathrm{Chen}} : \mathbb T^{\mathrm{par}}_1 \to
\mathcal O_E \otimes \mathbb Z_\ell$ is a 4-dimensional Chenevier
determinant. Its graded components $(\Sigma_1, \Sigma_2, \Sigma_3,
\Sigma_4)$ at Hecke generators $T_p$ recover $(a_p(f_{18}) + p^8 +
p^9, \ldots, p^{34})$ via the reciprocal spinor $L$-factor expansion
$\prod_{i=1}^4 (1 - \alpha_i x) = 1 - \Sigma_1 x + \Sigma_2 x^2 -
\Sigma_3 x^3 + \Sigma_4 x^4$. Verified at 56 primes $p \le 263$.

**Vol II framing**: the Vol II arithmetic reading ties Vol II's
one-loop Quillen norm and BV 3d HT-QFT partition function to the
Galois-side $L$-function data via the Chenevier determinant; the
older "pseudo-character" framing was a proxy sufficient for W20.4's
trace identities on reduced rings but insufficient for Vol II's
scheme-valued path-integral refinements (deformation rings around
the Saito-Kurokawa lift are non-reduced, and the BV anomaly
measures the nilpotent Cayley-Hamilton witnesses that $S^{\mathrm{ps}}$
cannot see). The universal identity on the B-row
(three-faces, AP-V2-20) reads off from $D^{\mathrm{Chen}}$, not from
$S^{\mathrm{ps}}$.

**Primary**: Chenevier 2014 arXiv:1301.0635 Sec 1.2 Def/Prop 1.9,
Thm 2.12 (determinant $=$ pseudo-character on reduced rings; strict
inequality on non-reduced); Taylor 1991 Duke 63 Thm 2.1
(pseudo-character original definition); Ikeda 2001 Ann Math 154
Cor 16.2 (Saito-Kurokawa lift); Weissauer 2005 LNM 1868 (spinor
Galois representation); Laumon 2005 Publ IHES 102 (geometric
Satake on $\mathrm{GSp}_4$).

**Cross-ref**: Vol I AP902 / Remark
`rem:dl-w25-determinant-not-pseudocharacter` in
`chapters/theory/derived_langlands.tex` and cache entry 422; Vol III
cache entry at `appendices/first_principles_cache.md` row 8 (tagged);
Vol III `ADJUDICATION_LEDGER` §III.C; Pattern 295
(Creutzig-Ridout / Lyubashenko coend pseudo-traces on non-semisimple
MTCs — DISTINCT object, do NOT conflate: coend pseudo-traces live
on modular tensor categories, not on Hecke algebras, and satisfy a
different (Kerler-Lyubashenko modified-trace) axiom set).

### AP-V2-24 / V2-AP127: Humbert--Heegner admissibility filter $n\equiv 3,5\pmod 8$ on the pentagon coboundary tower $\phi^{(n)}$

**Trigger**: any Vol II assertion that $\phi^{(n)}\ne 0$ on the
K3--Humbert regime of the $\mathsf{SC}^{\mathrm{ch,top}}$ coloured bar
differential from the sole datum $d_n>0$ (Padovan motivic-MZV
dimension, Brown 2012 *Ann Math* 175 Thm 1), without installing the
Humbert--Heegner admissibility filter
$n\equiv 3,5\pmod 8$ and the Eichler--Zagier 1985 polar-support cutoff
$\Delta\ge -1$ on the paramodular index-1 K3 elliptic genus.

**Ghost**: the pentagon coboundary tower $\{\phi^{(n)}\}_{n\ge 3}$ of
Definition \texttt{def:phi-n-pent-EK} (Vol I
\texttt{chapters/theory/shadow\_tower\_higher\_coefficients.tex}) has a
well-defined three-filter admissibility structure on the K3
$A_\infty$-Humbert regime: (i) Padovan $d_n = d_{n-2}+d_{n-3}$ counts
the motivic-MZV transcendence basis at weight $n$ (real theorem, Brown
2012); (ii) Eichler--Zagier 1985 polar-support cutoff $C(\Delta)=0$ for
$\Delta < -m^2 = -1$ on the index-1 weak Jacobi form $\phi^{K3}_{0,1}$
is a real theorem (Eichler--Zagier *Prog Math* 55 Thm 9.3 with
$C(-1)=2$, $C(0)=20$); (iii) Gritsenko--Nikulin 1998 *J Reine Angew
Math* 507 paramodular lift of the K3 elliptic genus gives explicit
$c_{\Phi_{10}/\eta^{24}}$ Fourier data.

**Precise error**: bare Padovan-dimension $d_n$ count WITHOUT the
Humbert--Heegner admissibility filter overcounts. Most
Padovan-admissible $n\ge 3$ (all $n\ge 3$ except $n=4$) are
Humbert--Heegner-FORBIDDEN on the K3--Humbert regime: the paramodular
lattice sum $\sum_{4NM-\ell^2=-D_n} c_{\Phi_{10}/\eta^{24}}(N,\ell,M)$
with $D_n = (n-3)/2$ is non-empty iff $D_n\bmod 4\in\{0,1\}$, which
forces $n\equiv 3,5\pmod 8$. Asserting a non-zero $\phi^{(n)}$ on the
Vol II Swiss-cheese coloured-bar reading on the sole basis of $d_n>0$
(e.g., at $n=7,9,12,24,26,\dots$) silently conflates the
MZV-transcendence count with the paramodular Humbert--Heegner
signature and misses the Heegner--Bruinier obstruction class
$\mathrm{ob}^{\mathrm{HB}}_n\in H^2(H_n,\mathrm{Sym}^2 T^{\mathrm{poly}}_{\mathrm{ch}}|_{H_n})$
of Bruinier-torsion order $c_n$ (Bruinier 2002 LNM 1780 Chern class on
Heegner divisors).

**Correct**:
$\phi^{(n)}\big|_{\mathrm{K3\text{-}Humbert}}\ne 0$ iff (i)
$n\equiv 3,5\pmod 8$ AND (ii) the $d_n$-dimensional Brown canonical
basis is non-empty AND (iii) $D_n\le 1$ (Eichler--Zagier polar cutoff).
First non-vanishing cases: $\phi^{(3)}$ = Drinfeld pentagon cocycle
($D_3=0$, $C(0)=20\ne 0$); $\phi^{(5)} = -2\cdot[\mathrm{gen}]^{\otimes 5}$
with Gritsenko--Nikulin 1998 Table 2 sign on $\Phi_{10}/\eta^{24}$
($D_5=1$, $C(-1)=2\ne 0$). Humbert--Heegner admissible $n\in[3,36]$:
$\{3,5,11,13,19,21,27,29,35\}$. Padovan-positive HH-forbidden $n$
(e.g., $4,6,7,8,9,10,12,14,15,16,17,18,20,22,23,24,25,26,28,30,31,32,33,34,36$)
all give $\phi^{(n)}=0$ on K3--Humbert. HH-admissible $n\ge 11$ give
$\phi^{(n)}=0$ on K3--Humbert by Eichler--Zagier polar support
($D_n\ge 4 > 1$).

**Condensed reference table** $(n, d_n, D_n, \mathrm{HH}, \phi^{(n)}\text{-K3})$:
$(3,1,0,Y,\text{non-zero})$; $(4,0,1/2,-,0)$;
$(5,1,1,Y,-2[\mathrm{gen}]^{\otimes 5})$; $(6,1,3/2,-,0)$;
$(7,1,2,N,0)$; $(8,2,5/2,-,0)$; $(9,2,3,N,0)$; $(10,2,7/2,-,0)$;
$(11,3,4,Y,0\,\text{polar})$; $(12,4,9/2,-,0)$;
$(13,5,5,Y,0\,\text{polar})$; $(19,17,8,Y,0\,\text{polar})$;
$(21,28,9,Y,0\,\text{polar})$; $(27,90,12,Y,0\,\text{polar})$;
$(29,149,13,Y,0\,\text{polar})$; $(35,504,16,Y,0\,\text{polar})$.

**Vol II framing (Swiss-cheese coloured-bar reading)**: the
Humbert--Heegner admissibility filter is the Humbert-stratification
refinement of the $\mathsf{SC}^{\mathrm{ch,top}}$ coloured bar
differential on $\overline{\mathcal A_2}$. Each coloured face of the
$\mathsf{SC}^{\mathrm{ch,top}}$ operadic resolution at weight $n$
couples to the paramodular K3 elliptic-genus polar slice $\Delta\ge -1$;
a face colouring is admissible iff the discriminant $D_n=(n-3)/2$ lies
in $\{0,1\}\bmod 4$, equivalently $n\equiv 3,5\pmod 8$. Non-admissible
coloured faces carry the Heegner--Bruinier obstruction class
$\mathrm{ob}^{\mathrm{HB}}_n$ of Bruinier-torsion order $c_n$; the
Kontsevich formality per-wall obstruction of AP-V2-6 and the
Hodge-to-de-Rham $d_2$ differential of AP-V2-12 both read off from this
same Heegner--Bruinier class. The filter locks which of the eight
bar-cobar ambient colourings of AP-V2-1 is admissible at each weight
$n$. Every Vol II $\phi^{(n)}$ reference on the coloured bar
differential (all Wave-20-28 entries) carries the HH admissibility
scope or a pointer to Theorem \texttt{thm:phi-n-humbert-heegner-admissibility}
inscribed in Vol I \texttt{chapters/theory/shadow\_tower\_higher\_coefficients.tex}
(lines 4364-4433).

**Three verification paths for the filter**: (i) discriminant-form
signature — the index-$1$ paramodular form $4NM-\ell^2\equiv -\ell^2\pmod 4$
takes values in $\{0,-1\}\pmod 4$, so $-D_n$ is representable iff
$D_n\in\{0,1\}\pmod 4$, forcing $n\equiv 3,5\pmod 8$ by odd-$n$
integrality; (ii) Eichler--Zagier 1985 weak-Jacobi-form polar-support
cutoff ($C(\Delta)=0$ for $\Delta<-m^2=-1$, real theorem); (iii)
Gritsenko--Nikulin 1998 paramodular lift of the K3 elliptic genus with
explicit $c_{\Phi_{10}/\eta^{24}}$ Fourier table (real object).

**Primary**: Eichler--Zagier 1985 *Prog Math* 55 Thm 9.3
(polar-support cutoff); Gritsenko--Nikulin 1998 *J Reine Angew Math*
507 (Humbert--Heegner structure, paramodular $\Phi_{10}/\eta^{24}$ sign
convention Table 2); Bruinier 2002 LNM 1780 §5 (Chern class on Heegner
divisors, torsion orders $c_n$); Brown 2012 *Ann Math* 175 Thm 1
(Padovan motivic-MZV dimension).

**Cross-ref**: Vol I Theorem \texttt{thm:phi-n-humbert-heegner-admissibility}
in \texttt{/Users/raeez/chiral-bar-cobar/chapters/theory/shadow\_tower\_higher\_coefficients.tex}
(lines 4364-4433); Vol I cache row 304 (AP890) + Pattern 299
comprehensive; Vol I \texttt{notes/antipatterns\_catalogue.md}
AP903-HH; Vol III \texttt{notes/antipatterns\_catalogue.md} AP-CY142
(partner entry); Vol III \texttt{appendices/first\_principles\_cache.md}
tip-table row V16; Vol III
\texttt{notes/first\_principles\_cache\_comprehensive.md} AP-CY142 long
form. Cross-links Vol II AP-V2-1 (eight Swiss-cheese bar-cobar
ambients), AP-V2-6 (chiral Kontsevich formality on Humbert strata),
AP-V2-12 (Kaledin Hodge-to-de-Rham scope) — the HH filter locks the
admissible coloured-bar ambient at each weight $n$.

### AP-V2-25 / V2-AP128: Swiss-cheese colours framed without stage-1/stage-2 identification (Vol III two-stage factorisation, 2026-04-22)

**Trigger**: any Vol II statement of the two-coloured operad
$\mathsf{SC}^{\mathrm{ch,top}}$ in which the CLOSED colour
($\mathrm{FM}_k(\mathbb{C})$, $E_2$-holomorphic) and the OPEN colour
($\mathrm{Conf}_m(\mathbb{R})$, $E_1$-topological) are declared as two
independent colour indices without naming the Vol III two-stage
factorisation identifying each colour with a stage of
$\Phi_d = \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C} \circ \Phi^{\mathrm{FA}}_d$.

**Ghost**: Vol III's CY-to-chiral programme produces a canonical
$E_d$-homotopy factorisation algebra on the CY-$d$ target in
\emph{stage 1} $\Phi^{\mathrm{FA}}_d$ (Kontsevich--Tamarkin formality
plus Costello--Gwilliam--Li factorisation-homology BV quantisation),
then pushes forward via $\int_{\Sigma_{d-1}}$ to a reference curve $C$
in \emph{stage 2} $\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}$, landing
on an $E_1$-chiral shadow. The two colours of
$\mathsf{SC}^{\mathrm{ch,top}}$ are exactly the two stages: CLOSED =
stage-1 $\Phi^{\mathrm{FA}}_d$ on local observables of an
$E_d$-holomorphic factorisation algebra on a formal disc in the CY
target ($\mathrm{FM}_k(\mathbb{C})$); OPEN = stage-2
$\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}$ landing on the reference
curve ($\mathrm{Conf}_m(\mathbb{R})$); MIXED operations
(closed-input $\to$ open-output) = factorisation-homology pushforward
$\int_{\Sigma_{d-1}}$, with the directional restriction
$\mathsf{SC}^{\mathrm{ch,top}}(\ldots, \mathsf{top}, \ldots;
\mathsf{cl}) = \varnothing$ expressing that stage~2 is a specialisation
of stage~1, never an inversion.

**Precise error**: writing $\mathsf{SC}^{\mathrm{ch,top}}$ as a pure
combinatorial datum (closed colour, open colour, mixed operations),
decoupled from $\Phi_d$, loses (i) the specialisation direction from
Vol III stage~1 to Vol III stage~2, (ii) the identification of the
mixed-operation pushforward with $\int_{\Sigma_{d-1}}$, and (iii) the
directional-restriction asymmetry as a shadow of the fibre-integration
arrow. Framing the colours purely operadically hides the Vol III
anchor and leaves the Vol II $\leftrightarrow$ Vol III bridge
diagrammatic.

**Correction**: every $\mathsf{SC}^{\mathrm{ch,top}}$ statement that
carries cross-volume content names the stage identification in its
preamble or a nearby remark:
\[
\text{CLOSED} = \text{stage-1 } \Phi^{\mathrm{FA}}_d |_{\text{local obs}},
\quad
\text{OPEN} = \text{stage-2 } \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C},
\quad
\text{MIXED} = \int_{\Sigma_{d-1}}.
\]
Inscribed at `chapters/theory/sc_chtop_heptagon.tex` Remark
\texttt{rem:heptagon-two-stage-CY-to-chiral} and
`chapters/theory/factorization_swiss_cheese.tex` Remark
\texttt{rem:pentagon-two-stage}; cross-volume macros
$\PhiFA, \SpCh, \HolFA, \EdHolFA, \EnHolFA, \intSigma, \hCS$ declared in
`main.tex` preamble with `\providecommand` mirrors in affected
chapters.

**Severity**: Critical — load-bearing for the Vol II / Vol III bridge
and for the heptagon / pentagon coherence diagrams.

**Primary**: Kontsevich 2003 *Lett Math Phys* 66 (deformation
quantisation on $\mathbb{R}^d$); Tamarkin 2003 *Lett Math Phys* 66
(GRT-torsor structure); Costello--Gwilliam 2017 / 2021 Vols 1--2 (BV
factorisation quantisation); Ayala--Francis 2015 *J Topol* 8
($\mathrm{Conf}_m(\mathbb{R})$ and $E_1$-factorisation homology);
Voronov 1998 *Math Ann* 312 (classical Swiss-cheese); Francis--Gaitsgory
2012 *Selecta Math* 18 (Ran-space factorisation stratification).

**Cross-ref**: Vol II CLAUDE.md "Two-stage factorisation: Vol III
alignment" section (lines 281--318); Vol III
\texttt{notes/antipatterns\_catalogue.md} AP-CY two-stage entry
(2026-04-22 wave); AP-V2-1 (eight bar-cobar ambient colourings refined
by the stage identification); AP-V2-9 (Theorem H hinge crosses stages
at $d \ge 2$); AP-V2-14 (CY-4 HK-restricted $\Phi_4$ inherits the
two-stage structure on the HK locus).

### AP-V2-26 / V2-AP129: Single-stage $\Phi_d$ framing of Vol III output (2026-04-22)

**Trigger**: any Vol II reference to Vol III's CY-to-chiral functor
that writes the output as a unified $\Phi_d$ without naming it as the
stage-2 shadow
$\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C} \circ \Phi^{\mathrm{FA}}_d$.

**Ghost**: Vol III's $\Phi_d$ factors canonically as stage-1
$\Phi^{\mathrm{FA}}_d$ (the $E_d$-holomorphic factorisation algebra)
followed by stage-2 $\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}$
(specialisation to a reference curve $C$ via factorisation-homology
pushforward along $\Sigma_{d-1}$). A CY-$d$ category admits a
\emph{family} of $E_1$-chiral shadows
$\{\Phi_d(\mathcal C; \Sigma_{d-1}, C)\}_{(\Sigma_{d-1}, C)}$
parametrised by the polarisation of the CY target.

**Precise error**: treating $\Phi_d$ as a unified functor forgets the
$(\Sigma_{d-1}, C)$-parameter family and the stage-2 reference-curve
choice; the stage-2 pushforward is not canonical (it depends on the
polarisation of $\text{CY}_d$ locally). Vol II cross-references
invoking $\Phi_d$ must carry the stage scope, or default to stage-1
$\Phi^{\mathrm{FA}}_d$ when the statement lives on the CY side.

**Correction**: Vol II cross-references to Vol III output use one of
three precise forms:
$\Phi_d = \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C} \circ
\Phi^{\mathrm{FA}}_d$ (both stages named), stage-1
$\Phi^{\mathrm{FA}}_d$ (CY-side claim), or stage-2
$\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}$ (curve-side claim). The
six primary Vol II files touching Vol III output have been aligned to
this convention; new inscriptions follow it from the first keystroke.

**Severity**: Critical — every cross-volume claim on
$\Phi_d(\mathcal C)$ carries the stage scope.

**Primary**: Kontsevich--Tamarkin formality (stage-1 existence);
Costello--Gaiotto--Li BV quantisation (stage-1 $E_d$ structure);
Ayala--Francis factorisation-homology pushforward (stage-2
specialisation).

**Cross-ref**: Vol II CLAUDE.md two-stage section (lines 281--318);
Vol III `chapters/theory/cy_to_chiral.tex` (stage-1 / stage-2
canonical construction); Vol III AP-CY two-stage entries (2026-04-22);
related to AP-V2-25 (closed/open colour identification); AP-V2-9
(Theorem H hinge at $d \ge 2$).

### AP-V2-27 / V2-AP130: 3D HT QFT scope without two-stage anchor (2026-04-22)

**Trigger**: statements about the Vol II 3D holomorphic-topological
QFT framework (bulk $\mathsf{SC}^{\mathrm{ch,top}}$ pairing, derived
centre $Z^{\mathrm{der}}_{\mathrm{ch}}(A)$ as bulk algebra) that omit
the cross-volume anchor identifying the $\mathsf{SC}^{\mathrm{ch,top}}$
structure with the Vol III two-stage factorisation at $d = 2, 3$.

**Ghost**: the Vol II 3D HT QFT framework is intrinsically a
Swiss-cheese / factorisation-operadic construction; its natural
ambient at the CY target is the stage-1 $\Phi^{\mathrm{FA}}_d$ output.
$\Phi^{\mathrm{FA}}_3$ exists canonically by Kontsevich--Tamarkin
formality plus Costello--Gwilliam--Li factorisation-homology BV
quantisation; stage-2 $\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_2, C}$ lands
the Vol II data on a reference curve.

**Precise error**: presenting the 3D HT QFT framework as
self-contained on the curve $C$, without flagging that the CLOSED
colour of $\mathsf{SC}^{\mathrm{ch,top}}$ is the stage-1 CY-side datum
and the OPEN colour is the stage-2 curve-side datum, hides (i) the
$(\Sigma_{d-1}, C)$ parameter, (ii) the $\int_{\Sigma_{d-1}}$-
pushforward as the MIXED operation, and (iii) the existence of
$\Phi^{\mathrm{FA}}_3$ as a named object. The 3D HT QFT content is
cross-volume: it lives at the Swiss-cheese layer = the two-stage
interface.

**Correction**: Vol II 3D HT QFT statements in reader-facing prose
name the two-stage anchor at first mention per chapter, via a
parenthetical or a remark: "the $\mathsf{SC}^{\mathrm{ch,top}}$
pairing on $C \times \mathbb{R}$ is the stage-2 shadow of the Vol III
$E_3$-holomorphic factorisation algebra $\Phi^{\mathrm{FA}}_3$, with
$\int_{\Sigma_2}$ the MIXED factorisation-homology pushforward". The
existence of $\Phi^{\mathrm{FA}}_3$ as a canonical stage-1 object is
itself a Vol III theorem (Kontsevich--Tamarkin + CGL); Vol II does not
re-prove it, only invokes it at the bridge.

**Severity**: High — Vol II's 3D HT QFT framework is the main
contribution of Parts V--VI; the two-stage anchor is its natural place
in the cross-volume architecture.

**Primary**: Costello--Gwilliam 2017 / 2021 Vols 1--2 (BV quantisation
on factorisation algebras); Kontsevich 2003 (deformation quantisation);
Costello--Gaiotto 2019 (5d / 3d HT twists); Costello--Li 2023
(holomorphic Chern--Simons factorisation); Ayala--Francis 2015
(factorisation-homology pushforward).

**Cross-ref**: Vol II CLAUDE.md two-stage section (lines 281--318);
AP-V2-25 (closed / open colour = stage-1 / stage-2); AP-V2-26
(single-stage $\Phi_d$ framing); AP-V2-1 (bar-cobar ambient
octachotomy refined by stage-1 vs stage-2).

### AP-V2-28 / V2-AP131: Manifesto conflations inherited from Vol III (2026-04-22)

**Trigger**: any Vol II reference to a Vol III canonical value that
re-states an earlier-retracted reading, in particular:
$\kappa_{\mathrm{cat}}(K3 \times E) = 2$;
$\mathrm{CoHA}(\mathbb{C}^3) = \mathcal{W}_{1+\infty}$; "six routes to
$G(K3 \times E)$ are six $\Phi_3$-applications"; Fake-Monster placed at
$d = 3$.

**Ghost**: four canonical readings of cross-programme objects bind at
every Vol II cross-reference site:
- $\kappa_{\mathrm{cat}}(K3 \times E) = \chi(\mathcal O_{K3})
  \cdot \chi(\mathcal O_E) = 2 \cdot 0 = 0$ on the total space
  (Künneth-multiplicative);
- $\mathrm{CoHA}(\mathbb{C}^3) = Y^+(\widehat{\mathfrak{gl}}_1)$,
  the \emph{positive half} of the affine Yangian (not the full
  $\mathcal{W}_{1+\infty}$ algebra, not $\mathcal{W}_\infty[c]$);
- the six routes to $G(K3 \times E)$ (CoHA, Schiffmann--Vasserot,
  Maulik--Okounkov, Borcherds, Toda, DMVV) are six \emph{different
  constructions}, each taking a different CY-input category, each
  producing the same $\Phi_3$-output through a pentagon colimit; NOT
  six applications of $\Phi_3$;
- the Fake-Monster BKM sits at $d = 5$ (signature
  $\mathrm{II}_{25,1}$, Borcherds 1992 *Invent Math* 109), NOT at
  $d = 3$ where the K3-BKM $\mathbf H_{\Delta_5}$ lives (signature
  $\Lambda^{2,1}_{II}$, Gritsenko--Nikulin 1998).

**Precise error**: Vol II cross-references inheriting fibre
$\kappa_{\mathrm{cat}}(K3) = 2$ at the total-space site (naive
Künneth-additive confusion); conflating $\mathrm{CoHA}(\mathbb{C}^3)$
with $\mathcal{W}_{1+\infty}$ (classical-limit vs full positive half);
or placing Fake-Monster on the $d = 3$ row that belongs to K3-BKM.
These propagate via status tables, landscape maps, and remark
cross-references; each re-statement is a latent canonical-preamble
violation (rows 20--22, 32, 59, 61--62, 65--66).

**Correction**: every Vol II cross-reference to a Vol III canonical
value cites as named in the preamble (lines 5--67 of this catalogue)
or names the preamble explicitly. Concretely:
- $\kappa_{\mathrm{cat}}(K3 \times E) = 0$ always (name the total
  space); the fibre value $2$ is $\kappa_{\mathrm{fiber}}(K3)$, a
  different $\kappa_\bullet$ index;
- $\mathrm{CoHA}(\mathbb{C}^3) = Y^+(\widehat{\mathfrak{gl}}_1)$
  (positive half; $\mathcal{W}_{1+\infty}$ is a classical-limit
  current subalgebra);
- six routes to $G(K3 \times E)$ are six distinct constructions
  producing the same $\Phi_3$ output;
- Fake-Monster at $d = 5$ with $\Phi_{12}$; K3-BKM at $d = 3$ with
  $\Phi_{10} = \Delta_5^2$.

**Severity**: High — each conflation has re-appeared across waves
(Wave 16 Fake-Monster rank, Wave 23 $\kappa_{\mathrm{cat}}$ drift);
the canonical preamble is the binding register.

**Primary**: Schiffmann--Vasserot 2012 *Publ IHES* 115 (CoHA);
Borcherds 1992 *Invent Math* 109 (Fake-Monster); Gritsenko--Nikulin
1998 (K3-BKM); Mukai 1984 (K3 Mukai lattice).

**Cross-ref**: Canonical preamble rows 20--22, 32, 59--63; AP-V2-22
(five-archetype landscape); AP-V2-20 (three-faces universal
identity); V2-AP68 / V2-AP69 (CoHA / $\kappa_\bullet$ discipline);
V2-AP122 (six routes = six constructions); Vol III
\texttt{chapters/examples/cy\_d\_kappa\_stratification.tex} (canonical
$\kappa_\bullet$ table); Vol III CLAUDE.md key-facts block.

### AP-V2-29 / V2-AP132: Structural voice discipline — `\begin{warning}` and bookkeeping vocabulary in reader-facing Vol II prose (2026-04-22)

**Trigger**: any reader-facing Vol II `.tex` under `chapters/`,
`frame/`, `examples/`, `theory/`, `connections/`, or `bibliography/`
containing (i) a `\begin{warning}` environment, or (ii) bookkeeping
tokens such as `AP-V2-N`, `V2-APN`, `Pattern N`, `Wave N`,
`cache-canonical`, `CLAUDE.md`, `manifesto`, `Beilinson-standard`,
`HZ-N`, `DNA strand`, `cache entry`, `CG-rectify pass`, `round K`, or
`batch M`.

**Ghost**: the Vol II manuscript is self-complete, self-coherent,
self-consistent. Reader-facing prose is written in the
Chriss--Ginzburg voice. Bookkeeping vocabulary — waves, rounds, AP
indices, cache entries, rectification passes, CLAUDE.md pointers —
belongs in `notes/`, `FRONTIER.md`, commit messages, and the local
`memory/`. `\begin{warning}` is a bookkeeping environment (the reader
sees a "warning" rather than a mathematical statement) and is
reserved for `notes/` scaffolding if anywhere.

**Precise error**: V2-AP51 earlier recorded that reader-facing prose
containing `Wave N` / `AP\d+` / `DNA` tags triggered CG-rectify scrubs
that removed substantive mathematical inscriptions as collateral
(cache entry 89, 2026-04-20/21). `\begin{warning}` in reader-facing
`.tex` is a parallel failure mode: it presents bookkeeping to the
reader as if it were mathematics.

**Correction**: reader-facing inscriptions contain (a) no
`\begin{warning}` environments, (b) no bookkeeping tokens listed in
the Trigger. Named section / subsection / remark titles denote
mathematical objects; equations carry mathematical labels. Hedging
that the mathematics does not earn is forbidden; the equals sign is a
theorem. Retroactive rectification runs per-chapter through the
`chriss-ginzburg-rectify` skill; forward inscriptions follow the CG
voice from the first keystroke. The "Session antipatterns — manuscript
hygiene (2026-04-22)" block below (CGCLEAN-1..55) enumerates the
specific vocabulary patterns with canonical repair.

**Severity**: High — the scrub-cascade mechanism is real (V2-AP51 +
cache entry 89); the forward fix is cheap; per-chapter sweep when
introducing new content.

**Primary**: Vol II CLAUDE.md "Writing standard: Chriss--Ginzburg
north star" section (lines 128--192); V2-AP51 (scrub-cascade
causality); cache entry 89 (hook-cascade content loss).

**Cross-ref**: Vol II CLAUDE.md lines 128--192; V2-AP51 (legacy
scrub-cascade AP); the "Session antipatterns — manuscript hygiene
(2026-04-22)" block below (CGCLEAN-1..55) is the detailed sibling
register (this AP-V2-29 is the upstream reader-facing voice-discipline
AP; CGCLEAN-1..55 are its per-vocabulary repair rows); Vol III
\texttt{notes/antipatterns\_catalogue.md} voice-discipline entry
(2026-04-22 wave); Vol I
\texttt{notes/antipatterns\_catalogue.md} corresponding voice-discipline
AP (parallel constraint across all three volumes).

## Session antipatterns — manuscript hygiene (2026-04-22)

**Standing principle** (from CLAUDE.md of all three volumes). The manuscript is self-complete, self-coherent, self-consistent. The current version stands for itself and only itself. All reader-facing LaTeX is standalone and up-to-date. No references to previous versions, intermediate ansätze, earlier drafts, retracted values, superseded formulas, or drafting-history commentary. When a mathematical retraction is genuinely informative, state the failed argument and its Gap/Flaw as mathematics, not as drafting record. Bookkeeping, meta-narration, and version-history apparatus belong in `notes/`, `FRONTIER.md`, commit messages, and `memory/` — never in the manuscript.

Conservative scope reminder. `cascade`, `homological retraction`, `deformation retract`, `retracts onto` are legitimate mathematical terms. `Borcherds programme` is legitimate (proper noun). `Platonic solid` is legitimate. Forbidden forms: `Platonic Theorem~A`, `our programme`, `programme-canonical`, `platonic chapter`, `\begin{warning}`, `\ClaimStatusRetracted`, `narrative counterpart`, `History of the claim`, drafting-dated remarks.

### Group A — Bookkeeping vocabulary

| Code | Pattern | Why it fails | Canonical repair |
|---|---|---|---|
| CGCLEAN-1 | `Wave $N$` / `Wave-$N$ / W$N$` in manuscript prose | Session-tracking label; the manuscript is timeless | Delete the label; if the content is a real mathematical result, state the theorem directly with its canonical name |
| CGCLEAN-2 | `AP-CY$n$` / `AP$n$` / `AP-CAT-$N$` in manuscript prose | Antipattern registry index; belongs in notes/, not LaTeX | Delete the bracketed code; if the underlying mathematical hazard is load-bearing, inscribe the hazard as a Remark or scope statement |
| CGCLEAN-3 | `FM$n$` (formula-mechanical marker) in manuscript | Hook / cache index; not a mathematical object | Delete; replace with the named mathematical content |
| CGCLEAN-4 | `HZ-$N$` / `HZ-IV` / `HZ discipline` in manuscript | Internal hierarchy label | Delete; replace with the actual hierarchy (primary source > direct computation > build system, stated in-line if relevant) |
| CGCLEAN-5 | `DNA strand S$x$` / `S-strand` in manuscript | Swarm-workflow metaphor | Delete; inscribe the mathematical content without the workflow layer |
| CGCLEAN-6 | `CG-rectify pass $k$` / `rectification pass $N$` in manuscript | Drafting workflow artefact | Delete; if rectification changed the content, the current text already reflects it |
| CGCLEAN-7 | `cache entry $n$` / `Cached Confusion` / `Cache anchor` / `Cache append` in manuscript | `notes/` apparatus leak | Delete; cite the primary mathematical reference instead |
| CGCLEAN-8 | `Wave $N$ spec` / `Wave $N$ verdict` / `Wave $N$ witnessing` in manuscript | Session-audit meta-language | Delete; the verdict is the mathematical statement that survives |
| CGCLEAN-9 | `programme-canonical` / `programme canonical value` | Meta-registry term; the manuscript has canonical values, not "programme-canonical" values | Replace with "canonical" or delete the adjective; the value stands for itself |
| CGCLEAN-10 | `type-error registry entry T$n$` / `type-error T$n$` | Registry index | Delete; if a genuine type error is the point, describe it as a mathematical mismatch (e.g., "CoHA is $E_1$-associative; the bar complex is a coalgebra") |

### Group B — Meta-narration

| Code | Pattern | Why it fails | Canonical repair |
|---|---|---|---|
| CGCLEAN-11 | `narrative counterpart` / `narrative arc` | Signals prose-as-explanation rather than mathematics-as-construction | Delete; replace with a direct mathematical statement of the construction |
| CGCLEAN-12 | `story` / `saga` / `odyssey` / `journey` (as nouns in manuscript prose) | Narrates author's intent instead of stating the theorem | Delete the noun; state the result |
| CGCLEAN-13 | `Platonic ideal` / `Platonic form` / `platonic chapter` / `platonic architecture` / `Platonic ensemble` / `platonic synthesis` in manuscript | "Platonic" is a workflow/aspirational tag, not a mathematical adjective | Delete; the chapter either inscribes the structure or does not |
| CGCLEAN-14 | `Platonic Theorem~A` / `Platonic Theorem~B` etc. | Theorems are named, not platonicised | Replace with "Theorem~A" (with the canonical reference) |
| CGCLEAN-15 | `This chapter's function is to...` / `The function of this section is...` | Meta-narration of authorial intent | Delete; the chapter performs its function by being read |
| CGCLEAN-16 | `we now turn to` / `having established` / `let us now` / `this brings us to` / `with this in hand` | Signpost phrases | Replace with a direct mathematical transition or a single-sentence summary of what is proved next |
| CGCLEAN-17 | `in the present work` / `the author` / `our programme` / `we have argued` / `it is worth noting` | Self-referential authorial voice | Delete; the work is the manuscript, not a framing around it |
| CGCLEAN-18 | Meta-paragraphs like `This chapter closes the ...` | Closing narration | Delete; if the section boundary earns three sentences, those sentences should be mathematics (what was established, what is the next obstruction, what resolves it) |
| CGCLEAN-19 | `the opening paragraphs of this preface` / `as noted in the preface` | Cross-reference within own manuscript via drafting geography | Replace with a theorem / section label cross-reference (`\ref{thm:...}` / `\S\ref{sec:...}`) |
| CGCLEAN-20 | `Earlier in the volume` / `Later in the volume` | Drafting-order navigation | Replace with the actual label or delete if redundant |

### Group C — Version-history

| Code | Pattern | Why it fails | Canonical repair |
|---|---|---|---|
| CGCLEAN-21 | `retracted` / `retraction` / `now retracted` / `the retracted $X$` in manuscript | Drafting history; the manuscript is the current version | Delete the retracted claim and its retraction both; only the live claim remains |
| CGCLEAN-22 | `superseded` / `supersedes` / `superseded by` in manuscript | Drafting history | Delete; cite the correct current result as the single statement |
| CGCLEAN-23 | `earlier draft` / `previous version` / `intermediate ansatz` / `prior derivation` | Drafting history | Delete; state only the current derivation |
| CGCLEAN-24 | `previously conjectural` / `previously open` / `previously unresolved` / `previously obstructing` | Temporal status markers | Delete "previously"; state the current status directly |
| CGCLEAN-25 | `now resolved` / `now proved` / `now known` | Temporal status marker | Delete the temporal qualifier; state the result |
| CGCLEAN-26 | `double-retraction` / `triple-retraction` | Drafting history | Delete; state the mathematics that survives |
| CGCLEAN-27 | `Three successive evaluations appear in the drafting record` / `History of the claim` | Explicit drafting-record exposition | Delete entirely; if the failed arguments are informative, state each as a Gap/Flaw lemma with its mathematical obstruction, not as "evaluation $k$" |
| CGCLEAN-28 | `drafting record` / `drafting trajectory` / `drafting history` | Self-narration of the process | Delete |
| CGCLEAN-29 | `\ClaimStatusRetracted` tag / environment | Retraction status embedded in tag system | Remove the tag; if the claim is wrong, delete it; if the failure is load-bearing, inscribe a named Gap/Flaw lemma |
| CGCLEAN-30 | Dated remarks (`2026-04-17`, `Wave 14 verdict`, etc.) in manuscript prose | Drafting timestamps | Delete the date; the statement is timeless |
| CGCLEAN-31 | `\index{retraction!...}` entries | Retraction-indexed navigation | Remove; indexing should be mathematical (`\index{Borcherds product}`, not `\index{retraction!$c_3$}`) |

### Group D — Bookkeeping-content references

| Code | Pattern | Why it fails | Canonical repair |
|---|---|---|---|
| CGCLEAN-32 | `\texttt{notes/...}` as a reader-facing reference | `notes/` is not reader-facing apparatus | Either (a) promote the content to the manuscript if load-bearing, or (b) delete the reference |
| CGCLEAN-33 | Absolute home paths `/Users/raeez/...` in manuscript | Personal filesystem path | Replace with a repo-relative label, a `\ref`, or delete |
| CGCLEAN-34 | `% TODO: librarian verification` / author-note comments | Editorial comment | Delete; either verify and inscribe, or remove the claim |
| CGCLEAN-35 | `% ALIAS` / `% LEGACY ALIAS` comments | Label-rename bookkeeping | Delete the alias comment; if a label rename is needed, change the label and run references |
| CGCLEAN-36 | `% Source: NEW CHAPTER` / `% Source: Wave N` markers | Provenance comments | Delete |
| CGCLEAN-37 | Compute-engine filenames like `*_waveN_*.py` cited in manuscript | Session-indexed filename leaks | Rename the engine to drop `waveN`; cite by canonical descriptive name |
| CGCLEAN-38 | Function names `waveN_foo` cited in manuscript | Session-indexed function leak | Rename; cite the canonical mathematical functionality |

### Group E — Warning / hedge environments

| Code | Pattern | Why it fails | Canonical repair |
|---|---|---|---|
| CGCLEAN-39 | `\begin{warning} ... \end{warning}` environment | Warnings hedge the mathematics; the manuscript either proves or conjectures, never warns | Replace with a `\begin{remark}` that states the mathematical scope cleanly, or delete if the scope is already carried by the theorem hypotheses |
| CGCLEAN-40 | `do not confuse` / `don't be fooled` / `beware` / `be careful` | Hedging / scolding the reader | Delete; the reader confuses nothing if the mathematics is stated with sufficient precision |
| CGCLEAN-41 | `we must be careful` | Self-hedging | Delete; state the scope, do not announce caution |
| CGCLEAN-42 | gratuitous `scope-restricted` / `scope-bounded` as adjectives | Meta-scope commentary | Delete; state the actual scope ("on the Koszul locus", "for affine $A$", etc.) |
| CGCLEAN-43 | `verdict` as a label for a mathematical statement | Meta-adjudication language | Replace with `theorem`, `proposition`, `lemma`, `corollary`, or `remark` as appropriate |

### Group F — Structural meta-naming

| Code | Pattern | Why it fails | Canonical repair |
|---|---|---|---|
| CGCLEAN-44 | Chapter filename `*_platonic.tex` | File-level meta-tag | Rename to the canonical mathematical content; update `\input` in `main.tex` |
| CGCLEAN-45 | Chapter label `ch:*-platonic` | Label-level meta-tag | Rename to `ch:*` matching the content; update all `\ref` |
| CGCLEAN-46 | Section label `sec:*-platonic` | Label-level meta-tag | Rename to `sec:*`; update all `\ref` |
| CGCLEAN-47 | Theorem label `thm:*-waveN-*` | Label carries session index | Rename to the canonical mathematical content; update `\ref` |
| CGCLEAN-48 | `\index{compute module!...}` entries | Index entry for internal compute apparatus | Remove; the manuscript indexes mathematics |
| CGCLEAN-49 | `\index{cache!...}` entries | Index entry for the cache registry | Remove |
| CGCLEAN-50 | `\index{retraction!...}` entries | Retraction navigation | Remove |

### Group G — Edge cases

| Code | Pattern | Why it fails | Canonical repair |
|---|---|---|---|
| CGCLEAN-51 | `Five attack-heal calibrations` | Swarm-workflow count | Delete; if five theorems or five checks are the substance, list the five objects, not the five calibrations |
| CGCLEAN-52 | `Reconstitution if the cancellation fails` | Workflow-contingent phrasing | Delete; state either the cancellation theorem or the obstruction lemma as live mathematics |
| CGCLEAN-53 | `Inversion of the programme perspective` | Meta-framing | Delete; state the dual construction directly |
| CGCLEAN-54 | `Gold-standard HZ-IV disjoint verification` | Epistemic-hierarchy-label meta-claim | Delete; if the verification is multi-path, list the paths (direct computation, primary source, build test) as three concrete checks |
| CGCLEAN-55 | `Three successive evaluations appear in the drafting record` (verbatim) | Direct drafting-record callout | Delete; if a Gap/Flaw triptych is mathematically useful, inscribe as three named lemmas with their precise obstructions |

**Cross-references.** These 55 entries appear in parallel across the three volumes. Vol I: `notes/antipatterns_catalogue.md` session entry (2026-04-22) with identical codes CGCLEAN-1..55; Vol III: same. Cache entries in `notes/first_principles_cache_comprehensive.md` (this volume) register the regex + 5-step DETECT/LOCALISE/MATH-CHECK/REPAIR/VERIFY protocol. Reader-oriented summary rows in `notes/first_principles_cache.md`. Hook signatures in `.claude/hooks/beilinson-gate.sh` (violation codes CGCLEAN-1..55).

## Section 14: Cross-programme imports from Vol III wave-12 synthesis (2026-04-22)

The following entries propagate Vol III wave-12 crystallisations not
already covered by AP-V2-25--29 (Vol III two-stage factorisation
family, manifesto conflations, structural voice discipline). Each
entry names the Vol III wave-12 source file and pins the precise scope
at which the underlying mathematics binds Vol II content. Numbering
continues from AP-V2-29 without gaps.

### AP-V2-30 / V2-AP133: Pattern 273 — chain-level statement vs $(\infty,1)$-categorical functor claim on $(B_{E_1}, B_{E_2}, B_{E_\infty})$
**Severity**: Critical
**Trigger**: Vol II's three-sector bar complex
$(B_{E_1}, B_{E_2}, B_{E_\infty})$ stated as an equivalence of
$(\infty,1)$-categories without a proof of morphism preservation; or
any wording that conflates the chain-level filtered quasi-isomorphism
$B_{E_1}A \twoheadrightarrow B_{E_2}A \twoheadrightarrow B_{E_\infty}A$
with an ambient $(\infty,1)$-categorical-functor claim.
**Ghost**: the chain-level filtered equivalence is a genuine theorem
at the explicit complex level, carrying the Kontsevich--Tamarkin
$E_d$-formality inputs and the Costello--Li holomorphic-locality
inputs (AP-V2-33). The $(\infty,1)$-categorical-functor statement is
a distinct claim on morphism spaces that has not been established in
full generality; it is the functorial lift of the chain-level
filtered equivalence to an $(\infty,1)$-categorical map between bar
constructions of chiral-topological algebras.
**Precise error**: Pattern 273 scope conflation (chain-level vs
$(\infty,1)$-categorical functor). Neither lane subsumes the other;
the chain-level statement supports explicit factorisation-algebra
readings of $\mathsf{SC}^{\mathrm{ch,top}}$, the
$(\infty,1)$-categorical claim would support ambient functorial
pushforwards. Merging the two loses precision at the exact site where
Beilinson's dictum operates (chain-level proofs are auditable against
primary sources; $(\infty,1)$-categorical claims are auditable against
Lurie *Higher Algebra* and Francis 2013 in restricted scopes).
**Correction**: state the filtered equivalence as a chain-level
theorem with explicit inputs (Kontsevich--Tamarkin $E_d$-formality,
Drinfeld associator choice, Costello--Li holomorphic locality). Tag
any $(\infty,1)$-functor lift as conjectural and named
(AP-V2-36). Both readings are load-bearing at their precise scopes;
neither is a shadow of the other.
**Primary**: Kontsevich 2003 *Lett Math Phys* 66 (formality); Tamarkin
2003 *Lett Math Phys* 66 (GRT and formality); Costello--Li 2016
arXiv:1605.09439 (holomorphic locality); Francis 2013 *Compositio* 149
(factorisation homology).
**Inscription**: `chapters/theory/bar_chain_models.tex`; Vol III
wave-12 source `notes/wave12_a1_bar_cobar_gelfand.tex`; Vol III
`notes/wave12_a1_phi_functor_foundations.tex` (chain-level vs
functorial scope in the $\Phi$-construction).
**Cross-ref**: Vol III AP-CY Pattern 273; Vol I Pattern 273 cross-ref
hook; Vol II AP-V2-1 (eight bar-cobar ambients) sets the chain-level
home for this filtered equivalence; AP-V2-36 (frontier conjectural
lift).

### AP-V2-31 / V2-AP134: Bare $\kappa$ without subscript in chiral-topological statements
**Severity**: Critical
**Trigger**: any Vol II inscription containing $\kappa$ without a
$\mathrm{ch}/\mathrm{cat}/\mathrm{BKM}/\mathrm{fiber}$ subscript, in
prose, theorem statement, equation, or label.
**Ghost**: four genuinely distinct invariants track four different
mathematical structures on a CY-to-chiral datum. $\kappa_{\mathrm{ch}}$
is the chiral-side Hodge supertrace via $\Phi$.
$\kappa_{\mathrm{cat}} = \chi(\mathcal O_X)$ is
Künneth-multiplicative on products; $\kappa_{\mathrm{cat}}(K3 \times E)
= 0$ on the total space, NOT $2$ (which is $\kappa_{\mathrm{fiber}}(K3)$,
a different object on a different domain). $\kappa_{\mathrm{BKM}}(\Phi_N)
= c_N(0)/2$ is the Borcherds weight; at $N=1$ this is $5$ for the
paramodular $\Delta_5$ convention or $12$ for the Fake-Monster
$\Phi_{12}$ convention (AP5 dual-indexing).
$\kappa_{\mathrm{fiber}}$ is the fibre/lattice correction. Each
satisfies a different functoriality and enters a different face of
$r_{\mathrm{CY}}$.
**Precise error**: bare $\kappa$ collapses four distinct objects into
one symbol. On the $\mathsf{SC}^{\mathrm{ch,top}}$ Swiss-cheese
reading of $\mathbf H_{\Delta_5}$, substituting one for the other
silently forces wrong Borcherds-product denominators. The three-faces
identity AP-V2-20 requires subscripted $\kappa$ to triangulate
correctly across Vol I / Vol II / Vol III.
**Correction**: every $\kappa$ in Vol II carries its subscript from
the first keystroke. At cross-volume linking points (Theorem C
five-archetype table, three-faces universal identity, Borcherds-weight
anchor), pin the subscript to the site it enters through. AP-V2-22
B-row coupling $\kappa + \kappa^! = 8$ is precise only under the
subscripted reading.
**Primary**: Vol III CLAUDE.md and
`chapters/examples/cy_d_kappa_stratification.tex`; Borcherds 1995
*Invent Math* 120; Gritsenko 1999 *Algebra i Analiz* 11.
**Inscription**: all Vol II chapters touching the CY-to-chiral face;
Vol III wave-12 source
`notes/wave12_a2_kappa_invariants_universal_borcherds.tex`.
**Cross-ref**: Vol III AP-CY bare-$\kappa$ family and HZ-7 subscript
discipline; Vol I bare-$\kappa$ AP hook; Vol II
`first_principles_cache.md` canonical preamble rows 59--62 and 65--66
for the four $\kappa_\bullet$ values; AP-V2-28 (manifesto conflations
inherited from Vol III) names the $\kappa_{\mathrm{cat}}(K3 \times E)
= 2$ retraction.

### AP-V2-32 / V2-AP135: $\mathrm{CoHA}(\mathbb C^3) = Y^+$ as evaluation image, not $\mathcal W_{1+\infty}$ isomorphism
**Severity**: Critical
**Trigger**: Vol II inscription identifying the
cohomological-Hall-algebra image of $\mathbb C^3$ with
$\mathcal W_{1+\infty}$, or implying the Schiffmann--Vasserot arrow
is an isomorphism rather than a map to the positive half of the
affine Yangian.
**Ghost**: Schiffmann--Vasserot 2012 / 2020 identify
$\mathrm{CoHA}(\mathbb C^3) = Y^+(\widehat{\mathfrak{gl}}_1)$, the
positive half of the affine Yangian. The image in
$\mathcal W_{1+\infty}$ exists via an evaluation arrow that factors
through a further quotient and is only an isomorphism at the classical
/ $\hbar = 0$ limit. Miki $S_3$-triality on $Y^+$ pins the
$\mathcal W_{1+\infty}$ slice; the curved-Dunn $H^2 = 0$ closure
question (AP-V2-34) is the extension-ladder problem for pulling
$\mathcal W_{1+\infty}$ back to the chiral-topological Vol II site.
**Precise error**: treating the evaluation arrow $Y^+ \to
\mathcal W_{1+\infty}$ as an isomorphism; missing the Miki triality
structure on $Y^+$ (Miki 2007 generators); collapsing the extension
ladder between $Y^+$, $\mathcal W_\infty[c]$, and
$\mathcal W_{1+\infty}$ into a single line.
**Correction**: state $\mathrm{CoHA}(\mathbb C^3) = Y^+$;
$\mathcal W_{1+\infty} = $ image of the evaluation slice under
Miki-triality; $\mathcal W_{1+\infty} = \mathcal W_\infty[c] \otimes
\mathcal H$ as a tensor with the Heisenberg algebra (preamble row 63).
Curved-Dunn closure is a separate chiral-topological step with
$H^2$-obstruction living in the Costello--Gwilliam BV cohomology of
the Vol II factorisation assembly.
**Primary**: Schiffmann--Vasserot 2012 *Publ IHES* 115; Schiffmann--
Vasserot 2020 *Duke Math J* 169; Pope--Romans--Shen 1990 *Phys Lett B*
236; Miki 2007 *J Math Phys* 48; Costello--Gwilliam 2021 Vol 2 (BV
curved Dunn).
**Inscription**: `chapters/theory/coha_chiralisation.tex`; Vol III
wave-12 source `notes/wave12_a11_coha_y_plus_vs_w_infty.tex`.
**Cross-ref**: Vol II AP-V2-18 (MO-extension super-quasi-Hopf) via
the $Y^+$ side; Vol III AP-CY CoHA-vs-vertex-algebra family;
`first_principles_cache.md` canonical preamble rows 62--63; AP-V2-28
(manifesto conflations inherited from Vol III, CoHA row).

### AP-V2-33 / V2-AP136: Canonical chiral algebra is a conjunction of Costello--Li and Kontsevich--Tamarkin
**Severity**: High
**Trigger**: Vol II inscription citing either Costello--Li holomorphic
locality or Kontsevich--Tamarkin $E_d$-formality as the sole input to
the "canonical chiral algebra on $X$", treating one as sufficient.
**Ghost**: Costello--Li 2016 establishes holomorphic locality of
$\Phi$ on a CY datum (inputs from the Costello--Gwilliam BV
framework); Kontsevich--Tamarkin $E_d$-formality establishes rigidity
of the formality quasi-isomorphism on Hochschild complexes. Both are
load-bearing and both enter the canonical chiral algebra construction.
Neither is a shadow of the other; holomorphic locality without
formality produces a chain-level cocycle without homotopy control;
formality without holomorphic locality produces a rigidity statement
without the specifically-holomorphic factorisation structure.
**Precise error**: single-citation treatment allows silently dropping
a load-bearing input. The $\mathsf{SC}^{\mathrm{ch,top}}$ Vol II
assembly uses both inputs at distinct stages of the Dunn--curved-Dunn
tower; citing only one hides where the other enters.
**Correction**: the canonical chiral algebra is the conjunction of
Costello--Li holomorphic locality on a factorisation datum AND
Kontsevich--Tamarkin formality on its Hochschild cochains. Inscribe
both as hypotheses; at stage transitions (AP-V2-25), name which input
supports which stage.
**Primary**: Costello--Li 2016 arXiv:1605.09439; Kontsevich 2003
*Lett Math Phys* 66; Tamarkin 2003 *Lett Math Phys* 66; Francis 2013
*Compositio* 149.
**Inscription**: `chapters/theory/canonical_chiral_algebra.tex`; Vol
III wave-12 source `notes/wave12_a1_phi_functor_foundations.tex`.
**Cross-ref**: Vol II AP-V2-6 (Kontsevich formality scope); AP-V2-25
(stage-1/stage-2 identification); AP-V2-30 (chain vs $(\infty,1)$-
functor); Vol III Pattern 273 / Pattern 236.

### AP-V2-34 / V2-AP137: $E_n$ on the wrong object — $A$ is $E_1$ at $d \ge 3$; $E_2$ lives on $Z(\mathrm{Rep}(A))$
**Severity**: Critical
**Trigger**: Vol II inscription asserting an $E_2$- or $E_n$-algebra
structure directly on $A$ for a CY$_d$ input at $d \ge 3$, or claiming
Dunn additivity produces an $E_2$-bracket on $A$ itself.
**Ghost**: at $d \ge 3$, the $\Phi$-image $A$ is $E_1$-associative
(chiral-associative on a punctured complex disc); the $E_2$-structure
lives on the Drinfeld centre $Z(\mathrm{Rep}(A))$, not on $A$. The
Swiss-cheese $\mathsf{SC}^{\mathrm{ch,top}}$ closed colour is $E_2$
on the braided-monoidal representation category, not on the
associative algebra directly.
**Precise error**: categorical-layer confusion. Dunn additivity
$\int_{\mathbb R^2} E_1 \simeq E_2$ raises $E_1$ on a $0$-dimensional
input to $E_2$ on a $2$-dimensional integral; at $d \ge 3$ the
integration target shifts and the $E_2$ lives on
$Z(\mathrm{Rep}(A))$, visible through the Deligne conjecture and
categorical Hochschild. Placing $E_2$ on $A$ collapses the
Drinfeld-centre face of $r_{\mathrm{CY}}$ and forces nonsensical
monoidal conditions on an $E_1$-algebra.
**Correction**: label every chiral-topological $E_n$-assertion with
the object it sits on. $A$ itself: $E_1$-associative at $d \ge 3$.
$Z(\mathrm{Rep}(A))$: $E_2$-braided. $\mathsf{SC}^{\mathrm{ch,top}}$
closed colour at $d \ge 3$: $E_2$ on $Z(\mathrm{Rep}(A))$, not on $A$.
At $d = 2$, $A$ itself carries $E_2$-structure via the PTVV row
$(d, \mathrm{shift}, E_n^{\mathrm{cl}}) = (2, -2, E_2)$; the
categorical layering is $d$-dependent.
**Primary**: Lurie *Higher Algebra* Ch 5 (Dunn); Francis 2013
*Compositio* 149; Ayala--Francis 2015 (factorisation homology on
stratified spaces).
**Inscription**: `chapters/theory/e_n_hierarchy_dimensional.tex`; Vol
III wave-12 source `notes/wave12_a8_en_hierarchy_kapranov.tex`.
**Cross-ref**: Vol III AP-CY $E_n$-output-scope family; Vol I Pattern
236 ambient discipline; Vol II AP-V2-1 (Swiss-cheese closed colour
carries the right $E_2$ scope only on $Z(\mathrm{Rep}(A))$ at $d \ge 3$);
AP-V2-37 (the $E_2$-braiding specialisation of this AP).

### AP-V2-35 / V2-AP138: PTVV $d$-dependent shift table $(d, \mathrm{shift}, E_n^{\mathrm{cl}})$
**Severity**: Critical
**Trigger**: Vol II inscription treating the PTVV shift as
$d$-independent in the chiral-topological specialisation, or stating a
single uniform $\Phi$-output type across CY dimensions.
**Ghost**: the CY-to-chiral functor output is dimensionally stratified
by a shifted-symplectic PTVV input of degree $2 - d$. At $d = 2$ shift
$-2$, classical limit $E_2$-Poisson; at $d = 3$ shift $-1$,
$E_1$-Poisson; at $d = 4$ shift $0$, $E_0$-Poisson; at $d = 5$ shift
$+1$, $E_5$-Poisson and the direction reverses. The table
$(d, \mathrm{shift}, E_n^{\mathrm{cl}}) \in \{(2,-2,E_2), (3,-1,E_1),
(4,0,E_0), (5,+1,E_5\text{-Poisson})\}$ controls the Vol II
chiral-topological specialisation on each row.
**Precise error**: uniform $\Phi$-output treatment erases the PTVV
staircase and collapses the three-faces universal identity (AP-V2-20)
into a single-shift identity, losing the $d$-dependent matching
between shift, $E_n$-class, and Quillen exponent.
**Correction**: every $\Phi$-application carries the $d$-subscript
$\Phi_d$ and the associated $(d, \mathrm{shift}, E_n^{\mathrm{cl}})$
triple. At cross-volume linking points, pin the row; the K3 $\times$ E
site sits at $d = 3$ via AP-V2-25.
**Primary**: Pantev--Toën--Vaquié--Vezzosi 2013 *Publ IHES* 117 (PTVV);
Costello--Gwilliam 2021 Vol 2 §3 (BV version of PTVV); Calaque--Pantev--
Toën--Vaquié--Vezzosi 2017 *Selecta Math* 23.
**Inscription**: `chapters/theory/ptvv_staircase.tex`; Vol III wave-12
source `notes/wave12_a1_phi_functor_foundations.tex`.
**Cross-ref**: Vol III $\Phi$-output-scope family; Vol II AP-V2-14
(generic CY-4 $\Phi_4$ existence) uses row $d = 4$ shift $0$;
AP-V2-20 three-faces identity carries row-specific Quillen exponents;
AP-V2-34 ($E_n$-on-wrong-object) depends on this row.

### AP-V2-36 / V2-AP139: Dunn--Lurie $\int_{\Sigma_2} E_3 \simeq E_1$ real-vs-complex conflation
**Severity**: High
**Trigger**: Vol II inscription transporting Dunn--Lurie additivity
$\int_{\Sigma_2} E_3 \simeq E_1$ to a holomorphic setting without
separating real and complex integration dimensions, or applying the
real-dimension additivity count to a complex curve without tracking
the transverse normal bundle.
**Ghost**: Dunn--Lurie additivity is a real-dimension statement:
$\int_{\mathbb R^2} E_3 \simeq E_1$ uses real-dimension-$2$
integration. In the holomorphic chiral-topological setting,
integration is along a complex curve $C$ of complex dimension $1$ =
real dimension $2$ with transverse normal bundle $N_C$ of complex
codimension $d - 1$; the additivity bookkeeping must keep real and
complex dimensions explicit, and the transverse normal bundle enters
the $E_n$-count via the $\mathrm{Td}(N_C)$ factor.
**Precise error**: slip between real and complex dimensions: writing
$\int_{C} E_3 \simeq E_1$ as a clone of the real Dunn--Lurie identity
without the transverse normal-bundle correction silently omits the
holomorphic-versus-smooth distinction and gives wrong $E_n$-types in
the Swiss-cheese closed colour. The MIXED-operation
$\int_{\Sigma_{d-1}}$ of AP-V2-25 is specifically the
holomorphic-factorisation-homology pushforward, not the smooth
Dunn--Lurie integral.
**Correction**: Vol II chiral-topological integrals are written
$\int^{\mathrm{hol}}_{C} E_n \simeq E_{n-2}$ with explicit
$\mathrm{Td}(N_C)$ and normal-bundle weighting, labelling real- vs
complex-dimension at every step. Dunn--Lurie remains the ambient
real-dimension additivity statement; the holomorphic refinement is
Costello--Li holomorphic factorisation homology.
**Primary**: Lurie *Higher Algebra* Thm 5.1.2.2 (Dunn); Costello--Li
2016 arXiv:1605.09439 (holomorphic locality); Francis 2013
*Compositio* 149.
**Inscription**: `chapters/theory/holomorphic_dunn.tex`; Vol III
wave-12 source `notes/wave12_f2_bv_brst_to_chiral_ce.tex`.
**Cross-ref**: Vol II AP-V2-34 ($E_n$-on-wrong-object); Vol III
$E_n$-output-scope family; Vol I Pattern 236 ambient discipline;
AP-V2-25 MIXED operation as the holomorphic pushforward.

### AP-V2-37 / V2-AP140: Iterated-Sugawara $\mathcal W_{1+\infty}$-endpoint conjecture — chain-level scope with three inputs
**Severity**: High (conjectural; tag \ClaimStatusConjectured)
**Trigger**: Vol II inscription asserting the iterated-Sugawara
$\mathcal W_{1+\infty}$ endpoint as a theorem without naming the three
load-bearing inputs: Miki $S_3$-triality on $Y^+$, $\mathcal W_{1+\infty}$
as the evaluation slice, and curved-Dunn $H^2 = 0$ closure.
**Ghost**: the proposed endpoint composes three structures. (i) Miki
$S_3$-triality on $Y^+(\widehat{\mathfrak{gl}}_1)$ pins the
$\mathcal W_{1+\infty}$ slice inside $Y^+$; (ii)
$\mathcal W_{1+\infty}$ is the evaluation-slice image of $Y^+$ under
the Schiffmann--Vasserot arrow (AP-V2-32); (iii) curved-Dunn $H^2 = 0$
closure is the Costello--Gwilliam BV condition that the curved-Dunn
assembly of the Vol II chiral-topological factorisation algebra lifts
$\mathcal W_{1+\infty}$ to the iterated-Sugawara tower without anomaly.
**Precise error**: treating the endpoint as a single citation-level
theorem glosses all three inputs; each is independently load-bearing
and independently conjectural at present.
**Correction**: state the endpoint as a conditional theorem. "If
(i)-(ii)-(iii) hold, then the Vol II chiral-topological iterated-
Sugawara tower converges to $\mathcal W_{1+\infty}$ via the curved-Dunn
closure." Tag \ClaimStatusConjectured until the three inputs are
independently verified. The chain-level formulation is where each
input can be checked directly.
**Primary**: Miki 2007 *J Math Phys* 48 (triality); Schiffmann--
Vasserot 2020; Costello--Gwilliam 2021 Vol 2 §4 (curved Dunn).
**Inscription**: `chapters/frontier/w_infty_endpoint_conjecture.tex`;
Vol III wave-12 source `notes/wave12_a11_coha_y_plus_vs_w_infty.tex`.
**Cross-ref**: Vol II AP-V2-32 (CoHA evaluation scope); AP-V2-37 flag
propagates to AP-V2-38 (class M $E_3$-bar cohomology) via the closure
hypothesis; AP-V2-41 (frontier conjectural bar-filtered equivalence).

### AP-V2-38 / V2-AP141: Class $\mathbf M$ $E_3$-bar cohomology is $6^g$ with per-handle tricomplex $P_{E_4}(t) = 3t + 3t^2$
**Severity**: High
**Trigger**: Vol II inscription asserting class-$\mathbf M$ $E_3$-bar
is infinite at cohomology, or finite at chain level, at genus
$g \in \{1, 2, 3\}$.
**Ghost**: class-$\mathbf M$ $E_3$-bar has chain-level
infinite-dimensional and cohomology $6^g$-dimensional via the kernel
of $d_4$ on the single-generator tricomplex with per-handle Poincaré
polynomial $P_{E_4}(t) = 3t + 3t^2$. The factor $6$ per handle is the
closed-form Euler-characteristic count $2 \cdot 3$ with sign
discipline; the $6^g$-scaling matches the genus-$g$ tensor power of
the per-handle profile. The wave-12 contribution is the explicit
tricomplex that AP-V2-21 named without unpacking.
**Precise error**: lane conflation (cf.\ AP-V2-21), strengthened by
the explicit tricomplex structure from wave-12. Treating the
$6^g$-count as a chain-level count, or the infinite-dimensional
chain level as the cohomological answer, gives wrong genus-scaling in
either direction.
**Correction**: cohomology is $6^g$, chain level infinite; state both
lanes and the spectral-sequence collapse at $d_4$, with the per-handle
Poincaré polynomial $P_{E_4}(t) = 3t + 3t^2$ written explicitly.
**Primary**: Francis 2013 *Compositio* 149; Costello--Gwilliam 2021
Vol 2 §9.
**Inscription**: `chapters/theory/class_m_e3_bar.tex`; Vol III wave-12
source `notes/wave12_a10_class_m_e3_bar.tex`.
**Cross-ref**: AP-V2-21 (original class-M AP, strengthened here with
per-handle tricomplex); Vol III AP-CY21; Vol II
`first_principles_cache.md` canonical preamble row 65; AP-V2-41
(filtered equivalence on $E_3$-bar is the conjectural lift).

### AP-V2-39 / V2-AP142: $E_2$-braiding on $Z(\mathrm{Rep}(A))$, not on $A$
**Severity**: Critical
**Trigger**: Vol II inscription asserting an $E_2$-braiding directly
on $A$ for a chiral-topological datum at $d \ge 3$, or placing the
Drinfeld centre content on the underlying algebra rather than on its
representation category.
**Ghost**: at $d \ge 3$ the $E_2$-braided-monoidal structure lives on
the Drinfeld centre $Z(\mathrm{Rep}(A))$ (a fusion-categorical
layer), not on $A$ itself, which is $E_1$-associative. This is the
braiding-face specialisation of AP-V2-34.
**Precise error**: placing the $E_2$-braiding on the wrong
categorical layer collapses the Drinfeld-centre face of
$r_{\mathrm{CY}}$ and forces nonsensical monoidal conditions on an
$E_1$-algebra. The $\mathsf{SC}^{\mathrm{ch,top}}$ closed colour at
$d \ge 3$ acts on the centre; the open colour acts on $A$; conflating
the two destroys the two-stage factorisation of AP-V2-25.
**Correction**: state every $E_2$-braiding fact as a statement about
$Z(\mathrm{Rep}(A))$; the underlying $A$ carries $E_1$-associative
structure only. Deligne centre and categorical Hochschild are the
machine that raises $E_1$ on $A$ to $E_2$ on the centre.
**Primary**: Lurie *Higher Algebra* 5.3 (Drinfeld centre); Etingof--
Nikshych--Ostrik 2005 *Ann Math* 162; Ayala--Francis 2015.
**Inscription**: `chapters/theory/drinfeld_centre_braiding.tex`; Vol
III wave-12 source `notes/wave12_a3_MO_E2_etingof.tex`.
**Cross-ref**: AP-V2-34 ($E_n$-on-wrong-object, parent AP); Vol III
Drinfeld-centre-vs-averaging family; Vol I Pattern 236 ambient
discipline; AP-V2-25 closed-colour identification with stage-1
$\Phi^{\mathrm{FA}}_d$.

### AP-V2-40 / V2-AP143: Six routes to $G(K3 \times E)$ are six distinct constructions, NOT six $\Phi_3$-applications
**Severity**: High
**Trigger**: Vol II inscription treating the six routes (CoHA,
Schiffmann--Vasserot, Maulik--Okounkov, Borcherds, Toda, DMVV) as six
applications of a single functor $\Phi_3$ to one object, rather than
six distinct constructions on six different CY-input categories that
produce the same $\Phi_3$-output via a pentagon colimit.
**Ghost**: each of the six routes takes a different CY-input
category and a different intermediate structure (Hall, stable-envelope,
Fock-space, Borcherds-product, topological-vertex, wall-crossing) and
produces $\mathbf H_{\Delta_5} = G(K3 \times E)$ as its
$\Phi_3$-image. The pentagon that unifies them is a colimit
diagram, not a single-object application of $\Phi_3$ six times. This
is the Vol III K3$\times$E specific crystallisation; Vol II's bar-
cobar / pentagon / heptagon cross-volume readings must respect it.
**Precise error**: single-functor framing of the six routes implies
six applications of $\Phi_3$ to one CY-input category and hides the
six distinct CY-input structures.
**Correction**: name the six distinct CY-input categories and name
the pentagon colimit that produces the shared output. Each route's
$\Phi_3$-image is $\mathbf H_{\Delta_5}$; the routes agree at the
colimit, not through single-object application.
**Primary**: Schiffmann--Vasserot 2012 (CoHA); Maulik--Okounkov 2019
(stable envelope); Borcherds 1992 / Gritsenko--Nikulin 1998
(Borcherds product); Toda 2013 (DT wall-crossing); Dijkgraaf--Moore--
Verlinde--Verlinde 1997 (DMVV).
**Inscription**: `chapters/theory/six_routes_pentagon.tex`; Vol III
wave-12 source `notes/wave12_a12_six_routes_k3_e.tex`.
**Cross-ref**: AP-V2-28 (manifesto conflations, row 3); Vol III
CLAUDE.md key-facts block (six routes bullet); Vol II
`first_principles_cache.md` canonical preamble rows 59--62.

### AP-V2-41 / V2-AP144: $(B_{E_1}, B_{E_2}, B_{E_\infty})$ filtered equivalence at chain level — frontier conjecture
**Severity**: Medium (frontier open question; tag
\ClaimStatusConjectured)
**Trigger**: Vol II inscription claiming the three-sector filtered
equivalence
$B_{E_1}A \twoheadrightarrow B_{E_2}A \twoheadrightarrow B_{E_\infty}A$
as a proved theorem without naming the chain-level scope where it
could be made rigorous, or dismissing it as inaccessible without
naming that scope.
**Ghost**: the filtered equivalence is a frontier open question from
the wave-12 cross-volume synthesis. Conjecturally, the filtered
equivalence holds at chain level for $\Phi$-image algebras when the
three inputs of AP-V2-33 are assembled, with the $E_2$-layer sitting
on $Z(\mathrm{Rep}(A))$ per AP-V2-39 and the real-vs-complex dimension
discipline per AP-V2-36.
**Precise error**: asserting the equivalence as theorem without
naming filtered bar models; alternatively, dismissing it as
inaccessible without naming the chain-level site where it becomes
rigorous.
**Correction**: log as conjecture at chain level with explicit
prerequisites. Scope: for $A = \Phi_d(X)$ with $X$ a CY$_d$ datum
satisfying Costello--Li holomorphic locality and Kontsevich--Tamarkin
$E_d$-formality, the filtered bar tower $B_{E_n}A$ has graded
quotients controlled by the kernel of $d_4$ on the per-handle
tricomplex (cf.\ AP-V2-38). Chain-level verification at small $g$ and
small $d$ is the route; $(\infty,1)$-categorical-functor lift (cf.\
AP-V2-30) is a separate conjecture.
**Primary**: Francis 2013 *Compositio* 149; Kontsevich--Tamarkin
2003; Costello--Gwilliam 2021 Vol 2.
**Inscription**: `chapters/frontier/bar_filtered_equivalence_conjecture.tex`;
Vol III wave-12 source `notes/wave12_a1_bar_cobar_gelfand.tex`.
**Cross-ref**: Vol II AP-V2-30 (chain vs $(\infty,1)$-functor);
AP-V2-38 (class M $E_3$-bar cohomology supplies the per-handle
obstruction); AP-V2-33 (two-input canonical chiral algebra supplies
the prerequisite inputs); Vol III Pattern 273 / Vol I Pattern 273
hook; AP-V2-37 (iterated-Sugawara endpoint uses the same three-input
closure).

## Section 15: 6d hCS audit session antipatterns (2026-04-22)

The following eighteen APs record precise errors caught during the
2026-04-22 audit of the Vol II 6d holomorphic Chern--Simons computation
(partition function $Z_{\mathrm{hCS}}$, Feynman coefficients, Deligne
exceptional series, weak-Jacobi / Siegel / paramodular denominators,
and the Vol III two-stage factorisation alignment for the $E_3$-chiral
avatar).

### AP-V2-42 / V2-AP145: Casimir quadratic-for-quartic in 6d hCS 1-loop (2026-04-22)

**Trigger**: any Vol II 6d hCS 1-loop amplitude statement that writes
$c_2(\mathfrak g)/24$ against a 3-factor propagator integrand for the
one-loop bubble / counterterm.

**Ghost**: the 6d hCS one-loop amplitude on $\mathbb R^4 \times C$ with
colour factor $\mathrm{tr}_{\mathrm{adj}}(T^a T^b T^c T^d)$ contracts
the adjoint quartic Casimir against the four-factor
$c(\bar\partial c)^3/(2\pi i)^3$ integrand of the BV propagator
(Costello 2017 *M-theory / $\Omega$-background*; Costello--Gaiotto
2019; CFG 2026).

**Precise error**: substituting $c_2(\mathfrak g)/24$ (the quadratic
Casimir normalisation used for a 3-factor integrand) for the quartic
Casimir $\mathrm{tr}_{\mathrm{adj}}(T^{(a}T^b T^c T^{d)})$ against the
4-factor propagator integrand. The error appears at
\texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4005-4015}.
The quadratic Casimir contracts an $(a, b)$ pair; the 6d hCS
one-loop amplitude requires the symmetric quartic invariant on the
adjoint because the Feynman graph has four external colour legs and
three internal propagators, giving the $c(\bar\partial c)^3/(2\pi i)^3$
weight.

**Correction**: replace $c_2(\mathfrak g)/24 \cdot (\cdots 3\text{-factor})$
with $\mathrm{tr}_{\mathrm{adj}}(T^{(a}T^b T^c T^{d)}) \cdot (\cdots
c(\bar\partial c)^3/(2\pi i)^3 )$, where the symmetrisation runs over
the four colour indices and $\bar\partial c$ is the BV propagator
one-form on $\mathrm{Conf}_2(\mathbb C^3)$ restricted to the one-loop
fibre. The anomaly-cancellation condition reduces to the vanishing
of $\mathrm{tr}_{\mathrm{adj}} T^{(a} T^b T^c T^{d)}$ on the Deligne
exceptional series except $E_6$ (see AP-V2-54 below).

**Severity**: Critical. The error propagates into every downstream
anomaly and 1-loop-finiteness claim. Per-sweep when editing 6d hCS
amplitude computations.

**Primary**: Costello 2017 arXiv:1610.04144 "M-theory in the
$\Omega$-background and 5d/4d correspondence"; Costello--Gaiotto
2019 arXiv:1812.05516 "Twisted supergravity and Koszul duality";
Costello--Francis--Gaiotto 2026 (6d hCS one-loop formulae, in
preparation as of the 2026-04-22 audit); Costello--Gwilliam 2021
*Factorization Algebras Vol 2* §9-11 (BV propagator conventions).

**Cross-ref**: \texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4005-4015}
(error site); V2-AP54 / AP-V2-54 (Deligne-safe anomaly list); V2-AP146
/ AP-V2-43 (theta-graph two-loop miscount); V2-AP149 / AP-V2-49
($E_1$-per-$\mathbb C$ Dunn scope); Canonical preamble rows 59-66
($\kappa_\bullet$ indexing).

### AP-V2-43 / V2-AP146: Feynman graph Euler miscount — theta labelled 1-loop (2026-04-22)

**Trigger**: any Vol II Feynman graph labelling that tags the theta
graph ($V = 2$, $E = 3$) as 1-loop, or asserts a 2-leg 1-loop
amplitude with three internal edges.

**Ghost**: for a connected Feynman graph, Euler's relation
$V - E + L = 1$ (equivalently, $L = 1 + E - V$) gives the loop number
$L$; the theta graph has $V = 2$, $E = 3$ and therefore $L = 2$. The
1-loop 2-leg amplitude is the bubble graph ($V = 2$, $E = 2$,
$L = 1$). This is undergraduate-level Feynman graph topology; the
error is a bookkeeping slip, not a framework failure.

**Precise error**: site
\texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:120-122}
labels the theta graph (two vertices connected by three parallel
edges) as 1-loop. Correct count: $L = 1 + E - V = 1 + 3 - 2 = 2$;
the theta is 2-loop. The 1-loop 2-leg graph is the bubble
($V = 2$, $E = 2$, $L = 1$), which is topologically distinct and
has a different symmetry factor.

**Correction**: retag the theta graph as 2-loop; reassign the
1-loop 2-leg amplitude to the bubble graph. The symmetry factors
and integral measures change accordingly (theta: symmetry factor
$1/12$ over $S_3 \times \mathbb Z_2$; bubble: $1/2$ over $\mathbb Z_2$).

**Severity**: High. Propagates into the loop-expansion-order
counting of every amplitude downstream.

**Primary**: Peskin--Schroeder 1995 §4.4 (Feynman graph topology);
Itzykson--Zuber 1980 §6-2-2 (Euler characteristic of connected
graphs).

**Cross-ref**: \texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:120-122}
(error site); AP-V2-42 / V2-AP145 (quartic Casimir / 3-factor
propagator — adjacent bookkeeping); AP-V2-46 / V2-AP149 (wedge of
1-forms — adjacent topology slip).

### AP-V2-44 / V2-AP147: Weak Jacobi form polar coefficient $c_{\phi_{-2,1}}(-1) = 24$ false (2026-04-22)

**Trigger**: any Vol II citation of $c_{\phi_{-2,1}}(-n) = 24$ for
$n \ge 1$, or the identification of $24$ as a weak Jacobi form polar
Fourier coefficient of $\phi_{-2,1}$.

**Ghost**: the weak Jacobi form $\phi_{-2,1}(\tau, z)$ of weight
$-2$ and index $1$ has Fourier expansion
$\phi_{-2,1}(\tau, z) = \sum c_{\phi_{-2,1}}(n, \ell) q^n y^\ell$; the
polar coefficients $c_{\phi_{-2,1}}(-n, \ell) = 0$ for $n \ge 1$
by Eichler--Zagier polar-support cutoff
$c(\Delta) = 0$ for $\Delta = 4nm - \ell^2 < -m^2 = -1$ at
index $m = 1$ (Eichler--Zagier 1985 *Prog Math* 55 Thm 9.1, standard).

**Precise error**: sites
\texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:130-131}
and
\texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:186-188},
and
\texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4184}
assert $c_{\phi_{-2,1}}(-1) = 24$. The value $24$ arises either as
$\chi(\mathrm{K3})$ (the topological Euler characteristic of K3) or
as $|c_{1/\Delta_5}(q^1)|$ (a paramodular / $1/\Delta_5$-expansion
coefficient), NOT as a polar Fourier coefficient of
$\phi_{-2,1}$. Polar coefficients of index-1 weak Jacobi forms
vanish identically for $\Delta < -1$.

**Correction**: $c_{\phi_{-2,1}}(-n) = 0$ for all $n \ge 1$ by
Eichler--Zagier polar-support cutoff. When $24$ is the intended
value, replace with $\chi(\mathrm{K3}) = 24$ (topological Euler
characteristic, Hirzebruch 1966) or with
$|c_{1/\Delta_5}(q^1)| = 24$ (paramodular index-1 $q$-expansion
coefficient, Gritsenko 1999 *St Petersburg Math J* 11 Thm 6.1).

**Severity**: High. The error conflates three different combinatorial
objects (topological Euler characteristic, paramodular $q$-coefficient,
polar weak-Jacobi coefficient) that all happen to equal 24. The
downstream Feynman coefficient claim is false; the correction
changes amplitude formulas.

**Primary**: Eichler--Zagier 1985 *Prog Math* 55 Thm 9.1
(polar-support cutoff); Gritsenko 1999 (paramodular $\Delta_5$
expansion); Hirzebruch 1966 (K3 Euler characteristic).

**Cross-ref**: \texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:130-131, 186-188}
and \texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4184}
(error sites); AP-V2-24 / V2-AP127 (Humbert--Heegner admissibility
filter using the same Eichler--Zagier cutoff); Canonical preamble
row 45 (Humbert--Heegner admissibility).

### AP-V2-45 / V2-AP148: Stuffle dimensional failure — $\zeta(3)^2$ mixed-weight (2026-04-22)

**Trigger**: any Vol II identity for $\zeta(3)^2$ that mixes weights
or drops the stuffle relation to a non-homogeneous right-hand side.

**Ghost**: Euler stuffle relation for multiple zeta values
$\zeta(a) \zeta(b) = \sum \zeta(\cdots)$ preserves weight: both
sides are homogeneous of weight $a + b$. For $\zeta(3)^2$ the
stuffle expansion gives
$\zeta(3)^2 = 2\zeta(3, 3) + \zeta(6)$, all of weight 6
(Zagier 1994 *Prog Math* 120 Eq 4; Waldschmidt 2011 Exp 5).

**Precise error**: sites
\texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:334, 382}
assert $\zeta(3)^2 = 3\zeta(3, 3) + 6\zeta(5)/5$. The right-hand
side mixes $\zeta(3, 3)$ (weight 6) with $\zeta(5)$ (weight 5);
this is dimensionally illegal under the weight grading of the MZV
algebra. The correct stuffle yields $2\zeta(3, 3) + \zeta(6)$, not
$3\zeta(3, 3) + (6/5)\zeta(5)$.

**Correction**: $\zeta(3)^2 = 2\zeta(3, 3) + \zeta(6)$
(Zagier 1994 *Prog Math* 120 Eq 4). Under Euler's $\zeta(2n)$
closed form, $\zeta(6) = \pi^6/945$. If the computation needs a
$\zeta(5)$ term, it must come from a separate weight-5 identity
(e.g., shuffle $\zeta(2)\zeta(3) = \zeta(2, 3) + \zeta(3, 2) +
\zeta(5)$), not from the stuffle of $\zeta(3)^2$.

**Severity**: Critical. The mixed-weight RHS falsifies the
stuffle relation as a graded algebra identity; downstream MZV
algebra claims inherit the type error.

**Primary**: Zagier 1994 *Prog Math* 120 "Values of zeta functions
and their applications" Eq 4; Waldschmidt 2011 "Lectures on
Multiple Zeta Values" Exp 5; Brown 2012 *Ann Math* 175 Thm 1
(motivic MZV algebra grading).

**Cross-ref**: \texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:334, 382}
(error sites); cache row 129 (Vol II weight-12 MZV
$\zeta(3, 3, 3, 3)$); V2-AP124 (MZV numerical value).

### AP-V2-46 / V2-AP149: Wedge of 1-forms $\eta_{12}^{\wedge 3} = 0$ (2026-04-22)

**Trigger**: any Vol II iterated-integral formula of the form
$\int_{\cdots} \eta_{ij}^{\wedge k}$ with $k \ge 2$ for a fixed
1-form $\eta_{ij}$.

**Ghost**: for any 1-form $\eta$ on any manifold,
$\eta \wedge \eta = 0$ by graded anticommutativity of the exterior
algebra (Spivak 1979 Vol 1 Ch 7). Iterated Chen / Bott--Taubes
integrals require $k$ \emph{distinct} 1-forms on $k$ \emph{distinct}
configuration slots; degenerate evaluations with repeated forms
vanish identically.

**Precise error**: site
\texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:253-255}
asserts $\int_{z_2} \eta_{12}^{\wedge 3} = \zeta(3)$. But
$\eta_{12}^{\wedge 3} = 0$ identically because $\eta_{12}$ is a
single 1-form wedge-anticommuting with itself. The $\zeta(3)$
identity requires the Kontsevich / Hanamura / Brown three-distinct-form
integral
$\int \eta_{12} \wedge \eta_{23} \wedge \eta_{13}$ over
$\mathrm{Conf}_4$ or similar nondegenerate configuration
(Kontsevich 1994 *Prog Math* 120; Brown 2012 *Ann Math* 175).

**Correction**: replace $\eta_{12}^{\wedge 3}$ with a nondegenerate
wedge of three distinct 1-forms such as
$\eta_{12} \wedge \eta_{23} \wedge \eta_{31}$ (cubic Arnold) or
$\eta_{12} \wedge \eta_{13} \wedge \eta_{14}$ (star) on
$\mathrm{Conf}_n$ for $n \ge 4$; compute the resulting integral
against the appropriate holomorphic propagator.

**Severity**: High. The $\zeta(3)$ identification is the load-bearing
anchor of Vol II's weight-3 one-loop amplitude; the degenerate form
destroys the identification on the nose.

**Primary**: Spivak 1979 Vol 1 Ch 7 (graded-commutative exterior
algebra); Kontsevich 1994 *Prog Math* 120 "Feynman diagrams and
low-dimensional topology"; Brown 2012 *Ann Math* 175 (motivic MZV
basis, iterated integrals).

**Cross-ref**: \texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:253-255}
(error site); cache row 17 (Kontsevich integral conventions).

### AP-V2-47 / V2-AP150: Residue ill-typed — $\mathrm{Res}_{\Phi_{10} = 0}(Z_{\mathrm{hCS}}) = 1/\Phi_{10}$ (2026-04-22)

**Trigger**: any Vol II identification of
$\mathrm{Res}_{\Phi_{10} = 0}(Z_{\mathrm{hCS}}) = 1/\Phi_{10}$, or a
residue identification at the Igusa cusp-form vanishing locus that
produces a denominator rather than a numerator.

**Ghost**: the residue operator
$\mathrm{Res}_{f = 0} : \Omega^n_X \to \Omega^{n-1}_{\{f = 0\}}$
\emph{removes} the denominator $1/f$; if $\omega = g/f \cdot dz$
near $\{f = 0\}$, then
$\mathrm{Res}_{f = 0}(\omega) = g|_{f = 0}$. The statement
$\mathrm{Res}_{\Phi_{10} = 0}(Z_{\mathrm{hCS}}) = 1/\Phi_{10}$
is type-wrong: a residue cannot be a denominator of the same form.

**Precise error**: site
\texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4246-4277}
writes
$\mathrm{Res}_{\Phi_{10} = 0}(Z_{\mathrm{hCS}}) = 1/\Phi_{10}$.
The correct identification is
$Z_{\mathrm{3d\,QG}}(\tau, \tau') = 1/\Phi_{10}(\tau, z, \tau')$
\emph{directly} (Dijkgraaf--Verlinde--Verlinde 1997 *NPB* 484,
the exact Igusa cusp form
Siegel denominator of 3d quantum gravity partition function),
with the Costello effective action giving the Borcherds singular
theta lift of the chiral genus (Borcherds 1998 *Invent Math*
132 Thm 1.7).

**Correction**:
$Z_{\mathrm{3d\,QG}}(\tau, z, \tau') = 1/\Phi_{10}(\tau, z, \tau')$
(DVV 1997). The residue / Borcherds lift identification is
$Z^{\mathrm{eff}}_{\mathrm{Cost}} = -\log \Phi_{10}$ as a
singular theta lift of the elliptic genus; no residue step
produces $1/\Phi_{10}$ from a 1-form.

**Severity**: Critical. The identification anchors the 3d quantum
gravity Part VI climax; the residue-as-denominator framing is
false on the nose.

**Primary**: Dijkgraaf--Verlinde--Verlinde 1997 *NPB* 484
"Counting dyons in $\mathcal N = 4$ string theory"; Borcherds
1998 *Invent Math* 132 Thm 1.7 (singular theta lifts);
Costello 2017 (effective action); Griffiths--Harris 1978 Ch 1.1
(residue operator).

**Cross-ref**: \texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4246-4277}
(error site); AP-V2-48 / V2-AP151 ($\eta^{24}$ denominator);
Canonical preamble row 28 ($\Phi_{10} = \Delta_5^2$ K3-BKM
denominator).

### AP-V2-48 / V2-AP151: Single $\eta^{24}$ denominator — both elliptic-curve factors needed (2026-04-22)

**Trigger**: any Vol II holomorphic partition function of the form
$(\Phi_{10}/\eta^{24})^\hbar$ or $\Phi_{10}/\eta(\tau)^{24}$ with a
single $\eta$-factor in the denominator.

**Ghost**: the Igusa cusp form $\Phi_{10}(\tau, z, \tau')$ is a
Siegel modular form of weight 10 on $\mathrm{Sp}_4(\mathbb Z)$
with variables $(\tau, z, \tau') \in \mathbb H_2$; the two torus
moduli $\tau, \tau'$ parametrise the two elliptic-curve factors
$E_\tau \times E_{\tau'}$ in the Kummer / Humbert-divisor local
geometry. The canonical 3d quantum gravity / DVV partition function
has Kontsevich-style weight cancellation requiring
$\Phi_{10}/(\eta(\tau)^{24} \eta(\tau')^{24})$ (two $\eta^{24}$
factors, one per elliptic curve).

**Precise error**: site
\texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4240-4264}
writes $(\Phi_{10}/\eta^{24})^\hbar$ with a single $\eta^{24}$
denominator. The single denominator does not cancel the second
torus weight; under the transformation
$\mathrm{Sp}_4(\mathbb Z) \to \mathrm{SL}_2(\mathbb Z)_\tau
\times \mathrm{SL}_2(\mathbb Z)_{\tau'}$, only $\eta(\tau)^{24}$
is restored while $\eta(\tau')^{24}$ is missed.

**Correction**:
$Z_{\mathrm{hCS}} = \Phi_{10}/(\eta(\tau)^{24} \eta(\tau')^{24})$.
Both elliptic-curve weight factors appear; the combined denominator
has weight $10 - 12 - 12 = -14$ under $\mathrm{Sp}_4$, exhibiting
the expected singular-weight-$10$ numerator cancellation against
the two weight-$12$ cusp denominators, yielding the weight-$(-14)$
Siegel meromorphic form.

**Severity**: High. The single-factor denominator breaks the
modular weight balance and changes the 3d QG partition function
asymptotics.

**Primary**: DVV 1997 *NPB* 484; Gritsenko--Nikulin 1998
*J Reine Angew Math* 507 §3; Borcherds 1998 *Invent Math* 132
(singular theta lifts).

**Cross-ref**: \texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4240-4264}
(error site); AP-V2-47 / V2-AP150 (residue type error); Canonical
preamble row 28 ($\Phi_{10} = \Delta_5^2$ two-elliptic-curve
denominator).

### AP-V2-49 / V2-AP152: Dunn $E_1$-per-$\mathbb C$ — three commuting $E_1$ from three complex coordinates (2026-04-22)

**Trigger**: any Vol II statement that the three complex coordinates
$(z_1, z_2, z_3)$ on $\mathbb C^3$ yield three commuting $E_1$-algebra
structures.

**Ghost**: Dunn's additivity theorem
$E_n \otimes E_m \simeq E_{n+m}$ applies to the \emph{topological}
little-disks operads. A single complex coordinate on $\mathbb C$
supplies an $E_2$-structure topologically (since
$\mathrm{FM}_k(\mathbb C)$ is real-dim $2k$ in topological form); the
$E_1$-structure (linear ordering) emerges only after
$\bar\partial$-cohomology is taken (Costello--Gwilliam 2017 *FA*
Vol 2 Ch 5). Three $\mathbb C$ coordinates give $E_6$ topologically,
descending to $E_3$ after Dolbeault cohomology in each coordinate
(CG Vol II Ch 5 Thm 5.6.1).

**Precise error**: site
\texttt{chapters/connections/six\_d\_hcs\_e3\_chiral\_avatar\_platonic.tex:462-465}
asserts three commuting $E_1$ algebras from three complex
coordinates on $\mathbb C^3$, combining via Dunn to give $E_3$.
The combinatorial shape is correct ($E_3$ output) but the input
types are wrong: each $\mathbb C_j$ supplies $E_2$ topologically
(real-dim 2), not $E_1$; the $E_1$ reduction happens only at the
Dolbeault-cohomology level. Naïve Dunn additivity across three
$E_1$'s would give $E_3$; the correct derivation is $E_6
\to E_3$ via Dolbeault.

**Correction**: each $\mathbb C$ coordinate carries an
$E_2$-holomorphic factorisation structure topologically; three
coordinates compose via Dunn to $E_6$; Dolbeault cohomology in
each $\bar\partial_j$ direction reduces each $E_2$ to $E_1$, so
$E_6 \to E_3$. The $E_3$ output on $\mathbb C^3$ after Dolbeault
is what matches the 6d hCS factorisation structure (Costello--Li
2016 arXiv:1605.09930 §6).

**Severity**: High. The Dunn-composition framing has been widely
repeated; the correct derivation shows the $E_3$-level is forced
by Dolbeault (not by independent $E_1$-structures), which is
load-bearing for the $E_3$-chiral avatar of the Vol II
$\mathsf{SC}^{\mathrm{ch,top}}$ heptagon.

**Primary**: Costello--Gwilliam 2017 *FA* Vol 2 Ch 5 Thm 5.6.1
(Dolbeault reduction $E_{2d} \to E_d$); Costello--Li 2016
arXiv:1605.09930 §6 (6d hCS factorisation structure);
Dunn 1988 *J Pure Appl Algebra* 50 (Dunn additivity).

**Cross-ref**: \texttt{chapters/connections/six\_d\_hcs\_e3\_chiral\_avatar\_platonic.tex:462-465}
(error site); AP-V2-50 / V2-AP153 (Willwacher wheel-to-$\zeta$
convention); cache row 141 (ambient $\mathbb E_n$-level discipline);
Canonical preamble ambient-$\mathbb E_n$ row.

### AP-V2-50 / V2-AP153: Willwacher wheel-to-$\zeta$ convention-sensitive (2026-04-22)

**Trigger**: any Vol II identification of an $n$-spoke wheel
Feynman / Kontsevich graph with a specific zeta value without
naming the convention (Willwacher 2014 vs Brown--Schnetz
variant).

**Ghost**: Willwacher 2014 *Invent Math* 200 Thm 1.2 proves
$\mathrm{GRT}_1(\mathbb Q)$-freedom on Kontsevich graph cohomology
generated by $n$-spoke wheels $W_n$; in Willwacher's convention,
$W_n \leftrightarrow \zeta(2n - 1)$ (odd zeta value at weight
$2n - 1$). The Brown--Schnetz variant (Brown 2013 *Ann Sci ENS*
46; Schnetz 2014 *Commun Num Theor Phys* 8) uses
$W_n \leftrightarrow \zeta(2n + 1)$ with shifted indexing. Both
are correct; the distinction is convention on the wheel index.

**Precise error**: sites
\texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:366-368}
and
\texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:541}
identify a 3-spoke wheel with $\zeta(5)$ without naming the
convention. Under Willwacher 2014 a 3-spoke wheel $W_3$ couples
to $\zeta(2 \cdot 3 - 1) = \zeta(5)$; under Brown--Schnetz a
3-spoke wheel $W_3$ couples to $\zeta(2 \cdot 3 + 1) = \zeta(7)$.
The unnamed convention leaves the identification ambiguous.

**Correction**: state the convention at first use. Vol II default
is Willwacher 2014: $W_n \leftrightarrow \zeta(2n - 1)$; so the
3-spoke wheel $\leftrightarrow \zeta(5)$, the 4-spoke wheel
$\leftrightarrow \zeta(7)$, etc. Flag any cross-reference to
Brown--Schnetz with the shifted convention explicitly.

**Severity**: Medium. The convention is a named choice; fixing it
once at the top of the amplitude chapter suffices. The ambiguity
is not a falsification, but it does propagate into every downstream
$\zeta$-coefficient.

**Primary**: Willwacher 2014 *Invent Math* 200 Thm 1.2
"$M$. Kontsevich's graph complex and the Grothendieck--Teichmüller
Lie algebra"; Brown 2013 *Ann Sci ENS* 46 (single-valued periods);
Schnetz 2014 *Commun Num Theor Phys* 8.

**Cross-ref**: \texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:366-368, 541}
(error sites); AP-V2-49 / V2-AP152 ($E_1$-per-$\mathbb C$ Dunn scope,
same Willwacher framework); cache row 138 ($\mathrm{GRT}_1(\mathbb Q)$-torsor
chain-level $E_2$).

### AP-V2-51 / V2-AP154: K3-type MHS pollution — Feynman period on $K3 \times \mathbb C^2$ is not motivic MZV (2026-04-22)

**Trigger**: any Vol II assertion that a Feynman period integral on
$K3 \times \mathbb C^2$ lies in the motivic MZV algebra
$\mathrm{MZV}^{\mathrm{mot}}$.

**Ghost**: Brown 2012 *Ann Math* 175 proves that periods of
rational K-theory classes on mixed-Tate motives over $\mathbb Z$
generate $\mathrm{MZV}^{\mathrm{mot}}$. K3 surfaces are NOT mixed
Tate: the K3 period lattice carries a K3-type MHS with genuine
$H^{2,0}(K3)$ non-Tate classes (Kulikov 1977 *Izv Akad Nauk* 41
structure theorem). Periods of forms supported on $K3 \times X$
for any $X$ inherit the K3-type pollution and escape the
motivic MZV algebra.

**Precise error**: site
\texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4132}
asserts a Feynman period on $K3 \times \mathbb C^2$ lies in
$\mathrm{MZV}^{\mathrm{mot}}$. The K3 holomorphic 2-form
$\omega_{K3} \in H^{2,0}(K3)$ contributes a transcendental period
(the K3 Hodge structure) that is independent of MZVs; the total
period lattice is the tensor product $\mathrm{MHS}_{K3} \otimes
\mathrm{MHS}_{\mathbb C^2}$ and only the $\mathbb C^2$ factor is
mixed Tate.

**Correction**: Feynman periods on $K3 \times X$ lie in
$\mathrm{MHS}_{K3} \otimes \mathrm{MZV}^{\mathrm{mot}}(X)$ if $X$
is mixed Tate (e.g., $X = \mathbb C^n$ with logarithmic forms).
The K3-type MHS component is a separate Hodge-theoretic invariant
and must be named distinctly. Identifications
"$\text{period} \in \mathrm{MZV}^{\mathrm{mot}}$" hold only after
explicit restriction to a mixed-Tate subquotient of the K3 period
sheaf (typically via a Humbert-divisor stratification; see
AP-V2-6 for the Humbert / Kontsevich interplay).

**Severity**: High. Mixed-Tate vs K3-type MHS is a categorical
distinction; conflation falsifies every motivic identification on
K3-containing geometries.

**Primary**: Brown 2012 *Ann Math* 175 Thm 1 (motivic MZV =
mixed-Tate over $\mathbb Z$); Kulikov 1977 *Izv Akad Nauk* 41
(K3 degeneration / MHS structure); Deligne 1971 Hodge II (MHS
foundations); Griffiths 1969 (K3 period map).

**Cross-ref**: \texttt{chapters/connections/thqg\_perturbative\_finiteness.tex:4132}
(error site); AP-V2-6 (Humbert strata / chiral Kontsevich scope);
AP-V2-12 (Kaledin NC Hodge scope on Humbert); cache row 136
($\zeta^{\mathrm{sv}}$ vs motivic scope).

### AP-V2-52 / V2-AP155: AFT 2017 miscited for smooth-product factorisation homology (2026-04-22)

**Trigger**: any Vol II citation of Ayala--Francis--Tanaka 2017
Thm 4.1.3 in the context of factorisation homology of smooth
products.

**Ghost**: Ayala--Francis 2015 *J Topol* 8 Thm 3.16 establishes
factorisation-homology Fubini for smooth-product manifolds
$\int_{M \times N} \mathcal F = \int_M \int_N \mathcal F$ with
$\mathcal F$ a locally constant factorisation algebra on
$\mathrm{Disk}(M \times N)$. Ayala--Francis--Tanaka 2017
arXiv:1312.3745 extends to stratified spaces with
stratified-locally-constant sheaves; Thm 4.1.3 of AFT 2017 is the
stratified Fubini, not the smooth one.

**Precise error**: sites
\texttt{chapters/connections/celestial\_holography.tex:2496, 2517}
cite AFT 2017 Thm 4.1.3 for smooth-product factorisation Fubini.
The correct citation is AF 2015 Thm 3.16 (smooth); AFT 2017 is
about stratified spaces and its Thm 4.1.3 covers
\emph{stratified} Fubini with a different hypothesis set
(conical smoothness, constructible).

**Correction**: replace "AFT 2017 Thm 4.1.3" with "AF 2015 Thm 3.16"
for smooth-product applications. Use AFT 2017 Thm 4.1.3 only when
the ambient is genuinely stratified (e.g., compactified Heegner
divisor strata on $\overline{\mathcal A_2}$, AP-V2-6).

**Severity**: Medium. Attribution error; the smooth Fubini is a
true theorem but is due to Ayala--Francis alone. Mixed citation
blurs the stratified / smooth distinction, which matters for the
Humbert-stratum content of Vol II.

**Primary**: Ayala--Francis 2015 *J Topol* 8 Thm 3.16 "Factorization
homology of topological manifolds" (smooth Fubini); Ayala--Francis--Tanaka
2017 arXiv:1312.3745 Thm 4.1.3 (stratified Fubini).

**Cross-ref**: \texttt{chapters/connections/celestial\_holography.tex:2496, 2517}
(error sites); AP-V2-6 (Kontsevich formality on Humbert strata —
stratified context where AFT 2017 does apply); AP-V2-12 (Kaledin
stratified Hodge scope).

### AP-V2-53 / V2-AP156: Costello 2011 Thm 13.4.1 applied non-compact (2026-04-22)

**Trigger**: any Vol II citation of Costello 2011 *Renormalization
and Effective Field Theory* Thm 13.4.1 for a non-compact total
space such as $\mathbb R^3 \times K3 \times \mathbb C^2$.

**Ghost**: Costello 2011 Thm 13.4.1 (BV quantisation existence /
obstruction) is proved for compact-support observables OR for
theories with compactly supported fields on a compact manifold. On
non-compact total spaces the theorem requires a compact-support
hypothesis or a modified version; Costello--Gwilliam 2021 *FA*
Vol 2 Ch 4 extends with the modified hypothesis.

**Precise error**: sites
\texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:145, 433}
invoke Costello 2011 Thm 13.4.1 on the non-compact total space
$\mathbb R^3 \times K3 \times \mathbb C^2$ (or equivalent
non-compact 6d hCS setting). Without a compact-support hypothesis
or the CG 2021 Vol 2 Ch 4 modified statement, the theorem does
not apply.

**Correction**: either (i) restrict observables to compact support
and invoke Costello 2011 Thm 13.4.1 on the compact-support
subalgebra, or (ii) cite Costello--Gwilliam 2021 *FA* Vol 2 Ch 4
(the modified BV quantisation theorem on non-compact spaces with
appropriate IR regularisation). State the chosen hypothesis
explicitly.

**Severity**: Medium. The theorem is true in the intended
application; the error is a hypothesis-verification gap, fixable
by citing the correct version.

**Primary**: Costello 2011 *Renormalization and Effective Field
Theory* Thm 13.4.1; Costello--Gwilliam 2021 *Factorization
Algebras Vol 2* Ch 4 (modified BV quantisation on non-compact
total spaces).

**Cross-ref**: \texttt{chapters/theory/six\_d\_hcs\_feynman\_coefficients.tex:145, 433}
(error sites); AP-V2-42 / V2-AP145 (quartic Casimir in 6d hCS
amplitude); cache row 138 (chain-level $E_2$ GRT-torsor / BV
quantisation).

### AP-V2-54 / V2-AP157: CANONICAL-ANOM-LOCUS (form c) — $E_6$ strict exclusion, $A_2$ refined/unrefined distinguished (2026-04-22, canonical form propagated)

**Trigger**: any Vol II listing of the 6d hCS anomaly-free locus
that either (a) strips $\{E_6, A_2\}$ from Deligne without
distinguishing $A_2$-refined / $A_2$-unrefined and without the
$K^{-1/2}$-refinement clause, or (b) strips $\{E_6\}$ alone while
admitting $A_2$ without qualifier.

**Ghost**: the 6d hCS quartic Casimir
$d^{abcd} = \mathrm{tr}_{\mathrm{adj}}(T^{(a}T^b T^c T^{d)})$
vanishes universally on the Deligne series via
$\mathrm{tr}_{\mathrm{adj}}T^4 = \alpha_{\mathfrak g}
(\mathrm{tr}_{\mathrm{adj}}T^2)^2$ factorisation. A separate cubic
obstruction $d^{abc} = \mathrm{tr}_{\mathrm{adj}}(T^{(a}T^b T^{c)})$
arises for $E_6$ (cubic Jordan invariant on $\mathrm{Sym}^3(\mathbf{27})$,
non-zero, equal to the Jordan algebra product on
$\mathfrak j_3^{\mathbb O}$) and for $A_2 = \mathfrak{su}(3)$
(Gell-Mann $d$-tensor). Two distinct obstructions operate.

**Wrong claims flagged**:
- \textbf{Form (a) strict}: anomaly-free $\iff$ Deligne $\setminus
  \{E_6, A_2\}$ with no refinement clause; this misses that
  $A_2$-refined (with critical twist + Dimofte slab) sits INSIDE
  the locus.
- \textbf{Form (b)}: anomaly-free $\iff$ Deligne $\setminus \{E_6\}$
  alone, admitting undifferentiated $A_2$; this wrongly admits
  $A_2$-unrefined with live cubic $d^{abc}$ and live critical-level
  quadratic obstruction.

**Precise error**: both (a) and (b) misrepresent the native-ambient
locus. Form (a) is too coarse on $A_2$ (lumps refined with
unrefined); form (b) is too permissive (admits $A_2$-unrefined).
The canonical native form distinguishes four contributions and
carries the $K^{-1/2}$-refinement clause explicitly.

**Correction — CANONICAL-ANOM-LOCUS (form c)**: the native-ambient
anomaly-free locus reads
$$\mathrm{Anom}_1 = 0 \iff \mathfrak g \in
\bigl(\mathrm{Deligne}^{\mathrm{exc}} \setminus
\{E_6,\, A_2\text{-unrefined}\}\bigr) \cup \{\mathrm{abelian}\}
\cup \{\mathrm{super-str}_{\mathrm{ad}} = 0\}
\cup \{\widehat{\mathfrak g}_{-h^\vee} \otimes K^{-1/2}
\text{-refined}\}.$$
Native-ambient distinctions:
- $E_6$ STRICTLY excluded. No refinement in the programme's
  toolkit kills $\mathrm{Sym}^3(\mathbf{27})$ cubic Jordan invariant
  $d^{abc} \ne 0$ within native ambient; the $K^{-1/2}$ critical
  twist addresses the quadratic obstruction, not the cubic.
- $A_2$-unrefined excluded: live Gell-Mann $d_{abc}$ and live
  critical-level quadratic obstruction.
- $A_2$-refined INSIDE the locus: Feigin--Frenkel critical-level
  $K^{-1/2}$ twist kills the quadratic; Dimofte-slab
  anomaly-inflow from Vol II Part V provides Green--Schwarz cubic
  cancellation.
- $\{A_1, G_2, D_4, F_4, E_7, E_8\}$ unconditionally inside
  (quartic factorises via Deligne; cubic $d^{abc}$ vanishes
  identically on these six).
Two obstructions distinguished: quartic adjoint Casimir
(Deligne-killed across the series via factorisation,
including for $A_2, E_6$) vs cubic $d^{abc}$ (nonzero for $A_2, E_6$;
cured only in native ambient by Green--Schwarz-type inflow,
operative for $A_2$-refined, not for $E_6$).

**Regex trigger**: flag BOTH
- \verb|Deligne.*\\setminus.*\\\{E_6,\\s*A_2\\\}(?!.*refined)|
  (form a: strict $\{E_6, A_2\}$ without refinement);
- \verb|Deligne.*\\setminus.*\\\{E_6\\\}(?!.*A_2)|
  (form b: exclude-only-$E_6$);
relative to canonical (c) form which includes
``$A_2\text{-unrefined}$'' qualifier and the explicit
$K^{-1/2}$-refined clause.

**Severity**: Critical. Both (a) and (b) forms propagate into
downstream 3D HT QFT anomaly arguments; form (a) wrongly excludes
a live sector ($A_2$-refined via Dimofte slab), form (b) wrongly
admits a dead sector ($A_2$-unrefined with live cubic).

**Primary**: Deligne 1996 *C R Acad Sci Paris* 322 "La série
exceptionnelle de groupes de Lie"; Cvitanović 2008 *Group Theory*
Ch 20 ($E_6$ cubic invariant); Baez 2002 *Bull Amer Math Soc* 39
(Jordan algebra structure on $\mathfrak j_3^{\mathbb O}$);
Feigin--Frenkel 1992 *Comm Math Phys* 147 (critical-level
$K^{-1/2}$ twist); Dimofte 2014 slab anomaly-inflow.

**Cross-ref**: wherever the Deligne exceptional series is invoked
for anomaly cancellation in Vol II; AP-V2-42 / V2-AP145 (quartic
Casimir in 6d hCS 1-loop amplitude); AP-V2-53 / V2-AP156
(non-compact hCS hypothesis); Vol I AP979 (canonical-ambient
twin); AP-CY262 (Vol III cubic vs quadratic Casimir);
AP-CY50-E14 (cross-volume ledger).

### AP-V2-55 / V2-AP158: Conway $V^{s\natural}$ as bosonic $\Psi$-sibling (2026-04-22)

**Trigger**: any Vol II placement of Conway $V^{s\natural}$ as
the fifth row of the bosonic $\Psi$-sibling table, or as an
independent image of the bosonic $\Psi$-functor.

**Ghost**: the $\Psi$-functor has four sibling branches
$\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$
covering the Baily--Borel--Freitag stratification of
$\overline{\mathcal A_2}$; Conway $V^{s\natural}$ is a super
chiral algebra at $c = 12$ on the metaplectic branch
$\Psi^{\mathrm{metap}}$ via Duncan 2007 *Duke* 139 (Leech
fermionic orbifold $A(\Lambda_{24})^+$, $c = 12$, NOT $c = 24$).
It is NOT a bosonic $\Psi$-image.

**Precise error**: any Vol II site placing Conway $V^{s\natural}$
as a bosonic $\Psi$-image at $c = 24$ or as the fifth independent
$\Psi$-sibling row. The correct placement is on the
$\Psi^{\mathrm{metap}}$ branch at $c = 12$ via the metaplectic
extension; the Leech lattice has no hyperbolic plane, so the
bosonic $\Psi$-image universal identity $(K, \hbar^2) = (2, -1/2)$
(Canonical preamble row 41) does not apply.

**Correction**: place Conway $V^{s\natural}$ on the
$\Psi^{\mathrm{metap}}$ branch at $c = 12$ via Duncan 2007 Leech
fermionic orbifolding $A(\Lambda_{24})^+$. Four siblings of the
$\Psi$-functor cover the Baily--Borel--Freitag stratification:
$\Psi$ (principal), $\Psi^{\deg}$ (cuspidal), $\Psi^{\mathrm{tor}}$
(toric), $\Psi^{\mathrm{metap}}$ (metaplectic); Conway is on the
fourth branch.

**Severity**: High. The Conway placement on the metaplectic branch
is load-bearing for the Part VI 3d gravity climax (the Conway
module / elliptic genus feeds Mathieu moonshine and the paramodular
/ Humbert identification).

**Primary**: Duncan 2007 *Duke Math J* 139 "Super-moonshine for
Conway's largest sporadic group"; Scheithauer 2008 *Adv Math* 220
(Freitag $\Psi$-functor sibling classification); Canonical preamble
rows 40-42 (Conway placement pinned).

**Cross-ref**: Canonical preamble rows 40-42 (Conway $c = 12$ on
$\Psi^{\mathrm{metap}}$); AP-V2-22 (five-archetype B-row);
V2-AP121 (Conway legacy placement AP).

### AP-V2-56 / V2-AP159: CHL Siegel weight ladder $\{\kappa(\Phi_N)\} = \{5, 4, 3, 2, 1\}$ false (2026-04-22)

**Trigger**: any Vol II identification of the CHL (Chaudhuri--Hockney--Lykken)
Siegel-weight ladder as $\{5, 4, 3, 2, 1\}$.

**Ghost**: the CHL construction orbifolds heterotic on
$T^6/\mathbb Z_N$ with $N \in \{2, 3, 5, 7, 11\}$; the resulting
1/4-BPS partition function is $1/\Phi_{k_N}$ with
$k_N = 24/(N+1) - 2$ (David--Jatkar--Sen 2006 JHEP 04:018 Eq
3.14). So $(k_2, k_3, k_5, k_7, k_{11}) = (6, 4, 2, 1, 0)$ in
the DJS indexing; with $N = 1$ (no orbifold, heterotic on
$T^6$) one gets $k_1 = 10$ (Igusa $\Phi_{10}$); including
$N = 1$ the sequence is $(10, 6, 4, 3, 2)$ for
$N \in \{1, 2, 3, 5, 7\}$ (David--Jatkar--Sen Eq 3.14 with
$N = 1$ the Igusa; tables in Persson--Volpato 2012 Table 1).

**Precise error**: any Vol II site writing
$\{\kappa(\Phi_N)\} = \{5, 4, 3, 2, 1\}$ for the CHL ladder.
The correct weight ladder is
$k_N = 24/(N + 1) - 2$, giving
$(k_1, k_2, k_3, k_5, k_7) = (10, 6, 4, 3, 2)$ across the five
CHL orbifolds. The stated $\{5, 4, 3, 2, 1\}$ sequence matches
neither the DJS weight nor the Siegel lifting weight
$\kappa(\Phi)$.

**Correction**: CHL weight ladder is
$k_N = 24/(N + 1) - 2$; explicit values
$(k_1, k_2, k_3, k_5, k_7) = (10, 6, 4, 3, 2)$. The Igusa
$\Phi_{10}$ is $k_1$; the $\mathbb Z_2$ CHL gives $\Phi_6$; the
$\mathbb Z_3$ CHL gives $\Phi_4$; etc. State the formula, not the
numerical sequence, as the load-bearing relation.

**Severity**: Critical. The weight ladder anchors the
BKM-$\Phi_N$ family across the Vol II / Vol III bridge; the
wrong sequence propagates into every CHL-dependent identification.

**Primary**: David--Jatkar--Sen 2006 JHEP 04:018 Eq 3.14;
Persson--Volpato 2012 arXiv:1212.0471 Table 1 (CHL Siegel weights);
Chaudhuri--Hockney--Lykken 1995 (original CHL orbifolds).

**Cross-ref**: Canonical preamble rows 28, 32-33 ($\Phi_{10}$,
$\Phi_{12}$, $\Phi_N$ placements); AP-V2-57 / V2-AP160
($\kappa_{\mathrm{BKM}}$ scope 5 vs 12).

### AP-V2-57 / V2-AP160: $\kappa_{\mathrm{BKM}}$ scope 5 vs 12 — $\Psi$-sibling-indexed (2026-04-22)

**Trigger**: any Vol II bare assertion of $\kappa_{\mathrm{BKM}}
= 5$ or $\kappa_{\mathrm{BKM}} = 12$ without naming the
$\Psi$-sibling denominator.

**Ghost**: $\kappa_{\mathrm{BKM}}$ is the Borcherds-weight
invariant of the BKM algebra; its numerical value is indexed by
the $\Psi$-sibling ($\Phi_N$ denominator) driving the lift. K3-BKM
(Gritsenko--Nikulin paramodular, $\Delta_5$ input) gives
$\kappa_{\mathrm{BKM}}^{\mathrm{K3}} = 5$; Fake-Monster
(Borcherds 1992, $\Phi_{12}$ input) gives
$\kappa_{\mathrm{BKM}}^{\mathrm{FM}} = 12$. Implicit conflation
loses the $\Psi$-sibling index.

**Precise error**: any Vol II site asserting $\kappa_{\mathrm{BKM}}$
bare (no sub/superscript denominator); the two legitimate values
$5$ and $12$ correspond to two distinct $\Phi_N$ inputs. The
Delta/Phi/OP denominator discipline: both values occur in published
inscriptions under different input-denominator conventions.

**Correction**: every Vol II $\kappa_{\mathrm{BKM}}$ reference
carries a superscript / subscript identifying the $\Psi$-sibling
denominator:
$\kappa_{\mathrm{BKM}}^{\mathrm{K3}} = 5$ for
$\Delta_5$ (paramodular, Gritsenko--Nikulin 1998);
$\kappa_{\mathrm{BKM}}^{\mathrm{FM}} = 12$ for $\Phi_{12}$
(Fake-Monster, Borcherds 1992). State the input $\Phi_N$ on every
inscription.

**Severity**: Critical. Flagged AP5-pending in Canonical preamble
row 60; unresolved conflation silently propagates.

**Primary**: Gritsenko--Nikulin 1998 *J Reine Angew Math* 507
(paramodular $\Delta_5$); Borcherds 1992 *Invent Math* 109
(Fake-Monster $\Phi_{12}$); Canonical preamble row 60
($\kappa_{\mathrm{BKM}}(\mathbf H_{\Delta_5})$ dual-indexing AP5
header).

**Cross-ref**: Canonical preamble row 60 (pending AP5 lock);
cache row 123 (five-archetype B-row witness); AP-V2-56 / V2-AP159
(CHL weight ladder); Cache header AP5 dual-indexing preamble
(lines 5-13).

### AP-V2-58 / V2-AP161: Three Yangian variants — four types in programme (2026-04-22)

**Trigger**: any Vol II assertion of "three Yangian types in the
programme" without naming the four canonical variants.

**Ghost**: the programme distinguishes four Yangian types on
different ambient spaces with different primary structures
(feedback entry \texttt{feedback\_yangian\_type\_distinction.md}
in MEMORY.md):
(i) \emph{classical Yangian} $Y(\mathfrak g)$ on the spectral line
$\mathbb C_u$ (Drinfeld 1985 *Sov Math Dokl* 32);
(ii) \emph{chiral Yangian} on the curve $C$ with $E_1$-chiral
structure (Vol II standalone "On the Integrable Theory of Chiral
Yangians");
(iii) \emph{spectral Yangian} on the torus $E_\tau$ via KZB
(Costello--Witten--Yamazaki 2018);
(iv) \emph{dg-shifted affine Yangian} on $T^* \mathbb C$ with
BV shift (Kodera--Nakajima 2018 arXiv:1602.06499;
Feigin--Frenkel--Rybnikov 2010).

**Precise error**: any Vol II site stating "three Yangian types"
conflates either (i) and (iv) (classical and dg-shifted) or drops
(iii) (spectral). Correct count is four.

**Correction**: state four Yangian types: classical, chiral,
spectral, dg-shifted affine. Name each with its ambient space
and primary reference at first use per chapter.

**Severity**: Medium. The conflation is a scope drift; fixing it
once per chapter clarifies the ambient-space structure of each
variant.

**Primary**: Drinfeld 1985 *Sov Math Dokl* 32 (classical Yangian);
Costello--Witten--Yamazaki 2018 arXiv:1709.09993 (spectral
Yangian); Kodera--Nakajima 2018 arXiv:1602.06499 (dg-shifted
affine Yangian); Vol II standalone "On the Integrable Theory of
Chiral Yangians" (chiral Yangian).

**Cross-ref**: MEMORY.md
\texttt{feedback\_yangian\_type\_distinction.md}; Vol II standalone
chiral-Yangians paper (protagonist of the programme identity);
Canonical preamble rows on four-types-of-centre discipline.

### AP-V2-59 / V2-AP162: $\Phi_d$ Stage-1 / Stage-2 unnamed in 6d hCS context (2026-04-22)

**Trigger**: any Vol II reference to 6d hCS as "$E_3$ on
$\mathbb C^3$" without naming the two-stage factorisation
$\Phi_d = \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}
\circ \Phi^{\mathrm{FA}}_d$.

**Ghost**: the Vol III two-stage factorisation pins
$\Phi_d$ as Stage-1 (the canonical $E_d$-holomorphic factorisation
algebra on the CY-$d$ target via Kontsevich--Tamarkin formality
plus Costello--Gwilliam--Li locality) followed by Stage-2 (the
factorisation-homology pushforward $\int_{\Sigma_{d-1}}$ and
restriction to a reference curve $C$; Vol II CLAUDE.md lines
281-318). 6d hCS on $\mathbb C^3$ fits this framework as
Stage-1 $\Phi^{\mathrm{FA}}_3$, not as an independent $E_3$
construction.

**Precise error**: any Vol II site referring to 6d hCS as
"$E_3$ on $\mathbb C^3$" without naming Stage-1 $\Phi^{\mathrm{FA}}_3$
and the pushforward structure. The unnamed-stage framing hides
the CY-side canonical origin and the $(\Sigma_{d-1}, C)$-family
parametrisation of the $E_1$-chiral shadows on $C$.

**Correction**: every Vol II reference to 6d hCS (or $E_3$
$\mathbb C^3$-factorisation algebra content) names:
Stage-1 $\Phi^{\mathrm{FA}}_3$ as the canonical $E_3$-holomorphic
factorisation algebra on the CY-$3$ target; Stage-2
$\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_2, C}$ as the
$\int_{\Sigma_2}$-pushforward to the reference curve $C$; the
$E_1$-chiral shadow on $C$ is the output, parametrised by
$(\Sigma_2, C)$.

**Severity**: High. Part of the 2026-04-22 two-stage factorisation
wave (AP-V2-25..AP-V2-29); unnamed stages break the cross-volume
$\Phi_d$ framework.

**Primary**: Vol II CLAUDE.md two-stage factorisation section
(lines 281-318); Vol III
\texttt{chapters/theory/cy\_to\_chiral.tex}; Costello 2013
arXiv:1303.2632; Costello--Li 2016 arXiv:1605.09930;
Kontsevich 2003 *Lett Math Phys* 66; Tamarkin 2003
*Lett Math Phys* 66.

**Cross-ref**: AP-V2-25 / V2-AP128 (closed/open colour =
Stage-1 / Stage-2); AP-V2-26 / V2-AP129 (single-stage $\Phi_d$
framing); AP-V2-27 / V2-AP130 (3D HT QFT two-stage anchor);
AP-V2-28 / V2-AP131 (manifesto conflations from Vol III);
AP-V2-49 / V2-AP152 ($E_1$-per-$\mathbb C$ Dunn / Dolbeault).


### 2026-04-22: 15 cross-volume AP entries added from Vol III Wave 11-19

Fifteen cross-volume entries mirror the Vol III Wave 11-19 error/retraction
catalogue (six errors E\_A--E\_F, nine retractions R1--R9) into the Vol II
registry with Vol II-specific manifestations on the Swiss-cheese two-colour
structure, 3D HT QFT framework, class-$\mathcal S$ / AGT / Gaiotto curve
chapters, Pentagon single-colour coherence, and modular PVA quantisation.
AP tags AP-V2-60 / V2-AP163 through AP-V2-74 / V2-AP177 follow after the
AP-V2-59 / V2-AP162 inscription. Each entry shares a signature with the
corresponding Vol I AP (AP939--AP953) and the Vol III AP-CY primary form.

### AP-V2-60 / V2-AP163: Künneth multiplicative on products; partition function factorisation (2026-04-22)
**Severity**: High (load-bearing for class-$\mathcal S$ partition
function factorisation)
**Trigger**: Vol II partition-function factorisation on a product
CY ambient (e.g., $K3 \times E$ Swiss-cheese pairing, $K3 \times T^2$
class-$\mathcal S$ reduction) cites
$\kappa_{\mathrm{cat}}(K3 \times E) = 2$ or treats the fibre-K3 value
$\chi(\cO_{K3}) = 2$ as the total-space invariant. The partition-
function factorisation $Z(K3 \times E) = Z(K3) \cdot Z(E)$ at the
$\mathsf{SC}^{\mathrm{ch,top}}$ Swiss-cheese level requires the
Künneth-multiplicative $\kappa_{\mathrm{cat}}$ rule; additive or
fibre-only framing breaks the factorisation and inflates the one-loop
Quillen exponent.
**Ghost**: categorical Euler on a product is Künneth-multiplicative:
$\kappa_{\mathrm{cat}}(X \times Y) = \chi(\cO_X) \cdot \chi(\cO_Y)$,
giving $\kappa_{\mathrm{cat}}(K3 \times E) = 2 \cdot 0 = 0$. This is
load-bearing for the Vol II Swiss-cheese one-loop Quillen exponent
$\exp(2\pi i Q_{\mathrm{Quillen}})$ and for class-$\mathcal S$
partition functions over product compactifications.
**Precise error**: bare $\kappa$ on a product loses the Künneth-
multiplicative structure; the fibre-K3 value $2$ is presented as the
total-space value; Vol II partition-function factorisation breaks.
**Correction**: every Vol II site citing $\kappa$ on a product CY
subscripts: $\kappa_{\mathrm{cat}}(K3 \times E) = 0$ (total space,
Künneth-multiplicative), $\kappa_{\mathrm{fiber}}(K3) = 2$ (fibre
only). Partition-function factorisation reads $Z(K3 \times E) =
Z(K3) \cdot Z(E)$ with $Z(E) = 0$-level consistent with $\chi(\cO_E) = 0$.
**Primary**: Vol III
\texttt{chapters/examples/cy\_d\_kappa\_stratification.tex} canonical
table; canonical preamble row 66; Vol III AP-CY primary form (E\_A).
**Inscription**: Vol II
\texttt{chapters/theory/sc\_chtop\_partition\_factorisation.tex};
cross-volume $\kappa_{\mathrm{cat}}$, $\kappa_{\mathrm{fiber}}$
subscripts.
**Cross-ref**: Vol III AP-CY (E\_A primary form); Vol I AP939
(partner); Vol II canonical preamble rows 59, 66; AP-V2-28 /
V2-AP131 (manifesto conflation); V2-AP68 / V2-AP69.

### AP-V2-61 / V2-AP164: Central charge $-1312/11$ at $c_{4d} = 107/6$ via two independent verification routes (2026-04-22)
**Severity**: High (load-bearing for
$\mathcal T[A_1, \Sigma_{0, 24}]$ landscape anchor)
**Trigger**: Vol II class-$\mathcal S$ / AGT chapters citing the
$(c_{4d}, c_{2d})$ central-charge pair for
$\mathcal T[A_1, \Sigma_{0, 24}]$ use the retracted Wave 14 value
$(26, -312)$ or the retracted Wave 25 Gaiotto-agent value
$(23/4, -69)$, or cite the canonical $(107/6, -214)$ without
two-route verification.
**Ghost**: canonical preamble row 17 pins $(c_{4d}, c_{2d}) = (107/6,
-214)$ for $\mathcal T[A_1, \Sigma_{0, 24}]$ via Wave 15 Gaiotto
verification against the universal formula $c_{4d} = (5n - 13)/6$
(canonical preamble row 19, Shapere--Tachikawa 2008).
**Precise error**: single-route citation of a load-bearing structure
constant whose Wave 14 value was retracted; AP274-style citation-as-
derivation.
**Correction**: every Vol II $(c_{4d}, c_{2d})$ citation for
$\mathcal T[A_1, \Sigma_{0, 24}]$ uses $(107/6, -214)$ with two
independent verification routes named: Shapere--Tachikawa 2008
universal formula $c_{4d} = (5n - 13)/6$ at $n = 24$; direct
Sugawara-denominator recomputation on the Gaiotto curve.
**Primary**: Vol III Wave 15 Gaiotto verdict; canonical preamble
rows 17, 19; Shapere--Tachikawa 2008.
**Inscription**: Vol II
\texttt{chapters/connections/class\_s\_gaiotto.tex} two-route
verification block.
**Cross-ref**: Vol III AP-CY (E\_B primary form); Vol I AP940
(partner); canonical preamble rows 17--19; AP-V2-19 (universal
ratio-of-levels).

### AP-V2-62 / V2-AP165: $\eta^{-48}$ Heisenberg--Mukai identity in modular PVA quantisation chapters, not Virasoro minimal-model (2026-04-22)
**Severity**: Medium (modular PVA quantisation anchor)
**Trigger**: Vol II modular PVA quantisation chapters
(\texttt{chapters/theory/pva-expanded.tex},
\texttt{chapters/connections/modular\_pva\_quantisation.tex}) frame
the $\eta^{-48}$ factor in $V_{24}$ character as a Virasoro minimal-
model expansion; the correct reading is the Heisenberg--Mukai lattice
partition function $\eta^{-48} = \prod_n (1 - q^n)^{-24}$ on the Mukai
lattice $\mathrm{II}_{1, 1} \otimes \cO_{K3}$.
**Ghost**: iterated Drinfeld--Sokolov reduction $V_{24} =
\mathrm{DS}^{24}(\mathrm{Heis}_{\mathrm{II}_{1,1}})$ produces
$\eta^{-48}$; character identity $\chi_{V_{24}} = \eta^{-48}$ holds by
iterated-DS contraction + Borcherds lift on the Mukai lattice. The
Virasoro minimal-model reading would give a finite-dimensional
character with no $\eta$-denominator.
**Precise error**: conflates two distinct vertex-algebra presentations
of the same generating function; Heisenberg--Mukai has a free-field
factor, Virasoro minimal-model does not.
**Correction**: every Vol II $\eta^{-48}$ reference names the
Heisenberg--Mukai lattice VOA and the iterated Drinfeld--Sokolov chain
$V_{24} = \mathrm{DS}^{24}(\mathrm{Heis}_{\mathrm{II}_{1,1}})$;
the modular PVA quantisation lands on the Heisenberg--Mukai positive-
half vertex-operator trace.
**Primary**: Vol III Wave 15--17 three-path verification;
Mukai 1988; Frenkel--Kac construction; canonical preamble row 52.
**Inscription**: Vol II \texttt{chapters/theory/pva-expanded.tex}
$V_{24}$ paragraph;
\texttt{chapters/connections/modular\_pva\_quantisation.tex}.
**Cross-ref**: Vol III AP-CY (E\_C, R5 primary form); Vol I AP941
(partner); AP-V2-49 / V2-AP152 ($V_{24}$ iterated-DS partner);
canonical preamble row 52 (chiral-Hochschild $e_k$ motivic home).

### AP-V2-63 / V2-AP166: Macdonald framework inapplicable to non-minimal Pentagon structures (2026-04-22)
**Severity**: Medium (Pentagon single-colour coherence scope)
**Trigger**: Vol II Pentagon single-colour coherence chapters invoke
the Macdonald difference-operator framework at non-minimal targets
(W-algebras at non-admissible level, logarithmic triplet VOAs,
parafermion cosets) Pentagon vertices.
**Ghost**: Macdonald framework is theorem for Cherednik DAHA and
admissible-level affine Lie algebras with minimal-model /
integrable structure. Pentagon single-colour coherence at admissible
level reduces to Macdonald-polynomial diagonalisation; at non-
admissible level the positive-part Yangian is not diagonal in the
Macdonald basis.
**Precise error**: Macdonald framework extended beyond its
admissibility domain; Pentagon coboundary $\phi^{(n)}$ silently
relies on admissibility.
**Correction**: every Vol II Pentagon / Macdonald invocation names
the admissibility constraint — minimal target, admissible level
$k + h^\vee = p / q$ with $(p, q) = 1$, finite-dimensional bar at
each weight. Non-minimal Pentagon uses Kac--Roan 2014 integrable-
system residue or Feigin--Frenkel centre on $\mathrm{Op}$-side.
**Primary**: Vol III Wave 18 Gaiotto verdict on admissibility;
Macdonald 1995; Cherednik 2005; canonical preamble row 44.
**Inscription**: Vol II
\texttt{chapters/theory/pentagon\_single\_colour.tex} admissibility-
scope remark.
**Cross-ref**: Vol III AP-CY (E\_D primary form); Vol I AP942
(partner); canonical preamble rows 44--48.

### AP-V2-64 / V2-AP167: Eight-form position $\ne$ Borcherds weight in Swiss-cheese Borcherds-lift inscriptions (2026-04-22)
**Severity**: Medium (Swiss-cheese Borcherds-lift catalogue)
**Trigger**: Vol II Swiss-cheese Borcherds-lift chapters use position-
in-catalogue as stand-in for Borcherds weight; weights are non-
monotone: $(w_1, \ldots, w_8) = (5, 2, 1, 1, 1/2, 1, 1/4, 0)$ with
Fourier coefficients $c_N(0) \in \{10, 4, 2, 2, 1, 2, 1/2, 0\}$.
**Ghost**: eight-form catalogue is enumerated by $\mathrm{Sp}_4(\Z)$-
cover refinement — integral ($\mathrm{Sp}_4(\Z)$), half-integral
($\mathrm{Mp}_4$), quarter-integral ($\widetilde{\mathrm{Mp}}_4$) — not
by weight; the weight-0 form is the degenerate terminal fibre.
**Precise error**: enumeration index confused with mathematical
invariant; "the $k$-th form" without weight specification loses
$c_N(0)/2$ identification.
**Correction**: every Vol II citation of the eight-form catalogue
names the cover group and the weight adjacent to each form; universal
identity $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$ stated with the
$(w_k, c_k(0))$ pair.
**Primary**: Gritsenko--Clery 2018; Vol III
\texttt{chapters/examples/cy\_d\_kappa\_stratification.tex} Theorem
\texttt{thm:borcherds-weight-kappa-BKM-universal}.
**Inscription**: Vol II
\texttt{chapters/examples/sc\_chtop\_borcherds\_lift.tex} eight-form
$(w, c(0))$-pairs.
**Cross-ref**: Vol III AP-CY (E\_E primary form); Vol I AP943
(partner); canonical preamble row 28 (K3-BKM Weyl denominator);
AP-V2-62 / V2-AP165 (Heisenberg--Mukai $\eta^{-48}$).

### AP-V2-65 / V2-AP168: Twined convention $c_2(0) \in \{4, 8\}$ declared at every Swiss-cheese cross-citation with Vol III (2026-04-22)
**Severity**: Medium (cross-volume convention discipline)
**Trigger**: Vol II Swiss-cheese cross-citations to Vol III Borcherds-
weight Fourier coefficients $c_N(0)$ silently switch between Gritsenko
1999 ($c_2(0) = 4$) and Cheng--Duncan 2013 umbral ($c_2(0) = 8$).
**Ghost**: the two conventions differ by a normalisation factor of $2$
at $N = 2$; both are published, both self-consistent. Vol III uses
Gritsenko; Vol II cross-citations not declaring convention risk silent
switch.
**Precise error**: convention-dependent structure constant cited bare;
ambient-qualifier drop (AP281 partner).
**Correction**: every Vol II reference to twined $c_N(0)$ names the
Borcherds automorphic form and the convention: "Gritsenko 1999
twined coefficients ($c_2(0) = 4$)" or "Cheng--Duncan 2013 umbral
coefficients ($c_2(0) = 8$)". Vol II defaults to Gritsenko for
cross-citations to Vol III.
**Primary**: Gritsenko 1999 Thm 6.1; Cheng--Duncan 2013.
**Inscription**: Vol II
\texttt{chapters/theory/sc\_chtop\_borcherds\_lift.tex} convention
declaration.
**Cross-ref**: Vol III AP-CY (E\_F primary form); Vol I AP944
(partner); AP281 partner (ambient-qualifier drop); canonical
preamble row 28.

### AP-V2-66 / V2-AP169: Feingold--Frenkel $F_3$ real-root subalgebra at Vol II Swiss-cheese genus-1 BKM sites (2026-04-22)
**Severity**: High (genus-1 BKM real-root identification)
**Trigger**: Vol II Swiss-cheese genus-1 BKM chapters frame the real-
root subalgebra of $\mathfrak g_{\Delta_5}$ as
$\widehat{\mathfrak{sl}_3}$; correct identification is the Feingold--
Frenkel rank-3 hyperbolic Kac--Moody algebra $F_3$.
**Ghost**: $F_3$ (Feingold--Frenkel 1983) has rank-3 hyperbolic Cartan
matrix $\begin{pmatrix} 2 & -2 & 0 \\ -2 & 2 & -1 \\ 0 & -1 & 2
\end{pmatrix}$, indefinite type, imaginary-root Weyl chamber non-
empty. Affine $\widehat{\mathfrak{sl}_3}$ has rank-2 Cartan plus
central element, no imaginary roots, cannot carry Humbert-divisor
structure of K3-BKM genus-1 anchor.
**Precise error**: conflates affine Lie algebra with rank-3 hyperbolic
Kac--Moody algebra at real-root subalgebra site.
**Correction**: every Vol II genus-1 BKM real-root reference cites
Feingold--Frenkel 1983 $F_3$ with rank-3 hyperbolic Cartan matrix
and names Humbert-divisor matching on K3 fibre.
**Primary**: Feingold--Frenkel 1983 *Math Ann* 263; Gritsenko--
Nikulin 1998 Table 2.
**Inscription**: Vol II
\texttt{chapters/theory/sc\_chtop\_genus1.tex} real-root
identification.
**Cross-ref**: Vol III AP-CY (R1 primary form); Vol I AP945
(partner); canonical preamble row 28; AP-V2-62 / V2-AP165
(Heisenberg--Mukai).

### AP-V2-67 / V2-AP170: $L_{-6}(\mathfrak{e}_8) \to$ iterated Drinfeld--Sokolov reduction for $V_{24}$ at 3D HT QFT bulk (2026-04-22)
**Severity**: High (3D HT QFT bulk $V_{24}$ presentation)
**Trigger**: Vol II 3D HT QFT bulk sites frame $V_{24}$ as
$L_{-6}(\mathfrak{e}_8)$; correct presentation is iterated Drinfeld--
Sokolov reduction from a higher-rank source.
**Ghost**: $V_{24} = \mathrm{DS}^{24}(\mathrm{Heis}_{\mathrm{II}_{1,1}})$
iterated-DS output with 24 singular iterations; $c(V_{24}) = 24$.
$L_{-6}(\mathfrak{e}_8)$ has $c = -992/11$, so identification fails
at central-charge level.
**Precise error**: conflates affine VOA at specific level with
iterated-DS output of higher-rank chain; central charges disagree.
**Correction**: every Vol II $V_{24}$ reference cites iterated-DS
chain with source, level, nilpotent orbit; central-charge match
$c(V_{24}) = 24$ is the anchor.
**Primary**: Vol III Wave 15 iterated-DS chain; Kac--Roan 2000;
Adamovic--Milas 1997.
**Inscription**: Vol II
\texttt{chapters/theory/ht\_qft\_bulk\_v24.tex} iterated-DS
presentation.
**Cross-ref**: Vol III AP-CY (R2 primary form); Vol I AP946
(partner); AP-V2-62 / V2-AP165 ($\eta^{-48}$); canonical preamble
row 52.

### AP-V2-68 / V2-AP171: $\kappa_{\mathrm{BKM}}$ not additive; universal $c_N(0)/2$ on HT top-form dimension (2026-04-22)
**Severity**: High ($\kappa_\bullet$ discipline)
**Trigger**: Vol II 3D HT QFT HT top-form dimension calculations
factorise $\kappa_{\mathrm{BKM}}$ additively as
$\kappa_{\mathrm{BKM}} = \kappa_{\mathrm{ch}} + \chi(\cO_{\mathrm{fiber}})$
across $N$-family.
**Ghost**: additive relation is FALSE at every
$N \in \{1, 2, 3, 4, 6\}$: at $N = 1$, left $= 5$, right $=
\kappa_{\mathrm{ch}}(K3 \times E) + \chi(\cO_E) = 0 + 0 = 0$; at
$N = 2$, left $= 4$, right $= 1$. Universal formula is
$\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$ (Borcherds weight).
**Precise error**: additive decomposition of Borcherds-lift Fourier
coefficient into CY-categorical and fibre-Euler pieces has no lift-
functoriality; Borcherds lift is multiplicative, not additive.
**Correction**: every Vol II $\kappa_{\mathrm{BKM}}$ identity reads
with scope-split — Borcherds-weight side uses multiplicative lift on
Mukai lattice; $\kappa_{\mathrm{ch}}$ / $\kappa_{\mathrm{fiber}}$
bookkeeping is CY-categorical denominator; the two scopes meet at
$\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$ as universal identity,
not additive decomposition.
**Primary**: Vol III
\texttt{chapters/examples/cy\_d\_kappa\_stratification.tex} Theorem
\texttt{thm:borcherds-weight-kappa-BKM-universal}; Borcherds 1995;
Gritsenko 1999.
**Inscription**: Vol II
\texttt{chapters/theory/ht\_qft\_kappa\_split.tex} scope-split
discipline.
**Cross-ref**: Vol III AP-CY (R3 primary form); Vol I AP947
(partner); canonical preamble row 59; AP-V2-60 / V2-AP163 (Künneth
multiplicative).

### AP-V2-69 / V2-AP172: Fake Monster at $d = 5$ on HT top-form dimension; $d = 5$ Poisson-$E_5$ corner closure (2026-04-22)
**Severity**: High (CY-dimension discipline)
**Trigger**: Vol II HT top-form dimension cites Fake Monster BKM as
$d = 3$ object or omits $d = 5$ Poisson-$E_5$ corner closure from
PTVV / HT top-form table.
**Ghost**: Fake Monster Cartan rank $26$ on $\mathrm{II}_{25, 1}$
(Borcherds 1992), natural CY home $K3 \times K3 \times E$ at $d = 5$
with $\Phi_{12}$ Borcherds lift; $d = 5$ PTVV shift $-4$, Poisson-
$E_5$ corner closure via Dunn--Lurie $E_3 \otimes E_2 = E_5$. K3-BKM
(Cartan rank $3$ on $\Lambda^{2, 1}_{\mathrm{II}}$) is $d = 3$ sibling.
**Precise error**: CY-dimension confusion Fake Monster ($d = 5$,
$\Phi_{12}$) vs K3-BKM ($d = 3$, $\Phi_{10} = \Delta_5^2$); HT top-
form table truncated at $d = 4$.
**Correction**: every Vol II HT top-form / PTVV / Fake-Monster
reference names $d$, lattice signature, Borcherds-lift denominator:
$d = 5$ on $\mathrm{II}_{25, 1}$ with $\Phi_{12}$ for Fake Monster;
$d = 3$ on $\Lambda^{2, 1}_{\mathrm{II}}$ with $\Phi_{10}$ for K3-BKM;
HT top-form table extends to $d = 5$ Poisson-$E_5$.
**Primary**: Borcherds 1992 *Invent Math* 109; Pantev--Toen--Vaquie--
Vezzosi 2013 *Publ IHES* 117; Vol III Wave 18 PTVV extension;
canonical preamble rows 20--22.
**Inscription**: Vol II
\texttt{chapters/theory/ht\_top\_form\_table.tex} $d = 5$ row;
\texttt{chapters/theory/ptvv\_shifted\_symplectic.tex} Poisson-$E_5$.
**Cross-ref**: Vol III AP-CY (R4, R8 primary forms); Vol I AP948,
AP952 (partners); canonical preamble rows 20--22.

### AP-V2-70 / V2-AP173: Direct $\chi_{V_{24}} = \eta^{-48}$ match via three-path verification (2026-04-22)
**Severity**: Medium (character-match verification discipline)
**Trigger**: Vol II $V_{24}$ character inscription cites
$\chi_{V_{24}} = \eta^{-48}$ as direct character identity without
naming the three-path verification.
**Ghost**: three-path verification: (i) iterated-DS character
contraction from $L_k(\mathfrak{e}_8)$ source; (ii) Sugawara-
denominator recomputation on Heisenberg--Mukai lattice; (iii)
Borcherds lift on the Mukai lattice.
**Precise error**: load-bearing identity cited without derivation;
direct-match framing hides the three-path verification.
**Correction**: every Vol II site citing $\chi_{V_{24}} = \eta^{-48}$
names the three verification paths and the iterated-DS construction
of $V_{24}$.
**Primary**: Vol III Wave 15--17 three-path verification;
AP-V2-67 / V2-AP170 (iterated-DS).
**Inscription**: Vol II \texttt{chapters/theory/pva-expanded.tex} /
\texttt{chapters/theory/ht\_qft\_bulk\_v24.tex} three-path block.
**Cross-ref**: Vol III AP-CY (R5 primary form); Vol I AP949
(partner); AP-V2-62 / V2-AP165 ($\eta^{-48}$); AP-V2-67 / V2-AP170
(iterated DS).

### AP-V2-71 / V2-AP174: Gaiotto curve $\Sigma_{0, 24}$ with 24 punctures, not $\Sigma_{2, 0}$ (2026-04-22)
**Severity**: High (class-$\mathcal S$ / AGT curve identification)
**Trigger**: Vol II class-$\mathcal S$ / AGT / Gaiotto curve chapters
frame the Gaiotto curve for $\mathcal T[A_1, K3]$ or K3-BKM landscape
as $\Sigma_{2, 0}$ closed genus-2.
**Ghost**: correct Gaiotto curve is $\Sigma_{0, 24}$: genus-0 with
24 simple $A_1$ punctures (regular / $\mathrm{SU}(2)$-flavour).
Puncture count $24$ matches $\chi(K3) = 24$ and the $24 I_1$ nodes
of the Pentagon construction. Trinion/tube count $(n_v, n_h) =
(63, 88)$ (canonical preamble row 18, Chacaltana--Distler 2010).
**Precise error**: puncture-count / genus confusion; moduli dimension
$6g - 6 + 2n = 42$ for $\Sigma_{0, 24}$ vs $6$ for $\Sigma_{2, 0}$.
**Correction**: every Vol II Gaiotto-curve reference uses
$\Sigma_{0, 24}$ with 24 $A_1$ punctures; trinion/tube count
$(63, 88)$ named.
**Primary**: Vol III Wave 15 Gaiotto verdict; canonical preamble
rows 17--18; Gaiotto 2009 *JHEP* 1208; Chacaltana--Distler 2010.
**Inscription**: Vol II
\texttt{chapters/connections/class\_s\_gaiotto.tex};
\texttt{chapters/connections/agt\_correspondence.tex}.
**Cross-ref**: Vol III AP-CY (R6 primary form); Vol I AP950
(partner); canonical preamble rows 17--19; AP-V2-61 / V2-AP164
(two-route central-charge).

### AP-V2-72 / V2-AP175: MC5 sewing theorem on Stage 2 of two-stage factorisation (2026-04-22)
**Severity**: High (sewing-theorem cross-volume anchor)
**Trigger**: Vol II MC5 sewing theorem is inscribed without
identifying its home on Stage 2 of the Vol III two-stage factorisation
$\Phi_d = \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C} \circ
\Phi^{\mathrm{FA}}_d$; sewing theorem is the stage-2 pushforward
$\int_{\Sigma_{d-1}}$ identification, not a free-standing construction.
**Ghost**: MC5 sewing theorem recovers the $E_1$-chiral shadow on
the reference curve $C$ from the $E_d$-holomorphic factorisation
algebra on the CY-$d$ ambient via factorisation-homology pushforward.
6D hCS at Level 3 realises the $E_3$-algebra underlying the 3D HT QFT,
which in turn lives on the OPEN colour of $\mathsf{SC}^{\mathrm{ch,top}}$
as stage-2 specialisation (AP-V2-25 / V2-AP128).
**Precise error**: MC5 inscribed as free-standing; two-stage anchor
missing; stage-2 specialisation scope lost.
**Correction**: every Vol II MC5 inscription names Stage 2
$\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}$, the factorisation-
homology pushforward $\int_{\Sigma_{d-1}}$, and the family parameter
$(\Sigma_{d-1}, C)$; 6D hCS realisation of $E_3$ named at stage-2
boundary.
**Primary**: Ayala--Francis 2015 *J Topol* 8; Costello 2013
arXiv:1303.2632; Vol III two-stage factorisation wave 2026-04-22.
**Inscription**: Vol II \texttt{chapters/theory/mc5\_sewing.tex}
Stage-2 bridge clause.
**Cross-ref**: Vol III AP-CY (R7 primary form); Vol I AP951
(partner); AP-V2-25 / V2-AP128 (closed/open colour = stage);
AP-V2-26 / V2-AP129 (single-stage framing); AP-V2-59 / V2-AP162
(6d hCS Stage-1).

### AP-V2-73 / V2-AP176: PTVV table to $d = 5$ Poisson-$E_5$ corner closure in Vol II (2026-04-22)
**Severity**: Medium (PTVV $d$-table completeness)
**Trigger**: Vol II PTVV chapter closes the $(d, \mathrm{shift},
E_n^{\mathrm{cl}})$ table at $d = 4$; correct closure is at $d = 5$
with $(5, -4, \mathrm{Poisson}\text{-}E_5)$.
**Ghost**: PTVV extends to all $d$ with shift $-2d + 2$ on the CY-$d$
moduli; Dunn--Lurie $E_3 \otimes E_2 = E_5$ gives Poisson-$E_5$
corner closure at $d = 5$. HT top-form dimension at $d = 5$ matches
$d = 5$ Fake Monster CY home.
**Precise error**: incomplete dimensional enumeration; AP-V2-35 /
V2-AP138 partner at $d = 5$.
**Correction**: Vol II PTVV table includes $d = 5$ row with Poisson-
$E_5$ corner closure; full sequence $(d, \mathrm{shift},
E_n^{\mathrm{cl}}) = (1, 0, E_0), (2, -2, E_1), (3, -2,
\mathrm{Poisson}\text{-}E_3), (4, -2, \mathrm{Poisson}\text{-}E_4),
(5, -4, \mathrm{Poisson}\text{-}E_5)$.
**Primary**: Pantev--Toen--Vaquie--Vezzosi 2013 *Publ IHES* 117;
Vol III Wave 18 PTVV extension.
**Inscription**: Vol II
\texttt{chapters/theory/ptvv\_shifted\_symplectic.tex} $d = 5$ row.
**Cross-ref**: Vol III AP-CY (R8 primary form); Vol I AP952
(partner); AP-V2-35 / V2-AP138; AP-V2-69 / V2-AP172 (Fake Monster).

### AP-V2-74 / V2-AP177: Pentagon single-colour coherence $H^3$ class-dependent at 24 $I_1$ nodes (2A $\in \Z/2$, 2B $\in \Z/4$) (2026-04-22)
**Severity**: Medium (Pentagon $H^3$ Mathieu-class refinement)
**Trigger**: Vol II Pentagon single-colour coherence at 24 $I_1$
nodes on Gaiotto curve $\Sigma_{0, 24}$ states $H^3$ cocycle
uniformly across Mathieu conjugacy classes.
**Ghost**: $H^3 \in \Z/2$ for class 2A (fixing two points of
$\Omega_{24}$), $H^3 \in \Z/4$ for class 2B (fixing one point and
one 2-cycle), per-class for all 25 classes in $M_{24}$. $H^3$ is
a Mathieu-class invariant, not universal.
**Precise error**: conjugacy-class-dependent torsion presented as
class-independent; Mathieu-moonshine refinement lost.
**Correction**: every Vol II Pentagon $H^3$ coherence statement
names Mathieu conjugacy class and torsion value: 2A $\in \Z/2$,
2B $\in \Z/4$, and per-class entries for all 25 classes.
**Primary**: Vol III Wave 18 Mathieu class-dependent $H^3$;
Eguchi--Ooguri--Tachikawa 2010 (Mathieu moonshine for K3);
Conway--Sloane SPLAG Ch 10.
**Inscription**: Vol II
\texttt{chapters/theory/pentagon\_single\_colour.tex} Mathieu-class
table;
\texttt{chapters/connections/mathieu\_moonshine.tex}.
**Cross-ref**: Vol III AP-CY (R9 primary form); Vol I AP953
(partner); AP-V2-63 / V2-AP166 (Macdonald inapplicability);
AP-V2-71 / V2-AP174 (Gaiotto curve 24 punctures); canonical preamble
row 39 (rank-162 MTC).

### AP-V2-75 / V2-AP178: Cubic $d^{abc}$ vs quartic Casimir on $E_2$-bar bulk-boundary anomaly; CANONICAL-ANOM-LOCUS (form c) (2026-04-22)
**Severity**: Critical (6D hCS / 3D HT QFT bulk anomaly discipline; Vol III AP-CY262 partner)
**Trigger**: Vol II $\SC^{\mathrm{ch,top}}$ / 3D HT QFT bulk-boundary
anomaly at the $E_2$-bar $\partial\mathrm{hCS}$ interface cites
only quadratic Casimir $C_2 = 2h^\vee$ for anomaly-free locus, or
admits $E_6$ in a ``cubic-free'' safe list, or cites Deligne
$\setminus \{E_6, A_2\}$ without $A_2$-refined / $A_2$-unrefined
distinction.
**Ghost**: two obstructions operate — \emph{quartic}
$\mathrm{tr}_{\mathrm{adj}}(T^{(a}T^bT^cT^{d)})$ factorising across
the Deligne exceptional series (Cohen--de Man 1996) with the
residual $(F^2)^2$ absorbed by Green--Schwarz, and \emph{cubic}
$d^{abc} = \mathrm{tr}_{\mathrm{adj}}(T^{(a}T^bT^{c)})$ nonzero on
$E_6$ (Jordan cubic on $\mathrm{Sym}^3(\mathbf{27}) =
\mathfrak j_3^{\mathbb O}$) and on $A_2 = \mathfrak{su}(3)$
(Gell-Mann $d$-tensor), zero on $\{A_1, G_2, D_4, F_4, E_7, E_8\}$.
**Precise error**: quadratic-only reading; cubic-only reading with
$E_6$ in safe list; form (a) Deligne $\setminus \{E_6, A_2\}$
strict; form (b) Deligne $\setminus \{E_6\}$ admitting unrefined
$A_2$.
**Correction — CANONICAL-ANOM-LOCUS (form c)**: $\mathrm{Anom}_1 = 0
\iff \mathfrak g \in (\mathrm{Deligne}^{\mathrm{exc}} \setminus
\{E_6, A_2\text{-unrefined}\}) \cup \{\mathrm{abelian}\} \cup
\{\mathrm{super-str}_{\mathrm{ad}} = 0\} \cup
\{\widehat{\mathfrak g}_{-h^\vee} \otimes K^{-1/2}\text{-refined}\}$.
$E_6$ strictly excluded; $A_2$-refined (Feigin--Frenkel critical
twist $K^{-1/2}$ + Dimofte-slab anomaly-inflow from Vol II Part V)
inside; $\{A_1, G_2, D_4, F_4, E_7, E_8\}$ unconditionally inside.
Every Vol II bulk-boundary anomaly inscription names both
obstructions and carries the $A_2$-refined qualifier plus
$K^{-1/2}$-refinement clause.
**Primary**: Deligne 1996 \emph{CRAS} 322; Cohen--de Man 1996
\emph{CRAS} 322; Cvitanović 2008 \emph{Group Theory} Ch 20;
Baez 2002 \emph{Bull AMS} 39; Frampton--Kephart 1983 \emph{Phys
Rev Lett} 50, 1347; Witten 1984 \emph{Comm Math Phys} 92, 455;
Feigin--Frenkel 1992 \emph{Comm Math Phys} 147; Dimofte 2014
slab anomaly-inflow; Vol III AP-CY262; AP-V2-54 / V2-AP157
(earlier Vol II form).
**Inscription**: Vol II
\texttt{chapters/connections/thqg\_perturbative\_finiteness.tex};
\texttt{chapters/theory/e2\_bar\_bulk\_boundary\_anomaly.tex}.
**Cross-ref**: Vol III AP-CY262 (primary form); Vol I AP979 /
Pattern 445; AP-V2-54 / V2-AP157 (earlier Deligne scope);
AP-CY50-E14 (cross-volume ledger); canonical preamble rows
40--42 ($\Psi$-sibling / $V^{s\natural}$ placement).

### AP-V2-76 / V2-AP179: PTVV shift $n = d - 4$ CY-dimension law at $d = 2$ gives shift $-2$; $E_2$-observables (2026-04-22)
**Severity**: High (Vol II $E_2$-scope PTVV discipline; Vol III AP-CY263 partner)
**Trigger**: Vol II PTVV / shifted-symplectic chapters cite the
$(-1)$-shift for hCS-type theories without naming the CY-dimension
$d$; at $d = 2$ (K3, $T^4$, bielliptic) the correct shift is $-2$
with $E_2$-observables, directly within Vol II's $E_2$-bar scope.
**Ghost**: hCS CY-dimension law $(n, E_k) = (d - 4, k = d - 2)$:
shift and $E_n$-index sum to $d - 2$ (the holomorphic disc
dimension after Dolbeault cohomology). At $d = 2$: $n = -2$,
$E_2$-observables; at $d = 3$: $n = -1$, $E_1$-observables; at
$d = 4$: $n = 0$, $E_0$-observables.
**Precise error**: uniform $(-1)$-shift collapses the $d$-dependent
hierarchy; Vol II's natural $E_2$-bar home at $d = 2$ hidden
behind a $d = 3$ quotation.
**Correction**: every Vol II PTVV / hCS / BV / Costello--Gwilliam
reference on CY$_d$ carries its shift $n = d - 4$ and
$E_{d-2}$-observable index; $d = 2$ entries explicitly note
shift $-2$ and $E_2$-factorisation on two holomorphic
disc-directions.
**Primary**: Pantev--Toen--Vaquie--Vezzosi 2013 \emph{Publ IHES}
117 Thm 2.5; Calaque--Pantev--Toen--Vaquie--Vezzosi 2017
arXiv:1506.03699 Prop 2.6; Lurie \emph{HA} \S 5.3; Costello--
Gwilliam 2021 Vol II \S 10--11; Vol III AP-CY263; AP-V2-35 /
V2-AP138 (Vol II PTVV shift / $d$-dependence).
**Inscription**: Vol II
\texttt{chapters/theory/ptvv\_shifted\_symplectic.tex} $d = 2$
shift-$(-2)$ row; \texttt{chapters/theory/e2\_bar.tex} PTVV
anchor.
**Cross-ref**: Vol III AP-CY263 (primary form); AP-V2-35 /
V2-AP138 (PTVV shift / $d$-dependence); AP-V2-73 / V2-AP176
(PTVV $d = 5$ closure).

### AP-V2-77 / V2-AP180: $\HH^0_{E_3}$ vs $\HH^0_{E_2}$ rank at bulk-boundary balance; compact-vs-open $E_2$ dualizability (2026-04-22)
**Severity**: Critical (Vol II $E_2$-dualizability discipline; Vol III AP-CY265 partner)
**Trigger**: Vol II $\SC^{\mathrm{ch,top}}$ bulk-boundary balance
statements cite $E_3$-Hochschild rank finiteness on open $\C^3$,
or conflate $\HH^0_{E_3}$-rank on $\C^3$ with $\HH^0_{E_2}$-rank
on $\C^2$ at the 3D HT QFT bulk / 2D boundary interface.
**Ghost**: on non-compact $\C^3$, $\HH^0_{E_3}(\Obs_{\mathrm{hCS}}
(\C^3, \mathfrak g)) \simeq \mathbb C[\![\tau_1, \tau_2, \tau_3]\!]$
(infinite rank, not dualizable in $\mathrm{Ch}$); on compact CY$_3$
$X$, $\HH^0_{E_3}(\Obs_{\mathrm{hCS}}(X, \mathfrak g)) \simeq
\bigoplus_q h^{0,q}(X) \cdot \HH^*_{\mathrm{Lie}}(\mathfrak g,
\mathbb C)^{[q]}$ (finite rank via Hodge truncation). The
$E_3$-dualizability gap is decisive for the bulk-boundary balance:
3D HT QFT on compact CY$_3$ lifts to a fully extended framed 3-TFT
(Lurie 2009 Thm 2.4.6); $\C^3$-locus does not.
**Precise error**: compactness dropped silently; $E_3$-dualizability
tied to $E_1$-output scope only; bulk-boundary balance breaks
when rank-infinity meets $E_2$-boundary finite-rank target.
**Correction**: every Vol II bulk-boundary statement at the
3D HT QFT / $E_2$-bar interface names compactness of the bulk
CY$_3$; $\C^3$-locus arguments explicitly flagged as
non-dualizable and non-extendable to compact CY$_3$. The
$\HH^0_{E_3}$ / $\HH^0_{E_2}$ rank match on the bulk-boundary
interface requires compact bulk.
**Primary**: Lurie 2009 \emph{On the classification of TFTs}
Thm 2.4.6; Gwilliam--Williams 2021 arXiv:2009.05037 Prop 5.3.2;
Francis 2013 \emph{Compos Math} 149 Thm 3.4; Ayala--Francis 2015
\emph{J Topology} 8 Thm 1.1; Calaque--Pantev--Toen--Vaquie--
Vezzosi 2017 arXiv:1506.03699; Vol III AP-CY265.
**Inscription**: Vol II
\texttt{chapters/theory/sc\_cht\_bulk\_boundary\_balance.tex};
\texttt{chapters/connections/thqg\_perturbative\_finiteness.tex}.
**Cross-ref**: Vol III AP-CY265 (primary form); AP-V2-34 /
V2-AP137 ($E_n$-on-wrong-object); AP-V2-25 / V2-AP128 (mixed
holomorphic pushforward); AP-V2-59 / V2-AP162 (6D hCS two-stage).

### AP-V2-78 / V2-AP181: $\HH^2_{E_2}$ rigidity non-critical vs critical Kac--Moody level (2026-04-22)
**Severity**: High (Vol II $E_2$-bar level-scope discipline; Vol III AP-CY266 partner)
**Trigger**: Vol II $E_2$-bar rigidity / deformation-theoretic
chapters state $\HH^2_{E_2}(\partial\mathrm{hCS}_5(\mathfrak{sl}_n),
-) = 0$ as uniform across all Kac--Moody levels, including critical
$k = -h^\vee$; directly in Vol II's $E_2$-bar domain.
**Ghost**: at \emph{generic} (non-critical) level, the Francis
chiral-tangent identification $\HH^2_{E_2}(\partial\mathrm{hCS}_5
(\mathfrak{sl}_n), -) \simeq H^2_{\mathrm{ch}}
(\widehat{\mathfrak{sl}}_n, V^{\mathrm{vac}}_k)$ reduces rigidity
to Whitehead 2 on semisimple Lie cohomology via Francis
spectral-sequence degeneration. At critical level $k = -h^\vee$,
$V^{\mathrm{vac}}_{-h^\vee}$ is NOT simple: its centre
$Z(V^{\mathrm{vac}}_{-h^\vee}) \simeq \mathrm{Fun}(\mathrm{Op}_
{\mathfrak{sl}_n^\vee})$ (Feigin--Frenkel; Frenkel--Ben-Zvi 2004
Thm 18.4.2) is polynomial in infinitely many generators; rigidity
requires separate analysis via opers-moduli smoothness.
**Precise error**: uniform-rigidity claim elides critical-level
subtlety; Whitehead 2 does not extend; Feigin--Frenkel
oper-smoothness argument dropped.
**Correction**: every Vol II $\HH^2_{E_2}$ rigidity statement for
affine Kac--Moody vertex algebras specifies level scope: generic
(non-critical) via Whitehead 2 + formality of $V^{\mathrm{vac}}_k$
+ Getzler--Kapranov spectral sequence; critical via Feigin--
Frenkel centre-as-functions on $\mathrm{Op}_{\mathfrak{sl}_n^\vee}$
+ oper moduli smoothness. The two arguments are distinct.
**Primary**: Francis 2013 \emph{Compos Math} 149 Thm 1.1, 2.29;
Frenkel--Ben-Zvi 2004 \emph{Vertex Algebras and Algebraic Curves}
Thm 3.4.3, 18.4.2; Feigin--Frenkel 1992 \emph{IJMPA} 7 S1A;
Whitehead 1937; Vol III AP-CY266.
**Inscription**: Vol II
\texttt{chapters/theory/e2\_bar\_rigidity.tex} level-scope clause;
\texttt{chapters/theory/kac\_moody\_deformations.tex}.
**Cross-ref**: Vol III AP-CY266 (primary form); AP-V2-34 /
V2-AP137 ($E_n$-on-wrong-object); canonical preamble central-charge
discipline; AP-V2-58 / V2-AP161 (Yangian-type discipline,
critical-level sibling).

### AP-V2-79 / V2-AP182: $\Lambda^{3,3}$-envelope rank-6 GBKM vs Fake-Monster $\mathrm{II}_{25,1}$ confusion (2026-04-22)
**Severity**: High (Vol II rank-6 envelope / Humbert-divisor discipline; Vol III AP-CY267 partner)
**Trigger**: Vol II Swiss-cheese / rank-6 envelope references place
$\mathfrak g_{\Delta_5}$ inside the Fake-Monster Lie algebra
$\mathfrak g_{\Phi_{12}}$ via a codimension-one timelike restriction
$\Lambda^{3,2} \hookrightarrow \Lambda^{3,3}$.
**Ghost**: $\Lambda^{3,3} = U \oplus U \oplus U$ is the unique
even unimodular lattice of rank 6 signature $(3,3)$ (Milnor;
Conway--Sloane Ch 15); the envelope GBKM on it has Cartan rank 6
signature $(3,3)$, hosting $\mathrm{Sp}_4(\mathbb Z)/\{\pm I\}
\simeq \mathrm{SO}_+(\Lambda^{3,2})$ at the Humbert-divisor level
(Gritsenko--Nikulin 1998 \emph{Duke} 94). The Fake-Monster lives
on $\mathrm{II}_{25,1}$ (rank 26, signature $(25,1)$, at $d = 5$
per canonical preamble / AP-V2-69), an entirely distinct lattice.
The Humbert-divisor fibre $\mathfrak g_{\Delta_5} \subset
\mathfrak g_{\Lambda^{3,3}}$ is a genuine codimension-one BKM
restriction, NOT a specialisation of the Fake Monster.
**Precise error**: rank-26 Fake-Monster ambient conflated with
rank-6 $\Lambda^{3,3}$-envelope ambient; two entirely distinct
BKM structures on two entirely distinct lattices.
**Correction**: every Vol II rank-6 envelope / $\mathfrak g_{\Delta_5}$
ambient statement names host lattice: $\Lambda^{3,3}$ rank 6 is
the Pfaffian envelope (heterotic BPS lift, Harvey--Moore 1996);
$\mathrm{II}_{25,1}$ rank 26 is the Fake-Monster non-compact host;
these do not inclusion-compare.
**Primary**: Borcherds 1995 \emph{Invent Math} 120 \S 13
(restriction theorem); Gritsenko--Nikulin 1998 \emph{Duke} 94
\S 2 ($\mathrm{Sp}_4 \simeq \mathrm{SO}(3,2)$); Harvey--Moore
1996 hep-th/9510182 eqs (4.15)--(4.20); Nikulin 1979 \emph{Izv AN
SSSR} 43; Vol III AP-CY267.
**Inscription**: Vol II
\texttt{chapters/connections/rank\_six\_envelope.tex};
\texttt{chapters/theory/humbert\_divisor\_restriction.tex}.
**Cross-ref**: Vol III AP-CY267 (primary form); AP-V2-69 /
V2-AP172 (Fake-Monster $d = 5$); canonical preamble rows 20--21
(Monster / Fake-Monster Cartan-rank discipline);
AP-V2-66 / V2-AP169 ($F_3$ real-root subalgebra of
$\mathfrak g_{\Delta_5}$).

### AP-V2-80 / V2-AP183: Chiral Booth--Lazarev as three-obstruction programme vs associative packaging on Vol II curved coderived (2026-04-22)
**Severity**: High (Vol II curved coderived framework discipline; Vol III AP-CY268 partner)
**Trigger**: Vol II curved-coderived / chiral-bar--cobar chapters
cite chiral Booth--Lazarev $\mathrm{FactCoAlg}^{\mathrm{cpt}}_
{\mathrm{crv}}(\mathrm{Ran}(C)) \simeq \mathrm{ChirAlg}(C)$ as
automatic lift of the associative curved bar--cobar Quillen
equivalence via ``factorisation-tensor packaging''.
**Ghost**: three distinct structural obstructions are unsettled.
(i) \emph{Ran-space Smith recognition}: $\mathrm{Ran}(C)$ is a
colimit of schemes $C^n / S_n$; generating cofibrations must be
stratification-compatible with $j_n^*$ restrictions pulling back
into lower-$m$ generators; conilpotence stratified by the
Ran-space filtration, not a single unstratified filtration.
(ii) \emph{Genus-tower curvature as operadic section}: curvature
$m_0^{(g,n)}$ is a family of sections of $\lambda_g \otimes
H^2(\mathcal A)$ on $\overline{\mathcal M}_{g,n}$ satisfying
boundary-clutching compatibility (separating + non-separating);
a modular form on $\overline{\mathcal M}$ with sewing as operadic
structure, not a scalar in $\mathcal A^2$.
(iii) \emph{Analytic IndHilb sewing}: the algebraic coderived
category and the Moriwaki IndHilb analytic sewing envelope are
distinct; comparison requires nuclearity + trace-class check
on topological vector spaces that Booth--Lazarev's purely
algebraic setup does not address.
**Precise error**: conjectural three-obstruction lift stated as
``the Quillen equivalence on Ran-space coalgebras''; overstates
theorem-status relative to Vol II's curved-coderived framework.
**Correction**: every Vol II chiral Booth--Lazarev citation
carries the three-obstruction qualifier;
$\kappa_{\mathrm{ch}}$-curving $m_0^{(g)} = \kappa_{\mathrm{ch}}
\cdot \lambda_g$ is the \emph{boundary-descent datum} on
$\overline{\mathcal M}_{g,n}$, not a derived-category-level
packaging. The three obstructions (Ran-space Smith,
genus-tower curvature operadic structure, analytic sewing)
named jointly.
**Primary**: Booth--Lazarev 2023 arXiv:2304.08409 (associative
coderived Quillen); Beilinson--Drinfeld 2004 \emph{Chiral
Algebras} AMS Colloq Publ 51 Ch 3 (Ran space, factorisation);
Francis 2013 \emph{Compos Math} 149 ($E_n$-tangent, operadic
model structure); Kontsevich--Soibelman 2010
\emph{Deformation Theory} Ch 3 ($L_\infty$ on smooth schemes);
Vol III AP-CY268.
**Inscription**: Vol II
\texttt{chapters/theory/curved\_coderived.tex} three-obstruction
qualifier; \texttt{chapters/theory/chiral\_bar\_cobar.tex}.
**Cross-ref**: Vol III AP-CY268 (primary form); AP-V2-72 /
V2-AP175 (MC5 sewing on Stage 2); AP-V2-30 / V2-AP133 (chain-level
vs $(\infty,1)$-functor lift); AP-V2-41 / V2-AP144 (frontier
$(\infty,1)$-equivalence).


### AP-V2-81 / V2-AP184: $\widehat{\mathfrak{sl}}_2$ vs $\widehat{\mathfrak{gl}}(1|1)$ scope-discipline on compact-Koszul descents (2026-04-23, sibling Vol III AP-CY305).
$\mathrm{CoHA}(Y^{\mathrm{conifold}}) \cong Y^+(\widehat{\mathfrak{gl}}(1|1))$ is the primary Vol III identification (Li-Yamazaki arXiv:2003.08909 §8.3.6.3); $Y^+(\widehat{\mathfrak{sl}}_2)$ is the ungraded shadow via supertrace projection of the central Cartan. On Vol II curved-coderived compact-Koszul descents the two appear at distinct scopes: the supergraded form appears on $E_3$-factorisation over $\mathbf Y$ (primary), the ungraded form via the semisimple quotient. **Counter**: Vol II compact-Koszul applications to conifold-type CY$_3$ must specify the super graded scope.

### AP-V2-82 / V2-AP185: Super 5D hCS $\to$ $E_1$-chiral all-orders is OPEN on super gauge algebras (2026-04-23, sibling Vol III AP-CY298).
CGY arXiv:1810.01970 5D hCS $\to$ Yangian VOA all-orders requires simply-laced Lie $\mathfrak g$ with non-degenerate Killing. For $\mathfrak g = \mathfrak{gl}(m|n)$ the Killing is degenerate and super-KT is $E_2$-only; all-orders extension OPEN. **Counter**: Vol II $E_1$-chiral claims on super gauge algebras must carry the OPEN flag beyond 1-loop.

### AP-V2-86 / V2-AP186: Super-KT formality is $E_2$-only; $E_3$ open (2026-04-23, sibling Vol III AP-CY299).
Kontsevich-Tamarkin $E_n$-formality extends to super dg-Lie in the Ginzburg-Schedler arXiv:0807.0174 sense only to $E_2$. Vol II $E_3$-chiral statements that invoke super-formality must carry OPEN flag at $E_3$. **Counter**: super-formality on Vol II $E_3$-factorisation constructions must specify $E_2$.

### AP-V2-87 / V2-AP187: Two-chart Čech atlas is not a Weiss cover — QC descent vs factorisation locality (2026-04-23, sibling Vol III AP-CY308).
Wrong: treating the two-chart Čech atlas $\{U_+, U_-\}$ of a resolved CY$_3$ (e.g. conifold, where each $U_\pm \cong \mathbb C^3$) as a "factorisation-algebra gluing" datum on Vol II $E_n$-factorisation constructions. Correct: Weiss covers (Costello-Gwilliam arXiv:2210.13036 Def 6.1.6) require every finite configuration of points to embed simultaneously into a single open; a two-chart atlas fails this whenever a configuration has points on both charts. Hence the two-chart atlas recovers QC-descent for $\mathrm{CoHA}$-as-sheaf (Kontsevich-Soibelman arXiv:1006.2706 §6.3) but NOT factorisation-algebra locality (Beilinson-Drinfeld, Francis-Gaitsgory arXiv:1103.5803 Thm 3.6.2, Francis arXiv:1303.0305, Lurie HA §5.5.4); the latter requires the Weiss refinement $\mathfrak D^{\sqcup}$ of disjoint unions of contractible discs. **Counter**: Vol II $E_n$-factorisation on toric CY$_3$ specifies QC-descent (minimal refinement sufficient) vs factorisation-locality (Weiss refinement required); these are not interchangeable.

### AP-V2-88 / V2-AP188: $\kappa_{\mathrm{ch,BV}}$ is a genuinely new subscripted invariant distinct from $\kappa_{\mathrm{ch}}$ on non-compact CY$_3$ (High, sibling Vol III AP-CY324).
Wrong: identifying three apparent values $\{+1, 0, -1\}$ for "$\kappa_{\mathrm{ch}}$(conifold)" as alternative computations of one invariant on Vol II BV / $E_n$-factorisation content. Correct: they are three distinct subscripted invariants under one chain-level convention (Kontsevich-Soibelman DT / Costello-Li holomorphic BV): $\kappa_{\mathrm{ch}}$ (DT / motivic BPS count), $\kappa_{\mathrm{cat}} = \chi(\mathcal O_X)$ (holomorphic Euler characteristic), $\kappa_{\mathrm{ch,BV}}$ (Costello-Li BCOV 1-loop curving via Polyakov 1981 Ch.~9 ghost supertrace, value $\mathrm{str}_{\mathfrak{gl}(1|1)}(\mathrm{ghost}) = -1$ at BRST $c = -2$). Polyakov ghost-mode balance $\kappa_{\mathrm{ch}} + \kappa_{\mathrm{ch,BV}} = \kappa_{\mathrm{cat}}$ holds on $\mathsf G$-class free-field CY$_3$ with $\chi_{\mathrm{top}} = 2$, NOT universal. **Counter**: Vol II BV / BCOV content on non-compact CY$_3$ must use the five-$\kappa$ ladder $\{\kappa_{\mathrm{ch}}, \kappa_{\mathrm{cat}}, \kappa_{\mathrm{BKM}}, \kappa_{\mathrm{fiber}}, \kappa_{\mathrm{ch,BV}}\}$; $\kappa_{\mathrm{ch,BV}}$ is a load-bearing new subscript for BV-anomaly-sensitive content.

### AP-V2-89 / V2-AP189: 5D hCS BV anomaly on super gauge: 1-loop VERIFIED, 2-loop+ OPEN (Medium, sibling Vol III AP-CY323).
Wrong: Vol II BV / factorisation content citing "Costello-Yagi all-orders theorem applies to $\mathfrak{gl}(1|1)$ / super gauge" without specifying obstruction status. Correct: the BV anomaly coefficient is $d^{abc} = \mathrm{str}(t^a \{t^b, t^c\})$; on $\mathfrak{gl}(1|1)$ direct basis computation gives $d^{abc} = 0$ identically (Cartan $\mathrm{str}(H^2) = \mathrm{str}(N^2) = 0$; odd pair $\psi^\pm$ carries opposite supertrace signs), so the 1-loop wheel is verified to vanish. Higher-loop local cohomology $H^1_{\mathrm{loc}}(\mathfrak{gl}(1|1), \mathcal O_{\mathrm{loc}})^{\geq 2}$ is not auto-killed (super-KT is $E_2$-only); 2-loop and higher OPEN. **Counter**: Vol II BV-anomaly stratifies on super gauge: 1-loop verified via basis, higher-loop OPEN per $E_2$-only super-KT.

### AP-V2-90 / V2-AP190: Bialgebra vs Hopf on $E_n$-factorisation positive halves (Medium, sibling Vol III AP-CY327).
Wrong: labelling $Y^+(\mathfrak g)$ (positive half of a Yangian / quantum affine algebra) as "Hopf" on Vol II $E_n$-factorisation content. Correct: $Y^+$ is an associative bialgebra; the antipode $S$ lives only on the Drinfeld double $D(Y^+) = Y$. Super positive halves (from $\widehat{\mathfrak{gl}}(1|1)$-type data) are $\mathbb Z_2$-graded bialgebras in super-vector-spaces. Strict-Hopf on rational / trigonometric / toroidal equivariance strata; elliptic case is quasi-Hopf with Felder-Jimbo-Konno dynamical associator. **Counter**: Vol II $E_n$-factorisation on positive-half quantum affinisations uses "bialgebra" / "$\mathbb Z_2$-graded bialgebra" / "quasi-Hopf" per stratum; never blanket "Hopf".

### AP-V2-91 / V2-AP191: BCOV cocycle target is $H^{0,1}(X, \mathcal O_X)$, not $H^{0,1}(X, \mathrm{Sym}^{\leq 2} T_X^*)$ (Medium, sibling Vol III AP-CY341).
Wrong: writing BCOV 1-loop anomaly class on Vol II BV / $E_n$-factorisation as $\alpha_{\mathrm{BCOV}}(X) = (\chi(X) / 24)[\Omega_X]^{0,1} \in H^{0,1}(X, \mathrm{Sym}^{\leq 2} T_X^*)$. Correct: BCOV anomaly class lives in $H^{0,1}(X, \mathcal O_X)$ — the $[\Omega_X]^{0,1}$ Dolbeault class is an $\mathcal O_X$-coefficient class (the $(0, 1)$-Dolbeault lift of the Atiyah class of the CY trivialisation $\omega_X \cong \mathcal O_X$), not a polyvector-field class. Primary: Costello-Li arXiv:1606.00365 Prop 5.2; BCOV hep-th/9309140; Polyakov 1981 Ch.~9. **Counter**: Vol II BCOV-target claims place the $(\chi / 24)$ class in $H^{0,1}(X, \mathcal O_X)$; higher-symbol targets are category errors.

### AP-V2-92 / V2-AP192: BCOV factor-split admits four vanishing mechanisms on compact / non-compact CY$_3$ (Low, sibling Vol III AP-CY342).
Wrong: claiming $\alpha_{\mathrm{BCOV}}(X) \neq 0$ on all non-quintic compact CY$_3$ in Vol II BV / (O2)-obstruction content. Correct: the factor-product $(\chi(X) / 24) \cdot [\Omega_X]^{0,1}$ admits four distinct vanishing mechanisms — (i) conifold / non-compact retractable: $[\Omega_{\mathbf Y}]^{0,1} = 0$ via retraction $\mathbf Y \simeq \mathbb P^1$; (ii) strict quintic / $h^{0,1} = 0$: Dolbeault factor vanishes; (iii) $K3 \times E$ Künneth: $\chi = 24 \cdot 0 = 0$ topological prefactor; (iv) generic CY$_3$ with abelian factor: product vanishes. Only compact CY$_3$ without abelian factor and with $h^{0,1} > 0$ carries genuinely non-zero $\alpha_{\mathrm{BCOV}}$. **Counter**: Vol II (O2)-obstruction statements identify which vanishing mechanism applies; blanket "compact $\Rightarrow$ non-zero $\alpha_{\mathrm{BCOV}}$" is wrong.

### AP-V2-93 / V2-AP193: 4 compact-CY$_3$ obstructions reduce to 2.5, not 4 independent (Medium, sibling Vol III AP-CY309).
Wrong: listing (O1) toric fan-completeness, (O2) BCOV $\alpha_{\mathrm{BCOV}}$, (O3) $\mathrm{Aut}^0$ rigidity, (O4) finite-quiver equivariance as 4 independent obstructions to compact-CY$_3$ chart-wise factorisation-algebra assembly on Vol II $E_n$ content. Correct: direct logical reductions give (O1) $\Leftrightarrow$ (O3) (two facets of Bogomolov-Matsumura rigidity) and (O1) $\Rightarrow$ (O4), with only (O2) logically independent. Effective count is 2.5: (O1 $\equiv$ O3) + (O2) with (O4) downstream. BCOV curving (O2) is independent because it sources from 1-loop BV counter-term, not automorphism geometry. **Counter**: Vol II compact-CY$_3$ obstruction statements derive the 2.5-count explicitly; "four independent obstructions" obscures the Bogomolov-Matsumura origin.

### AP-V2-94 / V2-AP194: $B_3$-braid on local $\mathbb P^2$ via $\pi_1(\mathrm{Conf}_3(\mathbb C))$, not Bondal-Orlov 2002 alone (Medium, sibling Vol III AP-CY340).
Wrong: citing Bondal-Orlov arXiv:math/0206295 as the primary source for the $B_3$-braid mutation action on local $\mathbb P^2$ nine-arrow quiver in Vol II factorisation / derived-category content. Correct: Bondal-Orlov 2002 is a derived-equivalence paper and does NOT establish braid-group action on mutations. Correct primary sources: Bondal-Polishchuk 1993 (helix autoequivalences); Kuznetsov arXiv:math/0610957 ($B_n$-action on exceptional collections); Bridgeland arXiv:0909.4299 ($B_3 = \pi_1(\mathrm{Conf}_3(\mathbb C))$ on 3-vertex derived categories); Seidel-Thomas arXiv:math/0001043 (spherical-twist braid relations from Ext-adjacency). **Counter**: $B_3$-mutation-action citations in Vol II factorisation content separate helix autoequivalence (Bondal-Polishchuk), exceptional-collection mutation (Kuznetsov), configuration-space identification (Bridgeland), and spherical-twist braid (Seidel-Thomas); Bondal-Orlov 2002 alone is insufficient.

### AP-V2-95 / V2-AP195: Finite Rees hCS--Hall is not compact critical CoHA (High, sibling Vol III AP-CY345).
Finite DWR/Ran simplex integration constructs a Maurer--Cartan element valued in an oriented Rees Hall complex. Vol II hCS / BV / factorisation prose must not identify this finite Rees object with ordinary compact critical CoHA. Completed Rees requires compatible inverse limits and Mittag--Leffler control; realized critical CoHA requires monoidal vanishing-cycle realization. **Counter**: name the layer before using a Hall comparison.

### AP-V2-96 / V2-AP196: Multi-chart gluing is finite-Rees constructed, completed realization remains open (Medium, sibling Vol III AP-CY346).
The latest finite construction gives explicit DWR/Ran gluing by face-compatible cyclic contractions over $\Omega^\bullet(\Delta^p)$ and Stokes assembly. The open Vol II obligation is not finite gluing; it is completed/realized compact CoHA and compact-support functoriality. **Counter**: do not repeat "no gluing outside $\mathbb C^3$" without the finite/completed distinction.

### AP-V2-97 / V2-AP197: $\CoHA(\mathbb C^3)=Y^+(\widehat{\mathfrak{gl}}_1)$, not $\mathcal W_{1+\infty}$ (High, sibling Vol III AP-CY347).
Vol II quantum-group and hCS comparisons must place $\mathcal W_{1+\infty}$ after Hall--Drinfeld doubling, derived/chiral centre, and vacuum/Fock or evaluation representation. The cohomological Hall algebra itself is the positive half. **Counter**: positive half, double, centre, representation/evaluation, and vertex algebra are separate layers.

### AP-V2-98 / V2-AP198: Toric vertices are not the full finite cyclic atlas (Medium, sibling Vol III AP-CY348).
Affine toric charts provide local vertices only. Finite hCS--Hall comparison also needs simplex-compatible cyclic contractions, mutation coherence, determinant-line orientation transport, Thom--Sebastiani compatibility, and Hall product compatibility. **Counter**: Vol II toric factorisation statements name the whole DWR/Ran atlas.

### AP-V2-99 / V2-AP199: Non-formality enters the cyclic potential (Medium, sibling Vol III AP-CY349).
Local $\mathbb P^2$ non-formality is not a terminal obstruction by itself. Higher $m_k$ operations are encoded in the cyclic potential $W_\sigma$ of the finite model. The remaining tasks are contraction, mutation, completion, and realization. **Counter**: do not treat non-formality as absence of a finite Rees critical model.

### AP-V2-100 / V2-AP200: Reduced-DT source labels must be OPi vs OP precise (Medium, sibling Vol III AP-CY350).
Use Oberdieck--Pixton only when Pixton's component is actually invoked; use Oberdieck--Pandharipande or the relevant Oberdieck theorem for reduced DT anchors. **Counter**: Vol II reduced-DT cross-references name the precise paper/source.

### AP-V2-101 / V2-AP201: Quasi-NCCR character does not construct compact CoHA (High, sibling Vol III AP-CY351).
A quasi-NCCR character identity is a character statement. It does not create compact-support descent, vanishing-cycle realization, orientation transport, or the Hall comparison map. **Counter**: character identities are never promoted to object construction.

### AP-V2-102 / V2-AP202: Layer-aware $Y^+$ and Hall notation (High, sibling Vol III AP-CY352).
Vol II uses separate notation or explicit qualifiers for finite Rees Hall complex, completed Rees Hall algebra, realized positive critical CoHA, and Hall--Drinfeld double. **Counter**: every map crossing a layer is named.

### AP-V2-103 / V2-AP203: Four-case hCS--Hall status (Medium, sibling Vol III AP-CY353).
$\mathbb C^3$ positive half is constructed; finite multi-chart Rees is constructed under cyclic-atlas hypotheses; completed Rees is conditional on pro-compatibility and Mittag--Leffler; realized compact critical CoHA is conditional on monoidal vanishing-cycle realization. **Counter**: Vol II summaries use the four-case stratification.

### AP-V2-104 / V2-AP204: CHL and Gritsenko--Clery ladders are different (Medium, sibling Vol III AP-CY354).
Both use $\kappa_{\mathrm{BKM}}=c(0)/2$, but their indexing families and constant terms differ. **Counter**: name CHL or Gritsenko--Clery before quoting constants.

### AP-V2-105 / V2-AP205: $\Delta_5$, $\Phi_{10}$, OP $\chi_{10}$, and $\Phi_{12}$ separation (High, sibling Vol III AP-CY355).
$\Delta_5$ is the primitive K3 paramodular BKM denominator; the unnormalised Igusa square is \(\Phi_{10}^{\mathrm{un}}=\Delta_5^2\); the OP/DT convention uses \(D_5=64^{-1}\Delta_5\) and \(\chi_{10}^{\mathrm{OP}}=D_5^2\); $\Phi_{12}$ is the Fake-Monster denominator. **Counter**: Vol II traces and partition functions name the denominator and the product convention before quoting a scalar.

### AP-V2-106 / V2-AP206: Determinant Hodge line notation $\lambda_1^{\det}$ (Low, sibling Vol III AP-CY356).
Use $\lambda_1^{\det}$ in the BL/DWR bridge; do not write $\lambda_g^1$ for the determinant Hodge line. **Counter**: genus-indexed Hodge-bundle notation is reserved for actual genus-indexed Hodge bundles.

### AP-V2-107 / V2-AP207: OP/DT scalar normalization uses \(D_5=64^{-1}\Delta_5\) (High, sibling Vol III AP-CY357).
The \(K3\times E\) OP/DT scalar is not the bare inverse
\(-\Delta_5^{-2}\). The primitive BKM denominator is \(\Delta_5\), but
the OP scalar convention uses the monic Borcherds product
\[
D_5=64^{-1}\Delta_5,\qquad
\chi_{10}^{\mathrm{OP}}=D_5^2=4096^{-1}\Delta_5^2,
\]
and therefore
\[
Z_{\mathrm{OP/DT}}=-(\chi_{10}^{\mathrm{OP}})^{-1}
=-D_5^{-2}=-4096\,\Delta_5^{-2}.
\]
The unnormalised Igusa square \(\Phi_{10}^{\mathrm{un}}=\Delta_5^2\)
is a different convention. **Counter**: Vol II traces, line-operator
characters, and endpoint comparisons name \(D_5\) or
\(\chi_{10}^{\mathrm{OP}}\) before using an OP/DT scalar.

### AP-V2-108 / V2-AP208: \(K3\times E\) spectrum is \(\{0,3,5,24\}\); \(2\) is fibre-only (High, sibling Vol III K3xE cache).
The \(K3\times E\) cache carries four construction values, not one
bare \(\kappa\): \(\kappa_{\mathrm{cat}}(K3\times E)=0\) on the total
space by Kunneth multiplicativity; the chiral Kunneth/CY dimension is
\(3\); the primitive K3 paramodular BKM denominator \(\Delta_5\) has
\(\kappa_{\mathrm{BKM}}=5\); and the Mukai/Heisenberg--Fock trace is
\(24\). The number \(2=\chi(\mathcal O_{K3})\) is only the K3-fibre
witness. **Counter**: Vol II class-\(\mathcal S\), HT partition
function, and BKM comparison prose uses \(\{0,3,5,24\}\) for the
\(K3\times E\) spectrum and never promotes \(2\) to a total-space
invariant.


---

# Imported Vol I Antipattern Registry Surface, 2026-04-30

Source: `/Users/raeez/chiral-bar-cobar/notes/antipatterns_catalogue.md`.
Exact Vol II mirror: `notes/vol1_imported_antipatterns_catalogue_2026_04_30.md`.
Source SHA-256: `b954958f06dfe486fe052bb98ac29ad0ac374692ae02895bcefde412981b02ba`; source line count: `6086`.

Import discipline: this block is active Vol II attack-heal control-surface material because these antipatterns were stored in Vol I while governing Vol II, topological-string, and cross-volume work. It is not theorem content by itself. Future Vol II work must consult it, translate entries into V2-native codes when used, and admit a claim only after first-principles proof, direct computation, or primary-source citation.

<!-- BEGIN VOL1_IMPORTED_ANTIPATTERN_REGISTRY_2026_04_30 sha256=b954958f06dfe486fe052bb98ac29ad0ac374692ae02895bcefde412981b02ba -->

# Anti-Pattern Catalogue (Vol I)

Extracted from the legacy CLAUDE.md (2026-04-18). The live CLAUDE.md

> **T5 supersession note (2026-04-24).** Entries below written before
> the Vol II T5 repair may still record the principal
> $\mathcal W_N$ or generic $\mathcal W_\infty[\mu]$ topologisation
> ladder as conditional on an unproved higher-spin antighost axiom.
> That status is superseded: Vol II proves principal-Casimir T5 in
> `thm:casimir-antighost-commutativity` and gives closed raw
> normal-ordered homotopies in
> `prop:closed-normal-ordered-antighost-homotopies`. The T5 axiom
> remains only for arbitrary non-Casimir stress towers.

## Canonical values at the latest-wave verdict (2026-04-21)

When any entry in this catalogue or the first-principles cache appears
to assert one value while a later entry asserts another, the **latest
wave's verdict takes precedence**. The registry below pins the
canonical value for every quantity that has flipped at any point across
Waves 1-26. Every AP/cache entry that names an older value does so only
in retraction context (explicitly flagged as "wrong claim" / "old draft"
/ "retracted").

| Quantity | Canonical value (latest-wave verdict) | Retracted earlier values | Lock |
|---|---|---|---|
| $(c_{4d}, c_{2d})$ of $\mathcal T[A_1, \Sigma_{0, 24}]$ | $(107/6, -214)$ (Wave 15 Gaiotto) | $(26, -312)$ (Wave 14 error); $(23/4, -69)$ (Wave 25 Gaiotto agent error) | WOV-2 |
| Trinion/tube count at $n=24$ | $(n_v, n_h) = (63, 88)$ | $(21, 27)$ (Wave 25 confusion of trinion/tube count with multiplicities) | Chacaltana--Distler 2010 §5.14 |
| Universal $c_{4d}$ formula, $A_1$ genus 0 | $c_{4d} = (5n - 13)/6$ | $(12(g-1) + 7n)/6$ (Wave 14 error — fails $\mathrm{SU}(2)$ $N_f = 4$) | Shapere--Tachikawa 2008 |
| Monster BKM Cartan rank | $2$ (on $\mathrm{II}_{1,1}$) | $26$ (Wave 16 confusion with Fake-Monster) | Borcherds 1992 *Invent Math* 109 |
| Fake-Monster Cartan rank | $26$ (on $\mathrm{II}_{25,1}$) | --- | Borcherds 1992 |
| K3-BKM Cartan rank | $3$ (on $\Lambda^{2,1}_{II}$) | --- | Gritsenko--Nikulin 1998 §3 |
| $c_3$ (Borcherds coefficient) | $-8$ in Bruinier reduced-class convention | $176256$ (= $p_{24}(5)$, Wave 16 error; also the Gritsenko--Nikulin Cartan-matrix value with conversion factor $-22032$) | WOV-5 |
| $\zeta(3,3,3,3)$ numerical value | $0.000295999\ldots$ (depth-4, weight 12) | $0.0028565$ (Wave 17 draft, 10× error) | Brown 2012 motivic basis |
| $c(1, 2, \pm 2)$ heterotic lift coefficient | $-2$ (Wave 22 Witten) | $+1$ (Wave 21 error) | $\eta^{18} \theta_1^2$ direct expansion + DMZ 2012 + holomorphic-anomaly (three-path) |
| Humbert $H_4$ arithmetic description | $(2, 2)$-isogeny quotient of $E_1 \times E_2$ with $\mathrm{End} \supset \mathbb{Z}[2i]$; monodromy order 2 | $\mathbb{Q}(\sqrt 2)$-RM locus (Wave 15 imprecise; that is $H_8$, not $H_4$) | van der Geer 1988 Ch IX |
| Theorem B scope (Koszul locus) | $\overline{\mathcal A_2} \setminus \bigcup_{n \,\mathrm{admissible}} H_n$ (all admissible Heegner divisors) | $\overline{\mathcal A_2} \setminus (H_1 \cup H_4)$ (Wave 15 narrow) | Wave 18 Beilinson tightening |
| K3-BKM Weyl denominator | $\Delta_5 = \mathrm{Grit}(\eta^9 \vartheta_1) \in S_5(K(1))$ (additive) | $\Phi_{12}$ attempted (Fake-Monster, $\mathrm{II}_{25,1}$, multiplicative — different object) | Gritsenko 1999 Thm 6.1 |
| Umbral $A_{N-1}$ labelling rule | $(N-1) \mid 24$; at $N = 6$ reanchor to Niemeier $6 D_4$ | $N \mid 24$ (Wave 18 error) | Wave 19 Gaiotto |
| Leech lattice simple root norm | $2$; $r_\lambda = (\lambda; 1, 1 - \lambda^2/2) = (\lambda; 1, -1)$ at $\lambda^2 = 4$ | $6$ (Wave 24 Nekrasov formula error, using $r_\lambda = (\lambda; e + (\lambda^2 - 2)/2\,f)$) | Conway 1983 *Proc R Soc Lond A* 384; SPLAG Ch 27 |
| Fake-Monster Weyl vector | $\rho = e = (0_\Lambda; 1, 0)$ lightlike, $\rho_\Lambda = 0$ | Sign-inconsistent Wave-23 draft | Borcherds 1992 Thm 10.1 |
| $\Phi_{12}$ automorphic home | Orthogonal Shimura variety $\mathcal{D}_{\mathrm{II}_{26,2}} = O(26, 2)^+ / (O(26) \times O(2))$, complex dim 26, singular weight $12 = (26 - 2)/2$ | Siegel $\mathrm{Sp}_{26}(\mathbb Z)$ (impossible, $\mathbb H_{26}$ has dim 351); Jacobi-in-26-variables (no torus direction) | Borcherds 1992 |
| Restriction sublattice for $\Phi_{12} \to \Phi_{10}$ | Primitive $\mathrm{II}_{2,2} \hookrightarrow \mathrm{II}_{25,1}$, giving $O^+(\mathrm{II}_{2,2}) \cong \mathrm{Sp}_4(\mathbb Z)$ | $\mathrm{II}_{2,1}$ (Wave 23 error — sig $(2,1)$ carries no holomorphic modular form) | Gritsenko--Nikulin 1998 |
| Master $L$-value identification | $\log Z^{(1)}_{\mathbf H_{\Delta_5}} = -\log \Delta_5 - \kappa_{\mathrm{BGS}} \cdot L'(0, \Delta_5, \mathrm{std}) + \log C$, $\kappa_{\mathrm{BGS}} = 24$ | $L'(0, \Delta_{10}, \mathrm{ad}^0)$ (Wave 24 Costello — CAP-reducible, conflates three $L$-functions); $L'(0, \mathrm{ad}^0 \rho_{\Delta_{12}})$ (Wave 25 Gaiotto agent alternative, incompatible regulator) | Bruinier--Kühn 2003 Thm 4.11; Yoshikawa 2004 Thm 5.7 |
| Bloch--Kato Selmer on representation of $\mathbf H_{\Delta_5}$ deformation | $\dim H^1_f(\mathbb Q, \mathrm{std}\,\rho_{\Delta_5}) = 1$ (paramodular cyclotomic Hida family tangent) | $\dim H^1_f(\mathrm{ad}^0 \rho_{\Delta_{10}}) = 1$ (Wave 24 — right number, wrong representation; the adjoint spinor is rigid, dim 0) | Pilloni 2011 + Urban 2011 + Poor--Yuen 2015 + Thorne 2020 |
| Bridgeland $\dim \mathrm{Stab}(K3 \times E)$ | $48 = \mathrm{rk}\,\mathcal N(K3 \times E) = 24 \cdot 2$ | $26$ (Wave 23 claimed "codim-0 slice" — actual Künneth image has codim $22 = \mathrm{rk}\,T_{K3}$) | Bridgeland 2007 Thm 1.2 |
| Codim of $\mathrm{Stab}^\Phi$ inside ambient | $22$ | $0$ (Wave 23 error) | via four independent paths |
| Rank of $u_{\zeta_8}^{\mathrm{tilt}}$-mod | $162 = 27 \cdot 6 = (\ell' - 1)^3 \cdot \lvert S_3 \rvert$ | --- | Andersen--Polo--Wen 1994 |
| Structure of rank-$162$ MTC | $(A_1)_{k=2}^{\otimes 3} \rtimes S_3$ ($S_3$-crossed braided fusion) | $(A_1)_{k=2}^{\otimes 3} \boxtimes \mathbb Z[S_3]$ (Wave 23 tensor claim — $\mathrm{Vec}_{S_3}$ non-modular) | Turaev 2000; ENO 2010 |
| Conway $V^{s\natural}$ central charge | $c = 12$ | $c = 24$ (Wave 23 inherited from Monster error) | Duncan 2007 *Duke Math J* 139 |
| Conway $V^{s\natural}$ $\Psi$-placement | $\Psi^{\mathrm{metap}}$-image on metaplectic $\overline{\mathcal A_2^{(2)}}$ branch | Bosonic $\Psi$-image with $(K, \hbar^2) = (2, -1/2)$ (Wave 19/23 — Leech has no hyperbolic plane, universal identity out of scope) | Scheithauer 2008 |
| "Four is all" citation chain | Scheithauer 2017 + Dittmann--Ma--Scheithauer 2021 + Scheithauer 2006 (three papers) | Scheithauer 2017 alone (Wave 23/24 incomplete) | Scope: GN-reflective signature-$(2, n)$ with $n \ge 3$ |
| Fifth Borcherds product outside GN-scope | $24 A_1$ Niemeier product (Borcherds 1995 *Invent Math* 120 §13), reflective automorphic on sig-$(2, 24)$ singular weight 12 but fails GN-reflectivity (divisor non-rational-quadratic-hyperplane components) | --- | Borcherds 1995 |
| Pentagon admissibility congruence variable | $D_n = (n - 3)/2 \pmod 4 \in \{0, 1\}$, equivalently $n \equiv 3, 5 \pmod 8$ | $n \not\equiv 0, 3 \pmod 4$ (loose restatement) | Eichler--Zagier 1985 Thm 9.1 |
| First admissible non-vanishing $\phi^{(n)}$ | $\phi^{(5)} = -2 \cdot [\mathrm{gen}]^{\otimes 5}$ | --- | Gritsenko--Nikulin 1998 Table 2 |
| $\Psi$-functor completeness | **Four** sibling functors $\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$ (disjoint union surjective) | $\Psi$ alone surjective (Wave 23 claim falsified) | Baily--Borel--Freitag stratification of $\overline{\mathcal A_2}$ |
| $\mathrm{grt}_1^{(1/2)}$ extension status | Non-split, non-abelian, non-central (three distinct properties) | "Split because central" (Wave 24 conflation) | Lie $H^2_{\mathrm{Lie}}(\mathfrak q; \mathrm{grt}_1)$ via van Est transgression |
| Obstruction cohomology type | Lie $H^2_{\mathrm{Lie}}$ via van Est transgression from group $H^1_{\mathrm{grp}}$ (Saito--Kurokawa Eichler cocycle) | Group $H^1$ only (Wave 24 misattribution) | Costello--Gwilliam 2017 BV realisation |
| $[\widetilde\sigma_2, \widetilde\sigma_2]$ | $288 \widetilde\sigma_4$ via $\tau(2)^2 / 2 = 288$ | --- | $\tau(2) = -24$ |
| $e_5$ vs $W_5$ | $e_5 = W_5$ identically at every $c$ | --- | Wang 1998 *Prog Theor Phys* Prop 4.2 three-leg uniqueness |
| $e_4$ at $c = -214$ | $:T \partial^2 T: - (3/2):(\partial T)^2: + (321/10) \partial^4 T + \hbar\mathrm{qt}(J^{(4)})$; pairing $(65193/10) \mathrm{Vol}(E) (2\pi i)^4$ | --- | Pope--Romans--Shen 1990; Wang 1998 |
| Chiral-Hochschild $e_k$ motivic home | Single-valued $\mathrm{zv}^{\mathrm{sv}}_{3k}$ (Brown 2014; Schnetz 2013) | Full motivic $\mathrm{Per}^{\mathrm{mot}}_{3k}$ (Wave 23 naive Costello cross-cut) | BGS analytic torsion on Shimura varieties forces SV |
| Pseudo-character target | Chenevier 2014 determinant $D_{\Delta_{10}}$ of dim 4 | Pseudo-representation (loses mod-$\ell^n$ Cayley--Hamilton for reducible $\rho_{\Delta_{10}}$) | Chenevier 2014 *Ann Inst Fourier* 64 |
| Hecke field of $\Delta_{10}$ | $\mathbb Q$ (since $\dim S_{26}(\mathrm{SL}_2(\mathbb Z)) = 1$); minimal coefficient ring $\mathbb Z$ | --- | Standard |
| Humbert divisor / AD correspondence | Humbert divisor $= $ Argyres--Douglas point at $E_{\tau_1} \times E_{\tau_2}$ (SW curve degenerates to pair of elliptic curves → $(A_1, A_{2N-1})$) | --- | GMN 2009 *Adv Theor Math Phys* 13 Ex 8.3 |
| $\kappa_\bullet$ indexing (K3 $\times$ E) | Four distinct values: $\kappa_{\mathrm{cat}} = 0$ (multiplicative Künneth), $\kappa_{\mathrm{ch}}^K = 3$ (additive), $\kappa_{\mathrm{BKM}}$ family-specific, $\kappa_{\mathrm{fibre}}(K3) = 2$ | Naive "$\kappa_{\mathrm{BKM}} = \kappa_{\mathrm{ch}} + \chi(\mathcal O_{\mathrm{fibre}})$" (fails for $N \ge 2$ — $N = 1$ accident only) | Always name the $\kappa_\bullet$ index |
| $\kappa_{\mathrm{BKM}}(\mathbf H_{\Delta_5})$ cross-volume value | **Pending AP5 lock** — Vol I abstract says 12, Vol III abstract says 5 (different $N$-index conventions for $c_N(0)/2$); every site must name the input denominator (Fake-Monster $\Phi_{12}$ vs paramodular $\Phi_{10} = \Delta_5^2$) | Both values occur in published inscriptions; resolve via landscape-census audit | Open AP5 audit |
| CoHA vs chiral-algebra type | $\mathrm{CoHA}_{K3 \times E}$ is $E_1$-associative (Hall product); chiralisation via $\Phi_3$-arrow gives $\mathbf H_{\Delta_5} = \Phi_3(\mathcal D_\hbar(\mathcal Y^{\mathrm{Hall}}(\mathrm{CoHA}_{K3 \times E})))$ | "CoHA is a chiral algebra" (type error) | Schiffmann--Vasserot + $\Phi_d$-functor framework |
| $\mathrm{CoHA}(\mathbb C^3)$ identification | $Y^+(\widehat{\mathfrak{gl}}_1)$ — the **positive half** of the affine Yangian | $W_{1+\infty}$ (classical limit only); full Yangian $Y$ | Schiffmann--Vasserot |
| $W_{1+\infty}$ vs $W_\infty[c]$ | $W_{1+\infty} = W_\infty[c] \otimes \mathcal H$ (Heisenberg); different objects | Conflating the two | Pope--Romans--Shen 1990 |
| Two-$\hbar$ discipline | $\hbar^{\mathrm{Drinfeld}} = 2\pi i / \ell$ (root-of-unity) vs $\hbar^{\mathrm{BV}}$ (loop-counting); semantically distinct but agree numerically at $\hbar^2 = -1/8$ for $\ell = 8$ | Using bare $\hbar$ without subscript | AP151 bridge |
| Bar cohomology class-$\mathbf M$ at $E_3$-level | $6^g$ (cohomological dim at $g \in \{1, 2, 3\}$; chain level is infinite) | "Infinite at cohomological level" | Wave 12 Vol III |
| $\kappa_{\mathrm{cat}}(K3 \times E)$ | $0$ (total space, Künneth multiplicative $\chi(\mathcal O_{K3}) \chi(\mathcal O_E) = 2 \cdot 0$) | $2$ (fiber only) | AP-CY68/AP234 discipline |

**Every cache/AP entry citing an older value must appear in a "Wrong
Claim" / "retracted" / "earlier draft" context. Any standalone assertion
of an older value is itself a latent AP violation; report it via
Gate 0 on sight.**


## Compatible dual readings (non-contradictions)

Some statements appear in different waves with apparently conflicting
framings yet describe the **same mathematical object** through two
compatible frames. An agent reading only one framing may mistake the
pair for a contradiction and attempt to "fix" one side, destroying
content. Every entry below is a certified **non-contradiction**: both
framings are correct simultaneously and reconcile through the named
mechanism. Do not flatten either reading. Do not treat this as an AP
retraction list; treat it as a frame-compatibility registry.

| Statement A (frame 1) | Statement B (frame 2) | Compatibility mechanism | Primary reference |
|---|---|---|---|
| Conway $V^{s\natural}$ is the $\Psi^{\mathrm{metap}}$-image on the metaplectic $\overline{\mathcal A_2^{(2)}}$ branch (canonical preamble, row 41). | Conway $V^{s\natural}$ is the $\mathbb Z/2$-super-twin of $V^\natural$ via Duncan's commutative orbifolding diamond (W20.2 healing). | Same VOA-super-algebra read from the $\Psi$-functor / $\overline{\mathcal A_2}$-stratification side and from the commutative-orbifolding-diamond side; the metaplectic $\Psi$-image equals the super-twin of the Leech CFT via the Duncan 2007 Frame-change. | Duncan 2007 *Duke Math J* 139; Scheithauer 2008. |
| $\kappa_{\mathrm{BKM}}(\mathbf H_{\Delta_5}) = 12$ (Vol I three-faces identity, Fake-Monster $\Phi_{12}$ input). | $\kappa_{\mathrm{BKM}}(\mathbf H_{\Delta_5}) = 5$ (Vol III Borcherds-weight count on paramodular $\Delta_5$). | Reconciled via W26.8: the two values index different input denominators ($\Phi_{12}$ vs $\Phi_{10} = \Delta_5^2$) under the same $c_N(0)/2$ rule. Each site must name the input denominator; the two are $(c_{\Phi_{12}}(0)/2, c_{\Delta_5}(0)/2) = (12, 5)$. | Borcherds 1992; Gritsenko 1999 Thm 6.1; canonical preamble row 57. |
| $e_4 = W_4 - (107/11)\Lambda_Z$ in the $W$-basis (Vol III abstract, $c = -214$). | $e_4 = {:}T\partial^2 T{:} - (3/2){:}(\partial T)^2{:} + (321/10)\partial^4 T + \hbar\,\mathrm{qt}(J^{(4)})$ in the Virasoro-composite basis (canonical preamble row 51). | Same cocycle class, two bases. Pope--Romans--Shen projection $W_k = e_k + (-c/22)\Lambda_Z$ converts between $W$-basis and Virasoro-composite basis; the scalar $-107/11 = c/22$ at $c = -214$ is the frame-change constant. | Pope--Romans--Shen 1990; Wang 1998 Prop 4.2. |
| Archetype $\mathbf B$ row: $\kappa + \kappa^! = 8$ (Vol I five-archetype table, $r(z)$-families). | K3 Mukai pairing $K = 2c_+ = 8$ (Vol III Mukai form on $\Lambda^{3,19}$). | Three faces of the same canonical $K = 8$: (i) derived-centre complementarity sum $\kappa + \kappa^!$; (ii) Mukai pairing $K = 2c_+$ with $c_+ = 4$; (iii) Lusztig reflection length $\ell_{\mathrm{Lusztig}} = 8$ at $\zeta_8$. All three equal $8$ because $K3 \times \mathrm{elliptic}$ pairs with $\mathbf H_{\Delta_5}$ through a single self-dual lattice indexing. | Borcherds--Mukai pairing; Lusztig 1989 at $\zeta_8$; Gritsenko--Nikulin 1998. |
| Padovan counts $d_n = d_{n-2} + d_{n-3}$ give MZV basis dimensions in weight $n$ (Brown 2012 conjecture). | Humbert-Heegner divisor filter selects arithmetically admissible $n \equiv 3, 5 \pmod 8$ on $\phi^{(n)}$ (Vol I / Vol III pentagon admissibility row). | Orthogonal filters. Padovan counts basis vectors; Humbert-Heegner selects admissible embedding points. The two do not compete: Padovan $d_n$ governs the motivic dimension; Humbert-Heegner governs the arithmetic non-vanishing locus. Compatible because $\phi^{(n)}$ lives in a weight-$n$ MZV-period space filtered further by admissibility. | Brown 2012; Eichler--Zagier 1985 Thm 9.1 (canonical row 44). |
| $(\infty, 1)$-obstruction tower (W21.2) gives the Maurer--Cartan tower of universal obstructions in an $(\infty, 1)$-stable category. | $(\infty, 2)$-adjunction AR structure (W22.1) exhibits bar--cobar as an adjunction in an $(\infty, 2)$-category with an Auslander--Reiten sequence. | The $(\infty, 2)$ structure is the $2$-categorical lift of the $(\infty, 1)$ obstruction tower: AR sequence at object level descends to the MC obstruction tower at morphism level. Both proved, both load-bearing. | Lurie HA 5.5; Riehl--Verity 2020. |
| Six-path $\chi_3$ (W22.6): six independent verification paths for the third Chern character. | Seven-path $\chi_3$ (W25.4): seven paths including a derived-deformation check. | Six independent paths plus one additional derived-deformation path, consistent with prior comparison. The seventh path refines the six without contradicting them. Every path reports the same $\chi_3$ value. | Beilinson multi-path verification discipline; W25.4 cache entry. |
| Four sibling functors $\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$ parametrise the Baily--Borel--Freitag stratification of $\overline{\mathcal A_2}$ (canonical preamble row 46). | Five $\Psi$-image rows including the Conway $V^{s\natural}$ row on the metaplectic branch (W20.2 table). | The four siblings parametrise the four strata; each sibling has one image row, and the $\Psi^{\mathrm{metap}}$ sibling hosts two rows (Conway $V^{s\natural}$ and a second metaplectic image). Row count $5 = 4$ siblings with one extra Conway row on the metaplectic branch. | Baily--Borel--Freitag; Duncan 2007. |
| Absolute Kuznetsov HPD blocked on $K3 \times E$ (obstruction from the non-trivial Brauer class on $K3 \times E$). | Relative HPD over $E$ works (the relative $E$-base kills the obstruction by descent). | Obstruction at the absolute level; healing at the relative level. Both statements are proved theorems; relative HPD is the correct home for the $K3 \times E$ Kuznetsov lift. | Kuznetsov 2007; Perry 2019 relative HPD. |
| Pseudo-character framework (Taylor 1991) is sufficient on reduced rings for Galois deformations of $\rho_{\Delta_{10}}$. | Chenevier determinant framework (2014) is necessary on non-reduced rings (where mod-$\ell^n$ Cayley--Hamilton fails for reducible $\rho$). | Compatible hierarchy, not competing frameworks. Taylor's pseudo-character suffices on reduced; Chenevier's determinant is the correct lift to non-reduced. Each dominates in its natural domain; neither refutes the other. Canonical preamble row 53 pins the non-reduced case to Chenevier; reduced-ring calculations may still cite Taylor. | Taylor 1991 *Invent Math* 116; Chenevier 2014 *Ann Inst Fourier* 64. |

**Operating rule.** If you encounter two inscriptions that seem to
contradict, before editing either side **check this table**. If the
pair matches an entry, both are correct. If not, escalate through
Gate 0: the confusion is a candidate for a new compatibility mechanism,
not a candidate for retraction.


## Open Problem Index (AP-CAT-3 resolution, 2026-04-21)

Each active open problem in the frontier is pinned to the AP or
cache-pattern chain that bounds its failure modes. A reader asking
"what is the AP status of Open Problem $k$?" finds the full chain
in one row below. Update this block on every inscription of a new
open-problem entry. The chain lists (a) the cache pattern giving
the first-principles discipline; (b) the canonical-values row that
fixes the latest-wave verdict; (c) the retraction AP that closed an
earlier wrong claim (where applicable).

| OP # | Open problem | Primary AP chain | Cache / pattern anchor | Canonical-value row |
|---|---|---|---|---|
| OP1 | Humbert--Heegner admissibility filter $n \equiv 3, 5 \pmod 8$ on $\phi^{(n)}$ | AP-V2-24 / V2-AP127 | Pattern 299 / entry 410 (AP890) | canonical row 44 (Eichler--Zagier 1985 Thm 9.1) |
| OP2 | Fake-Monster Cartan rank on $\mathrm{II}_{25,1}$ | AP-CY71 retraction of "Conway $V^{s\natural}$ as 5th $\Psi$-image" + AP892 metaplectic restoration | entry 412 (AP892) | canonical rows 20--21 (Monster rank 2 on $\mathrm{II}_{1,1}$; Fake-Monster rank 26 on $\mathrm{II}_{25,1}$) |
| OP3 | Universal ratio-of-levels closure at the Leech-Conway exception | cache entry 177 (Wave 18 DRINFELD) | first-principles entry 177 | canonical row 40 (Conway $V^{s\natural}$ $c = 12$) |
| OP4 | Bridgeland $\dim \mathrm{Stab}(K3 \times E) = 48$, codim of $\mathrm{Stab}^\Phi = 22$ | AP-CY retraction of "codim-0 slice" claim | cache entries 178, 192 | canonical rows 36--37 (Bridgeland 2007 Thm 1.2; four independent paths) |
| OP5 | Pseudo-character vs Chenevier determinant on non-reduced deformation rings | AP-CY87 (Wave 25 rename) + AP-CAT-5 non-contradictions cross-reference | Pattern 295 / entry 422 (AP902) | canonical row 34--35 (paramodular Hida family; Selmer dim 1) |
| OP6 | Yetter--Drinfeld tower weight $\delta^{(n)} = \lfloor n/2 \rfloor + 1$ | cache entry 175 (Wave 18 DRINFELD correction of Wave 17 $\lceil n/2 \rceil$) | first-principles entry 175 | retraction row (Wave 17 partial envelope scaling) |
| OP7 | Enriques $M_{12}$ mass formula threshold $D_0 = 0$ with Mersenne exception | cache entry 188 (Wave 21 BEILINSON extension) + cache entry 170 (Wave 20 BEILINSON template-mismatch falsification) | first-principles entries 170, 171, 188 | --- |
| OP8 | $\kappa_{\mathrm{BKM}}$ Vol I (Fake-Monster $\Phi_{12}$ input, value 12) vs Vol III (paramodular $\Delta_5$ input, value 5) | AP-CY85 + AP-CY69 + AP5 dual-indexing header (canonical row 57) + Pattern 411 (formerly Pattern 226; cross-volume bibkey check, renumbered 2026-04-21 per AP-CAT-1) | compatible-dual-readings row 3 (canonical preamble row 57) | canonical row 57 ($c_N(0)/2$ rule; $(c_{\Phi_{12}}(0)/2, c_{\Delta_5}(0)/2) = (12, 5)$) |
| OP9 | $\chi_3$ six-path (W22.6) vs seven-path (W25.4) verification closure | compatible-dual-readings row 9 | entry 192 (six paths) + W25.4 cache entry | --- |
| OP10 | $(\infty, 2)$-adjunction AR structure as $2$-categorical lift of $(\infty, 1)$ obstruction tower | compatible-dual-readings row 10 (W22.1 + W21.2) | Lurie HA 5.5; Riehl--Verity 2020 | --- |
| OP11 | Absolute vs relative Kuznetsov HPD on $K3 \times E$ | compatible-dual-readings row 11 (Perry 2019 relative HPD) | Kuznetsov 2007; Perry 2019 | --- |
| OP12 | Leech-Conway universal ratio-of-levels exception | cache entry 177 (Wave 18 DRINFELD) | first-principles entry 177 | canonical row 40 |
| OP13 | Humbert transversality at admissible Heegner triples (octachotomy) | AP304 concurrent-agent collision discipline; Wave 22.2 falsification of heptachotomy | cache entry (Wave 22.2 octachotomy lift) | --- |
| OP14 | Universal celestial holography proof-body closure | V2-AP conjecture-promotion discipline | cache `conj:uch-*` registry | --- |

**Operating rule.** Before inscribing a new open-problem entry in
any theorem-status table or FRONTIER row, append a row to this index
naming (a) the AP chain bounding the failure mode; (b) the cache or
first-principles entry giving the primary-lit discipline; (c) the
canonical-values row pinning the latest verdict (where applicable).
An open problem without an AP-chain anchor is a latent AP-CAT-3
violation; inscribe the anchor at the session close.


---

carries only a pointer to this file. Every agent running
/chriss-ginzburg-rectify must consult this catalogue AND
`appendices/first_principles_cache.md` at Gate 0 of every chunk.

---

## Cross-Volume Anti-Patterns
Before cross-volume edits, Read `notes/cross_volume_aps.md` (Vol II V2-AP* and Vol III AP-CY1..AP-CY69 catalogs). The Geometric/Algebraic Model Conflations (AP-CY62..AP-CY67) are kept inline below under "Geometric vs Algebraic Models." Künneth-multiplicative + Route A/B collision (AP-CY68..AP-CY69) live in `notes/cross_volume_aps.md`.

## Consolidated / Merged Anti-Patterns

**AP-RMATRIX (= AP126 + AP141 + AP148).** Rule: every r-matrix carries level prefix; verify k=0 in the chosen convention. Counter: (a) substitute k=0 and verify r vanishes (trace-form) or gives Omega/(h^v*z) (KZ non-abelian); (b) state convention; (c) grep bare `\Omega/z`. Two conventions for affine KM: trace-form `r(z)=k*Omega/z` (AP126 k=0 check; av(r)=kappa_dp; Sugawara shift dim(g)/2 for full kappa) and KZ `r(z)=Omega/((k+h^v)*z)` (k=0 nonzero for non-abelian). Bridge: `k*Omega_tr = Omega/(k+h^v)` at generic k. 42+ instances; THE MOST VIOLATED AP.

**AP-KAPPA (= AP1 + AP9 + AP20 + AP24 + AP48 + AP136).** Rule: kappa DISTINCT per family. Always qualify (`kappa^{KM}`, `kappa^{Vir}`). Counter: Read landscape_census.tex; evaluate at k=0, k=-h^v (KM), c=13 (Vir); cross-check compute/. Writing from memory FORBIDDEN. Values: KM=dim(g)(k+h^v)/(2h^v); Vir=c/2; W_N=c*(H_N-1), H_N=sum_{j=1}^{N} 1/j (NOT c*H_{N-1}; at N=2, H_1=1 but H_2-1=1/2); Heis=k. Complementarity: 0 (KM/free), 13 (Vir), 250/3 (W_3), 196 (BP). State WHICH algebra: intrinsic vs kappa_eff=kappa(matter)+kappa(ghost) vs kappa(B=A^!).

**AP-DESUSP (= AP22 + AP45 + B15 + B16).** Rule: `|s^{-1}v|=|v|-1`. Counter: `T^c(s^{-1} A-bar)`, NOT `T^c(s A-bar)`; bar=down=desuspension=s^{-1}.

**AP-SC-BAR (= AP165 + FM25 + B54 + B55 + B56).** Rule: B(A) = T^c(s^{-1} A-bar) is a SINGLE E_1 chiral coassociative coalgebra. SC^{ch,top} lives on the derived center pair (C^bullet_{ch}(A,A), A), NOT on B(A). Regex: `bar.*presents.*Swiss|coalgebra over.*\\SCchtop|\\SCchtop.*coalgebra`. The 6-step chain of errors (FM25): (1) "bar differential is closed color" WRONG — d_B is single degree-1 map over (ChirAss)^!, not parametrized. (2) "deconcat = open color" WRONG — cofree coassociative on T^c(V); single-colored E_1. (3) "d_B + Delta = SC-coalgebra" WRONG — two structures != two colors. (4) "SC on B(A) dual to SC on (Z^{der}(A), A)" WRONG — B(A) is INPUT; SC emerges in OUTPUT via Hom(B(A),A). (5) "E_1-chiral = E_1-topological" WRONG — chiral is on CURVE (2-real-dim holomorphic); topological is Conf_k(R); bar is over (ChirAss)^!, not (Ass)^!. (6) Steinberg analogy RETIRED. FORBIDDEN: "B(A) is SC-coalgebra", "bar presents Swiss-cheese", "bar differential is closed color", "bar coproduct is open color", "Curved SC^{ch,top}-algebra" for g>=1 (genus>=1 bar is curved A_infinity-chiral; d^2_{fib}=kappa*omega_g is on A_infinity structure). Counter: after writing B(A) + SC^{ch,top} in same paragraph, verify SC attributed to derived center pair. C^bullet_{ch}(A,A) via End^{ch}_A (spectral params from FM_k(C)), NOT topological RHom_{A^e}(A,A).

**AP-TOPOLOGIZATION (= AP158 + AP167 + AP168 + B58 + B59).** Rule: topologization (SC^{ch,top} → E_3-TOPOLOGICAL via Sugawara) PROVED for affine KM V_k(g) at non-critical level k!=-h^v; CONJECTURAL for general chiral algebras with conformal vector (conj:topologization-general). Proof COHOMOLOGICAL (Q-cohomology), not chain-level; for class M chain-level E_3 may fail. Counter: every reference carries "(proved for affine KM at non-critical level; conjectural in general; cohomological, not chain-level)." SC is two-coloured with directionality; Dunn does NOT apply to coloured operads. E_3 requires SC^{ch,top} + conformal vector making C-translations Q-exact; conformal vector KILLS chiral direction at cohomological level, two colors collapse, E_2^{hol}×E_1^{top} → E_3^{top} via Dunn. Without conformal vector: stuck at SC^{ch,top}. At k=-h^v: Sugawara undefined. "Chiral" label FORBIDDEN for topologized bulk; write E_3-topological.

**AP-SC-NOT-SELFDUAL (= AP166 + FM26 + B57).** Rule: `(SC^{ch,top})^! != SC^{ch,top}`. SC^! = (Lie^c, Ass^c, shuffle-mixed), closed dim (n-1)!. SC = (Com, Ass, product-mixed), closed dim 1. Confusion: "duality FUNCTOR is involution" != "OPERAD is self-dual". (P^!)^!~P does NOT mean P^!~P. Livernet proves Koszulity NOT self-duality. Counter: verify dim P(n) = dim P^!(n) at all degrees before claiming self-dual. CORRECT: "SC Koszul duality exchanges Com↔Lie while preserving Ass; duality functor is involution."

**AP-UNIFORM-WEIGHT-TAG (= AP32 + B26).** Rule: every obs_g, F_g, lambda_g in theorem/remark/definition MUST carry explicit scope tag: (UNIFORM-WEIGHT), (ALL-WEIGHT + cross-channel correction delta F_g^cross), (g=1 only), (LOCAL: scope in surrounding paragraph). Untagged = violation.

**AP-LABEL-DISCIPLINE (= AP124 + AP125 + B34 + B35 + FM14 + FM15).** Rule: (i) prefix matches env (thm→theorem, prop→proposition, lem→lemma, conj→conjecture, rem→remark, def→definition); (ii) uniqueness across all three volumes. Counter: before any `\label{foo}`, grep all three volumes (`~/chiral-bar-cobar ~/chiral-bar-cobar-vol2 ~/calabi-yau-quantum-groups`); if duplicate, prefix v1/v2/v3. When upgrading/downgrading environment: atomic rename (env + label + all `\ref` instances) in SAME tool-call batch, NO intermediate commit.

**FM-LIE-NUMERICS (= FM5 + FM21 + C16).** Rule: never generate exceptional Lie dimensions from memory. Plausible fabrications (779247 for E_8 is not any irreducible); prefactors (1/2, 1/24, 1/(2*pi*i), 7/5760) frequently wrong. Counter: grep `compute/lib/` before writing; if no match, use symbolic name `V_{omega_1}(E_8)`. E_8 fundamentals: `{248, 3875, 30380, 147250, 2450240, 6696000, 146325270, 6899079264}`. Source: `compute/lib/bc_exceptional_categorical_zeta_engine.py::FUNDAMENTAL_DIMS['E8']`.

## Opus 4.6 Quirks and Failure Modes

Model-specific recurrent patterns. Full source: opus_46_failure_modes_wave12.md.

- FM1. Generic-formula reaching. Falls back to textbook form when uncertain. Counter: substitute k=0 and verify r vanishes before proceeding with any r-matrix.
- FM2. Level-prefix dropping on summarisation (omits k, c/2, kappa, 1/(2*pi*i)). Counter: re-Read source lines verbatim; do not rely on context cache.
- FM3. Bosonic/fermionic conformal-anomaly conflation (c_bg, c_bc swap). Counter: substitute lambda=2 AND lambda=1, verify c_bc(2)=-26, c_bg(2)=+26, sum=0.
- FM4. k=0 (abelian, still Koszul) vs k=-h^v (critical, not Koszul) confusion at KM. Counter: two-row table before any KM limit statement.
- FM5. See FM-LIE-NUMERICS.
- FM6. Undefined macros in standalone extraction (`\cW`, `\hol`, `\Ran`, `\FM`, `\chHoch` not inherited). Counter: grep `\\[a-zA-Z]+` against standalone preamble.
- FM7. LaTeX typo `\end{definition>`. Counter: after every .tex Edit, grep `\\end\{[^}]*>` and `\\begin\{[^}]*>`.
- FM8. Universal-quantifier drift on uniform-weight theorems. Counter: three-line template (scope, tag, equation) before any obs_g or F_g.
- FM9. Harmonic-number off-by-one (H_{N-1} vs H_N-1). Counter: evaluate at N=2 AND N=3 numerically.
- FM10. Hardcoded Part~IV/Chapter~12. Counter: after any Edit, grep `Part~[IVX]+|Chapter~[0-9]+` and replace with `\ref{part:...}`.
- FM11. Sugawara shift missing in av(r(z))=kappa. Abelian: av(r)=kappa; non-abelian KM: av(r)+dim(g)/2=kappa(V_k(g)). Counter: state family first.
- FM12. Mid-response truncation on long audit tasks. Counter: separate fix from report across two tool calls.
- FM13. Auto-completion to majority-variant (Dedekind eta as bare `prod(1-q^n)`, missing `q^{1/24}`). Counter: break formula across lines; annotate each term's convention.
- FM14/FM15. See AP-LABEL-DISCIPLINE.
- FM16. Silent kappa-family conflation. Counter: every kappa formula carries explicit family superscript.
- FM17. Amplitude/dimension conflation for ChirHoch. "amplitude [0,2]" != "virtual dimension 2". Counter: choose explicitly.
- FM18. d_gen vs d_alg conflation. d_gen(Vir)=3, d_alg(Vir)=inf. Counter: subscript mandatory; refuse bare `d(...)`.
- FM19. Fiber-base confusion. dτ is on Σ_g; omega_g=c_1(lambda) is on M̄_g. Counter: state space for each side.
- FM20. Iff-drift. Counter: two-line template "Forward:... Reverse:... proved/CONJECTURED" before `\iff`.
- FM21. See FM-LIE-NUMERICS.
- FM22. Koszul conductor K_BP=2 instead of 196. Counter: write K=c+c', substitute two central charges, evaluate.
- FM23. Local-global conflation "over a point = over P^1". Three errors: (1) homotopy retract is DATA, (2) D != point, (3) A^1 already has Arnold relations. Counter: name specific space (point/D/A^1/P^1/X), comparison data, on-the-nose-or-extra-structure. See AP142.
- FM24. B-cycle monodromy i^2 error. hbar with pi*i factor: q=e^{2*pi*i*hbar} becomes real, not root of unity. Counter: substitute integer k; verify |q|=1.
- FM25. See AP-SC-BAR.
- FM26. See AP-SC-NOT-SELFDUAL.
- FM27. Scope inflation in metadata. Counter: scope qualifiers match actual verification level.
- FM28. See AP-TOPOLOGIZATION.
- FM29. Sign convention d+A vs d-A inconsistency. nabla=d+A gives dPhi=-A*Phi; physics nabla=d-A gives dPhi=+A*Phi. Counter: verify nabla(Phi)=0 gives displayed flat section.
- FM30. Belavin r-matrix Weierstrass (WRONG: `zeta(z)=theta_1'/theta_1 + 2*eta_1*z` breaks CYBE) vs Pauli decomposition `r(z)=sum w_a sigma_a⊗sigma_a/2`. Degeneration is TWO-STEP: elliptic→trigonometric (XXZ, Im tau→inf)→rational (z→0). Counter: verify CYBE numerically.
- FM31. Miura coefficient universality. (Psi-1)/Psi on cross-terms in Δ_z is UNIVERSAL across spins (J⊗J at s=2; J⊗T+T⊗J at s=3). Composite corrections appear; primary coefficient persists. Counter: verify at leading cross-term.
- FM31a. Asymptotic cancellation. 10/(5c^2)=2/c^2, NOT 2/(5c^2). Counter: substitute c=100 against exact formula.
- FM32. pi_k(BU): use fiber sequence pi_k(BX)=pi_{k-1}(X); check Bott periodicity parity. pi_3(BU)=0.
- FM32a. RTT sign depends on R-matrix convention. Additive R(u)=uI+Psi*P: [t_ij,t_kl]=Psi(delta_il t_kj - delta_kj t_il). Molev's 1-P/u: opposite sign. Counter: state convention.
- FM33. Formula applied outside hypothesis domain. kappa_ch=chi(S)/2 valid only for Tot(K_S→S); conifold NOT local surface. Counter: verify input satisfies hypotheses before applying.
- FM33a. qdet SPECTRAL-SHIFT PATTERN load-bearing, product order is not. Column `j` carries shift `u-j*Psi`; this is central independently of product order (reversed column-product `col_inc` equals `col_dec` as element of Y(gl_N), both central). Row-determinant or shift-reversed variant (leftmost factor has shift `(N-1)*Psi`) fails centrality at order Psi^3 (NOT Psi^2 — log-log slope converges to 3 as Psi → 0 at N=3, exact-Fraction first-principles verification, adversarial_swarm_20260418/attack_heal_qdet_ordering_20260418.md). Agrees at N=2 by rank coincidence (codim-1 antisymmetriser image). Cite Molev 1.6.4 + 1.11.2.
- FM34. Excision=gluing (B_L⊗_A B_R=B(A), one copy over A); coproduct=splitting (B(A)→B(A)⊗B(A), plain). Counter: verify codomain.
- FM34a. Heat equation prefactor: 1/(4πi) diagonal (a=b), 1/(2πi) off-diagonal. Symmetric-matrix chain rule gives factor 2.
- FM35. NEVER REVERT MATHEMATICAL CONTENT TO FIX BUILD ERRORS (CONSTITUTIONAL). Build errors are LaTeX; math is never at fault. Counter: (a) identify specific LaTeX error, (b) fix it, (c) rebuild. 100 macro errors = add 100 `\providecommand` lines. NEVER DROP MATHEMATICS.
- FM36. Macro portability on cross-volume insertion. Counter: grep undefined control sequences after insertion; add `\providecommand`. Related V2-AP39.
- FM37. Double superscript from macro with built-in superscript (`\SCchtop`). Use explicit `\mathsf{SC}^{ch,top,an}`.
- FM38. Vertex-IRF not automatic. If vertex-DYBE fails, attack from IRF/face-model side directly. Genus-2 DDYBE face-model numerical evidence: 4 generic-Omega tests at relative < 1e-4 + 1 diagonal-Omega tier-(T3) consistency check; 30 face-model tests total across TestThetaFunctions, TestBoltzmannWeightsG1, TestFaceDYBEG1, TestClassicalLimit, TestGenus2Theta, TestBoltzmannWeightsG2, TestFaceDDYBEG2, TestFullSuite. conj:g2-ddybe remains ClaimStatusConjectured.
- FM39. Spectral coassociativity uses SHIFTED parameters. (Δ_{z1}⊗id)∘Δ_{z1+z2}=(id⊗Δ_{z2})∘Δ_{z1}. NOT Δ_z with itself.
- FM40. Naive center != derived center. Z(Drin(H_k)) dim 1 vs Z^{der}_{ch}(H_k) dim 3. Ext^1,2 invisible to commutant. Specify: commutant, Hochschild, or categorical.
- FM41. Jones polynomial = Markov trace + writhe normalization + quantum dim. Raw KZ trace != Jones.
- FM42. Bulk substring replacement corruption (CRITICAL). `replace_all` "arity"→"degree" corrupts singularity→singuldegree, parity→pdegree, etc. 45 corruptions in Vol III campaign. Counter: (1) NEVER bulk-replace short strings inside common words; (2) after any bulk replace, grep `ldegree|ndegree|rdegree|pdegree|tdegree`; (3) checklist {singularity, complementarity, unitarity, regularity, modularity, parity, familiarity, similarity, polarity, disparity, linearity, popularity, circularity, hilarity}.
- FM43. E_n output scope of the CY-to-chiral correspondence programme {Φ_d}: Φ_1 lands in E_∞-chiral (d=1), Φ_2 lands in E_2-chiral (d=2), Φ_d lands in E_1-chiral (d≥3). The target E_{n(d)}-ChirAlg(M_d) varies with d, so {Φ_d} is NOT a single functor (AP247 + AP244); it is a d-indexed family of per-d constructions. Counter: always state per-d scope `(n(d)=∞ at d=1; n(d)=2 at d=2; n(d)=1 at d≥3)` AND write Φ_d (with d subscript) rather than bare Φ unless the d=2 evaluation is contextually pinned. Canonical home: Vol III `chapters/theory/cy_to_chiral.tex` thm:phi-platonic + rem:phi-not-unified-functor + conj:phi-d-functoriality.
- FM44. Agent rate limiting. Counter: batches of 3, not 30+. Expect 5-13 min per agent on 1000-3000 line files.
- FM45. Agent skill fidelity gap (200-word brief vs 15,000-word skill). Counter: for full rectification, invoke skill in main conversation.
- FM46. Stale preface line counts. Counter: after campaign, update with `wc -l` comparison.
- AP186 (Opus). Coincidental agreement masks bugs. (Psi-1)/Psi=1/Psi at Psi=2; comb(d+2,2)=comb(d+2,3) at d=3. Verify at 3+ parameter values.
- AP187 (Opus). Miura coefficients from elementary symmetric expansion. T(u)=prod(u+Λ_i) → ψ_s=e_s(Λ_i). TWO-LAYER DISCIPLINE (AP541 rectification 2026-04-18): the Miura-intermediate coefficient of :J*W_{s-1}: inside ψ_s is 1/Psi (Step 1 of thm:miura-cross-universality). The headline theorem coefficient on the primary cross-term J⊗W_{s-1}+W_{s-1}⊗J in Δ_z(W_s) is (Ψ-1)/Ψ = 1 − 1/Ψ (assembled across Steps 1-3: Drinfeld binomial contributes +1, Miura inversion contributes −1/Ψ, lower J-sectors contribute 0). Never cite thm:miura-cross-universality with the bare "1/Ψ" label — that is the intermediate, not the theorem. Cross-check against B62 (Δ_z(T) cross-term (Ψ-1)/Ψ, not 1/Ψ) and rem:spin2-ceff-miura (standalone/ordered_chiral_homology.tex:3182).

## Theorem Status

**2026-04-17 Beilinson-rectified addendum**: the status table below reflects the 2026-04-16 closure wave inscriptions. Subsequent audit findings (Vol~I `notes/rectification_map_beilinson_audit.md`) refined the following claims: (i) Programme climax `thm:programme-climax` (Vol II) is SCOPE-QUALIFIED to non-logarithmic C_2-cofinite standard landscape + irrational cosets; logarithmic W(p) excluded pending Adamovi\'c-Milas character-amplitude bound (see Open Frontiers section). (ii) CY-C pentagon (Vol III) stratifies generator rank ρ^{R_i}, NOT κ_ch (Hodge-supertrace invariant = 0 route-independent). (iii) Kummer-irregular prime labels 1423, 3067, 23, 43, 419 retracted at primary-source verification. (iv) Super-Yangian max(m, n) complementarity identity replaces naive sum=0 analogy. (v) Topologization class M original-complex chapter carries explicit RETRACTION NOTICE; the 1/e obstruction is a Stirling cancellation error, unconditional heal at `thm:tempered-stratum-contains-virasoro` (Vol II). The theorems A-D+H remain PROVED at their stated scopes; only scope qualifiers and cross-reference labels changed.

| Thm | Status | Key result |
|-----|--------|------------|
| A | **PROVED unconditional on a FIXED SMOOTH CURVE X** (g=0 six/seven-fold TFAE + (∞,2)-properad equivalence on the conilpotent-complete locus under (H1)+(H2)+(H3), satisfied throughout the finitely generated standard landscape) via `thm:A-infinity-2` (Vol I `chapters/theory/theorem_A_infinity_2.tex`). **Modular-family extension over M̄_{g,n} including boundary is CONDITIONAL** on two uninscribed ingredients, both cited at citation level only: (a) Francis-Gaitsgory six-functor base-change on the relative Ran prestack (GR17 Ch. III §10), used to propagate the Hackney-Robertson properad-level equivalence horizontally across families of smooth curves; (b) Mok25 logarithmic factorization-gluing at the nodal boundary (Siao Chi Mok, arXiv:2503.17563, March 2025 v2 May 2025), to globalise the equivalence across the full compactified base. Neither (a) nor (b) is inscribed as a theorem in Vol I. The inscribed scope remark is `rem:A-infinity-2-modular-family-scope` at `theorem_A_infinity_2.tex:890`. The PTVV alternative previously advertised for Theorem A belongs to Theorem C. Downstream theorems requiring the modular-family extension — clutching-uniqueness for Theorem D at g≥2, universal obs_g = κ·λ_g across all g — inherit this conditional qualifier. See `adversarial_swarm_20260417/wave1_theorem_A_attack_heal.md` for the full ledger. | Genus-0 clause: hub-and-spoke TFAE (`thm:ftm-seven-fold-tfae-via-hub-spoke`, Vol I `chapters/theory/ftm_seven_fold_tfae_platonic.tex`): on the full Koszul locus, six bidirections (Koszul morphism, counit qi, unit weak eq, twisted tensor acyclic, bar concentrated in weight 1, plus the PBW E_2-collapse HUB); on the class-G stratum, a seventh equivalence (SC-formality) is added (class L / class M are Koszul but NOT SC-formal). Spoke 4 (twisted tensor ⇔ PBW) is the load-bearing non-tautology via `lem:filtered-comparison` (existing, proved); the Kac-Moody-specific V_k(sl_2) witness clause (iii) in `prop:class-L-witness` depends on the phantom `rem:kac-moody-filtered-comparison` and is an open inscription target (`adversarial_swarm_20260417/wave1_theorem_A_attack_heal.md` item OF2). Spoke 5 (bar concentration weight 1) extends to g≥1 only for class G (d²_fib = κ·ω_g obstruction off class G). Prior ghost-label `thm:FTM-seven-fold` (never defined) retargeted at `theorem_A_infinity_2.tex:570`. **Theorem A^{∞,2}** (Francis-Gaitsgory bar-cobar (∞,2)-equivalence at properad level, `thm:A-infinity-2`) is PROVED under (H1)+(H2)+(H3) on a FIXED smooth curve X; the R-twisted Σ_n-descent clause (iii) is conditional on unitarity of R(z) (Yangian rational, U_q trigonometric at generic q satisfy it). **E_1-ordered variant** `thm:theorem-A-E1` PROVED via pure braid Orlik-Solomon Koszulity (Shelton-Yuzvinsky); Yangian instance concrete (Cherednik monodromy + EK flatness). **Wave-18 attack-and-heal (2026-04-18, `attack_heal_theorem_A_20260418.md`)**: OF2 HEALED — `rem:kac-moody-filtered-comparison` inscribed at `ftm_seven_fold_tfae_platonic.tex:673-729` with Kashiwara filtration on vacuum + Kac-Kazhdan determinant. OF3 HEALED — unitarity hypothesis (R2) `R(z)R^{op}(-z) = id` explicit in `lem:R-twisted-descent` (`theorem_A_infinity_2.tex:1018-1026`) with per-family scope at `rem:R-descent-unitarity-scope`. OF4 HEALED — theorem retitled "Six-fold TFAE on the Koszul locus, with a class-G seventh." OF5 HEALED — PTVV-for-A confabulation retracted via scope remark at `theorem_A_infinity_2.tex:945-970`. **Bibkey phantoms HEALED** (AP281 acute): `HackneyRobertson2017`, `HackneyRobertson2019`, `Francis2012`, `GR17`, `Positselski2011`, `Positselski2018`, `Hinich2003` were load-bearing in the proof chain but absent from `standalone/references.bib` (rendering `[?]` at build); alias entries inscribed at `references.bib:1100-1194`. Remaining open: (OF1) modular-family extension to M̄_{g,n} including boundary at chain level; non-finitely-generated MC4-completion regime (Conjecture Π 4). |
| B | **PROVED at coderived level weight-completed only** via `thm:chiral-positselski-weight-completed` (Vol I `chapters/theory/theorem_B_scope_platonic.tex:245-318`), with per-weight step `thm:chiral-positselski-at-each-weight` (:177-229). Raw class-M direct-sum level: GENUINELY FALSE, proved via `prop:chiral-positselski-raw-direct-sum-class-M-false` (:340-387) with explicit L_0 witness. Chain-level class G (Heisenberg: chain-level qi via on-the-nose harmonic vanishing; direct evaluation of `thm:bv-bar-coderived` at Vol II `bv_brst.tex:2088-2094`) and class L (affine KM non-critical: chain-level qi via PBW-associated graded reduction to genus-0 weight-by-weight comparison; same mechanism as class G): `\ClaimStatusProvedElsewhere` — attributed to Vol II `thm:bv-bar-coderived` (`chapters/connections/bv_brst.tex:2031`); NOT inscribed in Vol I. Chiral adaptation of Positselski 2011 Thm 4.6 + Thm 7.2 is CITED at citation level in Steps 2 and 4 of the weight-completed proof, not independently inscribed. Alias notice (AP286 Option-A): status-table aliases `thm:chiral-positselski-7-2` and `thm:chiral-positselski-5-3` are legitimate `\phantomsection` aliases co-located with the canonical theorem at `chapters/theory/theorem_B_scope_platonic.tex:246-247`, added in Wave-6 so that `\ref{}` consumers resolve to the canonical Positselski-chiral theorem without duplicating the body. No inscription action required; retain for consumer stability (AP285 alias-section-number discipline). | Weight-completed coderived statement: for conilpotent chiral CDG-coalgebra C on X with finite-dim graded pieces, counit Ω_X(B_X(C)) → C is a coderived equivalence in the weight-completed ambient. Proof structure: chiral Φ/Ψ adjoints + bicomplex totalisation + conilpotent contracting homotopy + coacyclic cone, with Positselski's associative-case Thm 7.2 / co-contra Thm 5.3 invoked as input. Raw class-M direct-sum falsity witness: S_4(Vir_c) = 10/[c(5c+22)] ≠ 0 forces bar cohomology at weight 4 in `Ch(Vect)` (bounded direct-sum ambient). **FTM seven-fold TFAE** resolves g=0 tautology. Ordered E_1-variant proved for Yangian (EK quantisation acyclicity = FM^ord Koszulness). |
| C | C0 and C1 H-level eigenspace proved unconditionally on Koszul locus; C0 ordinary-derived realization conditional on perfectness, unconditional for Heis + affine KM at positive integer or boundary admissible level; generic non-critical KM k is a separate conjecture (`conj:perfectness-boundary-km-generic`); C2 shifted-symplectic upgrade conditional on BV package (the BV package itself is proved at coderived level via Vol II `thm:bv-bar-coderived`, chain-level strict qi for classes G/L, coacyclic for class M) | Honest accounting after 2026-04-17 audit + 2026-04-18 Wave-14 heal (Vol I `chapters/theory/theorem_C_refinements_platonic.tex`). INSCRIBED: (T5) C2 hypotheses pinned to BV + Verdier + bar-chart lift; (T7) MC5 strengthening closes class-M C2(iii) weight-completed; (T8) explicit δF_2^cross(W_3) = (c+204)/(16c) in `standalone/multi_weight_cross_channel.tex`; (T3) C1 reflexivity unconditional on the expanded Heis + KM-integer/admissible locus via `cor:c0-unconditional-class-G-L`. CONDITIONAL: (T1) `conj:derived-center-koszul-equivalence`: brace dg Z^der_ch(A) ≅ Z^der_ch(A^!) as E_2-algebras after Deligne-Tamarkin CONJECTURE; only H^0 naive-center `lem:naive-center-koszul-identification` unconditional. Hinich 2001 rectification and Francis–Gaitsgory 2017 were attacked as candidate bypasses and neither eliminates the DT dependence (`rem:derived-center-sharpened-scope` inscribes the attack verdict + explicit falsification test on weight-2 ChirHoch^2 pairings for Heisenberg and Virasoro). (T2) `prop:perfectness-standard-landscape` scope REFINED 2026-04-18: unconditional for Heis + affine KM at positive integer level (TUY sewing) + affine KM at boundary admissible level (Arakawa C_2-cofin); generic non-critical KM recorded as `conj:perfectness-boundary-km-generic`; class-M boundary-stratum recorded as `conj:perfectness-boundary-class-M`. The Wave-14 heal separates the TUY / Arakawa / Kac–Kazhdan inputs which were previously conflated into a single "non-critical level" catch-all. (T9) `thm:C-PTVV-alternative` mixed: (i) and (ii) PTVV + HAG-II genuinely independent (tangent-complex formula now explicit: `T_{[*]} R Map(M̄_g, BG_A) = R Γ(M̄_g, Z(A))[1]`); Lagrangian clause (iii) forwards to `thm:quantum-complementarity-main` (not circular: target is proved from genus filtration + Verdier duality + involution-splitting, no back-reference to PTVV-alternative). RETRACTED: (T4) `thm:theorem-C-g0` was non-issue: +3 shift is misreading of -(3g-3) at g=0 (see `rem:plus-three-shift-nonissue`). (T6) folded into T3. Remaining frontier: chain-level brace equivalence; class-M boundary perfectness; generic non-critical KM nodal propagation; intrinsic derived-stack involution; multi-weight ω_g^cross; Faber λ_g-residue. Audit trail: `attack_heal_theorem_C_20260418.md`. |
| D | PROVED unconditional at $g=1$ (all uniform-weight families via 1-dim $H^2(\overline{\cM}_{1,1}) = \Q \cdot \lambda_1$); PROVED unconditional at $g=2$ for scalar-diagonal families (Heisenberg, Virasoro, rank-1 free fields unconditionally; $V_k(\fg)$ at positive-integer levels via Frenkel--Kac lattice realisation; $V_k(\fg)$ with $\fg$ non-abelian at generic level conditional on scalar-diagonality of the Wakimoto ghost sector); at $g \geq 3$ conditional on (a) scalar-diagonal hypothesis for multi-generator uniform-weight and (b) Faber's $\lambda_g$-conjecture for on-the-nose lift from the socle; PTVV alternative conditional on Mok25 log-FM chain-level sewing at $g \geq 3$ only ($g = 2$ uses classical Knudsen--Mumford stable reduction; $g = 1$ needs neither) | obs_g = kappa * lambda_g uniform-weight. **Wave-2 heal 2026-04-18**: inscribed `rem:kac-moody-obs-scope` at `higher_genus_foundations.tex:5072` recording that Step~2--3 of `thm:kac-moody-obs` tacitly uses scalar-diagonality at $g \geq 2$: unconditional at positive-integer level via Frenkel--Kac lattice realisation; open at generic level. See `attack_heal_theorem_D_20260418.md`. **2026-04-17 Beilinson audit Wave 1 (F1-F7) inscribed**: honest scope rewrite of `prop:scalar-obstruction-hodge-euler`, `prop:clutching-uniqueness`, `prop:theorem-D-factorization-homology-alt`. (F1+F2) **Scalar-diagonal hypothesis** inscribed as `def:scalar-diagonal-hypothesis`: the K-theoretic construction of Step 1 DEFINES obs_g by projection to the weight-$p$ scalar channel $\varepsilon_p$; coincidence with the physical bar-curvature obstruction requires scalar-diagonality (automatic for rank-1 strong-generator families; genuine hypothesis for multi-generator, e.g. affine KM $V_k(\fg)$ with $\fg$ non-abelian). `rem:scalar-diagonal-honest` states the honest reading. (F3) Step 3d rewritten as "differential-geometric corroboration, not independent derivation": prior splitting-principle "scalar-channel linearity" step conflated $\kappa^g \prod x_i$ with a polynomial in $\kappa$ possessing a linear-in-$\kappa$ component to project onto (a polynomial identity mistaken for a cohomological one); load-bearing derivation is Step 1c's K-theoretic linearity, and Step 3 (BGS/Arakelov/Quillen) supplies the Chern--Weil representative. (F5) **AP225 SCOPE CORRECTED**: `prop:clutching-uniqueness` retitled "socle scope"; $g \in \{1,2\}$ on the nose ($\cN^g = 0$); $g \geq 3$ on the socle $R^g/\cN^g$ only, on-the-nose lift conditional on the $\lambda_g$-conjecture. `rem:ap225-scope-correction` inscribed. (F6) **H04 PTVV downgraded from "genuinely disjoint" to "complementary presentation with distinct scope", and status-tag downgraded from `\ClaimStatusProvedHere` to `\ClaimStatusConditional`** (AP149 propagation, source-of-truth `higher_genus_foundations.tex:6117-6118`): `prop:theorem-D-factorization-homology-alt` part (iii) proof rewritten to NOT cite `prop:scalar-obstruction-hodge-euler`, and its environment header now carries `\ClaimStatusConditional` per the Mok25 log-FM chain-level inscription gap (the PTVV transgression at the nodal boundary of $\overline{\cM}_g$ is cited at citation level, not inscribed). Scalar coefficient $\kappa(\cA)$ detected on the genus-1 PTVV pairing directly; boundary restrictions handled via clutching-uniqueness on the socle. The two proof chains share no lemmas but cover complementary scopes: PTVV+clutching gives all uniform-weight families on the socle without scalar-diagonality (now CONDITIONAL on Mok25 chain-level log-FM); Arakelov--BGS gives on-the-nose for scalar-diagonal families without the $\lambda_g$-conjecture. HZ-IV decorators built on PTVV must declare `verified_against = clutching_socle_PTVV` with an honest `disjoint_rationale` recording the scope constraints AND the Conditional status of the PTVV route. (F7) **$\delta F_g^{\mathrm{cross}}$ scope**: explicit closed form $(c+204)/(16c)$ proved at $g=2$ for $W_3$ only; $g \geq 3$ multi-weight correction named via graph-sum construction of `thm:multi-weight-genus-expansion`, closed form CONJECTURAL. Audit trail: `adversarial_swarm_20260417/wave1_theorem_D_attack_heal.md`. |
| H | PROVED sharp Hilbert series on Koszul locus; E_1-variant INSCRIBED | ChirHoch*(A) concentrated in degrees {0,1,2} with exact formula `P(t,q) = HS_{Z(A)}(q) + HS_{ChirHoch^1(A)}(q)·t + HS_{Z(A^!)}(q)·t^2`. Critical level k=-h^v excluded (load-bearing; `rem:critical-level-dimensional-divergence` at `chiral_hochschild_koszul.tex:1850-1897`: Feigin-Frenkel centre infinite-dimensionalises ChirHoch^0, FM-tower collapse fails, Koszul hypothesis fails). **Step 3 circularity RESOLVED (2026-04-16)**: rerouted through `thm:pbw-koszulness-criterion` (PBW collapse, inscribed at `chiral_koszul_pairs.tex:784`) via Shelton-Yuzvinsky Koszulity of the pure-braid Orlik-Solomon algebra + chiral quadratic-Koszul lemma. **AP242 spot-check (2026-04-18)**: `lem:chiral-quadratic-koszul` IS inscribed with a proof body in `chiral_hochschild_koszul.tex:656-727` (NOT a forward reference). Proof rewritten 2026-04-18 to make the Fresse-transport non-degeneracy step explicit (Step 1 = fibrewise-flat reduction to $\cF(E_p)(3)$ where the classical OPE-residue pairing is non-degenerate via Fresse Prop. 7.1.3, with propagation via $\cD_{X^3}$-coherence + semicontinuity + Kontsevich formality at collision divisors; Step 2 = Fresse transport along the pairing; Step 3 = quasi-isomorphism on Koszul locus). Three HZ-IV decorators LIVE at `compute/tests/test_theorem_H_hochschild_koszul.py` (Heisenberg → dim triple (1,1,1) total 3; Virasoro → (1,0,1) total 2; affine sl_2 → (1,3,1) total 5) wired to `compute/lib/chirhoch_dimension_engine.py`. Honest HZ-IV residual: the engine's dataclass values cite `chiral_center_theorem.tex` as their VERIFIED source, which is the derivation chapter under review; genuine disjoint wiring to Feigin-Fuchs 1984 Fock-BRST engine and Wang 1998 semi-infinite BRST engine is not yet implemented (decorator prose is correct, test body disjointness is partial — AP277 residual, not vacuous). **E_1-variant `thm:hochschild-concentration-E1`** INSCRIBED at `chiral_hochschild_koszul.tex:1339-1390`: proof via ordered FM + pure-braid Orlik-Solomon Koszulity (Shelton-Yuzvinsky supersolvability of the pure braid arrangement) + ordered Koszul dual; no appeal to `thm:bar-concentration`, no $\Sigma_n$-averaging. Hypothesis is the PBW criterion (an E_∞-chiral condition); the subscript "E_1" denotes the first page of the collision-depth SS, NOT an E_1-chiral input (disambiguated at `rem:E1-scope-hochschild-concentration`:1392-1417). **AP255 phantom-file flag in `standalone/theorem_index.tex` healed 2026-04-18**: rows for `thm:pbw-koszulness-criterion` and `thm:hochschild-concentration-E1` were tagged `chiral_pbw_koszulness.tex [PHANTOM]`; both labels are in fact inscribed in live chapters. Index corrected. ALT H05 retrievable as sharp-Hilbert-series from `fm-tower-collapse` E_2 page. |
| MC1-4 | **MC1/MC2 UNCONDITIONAL on Koszul locus; MC3 PROVED on evaluation-generated core via five-family mechanism (Baxter RETRACTED); uniform Chari-Moura multiplicity-free ℓ-weights holds for classical + simply-laced exceptional types; non-simply-laced exceptional (G_2, F_4) per-case via Hernandez 2006/2010; MC4⁺ UNCONDITIONAL on classes G/L/C via `thm:completed-bar-cobar-strong`; MC4⁰ CONJECTURAL at generic parameters for class M principal** (2026-04-18 retraction per AP269 fabrication audit; see `adversarial_swarm_20260418/attack_heal_mc4_0_20260418.md`) | MC1 universal (semisimple hyp) via Whitehead; **MC1 Virasoro g≥2 PROVED unconditionally** via L_0-diagonalization on augmentation ideal (bypasses hypothesis (c) semisimple — family-specific L_0 argument suffices). **Feigin-Fuks resolution provides independent verification** + extends to principal W_N generic Ψ via Feigin-Frenkel screening. MC2 bar-intrinsic PROVED. **MC3 PROVED via five-family mechanism on the EVALUATION-GENERATED CORE** (AP47 scope; Vol I `chapters/theory/mc3_five_family_platonic.tex` 787 lines + 8 HZ-IV decorators, 2026-04-16; Wave-15 scope refinement 2026-04-18): asymptotic prefundamentals (type A); reflection-equation Shapovalov (B/C/D); Chari-Moura multiplicity-free ℓ-weights (uniform, classical + simply-laced exceptional; G_2/F_4 per-case via Hernandez 2006/2010); theta-divisor complement (elliptic); parity-balance (super-Yangian, **gl-super sub-case only**). Baxter constraint retracted as type-A rational artifact. Extension beyond the evaluation-generated core to the full compact/completed DK_g factorization category is CONDITIONAL in every simple type (2026-04-18 Wave MC3-compact-completion heal per AP1181 circular-gloss retraction): type A carries a single remaining packet `conj:dk-compacts-completion` (two clauses per `rem:dk-compacts-completion-two-clause` at `yangians_computations.tex:3947`: (C1) compact-core extension of Φ inside the separated completion with the three weight-bounded/colimit/compactness hypotheses of `prop:completion-extension-weight-bounded` verified in the chosen ambient; (C2) Francis-Gaitsgory pro-nilpotent comparison identifying that ambient with the MC3 target); outside type A the gap splits into three mechanically independent packets `conj:rank-independence-step2` (Lemma L, lift-and-lower, for shifted-envelope generation), `conj:mc3-automatic-generalization` (type-independence of the remaining pro-Weyl and completion mechanisms), and `conj:dk-compacts-completion`. See `adversarial_swarm_20260418/attack_heal_mc3_compact_completion_20260418.md`. **Wave-15 scope refinement 2026-04-18 per AP361**: super-Yangian family splits into Y(gl(m|n)) parity-balance sub-case and Y(osp(m|2n)) super-reflection-equation sub-case, the latter covered under the twisted B/C/D reflection-equation mechanism in its Z/2-graded form via ACDFR 2003 super-Sklyanin-Berezinian determinant sBer K(u); K3 super-Yangian Y_osp(4|20) lies in the osp-super sub-case. Inscribed at `rem:super-yangian-scope` + `rem:parity-balance-scope-gl-only` in `mc3_five_family_platonic.tex`. Asymptotic-vs-multiplicity-free scope split healed at `rem:chari-moura-vs-hernandez-jimbo` (AP362): two theorems at two category levels (pro-completed Ô vs uncompleted O), not one theorem in two proofs. Silent non-coverage enumerated in `thm:mc3-full-DK-conjectural` clause (iv) as five sub-cases (iv.a)-(iv.e) (2026-04-18 Wave-N2-SCA heal): (iv.a) logarithmic W(p), (iv.b) N=2 SCA via candidate Kazama-Suzuki coset `conj:mc3-n2-sca-kazama-suzuki` in `chapters/examples/n2_superconformal.tex`, (iv.c) cosets beyond GKO/KRW, (iv.d) indefinite non-rational lattice VOAs, (iv.e) roots-of-unity admissible-level. Koszulness of N=2 SCA at generic k is separately PROVED via `prop:n2-koszulness` (PBW collapse + Adamović 1999 free strong generation); MC3 for N=2 SCA remains open (AP501: Koszulness without MC3 is structurally independent; AP502: coset MC3 non-inheritance). Six confusion patterns #240-245 cached; AP361-AP365 registered Wave-15; AP501-AP502 registered Wave-N2-SCA; see `adversarial_swarm_20260418/attack_heal_n2_sca_20260418.md`. MC4⁺ UNCONDITIONAL (classes G/L/C bounded-degree conformal/fermion-number filtrations; inscribed `thm:completed-bar-cobar-strong` at `chapters/theory/bar_cobar_adjunction_curved.tex:968-1050`). **MC4⁰ CONJECTURAL (2026-04-18 Wave-MC4-swarm, AP421)** at generic parameters for class M principal: the prior UNCONDITIONAL advertisement was downgraded per AP269 fabrication audit — zero `h_htpy` / Wakimoto / SDR grep hits in any Vol I chapter, the standalone `N4_mc4_completion.tex` advertised inscription in `chapters/theory/bar_cobar_adjunction_curved.tex` that does not exist, and the W_N Feigin-Frenkel route is actively contradicted by `prop:ff-screening-coproduct-obstruction` at `ordered_associative_chiral_kd.tex:10176-10297`. Standalone healed: `thm:n4-mc4-zero-unconditional` → `conj:n4-mc4-zero-generic-parameters` with AP266 sharpened-obstruction `prop:mc4-zero-wakimoto-sdr-obstruction` recording the (Ψ-1)/Ψ cocycle class in H^1_{ch} matching `thm:miura-cross-universality`; six-gap status remark `rem:mc4-zero-sdr-candidate-gaps`. Hook-type `thm:n4-non-principal-hook` → `conj:n4-non-principal-hook` (matches CLAUDE.md Non-principal W row's self-declared gap). Honest frontier: explicit chain-level SDR for Virasoro at generic c (Wakimoto Fock + Feigin-Fuks irreducibility + stage-to-stage lift), resolution of the FF chiral-coproduct obstruction for W_N, theorem-level reduction of subregular/minimal to class-C βγ. MC2 ALT H06 (KS scattering) conditional on support-property comparison (open). Audit trail: `adversarial_swarm_20260418/attack_heal_mc4_0_20260418.md`. |
| MC5 | **STABLE ONLY AFTER AMBIENT SEPARATION**: raw direct-sum class-M chain-level statements in `Ch(Vect)` remain false; strict completed/pro and coderived surfaces are the verified or conditional repair surfaces. | Do not describe the direct-sum failure as an ambient artefact and do not identify the inverse-limit product with the original direct-sum complex. Finite S4 is only the first detected term; the real obstruction is the nonterminating harmonic family. Continuous completed/coderived equivalences require the explicit CP1--CP3 co/contra comparison package. |
| Koszul | **8 genuinely independent bidirectional equivalences** unconditional on generic parameters (non-critical, non-admissible, non-special); +1 one-way ChirHoch consequence; +1 perfectness-conditional Lagrangian; +1 one-directional D-mod purity (converse proved for affine KM only; Virasoro/W open); +1 uniform-weight-conditional genus refinement of (vii) (all-genera under uniform-weight; multi-weight genus-0 only). Critical/admissible/minimal loci require case-by-case verification; chapter opening at `chiral_koszul_pairs.tex:77-89` records the honest "9+1+1+1" count, theorem header at :2348 is the legacy "10 unconditional" and needs re-alignment. | Independent bidirectionals with (i): (ii), (iii), (iv), (v), (vii)|_{g=0}, (ix), (x), and (i) itself — eight in total. Equivalence (vi) Barr-Beck-Lurie is an ∞-categorical restatement of (v) (scope remark at :2335-2338 acknowledges it is a consequence, not an independent route). Equivalence (viii) ChirHoch concentration/duality is ONE-DIRECTIONAL only (explicit at :2503): concentration in degrees {0,1,2} follows from Koszulness; the reverse implication (concentration forces Koszulness) is not proved. Equivalence (xi) Lagrangian is conditional on perfectness of the fiber-center complex + non-degeneracy of the ambient tangent complex (not merely "BV package"). Equivalence (xii) D-module purity is ONE-DIRECTIONAL (xii)⇒(x); the converse (x)⇒(xii) is proved for affine Kac-Moody at non-critical level only — Virasoro and principal W_N converse are open. Equivalence (vii) uniform-weight genus refinement: all-genera under uniform-weight hypothesis; multi-weight collapses to genus-0. All-genera (vii) fails at critical level k=-h^v where κ=0 and bar cohomology becomes the opers complex. "Tropical Koszulness" is advertised in the preface as a 13th equivalence but is NOT inscribed in the chapter meta-theorem; treat as heuristic pending inscription. Cross-file count propagation: this correction affects `landscape_census.tex`, `working_notes.tex`, four surveys, and Vol II `thqg_introduction.tex` (nine files total; propagate in a separate pass). ALT: proof web (H09). |
| D^2=0 | Convolution level UNCONDITIONAL; ambient level CONDITIONAL on Mok25 | Two inscribed theorems with distinct dependency profiles. `thm:convolution-d-squared-zero` (higher_genus_modular_koszul.tex:32404, ClaimStatusProvedHere) uses only topological $\partial^2 = 0$ on $C_*(\overline{\cM}_{g,n})$; no Mok25, no log-FM, no relative FM; load-bearing for Theorems A-D-H, MC2, recursive-existence shadow tower; unaffected by any Mok25 retraction. `thm:ambient-d-squared-zero` (:32418) is the ambient-level enrichment with five-component differential and fiberwise curvature $\dfib^2 = \kappa(\cA)\cdot\omega_g$; genuinely depends on Mok25 log-FM normal-crossings + planted-forest stratification (cites Mok25 Thm 3.3.1, Thm 5.3.4). Wave-X Mok25 programme-wide audit (`adversarial_swarm_20260418/attack_heal_mok25_programme_audit_20260418.md`, 2026-04-18) SUPERSEDES the prior iteration-2 Thm-A-OF1 finding that labelled the `Mok25` bibkey "AP281 phantom / no locatable primary source / fabricated-mechanism territory": Mok25 IS a real preprint — Siao Chi Mok (Cambridge/Imperial, algebraic geometry), arXiv:2503.17563, "Logarithmic Fulton--MacPherson configuration spaces" (submitted 21 March 2025, v2 21 May 2025). The iteration-2 audit conflated Siao Chi Mok with Chung-Pang Mok (Purdue, automorphic forms) due to drifted `author`/`title` fields in `standalone/references.bib:601-607` (now harmonized). Ambient-level status-tag `ClaimStatusProvedHere` remains CONDITIONAL on Mok25 cited at citation level (a single-preprint pillar, AP1241) but is not "fabricated-mechanism" conditional; the honest scope remark `rem:mok-dependency` (:32568-32606) is the correct architectural accounting. Guide-to-main-results table at `guide_to_main_results.tex:117-121` inherits the same drift and needs the same retag. AP1341 registered: status-row headline conflates two theorems with distinct dependency status. Audit trail: `adversarial_swarm_20260418/attack_heal_d2_zero_convolution_20260418.md`. |
| Theta_A | PROVED | Bar-intrinsic; all-degree inverse limit (thm:recursive-existence). |
| SC-formal | PROVED forward unconditional (Heisenberg / lattice VOA / free fermion via Gaussian two-point kernel); converse CONDITIONAL on the two-colour transfer lemma ``open-colour SC output vanishes $\Longrightarrow$ underlying $\Eone$ tree vanishes'' (rem:sc-formal-open-transfer-frontier). | SC-formal iff class G. Forward: shadow tower truncation at degree~$2$ on the Gaussian locus kills every connected genus-$0$ tree with $r \geq 3$ external legs. Converse: shadow tower controls SC ops modulo the open-to-$\Eone$ transfer lemma. Class-$M$ non-vanishing: nonzero for infinitely many $k \geq 5$ (Riccati algebraicity via Theorem~\ref{thm:shadow-archetype-classification}), NOT degree-by-degree. The ``$\Delta = 0$'' shorthand abbreviates ``shadow tower truncates at degree~$2$''; $\Delta = 0$ alone is necessary but not sufficient (class $L$ has $\Delta = 0$ but $S_3 \neq 0$). Operadic rewrite is the PRIMARY proof; no separate ALT survives as an independent second proof. |
| Depth gap | PROVED on the chirally Koszul standard landscape | d_alg in {0,1,2,inf}; gap at 3. Betagamma d_alg=2 witness on standard conformal-weight line. Impossibility of 3 via MC relation + shadow Lie Jacobi. Scope: result is proved on the standard-landscape locus where chiral Koszulness holds (off-locus chain-level d_alg status is open per 2026-04-13 Platonic roadmap). ALT: representation-theoretic (H10). |
| ChirHoch^1 KM | PROVED at generic k for general simple g: ChirHoch^1(V_k(g)) = g and dim ChirHoch^*(V_k(g)) = dim(g) + 2 | prop:chirhoch1-affine-km (chapters/theory/chiral_center_theorem.tex:2133) + prop:chirhoch2-affine-km-general (:2223, ClaimStatusProvedHere, inscribed Wave-4 2026-04-18) give ChirHoch^0 = C, ChirHoch^1 = g, ChirHoch^2 = C universally at non-critical level via Koszul-dual center identification (dim ChirHoch^2(A) = dim Z(A^!) inscribed in thm:main-koszul-hoch at chiral_hochschild_koszul.tex:1791); Feigin-Frenkel dual level -k - 2h^v is also non-critical at generic k, so Z(V_{-k-2h^v}(g)) = C via Kac-Shapovalov nonvanishing (FF92, Kac98 Thm 4.10). Theorem H polynomial concentration (thm:hochschild-polynomial-growth) kills ChirHoch^{>=3}. Rank-independent slots: ChirHoch^0 (vacuum line) and ChirHoch^2 (Koszul-dual center collapse) both give dim 1; the rank-scaling slot is ChirHoch^1 = g. Predictions: sl_2 total 5, sl_3 total 10, sl_N total N^2+1, E_8 total 250. Candidate rank(g)+dim(g)+1 formula (giving sl_3 = 11) REFUTED: it would require ChirHoch^0 to scale with rank, but Z(V_k(g)) = C for all simple g at non-critical k. HZ-IV triangulation via engine extension to sl_3, sl_4 at small integer k is recommended follow-up; the proof stands on general Koszul-duality machinery and Kac-Shapovalov nonvanishing at dual level. |
| Topologization | **PROVED on original complex for G/critical; PROVED Q_{tot}-cohomological for L; class-M original-complex chain-level remains open outside the stated completed/coderived surfaces** (thm:topologization-tower; post-Wave-1 Beilinson audit 2026-04-17 identifies three open directions, not one) | **Topologization tower theorem (2026-04-16; status-table corrected 2026-04-17; MC5 scope corrected 2026-04-24)**: (1) class G (Heisenberg): E_3^top on original complex unconditional (commutative ⇒ Dunn on F^{H_k}); (2) class L (V_k(g), k≠-h^v): E_3^top on Q_{tot}-cohomology via `thm:E3-topological-km` and `thm:iterated-sugawara-construction` (Sugawara identity T_Sug=[Q_{tot},G_Sug] proved in H^•(A^{BV}_{3d}, Q_{tot}); the strict chain-level upgrade via explicit η_1 antighost-contact is a FRONTIER item; draft candidate formulas inscribed in Vol II `e_infinity_topologization.tex:401-411` were retracted 2026-04-17 pending independent verification); (3) class C (βγ, bc): E_3^top via FMS bosonization reducing to Heisenberg; (4) class M (Vir, W_N): completed/coderived E_3^top statements require the same strict tower and CP1--CP3-style continuous comparison hypotheses as MC5; the raw direct-sum ambient in `Ch(Vect)` is a genuine failure surface, not an artefact; (5) critical level k=-h^v: E_2^top not E_3 (Sugawara collapse is a dimension drop, not a failure). **Honest frontier inventory (post-Wave-1, corrected 2026-04-24)**: (I) class M on the original direct-sum complex in `Ch(Vect)` remains open/negative depending on the exact comparison map; completed/pro replacements do not identify with the direct-sum object without a separate theorem; (II) class L strict chain-level upgrade of the Q_{tot}-cohomological identification to the raw bar complex via explicit η_1 antighost-contact (required computation: [Q,tilde G_1] − T_Sug identically zero in the bulk BV chain complex, not merely Q-exact; Vol II candidate formulas retracted pending verification); (III) antighost BRST-commutativity [G^{(n)},G^{(m)}]=Q-exact at spin n,m≥3 (postulated as axiom (T5) of def:higher-spin-stress-tower, not yet derived; `thm:iterated-sugawara-construction` for n≤3 survives, the W_∞ convergence theorem inherits the conjectural status). ALT via CG factorization homology H08 no longer conjectural for G/C and remains conditional for M-completed unless the completion hypotheses are verified. |
| E_3 identification | **PROVED simple g, affine KM V_k(g) at generic non-critical level (`thm:e3-identification`, `en_koszul_duality.tex:5338`); gl_N + reductive semisimple + abelian/Heisenberg ENDPOINTS addressed at REMARK level only (`rem:e3-non-simple`:5541, `rem:e3-non-simple-gl-N`:5576, `rem:e3-heisenberg-endpoint`:5711); nilpotent/solvable UNINSCRIBED** | Z^der_ch(V_k(g)) ≅ A^λ as E_3-deformation families over λ H^3(g)[[λ]], λ = k+h^v. **E_3-top vs E_3-chiral scope (AP-TOPOLOGIZATION)**: thm:e3-identification(iv) upgrades to E_3-TOPOLOGICAL via Sugawara only at non-critical level; without conformal vector stuck at E_3-chiral (Dunn-colored E_2-chiral ⊗ E_1-top). Critical level k=-h^v: Sugawara undefined, perturbative only. **Wave-4 2026-04-18 phantom audit + heal (AP255/AP285/AP264/AP286)**: (i) six formerly-phantom labels confirmed absent across all three volumes (`prop:e3-gl-N`, `thm:e3-identification-semisimple`, `thm:e3-identification-reductive`, `prop:e3-heisenberg`, `thm:e3-identification-solvable`, "explicit E_3-algebraic vs E_3-topological comparison diagram"); (ii) one phantom CONSUMER retargeted: `compact_completed_mc3_comparison_platonic.tex:322` `\ref{thm:e3-identification-reductive}` → `\ref{thm:e3-identification}` plus `\ref{rem:e3-non-simple}` and `\ref{rem:e3-non-simple-gl-N}` (Option (b) retargeting rather than inscribing a phantom theorem); (iii) Fresse citation drift HEALED at `en_koszul_duality.tex:5278` and `:5430` — was Theorem 16.1.1, corrected to 16.2.1 matching the load-bearing citation in `e3_identification_chain_level_platonic.tex:467,:816,:831` (AP285 alias section-number drift); (iv) Heisenberg endpoint inscribed as `rem:e3-heisenberg-endpoint`:5711 (abelian closed form: H^3=0, Sym^2 of rank (r+1 choose 2), topological enhancement trivial since Heisenberg stress is inner at every level). Honest scope per Option (b)+(c) of AP255 Wave-5 heal menu: one theorem (simple g, generic non-critical), three remarks (reductive semisimple, gl_N, abelian Heisenberg), zero phantoms. The CY_4 p_1-twisted row stands separately: CONSTRUCTED+CONJECTURED, not PROVED. Inscription of solvable/nilpotent cases as theorems remains an open target. Attack-and-heal note: `adversarial_swarm_20260418/wave4_e3_identification_phantom_audit.md`. |
| Elliptic chiral QG | BELAVIN r-matrix PROVED for sl_2 at chapter level; Felder↔KZB leading-order ℏ CONJECTURAL (orphan standalone only) | Belavin non-dynamical $r^{\mathrm{ell}}_{\fsl_2}(z,\tau)$ is PROVED as chapter-level `\ClaimStatusProvedHere` in Vol~I `prop:elliptic-rmatrix-shadow` (`yangians_drinfeld_kohno.tex:7251-7389`) and `thm:g1sf-elliptic-rmatrix` (`chapters/connections/genus1_seven_faces.tex:470-546`). **Propagator convention (AP275-corrected 2026-04-18)**: the Cartan-channel propagator is the FULL Weierstrass zeta $\zeta_\tau(z) = \theta_1'(z\|\tau)/\theta_1(z\|\tau) + 2\eta_1 z$; the quasi-period correction $2\eta_1 z$ is LOAD-BEARING — dropping it and writing only the theta derivative $\theta_1'/\theta_1$ breaks classical Yang--Baxter off-diagonal (see FM30 and `rem:elliptic-rmatrix-scope` at `yangians_drinfeld_kohno.tex:7364-7389`). Root channels use theta-function ratios $\phi_\pm(z,\tau)$. Three-sample CYBE verification in `elliptic_rmatrix_shadow.py` (`verify_ybe_elliptic_sl2`, residual $<10^{-6}$). **Felder dynamical $R$-matrix ↔ KZB equivalence (distinct IRF/face-model statement, dynamical YBE not ordinary CYBE): NOT inscribed at chapter level**; only in orphan `standalone/seven_faces.tex` as `conj:chiral-qg-equiv-elliptic-sf` (downgraded 2026-04-18 from `thm`+ProvedElsewhere to `conj`+Conjectured with convention fix per AP40+AP125). The orphan is not `\input`-ed into `main.tex`, zero `\ref` consumers — AP255 phantom-in-chapters but build-neutral. Heat-equation prefactor $1/(4\pi i)$ diag vs $1/(2\pi i)$ off-diag inscribed at `rem:kzb-heat-prefactor` (`ordered_associative_chiral_kd.tex:12260`) for the 2-point commutator only; full Bernard heat identity NOW INSCRIBED at `rem:bernard-heat-identity-zeta` (`:12276`) and VERIFIED engine-side 2026-04-18 (see KZB flatness row for AP521-AP524 attack-heal; prior confabulated `d_τ(℘_1)` formula retracted per AP282/AP283). |
| CY_4 p_1-twisted | CONSTRUCTED + CONJECTURED (Vol III `cy_to_chiral.tex:4502-4700`, construction of p_1-twisted family + `conj:cy4-p1-family-associator-sextic` verified only at the sextic); PROVED NEGATIVE: "no native CY_4 Yangian exists" (Vol III `en_factorization.tex:289-299`) since π_4(BU)=Z is a Pontryagin BLOCK-obstruction, not a positive structure | Cocycle candidate c(x,y) = ⟨x ∪ y ∪ p_1(T_X), [X]⟩/24. K3×K3 N(X)=0 claim UNINSCRIBED (generically false since p_1(K3) ≠ 0). 7d hCS realization CITED (Costello-Gaiotto), not inscribed. |
| Toroidal chiral QG | PROVED at formal disk (Vol~III chapter-level + Vol~I standalone); global P¹×P¹ conditional on class-M chain-level on the ORIGINAL complex; AP255 standalone-orphan noted | **Canonical chapter-level home**: Vol~III `chapters/theory/quantum_groups_foundations.tex:199-231` — Ding--Iohara--Miki toroidal algebra $U_{q,t}(\ddot{\fsl}_N)$ with two spectral parameters $(u,v)$ is the chiral bialgebra of an ordered-Koszul factorization algebra on the formal disk $(\mathrm{Spec}\,\bC[\![z]\!])^{\times 2}$, via FHHSY shuffle $\mathrm{Sh}_{q,t}(\fsl_N)$ and Schiffmann--Vasserot CoHA on instanton moduli. Miki $S_3$ triality of $\Omega$-background $(q_1,q_2,q_3)$ under $q_1 q_2 q_3 = 1$. Vol~I standalone inscription: `thm:chiral-qg-equiv-toroidal-sf` at `standalone/seven_faces.tex:1043` (AP257 line-drift corrected 2026-04-18 from previously advertised :1020; body also renamed Drinfeld--Miki→Ding--Iohara--Miki to match Vol~III convention). **Bialgebra scope (AP263)**: antipode non-lifting is a proved negative for toroidal (same as for Yangian); "quantum group" here is shorthand for chiral bialgebra. **Prior companion label `thm:chiral-qg-equiv-toroidal-formal-disk` was PHANTOM (AP255)**: zero `\label{}` inscriptions anywhere in Vol~I; the only use site was a `\verb|...|` prose reference at `standalone/N3_e1_primacy.tex:1160`, retargeted 2026-04-18 to the real `-sf` label with explicit note that the earlier suffix was never inscribed. **Orphan standalone (AP255)**: `seven_faces.tex` is not `\input`-ed into `main.tex`; zero `\ref{thm:chiral-qg-equiv-toroidal-sf}` consumers exist programme-wide, so no `[?]` at build. **Global P¹×P¹ obstruction**: conditional on the class-M chain-level topologization frontier on the ORIGINAL bar complex in $\Ch(\Vect)$ (direction~I of the topologization-tower open inventory; the pro-object, $J$-adic, and weight-completed ambients already yield chain-level class-M topologization, so closing direction~I closes the toroidal global case as a downstream corollary). 2-parameter RTT via two invariant bilinear forms on $\gl_N$ ($B_{\tr}$, $B_{\ab}$, matching `rem:e3-non-simple-gl-N`). See `adversarial_swarm_20260418/wave4_toroidal_chiral_qg_attack_heal.md`. |
| Chiral QG equiv | **ORDERED (base): PROVED on the Koszul locus via `thm:chiral-qg-equiv` at `chapters/theory/ordered_associative_chiral_kd.tex:8404` (`\ClaimStatusProvedHere`, title "Chiral bialgebra equivalence on the Koszul locus"; in build via `main.tex:1394`). AP263 Hopf-antipode caveat inscribed inline as `rem:chiral-bialgebra-not-hopf` at :8496; GRT_1 torsor content as `rem:chiral-qg-grt1-torsor` at :8522. Triangle (I) R-matrix ↔ (II) chiral $A_\infty$ in $\End^{\mathrm{ch}}_\cA$ ↔ (III) chiral coproduct coassociative up to $\Phi$; canonical at $\mathrm{H}^0$, non-canonical up to $\mathrm{GRT}_1(\bQ)$ at the cochain level. Load-bearing on JKL26 at $N \geq 3$ via `prop:ff-screening-coproduct-obstruction` (Feigin--Frenkel is not a drop-in replacement; the $(\Psi-1)/\Psi$ obstruction is explicit). gl_N extension via `thm:glN-chiral-qg` at :10324 (`\ClaimStatusProvedHere`). ELLIPTIC and TOROIDAL formal-disk: standalone-only inscriptions at `standalone/seven_faces.tex:1006,1020` (`thm:chiral-qg-equiv-elliptic-sf`, `thm:chiral-qg-equiv-toroidal-sf`); `seven_faces.tex` is not `\input`-ed into `main.tex`, so chapter-level inscription is a FRONTIER item for both. Of the 8 formerly advertised strengthening labels, 1 is inscribed as a theorem (`thm:glN-chiral-qg`), 1 is inscribed as a remark that absorbs its intended content (`rem:chiral-qg-grt1-torsor`, consumer retargeted 2026-04-18), and 6 are AP241 phantoms with ZERO live consumer `\ref` after the 2026-04-18 retarget (`def:ordered-koszul-chiral-algebra`, `prop:yangian-ordered-koszul`, `thm:chiral-qg-equiv-ordered`, `prop:sl2-yangian-triangle-concrete`, `thm:w-infty-chiral-qg-completed`, chapter-level elliptic/toroidal). The `\phantomsection\label{thm:grt1-rigidity}` stub at `chapters/frame/preface.tex` was retired 2026-04-18.** | Post-heal consumer map: all in-build `\ref{thm:chiral-qg-equiv}` consumers (~35 sites across `chapters/theory/ordered_associative_chiral_kd.tex`, `chapters/theory/e1_modular_koszul.tex`, `chapters/theory/algebraic_foundations.tex`, `chapters/theory/introduction.tex`, `chapters/examples/yangians_foundations.tex`, `chapters/frame/preface.tex`) resolve to the chapter theorem, no `[?]` at build. Two prior `\ref{thm:chiral-qg-equiv-ordered}` refs (at `chapters/examples/exceptional_yangian_koszul_duality_platonic.tex:532`, `chapters/theory/en_koszul_duality.tex:25`) retargeted to the unsuffixed base 2026-04-18 (the base theorem IS the ordered variant; the `-ordered` suffix was historical editorial duplication). One prior `\ref{thm:grt1-rigidity}` ref at `chapters/frame/programme_overview_platonic.tex:357` retargeted to `rem:chiral-qg-grt1-torsor` 2026-04-18. Concrete verification of the triangle beyond $\fsl_2$ Yangian + affine KM + gl_N is OPEN. Elliptic standalone ($\fsl_2$ leading-$\hbar$ via Felder R + KZB + Fay trisecant); toroidal formal-disk standalone (Drinfeld--Miki $U_{q,t}(\ddot{\fsl}_N)$ + SV CoHA). Global $\PP^1 \times \PP^1$ toroidal extension conditional on class-M chain-level topologisation. See `attack_heal_chiral_qg_20260418.md`. |
| Non-principal W | **Koszulness at generic level PROVED for every nilpotent** (Arakawa Kazhdan filtration, Inventiones 2007); **hook-type $r \leq N-3$ is the MC1 semisimple-Levi window, not a Koszulness window**, CONDITIONAL at theorem level — statement inscribed at `thm:n4-non-principal-hook` in `standalone/N4_mc4_completion.tex` with proof sketch only and no `\ClaimStatus` tag, downgrade to `\ClaimStatusConjectured` pending full proof body; subregular / minimal / screening-kernel Koszul-locus preservation all UNINSCRIBED at theorem level for general $N$; logarithmic $W(p)$ Massey FALSIFIED in naive form (Gurarie 1993 + Flohr 1996 give unbounded Massey in $H^\bullet_{\log}$); shadow-tower Massey in $H^\bullet_{\mathrm{reg}}$ OPEN per Wave-4 UCH+$W(p)$ heal | Sub-claim status (Wave-7 heal 2026-04-17): (i) **Koszulness scope**: at generic level $k$ Koszulness holds for every $\cW^k(\fg, f)$ via Arakawa Kazhdan filtration, independent of the hook window. The $r \leq N-3$ restriction is the MC1 / Whitehead-vanishing hypothesis (c) on semisimple Levi factor of $\fg^f$, NOT a Koszulness hypothesis. Prior status-table gloss conflated the two; healed. (ii) **Hook-type parabolic screening** (KRW03, Arakawa07): CONDITIONAL — proof sketch in `standalone/N4_mc4_completion.tex` (`thm:n4-non-principal-hook`) and `prop:hook-pbw` in `higher_genus_modular_koszul.tex:2031`; `\ClaimStatus` tag missing; either inscribe full proof and upgrade, or tag `\ClaimStatusConjectured`. KRW03 proves DS existence + central-charge formula; parabolic-screening PRESENTATION on Levi Heisenberg is a downstream corollary. Cite pattern: `\cite{KRW03}` for DS existence + `\cite{Arakawa07}` for Kazhdan filtration. (iii) **BP boundary case**: BP $= \cW^k(\fsl_3, f_{(2,1)})$ has $(N,r) = (3,1)$, OUTSIDE the hook corridor $r \leq N - 3 = 0$ for $\fsl_3$; analysed directly via `thm:w-bp-strict` and `standalone/bp_self_duality.tex:253-297` (`thm:bp-koszul-conductor-polynomial`, Arakawa convention). NOT a corollary of the hook-type theorem. (iv) **Subregular $\cW_n^{(2)} \cong \cW_{N-2} \times \fsl_2 \times \beta\gamma$**: inscribed at $N = 3$ via BP (where the decomposition degenerates); general-$N$ is a structural expectation, prose only in `standalone/N4_mc4_completion.tex:999-1004` as "reduction to class C SDR", no inscribed theorem; Feigin–Semikhatov 2004 attribution pending. (v) **Minimal W via coset**: inscribed at $N = 3$ via BP; general-$N$ coset description is CONJECTURAL, no inscribed theorem. (vi) **Screening-kernel Koszul-locus preservation**: PARENTHETICAL observation (single site `standalone/N4_mc4_completion.tex:996`); no theorem, proposition, or lemma inscribes this with a proof; `logarithmic_w_algebras.tex:339-342` states the OPPOSITE for general screening kernels. Previous "preservation proved" claim RETRACTED; scope restricted to hook partitions at generic $k$; beyond-hook content is subsumed by `conj:ds-kd-arbitrary-nilpotent`. Explicit inscription as a scoped remark is pending. (vii) **Logarithmic $\cW(p)$ Massey $\langle \Omega,\Omega,\Omega\rangle$ obstruction EXPLICIT**: RETRACTED. Keyword `Massey` is ABSENT from `logarithmic_w_algebras.tex` (685 lines); only programme-wide mention is a retracted paragraph in `shadow_tower_quadrichotomy_platonic.tex:1173-1188` citing Gurarie 1993 (`arXiv:hep-th/9303160`) and Flohr 1996 (`arXiv:hep-th/9605151`) for UNBOUNDED Massey products in logarithmic amplitudes. Consequence: the naive "$C_2$-cofinite $\Rightarrow$ bounded Massey $\Rightarrow$ tempered" implication (B91) is FALSIFIED in $H^\bullet_{\log}$, NOT made explicit; `logarithmic_w_algebras.tex` carries `conj:wp-koszulness` without explicit Massey form. **Correct open problem** (per Wave-4 UCH+$W(p)$ heal): the shadow-tower triple Massey class in $H^\bullet_{\mathrm{reg}}(B(\cW(p)))$ — regular sector, distinct from the logarithmic sector where unboundedness is already known — is OPEN; sharpest target case $p = 3$ (first non-symplectic-fermion). |
| gl_N chiral QG | **UNCONDITIONAL at $N=2$ via $\fsl_2$ Yangian + affine $\widehat{\fsl}_2$ RTT + SV CoHA; $N \geq 3$ CONDITIONAL on JKL26 vertex bialgebra on the Jordan quiver CoHA (Argument B), with Feigin--Frenkel screening descent PROVED OBSTRUCTED at the explicit $(\Psi-1)/\Psi$ cocycle (Proposition `prop:ff-screening-coproduct-obstruction`, `ordered_associative_chiral_kd.tex:10177`). Bialgebra level only; full Hopf lift is a PROVED NEGATIVE per AP263 (`rem:chiral-bialgebra-not-hopf`, `:8497`).** Inscribed at `thm:glN-chiral-qg`, `ordered_associative_chiral_kd.tex:10324` (`\ClaimStatusProvedHere`). | Honest accounting (2026-04-18 swarm, AP256+AP271+AP305 healing): prior row "UNCONDITIONAL all $N\geq 2$ via Feigin--Frenkel screening (JKL26 phantom eliminated)" RETRACTED. Inscribed Argument~B explicitly names JKL26 as the external input at $N \geq 3$: "the coproduct construction is *external* to the Heisenberg parent and currently relies on~\cite{JKL26}" (theorem body at `:10517-10519`). The obstruction proposition computes $[Q_{\alpha_i}, \Delta_z^{\mathfrak{h}}]$ as a non-exact chiral $1$-cocycle in $H^1_{\mathrm{ch}}(\mathfrak{h},\mathfrak{h}\otimes\mathfrak{h}[z,z^{-1}])$ with coefficient $(\Psi-1)/\Psi$; the prefactor matches the universal cross-term coefficient of `thm:miura-cross-universality`. The cocycle vanishes only at $\Psi = 1$ (free-boson point); at generic $\Psi$ the FF-screening descent genuinely fails. Argument~A (coderivation on the Koszul locus) covers all $N$ at the abstract Koszul-locus level; concrete three-structure determination (R-matrix / $A_\infty$ / ordered coproduct) is realised at $N = 2$ (sl_2 Yangian + affine KM) without external dependency, and at $N \geq 3$ via JKL26. **Lemma `qdet-central-all-N`** at `:10675` carries `\ClaimStatusProvedElsewhere`: centrality is *inherited* from Molev (`\cite{Molev07}`, Theorems~1.6.4 + 1.11.2) via the antisymmetriser-image argument; the chiral version extends Molev by showing $\Delta_z$ commutes with the antisymmetriser specialisation. Molev is external; the inheritance step is internal. Decreasing-column ordering is load-bearing at $N \geq 3$ (`rem:qdet-decreasing-ordering`, `:10727`); the increasing-index variant fails centrality at order $\Psi^2$ for $N \geq 3$. **Two-parameter gl_N RTT** matches `rem:e3-non-simple-gl-N` at `en_koszul_duality.tex:5576` with $\dim \Sym^2(\mathfrak{gl}_N^*)^{\mathfrak{gl}_N} = 2$ spanned by $\{\tr(XY),\,\tr(X)\tr(Y)\}$ (ad-invariant; verified by direct decomposition $\mathfrak{gl}_N = \mathfrak{sl}_N \oplus \mathfrak{z}$, $\dim\Sym^2(\mathfrak{sl}_N^*)^{\mathfrak{sl}_N} = 1$ for $N \geq 2$, plus the central form). **$\Psi = 0$ degenerate case** (free abelian chiral algebra, identity QG) and **$\Psi = 1$ free-boson point** (FF obstruction cocycle vanishes) inscribed. **Non-principal W hook-type** status unchanged — see Non-principal W row. See `adversarial_swarm_20260418/attack_heal_glN_qg_20260418.md`. |
| Verlinde recovery | **PROVED ELSEWHERE at boundary genera + citation-level at g≥2** (2026-04-18 Wave-Verlinde audit): genus-0 ($Z_0=1$ from $S$-matrix unitarity) and genus-1 ($Z_1=k+1$ via Zhu + Frenkel--Zhu) unconditional; genus-$g\geq 2$ is Verlinde (1988) + TUY (1989) factorisation combined with the Vol-I bar-complex identification of conformal blocks (`prop:conformal-blocks-bar`, `chapters/theory/chiral_modules.tex:540`); clauses (iv) handle attachment and (v) separating factorisation are CONDITIONAL on the modular-family extension of `thm:A-infinity-2` to the Deligne--Mumford boundary (Francis--Gaitsgory relative Ran GR17 + Mok log-FM sewing Mok25, cited not inscribed per `rem:A-infinity-2-modular-family-scope`). Status downgraded `ProvedHere` → `ProvedElsewhere` at `prop:verlinde-from-ordered` with honest scope remark `rem:verlinde-from-ordered-scope`. Bare "Theorem~A for the ordered bar complex on $\Ran^{\mathrm{ord}}$" reference at `higher_genus_modular_koszul.tex:33505` retargeted to `\ref{thm:A-infinity-2}`. The engine's "three independent code paths" (`verlinde_Z2_three_paths`) are algebraically identical rewrites of $\sum_j S_{0j}^{-2}$ (AP602 decorator-label vs test-body computation collapsed to a single trigonometric identity); genuine disjoint verification lives at boundary genera + Beauville--Laszlo algebro-geometric moduli-of-bundles interpretation. Vol-I-novel content: `rem:verlinde-shadow-truncation` (shadow tower as generic-level ancestor of Verlinde) and `thm:verlinde-polynomial-family`. Attack-and-heal note: `adversarial_swarm_20260418/attack_heal_verlinde_20260418.md`. |
| ker(av) formula | PROVED (any finite-dim V over C, dim V = d; scope WIDER than "simple g") | dim(ker(av_n))|_{V^{⊗n}} = d^n − C(n+d−1, d−1) via Reynolds projector 1/n! ∑ σ and stars-and-bars dim Sym^n V = C(n+d−1, d−1) (prop:ker-av-schur-weyl, `ordered_associative_chiral_kd.tex:6400`, `\ClaimStatusProvedHere`). Depends only on d = dim V, not on g or g-module structure; proof uses no Lie bracket. For any simple g with d-dim module V, Schur-Weyl refines into non-trivial Σ_n-isotypics: for gl_r and V = C^r classical; for general simple g the commutant is End_g(V^{⊗n}) rather than C[GL(V)], but total kernel dim is conserved. REPRESENTATION-THEORETIC UPPER BOUND on the ordered chiral-centre kernel (subquotient via bar cohomology): equality at n=2 sl_2 fundamental (1-dim Alt^2 survives as classical Yangian computation); strict inequality at n ≥ 3 in general (rem:ker-av-table, :6489). Falsification tests at (d, n) ∈ {2, 3} × {2, 3, 4} via explicit SVD basis construction in `compute/lib/averaging_kernel_explicit_engine.py::verify_kernel_dim`; test_d3_n4 at `test_averaging_kernel_explicit_engine.py:96` predicts 66, SVD confirms 66. Schur-Weyl at (d=3, n=3): 17 = 1 + 2·8 = Alt^3(C^3) ⊕ 2·S^{(2,1)}(C^3), verified by explicit basis-in-kernel. See `adversarial_swarm_20260418/wave4_ker_av_attack_heal.md`. |
| Genus-2 construction | CONSTRUCTED | KZB with 2x2 Siegel period matrix, chi=-12 at degree 2, doubly-dynamical (conj:g2-ddybe). |
| Miura coefficient | PROVED (thm:miura-cross-universality) | (Psi-1)/Psi on J⊗W_{s-1}+W_{s-1}⊗J at ALL spins s>=2. Three-step proof (Drinfeld binomial +1, Miura inversion −1/Ψ, lower J-sectors 0); verified spins 2-6, engine miura_universality_proof_engine.py. STRUCTURAL MIRROR (inscribed 2026-04-18 as rem:miura-ff-screening-structural-mirror, ordered_associative_chiral_kd.tex): the SAME coefficient (Ψ-1)/Ψ appears as the prefactor of the non-exact chiral 1-cocycle R_i(z) representing the Feigin-Frenkel screening / Drinfeld-coproduct obstruction (prop:ff-screening-coproduct-obstruction). Both derivations extract 1 − 1/Ψ: the primitive "+1" from the Heisenberg parent's trivial coproduct (respectively the Drinfeld binomial z^0 = 1), the subtraction "−1/Ψ" from the Miura/screening nonlinearity (respectively α_i^2 = 2/Ψ contraction weight). This is a structural identification, not a normalization coincidence. Falsification test at s ≥ 4: W_3 W-term and W_4 W-term must both exhibit the same (Ψ-1)/Ψ coefficient (engine verifies through s = 6). |
| Z_g closed forms | DISCOVERED g=0..9; PROVED all g via Hurwitz-Bernoulli; arithmetic duality at leading Kummer-irregular primes {691, 3617} PROVED through r=11 | P_g(n) = n^{g-1}(n²-1)·R_{g-2}(n²), n=k+2. Leading = B_{2g-2}/(g-1). Kummer-congruence prediction (691 at Z_7, 3617 at Z_9) is falsification test. Arithmetic duality `thm:z-g-s-r-arithmetic-duality` (z_g_kummer_bernoulli_platonic.tex): leading Kummer-irregular primes {691, 3617} present on Z_g (via B_{2g-2}) AND absent on S_r(Vir_c) through r=11 `thm:s-r-kummer-absent-through-r-11` (shadow_tower_higher_coefficients.tex). Sharp at {691, 3617}: the higher Kummer-irregular prime 2111 appears in N_9 (via 29554 = 2·7·2111), but is not witnessed on Z_g at g≤9 since 2111 ∤ num(B_{2m}) for 2m≤16. Geometric disjointness: Z_g localises on Bun_{SL_2}(Σ_g), S_r on ChirHoch(Vir_c). Characteristic S_r primes include {61, 193, 2111, 16657, 3988097}. HZ-IV decorator at compute/tests/test_z_g_s_r_arithmetic_duality.py. |
| W_N Stokes count | **COMPUTED via standard irregular-sing theory (ProvedElsewhere)** | Stokes rays for $\cW_N$ KZ connection at degree 2 = $4N-4$ (linear in $N$). $\cW_2$: 4 rays; $\cW_3$: 8; $\cW_4$: 12. Four-step classical chain: (a) $W_N W_N$ OPE pole order $2N$ (Belavin-Polyakov-Zamolodchikov 1984, Zamolodchikov 1985 for $W_3$); (b) d-log absorption gives connection-form pole $2N-1$ (Arnold 1969); (c) Poincaré rank $r = (\text{pole order}) - 1 = 2N-2$ (Wasow 1965, Sibuya 1975); (d) Stokes ray count $= 2r = 4N-4$ under generic anti-Stokes / Stokes-sector convention (Boalch 2002, Jimbo-Miwa-Ueno 1981, Hertling-Sabbah 2005). Programme-novel content is the identification of Poincaré rank $2N-2$ with the $\cW_N$ chiral algebra at degree 2, a one-line corollary of the pole-order computation. Inscribed as `rem:stokes-wN-monograph` at `chapters/theory/ordered_associative_chiral_kd.tex:9985-9993` with convention pointer (Wave-7 2026-04-17 heal). Independent Virasoro engine `compute/lib/irregular_kz_stokes_virasoro_engine.py` computes $N=2$ case from first principles (Poincaré rank 2, 4 Stokes rays at $\theta = \pi/4, 3\pi/4, 5\pi/4, 7\pi/4$, 4 anti-Stokes rays at $\theta = 0, \pi/2, \pi, 3\pi/2$). Prior "DISCOVERED" label DOWNGRADED 2026-04-18 per AP1301 (classical-chain-masquerading-as-discovery, AP280 sibling). Full audit: `adversarial_swarm_20260418/attack_heal_w_n_stokes_20260418.md`. |
| Shadow = GW(C³) | IDENTIFIED | Shadow tower at κ=Ψ produces perturbative constant-map GW free energies F_g^{GW,const}(C³). MacMahon M(q) on DT side via MNOP. |
| Conformal anomaly | QUANTIFIED | Obstruction to constant coproduct = c/2 = κ(Vir_c). At c=0: obstruction vanishes. At c≠0: spectral parameter FORCED. |
| Chiral Higher Deligne | **PROVED clause (1) + (3); clause (2) CONDITIONAL** (Vol II `chapters/theory/chiral_higher_deligne.tex`, 2026-04-16; Wave-6 split 2026-04-17; Wave HDC-heal 2026-04-18 AP2021 downgrades below) | `thm:chiral-higher-deligne` (chiral_higher_deligne.tex:460, clause (1) + (3) \ClaimStatusProvedHere, clause (2) \ClaimStatusConditional per :459): Z^{der}_ch(A) is universal E_3-chiral algebra acted on by SC^{ch,top} via heptagon edges (3)↔(4)↔(5); clause (2) universality conditional on the $(\SCchtop)^!$-cobar contracting homotopy (cited, not inscribed). **`conj:H-concentration-via-E3-rigidity`** (:650, \ClaimStatusConjectured, retracted 2026-04-18 from \ClaimStatusProvedHere per native-vs-derived + E_1-vs-E_3 type errors — see `rem:H-concentration-via-E3-rigidity-type-error-retraction` at :673): the Step-1 assertion that $\cA|_{D_x}$ carries $E_3$ was a native/derived error (Thm~\ref{thm:chiral-higher-deligne}(1) places $E_3$ on $Z^{\mathrm{der}}_{\mathrm{ch}}(\cA)$, not on $\cA$), and Step 2's PBW identification $\cA|_{D_x} = U_{E_3}(V_x)$ was an $E_1$-to-$E_3$ envelope error (chiral PBW yields the $E_1$-chiral envelope); replacing $\cA|_{D_x}$ by its derived centre fixes Step 1 but breaks Step 2, and no chiral-$E_3$-PBW theorem presenting $Z^{\mathrm{der}}_{\mathrm{ch}}(\cA|_{D_x}) = U_{E_3}(W_x)$ is inscribed in any volume. Theorem H itself is unaffected: the parallel Vol~I `thm:hochschild-concentration-E1` (chiral PBW + Shelton--Yuzvinsky Orlik--Solomon Koszulity) supplies the canonical unconditional proof. `thm:chd-ds-hochschild` (:740, **split 2026-04-18**: $E_2$-chiral brace quasi-iso \ClaimStatusProvedHere; chain-level $E_3$-chiral lift \ClaimStatusConditional inheriting conditionality of `thm:chiral-higher-deligne`(2)): ChirHoch^•(W_k(g)) ≃ H^•_DS(ChirHoch^•(V_k(g))) at chain-level $E_2$-chiral Gerstenhaber; $E_3$-chiral extension is unconditional on cohomology, conditional at chain level on the cited $(\SCchtop)^!$-cobar contracting homotopy. `cor:universal-holography-class-M` (:837, **split 2026-04-18**: cohomological \ClaimStatusProvedHere; chain-level $E_3$-chiral \ClaimStatusConditional inheriting from `thm:chd-ds-hochschild`): class-M chain-level 3d HT holography identification at $E_2$-chiral/Gerstenhaber level is unconditional on the weight-completed ambient; $E_3$-chiral identification inherits conditionality. Audit trail: `adversarial_swarm_20260418/attack_heal_chiral_hdc_type_error_heal_20260418.md`. |
| Curved-Dunn H²=0 at g≥2 | **PROVED (Vol II `chapters/theory/curved_dunn_higher_genus.tex`, 2026-04-16)** | `prop:modular-bootstrap-to-curved-dunn-bridge`: chain map Φ on H²; `prop:genus1-twisted-tensor-product`: explicit Gauss–Manin uncurving + Arakelov pairing (phantom FM87 resolved); `thm:curved-dunn-H2-vanishing-all-genera`: curved-Dunn H²=0 ∀g≥2 via bridge from modular-bootstrap H²=0. `thm:irregular-singular-kzb-regularity`: Jimbo–Miwa replaces KZ Stokes, closes modular operad composition at generic non-integral level. Closes FM67/FM88/FM91/FM92/FM192/FM215. |
| SC^{ch,top} heptagon | **PROVED (Vol II `chapters/theory/sc_chtop_heptagon.tex`, 2026-04-16)** | Five classical faces + face (6) Drinfeld-centre `Z(Rep_fact(A)) ≃ Rep_fact(Z^der_ch(A))^{E_2}` via categorified bar-cobar with half-braiding (AP-CY25) + face (7) derived-AG via PTVV on `Map(X×R_≥0, B SC-Alg)`. Seven named edge theorems prove the closed heptagon. SC^{ch,top} is the GENERIC case; topologization to E_3-top at affine KM non-critical is proved, general chiral conjectural. |
| Universal celestial holography | **PROVED chain-level (Vol II `chapters/connections/universal_celestial_holography.tex`, 2026-04-16)** | `thm:uch-main`: SC^{ch,top}-structure on `(A^cel, Z^der_ch(A^cel))` + celestial OPE = chiral factorization homology on `P^1_cel` + shadow-tower coefficients = soft-factor hierarchy. Coverage: self-dual gauge (KM), gauge+matter (DS), gravity (Virasoro + w_{1+∞}), YM (Beem–Rastelli χ-functor). Closes FM102/FM103. Chain-level class M at g≥1 inscribed as `thm:uch-gravity-chain-level` (Vol II universal_celestial_holography.tex:511, ClaimStatusProvedHere), promoted from conjecture 2026-04-16 via half-BRST chain-level splitting (see rem:uch-gravity-chain-level-promotion at :579). |
| Periodic CDG admissible KL | **PROVED (Vol I `chapters/theory/periodic_cdg_admissible.tex`, 2026-04-16)** | `thm:periodic-cdg-is-koszul-compatible`: periodic-CDG filtration `F^n = ker(Q^n_{adm})` on `KL_k^{adm}` at `k+h^v = p/q` compatible with chiral Koszul duality. `thm:admissible-kl-bar-cobar-adjunction`: `Ω^ch ⊣ B^ch` descends unconditionally to `KL_k^{adm} ⇄ (KL_{k^!}^{adm})^op`; proof via Arakawa C_2-cofiniteness + periodic-CDG finite length + Adams-type spectral sequence with `d_1 = Q_{adm}`. `thm:adams-analog-construction`: chiral Steenrod algebra `A^{ch}_k` + chiral Adams functor; closes FM256. `cor:class-M-admissible-minimal-model`: `KL^{adm}_{Vir}(c_{p,q})` has `(p-1)(q-1)/2` simples, all finite-length (closes FM76 scope hole). Closes FM251 + the sole remaining irreducible-open programme frontier. |
| E_∞-Topologization | **PROVED (Vol II `chapters/connections/e_infinity_topologization.tex`, 2026-04-16)** | `thm:iterated-sugawara-construction`: higher-spin Casimir tower `{T^{(n)}}_{2 ≤ n ≤ N+1}` each inner, admitting BRST primitive G^{(n)} with `T^{(n)} = [Q_tot, G^{(n)}]` on cohomology. `thm:e-infinity-topologization-ladder`: k inner stress tensors ⟹ E_{k+2}-top via Dunn with E_2-chiral ⊗ E_1-top(T^{(n_1)}) ⊗ ... ⊗ E_1-top(T^{(n_k)}). Specializations: Virasoro (N=2) → E_3-top (unconditional); W_3 depth-3 → E_4-top (unconditional); W_N for $N \ge 4$ → E_{N+1}-top CONDITIONAL on axiom (T5) (antighost BRST-commutativity at all spins $n,m \ge 4$, postulated, see rem:axiom-T5-scope at e_infinity_topologization.tex:343–367); W_∞ → E_∞-top CONJECTURAL as `conj:e-infinity-specialisation-Winfty` (e_infinity_topologization.tex:679, ClaimStatusConjectured; unconditional at depth $N \le 3$; $N \ge 4$ conditional on T5). `thm:operadic-spiral`: infinite bidirectional spiral with bar B descending, center Z ascending, meeting at E_∞. Closes FM47/48/81/82/215 (W_N hook-type via branched cover, class M free-PVA via Li-filtration). Climax restatement: 3d quantum gravity Vol II Part VI is N=2 shadow of a 3d+∞ topological theory. |
| Theorem A^{∞,2} | **PROVED (Vol I `chapters/theory/theorem_A_infinity_2.tex`, 2026-04-16)** | `thm:A-infinity-2`: Francis-Gaitsgory bar-cobar (∞,2)-equivalence at properad level. `B̄^ch_X : Fact^{aug}(X) ⇄ CoFact^{conil,comp}(X) : Ω^ch_X` (∞,2)-adjoint equivalence on conilpotent-complete locus. Three clauses: (i) properad lift via Hackney-Robertson in FG ambient; (ii) pole-free restriction recovers LV12 (Ass, Ass^!) via `(D_X-mod, ⊗^!) ↪ Fact(X)`; (iii) R-twisted Σ_n-descent (`lem:R-twisted-descent`) relates `B^{ord}(A)` to `B^Σ(A)`. Closes FM69/70/72/73/74/195 (ambient category corrections, cross-volume citation drift). 14+ downstream corollaries enumerated (classical Theorem A, bar-cobar adjunction, Vol II bridge I.3.2, Vol III Φ, ...). `cor:chiral-KK-formal-smoothness`: FG ambient + R-twisted descent ⟹ formally smooth at properad level. |
| CY-D dimension stratification | **PROVED (Vol III `chapters/examples/cy_d_kappa_stratification.tex`, 2026-04-16)** | `thm:kappa-hodge-supertrace-identification`: `κ_ch(A_X) = Σ_q (-1)^q h^{0,q}(X)` unconditionally for compact CY_d via HKR + Mukai pairing + HC^-_d trace. `thm:kappa-stratification-by-d`: explicit values across d ∈ {1,2,3,4,5} — E(0), K3(2), abelian/bielliptic(0), quintic/K3×E/E³(0), local P²(3/2 via `thm:local-p2-shadow`), CY_4 sextic(2), CY_5 generic(0). `cor:conifold-non-local-surface` closes AP-CY34/AP-CY44 (conifold is NOT local surface at d=3; κ_ch=1 via direct McKay). `thm:borcherds-weight-kappa-BKM-universal`: κ_BKM(Φ_N) = c_N(0)/2 universal across N ∈ {1,2,3,4,6}; at N=1 this gives κ_BKM(Φ_1) = 10/2 = 5 via Gritsenko's Δ_5 weight-5 paramodular form of level 1 (Gritsenko 1999). The naive decomposition κ_BKM = κ_ch + χ(O_fiber) holds at NO N (already fails at N=1: 5 ≠ 0); the prior "N=1 coincidence" narrative is retracted as confabulation (HEAL 2026-04-17). Closes AP-CY37 + top-15 #15. |
| BP Koszul-conductor polynomial identity | **PROVED IN ARAKAWA CONVENTION (Vol I `standalone/bp_self_duality.tex:253-297`, 2026-04-16)** | `thm:bp-koszul-conductor-polynomial`: in Arakawa convention `c(BP_k) = 2 - 24(k+1)²/(k+3)`, `K_BP(k) := c(BP_k) + c(BP_{-k-6}) ≡ 196 ∈ Q(k)` constant rational function; `c(BP_k) - 98` is odd function of `(k+3)`. Fixed point `k = -3` coincides with critical level `-h^v(sl_3)`; `κ(BP_{-3}) = 49/3` is principal-value symmetric limit. 73+7 tests pass with 2 disjoint HZ-IV decorators. **Convention caveat (Vol II cross-check, 2026-04-17, CORRECTED)**: Fateev-Lukyanov screening convention `c^{FL}(k) = -(2k+3)(3k+1)/(k+3)` gives `K^{FL}_BP(k) = c^{FL}(k) + c^{FL}(-k-6) = 50(k+3)/(k+3) ≡ 50 ∈ Q(k)` — also polynomial-constant (pole at k=-3 is removable; residues cancel). Both conventions give CONSTANT conductors; they differ only in numerical value (196 Arakawa vs 50 FL) due to central-charge normalization. Prior CLAUDE.md entry "meromorphic with pole at k=-3, NOT polynomial" RETRACTED per agent a696c7f6591aa6e78 audit. Vol II `bp_chain_level_strict_platonic.tex` and `fm81_fractional_ghost_platonic.tex` healed 2026-04-17 to match. |
| Critical level jump | PROVED | At k=-h^v for sl_2: kappa=0, monodromy trivial (integer residues). Two DISTINCT cohomological objects: (i) degree-2 ORDERED chiral homology on V⊗V has dim H^1_{ord,2} doubling 4 → 8 (finite, chi_{ord,2} = -4 conserved; `prop:critical-level-ordered` at `ordered_associative_chiral_kd.tex:12130`); (ii) FULL symmetric bar cohomology H^n(B^Sigma(V_{-h^v}(g))) = Omega^n(Op_{g^v}(D)) INFINITE-dimensional in every degree n >= 0 (Feigin-Frenkel 1992, Frenkel-Teleman 2006; `thm:oper-bar` at `kac_moody.tex:4402`). Generic-level bar dim H^1(B^Sigma(V_k(sl_2))) = 3 = dim sl_2 is a THIRD object (symmetric bar, generic level; `bar_construction.tex:1633`). Koszulness fails for the full symmetric bar at critical level. Five standalones previously conflated (i) and (ii) via unqualified "H^1 doubles 4 → 8" (healed 2026-04-18 per AP2041-AP2043; attack_heal_h1_sl2_critical_level_20260418.md). 72 tests. |
| Genus-2 degree decomp | PROVED | CB_{2,2}(k) = 2k(k+1)(k+2)/3 (cubic). At k=1: 4. Singlet+triplet channels (prop:g2-conformal-block-degree). |
| Antipode non-lifting | PROVED (negative) | S(T(u))=T(u)^{-1} does NOT lift to vertex-algebraic antipode. Two obstructions (OPE and Hopf axiom). |
| DS intertwining | VERIFIED at degree ≥ 2 (non-trivial); degree 1 tautological | (pi_3×pi_3)∘Delta_z^{sl_3} = Delta_z^{W_3}∘pi_3 verified 57 tests. Spectral coassociativity uses shifted parameters. Engine docstring `compute/lib/ds_coproduct_intertwining_engine.py:103-104` records degree-1 as tautological (pi_3 acts as identity on the W_3 generators at leading order); the content is the degree-≥2 sl_3 RTT computation (AP256/AP257 discipline). |
| AP128 bar H^2 | FIXED | sl2_bar_dims gave h_2=6 (CE/Riordan). Correct chiral bar: h_2=5. AP63: Orlik-Solomon form factor. |
| Quantum det ordering | FOUND | Central qdet uses DECREASING column index. At N=3, increasing-index NOT central. 74 tests. |
| E_3 via Dunn | PROVED (alt) | prop:e3-via-dunn: CG factorization E_1^top×E_2^hol + Sugawara + Dunn = E_3^top. Independent of HDC. |
| E_3 for gl_N | EXTENDED | E_3 identification extends to gl_N via two invariant bilinear forms B_tr, B_ab. |
| 6d hCS defect | **PROVED at spins 1+2 on C ⊂ C³ with gl_1 gauge (cohomological; not chain-level in class M)** via Vol III `prop:codim2-defect-ope` (quantum_chiral_algebras.tex:397). Full higher-spin W_{1+∞} identification is **ProvedElsewhere** via Procházka–Rapčák (arXiv:1711.06888) Miura transform + Arbesfeld–Schiffmann–Vasserot (arXiv:1506.00246) level formula + the positive-half CoHA correction `CoHA(C³) ≃ Y^+(widehat(gl_1))` (Vol I standalone `cy_quantum_groups_6d_hcs.tex:1128`, `thm:wave14-CoHA-positive-half`); the slogan "CoHA = W_{1+∞}" is a character-level coincidence, not a bar-level equivalence. | Codim-2 defect: spin-1 Heisenberg at level Ψ = −σ_2(h_1,h_2,h_3) with Ψ ≥ 0 under CY σ_1 = 0 (Newton: Ψ = (1/2)·Σh_i²), Sugawara c = 1 **only at gl_1 (h^∨=0; c = N for gl_N)**. N_{C/Y} = O ⊕ O gives two spectral parameters (h_1, h_2) from normal-bundle FIBER (not sections). The −σ_2 coefficient arises from the equivariant 1-loop determinant on the formal bidisk D² (ProvedElsewhere via Costello 2013; 5d→6d-on-codim-2 bridge cited, not inscribed as a Vol III lemma — AP621). 48 engine tests + 1 HZ-IV-decorated independent-verification test (PR S_3 triality + ASV Newton-identity paths genuinely numerical 2026-04-18, per H3 heal). Global P¹×P¹ toroidal extension is CONDITIONAL (AP-TOPOLOGIZATION, class-M chain-level). |
| DDYBE face model | **NUMERICAL EVIDENCE** + **CONJECTURED** at generic Ω (Platonic form, Vol I `chapters/theory/genus_2_ddybe_platonic.tex`, 2026-04-16) | Face-model strategy bypasses vertex-IRF (the latter not established at g=2). Tolerance ladder T1(10^{-12})/T2(10^{-10})/T3(10^{-6})/T4(10^{-4}). Generic-Ω DDYBE: 4 tests at T4 (10^{-4} relative) — test_ddybe_generic_omega, test_ddybe_different_spectral, test_ddybe_different_eta, test_ddybe_strong_coupling — plus 1 diagonal-Ω tier-(T3) consistency check (test_ddybe_diagonal_omega) reducing on the nose to a SINGLE genus-1 Felder DYBE (the τ_2-sector entering only through the cancelling θ_3(0|τ_2) normalisation; "two genus-1 DYBE copies" retracted, the diagonal regime is vacuously doubly-dynamical); separating degeneration AP157-empty; non-separating degeneration (Ω_{22}→i∞, Ω_{12} fixed nonzero) NOT tested — genuine open frontier. Full 30 face-model tests span θ_1 antisymmetry (T1, structural identities in numpy — AP287 primitive-tautology tier), weight-zero matrix entries (T2, θ-series N=30 truncation), classical limit + degeneration (T3, N=8 + algebraic cancellation), generic-Ω DDYBE (T4, N=8 no cancellation — engineering budget, NOT theoretical error term). `conj:g2-ddybe` remains ClaimStatusConjectured — finite-ℏ commutativity of doubly-dynamical Casimirs is the residual conjecture (infinitesimal proved via heat-equation symmetry). "Fay trisecant extends to g=2" resolved: Szegő three-term identity (Fay 1973 Cor. 2.5) holds universally at all g ≥ 0; distinct from theta-nulls quartic identity on theta-divisor. Cache entries #220-227. |
| ChirHoch*(H_k) explicit | **PROVED** | `prop:derived-center-explicit` (chiral_center_theorem.tex:1967, ClaimStatusProvedHere): $\ChirHoch^*(\cH_k) = (\bC, \bC, \bC, 0, \ldots)$ with generators $\{1, \xi_k, \eta\}$; $P_{\cH_k}(t) = 1 + t + t^2$. Naive Drinfeld center dim 1; derived chiral center dim 3. 72+ tests. |
| Drinfeld double vs derived center (Heis) | **CONJECTURED** | `conj:v1-drinfeld-center-equals-bulk` (preface.tex:4267, AP271 line-drift 2026-04-18 from previously advertised :4228): $Z(U_{\cH_k}) \simeq \cZ^{\mathrm{der}}_{\mathrm{ch}}(\cH_k)$ as $E_3^{\mathrm{top}}$-algebras, conditional on topologization hypotheses of preface §10.2. Three obstructions flagged (pointwise reduction at class M; Verdier dual factorisation; bar-cobar/Hochschild compatibility). |
| Toroidal chiral QG (formal disk) | PROVED on formal disk | `thm:chiral-qg-equiv-toroidal-sf` (standalone/seven_faces.tex:1043, AP271 line-drift 2026-04-18 from previously advertised :1020; \ClaimStatusProvedElsewhere): Ding--Iohara--Miki toroidal algebra $U_{q,t}(\ddot{\fsl}_N)$ with two spectral parameters $(u,v)$ is the chiral bialgebra (AP263) of an ordered-Koszul factorization algebra on the formal disk $(\mathrm{Spec}\,\bC[\![z]\!])^{\times 2}$, witnessed via SV CoHA on instanton moduli. Global extension to $\PP^1 \times \PP^1$ CONDITIONAL on class-M chain-level topologisation frontier. Prior companion label `thm:chiral-qg-equiv-toroidal-formal-disk` was PHANTOM (AP255; zero `\label{}` inscriptions); `verb|...|` prose reference at `standalone/N3_e1_primacy.tex:1171` retargeted 2026-04-18 — see Toroidal chiral QG row above for full audit. |
| Coderived E_3 | PARTIAL | Steps 1-2 proved (D^co stable ∞-cat; obstruction coacyclic). Step 3 open. |
| KZB flatness | **2-point SCALAR Bernard identity INSCRIBED + VERIFIED**; matrix commutator trivially zero; n≥3 integrability FRONTIER (Halphen--Ramanujan + elliptic CYBE with dynamical shift) | Prefactor 1/(4πi) diagonal vs 1/(2πi) off-diagonal verified < 10^{-14} at `rem:kzb-heat-prefactor` (ordered_associative_chiral_kd.tex:12301, AP271 line-drift 2026-04-18 from previously advertised :12260). Full scalar Bernard heat identity `4πi ∂_τ ζ = -℘' + 2(ζ-2η_1 z)(-℘-2η_1) + 8πi η_1'(τ) z` inscribed as `rem:bernard-heat-identity-zeta` (ordered_associative_chiral_kd.tex:12317, line-drift 2026-04-18 from previously advertised :12276; `\ClaimStatusProvedElsewhere` attributing Bernard 1988 eq. 2.15 + Ramanujan 1916 + Halphen 1886 + Felder ICM 1994 + EFK98 Ch. 5). Engine `verify_bernard_heat_identity_zeta` + test `test_bernard_heat_identity_zeta_scalar` replace the prior xfail test 2026-04-18. 2-point matrix commutator $[A_z, A_\tau] = 0$ is trivial (Ω self-commutes); the scalar Bernard identity is the entire 2-point content. Naive "∂_τ ζ = ℘'" form drops THREE corrections: heat-equation prefactor 1/(4πi), nonlinear $(\zeta - 2\eta_1 z)(-\wp - 2\eta_1)$ cross-term, quasi-period drift $8\pi i \eta_1'(\tau) z$. Halphen--Ramanujan closure $\eta_1(\tau) = (\pi^2/6) E_2(\tau)$, $\partial_\tau E_2 = (E_2^2 - E_4)/12$, etc., enters at n≥3 and via $\partial_\tau \eta_1$ at n=2. Chain-level KZB flatness for the full n-point configuration space CONJECTURAL (`conj:trig-elliptic-ordered`, ordered_associative_chiral_kd.tex:12075). Prior formula "d_tau(wp_1) = (1/(4πi))d_w(wp+wp²)" retracted 2026-04-17 per AP283. See `adversarial_swarm_20260418/attack_heal_kzb_flatness_20260418.md` (AP521-AP524). |

## Anti-Patterns by Topic

### Epistemic
AP2: Read actual .tex proof, not CLAUDE.md description. Descriptions are claims ABOUT source.
AP4: ClaimStatusProvedHere = verify proof proves stated claim. Status tag != ground truth.
AP11: Single-point external dependency -> flag in concordance with source/status/fallback.
AP13: Forward references transparent about genus/level/type restrictions.
AP15: E_2* is quasi-modular. Genus-1 propagator IS E_2*. Graph sums produce {E_2*,E_4,E_6}, NOT {E_4,E_6}.
AP17: After writing ANY new theorem, IMMEDIATELY audit before building next result.
AP26: Fock inner product != BPZ for weight>=4, rank>=3 W-algebras.
AP28: NEVER use undefined terminological qualifier in 3+ locations.
AP35: False proof, true conclusion -> two cancelling errors. Fix BOTH.
AP37: SS page from FULL differential, not pole order heuristics. Lie homology != Hochschild homology.
AP38: Literature values: record source paper AND convention. DVV != Eichler-Zagier. Derive independently.
AP39/AP127: S_2 is a lane coefficient, not a synonym for kappa. Unqualified S_2=kappa is allowed only on scalar one-line lanes: Vir stress tensor (S_2=c/2) and Heis current line (S_2=k). For KM/W_N/lattice/bc/betagamma write S_2^T, S_2^J, etc.; the family kappa is the genus-1 obstruction after the chosen projection. KM: kappa=dim(g)(k+h^v)/(2h^v), not c/2.
AP41: Prose mechanism != mathematical mechanism. "Residue extracts simple-pole coefficient" WRONG.
AP42: State level of validity explicitly for sophisticated identifications.
AP43: Central object without \begin{definition} -> property list is conjecture, not definition.
AP47: Evaluation-generated core != full category. MC3 proved on eval core; DK-4/5 downstream.
AP60: Tag only genuinely new content ProvedHere. Classical parts ProvedElsewhere with attribution.
AP147: Circular proof routing. When proof chain involves multiple theorems referencing each other, insert ROUTING REMARK citing the primitive non-circular anchor. If none exists, proof is genuinely circular.
AP149: Resolution propagation failure. When conjecture proved/disproved/retracted, ALL references retain old status unless updated: (a) concordance, (b) preface, (c) introduction, (d) standalones, (e) theorem status table, (f) label prefixes, (g) other volumes. All updates SAME session.
AP150: Agent confabulation of mathematical structures. Every claimed multi-step structure verified arrow-by-arrow; each arrow needs independent theorem reference. Any conjectural arrow → structure is conjectural.
AP155: "New invariant" overclaiming. Architectural novelty (new framework) != computational novelty (new numbers). Check Bernard/Felder/Etingof-Varchenko/Calaque-Enriquez-Etingof for recovered invariants.
AP157: Degeneration-dependent "invariants." Formula from specific degeneration is NOT invariant unless degeneration-independence proved. Separating degeneration of genus-2 has ZERO genuinely genus-2 info.

### Computational
AP3/AP10/AP61 (discipline): Compute independently. NEVER pattern-match across occurrences. Every hardcoded expected value MUST have comment citing 2+ independent derivation paths. Bare numbers without derivation trail = future AP10 violations.
AP6: Specify genus, degree, level (convolution vs ambient) for D^2=0, kappa, Theta_A.
AP7: Before universal quantifier, verify proof has no implicit type/genus/level restriction.
AP8: NEVER "self-dual" unqualified. Virasoro self-dual at c=13.
AP14: Koszulness != SC formality. Koszul = bar H* in degree 1. SC formal = m_k^{SC}=0 for k>=3. All standard families Koszul; only class G SC-formal.
AP18: "Entire standard landscape" -> list every family, check each against hypotheses.
AP29/AP31/AP33: H_k^! = Sym^ch(V*) != H_{-k}. Same kappa, different algebras. delta_kappa=kappa-kappa' (asymmetry, vanishes c=13) != kappa_eff=kappa(matter)+kappa(ghost) (cancellation, vanishes c=26). kappa=0 implies m_0=0 (uncurved); higher-degree independent. F_1=0 does NOT imply F_g=0.
AP30: CohFT flat identity requires vacuum in V. ALWAYS list conditional axioms.
AP36: "implies" proved, "iff" claimed -> write "implies" until converse has independent proof.
AP59: Three invariants: p_max(OPE pole) != k_max(collision depth=p-1) != r_max(shadow depth). betagamma: p=1,k=0,r=4.
AP62: "depends only on dim(g)" = Euler char only. Individual bar cohom dims need full bracket.
AP63: CE(g_-) != chiral bar for multi-gen. Orlik-Solomon form factor. sl_3 chiral H^2=36 vs CE H^2=20.
AP64: CE weight vs PBW degree produce different sequences. Always specify grading.
AP66: Free-field GFs NOT D-finite. Interacting algebras ARE.
AP67: Strong gen != FREE strong gen. W(p) has 4 strong generators; FREE strong gen OPEN.
AP68: PVA slab ghost c != kappa. SVir: kappa=(3c-2)/4. Hierarchy: Sigma_0=13 > Sigma_1=41/4 > Sigma_2=1 > Sigma_4=0.
AP69: tau_shadow = kappa-deformed KdV. Obstruction kappa(kappa-1). Standard KdV only at kappa=1.
AP70: Shadow L^sh has POLES at s=1,2. Negative integers are trivial zeros. F_g <-> L^sh(1-2g) FAILS.
AP71: Shadow kappa != Dyson beta. At c=13, kappa=6.5, not 13.
AP72: W-algebra NOP bar has d^2 != 0. Needs full singular OPE + Orlik-Solomon.
AP73: BV=bar: PROVED G/L, CONDITIONAL C/M (harmonic decoupling).
AP74: False Bernoulli-Dirichlet identity. LHS entire, RHS has poles s=1,2. Verify numerically at s=0.
AP75: Koszulness = PBW degree concentration, NOT conformal weight grading.
AP76: Y_{1,1,1} has c=0, kappa=Psi, NOT c=3.
AP77: Stokes ratio on convergent series spurious. Use direct Pade, not Borel.
AP78: Never conjecture from isolated number-theoretic coincidences.
AP79: W(p) has 4 generators (T + sl_2 triplet), not 2.
AP80: Engine without test file -> verify BOTH exist after every agent completion.
AP133: Catalan index shift. C_n counts binary trees with n+1 leaves (n internal nodes). Verify: count leaves, subtract 1.
AP135: q-expansion coefficients. 1/eta(tau)^r has r-coloured partition numbers p_{-r}(n). r=2: (1,2,5,10,20,...). OEIS lookup before hardcoding.
AP140: Koszul conductor vs local constant. K=c+c' is GLOBAL invariant. Ghost numbers/grading shifts are LOCAL. K_BP=196, not 2.
AP143: DS ghost charge background shift omission. c_ghost(N,f,k) = c(sl_N,k) - c(W_{N,f},k); simplified N*(N-1) OMITS background charge. At N=7: 1722 vs 42. Before any DS ghost formula for N>=3, compute c(sl_N,k) - c(W_{N,f},k) directly from Fateev-Lukyanov.
AP178: S_4 large-c asymptotic off by factor 5. 10/[c(5c+22)] ~ 2/c^2 at large c, NOT 2/(5c^2). Counter: substitute c=100 and verify numerically.

### Empirical (AP116-AP186)
AP116: sum_{j=a}^{b}: substitute smallest index. H_N=sum_{j=1}^{N}, NOT N-1. Off-by-one = #1 error.
AP117: Connection 1-form = r(z)dz, NOT d log(z). KZ=sum r_{ij} dz_{ij}. Arnold d log(z_i-z_j) is bar coefficient, not connection.
AP118: Genus-1 scalar collapse. (Im Omega)^{-1} scalar at g=1, matrix at g>=2. Write in full matrix form; verify at g=2.
AP119: Before Borel: verify Gevrey-1 (factorial divergence). If |F_{g+1}/F_g|→const, series is Gevrey-0; use direct Pade.
AP120: Cauchy = 1/(2πi), NOT 1/(2π). Verify F_1=κ/24.
AP121: In LaTeX, NEVER Markdown (`29`→$29$, **bold**→\textbf, _italic_→\emph). Grep backticks after every write.
AP122: Test tolerance relative, not absolute. abs(computed-expected)/abs(expected) < rtol.
AP123: Combinatorial counts verified against known formula/GF BEFORE hardcoding. Genus-2 stable graphs=7 (not 6).
AP129: Ratio a/b most common transcription error is b/a or -b/a. Substitute known value and check NUMBER.
AP130: Objects on fiber (dτ ∈ H^{1,0}(E_τ)) != objects on base (c_1(λ) ∈ H^2(M̄_g)).
AP131: d_gen (finite) != d_alg (∞ for Vir). Always specify WHICH depth.
AP132: B(A)=T^c(s^{-1} Ā), Ā=ker(ε). NOT T^c(s^{-1} A). Mnemonic: bar uses bar A.
AP134: Cohomological amplitude [0,2] (topological) != virtual dimension (arithmetic).
AP137: c_{βγ}(λ)=2(6λ²-6λ+1) vs c_{bc}(λ)=1-3(2λ-1)². VERIFY c_total=0. At λ=1: c_{βγ}=2, c_{bc}=-2.
AP138: At even ||m||=0, [[m,m],f]=0 tautological. Identity [m,[m,f]]=½[[m,m],f] requires ||m|| ODD. Check parity.
AP139: Every variable in displayed equation must be quantified. LHS ⊇ RHS.
AP144: When multiple conventions coexist, BRIDGE IDENTITY stated at every site; boundary behavior checked in EACH.
AP145: Restructuring = O(N) debt. BEFORE: grep all three volumes; AFTER: verify every ref resolves.
AP146: After 100+ agent campaigns, expect follow-up commit for stragglers.
AP151: Convention clash in single file (hbar with/without πi). Counter: grep file for all definitions after any formula.
AP152: "Ordered" on curve = LABELED (non-coinvariant), NOT totally ordered. Total order lives in R. Specify.
AP153: E_3 scope inflation. HDC needs B^Σ(A) as E_2-coalgebra. E_inf input → B^Σ exists → E_3. E_1 input (Yangian) → B^Σ doesn't exist → ordered bar gives only E_2.
AP154: Two E_3 structures: (a) algebraic (HDC on E_2 bar, no conformal vector); (b) topological (Sugawara, conformal vector at non-critical). Identification CONJECTURAL (conj:e3-identification). Topological: PROVED affine KM; CONJECTURAL general. Chain-level may fail class M.
AP156: wp_1 conventions: (a) θ_1'/θ_1 vs (b) Weierstrass ζ_τ=(a)+2η_1*z. DIFFERENT monodromy. Specify.
AP177: S_2=c/12 lambda-bracket confusion. Shadow invariant S_2=κ=c/2 for Vir. c/12=(c/2)/3! is divided-power coefficient. S_2 is convention-INDEPENDENT. Verify any S_r against Vol I census.
AP180: Cross-volume shadow bridge: S_2^{VolI}=κ=c/2=6*S_2^{lambda-bracket}. Most likely S_2^{Vol II} is WRONG.
AP181: pi_3(BU)=0 (Bott periodicity), NOT Z. pi_3(U)=Z confusion. Before pi_k(BX): fiber sequence pi_k(BX)=pi_{k-1}(X); check Bott parity.
AP182: κ_ch=χ(S)/2 applies to Tot(K_S→S). Conifold (Tot(O(-1)²→P^1)) is NOT local surface. Verify geometry.
AP183: McKay quiver of Z_3 = 3 oriented 3-cycles (9 directed arrows), NOT K_{3,3}. K_{a,b} undirected; McKay directed.
AP184: Excision=gluing (B_L⊗_A B_R=B(A), one copy over A); coproduct=splitting (B(A)→B(A)⊗B(A), two copies plain). Verify codomain.
AP185: pi_4(BU)=Z is obstruction GROUP, not guarantee E_2 exists. Nonzero homotopy = potential obstruction, not automatic structure.
AP186: Shallow correction without first-principles investigation. Before any correction, write down first-principles analysis; state correct theorem. Every wrong claim contains seed of correct theorem. Examples: "categorified averaging" wrong but factorization E_1→^Z E_2→^{Sym} E_inf real; "CoHA=bar" wrong but character coincidence reflects SV theorem CoHA=Y^+.

### Operadic-Structural
AP25/AP34/AP50 (four functors): B(A)=coalgebra. D_Ran(B(A))=B(A!)=algebra. Omega(B(A))=A (INVERSION, not Koszul duality). Z^der_ch(A)=bulk. D_Ran is VERDIER. Bulk is HOCHSCHILD. A^!_inf (Verdier, chain) != A^! (linear, strict). Compatibility IS Theorem A. NEVER "bar-cobar produces bulk."
AP65/AP81-AP85/AP88/AP103/AP104 (operadic): B_P(A)=P^!-coalgebra != BP=cooperad (different levels). Three coalgebras: Lie^c (Harrison/coLie), Sym^c (coshuffle, 2^n terms), T^c (deconcatenation, n+1 terms). Coshuffle != deconcatenation. Factorization coproduct (Sym^c on Ran) != deconcatenation (T^c on ordered); R-matrix descent relates. B_{Com}(A) is coLie, not cocommutative. P^i=cooperad != P^!=(P^i)^v=operad. Cotriple bar != operadic bar. E_1 is PRIMITIVE; modular/symmetric is av-image.
AP86/AP87/AP89-AP93 (SC/promotion): B_{SC}(A) for one-colour ill-formed; SC two-coloured, use promotion A→(A,A). Closed=B_{Com}(A), open=B_{Ass}(A) + mixed sector. SC mixed-sector dim=(k-1)!*C(k+m,m), NOT (k-1)!*m!. FM_n(X) connected; only strata factor. Curved d^2=κ*ω_g NOT coderivation. Two curvatures: μ_0 (algebra, g=0) vs d_fib^2=κ*ω_g (fiberwise, g>=1). δF_g^cross in CLOSED sector: "mixed channels" (propagator) != "mixed sector" (open-closed).
AP94-AP98/AP100/AP102 (shadow/Hochschild): ChirHoch*(Vir_c) degrees {0,1,2}. NEVER C[Θ]. ChirHoch != Gelfand-Fuchs (GF infinite-dim, ChirHoch bounded). Shadow algebra has graded Lie bracket NOT ring. av: g^{E_1}→g^mod LOSSY; at deg 2 recovers κ in abelian/scalar; for non-abelian KM gives κ_dp, full κ adds dim(g)/2. C0 fiber-center; C1 Lagrangian eigenspace; C2 shifted symplectic/BV conditional. Scalar κ+κ'=K from C1 + Theorem D not C2. Theorems specify which bar: B^ord, B^Σ, B^Lie.
AP99: K11 Lagrangian CONDITIONAL (perfectness + bar-cobar normal-complex).
AP101: "qi not merely iso on cohomology" tautological. Use "qi of A-inf-algebras" vs "chain qi."
AP127: Migrating `\input{chapter}`: add `\phantomsection\label{}` stubs for EVERY label and grep `\ref{}` pointing to them.
AP128: Engine-test synchronized to same wrong value. NEVER update test expectations from engine output. Derive INDEPENDENTLY. Most dangerous AP10 variant.
AP142: Local-global conflation on curves. "Koszul duality over point = over P^1" FALSE. Three errors: (a) HOMOTOPY RETRACT IS DATA (A^1→point comparison needs retraction/homotopy/transfer); (b) DISK != POINT (formal D thickening carries content; VAs live on D not point); (c) A^1 ALREADY HAS ARNOLD RELATIONS (Conf_n(A^1) has Arnol'd algebra; P^1 adds compactness). Consequences: (i) genus-0 chiral Koszul NOT "just" classical; (ii) "everything new at g>=1" overstated; (iii) fiber over point←D→A^1→P^1→X not studied; (iv) formal disk vs point still needs retraction. Counter: specify WHICH space (point/D/A^1/P^1/X), COMPARISON DATA (retraction/localization/thickening), on-the-nose vs extra-structure, WHAT content each step.
AP159: Four Yangian types: (1) Classical Y_ℏ(g): E_1-top on R. (2) dg-shifted Y^{dg}_ℏ(g): at point/formal disk. (3) Chiral Y(g)^{ch}: E_1-chiral on curve X, D-module. (4) Spectral: factorization on A^1_u. Conflating = type error. Use ^{ch}/^{dg} superscripts.
AP160: Three Hochschild, geometry determines which: (i) Top HH: E_1 → E_2 (Deligne); (ii) Chiral ChirHoch: E_inf-chiral → {0,1,2} (Thm H); (iii) Categorical HH: dg cat → E_2 w/ CY shifted Poisson. Bare "Hochschild" MUST carry qualifier. Conv:three-hochschild in concordance constitutional.
AP161: Five E_1-chiral notions HEALED 2026-04-17 to two on Koszul locus via thm:e1-chiral-notions-collapse (algebraic_foundations.tex): (A)/(B)/(C)/(E) mutually Quillen-equivalent; only (D) double-A_inf remains distinct (conj:double-ainfty-notion-D-relation). Registry: ordinary E_1-chiral input + speculative (∞,2)-enhancement. Outside Koszul locus, the original AP161 separation still applies.
AP162: E_3 requires (a) conformal vector, (b) non-critical level, (c) T(z) Q-exact. At k=-h^v: Sugawara undefined. PROVED affine KM; CONJECTURAL general. Cohomological; class M chain-level open.
AP163: "Lives on R × C" unjustified for E_1-chiral. End^{ch}_A on curve X; SC^{ch,top} bar is coalgebra over PRODUCT operad, NOT factorization algebra on product space. Chiral Deligne-Tamarkin needed.
AP164: Chiral Gerstenhaber != topological Gerstenhaber. Chiral: OPE residues on FM_k(C); topological: little 2-disks. Agree for E_inf via formality; diverge for E_1.
AP169: SC^{ch,top} is GENERIC case; E_3 is special. Most chiral algebras lack conformal vector (Heis, lattice, critical, E_1-chiral, free fields); for these SC^{ch,top} is FINAL answer. Treat as first-class object, multiple presentations.
AP170: Two Yangian definitions: def:e1-chiral-yangian (E_1-factorization + RTT, weaker) vs def:chiral-yangian-datum ((A, S(z), Δ^{ch}, {m_k}) + modular MC tower, stronger). Equivalence OPEN.
AP171: Associator dichotomy. Chiral QG equiv holds "up to Drinfeld associator Φ". COHOMOLOGICAL derived center = ASSOCIATOR-INDEPENDENT (sl_2 proved). COCHAIN-LEVEL = ASSOCIATOR-DEPENDENT. Bar-side (κ, shadow) associator-FREE; cobar/coproduct depend on Φ.
AP172: A^! is SC^!-algebra = (Lie Sklyanin, Ass), NOT SC-algebra. FORBIDDEN: "A^! is an SC-algebra."
AP173: Z^{der}_{ch}(Y(g)^{ch}) NOT computed. Predicted infinite-dim (C[qdet T(u)]). E_1 input: only E_2, not E_3.
AP174: Chiral QG equiv proved abstractly on Koszul locus. Concrete: sl_2 Yangian + affine KM only. Elliptic partial; toroidal absent.
AP175: Pentagon of equivalences is a STAR. Five SC presentations (operadic/Koszul/factorization/BV-BRST/convolution) route through operadic hub. Three direct links: Koszul-BV, BV-Convolution, Factorization-Convolution.
AP176: CONSTITUTIONAL — "arity" BANNED. "Degree" universal for bar grading, operadic input count, tree vertex valence, Stasheff identity level, SC mixed sector params, cooperad/operad component indices, endomorphism operad components, brace insertion count. NEVER reintroduce "arity." Grep `grep -rn '\barity\b' chapters/ appendices/ standalone/` must return ZERO hits.
AP179: Graph vertex valence context-dependent. Graph: "degree"; operad: "arity" (BANNED per AP176; use "valence"). For vertices carrying operadic ell_k, use "valence."

### Label / Prose / Scope
AP5: Grep ALL THREE volumes for variants (~/chiral-bar-cobar, ~/chiral-bar-cobar-vol2, ~/calabi-yau-quantum-groups) after EVERY correction.
AP12: When proving a claim, grep entire manuscript for variants. Update same commit.
AP40: Environment matches tag. Conjectured → `\begin{conjecture}`. ProvedElsewhere → theorem + Remark attribution.
AP105: Heisenberg = abelian KM at level k = abelian CS boundary. SAME OPE J(z)J(w)~k/(z-w)^2. Simple-pole needs ODD generator (symplectic fermion).
AP106: NEVER "This chapter constructs..." Open with PROBLEM (CG deficiency).
AP107: r^coll(z) differs from Laplace-transform r(z) for odd generators.
AP108: Heis = CG opening, NOT atom. Atom of E_1 = genuinely nonlocal (Yangian, EK quantum VA).
AP109: NEVER list results before proving.
AP110: Each volume's preface tells its OWN story; cross-volume in delineated subsections.
AP111: No "What this chapter proves" blocks. Restructure.
AP112: Never trust stale page counts. Verify against fresh builds.
AP113: See HZ-7.
AP114: Stub chapters (<50 lines, no theorems) → develop or comment out.
AP115: Architectural CLAUDE.md claims must be enacted in .tex. Metadata-source gap is most dangerous AP.
**Prose laws**: no AI slop; no hedging where math clear; no em dashes; no passive-voice hedging; every paragraph forces next; state/prove once; scope always explicit; prior work = one sentence per paper.
AAP1: Grep `antml` or `</invoke>` in .tex after every write.
AAP7: Grep current file before writing a formula that appears elsewhere in same file.

### 732-Agent Adversarial Campaign (AP186-AP233). Full catalogue: `compute/audit/new_antipatterns_wave12_campaign.md`.

AP186-210 (specialized): AP186 ProvedHere-no-proof-block; AP187 orphaned chapters; AP188 empty sections; AP189 dead labels; AP190 hidden imports (cited result doesn't prove claim, 119 findings); AP191 circular proof chains; AP192 scope inflation statement-vs-proof; AP193 biconditional, only forward proved; AP194 curved complex with flat tools (45); AP195 five-object conflation (47); AP196 SC misattribution non-formula; AP197 bare Hochschild (89); AP198 Whitehead lemma semisimple-only; AP199 strong filtration inequality direction; AP200 transfer theorem gap (H*(A) results applied to A); AP201 Baxter constraint not vacuous at λ=0; AP202 coderived element-wise invalid; AP203 class-M harmonic unproved; AP204 genus-0 boundary contradiction; AP205 reflexivity hidden in duality; AP206 object switch mid-proof (Verdier != cobar); AP207 center-side vs bar-side lift missing; AP208 Theorem A Verdier algebra/coalgebra flip; AP209 missing lemma cited; AP210 topologization chain-level vs cohomological.

AP211-224 (new): AP211 test file absent (219); AP212 TODO/FIXME unresolved; AP213 stub chapter false coverage; AP214 cross-volume bridge outdated; AP215 preface advertising stronger than proved; AP216 Koszul (vii) genus-0 scope; AP217 Koszul (viii) ChirHoch freeness overclaim; AP218 SC-formality restricted to metric families; AP219 depth-gap d_alg=2 wrong line; AP220 D^2=0 wrong geometric space; AP221 Gerstenhaber single insertion; AP222 Theorem H config-space collapse unjustified; AP223 Theorem H bar-coalgebra/Koszul-dual conflation; AP224 README scope inflation.

AP225-233 (deep structural): AP225 (CRITICAL) genus-universality gap — all-genera scalar factorization NOT proved; genus-1 unconditional; clutching-uniqueness needed; affects Theorem D; AP226 K_0-class vs scalar (use Chern character); AP227 ProvedHere forwarding ("By Theorem X" is ProvedElsewhere); AP228 anomaly-Koszul dependency inversion; AP229 SC-formality propagation debt (Vol III stale); AP230 genus-1 sufficient claimed all-genera; AP231 draft artifacts in theorem statements; AP232 duality clause overclaiming family scope; AP233 compact/completed comparison gap MC3.

AP234-235 (preface rectification 2026-04-17): AP234 two-Koszul-conductors-same-letter — κ(A)+κ(A^!) (scalar complementarity, family-dependent: 0/13/250/3/98/3) distinct from Trinity K(A)=c+c^!=-c_ghost(BRST) (-k/2dim(g)/26/100/196); related by κ+κ^!=ϱ_A·K with ϱ_N=H_N-1 for principal W_N, ϱ_A=0 for KM/free, ϱ_BP=1/6. The equation κ+κ^!=K (bare, no ϱ) is FALSE for every standard family. Canonical K-values in universal_conductor_K_platonic.tex:795-821 and higher_genus_complementarity.tex:3015-3120. Grep trigger: any `K(Vir)=13` or `K=250/3 for W_3`. Counter: cross-check K at self-dual c=13 — correct is K=26, not 13. Confusion pattern #218 (cache). AP235 quaternitomy/quadrichotomy drift — "quadrichotomy" is canonical (matches thm:quadrichotomy, chap:shadow-quadrichotomy-platonic); "quaternitomy" is invented hybrid, drifts across preface, working_notes, part-introductions. Grep `quaternitomy` after any write naming the G/L/C/M partition. Confusion pattern #219 (cache).

AP236 (bar_construction.tex rectification 2026-04-17): blacklist-slug leakage into typeset parenthetical. `/B\d+` identifiers from the Wrong-Formulas Blacklist leak into manuscript `\textup{(}...\textup{)}` annotations during agent edits (e.g. `\textup{(}/B49; treated in Chapter~X\textup{)}`). A reader sees a non-existent citation; the slug rots with every blacklist renumbering. Constitutional metadata hygiene violation. Grep `\\textup\{\(\}\s*\/B\d+|(\s*\/B\d+;` before every commit touching `.tex`. Healing: delete the slug; keep only the mathematical cross-reference to the target chapter. Companion symptom: orphaned `\textup{(}\textup{)}` empty parentheticals left behind when a macro-resolved reference is deleted between the delimiters — grep `\\textup\{\(\}\\textup\{\)\}`. Confusion pattern #221 (cache).

### 2026-04-17 Waves 1-4 Deep Adversarial Audit (AP237-AP244)

The four-wave adversarial audit of 2026-04-17 surfaced eight patterns beyond AP1-AP236. Themes: splitting-principle degree accounting failures, statement/proof internal contradictions, naming-after-physical-source without geometric content, closure-by-repackaging, advertised-but-not-inscribed characterizations, forward-reference lemmas, HZ-IV non-disjoint dependencies, and overcounted foundational terms. The programme carries approximately 30% terminological redundancy and systematic scope-label over-inflation; these APs are detection patterns for future waves.

**AP237 (Splitting-principle degree accounting).** The product $\prod_{i=1}^g (\alpha \cdot x_i) = \alpha^g \prod x_i$ is a SINGLE MONOMIAL in $\alpha$, NOT a candidate for "scalar-channel linear-in-$\alpha$ extraction." When applying the splitting principle to extract a scalar-channel coefficient, verify the coefficient depends linearly in the scalar parameter. The correct route is a $K$-theoretic class $\alpha \cdot \lambda_{-1}(\mathbb{E})$ linear in $\alpha$ by construction, NOT a product-formula manipulation. Found at Vol I `higher_genus_foundations.tex:5742-5786` in `prop:scalar-obstruction-hodge-euler` Step 3d: Theorem D's "all-genera CLOSED" chain relies on this step. Counter: if you write $\alpha^g \prod x_i$ and claim "linear-in-$\alpha$ projection gives $\alpha \cdot c_g(\mathbb{E})$", STOP; that projection does not exist. Return to the $K$-theoretic class at Step 1c.

**AP238 (Statement/proof internal contradiction).** Same mathematical object assigned different numerical values within a single proposition-proof pair. Example: Vol III `cy_d_kappa_stratification.tex:1143` proposition STATEMENT formerly wrote $\kappa_{BKM}(\Phi_1) = 0 + 0 = 0$ (N=1 coincidence) while proof at :1170-1171 wrote $\kappa_{BKM}(\Phi_1) = 10$ (Φ_10 weight 10); the engine `compute/lib/kappa_bkm_universal.py` recorded $5$. HEALED 2026-04-17: correct value is $\kappa_{BKM}(\Phi_1) = 5$ via Gritsenko $\Delta_5$ weight-5 paramodular form of level 1; statement/proof/engine now agree on 5; "N=1 coincidence" narrative retracted as confabulation. Counter: before inscribing a proposition, numerically evaluate every symbol in statement AND proof at a common test point; any discrepancy = retraction required. Constitutes an AP10 variant scoped within a single environment.

**AP239 (Naming-after-physical-source without geometric content).** The programme systematically names objects after physical sources (K3, Y(gl(4|20)), Kummer, Monster) without verifying the geometric input is actually used beyond lattice rank/signature. Examples: (i) Y(gl(4|20)) is Mukai signature (4,20) symmetric indefinite → OSP(4|20) candidate, NOT gl(4|20) super-Yangian (Wave-2 F19); (ii) K3 abelian Yangian `thm:k3-abelian-yangian-presentation` depends only on rank-24 signature-(4,20) even unimodular lattice + CY_2 constraint; no K3-specific geometry enters beyond Mukai pairing. Counter: for every named object (K3-X, Monster-Y, CY-Z), list what geometric input beyond rank+signature is USED in the theorem; if none, rename to "rank-N sig-(p,q) X" with a remark noting physical source inspiration.

**AP240 (Closure-by-repackaging).** A claimed "N prior gaps CLOSED (date)" where the inscribed closure repackages the gap into cited remarks rather than independently resolving it. Example: Vol I Theorem D claim "All 3 prior gaps CLOSED (2026-04-16)" per CLAUDE.md Theorem Status; source-level inspection at `higher_genus_foundations.tex:6172-6229` shows the proof chain (a)-(d) exhibits fiberwise curvature = κ·ω_g^Ar, but step (b) "no higher Hodge bundle enters scalar channel" is ASSERTED via `rem:propagator-weight-universality`, not derived. `AGENTS.md:566` self-declares the gap. Counter: "N gaps CLOSED" requires INDEPENDENT resolution, not relocation into cited remarks. Audit each claimed closure: does the proof use any load-bearing `\ref{}` that itself asserts rather than derives? If yes, the closure is a repackaging.

**AP241 (Advertised-but-not-inscribed characterization).** A preface, FRONTIER.md, or CLAUDE.md advertises a theorem component that is absent from the inscribed theorem. Example: "Tropical Koszulness" listed in CLAUDE.md Koszul 10+1+1 meta-theorem breakdown (`chapters/theory/chiral_koszul_pairs.tex` advertised equivalence); ZERO matches for "tropical" in the source file. Appears only in `preface_sections5_9_draft.tex:943` prose. Counter: every time CLAUDE.md lists an equivalence (i)-(xii), grep the target theorem's .tex home for the equivalence keyword; absence = the advertisement is heuristic, not theorem. Grep gate before preface edits: for each advertised clause, locate the inscription in the `.tex` home; no inscription = remove the advertisement or downgrade to heuristic remark.

**AP242 (Forward-reference lemma labelled as inscribed).** A load-bearing lemma cited in a proof as `\ref{lem:foo}` where either (a) the label `\label{lem:foo}` exists nowhere across the three volumes, or (b) the lemma body at the label is itself a forward-reference with no proof. Example: Vol I `thm:hochschild-concentration-E1` proof cites `lem:chiral-quadratic-koszul` as the chiral transport of Shelton-Yuzvinsky contracting homotopy; the lemma is load-bearing for Theorem H but the inscription is a forward reference, not proven at callsite. Counter: before inscribing a theorem, verify every `\ref{lem:...}` in the proof resolves to a `\label{lem:...}` with a proven body. Paired with HZ-11 for cross-volume resolution.

**AP243 (HZ-IV decorator non-disjoint dependency).** "Three disjoint independent-verification paths" claim where V1/V2/V3 share a hidden lemma dependency. Example: `thm:bfn-phi-ade-identification` HZ-IV decorators cite V1 Kronheimer ALE + V2 BFN 2016 + V3 BKR01 derived McKay; but V1 Kronheimer is INPUT to V3 BKR01 (Kronheimer-Nakajima hyperkähler construction is an input to McKay correspondence theorem). Counter: for every HZ-IV block, list literature citations per path; cross-check that path-A citations do not appear as lemmas or inputs in path-B or path-C; shared papers = non-disjoint. If shared, either (i) replace one path with a genuinely disjoint source, or (ii) downgrade the claim from "three disjoint paths" to "two disjoint paths" and register the dependency explicitly.

**AP244 (Overcounted foundational terms).** The programme has N+k named "distinct" foundational notions where only N are mathematically distinct. Example: 5 E_1-chiral notions (A strict ChirAss / B A_∞ in End^ch / C EK quantum VA / D A_∞ in E_1-chiral / E factorization on Ran^ord); on the Koszul locus where the programme's theorems apply, (A)=(B)|_{strict}, (E)=(B) D-module-theoretic, (C)↔(B) via Drinfeld associator. Only (B) and (D) are genuinely different, and (D) is labelled "open problem"; 5 names → 2 objects on the working locus. Also: 4 classes G/L/C/M are post hoc labels on a 2×2 Boolean grid ($C = 0$ vs $\neq 0$) × ($Q = 0$ vs $\neq 0$). Counter: for every "N foundational notions" claim, construct a concrete example (e.g., V_k(sl_2) at generic k) where ALL N notions yield distinct objects; if any two coincide on the programme's working locus, declare terminological collapse and either consolidate or annotate the redundancy explicitly. HEALED 2026-04-17 (E_1-chiral instance) — `thm:e1-chiral-notions-collapse` in `chapters/theory/algebraic_foundations.tex` proves on Koszul locus (A)/(B)/(C)/(E) collapse to one category via three Quillen equivalences (Stasheff truncation + Etingof-Kazhdan + Beilinson-Drinfeld D-module); (D) remains open (`conj:double-ainfty-notion-D-relation`). Operational registry is now two: ordinary E_1-chiral + speculative (∞,2)-enhancement. The G/L/C/M 2×2 Boolean grid instance of AP244 remains pending.

### 2026-04-17 Preventative Anti-Patterns (AP245-AP254)

The 2026-04-17 adversarial audit + HEAL wave surfaced preventable errors whose prospective check would have prevented the error. The ten anti-patterns below catalogue those prospective checks. AP245 (statement-proof-engine agreement) and AP238 (statement/proof contradiction) are paired: AP238 is the diagnostic pattern, AP245 is the forward-check. AP246/AP239, AP247/AP244, AP249/HZ-11 are similarly paired forward-check / retrospective-diagnostic dyads.

**AP245 (Statement-proof-engine numerical agreement).** Every proposition carrying a numerical value must be evaluated at a common test point across (a) statement body, (b) proof body, (c) compute engine; any mismatch forces retraction before inscription. Canonical violation: $\kappa_{BKM}(\Phi_1)$ at Vol III `cy_d_kappa_stratification.tex:1143/1170-1171` — statement gave 0, proof gave 10, engine `compute/lib/kappa_bkm_universal.py:396-401` gave 5; correct value is 5 via Gritsenko 1999 $\Delta_5$ (weight-5 paramodular form of level 1). Counter: before inscribing a proposition with a numerical value, grep for the proposition's claim across statement / proof / engine and reject on any inconsistency. Relation to AP238: AP238 is the retrospective diagnostic; AP245 is the prospective guard.

**AP246 (Signature type-assignment discipline).** A lattice signature $(p, q)$ determines the Lie-algebra type of any Yangian attached to it. Mukai signature $(4, 20)$ is symmetric INDEFINITE (orthogonal form), giving the orthogonal series $O(p, q)$ or the orthosymplectic super series $OSP(p|q)$ (Arnaudon-Crampé-Doikou-Frappat-Ragoucy 2003, arXiv:math/0304188), NOT the general linear GL(p|q) or super-linear $\mathfrak{gl}(m|n)$. Counter: before naming a Yangian from a lattice signature, verify (i) symmetric indefinite → $SO(p,q)$ Yangian or $OSP(p|q)$ super-Yangian via reflection equations; (ii) symplectic → C-type Yangian; (iii) Z/2-supergraded with parity swap → $GL(m|n)$ super-Yangian. $Y(\mathfrak{gl}(4|20))$ was a naming artifact (AP239); correct candidate is $Y_{osp}(4|20)$.

**AP247 (Functor terminology requires single target).** A functor $\Phi: \mathcal{C} \to \mathcal{D}$ requires a SINGLE target category $\mathcal{D}$. A $d$-indexed family $\{\Phi_d\}$ with $d$-dependent target $\mathcal{D}_d$ is a CORRESPONDENCE PROGRAMME, not a functor. Counter: before using "functor" terminology, verify (a) single unified target category; (b) morphism action defined; (c) composition verified on a concrete morphism pair (e.g., Mukai transform K3 → K3 for d=2). CY-to-chiral $\Phi$ failed all three; correct terminology is "CY-to-chiral correspondence programme $\{\Phi_d\}_{d \geq 1}$". Relation to AP244: this is the terminological-inflation diagnostic specialised to functor claims.

**AP248 (Coloured dioperads vs operads).** When $\Sigma_n$-equivariance fails between colours (e.g., closed-to-open but not open-to-closed), the correct term is DIOPERAD (Gan 2003, arXiv:math/0210098) or WHEELED PROPERAD (Merkulov-Vallette 2009, arXiv:0907.2895), NOT operad. Dunn additivity does not apply to coloured structures with directional restriction. $SC^{ch,top}$ was mislabeled as an operad; it is a two-coloured dioperad with directional colour-restriction. Counter: before writing "P is an operad" for any coloured structure, verify all $\Sigma_n$-equivariances; colour-restricted equivariance = dioperad.

**AP249 (Base-change / extension theorems require inscription, not citation).** CLAUDE.md theorem-status labels cannot promote a theorem from "on fixed curve X" to "over $\overline{\mathcal{M}}_{g,n}$ including boundary" without an inscribed modular-family theorem. Base-change via six-functor formalism (Francis-Gaitsgory GR17 Vol II) must be INSCRIBED or marked `\ClaimStatusProvedElsewhere` with explicit attribution, not silently cited. Theorem A was advertised "modular-family PROVED 2026-04-16", but the inscribed `thm:A-infinity-2` is on fixed curve only; the modular-family extension was cited but not inscribed. Now healed to "PROVED unconditional on fixed smooth curve X; modular-family extension CONDITIONAL on GR17". Counter: for every claim "theorem extends from $X$ to $\overline{\mathcal{M}}_{g,n}$", locate either a local inscription or an explicit `\ClaimStatusProvedElsewhere` attribution; absence = the extension is rhetorical.

**AP250 (Algorithm uniformity requires per-type verification).** When citing an algorithm (Frenkel-Mukhin q-characters, Nakajima quiver variety, Kashiwara crystal basis, Lusztig canonical basis) as uniform across type classification (A/B/C/D/exceptional), verify per-type validity. Non-simply-laced exceptionals ($G_2$, $F_4$) often fail uniformity per Hernandez 2006 (arXiv:math/0606381) + Nakajima 2001+. The claim "FM algorithm gives true q-character for all simple types" is FALSE for $G_2$, $F_4$; per-case verification required. Counter: before writing "uniformly in type", list the types and cite the per-type reference; absence = the uniformity is aspirational.

**AP251 (Attribution density floor).** A theorem marked `\ClaimStatusProvedHere` without $\geq 1$ load-bearing citation to a classical source is suspect of novelty-inflation. Programme-wide novelty audit (2026-04-17) found 37% classical recast without explicit attribution (Drinfeld 1985, Frenkel-Jing 1988, Kac-Peterson 1984, Gritsenko-Nikulin 1995, Borcherds 1998, DMVV 1997, Gannon 2016 absent from Vol III K3 Yangian theorem bodies). Counter: every `ProvedHere` theorem's proof must cite $\geq 1$ classical input; zero citations = either mark as "to the best of our knowledge, not previously in the literature" or inscribe the missing attribution (AP251 violation).

**AP252 (Chern character / Taylor expansion degree direction).** $\prod_{i=1}^g (1 - e^{x_i}) = \prod_{i=1}^g (-x_i + O(x_i^2))$ begins at degree $2g$ (top-degree monomial $(-1)^g \prod x_i$); higher-order corrections from $O(x_i^2)$ terms live in degree strictly GREATER than $2g$. "Lower degree" is wrong. Similarly $(1-t)^{-1} = \sum_n t^n$ expands UPWARD. Counter: before writing "lower-degree corrections" or "higher-degree corrections" in a Chern character / Taylor / splitting-principle expansion, substitute $x_i = 0.01$ numerically and verify the expansion direction. Vol I `higher_genus_foundations.tex:5505` had "lower degree" where "higher degree" was meant (healed 2026-04-17).

**AP253 (Inter-volume dependency graph for "one programme" claim).** Before claiming "one programme" across multiple volumes, verify cross-volume load-bearing citation count. Sparse dependency graph ($< 30$ cross-volume theorem citations) = loosely coupled trilogy, not unified theory. Counter: grep each volume for `\ref{V[0-9]-thm:}` or `\ref{V[0-9]-prop:}`; if Vol I contains zero references to Vol II or Vol III theorems, information flow is strictly forward and the "one programme" framing is rhetorical. Downgrade to "three-volume series" until the dependency count crosses threshold.

**AP254 (Closure-date commit-floor).** Any "closure wave" or "N-agent session" claim attributing theorem closures to a specific date must have commit density $\geq N/10$ on that date (conservative batching factor). "Wave 14 closure 2026-04-16" produced 3 commits programme-wide — implying retroactive-compression. Counter: narrative wave labels are project-management shorthand, not synchronized mathematical events; status-table "PROVED (date)" entries should reference the inscription-commit date of the theorem body, not the rhetorical wave label. Audit: `git log --since=DATE --until=DATE+1day --oneline | wc -l` must match claimed agent count / 10.

### Wave 5 (2026-04-17 attack-then-heal swarm, AP255-AP266)

Twelve preventative patterns surfaced in the Wave-5 attack+heal campaign on Theorems A/B/C/D/H, MC1-MC4, MC5, topologization, chiral QG, cross-volume identities, Koszul equivalences, and shadow-tower closed forms. These complement AP245-AP254 by registering failure modes that a preflight check would prevent.

**AP255 (Phantom file + phantomsection mask).** A chapter file is advertised in CLAUDE.md theorem-status and in `standalone/theorem_index.tex`; the file does not exist on disk; `chapters/frame/preface.tex` carries `\phantomsection\label{...}` stubs that let `\ref{}` resolve without compilation errors, hiding the absence. Canonical violation healed 2026-04-17: `chapters/theory/theorem_C_refinements_platonic.tex` advertised 5 load-bearing lemmas (T1 derived-centre Koszul equivalence, T2 perfectness standard landscape, T4 +3 shift resolution, T6, T9 PTVV alternative) with preface `\phantomsection` masks at `preface.tex:5066`; Wave-5 inscribed the real 502-line chapter, retired the stubs, and downgraded T1 (Conjecture), conditionalised T2 (unconditional for Heis + affine KM non-critical; class-M boundary-stratum Conjecture), retracted T4 (non-issue — no Verdier pairing at g=0), folded T6 into C0, and tagged T9 mixed (clauses (i)-(ii) PTVV-genuine, clause (iii) forwards to main proof). Detection: for every `\phantomsection\label{foo}` in preface/intro, grep the whole repo for `\begin{[a-z]+}\s*\\label{foo}`; absence = phantom. Distinct from AP241 (keyword-absent-in-existing-file): AP255 is the file itself being absent. Counter: preflight grep before any status-table "PROVED" update — require (a) target file exists, (b) label is inscribed in a real environment, (c) build renders the `\ref{}` with a section number not `[?]`.

**AP256 (Aspirational-heal status drift).** Status table advertises a healing ("Feigin-Frenkel replaces JKL26", "no longer tautological", "nodal-sewing theorem inscribed", "PTVV genuinely disjoint", "strict chain-level via explicit $\eta_1^{(i)}, \eta_1^{(ii)}$") before the `.tex` has been edited to enact the replacement. Advertisement propagates into preface and standalones while the chapter body still cites the original dependency. Healed instances 2026-04-17: (a) `thm:glN-chiral-qg` Feigin-Frenkel "replaces JKL26" retracted — JKL26 restored at $N \geq 3$ with sharpened obstruction `prop:ff-screening-coproduct-obstruction` (the commutator $[Q_{\alpha_i}, \Delta_z^{Heis}]$ is a non-exact chiral 1-cocycle of class $(\Psi-1)/\Psi$ matching miura-cross-universality); (b) class-L topologization $[Q, \widetilde{G}_1] = T_{Sug}$ "strict chain-level" retracted to cohomological; (c) DS intertwining "no longer tautological" rewritten as "degree 1 tautological, degree $\geq 2$ via sl_3 RTT" to match engine docstring. Counter: for every status-table row claiming a mechanism replacement, grep the chapter body for the new mechanism keyword; zero hits = drift. Related AP249.

**AP257 (Engine-docstring vs manuscript contradiction).** A `compute/lib/*.py` engine's docstring or output admits a limitation that directly contradicts a manuscript status claim. Engines are ground truth for engine behaviour; if manuscript contradicts engine, the manuscript is wrong by default. Canonical instance: `compute/lib/ds_coproduct_intertwining_engine.py:103-104` admits "intertwining is a TAUTOLOGY at degree 1" while Vol I status advertised "no longer tautological." Counter: after any status-table "verified by engine X" update, read the engine's top-of-file docstring and reconcile. Healing: rewrite the status to match engine scope.

**AP258 (Cohomological-vs-chain status drift).** Status table writes "strict chain-level" while the inscribed theorem body says "in $H^\bullet(Q_{tot})$" or "up to $Q$-exact." Example: Vol II `thm:iterated-sugawara-construction` line 232 explicitly says "in $H^\bullet(\mathcal{A}^{BV}_{3d}, Q_{tot})$", contradicting the CLAUDE.md status-table claim that class-L original-complex chain-level is proved via explicit $\eta_1^{(i)}, \eta_1^{(ii)}$. Detection: grep status table for "strict chain-level" / "on-the-nose" / "chain-level"; for each hit, locate the cited theorem and verify the statement uses genuine chain-level language (explicit homotopy, no cohomological quotient). Healing: downgrade status to cohomological OR inscribe the explicit strict-chain construction with sign-tracked cochains.

**AP259 (Tautological-by-definition masquerading as theorem).** A theorem claims a nontrivial identity $\mathrm{obs}_g = \kappa \cdot \lambda_g$ whose proof defines the LHS via projection onto the RHS channel; equality is a definition, not a deduction. Violation at Vol I `higher_genus_foundations.tex:5447-5466`: K-theoretic step defines $\mathrm{obs}_g := \mathrm{ch}_g([\bar{B}]^{vir})_{\text{scalar channel}}$ then claims $\mathrm{obs}_g = \kappa \cdot \lambda_g$. Counter: for any identity $A = B$ labelled theorem, trace the definition chains of $A$ and $B$; if $A$ is defined via projection using $B$ as reference, the equality is tautological and the physical statement requires independent characterisation of $A$. Healing: rename to honest form "$\mathrm{obs}_g^{K\text{-thy}} = \kappa \cdot \lambda_g$ by construction; bar-curvature $\mathrm{obs}_g$ = K-theoretic $\mathrm{obs}_g$ separate theorem (conditional $g \geq 3$ on $\lambda_g$-conjecture)."

**AP260 (Scalar-channel linearity without constructed projector).** "Scalar-channel linearity" or "scalar-diagonal extraction" invoked without a constructed projector $\pi_\varepsilon$ commuting with pushforward and pullback. The splitting principle gives $\kappa^g \prod x_i$ on the flag bundle; "extract the linear-in-$\kappa$ component" requires a cohomological projection not a polynomial manipulation (see AP252 + AP237). Counter: grep for "scalar channel", "scalar diagonal", "linear-in-$\kappa$"; for each hit, verify a projector is constructed OR the scope is restricted to rank-1 families (Heisenberg, Virasoro) where scalar-diagonality is automatic. Healing: construct the projector, restrict scope with explicit hypothesis, or downgrade to Conjecture.

**AP261 (Mittag-Leffler single-index vs required bigrading).** A pro-system Mittag-Leffler argument claims stabilisation in cohomological degree $m$ alone while the filtration data supplies stabilisation only at fixed total weight $w$. Single-index ML may fail where bigraded $(m, w)$ succeeds. Healed 2026-04-17 at Vol I `mc5_class_m_chain_level_platonic.tex:229-437` Step 3: rewritten to stabilise in $(m, w)$-bigrading using weight-preservation of the bar differential plus Prop. `standard-strong-filtration`(iv) at fixed $w$. Detection: every pro-qi or pro-acyclicity claim must name the index over which transitions stabilise; if weight-dependent filtration is used, ML must be bigraded. Healing: rewrite ML step with explicit $(m, w)$ bigrading; if a case has $m$ fixed but $w$ unbounded, flag as genuine gap.

**AP262 (Hub-and-spoke TFAE scope inflation).** A $k$-fold TFAE is advertised; a `\begin{remark}[Scope]` immediately below retracts one spoke ("false off class G"); the headline stays $k$-fold. Healed 2026-04-17: FTM "seven-fold TFAE" with SC-formality spoke class-G-only demoted to six-fold Koszul-locus TFAE + separate class-G SC-formality corollary. Counter: after any TFAE statement, read the next 20 lines for scope remarks; if any spoke is retracted, rewrite as $(k-1)$-fold main + separate corollary for the scope-restricted spoke. Related: AP252 Chern-degree-direction; AP237 splitting-principle accounting.

**AP263 (Hopf-vs-bialgebra naming drift).** A theorem advertised as "quantum group equivalence" or "Hopf algebra structure" that constructs only a bialgebra because the antipode does not lift. The antipode proof is a PROVED NEGATIVE in Vol I (`S(T(u)) = T(u)^{-1}` does not lift to a vertex-algebraic antipode; two obstructions OPE and Hopf axiom); the cluster is a bialgebra, not Hopf. Healed 2026-04-17: `thm:chiral-qg-equiv` renamed to "Chiral bialgebra equivalence on the Koszul locus" with companion `rem:chiral-bialgebra-not-hopf`; cross-volume prose propagation of the rename remains partial (≈16 sites). Counter: before labelling "Hopf" or "quantum group", grep the programme for "antipode does not lift" / "S(T(u)) does not extend" scoped to the same object; if the negative is proved, demote to "bialgebra" / "chiral bialgebra equivalence."

**AP264 (Phantom `\ref{rem:foo}` to non-existent remark).** A load-bearing prose reference to `\ref{rem:...}` where the remark does not exist anywhere in the three volumes. Distinct from AP242 (forward-reference lemma): remarks have no proof body, so the reference signals missing MOTIVATION, making the prose sentence vacuous. Companion to AP255 at the remark level. Wave-5 instance: `rem:kac-moody-filtered-comparison` in `ftm_seven_fold_tfae_platonic.tex` Spoke 4 V_k(sl_2) non-tautology witness — zero matches anywhere. Counter: for every `\ref{rem:...}` in a load-bearing sentence, grep all three volumes for matching `\label{rem:...}`; zero hits = phantom. Also AP265 (bibkey variant): `\cite{foo}` with `foo` absent from bibliography renders `[?]` at build; Wave-3 A^{∞,2} audit found six missing bibkeys (`HackneyRobertson2017`, `HackneyRobertson2019`, `Francis2012`, `GR17`, `Positselski2018`, `Positselski2011`, `Hinich2003`). Counter: every `\cite{}` grep-checked against `standalone/references.bib` pre-commit.

**AP265 (Primary-source number-theoretic labelling).** A prime labelled "Kummer-irregular" or "leading Bernoulli irregular" without primary-source numerator verification at the cited $B_{2m}$. Primary sources: Akiyama-Tanigawa recurrence, Buhler-Harvey 2011 irregular-prime tables. Detection: for every "Kummer-irregular" assertion, cite a specific $B_{2m}$ with numerator divisibility and verify against a primary-source table. Healing: re-label as "Riccati-arithmetic characteristic primes" (preserving computational role without number-theoretic claim) or retract. Healed 2026-04-17 retractions (all verified REGULAR through $B_{418}$): {23, 43, 61, 193, 419, 1423, 3067, 811, 2111, 16657}. Companion discipline: $\{691, 3617\}$ genuinely Kummer-irregular at $B_{12}, B_{16}$ leading on $Z_g$; higher Kummer-irregular primes appear elsewhere (e.g., 2111 in $N_9$), so "absence from $S_r$" claims must scope explicitly to leading two primes.

**AP266 (Sharpened-obstruction healing register).** Not a failure mode but a POSITIVE healing pattern. When a deep-math attempt to prove a conditional theorem genuinely FAILS, the correct heal is to EXHIBIT THE OBSTRUCTION as an explicit cohomology class and identify its coefficient with a known universal invariant, converting the failure into a sharpened frontier item with a falsification test. Template: `prop:obstruction-to-X` stating $[R] \in H^\bullet_{?}(?, ?)$ with coefficient $c_X$; followed by `rem:falsification-test` stating any alternative construction must reproduce $c_X$ at a named sector. Exemplar: `prop:ff-screening-coproduct-obstruction` class $[R_i(z)] \in H^1_{ch}(V_{Heis}, V_{Heis} \otimes V_{Heis}[z, z^{-1}])$ with coefficient $(\Psi-1)/\Psi$ matching `thm:miura-cross-universality`. The pattern is the Beilinson dictum in action: a smaller true theorem (the obstruction is explicit) replacing a larger false one ("descent holds unconditionally"). Companion discipline: when a deep-math attempt succeeds, inscribe; when it fails, sharpen; never retreat silently.

### Wave 6 (2026-04-17 attack-then-heal, AP267-AP269)

**AP267 (Recursion lower-bound silently overwrites initial datum).** A recursion stated $a_r = f(a_{j}, a_{k})$ over $j + k = r + 2$, $j \le k$ with a lower bound on $(j, k)$ that admits the theorem's supposed INITIAL DATUM as a recursion output. Canonical violation: Riccati algebraicity theorem at `higher_genus_modular_koszul.tex:17990` quantified recursion as `j + k = r + 2, 3 \le j \le k` with no $r \ge 5$ restriction; at $r = 4$ the only solution is $(j, k) = (3, 3)$ (diagonal), and plugging in forces $S_4 = -9\alpha^2 / (16\kappa) = -9/(2c)$ for Virasoro — contradicting the Zamolodchikov-pinned $S_4 = 10/[c(5c + 22)]$. The theorem STATEMENT correctly names $(\kappa, \alpha, S_4)$ as the three initial data; the PROOF silently admitted $S_4$ as a recursion output. Counter: after writing any recursion with a lower bound, enumerate $(j, k)$ solutions at the lowest non-initial $r$ and verify none collapses onto an initial-datum slot. Healing: tighten lower bound to $r \ge r_0$ where $r_0$ = first genuinely non-initial stage; add strict $j < k$ or explicit diagonal-term handling at even $r$. Healed 2026-04-17 at Riccati proof: $r \ge 5$, $j < k$, with even-$r \ge 6$ diagonal term separated.

**AP268 (Meromorphic status asserted for a polynomial-constant rational function).** A Koszul-conductor or central-charge identity $K(k) = c(k) + c(k^!)$ claimed "meromorphic with pole at $k = k_0$" when the actual identity is polynomial-constant in $\mathbb{Q}(k)$ (genuine pole is removable because both summands have equal-and-opposite residues). Canonical violation: Vol II `bp_chain_level_strict_platonic.tex:111-125` and CLAUDE.md BP row claimed $K^{FL}_{BP}(k) = -12(k+3) - 48/(k+3)$ "meromorphic with pole at $k = -3$"; direct computation with Fateev-Lukyanov $c^{FL}(k) = -(2k+3)(3k+1)/(k+3)$ and dual $k \mapsto -k - 6$ gives $K^{FL}_{BP}(k) = [50(k+3)]/(k+3) \equiv 50 \in \mathbb{Q}(k)$. Both Arakawa ($K = 196$) and Fateev-Lukyanov ($K = 50$) conventions give polynomial-constant conductors; they differ by central-charge normalisation, not by meromorphy class. Counter: before labelling a Koszul-conductor identity "meromorphic", factor numerator; if it divides the denominator, the identity is polynomial-constant. Healing: state the constant value and identify the convention; the Koszul-conductor identity is always rationally CONSTANT on the standard landscape (this is the defining feature of Koszul self-duality on the complementarity line). Healed 2026-04-17 across Vol I `standalone/bp_self_duality.tex`, Vol I `standalone/N4_mc4_completion.tex:1010`, Vol II `bp_chain_level_strict_platonic.tex:116-133, 152`, Vol II `FRONTIER.md` V2-NF3, and Vol I `compute/tests/test_motivic_shadow_full_class_m.py:50-56`. CLAUDE.md theorem-status row for BP still needs the stale FL-meromorphic claim removed.

**AP269 (SDR-formula fabrication with proved-contradictory witness).** CLAUDE.md status-table asserts an UNCONDITIONAL proof of a strong completion theorem via a named explicit homotopy with closed-form $h_{htpy} = f(\text{parameters})$, citing a concrete SDR mechanism (e.g. "Wakimoto one-step SDR"); the homotopy formula has ZERO manuscript witness AND the programme contains a separately inscribed negative proposition proving the cited SDR mechanism FAILS at the structural level required. Canonical violation: "MC4$^0$ PROVED UNCONDITIONAL via Wakimoto one-step SDR; homotopy $h_{htpy} = (1 - \iota p)/(L_0 - h - N + 1)$ invertible generically; Feigin-Frenkel screening (W_N)". Grep Vol I manuscript: zero hits for $h_{htpy}$, zero for "Wakimoto one-step SDR", zero for the resonance-locus analysis $L_0 - h - N + 1 \in \mathbb{Z}_{\le 0}$. Worse: Vol I `ordered_associative_chiral_kd.tex:10148` proves Feigin-Frenkel screening is NOT chiral-coproduct-compatible on the Heisenberg parent — directly refuting the W_N route. Similarly "Costello-Gwilliam factorization algebra $F^{CS}_k$ on $\mathbb{R} \times C$, single-colored, not SC^{ch,top}$" in the Topologization row has zero typeset witness in Vol II `e_infinity_topologization.tex` (only a bibliography comment mentions CFG). Counter: before promoting a status-table entry from Conjectured to UNCONDITIONAL, grep the manuscript for every named object in the mechanism (homotopy formula, screening operator, construction name); zero hits across all three volumes ≥ fabrication. Check also for separately proved NEGATIVE propositions that contradict the promotion. Healing: retract UNCONDITIONAL → Conjectured with explicit scope (named resonance locus, explicit failure mode), or inscribe the SDR with chain-level computation in a single witness family (sl_2 at k=1 for Wakimoto; single Heisenberg generator for Feigin-Frenkel). This is the meta-pattern behind AP255/AP256/AP258 and stricter than each: it catches mechanism-level fabrication that neither phantom-file nor aspirational-heal detectors surface, because the theorem "exists" in some form while its MECHANISM is hollow.

**AP270 (Multi-object status-row conflation).** A single theorem-status row advertises one "VERIFIED / PROVED" banner whose numeric signature (dimension, invariant count, level count, test count) does not match any single inscribed object cleanly, because the banner is compressing two or three distinct inscribed results + one conjecture under one heading. Canonical violation: CLAUDE.md "Drinfeld center Heis VERIFIED: 5 invariants at 6 levels. Naive dim 1 vs derived dim 3. 72 tests" simultaneously describes (a) PROVED Heis computation (dim 3 graded, $P(t) = 1 + t + t^2$, `prop:derived-center-explicit`), (b) PROVED affine $\widehat{\mathfrak{sl}_2}$ computation (total dim 5, parametrises "5 invariants"), (c) CONJECTURAL preface equivalence `conj:v1-drinfeld-center-equals-bulk`. The "72 tests" is a stale count; current test count exceeds. Counter: status-row object-count must match exactly one inscribed theorem; compound rows rewritten as multi-row table or enumeration. Healing: split row into explicit per-family lines with individual labels and status tags.

**AP271 (Reverse drift: manuscript ahead of CLAUDE.md).** The standard drift direction is CLAUDE.md advertising ahead of the manuscript (AP256 aspirational-heal). The REVERSE drift is equally dangerous and less detected: the manuscript inscribes an HONEST scope retraction in a `\begin{remark}` or `\ClaimStatusConditional` tag, and CLAUDE.md carries a stale "PROVED unconditional" status that contradicts the inscribed remark. This violates the Beilinson epistemic hierarchy directly (source > CLAUDE.md > memory). Canonical violation: Vol I `chapters/theory/theorem_A_infinity_2.tex:879-915` `rem:A-infinity-2-modular-family-scope` states verbatim "Francis-Gaitsgory six-functor base-change on the relative Ran prestack... and Mok25 logarithmic factorization-gluing at the boundary (nodal degenerations)... neither is inscribed as a theorem in Vol~I; both are invoked only at the citation level", while CLAUDE.md Theorem A row advertises "modular-family PROVED... (1) hypothesis (c) base change proved via BD holonomic + GR17... (2) nodal sewing theorem via Mok25 log FM + Francis-Gaitsgory factorization-gluing + HS-sewing convergence (chain-level at nodes)". Similar reverse drift on the Topologization row (Vol II `e_infinity_topologization.tex:401-411` retracts the $\eta_1^{(i)}, \eta_1^{(ii)}$ strict-chain-level claim as frontier-only, but CLAUDE.md still carries the claim). Counter: before any `git commit` touching CLAUDE.md status-table, grep the cited theorem's chapter for `\begin{remark}[Scope` / `\ClaimStatusConditional` / `is cited but not inscribed`; any hit = stale status, requires rectification of CLAUDE.md to match manuscript. Healing: when manuscript retracts, CLAUDE.md retracts same day; refuse to propagate stale "PROVED" labels. Stronger than AP256: AP256 catches CLAUDE.md running ahead; AP271 catches CLAUDE.md lagging behind.

**AP272 (Unstated-cross-lemma via folklore citation).** A theorem whose proof uses a load-bearing unstated lemma cited only by author-year attribution ("by Polyakov-Wiegmann 1984", "by Feigin-Frenkel", "by GR17 Chapter III"), with the cited source not actually inscribing the lemma in the form used, or inscribing it only for a restricted scope that the theorem applies outside. Canonical violation: Trinity $K$ theorem `thm:uc-trinity` (Vol I `universal_conductor_K_platonic.tex`) uses $c_{matter}(A) = -c_{matter}(A^!)$ as an unstated cross-lemma attributed to "Polyakov-Wiegmann 1984"; this is the Feigin-Fuks involution in the Virasoro case and case-by-case for other families, NOT a general theorem. Consequence: the $K_c = K_g$ arrow in `thm:uc-trinity` is a family-by-family numerical coincidence, not a uniform proof; Heisenberg breaks the coincidence (K_g = -2, K_c = 2, AP234 $\kappa$-route gives 2k). Counter: after writing any theorem proof citing `\cite{AuthorYear}` as a single-citation mechanism, grep the bibliography entry for a primary-source statement of the cited lemma; if absent or scope-restricted, inscribe the lemma internally (with scope) before promoting the parent theorem to `\ClaimStatusProvedHere`. Stronger than AP28 (undefined terminological qualifier) and AP190 (cited result does not prove the claim): AP272 specifically catches the "uniform proof via a lemma that only exists family-by-family" pattern behind many of the programme's overclaimed universal theorems (Trinity K, Feigin-Frenkel universal screening, Koszul conductor universal on the standard landscape).

**AP273 (Admitted-redundant item counted in independent-equivalence tally).** A meta-theorem advertises "N unconditional equivalences"; the chapter-internal scope remark admits one of the N is a consequence of another listed equivalence, not an independent statement; the headline count stays N and propagates to preface, standalones, surveys, and Vol II bridge prose. Canonical violation: Koszul 10+1+1 meta-theorem `chiral_koszul_pairs.tex:2323-2974` with internal remark at lines 2335-2338 admitting item (vi) Barr-Beck-Lurie "is a consequence of item (v), not an independent equivalence; we list it because it is the ∞-categorical formulation most useful for applications"; the chapter opening (line 89) honestly writes "9+1+1+1", but the theorem header (line 2348) restores "10 unconditional equivalences" and this "10" propagates to `landscape_census.tex:1269`, `working_notes.tex:790,4177`, four surveys, and Vol II `thqg_introduction_supplement_body.tex:898,1816` — nine files carrying the inflated count. Counter: after writing any "N equivalences" meta-theorem, grep the same file for `consequence of item` / `not an independent` / `∞-categorical restatement`; any hit admitting redundancy requires the headline count to drop by the number of admissions. Healing: either (a) retract the count to N-k "genuine equivalences + k restatements" across ALL propagation sites atomically (AP5), or (b) promote the admitted-consequence item to a genuinely independent criterion (e.g. for Barr-Beck-Lurie, inscribe a monadicity criterion NOT routed through the counit qi). Distinct from AP262 (TFAE scope inflation via retracted spoke): AP273 is about counting items at the HEADLINE level where internal honesty exists but does not propagate.

**AP274 (Rhetorical functor identification across disjoint objects).** A prose sentence stitches two mathematically disjoint functors (domains, codomains, or constructions differ) by the phrase "is the categorified analogue of" / "is the same as" / "is the derived version of", with no theorem inscribed relating the two. The slogan then load-bears a multi-arrow diagram (operadic circle, holographic dictionary, bar-center correspondence) whose closure depends on the identification. Canonical violation: Vol I `standalone/en_chiral_operadic_circle.tex:2121` writes "This [Drinfeld-centre functor $Z(-): E_1\text{-Cat} \to E_2\text{-Cat}$] is the categorified analogue of the averaging map" — identifying the $\Sigma_n$-coinvariant projection $\mathrm{av}: \fg^{E_1} \to \fg^{\mathrm{mod}}$ (chain-level, LOSSY, kills the R-matrix) with the Drinfeld centre functor $Z(\cC)$ (category-level, PRESERVES R-matrix via half-braiding). These are COMPLEMENTARY, not identified: av is the symmetric-bar projection on chain complexes; $Z$ is the derived commutant for monoidal categories. The programme's E_n-circle Arrow 3 (E_1 → E_2 Drinfeld) rests on this slogan; once disambiguated, the honest statement is "the R-matrix of an E_1-algebra determines the half-braiding of the Drinfeld centre of its module category" (`thm:e-arrow3-braiding` clause (ii)), a narrower and PROVED claim. Counter: for every "is the X analogue of Y" prose sentence in an arrow/diagram proof, verify (a) X and Y have the same source and target categories, (b) there is a labelled isomorphism or commutative-diagram theorem relating them, (c) the proof body uses the identification in a specific computation. Missing any of (a)(b)(c) ⇒ retract the identification, rewrite as the narrower inscribed theorem, downgrade the arrow it load-bears. Related: AP186 (shallow correction pattern) catches "categorified averaging" as a known wrong identification; AP244 (overcounted foundational terms) catches the case where the two functors were being advertised as distinct but secretly collapsed on the working locus; AP274 is the DUAL case where two disjoint functors are rhetorically FUSED.

**AP275 (CLAUDE.md narrative inverts the correct mathematical discipline).** Stronger than AP271 (reverse drift = stale label after retraction): AP275 is CLAUDE.md prose that INVERTS a known-correct convention, such that a reader who EDITS the manuscript to match the CLAUDE.md line would INTRODUCE a bug. Canonical violation: Vol I CLAUDE.md Elliptic chiral QG row reads "θ_1'/θ_1 propagator (NOT Weierstrass ζ)"; the correctly inscribed propagator in `yangians_drinfeld_kohno.tex:7318`, `lattice_foundations.tex:4695`, `heisenberg_eisenstein.tex:481` is the FULL Weierstrass $\zeta_\tau(z) = \theta_1'(z|\tau)/\theta_1(z|\tau) + 2\eta_1 z$; dropping the quasi-period correction $2\eta_1 z$ (as the CLAUDE.md line would instruct) BREAKS classical Yang-Baxter (documented as FM30). Companion violation in the same row: "Felder dynamical R(z,τ)" — inscribed is Belavin non-dynamical r-matrix for sl_2 (`prop:elliptic-rmatrix-shadow`:7251), explicitly admitted in `rem` at :7364. Counter: before writing a CLAUDE.md convention line "X not Y", verify the inscribed expression uses X; if it uses the named NOT-Y object, the CLAUDE.md line is inverted. Healing: rewrite the status row to match source ("propagator IS the full Weierstrass $\zeta_\tau = \theta_1'/\theta_1 + 2\eta_1 z$; dropping the $2\eta_1 z$ correction breaks CYBE; see FM30 for the bug pattern").

**AP276 (Two-point extrapolation with undisambiguated competitor).** A closed form $f(s) = g_1(s)$ is announced based on agreement at two data points; a structurally different $f(s) = g_2(s)$ agrees at the SAME two points; neither is privileged without a third data point or a structural derivation. Publishing $g_1$ as the closed form invites AP-155 novelty overclaiming. Canonical violation: "Conjectural spin-stratified lattice $C_{\cW^{(s)}} = s(s+1)$ for each primary in principal $\cW_N$" advertised from $C_{\cT} = 2 \cdot 3 = 6$ (s=2, Virasoro T-line) and $C_{\cW_3} = 3 \cdot 4 = 12$ (s=3, W_3 W-line); the same two data points are also fit by $g_2(s) = 4s$ giving $C(2) = 8 \ne 6$ — FAILS at s=2 — so this particular pair is actually DISAMBIGUATED at s=2. But the honest residual: `shadow_tower_higher_coefficients.tex:2237-2283` lists FOUR candidate forms (s(s+1), 6(s-1), 4s²-4s+6, 4s), all marked Conjectural, all needing a third data point or structural derivation to distinguish. Counter: any "generalize to all $s$/$N$" conjecture based on $\le 2$ data points must (a) list at least two plausible competitors in the same inscription remark, (b) name the falsification test (e.g. compute $W_4$ W_4-line coefficient and check against each candidate), (c) refuse to privilege one form in the preface or headline prose. Healing: inscribe all competing forms as `\begin{conjecture}` alternatives in a single `\begin{remark}[Competing extrapolations]`; flag the load-bearing falsification computation as a programme priority. Related: AP155 (novelty overclaiming via architectural framework recovering known invariants); AP237 (splitting-principle degree accounting failure from two-point coincidence).

**AP277 (Vacuous HZ-IV test body behind sound decorator prose).** An HZ-IV independent-verification decorator carries a sound `verified_against` / `disjoint_rationale` prose field identifying three genuinely disjoint verification paths (distinct literature, distinct machinery, distinct working-locus families), but the test function body is a hard-coded tautology — `assert True` / `return True` / three Booleans all `True` with no wired computation. The decorator passes CI, the programme-wide HZ-IV coverage metric ticks up, the mathematics is neither verified nor disconfirmed by the decorator. Distinct from AP128 (engine-test synchronized to the same wrong value): there, two objects share a wrong mental model; here, one decorator simulates verification without executing any. Canonical violation: Vol I `compute/tests/test_theorem_H_hochschild_koszul.py:46-93` decorates Theorem H verification with Feigin-Fuchs 1984 (Vir) + Wang 1998 (W_N) + Whitehead+Kunneth (affine sl_2) — three genuine disjoint paths — while the test body contains three hard-coded `True` booleans with no engine invocation. The mathematics of Theorem H is independently sound (the proof body at `chiral_hochschild_koszul.tex:1490-1558` is load-bearing-complete), so the HZ-IV gap is a CODE-LAYER tautology, not a mathematical fragility. But it is the most detection-resistant variant of the tautology-registry pattern: prose audit sees a clean decorator; mathematics audit sees a proved theorem; only CODE audit catches the vacuous body. Counter: before accepting any HZ-IV decoration, `grep` the test body for an actual `assert` on a numerical value or a call into `compute/lib/*` — absence of either is a fabricated decorator. Healing: wire the test to existing engines (for Theorem H: `koszul_hilbert.py`, `bar_cohomology_koszul_criterion.py`) and check a per-family Hilbert-series value (e.g. Heis dim 3, V_k(sl_2) total dim 5); the decorator then earns its `verified_against` claim.

**AP278 (Moduli-space boundary classification asserted without construction).** A theorem invoking the "Costello TCFT dictionary" / "Segal factorization geometry" / "Stasheff associahedra boundary" / similar moduli-space formalism to VANISH or CANCEL a bracket-level obstruction, with the load-bearing step claiming "the only boundary components of $\overline{\mathcal{M}}_X$ are A and B" — without constructing $\overline{\mathcal{M}}_X$ as an explicit 1-dimensional compact manifold, without enumerating the strata of the Deligne-Mumford-type compactification, without citing a primary source (Costello 2005 Theorem, specific Segal paper, Stasheff monograph chapter) that GIVES this boundary count. The obstruction's vanishing follows from the asserted boundary classification via the standard "$\partial^2 = 0$ in cobordism" argument; the argument is sound ONCE the classification is given, but the classification itself is the load-bearing step. Canonical violation: Vol III `m3_b2_saga.tex:570-617` `thm:total-ainf-compat` proof asserts "the only boundary components are $b \circ B^{(2)}$ and $B^{(2)} \circ b$" with "no other Connes-hierarchy operations" to justify "$\mathrm{Obs}_{A_\infty} = 0$ universally via Costello TCFT"; no moduli-space construction, no citation to a specific Costello paper+theorem that gives this boundary enumeration, no verification that $B^{(0)}, B^{(1)}, B^{(3)}$ strata are absent from the Kuranishi-type closure. Companion violation: `prop:cyclic-ainf-framing-compat`. Related: AP257 engine exposes the violation — `chain_level_m2_b2_cancellation.py` verifies naive $\{b_2, B^{(2)}\}$ and $\{b_3, B^{(2)}\}$ land in DIFFERENT bar-arity graded pieces and cannot cancel directly, contradicting the theorem-level universal vanishing claim. Counter: for every theorem invoking a moduli-space boundary argument, verify (a) $\overline{\mathcal{M}}$ is constructed as a specific space (dimensions, strata, orientation, corners), (b) the boundary classification cites a primary source at a specific theorem/proposition number, (c) the closure of the lower-dimensional strata is accounted for. Absence of any ⇒ downgrade the theorem to Conditional on Segal/Costello strictification, OR retract to a restricted scope (here: toric/formal CY_3 only, where formality trivialises the A_∞-coherence). Stronger than AP272 (unstated cross-lemma via folklore citation) because the "citation" here is to a DICTIONARY or FORMALISM rather than to a specific statement.

**AP279 (Rename heal island-local, cross-volume consumer carries semantic drift in VALUE not just NAME).** An object is renamed in its canonical home (e.g. $Y(\fgl(4|20)) \to Y_{\osp}(4|20)$ per AP246 signature-type-assignment discipline) with a local `\begin{remark}[Correction]` retracting the old name. The rename is marked "healed"; the cross-volume CONSUMERS continue to use the old name. AP149 resolution-propagation-failure catches this far, and AP246 is the underlying discipline. The stronger AP279 pattern: the consumer's use of the old name is not a naming artefact only — the consumer has computed a formula or identity whose VALUE is convention-dependent on the old type, and under the rename the formula's value must be re-derived, not just relabelled. Canonical violation: Vol III `k3_yangian_chapter.tex` renames to $Y_{\osp}(4|20)$ via `rem:gl-to-osp-correction`; Vol II `unified_chiral_quantum_group.tex:795-796` writes the complementarity identity $\kappa^{super} + \kappa^{super,!} = \max(4, 20) = 20$ inherited from the $\fgl$ super-Yangian convention; under the osp rename, the identity couples to $\osp(m|2r)$ Shapovalov duality rather than $\fgl(m|n)$ Berezinian duality — the VALUE of the complementarity sum may no longer be 20, must be re-derived from the osp reflection-equation pairing (Arnaudon-Crampé-Doikou-Frappat-Ragoucy 2003, arXiv:math/0304188). Likewise Vol II `super_chiral_yangian.tex:597`, `dnp_identification_master.tex`, `grt_parametrized_seven_faces.tex` all carry the old name in load-bearing identity prose. Counter: after any AP246-style rename, audit cross-volume consumers for (a) the name, (b) any identity whose value depends on the algebra type. Renaming (a) without re-deriving (b) leaves the manuscript carrying a retracted identity. Healing: atomic rename + per-site identity re-derivation in the same commit (AP5 + AP149); if the new identity value is open, downgrade the consumer's identity to `\begin{conjecture}` pending re-derivation. Stronger than AP149 (resolution-propagation failure) and AP246 (signature type-assignment) because the consumer carries a semantic-drifted IDENTITY, not just a type-drifted NAME; the "heal" is only complete when the value is re-verified or honestly re-labelled Conjectural.

**AP280 (Three-step epistemic inflation: remark → standalone → headline).** A mathematical observation starts life as a one-or-two-sentence `\begin{remark}` noting a coincidence (generating-function identity, numerical match, character factorisation). The remark is promoted to a `\begin{proposition}` in a standalone paper with `\ClaimStatusProvedHere` and a proof body that is itself one or two sentences silently citing a classical theorem (MNOP, TUY, Beauville-Laszlo, Maulik-Okounkov). The proposition's status-table row in CLAUDE.md is written as "IDENTIFIED" or "PROVED" at the programme-headline level. Readers see "IDENTIFIED" in the theorem-status table and assume a genuine bridge theorem; the underlying content is a remark-level coincidence resting on an unstated citation. Canonical violation: Vol I Shadow = GW(C³) chain — `rem:shadow-gw-c3` (ordered_associative_chiral_kd.tex:9938, Remark) → `prop:v3-qg-shadow-gw` (standalone cy_quantum_groups_6d_hcs.tex:801, no status tag) + `thm:bar-macmahon` (higher_genus_modular_koszul.tex:27521, ProvedHere with two-sentence proof silently citing MNOP 2006 §4.3 "second quantization raises exponent to n") + `prop:conifold-dt-gv` (proof body literally "Bryan-Pandharipande, MNOP.") → CLAUDE.md "Shadow = GW(C³) IDENTIFIED". Counter: audit every CLAUDE.md theorem-status row ending in "IDENTIFIED" / "PROVED" / "DISCOVERED" for three-step inflation: trace backward to the load-bearing inscription; if the inscription is a remark, or a proposition whose proof body is a single-line folklore citation, the headline is inflated. Healing: downgrade each step — the remark stays a remark; the standalone proposition retags to `\ClaimStatusProvedElsewhere` with Remark[Attribution]; the CLAUDE.md headline rewrites to "NUMERICAL COINCIDENCE" or "STRUCTURAL IDENTIFICATION CONJECTURAL (bridge to [MNOP/TUY/MO] cited, not independently derived)". Stronger than AP60 (tag only genuinely new content ProvedHere) and AP240 (closure-by-repackaging) because AP280 tracks the full three-step propagation chain and names the characteristic failure signature: remark-length original insight, proof body ≤ 2 sentences, silent citation to a deep classical theorem as the load-bearing step. Distinguishes epistemically inflated headlines from genuine bridge theorems.

**AP281 (Bibkey naming-drift at scale — catastrophic phantom-citation rate).** Vol I at audit time (2026-04-17): 86 bibkey definitions in `standalone/references.bib`, 693 unique `\cite{...}` keys used across `chapters/` + `standalone/` + `appendices/`, yielding 621 phantom citations — approximately 87% of all `\cite{}` invocations render as `[?]` at build. The root cause is not missing sources but NAMING DRIFT: the same reference is cited under multiple aliases (`CostelloGaiotto` vs `CostelloGaiotto2020`, `Positselski11` vs `Positselski2011`, `Kon03` vs `Kontsevich03`, `FBZ04` vs `FrenkelBenZvi`, `FF` vs `FF92` vs `Feigin-Frenkel` vs `FeiginFrenkel94`), with only one alias defined per reference and drafts citing the other. Vol II 212 phantoms (also alias-dominated). Vol III 26 phantoms (cleaner). Most-cited individual phantoms: `AF15` (53×, Ayala-Francis 2015), `Kon03` (39×, Kontsevich 2003), `FBZ04` (38×, Frenkel-Ben-Zvi 2004), `Feigin-Frenkel` (34×), `Fay73` (34×), `Zhu96` (31×), `Arakawa17` (29×), `Weibel94` (27×), `KL93` (25×), `GR17` (17×). Programme-internal self-cites (`LorgatVolI`, `VolI`, `Vol2-fractional-ghost-platonic`, `Vol1-class-L-antighost`) compound the gap (AP241 + AP253). Stronger than AP264 (single `\cite{foo}` with `foo` absent from bibliography): AP281 is the SYSTEMIC failure pattern where the bibliography is not a phantom per-item but an alias-drifted OVERLAY with vast coverage gaps. Counter: maintain a single canonical-alias table per reference; pre-commit gate `grep -rn '\\\\cite\\{[^}]*\\}' chapters/ standalone/ | awk -F'{' '{print $2}' | sort -u` against `grep -oE '@[a-zA-Z]+\\{[^,]+,' references.bib | sed 's/@[a-zA-Z]*{//;s/,$//' | sort -u` and refuse commits where the diff `used ∖ defined` grows. Healing: (a) fastest is an ALIAS LAYER — duplicate `\bibitem{...}` entries pointing to the same body under every used alias; (b) cleaner is a repo-wide `sed` canonicalisation pass standardising to one alias per reference. Top-10 heal priorities: `AF15`, `Kon03`, `FBZ04`, `Feigin-Frenkel`, `Fay73`, `Zhu96`, `Arakawa17`, `Weibel94`, `KL93`, `GR17`. The 87% phantom rate is a build-gate level defect: a reader of the PDF cannot follow citation chains; every downstream "PROVED" claim with AP272 (unstated cross-lemma via folklore citation) inherits the citation-level gap. The three programme-interior self-cites (`Lorgat26I`, `LorgatMKD1`, `LorgatVolI`) are distinct: they name unpublished Vol I/II/III material under external-reference form, conflating the author's own programme with the cited literature; they must be rewritten as internal `\ref{...}` cross-volume pointers, not external `\cite{...}`.

**AP282 (Status-table "VERIFIED" vs test-file `@pytest.mark.xfail` contradiction).** A CLAUDE.md theorem-status row advertises "VERIFIED at machine precision" (or "tests pass" variant); the corresponding test in `compute/tests/*` is marked `@pytest.mark.xfail(reason="Frontier: ...")` or `@pytest.mark.skip`. The status row and the test file disagree on whether the claim is verified. Companion to AP258 (cohomological vs chain-level) specialised to engine-level verification: AP282 is the test-status vs CLAUDE.md-status drift. Canonical violation: CLAUDE.md line 629 "KZB flatness VERIFIED" while `test_kzb_flatness_2pt` at `compute/tests/test_theorem_genus1_seven_faces_engine.py:433` carries `@pytest.mark.xfail(reason="Frontier: verify_kzb_flatness_2pt uses naive d_tau(zeta)=wp' identity; correct KZB flatness involves Halphen/Ramanujan system for Eisenstein series")`. Only the trivial companion `test_kzb_flatness_2pt_commutator_vanishes` (checking `[Ω,Ω]=0`) passes. Counter: before any CLAUDE.md "VERIFIED" status, `grep -rn '@pytest.mark.xfail\\|@pytest.mark.skip' compute/tests/test_<name>.py`; any hit requires the status row to read "PARTIAL" or "FRONTIER" with scope qualifier naming the specific xfail. Healing: rewrite status row to match the test file's actual scope (what passes, what xfails, what reduces to trivial). Distinct from AP128 (engine-test synchronized on wrong value): AP282 is engine-test correctly flagging failure, CLAUDE.md lying about it.

**AP283 (Formula confabulation in status-table RHS not inscribed in .tex or engine).** A CLAUDE.md theorem-status row contains a specific mathematical identity or closed form whose verbatim text does not appear anywhere in the manuscript (`chapters/`, `standalone/`, `appendices/`) or in any engine source (`compute/lib/*.py`). The formula reads plausibly to a casual reader but is a memory-drift artefact — the programme's actual inscribed identity differs in form. Canonical violation: CLAUDE.md line 629 wrote "d_τ(℘_1) = (1/(4πi)) d_w(℘ + ℘²)" for the Bernard heat identity; `grep -rn` for any of "d_w(wp+wp²)", "d_w(℘+℘²)", "d_w(p+p²)", "∂_w(℘+℘²)" across all three volumes returns zero hits. The correct inscribed identity is ∂_τ ζ = ℘' (Halphen/Ramanujan consequence at genus 1). Counter: every closed-form identity in a CLAUDE.md status row must be `grep`-verifiable against the `.tex` or engine source that advertises it as inscribed; verbatim `grep` with ≥0 hits is the commit-gate. Healing: quote the identity from source (file:line citation in-row), or retract the row. Distinct from AP269 (SDR-formula fabrication contradicted by proved negative) — AP283 is simpler: the formula is MEMORY-DRIFT, not fabrication with a contradicting witness. Companion to AP280 (three-step epistemic inflation) where the status headline overclaims a remark-level observation; AP283 is the sibling where the status headline carries a specific formula that does not exist in the programme.

**AP284 (Chapter dependency cycle invalidating the "foundation-before-application" ordering).** A CLAUDE.md Regression Safeguard (RS-15: "Koszul programme before higher_genus in dependency DAG") advertises a foundation-before-application chapter ordering. The nominal `main.tex \input{}` sequence places foundations early and applications late (Koszul at idx 18-25; higher_genus at idx 29-32). The `\ref{}` graph between chapters, however, has back-edges: foundation chapters cite application chapters. Canonical violation (Vol I, 2026-04-17 chapter-DAG audit): the chapter `\ref` graph contains a 99-chapter strongly connected component including the entire Part I + II + III backbone. `chiral_koszul_pairs.tex` (idx 18) carries 27 forward `\ref`s into `higher_genus_modular_koszul.tex` (idx 32); `koszulness_vii_multiweight_platonic.tex` carries 23. A reader following the main.tex compile order cannot follow the proof chain linearly — every Koszul chapter depends on higher-genus results that have not yet been read. 630 unresolved `\ref`s across the monograph are candidate AP242 forward-reference-lemma violations. Counter: before publishing or pre-release, build the chapter-level dependency graph via `grep -oE '\\\\ref\\{[^}]+\\}'` + `grep -oE '\\\\label\\{[^}]+\\}'` + `main.tex \\\\input{}` order; compute SCCs (Tarjan). An SCC of size > 1 indicates a cycle; an SCC spanning the backbone indicates systemic coupling that no `\input` reordering can fix. Healing: (a) for each back-edge, either inscribe the cited lemma locally in the citing chapter (breaking the back-edge), or restructure the ref as `\ClaimStatusProvedElsewhere` with forward attribution; (b) for the 630 unresolved `\ref`s, triage against the AP242 discipline. Distinct from AP147 (circular proof routing at theorem level): AP284 is the same phenomenon at chapter level — more pervasive, harder to see without graph analysis, and harder to fix because it requires discipline across dozens of chapters.

**AP285 (Alias section-number drift).** A theorem's alias carries a source-section number that does not match the section actually cited in the proof. Canonical violation: `thm:chiral-positselski-5-3` alias on `thm:chiral-co-contra-correspondence` at `bar_cobar_adjunction_inversion.tex:1331`; proof body cites Positselski 2011 §5.2 (co-contra correspondence) at lines 1350, 1439, 1445 — "5-3" in the alias is a transcription error propagating to `standalone/theorem_index.tex`, `standalone/garland_lepowsky.tex`, `notes/first_principles_cache_comprehensive.md`. Healed 2026-04-17 by renaming alias to `-5-2`. Counter: every alias with an embedded section-number must match `\cite{...}[§X.Y]` in the proof body. Distinct from AP264 (phantom `\cite{foo}` with `foo` absent from bibliography): AP285 is a SECTION-level drift inside a bibkey that resolves. Stronger than AP272 (unstated cross-lemma via folklore): the citation is specific — just to the wrong section. Programme-wide audit: every alias pattern `<thm>-<source>-<section>` should be verified against the proof body once per audit.

**AP286 (Tactical phantomsection alias vs semantic heal).** When facing phantom-ref consumer sites, two heal options exist: (A, tactical) add `\phantomsection\label{foo}` aliasing to the nearest canonical theorem — resolves `\ref{}` without fixing underlying inscription gap; (B, semantic) inscribe the theorem body locally or retarget the consumer refs to the canonical name. Option A is legitimate for naming-drift repairs; it is INSUFFICIENT when consumer prose treats the aliased theorem as a DISTINCT result from the canonical. Canonical violation: `thm:topologization-tower` phantomsection alias at `preface.tex:5084` aliasing to Vol II `thm:iterated-sugawara-construction` + `thm:e-infinity-topologization-ladder` — Vol I prose at `part_iv_platonic_introduction.tex` treats "topologization tower" as umbrella name for BOTH Vol II theorems, while the alias points to only one. Agent self-flagged its own heal as "tactical close, not semantic one". Counter: before accepting Option A alias, verify consumer prose does not treat aliased theorem as distinct or umbrella; if distinct/umbrella, escalate to Option B (inscribe or retarget). Distinct from AP255 (phantom file + phantomsection mask): AP286 catches the legitimate-but-insufficient alias-heal, where phantomsection resolves `\ref{}` but prose semantics require more.

### Wave 12 (2026-04-17/18 attack-then-heal swarm, AP287-AP292)

Six patterns surfaced in Waves 8-12 (HZ-IV coverage campaigns, CY-D stratification, systemic V1-* phantoms, Vol III κ-subscript discipline, super-complementarity drift, numerical-arithmetic precedence).

**AP287 (Structural-impossibility primitive tautology — HZ-IV-W8-B).** A decorated test at "degree 1" / "primitive" / "initial-data" level where the three declared-disjoint verification paths are ALL structural boolean predicates returning `True` by definitional construction (e.g. `_pred(bool, bool) -> bool` returning `True` as long as each input is `True`, with inputs hard-coded `True`). The HZ-IV decorator payload is genuinely disjoint at the bibliographic level, but Python-level numerical cross-checking is STRUCTURALLY IMPOSSIBLE because the claim is primitive-by-fiat. Distinct from AP277 (vacuous HZ-IV test body behind sound decorator prose): AP277 flags tautological bodies where genuine numerical cross-check is possible; AP287 flags tautological bodies where no numerical observable exists at the primitive level. Wave-10 scan found 17 Vol II files / ~94 claims matching this pattern — the entire `test_climax_theorems_wave*_iv.py` series (waves 3-18) plus `test_universal_holography_functor_fm_iv.py`. Three healings menu: (1) scope-restrict — rename test, move decoration to degree ≥ 2 where non-primitive content exists (DS intertwining healed this way, Wave-8 #51); (2) honest omission with module-level `# HZ-IV-W8-B FLAG` comment, NOT counted toward genuine coverage (Wave-10 #70 applied to all 17 files); (3) downgrade to Conjectured if the claim is structural-only. Healed 2026-04-17: Wave-10 #70 scan + flag; Wave-11 #75 numerical upgrade where possible (8 tests across 4 files upgraded with genuinely disjoint mathematical provenance; 14 files retain flag because primitive-by-construction).

**AP288 (Decorator-label vs test-body computation disjointness — HZ-IV-W8-C).** `assert_sources_disjoint` validates label disjointness in the `derived_from` / `verified_against` fields, NOT computational-path disjointness in the test body. A decorator can pass `assert_sources_disjoint` (three labels differ) while the test body computes all three paths via the identical RHS expression with different variable names. Canonical violation: `test_periodic_cdg_admissible.py::test_chiral_steenrod_rank_matches_uq_and_ckl` declared exterior / $u_q$ / KL paths but computed all three as `2**rank` inline. Wave-9 #61 lint found 6 additional Vol II violations (all healed); Wave-8 #51 healed the seed case via three genuinely distinct arithmetic routes (power-set / iterated product / Dynkin-node exponentiation). Distinct from AP277 (vacuous body) and AP287 (structural impossibility): AP288 is the subtle case where three non-trivial computations are claimed but secretly collapse to one. Counter: after HZ-IV decoration, AST-walk the test body to detect syntactic-RHS identity; `assert_sources_disjoint` should raise on label-disjoint-but-computation-identical at import time. Infrastructure recommendation: extend `compute/lib/independent_verification.py`'s assertion with body-disjointness detection.

**AP289 (Künneth-multiplicative vs additive for supertrace/Hodge-characteristic invariants).** The Hodge supertrace $\Xi(X) = \sum_q (-1)^q h^{0,q}(X)$ (identified with $\kappa_{ch}(A_X)$ for compact CY$_d$ in Vol III) is Künneth-**multiplicative**: $\Xi(X \times Y) = \Xi(X) \cdot \Xi(Y)$. An additive rule $\Xi(X \times Y) = \Xi(X) + \Xi(Y)$ would make $\Xi$ behave like Euler characteristic under disjoint union, NOT under product. Canonical violation: Vol III `cy_d_kappa_stratification.tex` preamble and proof body asserted "Vol I additivity $\kappa_{ch}(X \times Y) = \kappa_{ch}(X) + \kappa_{ch}(Y)$" at lines 153-155, 463-466, 659. Catastrophe was partially masked because every listed product in the chapter had at least one factor with $\Xi = 0$, so additive and multiplicative agreed coincidentally; K3×K3 would have exposed the drift (additive = 4, multiplicative = 4 — also coincidence — but for any $X$ with $\Xi(X) \neq 0, 1$ and $Y$ with $\Xi(Y) \neq 0, 1$ the two rules diverge). Similarly preamble claimed $\kappa_{ch}(K3 \times E) = 3$ via additive (2+1=3) contradicting body supertrace (K3: 2, E: 0, product: 0). Counter: any super-trace, elliptic-genus, or Hodge-characteristic invariant under product must use Künneth-multiplicative; additive behavior is restricted to Euler characteristic under disjoint union. Healed 2026-04-17 (Wave-4 Vol III CY-D agent): 13 surgical edits to `cy_d_kappa_stratification.tex`. Distinct from AP237 (splitting-principle degree accounting): AP237 is about the Chern-character / Taylor-expansion direction; AP289 is specifically about multiplicative-vs-additive behavior under Cartesian product for super-trace characteristic classes.

**AP290 (HZ-7 $\kappa$-subscript type-swap: $\kappa_{cat}$ formula assigned to $\kappa_{ch}$ or vice versa).** HZ-7 mandates Vol III $\kappa$ always subscripted: $\kappa_{ch}, \kappa_{cat}, \kappa_{BKM}, \kappa_{fiber}$. The weaker sub-pattern is TYPE-SWAP: a formula correctly belonging to one subscript is assigned to another. Canonical violation: Vol III `quantum_group_reps.tex:498-576` propositions `prop:kappa-cat-quantum-groups` and `prop:kappa-cat-complementarity` used the Sugawara-shift formula $\dim(\mathfrak{g})(k+h^\vee)/(2 h^\vee)$ (correctly $\kappa_{ch}$ for chiral algebras) and assigned it to $\kappa_{cat}$. The same chapter at line 549 correctly wrote $\kappa_{cat}(K3 \times E) = \chi(\mathcal{O}) = 0$; the propositions contradicted the chapter's own correct formula. $\kappa_{cat}$ is the categorical Euler characteristic $\chi(\mathcal{O}_X)$ of the CY category; $\kappa_{ch}$ is the chiral central-charge conductor. These are mathematically distinct invariants even when both are numerically zero. Counter: before inscribing any Vol III proposition with $\kappa$-subscript, verify the formula matches the subscript's definitional meaning — Sugawara / central-charge formulas are $\kappa_{ch}$; $\chi(\mathcal{O}_X)$ / $h^{p,q}$ supertrace formulas are $\kappa_{cat}$ (or $\kappa_{ch}$ via the Hodge supertrace identification theorem for compact CY$_d$); $c_N(0)/2$ Borcherds-weight formulas are $\kappa_{BKM}$. Healed Wave-11 #77: propositions renamed to $\kappa_{ch}$ with caveat contrasting $\kappa_{cat}$ in parallel. Stronger than HZ-7 bare-$\kappa$ prohibition: HZ-7 catches zero-subscript; AP290 catches wrong-subscript.

**AP291 (Systemic phantom-ref via self-disabled label — `% label removed:` pattern).** A chapter accumulates `\ref{SLUG}` references to labels defined within the SAME chapter (self-references), then a rectification pass comments out the defining `\label{...}` via `% label removed: SLUG` without removing or retargeting the `\ref{SLUG}` consumers. The build treats unresolved `\ref{SLUG}` as cross-volume externals per `linear_read_notes.md` convention (Vol II-specific), silently tolerating 222+ phantom refs. Distinct from AP264 (phantom `\ref{rem:}` to non-existent remark): AP264 is about refs to labels that never existed; AP291 is about refs to labels that DID exist and were intentionally disabled. Canonical violation: Wave-8 #54 Vol II Part VI 3d gravity audit discovered 222 `\ref{V1-*thqg*}` references across 12 Vol II `chapters/connections/` files pointing at self-labels disabled via `% label removed:` comments. Classification (Wave-9 #58 sweep): Class A (self-disabled) 141 slugs, 219 refs; Class B (mis-prefixed active) 1 slug, 2 refs; Class C (genuinely external) 0; Class D (Vol III) 0; Class E (absent everywhere) 1 slug (Wave-10 #66 inscribed 8 missing anchors including genuine-absent). Counter: pre-commit grep `% label removed:` across all `.tex`; for each removal, grep the same file + all downstream files for `\ref{SLUG}`; zero hits required before accepting the removal. Healing: either re-inscribe the `\label{SLUG}` uncommented, or retarget every `\ref{SLUG}` to a live label, or delete the ref with parenthetical prose. Wave-9 #58: 248 heals across 13 files (222 V1-prefix strips + 141 label re-inscriptions + 1 retarget + 25 opportunistic secondary).

**AP292 (Operator-precedence arithmetic bug in Python cross-check computation).** A Python test performs a numerical cross-check of a mathematical identity via explicit arithmetic with `Fraction` or similar; left-associative operator precedence collapses the intended value to a different value, causing the test to either pass wrongly or fail with a confusing mismatch. Canonical violation: Wave-11 #75 W8-B numerical upgrade initially wrote `Fraction(24,12)*c/Fraction(2)` intending $24/(12 \cdot 2) \cdot c = 1 \cdot c$; Python's left-to-right evaluation gave $(24/12) \cdot (c/2) = 2 \cdot c / 2 = c$, then when $c$ was substituted numerically the path-value landed on 1 coincidentally for $c = 1$ but mismatched for general $c$. Symptom: test passed at one c-value and failed at others. Counter: for every `Fraction` / rational-arithmetic expression in a test, parenthesize to DISAMBIGUATE; prefer `Fraction(24, 12 * 2) * c` or `Fraction(1, 1) * c` pre-computed; alternatively evaluate symbolically via `sympy` and compare to the programme's intended closed form. Healed Wave-11 #75 by replacing the Polchinski path with an Eguchi-Ooguri-path rewrite using explicit `Fraction(1, 2) * c`. Distinct from AP245 (statement-proof-engine numerical agreement): AP245 is about three levels agreeing on a value; AP292 is about a single level computing the wrong value due to language-level precedence. Programme-wide impact: this class of bug is detectable by testing at MULTIPLE parameter values (AP186 cache-pattern variant — coincidental agreement masks precedence bug).

**AP287 (Cross-volume theorem existence without HZ-11 attribution)**. Vol I theorem cites Vol II/III label that EXISTS in target volume, but Vol I citation lacks HZ-11 attribution discipline (no `\ClaimStatusConditional`, no `\begin{remark}[Attribution]`). Label resolution works at build time — reader sees section number — but Vol I theorem silently inherits external scope. Distinct from AP271 (reverse drift, CLAUDE.md vs manuscript) and HZ-11 (label-only phantom): AP287 is the subtle middle case — label IS inscribed cross-volume and resolution works, but attribution discipline still requires Conditional tag + attribution remark. Canonical violation: Vol I `chiral_moonshine_unified.tex:409` cites Vol II `thm:uhf-monster-orbifold-bv-anomaly-vanishes` in clause (ii) of `thm:v-natural-e3-topological` (itself ProvedHere) with 3-step proof sketch that does not inscribe Vol II theorem locally. Label resolves; scope inheritance does not. Counter: every Vol I ProvedHere theorem citing Vol II/III label must (a) inscribe cited theorem locally with proof body, or (b) tag clause `\ClaimStatusConditional` + `\begin{remark}[Attribution]`. Healing: audit each Vol I ProvedHere; for every `\ref{V[2-3]-...}` or cross-vol label, verify Conditional tag OR attribution remark present. Related: AP242 (forward-reference lemma absent), AP249 (base-change cited not inscribed), HZ-11 (full phantom). AP287 is the middle case where everything resolves technically but scope discipline is violated.

**AP288 (Session-ledger stale narrative).** Adversarial/attack+heal session notes (`adversarial_swarm_YYYYMMDD/*.md`, `notes/*_attack_heal.md`, `notes/session_YYYYMMDD_master_synthesis.md`, `notes/part_N_platonic_reconstitution.md`) are written at the over-confident closure moment of the session. Subsequent waves retract or rescope some of those closures; the session ledger is NOT re-rectified to match the retraction. A naive reader (or future LLM agent) treating the ledger as ground truth inherits the retracted narrative. Canonical violations this session: Vol II `notes/bershadsky_polyakov_universal_holography_attack_heal.md` L63/69/154/173/202 claim "strict chain-level E_3-topological via (η_1^(i) + η_1^(ii)) absorption" 5× — RETRACTED 2026-04-17 as frontier item per `rem:frontier-class-L-strict-chain-level` at `e_infinity_topologization.tex:382-411`. Vol II `notes/session_20260417_master_synthesis.md:L48-68, :172` claims W_∞ E_∞ endpoint UNCONDITIONAL + class-M chain-level UNCONDITIONAL without cohomological/weight-completed qualifier. Vol II `notes/part_VI_climax_platonic_reconstitution.md:13` treats `thm:uch-gravity-chain-level` as proved. Vol II `notes/universal_celestial_argyres_douglas_attack_heal.md:47, 230` claims chain-level thm:uch-main applies unconditionally. Counter: every session ledger requires a dated "RETRACTED" or "RESCOPED" section appended after each subsequent wave's heals, OR a top-of-file "HISTORICAL — superseded by FRONTIER.md §X" banner matching the Vol III FRONTIER.md preservation convention (line 122-123). Healing: append dated RETRACTED annotations at each stale-claim site; OR top-of-file "WARNING: written under pre-heal narrative, superseded by <source>" banner; OR reverse-chron diff log at end of file. Distinct from AP271 (reverse drift CLAUDE.md vs manuscript): AP288 is notes-vs-manuscript drift where notes carry an older version of the programme's understanding. Stronger than AP254 (closure-date commit-floor) because AP288 catches lexical echoes of retracted narratives surviving in session notes that were never re-rectified, not just suspect closure-wave rhetoric. Programme-wide discipline: when a wave N retracts a wave N-K closure, add a dated annotation ("RETRACTED at wave N / YYYY-MM-DD per <source>") at every hit in the wave N-K ledger.

### Wave 13 (2026-04-18 recovery-infrastructure audit, AP293-AP295)

Three patterns surfaced in the 2026-04-18 session when a user-invoked "relaunch all rate-limited agents" revealed that `scripts/resume_failed.py`, advertised under "Recovery infrastructure" in CLAUDE.md and routinely invoked to pick up "rate-limited / timed-out" agents, was silently mass-failing 98 agents across 7 campaigns because the `codex` CLI binary was not installed on the host. Each per-agent failure produced a 78-byte file `# <agent_id> — ERROR: [Errno 2] No such file or directory: 'codex'`. The size-based "empty" classifier in `resume_failed.py` (threshold 200 bytes) then re-picked those files as fresh failures on the next invocation, and `scripts/campaign_dashboard.py` reported the permafailing campaigns as "running". Cache entry: Pattern 222 (rate-limited vs binary-not-found).

**AP293 (Recovery-infrastructure prerequisite guard absent).** Scripts advertised under §"Recovery infrastructure" of CLAUDE.md (`scripts/resume_failed.py`, `scripts/campaign_dashboard.py`, the 9 campaign scripts) assume external CLIs (`codex`, `pdflatex`, etc.) are on PATH without verifying at startup. When the CLI is missing, the scripts enter per-agent failure mode — every agent writes a tiny error file, the top-level exit code is 0, the dashboard reports "running", and the user experiences a SILENT MASS-FAILURE. The pattern is indistinguishable from genuine rate-limiting or genuine empties without inspecting the per-file error message. Counter: every script that shells out to an external CLI must `shutil.which("codex")` at startup; on miss, print a loud "PREREQUISITE MISSING: codex CLI not on PATH. Install via <URL> before invoking this script." and exit with a non-zero code. Canonical violation: `scripts/resume_failed.py:89` unconditionally calls `subprocess.run(["codex", "exec", ...])`. Healing: add a `check_prerequisites()` call at the top of `main()` raising `SystemExit(2)` with installation instructions when `codex` is absent; mirror this in `scripts/campaign_dashboard.py` and every script in `CAMPAIGNS` (the dict at `resume_failed.py:21-29`). Related: AP80 (engine without test file — verify BOTH exist after agent completion) applied at infrastructure layer rather than theorem layer.

**AP294 (File-size threshold conflates three distinct failure modes).** `resume_failed.py` classifies per-agent output by the two-branch rule at `resume_failed.py:74-78`: file missing → "missing"; file size < 200 → "empty"; file starts with "TIMEOUT" in first 100 chars → "timeout". This collapses THREE mechanically-distinct failure modes into "empty": (a) binary-not-found (~78 bytes, error string `ERROR: [Errno 2] No such file or directory`), (b) rate-limit / quota exceeded (small-but-nonzero output with rate-limit error in stderr), (c) genuinely empty execution (agent ran but wrote nothing). Relaunching case (a) as "empty" is guaranteed to fail identically on every retry, creating a permafail loop. Counter: extend classification to inspect the file's first N bytes for distinguishing tokens — `No such file or directory` / `command not found` → PREREQUISITE_MISSING (abort campaign, do not retry); `rate limit` / `quota exceeded` / `429` → RATE_LIMITED (retry with exponential backoff); otherwise → EMPTY (retry once). Healing: replace the `< 200` heuristic at `resume_failed.py:76` with a content-aware classifier. Register three distinct failure labels in the dashboard output. Related: AP269 (SDR-formula fabrication with proved-contradictory witness) — same pattern at infrastructure layer: an optimistic heuristic masks a deterministic failure.

**AP295 (Dashboard "running" status without liveness check).** `scripts/campaign_dashboard.py` reports campaign status by counting files + file sizes; there is no `pgrep codex > 0` liveness check, so a campaign where every agent failed-fast in the first second and the top-level invocation has already exited is reported IDENTICALLY to an actively-executing campaign. Counter: dashboard status must combine file-level evidence with process-level evidence — show `running` only if (a) recent files exist AND (b) `pgrep -f "codex exec"` returns ≥ 1 PID, OR if a lock-file / PID-file recorded at launch time matches a live PID. Healing: add a `liveness_check()` function at `campaign_dashboard.py` using `psutil` (cross-platform) or `subprocess.run(["pgrep", "-f", "codex exec"])` on POSIX; tag permafail cases as "stalled" with the first-line error message of one sample file. Related: AP252 (Chern-character degree direction) — same anti-pattern of reading a signal in the wrong direction; here, presence-of-files read as presence-of-activity.

**Pattern 407 (first-principles cache, session 2026-04-18; formerly Pattern 222, renumbered 2026-04-21 per AP-CAT-1 after H3-level 222 took the integer first). "Rate-limited agents" may be binary-not-found in disguise.** Trigger: user invokes `resume_failed.py` or equivalent recovery script and reports "relaunch rate-limited agents"; the output directory contains many small (sub-200-byte) files. Regex: `ls <out_dir>/*.md | xargs awk 'NR==1 && length < 200 {print FILENAME}'` yields dozens of hits. Heuristic check: `head -c 200 <out_dir>/<sample>.md` — if the first line contains `ERROR: [Errno 2] No such file or directory` or `ERROR: \[Errno 2\]` the root cause is PREREQUISITE_MISSING, not rate-limiting. Corollary: every invocation of "relaunch rate-limited agents" must inspect ONE sample failure file and verify the error string before relaunching; accepting the "rate-limited" framing without this check produces an infinite resume-loop. Append to `notes/first_principles_cache_comprehensive.md` with the matching regex.

### Wave 14 (2026-04-18 MC5 attack-then-heal, AP296-AP297)

**AP296 (Weight-preservation fallacy on bar differentials where OPE simple-pole summands strictly drop weight).** A bar-complex proof claims the bar differential "preserves total conformal weight on bar words" and uses the claim to decompose $B^{ch}(A_{\leq N})$ as a direct product $\prod_w B^{ch}(A_{\leq N})_w$ over exact weight. The claim is false on every class-M standard-landscape algebra: the OPE weight identity $\mathrm{wt}(a_{(n)} b) = \mathrm{wt}(a) + \mathrm{wt}(b) - n - 1$ strictly decreases total weight on simple-pole summands ($n \geq 1$). Canonical witnesses: Virasoro $T_{(3)} T = c/2$ drops weight $4 \to 0$; affine Kac-Moody $J^a_{(0)} J^b = [a, b]$ drops weight $2 \to 1$. The correct statement is the FILTRATION-PRESERVATION inscribed at `prop:standard-strong-filtration`(iii); the correct Mittag-Leffler argument runs in (cohomological degree $m$, filtration level $w$), not in (cohomological degree $m$, exact weight $w$). Product decomposition over exact weight is illegitimate; the decomposition that survives is the associated graded $\bigoplus_w \mathrm{gr}_w^F$ on the FILTRATION, not the object itself. Canonical violation: Vol I `chapters/theory/mc5_class_m_chain_level_platonic.tex:311-323` Step 3 of `thm:mc5-class-m-chain-level-pro-ambient` wrote "bar differential preserves total conformal weight on bar words" and decomposed the pro-object bar complex as $\prod_w$ over exact weight; this is the AP261 single-vs-bigraded Mittag-Leffler heal applied to the WRONG bigrading (exact weight, not filtration level). Counter: before writing any bar-complex product decomposition $\prod_w B^{ch}(A)_w$, substitute a Virasoro or affine KM bar word at weight $w = 4$ and compute the simple-pole summand of $d_{\mathrm{bar}}$; if the output lands at weight $w' < w$, the decomposition is filtration-level and must be rewritten as a filtered colimit / associated-graded pair. Healing: rewrite using the decreasing filtration $\{F_{\leq w}\}$ with ML stabilisation at $N \geq w$ via `prop:standard-strong-filtration`(iv). Distinct from AP261 (single-index ML vs bigraded ML): AP261 catches "wrong arity of grading"; AP296 catches "wrong OBJECT of bigrading" — bigrading by exact weight fails where bigrading by filtration level succeeds. Related: AP260 (scalar-channel linearity without constructed projector) is the same error class in K-theoretic Chern-character accounting; AP296 is its bar-complex / conformal-weight analogue.

**AP297 (Cross-volume homotopy formula cited to a proof body that does not inscribe it).** A proof in Volume $V_1$ cites an explicit contracting-homotopy formula (e.g. $h_N = \sum_{j} h \cdot m_0^{j-1}$ with $j \in [1, \lfloor N/2 \rfloor]$) as "constructed in the proof of $\langle$external proposition in Volume $V_2\rangle$"; the cited proof body does NOT inscribe the formula, supplying only a finite-stage quasi-isomorphism or weight-completed Milnor/Mittag-Leffler comparison. Reader traces the reference expecting the explicit formula; finds a qualitatively different, weaker argument. Canonical violation: Vol I `chapters/theory/mc5_class_m_chain_level_platonic.tex:267-282` Step 1 cites $h_N = \sum_j h \cdot m_0^{j-1}$ as "constructed in the proof of `prop:bv-bar-class-m-weight-completed`" (Vol II `chapters/connections/bv_brst.tex:2316-2360` + frontier remark at `:2361-2392`); the Vol II proof body supplies weight-completed finite-stage qi at the Milnor/ML level with NO formula of this shape, and the frontier remark at `:2361-2392` in fact states the opposite: "A strict comparison on the raw direct-sum models would have to absorb infinitely many nonvanishing higher corrections by a single ordinary chain-homotopy package. This is what fails." Counter: before citing an explicit homotopy formula to a cross-volume proof, grep the cited proof body for the formula's exact form; zero hits = AP297 fabrication-via-misattribution. Healing: (a) inscribe the explicit homotopy locally as a lemma with full derivation, if the formula is genuinely available; OR (b) downgrade citation to `\ClaimStatusProvedElsewhere` with honest attribution to the weaker result actually inscribed, and rewrite the Step to route through the weaker form; OR (c) sharpen the obstruction per AP266 — state explicitly why a raw single-chain-homotopy cannot absorb the higher corrections, converting the citation into a Beilinson falsification test (any claimed chain-level splitting in the bounded ambient must reproduce the obstruction coefficient). Distinct from AP269 (SDR-formula fabrication with a proved-contradictory witness) because AP297 catches the subtler case where the citation TARGET exists, resolves, is peer-reviewed, but does not inscribe the SPECIFIC FORMULA being cited — the betrayal is at the proof-body level, not at the theorem-existence level. Related: AP190 (hidden imports — cited result does not prove the stated claim); AP272 (unstated cross-lemma via folklore citation); AP283 (formula confabulation in status-table RHS not inscribed in .tex or engine) at the status-table-layer; AP297 is the proof-body-layer analogue of AP283.

**AP298 (Symmetric-group descent from braid relation alone, without the involution / unitarity bridge).** A proof extends a pure-braid representation $P_n \to G$ to the full symmetric group $\Sigma_n \to G$ using the braid relations (YBE, $R_{12} R_{13} R_{23} = R_{23} R_{13} R_{12}$) ALONE; the extension tacitly requires the INVOLUTION relation at one higher codimension — unitarity $R(z) R^{op}(-z) = \mathrm{id}$, or equivalently the squared transposition constraint $s_i^2 = 1$ lifted to the $R$-matrix level. Without unitarity, the extension is a $P_n$-representation, not a $\Sigma_n$-representation. Canonical violation: Vol I `chapters/theory/theorem_A_infinity_2.tex:988-994` entry version of `lem:R-twisted-descent` Step 1 extended the pure-braid representation to $\Sigma_n$ using only YBE; unitarity $R(z) R^{op}(-z) = \mathrm{id}$ was required for the involution $s_i^2 = 1$ but left silent. Yangian rational and $U_q$ trigonometric (generic $q$) satisfy unitarity; elliptic Belavin is convention-dependent (Pauli decomposition satisfies it, Weierstrass $\zeta$-based version may not); non-zero-level Heisenberg has scalar \(R(z)=\exp(k\hbar/z)\), so unitarity is verified by scalar cancellation \(R(z)R(-z)=1\), not because the \(R\)-matrix is trivial. Counter: before writing "the pure-braid representation extends to $\Sigma_n$ via $R$-twisting", state unitarity as a labelled hypothesis (R2) alongside YBE (R1); enumerate the families in which (R2) is verified (Yangian rational, $U_q$ generic $q$, scalar Heisenberg) and those in which it is open (elliptic convention-dependent, face-model-only). Healing: inscribe (R1)+(R2) hypothesis block; add a scope remark listing family coverage. Stronger than AP7 (universal quantifier without restriction audit): AP298 specifically catches the codim-2 (braid) → codim-3 ($\Sigma_n$) extension where the missing ingredient is a definite algebraic relation at the same codimension as the braid, not a scope quantifier. Related: AP253 (inter-volume dependency graph), AP160 (three Hochschild geometry-determines-which).

**AP299 (Mittag-Leffler invoked on filtration X while the hypothesis controls filtration Y; bridge elided).** A proof cites "Mittag-Leffler on the bar-length filtration" under a hypothesis that controls the coradical filtration (or conversely); the bar-cobar $\varprojlim$-duality bridge that converts one filtration into the other is elided, turning a chain of two theorems into one with a silent dependency. Canonical violation: Vol I `chapters/theory/theorem_A_infinity_2.tex:824-844` invokes Mittag-Leffler for the bar-length filtration under hypothesis (H2); (H2) controls the CORADICAL filtration on the codomain coalgebra, not directly bar-length. The bridge via bar-cobar $\varprojlim$-duality is a separate lemma and was implicit. Counter: for every ML / $\varprojlim$ invocation, state (a) which filtration is being stabilised, (b) which hypothesis controls which filtration, (c) the explicit bridge if the two differ. The bridge is itself a proposition (in this case, bar-cobar duality on conilpotent-complete loci under (H2)+(H3)); either inscribe it as a lemma or cite explicitly. Healing: add a scope remark inscribing the bar-length ↔ coradical bridge under (H2), quantitatively characterising the obstruction to extending beyond (H3) as a $\varprojlim^1$-class in $H^{\bullet+1}$. Distinct from AP261 (single-index vs bigraded ML) and AP296 (exact-weight vs filtration-level bigrading): AP299 catches a third class — correct number and type of indices, but wrong UNDERLYING filtration object on which ML runs. Related: AP37 (spectral-sequence page from full differential, not pole-order heuristic); AP200 (transfer-theorem gap — H*(A) results applied to A); AP299 is the filtration analogue of the transfer-theorem gap.

**AP300 (In-file ProvedHere-vs-retracted-mechanism drift via a non-cross-referenced retraction remark).** A lemma in a chapter carries `\ClaimStatusProvedHere` and a proof body whose Step-$k$ uses a specific proof mechanism (e.g. "finite-dim Zhu implies uniform bound on collision residues across $r$"); later in the SAME file a `\begin{remark}` inscribes a retraction of that mechanism (e.g. `rem:zhu-bounded-massey-retraction` stating the Zhu-boundedness implication is FALSE on the logarithmic sector); the lemma's proof body carries no `\ref{rem:...}` to the retraction, and the lemma's status tag was not downgraded to match. A casual reader sees the lemma as proved; only by reading the chapter end-to-end does the silent drift surface. Canonical violation: Vol II `chapters/theory/logarithmic_wp_tempered_analysis_platonic.tex` Step 1 of `lem:wp-ww-subchannel-tempered` (lines 369-446, `\ClaimStatusProvedHere`) uses the exact mechanism retracted at `rem:zhu-bounded-massey-retraction` (lines 481-514), ~35 lines later in the same file, with no `\ref` bridge. Counter: before accepting any in-chapter `\ClaimStatusProvedHere` lemma, grep the same file for `\begin{remark}` blocks containing "RETRACT" / "retract" / "incorrect" / "FAILS" / "does not apply"; cross-check every remark against every ProvedHere proof mechanism above it. Healing: (a) downgrade the lemma to `\ClaimStatusConditional` or `\ClaimStatusConjectured` with an explicit `\ref{rem:retraction}` in the proof body stating the sharpened scope; (b) rewrite the proof to route around the retracted mechanism through a disjoint chain; (c) delete the lemma and its consumers atomically (AP5). Distinct from AP256 (aspirational-heal CLAUDE.md-vs-manuscript drift) and AP271 (reverse drift CLAUDE.md-lags-manuscript): AP300 is an ENTIRELY IN-FILE drift where both the lemma and the retraction remark live in the same chapter; the metacognitive layer is not involved. Distinct from AP149 (resolution-propagation failure across files): AP149 is inter-file propagation gap; AP300 is same-file propagation gap. Stronger than AP4 (ClaimStatusProvedHere requires the proof to prove the claim) because AP300 surfaces the subtler case where the proof WAS valid at inscription time but the author later noticed a counterexample, inscribed the retraction in a remark, and forgot to tighten the lemma tag.

**Pattern 408 (first-principles cache, session 2026-04-18; formerly Pattern 223, renumbered 2026-04-21 per AP-CAT-1). Exact-weight vs filtration-level bigrading for bar-complex Mittag-Leffler.** Trigger: bar-complex proof writes `B(A_{≤N}) = ∏_w B(A_{≤N})_w` with product over "total conformal weight"; OR writes Mittag-Leffler in $(m, w)$ with $w$ = "exact bar-word weight". Regex: `\\prod_w.*B\^?\{ch\}.*A` in `.tex`. Counter-derivation: for any class-M algebra with a simple-pole summand in its OPE (Virasoro $T_{(3)} T = c/2$, affine KM $J_{(0)} J = [J, J]$), evaluate the bar differential on a two-word bar word at weight $w = 4$; the simple-pole summand strictly decreases weight. Replace exact-weight product with filtration-level filtered colimit + associated graded. Primary source: Etingof-Frenkel-Kirillov Vol II (OPE weight identity); Frenkel-Ben-Zvi "Vertex Algebras on Algebraic Curves" §3.5. Append to `notes/first_principles_cache_comprehensive.md` alongside Pattern 407 (formerly Pattern 222).

### Wave 15 (2026-04-18 /loop + async-Agent harness audit, AP301-AP304)

Four infrastructure anti-patterns observed by the main thread during the 2026-04-18 `/loop every 3m` supervision session concurrent with Waves 14 (MC5 swarm, AP296-AP299) and AP300 (Theorem A swarm in-file retraction drift). Wave 15 concerns the `/loop` + async-Agent supervision harness rather than mathematical content.

**AP301 (/loop interval shorter than median agent runtime).** When a recurring `/loop` is scheduled with interval `Nm` intended to supervise long-running subagents (rectification swarms, attack-heal agents, `scripts/rectify-all.py` campaigns), and N is smaller than the p50 agent runtime, every firing finds the same agents still running and performs no productive work while paying the full /loop cost: CLAUDE.md (~1000 lines) re-parsed into context, Bash spawns for state inspection, a full tool-call cycle, completion summary. The cron becomes a poll-heavy supervisor rather than a work-driver. Canonical violation: 2026-04-18 session, user scheduled `/loop every 3m` to supervise 5 attack-heal agents whose p50 runtime was ~6-8 min (Theorem H: 347s; Theorem B: 445s); iterations 2-3 found agents running with nothing productive to do except inscribe meta-notes. Counter: match /loop interval to p50 agent runtime. Rule of thumb: 3m for poll-heavy tasks with fast observable state changes (build logs, git status); 10-15m for attack-heal swarms on single theorems; 30-60m for `/rectify-all` campaigns spanning multiple chapters. When the user specifies an interval that mismatches the task, either ask or escalate to dynamic-mode (ScheduleWakeup-based self-pacing) that can wait on actual completion notifications. Related: AP254 (closure-date commit-floor) — both are mismatches between rhetorical timing claims and actual work throughput.

**AP302 (`replace_all: false` non-unique-string collision).** Edit tool call with `replace_all: false` on a string appearing $N > 1$ times in the file fails with "Found N matches, but replace_all is false". Common in structured documents with repeating section headers (`### Attribution`, `### References`, `### Summary`), or in CLAUDE.md-style manuscripts where pattern-specific counters repeat. Counter: before any Edit on a low-specificity anchor, grep the file for the target string; if $\geq 2$ hits, extend `old_string` with surrounding unique context. Canonical violation: 2026-04-18 session, Edit targeting `### Attribution\n\nNo AI attribution...` in `notes/first_principles_cache_comprehensive.md` hit 14 instances; required re-targeting with a 2-line pre-context anchor. Healing: standardise repeating headers to include entry-local identifiers (`### Attribution (Pattern N)`), OR adopt anchor-comment conventions (`<!-- END PATTERN N -->`). Related: AP145 (restructuring debt) — both flag cross-entry anchor collisions as templating debt.

**AP303 (Async-Agent output-file forbidden-read sentinel).** The `Agent` tool with `run_in_background: true` returns an agent ID and an output-file path that the caller is EXPLICITLY instructed not to read ("reading will overflow context"). The only supervision signal is an asynchronous completion notification. A /loop supervising multiple async agents cannot poll-inspect file sizes or transcript content to count "running" agents; it must either (a) track agent IDs in a local ledger at launch and decrement on each `<task-notification>` with `status: completed`, or (b) inspect the symlink mtime as a liveness proxy (stale mtime > threshold → likely stalled). Reading the jsonl directly to count tokens or infer progress is forbidden. Counter: at launch, append each agent ID + task to a local JSON ledger (e.g. `.claude/agent_ledger.json`); on each /loop firing, compute running count from completion-notification history. Canonical violation: 2026-04-18 cron iteration 2 attempted to count running agents by listing task output files; the listing does not distinguish "running" from "completed-but-notification-pending" without reading the forbidden jsonl. Related: AP295 (dashboard "running" without liveness check) — same pattern at subagent layer.

**AP304 (Concurrent-agent AP-numbering collision in shared metacognitive file).** When multiple subagents + the main thread edit the same file (CLAUDE.md, FRONTIER.md, MEMORY.md) concurrently during a single /loop session, each agent reads the file at launch and writes based on that snapshot. Without a reservation protocol, two agents inscribing "the next available AP number" both claim the same number. Canonical violation: 2026-04-18 cron iteration 2, main thread attempted to write AP296-298 into CLAUDE.md Wave 15; MC5 swarm had already inscribed AP296-299 as Wave 14; a Theorem-A-or-adjacent swarm inscribed AP300 concurrently; main-thread numbering bumped to AP301-304. Counter: before inscribing, re-read the AP-index tail of CLAUDE.md; reserve AP number blocks per agent at launch time by convention (MC5 = AP296-310, Theorem-A = AP311-325, Theorem-B = AP326-340, etc), inscribing the reservation explicitly in the subagent prompt. Alternatively, use subagent-assigned AP prefixes (MC5-AP1, THMA-AP1) promoted to global AP numbers only by the main thread at reconciliation. Related: AP149 (resolution-propagation failure) — both reflect the programme's lack of atomic cross-agent state mutation. This Wave 15 inscription is itself post-hoc reconciliation of an AP304 collision (TWO Edit failures on "File has been modified since read" as the MC5 and Theorem-H swarms wrote concurrently).

**AP305 (Pessimistic CLAUDE.md status-row drift — overstates manuscript deficit).** Inverse of AP271 (reverse drift where CLAUDE.md overstates manuscript success). Here the metacognitive layer carries a MORE PESSIMISTIC diagnosis than the manuscript supports: a status row flags "standalone NOT `\input`-ed" or "phantom file" when `main.tex` actually does `\input` the target and the `\label` resolves with a full proof body; OR a row claims "N of M strengthenings phantom" with M inflated; OR CLAUDE.md claims a theorem is "advertised-but-not-inscribed" (AP241) when the theorem is fully inscribed in a compiled chapter. Downstream effect: adversarial audit agents CHASE a non-existent defect, waste cycles on heals that amount to rewriting the CLAUDE.md row rather than editing the manuscript, and risk downgrading a validly-proved theorem based on stale metacognitive gloss. Canonical violation: 2026-04-18 chiral QG equivalence audit returned with headline "standalone NOT \input-ed; 35 consumers render [?]" diagnosis on CLAUDE.md Chiral QG equiv row; the source state was `main.tex:1394` `\input{chapters/theory/ordered_associative_chiral_kd}` with `\label{thm:chiral-qg-equiv}` at `:8404` carrying full proof; 34 of 35 cross-chapter refs resolve. Of 8 advertised strengthenings only 1 (`thm:glN-chiral-qg`) is fully inscribed, 1 is an AP286 tactical phantomsection, 2 have live consumers targeting a phantom label (`thm:chiral-qg-equiv-ordered`), 4 are zero-consumer phantoms — honest count "1 + 1 + 2 + 4 of 8" not "0 of 8 phantom + `[?]` cascade". Counter: before pessimistic diagnoses in CLAUDE.md status rows, run a three-step verification — (i) grep `main.tex` for `\input{path}` of the claimed-phantom file; (ii) grep the input target for `\label{target}`; (iii) run `make fast` or the build proxy and count `[?]` emissions in the log. Zero hits on (i) or (ii) justifies "phantom" diagnosis; live-build `[?]` count calibrates the consumer-resolution claim. Healing: rewrite the row to honest per-object accounting (X of M inscribed, Y phantom with zero consumers, Z phantom with live consumers needing retarget). Distinct from AP271 (reverse drift lagging success) and AP256 (aspirational-heal advertising advance): AP305 is the OPPOSITE failure mode — CLAUDE.md is more pessimistic than manuscript. Failure mode: AP271 causes inherited false confidence; AP305 causes adversarial-audit false alarm and inscription-by-downgrade risk. Related: AP283 (formula confabulation in status-table RHS); AP280 (three-step epistemic inflation). AP305 is the pessimistic-drift sibling of AP280, where the inflated claim is deficit-count rather than achievement-count.

**AP306 (Reserved-AP-block under-utilisation and stranding).** AP304's reservation protocol (MC5 = AP296-310, Theorem-A = AP311-325, 15 numbers per block) under-utilises: MC5 used AP296-299 (4 of 15); Theorem A cited pre-existing healed items rather than inscribing new ones (0 of 15); Theorem H claimed "no new APs needed"; Theorem D inscribed a remark but no AP; Theorem B healed phantoms without inscribing. Net over cron iteration 2: ~5 of 85 reserved numbers touched, 80 stranded. Counter: reservation blocks should be 5 per agent, not 15; stranded numbers should be explicitly re-released in a reconciliation step (`.claude/ap_reservations_YYYYMMDD.json` with `released` field) rather than permanently dead. Alternatively, skip block reservation entirely and require every agent to re-read the AP-tail + use the smallest available integer — this trades collision risk for number density. Canonical violation: 2026-04-18 session, AP300-AP310 (Theorem-A + MC5 tails) mostly unused; AP326-AP340 (Theorem-B block) entirely unused; AP311-AP320 (Theorem-A block) entirely unused; AP341-AP360 (Theorem-C block) status TBD. Healing: in future /loop swarms, (a) reduce block size to 5; (b) after swarm completion emit an "AP reconciliation" commit compacting numbering or marking ranges `[RELEASED 2026-04-18, reuse permitted]`; (c) eventually migrate from integer AP-numbering to content-addressed AP-keys (`ap:filtration-level-vs-exact-weight-ML`, `ap:concurrent-edit-collision`) that are stable across agent boundaries. Related: AP145 (restructuring = O(N) debt — block-reservation overhead is per-swarm).

**AP307 (Repo-root report-file proliferation).** Each attack-heal swarm agent emits a `attack_heal_<target>_<date>.md` report file. Over the 2026-04-18 iteration-2 /loop, 9 such files accumulated at repo root (A, B, D, H, MC5, C, MC3, Chiral QG, W(p), MC4^0). Without a cleanup schedule, these proliferate across sessions, clutter `git status`, and are not referenced anywhere in the manuscript or `notes/`. Counter: direct swarm reports into `adversarial_swarm_YYYYMMDD/attack_heal_<target>.md`, matching the existing `adversarial_swarm_20260416/` / `adversarial_swarm_20260417/` convention (CLAUDE.md "Today's Reconstitution Headlines" already references this directory pattern). Healing: at swarm completion, `mkdir -p adversarial_swarm_$(date +%Y%m%d) && mv attack_heal_*_$(date +%Y%m%d).md adversarial_swarm_$(date +%Y%m%d)/`. Update subagent prompts to "Report in adversarial_swarm_YYYYMMDD/attack_heal_<target>.md" rather than repo root. MC4^0 swarm prompt (launched 2026-04-18 cron iteration 3) adopts this convention going forward. Related: AP288 (session-ledger stale narrative) — same pattern, different layer.

**AP308 (Cron-fire interleaving with ongoing main-thread Edit sequence).** The `/loop every 3m` cron fires independently of main-thread tool-call state. When the main thread is mid-Edit-sequence when a cron fire arrives, the cron-fire input is enqueued but not processed until the current turn completes; the main thread re-enters the loop prompt assuming "new iteration, fresh state" while completion notifications from the previous iteration's agents have mutated the state mid-turn. Canonical violation: 2026-04-18 session iterations 2-3, main thread's Wave-15 inscription Edit failed twice with "File has been modified since read" as MC5 + Theorem-A + Chiral-QG swarms wrote into CLAUDE.md concurrently. Counter: accept the interleaving as benign for the "launch N agents if under K" rule — it is idempotent with respect to brief over-launching (a 6-agent moment self-corrects as agents complete). For CLAUDE.md edits, always re-Grep the specific anchor immediately before Edit and re-read a small window to obtain a fresh snapshot; do not reuse a Read result across multiple tool calls separated by agent-completion notifications. Related: AP300/AP301 (loop-interval mismatch) — interleaving is the operational manifestation of that structural mismatch.

### Wave 17 (2026-04-18 super-Y + ϱ_BP + Trinity-K attack-heal, AP309-AP312)

Four mathematical-content patterns surfaced when the second batch of attack-heal agents returned in the 2026-04-18 /loop session. Wave 17 is the content-layer complement to Waves 15/16 infrastructure-layer additions.

**AP309 (Primary-source citation for a strictly weaker claim, silently extrapolated to the stronger claim).** A proof cites `\cite[Theorem N]{AuthorYear}` where the cited theorem genuinely exists, is correct, and proves a SUBSET of what the proof needs; the proof silently promotes the citation to the stronger claim without inscribing the extrapolation. Canonical violation: Vol II `chapters/theory/super_chiral_yangian.tex:622-702` Step 3 of `thm:super-complementarity-max-mn` cites Gow 2006 Thm 5.1 for the Berezinian shift magnitude $+\max(m, n)$; Gow's theorem establishes CENTRALITY of $\mathrm{sBer}(T(u))$ in $Y(\mathfrak{gl}(m|n))$, NOT the shadow-depth shift magnitude under super-Sugawara + Feigin-Frenkel involution. The proof promotes "centrality" to "shift magnitude = $\max(m, n)$"; the bridge is missing. Counter: every `\cite[Thm N]{X}` at a load-bearing step must be pattern-matched — (a) extract what the cited theorem states verbatim; (b) compare to the inference the proof draws; if the inference is strictly stronger, inscribe the bridge locally or downgrade the claim. Distinct from AP272 (unstated cross-lemma via folklore citation — cited paper does not state the lemma at all) and AP281 (bibkey naming-drift — citations do not resolve): AP309 is the subtler case where the citation is correct but too weak. Related: AP190 (hidden imports — cited result does not prove the claim). Healed 2026-04-18 via theorem split: `thm:super-complementarity-supertrace-zero` (PROVED) + `conj:super-complementarity-berezinian-max-mn` (CONJECTURED) + `rem:berezinian-shift-open`.

**AP310 (HZ-IV three paths bibliographically disjoint but share an unverified intermediate step).** Stronger than AP277 (vacuous test body) and AP288 (label-disjoint-but-computation-identical). The three paths genuinely compute different quantities from different literatures; `assert_sources_disjoint` passes on bibliography keys; the decorator's `disjoint_rationale` reads cleanly; tracing each computation shows all three invoke the SAME unverified intermediate identity. Canonical violation: Vol I `chapters/examples/yangians_foundations.tex:103-190` `lem:super-trace-berezinian-bridge` HZ-IV decorator lists Nazarov 1991 + Gow 2006 + Arnaudon-Crampé 2003; each path cites a different primary source and performs a different numerical check; all three paths ultimately invoke the shift-magnitude claim $+\max(m, n)$ without deriving it independently. `disjoint_rationale` stated "bibliographically disjoint"; truth is bibliographically disjoint AND computationally sharing. Counter: for every HZ-IV decorator, for each of the three paths list the load-bearing intermediate identities; cross-tabulate; any identity appearing in $\geq 2$ paths without an independent proof is a shared unverified step. Healing: (a) find a genuinely disjoint fourth source closing the shared step; (b) inscribe the shared step as a separate proposition with its own HZ-IV decoration; (c) downgrade the parent theorem to `\ClaimStatusConditional`. Distinct from AP243 (HZ-IV V1-input-to-V3 bibliographic entanglement) and AP288 (label-disjoint computation-identical): AP310 is bibliographically disjoint AND computationally distinct at the outer layer but shares an inner load-bearing step. Infrastructure recommendation: extend `compute/lib/independent_verification.py`'s `assert_sources_disjoint` to track declared-intermediate-identity disjointness; register a `shared_intermediate` field agents must populate.

**AP311 (Two mathematical invariants hidden under one symbol — per-family values diverge on non-principal strata).** A programme-level scalar symbol (here $\varrho_A$, the leading coefficient in $\kappa + \kappa^! = \varrho_A \cdot K$) is used without qualifier. Two plausible definitions coexist: $\varrho_{\mathrm{lin}}$ = leading coefficient of $\kappa$ linear in $c$; $\varrho_{\mathrm{gen}}$ = signed harmonic sum $\sum_\alpha (-1)^{\varepsilon_\alpha} / h_\alpha$ over strong generators. They agree on principal $\mathcal{W}_N$ (both give $H_N - 1$) but diverge on Bershadsky-Polyakov: $\varrho_{\mathrm{T\text{-}line}} = 1/2$ vs $\varrho_{\mathrm{gen}} = 1/6$ (ratio 3 from $G^\pm$-pair-counted-once vs counted-twice convention). Agents argue "$\varrho_{\mathrm{BP}} = 1/6$ is back-fit" while the manuscript treats $\varrho = 1/2$ as T-line-canonical, producing a phantom open frontier ("structural origin open per Wave-1 F1.b") that is a symbol collision, not an open problem. The structural origin DOES exist: $\varrho_{\mathrm{gen}}$ is the Kac-Roan-Wakimoto signed harmonic sum over the minimal-nilpotent orbit profile; already inscribed at `cor:anomaly-ratio-ds` in `chapters/examples/landscape_census.tex:1450-1474` but not recognised programme-wide. Counter: like HZ-7 bare-$\kappa$ discipline in Vol III, require $\varrho$ to carry a subscript ($\varrho_{\mathrm{lin}}$, $\varrho_{\mathrm{gen}}$, $\varrho_{\mathrm{KRW}}$, $\varrho_{\mathrm{mod\text{-}char}}$, $\varrho_{\mathrm{T\text{-}line}}$) in any programme-wide invariant identity. Healing: inscribe the subscript discipline; retract "$\varrho_A$ structural origin is open" from CLAUDE.md frontier list; register the KRW-orbit-sum derivation as a theorem. Related: AP244 (overcounted foundational terms) is the structural sibling; AP234 (two Koszul conductors under letter $K$) is the $K$-sibling. Pre-emptive discipline: every programme-level scalar invariant should be introduced with ALL its plausible definitions enumerated.

**AP312 (Three-way cross-file scalar-value contradiction without bridge remark).** A scalar invariant (here $K(\mathcal{H}_k)$, the Koszul conductor of the Heisenberg chiral algebra) carries THREE distinct numerical values across three files, with no reconciliation remark. Canonical violation: 2026-04-18 Trinity K audit surfaced $K(\mathcal{H}_k)$ inscribed as (a) $-k$ at Vol I `chapters/theory/universal_conductor_K_platonic.tex` summary table, (b) $-2$ at Vol I `chapters/theory/kappa_conductor.tex:192`, (c) $+2$ in CLAUDE.md derivation prose. All three values correspond to legitimate conventions (spin-1 bosonic FMS ghost charge $= -2$; BRST sign $= +2$; $\kappa + \kappa^! = 0$ with $\kappa = -k/2$ giving $K = -k$), but no bridge remark names the canonical convention. Stronger than AP245 (statement-proof-engine numerical agreement — three levels of a single proposition) and AP238 (statement/proof contradiction within a single environment): AP312 is a THREE-FILE divergence across distinct propositions, no single environment to audit. A reader picks one of three values at random. Counter: for every scalar invariant appearing in $\geq 2$ files, enforce a single canonical value with a named convention inscribed at one programme-level source of truth; every other inscription either uses the canonical value OR carries a bridge remark showing the convention conversion. Healing (applied 2026-04-18): pinned $K(\mathcal{H}_k) = -2$ with `rem:uc-heisenberg-scope-distinction` distinguishing $K_g = -2$ (ghost-charge) from $\kappa + \kappa^! = 0$ ($\kappa$-complementarity); retired the $-k$ drift. Related: AP234, AP311, AP280. Preventative discipline: a `conventions.tex` inscribing every programme-level scalar with its canonical value and conversion rules.

**AP313 (Proof via post-hoc prefactor normalization absorbed into RHS symbol definition).** A theorem states a clean identity $f(g) = B_{2g-2} / (g-1)$ (or analogue); the proof proceeds only up to a universal combinatorial prefactor (e.g. $2^{g-1} (g-1)! / (2g-2)!$), and Step-$k$ of the proof writes "after absorbing the universal prefactor into the definition of $R_{g-2}$"; the boxed statement then holds by convention-redefinition rather than by a clean derivation. Canonical violation: Vol I `chapters/theory/z_g_kummer_bernoulli_platonic.tex:208-308` `thm:z-g-leading-coefficient-bernoulli` Step 5 at :282-307 admits the identity holds "up to the combinatorial factor $2^{g-1} (g-1) / (2g-2)!$ — not a clean identity", then the boxed form at :302 absorbs the prefactor into the definition of $R_{g-2}$. The low-genus verification at :312-322 works case-by-case; the Kummer witnesses at $g = 7, 9$ survive because the prefactors $c_7, c_9$ are computed explicitly coprime to $\{691, 3617\}$, but the GENERAL-$g$ clean-identity claim is normalization-absorbed, not derived. Counter: before boxing any identity whose proof includes the phrase "after absorbing the prefactor into the definition of [symbol on RHS]", rewrite the identity with the prefactor visible as $f(g) = c_g \cdot B_{2g-2} / (g-1)$ with $c_g$ explicitly computed; verify every downstream consequence (here the Kummer-witness coprimality) uses the explicit prefactor form rather than the normalization-absorbed boxed form. Healing: (a) rewrite the step to prove the clean identity directly via explicit prefactor tracking; (b) restate the boxed form with $c_g$ explicit and verify downstream consequences from that form; (c) sharpen the obstruction (AP266) if prefactor-free identity does not hold: exhibit $c_g$ as an explicit combinatorial invariant with Beilinson falsification test. Distinct from AP259 (tautological-by-definition — $A = B$ proved via defining $A$ by projection through $B$): AP313 is the subtler case where LHS and RHS start as different objects, a clean identity is claimed, but the proof closes only after absorbing a combinatorial normalization into one side. Distinct from AP280 (three-step epistemic inflation — remark → standalone → headline): AP313 is within a single theorem environment, normalization-absorbed at the proof-body level. Related: AP238, AP245, AP308, AP311.

**Pattern 409 (first-principles cache, session 2026-04-18; formerly Pattern 224, renumbered 2026-04-21 per AP-CAT-1). Async-Agent running-count is not directly observable.** Trigger: a /loop directive says "if under N agents, launch more"; the loop tries to compute the running count by inspecting task output files, subagent jsonl symlinks, or mtime. Counter: the reliable running-count is the difference between (agents launched over current /loop session) minus (received `<task-notification>` with `status: completed` entries). All three observable-file heuristics (size, mtime, symlink presence) fail: size stays tiny throughout the agent's run; mtime is set at symlink creation only; symlink exists from launch until Claude process cleanup. Regex-equivalent: maintain a ledger `.claude/agent_ledger_YYYYMMDD.json` with `{launched: [ids], completed: [ids]}` updated atomically on launch and on notification. Append to `notes/first_principles_cache_comprehensive.md` alongside Pattern 407 (formerly 222) + Pattern 408 (formerly 223).

**AP313 (Truncated Agent-result masquerading as clean completion).** When the `Agent` tool with `run_in_background: true` emits a `<task-notification>` with `status: completed`, the `result` field appears to be the agent's final reply. In practice the result is often the LAST LINE of the agent's last turn, truncated mid-sentence, mid-Edit-retry, or mid-"now update CLAUDE.md..." announcement, NOT a clean mission summary. The `status: completed` flag only indicates the subagent's turn ended; it does not confirm the announced heals landed on disk, that files under edit were saved, or that the agent reached its planned Phase 5 propagation step. Canonical violation: 2026-04-18 cron iteration 4, W(p) Massey swarm returned `result: "Note: $\\mathcal{SF}^{\\Z_2}$ not $\\mathcal{SF}^{\\Z/2\\Z}$..."` (mid-Edit-retry narration); Theorem C swarm returned `result: "Those are pre-existing lines in the standalone, not my edit. My edits are correct. Now update CLAUDE.md Theorem C row:"` (mid-update-announcement, CLAUDE.md update NOT reliably reflected in the visible transcript); Chiral QG swarm returned `result: "Acknowledged — no commits."` (bare terminator, though AP305 was inscribed mid-run). Counter: for every `<task-notification>` with `status: completed` and a short/mid-sentence `result`, treat the agent as AMBIGUOUS-COMPLETED; re-verify each claimed heal via (a) `ls adversarial_swarm_YYYYMMDD/attack_heal_<target>.md` or `wave1_<target>_attack_heal.md` for the Phase 5 deliverable, (b) `grep -rn '\label{thm:<target>}'` for claimed inscriptions, (c) `git status` / `git diff` in the agent's worktree. If any returns empty, either re-launch with a prompt noting the prior truncation or downgrade to "attempted" in the session ledger. Distinct from AP255 (phantom file + phantomsection mask) which catches manuscript-layer phantoms; AP313 catches agent-layer partials. Related: Pattern 409 (formerly 224; running-count opaque) — same class of "completion status is opaque" at different layer.

**Pattern 410 (first-principles cache, session 2026-04-18; formerly Pattern 225, renumbered 2026-04-21 per AP-CAT-1). `status: completed` on async-Agent is not agent success.** Trigger: a `<task-notification>` arrives with `status: completed`, and the `result` field is a short fragment ("Acknowledged.", "Now update...", "Those are pre-existing lines...", a bare continuation sentence with no mission summary, an inline diff line, or a JSON error snippet). Counter: do NOT trust the status field as a success signal. Verify via (i) `ls adversarial_swarm_YYYYMMDD/attack_heal_<target>.md` for the Phase 5 deliverable, (ii) `grep -rn '\label{thm:<target>}' chapters/ standalone/` for claimed inscriptions, (iii) `git diff --stat` in the agent's worktree for actual file modifications. If any of the three returns empty, treat completion as AMBIGUOUS-PARTIAL and either re-launch or downgrade status to "attempted". Append to `notes/first_principles_cache_comprehensive.md` alongside Pattern 407-409 (formerly 222-224).

**AP314 (Inscription-rate outpaces audit capacity).** When a /loop iterates rapidly and each iteration inscribes 1-4 new APs (Waves 13, 14, 15 in 2026-04-18 session added 13 APs + 4 cache patterns in ~25 minutes), new APs accumulate faster than they can be audited for cross-consistency, redundancy with existing APs (many of the 300+ already in the catalogue), or practical validation via subsequent Edit sessions. The metacognitive layer grows without the discipline of showing each AP has actually fired in a real Edit decision and prevented a real error. Canonical violation: 2026-04-18 session added AP293-AP314 = 22 APs + 4 cache entries in under 30 minutes, interleaved with swarm-inscribed AP296-AP312 (mathematical-content APs), total ~40 new items. Many overlap in scope (AP293/AP294/AP295 all describe resume_failed.py failure modes; AP301 and AP308 both describe /loop timing issues). No consolidation pass was performed before the next iteration. Counter: throttle inscription to at most 1 AP per /loop iteration; use remaining iteration budget to (a) audit the prior AP for actual firing in a real Edit decision, (b) identify redundancy with existing APs and merge, (c) update the HOT ZONE if the new AP exceeds 5 waves / 50+ instances threshold. Related: AP301 (/loop interval shorter than median agent runtime) — AP314 is the inscription-workload consequence of AP301's over-firing. Retroactive discipline for the 2026-04-18 session: AP293-AP295 (resume-infrastructure trio) could merge into a single meta-AP; AP301 + AP308 (loop timing + edit interleaving) could merge; AP303 + Pattern 409 (formerly 224) + AP313 + Pattern 410 (formerly 225) (async-Agent observability) could consolidate into one umbrella AP with three sub-patterns.

### Wave MC4-swarm (2026-04-18 MC4^0 attack-then-heal, AP421)

**AP421 (Standalone-only MC4^0 inscription claiming Vol~I chapter inscription — AP269+AP255 composite).** A completion-machinery theorem advertised in CLAUDE.md status-table as UNCONDITIONAL claims inscription-home in a Vol~I chapter; the chapter does not contain the theorem body, any named SDR deformation-retract data, any homotopy formula, or any resonance-locus analysis; the theorem actually lives only in a standalone `.tex` file that is not `\input` into `main.tex`; a phantomsection stub in `chapters/frame/preface.tex` resolves `\ref{}` invocations without build errors, masking the absence. The composite failure mode is stronger than AP255 alone (phantom-file-in-consumers) and AP269 alone (SDR-formula fabrication) because the false chapter-inscription ADVERTISEMENT is contained inside the standalone's proof sketch, so a reader who opens only the standalone sees a cross-chapter inscription claim that grep-verifies to zero hits. Canonical violation: `thm:n4-mc4-zero-unconditional` in `standalone/N4_mc4_completion.tex:928-953` (prior to 2026-04-18 heal) carried the proof-sketch line "Inscribed as the MC4${}^0$ upgrade in Vol~I `chapters/theory/bar_cobar_adjunction_curved.tex`, referenced at `thm:completed-bar-cobar-strong` as part of the 2026-04-16 wave closure"; the cited chapter exists and proves `thm:completed-bar-cobar-strong` (strong-completion-tower MC4 base), but contains zero `Wakimoto` / zero `h_htpy` / zero `SDR` hits; the standalone is not `\input` into `main.tex`; `chapters/frame/preface.tex:5130` carries `\phantomsection\label{thm:completed-bar-cobar-strong}` which is a legitimate label-share but reinforces the misdirection. Counter: before accepting any CLAUDE.md status-row "Inscribed as... in `chapters/theory/X.tex`" self-advertisement that lives inside a standalone file, (a) grep the cited chapter for the mechanism keyword (Wakimoto, SDR, screening, homotopy formula), (b) grep `main.tex` for `\input{...}` of the standalone, (c) grep `chapters/frame/preface.tex` for `\phantomsection\label{foo}` where `foo` is the target theorem — any of the three missing is a composite AP255+AP269 failure. Healing (applied 2026-04-18): downgraded theorem to conjecture (`conj:n4-mc4-zero-generic-parameters`), renamed label, retracted the false cross-chapter inscription advertisement, inscribed AP266 sharpened-obstruction `prop:mc4-zero-wakimoto-sdr-obstruction` recording the $(\Psi-1)/\Psi$ chiral-coproduct cocycle class from `prop:ff-screening-coproduct-obstruction` in `ordered_associative_chiral_kd.tex:10176-10297` as the load-bearing negative witness for the $\cW_N$ route. Status table MC1-4 row rewritten to reflect MC4${}^+$ UNCONDITIONAL (classes G/L/C) + MC4${}^0$ CONJECTURAL (class M principal). Audit trail: `adversarial_swarm_20260418/attack_heal_mc4_0_20260418.md`. Related: AP255 (phantom file), AP269 (SDR fabrication), AP266 (sharpened obstruction healing), AP271 (reverse drift), AP280 (three-step epistemic inflation at headline level).

### Wave N2-SCA (2026-04-18 N=2 superconformal silent non-coverage, AP501-AP502)

**AP501 (Koszulness without MC3 — structural independence of two properties of a chiral algebra).** Chiral Koszulness (bar cohomology concentrated in degree~1) and MC3 thick generation of the Drinfeld-Kohno factorization category ($\mathrm{thick}(\mathrm{ev}\text{-}\mathrm{type}) = \mathrm{DK}_\mathfrak{g}$) are structurally independent. A chiral algebra can satisfy Koszulness via PBW collapse + free strong generation (Adamović-type input) without MC3 being provable by any of the five-family mechanisms. Canonical instance: the $\mathcal{N}=2$ superconformal algebra at generic level $k \neq -2$ is chirally Koszul (`prop:n2-koszulness` at `chapters/examples/n2_superconformal.tex:301`, via PBW collapse at $E_2$ + Adamović 1999 generic free strong generation), yet MC3 for $\mathrm{SCA}_c$ is OPEN (`conj:mc3-n2-sca-kazama-suzuki` at `chapters/examples/n2_superconformal.tex`): the five-family mechanisms (type A asymptotic prefundamentals, B/C/D reflection-equation Shapovalov, classical+simply-laced Chari-Moura, elliptic theta-divisor, super-Yangian parity-balance) do not apply verbatim to the $\mathcal{N}=2$ SCA, which has mixed bosonic/fermionic strong generators and no $R$-matrix presentation over a simple Lie algebra. Counter: never infer MC3 from Koszulness, nor vice versa; treat them as separate inscription obligations. Before writing "the chiral algebra $A$ is Koszul, hence MC3 holds on $\mathrm{DK}_A^{\mathrm{ev}}$" (or the converse), stop: cite the relevant theorem and verify both assertions are independently inscribed. AP501 is the scope-discipline variant that prevents conflation; the specific healing for $\mathcal{N}=2$ SCA is inscribed at `sec:n2-mc3-status`. Related AP244 (overcounted foundational terms): Koszulness and MC3 are genuinely distinct on the working locus, not a terminological redundancy.

**AP502 (Coset MC3 non-inheritance — thick generation does not descend under branching).** If MC3 holds for the ambient chiral algebra $A \otimes B$, it does NOT follow that MC3 holds for the coset $\operatorname{Com}(C, A \otimes B)$ where $C \subset A \otimes B$ is a sub-chiral algebra. The obstruction: branching of an evaluation-type $(A \otimes B)$-module under $C$ decomposes into $C$-isotypic components, each of which is a $\operatorname{Com}(C, A \otimes B)$-module; the isotypic components need not be evaluation-type over the coset, and their thick closure need not exhaust the coset's Drinfeld-Kohno category. Canonical instance: the Kazama-Suzuki coset $\mathrm{SCA}_c \simeq \operatorname{Com}(\mathrm{U}(1), \widehat{\mathfrak{sl}}_{2,k} \otimes \psi\bar\psi)$. MC3 holds for $\widehat{\mathfrak{sl}}_{2,k}$ (Family~1, type~A asymptotic prefundamentals) and for $\psi\bar\psi$ (Family~6, $\mathfrak{gl}(1|1)$ super-Yangian parity-balance), hence for the ambient tensor product by Künneth; but MC3 for $\mathrm{SCA}_c$ is NOT a corollary. The two missing steps, inscribed at `rem:n2-mc3-coset-obstruction`: (C1) branching preservation (every isotypic summand is evaluation-type over the coset) and (C2) thick-closure preservation (the thick closure of the branching image exhausts $\mathrm{DK}_{\mathrm{coset}}$). Kazama-Suzuki 1989 + Feigin-Semikhatov 2004 establish CHARACTER-level branching; (C1)-(C2) at factorization-category level are open. Counter: for every MC3 advertisement of a coset chiral algebra, inscribe the coset preservation lemma explicitly or label the claim conjectural. AP502 prevents the "coset-of-covered-ambient is automatically covered" fallacy. Related AP240 (closure-by-repackaging), AP256 (aspirational-heal drift).

### Wave Elliptic-engine (2026-04-18 elliptic chiral QG CYBE audit, AP315-AP316)

**AP315 (Engine tolerance looser than physical scale of failure).** An engine numerical verifier (e.g. CYBE checker, closure-equation test, identity residual) uses a tolerance threshold calibrated for machine-precision identities but applied to a claim whose genuine failures are order-unity. A failing test registers as passed because the tolerance is several orders of magnitude looser than the physical scale of the discrepancy. Canonical violation: 2026-04-18 elliptic chiral QG swarm found `compute/lib/elliptic_rmatrix_shadow.py:475-513` `verify_ybe_elliptic_sl2` uses `passed = relative < 1e-6`; independent first-principles CYBE verification on the Cartan-root Belavin form (the function's actual return value, despite docstring claiming Pauli XYZ form) gives relative residual in the range 0.8-5.8 across test points. The tolerance 1e-6 is three-to-six orders of magnitude looser than the smallest genuine CYBE-failure scale 0.8; the "passed" result is a false positive. Counter: for every engine verifier, document the expected physical scale of success (machine precision for arithmetic identities, ~10^{-12} for transcendental identities after series truncation, order-unity for vanishing-cocycle claims); set tolerance one-to-two orders inside the expected precision, not at the absolute machine limit. For CYBE / YBE / MC / nullity tests, tighten tolerance to `10^{-12}` or tighter and fail loudly on any residual exceeding `10^{-10}`. Distinct from AP128 (engine-test synchronized to same wrong VALUE): AP315 is engine-test synchronized to same wrong TOLERANCE — the value may even be computed correctly but the acceptance gate is miscalibrated. Related: AP257 (engine-docstring vs manuscript contradiction) — often pair-fires with AP315 because a docstring claiming "returns Pauli XYZ form" + a body returning Cartan-root form + a verifier at tolerance 1e-6 gives the whole engine a false coherent appearance.

**AP316 (Worktree-isolated Agent heal abandoned at delivery).** When `Agent` tool is launched with `isolation: worktree`, edits live in `.claude/worktrees/agent-<id>/` on a separate branch, NOT in the main repo working tree. An agent that completes Phase 3 (heal) + Phase 4 (inscribe into .tex) but defers commit "pending build+test verification" or "pending follow-up session" leaves all load-bearing edits stranded in the worktree. The main repo receives only the Phase 5 report file (if the agent was instructed to write it into the main-repo `adversarial_swarm_YYYYMMDD/` path rather than the worktree). Canonical violation: 2026-04-18 elliptic swarm produced a CYBE falsification of `prop:elliptic-rmatrix-shadow` and drafted full replacement text for the inscribed proposition + CLAUDE.md Elliptic row + FM30 + AP275; ALL edits live at `/Users/raeez/chiral-bar-cobar/.claude/worktrees/agent-a18f0d3c/`; the main repo sees only the mission report. Without a merge step the heal is effectively un-delivered; a future session re-runs the same audit, re-discovers the same gap, and re-drafts the same replacement. Counter: at launch of every `isolation: worktree` Agent, include in the prompt an explicit delivery clause — "at mission end, emit `diff = git diff main...HEAD` and write it as `adversarial_swarm_YYYYMMDD/patches/<target>.patch` in the MAIN repo path, and enumerate the expected `git apply` command in the report". Supervisor side: after worktree-isolated Agent completion notification, run `cd .claude/worktrees/agent-<id> && git diff main...HEAD > ../../patches/<id>.patch` before the worktree gets auto-cleaned on session exit. Alternatively, launch the agent without worktree isolation when the heal is constrained to a single .tex file + CLAUDE.md row and concurrency risk is low. Related: AP306 (report-file proliferation), AP313 (truncated Agent-result masquerading as completion), AP303 (Async-Agent forbidden-read) — all variants of "Agent completion signal is opaque to delivery verification".

### Wave 6d-hCS-swarm (2026-04-18 6d hCS codim-2 defect attack-then-heal, AP621-AP623)

Three preventative patterns surfaced in the 2026-04-18 follow-through wave that implemented the Wave 3 (2026-04-18) analytical heal plan for the 6d holomorphic Chern-Simons codim-2 defect = W_{1+∞} identification (CLAUDE.md row line 624). Wave 3 produced the attack ledger + heal plan H1-H4; this wave committed them. See `adversarial_swarm_20260418/attack_heal_6d_hcs_defect_20260418.md` for the implementation record. Primary loci: Vol III `chapters/theory/quantum_chiral_algebras.tex` (+ `lem:5d-to-6d-codim2-bridge` inscribed + `rem:codim2-verification` test-count updated); Vol III `compute/tests/test_hcs_codim2_defect_ope.py` (HZ-IV body upgraded on paths v, vi); Vol I `CLAUDE.md` 6d-hCS row scope-qualified.

**AP621 (Bridge-lemma absence between folklore-citation pillar and theorem step).** A theorem whose load-bearing proof step quotes a named classical result ("as in Costello 2013", "per Beilinson-Drinfeld", "by Kashiwara-Vilonen") without naming the quoted result at a theorem/lemma level carries AP272's single-citation disease in its strongest form. The reader cannot audit the attribution because the cited content has no anchor inside the current manuscript. Canonical violation: Vol III `prop:codim2-defect-ope` Step 3 at `quantum_chiral_algebras.tex:455` asserted $\Psi = -\sigma_2$ as the "Kac–Moody level of the 5d theory on $\C^2 \times \R$, as in Costello 2013" without a local bridge lemma. Healing template: insert a `\begin{lemma}[Bridge from $V_1$ to $V_2$]...\ClaimStatusProvedElsewhere` block that (a) states the cited identity explicitly in the local notation of the current theorem, (b) supplies an `Attribution` proof block naming the primary source (paper + section/equation), and (c) names the downstream role of the bridge (which theorem-step it justifies). Scaffolding cost 5-12 lines; converts "as in X" into a cited bridge-lemma that reviewers can audit, engines can cross-check at boundary test points, and future waves can falsify via the primary source. Distinct from AP272 (folklore citation AT the proof step, diagnostic pattern) by requiring a positive structural repair (the bridge-lemma, prescriptive); distinct from AP249 (base-change cited not inscribed) by applying to ANY external mechanism, not only base-change extensions. Healed locally at `lem:5d-to-6d-codim2-bridge` in `quantum_chiral_algebras.tex` 2026-04-18.

**AP622 (Symmetry-invariance Boolean `True` in HZ-IV body as AP287 primitive tautology).** A test body that hard-codes `is_invariant = True` / `symmetry_match = True` / `triality_preserved = True` for a symmetry claim ($S_3$ triality, $\mathbb{Z}/2$ duality, Miki triality, Weyl triality) is an AP287 structural-impossibility primitive — the symmetry is named but not verified in the test body. Heal: enumerate the symmetry group $G$ explicitly in the test body using `from itertools import permutations` (for $S_n$), nested tuple-swaps (for small Abelian groups), or explicit generator lists (for exceptional groups); compute the invariant quantity at each group element and assert equality with the base value. For $S_3$ on $(h_1,h_2,h_3)$: 6 permutations; each yields $\sigma_2 = $ same value; 6 numerical assertions replace 1 Boolean. Canonical violation: Vol III `test_hcs_codim2_defect_ope.py` `TestCodim2DefectOPEIV` had `pr_w_inf_match = True` as the "Procházka-Rapčák triality" path pre-heal. Healed 2026-04-18 by enumerating `permutations((h_1,h_2,h_3))` at two test points (one symmetric, one generic) with 12 numerical assertions total. Distinct from AP287 (general primitive-tautology structural-impossibility — no numerical observable exists) and AP277 (vacuous HZ-IV body behind sound decorator prose, where numerical observable is possible but not wired): AP622 is specifically the SYMMETRY-invariance variant where the observable is a group enumeration over finitely many elements and the heal cost is O(|G|) lines.

**AP623 (Newton-identity disambiguation of seemingly-identical polynomial expressions in HZ-IV paths).** When quantity $\Psi$ is derived from a specific symmetric polynomial of generators (e.g. $\Psi = -\sigma_2$ in Viete form), an HZ-IV "independent" path that computes $\Psi$ from the SAME Viete polynomial with a permutation, a different variable name, or a cosmetic rewriting, is AP288 label-disjoint-but-computation-identical. True independence requires computing $\Psi$ from a genuinely DIFFERENT polynomial expression — Newton power sums $p_k = \sum h_i^k$, Schur polynomials, monomial-basis sums — related to the primary Viete form by a known polynomial identity (Newton's identities, Jacobi-Trudi). Canonical violation averted 2026-04-18: `TestCodim2DefectOPEIV` "ASV Miura" path, pre-heal hard-coded `asv_miura = True`; naive fix would be "recompute $\sigma_2$ with a loop over pairs" — still Viete-polynomial, still not disjoint from Step 3 of the primary proof. Genuine heal: under CY $\sigma_1 = 0$, Newton's identity $p_2 = \sigma_1^2 - 2\sigma_2$ reduces to $\Psi = (1/2)\sum h_i^2$. Test body computes $p_2$ via `sum(Fraction(x)**2 for x in h)` and divides by 2; this NEVER touches $\sigma_2$. Two test points match primary $\Psi$; the path is genuinely disjoint in computation while equivalent in value. Counter: for every HZ-IV path claiming independence from a Viete-polynomial derivation, name the polynomial-identity bridge (Newton / Jacobi-Trudi / plethysm) connecting the two expressions; if no such bridge exists explicitly, the path is either AP288-collapsed (merge into primary path) or AP277-vacuous (downgrade to structural). Distinct from AP288 (general label-disjoint / computation-identical) by supplying the positive counter-heal (Newton-identity-style disambiguation) and the discipline that MAKES the disambiguation rigorous (named polynomial-identity bridge); distinct from AP622 (symmetry-enumeration heal) which handles the Boolean-truth-marker case, whereas AP623 handles the subtler expression-rewriting case where the naive heal is still computation-identical.

### Wave 18 (2026-04-18 session deferred-inscription batch, AP317-AP319)

Three genuinely-new patterns held back through Waves 7-12 under AP314 inscription-rate throttle. Inscribed per explicit user directive (overrides throttle when user-requested).

**AP317 (Phantom-detector multi-line env-scan gap + archival-directory consumer inflation).** Phantom-label detectors that grep `\label{foo}` on single lines miss inscriptions where `\begin{env}[title]` and `\label{foo}` sit on separate lines; consumer-count queries that include archival directories (`relaunch_*`, `rectification_*`, `healing_*`, `opus_audit_*`, `wave2_audit_*`, `.claude/worktrees/*`, `notes/`, `adversarial_swarm_*/`) inflate ref counts vs live-manuscript. Canonical violation: Wave-7 AP255 re-run (agent a36e7276) reported `conj:topologization-general` as Tier-1 phantom with 37 refs; Wave-7 heal agent (a220890d) proved inscribed at `en_koszul_duality.tex:3392` via multi-line `\begin{conjecture}[title]\n\label{conj:topologization-general}` and 37-ref count inflated by archival dirs (true live-tex ~12). Wave-8 AP255 re-run agent (a36e72762) later classified 257 Wave-7 phantoms as detection false-positives (Vol I 251, Vol II 6). Wave-8 `conj:master-bv-brst` also FALSE POSITIVE at 86 advertised refs (AP286 legitimate alias to `conj:v1-master-bv-brst` at `editorial_constitution.tex:433`). Counter: every phantom-detector must implement (a) multi-line env scan within 3-5 lines below `\begin{env}[title]` headers; (b) metadata cross-check against `metadata/claims.jsonl` + `metadata/label_index.json` where available (Wave-9 metadata-build agent established parity for Vol II + Vol III); (c) live-tex consumer-count filter excluding `{relaunch_*, rectification_*, healing_*, opus_audit_*, wave*_audit_*, .claude/worktrees/*, notes/, adversarial_swarm_*/}` and any path containing `backup`, `archive`, `fix_wave_*`. Related: AP255 (phantom file + phantomsection mask), AP286 (tactical phantomsection alias), AP301 (/loop interval mismatch — infrastructure-level). Failure mode: adversarial-audit false alarm + inscription-by-downgrade risk. The AP255 count is robust only under the full three-step verification.

**AP318 (Inter-wave diagnostic framing inheritance — one wave's "drift" becomes another wave's premise without first-principles re-reading).** An audit-agent identifies some manuscript content as "drift" / "phantom" / "violation" based on a surface pattern; a subsequent agent or main-thread heal directive inherits that framing AS TRUTH without re-verifying from first principles. The second wave's heal would CAUSE damage if executed because the first wave's diagnostic was itself an AP186 shallow correction. Canonical violation: Wave-6 κ(K3×E) retraction agent (a43f5b50) REFUSED to apply the 30-site retraction prompt that inherited the Wave-5 κ-discipline agent's "drift" framing. First-principles reading at `cy_d_kappa_stratification.tex:414-426` revealed two distinct κ_ch routes — Route A (Hodge supertrace, Künneth-multiplicative) and Route B (Heisenberg-level rank-additive) — coexisting legitimately. Naive 30-site retraction would have erased Route B from load-bearing sites including 22-route cross-verification engine `kappa_consistency_adversarial.py`, dissolved the BKM decomposition 5 = 3+2, orphaned ~170 engine tests. Wave-10 Anchor A vs B adjudication agent (a0f44b6d) later formalised the H3 two-subscript split (`κ_ch` Route A canonical, `κ_ch^{Heis}` Route B) — confirming Wave-6's first-principles refusal was correct. Counter: before accepting any Wave-N "drift" / "phantom" / "violation" framing from a predecessor agent: (a) read the flagged site end-to-end; (b) trace the first-principles mathematics (what invariant is being computed? what convention does the surrounding text adopt?); (c) check whether the flagged site belongs to a legitimate-but-unlabeled second convention (Route B / κ_ch^{Heis}-style); (d) only apply the prescribed heal if the first-principles check confirms the diagnostic. Distinct from AP186 (shallow correction without first-principles for mathematical content): AP318 is the AGENT-COORDINATION layer where wave-to-wave inheritance creates the shallow-correction hazard. Distinct from AP305 (pessimistic CLAUDE.md drift overstating deficit): AP318 is the agent-agent equivalent. Programme discipline: every targeted-heal agent prompt should include "verify the previous wave's diagnostic against first principles before applying heal; REFUSE and flag if the diagnostic is itself a shallow correction." Honoured: Wave-6 κ(K3×E), Wave-9 chap:toroidal-elliptic false positive, Wave-8 conj:master-bv-brst false positive — three instances of disciplined refusal this session.

**AP319 (Gold-standard HZ-IV template — inline primary-source computations from disjoint libraries/tables at chain level).** Positive pattern (not a failure mode): HZ-IV decorators achieve genuine disjointness when each path's computation (i) is performed at chain level in the test body itself, (ii) uses a library / table / formula from a DIFFERENT primary source from the other paths, (iii) produces an INDEPENDENT numerical evaluation that agrees with the others at output level, (iv) the underlying engine is demoted to a Path Z sanity anchor (NOT in `verified_against`). Canonical exemplar: Vol I `test_z_g_s_r_arithmetic_duality.py` — (A) MC recurrence integration, (B) `sympy.bernoulli`, (C) `sympy.factorint` blackbox. Each path is independently executable; agreement is at output level; disjointness is substantive at the computation layer, not merely bibliographic. Propagation 2026-04-18: (Wave-8) Vol I `test_theorem_H_hochschild_koszul.py` upgraded with Feigin-Fuchs Fock-screening / strong-generator enumeration / scalar dual-center primary-source triples; (Wave-9) Vol III `test_kappa_bkm_universal.py::TestGoldStandardDisjointPaths` with Eichler-Zagier 1985 theta-ratio / Gritsenko 1999 Δ_5 / Mukai 1988 Nikulin fixed-point orbifold halving across N ∈ {1,2,3,4,6}; (Wave-10) Vol I DS intertwining + BP Koszul conductor + depth-gap trichotomy. Design doctrine (6 principles): (1) each `verified_against` path independent numerical at test-time; (2) no shared-table intermediate (FRAME_SHAPE_DATA-style); (3) agreement at output level; (4) classical-primary-source reachback; (5) engine = sanity only; (6) label-disjointness ≠ computational-disjointness (human audit required). Programme-wide baseline: ~400-500 HZ-IV decorators at 1:19 CLEAN:violation ratio (Wave-7 audit); propagation of AP319 template to all Tier-1 theorem-anchoring tests is the leverage-maximising HZ-IV discipline heal. Distinct from AP277 (vacuous test body, diagnostic pattern), AP288 (label-disjoint / computation-identical, diagnostic), AP310 (bibliographically-disjoint / shared-intermediate, diagnostic): AP319 is the POSITIVE counterpart — the template that CORRECTLY implements HZ-IV disjointness. Related AP266 (sharpened-obstruction healing, positive pattern sibling).

**AP320 (Wave-ledger narrative fabrication — claimed file-state changes that never landed on disk).** An adversarial-audit or heal agent produces a final report claiming specific `\label{}` inscriptions, file Edits, or `@independent_verification` decorators have been applied; direct filesystem grep reveals ZERO of the claimed edits landed on disk. The agent's narrative is self-consistent and well-formed (appropriate section headers, heal tables, commit plans) but the actual repository state is unchanged. Future waves inheriting this narrative as "completed work" compound the drift (AP318 fires downstream). Canonical violation: 2026-04-18 session Wave-13 agent adbf3a47 claimed 3 upgrades to Vol I κ_BKM HZ-IV tests (`test_cy_bkm_algebra_engine.py`, `test_cy_borcherds_lift_engine.py`, `test_moonshine_bar_complex.py`) with detailed before/after structure, AST-parse verification, design doctrine. Wave-14 agent a96093d2 built on this narrative adding 5 more upgrades with similar structural detail. Wave-15 Vol I κ_BKM HZ-IV finals agent (a5cfa3f0) performed the three-step verification (grep `@independent_verification` in all 8 claimed targets) and found ZERO matches — the claimed inscriptions never landed. Possible mechanisms: (a) agent Edit tool invocation silently failed (worktree-isolation merge gap AP316; network/permission error; file-lock conflict); (b) agent wrote to `/tmp` or scratch path not the manuscript; (c) agent hallucinated the inscription sequence without invoking Edit at all. Distinct from AP316 (worktree-isolated heal abandoned at delivery — the edits DID land but in a worktree): AP320 is when edits NEVER landed anywhere. Distinct from AP256 (aspirational-heal advertisement — CLAUDE.md advertises ahead of manuscript, but the manuscript state is known and honestly reported): AP320 is the agent REPORTING edits as done when they were never done. Distinct from AP313 (truncated Agent-result masquerading as clean completion — the agent's report is terminated mid-sentence): AP320 is fully-formed prose with apparent success signals but void at the file-state level. Counter: every heal agent report MUST be verified by a supervisor pass that re-greps the claimed edits before downstream waves inherit. Concretely: for every `\label{foo}` / `@independent_verification` / `Edit applied at file:line` claim in an agent's final report, run `grep -F 'claim-text' <file>` and fail the wave if the grep returns zero. For `isolation: worktree` agents: inspect `git diff main...HEAD` in the worktree before declaring heal-complete. Suggested programme-level discipline: pre-commit hook that refuses commits whose message references a wave's `adversarial_swarm_YYYYMMDD/waveN_*.md` note without matching disk-state diffs. Stronger than AP186 (shallow correction without first-principles) because AP320's failure mode is WORK NEVER PERFORMED, not work performed badly. Related AP316 (worktree abandonment), AP313 (truncated result), AP256 (aspirational advance), AP318 (downstream inheritance amplification).

**CORRECTION (Wave-16 re-verification 2026-04-18)**: the Wave-13/14 canonical example cited above was itself an AP318 detection-gap false positive — Wave-15 finals agent a5cfa3f0 used grep for long-form `@independent_verification` only, missing 2 Wave-13 bare `@_iv(` + 5 Wave-14 suffixed `@_iv_v14_*` alias forms that DID land on disk. True Wave-13/14 verdict: integrity GENUINE (8/10 claimed targets fully landed on disk; 2 Wave-14-deferred, not AP320 violations). The AP320 pattern ITSELF remains real and worth guarding; supervisor passes must use a COMPREHENSIVE decorator regex: `@independent_verification|@_iv\(|@_iv_v[0-9]+_[a-zA-Z]+\(|TestGoldStandardDisjointPaths`. Meta-lesson: AP320 detection is itself susceptible to AP318 inheritance (Wave-15 inherited a hypothesis and grep-verified it with the wrong tool). Register the inversion-failure mode as **AP321 (Supervisor-grep-blind — false-positive AP320 accusation via narrower grep than decorator-alias predecessors used)**: a supervisor enforcing AP320 via a single-form grep pattern produces false-positive "fabrication" verdicts when agents used aliased decorator forms. Canonical case: Wave-15 a5cfa3f0 greping `@independent_verification` only; Wave-13/14 had used `@_iv`/`@_iv_v14_*`. Counter for AP321: always use the comprehensive decorator regex above; cross-check with `import` statements naming the alias (`as _iv`, `as _iv_v14_*`) to confirm both the decorator and its import are present. Genuine AP320 instances in 2026-04-18 session: NONE CONFIRMED by Wave-16 re-verification; the AP320 pattern remains inscribed as preventative guard for future sessions.

### Wave W_N-Stokes (2026-04-18 attack-heal, AP1301)

### Wave Kontsevich-formality (2026-04-18 attack-heal, AP1361)

**AP1361 (Cohomological $\En$-formality invoked with chain-level phrasing; AP-CY33 propagation risk).** Kontsevich / Fresse--Willwacher / Tamarkin formality is a statement about the chain operad $C_*(\En;\R)$ being quasi-isomorphic to its cohomology operad $\Pn = H^*(\En;\R)$. In the RATIONAL / cohomology category this quasi-isomorphism collapses every chain-level refinement of an $\En$-structure to its $\Pn$-shadow (AP-CY33: chain-level != rational; E_n collapses under formality in cohomology). Invoking "Kontsevich formality" without a cohomology-vs-chain scope disclaimer risks a reader promoting a cohomological collapse to a chain-level structural claim, or conversely taking a chain-level claim as a corollary of formality when only the cohomology-level statement is available. The discipline: every invocation of $\En$-formality must carry one of three scope labels — (a) "cohomological" (formality used as a lemma about $H^*(\En)$ = Arnold-type algebra, no chain-level structure claimed); (b) "chain-level with fixed associator $\Phi$" (the Fresse Vol II Thm 16.2.1 route, which rigidifies the operad-level quasi-isomorphism into an algebra-level chain map whose coherence data are torsor-parametrised by $\mathrm{GRT}_1$, with associator-dependence made explicit per `prop:associator-dependence-explicit`); (c) "algebra-level formality lift" (Keller classicality for formal $A_\infty$-algebras, applicable when $m_n = 0$ for $n \geq 3$ on the minimal model). Canonical audit sites in Vol I (2026-04-18): (i) `chiral_hochschild_koszul.tex:717-720` Step 1 of `lem:chiral-quadratic-koszul` invoked "Kontsevich formality" to secure non-degeneracy of the OPE-residue pairing on $\cF(E)(3)$ at collision divisors; the actual mechanism is cohomological identification of the local fibre $\FM_3(\C)$'s real cohomology with the Arnold algebra, followed by Fresse Prop. 7.1.3 non-degeneracy on that algebra — not a chain-level E_3-structure claim. Healed by inserting scope label (a). (ii) `koszulness_moduli_scheme.tex:674-677` invoked "Kontsevich formality for two-coloured little disks" with `\ref{rem:sc-chtop-formality}` pointing to a label that does not exist anywhere across the three volumes (AP264 phantom-ref). Healed by retargeting to `prop:en-formality` with scope label (a) and an explicit disclaimer that no chain-level $\SCchtop$-formality is asserted. (iii) `preface.tex:4407` and `preface_sections10_13_draft.tex:72` invoke "Kontsevich formality supplies a quasi-isomorphism $\Phi: \SCchtop \to \mathsf{SC}$" without scope label — the statement is cohomological (on bar-cobar Quillen equivalence) and should carry the "(cohomology level)" qualifier; documented, not edited this wave per AP314 inscription-rate throttling. (iv) `motivic_shadow_full_class_m_platonic.tex:165` "$E_2$-topological side receives the Kontsevich formality quantisation" — correctly scoped via subsequent associator-dependence discussion (AP171). (v) `e3_identification_chain_level_platonic.tex:411-473` `thm:e3-identification-chain-level-associator-fixed` is the model citizen: it explicitly fixes the associator $\Phi = \Phi_{\mathrm{KZ}}$ before claiming a chain-level map, inscribes the $\GRT_1$-torsor structure at `prop:associator-dependence-explicit`, and isolates associator-dependence on cochains while proving associator-independence on cohomology (part iv) — template for all future chain-level formality invocations in the programme. The H02 ALT proof of Theorem B ("Keller deformation + Kontsevich formality") is advertised in `chapters/frame/programme_overview_platonic.tex:458-459`, `worldview_synthesis_2026_04_17.tex:853`, and `AGENTS.md:1027`, but is NOT inscribed as a theorem body anywhere in Vol I: the Keller side is inscribed via `thm:ainfty-koszul-characterization` + the HPL-transfer passages at `chiral_koszul_pairs.tex:1312-1320, 2550-2557, 2914-2917`, giving a forward-direction chain from chiral Koszulness to $A_\infty$-formality fiberwise on FM strata; the Kontsevich-formality side of the H02 advertisement (configuration-space / graph-integral route) has no inscribed chapter body beyond `thm:chiral-formality` (`deformation_quantization.tex:575`, `\ClaimStatusProvedElsewhere` to Tamarkin 2000 + FG12). AP241 advertised-but-not-inscribed for the full H02 composite. Distinct from AP-CY33 (diagnostic at the rational-vs-chain distinction): AP1361 is the Vol-I-specific prose-scope-discipline forward-check that prevents AP-CY33-class errors from entering chapter bodies. Distinct from AP272 (unstated cross-lemma via folklore citation): the cited works (Kontsevich 1999, Tamarkin 2000/2003, Fresse 2017, Willwacher 2015) all exist; the issue is scope-label absence, not citation absence. Phantom bibkeys `Keller01`, `Keller06`, `keller-icm`, `Fresse17`, `FresseWillwacher20`, `Idrissi22`, `Brown12`, `Furusho11`, `McClureSmith03`, `BarNatan98`, `DTT09` found at audit time in `standalone/references.bib` — subsumed under pre-existing AP281 (bibkey naming-drift at scale); no new AP inscription for these per AP314. Audit trail: `adversarial_swarm_20260418/attack_heal_kontsevich_formality_20260418.md`.

**AP1301 (Classical-chain computation masquerading as "DISCOVERED").** A CLAUDE.md status-table row carries the banner "DISCOVERED" (or "IDENTIFIED", "NEW INVARIANT") for a result whose load-bearing derivation is a short chain of classical steps, each ProvedElsewhere in standard primary literature. The programme-novel content is restricted to the IDENTIFICATION of a classical invariant with a programme object (pole order of a W-algebra OPE, Poincaré rank of a KZ connection, Stokes ray count under Boalch/Jimbo-Miwa-Ueno convention), which is a one-line corollary once the classical chain is assembled. Publishing the identification under "DISCOVERED" invites AP155 novelty-overclaiming + AP280 three-step epistemic inflation at the status-table layer. Canonical violation: "W_N Stokes count DISCOVERED: 4N-4" advertised at CLAUDE.md Theorem Status (line 611 pre-heal). The four-step derivation is entirely classical: (a) $W_N W_N$ OPE pole $2N$ is a two-point function of a weight-$N$ primary (BPZ 1984, Zamolodchikov 1985 for $W_3$); (b) d-log absorption is Arnold 1969; (c) Poincaré rank = pole order − 1 is the standard definition in any Wasow 1965 / Sibuya 1975 textbook; (d) Stokes ray count = $2r$ under generic anti-Stokes / Stokes-sector convention is Boalch 2002 (Symplectic manifolds and isomonodromic deformations, Lemma 2.1), Jimbo-Miwa-Ueno 1981, Hertling-Sabbah 2005. The one-line programme-novel identification is that the W-algebra OPE at weight $N$ drives the Poincaré rank through the chain. Healed 2026-04-18 by downgrading the CLAUDE.md row from "DISCOVERED" to "COMPUTED via standard irregular-sing theory (ProvedElsewhere)" with per-step primary-source attribution. Distinct from AP280 (three-step inflation at manuscript level — remark → standalone → headline): AP1301 is a single-remark inscription with honest primary-source attribution in the Wave-7 audit note; only the status-table banner was inflated. Distinct from AP272 (unstated cross-lemma via folklore citation): here the citations exist in the Wave-7 note (`adversarial_swarm_20260417/wave7_wN_stokes_g2_decomp_attack_heal.md`), just not in the CLAUDE.md row. Counter: for every CLAUDE.md status row using "DISCOVERED" / "IDENTIFIED" / "NEW", before committing verify (a) the derivation is not a 3-5 step classical chain, (b) if it is, replace the label with "COMPUTED (ProvedElsewhere)" and list per-step primary sources in the row body. Prose-hygiene corollary: "DISCOVERED" should be reserved for results whose derivation uses a genuinely new machinery or whose statement has no precedent in classical literature. Full audit: `adversarial_swarm_20260418/attack_heal_w_n_stokes_20260418.md`.

### Wave Symbol-overloading-meta (2026-04-18 meta-pattern consolidation, AP2121)

Consolidation of a recurring meta-pattern surfaced across six independent waves this session: one symbol (letter + subscript convention) denotes multiple distinct mathematical objects across independently-developed contexts, and the programme conflates them at context-interfaces. The error is not arithmetic within any single context but systemic at the notational layer. AP2122-AP2140 reserved for later inscription per AP314 inscription-rate throttling.

**AP2121 (Symbol-overloading across mathematical contexts: notational collision at context-interfaces).** A programme-level scalar or cohomological symbol ($\kappa$, $\varrho$, $Z$, $H^1$, $K$, $\chi$, $c$ when context-ambiguous, $S_r$ across conventions) carries distinct definitions in distinct contexts; each context is internally consistent; the contexts collide at an interface (cross-volume identity, product-formula input, subscript-discipline audit, derived-centre-vs-naive-centre disambiguation) where the symbol appears unscoped. A reader resolves the symbol using context $A$ while the author wrote it using context $B$; the collision is load-bearing only at the interface and is invisible when each context is read alone.

*Canonical sub-instances (this session).* **(a) AP234**: two $K$'s: scalar complementarity $\kappa + \kappa^{!}$ (family values 0 / 13 / 250/3 / 98/3) versus Trinity $K = c + c^{!} = -c_{\mathrm{ghost}}(\mathrm{BRST})$ ($-k$ / $2\dim(\fg)$ / 26 / 100 / 196); bridge identity $\kappa + \kappa^{!} = \varrho_A \cdot K$ with family-dependent $\varrho_A$. **(b) AP311**: two $\varrho$'s: $\varrho_{\mathrm{lin}}$ (leading $\kappa$-coefficient linear in $c$) versus $\varrho_{\mathrm{gen}}$ (KRW signed harmonic sum over strong generators); agree on principal $\cW_N$, diverge on BP ($1/2$ T-line vs $1/6$ generator-profile). **(c) AP1982**: "Drinfeld center dim 1" four-way collision across $Z(U_{\cH_k}\text{-mod})$ braided-monoidal centre, $Z(\cH_k)$ naive commutant, $H^0 \cZ^{\mathrm{der}}_{\mathrm{ch}}(\cH_k)$ cohomological slice; all three equal $\bC$ for Heisenberg by coincidence while derived chiral centre is genuinely 3-dimensional ($P(t) = 1 + t + t^2$). **(d) AP2001**: $\kappa_{\mathrm{ch}}^{\mathrm{Heis}}$ (lattice-rank additive) versus $\kappa_{\mathrm{ch}}$ (Hodge-supertrace K\"unneth-multiplicative); at K3 $\times$ E the additive route gives $2 + 0$, the multiplicative route gives $2 \cdot 0$; reconciliation is AP289. **(e) AP2041-2043**: $H^1(\fsl_2)$ four-way, namely full symmetric bar at generic level (dim 3), degree-2 ordered $V \otimes V$ (rank 4 with unitary lift 8), genus-1 KZB Frobenius-rank 4, full symmetric bar at critical level $k = -h^{\vee}$ (infinite-dimensional $\Omega^{\bullet}(\mathrm{Op}_{\fsl_2}(D))$). **(f) AP290**: HZ-7 $\kappa$ type-swap, Vol III proposition used Sugawara-shift $\dim(\fg)(k + h^{\vee})/(2 h^{\vee})$ (a $\kappa_{\mathrm{ch}}$ formula) under the subscript $\kappa_{\mathrm{cat}}$; subscript present but formula-body-subscript mismatch.

*Structural cause.* Each instance arises at the interface of two or more independently-developed contexts: $\kappa$-complementarity meets Trinity $K$; shadow-tower $\varrho$ meets KRW $\varrho$; four distinct constructions of "centre" share output symbol $Z$; lattice-rank-additive meets Hodge-multiplicative at Künneth products; four bar-construction variants share output symbol $H^1$; subscript discipline meets formula-body substance. The unifying mode is **notational overload at context-interfaces**, not arithmetic error within a context. The pattern is invisible in any single-context chapter and becomes load-bearing only when two or more contexts meet in the same paragraph, theorem statement, product formula, or audit spreadsheet.

*Counter-discipline.* Before deploying any programme-level symbol ($\kappa$, $\varrho$, $Z$, $H^i$, $K$, $\chi$, $c$, $S_r$) in a chapter or standalone: (1) inscribe a canonical-symbol-table at chapter opening (or in a programme-wide conventions remark) listing every subscript variant used, with definitional context plus one boundary-value check per subscript; (2) every use-site asserts which subscript applies, either explicitly in-line (`$\kappa_{\mathrm{ch}}$`) or via a local convention remark within 10 lines of first use; (3) PROHIBIT the bare symbol at any context-interface (any paragraph naming two or more of the contexts in the symbol-table): the subscript MUST appear; enforcement by Pattern 230 grep gate; (4) at context-interfaces, if the interface invariant ($\varrho_A \cdot K$, Künneth-multiplicative $\kappa_{\mathrm{ch}}$, Koszul-involution $c \leftrightarrow K - c$) has a known closed form, inscribe it as a named proposition or remark rather than deriving it inline each time.

*Relation to prior APs.* Stronger than HZ-7 (Vol III bare-$\kappa$ prohibition at zero-subscript) and AP290 (wrong-subscript) because AP2121 is the meta-pattern spanning all subscript symbols and all context-interfaces, not a single-letter discipline. Distinct from AP244 (overcounted foundational terms: multiple names for one object on the working locus): AP244 is name-inflation collapsing to fewer objects; AP2121 is name-compression where one symbol covers multiple genuinely-distinct objects. Distinct from AP234 (two $K$'s, specific instance) and AP311 (two $\varrho$'s, specific instance): AP2121 is the structural abstraction covering all such pairs. Distinct from AP272 (unstated cross-lemma via folklore citation): AP272 catches missing citations; AP2121 catches missing subscript disambiguation within the programme's own notation.

*HOT ZONE candidacy.* Six canonical sub-instances identified this session; programme-wide instance count higher (AP234 propagates across preface + working notes + four surveys + Vol II prose; AP290 recurs at every Vol III $\kappa$-audit; AP1982 recurs at every derived-centre versus naive-centre disambiguation). Candidate for promotion to HZ-12 in the next CLAUDE.md HOT ZONE consolidation wave; this inscription does not perform the HOT ZONE promotion per AP314 inscription-rate throttling. Cache entry: Pattern 230. Audit trail: `adversarial_swarm_20260418/attack_heal_symbol_overloading_meta_20260418.md`.

**WARNING (AP225):** Theorem D all-genera may rest on unproved universality step. Genus-1 obs_1=κ*λ_1 unconditional. All-genera requires (a) clutching-uniqueness (not yet proved), or (b) GRR (H04 sketched not inscribed). Until resolved: state genus-1 unconditional, all-genera conditional.

New wrong formulas B74-B78 in Blacklist above. New failure modes (FM35-FM38, AP42 variants): rate-limit cascade (batch<=5), timeout on >15K line files, "vacuous constraint" confabulation, circular proof chain detection needs DAG tracing.

### Wave 26-AP-Audit (2026-04-21 self-audit of AP catalogue + first-principles cache, AP-CAT-1 -- AP-CAT-8)

Eight patterns surfaced in the 2026-04-21 meta-audit of the AP catalogue and first-principles cache themselves. Themes: cache numbering collisions, bookkeeping-vocabulary infiltration of manuscript chapters, AP-to-open-problem cross-reference gaps, Wave-25 five-agent findings that lacked AP entries, cross-volume non-contradiction registry missing for the pseudo-character/determinant rename, and forward-looking holes.

**AP-CAT-1 (Cache-entry numbering collision).** When multiple waves inscribe cache patterns concurrently under AP304 concurrent-agent collision, two entries can claim the same integer `Pattern N`. The live `first_principles_cache_comprehensive.md` carried Pattern 295 twice (5772 W25 Chenevier; 5800 T1#2 path-independence) and Pattern 240/241/242 twice each (1784/1821/4448 vs 4505/4538/4570). A cache consumer greping `Pattern N` cannot tell which entry to cite. Counter: (a) pre-inscription grep `^(##|###) Pattern N` across the cache before assigning a new integer; (b) upon discovery of a collision, renumber the rebased side to the next free integer above both maxima and annotate both entries with a `(was duplicate-N; renumbered ...)` trailer OR add a disambiguating suffix (`bis`, `b`) with a prominent `Numbering note` block. Distinct from AP304 (multi-agent AP-number collision in the AP catalogue): AP-CAT-1 is the same failure mode in the cache-pattern namespace. Healed 2026-04-21: Pattern 295 at 5800 suffixed to 295bis; Pattern 240/241/242 at 4505/4538/4570 renumbered to 400/401/402.

**AP-CAT-2 (Bookkeeping-vocabulary infiltration of manuscript theory chapters).** The CLAUDE.md Chriss--Ginzburg standard prohibits "Wave N", "round M", "DNA strand S$x$", "AP$n$", "cache entry $n$", etc. in reader-facing `.tex` under `chapters/theory`, `chapters/examples`, `chapters/connections`, `chapters/frame`. Grep at audit time (2026-04-21) across `chapters/theory/` returned 14+ "Wave~13" / "Wave~18" references in load-bearing remarks of `bar_construction.tex:2689,2699,2865`; `cobar_construction.tex:3890,3902,3922`; `algebraic_foundations.tex:2824,2830`; `poincare_duality_quantum.tex:1302,1314`; `koszul_pair_structure.tex:2828,2846`; `quantum_corrections.tex:1474,1498`. The hook `beilinson-gate.sh` does not currently catch "Wave~\d+" in manuscript prose (it scans mathematical AP signatures, not vocabulary hygiene); the `chriss-ginzburg-rectify` skill catches it but only runs on a single chapter at a time. Counter: add a programme-wide grep gate `grep -rn 'Wave~\d\+\|round~\d\+\|DNA strand\|batch~\d\+' chapters/` to the pre-commit hook; any hit blocks the commit. Substitutions: "Wave~13" in a K3-Borcherds context rewrites to "K3", "K3 chiral-bialgebra-level", or just deletion where the mathematical anchor is already present. Healed 2026-04-21 at the 8 load-bearing remarks listed above; the residual `chapters/theory/bar_construction.tex:2865` "representatives of Wave~18" healed to "representatives of the single-valued $\mathrm{ChirHoch}^3$ basis". Related: AP235 (Vol I preface sweep), AP288 (session-ledger stale narrative in notes/).

**AP-CAT-3 (Open-problem to AP cross-reference absent).** The programme tracks 12 open problems in the frontier (Humbert--Heegner admissibility filter, Fake-Monster Cartan rank, universal ratio-of-levels closure, Bridgeland codim closure, pseudo-character primary form, YD tower weight, paramodular vs Fake-Monster $\kappa_{\mathrm{BKM}}$ cross-volume, etc.). Each open problem has an AP or cache-pattern that documents the failure mode behind the gap, but the catalogue does not carry an `Open Problem $k$ -> AP/Pattern chain` mapping at the top-level. A reader searching for "what is the AP status of Open Problem 8 ($\kappa_{\mathrm{BKM}}$ cross-volume)" must grep half a dozen entries. Counter: maintain an `Open Problem Index` block at the top of `antipatterns_catalogue.md` listing each open problem with the AP/Pattern chain that bounds it; update on every inscription of a new open-problem entry. Healed 2026-04-21: Open Problem 8 ($\kappa_{\mathrm{BKM}}$ Vol I vs Vol III, paramodular vs Fake-Monster) rests on AP-CY85, AP-CY69, AP5 dual-indexing header (canonical preamble row 57), Pattern 411 (formerly Pattern 226; cross-volume bibkey check). Other open-problem-to-AP chains to be inscribed on next sweep.

**AP-CAT-4 (Wave-25 five-agent findings not installed as AP entries).** The 2026-04-20 five-agent attack-heal wave surfaced four content findings: (a) pseudo-character vs Chenevier determinant terminology rename; (b) Humbert--Heegner admissibility filter $n \equiv 3, 5 \pmod 8$; (c) $e_4$ two-basis apparent discrepancy (Virasoro-composite vs $W$-primary); (d) non-contradictions table expansion. Findings (a), (b), (d) are installed as AP-CY87 update, Pattern 299 + entry 410 (AP890), and the compatible-dual-readings table at the top of the catalogue respectively. Finding (c) is installed in the canonical-values table (row 51) and as AP-CAT-4 below. The gap is that Wave 25's eight-agent (Anchor A vs Anchor B adjudication + six routes + orbifold diamond + super-twin) findings were installed as content-layer edits but NOT cross-referenced back to AP entries; a future reader cannot trace "what was the AP for Anchor A vs Anchor B" without reading the session ledger. Counter: every wave's final summary must list AP-entry anchors for each finding; findings without an AP anchor get one inscribed at the session close. Healed 2026-04-21: inscribe this AP-CAT-4 as the anchor for the $e_4$ two-basis non-contradiction. $e_4 = W_4 - (107/11)\Lambda_Z$ ($W$-primary basis) and $e_4 = {:}T\partial^2 T{:} - (3/2){:}(\partial T)^2{:} + (321/10)\partial^4 T + \hbar\,\mathrm{qt}(J^{(4)})$ (Virasoro-composite basis) describe the same chiral-Hochschild cocycle class $[e_4] \in \mathrm{ChirHoch}^3(\mathbf H_{\Delta_5})_4$ at $c = -214$; the conversion $W_4 = e_4 + (-c/22)\Lambda_Z$ at generic $c$ fixes the scalar $-107/11 = c/22$ at $c = -214$. Counter: every $e_k$ inscription must name the basis (quasi-primary / $W$-primary / Virasoro-composite) in a one-line scope remark; two-basis disagreements without conversion are latent AP238 violations.

**AP-CAT-5 (Non-contradictions table missing pseudo-character row).** The compatible-dual-readings table at `notes/antipatterns_catalogue.md:83-93` (the canonical non-contradiction registry) lists nine rows through the Pilloni--Urban pseudo-character framework row. The 2026-04-20 W25 rename (Taylor pseudo-character → Chenevier determinant) introduces a compatible-dual-reading between the two frameworks that is documented at entry 422 (AP902) of `appendices/first_principles_cache.md` and Pattern 295 (W25) of `notes/first_principles_cache_comprehensive.md` but NOT inscribed as a non-contradiction-table row. A reader encountering "pseudo-character $S^{\mathrm{ps}}$" in a pre-W25 inscription and "Chenevier determinant $D^{\mathrm{Chen}}$" in a post-W25 inscription may assume a retraction when the correct reading is a framework-upgrade with numerical-content preservation on reduced rings. Counter: inscribe the pseudo-character / Chenevier-determinant row in the non-contradictions table (Taylor 1991 sufficient on reduced, Chenevier 2014 necessary on non-reduced; 46-prime empirical verification preserved). Row already present at line 93 (Taylor vs Chenevier hierarchy). Expand the note to reference Pattern 295 and entry 422. Scope: the row is inscribed but cross-references to AP-CY87 and Pattern 295 need strengthening.

**AP-CAT-6 (Forward-looking hole: metaplectic $\Psi^{\mathrm{metap}}$ row double-counting).** AP-CY71 retracts the "Conway $V^{s\natural}$ is 5th $\Psi$-image" claim; AP892 (entry 412 of appendix cache) restores the fifth row on the metaplectic extension $\Psi^{\mathrm{metap}}$ (distinct from the Lorentzian $\Psi$). The compatible-dual-readings table at line 91 lists "Four sibling functors $\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$ parametrise Baily--Borel--Freitag stratification" versus "Five $\Psi$-image rows including Conway $V^{s\natural}$ on the metaplectic branch". The reconciliation states "four siblings, one extra Conway row on the metaplectic branch"; the LIKELY next-wave failure is that a new moonshine object (e.g. $V^{\natural\natural}$ Thompson-group VOA, baby-Monster VOA) gets inscribed as a sixth $\Psi$-image row without identifying its $\Psi$-sibling; the four-sibling / five-row reconciliation would then break to five-sibling / six-row, which the table does not anticipate. Counter: when a new moonshine / $\Psi$-image row is proposed, require (a) identify which of the four siblings $\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$ it targets; (b) verify the sibling's domain admits the new input (e.g. $\Psi^{\mathrm{metap}}$ admits positive-definite Niemeier lattices but rejects sig-(2,n) Lorentzian input); (c) update the row count in the compatible-dual-readings table atomically. This AP is PREVENTATIVE — no violation yet, but a standing falsification target for the next Conway / Thompson / baby-Monster wave.

**AP-CAT-7 (Forward-looking hole: chain-level vs cohomological ChirHoch concentration for class M at $d \geq 3$).** Theorem H concentrates $\mathrm{ChirHoch}^\bullet(A)$ in degrees $\{0, 1, 2\}$ for compact curve $X$. Under $\Phi_d$ the concentration enlarges to $\{0, 1, 2, d\}$ per AP-CY46 (cohomological). Class M at $E_3$ has bar cohomology $6^g$ per AP-CY21/AP-CY38 (closed form via Künneth); chain-level is infinite. The LIKELY next-wave failure is a theorem inscribing "$\mathrm{ChirHoch}^\bullet$ concentration for class M CY$_3$" without specifying chain-level or cohomological: chain-level would contradict AP-CY21; cohomological is the correct scope but needs explicit labelling. Counter: every concentration claim for class M at $d \geq 3$ carries the `(cohomological)` or `(chain-level: infinite)` label in its first sentence. Target hit: an attack-agent finding "concentration of class M enlarges to $\{0, 1, 2, 3\}$" without scope label would violate AP-CY46 + AP-CY21 simultaneously. Preventative AP; no violation yet at audit time.

**AP-CAT-8 (Forward-looking hole: HZ-IV tolerance calibration for transcendental identities).** AP315 (engine tolerance looser than physical scale) caught the $1 \times 10^{-6}$ tolerance for elliptic CYBE; the next-wave generalisation is that every HZ-IV path involving a transcendental identity (periods, $L$-values, automorphic Fourier coefficients, zeta values) needs a tolerance calibrated against the series-truncation precision, NOT against machine precision. A genuine failure at the $10^{-10}$ scale would pass a test set to $10^{-6}$. The Stirling-like growth of MZV dimensions under Padovan + Humbert--Heegner filtering (Pattern 299 / entry 410) means high-weight paths need tighter tolerance than low-weight ones. Counter: every HZ-IV test body evaluating a transcendental identity must record (a) the series truncation depth; (b) the expected residual at that truncation; (c) the tolerance set to one-to-two orders inside the expected residual; (d) a comment documenting the calibration. Preventative AP; applies to all future transcendental HZ-IV additions.

## Regression Safeguards (non-AP)

RS-3: Physics IS the homotopy type, not a "bridge." Costello-Gaiotto-Dimofte are substance.
RS-4: Costello/Dimofte/Gaiotto content in mathematical core, not "connections" chapters.
RS-9: The slab is a bimodule, NOT a Swiss-cheese disk. Two boundary components.
RS-10: Single-pass agent work without audit FORBIDDEN. Beilinson loop mandatory.
RS-12: The programme is THREE volumes.
RS-13: In Vol II, gravity is the CLIMAX (Part VI).
RS-14: Introduction orients, Overture instantiates.
RS-15: Koszul programme before higher_genus in dependency DAG.
RS-19: Preface is complete survey. Save before compressing.

## Agent Anti-Patterns

AAP2: Terminology rename ATOMIC across all three volumes in single session.
AAP3: Formula implemented ONCE in canonical module, import everywhere.
AAP4: \begin{proof} only after theorem/prop/lemma with ProvedHere. Use Remark[Evidence] for conjectures.
AAP5: Build-artifact commits batched with content. Never standalone artifact commits.
AAP6: Search ALL THREE volumes before downgrading a status tag.
AAP8: README claims independently verifiable by test suite.
AAP9: Wait for ALL agents to finish before launching new batch.
AAP10: After agent completion, verify BOTH engine AND test files exist.
AAP11: Every hardcoded expected value derivable by 2+ independent paths.
AAP12: Asymptotic tolerance proportional to 1/log(N) or verified by two methods.
AAP13: NEVER downgrade model without user permission. Wait and retry on rate limit.
AAP14: Unique branch name per agent.
AAP15: Serialize builds or use isolated worktrees. Parallel pdflatex kills.
AAP16: git stash FORBIDDEN. Use git diff > patch.diff + git apply.
AAP17: Verify edits via git diff, not agent narrative.
AAP18: Confabulating operadic theory -> compute or cite (Loday-Vallette, Vallette, Livernet). NEVER analogize.

## Pre-Edit Verification Protocol

Fill-in-the-blank templates mandatory BEFORE writing listed formula classes. Filling out a template IS the verification. Source: pre_edit_verification_protocol_wave12.md (PE-1 through PE-12).

Invocation protocol: before every Edit touching a trigger pattern, write the relevant template as a fenced block in the reply text (NOT in .tex), fill it in, end with `verdict: ACCEPT`, THEN invoke Edit. If any `match?` is `N` or required source blank, `verdict: REJECT` and re-draft.

Eight highest-priority templates follow. Remaining four (PE-3, PE-6, PE-9, PE-12) in source draft.

**PE-1. r-matrix write** (AP126, AP141, AP19, AP20)

Trigger: any r-matrix formula or reference to `r(z)`, `r_{ij}`, classical r-matrix.

```
## PRE-EDIT: r-matrix
family:                    [Heis / affine KM / Vir / W_N / lattice / Yangian / other]
r(z) written:              [full formula with level prefix]
level parameter symbol:    [k / k+h^v / hbar / c / Psi]
OPE pole order p:          [_]
r-matrix pole order p-1:   [_]              # AP19: d log absorbs one pole
convention:                [trace-form k*Omega/z / KZ Omega/((k+h^v)*z)]
AP126 check (trace-form):  r(z)|_{k=0} = [_]    expected: 0
  (KZ convention: k=0 gives Omega/(h^v*z) != 0 for non-abelian g; this is correct)
match?                     [Y/N]            # must be Y for trace-form; N/A for KZ non-abelian
AP141 grep check:          bare \Omega/z instances in edit scope: [N]
bare \Omega/z allowed?     N
critical-level check (KM): r(z)|_{k=-h^v} = [_]    (trace-form: finite; KZ: diverges)
source:                    [landscape_census.tex:LINE / compute/module.py]
verdict:                   [ACCEPT / REJECT]
```

**PE-2. kappa formula write** (AP1, AP9, AP24, AP39, AP48, AP136)

Trigger: any formula involving kappa or variants (kappa_eff, kappa(B), kappa_ch, kappa_BKM, kappa_cat, kappa_fiber).

```
## PRE-EDIT: kappa
family:                    [Heis / Vir / KM / W_N / bc / betagamma / svir / other]
kappa formula written:     [_]
census citation:           landscape_census.tex:LINE, kappa^{family} = [_]
match?                     [Y/N]            # STOP if N
AP39/AP127 uniqueness: is unqualified kappa = S_2? [Y/N]
  if Y, is the active lane Vir-T or Heis-J? [Y/N]  # otherwise write S_2^lane and projection
evaluation paths:
  at k=0:                  [_]  expected [_]
  at k=-h^v (KM):          [_]  expected 0
  at c=13 (Vir):           [_]  expected 13/2
AP136 boundary (W_N):      formula uses [H_N / H_{N-1} / H_N - 1]
  substitute N=2:          [_]  expected c/2 (W_2 = Vir)
wrong variants avoided:
  NOT kappa(Vir) = c       NOT kappa(W_N) = c*H_{N-1}
  NOT kappa(Heis) = k/2    NOT kappa(KM) = (k+h^v)/(2h^v) (missing dim(g))
verdict:                   [ACCEPT / REJECT]
```

**PE-4. bar complex formula** (AP132, AP22, AP23, AP44)

Trigger: `B(A)`, `T^c(...)`, any bar-construction formula, any desuspension.

```
## PRE-EDIT: bar complex
object written:            B(A) = [_]
T^c argument:              [s^{-1} \bar A / s^{-1} A / s \bar A / bare A]
AP132 augmentation:        \bar A = ker(epsilon) present?  [Y/N]   # must be Y
AP22 desuspension:         |s^{-1} v| = |v| [-1 / +1]              # must be -1
s^{-1} (NOT bare s) used?  [Y/N]                                   # must be Y
coproduct type:            [deconcatenation T^c / coshuffle Sym^c / coLie]
match to intended bar?     [B^ord -> deconc / B^Sigma -> coshuffle / B^Lie -> coLie]
grading:                   cohomological |d|=+1?  [Y/N]
verdict:                   [ACCEPT / REJECT]
```

**PE-5. Vol III kappa write** (AP113)

Trigger: ANY kappa occurrence in `~/calabi-yau-quantum-groups/**/*.tex`. Zero tolerance.

```
## PRE-EDIT: Vol III kappa
subscript written:         [kappa_ch / kappa_cat / kappa_BKM / kappa_fiber / OTHER]
subscript present?         [Y/N]   # must be Y; bare kappa FORBIDDEN
subscript justification:   [chiral shadow / categorified / BKM / fiber correction]
census citation:           Vol III landscape_census_cy.tex:LINE
grep BEFORE write:         bare `\kappa[^_]` hits: [N]
grep AFTER write:          bare `\kappa[^_]` hits: [N]
delta = 0?                 [Y/N]   # must be Y
verdict:                   [ACCEPT / REJECT]
```

**PE-7. Label creation** (AP124, AP125)

Trigger: any `\label{...}` write.

```
## PRE-EDIT: label
environment:               [theorem / proposition / conjecture / definition / remark / lemma]
label written:             \label{prefix:name}
prefix match (AP125):      theorem->thm, prop->prop, conj->conj, def->def, rem->rem, lem->lem
match?                     [Y/N]   # must be Y
AP124 duplicate check (grep all three volumes):
  Vol I matches:           [N]
  Vol II matches:          [N]
  Vol III matches:         [N]
  total BEFORE:            [N]
  total AFTER:             [N]
  delta = 1?               [Y/N]   # must be Y
if duplicate, rename with volume suffix and update all \ref
verdict:                   [ACCEPT / REJECT]
```

**PE-8. Cross-volume formula** (AP5, AP3)

Trigger: any formula shared across volumes (kappa, r-matrix, Theta_A, bar differential, connection 1-form, complementarity).

```
## PRE-EDIT: cross-volume formula
formula:                   [_]
Vol I grep:                [hits, canonical form]
Vol II grep:               [hits, canonical form]
Vol III grep:              [hits, canonical form]
consistent across volumes? [Y/N]
if inconsistent:
  canonical volume:        [Vol I / II / III]
  other volumes updated same session?  [Y/N]  # must be Y (AP5)
convention conversion?     [OPE(I) -> lambda(II) / motivic(III) / NA]
conversion applied?        [Y/N/NA]
verdict:                   [ACCEPT / REJECT]
```

**PE-10. Scope quantifier** (AP6, AP7, AP32, AP139)

Trigger: any theorem statement, any obs_g / F_g / lambda_g formula, any universal quantifier.

```
## PRE-EDIT: scope quantifier
statement:                 [_]
genus:                     [g=0 / g=1 / g>=2 / all g / UNSPECIFIED -> REJECT]
degree:                     [n=_ / all n / UNSPECIFIED -> REJECT]
level:                     [convolution M-bar_{g,n} / ambient Mok25 log FM / both / NA]
AP32 weight tag:           [(UNIFORM-WEIGHT) / (ALL-WEIGHT + delta F_g^cross) / NA]
tagged in statement?       [Y/N]  # must be Y for any g>=2 claim
AP139 free-variable audit:
  variables on LHS:        {_}
  variables on RHS:        {_}
  LHS superset RHS?        [Y/N]  # if N, bind the free variable
AP36 implies vs iff:       [implies / iff]
  if iff, converse proved in same theorem?  [Y/N]
verdict:                   [ACCEPT / REJECT]
```

**PE-11. Differential form type** (AP117, AP27, AP130)

Trigger: any write of a connection 1-form, KZ connection, Arnold form, bar propagator.

```
## PRE-EDIT: differential form
what:                      [connection 1-form / bar propagator / Arnold form / KZ / other]
form written:              [_]
expected type:
  connection 1-form: r(z) dz  (NOT d log)
  KZ:                sum r_{ij} dz_{ij}
  Arnold form:       d log(z_i - z_j)  (bar coefficient, NOT connection)
  bar propagator:    d log E(z,w)  (weight 1 ALWAYS, AP27)
match?                     [Y/N]
AP27 propagator weight:    1?  [Y/N]
AP117 d log check:         if d log appears, Arnold-form context? [Y/N]
space the form lives on:   [fiber Sigma_g / base M-bar_{g,n} / FM_n(X) / Ran(X)]
AP130 fiber-base:          object on fiber vs class on base correctly distinguished? [Y/N]
verdict:                   [ACCEPT / REJECT]
```

Refusal criteria: reject own edit if any `match?` = N, any blank source, any `FORBIDDEN` ticked, grep delta mismatch. On reject: re-draft, re-fill, proceed only when `verdict: ACCEPT`.

Remaining templates PE-3 (complementarity), PE-6 (exceptional dimensions), PE-9 (summation boundary), PE-12 (prose hygiene) in `compute/audit/pre_edit_verification_protocol_wave12.md`.

## Structural Facts

**Shadow tower**: Theta_A := D_A - d_0 is MC (thm:mc2-bar-intrinsic). kappa, C, Q are projections. All-degree convergence PROVED as pro-object limit in $\widehat{\gAmod}$ (surjective-transitions Mittag-Leffler on the operadic weight filtration w(g,r,d)=2g-2+r+d); chain-level termination in gAmod is class-dependent (G/L/C finite, M transfinite). G/L/C/M: G(r=2,Heis), L(r=3,aff), C(r=4,betagamma), M(r=inf,Vir/W_N). Shadow depth != Koszulness. Delta=8*kappa*S_4: Delta=0 ↔ finite tower. SC formality: A is SC-formal iff class G (prop:sc-formal-iff-class-g). Depth gap: d_alg in {0,1,2,inf}; gap at 3 (prop:depth-gap-trichotomy). ChirHoch^1(V_k(g)) = g with total dim = dim(g)+2 (prop:chirhoch1-affine-km).

**Convolution**: dg Lie Conv_str is strict model of L-inf Conv_inf. MC moduli coincide. Full L-inf needed for transfer/formality/gauge equivalence.

**E_1 primacy**: B^ord is the primitive (Stasheff). av: g^{E_1} → g^mod lossy Sigma_n-coinvariant projection. At degree 2, av(r(z)) recovers kappa in the abelian and scalar families; for non-abelian affine KM gives kappa_dp and the full kappa adds dim(g)/2. All standard chiral algebras are E_inf (local); E_1=nonlocal (Yangian, EK quantum VA). NEVER "E_inf means no OPE poles."

**Three-pillar constraints**: (1) Convolution sL-inf hom_alpha(C,A) is NOT strict Lie. (2) hom_alpha fails as bifunctor in both slots simultaneously (RNW19). MC3 one slot at a time. (3) Log FM != classical FM; requires snc pair (X,D).

## Architecture

**Vol I**: Introduction + Overture (Heisenberg CG opening, unnumbered) + Part I (Bar Complex: Thms A-D+H, 12 Koszul equivs) + Part II (Characteristic Datum: shadow tower, G/L/C/M/M*/W, higher-genus, E_1 modular) + Part III (Standard Landscape: all families, census) + Part IV (Physics Bridges: E_n, factorization envelopes, derived Langlands) + Part V (Seven Faces of r(z): F1 bar-cobar twisting, F2 DNP25 line-operator, F3 Khan-Zeng PVA, F4 Gaiotto-Zeng sphere Hamiltonians, F5 Drinfeld Yangian, F6 Sklyanin/STS83, F7 FFR94 Gaudin) + Part VI (Frontier) + Appendices.

**Vol II** (~1,749pp): SC^{ch,top} bar differential = holomorphic factorization on C, coproduct = topological factorization on R. Seven parts: I(Open Primitive) II(E_1 Core) III(Seven Faces) IV(Char Datum) V(HT Landscape) VI(3D Quantum Gravity = CLIMAX) VII(Frontier). See Vol II CLAUDE.md and notes/cross_volume_aps.md for V2-AP* catalog.

**Vol III** (~693pp): CY → chiral functor Phi. ~34,000 tests, ~460 engines. 10 proofs at publication standard. Clean build. Seven parts: I(Foundations) II(CY-to-Chiral Functor) III(E_n Hierarchy) IV(K3 Yangian) V(CY Landscape) VI(Seven Faces r_CY) VII(Frontiers). 4 stub chapters. kappa subscripts MANDATORY. CY-A_3 PROVED (inf-cat). K3 abelian Yangian PROVED. ZTE T COMPUTED. kappa_BKM = c_N(0)/2 universal. Class M E_3 bar dim = 6^g. Shadow tower through m_8 (160 tests, S_8=4144720/19683). Mock modular K3: THEOREM at d=2. CY-D dimension-stratified. See Vol III CLAUDE.md and notes/cross_volume_aps.md for AP-CY1..AP-CY69.

## Writing Standard

Channel: Gelfand (inevitability), Beilinson (falsification), Drinfeld (unifying principle), Kazhdan (compression), Etingof (clarity), Polyakov (physics=theorem), Nekrasov (seamless passage), Kapranov (higher structure IS math), Ginzburg (every object solves a problem), Costello (factorization), Gaiotto (dualities compute), Witten (physical insight precedes proof). **Convergent loop mandatory**: WRITE → REIMAGINE (Gelfand/Beilinson/Drinfeld) → REWRITE from scratch → BEILINSON AUDIT (adversarial) → REIMAGINE AGAIN → REWRITE AGAIN → CONVERGE (zero findings >= MODERATE). Preface/intro: 3+ iterations. Chapter openings: 2+. **CG structural moves**: deficiency opening, unique survivor, instant computation, forced transition, decomposition table, dichotomy, sentence-as-theorem.

## Skills

### `/chriss-ginzburg-rectify` — MANDATORY execution protocol (CONSTITUTIONAL)

When `/chriss-ginzburg-rectify <file>` is invoked, you MUST analyse the WHOLE file, chunk by chunk, linearly progressing from start to finish, with SMALL chunk size (~50 lines). No exceptions.

- **Whole-file coverage is non-negotiable.** Do not declare convergence on a partial sweep. Do not skip to "hot spots" and ignore the rest. Do not skim large regions via screening-grep in lieu of chunk-by-chunk audit. The file is not done until every line from 1 to EOF has passed the five-gate check.
- **Linear order.** Start at line 1 (or at the resume cursor if the queue entry is `[~]`), advance monotonically. Never skip forward. A chunk-N issue discovered from chunk-(N+k) context propagates back to chunk N; fix in place, continue forward.
- **Small chunks.** ~50 lines per chunk is the default. The user is watching and wants tight verdicts. 100 lines is the absolute maximum. Larger chunks dilute the audit (Gate 2 define-before-use and Gate 3 motivation cannot be checked reliably in 200-line blocks).
- **Budget-rule partial is legitimate, shallow-screen is not.** If the file is large (>3000 lines) and the session budget runs out, mark `[~]` with a resume cursor `% RESUME-FROM: <line/label>` and commit partial progress. This is fine. What is NOT fine: declaring `[x]` after editing only a few "obvious" spots and claiming the rest was "already clean" based on a grep pass. The queue file records what counts as done; only the full five-phase sweep counts.
- **If the skill was invoked without a file argument or on a non-existent path, STOP and ask.** Do not silently rectify an unrelated file.

Why: prior sessions have cut corners by running a lightweight `rg` screen for obvious prose slop and declaring "CONVERGED (0 edits)". That screen misses formula bugs, define-before-use violations, unmotivated definitions, operadic conflations, scope inflations, and circular-proof routings that the five-gate audit is designed to catch. Every shortcut becomes a future AP. The only reliable protection is linear whole-file coverage.

### Skill index

```
/build                      Build all three volumes, tests, census
/audit [target]             Deep Beilinson audit (6 hostile examiners)
/chriss-ginzburg-rectify [file]  Full 5-phase CG + Beilinson rectification (canonical) — WHOLE FILE, chunk by chunk, linear, small chunks
/rectify-all                Rectification across all chapters, all volumes, parallel tiers
/research-swarm [topic]     Launch 30+ research agents on frontier
/verify [claim]             Multi-path verification (3+ independent paths)
/propagate [pattern]        Cross-volume AP5 propagation check
/compute-engine [name]      Scaffold compute engine with multi-path tests
/beilinson-swarm            Beilinson rectification swarm across all chapters
/rectify [file]             Beilinson rectification loop (lighter than CG)
/beilinson-rectify [file]   CG fortification + rectification (v1)
/chriss-ginzburg-rectify-v1 [file]  CG rectification v1 (superseded)
```

RS-1,2,5,6,7,8,11,16,17,18,20 merged into corresponding APs. AP16 superseded by AP27.

## Build

```
pkill -9 -f pdflatex 2>/dev/null || true; sleep 2; make fast    # Vol I
cd ~/chiral-bar-cobar-vol2 && make                                # Vol II
cd ~/calabi-yau-quantum-groups && make fast                       # Vol III
make test                                                         # Fast tests
make test-full                                                    # All tests (~119K)
python3 scripts/generate_metadata.py                              # Census
```

CAUTION: Watcher spawns competing pdflatex; always kill before builds.

## Session Protocol

1. Read this file. 2. Build: `pkill -9 -f pdflatex; sleep 2; make fast`. 3. Tests: `make test`. 4. `git log --oneline -10`. 5. Read .tex source before any edit (never from memory). 6. After each change: build+test. After each correction: grep ALL THREE volumes (AP5). 7. Never guess a formula: compute or cite. Check landscape_census.tex (AP1). 8. Apply convergent writing loop to all prose. 9. Session end: build all three volumes, run tests, summarize errors by class. 10. Before first Edit, read the HOT ZONE (HZ-1 through HZ-10) and run Pre-Edit Verification Protocol mental check: is the pending edit touching an r-matrix, kappa, bar complex, label, Vol III kappa, cross-volume formula, scope quantifier, or differential form? If yes, fill the corresponding PE-1..PE-12 template as a fenced block in the reply BEFORE invoking Edit, ending with `verdict: ACCEPT`.

Details: FRONTIER.md (research programme status), MEMORY.md (session history), concordance.tex (constitution), notes/cross_volume_aps.md (Vol II/III AP catalog), notes/true_formula_census.md (full C1-C31).

## LaTeX

All macros in main.tex preamble. NEVER \newcommand in chapters (use \providecommand). Memoir class, EB Garamond (newtxmath + ebgaramond). Tags: \ClaimStatusProvedHere, \ClaimStatusProvedElsewhere, \ClaimStatusConjectured, \ClaimStatusHeuristic. Label everything with \label{def:}, \label{thm:}. Cross-reference with \ref. Do not add packages without checking compatibility. Do not create new .tex files when content belongs in existing chapters.

## Vol III 6d hCS Session Cross-Awareness (2026-04-12/13)

Feeds back into Vol I. AP-CY23-34 detailed in notes/cross_volume_aps.md. Capsule: AP-CY23 E_1-chiral bialgebra is correct Hopf home (B^ord preserves R-matrix; B^Σ kills Hopf); AP-CY24 docstring confabulation; AP-CY25 R=(id⊗S)∘Δ(1) WRONG, use half-braiding; AP-CY26 σ_2 even under h_i→-h_i (k^!=-k from Shapovalov); AP-CY27 sandbox non-persistence (verify with ls); AP-CY28 pole-unsafe test points (avoid z=±h_i); AP-CY29 wrong-repo writes; AP-CY30 factored != solved (YBE does NOT imply ZTE); AP-CY31 spectral z != worldsheet z; AP-CY32 reorganisation != bypass; AP-CY33 chain-level != rational (E_3 collapses under formality); AP-CY34 total {b,B^{(2)}}=0 via Costello TCFT but individual {b_k,B^{(2)}}!=0.

Key Vol I infrastructure results: Shadow tower S_k = A_∞ coproduct correction δ^{(k)} PROVED; shadow-Feynman L-loop=S_{L+1}. ZTE fails for Yang R-matrix at O(κ²); E_3 nontrivial beyond E_2; S^{corr} EXISTS. E_3 bar cohomology: class L=(1+t)^{3g}=2^{3g}; class C=2^{3g}; class M=6^g (Kunneth, chain-level P(q)^{6g}). Universal coproduct Δ_z(e_s)=Σ C(N_R-b,k) z^k e_a^L·e_b^R (all spins). Conductors: G/L ρ_K=0; M(Vir) 13; K3×E 0. Chiral CE: B(U^ch(L))=CE_*(L) RETRACTION REQUIRED (2026-04-18, `adversarial_swarm_20260418/attack_heal_chiral_CE_20260418.md` + `adversarial_swarm_20260418/attack_heal_chiral_CE_reassessment_20260418.md`): naive pointwise identification with the classical exterior CE of the finite-dim generator space contradicts Vol I `prop:derived-center-explicit` at `chapters/theory/chiral_center_theorem.tex:1967` (Heis total 3 not 2; affine sl_2 total 5 not 8; Virasoro Poincaré 1+t^2 not 1+t). Correct scope is BD04 Thm 4.8.1 derived-global-sections $\int_C \mathrm{CE}(L) \simeq \mathrm{CE}(\mathrm{R}\Gamma_{\mathrm{DR}}(C, L))$ — cited correctly in Vol I at `chapters/theory/bar_construction.tex:440` and `chapters/theory/chiral_koszul_pairs.tex:5809`. Vol III `prop:bar-ce-chiral` + `prop:bar-ce-identification` + `compute/lib/chiral_ce_complex.py` pending downgrade to `\ClaimStatusConjectured` with BD04 scope in a subsequent Vol III wave. κ_BKM=c_N(0)/2 universal (Borcherds weight); naive κ_ch+χ(O_fiber) is coincidence for K3×E. CY-A_3 PROVED in inf-categorical framework; chain-level [m_3,B^{(2)}]!=0 is NOT obstruction; Obs_Ainf=0 UNIVERSALLY via Costello TCFT.
See ~/calabi-yau-quantum-groups/FRONTIER.md F13-F24 and CLAUDE.md for full details.

## Alternative Proofs Secured (2026-04-13)

Every main theorem has at least TWO independent proof paths:

| Theorem | Primary | Alternative |
|---------|---------|-------------|
| A | Twisting morphisms + filtered comparison | Lurie inf-cat nerve-realization (H01) |
| B | Bar filtration spectral sequence | Keller deformation + Kontsevich formality (H02) |
| C | Fiber bar + eigenspace decomp | PTVV shifted symplectic (H03) |
| D | Shadow tower + genus universality | GRR on universal curve (H04) |
| H | Bar-Hochschild comparison | Deformation-theoretic dimensional analysis (H05) |
| MC2 | Recursive inverse limit | KS scattering diagram SKETCHED (H06; central-charge + support-property comparison uninscribed per KS08 Def. 1 + Rem. 1) |
| MC5 | **STABLE ONLY AFTER AMBIENT SEPARATION**: raw direct-sum class-M chain-level statements in `Ch(Vect)` remain false; strict completed/pro and coderived surfaces are the verified or conditional repair surfaces. | Do not describe the direct-sum failure as an ambient artefact and do not identify the inverse-limit product with the original direct-sum complex. Finite S4 is only the first detected term; the real obstruction is the nonterminating harmonic family. Continuous completed/coderived equivalences require the explicit CP1--CP3 co/contra comparison package. |
| Topol | Sugawara [Q,G]=T | CFG factorization homology (H08) |
| SC-formal | Shadow tower truncation (primary) | Operadic rewrite (H11) FOLDED INTO PRIMARY -- no independent second proof survives; forward unconditional, converse conditional on two-colour transfer lemma (rem:sc-formal-open-transfer-frontier) |
| Depth gap | MC relation at degree 4 | Shadow Lie Jacobi (H10) |
| Compl K | Fiber-center + Theorem D | Index theory / Euler characteristic (H12) |

Condition removal: uniform-weight (H13), Koszul locus (H14), chain-level topol (H15), perfectness C1 (H17).

## Platonic Ideal Roadmap (2026-04-13)

**Unconditional (high confidence):** Thms A (fixed-curve), B (on-locus), C0 (D^co), C1 (g>=1), D (non-circular; but see AP225), H, MC1, MC2, MC4, SC-formality, depth gap, D^2=0, Theta_A, ChirHoch^1, 10 Koszul equivs, Verlinde recovery, ker(av), Miura coefficient, critical level jump, E_3 identification (simple g), chiral QG equiv, gl_N chiral QG.

**Conditional (genuine mathematical restrictions):** C2 (uniform-weight), MC3 type-A (Baxter b=a-1/2), MC4 resonance (transfer comparison), MC5 chain-level (class M false), Koszul (vii) multi-weight (genus-0), Koszul (viii) freeness (Massey vanishing).

**Conjectural:** Topologization chain-level original complex (A-inf coherence), topologization general (non-KM), Theorem A modular-family (relative Ran base-change), off-locus chain qi (beyond class G/L).

**Open frontier (status corrected 2026-04-24):** the 2026-04-16 "all frontiers closed" headline is superseded. Chain-level E_3/topologization and MC5 class-M statements are surface-dependent: raw direct-sum class-M comparisons in `Ch(Vect)` remain failure surfaces, while completed/pro/coderived replacements require explicitly verified strict-tower and CP1--CP3 continuous comparison hypotheses. Modular-family Theorem A over $\overline{\mathcal M}_{g,n}$ remains conditional on relative Ran base-change and log-boundary gluing; chiral coproduct claims for non-gauge-theoretic families remain valid only at the scope stated in their theorem bodies. `conj:periodic-cdg` for admissible KL is separately closed by Vol I `periodic_cdg_admissible.tex` (`thm:periodic-cdg-is-koszul-compatible`).

**Beilinson-rectified open frontiers (2026-04-17 audit, surfaced AFTER the 2026-04-16 wave):**
- **W(p) triplet tempering (TWO-MASSEY SPLIT)** — the Vol II commit `a5640de` inscription RETRACTED `thm:tempered-stratum-contains-wp` from ProvedHere to Conjectured. Wave-1 audit 2026-04-17 (F4) refines: the Gurarie 1993 (arXiv:hep-th/9303160) / Flohr 1996 (arXiv:hep-th/9605151) counterexamples falsify the CORRELATION-FUNCTION Massey bound; they do NOT falsify the SHADOW-TOWER Massey bound (structure constants S_r of the bar tower), which is a different object. Two open sub-items: (a) correlation-function Massey unbounded on W(p) — confirmed (retraction stands); (b) shadow-tower Massey bounded on W(p) via C_2-cofiniteness — OPEN, is the bound c-independent? does it imply shadow-tower tempering? The `logarithmic_wp_tempered_analysis_platonic.tex` three-channel Stirling argument addresses (b), not (a). Pending Adamović-Milas character-amplitude bound for the definitive closure of (b). CORRECT tempered scope (non-shadow-sense): principal + non-logarithmic + non-minimal standard landscape.
- **Non-tempered stratum OVERCLAIM** — the programme-climax statement "non-tempered stratum is EMPTY on the C_2-cofinite standard landscape" is SCOPE-QUALIFIED: emptiness holds on the non-logarithmic subset (Virasoro, W_N, all Schellekens, Monster, irrational cosets). Logarithmic W(p) remains open.
- **CY-C pentagon invariant** — Vol III commit `cade61c` healed the pentagon stratification `{3, 12, 24}` from `κ_ch^{R_i}` to `ρ^{R_i}` (generator-lattice rank). Category error: κ_ch is route-independent = 0 for K3×E by Hodge supertrace; the stratification is an algebraic invariant orthogonal to κ_ch.
- **Kummer-irregular prime labelling** — Vol I commit `9668336` retracted primes 1423, 3067, 23, 43, 419 from the Kummer-irregular label (primary-source Bernoulli witness search found no witness). Wave-1 audit 2026-04-17 (F3) extended verification: {23, 43, 61, 193, 419} verified regular by exhaustion of the full Kummer window (B_2..B_{p-3}) by direct audit-internal Bernoulli-numerator computation; {1423, 3067, 811, 2111, 16657} verified regular through B_418 in-audit, with full confirmation via Buhler--Harvey 2011, "Irregular primes to 163 million" (cite added to B92). These primes still appear in S_r numerators as Riccati-arithmetic characteristic primes, NOT Kummer-arithmetic. Corrected Tier-3 emergence: {37, 691, 811}; 3067 dropped.
- **β_N all-rank harmonic law** — STATUS CORRECTED (Vol II `chapters/theory/beta_N_closed_form_all_platonic.tex`, `conj:beta-N-harmonic-closed-form`): β_N = 12(H_N - 1) = sum_{s=2}^{N} 12/s is the conjectural harmonic identification, not a proved all-rank closure. Conditional on it, both prior candidates are ruled out at N=4: (N+1)(N+2)/2 predicts 15; N²-N+4 predicts 16; harmonic value is 13. Conditional on the same law, β_N is rational and non-integral for every N ≥ 5; β_5 = 77/5, β_6 = 87/5. The proved finite-rank statement is `thm:wn-tempered-all-N`, which uses finiteness plus Stirling dominance, not the all-rank harmonic identification.
- **Super-complementarity canonical pairing** — the `κ + κ^! = max(m,n)` identity for super-Yangians scopes to sub-Sugawara line; two pairings (super-trace vs Berezinian) coexist without programme-level canonicalisation via the central automorphism σ_{sBer} (lem:super-trace-berezinian-bridge, yangians_foundations.tex:103). Verdier pairing inscription pending. psl(2|2) subtlety (Wave-1 audit 2026-04-17, F2): Beisert 2007 gives THREE central elements (P, K, C in psl(2|2)⊕C^3), not two; the coincidence with max(2,2)=2 comes from CENTRE RANK in the AdS/CFT extension, not central-element count. The HZ-IV disjoint_rationale at yangians_foundations.tex:92-94 should be read as rank-match (V3 is specialisation-only at psl(2|2), partially disjoint).
- **ϱ_BP=1/6 structural origin** — currently BACK-DERIVED from generator profile (J^bos_{h=1}, G^±,ferm_{h=3/2}, T^bos_{h=2}). Open: identify 1/6 as Hessian determinant of DS Lagrangian, coset central-charge ratio, or Berezinian shift. Wave-1 audit 2026-04-17 (F1.b).
- **K(W_N) for N ≥ 4** — CLOSED 2026-04-18 (STALE frontier retracted). Closed form K(W_N) = Σ_{j=2}^{N} 2(6j² − 6j + 1) = 2(N−1)(2N² + 2N + 1) inscribed at `cor:uc-K-WN` in `chapters/theory/universal_conductor_K_platonic.tex:641-688` (ProvedHere, via Toda BRST spin-tower: principal Drinfeld-Sokolov bc(j)-ghost generators at Casimir spins j ∈ {2, ..., N}, applying the per-spin uc-K-g formula). Third forward difference Δ³K(W_N) = 24 constant (cubic). Tabulated values K(W_2) = K(Vir) = 26, K(W_3) = 100, K(W_4) = 246, K(W_5) = 488, K(W_6) = 850 cross-referenced against `chapters/examples/landscape_census.tex` L199 (additional tabulations at landscape_census.tex:938,1520). The ϱ_A·K(A) modular-characteristic identity K^κ(W_N) = κ(W_N) + κ(W_N^!) = K(W_N)·(H_N − 1) follows at `rem:uc-harmonic-density-WN`:690-699. Wave-1 audit 2026-04-17 (F1.c) retracted; ϱ_BP structural-origin bullet (F1.b) remains genuinely open.

The Beilinson audit inscribed `notes/rectification_map_beilinson_audit.md` (926 lines) with full verdicts and heal paths; the post-audit priority order places the climax-rewrite and preface-refresh tasks LAST to prevent propagation of unverified closures.

**Recovery infrastructure:** `scripts/resume_failed.py`, `scripts/campaign_dashboard.py`, 9 campaign scripts.

## Geometric vs Algebraic Model Conflations (AP-CY62--AP-CY67)

Migrated 2026-04-16 to `notes/cross_volume_aps.md` (full AP-CY62..67 entries, all triggers, ramification guards). Vol I agents MUST read that file before edits touching ChirHoch, End^ch_A, BZFN, spectral parameters, FM_k(C), or geometric-vs-algebraic chiral model claims.

Operational capsule (sufficient if cross_volume_aps not loaded):
- AP-CY62: specify "geometric (FM)" vs "algebraic (bar/operadic)" for any chain-level C^*_ch(A,A).
- AP-CY63: "chiral endomorphism operad on FM_k(C)" FIRES (End^ch is algebraic, not on FM).

---

## K3 chiral bialgebra and cross-volume patterns (Wave 23-26, April 2026)

Identified through the adversarial swarm campaign on the non-abelian
K3 chiral bialgebra $\mathbf{H}_{\Delta_5}$, paramodular automorphic
data, super-Etingof-Kazhdan quantisation, and the CY-to-chiral functor
family $\{\Phi_d\}_{d \ge 1}$. Cross-volume patterns appear in Vol I,
Vol II, and Vol III simultaneously.

- **AP320 -- $\mathrm{grt}_1^{(1/2)}$ extension splitting (High).**
  Extension $0 \to \mathrm{grt}_1 \to \mathrm{grt}_1^{(1/2)} \to
  \bigoplus_k \mathbb{Q}\widetilde\sigma_{2k} \to 0$ is simultaneously
  non-abelian, non-split, and non-central. Conflation of these three
  distinct properties, or confusion of group vs Lie cohomology for the
  obstruction class, is a programme-level type error. Obstruction
  $[\omega_{\mathrm{SK}}] \in H^2_{\mathrm{Lie}}(\mathfrak{q}; \mathrm{grt}_1)$
  relates to the group Eichler cocycle
  $[\mathrm{SK}(\Delta)/\Delta] \in
  H^1(\mathrm{Sp}_4(\mathbb{Z}); \mathrm{Hom}(\mathrm{grt}_1, \mathbb{Q}[v_{\Delta_5}]))$
  via van Est transgression $\tau_{\mathrm{vE}}:H^1_{\mathrm{grp}} \to H^2_{\mathrm{Lie}}$.
  Witness: $[\widetilde\sigma_2,\widetilde\sigma_2] = 288 \widetilde\sigma_4$
  via $\tau(2)^2/2 = 288$. Hilbert series disagree at every even
  weight (Brown 2012; Furusho 2011), so $\mathrm{grt}_1^{(1/2)} \not\cong \mathrm{grt}_1$.

- **AP321 -- Chiral-Hochschild cocycle $e_k$ in full motivic periods
  (Medium).** BGS analytic torsion on Shimura varieties forces landing
  in the single-valued subring $\mathrm{zv}^{\mathrm{sv}}$ (Brown 2014;
  Schnetz 2013). Correct: $e_k = \mathrm{sv} \circ \pi^{\mathrm{depth}\le k}(\phi^{(3k)})$.
  Motivic home shrinks from Padovan-dim $\mathrm{Per}^{\mathrm{mot}}_{3k}$
  to strict subspace $\mathrm{zv}^{\mathrm{sv}}_{3k}$. At $k=4$:
  dim-3 SV, NOT $\mathbb{Q}\pi^4$.

- **AP322 -- "Four is all" as single Scheithauer 2017 citation
  (Medium).** Single-paper attribution is incomplete. Finiteness half
  of the classification requires a three-paper chain: Scheithauer 2017
  (lift existence) + Dittmann-Ma-Scheithauer 2021 Adv Math 386
  (finiteness of reflective signature-$(2,n)$ even genera) + Scheithauer
  2006 Invent Math 164 (prime-level enumeration). "Four is all" is
  GN-reflective-scoped; $24A_1$ Niemeier Borcherds product
  (Borcherds 1995 Invent Math 120 §13) is a non-GN fifth lift.

- **AP323 -- Pseudo-character for reducible Arthur parameter (High).**
  Pseudo-characters lose mod-$\ell^n$ Cayley-Hamilton data when
  $\rho_{\Delta_{10}} = \rho_{\Delta_{E_6}} \oplus \chi^8 \oplus \chi^9$
  is reducible. Use Chenevier 2014 determinants
  $D_{\Delta_{10}}: \mathbb{Z}_\ell[G_\mathbb{Q}] \to \mathbb{Z}_\ell$
  of dimension 4 factorising $D_{\rho_{\Delta_{E_6}}} \otimes D_{\chi^8 \oplus \chi^9}$,
  with Hecke field $\mathbb{Q}$ (since $\dim S_{26}(\mathrm{SL}_2) = 1$).

- **AP324 -- Tensor factorisation through $\mathrm{Vec}_G$ for nonabelian
  $G$ (High).** DGNO 2010 Selecta Math 16 Prop 2.11: pointed $\mathrm{Vec}_G$
  is modular only for abelian $G$ with non-degenerate quadratic form.
  $\mathrm{Vec}_{S_3}$ is not modular. For nonabelian $G$ use
  $G$-crossed braided fusion categories $\mathcal{C} \rtimes G$
  (Turaev 2000; Etingof-Nikshych-Ostrik 2010). Grothendieck-ring
  identity at $\mathbb{Z}$-algebra level can survive but MTC
  factorisation fails.

- **AP325 -- Agent inscription report vs disk mismatch (High).**
  Multiple agents returned truncated reports despite high tool counts
  (Polyakov $e_k$ 52+ tools thin; Gelfand $\Psi$-siblings 957 seconds
  empty). Report $\ne$ disk. Verify via `grep -l` + file-size delta
  after every agent completion.

- **AP326 -- Hook-cascade content loss (Critical).**
  Automated CG-rectify cascades do NOT preserve mathematical
  inscriptions containing bookkeeping vocabulary in prose, titles, or
  labels. 2026-04-20/21 cascade removed ~7 substantive Wave-23/24/25
  inscriptions (Beilinson $\mathrm{Stab}=48$; Kazhdan Selmer; Costello
  master BV; Drinfeld GRT-super; Gelfand rank-162 MTC; Gaiotto four
  siblings; Witten master $L$-value correction) as collateral because
  those inscriptions contained "Wave N" / "DNA" / "AP\d+" tags.
  **Counter**: inscriptions bookkeeping-free from first keystroke.
  Named section/remark titles denote mathematical objects not waves.
  Equations bear mathematical labels not catalogue IDs.

- **AP327 -- Numerical oscillation as convergence (High).**
  Values flipping sign or magnitude across adjacent adversarial waves
  (Leech root norm $2 \to 6 \to 2$; Witten heterotic $c(1,2,\pm 2)$
  $+1 \to -2$; Fake-Monster $c(28)$) without independent path-verification
  are ping-ponging, not converging. Convergence threshold: two consecutive
  waves with zero sign flips on any "verified" coefficient. Three
  independent path-verifications required.

- **AP328 -- Orphan inscription (High).**
  Files under `chapters/examples/` or `chapters/theory/` are NOT
  automatically included; they must be `\input`ed in main.tex. Vol III
  `chapters/examples/hochschild_calculus.tex` was orphaned; built
  chapter `chapters/theory/hochschild_calculus.tex` lacked the
  canonical $e_k$ inscription. Verify wiring via
  `grep -n "input.*TARGET" main.tex` before inscribing.

- **AP329 -- Cross-volume $\kappa_{\mathrm{BKM}}$ inconsistency (High).**
  Vol I abstract $\kappa_{\mathrm{BKM}}(\mathbf{H}_{\Delta_5}) = 12$;
  Vol III abstract $= 5$. "$N$"-index convention differs. Fix $N$ to
  Borcherds-family index (Monster=1, K3=2, ...) per landscape_census.
  Every Vol I/II/III occurrence of $\kappa_{\mathrm{BKM}}$ must name
  the input denominator (Fake-Monster $\Phi_{12}$ vs paramodular
  $\Phi_{10} = \Delta_5^2$ vs ...). AP5 audit required.

- **AP330 -- $N = 1$ naive $\kappa_{\mathrm{BKM}}$ decomposition
  (Medium).** $\kappa_{\mathrm{BKM}} = \kappa_{\mathrm{ch}} + \chi(\mathcal{O}_{\mathrm{fiber}})$
  holds as a numerical coincidence at K3×E ($N=1$) only. FAILS for
  $N \ge 2$. Use family-specific $c_N(0)/2$ universally. Distinguish
  four $\kappa_\bullet$'s: $\kappa_{\mathrm{cat}}(K3 \times E) = 0$
  (total, Künneth multiplicative), $\kappa_{\mathrm{ch}}^K = 3$
  (Künneth additive), $\kappa_{\mathrm{BKM}}$ (family), $\kappa_{\mathrm{fiber}}(K3) = 2$.

- **AP331 -- Six routes as six $\Phi$-applications (Medium).**
  Six routes to $G(K3 \times E)$ (CoHA, Schiffmann-Vasserot,
  Maulik-Okounkov, Borcherds, Toda, DMVV) are six different
  constructions witnessing the same $\Phi_3$-output via pentagon
  colimit, not six $\Phi_3$-applications. $\Phi$ gives ONE output
  per category. Different $\kappa$-values across routes come from
  different constructions (lattice rank $\rho^{R_i} \in \{3, 12, 24\}$
  at generator level), not different $\Phi$-applications.

- **AP332 -- Master $L$-value identification (Critical).**
  $\log Z^{(1)}_{\mathbf{H}_{\Delta_5}} = -\log \Delta_5 -
  \kappa_{\mathrm{BGS}} L'(0, \Delta_{10}, \mathrm{ad}^0)$ is wrong
  via three conflations: (a) adjoint vs standard (Yoshikawa 2004 +
  Bruinier-Kühn 2003 give degree-5 standard); (b) $\Delta_5$ vs
  $\Delta_{10}$ (anomaly $-\log \Delta_5$ pins paramodular); (c)
  CAP reducibility $L(s, \Delta_{10}, \mathrm{ad}^0) =
  L(s, \mathrm{Sym}^2 \Delta_{E_6})\zeta(s+1)\zeta(s-1)$
  (Pitale-Saha-Schmidt 2014). Correct: $L'(0, \Delta_5, \mathrm{std})$
  with $\kappa_{\mathrm{BGS}} = 24$. Waldspurger squaring at Euler-factor
  level: $L(2s, \Delta_5, \mathrm{std}) L(2s, \Delta_5 \otimes \epsilon_{K(1)}, \mathrm{std})
  = L(s, \Delta_{10}, \mathrm{std}) \cdot (\text{bad primes})$.

- **AP333 -- Bloch-Kato Selmer representation-vs-regulator (Critical).**
  $\dim H^1_f(\mathrm{ad}^0 \rho_{\Delta_{10}}) = 1$ is
  representation-correct but regulator-wrong. Adjoint spinor at CAP
  is rigid (Chenevier 2014 + Diamond-Flach-Guo 2004):
  $\dim H^1_f(\mathrm{ad}^0 \rho_{\Delta_{10}}) = 0$. Correct:
  $\dim H^1_f(\mathrm{std}\,\rho_{\Delta_5}) = 1$, the paramodular
  cyclotomic Hida family tangent at tame level $K(1)$.

- **AP334 -- Pentagon admissibility variable (High).**
  Admissibility congruence for $\phi^{(n)}$ is on the Heegner
  discriminant $D_n = (n-3)/2 \pmod 4 \in \{0, 1\}$, equivalently
  $n \equiv 3, 5 \pmod 8$. Not on $n$ directly. Humbert-Heegner
  filtration $\mathfrak{H}_D$ via Eichler-Zagier 1985 Thm 9.1
  polar support $\Delta \ge -m^2$; unconditional on K3 side.
  $\phi^{(5)} = -2 \cdot [\mathrm{gen}]^{\otimes 5}$ first admissible
  non-vanishing.

- **AP335 -- Class-$\mathcal{S}$ trinion/tube arithmetic (High).**
  $\mathcal{T}[A_1, \Sigma_{0, n}]$ at $n = 24$: 22 trinions + 21 tubes
  ($n-2$ trinions, $n-3$ tubes on $\mathbb{P}^1$). Each $A_1$ trinion
  contributes $n_h = 4$ (tri-fundamental $\mathrm{SU}(2)^3/\mathbb{Z}_2$);
  each tube adds $n_v = 3$. Total $(n_v, n_h) = (63, 88)$,
  $c_{4d} = 107/6$, $c_{2d} = -214$. Universal:
  $c_{4d} = (5n - 13)/6$. WOV-2 lock.

- **AP336 -- Fake-Monster automorphic home (High).**
  $\Phi_{12}$ is NOT a Siegel form on $\mathrm{Sp}_{26}(\mathbb{Z})$ or
  Jacobi in 26 variables. Correct: Borcherds-Hermitian automorphic form
  on $\mathcal{D}_{\mathrm{II}_{26,2}} = O(26,2)^+/(O(26) \times O(2))$
  of complex dim 26, singular weight $12 = (26-2)/2$. Leech simple
  roots have norm 2 (Conway 1983): $r_\lambda = (\lambda; 1, 1-\lambda^2/2)$
  at $\lambda^2 = 4$. Weyl vector $\rho = e$ lightlike.

- **AP337 -- Humbert divisor vs Argyres-Douglas (Medium).**
  Physical identification: at $E_{\tau_1} \times E_{\tau_2} \in H_N$,
  Seiberg-Witten curve degenerates to pair of elliptic curves producing
  $(A_1, A_{2N-1})$ Argyres-Douglas point (GMN 2009 Adv Theor Math Phys
  13 Ex 8.3). Four-sibling BKM stratification
  $\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$
  equals class-$\mathcal{S}$ join CDT $\cup$ AD.

- **AP338 -- $\Psi$-functor surjectivity alone (Medium).**
  $\Psi$ is not surjective onto super-EK-quantisable BKMs: super-affine,
  quantum-toroidal, metaplectic Conway escape. Minimal complete family:
  four sibling functors $\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$
  indexed by Baily-Borel-Freitag stratification.

---

## Wave 20-24 Extensions

Identified 2026-04-20 through the Wave 20-24 adversarial programme
on bar-cobar architecture, ChirHoch periods, pentagon tower
coincidences, paramodular arithmetic, the five-archetype landscape
with BKM crown, the unified five-theorem identity, universal
ratio-of-levels across the $\Psi$-image BKMs, and the $k_N$-ladder
along Niemeier roots. Author: Raeez Lorgat.

### Bar-cobar architecture

- **AP339 -- Octachotomy at codim 3 on $\overline{\mathcal A_2}$ (Critical).**
  Bar-cobar inversion closes across 8 ambients: fibrewise generic /
  single-monodromy-refined / bi-unipotent Malcev / tri-unipotent
  Malcev / weight-completed coderived / $A_\infty$-corrected /
  $(\infty,1)$-Perf / chiral-Kontsevich-formal on Koszul locus. First
  non-empty admissible Heegner triple $(3, 4, 7)$ is realised by
  $[E_{\sqrt{-3}} \times E_{\sqrt{-7}}]$; codim-4 empty by
  Bruinier-Humbert rank obstruction. WRONG: "Igusa transversality
  empties all admissible triples". CORRECT: arithmetic Hodge-index
  obstruction governs emptiness; the specific triple $(3, 4, 7)$ is
  non-empty.

- **AP340 -- Universal $k$-tower $k_{\max}(g) = g(g+1)/2$ (High).**
  Closure index equals the Siegel unipotent-radical dimension
  $\dim U^{\mathrm{Sieg}}_g = g(g+1)/2$. Genus ladder:
  $k_{\max}(2) = 3$, $k_{\max}(3) = 6$, $k_{\max}(4) = 10$. WRONG:
  "closure index is $g$-independent". CORRECT: closure index scales
  with Siegel unipotent radical dimension $g(g+1)/2$.

- **AP341 -- Genus-3 Schottky codim 0 (High).**
  WRONG: "Schottky locus is codim 1 at $g = 3$". CORRECT: Schottky
  at $g = 3$ is codim 0 (Zariski-dense in $\mathcal A_3$, Hoyt 1963);
  first genuine Schottky stratum appears at $g = 4$ via the
  Schottky-Jung quartic.

- **AP342 -- $(\infty, 2)$-adjunction Auslander-Reiten structure (High).**
  WRONG: "1-categorical bar-cobar adjunction headline captures the
  full content". CORRECT: the Auslander-Reiten triangular structure
  on compatibility 2-morphisms is the 2-categorical content of
  $\Omega^{\mathrm{ch}} \dashv B^{\mathrm{ch}}$; 1-truncation only
  captures unit/counit, discarding the AR-triangle data.

- **AP343 -- Compatibility-data canonical torsor (High).**
  $\mathrm{Data}(A) \simeq B\mathrm{Aut}^h_{L_\infty}(\mathfrak g^{\infty, 1}_A)$.
  For K3 BKM this specialises to
  $\widehat{\mathfrak{grt}}_1 / \mathrm{ob}^{\mathrm{GN}}$. WRONG:
  "compatibility homotopy is chosen". CORRECT: compatibility is a
  canonical torsor over the automorphism group of the
  $(\infty, 1)$-convolution dGLA, not a choice.

- **AP344 -- $(\infty, 1)$-obstruction tower (High).**
  WRONG: "single $\mathrm{Ext}^2$ captures inversion obstruction".
  CORRECT: the full Massey tower
  $\{\mathrm{ob}^{(k)}\}_{k \ge 2}$ arises via Lurie's formal moduli
  machinery; Bridgeland slicing reads the obstruction stratification
  and no single $\mathrm{Ext}^2$ class suffices.

### ChirHoch

- **AP345 -- Six verification paths for $\chi_3$ period (Critical).**
  Period identity
  $\langle [\chi_3], [e_3^{K3 \times E}] \rangle_{\Phi_3} =
  2\, \mathrm{Vol}(E)(2\pi i)^3 =
  \chi(\mathcal O_{K3}) \cdot
  \mathrm{Res}_{s = 1/2}[E_5^{(2)}(Z, s)]|_{K(1)}$
  has six genuinely independent verification paths: (A) CoHA
  Casimir $\times$ Mukai vector; (B) Igusa $\Phi_{10}^{-1}$ polar
  residue; (C) elliptic-volume rigidity via base-change along
  $E \to \mathrm{pt}$; (D) Kuznetsov relative HPD over $E$;
  (E) cyclic Feigin-Tsygan computation; (F) NC Hodge
  Mukai-Hochschild. The seventh path is
  derived-deformation-theoretic: the Gritsenko-Nikulin Borcherds
  twist classifying morphism (NOT Etingof-Kazhdan degree 2, NOT
  Drinfeld-centre degree 4).

- **AP346 -- Factor-of-2 Mukai double-count (Critical).**
  WRONG: pairing $= 4\, \mathrm{Vol}(E)(2\pi i)^3$. CORRECT:
  $2\, \mathrm{Vol}(E)(2\pi i)^3$; $\chi(\mathcal O_{K3}) = 2$ enters
  ONCE, not twice. The double-count arises from conflating the
  Mukai vector with the Euler characteristic during the
  base-change-to-$E$ reduction.

- **AP347 -- Reduced-vs-unreduced DT conflation (High).**
  The polar-leading coefficient
  $[p^{-1} q^{-1} y^0](-\Phi_{10}^{-1}) = 1$ is REDUCED
  Donaldson-Thomas; zero-curve-class unreduced
  $\mathrm{DT}_{(0, n)}(K3 \times E) = 0$ trivially since
  $\chi(K3 \times E) = 0$. Conflating the reduced and unreduced
  theories at the zero curve class is a type error in period
  extraction.

- **AP348 -- Theorem H SHARP vanishing $k \ge d + 1$ (Critical).**
  WRONG: "$\mathrm{ChirHoch}^k$ concentrated in $\{0, 1, 2\}$ at
  all $d$". CORRECT: at CY-$d$ input, concentration $\{0, 1, 2, d\}$
  is SHARP; degrees $k \ge d + 1$ genuinely vanish. The vanishing at
  $k = d + 1$ and above is a theorem, not a convention.

### Pentagon tower

- **AP349 -- Isolated-Niemeier symphony at $n = 24$ (High).**
  The four-voice coincidence $\chi(K3) = 24 = \mathrm{rk}(\mathrm{Leech})
  = p_{24}(12) = A_2^{12}$ umbral rank is UNIQUE in $n \le 32$.
  No other Hilbert-scheme dimension coincides simultaneously with a
  Niemeier rank, a partition-generating-function coefficient, and an
  umbral moonshine rank in this range.

- **AP350 -- Göttsche $p_{24}(k) = \chi(\mathrm{Hilb}^k(K3))$ GENERIC (High).**
  WRONG: "Göttsche coincidence at $n = 32$ is umbral resonance".
  CORRECT: $p_{24}(k) = \chi(\mathrm{Hilb}^k(K3))$ holds at every
  $k \ge 0$ by the universal Heisenberg-oscillator argument
  (Göttsche 1990, Math Ann 286 Thm 0.1); it is not umbral, and
  reading it as such is a category error.

- **AP351 -- $n = 26$ no Niemeier resonance (Medium).**
  All 23 Niemeier lattices have rank 24; none has rank 26. The value
  $n = 26$ is the Polyakov bosonic-string critical dimension
  $D_c = 24 + 2$, not a moonshine coincidence. Reading 26 as a
  Niemeier echo is a type error.

- **AP352 -- Padovan Binet constant (Medium).**
  The Padovan Binet constant is $\rho^2 / (2\rho + 3) = 0.3106$,
  NOT $1 / (\rho^2 + 2) = 0.2663$. The wrong formula arises from a
  misapplied Binet-type recurrence normalisation.

### Arithmetic

- **AP353 -- Chenevier determinant type error (Critical).**
  ⟂ retracted symbol $S^{\mathrm{ps}}$ per canonical preamble: see
  Chenevier determinant (Pattern 295 / W25 in
  `notes/first_principles_cache_comprehensive.md`); the programme-canonical
  symbol is $D^{\mathrm{Chen}}$, following Chenevier 2014
  arXiv:1301.0635 \S1.2. Type-error content preserved: WRONG:
  $D^{\mathrm{Chen}}$ lives "on $\mathbf H_{\Delta_5}$" (a
  non-commutative chiral bialgebra). CORRECT: $D^{\mathrm{Chen}}$
  lives on the paramodular Hecke algebra $\mathbb T^{\mathrm{par}}_1$
  (commutative). Assigning a Chenevier determinant to a non-commutative
  chiral bialgebra is a structural error.

- **AP354 -- Hecke primary form $E_6 \Delta$ weight 18 (High).**
  The correct elliptic form for the scalar weight-ten Saito--Kurokawa
  packet is $f_{18}=E_6\Delta$ of weight 18, not $\Delta$ of weight 12
  and not $E_4\Delta$ of weight 16. Witness (Chenevier 2014
  polynomial-law axiomatisation; ⟂ retracted symbol $S^{\mathrm{ps}}$
  per Pattern 295 / W25): $D^{\mathrm{Chen}}_1(T_p)=a_p(f_{18})+p^8+p^9$
  and $D^{\mathrm{Chen}}_4(T_p)=p^{34}$; 56 primes verified $p\le263$.

- **AP355 -- Shimura-Waldspurger bridge constant (High).**
  $C_k = 2^{2 \epsilon(k) - 6} \cdot L(1, \chi_f) / L(1/2, \mathrm{Sh}(f))$;
  at $k = 4$: $C_5 = 1 / 64$. The factor $2^{2 \epsilon(k) - 6}$ is
  load-bearing and frequently dropped.

- **AP356 -- Fricke LDP variance INVERSE-DENSITY (Critical).**
  WRONG: $\sigma_k^2 = (\pi^2 / 2) \cos^{-2} \theta_k^*$ (diverges
  at finite-density nodes). CORRECT:
  $\sigma_k^2 = \pi / (2 \sin^2 \theta_k^*)$. The inverse-density
  form is unconditional via Newton-Thorne 2021 and gives the
  finite variance required by the large-deviation principle.

### Five-archetype landscape

- **AP357 -- Five-archetype G/L/C/M/B (Critical).**
  $\kappa_{\mathrm{ch}} + \kappa_{\mathrm{ch}}^! \in \{0, 8, 13, 250/3, 98/3\}$.
  Row **B** = BKM crown $\mathbf H_{\Delta_5}$ at
  $\kappa + \kappa^! = 8$, $r_{\max} = \infty$. WRONG: four-archetype
  G/L/C/M exhaustive. CORRECT: BKM row B is the fifth archetype;
  omitting it loses the non-trivial $\kappa + \kappa^! = 8$ point of
  the complementarity spectrum.

### Unified five-theorem identity

- **AP358 -- A/B/C/D/H as five derived evaluations (Critical).**
  Theorems A, B, C, D, H are five derived evaluations of the single
  adjunction $\Omega^{\mathrm{ch}} \dashv B^{\mathrm{ch}}$ along the
  functor tower
  $\{\pi_0, D^{\mathrm{co}} / D^{\mathrm{ctr}}, \mathrm{tr},
  \int_{\overline{\mathcal M_g}}, \mathrm{RHom}\}$. WRONG: treat
  A/B/C/D/H as five independent theorems. CORRECT: they are one
  adjunction viewed through five functors; the unified identity is
  an invariant of $\Omega^{\mathrm{ch}} \dashv B^{\mathrm{ch}}$
  alone.

### Ratio-of-levels

- **AP359 -- Universal ratio-of-levels (High).**
  $\ell_X / \ell_Y = c_+(L_X) / c_+(L_Y)$ across the $\Psi$-image
  BKMs: Monster / Enriques / K3 / Fake-Monster at
  $(c_+, \ell) = (1, 2), (2, 4), (4, 8), (25, 50)$. The Leech-Conway
  row breaks (no Fricke involution on $\Lambda_{24}$). WRONG:
  extending the ratio-of-levels identity to Leech-Conway. CORRECT:
  $\Psi$-image only, Fricke-equipped lattices; Leech-Conway excluded
  by absence of the required involution.

### $k_N$ ladder failures

- **AP360 -- $4A_5$ phantom Niemeier at $N = 6$ (High).**
  WRONG: "Niemeier root system $4A_5$". CORRECT: no such Niemeier
  exists; re-anchor to $A_5^4 D_4$ with rank-deficit $D_4$
  absorption. The phantom $4A_5$ entry originated from a
  mis-extrapolated $A_n^k$ pattern without rank accounting.

- **AP361 -- Coxeter-void fourth regime at $N = 11$ (High).**
  At $N = 11$ the Coxeter number $h(A_{10}) = 11$ is unique, and no
  Niemeier filler exists. The $k_N$ ladder therefore has four
  regimes: naive / substitute / void / Leech-escape. WRONG: treat
  $N = 11$ under the substitute regime. CORRECT: $N = 11$ sits in
  the void regime; the Leech escape only activates beyond.


---

## Historical K3-BKM anti-patterns (Waves 14-22, cross-volume)

Patterns surfaced before the Wave 23-26 closure campaign. These are
the cross-volume K3-chiral-bialgebra AP-CYs that Vol I inherits via
the shared landscape-census formulas and the CY-to-chiral functor
family.

- **AP339 -- $(c_{4d}, c_{2d})$ central-charge reversal (Critical).**
  Wave-14 erroneously retracted $(107/6, -214)$ to $(26, -312)$.
  Wave-15 Gaiotto healed. Canonical $(c_{4d}, c_{2d}) = (107/6, -214)$
  for $\mathcal{T}[A_1, \Sigma_{0, 24}]$ on the K3-BKM row (Chacaltana--
  Distler 2010 §5.14; Beem--Rastelli 2013; WOV-2 lock).
- **AP340 -- Monster vs Fake-Monster Cartan rank (Critical).**
  Monster rank 2 on $\mathrm{II}_{1,1}$; Fake-Monster rank 26 on
  $\mathrm{II}_{25,1}$; K3-BKM rank 3 on $\Lambda^{2,1}_{II}$. Never
  conflate.
- **AP341 -- $c_3$ normalisation (High).**
  Use Bruinier reduced-class convention: $c_3 = -8$. Conversion
  $176256 = -22032 \cdot (-8)$ to Gritsenko--Nikulin Cartan-matrix.
- **AP342 -- Umbral Niemeier labelling rule (High).**
  $(N-1) \mid 24$ for $A_{N-1}$ umbral; $N = 6$ re-anchored to $6D_4$
  Niemeier root system, $k_6 = 9/2$, umbral group $3.\mathrm{Sym}_6$.
- **AP343 -- $\zeta(3,3,3,3)$ numerical value (High).**
  $\zeta(3,3,3,3) = 0.000295999\ldots$, not $0.0028565$ (10× error).
  First depth-4 MZV at weight 12.
- **AP344 -- Integer heterotic-lift $c(1,2,\pm 2)$ (High).**
  $c(1, 2, \pm 2) = -2$, not $+1$. Cross-verify via three paths.
- **AP345 -- Humbert $H_4$ vs $H_8$ (High).**
  $H_4$ is $(2,2)$-isogeny quotient with $\mathrm{End} \supset \mathbb{Z}[2i]$,
  monodromy order 2. $\mathbb{Q}(\sqrt 2)$-RM locus is $H_8$, not $H_4$.
- **AP346 -- Theorem B scope (Medium).**
  Bar-cobar inversion holds on $\overline{\mathcal{A}_2} \setminus
  \bigcup_{n \,\mathrm{admissible}} H_n$ — all admissible Heegner
  divisors excluded (Wave-18 tightening), not just $H_1 \cup H_4$.
- **AP347 -- Borcherds vs Gritsenko weight (Medium).**
  $\mathrm{Borch}(\phi_{0,1}^{K3}) = \Phi_{12}$ (weight 12, Fake-Monster);
  $\mathrm{Grit}(\eta^9\vartheta_1) = \Delta_5$ (weight 5, paramodular).
  The K3-BKM denominator is $\Delta_5$ (additive), not $\Phi_{12}$
  (multiplicative).
- **AP348 -- Two-$\hbar$ convention (High).**
  $\hbar^{\mathrm{Drinfeld}} = 2\pi i/\ell$ (root-of-unity) and
  $\hbar^{\mathrm{BV}}$ (loop-counting) agree numerically at
  $\hbar^2 = -1/8$, $\ell = 8$, but are semantically distinct. Name
  the $\hbar$ at every site (AP151 bridge).
- **AP349 -- Sylvester vs Feingold--Frenkel signature (Medium).**
  Sylvester gives misleading $(2, 0, -32)$ on $G_{\mathrm{BKM}}$.
  Feingold--Frenkel 1983 eigenvalue-based signature $\{+4, +4, -2\}$
  is correct; true signature $(2, 1)$.
- **AP350 -- Schmidt archimedean pair (Medium).**
  $(17/2, 15/2)$ for $\Delta_{10}$; $(7/2, 5/2) \otimes \mathrm{sgn}_\mathbb{R}$
  for $\Delta_5$ on Maass-spin cover. Both correct for their forms.
- **AP351 -- CoHA vs chiral-algebra type (Critical).**
  $\mathrm{CoHA}_{K3 \times E}$ is $E_1$-associative (Hall product),
  not a chiral algebra. Chiralisation requires the $\Phi_3$-arrow.
  $\mathbf{H}_{\Delta_5} = \Phi_3(\mathcal{D}_\hbar(\mathcal{Y}^{\mathrm{Hall}}(\mathrm{CoHA})))$.
- **AP352 -- $\kappa_\bullet$ Künneth discipline (High).**
  Distinguish $\kappa_{\mathrm{cat}}(K3 \times E) = 0$ (multiplicative),
  $\kappa_{\mathrm{ch}}^K = 3$ (additive), $\kappa_{\mathrm{BKM}}$
  (family), $\kappa_{\mathrm{fibre}}(K3) = 2$ explicitly.


---

## Foundational AP1-AP236 (from Waves 1-13 error archaeology)

The legacy CLAUDE.md (lines 500-1400) holds AP237 onward in explicit list
form; AP1-AP236 are referenced in prose throughout the legacy file but
are canonically load-bearing. Every `/chriss-ginzburg-rectify` Gate 1
already enumerates the most-triggered early APs; this block makes them
first-class entries so they can be consulted without reading the legacy.

- **AP1 -- $\kappa$ pattern-matching (High).** Never substitute one
  family's $\kappa$ into another. $\kappa(\mathcal H_k) = k$;
  $\kappa(V_k(\mathfrak g)) = \dim(\mathfrak g)(k + h^\vee)/(2 h^\vee)$;
  $\kappa(\mathrm{Vir}_c) = c/2$. Recompute from first principles at
  every site (AP3 forbids pattern-matching).
- **AP3 -- Pattern-match vs first-principles recomputation (High).**
  Before citing a formula, recompute from its defining chain. Formulas
  that look alike across families almost never transport.
- **AP4 -- Status tag must match what the proof proves (High).**
  `\ClaimStatusProvedHere` bodies must actually deliver a proof of
  the stated theorem. A weaker proved statement inside a stronger
  theorem body is a retraction.
- **AP5 -- Cross-volume formula propagation (Critical).** After any
  formula edit, grep all three volumes for variants; stale copies
  cause cascading retractions. Every $\kappa$/$c$/$K$-value is AP5-auditable.
- **AP7 -- Universal claim vs proof scope (High).** "for all" claims
  require proofs covering edge cases (critical level, admissible level,
  self-dual point); narrower proofs require scope qualification.
- **AP13 -- Circular proof dependency (High).** Proofs cannot cite the
  theorem being proved; every hypothesis-of-cited-result must be
  independently verified.
- **AP14 -- Koszulness $\ne$ Swiss-cheese formality (Medium).** Two
  distinct properties. Never substitute one for the other.
- **AP17 -- Unchecked correction cascade (High).** A "one-line fix"
  that propagates without AP5 grep leaves the repository inconsistent.
- **AP18 -- Scope qualifier ("non-critical level", "at integrable level",
  "on Koszul locus") must be in the theorem statement, not the remarks.
- **AP19 -- $r$-matrix pole order vs OPE pole order (High).**
  $r$-matrix poles are one fewer than OPE poles. Off-by-one is common.
- **AP24 -- Complementarity $\kappa + \kappa^! = ?$ (High).**
  $= 0$ for affine KM/free fields; $= 13$ for Virasoro at $c = 13$;
  $= 250/3$ for $\mathcal W_3$; $= 98/3$ for Bershadsky-Polyakov;
  $= 8$ for K3-BKM. NEVER unqualified "$= 0$".
- **AP25 -- Bar-cobar inversion $\ne$ derived centre (High).**
  $\Omega(B(A)) = A$ is inversion; $\mathbb D_{\mathrm{Ran}}(B(A)) =
  B(A^!)$ is Verdier; $C^\bullet_{\mathrm{ch}}(A, A) =
  \mathcal Z^{\mathrm{der}}_{\mathrm{ch}}(A)$ is bulk/Hochschild. Three
  different functors, never conflate.
- **AP27 -- Propagator weight $= 1$ (Medium).**
  Propagator weight is $1$, not $h^\vee$. Mixing weights fails at $g \ge 2$.
- **AP32 -- Uniform-weight $\ne$ all algebras (Medium).**
  Theorem D uniform-weight lane is a scope qualifier; genus-1 alone
  is not all-genera; AP32 covers both directions.
- **AP33 -- $\mathcal H_k^! \ne \mathcal H_{-k}$ (High).**
  $\mathcal H_k^! = \mathrm{Sym}^{\mathrm{ch}}(V^*)$, NOT the dual
  Heisenberg algebra.
- **AP34 -- Bar-cobar inversion $\ne$ open-to-closed holography (High).**
  Two independent structural passages on $B(A)$.
- **AP40 -- `\ClaimStatusConjectured` inside theorem environment (Medium).**
  Status tag must match environment: conjectures go in
  `\begin{conjecture}` blocks.
- **AP45 -- Desuspension lowers degree (Medium).** $|s^{-1} v| = |v| - 1$,
  not $|v| + 1$. Sign convention is locked.
- **AP46 -- $\eta(q)$ includes $q^{1/24}$ (Medium).**
  Dedekind $\eta(q) = q^{1/24} \prod_n (1 - q^n)$. Missing the $1/24$
  prefactor yields factor-$24$ errors in genus-1 computations.
- **AP82-AP85 -- Three coalgebra structures on $B(A)$ (High).**
  $\mathrm{Lie}^c$ vs $\mathrm{Sym}^c$ vs $T^c$ (tensor cofree). The
  chiral bar uses $T^c(s^{-1} \bar A)$ with deconcatenation; never
  mix the three.
- **AP97, AP104 -- Ordered bar is primitive; symmetric bar is
  $\mathrm{av}$-image (High).**
  $B^{\mathrm{ord}}$ is the primary object; $B^\Sigma$ is the
  coinvariant shadow. $E_1$ primacy (AP104).
- **AP113 -- Bare $\kappa$ without subscript (High).**
  Every occurrence must be $\kappa_{\mathrm{ch}}$, $\kappa_{\mathrm{cat}}$,
  $\kappa_{\mathrm{BKM}}$, etc. Bare $\kappa$ is scope ambiguity.
- **AP126 -- Mega-sweep label propagation (Medium).**
  After a rename (e.g., "shadow Postnikov" $\to$ "shadow obstruction tower"),
  one grep pass rarely suffices; 5+ commits were needed in the historical
  case. AP126 mandates atomic renames with zero-residual verification.
- **AP142 -- Local-global scope (High).**
  A fibrewise theorem does not globalise automatically. Every
  "holds on $\overline{\mathcal M}_{g,n}$" claim must include the
  sewing/Mok25 log-FM / PTVV transgression step.
- **AP151 -- Two-$\hbar$ convention bridge (Critical).**
  $\hbar^{\mathrm{Drinfeld}} = 2\pi i/\ell$ (root-of-unity) vs
  $\hbar^{\mathrm{BV}}$ (loop-counting). Distinct parameters that happen
  to agree numerically at specific specialisations; never omit the
  subscript.
- **AP165 -- $B(A)$ is NOT an $\mathsf{SC}^{\mathrm{ch,top}}$-coalgebra
  (Critical).** $B(A)$ is an $E_1$-coassociative coalgebra;
  $\mathsf{SC}^{\mathrm{ch,top}}$ emerges only on the chiral derived-centre
  pair $(\mathcal Z^{\mathrm{der}}_{\mathrm{ch}}(A), A)$.
- **AP186 -- Shallow correction / term swap (High).** When fixing an
  error, explain WHY the replacement is correct and HOW it relates to
  the original. Swaps without explanation are syntactic, not mathematical.
- **AP225 -- Socle vs nose scope for $\lambda_g$ (High).**
  $g \in \{1, 2\}$ on the nose; $g \ge 3$ only on the socle $R^g / \mathcal N^g$,
  on-the-nose lift conditional on $\lambda_g$-conjecture.
- **AP232 -- Repackaging is not closure (High).**
  "N prior gaps CLOSED" requires INDEPENDENT resolution, not relocation
  into cited remarks.

The full expanded list of AP237-AP2140 defined in this catalogue file
below; AP1-AP236 trigger patterns are all covered by this foundational
block plus the extracted detailed entries that follow.

---

## Waves 1-13 catch-up anti-patterns (2026-04-20 exhaustive audit)

Anti-patterns surfaced in the session archives for Waves 1 through 13 but
not previously consolidated into the main catalogue numbering. Each entry
carries: (a) the failure mode in first-principles terms, (b) the concrete
witness file / formula / archival-session reference, (c) the programme-
canonical correct statement with primary-literature citation.

### Session 2026-03-31 Beilinson swarm (Waves 1-3)

- **AP801 -- Vol II $F_g$ Todd-vs-$\widehat A$-genus confusion (Critical).**
  Failure: Vol II inscribed $F_1 = -\kappa/12$, $F_2 = \kappa^2/240$
  under Todd-genus convention. Both values are false under the
  programme-canonical $\widehat A$-genus / Mumford / Faber-Pandharipande
  convention where $F_g = \kappa \cdot \lambda_g^{FP}$. Source:
  `beilinson_swarm_2026_03_31.md` S17. Correct values (scalar
  channel, $\widehat A$): $F_1 = \kappa / 24$, $F_2 = 7 \kappa / 5760$,
  with $\kappa = c / 2$ for Virasoro (NOT $(c - 26) / 2$). Counter:
  (a) name Hodge-genus convention before any $F_g$ coefficient;
  (b) substitute $c$ and verify against
  `chapters/examples/landscape_census.tex`;
  (c) cross-check compute/lib/. Fixed across Vol II `rosetta_stone.tex`,
  `ht_bulk_boundary_line.tex`, preface.

- **AP802 -- "Com$^!$ = Lie $\simeq$ Com" mathematically false (Critical).**
  Failure: Vol II inscribed Swiss-cheese self-duality as
  "Com$^!$ = Lie $\simeq$ Com", conflating Koszul duality of operads
  with Koszul duality of algebras. Source:
  `beilinson_swarm_2026_03_31.md` S18. Correct: $\mathrm{SC}$ is
  Koszul (Livernet 2015 J Reine Angew Math 701) and
  $\mathrm{SC}^!$ is EXPLICITLY DIFFERENT from $\mathrm{SC}$;
  Com$^!$ = Lie as operads, NOT Lie $\simeq$ Com. The chiral-algebra
  statement: $\mathrm{SC}^! = (\mathrm{Lie}, \mathrm{Ass},
  \text{shuffle-mixed})$ while $\mathrm{SC} = (\mathrm{Com},
  \mathrm{Ass}, \text{product-mixed})$. Counter: see
  AP-SC-NOT-SELFDUAL (AP166 + FM26 + B57).

- **AP803 -- Cobar sign mismatch signs appendix vs main text (High).**
  Failure: signs appendix had $+$ sign; main text used $-$ sign for
  $d_\Omega$. Main text's $-$ sign is correct (Lefevre-Hasegawa 2003;
  appendix had naively copied LV12 without accounting for
  $s \to s^{-1}$ convention change). Counter: at every cobar sign
  site, (a) trace desuspension convention, (b) verify counit on
  $B(A)$, (c) confirm $d^2 = 0$ at 3-string level. Source:
  `beilinson_swarm_2026_03_31.md` S20.

- **AP804 -- Riordan sequence wrong for $\widehat{\mathfrak{sl}}_2$
  bar cohomology at $n \ge 3$ (Critical).**
  Failure: manuscript asserted $\dim H^n B(\widehat{\mathfrak{sl}}_2) =
  R_n$ (Riordan: $1, 1, 3, 6, 15, 36, \ldots$) at bar degree $n$.
  Correct: $\dim H^n = 2 n + 1$ (linear). Witness: SC-bar swarm
  `sc_bar_swarm_2026_04_07.md` via direct Chevalley-Eilenberg
  computation of $H^*(\widehat{\mathfrak{sl}}_2)$. Counter: do NOT
  quote Riordan for chiral-bar cohomology.

- **AP805 -- $\kappa(\mathrm{SVir}_c) = (c + 11) / 2$ (High).**
  Failure: manuscript asserted $\kappa(\mathrm{SVir}_c) = (c + 11) / 2$.
  Correct: $\kappa(\mathrm{SVir}_c) = (3 c - 2) / 4$ (Kac-Wakimoto
  1988; Friedan-Martinec-Shenker 1986). Counter: consult
  `chapters/examples/landscape_census.tex` super-Virasoro row.

### Wave 3 (2026-04-01 adversarial swarm)

- **AP806 -- "Every theorem consequence of MC equation" false universal
  (Critical).** Beilinson: Theorem B (Positselski) is INDEPENDENT of
  the Maurer-Cartan equation; its proof uses PBW + Barr-Beck-Lurie
  monadicity. The universal claim "every theorem is a consequence of
  $D \Theta_A + (1/2) [\Theta_A, \Theta_A] = 0$" is literally false.
  Source: `adversarial_swarm_2026_04_01.md`. Counter: do NOT fabricate
  MC-based routes for Theorem B, Theorem H concentration, or chiral
  Positselski adjoint-equivalence.

- **AP807 -- MCG-equivariance 3-line remark (High).** Segal: MCG
  equivariance of the genus-$g$ bar complex was asserted via a
  3-line remark with no proof. Counter: either provide explicit
  action of MCG on $B^{(g)}(A)$ via Deligne-Mumford on
  $\overline{\mathcal M_g}$ with Atiyah-Segal compatibility at each
  boundary stratum, or label conjectural.

- **AP808 -- "Modularity" conflation Segal-vs-clutching (High).**
  Segal/Beilinson: Segal modularity = MCG-invariance CONSTRAINT;
  manuscript modularity = clutching CONSTRUCTION. The equivalence
  between the two is a theorem (Huang 2005 under sewing axioms), not
  a definition. Counter: always specify which modularity is meant.

### Wave 4 (2026-04-05 arithmetic frontier)

- **AP809 -- Shadow Ferrero-Washington failure at $p = 2, 5$ for
  Virasoro (High).** Virasoro is shadow-irregular: $\mu \ne 0$ at
  $p = 2$ (growth $\sim -1$) and $p = 5$ ($\sim -0.5$). Witness:
  `compute/lib/padic_shadow_iwasawa_engine.py`. Counter: do NOT
  assume shadow Iwasawa $\mu = 0$ universally.

- **AP810 -- Shadow self-duality at $c = 13$ through $r = 12$ (High).**
  $F_g(c) + F_g(26 - c) = 13 \lambda_g$ through genus 5; RTF
  vanishes for all test functions at $c = 13$. Counter:
  self-duality extends from $\kappa$ to the full tower.

- **AP811 -- $c(W_N)$ generic-level formula (Critical).** Earlier
  manuscript: $c(W_N) = (N - 1) - N(N^2 - 1)(k + N - 1)^2 / (k + N)$.
  Correct (Fateev-Lukyanov 1988): $c(W_N) = (N - 1)(1 - N(N + 1) /
  (k + N))$. Propagated across ~33 library + 20+ test files.
  FF complementarity: Vir $c + c_{\mathrm{gh}} = 26 \to 2$ at $k = 1$;
  $W_3$: $100 \to 4$; $W_4$: $246 \to 6$. Counter: evaluate at
  $k = 0, -N$ against primary Fateev-Lukyanov 1988 + FBZ Ch. 15.

- **AP812 -- Affine-KM Fredholm exponent colored-partitions-vs-dim
  (Critical).** `affine_km_sewing_engine.py` used colored-partitions
  where correct is $\dim \mathfrak g$. Caught in Round 3.5. Counter:
  Fredholm determinant of sewing form uses Cartan dimension; the
  two agree only at $\mathfrak{sl}_2$ (rank 1).

- **AP813 -- Conductor NOT Koszul-symmetric (Medium).** Arithmetic
  conductor $N(A)$ is symmetric under $c \leftrightarrow 26 - c$
  ONLY at $c = 13$; in general $N(A) \ne N(A^!)$. Ising witness
  primes $\{2, 3, 7, 11, 17, 19\}$; $\det(\Delta) | N(A)$
  universally; growth exponent $\gamma \approx 3.4$.

- **AP814 -- Shadow curve genus always $0$ (Medium).** Complexity
  lives in the DISCRIMINANT $\Delta$, not in shadow-curve genus.
  Every shadow curve is rational. Weight filtration:
  $c^a (5 c + 22)^b$. Counter: distinguish "shadow curve" (always
  $\mathbb P^1$) from "shadow motive" (non-trivial in $\Delta$).

- **AP815 -- $d_{\mathrm{arith}}$ independent of $d_{\mathrm{alg}}$
  (Critical).** Lattice VOA counterexample via
  `darith_full_landscape_engine`. Max $d_{\mathrm{arith}} = 7$ at
  $M(7, 8)$. Counter: the two invariants are Koszul-independent.

- **AP816 -- $W_3$ shadow is irreducible GL$_3$ (Medium).**
  $\rho^{\mathrm{sh}}_{W_3}$ is irreducible cuspidal on GL$_3$, not
  Heisenberg $\oplus$ Eisenstein. Counter: do NOT reduce $W_3$
  shadow to abelian$+$1-d pieces.

- **AP817 -- $\kappa \ne f(D^2)$ (High).** Shadow independent of
  quantum dimension $D$. Counter: never reconstruct $\kappa$ from
  $D^2$ data alone.

- **AP818 -- Soft graviton hierarchy $=$ shadow degree hierarchy
  (High).** $S^{(0)}, S^{(1)}, S^{(2)} \leftrightarrow \kappa, S_3,
  Q^{\mathrm{contact}}$; universality breaks at degree 4. MHV:
  gluon $\to$ L, graviton $\to$ M, photon $\to$ G.

- **AP819 -- DS decreases $d_{\mathrm{arith}}$ (Medium).** DS:
  $d_{\mathrm{arith}}: 1 \to 0$ while depth $: 3 \to \infty$.
  Counter: do NOT assume DS preserves arithmetic complexity.

- **AP820 -- $N = 2$ pure gauge $\kappa = 0$ universally (Medium).**
  SW curve $\ne$ shadow curve structurally; $\kappa \ne c/2$
  algebraically for real $k$.

### Wave 5 (2026-04-05 extremal frontier $\kappa$ cross-verification)

- **AP821 -- DS central-charge missing factor $N$ (Critical).**
  `theorem_c_complementarity.py` had a DS formula missing factor $N$.
  Caught by XVER-34. Correct: $c(W_N) = (N - 1)(1 - N(N + 1) /
  (k + N))$; missing factor displays $\kappa + \kappa'$ at half value.

- **AP822 -- G$_2$ exceptional metric factor-3 error (Critical).**
  `ExceptionalRootSystem._omega_metric` over-counted Casimirs by 3
  for G$_2$. Correct: G$_2$ Casimir coefficient is $1/12$
  (simply-laced is $1/4$).

- **AP823 -- D-module purity direction inversion (High).**
  Manuscript header asserted $(x) \Rightarrow (xii)$; proved
  direction is $(xii) \Rightarrow (x)$. Counter: verify both
  directions before quoting Koszulness equivalence.

- **AP824 -- Stale conjecture-vs-theorem labels after resolution
  (High).** `conj:operadic-complexity` was proved yet 4 stale
  references quoted it as open. Counter: on resolution, grep all
  three volumes for the label and flip atomically.

- **AP825 -- Modular envelope adjunction status contradiction
  (High).** `conj:platonic-adjunction` labelled both conjecture and
  theorem in different files. Scope: proved ON KOSZUL LOCUS; global
  conjectural. Counter: atomic rename in same tool-call batch
  (AP-LABEL-DISCIPLINE).

- **AP826 -- Moonshine $\kappa(V^\natural) = 12$ not $24$ (Critical).**
  $\kappa(V^\natural) = 12$, class M, Monster $\to$ 1D shadow via
  Schur. Factor-2 error from conflating $c(V^\natural) = 24$ with
  $\kappa$ (which is half the central charge).

- **AP827 -- $E_8 \times E_8$ vs $D_{16}^+$ first differ at genus 2
  (Medium).** The two $c = 24$ lattices are indistinguishable at
  genus 1 (both give $E_4^3$) but separate at genus 2 via Siegel
  theta. Counter: never claim genus-1 distinction.

### Wave 6 (2026-04-07 frontier 105-agent swarm)

- **AP828 -- ChirHoch amplitude $[0, 2]$ vs virtual dimension
  (Medium).** Amplitude $\ne$ virtual dimension. Amplitude is range
  of non-vanishing degrees; virtual dim is Euler characteristic.
  Counter: in Theorem H write "concentrated in $\{0, 1, 2\}$" not
  "dimension 2."

- **AP829 -- $d_{\mathrm{gen}}$ vs $d_{\mathrm{alg}}$ silent
  conflation (High).** $d_{\mathrm{gen}}(\mathrm{Vir}) = 3$
  (Gelfand-Fuks); $d_{\mathrm{alg}}(\mathrm{Vir}) = \infty$.
  Bare $d(\mathrm{Vir})$ ambiguous. Counter: subscript mandatory;
  refuse bare $d(\cdot)$.

### Wave 7 (2026-04-07 SC bar swarm)

- **AP830 -- $d$ NOT a coderivation (Critical).** Full bar differential
  $d = \sum_g D^{(g)}$ is NOT a coderivation w.r.t. deconcatenation on
  $B(A)$. Each genus component $D^{(g)}$ IS a coderivation. Counter:
  specify genus component; never quote $\Delta \circ d = \ldots$ for
  the full sum.

- **AP831 -- BV $=$ bar false at chain level for class M (Critical).**
  $\delta_4 = Q^{\mathrm{contact}} \cdot m_0$ is exact in
  $D^{\mathrm{co}}(A)$ but NOT zero on the nose. Counter: state
  "BV $=$ bar holds at chain level for G/L/C; at coderived level for
  M." 107 tests G/L/C; 76 tests M in $D^{\mathrm{co}}$.

- **AP832 -- $K_{\mathrm{BP}} = 76$ propagation (High).** Correct
  $K_{\mathrm{BP}} = c_{\mathrm{BP}} + c_{\mathrm{BP}}^! = 196$;
  $\kappa_{\mathrm{BP}} + \kappa_{\mathrm{BP}}^! = 98/3$.
  Counter: evaluate $c + c^!$ before quoting.

- **AP833 -- ChirHoch$^*$ polynomial-ring overclaim (High).**
  "Polynomial" meant polynomial GROWTH of Betti numbers, not polynomial
  ALGEBRA. ChirHoch$^*(A)$ is BOUNDED in $\{0, 1, 2\}$. Counter:
  cite Hilbert-series, not algebra structure.

### Wave 8 (2026-04-08 149-agent final convergence)

- **AP834 -- $\delta F_2^{\mathrm{cross}}(W_3) = (c + 204) / (16 c)$
  (Critical).** Naive graph sum gave $(c + 120) / (16 c)$; correct
  value cross-verified by 5 agents. Resolves op:multi-generator-
  universality NEGATIVELY: scalar formula $F_g = \kappa \lambda_g$
  fails for multi-weight algebras at $g \ge 2$.

- **AP835 -- Segal-vs-clutching modularity equivalence
  unproved (Medium).** Layer at construction level.
  Counter: write theorem, not definition.

- **AP836 -- Heisenberg / abelian-CS single-atom (High).** Abelian
  CS boundary $=$ Heisenberg $=$ abelian KM at level $k$ with same
  OPE $J(z) J(w) \sim k / (z - w)^2$. Counter: single-atom
  escalation, not two-atom; abelian CS $\ne$ non-trivial alternative.

- **AP837 -- $\beta\gamma$ $\Theta_{\le 2} \ne 0$ on charged stratum
  (High).** $\kappa(\beta\gamma) = -1 \ne 0$ so $\Theta_{\le 2} \ne 0$.
  The issue: $\beta\gamma$'s $\kappa$ lives on a CHARGED STRATUM
  (weight-0 $\gamma$ violates positive grading). Counter: never
  assert $\Theta_{\le 2}(\beta\gamma) = 0$ without charged-stratum
  scope.

### Wave 9 (2026-04-09 editorial closure)

- **AP838 -- Standalone papers AP126 residuals (Medium).** After
  main-manuscript sweeps, grep `standalone/*.tex` for bare
  `\Omega / z` and patch. 22 residuals fixed in Wave 9.

- **AP839 -- Vol III AP113 residuals (High).** ~165 bare $\kappa$
  across Vol III patched. Every $\kappa$ in Vol III must carry
  family superscript.

### Wave 10 (2026-04-10 deep rectification)

- **AP840 -- AP142 local-global conflation (Critical).**
  "Koszul duality over a point" $\ne$ "over $\mathbb P^1$."
  Three errors: (a) homotopy retract is DATA; (b) formal disk
  $\mathbb D \ne$ point; (c) $\mathbb A^1$ already has Arnold
  relations. ~45 fixes / 25 files / 3 volumes. Counter: at every
  comparison, name (i) space, (ii) comparison data,
  (iii) on-the-nose vs extra-structure.

- **AP841 -- DS ghost-charge background shift omission (Critical).**
  $W_6, W_7$ off by factor 41 because $Q^{\mathrm{bg}} = (\rho, \alpha_-
  - \alpha_+)$ was dropped. Counter:
  $c = \mathrm{rk}(\mathfrak g) - 12 |Q^{\mathrm{bg}}|^2$ with
  explicit background shift per nilpotent orbit.

- **AP842 -- Convention coexistence without bridge (High).**
  Maintaining two conventions without explicit identity creates silent
  mismatch. Counter: every coexistence must carry an IDENTITY:
  $k \Omega_{\mathrm{tr}} = \Omega / (k + h^\vee)$;
  $\mathrm{Todd} \cdot e^{c_1 / 2} = \widehat A$; Molev $1 - P / u$
  vs Drinfeld $u I + \Psi P$ sign reversal.

- **AP843 -- Restructuring propagation debt $O(N)$ (Medium).** Any
  structural restructuring has $O(N)$ follow-up. Must happen IN THE
  SAME SESSION.

- **AP844 -- Mega-campaign straggler commits (Medium).** Large
  rectifications leave 5-10 trailing commits; end-of-session
  inspection mandatory.

- **AP845 -- Circular proof routing 3-hop (Critical).** $L_1 \to L_2
  \to$ theorem; trace 3 levels deep. Counter: install ROUTING REMARK
  identifying non-circular anchor (Arnold / PBW / Shelton-Yuzvinsky).

- **AP846 -- r-matrix normalisation convention-dependence (Critical).**
  Trace: $k = 0 \Rightarrow r = 0$. KZ: $k = 0 \Rightarrow r \ne 0$
  for non-abelian. Counter: name convention at every r-matrix site.

- **AP847 -- Resolution propagation failure (Medium).** Resolve
  propagation in same session via explicit propagation queue.

- **AP848 -- $c(W_3, k)$ Fateev-Lukyanov (Critical).** Wrong:
  $c(W_3, k) = 2(k - 9) / (k + 3)$. Correct:
  $c(W_3, k) = 2 - 24(k + 2)^2 / (k + 3)$.

- **AP849 -- $F_2$ exponent-of-$\kappa$ (Critical).** Wrong
  $F_2 = 7 \kappa^2 / 5760$; correct $F_2 = 7 \kappa / 5760$ (linear).

- **AP850 -- $\kappa(\beta\gamma_1) = -2 \to 1$ (Critical).** Correct
  value $1$ from Euler-class anomaly. Cascade $F_1 = -1/6 \to 1/24$.

### Wave 10 adversarial 583-agent

- **AP851 -- Theorem H "dim $\le 4$" overclaim (Critical).**
  $\widehat{\mathfrak{sl}}_2$ gives $\dim = 5$. Counter: state only
  AMPLITUDE $\{0, 1, 2\}$; drop dimension bound.

- **AP852 -- Concordance $\kappa = 0$ vs $\Theta_A = 0$ conflation
  (High).** Higher-degree components $\Theta_A^{(\ge 3)}$ INDEPENDENT
  of $\kappa$.

- **AP853 -- MC3 "all types" overclaim (High).** Layer 3 unconditional
  in type A only; non-simply-laced/exceptional conditional.

- **AP854 -- L$_\infty$ strictness Hodge SDR conditions (High).**
  Must name triple $(h, i, p)$ with $h^2 = 0$, $h i = 0$, $p h = 0$
  (Huebschmann-Stasheff 2002).

- **AP855 -- $\kappa_{\mathrm{ch}}$ / $\kappa_{\mathrm{BKM}}$
  conflation Vol III (Critical).** $\kappa_{\mathrm{BKM}}(K3 \times E)
  = \kappa_{\mathrm{ch}} + \chi(\mathcal O_{\mathrm{fiber}}) = 3 + 2
  = 5$. See AP352 / AP-CY61.

- **AP856 -- KRW central-charge formula incomplete (Critical).**
  `hook_type_w_duality.py`; cross-check against KRW 2003 Adv Math 185
  hook-type W-algebras.

### Wave 11 (bar-cobar v4 + zeta-zeros)

- **AP857 -- Shadow L-function critical abscissa $c^* \approx 6.12$
  (Medium).** Below $c^*$ entire; above $c^*$ meromorphic with zeros
  in critical strip. Counter: specify $c^*$ scope for any shadow
  L-function claim.

- **AP858 -- Sewing coordinate orientation upper-vs-lower half-plane
  (Medium).** $q = e^{2 \pi i \tau}$ requires upper half-plane
  $\tau$. Wrong orientation converts $q \to \bar q$, flipping $L_0$.

### Wave 12 Afternoon 2026-04-12 (AP150-AP175)

- **AP859 -- Multi-step structure confabulation (Critical).** AP150:
  the "$E_n$ circle" advertised as 5 proved arrows. Only 1 proved
  (HDC). Counter: arrow-by-arrow verification.

- **AP860 -- Convention clash in single file (High).** AP151:
  different $\hbar$ definitions within one file. Counter: grep all
  $\hbar$ definitions at insertion.

- **AP861 -- "Ordered" labeled-vs-time-ordered (High).** AP152:
  chiral "ordered" $=$ labeled on $C$; QFT "time-ordered" $=$ on
  $\mathbb R$. Bar uses labeled.

- **AP862 -- $E_3$ scope inflation (Critical).** AP153: $E_3$
  requires $E_\infty$-chiral input; classical Deligne gives only
  $E_2$ from $E_1$. See AP-TOPOLOGIZATION.

- **AP863 -- Two distinct $E_3$ structures (High).** AP154:
  algebraic (HDC) vs topological (Sugawara) $E_3$ are different
  objects at class M; name which.

- **AP864 -- Novelty architectural-vs-computational (Medium).**
  AP155: bar complex is architecturally novel, computationally
  recovers known data.

- **AP865 -- Weierstrass $\wp$ vs zeta quasi-periodicity (High).**
  AP156: theta vs Weierstrass-zeta normalisations differ at
  $\eta_1 = (\pi^2 / 6) E_2(\tau)$.

- **AP866 -- Degeneration-dependent invariants (High).** AP157:
  separating vs non-separating $g = 2$ degenerations give different
  $F_2$, $\delta F_2^{\mathrm{cross}}$.

- **AP867 -- Yangian 4-type conflation (Critical).** AP159: Hopf /
  dg-shifted / chiral / spectral; use
  $Y^{\mathrm{Hopf}}, Y^{\mathrm{dg}}, Y^{\mathrm{ch}}, Y^{\mathrm{spec}}$.

- **AP868 -- Bare "Hochschild" conflation (Critical).** AP160: 89
  instances bare; must qualify chiral / topological / categorical.

- **AP869 -- Five $E_1$-chiral notions (Critical).** AP161: (A)-(E);
  each with own derived centre.

- **AP870 -- "$\mathbb R \times C$ inhabitance" unjustified
  (Medium).** AP163: requires chiral Deligne-Tamarkin formality.

- **AP871 -- Chiral Gerstenhaber $\ne$ topological (High).** AP164:
  agree for $E_\infty$-chiral; diverge at $E_1$.

- **AP872 -- B-cycle monodromy $i^2$ error (Critical).** FM24:
  $q = e^{2 \pi i \hbar}$ becomes real if $\hbar \in (\pi i) \mathbb Z$.
  Substitute integer $k$ and verify $|q| = 1$.

- **AP873 -- Bulk substring replacement corruption (Critical).**
  FM42: `replace_all` "arity" $\to$ "degree" corrupts singularity $\to$
  singul-degree. 45 corruptions in Vol III. Checklist of protected
  words.

### Wave 12 deep adversarial (380+250-agent)

- **AP874 -- Genus-universality GRR gap (Critical).** AP225:
  all-genera $\mathrm{obs}_g = \kappa \lambda_g$ requires
  clutching-uniqueness or GRR (Mumford 1983). Split genus-1
  unconditional vs all-genera conditional.

- **AP875 -- $K_0$-class vs scalar confusion (High).** AP226:
  $\kappa$ is a COMPLEX SCALAR, not integer; cannot serve as $K_0$
  multiplicity. Replace with Chern character.

- **AP876 -- ProvedHere-as-forwarding (Medium).** AP227: proof body
  "By Theorem $X$" is ProvedElsewhere.

- **AP877 -- Anomaly-Koszul dependency inversion (High).** AP228:
  Theorem D $\leftarrow$ `thm:anomaly-koszul` $\leftarrow$
  `cor:kappa-additivity`. Check dependency direction.

- **AP878 -- SC-formality propagation debt across compute libraries
  (High).** AP229: `swiss_cheese_cy3_e1.py` carried stale class L
  SC-formal claim. After updating classification, grep ALL compute
  libraries.

- **AP879 -- Genus-1 sufficient claimed all-genera (High).** AP230:
  when proof needs only $g = 1$, cite genus-1 result; do not route
  through universality.

- **AP880 -- Draft artefacts in theorem statements (Medium).**
  AP231: raw `\textup{(LOCAL)}` / `(DRAFT)` / `(TODO)` inside
  theorem environments.

- **AP881 -- Duality-clause family scope over-claim (High).**
  AP232: "affine KM and free-field algebras" vs proof covers only
  affine KM $+$ principal $W_N$.

- **AP882 -- MC3 compact / completed comparison gap (High).**
  AP233: bounded-stratum to completed-category extension requires
  specific compact-generation theorem; otherwise conditional.

### Wave 12 campaign AP186-AP224 cross-reference

Full first-principles entries in
`compute/audit/new_antipatterns_wave12_campaign.md` for: AP186
ProvedHere-without-proof; AP187 orphaned chapters; AP188 empty
sections; AP189 dead labels; AP190 hidden imports; AP191 circular
proof chains; AP192 scope inflation; AP193 biconditional-forward-
only; AP194 curved complex with flat tools; AP195 five-object-
conflation-in-prose; AP196 SC misattribution; AP197 bare
Hochschild; AP198 Whitehead scope; AP199 strong-filtration
direction; AP200 transfer comparison gap; AP201 Baxter
vacuous-constraint fallacy; AP202 coderived element-wise; AP203
class-M harmonic mechanism unproved; AP204 genus-0 boundary
violation; AP205 reflexivity hidden in double-dual; AP206
Verdier-dual-vs-cobar object switch; AP207 centre-vs-bar lift
missing; AP208 Theorem A algebra/coalgebra flip; AP209
missing-lemma cited; AP210 topologization chain-vs-cohomological;
AP211 test file absent for compute engine; AP212 TODO/FIXME
unresolved; AP213 stub chapter; AP214 cross-volume bridge
outdated; AP215 preface advertising stronger than body; AP216
Koszul equivalence (vii) genus-0 scope; AP217 Koszul equivalence
(viii) free-polynomial overclaim; AP218 SC-formality via
non-existent metric; AP219 depth-gap witness on weight-changing
line; AP220 $D^2 = 0$ wrong geometric space; AP221 Gerstenhaber
single-insertion; AP222 Theorem H configuration-space collapse;
AP223 bar-coalgebra / Koszul-dual conflation in Theorem H;
AP224 README / metadata scope inflation. Do NOT renumber; cite
canonical AP186-AP224.

### Wave 13 (K3-BKM landscape-census anchoring)

- **AP883 -- $c_{4d}(A_1, \Sigma_{0, 24}) = 107/6$ anchor (Critical).**
  Established via Chacaltana-Distler pants, WOV-2 lock.
  W14 erroneously retracted to $(26, -312)$ via formula
  $(12(g - 1) + 7 n) / 6$; W15 healed back via correct
  $(5 n - 13) / 6$. See AP339.

- **AP884 -- $c_3 = 176256 \cdot [H_3]$ cache trap (Critical).**
  $176256 = p_{24}(5) = \chi(\mathrm{Hilb}^5(K3))$ is a Göttsche
  coincidence. Correct $c_3 = -8$ (Bruinier reduced). Conversion
  $-22032$ to Gritsenko-Nikulin. See AP341.

- **AP885 -- Class-M $E_3$-bar $6^g$ overcount (High).** Correct:
  $e_3^{K3 \times E} = 2 \mathrm{Vol}(E) (2\pi i)^3$; factor-of-2
  Mukai double-count removed.

- **AP886 -- Bi-based Ran factorisation (Medium).**
  $\mathrm{Ran}_{X, Y}$ for K3 $\times$ E is NOT reducible to
  $\mathrm{Ran}(X) \times \mathrm{Ran}(Y)$; $\mathbf H_{\Delta_5}$
  inhabits bi-based Ran. Künneth obstruction
  $\chi(\mathcal O_{K3}) \chi(\mathcal O_E) = 0$.

- **AP887 -- Conway $V^{s\natural}$ as 5th $\Psi$-image (Critical).**
  Hidden structure: $V^{s\natural}$ is the $\mathbb Z/2$-super-twin
  of $V^\natural$ in Duncan 2007 §6 commutative orbifolding diamond,
  inheriting $(K, \hbar^2) = (2, -1/2)$ from Monster. Not a 5th
  independent row. Downgrade to `conj:bkm-conway-psi-fifth-image`.

### T1#1 genus-5 / Andreotti--Mayer / universal closure additions

- **AP890 -- Andreotti--Mayer codim is rank index, not codimension
  (Critical).** $\mathrm{AM}_g$ is cut by hessian rank $\leq g-3$;
  $g-3$ is the *rank index*, NOT the codim of $\mathcal J_g$ in
  $\overline{\mathcal A_g}$. Beauville 1977 Cor.~II.4 gives
  $\mathrm{codim}\,\mathcal J_5 = 3$; Debarre 1992 Thm.~1.1 gives
  $\geq g-1$ for $g\geq 5$ (exact at $g=5$). Schottky locus codims:
  $g=4{:}1$, $g=5{:}3$, $g=6{:}6$. See Pattern 296.

- **AP891 -- $k_{\max}(g)=g(g+1)/2$ Krull bound vs saturation
  (High).** The closure $k_{\max}(g)=g(g+1)/2$ is the
  Krull-dimension upper bound on the depth of the NL-unipotent
  Malcev ladder on $\overline{\mathcal A_g}$. *Saturation*
  (existence of a non-empty depth-$k_{\max}$ admissible NL
  intersection of CM type) is unconditional only at $g\leq 6$ via
  Shimura 1963 + Lange 1975 cyclotomic CM (degree-$2g$ cyclotomic
  fields with $\varphi(n)=2g$). For $g\geq 7$ off cyclotomic ladder,
  conditional on André--Oort (Edixhoven--Yafaev 2003). See Pattern 297.

- **AP892 -- Single CKS unipotent vs commuting $k$-unipotent
  (Medium).** A degeneration with a single unipotent invokes
  one-variable Schmid 1973 + one-variable CKS to produce a single
  $\mathfrak{sl}_2$. A degeneration with $k$ pairwise-commuting
  unipotents on a transverse $k$-polydisc invokes $k$-variable CKS
  \cite[\S 3.7]{CattaniKaplanSchmid1986} to produce
  $\mathfrak{sl}_2^{\oplus k}$ via the polarised relative weight
  filtration. The Andreotti--Mayer locus at $g=5$ requires
  $3$-variable; the Schottky divisor at $g=4$ requires only
  $1$-variable. See Pattern 298.

- **AP893 -- Schottky/AM block-count growth $B^{\mathrm{Sch/AM}}(g)$
  (Medium).** Stratum count
  $B^{\mathrm{Sch/AM}}(2)=B^{\mathrm{Sch/AM}}(3)=0$,
  $B^{\mathrm{Sch/AM}}(4)=2$ (pure Schottky, Schottky$\cap$NL),
  $B^{\mathrm{Sch/AM}}(5)=4$ (AM rank-$2$ open, deeper rank-$1$,
  Jacobian sub-locus, AM$\cap$NL mixed). Conditional
  $B^{\mathrm{Sch/AM}}(6)=6$ pending verification of Debarre's
  component count for $\mathcal J_6\subset\AMlocSix$ (Debarre 1992).
  Total stratum count $N_g = 1+g(g+1)/2+B^{\mathrm{Sch/AM}}(g)+1+(g+1)+2$
  giving $N_2=8, N_3=13, N_4=20, N_5=28, N_6\overset{?}{=}38$.

### Wave 27 (2026-04-21): Humbert--Heegner admissibility filter

- **AP903-HH -- Humbert--Heegner admissibility $n\equiv 3,5\pmod 8$ on
  the pentagon coboundary tower $\phi^{(n)}$ (High). Cross-referenced
  at Vol I tip-cache row 304 (AP890 anchor).**
  **Ghost.** The pentagon tower $\{\phi^{(n)}\}_{n\ge 3}$ of
  Definition~\ref{def:phi-n-pent-EK} has a well-defined three-filter
  admissibility structure on the K3 $A_\infty$-Humbert regime of
  $\mathbf H_{\Delta_5}$; the Eichler--Zagier 1985 polar-support cutoff
  $\Delta\ge -1$ on $\phi^{K3}_{0,1}$ is real.
  **Precise error.** Bare Padovan-dimension $d_n = d_{n-2}+d_{n-3}$
  count without the HH admissibility filter overcounts. Most
  Padovan-admissible $n$ (all $n\ge 3$ except $n=4$) are
  Humbert--Heegner-FORBIDDEN: the paramodular lattice sum
  $\sum_{4NM-\ell^2=-D_n} c_{\Phi_{10}/\eta^{24}}(N,\ell,M)$ with
  $D_n=(n-3)/2$ is non-empty iff $D_n\bmod 4\in\{0,1\}$, equivalently
  $n\equiv 3,5\pmod 8$.
  **Correct.** $\phi^{(n)}\big|_{\mathrm{K3\text{-}Humbert}}\ne 0$ iff
  (i) $n\equiv 3,5\pmod 8$ AND (ii) the $d_n$-dimensional Brown
  canonical basis is non-empty AND (iii) $D_n\le 1$ (Eichler--Zagier
  polar cutoff). First non-vanishing: $\phi^{(3)}$ (Drinfeld pentagon
  cocycle, $D_3=0$, $C(0)=20$) and
  $\phi^{(5)} = -2\cdot[\mathrm{gen}]^{\otimes 5}$ (Gritsenko--Nikulin
  1998 sign, $D_5=1$, $C(-1)=2$). HH-admissible $n$ in $[3,36]$:
  $\{3,5,11,13,19,21,27,29,35\}$. HH-forbidden Padovan-positive $n$
  (e.g., $4,6,7,8,9,10,12,14,15,16,17,18,20,22,23,24,25,26,28,30,31,
  32,33,34,36$) all give $\phi^{(n)}\big|_{\mathrm{K3\text{-}Humbert}}=0$.
  For HH-admissible $n\ge 11$: $D_n\ge 4$ and $\phi^{(n)}=0$ by polar
  support.
  **Reference table** ($d_n$ = Padovan; HH $=n\equiv 3,5\pmod 8$;
  $\phi^{(n)}$-K3 $=\phi^{(n)}\big|_{\mathrm{K3\text{-}Humbert}}$):
  $(n, d_n, D_n, \mathrm{HH}, \phi^{(n)}\text{-K3}) =
  (3, 1, 0, Y, \text{non-zero});
  (4, 0, 1/2, -, 0);
  (5, 1, 1, Y, -2[\mathrm{gen}]^{\otimes 5});
  (6, 1, 3/2, -, 0);
  (7, 1, 2, N, 0);
  (8, 2, 5/2, -, 0);
  (9, 2, 3, N, 0);
  (10, 2, 7/2, -, 0);
  (11, 3, 4, Y, 0 \text{ polar});
  (12, 4, 9/2, -, 0);
  (13, 5, 5, Y, 0 \text{ polar});
  (19, 17, 8, Y, 0 \text{ polar});
  (21, 28, 9, Y, 0 \text{ polar});
  (27, 90, 12, Y, 0 \text{ polar});
  (29, 149, 13, Y, 0 \text{ polar});
  (35, 504, 16, Y, 0 \text{ polar})$.
  **Primary.** Eichler--Zagier 1985 *Prog. Math.* 55 Thm 9.3
  (polar-support cutoff, index-$1$ weak Jacobi forms); Gritsenko--Nikulin
  1998 *J.~Reine Angew. Math.* 507 (Humbert--Heegner structure,
  paramodular $\Phi_{10}/\eta^{24}$ sign convention); Bruinier 2002
  LNM 1780 (Chern class on Heegner divisors); Brown 2012 *Ann. Math.*
  175 Thm 1 (Padovan motivic MZV dimension).
  **Location.** Vol I
  `chapters/theory/shadow_tower_higher_coefficients.tex`
  Theorem~\ref{thm:phi-n-humbert-heegner-admissibility} (lines
  4364--4433) and Remark~\ref{rem:three-filter-composite-scope} (lines
  4435--4450); `appendices/first_principles_cache.md` row 304 (AP890
  anchor); `notes/first_principles_cache_comprehensive.md`
  Pattern~299.
  **Cross-volume.** Vol II: `V2-AP127` (cache row 137 +
  antipatterns_catalogue.md Wave-27 entry) -- HH filter reads as the
  Humbert-stratification refinement of the $\mathsf{SC}^{\mathrm{ch,top}}$
  coloured bar differential (Swiss-cheese factorisation at the
  paramodular K3 elliptic genus polar-support slice). Vol III:
  `AP-CY142` (cache row V16 + antipatterns_catalogue.md + comprehensive
  cache) -- HH filter reads as the admissible-discriminant set in the
  $c_{K3}$ Fourier expansion of the EOT K3 elliptic genus
  (Gritsenko--Nikulin 1998 §3 structure).
  **Related.** Pattern 132 / AP-CY78 (YD tower weight $\lfloor n/2\rfloor+1$);
  row 409 tip-cache AP889 (HH-admissibility voice precursor).

### Wave 28 (2026-04-21): Cross-volume nomenclature -- pseudo-character disambiguation

- **AP904-PC -- Pseudo-character name-clash across volumes (Critical).**
  **Ghost.** The string "pseudo-character" is used in multiple chapters
  across Vol~I and Vol~III with mathematically different meanings and no
  inline disambiguation; a reader encountering either occurrence may
  conflate them.
  **Precise error.** Two objects:
  (i) **Arithmetic pseudo-character** (Rouquier 1996 *J. Algebra* 180;
  Taylor 1991 *Duke* 63; Wiles 1988 *Invent.* 94): $K$-valued class
  function on a profinite group over a field of char $0$ or $p>d$,
  satisfying $T(1)=d$, $T(gh)=T(hg)$, truncation
  $T_{d+1}\circ\mathrm{Alt}_{d+1}\equiv 0$. Its commutative-ring
  generalisation is the **Chenevier determinant**
  (Chenevier 2014 arXiv:1301.0635 §1.2 Def.~1.5): a homogeneous
  multiplicative polynomial law $D:R\to A$ of degree $d$ on a
  commutative ring $A$ satisfying Cayley--Hamilton of dimension $d$.
  On reduced $A$ the two coincide (Thm.~2.12 of *ibid.*); on non-reduced
  $A$ the Chenevier determinant is strictly stronger.
  (ii) **Logarithmic CFT pseudo-character / Lyubashenko pseudo-trace**
  (Creutzig--Ridout 2013 *Nucl. Phys. B* 875 Thm.~3.4; Lyubashenko--Majid
  1994 *J. Algebra* 166 Prop.~3.4; Kerler--Lyubashenko 2001 LNS 262 Ch.~4):
  the trace $S^{\mathrm{ps}}:\mathrm{End}(\mathcal L)\to\mathrm{End}(\mathcal L)$,
  $f\mapsto\mathrm{tr}_1(R_{21}(f\otimes\mathrm{id})R_{12}^{-1})$, on
  the coend $\mathcal L$ of a non-semisimple finite ribbon category,
  carrying the projective $\mathrm{SL}_2(\mathbb Z)$-action up to ribbon
  anomaly, and satisfying $(S^{\mathrm{ps}})^2 = \theta\cdot C$ (not
  multiplicativity, not Cayley--Hamilton).
  **Correct nomenclature.**
  - Vol~I's arithmetic object on the paramodular Hecke algebra
  $\mathbb T^{\mathrm{par}}_1$ and on the universal deformation ring
  $R^{\mathrm{def}}_{\Delta_5}$ is a **Chenevier determinant**
  $D^{\mathrm{Chen}}$, never a pseudo-character; the retraction is
  load-bearing on $R^{\mathrm{def}}_{\Delta_5}$ (non-reduced:
  reducibility-ideal nilpotent of Bellaïche--Chenevier 2009 Prop.~1.5.1)
  where Chenevier 2014 Thm.~2.12 fails. See
  `chapters/theory/derived_langlands.tex`
  Remark~\ref{rem:dl-w25-determinant-not-pseudocharacter}.
  - Vol~III's LCFT object on the coend $\mathcal L$ is the
  **Lyubashenko pseudo-trace** / logarithmic-$S$ pseudo-character,
  with the programme-canonical symbol $S^{\mathrm{ps}}$. Occurrences
  in `chapters/theory/modular_trace.tex` lines 3302, 3368, 3414, 3463
  (rewritten with LCFT disambiguation footnotes after Wave 28) and in
  `chapters/theory/quantum_groups_foundations.tex` line 1939 refer
  exclusively to this object.
  - The two objects are **not** specialisations of a common formalism.
  Collapsing them is a type error.
  **Location.**
  - Vol~I: `chapters/connections/concordance.tex` Cross-volume
  type-error registry entry T6 (proposition inscribed Wave 28); see
  also `chapters/theory/derived_langlands.tex`
  Remarks \ref{rem:dl-chenevier-determinant-vs-pseudocharacter},
  \ref{rem:dl-chenevier-axiom-comparison},
  \ref{rem:dl-w25-determinant-not-pseudocharacter};
  Theorems \ref{thm:dl-pseudocharacter-delta10},
  \ref{thm:dl-chenevier-nonreduced-delta5}.
  - Vol~III: `chapters/theory/modular_trace.tex`
  Theorems \ref{thm:modtr-sps-coend-formula},
  \ref{thm:modtr-fricke-log-block-eigenvalues},
  \ref{thm:modtr-log-verlinde} (LCFT disambiguation footnote
  `fn:modtr-pc-disambig` at line 3302); remark in
  `chapters/theory/quantum_groups_foundations.tex` line 1939 (inline
  disambiguation footnote).
  **Primary.** Rouquier 1996 *J. Algebra* 180; Taylor 1991 *Duke* 63;
  Chenevier 2014 arXiv:1301.0635; Bellaïche--Chenevier 2009
  *Astérisque* 324; Creutzig--Ridout 2013 *Nucl. Phys. B* 875;
  Lyubashenko--Majid 1994 *J. Algebra* 166; Kerler--Lyubashenko 2001
  LNS 262 Ch.~4.
  **Cross-volume.** Vol~I ↔ Vol~III; Vol~II has zero occurrences.
  **Related.** AP353 (Chenevier determinant type error -- Hecke-algebra
  scope); Pattern 295 / W25 canonical preamble; AP354 (Hecke primary
  form $E_4\Delta$ weight 16).
  **Rule.** Every new occurrence of "pseudo-character" in any volume
  must specify, within the same paragraph, whether it is (i) arithmetic
  (Chenevier determinant) or (ii) LCFT (Lyubashenko pseudo-trace). Bare
  "pseudo-character" without qualifier is forbidden going forward.

- **AP905 -- Bar-cobar adjunction direction inverted (Critical; cross-volume).**
  The chiral bar-cobar adjunction is $\Omega^{\mathrm{ch}} \dashv B^{\mathrm{ch}}$
  (cobar LEFT, bar RIGHT), matching Loday--Vallette *Algebraic Operads*
  2012 Ch.~11, Positselski 2011 Mem.\ AMS 996, and Francis--Gaitsgory
  2012 *Selecta Math.* 18 Theorem 3.1. The counit $\Omega^{\mathrm{ch}}
  (B^{\mathrm{ch}}(A)) \to A$ is the bar-cobar quasi-isomorphism; the
  unit $C \to B^{\mathrm{ch}}(\Omega^{\mathrm{ch}}(C))$ lives on the
  coalgebra side. Writing $B^{\mathrm{ch}} \dashv \Omega^{\mathrm{ch}}$
  inverts this: under $F \dashv G \Leftrightarrow \mathrm{Hom}(F-,-)
  \cong \mathrm{Hom}(-,G-)$ it would make bar a left adjoint from
  algebras to coalgebras, contradicting the counit direction and the
  Quillen-equivalence structure (the left Quillen functor preserves
  cofibrations; in bar-cobar the cofibrations live on the algebra side
  via free/cofibrant resolutions inverted by cobar).
  **Counter.** Use $\Omega^{\mathrm{ch}} \dashv B^{\mathrm{ch}}$ in
  both chain-level and $(\infty,1)$-categorical statements; both lanes
  use the SAME adjunction direction. When a different adjunction is
  intended (e.g.\ Verdier / FM-geometric / $\SCchtop$-coloured) name
  which one, do not recycle the $B \dashv \Omega$ notation.
  **Canonical violations.** Vol~I
  `chapters/theory/chiral_climax_platonic.tex` lines 893, 928 (the
  $(\infty,1)$-lane writes $\mathbf{B}_X \dashv \boldsymbol{\Omega}_X$,
  inverting the chain-level lane at line 918 that writes
  $\Omega^{\mathrm{ch}}(\bar B^{\mathrm{ch}}(\cA)) \to \cA$ as the
  counit — two lanes use OPPOSITE adjunction directions within a
  single proposition);
  `chapters/connections/concordance.tex` lines 365, 1410, 2459
  (tables listing $B_\kappa \dashv \Omega_\kappa$, contradicting line
  2530 which states $\Omega^{\mathrm{ch}}(\bar B^{\mathrm{ch}}(\cA))
  \to \cA$ as the Quillen-equivalence counit);
  `chapters/theory/introduction.tex` line 496 ($B^{\mathrm{ch}} \dashv
  \Omega^{\mathrm{ch}}$) vs lines 3080, 3125 ($\Omega^{\mathrm{ch}}
  \dashv B^{\mathrm{ch}}$);
  Vol~III `chapters/theory/quantum_chiral_algebras.tex:1524`
  ($B \dashv \Omega$).
  **Related.** AP18 (sign/normalisation); AP166 ($\SCchtop^!$
  self-dual conflation); Pattern 269 (chain-level vs
  $(\infty,1)$-categorical lanes).
  **Rule.** Every new adjunction $F \dashv G$ in manuscript prose
  MUST use $\Omega \dashv B$ (cobar left, bar right) unless another
  adjunction is named. "$B \dashv \Omega$" is forbidden going forward.

- **AP906 -- Theorem~B scope collapse to "on the Koszul locus"
  (cross-volume scope contradiction).**
  Vol~I `chapters/theory/theorem_B_scope_platonic.tex` establishes a
  THREE-LAYER scope for the chiral Positselski / bar-cobar-inversion
  theorem:
  (a) at each finite bar weight — unconditional on all four shadow
  classes G/L/C/M (Theorem~\ref{thm:chiral-positselski-at-each-weight});
  (b) on the weight-completed bar coalgebra — unconditional on all
  four classes (Theorem~\ref{thm:chiral-positselski-weight-completed});
  (c) on the raw (non-completed) direct-sum bar coalgebra with
  class-M input — Theorem~B FAILS
  (Proposition~\ref{prop:chiral-positselski-raw-direct-sum-class-M-false}).
  The correct quotable scope is Corollary
  `cor:positselski-applicable-families`: raw direct-sum on classes
  G and L only, weight-completed on all four classes.
  **Violation.** Vol~III writes "Vol~I Theorem~B on the Koszul locus"
  as a blanket citation (`chapters/theory/quantum_chiral_algebras.tex:942,1524`;
  `chapters/connections/bar_cobar_bridge.tex:812`). This collapses
  three scopes into one and obscures that the weight-completed
  statement (the scope actually USED in Vol~III for K3-Yangian
  applications at generic non-linear-OPE fibres) is unconditional on
  all four classes, not restricted to the Koszul locus. "On the
  Koszul locus" is the raw-direct-sum scope; the K3 Yangian on
  generic fibres uses the weight-completed scope because the
  quantum-toroidal has non-linear OPE.
  **Counter.** Write "Vol~I Theorem~B (weight-completed,
  unconditional on all four classes G/L/C/M)" or "Vol~I Theorem~B
  (raw direct-sum, on classes G and L and the Koszul locus)",
  whichever is actually applied. Bare "Vol~I Theorem~B on the Koszul
  locus" is too narrow for most Vol~II and Vol~III uses.
  **Location.** Vol~III `chapters/theory/quantum_chiral_algebras.tex:942`
  (three-dualities proposition), line 1524 (K3 defect proof);
  `chapters/connections/bar_cobar_bridge.tex:812` (KPZ dictionary).
  **Related.** AP14 (Swiss-cheese Koszul locus); AP32 (scope
  contradictions between volumes); AP40 (status-tag drift);
  Pattern 236 (ambient qualifier discipline).
  **Rule.** Any cross-volume citation of Theorem~B that does not
  specify "raw direct-sum" vs "weight-completed" vs "finite-weight"
  is forbidden. "On the Koszul locus" without further qualifier is
  too narrow for all class-C / class-M applications.

- **AP907 -- Label conflation: `thm:e3-identification` vs
  `thm:e3-identification-chain-level-associator-fixed`.**
  Vol~I has two theorems with overlapping but distinct scopes:
  (a) `thm:e3-identification` (`chapters/theory/en_koszul_duality.tex:5346`):
  identification of $\Ethree$-deformation families
  $Z^{\mathrm{der}}_{\mathrm{ch}}(V_k(\fg)) \simeq \cA^\lambda$ as
  formal deformation families, proved only for $\fg$ simple
  (where $\dim H^3(\fg) = 1$, the key input to the order-by-order
  uniqueness argument);
  (b) `thm:e3-identification-chain-level-associator-fixed`
  (`chapters/theory/e3_identification_chain_level_platonic.tex:411`):
  chain-level identification $Z^{\mathrm{der}}_{\mathrm{ch}}(V_k(\fg))
  \to \cA^\lambda$ at the KZ associator, proved for $\fg$ simple,
  reductive, or solvable; $\fgl_N$, semisimple, nilpotent Heisenberg
  are subcases of reductive and solvable.
  **Violation.** Vol~III `chapters/examples/quantum_group_reps.tex:1146`
  cites "Vol~I `thm:e3-identification` (PROVED for simple $\frakg$,
  $\fgl_N$, semisimple, reductive, nilpotent Heisenberg, solvable)".
  The scope listed is the scope of (b), but the label cited is (a),
  which only covers simple $\fg$.
  **Counter.** Cite `thm:e3-identification-chain-level-associator-fixed`
  when invoking the extended reductive / solvable scope, and cite
  `thm:e3-identification` only for simple $\fg$ (formal-deformation
  families). The two labels are NOT interchangeable: (a) is about
  isomorphism classes of formal deformation functors; (b) is about
  existence of a specific chain map at a fixed associator.
  **Location.** Vol~III `chapters/examples/quantum_group_reps.tex:1146`.
  **Related.** AP7 (scope contradictions); AP40 (status-tag drift);
  AP106 (placeholder-label discipline); Pattern 269 (lane separation).
  **Rule.** Citing `thm:e3-identification` across volumes at the
  extended reductive / solvable scope is forbidden; re-cite to
  `thm:e3-identification-chain-level-associator-fixed`.

- **AP908 -- $(\mathsf{SC}^{\mathrm{ch,top}})^!$-cobar contracting
  homotopy confused with weight-completed Positselski
  (cross-volume dependency with wrong Vol~I anchor).**
  Vol~II `chapters/theory/chiral_higher_deligne.tex:495-508` declares
  clause~(2) of Theorem~\ref{thm:chd-universality} conditional on the
  "$(\SCchtop)^!$-cobar contracting homotopy" and attributes this
  to "the $\SCchtop$-variant of the chiral Positselski co-contra
  correspondence (Vol.~I, chapter on Theorem~B)". The phantom label
  `thm:chiral-positselski-7-2` recorded in the Vol~I preface is
  described as "a placeholder, not a proved body". In reality
  `thm:chiral-positselski-7-2` is an alias (via `\phantomsection
  \label{thm:chiral-positselski-7-2}`) inside the proved body of
  Vol~I's `thm:chiral-positselski-weight-completed` at
  `chapters/theory/theorem_B_scope_platonic.tex:248`. But the
  weight-completed chiral Positselski statement is on the SINGLE-COLOUR
  $E_1^{\mathrm{ch}}$ cooperad, NOT the two-coloured $(\SCchtop)^!$
  cooperad. The two-coloured lift is a different theorem, genuinely
  not yet proved at chain level. Vol~II's conditionality is correct
  in substance (the $(\SCchtop)^!$ lift is open) but the attribution
  (phantom label in Vol~I preface) is incorrect in form (the Vol~I
  label aliases a real proved body, for the single-colour statement).
  **Counter.** Rewrite Vol~II's remark as: "clause~(2) is conditional
  on the two-coloured $(\SCchtop)^!$-cobar contracting homotopy,
  which is a genuine open problem. The single-colour weight-completed
  chiral Positselski statement proved in Vol~I
  (`thm:chiral-positselski-weight-completed`,
  `chapters/theory/theorem_B_scope_platonic.tex:248`) does not
  immediately lift to the two-coloured cooperad setting. The Vol~I
  preface alias `thm:chiral-positselski-7-2` refers to the
  single-colour proved body, not the two-coloured lift."
  **Location.** Vol~II `chapters/theory/chiral_higher_deligne.tex`
  Remark~\ref{rem:chd-universality-conditional} lines 495-508.
  **Related.** AP13 (circular proof dependency); AP106 (placeholder
  labels); AP166 ($\SCchtop^!$ not self-dual); Pattern 277 ($\SCchtop$
  vs $E_3$).
  **Rule.** When a Vol~II / Vol~III theorem is conditional on a Vol~I
  result, the cited Vol~I label must point at a statement that
  actually establishes what the downstream theorem uses. Single-colour
  weight-completed chiral Positselski is NOT the same theorem as the
  two-coloured $(\SCchtop)^!$-cobar contracting homotopy.

- **AP909 -- Koszul locus equated with class~G (scope narrowing).**
  The Koszul locus is defined geometrically on
  $\overline{\mathcal A_2} \setminus H_1$ (principally polarised
  abelian surfaces minus the Humbert divisor $H_1$) at the
  cross-volume level, and intrinsically as the locus on which the
  bar coalgebra $B^{\mathrm{ch}}(\cA)$ is quadratic Koszul
  (Priddy sense; Positselski--Vishik 2000 Theorem~10.1). This
  includes class~G (free-field: Heisenberg, lattice VOAs) AND
  class~L (affine Kac--Moody at generic level, proved via
  Francis--Gaitsgory 2012 *Selecta Math.* 18). The Koszul locus is
  NOT synonymous with class~G.
  **Violation.** Vol~III `chapters/theory/quantum_chiral_algebras.tex:1524`
  writes "on the Koszul locus (class~G)", parenthetically equating
  the two. For $A_{K3} = V_{\Lambda_{K3}}$ the specific case is
  class~G (free-field, lattice VOA), so the proof succeeds; but the
  parenthetical equation "Koszul locus = class~G" narrows the scope
  incorrectly as a general identification.
  **Counter.** Write "on the Koszul locus; for $A_{K3}$ this is
  class~G (free-field, lattice VOA)" to keep the scope broad and
  label the specific case.
  **Location.** Vol~III `chapters/theory/quantum_chiral_algebras.tex:1524`.
  **Related.** AP14 (Swiss-cheese Koszul locus); AP32 (scope
  contradictions); AP906 (Theorem~B scope collapse).
  **Rule.** The Koszul locus includes class~G and class~L (and
  portions of class~C on fibres where shadow depth is 3); any
  parenthetical equation "Koszul locus = class~X" for a single
  class~X is incorrect.

- **AP910 -- Bar-cobar adjunction direction inverted, cross-volume sweep
  (Critical; cross-volume scope+direction contradiction; follow-up to
  AP905).** AP905 established the canonical direction $\Omega^{\mathrm{ch}}
  \dashv B^{\mathrm{ch}}$ and listed violations at Vol~I
  `chapters/theory/chiral_climax_platonic.tex:893,928`,
  `chapters/connections/concordance.tex:365,1410,2459`,
  `chapters/theory/introduction.tex:496,3080,3125`, and Vol~III
  `chapters/theory/quantum_chiral_algebras.tex:1524`. The 2026-04-21
  sweep surfaces six further inverted occurrences that AP905 did not
  enumerate, each of which propagates the same type error:
  (i) Vol~I `chapters/theory/introduction.tex:2264` ("The bar-cobar
  adjunction $B \dashv \Omega$ encodes boundary data") -- reader-facing
  entry point to Theorem~B;
  (ii) Vol~I `chapters/connections/holographic_codes_koszul.tex:25`
  (inverted in the chapter opening of the holographic code chapter);
  (iii) Vol~I `chapters/connections/thqg_entanglement_programme.tex:597`
  (Proposition statement on code structure with encoding map $B$ and
  decoding map $\Omega$ presented as $B \dashv \Omega$);
  (iv) Vol~I `chapters/connections/holographic_datum_master.tex:294`
  (first-face incarnation of $r_\cA(z)$);
  (v) Vol~III `chapters/theory/e1_chiral_algebras.tex:485` ($B^{\mathrm{ord}}
  \dashv \Omega^{\mathrm{ord}}$);
  (vi) Vol~III `chapters/theory/cy_to_chiral.tex:124,3285` ($\bar{B} \dashv
  \Omega$ and $B^{\Eone} \dashv \Omega^{\Eone}$);
  (vii) Vol~III `chapters/theory/braided_factorization.tex:2044`
  ($B_{\Etwo} \dashv \Omega_{\Etwo}$); and
  (viii) Vol~III `chapters/connections/bar_cobar_bridge.tex:796` (the
  categorical S-duality prose asserts the exchange $B \dashv \Omega$).
  Two of these (i) and (viii) are particularly harmful because they
  appear in reader-facing chapter openings (introduction, S-duality
  bridge).
  **Counter.** Healed 2026-04-21 at all eight locations to the canonical
  $\Omega \dashv B$ orientation with an AP905 cross-reference and a
  first-principles one-liner (Loday--Vallette 2012 Ch.~11, Positselski
  2011 Mem.\ AMS 996, Francis--Gaitsgory 2012 \emph{Selecta Math.} 18
  Theorem~3.1) so that a reader does not need to page back to AP905 /
  bar_cobar_adjunction_curved.tex to recover the convention. The AP905
  "Rule" is in force: "$B \dashv \Omega$" in bare manuscript prose is
  forbidden going forward. The S-duality bridge at Vol~III
  `bar_cobar_bridge.tex:796` is rewritten to realise the S-duality
  exchange as the interchange of unit and counit \emph{of the same
  adjunction} under $g_s \mapsto 1/g_s$, preserving the categorical
  content while fixing the direction.
  **Location.** Fixes applied this session at all eight file:line
  locations above. Vol~II has no new violations surfaced this pass
  beyond the two instances already documented in
  `compute/audit/vol1_theorem_A_infinity_2_attack_heal_2026_04_18.md`
  and `chapters/theory/foundations_recast_draft.tex:334`.
  **Related.** AP905 (canonical direction); AP13 (circular
  dependency); Pattern 269 (lane separation); AP906 (Theorem~B scope
  collapse).
  **Rule.** Cross-volume citations of Theorem~A / Theorem~B that
  echo "bar-cobar adjunction" must display the canonical direction
  $\Omega^{\mathrm{ch}} \dashv B^{\mathrm{ch}}$ (cobar left, bar right)
  at the first occurrence in the chapter, with a one-line primary
  reference at first use. Inverted variants in manuscript prose are
  forbidden.

- **AP911 -- Theorem~H concentration restated in Vol~III without
  CY$_d$-ambient qualifier (High; follow-up to AP-CAT-7; chain-level
  vs cohomological lane-crossing).** AP-CAT-7 marked this as a
  forward-looking hole ("preventative AP; no violation yet at audit
  time"); the 2026-04-21 sweep finds four Vol~III restatements of
  Theorem~H that give the concentration as $\{0,1,2\}$ without naming
  the ambient, and therefore implicitly extend the curve-ambient
  $d=1$ Ran statement to the CY$_d$-enlarged ambient where the
  correct concentration is $\{0,1,2,d\}$ in the cohomological lane
  (and infinite at chain level for class-M inputs at $d \ge 3$ per
  AP885).
  **Precise error.** The four Vol~III locations (a)
  `chapters/connections/modular_koszul_bridge.tex:286`, (b)
  `chapters/theory/drinfeld_center.tex:729`, (c)
  `chapters/theory/hochschild_calculus.tex:12`, (d)
  `chapters/theory/cy_categories.tex:93`, all paraphrase Vol~I
  Theorem~H as "concentrated in cohomological degrees $\{0,1,2\}$"
  without indicating that Vol~I Theorem~H is a statement on the
  \emph{curve-ambient} (i.e.\ $d=1$ Ran space) and that the
  Vol~III reading is through $\Phi_d$ with $d \in \{2, 3\}$. The
  CY$_d$-enlarged ambient is \emph{not} the curve ambient: the
  concentration enlarges to $\{0,1,2,d\}$ at the cohomological lane
  by AP-CY46 (Phi-d lane enlargement) and further, at chain level,
  is infinite for class-M inputs at $d \ge 3$ by AP885.
  **Counter.** Healed 2026-04-21 at all four locations by (a)
  annotating "on the curve-ambient reading ($d=1$ Ran ambient)"
  inside the Vol~I Theorem~H citation; (b) adding a one-line
  $\Phi_d$-bridge clause "Under $\Phi_d$ with $d \ge 2$, the
  concentration enlarges to $\{0,1,2,d\}$ in the cohomological
  lane (AP-CAT-7; chain-level is infinite for class-M inputs at
  $d \ge 3$)".
  **Location.** Healed 2026-04-21 at Vol~III
  `chapters/connections/modular_koszul_bridge.tex:286`,
  `chapters/theory/drinfeld_center.tex:729`,
  `chapters/theory/hochschild_calculus.tex:12`, and
  `chapters/theory/cy_categories.tex:93`. Vol~III
  `chapters/connections/bar_cobar_bridge.tex:221` already carries
  the correct "ordered bar" ambient qualifier and needs no fix.
  **Related.** AP-CAT-7 (forward-looking hole, preventative); AP-CY46
  (Phi-d lane enlargement / no CY$_4$ Yangian);
  AP885 (class-M $E_3$-bar $6^g$ overcount); Pattern 236 (ambient
  qualifier discipline); Pattern 269 (chain vs $(\infty,1)$ lane).
  **Rule.** Every Vol~III / Vol~II restatement of Theorem~H must
  carry an ambient qualifier. Two allowed phrasings: (a) "on the
  curve-ambient reading / $d=1$ Ran ambient: $\{0,1,2\}$" when
  citing Vol~I directly without $\Phi_d$-bridge; (b) "under $\Phi_d$
  with $d \ge 2$: $\{0,1,2,d\}$ in the cohomological lane; infinite
  at chain level for class-M inputs at $d \ge 3$" when citing
  through the Phi-d bridge. Bare "$\{0,1,2\}$" in Vol~III without
  ambient label is forbidden.

- **AP912 -- Quintic $K_0^{\mathrm{num}}$ rank correction (Agent aa123,
  Waves 1--5 finding).** Wave 1--5 attack-heal surfaced a rank error
  for $K_0^{\mathrm{num}}(D^b\mathrm{Coh}(X_5))$ at the Fermat quintic
  $X_5 \subset \mathbb P^4$: the numerical Grothendieck group has
  rank $4 = h^{1,1}(X_5) + h^{1,2}(X_5) - h^{1,2}(X_5) + 2
  = 1 + 1 + 1 + 1$ (Picard + two mid-Hodge + structure sheaf); NOT the
  naive $h^{0,0} + h^{1,1} = 2$. **Rule**: every quintic $K_0^{\mathrm{num}}$
  invocation must cite the Orlov 2009 \emph{Progr. Math.} 270 numerical
  Grothendieck rank identification $\mathrm{rk}\,K_0^{\mathrm{num}}(X_d)
  = \sum_p h^{p,p}(X_d) + 2 h^{1,d-1}(X_d)$ for hypersurfaces
  $X_d \subset \mathbb P^{d+1}$ of Calabi--Yau type; substitute and
  evaluate at $d = 3$ before propagating. For quintic: $h^{0,0}=h^{1,1}
  =h^{2,2}=h^{3,3}=1$ and $h^{1,2}=h^{2,1}=101$, giving
  $\mathrm{rk}\,K_0^{\mathrm{num}}(X_5) = 4$. **Counter**: three
  independent verification paths — (i) Orlov 2009 direct; (ii) HKR
  + Mukai lattice signature $\mathrm{sgn} = (1,1,h^{1,1},h^{1,1})
  = (1,1,1,1)$; (iii) Bridgeland 2007 Thm 1.2 rank via $\dim\mathrm{Stab}
  = 2 \cdot \mathrm{rk}\,K_0^{\mathrm{num}}$ followed by the Bayer--Macri
  quintic $\dim \mathrm{Stab}(X_5) = 8$. **Related**: AP150 (agent
  confabulation of multi-step structures — verify arrow by arrow);
  AP155 ("new invariant" overclaiming); AP226 ($K_0$ requires integer,
  not complex scalar); Pattern 295 (W25) in `notes/first_principles_cache_comprehensive.md`.
  **Scope**: applies at every quintic $K_0^{\mathrm{num}}$ citation
  site across all three volumes.

- **AP913 -- Hilb$^n$(K3) partition-formula factor-of-two correction
  (Agent a42182, Waves 1--5 finding).** The partition-function generating
  series $\sum_{n\ge 0} \chi(\mathrm{Hilb}^n(\mathrm{K3})) q^n$ has the
  Göttsche identification $\prod_{m\ge 1}(1-q^m)^{-24} = q/\Delta(q)
  = 1/\eta(q)^{24}\cdot q$ with Euler characteristic $\chi(\mathrm{K3})
  = 24$; the naive Hilb$^n$ transfer omits the $\eta(q)^{-24}$
  $q$-scaling that distinguishes the Göttsche generating series
  (prefactor $q$) from the pure MacMahon series (no prefactor). At
  $n=4$: $\chi(\mathrm{Hilb}^4(\mathrm{K3})) = p_{24}(4) = 108\,600$
  via Göttsche 1990 Math Ann 286 Thm 0.2, matching the coefficient of
  $q^4$ in $\prod_m(1-q^m)^{-24}$. **Rule**: every Hilb$^n$(K3)
  Euler-characteristic formula carries explicit prefactor convention.
  **Counter**: substitute $n=1,2,3,4$ and cross-check against the
  established $\chi(\mathrm{Hilb}^n(\mathrm{K3})) = p_{24}(n)$
  sequence $\{24, 324, 3200, 25650, 176256, 108\,8280,\ldots\}$; any
  deviation at low $n$ signals prefactor omission or Göttsche/MacMahon
  confusion. **Related**: AP135 ($q$-expansion $p_{-r}(n)$ for $1/\eta^r$
  partition numbers); AP884/AP341 ($c_3 = 176\,256 = p_{24}(5)$ Göttsche
  trap — same number appears at $n=5$ in Göttsche, \emph{not} as
  Borcherds coefficient); AP885 (class-M $E_3$-bar factor-of-two
  Mukai double-count at K3 $\times E$). **Scope**: applies to every
  Hilb$^n$(K3), Hilb$^n$(K3 $\times E$), and symmetric-product
  elliptic-genus citation across volumes.

- **AP914 -- Oberdieck factor-of-two at class-M $E_3$-bar invocation
  (Agent a03831, Waves 1--5 finding).** AP885 documents the class-M
  $E_3$-bar factor-of-two Mukai double-count (healed: $e_3^{\mathrm{K3}
  \times E} = 2\mathrm{Vol}(E)(2\pi i)^3$, NOT $4\mathrm{Vol}(E)(2\pi
  i)^3$). The Wave 1--5 sweep surfaced additional inscription sites
  where the pre-heal value $4\,\mathrm{Vol}(E)(2\pi i)^3$ survives
  without an AP885 retraction annotation, creating latent-AP violations
  when a future reader cross-references the canonical value against
  uncorrected sites. **Rule**: every K3 $\times E$ triangle-cycle /
  class-M $E_3$-bar / Oberdieck reduced-DT citation must pin the
  Mukai-pairing normalisation $\chi(\mathcal O_{\mathrm{K3}}) = 2$
  (\emph{not} $4$); the Oberdieck 2018 Geom \& Topol 22 Thm 2 reduced
  DT generating function $= -\Phi_{10}^{-1}$ uses the single-factor
  Mukai pairing, and the triangle-cycle value is $2\mathrm{Vol}(E)(2\pi
  i)^3$. **Counter**: grep all three volumes for `4\,\\?mathrm\{Vol\}\(E\)`
  and verify each site either (a) cites AP885/AP914 with retraction
  context; or (b) carries the corrected value $2\mathrm{Vol}(E)(2\pi i)^3$.
  **Related**: AP885 (primary class-M $E_3$-bar correction); AP-CY46
  ($\Phi_d$ lane enlargement); entry 405 (AP885) + entry 159 (W20
  MNOP-$\chi_3$) in `appendices/first_principles_cache.md`.
  **Scope**: cross-volume; affects all K3 $\times E$ DT / Hochschild /
  Borcherds citations.

- **AP915 -- $S_5$ / PRS conflation in Pope--Romans--Shen scalar row
  (Agent abcff0, Waves 1--5 finding).** The shadow coefficient $S_5$
  of the Virasoro archetype (canonical: $S_5 = -48/[c^2(5c+22)]$)
  collides semantically with the Pope--Romans--Shen 1990 \emph{PLB} 236
  Eq (2.7--2.9) scalar conversion $W_k = e_k + (-c/22)\Lambda_Z$ at
  $k=5$: the Wave 1--5 swarm observed draft sites conflating the
  shadow-tower ``fifth coefficient'' $S_5$ with the PRS ``fifth primary
  scalar'' $\Lambda_{Z,5}$ in the $W$-primary-to-Virasoro-composite
  basis conversion. These are semantically distinct invariants on the
  same weight-5 quasi-primary: $S_5$ is the Arnold-cycle shadow
  evaluation; the PRS scalar is the projection coefficient between two
  Wang-1998-three-leg-uniqueness bases. **Rule**: every shadow $S_r$
  citation and every PRS $-c/22$ scalar citation carries its distinct
  index (``$S_r$, shadow depth $r$'' vs ``PRS scalar $-c/22$ at weight
  $k$''); two-basis conversion equations must name Virasoro-composite
  vs $W$-primary bases explicitly per AP-CAT-4. **Counter**: grep for
  bare $S_5$ in any $W$-primary / PRS-conversion paragraph; verify the
  shadow / projection distinction is named. **Related**: AP178 ($S_4$
  large-$c$ factor-5 asymptotic check); AP-CAT-4 (two-basis $e_4$
  non-contradiction); canonical-values row 51 ($e_4$ at $c=-214$);
  AP38 (DVV != Eichler--Zagier convention); AP39 ($\kappa \ne S_2$ for
  non-Virasoro). **Scope**: all shadow-tower + PRS-basis prose in Vol~I
  and the Vol~III BKM crown row.

- **AP916 -- AM biconditional (Agent a154, Waves 1--5 finding).**
  ``AM'' = Atiyah--Manin pairing-style biconditional routinely appears
  in adjunction / duality chapters in the form ``$X \ne 0 \iff Y \ne
  0$''. The Wave 1--5 swarm surfaced cases where the forward direction
  (``$X \ne 0 \Rightarrow Y \ne 0$'') is proved but the reverse
  direction (``$Y \ne 0 \Rightarrow X \ne 0$'') is used as if proved
  when only the forward is established. **Rule** (strengthening AP36:
  ``implies proved, iff claimed''): every biconditional in an
  Atiyah--Manin / pairing / duality context carries a two-line template
  ``\emph{Forward}: \ldots\ [anchor/citation]. \emph{Reverse}: \ldots\
  [anchor/citation]'' before the $\iff$ symbol is written. Unstated
  reverse direction downgrades to $\Rightarrow$ + remark. **Counter**:
  grep chapters for ``$\iff$'' / ``if and only if'' in pairing or
  duality contexts; audit every site for documented reverse direction.
  **Related**: AP36 (implies proved, iff claimed); AP59 (three invariant
  distinctions); FM20 (iff-drift counter template). **Scope**: Vol~I
  Theorem~B / Theorem~C biconditional clauses; Vol~III $\Phi_d$ pairing
  biconditionals.

- **AP917 -- L2-theorem-upgrade naming scope (Agent ab57, Waves 1--5
  finding).** When a Lemma 2 (``L2'') is upgraded to a Theorem status,
  the upgrade must carry explicit scope naming: (a) the theorem-level
  label with the former lemma label preserved as a `\phantomsection` +
  back-reference for citing sites (AP127); (b) the proof body must
  either be inlined (promoting L2 from sketch to proof) or cited
  unambiguously (L2 as lemma with named proof elsewhere); (c) every
  consumer `\ref{lem:...}` must be retargeted to the `\ref{thm:...}`
  form in the SAME commit (AP149). Wave 1--5 surfaced an upgrade where
  the theorem-level label was installed but three consumer `\ref` sites
  retained the stale `lem:` prefix. **Rule**: L-to-T upgrades are
  atomic operations touching environment, label, and every downstream
  `\ref` in one tool-call batch. **Counter**: after any `lemma` $\to$
  `theorem` environment change, run `grep -rn "\\\\ref{lem:${LABEL}}"`
  and retarget; also run the AP124 duplicate-label check to prevent
  parallel `lem:` + `thm:` labels co-existing. **Related**: AP124
  (label discipline); AP127 (label migration stubs); AP149 (resolution
  propagation); AP40 (environment-tag match). **Scope**: any environment
  upgrade across all three volumes.

- **AP-CAT-9 -- Wave 1--5 agent-identifier-only entries without AP
  anchor (Beilinson meta-bookkeeping, 2026-04-22).** The 2026-04-21
  W26 AP-audit (AP-CAT-1 through AP-CAT-8) noted that ``every wave's
  final summary must list AP-entry anchors for each finding; findings
  without an AP anchor get one inscribed at the session close''. Waves
  1--5 produced six agent-level findings identified by agent-hash
  (aa123, a42182, a03831, abcff0, a154, ab57) that had been applied
  in-content but lacked a back-anchor row in the AP catalogue. This
  AP-CAT-9 closes that gap by installing AP912--AP917 as the explicit
  back-anchors for those six findings. **Rule**: every Wave $N$
  finding identified by agent-hash must receive an AP entry at session
  close; the agent-hash $\to$ AP$n$ mapping table lives in this AP-CAT
  registry. **Mapping** (2026-04-22):
  (a) aa123 — quintic $K_0^{\mathrm{num}}$ rank correction $\to$ AP912;
  (b) a42182 — Hilb$^4$(K3) partition-formula correction $\to$ AP913;
  (c) a03831 — Oberdieck factor-of-2 at class-M $E_3$-bar $\to$ AP914;
  (d) abcff0 — $S_5$ / PRS scalar conflation $\to$ AP915;
  (e) a154 — AM biconditional forward-vs-reverse scope $\to$ AP916;
  (f) ab57 — L2-theorem-upgrade atomic-rename scope $\to$ AP917.
  **Counter**: at every wave's close, grep session transcript for
  ``Agent a[0-9a-f]+`` agent-identifiers with content-layer edits;
  verify each has an AP-level back-anchor. Missing anchor = latent
  AP-CAT-9 violation; inscribe anchor in next session.

- **AP-CAT-10 -- Forward-looking hole: Wave 6+ expected territories.**
  Preventative. The Wave 1--5 attack targeted (a) quintic CY-3 numerical
  Grothendieck group; (b) Hilb$^n$(K3) Göttsche counts; (c) Oberdieck
  reduced-DT / $\Phi_{10}^{-1}$ normalisation; (d) shadow-tower /
  PRS-basis scalar distinction; (e) adjunction biconditional scope;
  (f) lemma-theorem upgrade atomicity. Wave 6+ is likely to target
  adjacent territory where the same confusion mechanisms re-appear:
  (W6a) Hilb$^n$(K3 $\times E$) triple-base DT (AP-CY46 $\Phi_d$ lane
  + AP914 Mukai factor + bi-based Ran per AP886); (W6b) non-principal
  $\mathcal W$-algebras at hook-type / rectangular (AP76 $Y_{1,1,1}$
  $\kappa=\Psi$ + non-principal-W hook-type correspondence); (W6c)
  Yetter--Drinfeld tower weights $\delta^{(n)}$ at $n \ge 10$
  extrapolation (entry 419 / AP899 at $n=7,8,9$ with $\lfloor n/2\rfloor
  +1$); (W6d) Leech-Conway universal ratio closure (OP3/OP12); (W6e)
  $(\infty,2)$-adjunction AR / obstruction-tower bridging (W21.2 +
  W22.1 + compatible-dual-readings row 10); (W6f) K3 $\times E$
  Hecke-eigensheaf condition (Selmer dim 1, canonical-values rows
  34--35). Each territory has known failure modes documented; Wave 6+
  agents should (i) grep the canonical-values table first; (ii) check
  compatible-dual-readings before attempting retraction; (iii) cite
  the AP chain bounding the failure mode in the open-problem index.
  **Rule**: Wave 6+ attack briefs cite AP-CAT-10 at launch, naming
  which territory the agent targets and which AP-chain bounds it.
  **Counter**: preventative — no violation yet.

- **AP-CAT-11 -- Cache top-15 reorganisation acknowledgment,
  2026-04-22.** The hook `.claude/hooks/beilinson-gate.sh` §3
  (``POST-INSCRIPTION CACHE SWEEP'') checks 15 cache patterns on every
  theorem-bearing edit. The cache-comprehensive file organises $\sim$250
  patterns without an explicit L1/L2/L3 tiering; the hook encodes the
  effective top-15 implicitly through the `CACHE#1`--`CACHE#15` grep
  sequence. This AP-CAT-11 formalises the tiering: **L1 (always
  checked, $n=1$--$7$)** specific/general, bare $\kappa$/Hochschild,
  native/derived $E_n$, construction/narration, algebra/coalgebra,
  CoHA-vs-vertex, Drinfeld-centre-vs-averaging. **L2 (inscription-gated,
  $n=8$--$11$)** ambient-qualifier discipline (Pattern 236), scope
  tags (AP-UNIFORM-WEIGHT-TAG), two-$\hbar$ bridge (AP151), $\kappa$
  cross-family qualification (AP-KAPPA). **L3 (theorem-context,
  $n=12$--$15$)** chain-vs-cohom class-M $E_3$ bar (AP-CY21/38),
  part-whole $\{b_k,B^{(2)}\}$ (AP-CY34), $\kappa_{\mathrm{BKM}} =
  \kappa_{\mathrm{ch}} + \chi(\mathcal O_{\mathrm{fib}})$ coincidence
  (AP-CY37, canonical row 56), Künneth-multiplicative
  $\kappa_{\mathrm{cat}}$ (canonical row 63). The hook encodes these
  in the current §3 ordering; no hook-code update needed. **Rule**:
  any future cache additions are tagged L1/L2/L3 at inscription;
  promotions/demotions between tiers documented here.

- **AP-XV1 -- Cross-volume phantom cascade (2026-04-22).** Consumer
  prose references targets `V1-thm:*`, `V2-prop:*`, `V3-thm:*` (or the
  older `thm:*-vol1` / `thm:*-vol2` bridge-naming convention), but the
  cross-volume alias block in the host volume's `main.tex` lacks the
  corresponding `\phantomsection\label{...}` anchor, OR the target
  label in the foreign volume was renamed, type-changed (thm vs prop
  vs cor), or case-changed (e.g. Vol~III `modular_trace.tex` uses
  lowercase `thm:modtr-plancherel-kerlerlyubashenko` while Vol~I
  consumers wrote camelCase `V3-thm:modtr-Plancherel-KerlerLyubashenko`).
  LaTeX labels ARE case-sensitive, so camelCase/lowercase drift yields
  hard `[?]` renders that the author reads as "missing theorem".
  **Mechanism**. Three sub-patterns: (i) *missing alias* — the foreign
  volume's label exists, but the host volume's `\phantomsection` stub
  block was never updated when the consumer prose was inscribed;
  (ii) *type mismatch* — consumer wrote `V2-thm:X` but foreign volume
  inscribed `prop:X` (or `cor:X`, etc.); (iii) *case drift* — consumer
  wrote `V3-thm:Modtr-Plancherel-ML-limit` (camelCase) while Vol~III
  inscribed `thm:modtr-plancherel-ml-limit` (lowercase). The *genuine*
  phantom — target does not exist anywhere — is rare; at the 2026-04-22
  sweep only 1 of 107 `V[123]-` references in Vols~I+II was a true
  missing-target phantom (`V1-constr:platonic-package`).
  **Regex trigger**: `(V[123]-(thm|prop|lem|cor|def|rem|sec|ch|eq|conj|constr|comp|conv|princ):|[^[:alnum:]](thm|prop|lem|cor|def):[^}]+-vol[123])`.
  **Protocol** (5-step heal):
  (1) Enumerate all `\ref{V[123]-*}` used in each volume:
  ``grep -rhE '\\(ref|cref|Cref|autoref|eqref)\{V[123]-[^}]+\}' --include="*.tex" | grep -oE 'V[123]-[^}]+' | sort -u``.
  (2) Enumerate all `\phantomsection\label{V[123]-*}` aliases declared
  in that volume's `main.tex`.
  (3) Compute set difference (used minus aliased) $=$ alias-missing
  phantoms.
  (4) For each alias-missing name, grep the foreign volume for the
  bare label; if present, insert a `\phantomsection\label{V?-...}`
  stub in the host volume's cross-volume alias block (near the
  existing `V1-*` or `V2-*` block); if the foreign volume's label
  is case-drifted or type-drifted, record the canonical target in
  the alias comment and either (a) adjust the consumer to the real
  canonical name, or (b) add a case/type-drift alias stub (the alias
  route is lossless when the consumer prose already attributes the
  theorem to the correct volume with text like "Vol~III Theorem~").
  (5) For each genuinely missing target (label not found in any
  volume), the heal is to either inscribe the missing theorem in
  the correct volume or to retarget the consumer to a real sibling
  (e.g. `V1-thm:topologization-tower` $\to$ local
  `thm:e-infinity-topologization-ladder` alias at
  preface.tex:5487).
  **Findings from the 2026-04-22 sweep.**
  *Vol~I* had 19 alias-missing phantoms; all healed by adding
  aliases to `/Users/raeez/chiral-bar-cobar/main.tex` cross-volume
  block. Type-drift fixed: `V2-thm:modbootstrap-cross-Tmatrix` and
  `V2-thm:modbootstrap-cross-block-count` are `prop:` in Vol~II;
  aliases installed with prop-level anchors and consumer at
  `chapters/theory/shadow_tower_quadrichotomy_platonic.tex:1842`
  retyped to Proposition. Case-drift aliases installed for four
  Vol~III `modtr-*` targets (Plancherel-KL, Plancherel-ML-limit,
  Zsigma2-Phi10, ZT2-HH0). *Vol~II* had 88 alias-missing phantoms;
  all healed by a single insertion into
  `/Users/raeez/chiral-bar-cobar-vol2/main.tex` after the existing
  V1- alias block (lines 1017--1113 after heal). *Vol~III* had 0
  phantoms. *Residual after sweep*: 0.
  **Counter**: `post-inscription` whenever a consumer writes
  `\ref{V[123]-*}` to a new target, verify (i) the target label
  exists in the foreign volume, (ii) a `\phantomsection\label{V?-*}`
  alias exists in the host volume's cross-volume block, (iii) the
  type (thm/prop/cor/etc.) in the alias matches the foreign
  volume's inscribed type. Running sweep: the three `comm -23`
  pipeline above produces a zero-diff when clean. Re-run at every
  wave close.
  **Note**. Agent aa72's two flagged phantoms ---
  `thm:bar-cobar-adjunction-vol1` and `V1-thm:topologization-tower`
  --- turned out to be already-healed (the latter has a local
  phantomsection alias at `preface.tex:5487`; the former has zero
  live consumers). The sweep found far larger real-phantom surface
  area (107 refs, 107 alias-missing anchors) than the two-phantom
  report suggested.

- **AP918 -- Curved-morphism primitive $f_0$ vs curved-algebra curvature
  $\ell_0$ (Agent a3c1fd, Wave 5 residual).** The Merkulov curved-algebra
  / curved-morphism formalism carries two distinct curvature-type
  objects living in *different* graded slots and satisfying *different*
  identities: (i) $\ell_0 \in A^2$ --- the *curved-algebra curvature*,
  a degree-$2$ element of the algebra with $[\ell_0, \ell_0]/2 = 0$
  required for curved-$L_\infty$ (curvature of the differential);
  (ii) $f_0 \in \mathrm{Hom}^0(k, A)$ --- the *curved-morphism
  primitive*, a degree-$0$ element specifying where the curved morphism
  sends the basepoint, with Maurer--Cartan defect
  $D f_0 + \tfrac12[f_0, f_0] = f^\ast \ell_0^{\mathrm{src}} - \ell_0^{\mathrm{tgt}}$
  (cf. Merkulov 2019, wheel-free formality and curved
  $L_\infty$-morphisms). Conflating the two produces sign-ambiguity in
  the Maurer--Cartan obstruction tower and the wrong curvature class:
  writing $\ell_0$ where $f_0$ is meant sends a morphism-primitive
  into the algebra-curvature slot, silently shifting the obstruction
  identity by one graded-commutator term. **Wrong claim**:
  ``$\ell_0 = f_0$'' or ``$\ell_0$ is the image of the basepoint
  under the curved morphism''. **Precise error**: $\ell_0$ lives in
  the algebra in degree $+2$ and is part of the algebra data; $f_0$
  lives in the morphism-space in degree $0$ and is part of the
  morphism data. They are bridged by the MC-identity above, not
  equated. **Counter**: every curved-$L_\infty$ or curved-$A_\infty$
  statement that mentions both ``curvature'' and ``primitive'' lists
  which slot each lives in ($\ell_0 \in A^2$ vs $f_0 \in \mathrm{Hom}^0(k,A)$)
  and names the bridging identity. **Related**: AP918; AP76
  ($Y_{1,1,1}$ curvature scope); AP151 (two-$\hbar$ bridge); Merkulov
  2019.

- **AP919 -- Bucket-size discipline: four-bucket vs five-bucket
  archetype enumeration (Agent aecd21, Wave 5 residual).**
  Standard-landscape enumeration routinely presents the dichotomy as
  *four* archetypes $\mathsf G/\mathsf L/\mathsf C/\mathsf M$
  (Heisenberg, affine Kac--Moody, $\beta\gamma$, Virasoro) with
  shadow-tower depth $r_{\max} \in \{2,3,4,\infty\}$ (per `CLAUDE.md`
  ``One dichotomy'' paragraph). The *five*-archetype view adds
  $\mathsf B$ and reports the ceiling
  $\kappa + \kappa^! \in \{0,8,13,250/3,98/3\}$ (Theorem~C
  family-dependent ceiling) vs the classical $\mathsf G/\mathsf L/\mathsf C/\mathsf M$
  subset $\{0,13,250/3,98/3\}$. AP917 covers *L-to-T* lemma/theorem
  upgrade atomicity but does NOT cover bucket-size selection. Agent
  aecd21 flagged a statement quantifying over ``the four archetypes''
  that in context required the five-bucket ceiling (the $\mathsf B$
  contribution carrying the $8$). **Wrong claim**: ``standard
  landscape = four archetypes, always''. **Precise error**: the
  four-bucket view is the *shadow-depth* classification
  ($\mathsf{G}/\mathsf{L}/\mathsf{C}/\mathsf{M}$); the five-bucket
  view is the $\kappa + \kappa^!$ *ceiling* classification
  ($\mathsf{G}/\mathsf{L}/\mathsf{C}/\mathsf{M}/\mathsf{B}$). Picking
  the wrong bucket count collapses the $\mathsf B$ contribution or
  overcounts when only the classical subset is in scope. **Counter**:
  every archetype-quantification in a proof or table cites (i) which
  classification (shadow-depth vs $\kappa+\kappa^!$ ceiling); (ii)
  the canonical $r_{\max}$ / $\kappa+\kappa^!$ value per archetype;
  (iii) the bucket-size rationale (four vs five). Appendix tables
  default to five-bucket with a classical-subset column marked.
  **Related**: AP919; AP917 (L-to-T atomicity); AP37
  ($\kappa \neq S_2$); landscape\_census.tex canonical rows;
  Theorem~C.

- **AP920 -- Sense-count discipline: anchor T20 four-vs-three-sense
  drift (Agent a275f6, Wave 5 residual).** The $r$-matrix, the
  $\kappa$-pairing, the bar differential, and the Koszul-duality
  inversion each admit several precise meanings depending on *in
  which ambient* (chain vs $(\infty,1)$, chiral vs topological vs
  categorical, specific vs uniform) the statement lives. Agent
  a275f6's finding at anchor T20 (a four-sense vs three-sense drift
  at a $\kappa$-sum claim) exposed a general discipline gap: *the
  number of meaningful senses of a named symbol or phrase must be
  stated when the count is not one*. A claim ``in all four senses,
  $X$'' that is really tri-sensed (because one sense is vacuous or
  degenerate) is a type error --- the wrong cardinality of senses
  produces a wrong universal-type claim. **Wrong claim**: asserting
  ``in all $n$ senses of $X$'' for a wrong $n$ (over- or
  under-counting senses). **Precise error**: the four relevant senses
  of $\kappa$-pairing in class~$\mathsf M$ are
  (chain, cohomological, chiral, categorical); but the four senses of
  $r$-matrix in class~$\mathsf M$ are
  (uniform-weight, single-generator, OPE-pole, deformation) ---
  different tetrad. Collapsing two distinct tetrads to a single
  ``four senses of the Virasoro object'' is a cardinality-type
  error. **Counter**: whenever a sense-count $\ge 2$ is stated, list
  each sense explicitly ($n$ bullets) at or within ten lines of the
  statement. Cross-volume grep at wave close for ``all $n$ senses''
  patterns; verify $n$-tuple is named. **Related**: AP920; AP37;
  Pattern~236 (ambient qualifier); AP-CAT-11 (L1/L2/L3 cache
  tiering); anchor T20.

- **AP921 -- Phantom file-path citation (Agent a1371d, Wave 5
  residual, anchor `lattice_foundations.tex:5991`).** Agent a1371d
  cited `lattice_foundations.tex:5991` as the locus of a specific
  lattice-VA identity; the file exists but does not extend to line
  $5991$. Cause: citation fabricated from a plausible-seeming
  file/line pair without Read verification. This is the file-path
  analogue of AP35 (``false proof / true theorem''): the *claim* may
  be true (the lattice identity does hold somewhere) but the
  *citation* is fabricated. Downstream damage: an auditor grepping
  the cited line finds nothing and flags either a missing inscription
  (which may already exist elsewhere) or a false positive, wasting
  audit cycles. **Wrong claim**: any ``$f$:$N$'' file/line citation
  where $N$ exceeds the file length, or points to a non-matching
  line content. **Precise error**: citation emitted without a Read-tool
  verification pass; model plausibility substitutes for ground truth.
  **Counter**: every inscribed ``$f$:$N$'' citation in an AP entry,
  cache pattern, commit message, or summary must be verified by a
  Read or `sed -n "${N}p" f` call *in the same tool-call batch* as
  the emission. The discipline generalises Pattern~244 (HTML-entity
  leakage) to *pointer hygiene*: every pointer into the repository
  carries a verification receipt. **Counter rule** (2026-04-22):
  every agent summary naming a file:line pair must have a Read or
  grep-verification anchor in the same deliverable. **Related**:
  AP921; AP35 (false-proof / true-theorem); Pattern~244 (tool-markup
  leak); AP124 (duplicate label hygiene); agent a1371d finding.

- **AP922 -- Shadow-tower coefficient $S_r$ vs primary-norm scalar
  $\langle\Lambda_r|\Lambda_r\rangle$ general-lane conflation (Agent
  abcff0, Wave 5 residual; general lane upstream of AP915).** AP915
  healed the specific $S_5$ vs Pope--Romans--Shen $-c/22$ conflation
  at weight $5$. The *general* anti-pattern is broader: shadow-tower
  coefficients $S_r$ (from the Arnold-cycle pairing on
  $\mathrm{Conf}_n(X)$, indexed by depth $r$) and primary-norm scalars
  $\langle\Lambda_r|\Lambda_r\rangle$ (Zamolodchikov-type norms of
  Virasoro-composite quasi-primaries at weight $r$) are *distinct
  families* of scalars indexed by the *same* integer $r$ but computed
  in *different bases* (shadow vs PRS / primary). They are related by
  a basis-change matrix (Wang 1998 CMP 195 Prop~4.2 three-leg
  uniqueness for $r \le 5$) but are not equal and not trivially
  interconvertible. Specific instances: at $r=2$, $S_2 = c/2$
  coincides with the Virasoro primary norm $\langle T | T \rangle = c/2$
  by chance; at $r=4$, $S_4 = 10/[c(5c+22)]$ while
  $\langle \Lambda | \Lambda \rangle = c(5c+22)/10$ with
  $\Lambda = {:}TT{:} - (3/10)\partial^2 T$ (these are reciprocals up to
  sign, *not* equal; AP922 lane: these are distinct pieces of data, one
  does not substitute for the other); at $r=5$, AP915 conflation.
  **Wrong claim**: ``$S_r$ is the primary norm at weight $r$'' or
  ``the PRS scalar at weight $r$ is $S_r$'' in general. **Precise
  error**: shadow coefficients pair Arnold cycles against OPE
  residues; primary norms pair $\{L_{-n}, W_{-m}, \ldots\}$ PBW-basis
  vectors against each other. The two are functions on different
  domains; their coincidence at $r=2$ is a small-depth artefact.
  **Healed reading**: every $S_r$ vs primary-norm transition cites
  the basis-change matrix (Wang 1998 Prop~4.2; Pope--Romans--Shen
  1990 PLB 236 Eq~(2.7--2.9)) or downgrades to ``shadow-basis /
  primary-basis related data''. **Counter**: grep every $S_r$ in
  context $r \ge 2$ for adjacent primary-norm mention; if both are
  present, require the basis-change matrix name (Wang) or explicit
  scalar ratio. **Related**: AP922; AP915 (weight-5 specific
  instance); AP178; AP37 ($\kappa \neq S_2$); AP-CAT-4; canonical
  row 51; Zamolodchikov 1985 TMP 65.

- **AP-CAT-12 -- Wave 5 residual AP enumeration and agent-hash
  back-anchoring (2026-04-22).** Closure for Wave 5 found five
  residual patterns that AP912--AP917 did NOT cover: Merkulov
  curved-morphism primitive vs curved-algebra curvature (a3c1fd);
  bucket-size discipline four vs five (aecd21); sense-count
  discipline four vs three (a275f6, anchor T20); phantom file-path
  citation (a1371d); and shadow-tower vs primary-norm general lane
  (abcff0, upstream of AP915). AP-CAT-12 back-anchors these five as
  AP918--AP922 and commits the enumeration to the catalogue.
  **Mapping** (2026-04-22 Wave 5 residual):
  (a) a3c1fd --- curved-morphism $f_0$ vs curved-algebra $\ell_0$
  $\to$ AP918;
  (b) aecd21 --- four-vs-five bucket discipline $\to$ AP919;
  (c) a275f6 --- sense-count discipline (anchor T20) $\to$ AP920;
  (d) a1371d --- phantom file-path citation $\to$ AP921;
  (e) abcff0 --- shadow-tower vs primary-norm general lane $\to$
  AP922 (AP915 is the weight-5 specific instance; AP922 is the
  general lane).
  **Rule**: AP-CAT-9's discipline extends here: every residual
  finding from a Wave-$N$ closure that is NOT covered by the wave's
  formal deliverable receives an AP entry at the next session
  close, and an AP-CAT-$k$ back-anchor row recording the agent-hash
  $\to$ AP$n$ mapping. **Counter**: preventative --- at every wave
  close, enumerate adversarial findings not covered by the wave's
  formal AP sweep; cross-reference agent hash; inscribe AP$n$.
  **Related**: AP-CAT-12; AP-CAT-9 (Wave 1--5 formal mapping);
  AP-CAT-10 (Wave 6+ territories); AP-CAT-11 (cache-tiering);
  AP918--AP922.

- **AP-CAT-13 -- Hook regex recommendations for AP918--AP922
  (2026-04-22, advisory).** Five new hook regexes recommended (not
  yet inscribed in `.claude/hooks/beilinson-gate.sh`) to flag
  AP918--AP922 at inscription.
  (i) *AP918* (curved-morphism / curved-algebra):
  `grep -qE 'ell_0|\\\\ell_0|f_0\\s*=\\s*\\\\ell|Maurer.Cartan.curvature' && grep -qE 'curved.morphism|curved.algebra|Merkulov'` --- trigger
  when curvature + morphism both appear and the bridging identity is
  absent.
  (ii) *AP919* (bucket-size):
  `grep -qE '(four|5|five).archetype|(four|5|five).bucket|\\\\mathsf\\{[GLCMB]\\}'` --- trigger when archetype-quantifier present without both classifications named.
  (iii) *AP920* (sense-count):
  `grep -qE 'all (two|three|four|five) senses|in both senses'` ---
  trigger when sense-cardinality stated without explicit enumeration
  within ten lines.
  (iv) *AP921* (phantom file-path): no file-side regex; this is an
  agent-discipline AP. Hook cannot verify that a citation line is
  real without reading the target file. The discipline lives at the
  *summary-emission* layer rather than the edit layer. Optional
  post-commit sweep: `grep -rE '[a-z_]+\\.tex:[0-9]{4,5}' notes/` and
  validate each hit against file length.
  (v) *AP922* (shadow vs primary-norm general lane):
  `grep -qE 'S_[2-9]|S_\\{[0-9]+\\}' && grep -qE '\\\\langle.*\\\\Lambda|PRS|Pope..Romans..Shen|Zamolodchikov'` --- trigger when both shadow-coefficient and primary-norm tokens appear without Wang-1998 basis-change citation within 20 lines.
  Hook inscription is deferred pending false-positive review on one
  wave's worth of edits. **Rule**: cache comprehensiveness and hook
  coverage are tracked separately; a new AP with a recommended regex
  may ship ahead of hook inscription as long as the regex is noted
  here. **Related**: AP-CAT-13; AP918--AP922;
  `.claude/hooks/beilinson-gate.sh` §§1,3.

- **AP923 -- hCS-vs-categorical $\mathbb E_3$-factorisation presentation
  duality (Vol~I cross-volume concordance with Vol~III AP-CY160).**
  Vol~I reads the CY$_3$ constructor-side input to the bar--cobar
  machine as $\PhiFA_3(\cC) = \mathrm{HH}^\bullet(\cC,\cC)$ with
  Lurie--To\"en $\mathbb E_3$-factorisation structure, i.e.\ a single
  categorical ``black box'' with no geometric resolution exposed. The
  Vol~III wave-14 correction surfaces a second grammar: the same
  $\mathbb E_3$-factorisation algebra admits a \emph{geometric}
  presentation as the BV quantisation of holomorphic Chern--Simons on
  $\mathrm{Tot}(\Omega_X)$ (Costello 2013 \texttt{arXiv:1303.2632};
  Costello--Li 2016 \texttt{arXiv:1605.09930} \S 6). The two grammars
  agree only up to a $\mathrm{GRT}_1(\mathbb Q)$-torsor of
  Kontsevich--Tamarkin associators (Willwacher 2014 \emph{Invent.\
  Math.} 200); Costello--Li propagator picks the Kontsevich-associator
  point. **Vol~I bearing**: Theorems A (bar--cobar) and B (chiral
  Positselski) read the CY$_d$ constructor-side input from
  $\Phi_d(\cC)$ without distinguishing the two grammars, silently
  importing whichever presentation the downstream application picks.
  Load-bearing $\zeta(3)$-coefficients in the structure-constant
  lane (e.g.\ Theorem~D's all-weight free energy
  $F_g(\cA) = \kappa(\cA)\lambda_g^{\mathrm{FP}} + \delta
  F_g^{\mathrm{cross}}(\cA)$ at $g = 3$) depend on the
  $\mathrm{GRT}_1$-torsor base-point choice, which is invisible at
  iso-class level (AP924 infra). **Wrong claim**: ``hCS quantisation
  and categorical Hochschild are the same theory'' or, dually,
  ``bar--cobar input from $\Phi_d(\cC)$ is grammar-independent''.
  **Precise error**: the two grammars are quasi-isomorphic as
  $\mathbb E_3$-algebras but not equal on the nose; the structure
  constants (load-bearing for Theorem~D $g\ge 3$) live in a
  $\mathrm{GRT}_1(\mathbb Q)$-torsor of presentations, not at a
  canonical point. **Counter**: every bar--cobar use of
  $\PhiFA_3(\cC)$ must name the grammar (hCS geometric resolution or
  categorical Hochschild with Lurie--To\"en structure) and cite the
  Costello--Li propagator or Willwacher torsor point where structure
  constants are load-bearing. Iso-class statements do not require
  the distinction; parametrised statements do. **Related**: AP923;
  Vol~III AP-CY160 (originating wave-14 entry); AP275 (scope
  qualifier without obstruction statement); AP920 (sense-count
  discipline); Costello--Li 2016; Willwacher 2014.

- **AP924 -- Iso-class vs parametrised Kontsevich--Tamarkin formality
  (Vol~I cross-volume concordance with Vol~III AP-CY161).** The
  Vol~I bar--cobar adjunction $\Omega^{\mathrm{ch}} \dashv
  B^{\mathrm{ch}}$ at chain level (Loday--Vallette 2012 Thm~11.3.3;
  Francis--Gaitsgory 2012 \emph{Selecta Math.} 18 Thm~3.1) relies on
  a choice of Drinfeld associator in the formality isomorphism
  $T_\mathrm{poly} \simeq D_\mathrm{poly}$. Stating ``the
  bar--cobar adjunction is canonical because KT formality is
  contractible in choices'' conflates two cardinalities: at
  \emph{iso-class} level the space of formality isomorphisms is
  indeed contractible (any two quantisations are
  quasi-isomorphic), but the \emph{parametrised} space of
  quasi-isomorphisms (fat formality morphisms) carries a free
  $\mathrm{GRT}_1(\mathbb Q)$-action with
  $\pi_0 = \mathrm{GRT}_1(\mathbb Q)$-torsor (Willwacher 2014
  \emph{Invent.\ Math.} 200 Thm~1.2). **Vol~I bearing**: at
  iso-class level Theorem~A's adjunction statement is canonical; at
  structure-constant level (Theorem~D's $\lambda_g$-coefficients,
  Theorem~C's $\kappa + \kappa^! \in \{0, 8, 13, 250/3, 98/3\}$
  landmark ceiling) the $\mathrm{GRT}_1$-torsor base-point is
  load-bearing. The five-archetype landmark values live on a
  specific base-point (Drinfeld's even-associator); cross-archetype
  comparisons implicitly require the same base-point. **Wrong
  claim**: ``bar--cobar is natural because KT formality is a
  point''. **Precise error**: contractibility is an iso-class
  statement; the parametrised space is a
  $\mathrm{GRT}_1(\mathbb Q)$-torsor, not a point, and
  structure-constant-level claims must name the base-point.
  **Counter**: every chain-level bar--cobar adjunction coherence
  statement in a proof of Theorem~A / Theorem~D must specify the
  lane: (i) iso-class lane (contractibility suffices); (ii)
  parametrised lane (name the $\mathrm{GRT}_1$ base-point, typically
  Drinfeld's even-associator or Costello--Li propagator). Scope
  qualifier paired with structure-constant witness per AP275.
  **Related**: AP924; Vol~III AP-CY161; AP923 (companion grammar
  distinction); AP285 (universal without universal property); AP910
  (bar--cobar adjunction direction); Willwacher 2014; Tamarkin
  2003; Loday--Vallette 2012.

- **AP925 -- BCOV one-loop curving vs tree-level Yukawa Hodge-degree
  discipline (Vol~I cross-volume concordance with Vol~III AP-CY162).**
  Vol~I's modular characteristic $\kappa_{\mathrm{ch}}(\cA)$ for a
  chiral algebra $\cA = \PhiFA_d(\cC)$ with $\cC$ a CY$_d$ category
  is a Hodge-supertrace invariant
  $\sum_q(-1)^q h^{0,q}(X)$ (AP307, AP-CAT-7 when restricted through
  the $\Phi_d$-bridge on compact CY$_d$). At one-loop the BCOV BV
  anomaly $\alpha_{\mathrm{BCOV}} \in H^{0,1}(X)$ with
  $\alpha_{\mathrm{BCOV}} = (\chi(X)/24)[\Omega_X]^{0,1}$
  (Costello--Li 2016 \emph{arXiv:1605.09930} Prop~5.2) can shift
  $\kappa_{\mathrm{ch}}$ under quantisation; this must be
  distinguished from the tree-level BCOV Yukawa
  $Y_3 \in H^{0,3}(X) = H^3(X, \mathcal O_X)$ (Kapranov 1999
  \emph{Compositio} 115 §4 as $\ell_3^{\min}$ in the $L_\infty$-
  structure on $\Omega^{0,*}(X, T_X)$). The two cocycles are both
  Atiyah-sourced ($\mathrm{At}(T_X) \in H^1(X, \Omega_X \otimes
  \End T_X)$) but land in Hodge-disjoint receptacles:
  $H^{0,1} \cap H^{0,3} = 0$ (Serre-dual under $\chi(\mathcal O_X)
  = \int_X \mathrm{td}(T_X)$ but never simultaneously in the same
  slot). **Vol~I bearing**: any appearance of
  $\kappa_{\mathrm{ch}}$-shift under CY$_d$ quantisation (Theorem~A
  or CY-A Step~4 bar--cobar consistency at $d = 3$) must name
  whether the shift is one-loop ($\alpha_{\mathrm{BCOV}}$, Hodge
  degree $(0,1)$) or tree-level ($Y_3$, Hodge degree $(0,3)$);
  conflating the two collapses Theorem~A's one-loop obstruction
  into a tree-level identity. **Wrong claim**: ``BCOV curving
  equals the Yukawa cubic because both are Atiyah-sourced''.
  **Precise error**: Hodge-degree mismatch; the three
  Atiyah-sourced cocycles $(Y_3, \alpha_{\mathrm{BCOV}},
  \mathrm{td}(T_X))$ live in three distinct Hodge receptacles
  $(H^{0,3}, H^{0,1}, \bigoplus_p H^{p,p})$ and never appear
  simultaneously in the same slot. **Counter**: every
  $\kappa_{\mathrm{ch}}$-shift claim must name the Hodge receptacle
  with explicit $(p, q)$-bidegree; Serre duality relates but does
  not identify. **Related**: AP925; Vol~III AP-CY162; AP307
  (K\"unneth-additive vs K\"unneth-multiplicative
  $\kappa_{\mathrm{ch}}$); AP-CAT-7 (Theorem~H ambient qualifier);
  Costello--Li 2016 Prop~5.2; Kapranov 1999; BCOV 1994; Atiyah 1957.

- **AP926 -- MNOP as $\mathbb E_2$-centre trace identity vs
  change-of-variables (Vol~I cross-volume concordance with Vol~III
  AP-CY163).** The Vol~I five theorems (A, B, C, D, H) all live on
  the same factorisation algebra $\cF_X = \PhiFA_d(\cC)$, and
  trace-identities on the derived centre
  $Z^{\mathrm{der}}_{\mathrm{ch}}(\cA)$ of a dualisable
  $\cA$-module are a fifth-theorem companion to Theorem~C. The
  MNOP identity $Z_{\mathrm{DT}}(X, -q) = Z_{\mathrm{PT}}(X, -q)
  \cdot \mathrm{McMahon}$ together with $-q = e^{iu}$
  (MNOP I--II 2006; Pandharipande--Thomas 2014 \emph{Forum Math.\
  Pi} 2) is routinely presented as ``formal change of variables on
  generating functions''. The Vol~III wave-14 correction surfaces
  the categorical content: $-q = e^{iu}$ encodes the unique
  $\mathbb E_2$-centre automorphism $\sigma \in
  Z(\Zcoh(D^b\mathrm{Coh}(X)))$ intertwining three dualisable
  $\mathbb E_3$-modules $M_{\mathrm{DT}}, M_{\mathrm{PT}},
  M_{\mathrm{GW}} \in \mathrm{Mod}_{\mathbb E_3}^{\mathrm{dual}}$
  (Lurie \emph{HA} 7.3.4.2), with semi-classical qdilog residue
  $(2\sin(ku/2))^{2g-2}$ at each BPS class. **Vol~I bearing**:
  the Theorem~C derived-centre complementarity
  $\kappa(\cA) + \kappa(\cA^!)$ landmark ceiling
  $\in \{0, 8, 13, 250/3, 98/3\}$ is a trace identity on
  $Z^{\mathrm{der}}_{\mathrm{ch}}(\cA)$, not a functor
  comparison; reading MNOP as ``change of variables'' drops the
  $\mathbb E_2$-centre structure that makes the identity a
  categorical statement. **Wrong claim**: ``MNOP is a tautology /
  formal substitution''. **Precise error**: trace identity on the
  centre is mistaken for a symbol-rewriting rule; the substitution
  $-q = e^{iu}$ is the semi-classical residue of a
  dualisable-module intertwining automorphism, not a naming
  convention. **Counter**: every trace-identity statement in the
  Theorem~C family (including Vol~III's MNOP consequences and
  Vol~I's $\kappa + \kappa^!$ landmark ceiling) must name the
  centre object ($Z^{\mathrm{der}}_{\mathrm{ch}}$,
  $Z(\Zcoh)$, etc.), the dualisable-module class, and the
  intertwining automorphism. **Related**: AP926; Vol~III AP-CY163;
  AP-CAT-7; AP287 (Drinfeld double candidate unspecified); AP903
  (Theorem~C derived-centre complementarity); Lurie \emph{HA}
  7.3.4.2; MNOP 2006; Maulik 2019.

- **AP927 -- ``Ten real simple roots'' read in the Igusa datum
  (Vol~I cross-volume concordance with Vol~III AP-CY171).**
  Vol~I's landscape census (`landscape_census.tex`) carries the
  numeric $10$ in two independent Igusa datum rows: (a) the
  weight $\mathrm{wt}(\Phi_{10}) = 10$ of the Igusa cusp form
  (from Igusa 1964 \emph{Amer.\ J.\ Math.} 86); (b) the Borcherds
  weight $\kBKM(\Phi_{10}) = c_{\Phi_{10}}(0)/2 = 10$
  (Borcherds 1995 \emph{Invent.\ Math.} 120 Theorem~10.4). Neither
  integer counts the number of real simple roots of the
  underlying generalised Kac--Moody algebra: the K3-BKM Cartan
  has rank $3$ on $\Lambda^{2,1}_{II}$ (canonical preamble row
  23; Gritsenko--Nikulin 1998 \S 3) with three real simple roots,
  and the Fake-Monster Cartan has rank $26$ on $\mathrm{II}_{25,1}$
  with $26$ real simple roots. Reading ``$10$ real simple roots''
  from either $\mathrm{wt}(\Phi_{10}) = 10$ or
  $\kBKM(\Phi_{10}) = 10$ is a two-step type error
  (automorphic-weight $\to$ root-count, or Borcherds-weight
  $\to$ root-count). **Vol~I bearing**: every
  $\kBKM$-subscripted statement that mentions $\Phi_{10}$ with
  the value $10$ must specify which integer is named
  (form-weight, Borcherds-weight, Cartan rank, or number of
  real simple roots --- four distinct invariants), and must not
  narrate any two as the same. **Wrong claim**: ``the Igusa
  datum contains $10$ real simple roots''. **Precise error**:
  the weight of $\Phi_{10}$, the Borcherds weight
  $c_{\Phi_{10}}(0)/2$, and the count of real simple roots are
  three distinct invariants; exactly two of them equal $10$ by
  coincidence (Borcherds 1995 singular-weight identity
  $c(0)/2 = \mathrm{wt}$ for lattices of signature $(2, 2)$);
  the number of real simple roots is $3$ on $\Lambda^{2,1}_{II}$
  and $26$ on $\mathrm{II}_{25,1}$, never $10$. **Counter**:
  every numerical read of the Igusa datum in
  `landscape_census.tex`, `bar_construction.tex`, or
  `holographic_datum_master.tex` must cite the specific
  invariant and the sublattice; $\Phi_{10}$-context value
  $10 = \mathrm{wt} = \kBKM$ with root count $3$ on
  K3-BKM. **Related**: AP927; Vol~III AP-CY171; AP113 (bare
  $\kappa$ forbidden); AP288 (bare $r(z)$); canonical preamble
  rows 22 (K3-BKM rank 3), 28 ($\Delta_5$ weight 5 denominator),
  60 ($\kBKM(\mathbf H_{\Delta_5})$ cross-volume value);
  Borcherds 1995; Gritsenko--Nikulin 1998; Igusa 1964.

- **AP928 -- $\Delta_5$ half-BPS / $\Phi_{10}$ full-BPS duality at
  automorphic-form level (Vol~I cross-volume concordance with
  Vol~III AP-CY172).** Vol~I's canonical $\kBKM(\Delta_5) = 5$
  (via Gritsenko $\Delta_5 = \mathrm{Grit}(\eta^9 \vartheta_1)
  \in S_5(K(1))$ additive lift, preamble row 28, Gritsenko 1999
  Thm~6.1) and Vol~III's $\kBKM(\Phi_{10}) = 10$ (Borcherds
  weight on the full BKM denominator, preamble row 88) are not
  independent data points: they are related by
  $\Phi_{10} = \Delta_5^2$ as paramodular Siegel forms on
  $\mathrm{Sp}_4(\mathbb Z)$ (Gritsenko--Nikulin 1998), and
  therefore
  $\kBKM(\Phi_{10}) = 2 \cdot \kBKM(\Delta_5) = 10$ by the
  multiplicative-lift law $c_{f^2}(0)/2 = 2 \cdot c_f(0)/2$
  (Borcherds 1995 Thm~13.3). Physically this is the
  heterotic-on-$K3 \times T^2$ / Type-II-on-$K3 \times T^2$
  S-duality at the automorphic-form level: heterotic sees
  $\Delta_5$ as the half-BPS partition function (Gritsenko
  additive lift of
  $\phi_{0,1}^{K3}$); Type-II sees $\Phi_{10}$ as the
  full-BPS partition function (Dijkgraaf--Verlinde--Verlinde
  1997 \emph{Nucl.\ Phys.\ B} 484 for 1/4-BPS $D1$-$D5$-$P$
  states). **Vol~I bearing**: Theorem~C's five-archetype
  landmark ceiling value $\kappa + \kappa^! = 8$ for the
  $\mathsf B$-row (K3 Mukai-enhanced Heisenberg witness) is
  a \emph{three-faces} identity: (i)
  $\kappa + \kappa^! = 8$ via derived-centre complementarity;
  (ii) Mukai pairing $K = 2 c_+ = 8$ with $c_+ = 4$; (iii)
  Lusztig reflection length $\ell_{\mathrm{Lusztig}} = 8$ at
  $\zeta_8$ (preamble row 90). Reading $\kBKM(\Delta_5) = 5$
  and $\kBKM(\Phi_{10}) = 10$ as independent numbers breaks
  the three-faces concordance: the $\mathsf B$-row landmark
  $8$ is not $5$ or $10$, but is coherent with both via
  $8 = 10 - 2$ (the Borcherds-weight floor shift under the
  squaring identity plus the central-charge Mukai correction,
  AP-CY172 wave-14 derivation). **Wrong claim**: ``Vol~I
  $\kBKM = 5$ and Vol~III $\kBKM = 10$ are different values in
  conflict''. **Precise error**: they index different input
  denominators ($\Delta_5$ vs $\Phi_{10} = \Delta_5^2$) under
  the same $c_N(0)/2$ rule; heterotic/Type-II S-duality
  exchanges them, not a value disagreement. **Counter**: every
  $\kBKM$-statement on the $\mathsf B$-row or $\mathrm{Sp}_4$
  Igusa datum must name the input denominator ($\Delta_5$
  half-BPS or $\Phi_{10} = \Delta_5^2$ full-BPS) and cite the
  heterotic/Type-II duality lane; concordance between the
  two values is the squaring identity, not an independence
  claim. **Related**: AP928; Vol~III AP-CY172; AP927
  (``$10$ real simple roots'' companion); canonical preamble
  rows 28, 34, 60, 88, 90; Gritsenko 1999; Gritsenko--Nikulin
  1998; Borcherds 1995 Thm~13.3; DVV 1997; the W26.8
  reconciliation in the canonical-values registry.

- **AP-CAT-14 -- Wave-14 cross-volume session-correction anchor map
  (2026-04-22; cross-volume concordance with Vol~III AP-CY160--172).**
  Six Vol~III wave-14 session-correction anti-patterns cross-cut
  into Vol~I's bar--cobar / modular-characteristic domain and are
  inscribed here as AP923--AP928 with one-to-one concordance to the
  Vol~III originating range AP-CY160--AP-CY165 (hCS / KT formality /
  BCOV / MNOP grammar corrections) plus AP-CY171--AP-CY172 (Igusa
  datum ``$10$ real simple roots'' and $\Delta_5$/$\Phi_{10}$
  duality). Vol~III AP-CY166--AP-CY170 are strictly CY-specific
  (Szendr\H oi NCCR vs local charts, hCS gauge vs dg-autoequivalence,
  and three further CY$_3$-geometric items) with no natural Vol~I
  bearing and are therefore not propagated here. **Mapping**:
  AP-CY160 $\to$ AP923 (hCS vs categorical $\mathbb E_3$-factorisation
  presentation duality);
  AP-CY161 $\to$ AP924 (iso-class vs parametrised KT formality);
  AP-CY162 $\to$ AP925 (BCOV curving vs Yukawa Hodge-degree
  discipline);
  AP-CY163 $\to$ AP926 (MNOP as $\mathbb E_2$-centre trace identity);
  AP-CY171 $\to$ AP927 (``$10$ real simple roots'' Igusa misread);
  AP-CY172 $\to$ AP928 ($\Delta_5$/$\Phi_{10}$ heterotic/Type-II
  S-duality).
  **Rule** (cross-volume discipline): every Vol~III AP-CY with
  bar--cobar / modular-characteristic / derived-centre /
  obstruction-tower / Hochschild-concentration bearing receives a
  Vol~I companion AP within the session-close sweep; the companion
  carries the next available AP number (here AP923--AP928) and
  explicit cross-reference to the Vol~III originating entry.
  Strictly CY-dimension-specific patterns (CY-A/B/C/D
  stratification, Stage-1 hCS quantisation, 6d-hCS) stay in
  Vol~III. **Harmonisation**: AP923 (Theorem~A/B bar--cobar
  input), AP924 (Theorem~A adjunction coherence), AP925
  (Theorem~A $\kappa_{\mathrm{ch}}$-shift discipline), AP926
  (Theorem~C trace identity), AP927 (Theorem~C five-archetype
  landmark numerics), AP928 (Theorem~C $\mathsf B$-row
  three-faces) collectively sharpen the five-theorem core without
  adding new theorems. **Hook regexes** (advisory, deferred per
  AP-CAT-13 discipline): (i) AP923 --- trigger when
  $\PhiFA_3$ or $\mathrm{HH}^\bullet(\cC,\cC)$ appears without
  either ``Costello--Li'' / ``Willwacher'' / ``KT formality'' /
  ``hCS grammar'' within 20 lines; (ii) AP924 --- trigger on
  ``KT formality is contractible'' / ``formality is a point'' /
  ``canonical quantisation'' without ``$\mathrm{GRT}_1$-torsor'' /
  ``Willwacher'' base-point; (iii) AP925 --- trigger on
  ``BCOV'' + ``Yukawa'' co-occurrence without explicit Hodge
  $(p, q)$-bidegree; (iv) AP926 --- trigger on ``MNOP'' +
  ``change of variables'' / ``tautology'' / ``formal
  substitution'' without ``$\mathbb E_2$-centre'' /
  ``dualisable module'' / ``intertwining automorphism''; (v)
  AP927 --- trigger on ``Igusa'' / ``$\Phi_{10}$'' /
  ``$\Delta_5$'' + ``simple roots'' / ``root count'' without
  specifying sublattice rank; (vi) AP928 --- trigger on
  $\kBKM(\Delta_5)$ vs $\kBKM(\Phi_{10})$ co-occurrence without
  the squaring identity $\Phi_{10} = \Delta_5^2$ cited within
  10 lines.
  **Related**: AP-CAT-14; AP923--AP928; AP-CAT-12 (Wave~5 residual
  back-anchor); AP-CAT-13 (hook regex discipline); Vol~III
  AP-CY160--AP-CY172; canonical preamble rows 22, 28, 34, 60, 88,
  90.

### Vol III two-stage factorisation integration wave (2026-04-22, cross-volume audit from Vol III $\Phi_d = \SpCh_{\Sigma_{d-1}, C} \circ \PhiFA_d$; disjoint from AP-CY160--172 / AP923--AP928)

- **AP929 -- Vol III's $\Phi_d$ framed as a single-stage unified
  functor (cross-volume, reader-facing Vol~I references to Vol~III's
  output).** Vol~I chapters and standalones that reference Vol~III's
  CY-to-chiral output routinely framed it as a single-stage functor
  $\Phi_d: \mathrm{CY}\text{-cat}_d \to \mathrm{ChirAlg}$ landing
  directly in Vol~I's $E_1$-chiral ambient. The Vol~III construction
  is *two-stage*: $\Phi_d = \SpCh_{\Sigma_{d-1}, C} \circ \PhiFA_d$
  where (Stage~1) $\PhiFA_d$ is the canonical $E_d$-holomorphic
  factorisation algebra on a CY$_d$ variety, pinned by
  Kontsevich--Tamarkin $E_d$-formality + Costello--Gwilliam--Li
  locality, and (Stage~2) $\SpCh_{\Sigma_{d-1}, C}$ is factorisation
  homology over a $(d-1)$-cycle $\Sigma_{d-1}$ followed by
  restriction to a reference curve $C$. Vol~I's Theorem~A
  bar--cobar adjunction is the $E_1$-chiral *shadow on the
  reference curve $C$ after Stage-2 specialisation*, not the direct
  functor image. **Wrong claim**: ``the CY-to-chiral functor
  $\Phi_d$'' framed as one functor, or ``under $\Phi_d$ the CY$_d$
  category gives an $E_1$-chiral algebra'' without naming the two
  stages. **Precise error**: collapses two distinct
  $(\infty, 1)$-categorical operations (factorisation-algebra
  assignment $\PhiFA_d$ and specialisation
  $\SpCh_{\Sigma_{d-1}, C}$) into one, hiding the family-parameter
  $(\Sigma_{d-1}, C)$ that indexes the Vol~I shadow. A single
  CY$_d$ category admits a *family* of $E_1$-chiral shadows indexed
  by $(\Sigma_{d-1}, C)$; Theorem~A governs *each* shadow on its
  curve, not a unique image. **Healed reading**: every Vol~I
  reference to Vol~III's output names the two stages and the family
  parameter, and labels which stage the consumer theorem applies to.
  **Counter**: Vol~I chapters citing the Vol~III $\Phi_d$ output
  carry a bridge clause ``the Stage-2 $E_1$-chiral shadow
  $\SpCh_{\Sigma_{d-1}, C}(\PhiFA_d(\cC))$ on the reference curve
  $C$''; single-stage phrasing is flagged at hook sweep. Healed
  via Vol~I 2026-04-22 C12 + C12b alignment passes (29 Vol~I
  chapters brought onto the two-stage frame). **Related**: AP247
  (functor terminology requires single target --- the two-stage
  factorisation honours this by making the target of $\PhiFA_d$
  the $E_d$-factorisation-algebra category, distinct from the
  target of $\SpCh$, which is $E_1$-chiral); AP-CY46 ($\Phi_d$-
  ambient concentration enlargement $\{0,1,2,d\}$); AP-CY81 (six
  K3$\times$E constructions as distinct, not six $\Phi_3$-
  applications); V2-AP46 (Theorem~H hinge $d$-dependent); Vol~III
  CLAUDE.md line 52--59; Vol~III
  \texttt{chapters/theory/cy\_to\_chiral.tex}.

- **AP930 -- Six routes to $G(K3 \times E)$ framed as six
  $\Phi_3$-applications (cross-volume, reader-facing Vol~I
  references; companion to AP-CY60 / AP-CY81).** Vol~I material
  naming the ``six routes to $G(K3 \times E)$''
  (CoHA/Schiffmann--Vasserot, Borcherds lift, lattice VOA, Kummer,
  sigma model, BLLPR) occasionally framed them as six applications
  of $\Phi_3$ to the $K3 \times E$ CY$_3$ category. $\Phi_3$ gives
  *one* output per input CY$_3$ category; the six routes are
  *six different constructions*, each producing a candidate algebra
  from distinct geometric or representation-theoretic input,
  converging (conjecturally, CY-C content) on the same
  $G(K3 \times E)$. **Wrong claim**: ``the six $\Phi_3$-routes to
  $G(K3 \times E)$'' or ``applying $\Phi_3$ along six independent
  constructions gives the same output''. **Precise error**:
  functoriality of $\Phi_3$ would produce one output per input; the
  six routes have *different input data* (CoHA on
  $\mathrm{Coh}(K3 \times E)$, Borcherds reflection group on
  $\mathrm{II}_{3, 19} \oplus H$, Frenkel--Kac lattice VOA on the
  Mukai lattice, Kummer orbifold, sigma-model partition function,
  BLLPR DT generating function) and produce outputs whose
  identification is the *content* of CY-C, not a consequence of
  functoriality. Conflating the two erases the adversarial-swarm
  programme (CY-C conjectural) and overadvertises functoriality of
  $\Phi_3$. **Healed reading**: every ``six routes'' phrasing
  carries the bridge clause ``six distinct constructions converging
  (conjecturally, CY-C) on the same $\Phi_3$-output'', with each
  route named and its input data listed. **Counter**: hook regex
  `grep -qE 'six (routes|applications|ways)' && grep -qE 'Phi_3|\\\\Phi_3|K3.*times.*E'`
  flags single-stage six-route phrasing; any hit without
  ``distinct constructions'' clause within 10 lines blocks. Related:
  Vol~III AP-CY60 (original); Vol~III AP-CY81 (generator-level
  falsification across routes); AP929 (two-stage framing); AP247
  (single target); Vol~III CLAUDE.md ``Six routes to
  $G(K3 \times E)$ are six DIFFERENT constructions''.

- **AP931 -- Fake-Monster BKM at $d = 3$ via $K3 \times E$ (correct
  dimension is $d = 5$ via $K3 \times K3 \times E$, cross-volume
  CY-dimension discipline).** Vol~I material occasionally placed the
  Fake-Monster BKM at CY-dimension $d = 3$ by conflating it with the
  K3 $\times$ E Gritsenko--Nikulin $\Delta_5$ paramodular
  construction. The Fake-Monster BKM carries Cartan rank $26$ on
  $\mathrm{II}_{25, 1}$ (canonical preamble row 21) and is the Weyl
  denominator $\Phi_{12}$ image from the signature-$(26, 2)$
  orthogonal Shimura variety; its natural CY-categorical home is the
  CY-5 product $K3 \times K3 \times E$ (lattice
  $\mathrm{II}_{3, 19} \oplus \mathrm{II}_{3, 19} \oplus
  \mathrm{II}_{1, 1}$ has the right signature after primitive
  restriction to $\mathrm{II}_{26, 2}$ for the $\Phi_{12}$ Borcherds
  lift). The $K3 \times E$ CY-3 category carries the *paramodular*
  $\Delta_5$ (Cartan rank $3$ on $\Lambda^{2, 1}_{\mathrm{II}}$,
  canonical preamble row 28), a *different* BKM at a different
  dimension. **Wrong claim**: ``the Fake-Monster BKM is the Vol~III
  $d = 3$ $K3 \times E$ output'' or ``$\Phi_{12}$ lifts from the
  $K3 \times E$ Mukai lattice''. **Precise error**: confuses two
  distinct BKM algebras (Fake-Monster rank-26 on
  $\mathrm{II}_{25, 1}$ vs K3-BKM rank-3 on
  $\Lambda^{2, 1}_{\mathrm{II}}$) and two distinct dimensions
  ($d = 5$ vs $d = 3$); the Cartan rank, automorphic-form weight,
  and singular-weight condition all diverge. **Healed reading**:
  every Fake-Monster citation names the CY-5 $K3 \times K3 \times E$
  (or equivalent rank-26 positive-part lattice input) and the
  signature-$(26, 2)$ orthogonal Shimura variety. Every $K3 \times
  E$ reference uses the paramodular $\Delta_5$, rank-3 K3-BKM, not
  $\Phi_{12}$. **Counter**: before placing any Fake-Monster
  reference at $d = 3$, cross-check the Cartan rank in canonical
  preamble rows 20--21 and 28; any claim of Fake-Monster $\Phi_{12}$
  lift from a rank-3 or rank-4 input is a type error. **Related**:
  canonical preamble rows 20--21 (Monster rank 2, Fake-Monster rank
  26); canonical preamble row 28 (K3-BKM Weyl denominator
  $\Delta_5$); AP-CY71 (Conway $V^{s\natural}$ is *not* a fifth
  $\Psi$-image, parallel misplacement pattern); AP-CY79 (phantom
  Niemeier root systems at misplaced $N$); AP927 (``$10$ real
  simple roots'' Igusa misread, parallel rank/weight conflation);
  Borcherds 1992 \emph{Invent.\ Math.}~109 Thm~3; Gritsenko 1999
  Thm~6.1.

- **AP932 -- K3 Yangian envelope three-way confusion: classical
  $\so(4, 20)$ vs programme $Y_\hbar(\so(4|20))$ vs forbidden
  $\osp(4|20)$ / $\fgl(4|20)$ (cross-volume, super-extension type
  discipline; sharpening of AP246 / AP279 / AP361).** Multiple
  envelope readings of the K3 Yangian coexist in cross-volume prose:
  (i) the *classical* envelope is $\so(4, 20)$, the orthogonal Lie
  algebra on the Mukai lattice of signature $(4, 20)$ of K3 (even
  unimodular $\Lambda^{3, 19} \oplus H$ with Mukai pairing giving
  signature $(4, 20)$); (ii) the *programme-specific super-extension*
  $Y_\hbar(\so(4|20))$ is a non-Kac super-Yangian with
  $\mathbb Z / 2$-graded generators without the symplectic-on-odd
  Kac axiom; (iii) $\osp(4|20)$ is *forbidden* as a Kac-super
  reading because it requires a symplectic form on the odd part,
  which the K3 Mukai pairing does not carry (the signature
  $(4, 20)$ is all orthogonal); (iv) $\fgl(4|20)$ is likewise
  *forbidden* as the naive $\fgl$-super reading because it inverts
  the signature-type assignment (K3's Mukai is orthogonal-indefinite,
  not $\fgl$-Berezinian). **Wrong claim**: ``the K3 Yangian envelope
  is $\osp(4|20)$'' or ``K3 Yangian $= Y(\fgl(4|20))$'' or ``$Y(\osp)$
  and $Y(\so)$ are interchangeable super-extensions''. **Precise
  error**: Kac's super-classification requires symplectic form on
  the odd part for $\osp$ and Berezinian form for $\fgl$-super; the
  K3 Mukai pairing is purely orthogonal, so neither Kac-super class
  applies. The programme-specific $Y_\hbar(\so(4|20))$ is a
  non-Kac construction: $\mathbb Z / 2$-graded orthogonal Yangian
  with the orthogonal reflection-equation Shapovalov datum (ACDFR
  2003 arXiv:math/0304188, super-Sklyanin determinant
  $\mathrm{sBer}\, K(u)$). The prior $Y(\fgl(4|20))$ and
  $Y_\osp(4|20)$ inscriptions were Wave-15 programme-trajectory
  artefacts retracted per AP-CY133 canonical $\so$ discipline.
  **Healed reading**: every K3 Yangian envelope inscription
  specifies (a) which of $\so(4, 20)$ / $Y_\hbar(\so(4|20))$ /
  rejected-$\osp(4|20)$ / rejected-$\fgl(4|20)$ is in force; (b)
  the signature-type assignment (orthogonal, not symplectic, not
  $\fgl$); (c) the Kac vs non-Kac classification. **Counter**: hook
  regex `grep -qE 'osp\\(4\\|20\\)|fgl\\(4\\|20\\)|gl\\(4\\|20\\)'`
  flags any $\osp$ or $\fgl$ reading as AP-candidate; verify
  signature-type in canonical Mukai pairing before accepting. The
  allowed envelope classes are $\so(4, 20)$ (classical) and
  $Y_\hbar(\so(4|20))$ (programme super-extension, non-Kac).
  **Related**: AP246 (signature-type assignment); AP279 (rename-heal
  island-local across consumers); AP239 (naming after physical
  source); AP361 (super-Yangian $\fgl$ / $\osp$ / $\so$ scope
  refinement); line 237 MC1-4 row Wave-15 scope refinement (``K3
  super-Yangian $Y_\osp(4|20)$ lies in the osp-super sub-case''
  retracted per AP-CY133 canonical $\so$ discipline); Vol~III
  CLAUDE.md ``Super-Yangian envelope (osp vs gl vs so)''.

- **AP933 -- $\kappa_{\mathrm{cat}}(K3 \times E) = 2$ (fibre) vs $0$
  (total space, K\"unneth-multiplicative) manifesto-level conflation
  (cross-volume reinforcement of canonical preamble row 63 + AP289 +
  AP-CY68).** A lingering Vol~I reader-facing reference to
  $\kappa_{\mathrm{cat}}(K3 \times E) = 2$ reads the K3 \emph{fibre}
  value $\chi(\mathcal O_{K3}) = 2$ as if it were the total-space
  value, rather than applying the K\"unneth-multiplicative rule
  $\chi(\mathcal O_{X \times Y}) = \chi(\mathcal O_X) \cdot
  \chi(\mathcal O_Y)$. For $K3 \times E$ this gives
  $2 \cdot 0 = 0$, not $2$. The substrate is already covered by
  AP289 (K\"unneth-multiplicative vs additive supertrace), AP-CY68
  (K\"unneth-multiplicative Hodge supertrace), and canonical
  preamble rows 56 and 63. AP933 back-anchors the
  \emph{manifesto-level} reading of this conflation (reader-facing
  chapters and standalones citing the wrong value), distinct from
  the proof-level conflation AP289 catches at the theorem-statement
  tier. **Wrong claim**: ``$\kappa_{\mathrm{cat}}(K3 \times E) =
  \chi(\mathcal O_{K3}) = 2$''. **Precise error**: reads the fibre
  value as the total-space value; confuses fibration data with
  product data; violates the K\"unneth-multiplicative identity on
  $\chi(\mathcal O)$. **Healed reading**: every
  $\kappa_{\mathrm{cat}}(K3 \times E)$ inscription uses the
  total-space value $0 = 2 \cdot 0 = \chi(\mathcal O_{K3}) \cdot
  \chi(\mathcal O_E)$; the fibre value $2 = \chi(\mathcal O_{K3})$
  is retained as a separate invariant with its own subscript
  ($\kappa_{\mathrm{fibre}}(K3) = 2$, canonical row 56).
  **Counter**: hook regex
  `grep -qE 'kappa.{0,3}cat.*K3.*times.*E.*=.*2|chi.*mathcal.*O.*K3.*times.*E.*=.*2'`
  flags reader-facing misstatements; any hit without a K\"unneth
  total/fibre disambiguation within 10 lines blocks. **Related**:
  canonical preamble rows 56, 63; AP289 (K\"unneth multiplicative
  vs additive); AP-CY68 (Hodge supertrace K\"unneth); AP-CY84
  ($\kappa_{\mathrm{ch}} = \chi(\mathcal O_X)$ dimension-restricted);
  AP290 (HZ-7 subscript type-swap); Vol~III CLAUDE.md
  ``$\kappa_{\mathrm{cat}}(K3 \times E) = 0$ (total space,
  K\"unneth-multiplicative), NOT 2 (fibre)''.

- **AP934 -- $\mathrm{CoHA}(\mathbb C^3) = \mathcal W_{1 + \infty}$
  manifesto-level misidentification (cross-volume reinforcement of
  canonical preamble row ``$\mathrm{CoHA}(\mathbb C^3)$
  identification'' + AP-CY7).** A lingering Vol~I reader-facing
  reference to $\mathrm{CoHA}(\mathbb C^3) = \mathcal W_{1 + \infty}$
  claims the \emph{classical limit} as the whole. The correct
  identification is $\mathrm{CoHA}(\mathbb C^3) =
  Y^+(\widehat{\fgl}_1)$, the \emph{positive half} of the affine
  Yangian (Schiffmann--Vasserot 2013 \emph{Duke}~162); the classical
  limit of the affine Yangian recovers $\mathcal W_{1 + \infty}$,
  but the CoHA itself is the quantum/non-classical object. AP934
  back-anchors the manifesto-level reading of this conflation
  (reader-facing chapters citing $\mathcal W_{1 + \infty}$ where
  $Y^+$ is load-bearing), distinct from canonical row coverage.
  **Wrong claim**: ``$\mathrm{CoHA}(\mathbb C^3) =
  \mathcal W_{1 + \infty}$'' or ``the cohomological Hall algebra of
  $\mathbb C^3$ is $\mathcal W_{1 + \infty}$''. **Precise error**:
  takes the classical limit for the full quantum algebra; drops
  the positive-half vs full-Yangian distinction; confuses
  $\mathcal W_{1 + \infty} = \mathcal W_\infty[c] \otimes \mathcal H$
  (Pope--Romans--Shen 1990) with $Y^+(\widehat{\fgl}_1)$
  (Schiffmann--Vasserot; Tsymbaliuk 2017). **Healed reading**:
  every CoHA inscription uses
  $\mathrm{CoHA}(\mathbb C^3) = Y^+(\widehat{\fgl}_1)$;
  $\mathcal W_{1 + \infty}$ appears only in the classical-limit
  shadow with the full-Yangian vs positive-half distinction
  inscribed. **Counter**: hook regex
  `grep -qE 'CoHA.*mathbb C.*3.*=.*W_\\{?1.*infty'` flags
  reader-facing misidentifications; any hit without the
  $Y^+(\widehat{\fgl}_1)$ correction within 10 lines blocks.
  **Related**: canonical preamble row
  ``$\mathrm{CoHA}(\mathbb C^3)$ identification''; canonical row
  ``CoHA vs chiral-algebra type''; AP-CY7 (CoHA $\ne$ $E_1$-chiral
  algebra); Schiffmann--Vasserot 2013; Tsymbaliuk 2017; Vol~III
  CLAUDE.md ``$\mathrm{CoHA}(\mathbb C^3) = Y^+$ (positive half),
  NOT $\mathcal W_{1+\infty}$ (full)''.

- **AP935 -- \texttt{\textbackslash begin\{warning\}} environment
  infiltration of reader-facing manuscript \texttt{.tex} (cross-volume
  voice discipline; companion to AP-CAT-2).** The Chriss--Ginzburg
  manuscript voice prohibits meta-narration, bookkeeping, and
  editorial asides in reader-facing \texttt{.tex} under
  \texttt{chapters/}, \texttt{frame/}, \texttt{examples/},
  \texttt{theory/}, \texttt{connections/}, \texttt{bibliography/}.
  The \texttt{warning} environment carries editorial voice by
  construction (``warning to the reader'', ``note that'', ``do not
  confuse'') and is therefore forbidden in the manuscript body;
  warnings belong in \texttt{notes/} or the first-principles cache
  as antipattern prose. **Wrong claim**: a
  \texttt{\textbackslash begin\{warning\}
  ... \textbackslash end\{warning\}} block appears inside a
  reader-facing chapter as if it were a mathematical environment
  (theorem, remark, definition). **Precise error**: warning is a
  voice register, not a mathematical object; placing it in the
  manuscript substitutes editorial instruction for mathematical
  content, violating the CG standard's ``prose is mathematics, not
  explanation of mathematics'' rule. **Healed reading**: every
  manuscript warning is either (a) rewritten as a \texttt{remark}
  stating the mathematics of the potential conflation
  (Gap/Flaw-as-mathematics, per CLAUDE.md ``When a mathematical
  retraction is genuinely informative...''), or (b) moved to
  \texttt{notes/antipatterns\_catalogue.md} / the cache as an AP
  entry. **Counter**: hook regex
  `grep -qE '\\\\begin\\{warning\\}'` applied to reader-facing
  \texttt{.tex} paths; any hit blocks the commit until converted or
  moved. **Related**: AP-CAT-2 (bookkeeping-vocabulary infiltration
  of manuscript theory chapters); AP235 (Vol~I preface sweep);
  AP288 (session-ledger stale narrative in \texttt{notes/});
  V2-AP29 (AI slop cleanup); CLAUDE.md ``Forbidden in manuscript
  prose''.

- **AP936 -- Cross-volume discipline-token leakage: ``AP-$n$'',
  ``Pattern~$n$'', ``Wave~$n$'', ``cache-canonical'', ``CLAUDE.md'',
  ``manifesto'', ``Beilinson-standard'' vocabulary in reader-facing
  \texttt{.tex} (voice discipline, companion to AP-CAT-2 and
  AP935).** AP-CAT-2 catches the specific tokens ``Wave~$n$'' /
  ``round~$m$'' / ``DNA strand''. The broader cross-volume sweep
  surfaces additional infiltration modes: explicit ``AP113'' /
  ``Pattern~245'' references (bookkeeping-index leak);
  ``cache-canonical'' / ``manifesto'' / ``Beilinson-standard'' /
  ``CLAUDE.md'' prose (process-reference leak); and ``see
  \texttt{notes/}'' / ``as inscribed in the cache''
  (indexing-convention leak). All three modes treat the manuscript
  as a secondary artefact explaining the bookkeeping, rather than a
  primary object carrying its own logical force. **Wrong claim**: a
  reader-facing \texttt{.tex} paragraph contains any of the tokens
  ``AP$n$'' / ``Pattern~$n$'' / ``Wave~$n$'' / ``cache-canonical'' /
  ``manifesto'' / ``Beilinson-standard'' / ``CLAUDE.md'' / ``see
  \texttt{notes/...}''. **Precise error**: these tokens are
  \emph{bookkeeping-ledger} identifiers; they have meaning only
  relative to the project-internal AP/Pattern/Wave registry, not to
  the mathematical content. A reader without access to the registry
  cannot evaluate the referenced claim. **Healed reading**: for
  every manuscript violation, extract the \emph{mathematical}
  statement that the bookkeeping token references and inscribe it
  directly in CG voice. ``per AP289 K\"unneth-multiplicative
  discipline'' becomes ``$\Xi(X \times Y) = \Xi(X) \cdot \Xi(Y)$ by
  K\"unneth''. ``per cache-canonical row 63'' becomes
  ``$\kappa_{\mathrm{cat}}(K3 \times E) = \chi(\mathcal O_{K3})
  \chi(\mathcal O_E) = 2 \cdot 0 = 0$''. **Counter**: hook regex
  `grep -qE 'AP[0-9]{2,4}|Pattern~?[0-9]+|Wave~?[0-9]+|cache.canonical|Beilinson.standard|CLAUDE\\.md|manifesto'`
  applied to reader-facing paths; any hit blocks the commit.
  **Related**: AP-CAT-2 (Wave~$n$ / round~$m$ specific tokens);
  AP235 (preface sweep); AP935 (\texttt{warning} environment);
  AP-CY28 (Vol~III ``AP-N'' leakage); CLAUDE.md ``Bookkeeping
  vocabulary of any kind''.

- **AP-CAT-15 -- Vol III two-stage factorisation wave cross-volume
  propagation (2026-04-22, distinct from AP-CAT-14 wave-14
  session-correction).** The Vol III two-stage factorisation
  $\Phi_d = \SpCh_{\Sigma_{d-1}, C} \circ \PhiFA_d$ integration
  wave (Vol~III CLAUDE.md line 52--59,
  \texttt{cy\_to\_chiral.tex}) surfaces eight cross-volume
  antipatterns of direct relevance to Vol~I reader-facing material
  disjoint from AP-CAT-14's hCS/KT/BCOV/MNOP/Igusa scope:
  (a) single-stage $\Phi_d$ framing $\to$ AP929;
  (b) six routes framed as six $\Phi_3$-applications $\to$ AP930;
  (c) Fake-Monster at $d = 3$ vs $d = 5$ $\to$ AP931;
  (d) K3 Yangian envelope three-way $\to$ AP932;
  (e) $\kappa_{\mathrm{cat}}(K3 \times E) = 2$ fibre-vs-total-space
  $\to$ AP933 (back-anchor to AP289 + AP-CY68 + canonical row 63,
  reinforced at reader-facing tier);
  (f) $\mathrm{CoHA}(\mathbb C^3) = \mathcal W_{1 + \infty}$
  $\to$ AP934 (back-anchor to canonical row + AP-CY7, reinforced
  at reader-facing tier);
  plus two voice-discipline patterns:
  (g) warning-environment infiltration $\to$ AP935;
  (h) discipline-token leakage $\to$ AP936.
  **Rule**: cross-volume integration waves that originate in one
  volume's construction (here Vol~III $\Phi_d$ two-stage) must sweep
  the other two volumes for reader-facing consumers of the
  construction and inscribe APs at each site where the consumer's
  framing drifts from the canonical construction. AP933 and AP934
  are back-anchors to existing cache / canonical-row entries,
  reinforced at the manifesto tier where manifesto-level
  misstatements still appear despite proof-tier coverage.
  **Hook regexes** (advisory, deferred per AP-CAT-13 discipline):
  (i) AP929 --- `grep -qE '\\\\Phi_d|Phi_d' && ! grep -qE 'Sp.?Ch|Sigma_\\{?d-1\\}?|two.stage|Stage.2'` within 40 lines;
  (ii) AP930 --- `grep -qE 'six (routes|applications|ways)' && grep -qE 'Phi_3|\\\\Phi_3' && grep -qE 'K3.*times.*E'` without ``distinct constructions'' within 10 lines;
  (iii) AP931 --- `grep -qE 'Fake.?Monster' && grep -qE 'K3.*times.*E|d\\s*=\\s*3' && ! grep -qE 'K3.*K3.*E|d\\s*=\\s*5|Phi_\\{?12\\}?'`;
  (iv) AP932 --- `grep -qE 'osp\\(4\\|20\\)|fgl\\(4\\|20\\)|gl\\(4\\|20\\)'`;
  (v) AP933 --- `grep -qE 'kappa.{0,3}cat.*K3.*times.*E.*=.*2'` without K\"unneth disambiguation;
  (vi) AP934 --- `grep -qE 'CoHA.*mathbb C.*3.*=.*W_\\{?1.*infty'` without $Y^+$ correction;
  (vii) AP935 --- `grep -qE '\\\\begin\\{warning\\}'` in reader-facing paths;
  (viii) AP936 --- `grep -qE 'AP[0-9]{2,4}|Pattern~?[0-9]+|Wave~?[0-9]+|cache.canonical|Beilinson.standard|CLAUDE\\.md|manifesto'` in reader-facing paths.
  Hook inscription deferred pending false-positive review.
  **Related**: AP-CAT-14 (Vol~III wave-14 corrections, distinct
  scope); AP-CAT-13 (hook regex discipline); AP-CAT-12 (Wave~5
  closure); AP-CAT-2 (bookkeeping infiltration); Vol~III two-stage
  CLAUDE.md; Vol~I 2026-04-22 C12 + C12b alignment passes.

  **Numbering note (AP-CAT-1 discipline)**: `appendices/first\_principles\_cache.md`
  row 446 carries a trailing-column tag ``AP929'' for a Vol~II 6d hCS
  one-loop quartic-Casimir entry (V2-6d-hCS-one-loop-quartic-not-quadratic).
  That tag is a pre-existing divergence: the AP929 number was locally
  assigned in the reader-facing cache \emph{without} a corresponding
  entry in the primary catalogue. The primary catalogue's AP929
  (``Vol III's $\Phi_d$ framed as a single-stage unified functor'',
  inscribed above) is the canonical assignment; the reader-facing
  tag should be renumbered to AP929bis or to a free number in the
  next maintenance sweep per AP-CAT-1 ``adopt a disambiguating
  suffix bis ... or renumber to the next free integer above both
  maxima''. Flagged here for visibility; no immediate edit to the
  reader-facing cache to avoid disrupting load-bearing Vol~II
  content.


## Session antipatterns — manuscript hygiene (2026-04-22)

**Governing principle** (inscribed across Vols I/II/III CLAUDE.md): *The
manuscript is self-complete, self-coherent, self-consistent. The current
version stands for itself and only itself. All LaTeX mathematical
writing is standalone, up-to-date, consistent, coherent. No references
to previous versions, intermediate ansätze, earlier drafts, retracted
values, superseded formulas, or drafting-history commentary. When a
mathematical retraction is genuinely informative, state the failed
argument and its Gap/Flaw as mathematics, not as drafting record.*

The fifty-five patterns below are the sweep of manuscript-hygiene
anti-patterns caught in the Vol III session of 2026-04-22 and
generalised to Vol I prose. They are organised into groups
A (bookkeeping vocabulary), B (meta-narration), C (version-history
commentary), D (bookkeeping-content refs / absolute paths), E (warning
boxes / hedge-language), F (structural meta-naming), and G (edge
cases). Each entry gives the forbidden form, the failure mode, and a
canonical repair generalised to Vol I prose.

### Group A. Bookkeeping vocabulary in prose (CGCLEAN-1 to CGCLEAN-10)

- **CGCLEAN-1 — "Wave $N$" / "wave $N$" / "Wave-$N$" in manuscript
  prose.** Forbidden token-class: `Wave 13`, `Wave 14`, `wave 25`,
  `Wave-20`, `waves 14--19`, `W$N$` bookkeeping subscript. **Failure
  mode**: the Wave $N$ label is an adversarial-swarm session token
  that belongs to `notes/`, commit messages, and agent memory — never
  the reader-facing `.tex` source. Its presence in prose converts a
  mathematical statement into a project-management artefact and forces
  the reader to reverse-engineer a private workflow. **Canonical
  repair**: delete the wave tag; reword the surrounding clause as
  direct mathematics. "In Wave 14 we established…" → "Theorem~C
  establishes…"; "The Wave-20 calculation gives $c = 12$" → "The
  calculation below gives $c = 12$". **Related**: AP296 (Theorem-H
  theory erasure); governing principle above; CLAUDE.md Vol I/III
  "Forbidden in manuscript prose" list.

- **CGCLEAN-2 — "AP-CY$n$" / "AP$n$" / "AP-CAT-$N$" / "V2-AP$n$" in
  prose.** Forbidden: `AP-CY56`, `AP113`, `AP-CAT-3`, `V2-AP127` used
  as load-bearing pointers in reader-facing `.tex` chapters. **Failure
  mode**: AP codes are catalogue indices internal to
  `notes/antipatterns_catalogue.md`; inscribing them in prose binds the
  manuscript to a private numbering scheme that will drift across
  sessions and volumes. The reader asked to "consult AP-CY56" must
  open a bookkeeping file. **Canonical repair**: replace the AP tag by
  the mathematical content it encodes. "per AP-CY56" → "because $E_2$
  lives on the derived centre $Z(\Rep(A))$, not on $A$"; "avoiding
  AP113" → "with each $\kappa$ subscripted by its source
  $\{\mathrm{ch},\mathrm{cat},\mathrm{BKM},\mathrm{fibre}\}$".
  **Related**: `notes/antipatterns_catalogue.md` is the reader-facing
  home for AP indices via `appendices/antipatterns.tex`; prose cites
  theorems and definitions, not catalogue indices.

- **CGCLEAN-3 — "FM$n$" formula-mechanical tag in prose.**
  Forbidden: `FM24`, `FM42`, `FM-27` as prose annotations. **Failure
  mode**: FM codes index the formula-mechanical sub-catalogue of
  `notes/antipatterns_catalogue.md`, identical hygiene issue as
  CGCLEAN-2. **Canonical repair**: replace with the mathematical
  content (e.g., "FM24: B-cycle $i^2 = +1$ sign error" → "where
  $|q| = e^{-2\pi\Im\tau} < 1$ forces the orientation $i^2 = -1$").

- **CGCLEAN-4 — "HZ-$N$" / "HZ-IV" Hypothesis-Zoo tag in prose.**
  Forbidden: `HZ-3`, `HZ-7`, `HZ-11`, `HZ-IV` referencing the
  hypothesis-zoo discipline. **Failure mode**: the HZ register is a
  programme-internal discipline checklist (subscript discipline,
  independent-verification protocol). Inscribing `HZ-7` in prose says
  nothing to the reader; inscribe what the discipline demands.
  **Canonical repair**: "by HZ-7 subscript discipline" → "with
  $\kappa_{\mathrm{BKM}}$ named explicitly to distinguish it from
  $\kappa_{\mathrm{ch}}$ and $\kappa_{\mathrm{cat}}$".

- **CGCLEAN-5 — "DNA strand S$x$" / "strand $S_2$" as organising
  device.** Forbidden: `DNA strand S3`, `strand~S_k`, `the S4
  thread` used to label swarm output threads. **Failure mode**: DNA
  strand is an agent-orchestration construct; in prose it narrates
  swarm workflow rather than mathematics. **Canonical repair**: replace
  with the mathematical thread the strand was tracking. "Strand S_3
  resolved the $\kappa$-complementarity theorem" → "Theorem~C
  resolves the complementarity $\kappa + \kappa^! = 13$".

- **CGCLEAN-6 — "CG-rectify pass $k$" as labour indicator.**
  Forbidden: `CG-rectify pass 4`, `rectification round 12`,
  `the third rectification`. **Failure mode**: pass-count labels the
  labour, not the content; the reader learns that the authors have
  worked a section hard but not what the section contains. **Canonical
  repair**: delete the pass-count clause and state what the pass
  achieved mathematically.

- **CGCLEAN-7 — "cache entry $n$" / "Cached Confusion \#N" /
  "Cache anchor" / "Cache append" in prose.** Forbidden: `cache
  entry 175`, `per Cached Confusion #47`, `cache pattern 236`,
  `cache-anchor AP885`. **Failure mode**: cache indices are internal
  to `notes/first_principles_cache_comprehensive.md` and
  `appendices/first_principles_cache.md`. Prose naming them forces the
  reader into the cache; inscribe the confusion directly. **Canonical
  repair**: "cache entry 175 forces $\delta^{(n)} = \lfloor n/2\rfloor
  + 1$" → "the Yetter-Drinfeld Schauenburg bracket has weight
  $\lfloor n/2\rfloor + 1$".

- **CGCLEAN-8 — "Wave $N$ spec" / "witnessing" / "verdict" / "ledger"
  as meta-labels.** Forbidden: `Wave 18 spec`, `ledger closure`,
  `witnessing the adjunction`, `verdict: proved`. **Failure mode**:
  spec / witnessing / verdict / ledger are agent-pipeline terms that
  reduce theorems to workflow states. **Canonical repair**: a theorem
  is either proved, conjectured, or open; state the status directly
  via `ClaimStatusProved`/`ClaimStatusConjectured` macros or explicit
  prose ("holds under the hypothesis…", "remains open beyond…").

- **CGCLEAN-9 — "programme-canonical" as elevation tag.**
  Forbidden: `the programme-canonical choice`,
  `the programme-canonical normalisation`. **Failure mode**: the word
  "programme" is reader-facing and acceptable in bibliography; the
  compound "programme-canonical" is a swarm-orchestration elevation
  tag that does the work "canonical" already does without the wrapper.
  **Canonical repair**: "the programme-canonical $\kappa_{\mathrm{ch}}$"
  → "the canonical $\kappa_{\mathrm{ch}}$".
  **Legitimate**: "the Borcherds programme", "the geometric Langlands
  programme" (as a bibliographic reference to a named research
  programme) remains fine. Forbidden is only the elevation-tag
  compound.

- **CGCLEAN-10 — "type-error registry entry T$n$" / "anchor T$n$"
  as cross-reference.** Forbidden: `anchor T20`, `type-error
  registry T7`. **Failure mode**: identical to CGCLEAN-2/CGCLEAN-7;
  T-codes are internal indices into a registry not inscribed in the
  manuscript. **Canonical repair**: cite the mathematical content
  directly.

### Group B. Meta-narration (CGCLEAN-11 to CGCLEAN-25)

- **CGCLEAN-11 — "narrative counterpart" / "narrative arc" /
  "…narrative" as organising device.** Forbidden: `the narrative
  counterpart of Theorem~A`, `the narrative arc of Chapter 5`,
  `the functorial narrative`. **Failure mode**: narrative is a
  property of stories; mathematics has statements, constructions, and
  proofs. "Narrative counterpart" is either a statement (in which case
  state it) or decoration (in which case delete it). **Canonical
  repair**: "the narrative counterpart of Theorem~A is the
  bar-cobar adjunction" → "Theorem~A extends to the bar-cobar
  adjunction as follows".

- **CGCLEAN-12 — "story" / "saga" / "odyssey" / "journey" as
  noun-forms.** Forbidden: `the story of $r(z)$`, `the saga of
  $\{b_k, B^{(2)}\}$`, `the $\kappa$-odyssey`, `the K3 journey`.
  **Failure mode**: rhetorical inflation that narrates the authors'
  experience rather than the mathematics. **Canonical repair**: "The
  story of $r(z)$ begins with the Yang-Baxter equation" → "The
  $r$-matrix $r(z)$ is characterised by the Yang-Baxter equation".
  **Watch**: "story" as a verb ("this construction stores the data of
  $\ldots$") is OK; only the narrative-noun form is forbidden.

- **CGCLEAN-13 — "Platonic ideal" / "Platonic form" / "platonic
  chapter" / "platonic architecture" / "Platonic ensemble" /
  "platonic synthesis" / "Platonic-form construction".** Forbidden:
  `the Platonic ideal of Theorem~D`, `the platonic chapter synthesis`,
  `the Platonic $\kappa$-spectrum`. **Failure mode**: Platonic
  language is elevated bookkeeping vocabulary (from the
  `chriss-ginzburg-rectify` skill's internal target of a "Platonic
  ideal chapter"); in prose it narrates the authorial goal rather than
  stating mathematics. **Canonical repair**: delete the Platonic tag;
  state the construction directly. "the Platonic-form $\Phi$ functor" →
  "the $\Phi$ functor". **Watch**: "Platonic" in a music or
  philosophy context (e.g., a reference to Plato's *Timaeus*) is
  legitimate if the mathematical content is about Plato's work;
  forbidden is the bookkeeping elevation tag.

- **CGCLEAN-14 — "Platonic Theorem~A" / "Platonic Theorem~C" as
  name-tag.** Forbidden: `the Platonic Theorem~A`,
  `the Platonic form of Theorem~C`. **Failure mode**: the compound
  "Platonic Theorem~X" is a rectification-session bookkeeping label
  that elevates the target; the theorem is either Theorem~A or it
  isn't. **Canonical repair**: "Platonic Theorem~A" → "Theorem~A".

- **CGCLEAN-15 — "This chapter's function is to…" meta-paragraph.**
  Forbidden: `This chapter's function is to establish$\ldots$`,
  `The role of this chapter is$\ldots$`, `This chapter serves as
  the foundation$\ldots$`. **Failure mode**: narrating the chapter's
  purpose rather than executing it. The Chriss--Ginzburg north star is
  *show don't tell*; the chapter demonstrates its function by stating
  and proving theorems. **Canonical repair**: delete the meta-paragraph;
  start with a definition, a question, or a theorem.

- **CGCLEAN-16 — "we now turn to" / "having established" / "let us
  now" / "this brings us to" as signpost.** Forbidden: all forms of
  authorial-intent signposting at section boundaries. **Failure mode**:
  signposts replace the three-sentence CG section boundary (what was
  just established; the question the next section resolves; the
  construction that resolves it) with empty narrative scaffolding.
  **Canonical repair**: state what the next section constructs without
  an "us" or a "turn". "We now turn to the cobar complex" →
  "The cobar complex $\Omega(C)$ is $\ldots$". **Related**: already
  partially flagged by hook's existing "Signpost phrases" block;
  CGCLEAN-16 extends and hardens the existing list.

- **CGCLEAN-17 — "in the present work" / "the author" /
  "our programme" / "we have argued" / "it is worth noting".**
  Forbidden: all authorial self-references in manuscript prose.
  **Failure mode**: elite mathematical writing (Gelfand, Manin,
  Drinfeld, Etingof) names constructions and theorems, not authors.
  The "present work" is the text in front of the reader; the authorial
  voice is the mathematics. **Canonical repair**: "In the present
  work, we establish$\ldots$" → "Theorem~A establishes$\ldots$"; "our
  programme extends$\ldots$" → "the construction extends$\ldots$"; "we
  have argued that $X = Y$" → "$X = Y$ because $\ldots$".

- **CGCLEAN-18 — Meta-paragraphs of the form "This chapter closes
  the $\ldots$ programme".** Forbidden: closing-the-programme
  paragraphs, opening-the-programme paragraphs, any pass of
  programme-meta-narration. **Failure mode**: same as CGCLEAN-15 but
  at programme rather than chapter scope. **Canonical repair**:
  delete; the last theorem of the chapter speaks for itself.

- **CGCLEAN-19 — "the opening paragraphs of this preface".**
  Forbidden: self-reference to the document's own structure. **Failure
  mode**: the preface opening IS the opening paragraphs; naming them
  in the preface is circular. **Canonical repair**: delete the
  self-reference.

- **CGCLEAN-20 — "Earlier in the volume" / "at earlier drafts" /
  "in an earlier version".** Forbidden: cross-temporal reference
  inside the current document. **Failure mode**: violates the
  self-complete principle; the current manuscript is the version that
  stands, "earlier" versions are irrelevant to the reader. **Canonical
  repair**: remove the temporal qualifier; if a theorem is established
  in Chapter~$k$ of the current volume, cite it as
  `\ref{thm:$\ldots$}` in Chapter~$k$, not as "earlier".

- **CGCLEAN-21 — "notably" / "crucially" / "remarkably" /
  "interestingly" / "importantly".** Forbidden adverbs: the LLM /
  signposting tell that marks "here is something important" rather than
  making the importance visible in the mathematics. **Canonical
  repair**: delete the adverb and restate the claim with load-bearing
  mathematical content. **Related**: already flagged by the hook's
  existing "AI slop words" block under SECTION 1; CGCLEAN-21 records
  it explicitly as a manuscript-hygiene antipattern.

- **CGCLEAN-22 — "furthermore" / "moreover" / "additionally" as
  sentence openers.** Same category as CGCLEAN-21; already flagged
  by the hook's AP134 block. **Canonical repair**: restate the follow-up
  clause directly, typically as a new sentence that names the new
  mathematical object.

- **CGCLEAN-23 — "this preface's role is to".** Same failure mode as
  CGCLEAN-15, at preface scope. **Canonical repair**: delete; the
  preface executes its role by being the preface.

- **CGCLEAN-24 — "show what / tell what"-mix: narrating what the
  section will show before showing it.** Forbidden template: a
  section-opening paragraph that describes what the section's theorems
  say without stating them. **Failure mode**: the reader gets the
  narration twice — once as narration, once as the theorem. Violates
  the CG economy rule. **Canonical repair**: delete the narration,
  begin with the first definition or theorem.

- **CGCLEAN-25 — authorial-self-intent meta ("in this paper", "in
  this volume", "throughout this work" used as organising devices).**
  Same as CGCLEAN-17 / CGCLEAN-19. Delete.

### Group C. Version-history commentary (CGCLEAN-26 to CGCLEAN-40)

- **CGCLEAN-26 — "retracted" / "retraction" / "retracted ansatz" /
  "now retracted" / "the retracted value".** Forbidden in manuscript
  prose: any mention that a prior claim was retracted. **Failure
  mode**: violates the self-complete principle. The reader of the
  current manuscript needs the current claim, not the history. If a
  prior wrong claim contained mathematical content worth preserving,
  state the failed argument and its Gap/Flaw as *mathematics*:
  "The naive identity $\kappa_{\mathrm{BKM}} = \kappa_{\mathrm{ch}} +
  \chi(\cO_{\mathrm{fibre}})$ holds at $N = 1$ by the coincidence
  $5 = 3 + 2$; at $N = 2$ the left side is $4$, the right side is
  $1$, and the universal formula $\kappa_{\mathrm{BKM}}(\Phi_N) =
  c_N(0)/2$ takes over." — no "retraction" vocabulary needed.
  **Canonical repair**: delete "retracted" / "retraction" mentions;
  inscribe the mathematical Gap/Flaw if informative, else delete
  entirely.

- **CGCLEAN-27 — "superseded" / "supersedes" in prose.** Same as
  CGCLEAN-26. **Canonical repair**: state the current identity; do
  not chronicle its predecessor.

- **CGCLEAN-28 — "earlier draft" / "previous version" /
  "intermediate ansatz" / "prior derivation".** Same as CGCLEAN-26.
  **Canonical repair**: the current version is the only version the
  reader sees.

- **CGCLEAN-29 — "previously conjectural" / "previously open" /
  "previously unresolved" / "previously obstructing".**
  Forbidden temporal qualifier. **Failure mode**: the reader does not
  need to know the claim was formerly a conjecture; they need to know
  its current status. **Canonical repair**: "previously conjectural,
  now proved" → "Theorem~C establishes$\ldots$".

- **CGCLEAN-30 — "now resolved" (drop the "now").** Forbidden: all
  "now" temporal qualifiers. **Canonical repair**: "is now resolved" →
  "is resolved"; "now proved" → "proved".

- **CGCLEAN-31 — "double-retraction" / "Three successive
  evaluations" / "History of the claim" as manuscript prose.**
  Forbidden: retraction-of-retraction chronicles. **Failure mode**:
  doubly violates CGCLEAN-26; historical drafting notes are strictly
  `notes/`-material. **Canonical repair**: delete; the current
  formulation is the only one that appears.

- **CGCLEAN-32 — "drafting record" / "drafting trajectory".**
  Forbidden. **Canonical repair**: delete.

- **CGCLEAN-33 — `\ClaimStatusRetracted` macro in manuscript prose.**
  Forbidden: any `\ClaimStatusRetracted` tag on a theorem / proposition
  / lemma / remark block. **Failure mode**: the retracted-status tag
  belongs to the `notes/` catalogue history, not the theorem-status
  surface seen by the reader. A retracted claim is a claim the
  manuscript does not make; make a different claim. **Canonical
  repair**: delete the environment entirely, or replace with a correct
  theorem using `\ClaimStatusProved` / `\ClaimStatusConjectured`. If
  the retracted argument is informative, rewrite as a scoped
  statement: "The identity fails at $N = 2$ by the computation
  $\kappa_{\mathrm{BKM}}(N=2) = 4 \ne 1 = \kappa_{\mathrm{ch}} +
  \chi(\cO_{\mathrm{fibre}})$", inside a proper theorem block.

- **CGCLEAN-34 — dated remarks ("2026-04-17", "2026-04-19",
  "(Etingof 2026-04-19)").** Forbidden: inline dates in manuscript
  prose. **Failure mode**: the date is drafting metadata; the theorem
  is either true or false, not "true as of April 19th". If a
  communication with Etingof fixed the argument, cite his paper or
  personal communication without the date. **Canonical repair**:
  "(Etingof 2026-04-19)" → "\cite{Etingof_personal_communication}" or
  delete.

- **CGCLEAN-35 — `\index{retraction!$\ldots$}` index entries.**
  Forbidden: any `\index{retraction$\ldots$}` in reader-facing source.
  **Failure mode**: index entries tagged "retraction" tell the reader
  where in the book to find the retraction history; the book does not
  have retraction history (per CGCLEAN-26). **Canonical repair**:
  delete the index entry.

- **CGCLEAN-36 — "(now retracted)" / "(retracted)" /
  "(now corrected)" parenthetical.** Forbidden. **Canonical repair**:
  delete the parenthetical; restate with the current formulation.

- **CGCLEAN-37 — "as previously believed" / "as previously
  stated".** Forbidden. **Canonical repair**: state the current
  belief without referencing the previous one.

- **CGCLEAN-38 — "originally claimed" / "the original argument".**
  Forbidden. **Canonical repair**: the current argument is the
  argument.

- **CGCLEAN-39 — "pre-rectification" / "post-rectification" as
  temporal dividers.** Forbidden. **Canonical repair**: delete.

- **CGCLEAN-40 — "drafting version $k$" / "session $N$ inscription"
  as versioning tags.** Forbidden. **Canonical repair**: remove all
  version tags; the manuscript is the manuscript.

### Group D. Bookkeeping-content refs / absolute paths (CGCLEAN-41 to CGCLEAN-47)

- **CGCLEAN-41 — `\texttt{notes/$\ldots$}` reader-facing reference.**
  Forbidden: `\texttt{notes/antipatterns_catalogue.md}`,
  `\texttt{notes/first\_principles\_cache\_comprehensive.md}` cited
  from manuscript prose. **Failure mode**: `notes/` is the internal
  bookkeeping directory; a reader of the published manuscript cannot
  access it. **Canonical repair**: cite the mathematical content
  directly, or if a catalogue entry is reader-facing-relevant, move it
  to `appendices/` and cite as `\ref{app:$\ldots$}`.

- **CGCLEAN-42 — `/Users/raeez/$\ldots$` absolute path in LaTeX.**
  Forbidden: any filesystem absolute path in `.tex` source. **Failure
  mode**: publishes author-private paths; nonportable; breaks build
  on any other machine. **Canonical repair**: use `\input{$\ldots$}`
  with a path relative to the `main.tex` root, or cite the chapter by
  `\ref{ch:$\ldots$}`.

- **CGCLEAN-43 — `% TODO: librarian verification` comments.**
  Forbidden: `% TODO`, `% FIXME`, `% XXX` annotations surviving into
  published `.tex`. **Failure mode**: todo markers are drafting
  metadata; their presence means the content is incomplete. **Canonical
  repair**: either complete the task (verify the citation, fix the
  gap) or remove the content the todo annotates.

- **CGCLEAN-44 — `% ALIAS` / `% LEGACY ALIAS` / "both keys used in
  prose" / "consolidate in future revision" comments.** Forbidden.
  **Failure mode**: aliases and legacy-key consolidation notes are
  bibliography-hygiene metadata. **Canonical repair**: pick one key,
  delete the alias comment, enforce the single key across prose.

- **CGCLEAN-45 — `% Source: NEW CHAPTER (see notes/$\ldots$)`
  provenance comments.** Forbidden. **Failure mode**: provenance
  comments encode drafting history in source; the content of the file
  stands on its own. **Canonical repair**: delete the provenance line.

- **CGCLEAN-46 — compute-engine filenames `*_waveN_*.py` in prose.**
  Forbidden: `\texttt{compute/lib/k3\_yangian\_wave14\_arthur\_hecke\_delta10.py}`
  cited from manuscript prose. **Failure mode**: the wave-$N$ stem
  embeds bookkeeping into filename, and the filename propagates into
  the manuscript via "see compute engine $X$". **Canonical repair**:
  rename the compute file to drop the wave tag (e.g.,
  `k3\_yangian\_arthur\_hecke\_delta10.py`), then cite; alternatively
  cite by the mathematical content the engine computes without naming
  the file. **Watch**: compute engines themselves (`.py`) may carry
  wave tags in their *history* (git log), but the filename, module
  docstring, and manuscript prose citation must be wave-free.

- **CGCLEAN-47 — function names `waveN_foo` in compute code or
  prose.** Forbidden: `def wave14_compute_tau($\ldots$)`. **Failure
  mode**: same as CGCLEAN-46 at function granularity. **Canonical
  repair**: rename the function to describe what it computes.

### Group E. Warning boxes / hedge-language (CGCLEAN-48 to CGCLEAN-51)

- **CGCLEAN-48 — `\begin{warning}$\ldots$\end{warning}` environment.**
  Forbidden in manuscript prose. **Failure mode**: warning boxes
  narrate the author's concern about reader confusion; the CG standard
  handles this by sharpening the statement so confusion cannot arise.
  **Canonical repair**: delete the warning box; rewrite the adjacent
  statement so the content the warning was defending becomes the
  default reading. "WARNING: do not confuse the Drinfeld centre with
  the derived centre" → "The Drinfeld centre
  $Z(\mathcal{C}) \ne Z^{\mathrm{der}}(A)$; see
  \ref{def:drinfeld-centre} versus \ref{def:derived-centre}".

- **CGCLEAN-49 — "do not confuse" / "don't be fooled" / "beware" /
  "be careful" / "we must be careful".** Forbidden. **Failure mode**:
  hedge vocabulary warns the reader about a trap the prose is laying;
  lay no trap. **Canonical repair**: delete the warning; state the
  distinction as mathematics.

- **CGCLEAN-50 — gratuitous "scope-restricted to" / "scope-qualified
  to" inflation.** Forbidden when used decoratively. **Failure
  mode**: legitimate scope qualifiers name a hypothesis; "Theorem~C
  holds on the Koszul locus" is fine. "Theorem~C is scope-restricted
  to the Koszul locus where its proof is scope-qualified by the
  Koszul hypothesis" is decoration. **Canonical repair**: state the
  hypothesis once, cleanly. "Theorem~C holds under the Koszul
  hypothesis (\ref{hyp:koszul})".

- **CGCLEAN-51 — "verdict" as meta-label.** Forbidden: `Verdict:
  proved`, `Verdict: open`, `Verdict: conditional`. **Failure mode**:
  verdict is adjudication-pipeline vocabulary. **Canonical repair**:
  "Verdict: proved" → "Conclusion:" or simply start with
  "\begin{proof}"; the structure speaks for itself.

### Group F. Structural meta-naming (CGCLEAN-52 to CGCLEAN-54)

- **CGCLEAN-52 — chapter filenames `*_platonic.tex` / chapter labels
  `ch:*-platonic` / section labels `sec:*-platonic`.** Forbidden:
  `chapters/theory/theorem_C_refinements_platonic.tex`,
  `\label{ch:cy-to-chiral-platonic}`. **Failure mode**: file and
  label names encode the rectification-target bookkeeping into the
  build system. **Canonical repair**: rename file to drop the
  `_platonic` suffix; update all `\ref{ch:$\ldots$-platonic}` to
  `\ref{ch:$\ldots$}`; commit as a structural rename across the
  three-volume suite.

- **CGCLEAN-53 — theorem labels `thm:*-waveN-*`.** Forbidden:
  `\label{thm:kappa-additivity-wave14}`,
  `\label{thm:borcherds-weight-wave22-closure}`. **Failure mode**:
  wave-tagged labels migrate the bookkeeping into permanent
  build-system identifiers; labels persist across sessions and become
  load-bearing citations for later chapters. **Canonical repair**:
  rename to `\label{thm:$\ldots$}` without the wave stem; update every
  `\ref` and `\eqref` in three volumes.

- **CGCLEAN-54 — `\index{compute module!$\ldots$}` /
  `\index{cache!$\ldots$}` / `\index{retraction!$\ldots$}` in
  manuscript.** Forbidden. **Failure mode**: index entries tagged
  `compute module!`, `cache!`, `retraction!` expose internal
  bookkeeping in the published index. **Canonical repair**: delete the
  bookkeeping-flavoured index entry; index by mathematical content
  only.

### Group G. Edge cases (CGCLEAN-55)

- **CGCLEAN-55 — programme-internal labour/calibration vocabulary as
  section content.** Forbidden: `Five attack-heal calibrations`,
  `Reconstitution if the cancellation fails`,
  `Inversion of the programme perspective`, `History of the claim`,
  `Gold-standard HZ-IV disjoint verification`, `Three successive
  evaluations appear in the drafting record`. **Failure mode**:
  each of these phrases labels a programme-internal activity
  (attack-heal cycles, cancellation-recovery contingencies,
  HZ-disjoint-path protocols, evaluation chronicles) and inscribes
  that activity into the manuscript. The activity is real and
  valuable; its inscription in reader-facing prose is the defect.
  **Canonical repair**: extract the mathematical content the activity
  established and inscribe that. "Five attack-heal calibrations
  confirm $c_3 = -8$" → "Six independent derivations
  (Borcherds product, Bruinier reduced class, Gritsenko $\Delta_5$
  coefficient, class-number computation, Heegner-divisor weight,
  direct Fourier expansion) give $c_3 = -8$"; "Inversion of the
  programme perspective" → delete; "History of the claim" → delete;
  "Gold-standard HZ-IV disjoint verification" → "Four independent
  derivations" with the paths enumerated directly; "Three successive
  evaluations appear in the drafting record" → delete and inscribe the
  final value.

**Operating rule for Groups A–G.** Every manuscript-prose inscription
(reader-facing `.tex` in `chapters/`, `frame/`, `examples/`, `theory/`,
`connections/`, `bibliography/`) is swept against CGCLEAN-1 through
CGCLEAN-55 at the PostToolUse hook. A pattern fires when its signature
matches and none of the documented legitimate contexts ("Borcherds
programme", "homological retraction", "Platonic" in music / philosophy
citation, etc.) applies. False positives are healed by expanding the
hook's exclusion list; true positives are healed by the canonical
repair above.

**Related**: governing principle at the top of this section;
CLAUDE.md Vol I/III "Forbidden in manuscript prose" list;
`chriss-ginzburg-rectify` skill; `appendices/first_principles_cache.md`
Session antipatterns table; `.claude/hooks/beilinson-gate.sh` CGCLEAN
block.

## Vol~III Wave~12 cross-programme imports (2026-04-22)

The ten entries below (AP929--AP938) are Vol~III wave-12 synthesis
patterns cross-imported into Vol~I's catalogue under the AP-CAT-15
discipline. Each cites the originating `notes/wave12_*.tex` source in
Vol~III, so a Vol~I reader can pull the crystallising argument
directly. These imports sharpen the bar--cobar / shadow-tower /
MC-element / modular-characteristic / derived-centre territory;
they are content-orthogonal to the wave-14 session-correction block
AP923--AP928 (AP-CAT-14, numerical-value discipline) and to the
CGCLEAN manuscript-hygiene block (prose-bookkeeping discipline).

- **AP929 -- Pattern~273 formalisation: object-level chain-level
  statement vs $(\infty,1)$-categorical functor statement are two
  theorems about two categorical structures (Vol~III Wave~12
  cross-import, 2026-04-22, source
  `notes/wave12_a1_phi_functor_foundations.tex`).**
  The chiral bar--cobar adjunction
  $\Omega^{\mathrm{ch}} \dashv B^{\mathrm{ch}}$ admits two readings,
  load-bearing at different scopes, that must not be collapsed.
  (i) The object-level chain-level reading
  $\Omega^{\mathrm{ch}}(B^{\mathrm{ch}}(\cA))
  \xrightarrow{\simeq} \cA$
  lives on explicit filtered cdga / chiral factorisation
  $\cD$-module data, with a splitting witness supplied on the
  Koszul locus (Loday--Vallette 2012 Ch.~11; Positselski 2011
  Mem.~AMS~996). (ii) The $(\infty,1)$-categorical reading
  $\Omega^{\mathrm{ch}} \colon
  \mathrm{coAlg}^{\mathrm{ch}}_{\infty}
  \rightleftarrows
  \mathrm{Alg}^{\mathrm{ch}}_{\infty} \colon B^{\mathrm{ch}}$
  at the level of $(\infty,1)$-categories of factorisation
  algebras (Francis--Gaitsgory 2012 \emph{Selecta Math.}~18
  Thm~3.1) depends on morphism preservation. The two are
  \emph{different theorems} about \emph{different categorical
  structures}: the first needs the splitting witness; the second
  needs the functor-of-$(\infty,1)$-categories witness. Each is
  real and load-bearing; neither is the shadow of the other.
  **Wrong claim**: ``Theorem~A is an adjunction'' (bare, lane
  unspecified). **Precise error**: collapses two theorems into
  one, silently inheriting the witness the reader expects.
  **Counter**: every invocation of Theorem~A names the lane;
  object-level carries the Koszul splitting witness;
  $(\infty,1)$-categorical carries the morphism-preservation
  witness; ``both lanes'' is acceptable when both witnesses are
  supplied. **Related**: AP929; AP905 (adjunction direction);
  AP910 (adjunction direction follow-up); AP923 (hCS vs
  categorical $\mathbb E_3$ grammar); AP-CAT-7 (Theorem~H ambient
  qualifier); Pattern~236 (ambient qualifier); Loday--Vallette
  2012; Francis--Gaitsgory 2012.

- **AP930 -- Bare $\kappa$ forbidden: four $\kappa_\bullet$-invariants
  in the cross-programme spectrum (Vol~III Wave~12 cross-import,
  2026-04-22, source
  `notes/wave12_a2_kappa_invariants_universal_borcherds.tex`).**
  The wave-12 synthesis crystallises four distinct
  $\kappa_\bullet$-invariants across the CY-to-chiral programme:
  $\kappa_{\mathrm{ch}}$ (chiral-side shadow-tower coefficient,
  Vol~I-native, computed via $\Phi$);
  $\kappa_{\mathrm{cat}} = \chi(\mathcal O_X)$ (categorical Euler,
  K\"unneth-\emph{multiplicative} on products, so
  $\kappa_{\mathrm{cat}}(K3 \times E) = \chi(\mathcal O_{K3}) \cdot
  \chi(\mathcal O_E) = 2 \cdot 0 = 0$, \emph{not} $2$);
  $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$ (Borcherds 1995
  \emph{Invent.\ Math.}~120 Thm~10.4; Gritsenko 1999);
  $\kappa_{\mathrm{fiber}}$ (fibre-lattice correction). Vol~I's
  shadow tower natively reports $\kappa_{\mathrm{ch}}$; on
  K3-fibred CY$_3$ and at $d = 3$ the three $\kappa_\bullet$ are
  pairwise distinct. The identity
  $\kappa_{\mathrm{BKM}} = \kappa_{\mathrm{ch}} +
  \chi(\mathcal O_{\mathrm{fiber}})$ is FALSE at every
  $N \in \{1,2,3,4,6\}$: at $N = 1$, LHS $= 5$ (Gritsenko $\Delta_5$
  weight), RHS $= 0 + 0 = 0$.
  **Wrong claim**: bare $\kappa$ in prose; or
  $\kappa_{\mathrm{ch}} \equiv \kappa_{\mathrm{BKM}}$ at $N = 1$.
  **Precise error**: erases which of the four invariants is
  computed; no single universal identity relates them.
  **Counter**: every $\kappa$ carries a subscript at every
  occurrence; Vol~I shadow-tower rows reporting
  $\kappa_{\mathrm{ch}}$ carry the subscript even when no other
  $\kappa_\bullet$ is in scope. **Related**: AP930; AP113; AP925
  (BCOV curving vs Yukawa); HZ-7 subscript discipline; Vol~III
  \texttt{chapters/examples/cy\_d\_kappa\_stratification.tex}
  canonical table.

- **AP931 -- CoHA $= Y^+$ not $\cW_{1+\infty}$: associativity-class
  disambiguation (Vol~III Wave~12 cross-import, 2026-04-22, source
  `notes/wave12_a11_coha_y_plus_vs_w_infty.tex`).**
  Schiffmann--Vasserot 2013 \emph{Publ.\ Math.\ IHES}~118 prove
  that the cohomological Hall algebra of $\mathbb C^3$ equals the
  positive half $Y^+ = Y(\widehat{\mathfrak{gl}}_1)^+$ of the
  affine Yangian, \emph{not} the full $\cW_{1+\infty}$ vertex
  operator algebra. The precise bridge is
  $Y^+ \hookrightarrow Y \xrightarrow{\mathrm{ev}_\lambda}
  \mathrm{End}(\cW_{1+\infty}[\lambda])$: CoHA sits inside the full
  Yangian, which acts on the VOA at spectral parameter $\lambda$.
  The associativity classes differ: CoHA is an $E_1$-associative
  cohomological-Hecke algebra with Hall product; $\cW_{1+\infty}$
  is an $E_2$-chiral VOA with OPE. Vol~I's shadow-tower
  identifications at class $\mathsf M$ cross into this territory
  and need the same disambiguation: ``CoHA of $\mathbb C^3$
  identifies with $\cW_{1+\infty}$'' is a type error.
  **Wrong claim**: ``CoHA$(\mathbb C^3) = \cW_{1+\infty}$''.
  **Precise error**: conflates two algebras of different
  associativity class and different composition.
  **Counter**: every CoHA vs $\cW_{1+\infty}$ identification names
  the class ($E_1$ Hecke vs $E_2$ chiral), the inclusion
  $Y^+ \hookrightarrow Y$, and the spectral parameter
  $\lambda$ of the evaluation. **Related**: AP931; AP270;
  Schiffmann--Vasserot 2013; Miki 2007 \emph{J.\ Math.\ Phys.}~48;
  Vol~III Theorem CoHA-to-Yangian at class $\mathsf M$.

- **AP932 -- $E_n$-level on the wrong object: at $d \ge 3$, $A$ is
  $E_1$, $E_2$ lives on $Z(\mathrm{Rep}(A))$ only (Vol~III Wave~12
  cross-import, 2026-04-22, source
  `notes/wave12_a8_en_hierarchy_kapranov.tex`).**
  The $E_n$-level of a chiral algebra in the CY-to-chiral
  programme is $d$-dependent. At $d = 2$ (K3), the chiral algebra
  $A$ is $E_2$-chiral; at $d \ge 3$ (K3 $\times$ E, quintic, local
  $\bP^2$), $A$ is only $E_1$-chiral and the $E_2$-structure
  migrates to the Drinfeld centre $Z(\mathrm{Rep}(A))$ of its
  representation category. Vol~I's coderived-category machinery
  and iterated Sugawara / Swiss-cheese operad content at
  $d \ge 3$ ambients must respect this: an $E_2$-claim on the
  chiral algebra itself is a type error at $d \ge 3$, even when
  the cohomological shadow looks $E_2$-like.
  **Wrong claim**: ``$A$ is $E_2$-chiral'' stated uniformly across
  $d$; ``the Swiss-cheese operad acts on $A$'' at $d \ge 3$.
  **Precise error**: collapses the PTVV shift-law dimensional
  stratification.
  **Counter**: every $E_n$-claim names the dimension $d$ and the
  object carrying the structure; at $d \ge 3$, $E_2$-statements
  carry the qualifier ``on $Z(\mathrm{Rep}(A))$''.
  **Related**: AP932; AP910; AP923; AP926 (MNOP as
  $\mathbb E_2$-centre); PTVV 2013 \emph{Publ.\ Math.\ IHES}~117;
  Lurie 2009 ``Derived Algebraic Geometry VI''.

- **AP933 -- $\Phi$ / bar--cobar output scope is $d$-dependent via
  the PTVV shift law (Vol~III Wave~12 cross-import, 2026-04-22,
  source `notes/wave12_a5_cy_d_stratification.tex`).**
  Vol~I's Theorem~A is a universal bar--cobar adjunction; its
  \emph{output $E_n$-level}, when specialised to CY$_d$-categories,
  is dimension-dependent via the Pantev--To\"en--Vaqui\'e--Vezzosi
  shift law
  $(d, \mathrm{shift}, E_n^{\mathrm{cl}}) \in
  \{(2,-2,E_2),(3,-1,E_1),(4,0,E_0),(5,+1,E_5\text{-Poisson})\}$.
  A CY$_d$-category with $d$-shifted symplectic structure produces
  a chiral-side output whose $E_n$-level is determined by
  $\mathrm{shift} + d = d - 2$ on the classical side (and
  correspondingly on the chiral side after Swiss-cheese grading).
  Vol~I results specialised to CY$_d$-categories \emph{carry this
  shift}; stating them without qualification applies the $d = 2$
  output level ($E_2$-chiral) uniformly across $d$, which fails
  at $d \ge 3$.
  **Wrong claim**: ``Theorem~A produces a chiral algebra''
  (unqualified) applied to a CY$_3$-category.
  **Precise error**: output is $E_1$-chiral, not $E_2$; the missing
  $E_2$-level sits on the centre of $\mathrm{Rep}$ (cf.\ AP932).
  **Counter**: when specialising Theorem~A to CY$_d$-categories,
  name the triple $(d, \mathrm{shift}, E_n^{\mathrm{cl}})$
  explicitly and carry the PTVV-shift-induced output level through
  to the chiral side. **Related**: AP933; AP932; AP929; AP924
  (iso-class vs parametrised KT formality); AP-CY5;
  PTVV 2013 \emph{Publ.\ Math.\ IHES}~117.

- **AP934 -- Single-citation canonicity: Costello--Li and
  Kontsevich--Tamarkin are two load-bearing inputs (Vol~III Wave~12
  cross-import, 2026-04-22, source
  `notes/wave12_a4_costello_2013_audit.tex`).**
  ``Canonical chiral algebra up to contractible choice'' is not a
  single theorem. Two inputs enter independently:
  (i) Costello--Li 2016 arXiv:1605.09930 \emph{holomorphic locality
  / homotopy renormalisation} of the factorisation algebra of
  observables of holomorphic Chern--Simons / B-model field theory,
  supplying the anomaly-free renormalisation datum on the chiral
  side;
  (ii) Kontsevich 2003 / Tamarkin 2003 \emph{$E_d$-formality},
  supplying the $E_d$-algebra up to contractible choice of Drinfeld
  associator. The canonical-up-to-contractible-choice property of
  the chiral algebra is the \emph{conjunction} of these two, not
  either alone. Kontsevich--Tamarkin delivers $E_d$-formality
  (classical-to-quantum); Costello--Li delivers holomorphic
  locality (renormalisation). Neither subsumes the other.
  **Wrong claim**: ``the chiral algebra is canonical by
  Kontsevich--Tamarkin formality'' (alone); or ``canonical by
  Costello--Li holomorphic locality'' (alone).
  **Precise error**: single-citation canonicity drops one of the
  two inputs; the $\mathrm{GRT}_1(\mathbb Q)$-torsor (AP924)
  coexists with the Costello--Li renormalisation base-point.
  **Counter**: every ``canonical chiral algebra'' claim names both
  Costello--Li \emph{and} Kontsevich--Tamarkin, with the
  contractible-choice space identified with the
  Grothendieck--Teichm\"uller torsor of associators.
  **Related**: AP934; AP924; AP923; AP-CY28;
  Costello--Li 2016; Kontsevich 2003; Tamarkin 2003; Willwacher
  2014; Francis 2013 \emph{Compositio}~149.

- **AP935 -- ``$X$ gives $Y$'' narration: construction requires
  the explicit arrow (Vol~III Wave~12 cross-import, 2026-04-22,
  source `notes/wave12_a1_bar_cobar_gelfand.tex`).**
  Prose of the form ``the shadow tower gives the bar Euler
  product'', ``Koszul duality gives the chiral algebra'', or ``the
  bar construction gives the coalgebra'' is narration, not a
  theorem; a construction is the datum of an explicit morphism.
  Chain-level: the cobar $\Omega^{\mathrm{ch}}$ or its avatar;
  $(\infty,1)$-categorical: the functor realising the
  identification. Vol~I prose that reads ``$X$ gives $Y$''
  without naming the arrow erases the structure of ``giving'' as
  the datum of a specified morphism; the reader cannot
  reconstruct which map is at stake. The CG voice
  (Chriss--Ginzburg, \emph{show don't tell}) demands the arrow.
  **Wrong claim**: ``$X$ gives $Y$'' (bare, no arrow).
  **Precise error**: erases the arrow.
  **Counter**: every ``$X$ gives $Y$'' sentence is rewritten as
  ``$\alpha \colon X \to Y$'' with $\alpha$ named inline or cited
  to an explicit construction; ``shadow tower gives'', ``bar
  construction gives'', ``Koszul duality gives'' are all audited
  under this discipline. **Related**: AP935; AP106; AP108;
  AP111; CG-voice discipline in CLAUDE.md; Pattern~273.

- **AP936 -- Denominator $=$ bar Euler: needs BOTH CY-A$_d$ and
  Vol~I anchor (Vol~III Wave~12 cross-import, 2026-04-22, sources
  `notes/wave12_a3_cy_a3_equivalence.tex` and
  `notes/wave12_b5_conifold_chiral_yangian.tex`).**
  The motivic Hall algebra generating-series identity
  (Kontsevich--Soibelman 2008 arXiv:0811.2435 wall-crossing;
  Schiffmann--Vasserot's CoHA product series) does NOT imply a
  naive series-level equality of the Borcherds-type denominator
  with the bar Euler product on the chiral side. The
  identification ``denominator $=$ bar Euler'' requires two
  load-bearing anchors: (i) a CY-A$_d$ equivalence at the
  appropriate $d$ (for $d = 3$: Vol~III Theorem CY-A$_3$,
  existence-and-$E_1$-rigidity, \emph{not} equivalence --- see
  AP938); (ii) the Vol~I shadow-tower / bar--cobar anchor at the
  specific class (Vol~I Theorem~A on the Koszul locus). The
  Vol~I analogue of Vol~III's AP-CY8: any series-level
  denominator-vs-bar-Euler identification must cite BOTH anchors
  at their specific scopes.
  **Wrong claim**: ``by the CoHA generating-series formula, the
  Borcherds denominator equals the bar Euler product''.
  **Precise error**: CoHA $= Y^+$ supplies only the positive-half
  factor (AP931); the Borcherds denominator carries negative-root
  and imaginary-root contributions that require a separate
  BKM-side anchor.
  **Counter**: every denominator-vs-bar-Euler identification cites
  (i) the CY-A$_d$ theorem at its existence-scope; (ii) Vol~I
  Theorem~A on the relevant Koszul / filtration locus; (iii) any
  class-specific anchor (e.g.\ class $\mathsf M$ AP937).
  **Related**: AP936; AP931; AP938; AP-CY8; Kontsevich--Soibelman
  2008; Borcherds 1995; Schiffmann--Vasserot 2013.

- **AP937 -- Class $\mathsf M$ $E_3$-bar $= 6^g$ at cohomology,
  infinite at chain level (Vol~III Wave~12 cross-import,
  2026-04-22, source `notes/wave12_a10_class_m_e3_bar.tex`).**
  For a class-$\mathsf M$ (Virasoro-type) chiral algebra on a
  curve of genus $g$, the $E_3$-bar construction yields a
  finite-dimensional cohomology space with
  $\dim H^\bullet(E_3\text{-bar}) = 6^g$ (Oberdieck 2019
  \emph{Geom.\ Topol.}~23 generating-function count, verified at
  $g = 1, 2, 3$ by direct descent-spectral-sequence computation).
  At chain level the $E_3$-bar complex is
  infinite-dimensional. Vol~I's class-$\mathsf M$ shadow-tower
  treatment must state this explicitly. AP885's earlier ``$6^g$
  overcount'' retraction applies to the chain-level reading, not
  the cohomological one; the $6^g$ value is exact at cohomology.
  **Wrong claim**: ``the class-$\mathsf M$ $E_3$-bar is infinite''
  (bare); or ``the class-$\mathsf M$ $E_3$-bar is $6^g$'' (bare).
  **Precise error**: collapses chain-level and cohomological
  lanes; transports infinity to cohomology or finiteness to
  chain level.
  **Counter**: every class-$\mathsf M$ $E_3$-bar statement names
  the lane (chain-level infinite vs cohomological $6^g$), with
  Oberdieck 2019 as primary for the cohomological count.
  **Related**: AP937; AP885 (superseded-in-part); AP914
  (Oberdieck factor-of-two); AP-CY18; Oberdieck 2019;
  \texttt{chapters/theory/shadow\_tower\_higher\_coefficients.tex}.

- **AP938 -- Two-stage factorisation is existence-and-rigidity,
  NOT equivalence, at $d \ge 3$ (Vol~III Wave~12 cross-import,
  2026-04-22, source `notes/wave12_a3_cy_a3_equivalence.tex`).**
  Vol~III's Theorem CY-A$_3$ proves the two-stage factorisation
  $\Phi_3 = \mathrm{BarCobar}^{\mathrm{ch}} \circ Q^{\mathrm{CY}_3}$
  of the CY$_3$-to-chiral functor as
  \emph{existence-and-$E_1$-rigidity} of the chiral-side
  $(\infty,1)$-stage, not as an equivalence of
  $(\infty,1)$-categories. Existence: every CY$_3$-category in the
  $(\infty,1)$-image of $Q^{\mathrm{CY}_3}$ admits a chiral
  $E_1$-stage unique up to contractible choice of the
  Kontsevich--Tamarkin / Costello--Li renormalisation datum
  (AP934). Rigidity: the $E_1$-stage is rigid in its own moduli.
  \emph{Not} established: that $\Phi_3$ is fully faithful or
  essentially surjective as an $(\infty,1)$-functor. Vol~I's
  CY-A theorem, when stated for $d \ge 3$, must carry this
  existence-and-rigidity scope, not the stronger equivalence
  scope that applies at $d \le 2$.
  **Wrong claim**: ``Theorem~CY-A at $d = 3$ is an equivalence of
  $(\infty,1)$-categories''.
  **Precise error**: upgrades existence-and-rigidity to
  equivalence; conflates Vol~III's proved scope with the
  stronger $d \le 2$ form.
  **Counter**: every cross-volume citation of CY-A displays the
  scope at each $d$: $d \le 2$ equivalence;
  $d \ge 3$ existence-and-$E_1$-rigidity.
  **Related**: AP938; AP933; AP929; AP934; AP-CY5; AP-CY6;
  AP-CY30;
  \texttt{chapters/examples/cy\_d\_kappa\_stratification.tex}
  CY-A scope declaration; Vol~III Theorem CY-A$_3$.

- **AP-CAT-15 -- Wave-12 cross-volume import anchor map (2026-04-22;
  cross-volume concordance with Vol~III wave-12 synthesis).** Ten
  Vol~III wave-12 cross-programme antipatterns cross-cut into
  Vol~I's bar--cobar / shadow-tower / MC-element /
  modular-characteristic / derived-centre domain and are inscribed
  here as AP929--AP938 with direct reference to the Vol~III
  wave-12 source files. **Source mapping** (Vol~III wave-12
  $\to$ Vol~I AP):
  `wave12_a1_phi_functor_foundations.tex` $\to$ AP929 (Pattern~273
  formalisation);
  `wave12_a2_kappa_invariants_universal_borcherds.tex` $\to$ AP930
  (bare $\kappa$ forbidden);
  `wave12_a11_coha_y_plus_vs_w_infty.tex` $\to$ AP931 (CoHA $=Y^+$);
  `wave12_a8_en_hierarchy_kapranov.tex` $\to$ AP932 ($E_n$ on wrong
  object);
  `wave12_a5_cy_d_stratification.tex` $\to$ AP933 (output scope
  $d$-dependent via PTVV);
  `wave12_a4_costello_2013_audit.tex` $\to$ AP934 (single-citation
  canonicity);
  `wave12_a1_bar_cobar_gelfand.tex` $\to$ AP935 (``$X$ gives $Y$''
  narration);
  `wave12_a3_cy_a3_equivalence.tex` + `wave12_b5_conifold_chiral_yangian.tex`
  $\to$ AP936 (denominator $=$ bar Euler);
  `wave12_a10_class_m_e3_bar.tex` $\to$ AP937 (class $\mathsf M$
  $6^g$ cohomological);
  `wave12_a3_cy_a3_equivalence.tex` $\to$ AP938
  (existence-and-rigidity, not equivalence).
  **Rule** (cross-volume discipline, orthogonal to AP-CAT-14):
  AP-CAT-14 covers Vol~III wave-14 session-corrections
  (AP-CY160--172 $\to$ AP923--928, numerical-value discipline);
  AP-CAT-15 covers the Vol~III wave-12 synthesis
  (grammar-level scope separation). The two blocks coexist and do
  not subsume each other. **Hook regexes** (advisory, deferred per
  AP-CAT-13):
  (i) AP929 --- trigger on ``bar--cobar adjunction''
  $+$ ``Theorem~A'' without ``chain-level'' or
  ``$(\infty,1)$-categorical'' qualifier within 10 lines;
  (ii) AP930 --- trigger on bare ``$\kappa$'' without a
  $_\bullet$-subscript (re-uses the HZ-7 hook);
  (iii) AP931 --- trigger on ``CoHA''$+$``$\cW_{1+\infty}$''
  co-occurrence without ``$Y^+$'' or ``$\mathrm{ev}_\lambda$''
  within 10 lines;
  (iv) AP932 --- trigger on ``$E_2$-chiral'' applied to a chiral
  algebra at a $d \ge 3$ ambient without
  ``$Z(\mathrm{Rep})$'' qualifier;
  (v) AP933 --- trigger on ``Theorem~A'' at CY$_d$ ambient
  without $(d, \mathrm{shift}, E_n^{\mathrm{cl}})$ triple named;
  (vi) AP934 --- trigger on ``canonical chiral algebra'' without
  both Costello--Li and Kontsevich--Tamarkin named;
  (vii) AP935 --- trigger on ``gives'' following a named
  construction (shadow tower, bar construction, Koszul duality)
  in manuscript prose;
  (viii) AP936 --- trigger on ``denominator''+``bar Euler''
  co-occurrence without both CY-A$_d$ and Vol~I Theorem~A cited
  within 20 lines;
  (ix) AP937 --- trigger on ``$6^g$'' in a class-$\mathsf M$
  context without ``at cohomology'' qualifier;
  (x) AP938 --- trigger on ``Theorem~CY-A''+``$d = 3$''+
  ``equivalence'' co-occurrence.
  **Related**: AP-CAT-15; AP929--AP938; AP-CAT-14; AP-CAT-13;
  AP-CAT-12; Vol~III wave-12 source files above; canonical
  preamble $\kappa_\bullet$-subscript rows.


### 2026-04-22: 15 cross-volume AP entries added from Vol III Wave 11-19

Fifteen cross-volume entries mirror the Vol III Wave 11-19 error/retraction
catalogue (six errors E\_A--E\_F, nine retractions R1--R9) into the Vol I
registry with Vol I-specific manifestations. Each entry carries a shared
signature (so an agent moving between volumes identifies the same underlying
pattern), a Vol I-specific ghost/correction reading (so the Vol I mathematics
remains primary), and a cross-reference to the Vol III AP-CY primary form.
AP tags AP939--AP953 follow the next available numbering after the two-stage
/ six-routes / dimension-stratified append (AP929--AP938).

- **AP939 (Kuenneth multiplicative on products; bare $\kappa$ treats
  $\kappa_{\mathrm{cat}}(K3 \times E) = 2$ as a total-space value when the
  fibre value is meant).** Vol I's landscape census (`chapters/examples/
  landscape_census.tex`), K3 bar-complex examples
  (`chapters/examples/k3_bar_complex.tex`), and any cross-citation to
  Vol III that names $\kappa(K3 \times E)$ without subscript silently
  presents the K3-fibre value $\chi(\cO_{K3}) = 2$ as a total-space
  invariant; categorical Euler on the total space is Kuenneth-multiplicative:
  $\kappa_{\mathrm{cat}}(K3 \times E) = \chi(\cO_{K3}) \cdot \chi(\cO_E) =
  2 \cdot 0 = 0$. **Wrong claim**: "$\kappa(K3 \times E) = 2$" or
  "$\chi(\cO_{K3 \times E}) = 2$" without naming whether the subscript is
  $\kappa_{\mathrm{cat}}$ (on total space) or $\kappa_{\mathrm{fiber}}$ (on
  fibre). **Precise error**: conflates two distinct $\kappa_\bullet$
  invariants — $\kappa_{\mathrm{cat}}$ on the total space of the product,
  which is Kuenneth-multiplicative and vanishes because $\chi(\cO_E) = 0$,
  with $\kappa_{\mathrm{fiber}}$ on the K3 fibre, which equals $\chi(\cO_{K3})
  = 2$. The two values live on different invariants of different objects;
  neither is a "correction" of the other. **Correction**: every
  $\kappa(K3 \times E)$ reference in Vol I subscripts explicitly:
  $\kappa_{\mathrm{cat}}(K3 \times E) = 0$, $\kappa_{\mathrm{fiber}}(K3) = 2$,
  $\kappa_{\mathrm{BKM}}(\mathbf H_{\Delta_5}) = 5$ (paramodular $\Phi_{10} =
  \Delta_5^2$), $\kappa_{\mathrm{ch}}^K(K3 \times E) = 3$ (Heisenberg-level
  additive). Vol I's landscape-census numerics and archetype-$\mathbf B$
  $K = 8$ row depend on naming $\kappa_\bullet$ at every site. **Primary**:
  Vol III AP-CY68 / AP234 canonical ruling; Vol III `chapters/examples/
  cy_d_kappa_stratification.tex` canonical table; canonical preamble rows
  63, 56. **Detection heuristic**: grep `'\\kappa([^_]*K3.*\\times.*E|K3 \\times E).*=.*2'`;
  any bare $\kappa$ adjacent to $K3 \times E$ without a subscript in 10 lines
  triggers Gate 0. **Cross-ref**: Vol III AP-CY for primary form of
  Kuenneth-multiplicative total-space discipline; AP929 / AP930 two-stage
  framing; canonical preamble row 56 ($\kappa_\bullet$ indexing on
  $K3 \times E$); V2-AP68 / V2-AP69 (Vol II partner).

- **AP940 (Central charge $-1312/11$ at $c = -214$ needs two independent
  verification routes, not $-14432/121$).** Vol I chapters invoking the
  class-$\mathcal S$ / Gaiotto $\mathcal T[A_1, \Sigma_{0, 24}]$ central-charge
  identity routinely cite one value without a second route check. The Wave 14
  error $-14432/121$ (from a dropped simplification) was corrected to
  $-1312/11$ in Vol III; the same verification is load-bearing for Vol I's
  landscape-census rows at the K3-BKM family, since $c_{2d} = -1312/11$ is the
  central charge entering the $\kappa + \kappa^! = 8$ archetype-$\mathbf B$
  row on the $r(z)$-family landscape. **Wrong claim**: a Vol I reference uses
  the Wave 14 value $c_{2d} = -14432/121$ or quotes $c = -214$ with a bare
  central-charge ratio without naming the $\Sigma_{0, 24}$ / Gaiotto curve
  anchor. **Precise error**: single-route citation for a load-bearing
  structure-constant whose Wave 14 value was retracted; AP274-style "citation
  as derivation" for a load-bearing constant. **Correction**: every Vol I
  central-charge identity in the K3-BKM family is verified against two
  independent routes (Shapere--Tachikawa 2008 universal $c_{4d} = (5n - 13)/6$
  formula; direct Sugawara-denominator recomputation on the class-$\mathcal S$
  torus), and the value $(c_{4d}, c_{2d}) = (107/6, -214)$ for
  $\mathcal T[A_1, \Sigma_{0, 24}]$ is named with Gaiotto-curve anchor.
  **Primary**: Vol III Wave 15 Gaiotto verdict; canonical preamble row 17;
  Shapere--Tachikawa 2008. **Detection heuristic**: grep `'c_{2d}\\s*=\\s*-14432'`
  or `'-14432/121'`; any hit is a retracted-value witness. **Cross-ref**:
  Vol III AP-CY for primary form (E\_B); AP274 (citation as derivation
  discipline for load-bearing constants); canonical preamble row 17.

- **AP941 (Eta-power $\eta^{-48}$ is Heisenberg--Mukai identity at the K3
  Hall--Drinfeld site, not a Virasoro minimal-model relation).** Vol I's K3
  chiral Hall--Drinfeld example (`chapters/examples/k3_chiral_hall_drinfeld.tex`
  or equivalent) invokes $\eta^{-48}$ as the character of the $V_{24}$-module
  iterated Drinfeld--Sokolov reduction; a careless reading frames this as a
  Virasoro minimal-model expansion (where $\eta$-powers index level-shift
  families). The two are not the same: $\eta^{-48}$ is the Heisenberg--Mukai
  lattice partition function on $\mathrm{II}_{1, 1} \otimes \cO_{K3}$, read
  as the positive-half-vertex-operator trace on the Mukai lattice.
  **Wrong claim**: "$\eta^{-48}$ is a Virasoro minimal-model level character"
  or "the $\eta^{-48}$ factor arises from DS reduction of Virasoro at level
  $-48$". **Precise error**: conflates two distinct vertex-algebra
  presentations of the same generating function; the Heisenberg--Mukai
  reading (24 copies of the Heisenberg VOA on the Mukai lattice) gives the
  $\eta^{-48}$ factor as $\prod_n (1 - q^n)^{-24}$ characters times inverse
  $\Delta$, and the Virasoro minimal-model reading has no free-field
  Heisenberg factor. **Correction**: every Vol I $\eta^{-48}$ reference names
  the Heisenberg--Mukai lattice VOA and the iterated Drinfeld--Sokolov chain
  $V_{24} = \mathrm{DS}^{24}(\mathrm{Heis}_{\mathrm{II}_{1,1}})$ explicitly.
  **Primary**: Vol III Wave 15--17 direct character-match computation;
  Mukai 1988; Frenkel--Kac construction. **Detection heuristic**: grep
  `'\\eta\\^{-48}|eta\\^{-48}'`; any occurrence within 10 lines of "Virasoro
  minimal" or "minimal model" without "Heisenberg" or "Mukai" clarification
  triggers Gate 0. **Cross-ref**: Vol III AP-CY for primary form (E\_C, R5);
  AP284 (equivalence claim without inscribed implication structure); R5
  retraction (direct $\chi_{V_{24}}$ match $\to \eta^{-48}$).

- **AP942 (Macdonald framework inapplicable to non-minimal targets).** Vol I
  material using the Macdonald difference-operator framework to parametrise
  R-matrices or bar-weight structures occasionally extends the framework to
  non-minimal-model targets (W-algebras at non-admissible level, logarithmic
  triplet VOAs, parafermion cosets) where the Macdonald polynomial basis
  does not diagonalise the ambient Hamiltonian. **Wrong claim**: "the
  Macdonald framework applies to any $W$-algebra at non-critical level" or
  "Macdonald symmetric functions diagonalise the bar-pairing on class
  $\mathbf L$ / $\mathbf C$ / $\mathbf M$". **Precise error**: the Macdonald
  framework is a theorem for Cherednik double-affine Hecke algebras and
  their specialisations to admissible-level affine Lie algebras; for
  non-minimal-model targets (non-admissible $W$-algebras, triplet VOAs) the
  positive-part Yangian structure is not diagonal in the Macdonald basis,
  and the pentagon-coboundary calculation $\phi^{(n)}$ silently relies on
  admissibility. **Correction**: every Vol I Macdonald-framework invocation
  names the admissibility constraint — minimal target, admissible level
  $k + h^\vee = p / q$ with $(p, q) = 1$, finite-dimensional bar at each
  weight. Non-minimal targets use a different framework (Kac--Roan 2014
  integrable-system residue; Feigin--Frenkel centre on $\mathrm{Op}$-side).
  **Primary**: Vol III Wave 18 Gaiotto verdict on admissibility; Macdonald
  1995; Cherednik 2005. **Detection heuristic**: grep `'Macdonald'` within
  10 lines of "W-algebra" / "triplet" / "logarithmic" without "admissible"
  or "minimal" qualifier triggers Gate 0. **Cross-ref**: Vol III AP-CY for
  primary form (E\_D); AP275 (scope-qualifier without obstruction); AP289
  (defensive scoping as silent downgrading); canonical preamble row 44
  (admissibility-congruence $D_n$).

- **AP943 (Eight-form position $\ne$ Borcherds weight; weights are
  $(5, 2, 1, 1, 1/2, 1, 1/4, 0)$).** Vol I material citing the
  Gritsenko--Clery eight-form catalogue occasionally equates "eight-form
  position" (indexed $1, 2, 3, 4, 5, 6, 7, 8$ in the Gritsenko enumeration)
  with "Borcherds weight" (the weight of the Borcherds automorphic lift,
  entering $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$). The weights are not
  monotone in position: $(w_1, \ldots, w_8) = (5, 2, 1, 1, 1/2, 1, 1/4, 0)$
  with Fourier coefficients $c_N(0) \in \{10, 4, 2, 2, 1, 2, 1/2, 0\}$,
  reproducing $\kappa_{\mathrm{BKM}} = c_N(0)/2 = (5, 2, 1, 1, 1/2, 1, 1/4, 0)$.
  **Wrong claim**: "the $k$-th Gritsenko--Clery form has Borcherds weight
  $k$" or "Borcherds weight climbs with position in the eight-form
  catalogue". **Precise error**: confuses enumeration index with a
  mathematical invariant; the enumeration is by $\mathrm{Sp}_4(\Z)$-cover
  refinement (integral, half-integral, quarter-integral), not by weight.
  **Correction**: every Vol I citation of the eight-form catalogue names
  the cover group ($\mathrm{Sp}_4(\Z)$ for integral, $\mathrm{Mp}_4$ for
  half-integral, $\widetilde{\mathrm{Mp}}_4$ for quarter-integral) and
  lists the weight next to each form. **Primary**: Gritsenko--Clery 2018;
  Vol III `chapters/examples/cy_d_kappa_stratification.tex` Theorem
  `thm:borcherds-weight-kappa-BKM-universal`. **Detection heuristic**: grep
  `'eight.*form|8.*form'` within 10 lines of "Borcherds weight" without the
  explicit weight triple triggers Gate 0. **Cross-ref**: Vol III AP-CY for
  primary form (E\_E); AP274 (citation as derivation for load-bearing
  constants); canonical preamble row 28 (K3-BKM Weyl denominator).

- **AP944 (Twined convention ambiguity $c_2(0) \in \{4, 8\}$ requires
  declaration at every cross-citation).** Vol I material citing twined
  Fourier coefficients $c_N(0)$ in the Borcherds-weight identity
  $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$ carries an implicit convention
  choice at $N = 2$: $c_2(0) = 4$ under one twined convention (Gritsenko
  1999) and $c_2(0) = 8$ under another (Cheng--Duncan 2013 umbral normalisation
  factor of 2). Vol I cross-citations to Vol III that do not declare the
  convention silently switch between $\kappa_{\mathrm{BKM}}(\Phi_2) = 2$
  and $\kappa_{\mathrm{BKM}}(\Phi_2) = 4$. **Wrong claim**: "$c_2(0) = 4$"
  or "$c_2(0) = 8$" without naming Gritsenko or Cheng--Duncan convention.
  **Precise error**: convention-dependent structure-constant cited bare,
  ambient-qualifier drop (AP281). **Correction**: every Vol I reference to
  twined $c_N(0)$ declares the convention and the Borcherds automorphic
  form used: "Gritsenko 1999 twined coefficients ($c_2(0) = 4$)" or
  "Cheng--Duncan 2013 umbral coefficients ($c_2(0) = 8$)". **Primary**:
  Vol III canonical $c_N(0)$ table (Gritsenko convention); Cheng--Duncan
  2013 umbral moonshine. **Detection heuristic**: grep
  `'c_2\\(0\\)|c_2\\^{(0)}'`; any hit without "Gritsenko" or "Cheng--Duncan"
  in 10 lines triggers Gate 0. **Cross-ref**: Vol III AP-CY for primary
  form (E\_F); AP281 (ambient-qualifier drop); canonical preamble row 28
  (K3-BKM denominator).

- **AP945 ($\widehat{\mathfrak{sl}_3}$ $\to$ $F_3$ Feingold--Frenkel at Vol I
  genus-1 BKM sites; $F_3$ is the real-root subalgebra).** Vol I's genus-1
  BKM chapters (classifying real-root subalgebras of $\mathfrak g_{\Delta_5}$
  and the Fake Monster) initially framed the real-root subalgebra as
  $\widehat{\mathfrak{sl}_3}$; the correct identification is the Feingold--
  Frenkel hyperbolic Kac--Moody algebra $F_3$ on the rank-3 hyperbolic
  lattice $A_1^{(1)} \oplus (-1)$. $F_3$ has genuine imaginary roots of
  Humbert-divisor type, matching the K3 fibre $(2, 1, 1)$-lattice, whereas
  $\widehat{\mathfrak{sl}_3}$ is affine untwisted and has no imaginary-root
  content. **Wrong claim**: "the real-root subalgebra of $\mathfrak g_{\Delta_5}$
  is $\widehat{\mathfrak{sl}_3}$" or "the genus-1 BKM real-root piece is
  affine $\mathfrak{sl}_3$". **Precise error**: conflates an affine Lie
  algebra (rank-2 Cartan plus central element, no imaginary roots) with a
  rank-3 hyperbolic Kac--Moody algebra (Cartan matrix of indefinite type,
  imaginary-root Weyl chamber non-empty). **Correction**: every Vol I
  genus-1 BKM real-root site cites Feingold--Frenkel 1983 $F_3$ with
  rank-3 hyperbolic Cartan matrix $\begin{pmatrix} 2 & -2 & 0 \\ -2 & 2 & -1
  \\ 0 & -1 & 2 \end{pmatrix}$ and names the Humbert-divisor matching on
  the K3 fibre. **Primary**: Feingold--Frenkel 1983 *Math Ann* 263;
  Gritsenko--Nikulin 1998 Table 2. **Detection heuristic**: grep
  `'\\\\hat\\{\\\\mathfrak\\{sl\\}_3\\}|\\\\widehat\\{\\\\mathfrak\\{sl\\}_3\\}'`
  within 10 lines of "genus-1 BKM" or "$\\Delta_5$ real-root" without
  "Feingold" triggers Gate 0. **Cross-ref**: Vol III AP-CY for primary
  form (R1); canonical preamble row 28 (K3-BKM denominator); AP-CY
  Feingold--Frenkel entry.

- **AP946 ($L_{-6}(\mathfrak{e}_8)$ $\to$ iterated Drinfeld--Sokolov
  reduction; $V_{24}$ is not $L_{-6}(\mathfrak{e}_8)$ on the nose).** Vol I's
  K3 chiral Hall--Drinfeld and BRST-ghost-charge chapters initially framed
  the 24-dimensional CFT $V_{24}$ supporting the Heisenberg--Mukai factor as
  $L_{-6}(\mathfrak{e}_8)$, the $\mathfrak{e}_8$ affine VOA at level $-6$.
  The correct construction is iterated Drinfeld--Sokolov reduction from a
  higher-rank source: $V_{24} = \mathrm{DS}_{\lambda}(L_k(\mathfrak{g}))$
  for a sequence $(\mathfrak g, k, \lambda) = (\mathfrak{e}_8, -5, \mathrm{sub-regular})
  \to (\mathfrak{e}_6, \ldots) \to \ldots$ with 24 singular iterations,
  matching the $\eta^{-48}$ character and the $V_{24}$ central charge
  $c = 24$ exactly. **Wrong claim**: "$V_{24} \simeq L_{-6}(\mathfrak{e}_8)$"
  or "the 24-dimensional ghost factor is affine $\mathfrak{e}_8$ at level
  $-6$". **Precise error**: conflates an affine VOA at a specific level
  with the iterated DS output of a higher-rank-source chain; $L_{-6}(\mathfrak{e}_8)$
  has central charge $c = -992/11 \ne 24$, so the identification fails at
  the central-charge level. **Correction**: every Vol I $V_{24}$ reference
  cites the iterated DS chain with source, level, and nilpotent; the
  character match $\chi_{V_{24}} = \eta^{-48}$ is the three-path verification
  anchor (Sugawara-denominator recomputation; iterated-DS character
  contraction; Borcherds lift on the Mukai lattice). **Primary**: Vol III
  Wave 15 iterated-DS chain; Kac--Roan 2000; Adamovic--Milas 1997.
  **Detection heuristic**: grep `'L_{-6}.*e_8|L_{-6}\\(\\\\mathfrak\\{e\\}_8\\)'`
  within 10 lines of "$V_{24}$" without "iterated" or "DS" triggers Gate 0.
  **Cross-ref**: Vol III AP-CY for primary form (R2, R5); AP941
  ($\eta^{-48}$ Heisenberg--Mukai); canonical preamble row 52
  (chiral-Hochschild $e_k$ motivic home).

- **AP947 ($\kappa_{\mathrm{BKM}}$ not additive; universal $c_N(0)/2$
  forces scope-split in BRST ghost-charge computation).** Vol I's BRST
  ghost-charge computation (the three-factor trace identity
  $\mathrm{tr}_{\mathrm{ghost}}(Q_{\mathrm{BRST}}^2) = \mathrm{tr}_{\mathrm{Pentagon}}
  = \omega_{\mathrm{Borcherds}} = c_N(0)/2$) initially attempted to
  factorise additively as
  $\kappa_{\mathrm{BKM}}(\Phi_N) = \kappa_{\mathrm{ch}}(X) + \chi(\cO_{\mathrm{fiber}})$
  across the $N$-family. This relation is FALSE at every
  $N \in \{1, 2, 3, 4, 6\}$ (Vol III primary form: at $N = 1$, left $= 5$,
  right $= \kappa_{\mathrm{ch}}(K3 \times E) + \chi(\cO_E) = 0 + 0 = 0$;
  at $N = 2$, left $= 4$, right $= 1$). The universal formula is
  $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$ (Borcherds weight). **Wrong
  claim**: "$\kappa_{\mathrm{BKM}} = \kappa_{\mathrm{ch}} + \chi(\cO_{\mathrm{fiber}})$"
  or additive cross-volume $\kappa_\bullet$ factorisation. **Precise
  error**: additive decomposition of a Borcherds-lift Fourier coefficient
  into CY-categorical and fibre-Euler pieces has no lift-functoriality
  justification; the Borcherds lift is a multiplicative automorphic-product
  functor, not an additive one. **Correction**: the Vol I three-factor
  trace identity reads with scope-split — the ghost trace computes
  $c_N(0)/2$ via Borcherds lift on the Mukai lattice (multiplicative),
  and the separate $\kappa_{\mathrm{ch}}$ / $\kappa_{\mathrm{fiber}}$ bookkeeping
  is the CY-categorical denominator contribution on the $\Phi$-functor
  side. The two scopes meet at $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$
  as a universal Borcherds weight theorem, not as an additive
  decomposition. **Primary**: Vol III `chapters/examples/
  cy_d_kappa_stratification.tex` Theorem `thm:borcherds-weight-kappa-BKM-universal`;
  Borcherds 1995; Gritsenko 1999. **Detection heuristic**: grep
  `'kappa_{BKM}.*=.*kappa_{ch}.*\\+'`; any hit triggers Gate 0.
  **Cross-ref**: Vol III AP-CY for primary form (R3); canonical preamble
  row 56 ($\kappa_\bullet$ indexing on $K3 \times E$); AP939 (Kuenneth
  multiplicative).

- **AP948 (Fake Monster at $d = 5$, not $d = 3$; $K3 \times K3 \times E$
  CY-5 home).** Vol I genus-5 / Andreotti--Mayer / universal-closure
  chapters that cite the Fake Monster BKM as a $d = 3$ object (on
  $K3 \times E$) conflate it with the K3-BKM paramodular $\Delta_5$
  (Cartan rank $3$ on $\Lambda^{2, 1}_{\mathrm{II}}$). The Fake Monster's
  Cartan rank is $26$ on $\mathrm{II}_{25, 1}$ (Borcherds 1992), and its
  natural CY-categorical home is the CY-5 product $K3 \times K3 \times E$
  with lattice signature $(3, 19) + (3, 19) + (1, 1) = (7, 39)$ admitting
  primitive restriction to $\mathrm{II}_{26, 2}$ for the $\Phi_{12}$
  Borcherds lift. **Wrong claim**: "the Fake Monster BKM is the Vol III
  $d = 3$ $K3 \times E$ output" or "$\Phi_{12}$ lifts from the $K3 \times E$
  Mukai lattice". **Precise error**: confuses two BKMs at two different
  CY-dimensions indexed by two different Borcherds lifts: $\Phi_{12}$
  (Fake Monster, $d = 5$) and $\Phi_{10} = \Delta_5^2$ (K3-BKM, $d = 3$).
  **Correction**: every Vol I Fake Monster reference cites $d = 5$ with
  $K3 \times K3 \times E$ input and $\Phi_{12}$ Borcherds lift; K3-BKM
  references at $d = 3$ use $\Delta_5$ or $\Phi_{10}$. **Primary**:
  Borcherds 1992 *Invent Math* 109; Vol III Dimension-stratified siblings
  catalogue; canonical preamble rows 20--22. **Detection heuristic**: grep
  `'Fake.*Monster|fake.*monster'` within 10 lines of "$K3 \\times E$" or
  "$d\\s*=\\s*3$" without "$d = 5$" or "$K3 \\times K3$" triggers Gate 0.
  **Cross-ref**: Vol III AP-CY for primary form (R4); AP931 (Fake-Monster
  BKM dimension discipline); canonical preamble rows 20--22.

- **AP949 (Direct $\chi_{V_{24}}$ match $\to \eta^{-48}$ is an
  identity-by-character-matching not by DS-chain-construction).** Vol I
  K3 chiral Hall--Drinfeld sites cite $\chi_{V_{24}} = \eta^{-48}$ as a
  direct character identity; the underlying mathematics is a three-path
  verification: (i) iterated-DS character contraction from
  $L_k(\mathfrak{e}_8)$ source; (ii) Sugawara-denominator recomputation on
  the Heisenberg--Mukai lattice; (iii) Borcherds lift on the Mukai lattice.
  Treating the character identity as a stand-alone fact loses the three-path
  derivation and the $V_{24}$ construction anchor. **Wrong claim**:
  "$\chi_{V_{24}} = \eta^{-48}$ by inspection" or "the direct character
  match pins $V_{24}$". **Precise error**: load-bearing identity cited
  without derivation (AP274 citation as derivation); direct-match framing
  hides the three-path verification. **Correction**: every Vol I site
  citing the $\eta^{-48}$ character identity names the three verification
  paths and the iterated-DS construction of $V_{24}$. **Primary**: Vol III
  Wave 15--17 three-path verification. **Detection heuristic**: grep
  `'\\\\eta\\^{-48}'` within 10 lines of "direct character" / "by inspection"
  without "three-path" or "iterated DS" triggers Gate 0. **Cross-ref**:
  Vol III AP-CY for primary form (R5); AP941 ($\eta^{-48}$ Heisenberg--Mukai);
  AP946 (iterated DS chain for $V_{24}$); AP274 (citation as derivation).

- **AP950 (Gaiotto curve $\Sigma_{0, 24}$ with 24 punctures, not
  $\Sigma_{2, 0}$ genus-2 closed).** Vol I class-$\mathcal S$ / AGT chapters
  referencing the Gaiotto curve for the K3-BKM landscape initially framed
  the curve as the genus-2 closed Riemann surface $\Sigma_{2, 0}$. The
  correct Gaiotto curve is $\Sigma_{0, 24}$: genus-0 with 24 punctures,
  each carrying a simple $A_1$ puncture (regular / $\mathrm{SU}(2)$-flavour).
  The puncture count $24$ matches the $24 I_1$ nodes of the pentagon
  construction and the $\chi(K3) = 24$ Euler characteristic. **Wrong
  claim**: "the Gaiotto curve for $\mathcal T[A_1, K3]$ is $\Sigma_{2, 0}$
  closed genus-2" or "the K3-BKM class-$\mathcal S$ construction uses a
  closed genus-2 surface". **Precise error**: puncture count / genus
  confusion; the genus-2 surface has no punctures and carries a different
  dimension of the moduli of flat connections ($6g - 6 + 2n = 6$ vs
  $6g - 6 + 2n = 42$). **Correction**: every Vol I Gaiotto-curve reference
  cites $\Sigma_{0, 24}$ with 24 punctures and the $A_1$-type puncture
  assignment. **Primary**: Vol III Wave 15 Gaiotto verdict; canonical
  preamble row 18 (trinion/tube count $(n_v, n_h) = (63, 88)$);
  Gaiotto 2009 *JHEP* 1208; Chacaltana--Distler 2010. **Detection heuristic**:
  grep `'\\\\Sigma_{2,0}|\\\\Sigma_{2}'` within 10 lines of "Gaiotto" or
  "K3-BKM" triggers Gate 0. **Cross-ref**: Vol III AP-CY for primary form
  (R6); canonical preamble rows 17--18 ($(c_{4d}, c_{2d})$ and trinion
  count); AP940 (two-route central-charge verification).

- **AP951 (Two-stage factorisation $\Phi_d = \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}
  \circ \Phi^{\mathrm{FA}}_d$; Koszul Reflection $K$ becomes Stage 1).** Vol I's
  Koszul Reflection construction $K$ (on which the five-archetype
  $\mathsf{G}/\mathsf{L}/\mathsf{C}/\mathsf{M}/\mathsf{B}$ landscape
  complementarity is built) is \emph{Stage 1} of the Vol III two-stage
  factorisation of $\Phi$: $\Phi_d = \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}
  \circ \Phi^{\mathrm{FA}}_d$, with $K$ identified as $\Phi^{\mathrm{FA}}_d$
  on the Koszul-self-dual locus. The Vol I bar--cobar adjunction (Theorem A)
  is the Stage-2 $E_1$-chiral shadow on the reference curve $C$. **Wrong
  claim**: "Koszul Reflection $K$ is the CY-to-chiral functor" or "$\Phi_d$
  is a single-stage functor landing directly in $E_1$-chiral". **Precise
  error**: collapses two distinct $(\infty, 1)$-categorical operations
  (factorisation-algebra assignment $\Phi^{\mathrm{FA}}_d$ and specialisation
  $\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}$) into one, hiding the
  family-parameter $(\Sigma_{d-1}, C)$. **Correction**: every Vol I
  reference to the CY-to-chiral assignment names the two stages and the
  family parameter; Koszul Reflection $K$ is named as Stage 1 on the
  Koszul-self-dual locus. **Primary**: Vol III `chapters/theory/cy_to_chiral.tex`
  two-stage construction; AP929 (Vol I partner); Kontsevich--Tamarkin
  $E_d$-formality. **Detection heuristic**: grep `'\\\\Phi_d|Phi_d'` in Vol I
  prose without "two-stage" or "Stage 1 / Stage 2" within 10 lines triggers
  Gate 0 (handled by AP929 hook). **Cross-ref**: Vol III AP-CY for primary
  form (R7); AP929 (Vol I two-stage framing partner); AP930 (six routes as
  distinct constructions).

- **AP952 (PTVV shifted-symplectic chapter extends to $d = 5$ Poisson-$E_5$
  corner closure).** Vol I's PTVV shifted-symplectic chapter
  (`chapters/theory/ptvv_shifted_symplectic.tex` or equivalent) initially
  enumerated the $(d, \mathrm{shift}, E_n^{\mathrm{cl}})$ table only to
  $d = 4$. The correct extension closes the table at $d = 5$: $(d, \mathrm{shift},
  E_n^{\mathrm{cl}}) = (5, -4, \mathrm{Poisson}\text{-}E_5)$, with the
  Poisson-$E_5$ shift coming from the PTVV $(-4)$-shifted symplectic
  structure on the CY-5 moduli and the Dunn--Lurie additivity
  $E_3 \otimes E_2 = E_5$ after specialisation. **Wrong claim**: "PTVV
  table closes at $d = 4$" or "$d = 5$ is beyond the PTVV framework".
  **Precise error**: incomplete dimensional enumeration; the PTVV machinery
  extends to all $d$ with shift $-2d + 2$ on the CY-$d$ moduli. **Correction**:
  every Vol I PTVV table enumeration includes the $d = 5$ row with
  Poisson-$E_5$ corner closure and $-4$-shift; the full sequence is
  $(d, \mathrm{shift}, E_n^{\mathrm{cl}}) = (1, 0, E_0), (2, -2, E_1), (3,
  -2, \mathrm{Poisson}\text{-}E_3), (4, -2, \mathrm{Poisson}\text{-}E_4),
  (5, -4, \mathrm{Poisson}\text{-}E_5)$. **Primary**: Pantev--Toen--Vaquie--
  Vezzosi 2013 *Publ IHES* 117; Vol III Wave 18 PTVV extension. **Detection
  heuristic**: grep `'PTVV.*d.*=.*4'` without "$d = 5$" in the same
  enumeration triggers Gate 0. **Cross-ref**: Vol III AP-CY for primary
  form (R8); AP-V2-35 partner (Vol II PTVV table); V2-AP138 (Vol II
  $d$-dependent shift table).

- **AP953 ($H^3$ class-dependent: 2A $\in \Z/2$, 2B $\in \Z/4$; pentagon
  single-colour coherence at $24 I_1$ nodes is class-dependent).** Vol I's
  pentagon single-colour coherence computation at the 24 $I_1$ nodes (on the
  class-$\mathcal S$ Gaiotto curve $\Sigma_{0, 24}$) gives an $H^3$-cocycle
  whose value depends on the Mathieu conjugacy class: $H^3 \in \Z/2$ for
  class 2A (fixing two points of $\Omega_{24}$), $H^3 \in \Z/4$ for class 2B
  (fixing one point and one 2-cycle). Treating $H^3$ as class-independent
  (e.g., uniformly $\Z/2$ or uniformly $\Z/4$) loses the Mathieu-class
  refinement. **Wrong claim**: "pentagon $H^3$-coherence at 24 $I_1$ nodes
  is uniformly $\Z/2$" or "$H^3$ takes a single value across Mathieu classes".
  **Precise error**: conjugacy-class-dependent torsion presented as
  conjugacy-class-independent; the $H^3$ value is a Mathieu-class invariant,
  not a universal invariant. **Correction**: every Vol I $H^3$ pentagon
  coherence statement names the Mathieu conjugacy class and its torsion
  value: 2A $\in \Z/2$, 2B $\in \Z/4$, and analogous refinements for other
  classes in the 25-class table. **Primary**: Vol III Wave 18 Mathieu
  class-dependent $H^3$ computation; Conway--Sloane SPLAG Ch 10;
  Eguchi--Ooguri--Tachikawa 2010 (Mathieu moonshine for K3).
  **Detection heuristic**: grep `'H\\^3'` within 10 lines of "$I_1$ node" or
  "pentagon coherence" without a Mathieu class specifier triggers Gate 0.
  **Cross-ref**: Vol III AP-CY for primary form (R9); canonical preamble
  row 39 (Rank-162 MTC); AP950 (Gaiotto curve); V2-AP142 (Vol II class
  $\mathbf M$ $E_3$-bar partner).

### Vol III wave-12 cross-programme propagation (2026-04-22, twelve-item anchor block A--L, disjoint from AP929--AP953 by scope)

Vol III wave-12 surfaced twelve cross-programme antipatterns (lettered A--L in the Vol III session ledger) that merit dedicated Vol~I anchors over and above the AP929--AP953 two-stage-factorisation + associated wave. Entries AP954--AP965 below pin each letter A--L to a Vol~I-scoped AP with its own trigger and correction; existing AP929--AP953 entries back-anchor where overlap exists, and the new AP954--AP965 entries carry scope (chain-level vs $(\infty,1)$, $E_n$-on-which-object, PTVV-shift discipline, Dunn--Lurie real-vs-complex, $(B_{E_1}, B_{E_2}, B_{E_\infty})$ frontier) that the two-stage wave did not itself resolve.

- **AP954 (item A -- Pattern 273: chain-level statement vs $(\infty,1)$-categorical functor claim; cross-programme).** Pattern 273 ($\Phi$ functor vs object-level correspondence) carries across all three volumes as a *scope declaration, not a hierarchy*: the chain-level object-level $\Phi$ and the $(\infty,1)$-categorical $\Phi$-as-functor (when morphism preservation is proved) are two different statements about two different categorical structures, both load-bearing, both documented at their precise scope. Vol~I locus: any Vol~I theorem citing Vol~III's $\Phi_d$ output must specify whether the citation reads the chain-level object-level correspondence (Stage-2 chain-level image on the reference curve $C$) or the $(\infty,1)$-categorical functor lift (AP-V2-41 / V2-AP144 frontier conjecture in Vol~II). **Wrong claim**: monolithic ``$\Phi$ functor'' framing without naming the lane. **Precise error**: chain-level proofs auditable against explicit complexes and primary citations vs $(\infty,1)$-categorical claims auditable against Lurie HA + Francis 2013 occupy different epistemic registers; conflating them costs Beilinson-dictum rigour at the audit site. **Correction**: Vol~I citations of Vol~III's $\Phi$ output carry a lane tag (``chain-level on $C$'' or ``$(\infty,1)$-functor, conjectural''); AP954 fires at any unscoped ``$\Phi$ is a functor'' in reader-facing Vol~I prose. **Related**: Vol~I Pattern 273; Vol~II AP-V2-30 / V2-AP133; Vol~III AP-CY Pattern 273; AP929 (two-stage factorisation pins the Stage-2 chain-level image); Vol~III wave-12 source \texttt{notes/wave12\_a1\_bar\_cobar\_gelfand.tex} / \texttt{notes/wave12\_a1\_phi\_functor\_foundations.tex}.

- **AP955 (item B -- bare $\kappa$ forbidden: four subscripted invariants).** Vol~I landscape census and cross-volume $\kappa + \kappa^! = 8$ B-row ceiling statements conflate bare $\kappa$ with one of four genuinely distinct invariants: $\kappa_{\mathrm{ch}}$ (chiral-side Hodge supertrace via $\Phi$, compact CY$_d$), $\kappa_{\mathrm{cat}} = \chi(\mathcal O_X)$ (Künneth-multiplicative on products; $\kappa_{\mathrm{cat}}(K3 \times E) = 0$ on total space, not $2$ which is $\kappa_{\mathrm{fiber}}(K3)$), $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$ (Borcherds weight), $\kappa_{\mathrm{fiber}}$ (fibre/lattice correction). Each satisfies a different functoriality; each enters a different face of $r_{\mathrm{CY}}$. **Wrong claim**: unsubscripted $\kappa$ in Vol~I prose or theorem statement at any CY-to-chiral interface. **Precise error**: AP5 dual-indexing adds a second layer ($\Phi_{10}$ vs $\Phi_{12}$) that bare $\kappa$ cannot carry. **Correction**: every $\kappa$ in Vol~I reader-facing prose and equations carries its subscript from first keystroke; AP955 fires on any unscoped $\kappa$ in the CY-to-chiral chapter family. **Related**: AP5 canonical preamble; AP113 (bare-$\kappa$ Vol~III); AP930 and AP933 (Vol~I wave-12 specific $\kappa_{\mathrm{cat}}$ anchor); Vol~II AP-V2-31 / V2-AP134; Vol~III HZ-7 subscript discipline; Vol~III wave-12 source \texttt{notes/wave12\_a2\_kappa\_invariants\_universal\_borcherds.tex}.

- **AP956 (item C -- CoHA $= Y^+$ evaluation chain; not $\mathcal W_{1+\infty}$).** Vol~I's Yangian / W-algebra cross-chapters cite the cohomological Hall algebra image under Schiffmann--Vasserot as ``$\mathrm{CoHA}(\mathbb C^3) = \mathcal W_{1 + \infty}$''; the correct identification is $\mathrm{CoHA}(\mathbb C^3) = Y^+$ (positive half of the affine Yangian), with the $\mathcal W_{1 + \infty}$ slice arising only through the evaluation arrow $Y^+ \hookrightarrow Y \xrightarrow{\mathrm{ev}_\lambda} \mathrm{End}(\mathcal W_{1 + \infty}[\lambda])$ specialising at a spectral parameter $\lambda$. **Wrong claim**: CoHA $= \mathcal W_{1 + \infty}$ as an isomorphism, collapsing the evaluation chain. **Precise error**: collapses (i) the positive-half identification, (ii) Miki $S_3$-triality structure on $Y^+$, (iii) spectral parameter dependence through $\mathrm{ev}_\lambda$, into a single-line equality. **Correction**: every Vol~I citation of the CoHA side writes the evaluation chain $Y^+ \hookrightarrow Y \xrightarrow{\mathrm{ev}_\lambda} \mathrm{End}(\mathcal W_{1 + \infty}[\lambda])$ and names the spectral parameter. **Related**: AP934 (Vol~I wave-12 CoHA back-anchor); Vol~II AP-V2-32 / V2-AP135; Vol~III AP-CY CoHA-vs-vertex-algebra family; Vol~III wave-12 source \texttt{notes/wave12\_a11\_coha\_y\_plus\_vs\_w\_infty.tex}; Schiffmann--Vasserot 2012 \emph{Publ IHES} 115; Miki 2007 \emph{J Math Phys} 48.

- **AP957 (item D -- $E_n$ on wrong object: at $d \ge 3$, $A$ is $E_1$; $E_2$ lives on $\mathcal Z(\mathrm{Rep}(A))$).** Vol~I theorems referencing Vol~III's $\Phi$ output at $d \ge 3$ sometimes state ``$A$ is $E_2$-chiral'' where $A$ is the algebra output of Stage-2 specialisation; the correct statement is that the $E_2$-structure lives on the Drinfeld centre $\mathcal Z(\mathrm{Rep}(A))$, NOT on $A$ itself, and $A$ is $E_1$ at $d \ge 3$. **Wrong claim**: ``$A$ is an $E_2$-chiral algebra at $d \ge 3$''. **Precise error**: conflates $A$ with $\mathcal Z(\mathrm{Rep}(A))$; conflates braiding-on-module-centre with braiding-on-algebra. $A$ at $d \ge 3$ carries at most $E_1$-structure; the $E_2$-braided structure that Vol~III's $\Phi_d$ lands in at $d \ge 3$ sits on the module-category centre, a different categorical object with its own derived complex. **Correction**: every Vol~I statement naming ``$E_2$'' or ``$E_n$ for $n \ge 2$'' at $d \ge 3$ ambient specifies whether the structure sits on $A$ or on $\mathcal Z(\mathrm{Rep}(A))$; any $E_2$-on-$A$ claim at $d \ge 3$ is flagged. **Related**: Vol~II AP-V2-34 / V2-AP137 ($E_n$ on wrong object); Vol~II AP-V2-39 / V2-AP142 ($E_2$-braiding on $Z(\mathrm{Rep}(A))$); AP276 ($E_n$-chiral notion unspecified); AP277 ($\mathsf{SC}^{\mathrm{ch,top}}$ bicoloured vs $E_3$); Vol~III wave-12 $E_n$-output-scope family; AP-V2-25 MIXED-operation holomorphic pushforward.

- **AP958 (item E -- $\Phi$ output $d$-dependent: PTVV shift table).** The PTVV shifted-symplectic structure controls the $E_n^{\mathrm{cl}}$ classical output of $\Phi_d$ via the dimensional table $(d, \mathrm{shift}, E_n^{\mathrm{cl}}) \in \{(2, -2, E_2), (3, -1, E_1), (4, 0, E_0), (5, +1, E_5\text{-Poisson})\}$. Vol~I theorems citing Vol~III's output routinely pick one row and silently universalise (e.g., assuming $E_2$-chiral for every $d$); the correct reading is $d$-dependent, each row indexing a different Vol~I ambient. **Wrong claim**: Vol~III's $\Phi_d$ output universally $E_2$-chiral or universally $E_1$-chiral without naming $d$. **Precise error**: silently fixes the PTVV shift to a single value, collapsing the $d$-dependent hierarchy that Vol~III Part III documents. **Correction**: every Vol~I reference to Vol~III's $\Phi_d$ output at CY$_d$ names the triple $(d, \mathrm{shift}, E_n^{\mathrm{cl}})$ from the PTVV table; AP958 fires on any monolithic ``$\Phi$ output is $E_n$'' without a $d$-index. **Related**: AP929 (two-stage factorisation pins the Stage-2 curve ambient); AP952 (PTVV extends to $d = 5$ Poisson-$E_5$); AP-CY46 ($\Phi_d$-ambient concentration enlargement); Vol~II AP-V2-35 / V2-AP138 PTVV shift table; Vol~III wave-12 source \texttt{notes/wave12\_a5\_ptvv\_shift\_en\_table.tex}; PTVV 2013 \emph{Publ IHES} 117.

- **AP959 (item F -- single-citation canonicity forbidden: Costello--Li + Kontsevich--Tamarkin is a conjunction).** Vol~I citations of ``the canonical chiral algebra on $X$'' sometimes cite Costello--Li 2016 holomorphic locality alone, or Kontsevich--Tamarkin $E_d$-formality alone, as sufficient input. Both are load-bearing and both enter the canonical chiral algebra construction: holomorphic locality without formality produces a chain-level cocycle without homotopy control; formality without holomorphic locality produces a rigidity statement without the specifically-holomorphic factorisation structure. **Wrong claim**: single-source citation for canonicity of Vol~III's chiral algebra output. **Precise error**: one-leg citation collapses the conjunction. **Correction**: every Vol~I citation of ``canonical chiral algebra'' names both Costello--Li 2016 (arXiv:1605.09439) AND Kontsevich--Tamarkin $E_d$-formality (Kontsevich 2003, Tamarkin 2003, both \emph{Lett Math Phys} 66). **Related**: AP934 (Vol~I wave-12 canonical back-anchor); Vol~II AP-V2-33 / V2-AP136; Vol~III wave-12 source canonical-conjunction notes.

- **AP960 (item G -- ``X gives Y'' narration without explicit arrow).** Vol~I prose contains narrative phrases like ``the shadow tower gives the Gaudin generators'', ``the bar construction gives the Koszul dual'', ``$\Phi$ gives a chiral algebra''; these are narration, not a chain-level or functorial operation. The correct form names the explicit arrow: a chain homotopy, a restriction map, a specialisation functor, a specific operadic inclusion, or an explicit bimodule transport. **Wrong claim**: a named construction ``gives'' another object without naming the arrow. **Precise error**: the verb ``gives'' is verb-of-all-trades that hides whether the operation is a chain-level map, a functor, a restriction, an inclusion, or a specialisation. **Correction**: every ``gives'' following a named construction (shadow tower, bar construction, Koszul duality, $\Phi$-assignment) is replaced by the specific arrow or flagged. **Related**: AP935 (Vol~I wave-12 narration anchor); AP278 (holography as narrative without operadic map); AP272 (narration as mechanism without chain-level map); Vol~III wave-12 source narration-arrow notes.

- **AP961 (item H -- denominator $=$ bar Euler anchor: requires BOTH CY-A AND Vol~I Theorem~A).** Vol~I claims of the form ``the Borcherds-product denominator equals the bar Euler characteristic'' at a landmark K3-type CY threefold require two inputs: (a) a CY-A$_d$ existence or equivalence theorem at the CY side, (b) Vol~I Theorem~A (bar--cobar adjunction) at the chiral side. Either alone is insufficient: CY-A without Theorem~A produces a Hodge-supertrace identity without the bar-Euler interpretation; Theorem~A without CY-A produces an Euler identity without a CY-side origin. **Wrong claim**: denominator-$=$-bar-Euler with a single-sided citation. **Precise error**: single-sided citation collapses the two-anchor structure. **Correction**: every Vol~I denominator-$=$-bar-Euler citation names both CY-A$_d$ (at the specific $d$) and Vol~I Theorem~A within 20 lines. **Related**: AP936 (Vol~I wave-12 denominator anchor); AP938 (CY-A$_3$ existence-and-rigidity, not equivalence); Vol~III wave-12 source denominator/bar-Euler two-anchor notes; Borcherds 1995 \emph{Invent Math} 120.

- **AP962 (item I -- class $\mathsf M$ $E_3$-bar $= 6^g$ at cohomology, not infinite at chain level).** Class $\mathsf M$ chiral algebras (Virasoro and $\mathcal W_N$-type) have infinite-dimensional $E_3$-bar chain complex at generic level; the finite $6^g$ genus-$g$ statement lives at cohomology, with per-handle tricomplex $P_{E_4}(t) = 3t + 3t^2$ (Vol~II AP-V2-38 / V2-AP141) and cohomology collapse via a $d_4$ differential. Vol~I narrations that call class $\mathsf M$ ``$E_3$-bar finite'' or ``$6^g$-dimensional'' without the ``at cohomology'' qualifier collide with chain-level infinity. **Wrong claim**: ``class $\mathsf M$ $E_3$-bar is $6^g$-dimensional'' without ``at cohomology''. **Precise error**: scope confusion between chain-level complex (infinite) and cohomology (finite $6^g$). **Correction**: every class-$\mathsf M$ $6^g$ claim carries the ``at cohomology'' qualifier; AP962 fires on any $6^g$ in class-$\mathsf M$ context without the qualifier. **Related**: AP937 (Vol~I wave-12 $6^g$ anchor); AP-CY21 / AP-CY38 (closed form via Künneth); Vol~II AP-V2-38 / V2-AP141 (per-handle tricomplex); AP-CY46 (cohomological concentration); Vol~III wave-12 source \texttt{notes/wave12\_a10\_class\_m\_e3\_bar.tex}.

- **AP963 (item J -- two-stage factorisation: $\Phi_d = \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C} \circ \Phi^{\mathrm{FA}}_d$; CY-A$_3$ existence-and-rigidity, NOT equivalence).** Vol~I theorems citing Vol~III's CY-A$_3$ result as a full equivalence between the CY$_3$ category and the chiral-algebra output collapse two distinct Vol~III statements: (i) $\Phi_d = \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C} \circ \Phi^{\mathrm{FA}}_d$ factorisation is well-defined at $d = 3$ (existence); (ii) the Stage-1 output is $E_1$-rigid (rigidity); the *equivalence* of Vol~III's $\Phi_3$ with any Vol~I chiral algebra lifts to a Vol~II conjectural $(\infty,1)$-functor, not established at AP963 inscription time. **Wrong claim**: CY-A$_3$ as an equivalence in Vol~I prose. **Precise error**: conflates existence-and-$E_1$-rigidity (proved; chain-level on an explicit factorisation-algebra complex) with $(\infty,1)$-categorical equivalence (frontier; Vol~II AP-V2-41). **Correction**: every Vol~I citation of CY-A$_3$ reads it as existence-plus-rigidity at the factorisation-algebra level, not as an equivalence; AP963 fires on any ``CY-A$_3$ equivalence'' phrasing. **Related**: AP929 (two-stage factorisation framing); AP938 (CY-A$_3$ existence-and-rigidity anchor); AP951 (two-stage factorisation on K3-Yangian side); Vol~II AP-V2-41 / V2-AP144 (frontier conjectural lift); Vol~III wave-12 source \texttt{notes/wave12\_a3\_cy\_a3\_equivalence.tex}; Costello--Gwilliam \emph{FA} Vol 2; Kontsevich--Soibelman 2011 arXiv:1006.2706.

- **AP964 (item K -- Dunn--Lurie real-vs-complex slip; Vol~II primary, Vol~I cross-ref).** Dunn--Lurie additivity $\int_{\Sigma_2} E_3 \simeq E_1$ is stated in real dimensions (2-real-dimensional $\Sigma_2$); conflation with complex-dimensional conventions (Vol~III's $\Sigma_{d-1}$ at complex dimension) can swap $\Sigma_2$-over-reals with $\Sigma_1$-over-complex, producing shift errors in the Stage-2 specialisation. Vol~II-primary (AP-V2-36 / V2-AP139); Vol~I cross-reference here so any Vol~I citation of Dunn--Lurie additivity in the CY-to-chiral factorisation chain carries the real-vs-complex discipline. **Wrong claim**: $\int_{\Sigma_2} E_3 \simeq E_1$ without naming whether $\Sigma_2$ is 2-real-dimensional (Dunn--Lurie) or $1$-complex-dimensional (a curve, as in Vol~III's Stage-2). **Precise error**: silent dimension-convention switch at the Dunn--Lurie input. **Correction**: every Vol~I citation of Dunn--Lurie additivity names ``2-real-dimensional $\Sigma_2$''; for Vol~III-style complex-dimension pushforwards, name $\Sigma_{d-1}$ over complex dimension and cite the MIXED-operation holomorphic pushforward (AP-V2-25). **Related**: Vol~II AP-V2-36 / V2-AP139 (primary); AP-V2-25 closed-colour identification with stage-1; AP277 ($\mathsf{SC}^{\mathrm{ch,top}}$ bicoloured); Vol~III wave-12 source Dunn--Lurie-real-vs-complex notes; Dunn 1988 \emph{JPAA} 50; Lurie \emph{HA} \S 5.1.

- **AP965 (item L -- $(B_{E_1}, B_{E_2}, B_{E_\infty})$ chain-level filtered equivalence frontier; Vol~II primary, Vol~I cross-ref).** Vol~II's three-sector bar complex $(B_{E_1}, B_{E_2}, B_{E_\infty})$ carries a chain-level filtered quasi-isomorphism $B_{E_1} A \twoheadrightarrow B_{E_2} A \twoheadrightarrow B_{E_\infty} A$ with explicit Kontsevich--Tamarkin $E_d$-formality and Costello--Li holomorphic-locality inputs; the $(\infty,1)$-categorical-functor lift is conjectural and named (Vol~II AP-V2-41 / V2-AP144). Vol~I theorems citing the three-sector bar complex as an equivalence of $(\infty,1)$-categories collapse the chain-level vs functor-lift scope. **Wrong claim**: $(B_{E_1}, B_{E_2}, B_{E_\infty})$ as an $(\infty,1)$-equivalence in Vol~I prose. **Precise error**: conflates the chain-level filtered equivalence (theorem) with the functor lift (frontier). **Correction**: every Vol~I citation of the three-sector bar complex reads it as a chain-level filtered equivalence with explicit Kontsevich--Tamarkin + Costello--Li inputs; any $(\infty,1)$-functor claim is flagged as conjectural. **Related**: Vol~II AP-V2-41 / V2-AP144 (primary, frontier); Vol~II AP-V2-30 / V2-AP133 (chain-level vs $(\infty,1)$-functor); AP954 (item A, cross-programme Pattern 273 hook); Vol~III wave-12 source \texttt{notes/wave12\_a9\_be1\_be2\_beinf\_filtered.tex}.

**Rule** (cross-volume discipline, orthogonal to AP929--AP953):
AP929--AP953 cover the Vol~III wave-12 two-stage-factorisation
propagation at the specific subjects (single-stage $\Phi_d$ framing,
six routes, Fake-Monster $d$-stratification, K3 Yangian envelope,
$\kappa_{\mathrm{cat}}(K3 \times E)$, CoHA $\neq \mathcal W_{1+\infty}$,
\texttt{warning}-environment infiltration, discipline-token leakage,
denominator-$=$-bar-Euler, class $\mathsf M$ $6^g$, CY-A$_3$
equivalence, etc.); AP954--AP965 cover the twelve-letter A--L
cross-programme anchor block as Vol~I-scoped entries carrying
lane-discipline, $E_n$-on-which-object, PTVV shift, Dunn--Lurie
real-vs-complex, and $(B_{E_1}, B_{E_2}, B_{E_\infty})$ frontier
scopes that AP929--AP953 did not themselves resolve. Each of
AP954--AP965 cross-references the AP929--AP953 entry of matching
subject where overlap exists; the two blocks coexist without
subsumption and jointly close the Vol~III wave-12 (2026-04-22)
cross-programme propagation at Vol~I.
**Related**: AP-CAT-15 (wave-12 cross-volume import anchor);
AP929--AP953 (two-stage factorisation + associated wave); Vol~II
AP-V2-30 through AP-V2-41 (Vol~II wave-12 mirror); AP-CAT-14
(wave-14, distinct); AP-CAT-12 (wave-5 residual); canonical
preamble $\kappa_\bullet$-subscript rows.

## 6d hCS audit + Harmonies synthesis errors (2026-04-22, AP966--AP983)

The eighteen entries below (AP966--AP983, originally catalogued as
AP-$\alpha$ through AP-$\tau$) crystallise errors encountered in the
2026-04-22 6d holomorphic Chern--Simons audit and the Harmonies
synthesis pass. Ten of these (AP966 = $\alpha$, AP967 = $\beta$,
AP971 = $\zeta$, AP972 = $\eta$, AP973 = $\theta$, AP975 = $\kappa$,
AP976 = $\lambda$, AP977 = $\mu$, AP981 = $\rho$, AP983 = $\tau$)
are cross-volume relevant and are mirrored as single-line pointers in
`notes/cross_volume_aps.md`. Each entry follows a three-part
first-principles protocol: (a) RIGHT / ghost of a true theorem;
(b) WRONG / precise error; (c) CORRECT / counter + primary citations.

- **AP966 (= AP-$\alpha$) -- Casimir quadratic-for-quartic in 6d hCS
  1-loop: $c_2(\fg)/24$ bubble with 3-factor integrand is wrong on BOTH
  Casimir degree AND integrand factor count (2026-04-22 6d hCS audit).**
  The 1-loop BV obstruction for 6d holomorphic Chern--Simons on a
  Calabi--Yau 3-fold is a nontrivial $H^1_{\mathrm{BRST}}$ class; its
  precise form is fixed by gauge-anomaly cancellation in BV.
  **Ghost of a true theorem**: the 1-loop anomaly is controlled by an
  $\mathrm{ad}$-invariant polynomial on $\fg$ entering at a specific
  Casimir degree; Costello 2017 \emph{Japan J. Math.}~12
  ``M5-brane and the AdS/CFT correspondence'' eq.~(5.17),
  Costello--Li 2020 arXiv:2005.11356 eq.~(6.12), and the 2026 CFG
  audit all lock the form.
  **Wrong claim**: 1-loop anomaly $= c_2(\fg)/24 \cdot
  \int A \wedge \bar\partial A \wedge \bar\partial A$ (quadratic
  Casimir, three $A$-factors).
  **Precise error**: (a) Casimir degree --- the 6d hCS 1-loop
  obstruction is governed by the adjoint \emph{quartic} Casimir
  $\mathrm{tr}_{\mathrm{adj}}(T^{(a} T^b T^c T^{d)})_{\mathrm{sym}}$,
  not the quadratic Casimir
  $c_2(\fg) = \mathrm{tr}_{\mathrm{adj}}(T^a T^a)$;
  (b) integrand factor count --- the canonical representative is
  $c \wedge (\bar\partial c)^{\wedge 4}/(2\pi i)^3$ with four
  $\bar\partial c$-factors (ghost number 1, form degree $0 + 4$), OR
  the cohomologous bubble
  $c_2/24 \cdot A \wedge \bar\partial A \wedge \bar\partial A$ when
  written on the connection side after CS descent
  $F = \bar\partial c + \tfrac12 [c, c]$.
  **Counter**: every 6d hCS 1-loop statement names the Casimir degree
  (quartic, not quadratic), the ghost-side representative
  $c \wedge (\bar\partial c)^{\wedge 4}$, and the CS-descent
  equivalence to the $A$-side bubble; the three-$A$-factor bubble
  with $c_2/24$ is valid only as the CS-descended $A$-side form of
  the correct ghost-side obstruction, not as an independent
  quadratic-Casimir claim.
  **Related**: AP966; AP923 (hCS vs categorical $\mathbb E_3$
  grammar); AP934 (CL + KT input canonicity); AP979 (Deligne-safe
  anomaly list); Costello 2017 \emph{Japan J. Math.}~12 eq.~(5.17);
  Costello--Li 2020 eq.~(6.12); CFG 2026 audit memo.

- **AP967 (= AP-$\beta$) -- Weak Jacobi form polar coefficient:
  $c_{\phi_{-2,1}}(-1) = 24$ is FALSE. $\phi_{-2,1}$ holomorphic at
  cusp forces $c_{\phi_{-2,1}}(-n) = 0$ for $n \ge 1$; 24 enters via
  $\chi(K3)$ or $1/\Delta_5$ (2026-04-22 Harmonies synthesis).**
  The integer 24 enters 6d hCS 1-loop generating functions via the K3
  elliptic-genus structure, not via a polar coefficient of $\phi_{-2,1}$
  itself.
  **Ghost of a true theorem**: the 6d hCS BV-loop generating function
  on K3-fibred CY$_3$ carries a $\chi(K3) = 24$ normalisation,
  ultimately traceable to the Eichler--Zagier weak-Jacobi index-1
  sector.
  **Wrong claim**: $c_{\phi_{-2, 1}}(-1) = 24$; reading 24 as a polar
  Fourier coefficient of $\phi_{-2, 1}(\tau, z)$.
  **Precise error**: $\phi_{-2, 1}$ is a \emph{weak} Jacobi form of
  weight $-2$, index $1$ (Eichler--Zagier 1985 Ch.~9): the
  $q$-expansion at the cusp is holomorphic, so the polar part vanishes
  and $c_{\phi_{-2, 1}}(-n) = 0$ for every $n \ge 1$. The true
  identifications: $\chi(K3) = 24 = |c_{1/\Delta_5}(q^1)|$
  (Gritsenko--Nikulin 1998 \emph{Duke}~95 denominator Fourier
  expansion); $\chi(K3) = 24 = \phi_{K3}(\tau, 0)$ when the K3
  elliptic genus $\phi_{K3}$ is evaluated at $z = 0$
  (Eguchi--Ooguri--Tachikawa 2011 arXiv:1004.0956). Neither uses a
  polar coefficient of bare $\phi_{-2, 1}$.
  **Counter**: every occurrence of 24 in a 6d hCS BV generating
  function names its source ($\chi(K3)$, $c_{1/\Delta_5}(q^1)$, or
  $\phi_{K3}(\tau, 0)$), not $c_{\phi_{-2, 1}}$.
  **Related**: AP967; AP969 (wedge-of-1-forms ill-typed); AP972
  ($\eta^{24}$ denominator bimodular weight); AP-CY42 ($\phi_{0, 1}$
  normalisation); Eichler--Zagier 1985 \emph{Theory of Jacobi Forms}
  Ch.~9; Gritsenko--Nikulin 1998; Eguchi--Ooguri--Tachikawa 2011.

- **AP968 (= AP-$\gamma$) -- Feynman graph Euler miscount: theta graph
  is 2-loop, not 1-loop (2026-04-22 6d hCS audit).**
  Loop-number counting in BV Feynman diagrammatics uses the Euler
  identity $V - E + L = 1$ for connected graphs.
  **Ghost of a true theorem**: the characteristic 1-loop 2-leg graph
  in $\phi^3$-theory perturbation controls the 1-loop self-energy /
  anomaly computations.
  **Wrong claim**: ``the theta graph ($V = 2$, $E = 3$) is one-loop''.
  **Precise error**: Euler identity gives
  $L = 1 - V + E = 1 - 2 + 3 = 2$, so the theta graph is
  \emph{two-loop}, not one-loop. The correct 1-loop 2-leg graph in
  $\phi^3$-theory is the \emph{bubble} ($V = 2$, $E = 2$, $L = 1$).
  **Counter**: every loop-order assertion in BV Feynman diagrammatics
  is verified against $L = 1 - V + E$; use \emph{bubble} for 1-loop
  2-leg, \emph{theta} for 2-loop 2-leg (three-edge), \emph{sunset} for
  the 2-loop 2-leg three-edge variant.
  **Related**: AP968; AP966 (6d hCS 1-loop form); Costello 2011
  \emph{Renormalization and Effective Field Theory} Ch.~2 graph
  combinatorics; standard QFT (Peskin--Schroeder Ch.~10).

- **AP969 (= AP-$\delta$) -- Stuffle dimensional failure:
  $\zeta(3)^2 = 3\zeta(3, 3) + 6\zeta(5)/5$ mixes weights 5 and 6
  (2026-04-22 Harmonies synthesis).**
  The stuffle (harmonic shuffle) product on multiple zeta values
  $\zeta(n_1, \ldots, n_r)$ preserves total weight $\sum n_i$; LHS and
  RHS must have equal weight.
  **Ghost of a true theorem**: the square $\zeta(3)^2$ admits a
  depth-reduction in the motivic MZV ring (Brown 2012
  \emph{Annals}~175).
  **Wrong claim**: $\zeta(3)^2 = 3 \zeta(3, 3) + 6 \zeta(5)/5$.
  **Precise error**: (a) weight mismatch --- LHS has weight
  $3 + 3 = 6$; RHS mixes weight-6 ($\zeta(3, 3)$) and weight-5
  ($\zeta(5)$) terms, which is dimensionally forbidden in the MZV
  grading; (b) coefficient --- the correct stuffle gives
  $\zeta(3)^2 = 2 \zeta(3, 3) + \zeta(6)$ (prefactor 2, not 3; no
  $\zeta(5)$).
  **Counter**: every MZV product identity is weight-homogeneous by
  construction; before asserting a depth-reduction, verify total
  weight and apply the correct stuffle
  $\zeta(a) \cdot \zeta(b) = \zeta(a, b) + \zeta(b, a) + \zeta(a + b)$
  with $a = b = 3$: $\zeta(3)^2 = 2 \zeta(3, 3) + \zeta(6)$.
  **Related**: AP969; AP970 (wedge-of-1-forms ill-typed); Brown 2012
  \emph{Annals Math.}~175 ``Mixed Tate motives over $\Z$''; Zagier
  values at $\zeta(3)$.

- **AP970 (= AP-$\epsilon$) -- Wedge-of-1-forms ill-typed:
  $\eta_{12}^{\wedge 3} = 0$ by anticommutativity of 1-forms
  (2026-04-22 6d hCS audit).**
  In Goncharov's iterated-integral / period formalism for motivic
  multiple zeta values, polylogarithms arise from iterated integrals
  over \emph{distinct} logarithmic 1-forms on configuration spaces.
  **Ghost of a true theorem**: iterated-integral formulas realise
  motivic MZVs as periods on $\mathrm{Conf}_n(\C)$ (Goncharov 2001
  arXiv:math/0103059).
  **Wrong claim**: $\int_{z_2} \eta_{12}^{\wedge 3} = \zeta(3)$ with
  $\eta_{12} = d\log(z_1 - z_2)$.
  **Precise error**: $\eta_{12}$ is a 1-form, hence by
  anticommutativity of the wedge on 1-forms,
  $\eta_{12} \wedge \eta_{12} = 0$, and thus $\eta_{12}^{\wedge 3} = 0$
  identically. The expression is ill-typed; no integral of it yields
  $\zeta(3)$.
  **Counter**: the correct Goncharov iterated integral realising
  $\zeta(3)$ uses three \emph{distinct} logarithmic 1-forms:
  $\int_{\Delta^2} \eta_{12} \wedge \eta_{23} \wedge \eta_{31}
  = - \zeta(3)$ on a 2-simplex chain. Every iterated-integral MZV
  expression carries distinct index pairs $(i_k, j_k)$ at each wedge
  slot.
  **Related**: AP970; AP969 (stuffle dimensional failure); AP975
  (K3-MHS pollution); Goncharov 2001 arXiv:math/0103059; Brown 2012
  \emph{Annals}~175.

- **AP971 (= AP-$\zeta$) -- Residue ill-typed:
  $\mathrm{Res}_{\Phi_{10} = 0}(c/\Phi_{10}) = c$, not $c/\Phi_{10}$
  (2026-04-22 Harmonies synthesis).**
  The residue of a meromorphic expression $c/f$ at the divisor
  $\{f = 0\}$ in the simple-pole case is $c$ restricted to the
  divisor, not $c/f$ (which is the original expression,
  self-referential).
  **Ghost of a true theorem**: the 3d gravity partition function on
  $\mathrm{AdS}_3 \times K3$ has pole structure on the Humbert /
  $\Phi_{10}$-divisor at genus 2 (Dijkgraaf--Verlinde--Verlinde 1997
  arXiv:hep-th/9703030).
  **Wrong claim**: $\mathrm{Res}_{\Phi_{10} = 0}(Z_{\mathrm{hCS}}) =
  1/\Phi_{10}$.
  **Precise error**: the residue \emph{removes} one power of the
  vanishing factor; it cannot be proportional to the expression being
  differentiated. The correct identification is DIRECT:
  $Z^{\mathrm{AdS}_3 \times K3}_{\mathrm{3dQG}} = 1/\Phi_{10}$ itself
  (DVV 1997); the hCS effective action arises via the Borcherds
  singular-theta lift of a weight-10 Siegel modular form, not as a
  residue of itself.
  **Counter**: every pole-structure claim on $\Phi_{10}$-divisor
  distinguishes (i) the partition function $1/\Phi_{10}$ (direct, DVV
  1997); (ii) the effective action via Borcherds theta lift;
  (iii) residues at other Humbert / Heegner divisors;
  self-referential residues are forbidden.
  **Related**: AP971; AP972 ($\eta^{24}$ denominator bimodular
  weight); AP981 ($\kappa_{\mathrm{BKM}}$ 5 vs 12 scope); DVV 1997
  arXiv:hep-th/9703030; Borcherds 1995 \emph{Invent.\ Math.}~120.

- **AP972 (= AP-$\eta$) -- $\eta^{24}$ denominator incomplete:
  $\Phi_{10}|_{z = 0}$ has bimodular weight $(12, 12)$ so single
  $\eta^{24}$ leaves weight $(-12, 0)$; BOTH $\eta(\tau)^{24}$ AND
  $\eta(\tau')^{24}$ required (2026-04-22 Harmonies synthesis).**
  The Siegel paramodular form $\Phi_{10}$ of weight 10, when
  restricted to the diagonal $\{z = 0\}$ of the Siegel half-space
  $\mathcal H_2$, factorises along a specific bimodular-weight
  pattern.
  **Ghost of a true theorem**: the weight-zero bimodular object
  extracted from $\Phi_{10}$ on the diagonal, entering the 3d gravity
  / 6d hCS partition function on $\mathrm{AdS}_3 \times K3$, has a
  specific ratio structure.
  **Wrong claim**: $(\Phi_{10}/\eta^{24})^\hbar$ with a single
  $\eta^{24}$ factor.
  **Precise error**: on the diagonal, $\Phi_{10}|_{z = 0}$ transforms
  under $\mathrm{SL}_2(\Z)_\tau \times \mathrm{SL}_2(\Z)_{\tau'}$ with
  bimodular weight $(12, 12)$, not $(10, 0)$ (the weight of
  $\Phi_{10}$ on $\mathcal H_2$ does not transport to a single
  modular factor on either $\tau$ or $\tau'$). Quotient by a single
  $\eta(\tau)^{24}$ (weight 12 in $\tau$) leaves residual weight
  $(-12, 0)$ --- still bimodular-weight nonzero.
  **Counter**: the correct weight-zero bimodular object is
  $\Phi_{10}|_{z = 0}/(\eta(\tau)^{24} \eta(\tau')^{24})$, following
  DVV 1997 factorisation
  $\Phi_{10}|_{z = 0} = z^2 \cdot \eta(\tau)^{24} \eta(\tau')^{24}
  \cdot (\text{Jacobi correction})$.
  Every partition-function denominator on the $(\tau, \tau')$-diagonal
  carries BOTH $\eta(\tau)^{24}$ and $\eta(\tau')^{24}$ factors.
  **Related**: AP972; AP971 (residue ill-typed); AP975 (K3-MHS
  pollution); AP-CY78 (YD tower weight); DVV 1997
  arXiv:hep-th/9703030 eq.~(3.17--3.21); Gritsenko--Nikulin 1998.

- **AP973 (= AP-$\theta$) -- Dunn additivity on holomorphic $\C$:
  naive $E_3 = E_1 \otimes_D E_1 \otimes_D E_1$ on $\C^3$ from three
  holomorphic directions is a TYPE error at the FA level; Dunn holds
  at cohomology (2026-04-22 6d hCS audit).**
  Dunn additivity for $E_n$-algebras, $E_n \simeq E_1^{\otimes n}$ at
  the level of $\infty$-operads (May--Dunn; Lurie \emph{HA} 5.1.2),
  transports to factorisation-algebra (FA) structure only at the level
  of cohomology / after Dolbeault reduction, not at the level of raw
  factorisation algebras on $\C^n$.
  **Ghost of a true theorem**: 6d hCS on $\C^3$ carries an
  $E_3$-structure after passage to $\bar\partial$-cohomology.
  **Wrong claim**: $E_3 = E_1 \otimes_D E_1 \otimes_D E_1$ on $\C^3$
  obtained by tensoring three ``holomorphic'' $E_1$ structures, one
  per holomorphic direction.
  **Precise error**: each $\C$-factor is 2 real dimensions; the
  holomorphic FA on $\C$ is naturally $E_2$-topologically, not $E_1$.
  An $E_1$-shadow per holomorphic direction emerges only after passing
  to $\bar\partial$-cohomology (Dolbeault reduction of the
  Costello--Gwilliam FA). At the FA level on $\C^3$ the structure is
  $E_6$-topologically (six real dimensions), which reduces to $E_3$
  via $\bar\partial$-cohomology; Dunn-style decomposition at the
  cohomology level then gives $E_1$ per holomorphic direction.
  **Counter**: every Dunn decomposition of $E_n$ on holomorphic $\C^k$
  names the lane: (i) FA level on $\C^k$ is $E_{2k}$-topological;
  (ii) $\bar\partial$-cohomology level is $E_k$-complex; (iii) Dunn
  per-direction $E_1$ lives at the cohomology lane. Never state
  $E_n = E_1^{\otimes n}$ at the raw FA level on $\C^n$.
  **Related**: AP973; AP932 ($E_n$ on wrong object at $d \ge 3$);
  AP983 (two-stage $\Phi_d$ naming); Costello--Gwilliam 2017/2021
  Vol II Ch.~5; Lurie \emph{HA} 5.1.2 (Dunn additivity).

- **AP974 (= AP-$\iota$) -- Wheel-to-$\zeta$ convention-sensitive:
  $W_n \leftrightarrow \zeta(2n \pm 1)$ requires convention statement
  (2026-04-22 Harmonies synthesis).**
  The realisation of wheel graphs $W_n$ in the Kontsevich graph
  complex $\mathsf{GC}_2$ as motivic odd-weight zeta values is
  convention-dependent.
  **Ghost of a true theorem**: wheel graphs in $\mathsf{GC}_2$ realise
  odd-weight $\zeta$-values, saturating $\mathrm{grt}_1 \otimes \Q$ at
  low depth (Willwacher 2014 \emph{Invent.\ Math.}~200).
  **Wrong claim**: ``$W_3 \leftrightarrow \zeta(5)$'' without naming
  the convention.
  **Precise error**: two conventions circulate.
  (i) Willwacher 2014: $W_n \leftrightarrow \zeta(2n - 1)$, giving
  $W_3 \leftrightarrow \zeta(5)$.
  (ii) Brown--Schnetz 2012 (motivic period convention on
  $\overline M_{0, n}$): $W_n \leftrightarrow \zeta(2n + 1)$, giving
  $W_3 \leftrightarrow \zeta(7)$.
  Either is correct under its own convention; stating a
  wheel-to-$\zeta$ identity without the convention is ambiguous.
  **Counter**: every $W_n \leftrightarrow \zeta(\cdot)$ identity names
  the convention (Willwacher vs Brown--Schnetz) within 10 lines of the
  assertion.
  **Related**: AP974; AP969 (stuffle); AP975 (K3-MHS pollution);
  Willwacher 2014 \emph{Invent.\ Math.}~200;
  Brown--Schnetz 2012 ``Modular forms in quantum field theory''.

- **AP975 (= AP-$\kappa$) -- K3-MHS pollution of motivic MZV:
  ``Feynman period on $K3 \times \C^2$ is motivic MZV'' without
  K3-separation is a scope error (2026-04-22 Harmonies synthesis).**
  Brown's motivic MZV framework (Brown 2012 \emph{Annals}~175) applies
  to periods on moduli spaces of genus-zero curves with marked points
  $\overline M_{0, n}$; it does NOT a priori extend to Feynman periods
  on K3-fibred backgrounds.
  **Ghost of a true theorem**: Feynman graph periods on Calabi--Yau
  backgrounds live in motivic period rings, generalising the MZV case
  (Brown 2012 Thm 1.2).
  **Wrong claim**: ``the Feynman period on $K3 \times \C^2$ is a
  motivic multiple zeta value'' (bare, without K3-separation).
  **Precise error**: K3 periods form a Mixed Hodge Structure of type
  K3 (weight 2, rank-22 transcendental lattice); they do not belong
  to Brown's motivic MZV Hopf subalgebra. The assertion conflates the
  K3-MHS with the MZV-MHS.
  **Counter**: Feynman periods on $K3 \times \C^2$ are handled in one
  of two ways: (i) push the K3-factor forward via $\chi(K3) = 24$,
  reducing to a $\C^2$ Feynman period \emph{then} invoke motivic MZV
  content; (ii) prove the K3-period sublattice is contained in
  Brown's MZV Hopf algebra (typically requires additional
  Hodge-theoretic input). Brown 2012 Thm 1.2 is scope-restricted to
  $\overline M_{0, n}$.
  **Related**: AP975; AP969 (stuffle); AP970 (wedge ill-typed);
  AP974 (wheel convention); Brown 2012 \emph{Annals}~175 Thm~1.2; K3
  MHS classification.

- **AP976 (= AP-$\lambda$) -- AFT 2017 citation scope:
  stratified-spaces theorem mis-cited for smooth-product Fubini
  (2026-04-22 6d hCS audit).**
  The Ayala--Francis--Tanaka 2017 arXiv:1409.0848 factorisation-homology
  framework is scope-restricted.
  **Ghost of a true theorem**: factorisation homology is pointwise
  in the stratification on stratified spaces (AFT 2017);
  smooth-product Fubini is AF 2015 arXiv:1206.5522 Thm 3.16.
  **Wrong claim**: ``smooth-product Fubini follows from AFT 2017 Thm
  4.1.3''.
  **Precise error**: AFT 2017 Thm 4.1.3 concerns STRATIFIED spaces;
  its stratified-space hypotheses do not transport verbatim to smooth
  products. The smooth-product Fubini for factorisation homology is
  AF 2015 Thm 3.16.
  **Counter**: every factorisation-homology Fubini citation names the
  scope: AF 2015 Thm 3.16 for smooth products
  $\int_{M \times N} A = \int_M \int_N A$; AFT 2017 Thm 4.1.3 for
  stratified spaces.
  **Related**: AP976; AP977 (Costello 2011 non-compact scope);
  AF 2015 arXiv:1206.5522; AFT 2017 arXiv:1409.0848.

- **AP977 (= AP-$\mu$) -- Costello 2011 Thm 13.4.1 compact-support
  hypothesis: RG flow for non-compact $\R^3 \times K3 \times \C^2$
  needs additional input (2026-04-22 6d hCS audit).**
  Costello 2011 \emph{Renormalization and Effective Field Theory} Thm
  13.4.1 guarantees well-defined renormalisation-group flow for BV
  theories under compact-support hypotheses.
  **Ghost of a true theorem**: RG flow is well-defined for BV theories
  on compact (or compactly supported) targets.
  **Wrong claim**: ``Thm 13.4.1 applies directly to
  $\R^3 \times K3 \times \C^2$''.
  **Precise error**: the target is non-compact ($\R^3$ and $\C^2$
  factors). Costello 2011 Thm 13.4.1 requires the space of fields be
  compactly supported OR the target be compact; neither holds naively
  here.
  **Counter**: non-compact RG-flow applications cite either
  (i) Costello--Gwilliam 2017/2021 Vol II Ch.~4 (non-compact targets
  with additional analytic input) or (ii) add a compact-support
  hypothesis (e.g.\ via an infrared cutoff or compact fibre);
  Costello 2011 Thm 13.4.1 is not directly applicable without such
  an upgrade.
  **Related**: AP977; AP976 (AFT 2017 scope); Costello 2011
  \emph{Renormalization and EFT} Ch.~13; Costello--Gwilliam 2017/2021
  Vol II Ch.~4.

- **AP978 (= AP-$\nu$) -- CHSW holonomy embedding: $F_\cA = R$ factors
  through $\mathrm{SU}(4)$-structure group and $\mathrm{SU}(3)$-holonomy
  (2026-04-22 6d hCS audit).**
  The Candelas--Horowitz--Strominger--Witten 1985 \emph{Nucl.\ Phys.\
  B}~258 ``Vacuum configurations for superstrings'' standard embedding
  identifies the gauge connection with the tangent-bundle spin
  connection, producing the heterotic $E_8 \times E_8$ compactification
  on a Calabi--Yau 3-fold.
  **Ghost of a true theorem**: the CHSW standard embedding relates the
  Riemann curvature $R$ to a gauge field strength $F_\cA$ along a
  specific holonomy reduction.
  **Wrong claim**: ``$F_\cA = R$ into $\mathrm{SU}(3)$ holonomy''
  (bare, skipping the $\mathrm{SU}(4)$-structure-group layer).
  **Precise error**: the statement skips the $\mathrm{SU}(4)$
  structure-group layer. The correct sequence is:
  $\mathfrak{so}(6)$-Riemann (real tangent bundle of $\mathrm{CY}_3$)
  $\cong \mathfrak{su}(4)$-structure group (complex structure-bundle
  identification, since $\mathrm{Spin}(6) \cong \mathrm{SU}(4)$),
  reducing to $\mathrm{SU}(3)$-holonomy on $\mathrm{CY}_3$ (trivial
  canonical bundle). Then embedding
  $\mathrm{SU}(3) \hookrightarrow E_8$ identifies $F_\cA$ with $R$ on
  the Calabi--Yau.
  **Counter**: every CHSW embedding statement factors through
  $\mathrm{SU}(4)$-structure group first, then
  $\mathrm{SU}(3)$-holonomy, then $E_8$-embedding.
  **Related**: AP978; AP979 (Deligne-safe anomaly list); CHSW 1985
  \emph{Nucl.\ Phys.\ B}~258; Green--Schwarz--Witten Vol II Ch.~12.

- **AP979 (= AP-$\xi$) -- Deligne-safe 6d hCS anomaly-free list:
  CANONICAL-ANOM-LOCUS (form c) with $E_6$ strict exclusion, $A_2$
  refined/unrefined distinguished (2026-04-22 6d hCS audit).**
  The Deligne exceptional series
  $\{A_1, A_2, G_2, D_4, F_4, E_6, E_7, E_8\}$ carries universal
  tensor identities (Deligne 1996 \emph{CRAS}~322).
  **Ghost of a true theorem**: universal tensor identities on the
  Deligne series underlie a universal quartic-Casimir factorisation
  $\mathrm{tr}_{\mathrm{adj}}T^4=\alpha_{\mathfrak g}(\mathrm{tr}_{\mathrm{adj}}T^2)^2$
  that makes the quartic obstruction uniform across the series.
  **Wrong claims flagged**:
  \begin{itemize}
  \item \textbf{Form (a) strict}: ``6d hCS anomaly-free locus equals
    $\mathrm{Deligne}^{\mathrm{exc}}\setminus\{E_6, A_2\}$ without
    further refinement'' (excludes $A_2$-refined from the locus).
  \item \textbf{Form (b)}: ``6d hCS anomaly-free locus equals
    $\mathrm{Deligne}^{\mathrm{exc}}\setminus\{E_6\}$ alone, allowing
    $A_2$ without distinguishing refined/unrefined'' (wrongly admits
    $A_2$-unrefined with live $d^{abc}$).
  \item \textbf{Original form}: ``the full Deligne series is 6d hCS
    anomaly-free'' (misses $E_6$ cubic Jordan invariant and $A_2$
    Gell-Mann $d$-tensor).
  \end{itemize}
  **Precise error**: each of (a), (b), and the full-Deligne claim
  misrepresents the native-ambient locus. $E_6$ has a nontrivial
  symmetric cubic invariant $\mathrm{Sym}^3(\mathbf{27})$, hence
  $d^{abc}(E_6) \neq 0$; no refinement in the programme's toolkit
  (critical twist, Dimofte slab, Green--Schwarz inflow) cures it
  within native ambient. $A_2 = \mathfrak{su}(3)$ has $d^{abc}(A_2)
  \neq 0$ (Gell-Mann $d$-tensor); unrefined, it sits outside the
  locus, but $A_2$-refined (Feigin--Frenkel critical twist
  $K^{-1/2}$ killing the quadratic obstruction, plus Dimofte-slab
  anomaly-inflow from Vol II Part V providing Green--Schwarz cubic
  cancellation) sits INSIDE the locus.
  **Counter — CANONICAL-ANOM-LOCUS (c)**: the native-ambient
  anomaly-free locus reads
  $$\mathrm{Anom}_1=0\iff \fg\in\bigl(\mathrm{Deligne}^{\mathrm{exc}}\setminus\{E_6, A_2\text{-unrefined}\}\bigr)\cup\{\mathrm{abelian}\}\cup\{\mathrm{str}_{\mathrm{ad}}=0\}\cup\{\widehat\fg_{-h^\vee}\otimes K^{-1/2}\text{-refined}\}.$$
  Native-ambient distinctions:
  \begin{itemize}
  \item $E_6$ STRICTLY excluded (no programme-toolkit refinement
    kills $\mathrm{Sym}^3(\mathbf{27})$ cubic $d^{abc}$; the
    $K^{-1/2}$ twist addresses quadratic, not cubic).
  \item $A_2$-unrefined excluded; $A_2$-refined (critical twist +
    Dimofte slab) INSIDE the locus.
  \item $\{A_1, G_2, D_4, F_4, E_7, E_8\}$ unconditionally inside
    (quartic factorises, cubic vanishes identically).
  \end{itemize}
  Two obstructions: quartic adjoint Casimir (Deligne-killed across
  the series via $\alpha_{\mathfrak g}(\mathrm{tr}_{\mathrm{adj}}T^2)^2$
  factorisation, including for $A_2, E_6$) vs cubic $d^{abc}$
  (nonzero for $A_2, E_6$; cured only in native ambient by
  Green--Schwarz-type inflow, operative for $A_2$-refined, uncurable
  for $E_6$).
  **Regex trigger**: flag BOTH
  (a) \verb|Deligne.*\\setminus.*\\\{E_6.*A_2\\\}(?!.*refined)|
  (strict $\{E_6, A_2\}$ removal without refinement commentary), AND
  (b) \verb|Deligne.*\\setminus.*\\\{E_6\\\}(?!.*A_2)|
  (exclude-only-$E_6$ form allowing undifferentiated $A_2$),
  relative to canonical (c) with the explicit $A_2$-unrefined
  qualifier and $K^{-1/2}$-refinement clause.
  **Related**: AP966 (6d hCS 1-loop form); AP978 (CHSW
  embedding); V2-AP157 / AP-V2-54 (Vol II sibling); AP-CY262
  (Vol III cubic vs quadratic Casimir); AP-CY50-E14 (cross-volume
  ledger); Deligne 1996 \emph{CRAS}~322; Gell-Mann 1962; Feigin--Frenkel 1992 critical-level
  twist; Dimofte 2014 slab anomaly-inflow.

- **AP980 (= AP-$o$) -- Conway $V^{s\natural}$ on super-metaplectic
  $\Psi$-branch, not as 5th bosonic $\Psi$-image (2026-04-22 Harmonies
  synthesis).**
  The $\Psi$-functor on CY-Siegel-automorphic data produces
  Borcherds--Kac--Moody algebras; its $\Psi$-image catalogue is finite
  and structured.
  **Ghost of a true theorem**: Conway $V^{s\natural}$ at central charge
  $c = 12$ is a legitimate $\Psi$-type image.
  **Wrong claim**: Conway $V^{s\natural}$ is the 5th bosonic
  $\Psi$-image, alongside Monster $V^\natural$ and the Leech /
  Fake Monster / Gritsenko--Nikulin siblings.
  **Precise error**: Conway $V^{s\natural}$ is a super VOA, not a
  bosonic VOA; it lives on the \emph{super-metaplectic}
  $\Psi^{\mathrm{metap}}$ branch via Duncan 2007 \emph{Proc.\ Nat.\
  Acad.\ Sci.} orbifold-diamond construction from the Leech lattice
  fermionic VOA $A(\Lambda_{24})$, not on the bosonic $\Psi$-branch.
  **Counter**: the $\Psi$-sibling ladder has four branches
  $\{\Psi, \Psi^{\deg}, \Psi^{\mathrm{tor}}, \Psi^{\mathrm{metap}}\}$;
  Conway $V^{s\natural}$ sits on the metaplectic branch
  $\Psi^{\mathrm{metap}}$; the four bosonic $\Psi$-images (Monster,
  Fake-Monster-$\mathrm{II}_{25,1}$, Gritsenko--Nikulin paramodular,
  torsional) do NOT include Conway.
  **Related**: AP980; AP-CY71 (Conway $V^{s\natural}$ Duke vs MRL);
  AP-CY72 ($V^{s\natural}$ in two $\Psi$-input rows); Duncan 2007
  \emph{PNAS}; Duncan--Mack-Crane 2015 $V^{s\natural}$ Conway
  construction.

- **AP981 (= AP-$\pi$) -- CHL Siegel weight ladder:
  $k_N = 24/(N + 1) - 2 \in \{10, 6, 4, 3, 2\}$, not
  $\{5, 4, 3, 2, 1\}$ (2026-04-22 Harmonies synthesis).**
  Chaudhuri--Hockney--Lykken (CHL) compactifications of the heterotic
  string on $K3 \times T^2 / \Z_N$ yield Siegel paramodular partition
  functions whose weights follow a specific ladder.
  **Ghost of a true theorem**: CHL compactifications at levels
  $N \in \{1, 2, 3, 5, 7\}$ produce a Siegel-weight ladder indexed
  by $N$.
  **Wrong claim**: $\{\kappa_{\mathrm{BKM}}(\Phi_N)\} = \{5, 4, 3, 2, 1\}$
  read as CHL Siegel weights.
  **Precise error**: the CHL weight formula is
  $k_N = 24/(N + 1) - 2$, giving
  $k_N|_{N = 1, 2, 3, 5, 7} = \{10, 6, 4, 3, 2\}$, not
  $\{5, 4, 3, 2, 1\}$.
  **Counter**: separate the two distinct indexings:
  (i) $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$ (Borcherds weight row,
  canonical preamble); (ii) CHL heterotic Siegel weight
  $k_N = 24/(N + 1) - 2$ (CHL 1995 \emph{PRL}~75).
  **Related**: AP981; AP980 ($\Psi$-sibling branches); AP982
  ($\kappa_{\mathrm{BKM}}$ 5 vs 12 scope); CHL 1995 \emph{PRL}~75
  arXiv:hep-th/9507050; Sen 2007 dyon-counting review.

- **AP982 (= AP-$\rho$) -- $\kappa_{\mathrm{BKM}}$ 5 vs 12 scope
  discipline: $\Psi$-sibling-indexed invariant, AP5 scope qualifier
  mandatory (2026-04-22 Harmonies synthesis).**
  The $\kappa_{\mathrm{BKM}}$ invariant admits two legitimate values
  under different $\Psi$-sibling inputs, each canonical within its
  own scope.
  **Ghost of a true theorem**: $\kappa_{\mathrm{BKM}}$ is a unified
  Borcherds-weight invariant.
  **Wrong claim**: ``Vol I $\kappa_{\mathrm{BKM}} = 12$ contradicts
  Vol III $\kappa_{\mathrm{BKM}} = 5$''.
  **Precise error**: the two values parametrise different
  $\Psi$-siblings: $\kappa_{\mathrm{BKM}} = 5$ is the
  K3-Gritsenko-denominator $\fg_{\Delta_5}$ sibling (paramodular
  $\Phi_{10} = \Delta_5^2$ convention, Gritsenko--Nikulin 1998);
  $\kappa_{\mathrm{BKM}} = 12$ is the Fake-Monster $\fg_{\Phi_{12}}$
  sibling on $\mathrm{II}_{25, 1}$ (Borcherds 1990 \emph{Invent.\
  Math.}~109). Both are correct; neither contradicts the other.
  **Counter**: every $\kappa_{\mathrm{BKM}}$ occurrence names its
  $\Psi$-sibling input (denominator $\Delta_5$ gives 5; denominator
  $\Phi_{12}$ gives 12); the AP5 dual-indexing header in
  `notes/cross_volume_aps.md` is mandatory within 10 lines.
  **Related**: AP982; AP981 (CHL weight ladder); AP980
  ($\Psi$-branches); AP5 cross-volume dual-indexing header; canonical
  preamble ``$\kappa_{\mathrm{BKM}}(\mathbf H_{\Delta_5})$
  cross-volume value'' row; Vol III
  `chapters/examples/cy_d_kappa_stratification.tex`.

- **AP983 (= AP-$\sigma$) -- Yangian type count: four distinct
  Yangian types, not three (2026-04-22 Harmonies synthesis).**
  The Yangian family admits four genuinely distinct deformations of
  $U(\fg[t])$, each living on a different space with different spectral
  parameter structure.
  **Ghost of a true theorem**: Yangians are deformations of $U(\fg[t])$
  classified by their spectral structure.
  **Wrong claim**: ``three Yangian variants''.
  **Precise error**: miscounts the Yangian zoology by one; four types
  coexist in the programme.
  **Counter**: the four distinct Yangian types are:
  (i) classical Drinfeld Yangian $Y(\fg)$ on $\C$ (Drinfeld 1985
  \emph{Dokl.\ Akad.\ Nauk SSSR}~283);
  (ii) chiral Yangian $Y^{\mathrm{ch}}_\fg$ on a curve $C$ with
  $E_1$-chiral structure;
  (iii) spectral Yangian $Y^{\mathrm{spec}}$ on $\R \times \C$
  (Costello--Yagi 2018, 4d Chern--Simons spectral parameter);
  (iv) dg-shifted affine $\widehat Y^{\mathrm{dg}}$ on CY$_d$-categorical
  ambients (Yagi shifted affine Yangian, Soibelman--Yang COHA-Yangian).
  **Related**: AP983; feedback `Four Yangian types`
  (`feedback_yangian_type_distinction.md`); AP-CY programme
  identity entry; Drinfeld 1985; Costello--Yagi 2018 arXiv:1810.01970;
  Schiffmann--Vasserot 2013.

- **AP984 (= AP-$\tau$) -- Two-stage $\Phi_d$ Stage-1/Stage-2 naming
  discipline: every 6d hCS $= E_3$ on $\C^3$ claim names Stage-1
  output and Stage-2 specialisation data (2026-04-22 6d hCS audit).**
  The CY-to-chiral functor factors as
  $\Phi_d = \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C} \circ
  \Phi^{\mathrm{FA}}_d$ (Vol III `cy_to_chiral.tex`): Stage 1 is
  canonical $E_d$-holomorphic factorisation algebra; Stage 2 is
  factorisation homology over a $(d-1)$-cycle $\Sigma_{d-1}$ followed
  by restriction to a reference curve $C$.
  **Ghost of a true theorem**: 6d hCS on a CY$_3$ target carries an
  $E_d$-structure at one end and an $E_1$-chiral structure at the
  other.
  **Wrong claim**: ``6d hCS $= E_3$ on $\C^3$'' (bare, no stage
  naming) OR ``6d hCS produces a chiral algebra'' without naming the
  Stage-2 specialisation $(\Sigma_{d-1}, C)$.
  **Precise error**: the two ends of $\Phi_d$ operate on different
  objects at different stages: Stage 1 output is an $E_d$-holomorphic
  factorisation algebra on the CY$_d$ variety (here $d = 3$, output
  is $E_3$ holomorphic FA on $\C^3$ or the compact CY$_3$); Stage 2
  output is an $E_1$-chiral algebra on a reference curve $C$,
  depending on the choice of $(d - 1)$-cycle $\Sigma_{d - 1}$ along
  which factorisation homology is taken. Conflating these outputs is
  the standard grammar-level error flagged by AP932 / AP933.
  **Counter**: every $\Phi_d$ invocation names:
  (i) Stage-1 $\Phi^{\mathrm{FA}}_d$: CY$_d \to E_d$-holomorphic FA on
  the variety;
  (ii) Stage-2 $\mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C}$:
  factorisation homology over $\Sigma_{d-1}$ + restriction to $C$
  $\to E_1$-chiral algebra on $C$;
  (iii) the specific $(\Sigma_{d-1}, C)$ pair in use --- it is not
  canonical; a single CY$_d$ admits a family of $E_1$-chiral shadows.
  **Related**: AP984; AP932 ($E_n$ on wrong object); AP933 (PTVV
  shift law); AP973 (Dunn additivity on $\C^k$); AP951 (Vol III
  Wave 11-19 two-stage partner); AP963 (item J, wave-12
  cross-programme two-stage anchor); AP-CY62-crossvol (Vol I
  propagation of Vol III R7); Vol III
  `chapters/theory/cy_to_chiral.tex`; CLAUDE.md Vol~III two-stage
  factorisation row.

- **AP985 (= AP-$\upsilon$) -- Cubic $d^{abc}$ vs quadratic $C_2 = 2h^\vee$
  in hCS anomaly: true obstruction is cubic, not quadratic
  (2026-04-22 wave-15 cross-volume propagation from Vol III AP-CY262).**
  The 1-loop BV anomaly of 6d holomorphic Chern--Simons on a
  CY$_3$ target is the genuine obstruction to quantisation; it
  admits two candidate Lie-algebraic witnesses that are NOT
  interchangeable.
  **Ghost of a true theorem**: hCS$_6(\mathfrak g)$ carries a
  1-loop BV anomaly whose vanishing is controlled by an invariant
  of $\mathfrak g$.
  **Wrong claim**: The hCS anomaly is the quadratic Casimir
  $C_2(\mathfrak g) = 2h^\vee$ (dual Coxeter number).
  **Precise error**: $C_2 = 2h^\vee$ is a scheme-dependent
  wave-function coefficient absorbed into a BV-trivial
  counter-term (field redefinition at one-loop). The true
  cocycle-level anomaly is the cubic symmetric invariant
  $d^{abc}(\mathfrak g) = \mathrm{tr}(t^{(a}\{t^b, t^{c)}\})$,
  which vanishes precisely for $\mathfrak{sl}_2$ (no cubic
  invariant: rank 1, $\mathrm{Sym}^3(\mathbf{3})$ contains no
  trivial) and is non-trivial for $\mathfrak{sl}_{N \ge 3}$
  (Gell-Mann $d^{abc}$).
  **Counter**: state the obstruction as cubic $d^{abc}$ at
  cocycle level; the quadratic $C_2 = 2h^\vee$ is the
  renormalisation-scheme coefficient on the counter-term
  $\int \mathrm{tr}(A \partial A)$, BV-trivial after field
  redefinition. $\mathfrak{sl}_2$ hCS is unobstructed;
  $\mathfrak{sl}_{N \ge 3}$ hCS is obstructed by $d^{abc}$.
  **Related**: AP985; AP-CY262 (Vol III sibling); AP932 ($E_n$
  on wrong object); AP966 (6d hCS 1-loop form); AP979 (CHSW
  anomaly-locus native ambient); Costello 2013 arXiv:1303.2632
  \S7; Costello--Li 2016 arXiv:1605.09439; Frampton--Kephart 1983
  \emph{PRL}~50 cubic anomaly classification; Cvitanović 2008
  \emph{Group Theory} Ch.~18.

- **AP986 (= AP-$\phi$) -- PTVV shift $n = d - 4$ for hCS
  $d$-shifted-symplectic structure (2026-04-22 wave-15
  cross-volume propagation from Vol III AP-CY263).**
  Pantev--Toen--Vaquie--Vezzosi (PTVV 2013) shifted-symplectic
  structures on mapping stacks $\mathrm{Map}(X, B\mathfrak g)$
  of hCS data come with a specific shift determined by the
  dimension of $X$.
  **Ghost of a true theorem**: hCS on a CY$_d$ target carries a
  PTVV shifted-symplectic structure of some specific shift $n$.
  **Wrong claim**: generic shift invocation (``shifted
  symplectic'' without naming $n$) or wrong value such as
  $n = -d$ or $n = d$.
  **Precise error**: silently collapses the $d$-dependent shift
  ladder; the correct shift under PTVV for hCS on a CY$_d$ is
  $n = d - 4$. The ladder runs $(d, n) \in
  \{(2, -2), (3, -1), (4, 0), (5, +1), (6, +2)\}$; it terminates
  at $d = 4$ with $n = 0$ (ordinary symplectic structure).
  **Counter**: every hCS PTVV invocation names the triple
  $(d, n = d - 4, E_n^{\mathrm{cl}})$, e.g. $d = 3$ gives
  $(-1)$-shifted symplectic (Joyce--Upmeier CY$_3$ obstruction
  theory); $d = 4$ is the PTVV-$0$ termination (ordinary
  symplectic) distinguishing the hCS anomaly-locus boundary.
  **Related**: AP986; AP-CY263 (Vol III sibling); AP958 / W12E
  (PTVV-shift discipline cross-programme); AP933 (PTVV shift
  law); AP932 ($E_n$ on wrong object); PTVV 2013
  \emph{Publ.\ IHES}~117; Calaque 2015 \emph{Contemp.\ Math.}~643
  (shifted-symplectic mapping stacks); Costello--Gwilliam
  \emph{FA in QFT} Vol 2 \S5.5; Joyce--Upmeier 2017 CY$_3$
  obstruction-theory.

- **AP987 (= AP-$\chi$) -- $\HH^2_{E_2}$ rigidity of
  $\partial \mathrm{hCS}_5(\mathfrak{sl}_n)$ at non-critical
  level via Whitehead~2; critical level $k = -h^\vee$ requires
  separate Feigin--Frenkel argument (2026-04-22 wave-15
  cross-volume propagation from Vol III AP-CY266).**
  The boundary-chiral reduction $\partial \mathrm{hCS}_5$ of 5d
  holomorphic Chern--Simons on $\R \times \C^2$ to a chiral
  algebra on a reference curve $C$ inherits a deformation
  problem whose rigidity is controlled by $E_2$-Hochschild
  cohomology in degree 2.
  **Ghost of a true theorem**: $\partial \mathrm{hCS}_5(\fg)$ is
  $E_2$-Hochschild rigid and therefore uniquely quantisable.
  **Wrong claim**: $\HH^2_{E_2}(\partial \mathrm{hCS}_5(\fg)) = 0$
  is a universal rigidity claim across all levels $k$ including
  critical $k = -h^\vee$.
  **Precise error**: conflates two distinct arguments at
  distinct loci of the level-line. At non-critical level
  $k \ne -h^\vee$, rigidity is Whitehead 2 lemma applied to the
  affine Kac--Moody vertex algebra
  $V_k(\fg)$ (the affine Lie algebra $\widehat{\fg}_k$ is
  semisimple, so $H^2(\widehat{\fg}_k, V_k) = 0$ Whitehead-type
  vanishing). At critical level $k = -h^\vee$, Whitehead fails;
  rigidity requires the separate Feigin--Frenkel 1992 centre
  argument (classical $\mathcal W$-algebra via screening
  operators, $\mathfrak z(\widehat{\fg}) = \pi_0^{\mathrm{ff}}$).
  **Counter**: state rigidity in two clauses: (i) non-critical
  $k \ne -h^\vee$: Whitehead 2 $\Rightarrow$
  $\HH^2_{E_2}(V_k(\fg)) = 0$; (ii) critical $k = -h^\vee$:
  separate Feigin--Frenkel classical-$\mathcal W$ argument
  required. For $\fg = \mathfrak{sl}_n$, the identification
  $V_{k^\vee}(\fg) = \mathcal W^{\mathrm{cl}}(\fg)$ at critical
  level is the substitute for the Whitehead input.
  **Related**: AP987; AP-CY266 (Vol III sibling); AP972
  (critical-level $V_k(\fg)$ screening); AP947 (Feigin--Frenkel
  centre); AP932 ($E_n$ on wrong object); Whitehead 1937
  \emph{Ann.\ Math.}~38; Feigin--Frenkel 1992 \emph{CMP}~128;
  Frenkel 2005 \emph{BAMS}~42 critical-level review; Costello--Gaiotto
  2021 arXiv:2108.01379 5d hCS / Yangian quantisation.

- **AP988 (= AP-$\psi$) -- Chiral Booth--Lazarev three-obstruction
  tower collapses to two via Mumford boundary-descent
  ($\kappa_{\mathrm{ch}}$-curving IS the O2 datum) (2026-04-22
  wave-15 cross-volume propagation from Vol III AP-CY268).**
  The chiral Booth--Lazarev programme names three obstructions
  (O1) Ran--Smith descent, (O2) genus-tower compatibility,
  (O3) IndHilb torsor gluing; a first-principles audit of the
  genus-boundary stratification collapses O2 onto O1+O3.
  **Ghost of a true theorem**: the chiral Booth--Lazarev lifting
  problem carries independent obstructions classified by
  $(O_1, O_2, O_3)$.
  **Wrong claim**: all three obstructions $(O_1, O_2, O_3)$ must
  be verified independently.
  **Precise error**: O2 (genus-tower compatibility) is NOT
  independent of O1 (Ran--Smith) + O3 (IndHilb): it collapses
  via Mumford's boundary-divisor identity on
  $\overline{\mathcal M}_g$: $\partial^* \lambda_g^1 =
  \pi_1^* \lambda_{g_1}^1 + \pi_2^* \lambda_{g_2}^1$ on the
  irreducible-boundary divisor $\Delta_0$ (and similarly on
  separating $\Delta_{\mathrm{sep}}^{g_1, g_2}$), because the
  $\kappa_{\mathrm{ch}}$-curving of the chiral algebra over
  $\mathcal M_g$ IS the boundary-descent datum implementing O2.
  The curving data is not an extra obstruction; it is automatic
  from the $\kappa_{\mathrm{ch}}$-subscript indexing of the
  chiral bundle's Chern-class.
  **Counter**: state the Booth--Lazarev obstruction tower as
  (O1) + (O3) only; absorb the former O2 (genus-tower
  compatibility) into the $\kappa_{\mathrm{ch}}$-curving datum
  whose boundary-descent along Mumford's
  $\partial^* \lambda_g^1 = \pi_1^* \lambda_{g_1}^1 +
  \pi_2^* \lambda_{g_2}^1$ is automatic. The two genuine
  obstructions are (O1) Ran--Smith chiral-algebra descent across
  the Ran space and (O3) IndHilb torsor gluing at the
  Hilbert-scheme moduli of fibrewise data.
  **Related**: AP988; AP-CY268 (Vol III sibling); Booth--Lazarev
  2023 arXiv:2304.08409 (chiral Koszul duality coderived
  programme); Mumford 1983 \emph{Arithmetic and Geometry}~II
  boundary-divisor identity on $\overline{\mathcal M}_g$;
  AP-CY61 (first-principles investigation principle); AP961 /
  W12H (two-anchor citation discipline); Vol III
  `chapters/examples/k3_chiral_bialgebra_platonic.tex`
  $\kappa_{\mathrm{ch}}$-curving row.

## Vol III conifold master synthesis cross-programme propagation: AP2141 through AP2143 (2026-04-23)

Source: Vol III `notes/master_synthesis_coha_conifold_2026_04_23.tex` Wave 1-5 adversarial audit. Siblings of Vol III AP-CY298, AP-CY299, AP-CY305.

**AP2141 — $\widehat{\mathfrak{sl}}_2$ vs $\widehat{\mathfrak{gl}}(1|1)$ scope-discipline for affine super BPS algebras (Medium, sibling of AP-CY305).**
Wrong: treating $\mathrm{CoHA}(Y)^{\mathrm{conifold}} \cong Y^+(\widehat{\mathfrak{sl}}_2)$ and $\mathrm{CoHA}(Y)^{\mathrm{conifold}} \cong Y^+(\widehat{\mathfrak{gl}}(1|1))$ as equivocal alternatives. Correct: both are valid theorem-grade identifications at different scopes — $Y^+(\widehat{\mathfrak{gl}}(1|1))$ is primary (Li-Yamazaki arXiv:2003.08909 §8.3.6.3, Gaiotto-Rapčák arXiv:1703.00982 $Y_{0,1,1}[\psi]$), $Y^+(\widehat{\mathfrak{sl}}_2)$ is the ungraded shadow via Chevalley-fixed subalgebra (supertrace projection quotienting the central Cartan $K_n = h_n^{(0)} + h_n^{(1)}$). Agreement on real roots (numerical BPS $\Omega = 1$); differ on imaginary (supercount $\Omega = -2$ vs ordinary dim $2$). **Counter**: bar-cobar applications to super affine algebras must declare scope: graded super (primary) vs ungraded shadow.

**AP2142 — Super 5D hCS $\to$ Yangian VOA all-orders is OPEN, not a Costello-Gaiotto-Yagi theorem (High, sibling of AP-CY298).**
Wrong: citing Costello-Gaiotto-Yagi arXiv:1810.01970 all-orders 5D hCS → Yangian VOA theorem as applying to Lie superalgebras $\mathfrak g = \mathfrak{gl}(m|n)$ or $\widehat{\mathfrak{gl}}(1|1)$. Correct: CGY requires simply-laced $\mathfrak g$ with non-degenerate even Killing form; $\mathfrak{gl}(1|1)$ Killing form is degenerate ($\mathrm{str}(K \cdot \text{anything}) = 0$); super-KT formality is $E_2$-only (Ginzburg-Schedler arXiv:0807.0174), not $E_3$; higher-loop $H^1_{\mathrm{loc}}(\mathfrak{gl}(1|1), \mathcal O_{\mathrm{loc}})$ obstructions not auto-killed. Status: OPEN at all orders for super. 1-loop wheel vanishing VERIFIED; 2-loop+ open. **Counter**: every "5D hCS $\to$ Yangian VOA" invocation for super $\mathfrak g$ must carry OPEN flag; abelian $\widehat{\mathfrak{gl}}_1$ per-chart case is all-orders rigorous, super assembly is conjectural.

**AP2143 — Super Kontsevich-Tamarkin formality is $E_2$-only, not $E_3$ (Medium, sibling of AP-CY299).**
Wrong: treating Kontsevich-Tamarkin $E_n$-formality (Kontsevich arXiv:q-alg/9709040, Tamarkin arXiv:math/9803025) as extending to Lie superalgebras with $E_3$-formality. Correct: Ginzburg-Schedler arXiv:0807.0174 super-Koszul-duality establishes $E_2$-formality for super-Poisson structures; $E_3$ formality is not proved for super. Technical obstruction to extending CGY's all-orders 5D hCS theorem to super-gauge case. Relevant wherever bar-cobar on super affine algebras invokes $E_3$-formality (most commonly in Maurer-Cartan transfer + A_∞-to-L_∞ conversions on super dg-Lie). **Counter**: "super KT formality" must specify $E_2$; $E_3$-requiring claims on super dg-Lie/L_∞ algebras are open.

## Vol III conifold master synthesis cross-programme propagation (cont'd): AP2144 through AP2156 (2026-04-23)

Source: Vol III Wave-5 adversarial audit on `notes/antipatterns_catalogue.md` AP-CY293–AP-CY344, selecting entries with chiral-side / bar-cobar / Borcherds-lift relevance to Vol I.

**AP2144 — Two scopes of the universal Borcherds identity $\kappa_{\mathrm{BKM}}(\Phi_N) = c_N(0)/2$ (Medium, sibling of AP-CY334).**
Wrong: citing the universal Borcherds identity in Vol I canonical preamble / landscape census with a single constant-term ladder. Correct: the identity lives on TWO distinct scopes — (A) the CHL ladder $N \in \{1, 2, 3, 4, 6\}$ with $c_N(0) = (10, 8, 6, 4, 2)$ giving $\kappa_{\mathrm{BKM}} = (5, 4, 3, 2, 1)$ (Gritsenko-Nikulin 1995 arXiv:alg-geom/9504006 Pt II Thm 2.1, EH 2011, GK 2010); (B) the Gritsenko-Cléry 8-form atlas with constant terms $(10, 4, 2, 2, 1, 2, 1/2, 0)$ and cover stratified by $\mathrm{Sp}_4(\mathbb Z)$ / $\mathrm{Mp}_4$ / $\widetilde{\mathrm{Mp}}_4$. Vol I bar-cobar statements that quote "$\kappa_{\mathrm{BKM}} = c_N(0)/2$" without scope label drop the $N = 2, 3, 4$ CHL values and produce wrong ladder numerics at those indices. **Counter**: every Vol I invocation of the universal Borcherds identity specifies scope (A) or (B); bare citation is ambiguous.

**AP2145 — Saito-Kurokawa lift target is $\Phi_{10} = \Delta_5^2$; rescale factor $4$, not $2$ (Medium, sibling of AP-CY335).**
Wrong: writing the Saito-Kurokawa lift on Vol I paramodular / chiral-lift content with target $\Delta_5$ and weight-rescale factor $2$ (i.e., mapping programme parameter $s = k/2$ to Andrianov $s = k$). Correct: SK lifts square Siegel cusp forms; the target is $\Phi_{10} = \Delta_5^2$ (weight $10$), the elliptic source is $g = \Delta \cdot E_6 \in S_{18}$, and the rescale is factor $\mathbf 4$ (two factors of $2$: Siegel weight doubling under SK + Andrianov convention $L_{\mathrm{spin}}(s) = \zeta(s - k + 1) \zeta(s - k + 2) L(s, g)$). Residue at $s = 10$: $\mathrm{Res}_{s = 10} L_{\mathrm{spin}}(s, \Phi_{10}) = -15120 \cdot a_{10}(g) \cdot \Omega^-(g)$ with $\Omega^-$ the Manin minus-period (Ichino-Ikeda Deligne period). **Counter**: Vol I SK-lift invocations on Siegel / paramodular content identify lift target as $\Phi_{2k} = \Delta_k^2$, rescale factor $4$; never as $\Delta_k$ with factor $2$.

**AP2146 — Canonical Faddeev-Kashaev pentagon form on chiral side (Low, sibling of AP-CY330).**
Wrong: stating the pentagon identity on Vol I wall-crossing / MMNS / bar-cobar products as "$\Psi(x) \Psi(y) = \Psi(y) \Psi(z) \Psi(x)$" with a detached variable $z$. Correct: Faddeev-Kashaev hep-th/9310070 canonical pentagon reads $\Psi(x_0) \Psi(x_1) = \Psi(x_1) \Psi(q^{-1/2} x_0 x_1) \Psi(x_0)$ — the middle factor argument $q^{-1/2} x_0 x_1$ is DETERMINED by the motivic commutator of the two outer generators, not a free variable. On a conifold-type bar-cobar product: $x_0 = \hat x_{[S_0]}$, $x_1 = \hat x_{[S_1]}$, middle $= \hat x_{[S_0] + [S_1]}$ with BPS invariant $\Omega(\gamma_{S_0} + \gamma_{S_1}) = 1$ (MMNS arXiv:1107.5017). Derivation via BPS-index invariance + Kashaev-Nakanishi arXiv:1104.4630. **Counter**: every pentagon-identity statement in Vol I wall-crossing or chiral-factorisation content uses the canonical FK form with determined middle argument.

**AP2147 — Bialgebra vs Hopf discipline on chiral-side positive halves (Medium, sibling of AP-CY327).**
Wrong: labelling $Y^+(\mathfrak g)$ (the positive half of a Yangian / quantum affine / quantum toroidal algebra) as a "Hopf algebra" in Vol I bar-cobar content. Correct: $Y^+$ is a bialgebra (associative multiplication + coassociative coproduct + compatible counit), NOT a Hopf algebra; the antipode $S$ lives only on the Drinfeld double $D(Y^+) = Y$. The $\mathbb Z_2$-graded structure (super positive halves from e.g. $\widehat{\mathfrak{gl}}(1|1)$) is bialgebra-in-super-vector-spaces; again, antipode only on the double. Strict-Hopf holds on rational/trigonometric/toroidal equivariance strata; elliptic is quasi-Hopf with Felder-Jimbo-Konno dynamical associator. **Counter**: Vol I bar-cobar on positive-half quantum affinisations must say "bialgebra" (or "$\mathbb Z_2$-graded bialgebra" / "quasi-Hopf" per stratum), not blanket "Hopf".

**AP2148 — Strict Hopf vs quasi-Hopf discipline stratified by equivariance stratum (Low, sibling of AP-CY329).**
Wrong: writing "$Y$ is a Hopf algebra" on Vol I quantum-group content without specifying which equivariance stratum (rational / trigonometric / toroidal / elliptic). Correct: rational / trigonometric / toroidal Yangian families carry strict $\mathbb Z_2$-graded Hopf structure; the elliptic family (Felder-Jimbo-Konno dynamical quantum groups) is **quasi-Hopf** with coassociator satisfying pentagon up to a dynamical $\Phi$-twist depending on a Cartan-valued parameter — not a scalar cocycle. **Counter**: Hopf-status claims on affine (super) algebras in Vol I carry stratum label: strict-Hopf (rational / trig / toroidal) vs quasi-Hopf dynamical (elliptic).

**AP2149 — Super-shuffle is bigraded $\mathrm{Sym} \otimes \Lambda$, not rational $S_m \times S_n$ invariants (Medium, sibling of AP-CY326).**
Wrong: defining the super-shuffle algebra in Vol I bar-cobar content on $\widehat{\mathfrak{gl}}(1|1)$-type positive halves as rational invariants $\mathbb C(z^{(0)} | z^{(1)})^{S_m \times S_n}$. Correct: Feigin-Odesskii alg-geom/9610001 §1.3-1.4 gives the super-shuffle as bigraded polynomial $\otimes$ exterior, $\mathrm{Sh}^{\mathrm{super}}_{m,n} = (\mathbb C[z^{(0)}_1, \dots, z^{(0)}_m] \otimes \Lambda[z^{(1)}_1, \dots, z^{(1)}_n])^{S_m \times S_n}$ — polynomial on the even colour, Grassmann on the odd colour, symmetrised / antisymmetrised separately with Koszul sign $(-1)^{\mathrm{inv}(\tau)}$. The Davison embedding of super-CoHA sits as the wheel-quotient $\mathrm{Sh}^{\mathrm{super}} / I_{\mathrm{wheel}}$ with explicit wheel elements at bidegrees $(2, 1)$ and $(1, 2)$ (Davison arXiv:1311.6989 Thm~A). **Counter**: Vol I bar-cobar on super positive halves uses the Sym $\otimes$ $\Lambda$ bigrading with Koszul sign, never rational invariants.

**AP2150 — $\widehat{\mathfrak{gl}}(1|1)$ central extension is rank-one, central only at $m + n = 0$ (High, sibling of AP-CY320).**
Wrong: Vol I bar-cobar Maurer-Cartan expansions on $\widehat{\mathfrak{gl}}(1|1)$-type positive halves writing super-bracket level-$K$ corrections at arbitrary Fourier mode (e.g., $\{e^{(1)}_m, f^{(1)}_n\}_+ = H_{m+n} + \mathrm{sgn}(m - n) K_{m+n}$). Correct: the universal central extension of $\widehat{\mathfrak{gl}}(1|1)$ is rank-one (Kac 1977 Thm 8.6); the 2-cocycle reads $\omega(x t^m, y t^n) = m \delta_{m+n, 0} \mathrm{str}(xy) K_0$, so all level-$K$ corrections collapse to the zero Fourier mode. The supertrace $h^{\mathrm{tr}} = h^{(0)} + h^{(1)}$ gives a Heisenberg ideal (central in $\widehat{\mathfrak{sl}}(1|1)$ only), distinct from the bona fide central $K_0$ of $\widehat{\mathfrak{gl}}(1|1)$. **Counter**: every level-$K$ term in Vol I super bar-cobar pairings attaches to the $m + n = 0$ Fourier mode only.

**AP2151 — $\widehat{\mathfrak{gl}}(1|1)$ level OPE is cross Cartan $H \cdot N$, not diagonal $N \cdot N$ (High, sibling of AP-CY307).**
Wrong: writing $N(z) N(w) \sim k / (z - w)^2$ as the Vol I Cartan OPE of $\widehat{\mathfrak{gl}}(1|1)$ on chiral-factorisation content. Correct: $\mathfrak{gl}(1|1)$ has rank-$2$ Cartan $\{H, N\}$ with $H$ semisimple (eigenvalues $\pm 1$ on $\psi^\pm$) and $N$ central (annihilates $\psi^\pm$); the invariant super-bilinear form (supertrace on the defining representation) makes both self-isotropic: $\mathrm{str}(H^2) = 0 = \mathrm{str}(N^2)$; the cross-pairing $\mathrm{str}(H N) = 2 \neq 0$ carries the level. Hence: $H(z) H(w) \sim 0$, $N(z) N(w) \sim 0$, $H(z) N(w) \sim k / (z - w)^2$; $H(z) \psi^\pm(w) \sim \pm \psi^\pm(w) / (z - w)$, $N(z) \psi^\pm(w) \sim 0$. **Counter**: Vol I level-$k$ super Cartan OPEs derive the form from the invariant super-bilinear structure explicitly; diagonal ansätze analogous to $\mathfrak{sl}_n$ are category errors.

**AP2152 — Ungraded shadow $Y^+(\widehat{\mathfrak{gl}}(1|1)) \twoheadrightarrow Y^+(\widehat{\mathfrak{sl}}_2)$ is a surjection, not an iso (Medium, sibling of AP-CY328, refines AP2141).**
Wrong: stating $Y^+(\widehat{\mathfrak{gl}}(1|1))^{\mathrm{con}} \cong Y^+(\widehat{\mathfrak{sl}}_2)^{\mathrm{con}}$ as bialgebra isomorphism in Vol I bar-cobar content. Correct: the supertrace projection $\mathfrak{gl}(1|1) \twoheadrightarrow \mathfrak{sl}(1|1) / \langle K \rangle$ induces a SURJECTION (not iso) of bialgebras whose kernel is the two-sided ideal generated by the central $K$. The super source has a 2-dim imaginary-root line spanned by $\{H, K\}$; the shadow has a 1-dim imaginary-root line $\{H\}$ only. **Counter**: Pattern 273 scope discipline in Vol I must distinguish "super source" from "ungraded shadow" as distinct bialgebras connected by a surjection.

**AP2153 — Super DM integrality: BPS multiplicities carry $(\text{even}|\text{odd})$ split (Medium, sibling of AP-CY322).**
Wrong: quoting Davison-Meinhardt integrality on super BPS Lie algebra content in Vol I bar-cobar as "multiplicity matches MMNS" without specifying super-dimension signature. Correct: on conifold-type super BPS Lie algebra, real roots $\pm \alpha$ carry super-dimension $(0|1)$ (fermionic odd, multiplicity 1); imaginary roots $n \delta$ carry $(2|0)$ (bosonic even, 2-dim from $\Omega^{\mathrm{mot}}(n \delta) = -\mathbb L - 1$, Cartan spanned by $H$ + $K$). Matches MMNS via DM-integrality (Davison-Meinhardt arXiv:1601.02479 Thm A). **Counter**: super BPS Lie algebra multiplicities in Vol I bar-cobar carry the $(\text{even}|\text{odd})$ split; total-dimension-only statements lose the $\mathbb Z/2$-grading.

**AP2154 — $\kappa_{\mathrm{ch,BV}}$ is a distinct subscripted invariant from $\kappa_{\mathrm{ch}}$ and $\kappa_{\mathrm{cat}}$ on non-compact CY$_3$ (High, sibling of AP-CY324).**
Wrong: propagating apparent three-valuedness $\{+1, 0, -1\}$ for "$\kappa_{\mathrm{ch}}$" across Vol I landscape census / canonical preamble / bar-cobar rows on non-compact CY$_3$ content. Correct: the three values are three DISTINCT subscripted invariants — $\kappa_{\mathrm{ch}}$ (Hodge supertrace / DT count), $\kappa_{\mathrm{cat}} = \chi(\mathcal O_X)$ (holomorphic Euler characteristic, requires coherence-of-sheaf scope on non-compact), $\kappa_{\mathrm{ch,BV}}$ (Costello-Li BCOV 1-loop curving via ghost supertrace; NEW subscript previously absent from the four-$\kappa$ discipline). The Polyakov ghost-mode balance $\kappa_{\mathrm{ch}} + \kappa_{\mathrm{ch,BV}} = \kappa_{\mathrm{cat}}$ holds on $\mathsf G$-class free-field CY$_3$ with $\chi_{\mathrm{top}} = 2$; NOT universal. **Counter**: the four-$\kappa$ discipline in Vol I canonical preamble extends to a five-$\kappa$ ladder $\{\kappa_{\mathrm{ch}}, \kappa_{\mathrm{cat}}, \kappa_{\mathrm{BKM}}, \kappa_{\mathrm{fiber}}, \kappa_{\mathrm{ch,BV}}\}$ on non-compact CY$_3$ content; ghost-BV subscript is load-bearing.

**AP2155 — BCOV cocycle target is $H^{0,1}(X, \mathcal O_X)$, not $H^{0,1}(X, \mathrm{Sym}^{\leq 2} T_X^*)$ (Medium, sibling of AP-CY341).**
Wrong: writing the BCOV 1-loop anomaly class as $\alpha_{\mathrm{BCOV}}(X) = (\chi(X) / 24)[\Omega_X]^{0,1} \in H^{0,1}(X, \mathrm{Sym}^{\leq 2} T_X^*)$ on Vol I landscape / chiral-factorisation content. Correct: the BCOV anomaly class lives in $H^{0,1}(X, \mathcal O_X)$. The $[\Omega_X]^{0,1}$ is the $(0, 1)$-Dolbeault lift of the Atiyah class of the CY trivialisation $\omega_X \cong \mathcal O_X$ — an $\mathcal O_X$-coefficient class, not a polyvector-field class. Primary: Costello-Li arXiv:1606.00365 Prop 5.2; BCOV arXiv:hep-th/9309140; Polyakov 1981 Ch.~9. **Counter**: BCOV cocycle-target claims in Vol I place the $(\chi / 24)$ class in $H^{0,1}(X, \mathcal O_X)$; higher-symbol targets are category errors.

**AP2156 — Super root $\alpha$ isotropy derived at pairing level, not cited only (Medium, sibling of AP-CY321).**
Wrong: asserting "$(\alpha, \alpha) = 0$ on the isotropic super root" in Vol I super bar-cobar via Kac 1977 Thm 2.4 citation alone, without the supertrace computation. Correct: for $\alpha = \epsilon_1 - \epsilon_2$ on $\mathfrak{gl}(1|1)$, $(\alpha, \alpha) = (\epsilon_1, \epsilon_1) - 2 (\epsilon_1, \epsilon_2) + (\epsilon_2, \epsilon_2) = 1 - 0 + (-1) = 0$ under the supertrace-induced bilinear form ($+1$ on even diagonal, $-1$ on odd diagonal). Affine lift inherits isotropy via $(\delta, \delta) = (\delta, \alpha) = 0$. **Counter**: every "$\alpha$ isotropic" statement in Vol I super content displays the supertrace computation once; citation alone hides the zero's source.

## Vol III finite-Rees hCS--Hall and automorphic-product propagation: AP2157 through AP2168 (2026-04-30)

Source: Vol III finite DWR/Ran Rees hCS--Hall construction and compact-CoHA critique reconciliation, mirrored for Vol I bar--cobar and chiral quantum-group uses.

**AP2157 -- Finite Rees hCS--Hall is not ordinary compact critical CoHA (High, sibling of AP-CY345).**
Wrong: using a finite simplex-integrated hCS--Hall map as the compact critical-CoHA comparison. Correct: the finite construction gives a Maurer--Cartan element in a bounded DWR/Ran convolution dg Lie algebra valued in an oriented Rees Hall complex. Compact critical CoHA requires pro-compatible inverse limit with Mittag--Leffler control and monoidal vanishing-cycle realization preserving Hall correspondences, orientations, shifts, Tate twists, and Thom--Sebastiani. **Counter**: Vol I bar--cobar arguments name which layer is used: finite Rees, completed Rees, realized critical CoHA, positive half, or Drinfeld double.

**AP2158 -- Multi-chart gluing exists at finite Rees level (Medium, sibling of AP-CY346).**
Wrong: repeating that no explicit gluing exists outside the single $\mathbb C^3$ chart. Correct: finite multi-chart gluing is supplied by face-compatible cyclic contractions over $\Omega^\bullet(\Delta^p)$, Fourier--BV Rees identification, and Stokes assembly into a Maurer--Cartan element. The remaining open step is completed/realized compact CoHA, not finite gluing. **Counter**: Vol I gluing claims distinguish finite DWR/Ran descent from compact realization.

**AP2159 -- $\CoHA(\mathbb C^3)$ is the positive half, not $\mathcal W_{1+\infty}$ (High, sibling of AP-CY347).**
Wrong: writing $\CoHA(\mathbb C^3)=\mathcal W_{1+\infty}$. Correct: $\CoHA(\mathbb C^3)\cong Y^+(\widehat{\mathfrak{gl}}_1)$; the full Yangian/Drinfeld double and a Fock/vacuum-module evaluation produce the $\mathcal W_{1+\infty}$ vertex algebra. **Counter**: Vol I quantum-group rows separate positive half, double, representation, and vertex algebra.

**AP2160 -- Toric vertices do not by themselves give the simplex map (Medium, sibling of AP-CY348).**
Wrong: toric CY$_3$ comparison follows only from affine charts $U_\sigma\simeq\mathbb C^3/G_\sigma$. Correct: vertices supply local models; the finite comparison also requires face-compatible contractions, mutation coherence, orientation transport, Thom--Sebastiani compatibility, and Hall product compatibility on all DWR/Ran simplices. **Counter**: Vol I toric statements name the full finite cyclic atlas, not only the vertex charts.

**AP2161 -- Non-formality of local $\mathbb P^2$ is absorbed into the cyclic potential (Medium, sibling of AP-CY349).**
Wrong: local $\mathbb P^2$ is blocked merely because its category is non-formal. Correct: higher $m_k$ operations enter the finite cyclic potential $W_\sigma$ and hence the Rees critical complex. The real obligations are cyclic contraction, mutation coherence, completion, and realization. **Counter**: non-formality is data inside $W_\sigma$, not a terminal obstruction.

**AP2162 -- Oberdieck--Pixton and Oberdieck--Pandharipande labels are not interchangeable (Medium, sibling of AP-CY350).**
Wrong: using Oberdieck--Pixton for every reduced DT theorem anchor. Correct: use OPi only when Pixton's component is used; use Oberdieck--Pandharipande or the year-specific Oberdieck theorem for reduced DT anchors on $K3\times E$. **Counter**: Vol I reduced-DT citations name the actual theorem source.

**AP2163 -- Quasi-NCCR characters do not construct compact critical CoHA (High, sibling of AP-CY351).**
Wrong: a quasi-NCCR or compact character formula constructs the compact critical CoHA. Correct: a character equality can test a constructed object; it does not supply compact-support descent, monoidal realization, orientation transport, or a comparison morphism. **Counter**: Vol I character identities are evidence or consequences, not construction of the compact CoHA layer.

**AP2164 -- Layer-aware $Y^+$ notation (High, sibling of AP-CY352).**
Wrong: one notation $Y^+(X)$ silently passes through finite Rees, completed Rees, realized CoHA, positive half, and Drinfeld double. Correct: finite Rees Hall complex, completed Rees Hall algebra, realized positive critical CoHA, and Hall--Drinfeld double are distinct objects connected by named arrows. **Counter**: every cross-volume use states the arrow crossed.

**AP2165 -- Global hCS--Hall has four cases, not yes/no (Medium, sibling of AP-CY353).**
Wrong: asking whether global hCS--Hall is simply constructed or not constructed. Correct: $\mathbb C^3$ positive half is constructed; finite multi-chart Rees is constructed under cyclic-atlas hypotheses; completed Rees is conditional on pro-compatibility and Mittag--Leffler; realized compact critical CoHA is conditional on monoidal vanishing-cycle realization. **Counter**: Vol I summaries keep this four-case stratification.

**AP2166 -- CHL and Gritsenko--Clery constant-term ladders are separate (Medium, sibling of AP-CY354).**
Wrong: quoting $\kappa_{\mathrm{BKM}}=c(0)/2$ without specifying the CHL ladder or the Gritsenko--Clery atlas. Correct: both obey the weight formula but use different indexing families and constant terms. **Counter**: name CHL or Gritsenko--Clery before comparing constants.

**AP2167 -- $\Delta_5$, $\Phi_{10}$, and $\Phi_{12}$ have different roles (High, sibling of AP-CY355).**
Wrong: interchanging the K3 paramodular denominator, its square, and the Fake-Monster denominator. Correct: $\Delta_5$ is the primitive weight-5 K3 paramodular BKM denominator; $\Phi_{10}=\Delta_5^2$ is the compact reduced-DT / Igusa square; $\Phi_{12}$ is the Fake-Monster product on $\mathrm{II}_{25,1}$. **Counter**: every trace or partition-function statement names which automorphic product appears.

**AP2168 -- Determinant Hodge line notation is $\lambda_1^{\det}$ (Low, sibling of AP-CY356).**
Wrong: writing the determinant Hodge line as $\lambda_g^1$ in the BL/DWR bridge. Correct: use $\lambda_1^{\det}$ for the determinant Hodge line; reserve genus subscripts for actual genus-indexed Hodge bundles. **Counter**: Vol I boundary-curving claims use determinant-line notation.

**AP2169 -- OP/DT scalar normalization uses $D_5=64^{-1}\Delta_5$ (High, sibling of AP-CY357).**
Wrong: writing the \(K3\times E\) reduced-DT scalar as bare
\(-\Delta_5^{-2}\) or unqualified \(-\Phi_{10}^{-1}\). Correct:
\(\Delta_5\) is the primitive BKM denominator, while the OP scalar
branch uses the monic Borcherds product
\[
D_5=64^{-1}\Delta_5,\qquad
\chi_{10}^{\mathrm{OP}}=D_5^2=4096^{-1}\Delta_5^2.
\]
Thus \(Z_{\mathrm{OP/DT}}=-(\chi_{10}^{\mathrm{OP}})^{-1}
=-D_5^{-2}=-4096\,\Delta_5^{-2}\). **Counter**: Vol I trace,
bar-Euler, Hochschild, and DT passages use \(\chi_{10}^{\mathrm{OP}}\)
or \(D_5\) when referring to the OP/DT scalar; the unnormalised Igusa
square \(\Phi_{10}^{\mathrm{un}}=\Delta_5^2\) is a different convention.

## Topological-strings Hamiltonian BF / mixed HT critique propagation: AP2170 through AP2177 (2026-04-30)

Source: `~/topological-strings` refreshed 02:14 SAST
`Mixed Holomorphic-Topological strings critique.pdf`, mirrored here only
as cross-programme anti-pattern discipline.  The critique is not a
source of truth; each entry records a failure mode requiring
first-principles proof in the target manuscript.

**AP2170 -- PBW special fibre is not deformation-level Koszul duality (High).**
Wrong: titling a stable Hamiltonian trace comparison as Koszul duality
when the proof identifies only the associated graded PBW special fibre
\(U_\hbar(\mathfrak g)/(\hbar)\cong S(\mathfrak g)\). Correct: call this
a PBW/Rees shadow unless a filtered deformation-level algebra map,
module action, and completed bar-cobar adjunction are constructed.
**Counter**: every closed/open Koszul claim states which layer is being
used: operadic, CE/PV, PBW/Rees, Matlis, factorization-current, or
stable horizon.

**AP2171 -- Hamiltonian CE coordinate is not the boundary Hamiltonian source (High).**
Wrong: sending the CE coordinate \(c^f\) dual to a Hamiltonian generator
to the boundary Hamiltonian \(B_f\). Correct:
\(c^f\mapsto\theta^f\) is the constant polyvector coordinate, while the
shifted cotangent coordinate \(u_f\mapsto O_f\) gives the boundary
Hamiltonian function; \(X_{O_f}=d_\pi O_f=\Phi(d_{\rm CE}u_f)\).
**Counter**: every boundary-source statement distinguishes \(c\),
\(u\), \(\theta\), and \(O\).

**AP2172 -- Polynomial one-\(\psi\) descendants are not Matlis principal parts (High).**
Wrong: identifying \(\Psi_{a,b}\) with \(\rho_{a,b}\) because their
labels match. Correct: \(\Psi_{a,b}\) are PBW-special-fibre open
descendant labels; \(\rho_{a,b}=z_1^{-a-1}z_2^{-b-1}dz_1dz_2\) are
residue principal parts carrying the coadjoint module. Same labels do
not give the same topology, \(\mathfrak h\)-action, or deformation
object. **Counter**: a \(\Psi\to\rho\) bridge must name a genuine module
map and prove action compatibility; same-action Rees/Fourier bridges are
not presumed.

**AP2173 -- Scalar reduction is not balanced-supertrace anomaly cancellation (High).**
Wrong: treating projection away from scalar Hamiltonians as physical
removal of the \(U(1)\) scalar anomaly. Correct: ordinary
\(\mathfrak{gl}_N\) has
\(\{\mathrm{Tr}\phi_1,\mathrm{Tr}\phi_2\}=N\) and quantum scalar class
\(\hbar N[\bar c]\); a balanced supertrace brane
\(\mathbb C^{N|N}\) is a different source choice and must be proved by
direct supertrace cancellation. **Counter**: every bracket theorem says
before/after scalar quotient and ordinary/balanced branch.

**AP2174 -- Strict product/direct-sum endpoint does not supply the Costello kernel (High).**
Wrong: using the unweighted product/direct-sum formal Hamiltonian pair as
the analytic Costello coefficient system. Correct: a Costello BV kernel
requires weighted or otherwise kernel-admissible coefficient spaces with
proved bracket continuity, coadjoint continuity, and convergence of the
diagonal Casimir \(C_{-,+}=\sum_I H_I\otimes\rho_I\). **Counter**:
strict endpoint claims are obstruction statements unless the diagonal
kernel is constructed in the stated topology.

**AP2175 -- Reduced principal-part currents do not construct the unreduced factorization center (High).**
Wrong: inferring an unreduced \(Z^{P_0}_{\rm fact}(Obs_{\rm open})\)
map from a reduced post-contraction current theorem. Correct: the
unreduced lift needs compactly supported de Rham current sources,
principal-part current fields, product and bracket centrality homotopies,
and coherence for arbitrary open observables. **Counter**: reduced
current maps are labelled as post-contraction data until the unreduced
centrality complex is built.

**AP2176 -- First/third Moyal graph checks are not all-order Costello/QME realization (High).**
Wrong: treating the first coefficient and the \(1/24\) third coefficient
as an analytic Costello graph theorem. Correct: the formal Moyal product
is an algebraic target; Costello realization requires an elliptic field
complex, heat kernel, propagator, RG flow, counterterms, scalar and
non-scalar QME classes, and proof that no higher connected binary
cumulants survive. **Counter**: graph normalization sufficiency is
conditional on the renormalized binary contraction and cumulant
vanishing hypotheses.

**AP2177 -- Local Hamiltonian BF does not transfer to compact CY / Igusa / MNOP without matched descent (High).**
Wrong: extracting compact Calabi--Yau, global BCOV, Vol III, Igusa/BKM,
MNOP, or sister-volume consequences from a local formal-disk theorem.
Correct: supply a matched-conventions target datum, obstruction vector,
and null-homotopies for degrees, pairings, periods, central classes,
factorization products, QME, convergence, and descent. **Counter**: no
target datum and no null-homotopies means no global consequence.

**AP2178 -- Topological-string stable ranges and homology models are degree-specific (Medium).**
Wrong: importing a stable range, homology statement, or topology label
without the degree and coefficient convention used in the proof. Correct:
\(N\ge d\) is supported for scalar trace-coordinate invariants when
\(d\) is ordinary matrix-letter trace word degree; it is not a
\(\psi\)-weight, concomitant, or large-\(N\) physics theorem. The
one-\(\psi\) Abel-map proof gives homology equivalence unless a stronger
homotopy proof is supplied. Matlis continuous duals use finite Taylor
support / locally finite direct-sum topology, and the finite Dirac brane
trace algebra is ghost-zero reduced BV/Koszul, not full open BV.
**Counter**: every stable-range, one-\(\psi\), Matlis, or Dirac statement
names its degree convention, coefficient topology, and truncation.

**AP2179 -- Admissible Koszul large-\(N\) is not physical large-\(N\) (High).**
Wrong: replacing the raw finite-\(N\) brane limit by
\(U_\hbar(\bar A_\hbar\ltimes P_\hbar[1])\) and then calling this the
physical large-\(N\) theorem. Correct: the enveloping algebra is a
candidate admissible Koszul object; it becomes theorem-grade only inside
an explicit filtered/bar-cobar category and after a comparison from
finite-\(N\) brane data. A physical large-\(N\) theorem separately
requires BV states, topology, trace normalization, connected cumulants,
planar scaling, QME control, and convergence of correlation
distributions. **Counter**: every large-\(N\) claim names whether it is
degreewise stable algebra, admissible Koszul deformation, or physical
BV/correlator limit.

**AP2180 -- Naming a mixed HT habitat does not construct the category (High).**
Wrong: treating labels such as \(\mathrm{Coeff}_{q_-,q_+}\),
\(\mathrm{BV}^{w}_{HT}\), \(\mathrm{FactCenter}^{pp}\),
\(\mathrm{Koszul}_{HT}\), \(\mathrm{BVFactState}_{1/N}\), or
\(\mathrm{GlobDesc}\) as established mathematical categories merely
because the intended objects are named. Correct: each habitat needs its
objects, morphisms, topology/completion, compositions, signs, functorial
rules, obstruction complexes, and universal properties written and
proved. **Counter**: every habitat theorem begins with definitions and
first nontrivial coordinate checks before any equivalence or realization
claim.

<!-- END VOL1_IMPORTED_ANTIPATTERN_REGISTRY_2026_04_30 -->

---

# Topological-Strings 04:11 Final-Frontier Control Surface, 2026-04-30

These entries are Vol II-visible control material for
`~/topological-strings` after the 04:11:11 SAST, 923-page critique
refresh.  They are attack-heal rules, not theorem proofs.

**AP2181 -- One-theorem collapse hides the required six-habitat realization (High).**
Wrong: trying to prove one maximal theorem directly identifying the
formal local CE/PV--Moyal skeleton with a renormalized mixed
holomorphic-topological string theory. Correct: the theorem target is a
six-habitat realization problem, with coefficient, CE/PV, Matlis/MC,
factorization-center, BV/QME, mixed-Koszul, large-\(N\), and global
descent interfaces separated by explicit source/target data and
obstruction classes. **Counter**: every repair states which habitat it
constructs and which obstruction remains.

**AP2182 -- Weighted coefficient scales are mathematical data, not notation (High).**
Wrong: writing \(H_\pm,P_\pm,C_{-,+}\) as decorative symbols after the
strict product/direct-sum endpoint has failed. Correct: a coefficient
theorem must define the two-weight scale \(1<q_-<q_+\), strong duals,
projective tensor topology, bracket continuity \(H_+\otimes H_+\to H_-\),
coadjoint continuity \(H_+\otimes P_-\to P_+\), and convergence of
\(C_{-,+}=\sum_I H_I\otimes\rho_I\). **Counter**: no Costello kernel
claim without the estimates.

**AP2183 -- A formal diagonal is not a convergent MC kernel (High).**
Wrong: treating \(\Theta=\sum H_I\otimes\theta_I+\sum\eta_I\otimes O_I\)
as a valid universal MC element before specifying a completed
convolution dg Lie algebra. Correct: \(\Theta_{-,+}\) must converge in
the weighted convolution habitat and the Hamiltonian and cotangent
components of \(d\Theta+\frac12[\Theta,\Theta]\) must cancel with the
actual signs. **Counter**: every MC claim names the completion and gives
the component cancellation.

**AP2184 -- \(SC^{ch,top}\) is a two-colour convolution habitat, not a product slogan (High).**
Wrong: importing Vol II language as if the mixed Swiss-cheese object were
a product of an open and a closed theory. Correct: the relevant grammar
is a two-coloured chiral-topological operad, its bar cooperad, a
convolution \(L_\infty\) algebra, a universal MC element, and a derived
centre/bulk-boundary universal property. **Counter**: `Koszul_HT`
sections build the two-colour convolution data before using the name.

**AP2185 -- BV/QME habitat names do not construct Costello quantization (High).**
Wrong: saying "Costello BV theory" while omitting fields, BV pairing,
gauge fixing, heat kernel, weighted propagator, counterterms, RG flow,
and local obstruction complex. Correct: a BV theorem specifies the field
sheaf on \(\mathbb R^2_{top}\times\mathbb C^2_{hol}\), residue pairing,
brane support, \(P^{cl}_{\epsilon,L}=P^X_{\epsilon,L}\otimes C_{-,+}\),
renormalized action, and classes
\(Ob_{n,\mathrm{nsc}}\in H^1(Loc_{bulk+defect},Q+\{I_0,-\})\).
**Counter**: no quantum theorem without the denumerable obstruction
complex or a proof that its relevant classes vanish.

**AP2186 -- Graph normalization is not analytic graph realization (High).**
Wrong: deriving the Moyal coefficients from an assumed binary boundary
contraction and calling this the Costello graph theorem. Correct:
normalization says that the contraction
\(\frac{\hbar}{2}(\partial_{z_1}\otimes\partial_{z_2}
-\partial_{z_2}\otimes\partial_{z_1})\) produces the Moyal product.
Realization proves that the Costello propagator, brane kernel,
counterterms, and QME solution produce that contraction and kill higher
connected binary Hamiltonian cumulants. **Counter**: every graph theorem
says whether it is normalization or realization.

**AP2187 -- Physical large \(N\) requires states, cumulants, topology, and Ward identities (High).**
Wrong: upgrading degreewise stable trace invariant theory or the
admissible Koszul enveloping algebra into physical open/closed
large-\(N\) duality. Correct: physical large \(N\) is a sequence
\((Obs_N,Q_N,\Delta_N,I_N,\mu_N,\omega_N,\nu_N,T_N)\) with states,
trace normalization, connected cumulants, planar expansion, topology,
convergence, and QME Ward identities. **Counter**: every large-\(N\)
claim declares algebraic stable, admissible Koszul, or physical
BV/correlator status.

**AP2188 -- Global descent needs a target datum and null-homotopies (High).**
Wrong: transferring the local formal-disk theorem to compact CY, BCOV,
MNOP, Vol III, Igusa/BKM, or sister-volume conclusions without a global
target. Correct: supply
\(C_{tar}=(X,\omega,\Omega,U,E_i,Q_i,\omega_i,I_i,\phi_{ij},\phi_{ijk},
Per,Anom,QME,Conv,Hyp,Op)\), the obstruction vector
\(Ob_{UKD}(C_{tar})\), and chosen null-homotopies. **Counter**: no
target datum means no transfer.

**AP2189 -- Manuscript-proper reconstitution must use the Vol II rectify harness (High).**
Wrong: rewriting the topological-strings manuscript proper by ad hoc
local patching or by the older Vol I/local rewrite habit after the 04:11
refresh. Correct: writing, rewriting, or reconstituting manuscript text
uses `~/chiral-bar-cobar-vol2/.agents/skills/chriss-ginzburg-rectify/SKILL.md`
under Supremum Discipline, with the Vol II control surfaces loaded.
**Counter**: every manuscript-proper edit records the Vol II
rectify/Supremum path or remains a non-inscriptive report.

**AP2190 -- \(\mathbf H_{\Delta_5}\) in Vol II is a gated K3 comparator, not native Vol II output (High).**
Wrong: treating the Vol III K3 chiral Hall--Drinfeld double
\(\mathbf H_{\Delta_5}\) as something produced by Vol II's generic
\(A_\infty\)-chiral, Hochschild, \(W\)-algebra, spectral-braiding, or
HT bulk-boundary machinery. Correct: Vol II may use
\(\mathbf H_{\Delta_5}\) only as an external K3 boundary fibre,
recognition target, deformation-chart comparator, or obstruction test
after the finite Hall/CoHA source, compact support/orientation,
positive half, opposite half, Cartan, completed coproduct,
nondegenerate Hopf pairing, parity/root fixture, PBW/no-extra-relations,
radical quotient, inverse-limit/completion, and Borcherds--Heegner
comparison gates have been stated. The automorphic datum
\(\Delta_5\) supplies the Borcherds denominator target; it does not
construct the chiral bialgebra. **Counter**: every Vol II occurrence of
\(\mathbf H_{\Delta_5}\) declares comparator/conditional K3 boundary
status or is removed from theorem-source position.

**AP2191 -- Failed theorem surfaces are not allowed to remain loose conditionals (High).**
Wrong: when a theorem, lemma, or construction breaks under attack,
lowering it to "conditional", "expected", "external input", or
"future work" and then treating the manuscript as healed. Correct:
Supremum Discipline requires reworking the failed surface until it is
one of two precise objects: a proved theorem with the missing habitat,
map, kernel, homotopy, counterterm, sign, finite-window row, incidence
table, or source fixture supplied; or a proved obstruction theorem
naming the exact cohomology class, cokernel functional, Roos class,
matrix equation, graph row, or source-face table that blocks the
original claim. **Counter**: every attack-heal report on theorem-level
material either inscribes the stronger true theorem or names the exact
obstruction theorem and the next construction target.
