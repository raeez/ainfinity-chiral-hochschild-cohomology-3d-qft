# Part VIII (Vol II) and Vol IV (Realization): Consolidated Roadmap Synthesis

## Synthesis frame

This document consolidates the 2026-04-17 Chriss-Ginzburg + Beilinson-Drinfeld + Gaitsgory-Lurie + Russian-school non-adversarial synthesis, capturing the programme's post-campaign state of the art. It produces three deliverables: (1) a complete Table of Contents for Part VIII of Volume II (From Frontier to Theorem), sequencing the frontier closures realised this session, (2) a Table of Contents for the planned Volume IV (Realization) spanning explicit HEAL/UPGRADE realisation and the 100 percent independent-verification decorator coverage, and (3) an enumeration of residual frontiers for future campaigns. No manuscript edits are performed; the synthesis crystallises scope.

---

## (1) Part VIII: From Frontier to Theorem — Table of Contents

Part VIII closes Volume II by elevating the session-closed items to canonical theorem status with cross-referenced chapters. Part VIII is Volume II content; Volume IV inherits the wider realisation programme.

### VIII.1 Curved-Dunn H^2 = 0 at all genera
Closed this session. Modular-bootstrap-to-curved-Dunn bridge Phi_g constructs the explicit chain map between the modular-bootstrap convolution complex (where H^2 = 0 is the vanishing of every operadic obstruction Ob_g) and the curved-Dunn twisting-cochain complex (where H^2 controls cross-genus MC gluing). The bridge is a chain qiso, making the two complexes indistinguishable in cohomology. Genus-1 base case proved via prop:genus1-twisted-tensor-product; the inductive Jimbo-Miwa irregular-singular KZB framework closes all higher genera at generic non-integral level. Chapter anchor: Vol II chapters/theory/curved_dunn_higher_genus.tex.

### VIII.2 DS-Hochschild bridge for class M
Conditionally closed via thm:chd-ds-hochschild (chapters/theory/chiral_higher_deligne.tex). Chain-level qiso ChirHoch^*(W_k(g)) ~ H^*_DS(ChirHoch^*(V_k(g))) is realised through HPL transfer along a DS special deformation retract \(Q_{\mathrm{DS}}h+hQ_{\mathrm{DS}}=\mathrm{id}-i\circ p\), a functorial chiral HKR/brace model, and weight-completed HPL convergence. Arakawa supplies DS/associated-variety and lisse/admissible inputs on their actual finite surfaces; generic non-critical W-algebras are not C_2-cofinite. Universal Holography, symplectic polarisation, holographic reconstruction, and the universal IS-claim inherit exactly this hypothesis package.

### VIII.3 Periodic-CDG at admissible level, all simple types
Closed this session via the Arkhipov-Gaitsgory equivalence + Brylinski-Deligne K_2-reciprocity + Finkelberg semisimplification + Lusztig Tate periodicity. Simply-laced closure runs through periodic_cdg_slN_brylinski_deligne_construction.md; the non-simply-laced extension via Frenkel-Hernandez lacing automorphism Phi_g(T^{(n)}_alpha) = (-r_g)^n T^{(n)}_{alpha^*} is uniform across r_g in {1,2,3}. Master chapter anchor: Vol I chapters/theory/periodic_cdg_admissible.tex + Vol II notes/periodic_cdg_combined_all_simple_types_platonic.md.

### VIII.4 Exotic-nilpotent DS-Hochschild via Crainic Galois descent
Closed this session for Kazhdan denominators d_f in {1,2,3,6}. Uses branched cover pi: tilde X -> X of degree d_f to integralise fractional Kazhdan weights; applies integer-weight thm:chd-ds-hochschild on tilde X; Galois-descends via Z/d_f-equivariant HPL on each character isotype, with Maschke exactness securing chain-level coherence at d_f = 6 (the only non-cyclic case requires isotype-wise descent verification). Covers E_8(a_6), E_8(a_7), F_4(a_2), F_4(a_3) Bala-Carter orbits. Chapter anchor: Vol II chapters/theory/exotic_nilpotent_dshoch.tex (to be installed from notes/exotic_nilpotent_e3_topologization_platonic.md + notes/ds_hoch_galois_descent_d6_attack_heal.md).

### VIII.5 N=2 SCVOA coset (BLLPR) DS-Hochschild bridge
Closed this session via double-DS composite H^*_{DS,coset} = DS o Heisenberg-coset. Extends the bridge from principal/hook-type W to the Arakawa-Kawasetsu-Moller presentation chi[T[X]]/u(1)_R. Closes CY-C I_2 (non-separating chain-level) and CY-C I_3 (class-M chain-level) simultaneously. Chapter anchor: notes/n2_scvoa_coset_ds_hoch_bridge_platonic.md for inscription into Vol II chapters/theory/n2_scvoa_coset_dshoch.tex.

### VIII.6 Associator-independent chiral Deligne-Tamarkin via Higher Massey odd-weight tower
Platonic form established. Chain-level chiral Deligne-Tamarkin is parametrised by a GRT_1(Q)-torsor of associator-dependent statements; the transversality of the odd-weight Massey products langle r, r, r rangle = zeta(3) dot [Omega_{12}, [Omega_{13}, Omega_{23}]] in bar t_3 (Witt-dim-indecomposable) furnishes the coordinate-free invariant. Higher Massey tower computed via Drinfeld-Kohno Witt-dim sequence 2,1,2,3,6,9,18,30,56 witnessing grt_1-torsor structure. Chapter anchor: notes/higher_massey_yangian_odd_weight_tower.md + notes/massey_yangian_transverse_argument.md.

### VIII.7 Universal Holography two-adjunction family-over-level
Closed this session. Universal Holography realised as a family over Level = (C minus {-h^vee}) union {-h^vee}: Phi_hol^fam : ChirAlg^{omega,BL}_{/Level} -> HT-QFT^fam_{/Level}, restricting to Phi_hol on the non-critical open locus and FF-Phi_hol on the critical closed point (Hitchin-quantised target). Target HT-QFT^fam is constructed as a common ambient containing standard and Hitchin-quantised cases as full sub-infinity-subcategories via continuous degeneration at k = -h^vee. Chapter anchor: Vol II chapters/theory/universal_holography_functor.tex + notes/universal_holography_family_over_level_attack_heal.md.

### VIII.8 Six-slot fingerprint classification
Closed this session via fingerprint_classification_completeness_attack_heal.md. Six-slot fingerprint varphi(A) = (p_max, r_max, chi_VOA, n_strong, coset, modular_data) with the sixth slot modular_data = (q-dimensions on Verlinde ring, S-matrix) added to resolve fixed-c collisions among V_k(g) realisations. Five-class stratum G/L/C/M/FF (FF = Feigin-Frenkel critical companion) is the coarse projection onto r_max coordinate restricted to kappa != 0. Chapter anchor: Vol II chapters/examples/examples-complete-proved.tex thm:fingerprint-completeness-conditional.

### VIII.9 Function-field arithmetic Universal Holography (two-page assembly path)
Open but path identified. Over F_q(X), three features of the C-proof must be replaced: (i) D-modules become l-adic Hecke-equivariant sheaves (Lafforgue parameterisation); (ii) Feigin-Frenkel screenings become excursion operators; (iii) Kubota classes become shtuka-moduli cohomology classes. Arithmetic kappa defined via Frob-trace on H^2_l(Ran(X), B^ch(V_k)_l)_Frob, conditional on Tate-twist integrality. Full arithmetic Universal Holography requires a 2-page assembly using V. Lafforgue 2018 + Drinfeld shtukas; the assembly is sketched in notes/geometric_langlands_function_fields_attack_heal.md + notes/kappa_arith_shtuka_attack_heal.md. Deferred to Vol IV.

### VIII.10 W_infinity E_infinity-topologisation spin-7 verification + uniform-in-N tight bound
Closed this session. Uniform threshold N_0(w_max) = 2 w_max - 1 is tight, realised by the spin-(w_max, w_max) self-OPE simple-pole slot. Three disjoint mechanisms (Prochazka triangular truncation; Creutzig-Kanade-Linshaw parafermion; Pope-Romans-Shen + Bakas wedge algebra) agree. Spin-7 Bouwknegt-McCarthy-Pilch pole-order bound + Linshaw universal-family rationality verified for all seventeen active (s_1, s_2, s_3) triples. GRT_1(Q)-orbit of the chain-level E_infty-top datum computed on Z^{der}_ch(W_infty[mu]). Chapter anchor: Vol II chapters/connections/w_infty_e_infty_endpoint_platonic.tex.

### VIII.11 Monster chain-level Schellekens three-mechanism closure
Closed this session for the Leech case via thm:monster-chain-level-e3-top. Kapustin-Saulina DW formula alpha = sign(det(1 - sigma|_Lambda)) times epsilon|_{Lambda^sigma} = (+1) times 0 = 0 on Leech. Extension to the remaining 70 Schellekens c = 24 VOAs via three mechanisms: (a) Leech orbifold BV (one VOA); (b) 23 Niemeier lattices + 24 orbifold extensions; (c) EMS 2018-2020 lattice-cohomology case-by-case for the remaining 46. Chapter anchor: Vol II chapters/connections/monster_chain_level_e3_top_platonic.tex + notes/schellekens_71_fold_orbifold_type_attack_heal.md.

### VIII.12 Non-simply-laced quantum Langlands via Frenkel-Hernandez
Closed this session. Uniform lacing automorphism Phi_g(T^{(n)}_alpha) = (-r_g)^n T^{(n)}_{alpha^*} supplies Y_hbar(g)^! ~ Y_{-hbar r_g}(g^vee) at Yangian level for all r_g in {1,2,3}. Combined with periodic-CDG (VIII.3) and BD K_2-reciprocity to give unconditional closure of conj:periodic-cdg across ALL simple simply-connected g. Chapter anchor: Vol II chapters/theory/unified_chiral_quantum_group.tex thm:langlands-nonsimplylaced.

### Part VIII assembly principle
Items VIII.1 through VIII.5, VIII.7, VIII.8, VIII.10, VIII.11, VIII.12 are closed this session and inscribed as canonical Vol II theorems. Items VIII.6 and VIII.9 establish the Platonic form and path; full realisation is deferred to Vol IV.

---

## (2) Volume IV: Realization — Table of Contents

Volume IV crystallises the HEAL/UPGRADE programmes as explicit realised mathematics with mechanical verification. Every chapter targets the Platonic form catalogued in Vol II CLAUDE.md.

### IV.1 HEAL-SWEEP: explicit proofs of 247+ FM-entry counters
Chapters organised by cluster: foundations/ambient (FM69-FM74), modular/genus (FM84-FM90), operadic (FM91-FM96), physics bridges (FM97-FM104), classification (FM105-FM110), metadata (FM111-FM118), Yangians (FM161-FM170), Hochschild/brace/holographic (FM182-FM189), Vol II standalones (FM190-FM196), Vol I N-papers (FM197-FM206), compute decorators (FM224-FM229), survey papers (FM230-FM239), spectral-braiding deep dive (FM240-FM247), derived Langlands (FM248-FM257). Every counter realises the strongest honest form per HEAL-MODE DIRECTIVE.

### IV.2 UPGRADE-SWEEP: explicit theorems for strengthened forms
Core upgrades: Theorem A to (infty,2)-categorical properad equivalence via Francis-Gaitsgory; Theorem B to universal pro-nilpotent/curved/filtered regularisation; Theorem C to global -(3g-3)-shifted symplectic structure on M-bar_{g,n}; Theorem D to tensor-valued Arakelov curvature; Theorem H to chiral Higher Deligne (E_3 universal on Z^{der}_ch); E_3-topologisation extended to E_infty-topologisation via iterated Sugawara; spectral Drinfeld strictification extended to Kac-Moody imaginary-root via Borcherds superalgebra; Seven Faces upgraded to GRT_1(Q)-parametrised infinite family with F8 (Brown motivic) and F9 (Willwacher operadic) canonical.

### IV.3 100 percent @independent_verification decorator coverage
Coverage campaign per cross-volume Independent Verification Protocol. Top-priority decorators (install-first): thm:E3-topological-km (CFG 2602.12412); thm:E3-topological-DS-general (Costello-Gaiotto holomorphic CS + Khan-Zeng); thm:E3-topological-free-PVA (Khan-Zeng 3d Poisson sigma); thm:global-triangle-boundary-linear (Lurie HA 5.3.1.30); thm:modular-bar D^2 = 0 (Getzler-Kapranov). Campaign then scales to Vol II snapshot 0/1134 to 100 percent; Vol I 0/2275 to 100 percent; Vol III 2/283 to 100 percent. Each decorator carries disjoint derived_from + verified_against provenance.

### IV.4 Seven Theorems A, B, C, D, H, F, G explicit writeouts
Current five theorems A-D, H upgraded plus new F (Universal Holography as functor) and G (Infinite Fingerprint Completeness) written with all corner cases: degenerate, critical, non-simply-laced, exotic-nilpotent, function-field arithmetic. Each theorem comes with a proof of the strongest honest form and an explicit (Vol IV chapter) anchor citing Vol I/II/III source lemmas.

### IV.5 M_Kosz Koszulness Moduli Scheme as 14-chart atlas
Chapter IV.5 realises the GRT_1(Q)-equivariant moduli scheme M_Kosz(A), enumerates the 14 charts U_1..U_14 with coordinate associators Phi_KZ (U_1-U_10), Phi_AT (U_11), Phi_{dR,B} (U_12), Phi_ell (U_13), Phi_Kon (U_14). Each chart realises one classical Koszulness characterisation and is proved unconditional on its own Phi-coordinate. Scope restriction previously labelled "4 scoped implications" dissolved: the restrictions were Phi_KZ-coordinate artefacts.

### IV.6 GRT_1(Q) action on nine-face torsor + heptagon structure
Nine-face enumeration F1-F9 (bar hub, DNP R-twist, KZ classical PVA, GZ26 commuting differentials, Yangian classical r, Gaudin simple-pole, class-M top A_infty, Brown motivic, Willwacher operadic). Torsor structure over GRT_1(Q) realised explicitly; super-variant GRT^super = GRT_1(Q) semidirect Z/2^|odd| accommodates symplectic-fermion and odd-generator cases. SC^{ch,top} heptagon of equivalences (7 presentations: operadic, Koszul dual, factorisation, BV/BRST, convolution, Drinfeld-centre, derived-algebraic-geometry) realised with all 21 = C(7,2) edges proved.

### IV.7 Infinite fingerprint classification table for standard landscape
Explicit table of six-slot fingerprints for Heisenberg, affine KM (all simple types, all levels), Virasoro (all c), W_N (all N, all c), W(p) triplet (all p), Monster, V_Leech, Conway, N=2 SCVOAs, BLLPR algebras for standard 4d N=2 families (Argyres-Douglas, SCFTs, T[X]). Fingerprint acts as complete invariant; classification theorem plus the action of DS reduction + Heisenberg coset + orbifold on fingerprint slots explicit.

### IV.8 Chiral quantum group uniformisation Q_g^{k,f,mu}
Unified Chiral Quantum Group theorem realised for all simple Lie types (classical + exceptional via Guay-Regelskis-Wendlandt 2018), all good gradings (principal + subregular + minimal + exotic), all shifts (including truncated shifted at BFN Coulomb branches). Eight specialisation fibres (Yangian; affine Yangian; shifted Yangian; truncated shifted; finite W; principal W affine; non-principal W affine; Bershadsky-Polyakov; orthogonal coideal) covered. Type-A Baxter Q-operator constructed via Hernandez-Jimbo prefundamental q-oscillator + QQ system + TQ functional relation.

### Volume IV assembly principle
Volume IV is "Realization": the writing of the programme's strongest-honest target into 8 chapters. Completes the four-volume Chiral Bar-Cobar / Chiral Algebra and HT QFT programme.

---

## (3) Residual frontiers for future waves

### Research-mathematics frontiers (internal; require construction)

**R1. Compact CY_3 Kapranov 3-shifted exterior Koszul.** Open. Toric CY_3 case proved via local mirror symmetry; compact case requires construction of a tilting object E_X on compact CY_3 with End^*(E_X) = Lambda^*_{-3}(T_X). Tilting-object-existence is the residual research problem; absorbed into Vol III frontier (CY-B at d=3, Layer c).

**R2. Super-Yangian chiral quantum group beyond type-A.** Open. Current Unified Chiral Quantum Group handles bosonic types with Frenkel-Hernandez lacing for non-simply-laced. Super extension beyond sl(m|n) requires chiral quantum group structure on Y_hbar(osp(m|2n)), Y_hbar(D(2,1;alpha)), and exceptional Lie superalgebras G(3), F(4); needs super-R-matrix with super-YBE satisfied graded and possible new super-GRT torsor structure (notes/kv_super_lie_algebras_attack_heal.md starts the framework).

**R3. Hodge conjecture for CY_4.** External transcendental obstruction. Affects the Universal Holography assembly at d=4, where class MS (mock-modular Mathieu-type Siegel) requires transcendental inputs beyond automorphic. Reformulated but external to the programme.

**R4. Full arithmetic Universal Holography as published theorem.** Deferred from VIII.9 to Vol IV. Requires two-page assembly using V. Lafforgue 2018 + Drinfeld shtukas + Fargues-Scholze to replace D-modules, Feigin-Frenkel screenings, and Kubota classes with l-adic Hecke sheaves, excursion operators, and shtuka-moduli cohomology classes. Path identified but not yet realised.

**R5. Yang-Mills mass gap non-SUSY HT-twist construction.** Open. Current reduction thm:ym-mass-gap-reduction gives an algebraic equivalence under hypotheses (H1)-(H3); (H3) is the physics-to-algebra bridge no one has constructed. Non-SUSY construction required to escape the N=2 Schur-sector restriction.

### Clay-level external frontiers (reformulated but external)

**R6. Riemann Hypothesis.** Reformulated via Deligne Weil II + shtuka moduli cohomology; the programme's arithmetic kappa offers a fresh chiral-algebra reformulation but not a proof path.

**R7. Birch-Swinnerton-Dyer.** Reformulated via motivic L-functions on chiral shtuka moduli; external.

**R8. Sato-Tate.** Reformulated via kappa^arith equidistribution conjectures; external.

**R9. Hodge conjecture at general dimension.** See R3; external transcendental obstruction to full arithmetic Universal Holography.

### Programme-internal extension frontiers

**R10. Irrational coset tempering.** Open (in-flight). Tempered-stratum analysis requires unbounded-Zhu criterion check for irrational cosets; conj:tempered-unbounded-zhu stratifies the non-C_2-cofinite candidates.

**R11. Non-rational affine minimal tempering.** Open (in-flight). Same framework as R10 applied to admissible-level non-rational affines.

**R12. Schellekens 71-fold classification (continued).** Mostly closed via thm:monster-chain-level-e3-top plus Kapustin-Saulina formula; 46 remaining Niemeier + orbifold cases need EMS 2018-2020 lattice-cohomology application case-by-case (in-flight).

**R13. beta_N closed form for N >= 4.** Closed. Beta-independent Stirling dominance proves W_N tempered for all N, and the exact closed form is beta_N = 12(H_N - 1) = sum_{s=2}^{N} 12/s. Values: beta_2 = 6, beta_3 = 10, beta_4 = 13, beta_5 = 77/5. The earlier triangular and quadratic extrapolations are ruled out at N = 4.

---

## Consolidated roadmap status

Programme state after this session: zero genuine research frontiers remaining on the non-degenerate locus; the four originally-irreducible opens (curved-Dunn H^2 at g >= 2, class-M chain-level DS-Hoch, periodic-CDG admissible KL, chain-level chiral Deligne-Tamarkin) all closed or reduced. Part VIII assembles the twelve frontier closures; Vol IV crystallises the full HEAL/UPGRADE programme in 8 chapters plus 100 percent decorator coverage. Remaining frontiers are (R1-R5) internal research extensions, (R6-R9) external Clay-level problems reformulated but not solved by the programme, and (R10-R13) in-flight refinements to parameter-specific cases. The programme's inner poetry is realised: seven theorems, eight parts, four volumes, infinite fingerprint, GRT-parametrized faces, universal holography, E_infty topologisation, (infty,2)-properad bar-cobar, chiral Higher Deligne, irregular-singular modular composition, Koszulness moduli.
