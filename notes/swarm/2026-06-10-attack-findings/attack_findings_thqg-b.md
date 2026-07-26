# THQG-B attack findings — Wed Jun 10 14:46:09 WIB 2026

## F1: Virasoro shadow table — false numeric column + sign-inconsistent symbolic rows
- severity: major; verdict: false (table entries), closed form itself verified correct
- anchors: chapters/connections/thqg_holographic_reconstruction.tex:1550-1591 (table), :1422-1438 (recursion+closed form), :1449-1486 (proof)
- evidence: closed form S_r^cub = (-1)^(r-4) 6^(r-5) 240/(c^(r-3)(5c+22) r) VERIFIED consistent with recursion S_{r+1}=-6r/(c(r+1))S_r from S_5=-48/(c^2(5c+22)). Exact values at c=26 (5c+22=152):
  r=6: 15/166972 = 8.98e-5  (table says 9.01e-5 — false 3rd digit)
  r=7: -135/7597226 = -1.78e-5 (table says -1.48e-5 — false, ratio 0.83)
  r=8: +405/112873072 = +3.59e-6 (table says +2.14e-6 — false, ratio 0.60); row-8 SYMBOLIC entry (-6)^3·240/(c^5(5c+22)·8) is NEGATIVE — wrong sign vs closed form (+) and vs its own numeric (+)
  r=9: -135/183418742 = -7.36e-7 (table says -2.77e-7 — false, ratio 0.38); row-9 symbolic (-6)^4 form is POSITIVE — wrong sign
  r=10: +729/4768887292 = +1.53e-7 (table says +3.26e-8 — false, ratio 0.21); row-10 symbolic (-6)^5 form is NEGATIVE — wrong sign
  Mechanism of numeric corruption ~ recursion factor frozen at r=5 (-30/(c(r+1))) reproduces -1.481e-5, 2.136e-6 for rows 7,8. Rows 8-10 symbolic use (-6)^(r-5) where closed form requires (-1)^(r-4)6^(r-5) = -(-6)^(r-5).
- heal: replace symbolic rows 8-10 by (-1)^(r-4)6^(r-5)·240/(c^(r-3)(5c+22)r); replace numeric column r>=6 by 8.98e-5, -1.78e-5, +3.59e-6, -7.36e-7, +1.53e-7.

## F2: Jet-dimension internal contradiction: r(r+1)/2 vs r(r+3)/2 (the N(N+1)/2 lead)
- severity: major; verdict: inconsistent
- anchors: thqg_holographic_reconstruction.tex:1257-1259 + table :1263-1286 (uses 1+r+r(r+1)/2: 6,10,15,21,28,36,45) VS :1352-1357 (proof derives 2r+C(r,2)=r(r+3)/2 using p(2)=2 states at weight 4)
- evidence: excess-2 extended jet dim: Computation says r + C(r,2) = r(r+1)/2 (1 state at weight 4 per factor); growth-theorem proof says 2r + C(r,2) = r(r+3)/2 (2 states at weight 4: correct, vacuum module has L_{-4}|0>, L_{-2}^2|0>). Both cannot hold. Proof count is correct given "descendants and quasi-primary composites" (line 1249-1251). Also proof names the weight-4 states as "L_{-2}|0> and L_{-1}^2|0>" (line 1353) — both have weight 2, and L_{-1}|0>=0 in the vacuum module; correct states are L_{-4}|0>, L_{-2}^2|0>.
- heal: table column should read 1+r+r(r+3)/2: 8,13,19,26,34,43,51 for r=2..8; fix state labels.

## F3: "states at weight h = p(h-2)" false for h>=5; extended-dimension formula overcounts
- severity: major; verdict: false
- anchors: thqg_holographic_reconstruction.tex:1239-1241, :1363-1372 (eq Π p(λ_i)), :1267-1286 (Σ_{d=0}^{2r} column)
- evidence: vacuum-module graded dim at weight h = #partitions of h into parts>=2: h=5: 2 (claimed p(3)=3); h=6: 4 (claimed 5); h=7: 4 (claimed 7); h=8: 7 (claimed 11). Computed by sympy. Hence dim J^r_ext(d) = Σ_{λ1+..+λr=d} Π p(λ_i) overcounts from excess d=3 on. O(r^D) conclusion survives (leading term from d factors of excess 1) but constant c_d and the formula are wrong; also c_d defined as "Σ_{|λ|=d} Π p(λ_i)/d!" conflates partitions with compositions. The Σ_{d=0}^{2r} column values 7,15,26,42,63,92,131 reproduce under NO natural counting (neither p, vacuum-module, nor L_{-1}-only); no formula given.
- heal: replace p(h-2) by partitions-into-parts>=2 count (equivalently P(q)=Π_{n>=2}(1-q^n)^{-1}); recompute table; delete or derive the Σ column.

## F4: Kac-table parametrization internally mismatched
- severity: minor; verdict: false (remark-level)
- anchors: thqg_holographic_reconstruction.tex:1303-1305
- evidence: pairing c = 13-6(t+1/t) with h_{r,s} = ((rt-s/t)^2-(t-1/t)^2)/4 fails: t=2 gives c=-2 and h_{1,2}=-5/16, but true h_{1,2}(c=-2)=-1/8 (symplectic fermion). The h-formula as written requires c = 13-6(t^2+t^{-2}) (verified at t=sqrt(2)).
- heal: use c = 13-6(t+1/t), h_{r,s} = ((rt-s)^2-(t-1)^2)/(4t).

## S1 (sound): closed form S_r^cub and recursion mutually consistent; Step2-5 arithmetic verified (480/(c^2(5c+22)), -2880/(c^3(5c+22)), -6r/(c(r+1)) factor). Weight-4 Gram det c^2(5c+22)/2 verified by direct Virasoro computation.

## F5: D3/N=4-SYM kappa: proposition vs summary table disagree by (N-1)/(2N)
- severity: major; verdict: inconsistent
- anchors: thqg_holographic_reconstruction.tex:2373-2375 (prop:twisted-n4-datum: kappa(gl_N,1) = (N^2-1)(N+1)/(2N) + 1) vs ~:2906 WIP table row "D3 (N=4 SYM): N(N+1)/2"; Vol I convention anchor: chiral-bar-cobar/chapters/examples/landscape_census.tex:2414 (kappa(g_k) = dim g (k+h^v)/(2h^v)), :1369 (kappa = rank for lattice/Heisenberg)
- evidence: (N^2-1)(N+1)/(2N) + 1 - N(N+1)/2 = (N-1)/(2N) != 0 for N>=2 (N=2: 13/4 vs 3). Prop value = Vol I census convention (sl_N part dim(N^2-1), h^v=N, k=1, plus Heisenberg u(1) kappa=1). Table value = treating gl_N as simple with dim N^2, h^v=N. Same file asserts both.
- heal: table row should read (N^2-1)(N+1)/(2N) + 1 (equivalently (N^3+N^2+N-1)/(2N)); or declare the u(1) convention and use it consistently.

## F6: Class-M bulk corollary invokes "Arakawa C_2-cofiniteness" outside its scope
- severity: major; verdict: unlicensed (inapplicable citation in proof body)
- anchors: thqg_holographic_reconstruction.tex:3161-3207 (cor thm:class-M-chain-bulk proof, :3196-3197); bibitem main.tex:3716 (Arakawa15 = Rationality of W-algebras, Ann. Math 182 (2015))
- evidence: corollary states A = W^k(g) UNIVERSAL principal W-algebra at arbitrary non-critical level; proof asserts HKR degeneration "cohomologically by Arakawa C_2-cofiniteness". Universal W^k(g) is NEVER C_2-cofinite: R_V = W^k/C_2 = polynomial ring on rank(g) generators, infinite-dimensional. Arakawa's C_2-cofiniteness (IMRN 2015, arXiv:1004.1554) and rationality (Annals 2015, the cited item) apply to SIMPLE quotients W_k(g) at admissible (nondegenerate) levels only. Also Arakawa15 bibitem gives arXiv:1505.09016; the real arXiv ID of that Annals paper is 1211.7124.
- heal: replace mechanism with Arakawa's associated-variety theorem X_{W^k(g)} = Slodowy slice / arc-space freeness of gr W^k(g) (Feigin-Frenkel grading), which is what gives HKR degeneration for the universal algebra; fix arXiv ID; or restrict corollary to C_2-cofinite simple quotients.

## F7: "Costello--Gaiotto identify the DS boundary condition for holomorphic Chern-Simons" - dubious attribution to 1812.09257
- severity: major; verdict: unlicensed (suspected misattribution; cannot fetch paper offline)
- anchors: thqg_holographic_reconstruction.tex:1834-1836, 1895-1899, 3099-3100 (cite [CostelloGaiotto18 S4]), 3136-3139 (cite [S5]: "DS ... commutes with bulk BV quantisation"); bibitem main.tex:2784 (Twisted holography, arXiv:1812.09257)
- evidence: 1812.09257 constructs the B-model on the deformed conifold dual to the N D1-D5... chiral algebra of gl_K type; its sections do not, to the auditor's knowledge, contain a principal Drinfeld-Sokolov boundary condition for holomorphic Chern-Simons nor a theorem that DS commutes with bulk BV quantisation. W-algebra boundary constructions live in Gaiotto-Rapcak (corner VOAs, 1703.00982) and Costello-Yamazaki-type 4d/5d CS. Load-bearing: this is hypothesis (H1) class-M coverage of Universal Holography.
- heal: either cite the actual construction (Gaiotto-Rapcak corner; or Costello-Gaiotto VOA/3d N=4 1804.06460 if intended) with verified section, or downgrade (H1) DS sector to a named hypothesis.

## F8: Tautology cluster tagged ProvedHere (proof = restatement of hypothesis)
- severity: critical; verdict: unproved (content deferred to assumption)
- anchors: thqg_holographic_reconstruction.tex:3343-3383 (thm:chiral-hochschild-trinity), :3235-3276 (thm:hc-verdier-distance), :3295-3318 (BRST lower bound), :3385-3407 (kappa-conductor trace)
- evidence: (a) Trinity theorem: hypothesis "Suppose Phi_GA, Phi_AB are filtered quasi-isomorphisms" and conclusion is the boxed zigzag of those same quasi-isomorphisms; proof says "the assertion is the two-out-of-three statement" - but both maps are hypothesised, nothing is derived; construction of Phi_GA/Phi_AB (the actual FM/enveloping/factorization comparison) appears nowhere. (b) Verdier distance: hypothesis identifies gr(undetectable ideal) = Rad(ev) weight-preservingly; conclusion d_QECC = d_Verdier is then definitional. (c) BRST bound: hypothesis "sends every nonzero radical class to bar weight >= twice its ghost spin"; conclusion d >= 2*lambda_min is the same inequality. (d) kappa-conductor: hypotheses include the identity K(A) = -c_ghost; conclusion asserts K(A) = tr(K_A) = -c_ghost, but NO step connects tr_{Z(A)}(K_A) to K(A): the proof says the Vol I identity "identifies this trace" with -c_ghost - non sequitur, the Vol I identity as quoted concerns K(A), not the trace.
- heal: retag all four as ClaimStatusConditional with the comparison morphisms as named hypotheses, or supply the actual constructions (for (d): prove tr_{Z(A)}(K_A) = K(A) or add it as hypothesis).

## F9: Verdier-distance calibrations unanchored; [[5,1,3]] vs HaPPY conflation
- severity: minor; verdict: missing-artifact
- anchors: thqg_holographic_reconstruction.tex:3278-3293
- evidence: d_Verdier(Vir_{-22/5}) = 3 cited to "genus-zero computations" with no computation anywhere in volume; the natural Lee-Yang degeneracy is the level-4 vacuum null (5c+22=0 kills the weight-4 Gram det), suggesting 4, not 3, absent a stated weight dictionary. "HaPPY pentagon code: d_Verdier = 3 matching the published QECC distance 3 of the [[5,1,3]] pentagon code" conflates the single-tile seed code [[5,1,3]] (d=3, correct) with the HaPPY holographic code (distance grows with radius).
- heal: supply the genus-zero computation with the weight dictionary, or strike the calibration values; say "single-pentagon seed code".

## S2 (sound): K3xE WIP healing verified correct: HH(K3xE) = (1,2,23,44,23,2,1) (Kunneth convolution (1,0,22,0,1)*(1,2,1), total 96); F_g = 24*lambda_g with lambda_g = (2^(2g-1)-1)|B_2g|/(2^(2g-1)(2g)!): F_1=1, F_2=7/240, F_3=31/40320 all verified; Br(K3xE)_tors = (Q/Z)^(22-rho) verified (old (Q/Z)^22+(Q/Z) was false, Br(E)=0); HH^2 = 21+1+1 = 23 with Poisson/gerby labels now correct; sl2 current bracket f^{ab}_c J^c + k delta^{ab} (level only on central term) now correct; Q_L = (2*24)^2 = 2304 correct.

## F10: Degree-6 proof body: double-counted brackets + false arithmetic step (640/12 = 80/3)
- severity: major; verdict: inconsistent (table correct, proof broken)
- anchors: thqg_gravitational_complexity.tex:800-820 (proof of thm:thqg-virasoro-tower-explicit)
- evidence: proof lists THREE contributions (a) {C,Sh5}, (b) {Sh4,Sh4}, (c) {Sh5,C} and sums without symmetry factor: o^(6) = -640(45c+193)/(c^3(5c+22)^2); then writes Sh_6 = 640(...)/12 = 80(...)/3 -- but 640/12 = 160/3, NOT 80/3 (factor 2). The table value 80(45c+193)/(3c^3(5c+22)^2) is correct under the unordered-pair (or ordered-with-1/2) convention: sympy: unordered+half gives exactly the table; ordered-no-half gives 2x. The phrase "Restoring the standard normalisation gives the stated formula" papers over the factor 2. S7 = -2880(15c+61)/(7c^4(5c+22)^2) VERIFIED exactly under the correct convention (o^(7) = (2*3*6/c)S3S6 + (2*4*5/c)S4S5).
- heal: state the bracket convention o^(r+1) = (1/2) sum_{ordered j+k=r+3} (2jk/c) S_j S_k (or unordered with half on diagonal), delete contribution (c) as separate item, and remove the "restoring normalisation" sentence; arithmetic then closes exactly.

## F11: Growth-bound paragraph false: claims divergence where its own product converges to 0
- severity: major; verdict: false (within ProvedHere proof)
- anchors: thqg_gravitational_complexity.tex:713-730 (proof of thm:thqg-virasoro-infinite); contradicts thqg_holographic_reconstruction.tex:1499-1512 (trichotomy: geometric decay for |c|>6)
- evidence: |S_5| prod_{j=5}^{r-1} 6j/(|c|(j+1)) = |S_5| 6^{r-5} 5/(r|c|^{r-5}) -> 0 for |c| > 6 (e.g. c=26: ~(3/13)^r/r). The text bounds it by C^r r!/|c|^r and concludes "This is O(r!/|c|^r), which diverges as a sequence in r" -- an upper bound diverging implies nothing, and the actual sequence converges to 0 for |c|>6, contradicting the companion file's trichotomy (geometric decay at c=26).
- heal: delete the factorial bound; state |S_r^cub| = 5|S_5| (6/|c|)^{r-5}/r and route divergence discussion to |c|<6 only, as in the trichotomy.

## F12: Inductive step of "tower is infinite" false as stated; counterexample from its own table
- severity: major; verdict: unproved (statement survives for generic c; proof step false pointwise)
- anchors: thqg_gravitational_complexity.tex:702-711 (proof), :843-847 (prop part iii: numerator zeros)
- evidence: step "if Sh_r != 0 then {C,Sh_r}_H != 0 so Sh_{r+1} != 0" ignores cross-term cancellation: at c = -193/45, S_5 = -174960/37249 != 0 but S_6 = 0 (numerator 45c+193 vanishes) -- the manuscript's own Proposition (iii) exhibits this. The cubic source {C,Sh_r} is only ONE summand of o^(r+1); concluding o != 0 requires the all-degree non-cancellation theorem (correctly cited as V1-thm:virasoro-rmax-infinity in thqg_holographic_reconstruction.tex:1190-1193 but omitted here).
- heal: argue at the level of rational functions in c: each S_r in Q(c) is nonzero by the cubic-source leading large-c asymptotics (S_r ~ cubic-source closed form, verified at r=6,7), hence S_r(c) != 0 off a finite set for each r; cite the non-cancellation theorem for the leading-term argument.

## S3 (sound): exact tower values S_6 = 80(45c+193)/(3c^3(5c+22)^2), S_7 = -2880(15c+61)/(7c^4(5c+22)^2) verified by direct recursion (unordered convention); potential coefficients c/4, 1/3, 5/12c(5c+22), -2/5c^2(5c+22), (45c+193)/27c^3(5c+22)^2, -4(15c+61)/49c^4(5c+22)^2 ALL = s_r/r! exactly; WIP 1/r-vs-1/r! normalisation repair is consistent at every site checked; genus-1 Hessian 120/(c^2(5c+22)) and rho = 480/(c^4(5c+22)) verified; denominator pattern and sign alternation verified r<=7.
