# Research Paper Refinement Vol II: Extracted Proof Obligations

Source: `materials/raw/2026-06-08-research-paper-refinement-vol2.pdf`.
The PDF is a 124-page browser export created on 2026-06-08 at 05:31
Asia/Ho_Chi_Minh. Its useful mathematical content begins after the
initial blank browser pages.

## Primary Obligations

The document's strongest criticisms are formula-level obligations.
They are not requests for style-only rewriting.

1. The two-coloured chain object must be written as operations with the
   no-open-to-closed component equal to the empty component, not as
   prose.
2. The \(A_\infty\)-chiral algebra and the chiral Hochschild complex
   must be written as coderivations with fixed suspension signs.
3. Sesquilinearity must be derived from translation equivariance and
   then used uniformly.
4. Fulton--MacPherson residue signs must control the operadic boundary
   calculus, including the edge-ordering orientation convention.
5. The Arnold relation must be the printed Jacobi proof at arity four,
   with the arity-five obstruction computed separately.
6. The missing higher homotopy coherence must be a tree formula for the
   contracting homotopy, with internal-edge orderings and signs.
7. The \(E_3\)-PBW route must be a theorem with an associated graded
   identification, not a hope.
8. The BV construction must specify fields, kernels, and the QME.
9. The Heisenberg, affine, Virasoro, and \(W_3\) tests must be fully
   normalized before they serve as benchmarks.
10. Virasoro must be treated as an infinite HPL tower in a
    completed/pro ambient.
11. DS--Hochschild must be an HPL transfer theorem of brace
    \(E_2^{\mathrm{ch}}\)-algebras, and an SC pair only after the
    mixed estimates and higher homotopies are supplied.
12. The modular graph complex, genus-one Arakelov kernel, and
    tensor-Arakelov curvature must replace scalar shortcuts at higher
    genus.
13. Anomaly completion must state the curved Maurer--Cartan convention
    \(d\eta+\eta^2+\Theta=0\) and distinguish it from the strict
    central transgression convention \(d\eta=\Theta\), \(u=\eta^2\).
14. The K3/Borcherds sector needs an operator-level Pfaffian/Hall
    square, not only a scalar modular identity.
15. The celestial bridge must be a residue theorem, not an analogy with
    scattering.
16. Physical gravity claims must be separated into exact algebraic
    statements, conditional analytic statements, and heuristic physics.

## Harvest Status Through A470

The primary-obligation list is historical.  It remains useful as an
index of pressure points, but a correction is counted as harvested only
when the live manuscript or compute layer has either proved it,
weakened it to the correct conditional statement, or installed a guard
against the false stronger form.

1. Two-coloured operations and the no-open-to-closed component have
   been harvested on the active \(\SCchtop\) surfaces and guarded by
   the chain-dioperad tests.
2. \(A_\infty\)-chiral coderivation and chiral Hochschild sign models
   have been harvested into the coderivation/Hochschild guards.
3. Sesquilinearity is treated as translation-equivariance data on the
   repaired chiral algebra surfaces; surviving uses are scoped by the
   local tests rather than by prose analogy.
4. FM residue signs and edge-ordering orientation were harvested in
   A433, with the translation-reduced, scale-reduced, and
   full-affine-reduced spaces separated.
5. The Arnold arity-four Jacobi lane is guarded together with the FM
   residue-sign tests; arity-five remains a separate obstruction
   computation where it is not locally proved.
6. Higher homotopy coherence has been harvested on the main live
   surfaces that used it.  A318 prints the signed two-coloured
   Boardman--Vogt cobar homotopy and its obstruction class for Chiral
   Higher Deligne.  A435 propagates a suspended Kadeishvili transfer
   formula to the celestial boundary-transfer theorem: internal edges
   are ordered, the orientation line is fixed, SDR side conditions are
   stated, and the old free \(\pm\)-sign placeholder is removed.  A437
   prints the corresponding suspended HPL transfer formula on the
   DS--Virasoro gravity surface, with source/target suspensions,
   left-to-right internal-edge order, orientation line, and no free
   sign choice.  A443 removes a remaining low-degree ordered-bar sign
   placeholder in the \(\widehat{\mathfrak{sl}}_2\) computation: the
   differential now displays the \((-1)^{i-1}\) ordered-bar sign, the
   \(\mathfrak{sl}_2\) normalisation, and the desuspended
   \(\Sigma_2\)-coinvariant sign explicitly.  Remaining claims that
   need a different tree homotopy must either display the same level of
   data or stay conditional.  A444 extends the fixed-sign convention
   back to the general low-degree ordered-bar construction and the
   Heisenberg ordered-bar computation, and A445 separates the
   associative open product from the affine \(0\)-mode Lie bracket:
   the manuscript now prints the face sign
   \(d_i=(-1)^{i-1}\partial_i\), the degree-\(2\) and degree-\(3\)
   \(\mu_2\)-formulas without free signs, the associativity
   cancellation with computed coefficients, and
   \(d(s^{-1}\alpha\otimes s^{-1}\alpha)=k\,s^{-1}(1)\).  The
   constants \(c^K_{IJ;n}\) now belong to the vertex-shadow
   symmetric/Arnold calculation, not to the open \(E_1\) proof.  A446
   applies the same separation to the Yang--Baxter-from-\(d^2=0\)
   boundary-face proof, so its two algebraic contraction paths are
   \(\mu_2\)-associativity paths and its \(R\)-matrix part is the
   monodromy/Yang--Baxter datum.  A447 applies the separation to the
   bar-cohomology concentration theorem and \(R\)-matrix leading term:
   \(d_\mu\) is the open-colour first differential,
   \(Q_{\mu_2}(\bar\cA)\) is the degree-one quotient, and the
   Lie-homology pole-depth spectral sequence is only the
   vertex-shadow symmetric chiral/Francis--Gaitsgory sequence.  A448
   applies the same correction to the Yangian double-bar claim:
   \(U(\fg[t])\) is the PBW associated graded/current shadow of
   \(Y_\hbar(\fg)\), not the full ordered bar object
   \(\Barchord(Y_\hbar(\fg))\).  A449 propagates the same correction
   into the executable sl2/sl3 oracles: the Yangian \(E_1\) product is
   the RTT associative multiplication \(\mu_2\), and the current
   Lie bracket is only the primitive collision-residue shadow.  A450
   applies the same scope to the quantum-lattice comparison API: the
   full Yangian is \(\hbar\)-visible, and only the PBW primitive
   collision-residue shadow loses \(\hbar\).  A451 scopes the
   q-stability oracle as an explicit \(A_1\) finite-window benchmark,
   not a firstness theorem or a complete all-higher-operations theorem.
   A452 separates the bare \(q\to1\) affine/loop limit from the
   additive rational scaling degeneration that produces the Yangian.
7. The \(E_3\)-PBW route is not treated as a hope where it is used.
   A319/A426 force conjectural status, Rees-flatness, associated
   graded identification, and spectral-sequence support.  A436 sharpens
   the certificate further: the object \(W_x\) must be the shifted
   \(E_3\)-tangent datum
   \(\mathbb T^{E_3}_{R_x,x}[-2]\), and
   \(H^\bullet(W_x[-2])\) must identify with the
   \(E_3\)-indecomposables of \(\operatorname{gr}_F R_x\).  Thus the
   PBW package cannot be satisfied by a tautological choice
   \(W_x=R_x\).  General extensions remain conditional.
8. BV fields, kernels, QME, and boundary data were harvested first in
   A181/A320 and then sharpened in A434 by the standard-family BV
   realisation datum table and comparison maps.
9. Heisenberg, affine, Virasoro, and \(W_3\) normalization gates were
   harvested across the A182--A188 sequence and remain guarded by the
   family-specific compute tests.  A454 adds a current input-harvest
   guard for the Desktop-PDF \(W_3\) \(\Lambda\)-channel obligation:
   the live W3 surface must retain the literal \(P_\Lambda\)
   projection, coefficient-one \(WW\Lambda\) channel, FM pair-corner
   \(m_{4,\Lambda}\) kernel, primary-line \(\Lambda_0\) eigenvalue,
   matrix-valued \(K^{W_3}\), and the statement that Gaudin
   commutativity requires a supplied \(\mathcal W_3\) Gaudin/RLL
   algebra rather than PVA Jacobi alone.
   A470 removes the remaining Heisenberg scalar/full overclaim from
   Rosetta, modular bootstrap, and the preface: the rigorous
   Faber--Pandharipande statement is
   \(F_g^{\mathrm{sc}}(\cH_k)=k\lambda_g^{\mathrm{FP}}\) on the scalar
   Hodge lane, not a bare full \(F_g(\cH_k)\) coefficient.  The same
   preface rung now separates the algebraic derived-centre bulk
   candidate from the physical HT bulk reading.
10. Virasoro is treated on completed/pro HPL surfaces where the tower
    is infinite; raw finite-type claims are not accepted without the
    completed ambient.
11. DS--Hochschild is now an HPL/coderivation-transfer theorem on the
    completed chiral bar surface; A189 and A430 remove the raw
    functoriality shortcut.  A455 pins the BRST complex, Hochschild
    SDR, conditional status, bounded-shift/weightwise convergence,
    and planar-rooted-tree brace formula in the input-harvest guard.
12. Higher-genus modular graph, genus-one Arakelov, and
    tensor-Arakelov curvature shortcuts have been converted into
    scoped scalar-shadow or conditional comparison statements where
    the local proof surface does not supply a theorem.  A455 guards
    the explicit relative modular graph theorem, elliptic Arakelov
    Green kernel and scalar tower, and tensor-Arakelov stable-graph
    component formula with its \(W_3\) channel matrix.
13. The anomaly-completion convention was repaired on 2026-06-08: the
    curved Maurer--Cartan equation \(d\eta+\eta^2+\Theta=0\) is
    separated from the strict central transgression convention.  A455
    guards the curved-MC anomaly-completion definition and Virasoro
    class-\(\mathbf M\) example.
14. The K3/Borcherds sector is not treated as a scalar modular
    identity.  The operator Pfaffian/Hall square is an explicit
    construction problem unless the relevant recognition hypotheses
    are present.  A438 makes the conditional square type-correct: the
    Hall source is the K3 fibre, but the chiral target is the
    \(K3\times E\) algebra \(A_X\).  The theorem now inserts
    \((-) \boxtimes \mathcal O_E\), \(\PhiFA_3\), and
    \(\SpCh_{E,C}\) before \(\Zderch(A_X)\).  A442 propagates the
    same typed route to the main K3\(\times E\) bridge summary, so the
    high-visibility operator paragraph no longer bypasses the
    total-space lift.  A455 guards the K3/Borcherds theorem as
    conditional on the P1 datum, finite Hall gates,
    \(\hypAmbientWtCpl\), and
    \(\effPfaffOrient+\effPBWnoExtra\), not as a scalar-character
    theorem.
15. Celestial statements are required to be residue theorems or
    labelled heuristic/conditional bridges; analogy-only statements
    are not counted as proofs.  A441 repairs the four-point MHV
    surface by replacing the unqualified bar-integral reproduction
    claim with a conditional cyclic FM boundary-residue identity:
    colour-order projection, spinor-helicity residue characters,
    negative-helicity pairing, and orientation are now explicit
    comparison data.
16. Gravity claims have been repeatedly separated into exact chiral
    algebraic statements, conditional analytic/non-chiral path-integral
    statements, and heuristic physics; A432 and A434 add the latest
    Brown--Henneaux and full non-chiral data gates.

The remaining unsolved work is therefore research-level construction:
global all-loop BV realisations beyond the stated standard-family
input, full non-chiral 3d gravity path-integral data, complete
higher-homotopy tree formulas where not printed, and operator-level
K3/Borcherds and celestial construction problems.  These are not
unharvested PDF corrections; they are now explicit mathematical
obligations.  A499 sharpens the class-\(\mathsf C\) all-loop BV/bar
piece of this list: the remaining condition is the simultaneous
vanishing and compatible lifting of the obstruction classes
\(\mathfrak h_{g,r}^{\mathsf C}\), not an undefined
"harmonic decoupling" slogan.  A439 further scopes the half-space reflected-weight
theorem as a conditional reduction: reflected Stokes weights verify
already constructed renormalized BV data, but do not by themselves
construct the Costello effective action, counterterms, or all-loop
QME.  A440 further scopes the CFG comparison surface: the spectral
remark now compares CFG and the Vol~II annular bar tower only through
a supplied trace-comparison datum and a common perturbative trace
shadow on standard lanes, not through a global comparison functor
between ordinary Chern--Simons, hCS, and the chiral avatar.  A500
replaces the remaining loose "global comparison functor" obligation by
the typed datum
\(\mathfrak C_{\mathrm{CFG}\to\hCS}\) and obstruction profile
\(\mathfrak o_{\mathrm{CFG}\to\hCS}\): the missing theorem is the
existence and compatibility of the elliptic replacement, BV
comparison, factorisation-algebra comparison, trace component
\(\tau_{\mathrm{tr}}\), regulator homotopies/counterterms,
Mellin-residue normalisation, and optional K3/Hall--BKM recognition
component.  A441
further scopes the celestial amplitude bridge: the four-point
Parke--Taylor formula is harvested as a conditional residue theorem,
while the Lorentzian scattering map, LSZ map, Mellin/celestial
transform, and general-\(n\) BCFW recursion remain explicit input or
construction data.  A442 further propagates the K3/Borcherds
operator-square typing to `main.tex`: the Hall source remains the K3
fibre, but the chiral branch passes through
\((-) \boxtimes\mathcal O_E\), \(D^b\mathrm{Coh}(K3\times E)\),
\(\PhiFA_3\), \(\SpCh_{E,C}\), and \(\Zderch(A_X)\).  A444--A452
further propagate the fixed ordered-bar sign convention from the
\(\widehat{\mathfrak{sl}}_2\) worked block back to the general
low-degree construction, the Yang--Baxter boundary-face proof, and the
bar-concentration depth spectral sequence; they also retarget the
Yangian double-bar equality to its PBW/current associated-graded
shadow, separating the open associative product and full
\(\hbar\)-dependent Yangian object from affine current Lie-bracket
shadows, and they keep the quantum-lattice q-stability evidence at its
proved \(A_1\) finite-window scope with the correct
quantum-affine/Yangian degeneration convention.

A501 applies the same discipline to the six-slot fingerprint theorem:
the rigorous statement is conditional branchwise reconstruction via
\(\mathfrak R_c\), not an unconditional assertion that six scalar slots
alone classify standard-landscape chiral algebras up to isomorphism.
The \(M_{\mathrm{Kosz}}\) and summary layers now carry the same
branchwise reconstruction datum.

A502 records the build-surface cleanup forced by the A501 verification
run: the active bar-cobar projection table now uses \(\Barch_X\) rather
than the undefined `\Bbarch_X`, and the umbral moonshine citation
`DuncanGriffinOno15` has a matching bibliography alias.

A504 closes a remaining scalar/full Hodge-lane residue in the \(M2\)
zeta benchmark.  The proposition now writes the genus-one scalar as
\[
F_1^\zeta=\kappa_{M2}^{\zeta}/24,
\]
and keeps the actual DDCA trace
\(\kappa_{M2}=\kappaChHodge(A_{M2,\infty})\) separated from
\(\kappa_{M2}^{\zeta}\) by the zeta-comparison conjecture.  The bare
form \(F_1=\kappa/24\) is now a rejected stale form.

A455 extends the input-harvest source guard to the mid-PDF
A188--A194 block; the focused guard passes `21` checks and the paired
ownership bundle passes `267` checks.

## Desktop 05:57 Extended Harvest Layer

The 05:57 Desktop export contains an additional local-theorem layer
after the primary sixteen obligations.  The rigorous harvestable part
is now indexed as A453 and guarded by
`compute/tests/test_vol2_input_harvest_guards.py`.

- Spectral discriminants and Casimir recurrence are treated as cyclic
  transport modules and reduced \(H^{1,\alg}\) data in
  `chapters/connections/casimir_divisor_core_transport.tex`.
- Relative Hochschild duality, instanton completion, and screening
  estimates are harvested through the \v Cech, Novikov, Hodge, and gap
  packages in the Yang--Mills higher-body and instanton-screening
  chapters.
- Celestial transfer is guarded as an obstruction tower and a
  conditional residue identity, not as a dictionary with Lorentzian
  scattering.
- The chiral CE, 6d hCS, Kontsevich-integral, non-principal
  \(W\)-algebra, annular supersymmetric-index, chiral Springer, and
  Yang--Mills central-formality fragments are each tied to active
  theorem/definition labels in the live `main.tex` input graph.
- The final false-pattern layer is harvested as scope discipline:
  product-formal \(\SCchtop\) recognition, completed class-\(\mathbf M\)
  ambient, operator-valued Virasoro \(R\)-matrix corrections,
  conditional Page/Stokes data, pro-window \(W_\infty\), conjectural
  harmonic \(\beta_N\), arity-four \(m_4\) scope, bulk--boundary cone
  reconstruction, one-wheel genus-one modular class, and bar-valued
  logarithmic superconnection are guarded against their stronger false
  forms.
- A466 closes the remaining \(W_3\) quartic odd-sector correction from
  the final PDF harvest pass: the \(W\)-odd part of
  \(m_4(W,W,W,W)\) vanishes for a \(\sigma\)-equivariant contraction,
  while the old weight-nine display is only a support envelope.

## Current Verification

After A457, the active theorem/status surface no longer contains the
undefined fused `\ClaimStatusProvedHereConditional` token except in
negative source-guard assertions.  A458--A460 then close three adjacent
scope defects on the harvested surface: the Coxeter--Todd Type C row
does not compute the BV/DW anomaly class from fixed-lattice data alone;
the Three-Hochschild theorem states the positive comparison rather than
reader-facing audit prose; and the DS--Hochschild class dichotomy is a
proved class-stratified criterion, with class-\(\mathsf M\) chain-level
promotion still tied to the bounded-shift HPL convergence obligation.
Current reruns: the consolidated input-harvest guard cluster passes
with `490 passed`, and the adjacent post-harvest status/scope cluster
passes with `392 passed, 1 skipped`.  A461 then closes the
finite-Zhu/Massey scope defect: finite-dimensional Zhu algebra is not
used as a proof of bounded Massey growth; generic Virasoro and generic
affine rows are separated from \(C_2\)-cofinite finite-Zhu rows; and
the logarithmic \(W(p)\) frontier is a regular/logarithmic amplitude
problem.  A462 closes the adjacent source-truth cleanup by replacing
residual protocol/all-caps audit markers with scoped GCFT/Pic, oper,
Weyl/Gel'fand--Fuchs, finite-envelope \(W(p)\), and finite-Zhu/BGG
provenance statements.  The focused finite-Zhu/logarithmic/programme
bundle passes with `41 passed`; the focused A462 suite passes with `28
passed`.  A463 then removes the remaining gravity-scope overclaim that
logarithmic \(\mathcal W(p)\) is inherited \(E_3\)-topological through
a finite-orbifold/coset lane; the focused A463 gravity guard passes
with `211 passed, 1 skipped`.  A464 then restricts the high-genus
gravity trace to the scalar-shadow comparison in the
weight-completed ambient and leaves the original-complex upgrade at
Frontier F5; the focused programme/gravity/curved-Dunn bundle passes
with `228 passed, 1 skipped`.  A466 then closes the final PDF
\(W_3\) odd-sector parity item, with its focused W3/harvest bundle
passing `27 passed`.  A470 closes the remaining Heisenberg scalar/full
export drift, with its focused scalar-linearity suite passing
`10 passed`.  A504 closes the M2 zeta bare-\(\kappa\) residue: the
focused BV/HT, M2, and scalar-lane guards pass with `13 passed`, and
the consolidated deletion/input-harvest/BV bundle passes with `109
passed`.  `make verify-licensing` reports zero blockers and zero
warnings, and the final `make fast` converges after two passes with
zero undefined citations, zero undefined references, and zero rerun
requests; direct log scan is clean for fatal errors, undefined
controls, unresolved citation/reference warnings, and rerun requests.

## Progress Recorded On 2026-06-08

Obligation 13 was repaired locally:

- `chapters/connections/anomaly_completed_core.tex` now contains the
  curved-MC identity
  \(d_\eta^2=[\Theta+d\eta+\eta^2,-]\) and states the flatness equation
  \(d\eta+\eta^2+\Theta=0\).
- `chapters/connections/anomaly_completed_topological_holography.tex`
  now states that the Virasoro bridge uses the strict central-shadow
  Ore convention, while the curved class-\(\mathsf M\) bar complex uses
  the curved-MC equation.
- `compute/lib/anomaly_completed_engine.py` now treats \(B_\Theta\) as
  infinite-rank over \(B\) unless an explicit \(\eta\)-power cutoff is
  supplied.
- `compute/tests/test_anomaly_completed_engine.py` verifies the
  infinite-rank default, finite \(\eta\)-truncations, and the
  curved-MC convention.

Verification:

- `compute/.venv/bin/python -m pytest compute/tests/test_anomaly_completed_engine.py -q -ra`
  passed: 34 tests.
- `make fast` converged after two passes with zero undefined citations
  and zero undefined references on pass 2.
