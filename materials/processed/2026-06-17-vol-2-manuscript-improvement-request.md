# Manuscript Improvement Request: Processed Triage

Source: `/Users/raeez/Desktop/vol 2Manuscript Improvement Request.pdf`.

Raw copy: `materials/raw/2026-06-17-vol-2-manuscript-improvement-request.pdf`.

Extraction date: 2026-06-17.

## Adopted

A432 adopts the PDF's rigorous first-page typing corrections.  The
abstract now defines the chiral derived centre as
\[
\cZ^{\mathrm{der}}_{\mathrm{ch}}(\cA)
:=\mathrm{REnd}^{\mathrm{ch}}_{\cA-\cA}(\cA)
\simeq C^\bullet_{\mathrm{ch}}(\cA,\cA)
\]
in the chiral \(\cA\)--\(\cA\)-bimodule category, with the Hochschild
presentation licensed by the chiral bar resolution on the Koszul locus.
The physical functor now has source
\(\mathrm{HTReal}^{\mathrm{Kosz},\omega,\mathrm{BL},\Xi}_X\), so the
BV/HT realisation datum is an input, not a consequence of a bare chiral
algebra.

A432 also adopts the PDF's \(\SCchtop\) typing correction on active
first-page and theorem-near surfaces: absent an explicitly constructed
operadic envelope, \(\SCchtop\) is written as the directed
two-coloured dioperad.  The \(E_3^{\mathrm{top}}\) statement is
post-topologisation, after \(Q\)-exact holomorphic translations.

A432 adopts the PDF's \(C_2\)-cofiniteness correction on the repaired
surfaces.  The standard landscape now separates rational/lisse
\(C_2\)-cofinite loci from generic non-critical tempered Virasoro and
\(\mathcal W\)-loci.

A432 adopts the verified Virasoro TQFT bibliography correction.  arXiv
metadata verifies Collier--Eberhardt--Zhang,
\emph{Solving 3d gravity with Virasoro TQFT}, arXiv:2304.13650, and
\emph{3d gravity from Virasoro TQFT: Holography, wormholes and knots},
arXiv:2401.13900.  The stale "Virasoro TQFT and Teichmuller theory"
entry at arXiv:2311.15063 was removed from the live bibliography.

A433 adopts the PDF's Fulton--MacPherson correction.  The FM calculus
now separates unreduced \(\FM_X(k)\), translation-reduced
\(\FM_k^{\mathrm{tr}}(\C)\), scale-reduced
\(\FM_k^{\mathrm{red}}(\C)\), and full affine-reduced
\(\FM_k^{\mathrm{aff}}(\C)\) spaces.  Stokes calculations use
translation-reduced coordinates \(w_i=z_i-z_k\); collision divisors
have scale-reduced inner factor
\(\FM_{|S|}^{\mathrm{red}}(\C)\).  The full affine quotient is not used
in the residue calculation because it erases the angular residue.

A433 also adopts the PDF's closed/open collision-stratum correction.
Non-consecutive divisors such as \(D_{\{1,3\}}\) are genuine closed
chiral FM divisors.  They vanish only after the open planar
\(E_1\)-ordered projection for tree-level planar weight forms.  The
dressed PVA descent and the closed chiral collision calculus retain
the non-consecutive divisors.

A434 harvests the remaining rigorous items as scoped theorem data.
The construction chapter now contains
`prop:standard-family-bv-realisation-data`, which says that a
standard-family HT realisation starts from a field complex, BRST
operator, interaction, BV pairing, gauge fixing, boundary Lagrangian,
scale-\(L\) QME solution, and bulk/boundary comparison maps.  The rows
for Heisenberg/lattice, affine \(V_k(\fg)\), beta-gamma/symplectic
boson, Virasoro, and \(\mathcal W_N\) are conditional where the BV
construction is transported from DS/KZ or harmonic-decoupling data.
The algebra \(A\) alone is explicitly not enough to construct its BV
theory.

A434 also harvests the quantum-group lane.  The spectral-braiding core
now states the KZ--Drinfeld--Kohno--Kazhdan--Lusztig proof path with
\(\hbar_{\mathrm{KZ}}=(k+h^\vee)^{-1}\), the KZ connection, the
trace-form/KZ-kernel separation, \(q_{\mathrm{KL}}=\exp(\pi
i/(k+h^\vee))\), and \(q_{\mathrm{DK}}=q_{\mathrm{KL}}^2\).  The
same pass propagated the typed FM reductions to the active spectral
surfaces: Stokes runs on \(\FM^{\mathrm{tr}}\), collision fibres use
\(\FM^{\mathrm{red}}\), and total collision records the associator
rather than disappearing.

A434 formulates the boundary/defect correction as
`cor:boundary-defect-slab-package`.  The slab package is the four-term
datum
\[
(B_\partial,\ A_{\mathrm{bulk}}\simeq \cZ^{\mathrm{der}}_{\mathrm{ch}}
(B_\partial),\ \widehat{\cC}_{\mathrm{line}},
(\beta_z^{\mathrm{hol}},\beta^{\mathrm{cat}})).
\]
The meromorphic spectral exchange, categorical braiding, line-module
comparison, and bulk centre are separate pieces of data; no one object
determines the other three by itself.

A434 records the CFG separation as a regression guard on
`appendices/cfg_side_by_side.tex`.  CFG supplies the ordinary
topological Chern--Simons factorisation lane.  Vol II hCS and
chiral-avatar statements use Dolbeault BV fields, holomorphic
topological reduction, and chiralisation data.  Equality between these
lanes requires a named comparison functor; it is not inferred from the
formal similarity of the BV equations.

A434 tightens the gravity scope.  The Universal Holography functor
surface retains the chiral Brown--Henneaux exact sector for
\(\mathrm{Vir}_c\).  The full non-chiral physical gravity partition
function now explicitly requires an anti-holomorphic Virasoro copy,
a real form relating \((c,\bar c)\) to the two Chern--Simons levels,
a modular invariant pairing of chiral blocks, an integration cycle in
real metrics or flat connections, and a saddle prescription.

A439 propagates the same discipline to the half-space reflected-weight
theorem.  Reflected Stokes weights verify already supplied
renormalized BV data; they do not construct the Costello effective
action, counterterms, or all-loop QME.

A498 closes the adjacent BV/bar comparison layer left implicit in the
same datum package.  The graph-residue boundary operations
\(m_k^{\mathrm{BV}}\) are now compared to the chiral bar coderivation
only after a boundary SDR \(i_\partial,p_\partial,h_\partial\) transfers
the BV operations to \(T^c(s^{-1}\bar A_\partial)\).  The live
construction prints \(D_A^{\mathrm{BV}}\) and the equality
\[
\pi_1D_A^{\mathrm{BV}}\big|_{(s^{-1}A_\partial)^{\otimes k}}
=
\pi_1D_A\big|_{(s^{-1}A_\partial)^{\otimes k}},
\]
and explicitly rejects bare pointwise equality between raw graph
residues and the bar coderivation.

A440 propagates the CFG/hCS lane separation to the live
spectral-braiding comparison remark.  CFG and the Vol II annular bar
tower are compared through a supplied trace-comparison datum and a
common perturbative trace shadow on standard lanes, not through a
global comparison functor between ordinary Chern--Simons, hCS, and the
chiral avatar.

A466 closes the final PDF \(W_3\) quartic odd-sector item.  The live
\(m_4(W,W,W,W)\) theorem now separates the weight-nine support
envelope from the actual odd-sector value: for a \(\sigma\)-equivariant
weight-preserving contraction with \(\sigma(T)=T\) and
\(\sigma(W)=-W\), the odd component vanishes.  The old proof-sketch
wording and the suggestion that the weight projection fixes a non-zero
\(W\)-valued coefficient are absent from the live \(W_3\) source.

## Remaining Mathematical Work

The source PDFs do not supply proofs of the full non-chiral gravity
path integral, an all-loop BV construction for every class-\(\mathsf C\)
or class-\(\mathsf M\) family, or the existence of all components of a
global CFG--hCS comparison datum.  After A499, the class-\(\mathsf C\)
all-loop BV/bar gap is no
longer a vague harmonic-decoupling slogan: it is the vanishing and
compatible lifting of the obstruction classes
\(\mathfrak h_{g,r}^{\mathsf C}\) for all \(g\ge1,r\ge4\).  The
manuscript now treats these as explicit hypotheses, conditional
constructions, or open construction problems rather than as silently
proved consequences.  After A500, the CFG/hCS comparison gap is no
longer an unnamed missing functor: it is the comparison datum
\[
\mathfrak C_{\mathrm{CFG}\to\hCS}
=
(\mathcal R_{\mathrm{ell}},\rho_{\mathrm{BV}},\rho_{\mathrm{fac}},
\tau_{\mathrm{tr}},\{h_n,c_n\}_{n\ge1},
\nu_{\mathrm{Mell}},\eta_{\mathrm{K3}})
\]
together with the obstruction profile
\(\mathfrak o_{\mathrm{CFG}\to\hCS}\).  A global comparison theorem is
the existence and compatibility of the required components, not a
consequence of the row-wise CFG comparison table.

After A501, the six-slot fingerprint theorem is likewise scoped by its
missing reconstruction input: equality of the scalar slots is
reconstructive only after the classwise datum \(\mathfrak R_c\) is
installed, and the \(M_{\mathrm{Kosz}}\) and summary layers no longer
state unconditional fingerprint completeness.

After A503, the earlier physical-origin boundary theorem also carries
the same BV/factorization discipline as the later BV statement.  The
manuscript no longer says that an HT boundary condition alone supports
a chiral algebra; the theorem now requires a chosen boundary chart,
renormalized boundary factorization algebra, bulk-to-boundary OPE
comparison, level normalization, and the
Costello--Dimofte--Gaiotto factorization-to-chiral comparison.

## Verification

After A432, the focused Hochschild bundle passed with `58 passed`.
After A433, the focused FM/Arnold bundle passed with `7 passed`.
After A434, the focused harvest bundle passed with `44 passed` across
the BV-datum, spectral-braiding, slab-package, CFG-separation,
universal-holography, and FM guards.  A439 passed `44` affine
half-space BV guards.  A440 passed `2` CFG/hCS lane-separation guards.
A498 passed the focused BV bundle with `11 passed` and the broader
harvest-adjacent cluster with `221 passed`.
A499 passed the focused class-C/BV bundle with `15 passed` and the
broader harvest-adjacent cluster with `225 passed`.
A500 passed the focused CFG/hCS guard with `3 passed`, the adjacent
harvest bundle with `82 passed`, and the half-space/class-C/BV/CFG
bundle with `60 passed`; `make verify-licensing` again reported zero
blockers and zero warnings.
A501 passed the focused fingerprint guard with `4 passed` and the
adjacent fingerprint/input-harvest bundle with `47 passed`; `make
verify-licensing` again reported zero blockers and zero warnings.
A502 passed the combined bar-cobar/moonshine/fingerprint/input bundle
with `60 passed`; the final `make fast` converged after two passes with
zero undefined citations, zero undefined references, and zero rerun
requests, and the final log scan found no fatal errors, undefined
controls, unresolved citation/reference warnings, or rerun requests.
A503 passed the focused BV/HT guard with `11 passed`; the consolidated
deletion/input-harvest/BV bundle passed with `109 passed`.  `make
verify-licensing` reported zero blockers and zero warnings, and the
final `make fast` converged after two passes with zero undefined
citations, zero undefined references, and zero rerun requests; direct
log scan found no fatal errors, undefined controls, unresolved
citation/reference warnings, or rerun requests.
Current rerun: the June 17 Manuscript Improvement guard bundle passes
`96` checks across the FM, BV-datum, spectral-braiding, slab-package,
CFG/hCS, half-space BV, and Universal Holography surfaces.
After A457, the active theorem/status surface no longer contains the
undefined fused `\ClaimStatusProvedHereConditional` token except in
negative source-guard assertions.  After A458--A460, the active
post-harvest surface also separates Schellekens Type C lattice
compatibility from BV/DW anomaly construction, removes audit prose from
the Three-Hochschild theorem, and titles the DS--Hochschild dichotomy
as a proved class-stratified criterion rather than a conditional-status
fusion.  A461--A462 then remove the finite-Zhu/Massey overclaim and
the remaining theorem-facing audit markers, replacing them by scoped
finite-Zhu amplitude, GCFT/Pic, oper, and Weyl/Gel'fand--Fuchs
statements.  A463 keeps logarithmic \(\mathcal W(p)\) out of the
inherited \(E_3\)-topological finite-orbifold/coset lane until its
logarithmic HT datum, antighost/topological lift, descent data, and
regular/logarithmic amplitude bounds are supplied.  A464 keeps the
high-genus gravity partition-function remark at scalar-shadow
weight-completed scope and leaves the original-complex upgrade at
Frontier F5.  A466 closes the final PDF \(W_3\) odd-sector
parity-vanishing item.  Current reruns: the consolidated input-harvest guard
cluster passes with `490 passed`, and the adjacent post-harvest
status/scope cluster passes with `392 passed, 1 skipped`; the focused
A462 suite passes with `28 passed`, and the focused A463 gravity guard
passes with `211 passed, 1 skipped`; the focused A464
programme/gravity/curved-Dunn bundle passes with `228 passed, 1
skipped`.
`make verify-licensing` reported zero blockers and zero warnings.
Scoped `git diff --check` passed.  After A434, `make fast` converged
after two passes at 2519 pages with zero undefined citations, zero
undefined references, and zero rerun requests; direct log scan found no
fatal TeX errors, undefined controls, unresolved citation/reference
warnings, rerun requests, or label-change warnings.
