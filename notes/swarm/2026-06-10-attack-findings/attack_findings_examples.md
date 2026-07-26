# Attack findings: EXAMPLES LANDSCAPE AND beta_N axis
# Auditor session 2026-06-10

## F1. Denominator formula in cor:beta-N-rational-not-integer is FALSE at N=18 (and 852 of 996 N in [5,1000])
severity: major | verdict: false
anchors: chapters/theory/beta_N_closed_form_all_platonic.tex:249-251 (claim), :254-264 (proof)
evidence: Claim: den(beta_N) = lcm(1..N)/12 for every N>=5 "where this quotient is a non-unit rational".
Computed (sympy/Fraction): den(beta_18)=340340 but lcm(1..18)/12=1021020 (=3*340340). 852 failures for N in 5..1000.
Proof body premise "H_N - 1 ... has denominator equal to lcm(2,...,N) after reduction" already false at N=6: den(H_6-1)=20, lcm(2..6)=60.
Non-integrality itself SURVIVES: no integral beta_N for 5<=N<=1000.
heal: replace clause by den(beta_N)=den(H_N)/gcd(12,den(H_N)); prove non-integrality all N>=5 by Bertrand (prime p in (N/2,N], p>=5, v_p(H_N-1)=-1, p coprime to 12 => v_p(beta_N)=-1).

## F2. Bertrand strengthening ABSENT; all-N>=5 non-integrality unproved in manuscript (proof-by-table to N=12)
severity: major | verdict: missing-artifact
anchors: chapters/theory/beta_N_closed_form_all_platonic.tex:243-264, :577-578; grep "Bertrand" over chapters/ appendices/ compute/ FRONTIER.md = no hits
evidence: cor:beta-N-rational-not-integer proof says "Subsequent N values follow the same pattern; see Table" (table stops at N=12). FRONTIER.md:29 asserts "rational for N >= 5" with no proof anywhere.
heal: insert proof: for N>=5 Bertrand gives prime p, N/2<p<=N, p>=5 (direct check N=5..8; automatic N>=9 since p>N/2>=4.5); unique multiple of p in {2..N} is p itself, so v_p(H_N-1)=-1; gcd(p,12)=1 so v_p(12(H_N-1))=-1<0; hence beta_N not integer. Integral exactly at N in {2,3,4}: 6,10,13.

## F3. WIP edit broke LaTeX: \begin{itemize}[nosep] at :570 closed by \end{enumerate} at :586
severity: major | verdict: inconsistent
anchors: chapters/theory/beta_N_closed_form_all_platonic.tex:570,586; main.tex:2340 inputs the file
evidence: git diff shows WIP replaced "-\end{itemize}" with "+\end{enumerate}" in the Summary section while the opening "\begin{itemize}[nosep]" (line 570) was untouched. pdflatex error guaranteed: "\begin{itemize} ended by \end{enumerate}".
heal: restore \end{itemize} at line 586.

## F4. Three-way status inconsistency for the harmonic closed form (CONFIRMED, extends main-thread item ii)
severity: major | verdict: inconsistent
anchors: chapters/theory/beta_N_closed_form_all_platonic.tex:153-155 (Conjecture + ClaimStatusConjectured), label thm:beta-N-closed-form-proved-all-N at :154; FRONTIER.md:21+29 (listed under "Vol II-specific load-bearing closures")
evidence: chapter = Conjecture; label = "proved-all-N"; FRONTIER = closure. Also FRONTIER.md:29 states "rational for N >= 5" as if unconditional (it is conditional on the conjecture AND unproved for all N, cf. F2).
heal: rename label to conj:beta-N-harmonic-closed-form (fix all \ref), move FRONTIER entry from closures to open problems, mark "rational for N>=5" as conditional-on-conjecture (or prove via F2 heal and keep).

## F5. Proof body of thm:wn-tempered-all-N asserts a provenance that does not exist (voice-rectify inverted a disclaimer)
severity: major | verdict: overstated
anchors: chapters/theory/wn_tempered_closure_platonic.tex:150-156 (current); git dffb667~1 version lines 228-235 (original)
evidence: Original (pre-dffb667): "For N>=4, this chapter does not derive it from Fateev--Lukyanov structure constants. The MISSING INPUT is a full Miura/OPE computation... The theorem therefore takes the finite envelope as an explicit assumption."
Current: "For N>=4 the envelope COMES FROM a full Miura/OPE computation of the shadow Riccati coefficients, visible at the W_4 bridge A_5/A_4=-52/5."
No such computation exists; -52/5 is the value REQUIRED by the conjecture (same chapter :421-427 lists it as residual open; beta chapter :145-147 says A_5^{W4} "is absent"). The committed voice rectification flipped "missing input is X" into "comes from X".
heal: restore disclaimer semantics: "For N>=4 no full Miura/OPE computation of the shadow Riccati coefficients exists; the required bridge is A_5^{W4}/A_4^{W4} = -52/5; the theorem takes the finite envelope as an explicit hypothesis."

## F6. test_w4_beta_direct.py failure reconciliation: 3 stale prose pins; 1 of them was a live canary for F5
severity: minor | verdict: fragile
anchors: compute/tests/test_w4_beta_direct.py:80,87,88; chapters/theory/beta_N_closed_form_all_platonic.tex:51-56; chapters/theory/wn_tempered_closure_platonic.tex:155-156
evidence: Pin "They do not prove it." removed at 957be0d but replaced by equivalent "The all-rank proof obligation is the direct computation..." (content preserved -> test stale). Pin "finite envelope as an explicit assumption" -> now "explicit hypothesis" (synonym, test stale). Pin "does not derive it from Fateev--Lukyanov" -> sentence REPLACED by false provenance (F5): here the manuscript regressed, the test was right.
heal: after F5 heal, repin test to content-bearing strings: "proof obligation", "is absent", "explicit hypothesis", "no full Miura/OPE computation" (or similar invariant content), not exact prose sentences.

## F7. Scaling law misattributed to the closed-form conjecture (circular cross-reference)
severity: minor | verdict: inconsistent
anchors: chapters/theory/wn_tempered_closure_platonic.tex:217-218
evidence: "the kappa-ratio scaling law conjectured in Conjecture~\ref{thm:beta-N-closed-form-proved-all-N}" -- that label is the harmonic CLOSED FORM, not the scaling law; the scaling law is prop:kappa-ratio-scaling-law (beta chapter :68-85). As written the remark derives beta_N=12(H_N-1) from the conjecture asserting beta_N=12(H_N-1).
heal: cite Conjecture prop:kappa-ratio-scaling-law.

## F8. Ill-formed quantifier in Conjecture prop:kappa-ratio-scaling-law
severity: minor | verdict: fragile
anchors: chapters/theory/beta_N_closed_form_all_platonic.tex:71-72
evidence: "for every ... $r \ge 3$ with $r \ne 4$ mod parity channels vanishing" is not a well-formed mathematical condition (word salad: "r != 4 mod parity channels vanishing").
heal: state the intended scope: "for every r >= 4" (the Vol I closed form A_r(Vir)=8(-6)^{r-4}/r is established for r >= 4; odd-spin channels vanish identically), or define the parity-channel exclusion precisely.

## Also verified SOUND in these two chapters:
- prop:beta-stirling-dominance (ProvedHere): Stirling argument checked line-by-line; each factor limit correct; conclusion (beta e)/(r|c|) -> 0 valid.
- kappa-ratio arithmetic: ratio 2(H_N-1); A_4^{W3}=10/3 = (5/3)*2; A_4(Vir)=2 from 8(-6)^0/4; A_5^{W4}=-676/15=-52/5*13/3; A_{r+1}/A_r(Vir)=-6r/(r+1) from 8(-6)^{r-4}/r -- all reproduced exactly by hand.
- beta_2=6, beta_3=10, beta_4=13 integral; no integral beta_N for 5<=N<=1000 (computed).
- Table tab:beta-N-values-extended N=2..12 values verified by direct Fraction computation (pending: done below in F-table check).

## F9. rem:beta-asymptotic: beta_100 ~ 55.3 and "factor ~93" are wrong; remark contradicts its own displayed asymptotic
severity: major | verdict: false
anchors: chapters/theory/beta_N_closed_form_all_platonic.tex:342-346 (cf. :331-338 for the correct formula)
evidence: Remark derives beta_N ~ 12 log N - 5.07, then evaluates "At N = 100, beta_100 ~ 12*4.61 ~ 55.3" -- this is 12 log(100), DROPPING the -5.07 shift just derived. Exact: beta_100 = 12(H_100 - 1) = 50.2485. Triangular beta^(A)_100 = 5151; ratio 5151/50.25 = 102.5, not "~93" (93.1 = 5151/55.3 inherits the error).
heal: replace by beta_100 = 12(H_100-1) ~ 50.25 (or "12*4.19"), factor ~ 102.5 larger.

## Table check: tab:beta-N-values-extended (N=2..12, H_N, exact, decimals) verified correct by exact Fraction computation. SOUND.

## F10. (d) deletion-ledger failure reconciled: untracked test pins pre-upgrade scalar formula; WIP manuscript correctly upgraded F_g -> F_g^{sc}
severity: minor | verdict: fragile
anchors: compute/tests/test_deletion_ledger.py:2379; chapters/connections/thqg_perturbative_finiteness.tex:483-505 (WIP); theorems_C_D_native_vol2_platonic.tex:538 (delta F_2^cross(W_3)=(c+204)/(16c))
evidence: test (untracked, new) requires "F_g(\cA)=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}" inside prop:thqg-I-Fg-closed-form. WIP layer renamed to F_g^{\mathrm{sc}}(\cA)=... and added F_g = F_g^{sc} + delta F_g^{cross}. The upgrade is mathematically forced: delta F_2^cross(W_3) = (c+204)/(16c) != 0, so full F_g != kappa lambda_g^FP for multi-weight algebras. The sibling test test_wN_tensor_arakelov_weight_distribution.py:245-254 FORBIDS the old string in the native Theorem D file -- the two tests pin opposite conventions; the WIP manuscript agrees with the newer one.
HEAD remark even used BARE kappa ("F_g = \kappa \cdot \lambda_g^{FP}"), violating voice-table #6; WIP fixed it.
heal: update deletion-ledger pin to "F_g^{\mathrm{sc}}(\cA)=\kappaChHodge(\cA)\lambda_g^{\mathrm{FP}}" (manuscript side is correct).

## F11. Stale [LT] line anchors in test_wn_tempered_closure.py docstrings
severity: minor | verdict: fragile
anchors: compute/tests/test_wn_tempered_closure.py:79,90 ("tempered_stratum_characterization_platonic.tex:687 (beta_2=6)", ":698 (beta_3=10)")
evidence: lines 687/698 of the current chapter contain the wn-tempered proof/remark text, not the beta pins.
heal: re-anchor or drop line numbers from VERIFIED comments.

## Verified SOUND (tests vs manuscript):
- test_w3_pva.py: all four lambda-brackets match Zamolodchikov 1985 normalisation: {T T}=dT+2T lam+(c/12)lam^3; W primary wt 3; W_(5)W=c/3, W_(3)W=2T, W_(2)W=dT, W_(1)W=(3/10)d2T+(32/(22+5c))Lambda, W_(0)W=(1/15)d3T+(16/(22+5c))dLambda. Pole-2 coefficient 2*beta with beta=16/(22+5c) (Zamolodchikov C-structure) checked by hand; skew {W T}=2dW+3W lam derived from {T W} by hand.
- lambda_g^FP = (2^{2g-1}-1)/2^{2g-1} |B_{2g}|/(2g)! matches Faber-Pandharipande lambda_g theorem constant (int_{Mbar_{g,1}} psi^{2g-2} lambda_g).
- x/2 / sin(x/2) Taylor coefficients (1/24 at g=1, 7/5760 at g=2) verified by hand against generating-function claim in thm:thqg-I-generating-function.

## F12. FALSE IDENTITY: Lambda_0|h> = (h^2 - 3h/5)|h> ; correct value is (h^2 + h/5)|h>
severity: critical | verdict: false
anchors: chapters/examples/w-algebras-w3.tex:293-301 (eq:w3-primary-lambda-eigenvalue, WIP); chapters/connections/dnp_identification_master.tex:711-717; compute/lib/examples/w3_algebra.py:273-277; compute/tests/test_w3_pva.py:224-234
evidence: With the chapter's own definition Lambda = :TT: - (3/10) d^2 T (Notation block of the same file) and the standard mode conventions verified by the same test suite:
 (TT)_0|h> = (L_1 L_{-1} + L_0^2)|h> = (h^2 + 2h)|h>   [checked in BOTH physics ordering and Borcherds (a_{(-1)}b)_{(3)} expansion]
 (d^2 T)_0 = (0+2)(0+3) L_0 = 6 L_0
 => Lambda_0|h> = (h^2 + 2h - (3/10)*6h)|h> = (h^2 + h/5)|h>, roots {0, -1/5}.
Equivalently Zamolodchikov/Bouwknegt-Schoutens Lambda_m = :LL:_m + (x_m/5) L_m with x_0 = 1. Manuscript claims h^2 - 3h/5 with roots {0, 3/5} "independent of c"; difference = 4h/5. "Primary-line normalisation used by the DNP/Gaudin comparison" is DEFINED NOWHERE (grep over chapters/ appendices/ finds only the two assertions). The corollary "normalized Lambda-channel vanishes on primary lines with h=0 and h=3/5" is false: vanishing locus is {0,-1/5}.
test_w3_pva.py::test_primary_lambda_eigenvalue passes only because compute/lib/examples/w3_algebra.py:277 HARDCODES the same expression (circular pin, no independent route).
heal: replace eigenvalue by h^2 + h/5 = h(5h+1)/5 with roots {0, -1/5} in both chapters, engine, and test; or, if a nonstandard normalisation is truly intended, define it explicitly and derive the eigenvalue within it (currently absent).

## F13. Citation fraud: thm:w3-CYBE (ClaimStatusProvedElsewhere) licensed by garbled bibitem BouFeiMal90 pointing to an unrelated N=2 paper
severity: major | verdict: unlicensed
anchors: chapters/examples/w-algebras-w3.tex:920-930 (WIP); chapters/examples/w-algebras.tex:869,876; chapters/examples/w-algebras-conditional.tex:219,226; main.tex:3065-3066
evidence: bibitem reads "A. Bouwknegt, J. Feigin, and J. McCarthy, Resolutions and characters of irreducible representations of the N=2 superconformal algebra, Lett. Math. Phys. 23 (1991), 193-202." Problems: (1) Peter Bouwknegt's initial is P., Boris Feigin's is B. ("J. Feigin" does not exist); (2) the title is Feigin-Semikhatov-Tipunin, Nucl. Phys. B 536 (1999) 617 (N=2 superconformal), not an LMP 1991 paper by these authors; (3) the cited content (N=2 resolutions) has nothing to do with the claim it licenses: W_3 PVA Jacobi identity / Zamolodchikov associativity "in the standard finite-pole normalisation". The bibitem sits under the comment "Entries added to resolve undefined citations" (main.tex:3062). The mathematical claim itself is TRUE (Zamolodchikov 1985).
heal: cite Zamolodchikov, Teor. Mat. Fiz. 65 (1985) 347 [TMP 65 (1985) 1205]; Bouwknegt-Schoutens, Phys. Rep. 223 (1993) 183; for the PVA formulation De Sole-Kac, Jpn. J. Math. 1 (2006). Delete or repair the BouFeiMal90 entry at all four call sites.

## F14. Explicit m_3/m_4 W_3 formulas deleted (they were wrong); replacement points to a NONEXISTENT appendix section
severity: major | verdict: missing-artifact
anchors: chapters/examples/w-algebras-w3.tex:882 ("see Section~\ref{sec:w-algebras-explicit-appendix}"); WIP diff hunks deleting eq:m3-WWW explicit block and eq:m4-WWWW-explicit sigma-polynomial block
evidence: WIP replaced the (RECTIFICATION-FLAGged, mis-graded) explicit m_3(W,W,W) and m_4^tree formulas by weight-support statements (Pi_7, Pi_9). Honest retreat, BUT: label sec:w-algebras-explicit-appendix is defined nowhere in chapters/ or appendices/ (grep: zero hits for \label{sec:w-algebras-explicit-appendix}) => dangling \ref (?? in PDF). Net state: the chapter titled "Detailed Computation of m_4(W,W,W,W)" contains no explicit computation and points to a void.
heal: either supply the corrected explicit m_3 (weight-7) and m_4 (weight-9) coefficient tables in an actual appendix section carrying that label, or delete the pointer sentence; the support statements themselves verified consistent (weights 9-3-j + (3+j) = 9 etc.).

## F15. False degree-bound justification: "minimum output weight is h=1"
severity: minor | verdict: false
anchors: chapters/examples/w-algebras-w3.tex:679-680
evidence: W_3 vacuum module has no weight-1 state (generators T wt 2, W wt 3; scalar wt 0; chapter itself proves scalar output vanishes). Minimum nonscalar output weight is 2, so a+b+c <= 7, not 8.
heal: replace by "minimum output weight is h=2 (T-sector), giving a+b+c <= 7".

## F16. prop:m4-box "Proof sketch" promoted to "Proof" with no new content; proof contains a false dimension claim
severity: major | verdict: unproved
anchors: chapters/examples/w-algebras-w3.tex:792-830 (prop:m4-box + proof); WIP diff hunk "-\begin{proof}[Proof sketch] +\begin{proof}"
evidence: Proof asserts FM_4(C) has "real dimension 2x3=6" then "the integral is over a 6-4=2-dimensional family" (a 4-form is not integrated over a 6-fold by passing to a 2-dim family without naming the fibration), and "absence of codimension-0 boundary contributions (all boundary strata have dimension < 2)" -- boundary strata of a 6-dim (or 5-dim) FM compactification have dimension 5 (codim 1), not < 2. The four claimed properties (polynomiality, weight 9, (c/360)^2 leading, c=0 behaviour) are plausible but NOT proved by this text.
heal: restore "Proof sketch" or supply the real argument: fibre-integrate the 4 dlog forms along the FM_4 fibres of the forgetful map, show boundary-strata contributions vanish by the Theta-ordering support of K(z,t)=Theta(t)/(2 pi z), and exhibit the wheel-integral evaluation; name the 2-dim family explicitly.

## F17. thm:central-charge-shift: untagged theorem, physics-heuristic proof, sketch->Proof promotion
severity: major | verdict: unlicensed
anchors: chapters/examples/w-algebras-w3.tex:1038-1054 (statement+proof); WIP hunk "-\begin{proof}[Proof sketch] +\begin{proof}"
evidence: Theorem asserts the algebraic duality A_c^! = A_{alpha-c}, alpha = sum 2(6 s_i^2 - 6 s_i + 1), self-dual c = alpha/2. Proof body is a 1-loop ghost vacuum-bubble heuristic with a factor-1/2 "chiral half of the ghost determinant" argument -- a physics metaphor, not a proof of the !-duality statement. No ClaimStatus tag on the theorem (only :907 in the file has one). NOTE verified: alpha-arithmetic is sound (26/100/246 = 2(N-1)+4N(N^2-1), cross-checked vs kappa+kappa^!=250/3 at W_3, c+c'=100); and this complementarity is NOT Feigin-Frenkel duality (FF-dual levels (k+3)(k'+3)=1 give c(k)=c(k'), checked symbolically) -- the corpus-internal Theorem C is the only available license.
heal: tag \ClaimStatusProvedElsewhere{Vol I Theorem C} (or the in-volume complementarity theorem), demote the ghost argument to a labelled heuristic remark, restore "Proof sketch" if kept.
