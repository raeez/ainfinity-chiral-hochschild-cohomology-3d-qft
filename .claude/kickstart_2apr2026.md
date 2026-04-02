# Kickstart Prompt: 2 April 2026 Session

Read this BEFORE touching any file. It captures the full state of the 2 April 2026 session.

---

## A. SESSION SUMMARY

The 2 April session was a 12-hour deep dive with ~43 agent worktrees running adversarial audits. Three major threads:

1. **The Polyakov programme revisited** -- 23 adversarial agents attacked the entire bar-cobar programme from first principles. Results are documented in working_notes.tex sections 15447-15828 (session notes) and 15830-16257 (second-wave findings, Wick rotation architecture, chromaticity). Key outcomes: arithmetic corrections to Brown-Henneaux specialisation, BP complementarity fixed (c_sum = 2 not 76), genus-g curvature-braiding uniformity proved, modular covariance from MC proved, full double series convergence proved for c > 6.125.

2. **dg-shifted Yangian computational frontier** -- Explicit presentations of Y(sl_2), Y(sl_3), Y(sl_4) from the ordered bar complex. Root quadrilateral at filtration 3. Three Koszul duality operations table. Two-colour architecture. Curvature of the Yangian. All in working_notes.tex sections 14807-15268.

3. **E_1/E_infty convention crisis** -- An LLM agent incorrectly concluded "vertex algebras are NOT E_infty" based on reading BD Chapter 4. This caused cascading edits to ordered_associative_chiral_kd_core.tex with contaminated parentheticals. The author corrected the convention, the contaminated file was reverted, and 24 new anti-patterns (AP35-AP58) were added to CLAUDE.md to prevent recurrence.

Net result: ~22,600 lines of additions across 36 files, but one critical revert (ordered_associative_chiral_kd_core.tex) destroyed ~200 lines of GOOD content alongside the bad parentheticals.

---

## B. THE E_1/E_infty FRAMEWORK (Author's Settled Convention)

**This is the single most important thing to internalise before editing.**

- **E_infty-chiral = LOCAL.** ALL vertex algebras are E_infty-chiral: Kac-Moody, Virasoro, Heisenberg, W-algebras, lattice VOAs. The factorization structure is on UNORDERED configurations with Sigma_n-equivariance. OPE poles are compatible with E_infty.

- **E_1-chiral = NONLOCAL.** A genuinely new concept introduced in Vol II. Factorization on ORDERED configurations, no Sigma_n symmetry. Examples: Etingof-Kazhdan quantum vertex algebras, Yangians, quantum groups.

- **The discriminant is PROVENANCE of R(z), not its value.** Both E_infty algebras with poles and E_1 algebras can have R(z) != tau. For E_infty, R(z) is DERIVED from the local OPE (monodromy of a flat connection). For E_1, R(z) is INDEPENDENT INPUT.

**Three-tier R-matrix picture (all within or adjacent to E_infty):**

| Tier | R(z) | E_infty? | Examples |
|------|-------|----------|----------|
| (i) Pole-free commutative | R(z) = tau (flip) | Yes (subclass) | Polynomial algebras, D-schemes |
| (ii) Vertex algebras with poles | R(z) != tau, derived from OPE | Yes | H_k, V_k(g), Vir_c |
| (iii) Genuinely E_1 | R(z) != tau, independent input | No | EK quantum vertex algebras, Yangians |

**Critical negatives:**
- NEVER say "E_infty means no OPE poles"
- NEVER say "Heisenberg/Virasoro/KM is not E_infty"
- NEVER add parentheticals like "(commutative chiral algebras in the sense of BD -- those whose chiral product has no OPE singularities)" -- this restricts E_infty to the pole-free subclass
- NEVER edit E_1/E_infty language without explicit author confirmation (AP50)

Read CLAUDE.md AP35-AP58 in full before touching any E_1/E_infty text.

---

## C. WHAT SURVIVES (Manuscript Insertions That Passed Audit)

All of these are in the unstaged working tree diff (git diff HEAD). They should be KEPT:

1. **spectral-braiding-core.tex** (+1164 lines): Arnold/IHX remark, explicit CYBE verification for Virasoro r-matrix, Virasoro r-matrix antisymmetry explanation (d log absorption), SC^{ch,top} operations on Y(sl_2) with full proof (explicit FM-integral formulas, collision residue, ordered bar structure, spectral Kohno connection).

2. **ordered_associative_chiral_kd_frontier.tex** (+546 lines): conj:FG-shadow PROMOTED TO THEOREM with complete proof via stratification spectral sequence + pole non-increase theorem. Pole non-increase theorem proved in full (Theorem thm:pole-non-increase). FG-shadow convergence corollary.

3. **3d_gravity.tex** (+1989/-~200 lines): Brown-Henneaux arithmetic correction (factor-of-2 fix in Q), Beilinson-fortified proofs.

4. **thqg_celestial_holography_extensions.tex** (+564 lines): Celestial holography extension material.

5. **thqg_perturbative_finiteness.tex** (+524 lines): Partial fraction sign correction ((-1)^{n+1} -> (-1)^n), residue correction, intermediate expansion coefficient fix.

6. **working_notes.tex** (+12085 lines): All session notes including the dg-shifted Yangian presentations, Polyakov synthesis, Wick rotation architecture, chromaticity section.

7. **rosetta_stone.tex** (+218 lines), **w-algebras-stable.tex** (+238 lines), **examples-worked.tex** (+78 lines), **examples-complete-conditional.tex** (+204 lines): Example chapter expansions.

8. **CLAUDE.md** (+61 lines): AP35-AP58 anti-patterns.

9. **Various core/frontier files**: anomaly_completed_core (+40), celestial_boundary_transfer_core (+68), log_ht_monodromy_core (+39), modular_pva_quantization_core (+188), and their frontier counterparts.

10. **compute/audit/linear_read_notes.md** (+3544 lines): Findings register from adversarial audits.

---

## D. WHAT WAS REVERTED AND NEEDS RE-INSERTION

**ordered_associative_chiral_kd_core.tex** was reverted to HEAD via `git checkout HEAD -- file`, destroying BOTH contaminated parentheticals AND good mathematical content. The file currently has ZERO diff from the last commit.

**Good content that was destroyed and must be re-inserted:**

1. **sl_3 computation** -- Explicit adjacent-root structure, triangle sectors, the coefficient 1/2 from FM-integral. This is mathematically correct and should be restored from working_notes.tex section 14906-14925.

2. **Two-colour Koszul duality architecture Remark** -- The table showing closed-colour (Verdier/Ran -> chiral algebra A^!_ch), open-colour (linear duality -> Yangian A^!_line), pure E_1 (classical -> bare A_infty). Correct content, documented in working_notes.tex section 15165-15253.

3. **Small vs homotopy Koszul distinction** -- The distinction between quadratic Koszulity (combinatorial, on generators/relations) and homotopy Koszulity (homotopical, on the operad). Important for the recognition theorem narrative.

4. **Gravitational Yangian Remark** -- The gauge-gravity complexity dichotomy table (gauge: m_k = 0 for k >= 3, nontrivial Yangian; gravity: m_k != 0 for all k, strictly primitive). From working_notes.tex section 14942-14958.

5. **Annular bar differential** -- The bar differential on the annulus and its relation to genus-1 curvature.

**How to re-insert:** Use the working_notes.tex content (sections 14807-15268) as the source. Insert into ordered_associative_chiral_kd_core.tex using targeted Edit operations. Do NOT add any parenthetical glosses on E_infty (the pre-existing text in the file is correct). Do NOT batch-propagate -- insert one piece at a time and verify.

---

## E. WHAT WAS PROVED

1. **conj:FG-shadow promoted to Theorem** (in ordered_associative_chiral_kd_frontier.tex). The proof is complete: stratification spectral sequence (Construction constr:commutator-stratification) + E_1 identification (Proposition prop:stratification-E1) + convergence (Corollary cor:FG-shadow-convergence) + pole non-increase (Theorem thm:pole-non-increase). This is a genuine mathematical result.

2. **Explicit SC^{ch,top} operations on Y(sl_2)** (in spectral-braiding-core.tex). The ordered bar complex of V_k(sl_2) presented in full: degree 1 (generators), degree 2 (9 generators, bar differential extracting Lie bracket), degree 3 (d^2 = 0 is YBE). The collision residue r(z) = hbar*Omega/z. The spectral Kohno connection. The RTT presentation derived from the bar complex.

3. **Gravitational CYBE verification** (in spectral-braiding-core.tex, Proposition prop:virasoro-cybe). Full term-by-term verification that the Virasoro r-matrix satisfies the CYBE, via reduction to the PVA Jacobi identity.

4. **BP complementarity: c_sum = 2** (in working notes). The Bershadsky-Polyakov algebra has c + c' = 2, not 76. The ghost formula (28) is wrong for non-principal W-algebras.

5. **Genus-g curvature-braiding uniformity** (in working notes). The Coisson bracket c_0 is genus-independent; the curvature-braiding dichotomy persists at all genera via the Bergman kernel.

6. **Full double series convergence** (in working notes). For c > 6.125, the bidisk convergence is proved with explicit radius bounds.

---

## F. AGENT OUTPUTS WORTH MINING

There are 43 agent worktrees in `.claude/worktrees/`. The most recent (2 April) contain:

- **Double-bar computation** (partially correct): The open-colour double bar B^{ord}(Y_hbar(sl_2)) recovers the current algebra sl_2[t] WITHOUT central extension. The Lie bracket recovery is RIGHT. The claim that the central extension is invisible to the open-colour double bar is RIGHT. But any claim that the full double bar recovers the affine algebra including central extension via the OPEN colour is WRONG -- the central extension requires the CLOSED colour (Verdier/Ran duality).

- **FG-shadow proof** (correct): The stratification spectral sequence approach to proving conj:FG-shadow. This has been incorporated into ordered_associative_chiral_kd_frontier.tex.

- **Explicit FM-integral formulas** (correct): The FM_k(C) integrals for collision residues at k = 2, 3, 4. The 1/2 coefficient at k = 3 (sl_3 triangle) from the beta integral. The 1/3 coefficient at k = 4 (sl_4 quadrilateral).

- **Heisenberg R-matrix** (corrected): R(z) = exp(k*hbar/z), NOT R = 1. The initial agent computation was wrong (forgot AP19: d log absorbs a pole). Corrected computation is in AP41.

- **Chromatic frontier** (quarantined): 1710 lines of material with deep mathematical errors flagged. The errors are documented in working_notes.tex section 15349-15445 (five specific anti-patterns). The correct parts (smashing localisation, four-stage pipeline, formal group identification) are extracted into working_notes.tex section 15288-15445.

---

## G. THE DANGER ZONES

**Read CLAUDE.md lines 130-178 (AP35-AP58) before ANY editing session.** These encode every error pattern discovered on 2 April. The most dangerous:

- **AP39**: NEVER equate "E_infty-chiral" with "no OPE poles"
- **AP42**: NEVER introduce parenthetical glosses that RESTRICT a defined term to a subclass
- **AP47**: NEVER treat an agent's literature claim over the author's explicit statement
- **AP49**: NEVER oscillate between conventions within a session
- **AP50**: NEVER edit E_1/E_infty language without explicit author confirmation
- **AP52**: NEVER revert an entire file when only specific lines are contaminated
- **AP54**: NEVER propagate a correction to multiple files before verifying it

**The Heisenberg R-matrix trap (AP41):** Multiple agents on 2 April computed R = 1 for Heisenberg by forgetting AP19 (the d log kernel absorption). The correct R-matrix is R(z) = exp(k*hbar/z). Always apply d log absorption before computing collision residues.

---

## H. RECOMMENDED FIRST ACTIONS FOR NEXT SESSION

1. **Read CLAUDE.md AP35-AP58 in full.** Non-negotiable.

2. **Verify the surviving insertions build.** Run `pkill -9 -f pdflatex 2>/dev/null || true; sleep 2; make` and check for compilation errors.

3. **Re-insert the five pieces of good content into ordered_associative_chiral_kd_core.tex.** Source: working_notes.tex sections 14807-15268. Insert surgically, one at a time. Do NOT add parenthetical glosses on E_1/E_infty. The existing text in that file is correct; only ADD new content.

   Priority order:
   - (a) Two-colour Koszul duality architecture Remark (working_notes 15165-15253)
   - (b) Gravitational Yangian Remark / gauge-gravity dichotomy table (working_notes 14942-14958)
   - (c) sl_3 adjacent-root computation (working_notes 14906-14925)
   - (d) Small vs homotopy Koszul distinction
   - (e) Annular bar differential

4. **Mine the agent worktrees** for any results not yet incorporated. Key targets:
   - The double-bar analysis (partially correct, needs careful extraction)
   - Any remaining FM-integral verifications
   - The chromatic quarantine material (correct parts only)

5. **Commit the surviving changes.** The working tree has ~22,600 lines of good additions that should be committed before any further editing.

6. **The open frontier from 2 April** (documented in working_notes 15796-15828):
   - Full Virasoro convergence beyond scalar lane
   - BP linear coproduct (first matter-coupled gravity example?)
   - Soft graviton p >= 3 (quintic shadow S_5 = -48/[c^2(5c+22)])
   - Closed-form P_k or spectral curve for the graviton self-energy resolvent

---

## FILE STATE SUMMARY

| File | Status | Notes |
|------|--------|-------|
| ordered_associative_chiral_kd_core.tex | COMMITTED (3e0587c) | All 5 pieces re-inserted: (a) two-colour architecture Remark, (b) gravitational Yangian / gauge-gravity dichotomy, (c) sl_3 adjacent-root computation, (d) small vs homotopy Koszul distinction, (e) annular bar differential. Plus double-bar proposition (B^{ord,ch}(Y) recovers current algebra without central extension). Correct E₁/E_∞ language throughout. |
| ordered_associative_chiral_kd_frontier.tex | +546 lines | conj:FG-shadow promoted to theorem. KEEP. |
| spectral-braiding-core.tex | +1164 lines | SC operations on Y(sl_2), CYBE, Arnold/IHX. KEEP. |
| 3d_gravity.tex | +1989 lines | Brown-Henneaux fix. KEEP. |
| thqg_perturbative_finiteness.tex | +524 lines | Sign corrections. KEEP. |
| working_notes.tex | +12085 lines | All session notes. KEEP. |
| CLAUDE.md | +61 lines | AP35-AP58. KEEP. |
| compute/audit/linear_read_notes.md | +3544 lines | Audit register. KEEP. |
| 28 other files | Various small additions | Core/frontier expansions. KEEP. |

---

## FULL CATALOG OF AGENT OUTPUTS AND UNWRITTEN RESULTS

### Computations produced by the 24-agent dg-shifted Yangian swarm (not yet in manuscript):

1. **Complete sl₂ Yangian from bar complex** (C1): All generators, RTT 16-component expansion, Gauss decomposition, twisted coproduct Δ_z, evaluation modules, PBW basis with Hilbert series ∏(1-q^d)^{-3}. IN WORKING NOTES (lines 14812-14904).

2. **sl₃ adjacent-root computation** (C2): Triangle sectors T^L_{12} = E_{13}^(1)E_{21}^(2)E_{32}^(3) in matrix units, coefficient 1/2 verified, Serre from g_{2α₁+α₂}=0. IN MANUSCRIPT (rosetta_stone.tex) AND WORKING NOTES.

3. **A₃ root quadrilateral** (C3): Coefficient 1/3 = ∫₀¹(1-t)²dt, BCH verification, root-space one-dimensionality for sl₄. IN WORKING NOTES (lines 14927-14939).

4. **RTT presentation from d²=0** (C4): Complete derivation showing RTT = bar d²=0 on degree-3. Explicit 16 RTT relations for sl₂. IN WORKING NOTES.

5. **Drinfeld new realization** (C5): All relations (D1)-(D6) traced to bar complex identities. Explicit for sl₂ and sl₃. IN WORKING NOTES.

6. **Explicit R-matrices** (C6): Yang R-matrix R(u)=uI+P verified for sl₂ (4×4) and sl₃ (9×9). YBE verified by Python. IN WORKING NOTES.

7. **Free-field realization** (C7): βγ nilpotent kernel Θ²=0, exact transport T(u,u₀) = 1+Θ·log(...), Wakimoto embedding. IN WORKING NOTES.

8. **PBW basis** (C8): gr(Y_dg(g)) ≅ U(g[t]), Hilbert series, shuffle algebra connection. IN WORKING NOTES.

9. **Twisted coproduct** (C9): Explicit Δ_z on Drinfeld generators, infinitesimal coassociativity verified. IN WORKING NOTES.

10. **Evaluation modules** (C10): ev_a construction, Clebsch-Gordan, Drinfeld polynomials, C_line ≃ A!-mod. IN WORKING NOTES.

11. **Kac-Moody frontier** (C11): Three-tier stratification (sl₂-hat trivial / affine rank≥2 conjectural / hyperbolic failure). CONJECTURE IN MANUSCRIPT (dg_shifted_bridge.tex).

12. **Spectral Drinfeld class** (C12): Explicit vanishing at filtrations 2,3,4 for types A, B₂, G₂. IN WORKING NOTES.

### Theoretical results from the 12-agent theory swarm:

13. **Koszul duality → Yangian mechanism** (T1): Three functors on B(A), E₁ vs E_∞ distinction, the recognition theorem. IN WORKING NOTES.

14. **Swiss-cheese E₁ wing** (T2): Why E₁ on open side (topology of ℝ forces it), propagator factorization. IN WORKING NOTES.

15. **Configuration space geometry** (T3): FM_k(ℂ)×Conf_k(ℝ) boundary strata, associahedron, transfer via HPL. IN WORKING NOTES.

16. **Strictification mechanism** (T4): Root-space 1-dim + Jacobi collapse + BCH=1/n = no obstruction. IN WORKING NOTES.

17. **KZ monodromy and quantum groups** (T5): Derived additive KZ → classical KZ via one-loop collapse, Laplace transform as spectral duality. IN WORKING NOTES.

18. **Gravitational Yangian** (T6): Gauge-gravity (m,Δ) dichotomy, opposite corners of complexity plane, DS as transport. IN MANUSCRIPT (3d_gravity.tex) AND WORKING NOTES.

19. **Spectral braiding** (T7): Arnold cancellation = IHX, Jones polynomial from bar complex (proved for affine lineage). IN MANUSCRIPT (spectral-braiding-core.tex).

20. **DS reduction transport** (T8): DS preserves Koszulness, destroys formality. Resolvent tree formula. IN MANUSCRIPT (w-algebras-stable.tex) AND WORKING NOTES.

21. **Shifted quantum groups** (T9): FT/BFN/COZZ comparison, three-language spine. IN WORKING NOTES.

22. **Factorization quantum groups** (T10): Latyntsev comparison, the three-stage pipeline. IN WORKING NOTES.

23. **Categorical enhancement** (T11): C_op as primitive, Morita invariance, 2-categorical structure. IN WORKING NOTES.

24. **Genus ≥ 1 curved Yangian** (T12): Curved A_∞ from κ·ω_g, genus-1 KZB, Fay trisecant, curvature-braiding dichotomy. IN WORKING NOTES.

### Two-colour Koszul duality architecture (from the 14-agent resolution swarm):

25. **Three bar complexes**: FG bar (zeroth product) ≠ full symmetric bar (all products, Σ_n-coinvariants) ≠ ordered bar (all products, no Σ_n). IN CLAUDE.md (AP37) AND WORKING NOTES.

26. **Verdier vs linear duality**: Closed colour uses D_Ran (produces chiral algebra), open colour uses H*(-)^∨ (produces Yangian). IN MANUSCRIPT (rem:two-colour-architecture).

27. **R-matrix as cross-colour datum**: B^ch = (B^ord)^{R-Σ_n}. For pole-free E_∞: R=τ. For E_∞ with poles: R≠τ but derived from local OPE. For E₁: R independent input. IN MANUSCRIPT AND WORKING NOTES.

28. **FG-shadow theorem PROVED**: Commutator filtration compatibility. Key insight: pole-order filtration fails for Virasoro but PBW filtration works. IN MANUSCRIPT (frontier.tex).

29. **Explicit SC^{ch,top} operations**: μ_{0,2} (Yangian product), μ_{2,0} (Sklyanin bracket), μ_{1,1} (Casimir adjoint), μ_{1,2} (RTT), general kernel formula. IN MANUSCRIPT (spectral-braiding-core.tex).

30. **Double bar computation**: B^{ord,ch}(Y_ℏ(sl₂)) recovers sl₂[t] (current algebra WITHOUT central extension). Central extension is closed-colour datum, invisible to open-colour double bar. IN WORKING NOTES.

### Results that STILL NEED to be written into the manuscript:

- The gravitational Yangian Remark (was in ordered_kd_core, now only in dg_shifted_bridge — consider whether it belongs in both)
- The annular bar differential construction (was in ordered_kd_core, now lost — needs re-insertion from working notes)
- The double-bar computation as a Computation environment (currently only in working notes)
- Mining the 100+ agent outputs for additional results

### Anti-patterns (AP35-AP58) — READ BEFORE ANY EDIT:

The E₁/E_∞ crisis produced 24 anti-patterns. The critical ones:
- AP35: E_∞ = LOCAL = ALL vertex algebras. E₁ = NONLOCAL.
- AP36: R(z)≠τ does NOT imply E₁. Discriminant is provenance.
- AP41: Heisenberg R(z) = exp(kℏ/z), NOT trivial. d log absorption!
- AP50: NEVER edit E₁/E_∞ without author confirmation.
