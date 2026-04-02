# Adversarial Audit: 3d_gravity.tex, Corollary cor:gauge-gravity-dichotomy

**Date:** 2026-04-02  
**File:** `/Users/raeez/chiral-bar-cobar-vol2/chapters/connections/3d_gravity.tex`  
**Target:** Corollary `\label{cor:gauge-gravity-dichotomy}` (lines 157-194)  
**Status:** AUDIT COMPLETE. ISSUES FOUND: 3 (1 HIGH, 2 MEDIUM).

---

## Examined Claims

1. **Table claim (line 171):** "Gravity: m_k ≠ 0 ∀k (non-formal)"
2. **Corollary body (line 172):** "strictly primitive Δ_z"
3. **Corollary body (lines 180-182):** "ghost-number obstruction (|h|=-1 in ghost degree) kills all non-primitive coproduct terms"
4. **Proof chain:** Claims sourced from `prop:pole-order-classification`, `thm:ds-hpl-transfer`, `thm:gravitational-primitivity`

---

## ISSUE 1: HIGH SEVERITY — CLAIM STATUS MISMATCH ON INFINITUDE OF m_k TOWER

### Description

The table entry "Gravity: m_k ≠ 0 ∀k (non-formal)" is claimed to follow from `Proposition~\ref{prop:pole-order-classification}` (line 187).

However, `prop:pole-order-classification` (lines 135-152) does not rigorously prove m_k ≠ 0 for ALL k. The proof only says:

> "the Stasheff identity then forces m_4, m_5, ... in perpetuity" (line 152)

This is qualitatively correct but not a rigorous proof of infinitude.

### Where the Actual Proof Lives

Lines 517-521 and line 681 cite `Proposition~\ref{prop:vir-truncation}` to support the infinitude claim:

> "The Ainf structure does not truncate: m_k ≠ 0 for all k ≥ 3 and generic c. ... (Proposition~\ref{prop:vir-truncation})."

**grep result:** `\label{prop:vir-truncation}` is NOT in 3d_gravity.tex. It is defined in `/Users/raeez/chiral-bar-cobar-vol2/chapters/examples/w-algebras-virasoro.tex`.

### Status Problem

The proposition in w-algebras-virasoro.tex is marked:

```
\begin{proposition}[Structure of Virasoro A_∞ Operations; \ClaimStatusConditional]
```

It begins: "Assume the Khan–Zeng Virasoro realization satisfies Theorem~\ref{thm:physics-bridge}."

**Contradiction:** Corollary `cor:gauge-gravity-dichotomy` is marked `\ClaimStatusProvedHere` (line 157), yet its proof depends on a **conditional** proposition (marked `\ClaimStatusConditional`). A proved result cannot depend on a conditional ingredient.

### Action Required (Choose One)

**(A)** Promote `cor:gauge-gravity-dichotomy` to `\ClaimStatusConditional` with explicit statement of the condition (that Khan–Zeng Virasoro realization satisfies `thm:physics-bridge`).

**(B)** Provide a self-contained unconditional proof of m_k ≠ 0 ∀k in 3d_gravity.tex (using the wheel-diagram argument without dependence on the conditional proposition).

**(C)** Move `prop:vir-truncation` to 3d_gravity.tex and upgrade it from conditional to proved.

---

## ISSUE 2: MEDIUM SEVERITY — GHOST-NUMBER OBSTRUCTION ATTRIBUTION ERROR

### Description

Corollary body (lines 180-182) claims:

> "the ghost-number obstruction (|h|=-1 in ghost degree) kills all non-primitive coproduct terms"

Corollary proof (line 192-193) attributes this to:

> "for gravity Δ_z is primitive (Theorem~\ref{thm:gravitational-primitivity}). The DS transport is Theorem~\ref{thm:ds-hpl-transfer}."

The wording "ghost-number obstruction... kills all non-primitive coproduct terms" should refer to the detailed argument, but the proof only cites the high-level existence statement via `thm:ds-hpl-transfer`.

### Where the Detailed Argument Lives

Theorem `thm:ds-hpl-transfer` (lines 1241-1320):
- States that HPL transfer produces transferred coproduct that is an A∞-algebra morphism.
- Does NOT explain the ghost-number obstruction mechanism.

Theorem `thm:gravitational-primitivity` (lines 1576-1601):
- Contains the complete ghost-number obstruction proof in part (iii) of the proof (lines 1663-1669).
- Rigorously shows that source-side trees have ghost number +1, target trees have ghost number -1, and projection `p ⊗ p` annihilates both.

### Problem

The ghost-number obstruction argument is in `thm:gravitational-primitivity`, not `thm:ds-hpl-transfer`. The Corollary proof is imprecise about attribution.

### Action Required

Revise proof text to clarify:

> "for gravity Δ_z(x) = τ_z(x) ⊗ 1 + 1 ⊗ x is strictly primitive (Theorem~\ref{thm:gravitational-primitivity}), where primitivity is enforced by the ghost-number obstruction (|h| = -1) that kills all non-primitive terms."

This makes clear that the obstruction argument is part of the primitivity theorem, not DS transfer.

---

## ISSUE 3: MEDIUM SEVERITY — UNDERLYING THEOREM HAS DISCLOSED GAPS

### Description

Corollary proof cites `Theorem~\ref{thm:ds-hpl-transfer}` (line 193) for the DS transport claim.

Theorem `thm:ds-hpl-transfer` (lines 1241-1320) states HPL transfer produces transferred A∞ products, coproduct, and r-matrix.

However, the theorem includes a remark (lines 1322-1334, `rem:ds-hpl-honest-gaps`) disclosing:

**(a)** Full sign verification of transferred products at arity ≥ 3 is beyond the proved core.

**(b)** Verification that transferred structure satisfies dg-shifted Yangian axioms (translation compatibility, rationality of r-matrix) as opposed to general A∞ structure is beyond the proved core.

The remark also notes: "The coproduct vanishing, originally listed as a separate gap, is resolved at all arities by Theorem~\ref{thm:gravitational-primitivity}."

### Assessment

The gaps are honestly disclosed. The coproduct part (the piece used in Corollary) is resolved by `thm:gravitational-primitivity`. However, the sign verification gap (a) affects the full transferred product structure at arity ≥ 3.

Since Corollary `cor:gauge-gravity-dichotomy` claims the full A∞ structure is non-trivial (m_k ≠ 0 ∀k), the sign verification gap is relevant to the claim's status.

### Action Required

No immediate fix needed if the signs are expected to be verified elsewhere or if the Corollary is promoted to `\ClaimStatusConditional` (as suggested in Issue 1). However, document the dependency: if sign verification is incomplete in `thm:ds-hpl-transfer`, the Corollary's claim about the full tower should flag this dependency.

---

## ISSUE 4: LOW-MEDIUM SEVERITY — ENVIRONMENT CLASSIFICATION UNCLEAR

### Description

The result is marked `\begin{corollary}` (line 157), suggesting it is a simple consequence of the preceding proposition. However:

- The dichotomy is presented as a major structural result ("Gauge-gravity complexity dichotomy"; caption and index).
- The proof depends on THREE independent major theorems: `prop:pole-order-classification`, `thm:ds-hpl-transfer`, and `thm:gravitational-primitivity`.
- It is a primary result of Part VI (Three-Dimensional Quantum Gravity), not a subsidiary consequence.
- It is heavily cited and indexed.

### Observation

While technically a corollary (it follows logically from the preceding proposition + two additional theorems), the breadth and significance suggest Theorem status would be clearer.

### Recommendation (Non-Critical)

Either:
- Promote to `\begin{theorem}[Gauge-gravity complexity dichotomy; ...]`, or
- Retain as `\begin{corollary}` but explicitly state in the proof that it combines three major results: "The following corollary assembles Propositions~\ref{prop:pole-order-classification}, Theorems~\ref{thm:ds-hpl-transfer} and~\ref{thm:gravitational-primitivity}..."

---

## ISSUE 5: VERIFIED — NO ISSUE FOUND

### Claim (2): "strictly primitive Δ_z"

**Status:** Theorem `thm:gravitational-primitivity` exists (line 1578) and rigorously proves strict primitivity.

**Result:** ✓ NO ISSUE.

---

## SUMMARY TABLE

| Issue | Severity | Type | Locations | Status |
|-------|----------|------|-----------|--------|
| 1. m_k ≠ 0 ∀k infinitude claim depends on conditional proposition | HIGH | Proof dependency mismatch | Lines 157, 187, 517-521, 681 | **REQUIRES ACTION** |
| 2. Ghost-number obstruction attribution | MEDIUM | Incomplete proof reference | Lines 180-182, 192-193 | **REQUIRES CLARIFICATION** |
| 3. DS transfer theorem has disclosed gaps | MEDIUM | Underlying incompleteness | Lines 1322-1334 | **DOCUMENT DEPENDENCY** |
| 4. Corollary vs Theorem classification | LOW-MEDIUM | Organizational clarity | Line 157 | Optional improvement |
| 5. Theorem `thm:gravitational-primitivity` existence | — | Verification | Line 1578 | ✓ VERIFIED |

---

## RECOMMENDATIONS (Priority Order)

### 1. URGENT: Resolve m_k ≠ 0 ∀k Status

Choose and implement one of options (A), (B), or (C) above. The status mismatch (Proved claim depending on Conditional proposition) is load-bearing.

### 2. HIGH: Clarify Ghost-Number Obstruction Attribution

Revise Corollary proof to make clear that the ghost-number obstruction argument is part of `thm:gravitational-primitivity`, not `thm:ds-hpl-transfer`. Correct the attribution in the proof (line 192-193).

### 3. MEDIUM: Document Dependencies if Issue 1 is Resolved via Promotion to Conditional

If Corollary is promoted to `\ClaimStatusConditional`, explicitly state the condition: "Assume Khan–Zeng Virasoro realization satisfies Theorem~\ref{thm:physics-bridge}."

### 4. OPTIONAL: Promote to Theorem for Clarity

If appropriate, promote `cor:gauge-gravity-dichotomy` to `\begin{theorem}[...]` to reflect its structural importance.

---

## Files Referenced

- Primary: `/Users/raeez/chiral-bar-cobar-vol2/chapters/connections/3d_gravity.tex` (lines 157-194, 1576-1669)
- Secondary: `/Users/raeez/chiral-bar-cobar-vol2/chapters/examples/w-algebras-virasoro.tex` (contains `prop:vir-truncation`, marked conditional)
- Main build: `/Users/raeez/chiral-bar-cobar-vol2/main.tex`

---

**Audit conducted:** 2026-04-02 by adversarial deep read.  
**All labels verified via grep and Read tool.**  
**No claims taken on faith; all assertions traced to source.**
