---
name: vol2-deep-audit
description: Use when the user asks to audit, review, red-team, falsify, pressure-test, or run a findings-first mathematical pass on a Vol II theorem, chapter, proof, or live surface. Do not use for straightforward fix-only requests where the rectification skill is the better fit.
---

# Vol II Deep Audit

This is the findings-first counterpart to rectification.

## Inputs

- exact target file, theorem, definition, or live surface if given
- otherwise the smallest live surface implied by the prompt

## Audit posture

- Read before editing.
- Default output is prioritized findings, not prose reassurance.
- Treat every strong claim as false until a local proof, citation, or computation survives hostile checking.
- If the user asks for both audit and repair, start with the audit pass, then hand off to `$vol2-beilinson-rectification` for fixes.

## Mandatory loop

1. Start with a short progress update naming the target and first verification step.
2. Track substantial work with an available planning tool or a concise written plan.
3. Read:
   - target file
   - neighboring context
   - selected active `\input` graph from `platonic/main.tex` (verify `Makefile` and the integrated ledger)
   - current dirty diff
   - relevant compute/tests or build logs
4. Select relevant checks in one audit. Repeat only where new evidence warrants:
   - `RED`: proof logic, hypotheses, signs, formulas, scope, hidden conditionality
   - `BLUE`: theorem/proof/status drift, label collisions, duplicate formulations, build/test collisions, cross-volume contradictions
   - `GREEN`: missing definitions, dangling references, unsupported objects, missing lemmas, places where the text overstates what is actually built
5. Record actionable findings in `compute/audit/linear_read_notes.md`.
6. Report findings ordered by severity and anchored to exact locations.
7. If no findings remain at `MODERATE+`, say so explicitly and state residual verification gaps.

## Reporting standard

- Findings first.
- Prioritize bugs, risks, regressions, and missing tests over summaries.
- If no finding is strong enough to survive reread, say "no findings" rather than inventing friction.
- If a claim is only conditionally or computationally supported, say that directly.

## Tool discipline

- Batch independent reads, searches, and checks with the available host tools.
- Use `exec_command` for `git diff`, `pytest`, `make fast`, and focused shell inspection.
- Do not start patching until the audit surface is understood.
