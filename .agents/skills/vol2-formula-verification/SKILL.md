---
name: vol2-formula-verification
description: Use when verifying or repairing a formula, coefficient, invariant, table entry, generating function, spectral sequence page, test oracle, or convention bridge in Vol II. Do not use for purely prose edits with no mathematical payload.
---

# Vol II Formula Verification

Every load-bearing formula should survive at least three independent attacks.

## Verification paths

Use at least three when feasible:

1. direct computation from the definition
2. alternative but equivalent formula
3. limiting or special case
4. symmetry, duality, or reduction
5. cross-family or functorial consistency
6. literature comparison with normalization check
7. degree, weight, sign, or units analysis
8. numerical spot check
9. operadic or factorization consistency
10. descent consistency to the PVA or classical shadow

## Vol II convention traps

- The bar kernel absorbs one pole: `r(z)` has pole order one less than the OPE.
- Vol II lambda-brackets use divided powers: the order-`n` coefficient is `a_(n)b / n!`.
- Desuspension lowers degree.
- `kappa` is not automatically `S_2`.
- `eta(q)` includes `q^(1/24)`.

## Workflow

1. State the exact formula and its scope.
2. Lock conventions before computing.
3. Recompute independently; do not repair by pattern.
4. Map the claim to compute code or tests if executable.
5. Add or strengthen the narrowest test that would catch the specific failure mode.
6. Update verified affected copies within scope. Report exact dependencies outside scope.

## Output standard

- If the claim is verified, say by which independent paths.
- If only part is verified, report that part and the exact unresolved obligation. Preserve the theorem target.
- If the formula differs across volumes, explicitly perform the OPE-mode versus lambda-bracket conversion before comparing.
