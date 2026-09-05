---
name: vol2-cross-volume-propagation
description: Use after any mathematical wording, status, or formula change that may appear elsewhere in Vol II, Vol I, Vol III, superseded split files, notes, or compute layers. Do not use for isolated edits that cannot plausibly propagate.
---

# Vol II Cross-Volume Propagation

This skill exists to stop local truth from coexisting with global drift.

## Search surface

After a load-bearing change, inspect relevant live occurrences:

- active Vol II chapter files
- active Vol II appendices
- historical references needed to locate active duplicates
- `~/chiral-bar-cobar`
- `~/calabi-yau-quantum-groups` when the bridge is genuinely cross-volume
- `compute/` and `compute/tests/`

## Propagation checklist

1. Search exact phrase variants and symbolic variants.
2. Check theorem statements, summaries, examples, remarks, introductions, and appendices.
3. Check build-facing surfaces:
   - label names
   - references
   - citation claims
4. Check compute-facing surfaces:
   - function docstrings
   - tests
   - hardcoded expected values
5. Update verified live duplicates within assigned paths. Report other owners' exact dependent paths.

## Convention alert

- Vol I uses OPE modes.
- Vol II uses lambda-brackets with divided powers.
- Vol III may use motivic or categorical normalizations.

Never compare formulas across volumes without converting conventions first.

## Stop condition

Do not call the propagation done until:
- the active surface within assigned paths is consistent;
- affected paths outside scope are reported to their integration owner;
- historical evidence remains unchanged and is distinguished from active claims;
- tests or log checks no longer contradict the new statement.

Read other repositories only where the changed claim depends on them. This skill grants no cross-repository write authority.
Preserve archive snapshots. Report active misleading references to their accountable owner.
