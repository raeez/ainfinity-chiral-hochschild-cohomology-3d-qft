---
name: vol2-build-surface
description: Run isolated local LaTeX builds, targeted tests, and log checks for affected mathematical changes.
---

# Build and test evidence

Read `Makefile` and `platonic/PLATONIC_LEDGER.md` to select the requested source.
The current integrated manuscript is `platonic/main.tex`, built by `make platonic` into `out/platonic.pdf`.
Run that target only in the assigned isolated worktree. It writes auxiliaries inside that worktree's `platonic/` directory.
Do not run concurrent builds in the same worktree. `MKD_BUILD_NS` does not isolate the platonic target.

```bash
make platonic
```

Inspect source freshness, decisive logs, output text, and changed rendered pages. Existing PDFs and a zero exit do not prove fresh output.
`make fast` and `scripts/build.sh` compile the legacy root `main.tex`. They do not verify the current integrated manuscript.
Use them only for explicit legacy-source work, with a unique `MKD_BUILD_NS` for that task.
Do not substitute a legacy build when a current source is missing or fails.

For compute changes, use the affected `python3 -m pytest` slice with the repository's available dependencies.
Metadata regeneration applies only when the change affects generated claim indexes.

Local builds after coherent changes need no repeated approval. Broader checks follow the changed dependencies and observed failures.
Release, iCloud, and publication targets have external effects. Run them only with the relevant existing authorization.
Never terminate processes by executable name. Stop only a process handle or PID launched and still owned by this task.
Before signaling a PID, confirm its identity and ownership. Prefer graceful termination.

Classify failures before repair:

- Fatal LaTeX errors require source or dependency diagnosis.
- Undefined references matter after stable reruns. Check external references against their intended source.
- First-pass warnings are provisional until auxiliary files stabilize.
- An interrupted build needs a fresh task namespace or repair of its own artifacts before comparison.
- Oracle disagreement needs independent mathematical and convention checks. Do not copy engine outputs into tests.

Report the command, source revision, result, and remaining warnings. A clean build verifies rendering and references, not mathematical truth.

## Changed ProvedHere claims

Preserve or add the disjoint-route witness in `compute/tests/test_*_iv.py` when changing a `ProvedHere` claim.
Tests that claim independent theorem verification use `@independent_verification` from `compute/lib/independent_verification.py`.
Name the claim label, derivation sources, verification sources, and the rationale for their independence.
Run the affected witness tests, then `make verify-independence` for registry tautology and orphan checks.
Read `compute/scripts/audit_independent_verification.py` when interpreting its result.
Its scanner covers legacy `chapters/`, `appendices/`, and specified note sources. It does not scan `platonic/` automatically.
For current integrated claims, explicitly match the witness to the actual theorem label and proof in the selected input graph.
A zero audit exit permits coverage gaps. It cannot establish complete witness coverage or prove a theorem.
Report uncovered or unscanned claims precisely, and retain the missing witness as an obligation rather than changing theorem status to pass the gate.
