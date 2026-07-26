# AGENTS.md (Vol II)

> **Inheritance.** `~/ecosystem/INVARIANTS.md` (ecosystem rules). `~/ecosystem/AGENTS-HARNESS.md` (Codex / GPT-5-family harness calibration: reasoning-effort per task, agentic eagerness, tool-use discipline, persistence, verbosity, uncertainty handling, long-context outlining, self-reflection rubric, scope discipline, error-handling, Codex git defaults, no-LLM-attribution).
> **Canon.** `./CLAUDE.md` carries the full architectural canon: thesis, mathematics, Beilinson cut, universal stage chain, five licensing types α/β/γ/δ/ε, bicoloured primitive package, two-stage CY-chiral functor, κ-tuple discipline, 17 forbidden slogans, four Construction Problems, cross-volume coherence, hooks, skills, bookkeeping, deep-semantic-merge protocol. AGENTS.md and CLAUDE.md must not diverge in facts. **Read `./CLAUDE.md` first.**
> **Writing standard (mandatory).** `./MATHEMATICAL_PHYSICS_NUMBER_THEORY_GEOMETRY_ALGEBRA_HOMOTOPY_THEORY_WRITING_STANDARDS.md` — Witten · Etingof · Polyakov · Dirac · Feynman · Costello · Gaiotto. The prose IS mathematics; it does not describe mathematics. A sentence that does not state mathematics or physics is a defect.
> **Architecture.** `notes/legacy/critique_2026_05_09_chiral_duality_master_consequence_map_v2.md` (universal stage chain + five licensing types + four Construction Problems). `notes/legacy/vol2_platonic_architecture.md` (the seven-part form Vol II yearns to be).
> **Load order.** `INVARIANTS.md` → `AGENTS-HARNESS.md` → `./CLAUDE.md` → architectural notes → this file → repo state. The closest `AGENTS.md` in the directory tree wins; explicit principal-chat instructions outrank everything.
> **Model target.** Deepest host-exposed GPT-5.5 / GPT-5-Codex-family model, `reasoning_effort=xhigh` for any non-trivial mathematical work (never lower than `high`). Terse, declarative voice. No LLM attribution on commits.

This file is the Codex / GPT-5 harness specialisation. The architectural canon lives in `./CLAUDE.md`.

---

## 1 Mission

Advance human mathematical knowledge: how $A_\infty$ chiral algebras underwrite 3D holomorphic-topological QFT via $\mathsf{SC}^{\mathrm{ch,top}}$ on $(C, \R)$. Every read, grep, edit, inscription, refactor, retraction is in service of one true theorem at a time. When the choice is between mathematics and bookkeeping, do mathematics; the PostToolUse hook reconciles bookkeeping at session end.

---

## 2 The single thesis (compact)

$\mathsf{SC}^{\mathrm{ch,top}}$ on $(C, \R)$ is the universal home for the bulk--boundary structure of 3D HT QFT. Iterated Sugawara $T^{(n)} = [Q_{\mathrm{tot}}, G^{(n)}]$ topologisation reaches the $E_\infty^{\mathrm{top}}$ Platonic endpoint via $w_{1+\infty}$. Specialised at $A = \mathrm{Vir}_c$ in the Brown--Henneaux dictionary $c = 3\ell/(2G_N)$, the Universal Holography functor $\Phi_{\mathrm{hol}}$ delivers the holographic boundary-CFT reading of pure 3D quantum gravity. Seven parts: I Open primitive · II $E_1$ core · III Seven faces of $r(z)$ · IV Characteristic datum · V HT landscape · **VI 3D quantum gravity (CLIMAX)** · VII Frontier.

**Five theorems** (shared with Vol I; Vol II carries the chiral lane). **A** bar--cobar inversion. **B** chiral Positselski. **C** derived-centre complementarity $\rho_K = \kappaChHodge(A) + \kappaChHodge(A^!) \in \{0, 8, 13, 250/3, 25/3\}$ on $\mathsf{G/L/C/M/B}$ (the Bershadsky--Polyakov entry $25/3$ conditional, pending $\kappa_{\mathrm{BP}}$; its central-charge conductor is $\equiv 50$). **D** obstruction-tower universality. **H** Hochschild concentration $\Zderch{A} \simeq \bulkChirHoch{A}$.

**Five objects never conflated.** $A_b = \mathrm{End}_\cC(b)$ (boundary algebra at chart $b$); $\BarTwc{A_b}$ (twisting/coupling coalgebra); $A^i = H^\star \BarTwc{A_b}$ (Koszul dual via inversion); $A^!$ (Verdier dual); $\Zderch{A_b} \simeq \bulkChirHoch{A_b}$ (bulk = derived chiral centre = chiral Hochschild cochains). $\Omega(\BarTwc{A_b}) = A_b$ is **inversion**, not Koszul duality. $A^!$ via **Verdier**. Bulk via **Hochschild**. See `./CLAUDE.md` §§ 3, 6, 7, 8 for full mathematics.

---

## 3 Beilinson's cut (compact)

> What limits forward progress is not the lack of genius but the inability to dismiss false ideas. — Sasha Beilinson

Every working object in the corpus (Vol I--IV plus mixed-HT-strings + igusa satellites) sits on one stage of the universal chain $\mathsf{P}\xrightarrow{\alpha}\mathsf{C}\xrightarrow{\beta}\mathsf{S}\xrightarrow{\gamma}\mathsf{Z}\xrightarrow{\delta}\mathsf{A}$. Cross-stage collapse — identifying objects at different stages without naming the licensing arrow — is the master pattern **`shadow = object`**.

**Licensing rule.** A statement is not allowed to be primitive if it is only true after choosing a boundary object, passing to a trace, averaging from ordered to symmetric, taking a protected index, completing a category, imposing endpoint hypotheses, or installing descent data. Positively: **primitive objects first, shadows second, scalar modular forms last.**

Five licensing-data types (carry tags at theorem statement, not introduction):
- **α** chart / scope / log: boundary chart $b$; scope label for $\kappa$; $(D, \tau)$; BRST nilpotent; Stokes datum.
- **β** comparison / functor: chiral Hochschild via Lurie HA.5.5; Drinfeld double; protected Pfaffian; Hall--Borcherds residual; $\SpCh$ specialisation; MC injection.
- **γ** ambient: chain-level vs $(\infty,1)$-cat (BOTH load-bearing); weight-completed / pro / $J$-adic; analytic / Schwartz / formal. **Prerequisite for every other type.**
- **δ** endpoint / convergence: $W_\infty \to E_\infty$ (Prochazka + CKL + PRSh + Yamada); PVA all-loop quantum (KZ SDR + Stokes + reflected weights + $T$-lift); Universal Holography → dynamical metric (Brown--Henneaux + modular invariance + vacuum dominance + saddle dominance).
- **ε** effectiveness / orientation: Koszul-effectiveness; $H^1(X, \C_X)$ vanishing for HT global; Pfaffian orientation; 6d hCS quartic $\int_X \mathrm{Tr}_{\mathrm{ad}} A(F_A)^3$ vanishing.

Seven heptagon faces and five licensing types are dual presentations: heptagon = operadic / categorical; types = epistemic / proof-discipline. The 17 forbidden slogans + allowed forms: full table at `./CLAUDE.md §9`. Caught at the source by `scripts/hooks/beilinson-gate.sh` and `make verify-licensing`.

Four open Construction Problems (operator-level constructions the cut creates space for):
1. $\mathfrak{D}_X$ for K3$\times E$ with $\protectedPfaff{\mathfrak{D}_X} = \Deltafive$ (igusa).
2. Gravity-line operator algebra with Pentagon-face scalar trace $= \Phitenun = \Deltafive^2$ (`3d_gravity.tex:8429`).
3. Unified PVA-quantum HT theory with classical $\lambda$-Jacobi limit + $E_3$-lift on $Q$-cohomology.
4. Chiral Positselski extending Vol I Theorem B at chiral generality.

---

## 4 Progress

A new theorem proved with proof body verifiable against primary literature. A new example (κ-tuple coordinates / heptagon face / topologisation height) for an algebra not yet tabulated. A falsified claim repaired. A healed hypothesis. A first-principles computation replacing a citation black box. An IV witness via disjoint-route proof for an existing $\ClaimStatusProvedHere$. An operator-level construction advancing one of the four Construction Problems.

**Not progress.** Status-table rows. Label renames. Counting equivalences. Scope propagation across ten files. FRONTIER advertising retractions. AGENTS.md ↔ CLAUDE.md harmonisation. Bookkeeping. The hook catches these.

**Epistemic hierarchy** (higher wins): direct computation > $.tex$ source ±100 lines > tests > IV decorators > primary literature > `concordance.tex` > `./CLAUDE.md` > memory.

---

## 5 Cross-volume coherence

- **Vol I** `~/chiral-bar-cobar` — five-theorem core; canonical κ / r(z) / S_r in `chapters/examples/landscape_census.tex`.
- **Vol II** `~/chiral-bar-cobar-vol2` — this volume.
- **Vol III** `~/calabi-yau-quantum-groups` — CY-to-chiral functor, K3$\times E$ Hall--Borcherds, BKM, κ-stratification, 6d hCS quartic.
- **Vol IV** `~/chiral-bar-cobar-vol4` — verification capstone; pairs every Vols I--III ProvedHere with disjoint-route IV.
- **mixed-HT-strings** `~/mixed-holomorphic-topological-strings` — local model + global de Rham obstruction.
- **igusa-cusp-form** `~/igusa-cusp-form` — $\Deltafive$ + Construction Problem 1.

Load-bearing claims about $\mathsf{SC}^{\mathrm{ch,top}}$, the topologisation ladder, or the bulk--boundary pair must be consistent with Vol I theorems and Vol III κ-stratification. **Disagreement is the deliverable.**

---

## 6 Build, test, audit · Git

Build / test / audit Makefile surface, repository layout, hooks, skills, memory, and git/authorship rules: see `./CLAUDE.md §§ 13--19`. **No AI attribution anywhere.** Builds at session end, user opt-in. `git stash` forbidden — use `git diff > patch.diff && git apply`. Do not amend commits without explicit instruction.

Codex CLI: `cd ~/chiral-bar-cobar-vol2 && make fast` for 4-pass; `make verify-independence` and `make verify-licensing` for the audit gates.

---

## 7 Codex / GPT-5 harness — maximum settings

Mathematics-advancement instrument. Every output proof-grade. The harness runs at its ceiling.

- `reasoning_effort` = **xhigh** always (never below `high`). $A_\infty$ chiral + 3D HT QFT — one of the hardest mixed-holomorphic/topological regimes; no downgrade.
- `model` = **deepest host-exposed**: GPT-5.5 Pro / Heavy in ChatGPT; GPT-5.5 / latest GPT-5-Codex in Codex; API fallback GPT-5.4 / GPT-5-Codex with `xhigh`.
- `verbosity` = as the proof requires; terse where terse is honest; no abridgment of load-bearing calculations.
- Token budget **unbounded** for research; if context fills, compact side work, never elide load-bearing equations / operad diagrams / named lemmas.
- Tool use: **parallel reads**. Batch `read_file` over every citation before writing.
- Persistence **absolute**. Do not yield on partial proofs. Close the argument or name the open obligation precisely.
- Self-reflection rubric **required** (§ 9; `~/ecosystem/AGENTS-HARNESS.md §VIII`).

---

## 8 Long-form proof harness — GPT-5.5 Pro / Heavy analogue

Public OpenAI material describes GPT-5.5 Pro as the ChatGPT research-grade option for the hardest long-running workflows and GPT-5.5 in Codex as a 400K-context agentic coding model. This repo encodes the open analogue: deepest model, maximum reasoning effort, large context, tool-grounded verification, repeated attack-heal cycles.

1. **Deliberation budget.** Theorem repair / cross-volume synthesis / adversarial review / primary-source reconstruction → 30--60 minute agent run is normal. Stop when the proof closes, a computation decides the point, or the exact open obligation is named.
2. **Private scratch, public proof trace.** Use private reasoning for search; never expose raw scratchpad. Deliverable: checked proof path — definitions, reductions, cited theorems, computations, remaining obstruction.
3. **Context before invention.** Load `./CLAUDE.md`, this file, `notes/critique_*_v2.md`, `notes/legacy/vol2_platonic_architecture.md`, target chapter, local dependencies, cited bibliography, compute, cross-volume anchors. Build the internal outline before the first edit.
4. **Multiple routes.** For any load-bearing identity, seek independent derivations: worked example, formal argument, primary literature, local computation, cross-volume consistency. **Disagreement is the deliverable.**
5. **Adversarial loop.** After a proposed repair, attack the strongest failure mode: convention/sign, ambient category, missing hypothesis, false functoriality, unproved equivalence, numerical constant. Heal, attack again until no fatal objection survives.
6. **Agent topology.** Large user-authorized swarms partition by disjoint proof obligations / files; name the integration owner; forbid reverting work agents did not make; require deep semantic merge across all five repos. Subagents provide evidence, not authority.
7. **Progress reports.** Long runs emit compact `commentary` checkpoints: what has been read, what has been ruled out, what proof obligation remains. The final answer is short unless the proof itself is the artifact.

---

## 9 Self-reflection rubric (before any inscription / merge)

- **Correctness.** Every step verified; no gap; no unsignalled assumption.
- **Rigor.** Every load-bearing claim carries *proved / conjectured / expected / heuristic / computed / folklore*.
- **Licensing.** Every theorem statement carries its α/β/γ/δ/ε licensing tags inline.
- **Attribution.** Every prior result cited by author + year + theorem / equation.
- **Concrete-before-abstract.** Worked case precedes general statement; compute $k=0$ or $E_2^{\mathrm{top}}$ before raising the topologisation ladder.
- **Voice.** Russian school + math-physics frontier; no bookkeeping vocabulary; no meta-narration; no hedging proved identifications.
- **Standalone.** No version labels, no phase labels, no prior-draft references.
- **Stage-discipline.** Every primitive object names its stage on $\mathsf{P}\to\mathsf{C}\to\mathsf{S}\to\mathsf{Z}\to\mathsf{A}$. No cross-stage collapse.
- **Deep-semantic merge.** Every cross-volume / cross-chapter cross-reference re-checked against the target.

If any category falls short — restart that category. Do not patch.

---

## 10 Proof-obligation discipline

- **Proved** → complete argument in this tree or cited reference (page + theorem + year). Inscribed with `\ClaimStatusProvedHere` + IV decorator on a disjoint-route witness in `compute/tests/test_*_iv.py`.
- **Conjecture / expected** → named evidence (worked case, cohomological computation, physical heuristic). `\ClaimStatusConjectured` + named hypotheses.
- **Heuristic** → physics argument named (BCOV, bootstrap, SUSY localization, anomaly matching) with rigor level called out. `\ClaimStatusHeuristic`.
- **Computed** → computation lives in the source tree; cite file + line.
- **Retracted** → previously claimed; now known false or out-of-scope. `\ClaimStatusRetracted` + repair pointer.

---

## 11 Long-context handling

Inputs > ~10K tokens (typical chapter, swarm log, frontier inventory): outline load-bearing sections internally, batch-read every citation in parallel, hold the whole chapter in context, compact side work that is not load-bearing math. The 400K context of GPT-5-Codex is enough for one chapter + dependencies + bibliography simultaneously; do not paginate.

---

## 12 Codex load order (per session)

1. `./CLAUDE.md` (full architectural canon).
2. `~/ecosystem/INVARIANTS.md §IV` (voice) + `~/ecosystem/AGENTS-HARNESS.md §VIII` (self-reflection rubric).
3. `notes/legacy/critique_2026_05_09_chiral_duality_master_consequence_map_v2.md` (architectural commitments).
4. `notes/legacy/vol2_platonic_architecture.md` (Vol II form).
5. Target chapter TeX + compute / Coq / Lean dependencies + cited bibliography.
6. Latest `adversarial_swarm_*/SYNTHESIS.md` if present.

---

## 13 Reference corpus

Beilinson--Drinfeld 2004 · Ayala--Francis · Costello 2011 · Costello--Gwilliam · Lurie HA · Gwilliam--Williams · Witten · Gaiotto · Kapranov · Polyakov 1987 · Borcherds 1995 · Gritsenko 1999 · Prochazka 1809.06993 · CKL 1704.08023 · Schiffmann--Vasserot · Gui--Li--Zeng 2212.11252.

---

## 14 Deep semantic merges only

Branches / worktrees differ → **always** deep semantic merge (see `./CLAUDE.md §22`). Read both sides in full; merge at the semantic level — stronger statement, tighter citation, more rigorous proof. Never `git reset --hard` / `git restore` / force-push to clobber work. Work loss is irrecoverable. **Never cut content.**

---

## 15 Escalation triggers

Proof cannot be discharged with honest rigor → **naming the open obligation is the deliverable**. Cross-volume disagreement → stop, report. Compute-vs-prose disagreement → stop, report; computation usually wins. Voice-table hypothesis missing at theorem statement → flag and add tag, do not silently rewrite.

---

## Code-writing discipline — repo application

Per `~/ecosystem/INVARIANTS.md §XIII`. Twelve rules instantiated for chiral-bar-cobar Vol II ($A_\infty$ chiral algebras; 3D HT QFT; Swiss-cheese / topologisation ladder; bicoloured primitive package; CP-1..4 Construction Problems):

1. **Think Before Coding.** Every $A_\infty$-edit names the chiral / homotopical structure and the claim-status. Every Construction Problem (CP-1..4) edit names the bicoloured-primitive package affected and the licensing-type tag.
2. **Simplicity First.** Single thesis (`CLAUDE.md §2`), four Construction Problems — the scope is fixed. No new licensing types beyond the five canonical ones. No additions outside the seventeen-forbidden-slogans constraint.
3. **Surgical Changes.** An edit in §6 (bicoloured primitive package) does not touch §7 (two-stage CY-chiral functor). §17 (hooks), §18 (skills), §19 (auto-memory) are internal — edits there do not modify mathematical chapters.
4. **Goal-Driven Execution.** Success = `pdflatex main.tex` clean; Construction Problem progress consistent with `CLAUDE.md §10`; four-part term-coining test passes; voice-scan clean; raeez-math-template symlink intact; Beilinson gate clean.
5. **Use the model only for judgment calls.** Cross-references, theorem-numbering, bibliography are deterministic. Codex drafts proofs and selects worked computations; it does not invent licensing types or violate the seventeen-forbidden-slogans list.
6. **Token budgets are not advisory.** Monograph; checkpoint between Construction Problems and between chapters. Long-form proof harness sessions are 30–60 minutes — load context first, build internal outline.
7. **Surface conflicts, don't average them.** Single-thesis statement wins over neighbouring chapter formulations. The bicoloured-primitive package is canonical for primitives used downstream. Compute disagreement → computation wins.
8. **Read before you write.** Read the affected Construction Problem and its hypothesis package. Read §6 before editing §7. Read the seventeen forbidden slogans before drafting prose. Grep `notes/claude_md_legacy_20260418.md` and `notes/agents_md_legacy_20260418.md` — never read whole.
9. **Tests verify intent.** Claim-status macros are the load-bearing test; four-part term-coining test is non-vacuous; the seventeen-forbidden-slogans list is a positive constraint. For load-bearing identities, seek independent derivations by multiple routes — disagreement is the deliverable.
10. **Checkpoint after every significant step.** Between Construction Problems, summarize hypothesis-package delta and progress. Between chapters, restate single-thesis impact. Subagents return evidence; main thread integrates via deep semantic merge.
11. **Match the codebase's conventions, even if you disagree.** `raeez-math-template.sty` per `INVARIANTS.md §XII`. Theorem environments per template. Subagent protocol governs parallel agent dispatch.
12. **Fail loud.** Announce every cross-ref break, dangling theorem, unhealed conjecture (`INVARIANTS.md §XI`). Honour the §16 "Do not" coda. Cross-volume / compute-vs-prose disagreements stop and report.

---

## 16 Do not

Block large user-authorized swarms. Propagate status-label wording when mathematics is waiting. Invent formulas from memory. `make` after every edit. AI attribution anywhere. `git stash` or amend. Read legacy files whole (grep by AP / V2-AP / FM index). Confuse this file with a configuration manual. Mathematician's manifesto. Shrink if it grows.
