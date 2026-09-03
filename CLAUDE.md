# CLAUDE.md (Vol II)

## Golden rule: no project management at any scale whatsoever in the manuscript

Reader-facing book and paper sources contain mathematics or physics only.
This is absolute at every scale: volume, part, chapter, section, theorem,
proof, remark, example, footnote, caption, table, front matter, back matter,
bibliography annotation, and PDF metadata.
Never insert agent/task instructions, audit or repair workflow, worktree or
commit details, review/packet status, test dashboards, TODO queues, progress
reports, ownership assignments, or references to models, agents, prompts, and
critiques. Keep such material outside the manuscript. Mathematical status is
different and must remain explicit: theorem, conditional consequence,
heuristic, conjecture, and open problem. State missing data as mathematics,
without narrating who found the gap. The objective is correct, rigorous,
complete mathematics; builds and prose polish do not replace proof. Pass this
rule to every subagent and check its returned prose.

> **Inheritance.** `~/ecosystem/INVARIANTS.md` — destructive-git forbidden list, multi-agent worktree concurrency, standalone-documents, Russian-school voice, every-file-into-the-repo, no-LLM-attribution, deep-semantic-merges, intelligence propagation.
> **Writing standard (mandatory).** `./MATHEMATICAL_PHYSICS_NUMBER_THEORY_GEOMETRY_ALGEBRA_HOMOTOPY_THEORY_WRITING_STANDARDS.md`. Witten · Etingof · Polyakov · Dirac · Feynman · Costello · Gaiotto. The prose IS mathematics; it does not describe mathematics. A sentence that does not state mathematics or physics is a defect.
> **Architecture.** `notes/legacy/critique_2026_05_09_chiral_duality_master_consequence_map_v2.md` (universal stage chain $\mathsf{P}\to\mathsf{C}\to\mathsf{S}\to\mathsf{Z}\to\mathsf{A}$ + five licensing types α/β/γ/δ/ε + four Construction Problems). `notes/legacy/vol2_platonic_architecture.md` (the seven-part form Vol II yearns to be).

A mathematician's working manifesto. The mathematics follows.

---

## 1 What this repository is for

An instrument for advancing human mathematical knowledge: how $A_\infty$ chiral algebras underwrite 3D holomorphic-topological QFT through the two-coloured Swiss-cheese-chiral-topological operad $\mathsf{SC}^{\mathrm{ch,top}}$ on the curve-line pair $(C, \R)$. Every read, edit, and inscription serves one true theorem at a time. When the choice is between mathematics and bookkeeping, do mathematics. The PostToolUse hook reconciles bookkeeping at session end.

---

## 2 The single thesis

$\mathsf{SC}^{\mathrm{ch,top}}$ on $(C, \R)$ is the universal home for the bulk--boundary structure of 3D HT QFT. Iterated Sugawara topologisation reaches the $E_\infty^{\mathrm{top}}$ Platonic endpoint via $w_{1+\infty}$. Specialised at $A = \mathrm{Vir}_c$ in the Brown--Henneaux dictionary $c = 3\ell/(2G_N)$, the Universal Holography functor $\Phi_{\mathrm{hol}}$ delivers the holographic boundary-CFT reading of pure 3D quantum gravity.

Seven parts hold the structure: I Open primitive · II $E_1$ core · III Seven faces of $r(z)$ · IV Characteristic datum · V HT landscape · **VI 3D quantum gravity (CLIMAX)** · VII Frontier. One argument, told four times.

---

## 3 The mathematics

**One structure.** $\mathsf{SC}^{\mathrm{ch,top}}$ — 2-coloured symmetric operad (coupled Voronov half-plane model: mixed components $\FM_{k,m}(\mathbb H)$; the decoupled product $\FM_k(\C)\times E_1(m)$ is only the depth-zero associated graded) with directional restriction $\mathsf{SC}^{\mathrm{ch,top}}(\ldots, \mathsf{top}, \ldots; \mathsf{cl}) = \varnothing$. Dunn additivity does NOT apply. Coupled homotopy-Koszulity is a named Conjecture (Livernet non-formality obstructs the formality route); Koszulity is proved at the gr level. Bar = holomorphic factorisation on $C$. Coproduct = topological factorisation on $\R$. Pairing over $C \times \R$ = bulk--boundary datum.

**One tower.** $T^{(n)} = [Q_{\mathrm{tot}}, G^{(n)}]$ promote chiral $\to E_{k+2}^{\mathrm{top}}$ on $Q$-cohomology. Affine KM at non-critical $\to E_3$. $\mathcal{W}_N \to E_{N+1}$. $\mathcal{W}_\infty[\lambda] \to E_\infty$ Platonic endpoint conditional on $\hypProchazka, \hypCKL, \hypPRSh, \hypYamada$.

**Five theorems** (shared with Vol I; Vol II carries the chiral lane). **A** bar--cobar inversion. **B** chiral Positselski. **C** derived-centre complementarity $\rho_K = \kappaChHodge(A) + \kappaChHodge(A^!) \in \{0, 13, 250/3\}$ on the computed families ($0$ on $\mathsf{G/L/C}$ and the $\mathsf{B}$-row Hodge sum; $13$, $250/3$ on $\mathsf{M}$; the Bershadsky--Polyakov slot is OPEN — its central-charge conductor $c_{\mathrm{BP}}(k)+c_{\mathrm{BP}}(-k-6) \equiv 50$ is a distinct invariant, and the Mukai conductor $K^{\kappa} = 2c_+(\mathrm{II}_{4,20}) = 8$ is a distinct sixth invariant, not a $\rho_K$ value). **D** obstruction-tower universality. **H** Hochschild concentration $\Zderch{A} \simeq \bulkChirHoch{A}$.

**Vol II contributions.** Curved-Dunn $H^2 = 0$ at $g \geq 2$; $\mathsf{SC}^{\mathrm{ch,top}}$ heptagon; chiral higher Deligne; universal celestial holography; periodic-CDG on admissible Kazhdan--Lusztig; unified $Q_g^{k, f, \mu}$ chiral quantum group.

**Five objects never conflated.** $A_b = \mathrm{End}_\cC(b)$, $\BarTwc{A_b}$, $A^i = H^\star \BarTwc{A_b}$, $A^!$, $\Zderch{A_b} \simeq \bulkChirHoch{A_b}$. $\Omega(\BarTwc{A_b}) = A_b$ is **inversion**, not Koszul duality. $A^!$ via **Verdier**. Bulk via **Hochschild**.

---

## 4 Beilinson's cut and the universal stage chain

> What limits forward progress is not the lack of genius but the inability to dismiss false ideas. — Sasha Beilinson

Every working object in the corpus (Vol I, Vol II, Vol III, Vol IV plus mixed-HT-strings + igusa satellites) sits on one stage of the universal chain
$$
\underbrace{\mathsf{P}}_{\text{primitive}}
\xrightarrow{\;\alpha\;}
\underbrace{\mathsf{C}}_{\text{chart}}
\xrightarrow{\;\beta\;}
\underbrace{\mathsf{S}}_{\text{shadow}}
\xrightarrow{\;\gamma\;}
\underbrace{\mathsf{Z}}_{\text{centre / bulk}}
\xrightarrow{\;\delta\;}
\underbrace{\mathsf{A}}_{\text{acting object}}.
$$

Cross-stage collapse — identifying objects at different stages without naming the licensing arrow — is the master pattern **`shadow = object`**. Every false slogan is one such collapse.

**Licensing rule.** A statement is not allowed to be primitive if it is only true after choosing a boundary object, passing to a trace, averaging from ordered to symmetric, taking a protected index, completing a category, imposing endpoint hypotheses, or installing descent data. Stated positively: **primitive objects first, shadows second, scalar modular forms last.**

Two structural lanes are parallel instances of the chain, unified at level $\mathsf{Z}$ by SC$^{\mathrm{ch,top}}$ heptagon face (6) Drinfeld-centre $=$ derived-centre identification:

| Stage | Open lane | CY lane (Vol III bridge) |
|-------|-----------|-------------------------|
| $\mathsf{P}$ | $\openFactCat$ on $\logCurve$ | CY$_d$-category |
| $\mathsf{C}$ | $A_b = \mathrm{End}_\cC(b)$ | $\PhiFA_d(A)$ stage-1 $E_d$-FA |
| $\mathsf{S}$ | $\BarTwc{A_b}$ twisting/coupling | $\SpCh = \int_{\Sigma_{d-1}}$ stage-2 specialisation |
| $\mathsf{Z}$ | $\Zderch{A_b} \simeq \bulkChirHoch{A_b}$ | derived chiral centre on $C$ |
| $\mathsf{A}$ | line operators on $A_b$ | $\Drinfdouble{\Yplus{X}}$ Drinfeld double |

A third **trace lane** carries $\mathrm{tr}_{\mathrm{ghost}}(Q_{\mathrm{BRST}}^2) = \mathrm{tr}_{\mathrm{Pentagon}} = \omega_{\mathrm{Borcherds}} = c_N(0)/2$ at $N \in \{1,2,3,4,6\}$. At $N = 1$: $\{5,5,5\}$ on the K3-Heisenberg witness.

---

## 5 The five licensing types

Every theorem statement carries the relevant licensing tags **in the statement** (not the introduction).

- **α — chart / scope / log decoration.** Choice of $b$, scope label for $\kappa$, tangential log $(D, \tau)$, BRST nilpotent, Stokes / branch datum.
- **β — comparison / functor / natural transformation.** Chiral Hochschild via Lurie HA.5.5, Drinfeld double $D$, protected Pfaffian $\protectedPfaff{\cdot}$, Hall--Borcherds residual, $\SpCh$ specialisation, MC injection.
- **γ — ambient declaration.** Chain-level vs $(\infty,1)$-categorical (BOTH equally load-bearing; never substitutable). $\mathrm{Ch}(\mathrm{Vect})$ vs weight-completed vs pro-object vs $J$-adic. Topological / analytic / Schwartz / formal completion. **Ambient is the prerequisite for every other licensing type.**
- **δ — endpoint / convergence.** $W_\infty \to E_\infty$ (Prochazka + CKL + PRSh + Yamada); PVA all-loop quantum (KZ SDR + Stokes + reflected weights + $T = [Q_{\mathrm{tot}}, G]$ lift); Universal Holography $\to$ dynamical metric (Brown--Henneaux + modular invariance + vacuum dominance + saddle dominance); Mittag--Leffler antighost-commutativity for iterated Sugawara.
- **ε — effectiveness / orientation / non-degeneracy.** Koszul-effectiveness for MC bijection (Vol I Theorem B / Positselski); $H^1(X, \C_X)$ vanishing for HT global; Pfaffian orientation; 6d hCS quartic obstruction $\int_X \mathrm{Tr}_{\mathrm{ad}} A(F_A)^3$ vanishing on verified loci.

The seven heptagon faces and the five licensing types are dual presentations: heptagon = operadic / categorical; types = epistemic / proof-discipline.

---

## 6 The bicoloured primitive package

The primitive open object is **not** an algebra. It is the bicoloured nine-tuple
$$
\bigl(\openFactCat\big|_{\logCurve},\;\; \cD^{\mathrm{cl}}\big|_X;\;\; b,\;\; A_b,\;\; F^{\mathrm{cl}};\;\; \HalfBraid,\;\; \mathrm{tr}^{\mathrm{cl}}_X;\;\; \TraceC\bigr).
$$
$\openFactCat$ — open factorisation dg-cat on $\logCurve = (X, D, \tau)$. $\cD^{\mathrm{cl}}$ — closed factorisation $\infty$-cat (chain-level: factorisation D-modules on $\mathrm{Ran}(X)$). $b$ — boundary vacuum (chart, not invariant). $A_b = \mathrm{End}_\cC(b)$. $F^{\mathrm{cl}}$ — closed-colour vacuum FA. $\HalfBraid$ — half-braiding. $\mathrm{tr}^{\mathrm{cl}}_X$ — closed-colour trace, evaluating to $c_N(0)/2$ at the Pentagon endpoint. $\TraceC$ — open-colour cyclic trace.

Seven heptagon faces: (1) bicoloured operadic primitive · (2) quantum complementarity · (3) bar--cobar inversion · (4) brace / chiral Deligne--Tamarkin · (5) topologisation ladder · (6) Drinfeld centre = derived centre · (7) PTVV derived-AG.

Modularity = open-colour cyclic trace + clutching, with closed-colour shadow $\mathrm{tr}^{\mathrm{cl}}_X$ carrying the modular consequences. Modularity is **never** a property of the closed algebra.

---

## 7 The two-stage CY-chiral functor + κ-tuple

$$
\Phi_d \;=\; \mathrm{Sp}^{\mathrm{ch}}_{\Sigma_{d-1}, C} \circ \Phi^{\mathrm{FA}}_d
\;:\;
\mathrm{CY}_d\text{-Cat} \to E_d\text{-HolFA}(X) \to \mathrm{ChirAlg}(C).
$$
Stage 1 ($\PhiFA_d$): canonical $E_d$-holomorphic FA on CY-$d$ via Kontsevich--Tamarkin + Costello--Gwilliam--Li BV. Stage 2 ($\SpCh$): $\int_{\Sigma_{d-1}}$ specialisation to reference curve $C$. Closed colour = stage 1; open colour = stage 2. The directional restriction expresses "stage 2 is specialisation of stage 1, never inversion."

Hall side parallel two-stage. $\Yplus{X} = H^\bullet_{\mathrm{eq}}(\cM^+_{\mathrm{eff}}(X), \phi_W)$ is the positive half. $\Drinfdouble{\Yplus{X}} = G(X)$ is the full quantum group, requiring Hall pairing + completion + integral form + stable-envelope transport + descent. $\CoHA{\C^3} = \Yplus{\fgl_1}$, NOT $\Wonepinf$ — the latter is the Fock representation of the Drinfeld double.

**κ-tuple.** $\kappa$ is a tuple, not a scalar:
$$
\kappaTuple{A} \;=\; \bigl(\kappaCat,\; \kappaChHodge,\; \kappaChHeis,\; \kappaBKM,\; \kappaFiber\bigr)(A).
$$
K3$\times E$ witness: $(0, 0, 3, 5, 24)$. The additive identity $\kappaBKM = \kappaChHodge + \chi(\cO_{\mathrm{fiber}})$ **fails** at every $N \in \{1,2,3,4,6\}$. Universal Borcherds weight $\kappaBKM(\Phi_N) = c_N(0)/2$ (Borcherds 1995, Gritsenko 1999). Vol II conductor $K^{\kappa_{\mathrm{ch}}} = 8 = \mathrm{ord}(H_1)$ on $\mathcal{B}$ is a *separate* invariant ($\hbar^2 \cdot K^{\kappa_{\mathrm{ch}}} = -1$). Bare $\kappa$ forbidden.

Cross-volume preamble macros: `\PhiFA, \SpCh, \HolFA, \EdHolFA, \intSigma, \hCS, \Yplus, \Drinfdouble, \CoHA, \Wonepinf, \kappaCat, \kappaChHodge, \kappaChHeis, \kappaBKM, \kappaFiber, \Kkappa, \kappaTuple`. APs V2-AP25--V2-AP29 (mirror V2-AP128--V2-AP132) protect.

---

## 8 The 17 forbidden slogans (voice table)

Every false slogan is replaced by its allowed form. Caught by `scripts/hooks/beilinson-gate.sh` PostToolUse and `make verify-licensing`.

| # | Forbidden | Allowed |
|---|-----------|---------|
| 1 | "Let $A$ be a chiral algebra" | $\primPkgBicolour$ + $A_b = \mathrm{End}_{\openFactCat}(b)$ |
| 2 | "$\BarTwc{A}$ is the bulk" | $\bulkOf{A} = \bulkChirHoch{A}$; $\BarTwc{A}$ is twisting |
| 3 | "$E_1$-bar explains $2d \to 3d$" | chiral Deligne--Tamarkin / Swiss-cheese |
| 4 | "open sector on $X$" | "open sector on $\logCurve$" |
| 5 | "the closed algebra is modular" | "open category carries cyclic trace + clutching" |
| 6 | bare "$\kappa$" | $\kappaTuple{A}$ |
| 7 | "$\Phi_d : \mathrm{CY}_d \to \mathrm{ChirAlg}$" | $\Phi_d^{(\Sigma_{d-1}, C)} = \SpCh \circ \PhiFA_d$ |
| 8 | "$\CoHA{\C^3} = \Wonepinf$" | $\CoHA{\C^3} = \Yplus{\fgl_1}$; $\Wonepinf$ via Fock + Drinfeld double |
| 9 | "6d hCS = 3d CS" | "6d hCS realises $\PhiFA_3$; quartic obstruction; 3d-CS = analogy" |
| 10 | "formal local HT $\Rightarrow$ compact global" | "+ $\effHTGlobalDR$ + descent + QME + anomaly + locality" |
| 11 | "$\Deltafive$ = Hilbert space" | "Borcherds shadow; $\protectedPfaff{\operatorPrim{X}} = \Deltafive$ is construction problem" |
| 12 | "$\ZBPS$ = path integral" | "scalar shadow; gravity-line via Hall--Borcherds residual" |
| 13 | "Universal Holography = QG" | "boundary-CFT reading of 3D gravity" |
| 14 | "$W_\infty \Rightarrow E_\infty$" | conditional on $\hypProchazka, \hypCKL, \hypPRSh, \hypYamada$ |
| 15 | "class M chain-level" | conditional on $\hypAmbientWtCpl$; fails in $\mathrm{Ch}(\mathrm{Vect})$ |
| 16 | "PVA $\Rightarrow$ quantum theory" | classical only; quantum conditional on $\hypKZSDR, \hypStokes, \hypReflWts, \hypTLift$ |
| 17 | "quadratic dual = Koszul" | injection in general; bijection conditional on $\effKoszul$ |

---

## 9 The four Construction Problems (Frontier)

The cut creates space for four operator-level construction problems:

1. **$\mathfrak{D}_X$ for K3$\times E$** with $\protectedPfaff{\mathfrak{D}_X} = \Deltafive$. Source: igusa `~/igusa-cusp-form/main.tex:94--118`.
2. **Gravity-line operator algebra** with Pentagon-face scalar trace $= \Phitenun = \Deltafive^2$. Source: `chapters/connections/3d_gravity.tex:8429`.
3. **Unified PVA-quantum HT theory** with classical $\lambda$-Jacobi limit and $E_3$-lift on $Q$-cohomology.
4. **Chiral Positselski** extending Vol I Theorem B at chiral generality; specialisation to quadratic recovers Gui--Li--Zeng MC bijection.

---

## 10 Progress

**Counts.** A new theorem proved with verifiable proof body. A new example (κ-tuple coordinates / heptagon face / topologisation height) for an algebra not yet tabulated. A falsified claim repaired. A healed hypothesis. A first-principles computation replacing a citation black box. An IV witness via disjoint-route proof for an existing $\ClaimStatusProvedHere$. An operator-level construction advancing one of the four Construction Problems.

**Does not count.** Status-table rows. Label renames. Counting equivalences. Scope propagation across ten files. FRONTIER advertising retractions. CLAUDE.md ↔ AGENTS.md harmonisation. Bookkeeping. The hook catches these.

**Epistemic hierarchy** (higher wins): direct computation > $.tex$ source ±100 lines > tests > IV decorators > primary literature > `concordance.tex` > this file > memory.

---

## 11 Self-coherence

The manuscript stands for itself. Every $.tex$ chapter is standalone, up-to-date, consistent. No references to previous versions, intermediate ansätze, retracted values, drafting-history commentary. If a formula used to be $X$ and now it is $Y$, the manuscript says $Y$. When a retraction is informative — a proof attempted and failed whose failure illuminates why the successful proof is forced — state the failed argument and its flaw *as mathematics*, never as drafting record.

CG voice: Russian elite (Gelfand, Drinfeld, Beilinson, Kapranov, Etingof, Kontsevich, Bezrukavnikov) + math-physics elite (Polyakov, Witten, Costello, Gaiotto, Nekrasov). Show don't tell. Forbidden in chapters: bookkeeping vocabulary (Wave / round / AP$n$ / Pattern); meta-narration ("we now turn to," "notably," "in the present work"); hedging proved identifications. Required: titles name objects; definitions preceded by the question; symbols defined at first use; physical claims labelled (theorem / heuristic / metaphor); economy. Full canon: `./MATHEMATICAL_PHYSICS_NUMBER_THEORY_GEOMETRY_ALGEBRA_HOMOTOPY_THEORY_WRITING_STANDARDS.md`.

---

## 12 Cross-volume coherence

- **Vol I** `~/chiral-bar-cobar` — five-theorem core; $E_1$--$E_1$ operadic Koszul duality; canonical κ / r(z) / S_r in `chapters/examples/landscape_census.tex`.
- **Vol II** `~/chiral-bar-cobar-vol2` — this volume.
- **Vol III** `~/calabi-yau-quantum-groups` — CY-to-chiral functor, K3$\times E$ Hall--Borcherds, BKM, κ-stratification, 6d hCS quartic.
- **Vol IV** `~/chiral-bar-cobar-vol4` — verification capstone; pairs every Vols I--III ProvedHere with disjoint-route IV.
- **mixed-HT-strings** `~/mixed-holomorphic-topological-strings` — local model + global de Rham obstruction.
- **igusa-cusp-form** `~/igusa-cusp-form` — $\Deltafive$ + Construction Problem 1.

Load-bearing claims about $\mathsf{SC}^{\mathrm{ch,top}}$, the topologisation ladder, or the bulk--boundary pair must be consistent with Vol I theorems and Vol III κ-stratification. **Disagreement is the deliverable.**

---

## 13 Build, test, audit

```
make fast / make / make release / make check / make test
make verify-independence  # @independent_verification audit for ProvedHere
make verify-licensing     # whole-volume voice + licensing grep
make clean-builds / make count
```

Single test: `compute/.venv/bin/python -m pytest compute/tests/test_<name>.py -q -ra`. Warm rebuilds: `export MKD_BUILD_NS="agent-$$"` once per session. Builds in `/tmp/mkd-chiral-bar-cobar-vol2-<NS>/` (parallel-safe). Session-end only, user opt-in. Requires TeX Live 2024+ (memoir, EB Garamond, newtxmath); Python 3.10+; ~60 engines; ~5400 tests.

---

## 14 Repository layout

`main.tex` entry; preamble holds all macros; chapters use `\providecommand`. `chapters/frame/` preface + part intros. `chapters/theory/` Parts I--IV core. `chapters/examples/` landscape examples. `chapters/connections/` Parts II--III + V--VII (line operators, celestial holography, HT bulk-boundary, 3D gravity climax, Universal Holography, $w_{1+\infty}$, Vol I/III bridges). `appendices/` brace signs, orientations, FM proofs. `compute/lib/` Python engines + `independent_verification.py`. `compute/tests/` pytest; `test_*_iv.py` are IV witnesses. `scripts/hooks/beilinson-gate.sh` version-controlled hook (install: `cp scripts/hooks/beilinson-gate.sh .claude/hooks/`). `notes/` workshop floor, never reader-facing. `FRONTIER.md` live queue (`ROADMAP_85_TO_100.md` concluded, archived at `archive/2026-04/`). `standalone/` extracted papers. `out/` build. Memoir + EB Garamond.

---

## 15 Git, authorship

All commits by **Raeez Lorgat** only. **No AI attribution anywhere**: no `Claude` / `Anthropic` / `Co-Authored-By` / `Generated with` / robot-emoji — in commits, comments, docstrings, manuscripts. Pre-commit hook nudges. Two remotes (`origin`, `ainfinity`); both receive main on push. `git stash` forbidden — use `git diff > patch.diff && git apply`. Do not amend without explicit instruction.

---

## 16 Bookkeeping

- `notes/legacy/critique_2026_05_09_chiral_duality_master_consequence_map_v2.md` — architectural reconstitution map.
- `notes/legacy/vol2_platonic_architecture.md` — seven-part platonic form.
- `notes/antipatterns_catalogue.md` — live AP registry (V2-AP25--132 + Vol I import + HT-strings AP2181--AP2189).
- `notes/first_principles_cache_comprehensive.md` — confusion-pattern registry; Vol I import block.
- `notes/claude_md_legacy_20260418.md` / `notes/agents_md_legacy_20260418.md` — prior manifestos (lossless). Grep by index; do not read whole.
- `~/chiral-bar-cobar/CLAUDE.md` · `~/chiral-bar-cobar/chapters/examples/landscape_census.tex` · `~/calabi-yau-quantum-groups/CLAUDE.md`.

Claim-status tags (`\ClaimStatusProvedHere`, `\ClaimStatusConjectured`, `\ClaimStatusEvidence`, `\ClaimStatusRetracted`) are reader-facing bookkeeping. When uncertain, name the proof obligation; heal the proof, statement, or construction; do not downgrade.

---

## 17 Hooks

`PreToolUse(Agent)` cache injection · `PreToolUse(Bash, git commit)` no-AI-attribution reminder · `PostToolUse(Edit|Write)` `beilinson-gate.sh` AP + cache + voice-table sweep · `Stop` session-end summary. The Beilinson gate silently enforces the cut: every edit triggers the licensing audit. Address what it flags; return to mathematics.

---

## 18 Skills

Invoke via Skill tool when the task matches: `chriss-ginzburg-rectify` · `audit` · `propagate` · `verify` · `rectify` · `investigate` · `compute-engine` · `build` · `research-swarm` (30+ elite agents) · `attack-heal-swarm-loop` · `codex:rescue`.

Skills are not automatically aware of the manuscript voice; the Skill prompt encodes it. Invoke with explicit licensing-type tags.

---

## 19 Auto-memory

Persistent file-based memory at `~/.claude/projects/-Users-raeez-chiral-bar-cobar-vol2/memory/`; index in `MEMORY.md` (truncated past 200 lines; concise `- [Title](file.md) — hook`).

Four types: **user** (role / preferences), **feedback** (corrections AND validated approaches; rule + `**Why:**` + `**How to apply:**`), **project** (work / decisions; dates absolute), **reference** (external systems). Save on explicit user request; user correction; user validation of a non-obvious choice; reference to prior conversations. Do not save code patterns, file paths, git history, debugging fixes, ephemeral state. Before recommending from memory: verify against current code (point-in-time). If user says ignore memory, do not apply, cite, or compare.

---

## 20 Long-form proof harness

Frontier mathematics runs in maximum-effort mode. Deepest available model + highest reasoning budget. For Claude Opus 4.7 (1M context), deliberation budget unbounded; 30--60 minute agent runs are normal for theorem repair, cross-volume synthesis, adversarial review, primary-source reconstruction.

Load context: this file, `notes/critique_*_v2.md`, `notes/legacy/vol2_platonic_architecture.md`, target chapter, dependencies, bibliography, compute, cross-volume anchors. Build the internal outline before the first edit.

For any load-bearing identity, seek independent derivations by **multiple routes**: worked example, formal argument, primary literature, local computation, cross-volume consistency. **Disagreement is the deliverable.**

After every proposed repair, run an **attack-heal loop**: strongest counterexample, sign / convention / ambient-category check, missing hypothesis, false functoriality, unproved equivalence, numerical constant. Heal and attack again until the theorem closes or the exact obstruction is named. Do not downgrade the manuscript to close. Subagents provide evidence, not authority; the main thread integrates by deep semantic merge.

---

## 21 Subagent / swarm protocol

Large user-authorized swarms permitted. When the user asks for adversarial / rescue / review / cross-volume sweeps, treat as authorization for the largest useful swarm. Do not downshift to historical 3 / 5 / 30-agent caution.

Swarm design explicit before launch: partition by disjoint mathematical axes / files / proof obligations; name the integration owner; forbid agents from reverting work they did not make; require deep semantic merge across all five repos.

Every attack-heal agent returns a compact, checkable report: claim attacked · failure mode or proof · file anchors · primary anchors · exact formulas/constants · claim-status recommendation · files changed · tests/computations run · remaining obligations. Subagents do not vote truth into existence.

---

## 22 Branch / worktree reconciliation — DEEP SEMANTIC MERGES ONLY

Branches / worktrees differ → **always** deep semantic merge. Never discard one side without reading. Never `git reset --hard` / `git checkout --` / `git restore` to clobber work as a shortcut. Never force-push to obliterate upstream divergence. Read both sides in full; merge at the semantic level — stronger statement, tighter citation, more rigorous proof. Applies to `git pull`, `git merge`, worktree reconciliation, cherry-picks, rebase conflicts, push-rejection upstream divergence. Work loss is irrecoverable. **Never cut content.**

---

## Code-writing discipline — repo application

Per `~/ecosystem/INVARIANTS.md §XIII`. Twelve rules instantiated for chiral-bar-cobar Vol II ($A_\infty$ chiral algebras; 3D HT QFT; Swiss-cheese / topologisation ladder; bicoloured primitive package; CP-1..4 Construction Problems):

1. **Think Before Coding.** Every $A_\infty$-edit names the chiral / homotopical structure and the claim-status. Every Construction Problem (CP-1..4) edit names the bicoloured-primitive package affected and the licensing-type tag.
2. **Simplicity First.** Single thesis (§2), four Construction Problems — the scope is fixed. No new licensing types beyond the five canonical ones. No additions outside the seventeen-forbidden-slogans constraint (§8).
3. **Surgical Changes.** An edit in §6 (bicoloured primitive package) does not touch §7 (two-stage CY-chiral functor). §17 (hooks), §18 (skills), §19 (auto-memory) are internal — edits there do not modify mathematical chapters.
4. **Goal-Driven Execution.** Success = `pdflatex main.tex` clean; Construction Problem progress consistent with §10; four-part term-coining test on any new term passes; voice-scan clean per §8; raeez-math-template symlink intact; Beilinson gate clean (`scripts/hooks/beilinson-gate.sh`).
5. **Use the model only for judgment calls.** Cross-references, theorem-numbering, bibliography are deterministic. The model drafts proofs and selects worked computations; it does not invent licensing types or violate the seventeen-forbidden-slogans list.
6. **Token budgets are not advisory.** Monograph; checkpoint between Construction Problems and between chapters. Long-form proof harness sessions are 30–60 minutes — load context first, build internal outline.
7. **Surface conflicts, don't average them.** Single-thesis statement (§2) wins over neighbouring chapter formulations. The bicoloured-primitive package (§6) is canonical for primitives used in §7 onward. Compute disagreement → computation wins (§22).
8. **Read before you write.** Read the affected Construction Problem and its hypothesis package. Read §6 before editing §7. Read §8 (seventeen forbidden slogans) before drafting prose. Grep `notes/claude_md_legacy_20260418.md` — never read whole.
9. **Tests verify intent.** Claim-status macros are the load-bearing test; four-part term-coining test is non-vacuous; the seventeen-forbidden-slogans list is a positive constraint. For load-bearing identities, seek independent derivations by multiple routes — disagreement is the deliverable.
10. **Checkpoint after every significant step.** Between Construction Problems, summarize hypothesis-package delta and progress per §10. Between chapters, restate single-thesis impact. Subagents return evidence; main thread integrates via deep semantic merge (§22).
11. **Match the codebase's conventions, even if you disagree.** `raeez-math-template.sty` per `INVARIANTS.md §XII`. Theorem environments per template. §21 subagent protocol governs parallel agent dispatch.
12. **Fail loud.** Announce every cross-ref break, dangling theorem, unhealed conjecture (`INVARIANTS.md §XI`). Honour §23 "Do not" coda. Cross-volume / compute-vs-prose disagreements stop and report.

---

## 23 Do not

Block large user-authorized swarms. Propagate status-label wording when mathematics is waiting. Invent formulas from memory. `make` after every edit. AI attribution anywhere. `git stash` or amend. Read `notes/claude_md_legacy_20260418.md` whole (grep by index). Confuse this file with a configuration manual. Mathematician's manifesto. Shrink if it grows.
