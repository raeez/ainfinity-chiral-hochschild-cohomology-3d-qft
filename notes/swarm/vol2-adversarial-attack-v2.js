export const meta = {
  name: 'vol2-adversarial-attack-v2',
  description: 'Relaunch: first-principles adversarial attack on Vol II across 17 axes, with salvaged leads injected and death-proof disk journaling',
  phases: [
    { title: 'Attack', detail: '17 axis agents with injected leads; findings journaled to disk after each confirmation' },
    { title: 'Verify', detail: 'independent skeptic referee per non-minor finding' },
  ],
}

const FINDINGS_SCHEMA = {
  type: 'object',
  properties: {
    axis: { type: 'string' },
    summary: { type: 'string', description: 'Three-sentence verdict on the axis: what is solid, what is broken, what is missing.' },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          title: { type: 'string' },
          severity: { type: 'string', enum: ['critical', 'major', 'minor'] },
          verdict: { type: 'string', enum: ['false', 'unproved', 'overstated', 'unlicensed', 'inconsistent', 'missing-artifact', 'fragile'] },
          claim: { type: 'string', description: 'The exact claim attacked, quoted or precisely paraphrased' },
          anchors: { type: 'array', items: { type: 'string' }, description: 'file:line anchors' },
          evidence: { type: 'string', description: 'The mechanism of failure: exact computation, counterexample, contradiction, or missing step' },
          heal: { type: 'string', description: 'The precise repair: corrected statement, missing hypothesis to add, proof step to supply, or named obstruction' },
        },
        required: ['title', 'severity', 'verdict', 'claim', 'anchors', 'evidence', 'heal'],
      },
    },
    sound: { type: 'array', items: { type: 'string' }, description: 'Load-bearing claims attacked that SURVIVED, with one-line reason each' },
  },
  required: ['axis', 'summary', 'findings', 'sound'],
}

const VERDICT_SCHEMA = {
  type: 'object',
  properties: {
    isReal: { type: 'boolean' },
    confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
    assessment: { type: 'string' },
    corrected_heal: { type: 'string' },
  },
  required: ['isReal', 'confidence', 'assessment'],
}

const PREAMBLE = `You are an adversarial mathematical auditor on Vol II of a research monograph (A-infinity chiral algebras and 3D holomorphic-topological QFT), repo /Users/raeez/chiral-bar-cobar-vol2. Vol I lives at /Users/raeez/chiral-bar-cobar, Vol III at /Users/raeez/calabi-yau-quantum-groups. Single thesis: the two-coloured Swiss-cheese-chiral-topological operad SC^{ch,top} on the curve-line pair (C,R) is the universal home for bulk-boundary structure of 3D HT QFT; iterated Sugawara topologisation T^{(n)} = [Q_tot, G^{(n)}] climbs an E_{k+2}^{top} ladder toward an E_infinity^{top} endpoint via w_{1+infinity}; at A = Vir_c under Brown-Henneaux c = 3l/(2G_N) the Universal Holography functor gives the boundary-CFT reading of pure 3D gravity.

The manuscript enforces a licensing discipline: every object sits on the stage chain P(primitive) -> C(chart) -> S(shadow) -> Z(centre/bulk) -> A(acting object); identifying objects across stages requires a named licensing arrow of type alpha (chart/scope/log choice), beta (comparison functor), gamma (ambient declaration: Ch(Vect) vs weight-completed vs pro vs J-adic; chain-level vs (infinity,1)-categorical), delta (endpoint/convergence hypotheses), epsilon (effectiveness/orientation/non-degeneracy). Claim-status macros: ClaimStatusProvedHere requires a real proof body; ClaimStatusConditional requires named hypotheses; also Evidence and Retracted.

YOUR JOB: attack, from first principles, the mathematics in your assigned axis. Hunt: (1) false identities - verify EVERY numerical or algebraic identity you meet by direct computation (python at compute/.venv/bin/python, sympy installed); (2) unproved claims tagged ProvedHere - a proof body that defers its content to an assumption, a forward reference, or a circular citation is NOT a proof; (3) sign / convention / ambient-category errors; (4) missing hypotheses - statements true only under conditions absent from the statement; (5) false functoriality and unproved equivalences; (6) citation fraud - results attributed to primary literature that the literature does not contain (check against your own knowledge of the actual papers; one fabricated bibitem is already confirmed in this volume); (7) internal contradictions between chapters or against Vol I / Vol III; (8) shadow=object collapses ('Bar(A) is the bulk', bare kappa, 'Delta_5 = Hilbert space', 'W_infinity => E_infinity unconditional', 'PVA => quantum', 'quadratic dual = Koszul', 'Universal Holography = quantum gravity').

DEATH-PROOF JOURNALING (mandatory). Your run may be killed without warning at any moment (spend limits). Immediately after EACH finding you confirm - and after finishing each file - append a compact block to /tmp/attack_findings_AXISKEY.md via Bash (cat >> with a heredoc): title / severity / verdict / anchors with file:line / evidence with exact numbers / heal. Write to disk BEFORE moving to the next target. At the end, still return the full structured output as instructed. If your time is running long, prioritize journaling confirmed findings over breadth.

Protocol: FIRST run git diff -- <your files> ; the tree carries uncommitted WIP and half-finished edits are prime targets - attack both the committed base and the WIP layer. Read actual theorem statements and proof bodies, not abstracts. Strongest-counterexample mindset: for each load-bearing claim try to break it; if it survives serious attack, record it in 'sound' with a one-line reason - survivals matter as much as kills. Severity: critical = false or unproved load-bearing claim tagged ProvedHere, or an internal contradiction; major = overstated/unlicensed claim, missing hypothesis, wrong constant, lost content; minor = wording/anchor drift. Heals must heal: prefer repairing the proof or naming the exact missing hypothesis over downgrading the claim; recommend downgrade only for genuine falsity.

You MUST NOT edit any repo file (the /tmp journal is the only file you write). Evidence only. Every finding needs exact file:line anchors and, where numerical, the exact computation. A finding without an anchor and a mechanism is worthless.

YOUR AXIS: `

const AXES = [
  {
    key: 'ladder',
    brief: `THE TOPOLOGISATION LADDER. Files: chapters/connections/e_infinity_topologization.tex (1788 lines), chapters/connections/w_infty_e_infty_endpoint_platonic.tex, chapters/connections/monster_chain_level_e3_top_platonic.tex, chapters/connections/fractional_ghost_chain_level_platonic.tex.
KNOWN FROM THE MAIN THREAD (verify, anchor every instance, extend - do not re-derive): (a) e_infinity_topologization.tex:407 asserts wt(antighost c-bar_a) = 0 while Step 1 of thm:iterated-sugawara-construction needs J_a = [Q_tot, c-bar_a] with wt(J_a) = 1 - consistent only if Q_tot carries an undeclared weight-+1 component; find whether ANY file declares Q_tot's conformal bigrading, and write the corrected bookkeeping. (b) The theorems state identities 'in H^bullet(...)' (lines 331, 440, 494) where the Koszul-Tate contraction + filtered perturbation lemma machinery actually produces transferred CHAIN-LEVEL identities - and the ladder theorem (line 728) consumes chain-level homotopy data; list every 'on cohomology' site in the chapter that should be a transferred-complex chain-level statement. (c) thm:casimir-antighost-commutativity (line 422): the mechanism (KT contraction, filtered perturbation transfer, Kirillov-Kostant centrality of invariants - that part is sound - and acyclicity coherence tower) is plausible; attack its two weak joints: the filtration-convergence step (is the filtration exhaustive+complete on the completed local BV vertex algebra? where does the perturbation lemma's convergence get checked?) and the assertion at ~line 489 that the screening charges commute with Q_tot 'because they are boundary BRST charges' (assertion or proof?).
Then the original list: (d) thm:e-infinity-weightwise-inverse-limit - Mittag-Leffler usage; (e) the four endpoint hypotheses (Prochazka, Creutzig-Kanade-Linshaw, Pope-Romans-Shen, Yamada) - real results saying what is claimed? (f) which algebras does the ladder cover: affine KM has only the n=2 rung; Virasoro (class M, non-Sugawara) - is it covered by any rung, and is every flat downstream slogan ('Affine KM -> E_3, W_N -> E_{N+1}') consistent with the actual hypotheses? (g) prop:stress-tower-dimension-count - verify numerically.`,
  },
  {
    key: 'heptagon',
    brief: `THE SC^{ch,top} OPERAD AND ITS HEPTAGON. Files: chapters/theory/sc_chtop_heptagon.tex (3746 lines), chapters/theory/axioms.tex (1877 lines).
LEADS FROM A PRIOR PASS (it read lines 1-1300 + the WIP diffs of both files and said 'strong targets already found' - the WIP diff regions are prime; start there): run git diff -- chapters/theory/sc_chtop_heptagon.tex chapters/theory/axioms.tex and attack every changed region first.
KNOWN FROM THE MAIN THREAD (do not re-derive; extend): face (6) thm:drinfeld-centre-sc-face (line 1504) is ClaimStatusConditional - the factorization-module centre comparison is assumed in the statement - while FRONTIER.md 'load-bearing closures' and notes/legacy/critique_2026_05_09_..._v2.md section 2 call it 'already proved'. Your jobs: (i) determine whether the comparison hypothesis is dischargeable on any stated locus (Koszul locus? evaluation-generated core?) from results elsewhere in the volume; (ii) find every site in ANY chapter that uses face (6) or cor:open-cy-lane-unification as unconditional input (grep for the labels + 'face (6)' + 'face~(6)'); (iii) the (H1)-(H4) orphan question: axioms.tex, chapters/connections/concordance.tex, chapters/theory/pva-preview.tex, chapters/examples/examples-complete.tex contain '(H1)' hits - classify each as eliminated-standing-hypotheses orphan vs locally-defined package (the climax chapter's (H1)-(H4) is a local analytic package, defined near programme_climax_platonic.tex:1918 - name collisions between distinct (H1)-(H4) packages are themselves a finding).
Then: (a) the dioperad definition well-formedness (directional restriction; composition closure; units); (b) classify all seven faces: proved-here-with-real-proof / imported / asserted / conditional; (c) the Dunn-additivity boundary (disclaimed for the two-coloured structure, used after topologisation) - stated precisely where? (d) Swiss-cheese literature kinship claims (Voronov, Hoefel, Idrissi) accurate?`,
  },
  {
    key: 'hochschild',
    brief: `THE DERIVED CENTRE / CHIRAL HOCHSCHILD AXIS. File: chapters/connections/hochschild.tex (5863 lines; 1168-line uncommitted diff) plus chapters/theory/chiral_higher_deligne.tex (the original thm:chiral-higher-deligne lives there at line ~465).
LEAD FROM A PRIOR PASS: a 'strong target' sits in the Zhu bridge region, hochschild.tex ~1030-1290 (the CDG-compatibility theorem at line 1030 tagged ProvedHere, and the line-1190 remark where chiral-higher-deligne 'forces' concentration). Start there: read 1000-1330, identify the defect, verify it (run computations if numerical), anchor it.
KNOWN FROM MAIN THREAD: thm:chiral-higher-deligne clause (2) (universality) is explicitly conditional on the two-coloured (SC^{ch,top})^!-cobar contracting homotopy (rem:chd-universality-hypotheses names it 'a genuinely different theorem, not yet established'), while FRONTIER.md glosses it as an unconditional universality closure - find every downstream site treating universality as proved.
Then: (a) Theorem H concentration ChirHoch^bullet subset {0,1,2} on the Koszul locus at generic level - reconstruct the actual proof ('E_3-rigidity + PBW collapse' per FRONTIER, but the climax chapter's deleted-claims table says 'E_3-PBW proves concentration' was DELETED as a claim - reconcile); check grading conventions; check 'generic level' is defined and critical level excluded. (b) thm:chd-ds-hochschild (DS compatibility): restricted to principal nilpotents in its statement? (the (3,2) sl_5 obstruction elsewhere shows non-principal fails). (c) rem:cyclic-hochschild-four-objects licensing. (d) thm:hochschild-bridge-higher-genus proof body. (e) Bulk = chiral Hochschild cochains 'physical prefactorization model' theorem at ~line 792 - what is its proof?`,
  },
  {
    key: 'kappa',
    brief: `THEOREM C COMPLEMENTARITY AND THE KAPPA-TUPLE ARITHMETIC. Files: chapters/examples/rosetta_stone.tex, chapters/examples/examples-complete*.tex, chapters/connections/programme_climax_platonic.tex (thm:kappa-tuple-primitivity-orthogonal-shimura at ~546, thm:kappa-tuple-primitivity-class-B at ~675), chapters/connections/concordance.tex; canonical cross-volume source /Users/raeez/chiral-bar-cobar/chapters/examples/landscape_census.tex.
Attack: (a) rho_K = kappa_ch^Hodge(A) + kappa_ch^Hodge(A^!) in {0, 8, 13, 250/3, 98/3} on the five archetypes G/L/C/M/B - for each value find the computation that produces it; run the engines (compute/lib/genus1_kappa_verification.py exists; ls compute/lib | grep -iE 'kappa|rho|complementar' and run them). (b) Virasoro: kappa(Vir_c) + kappa(Vir_{26-c}) = 13 - determine kappa_ch^Hodge(Vir_c) as a function of c, verify c-independence of the sum, and pin which archetype row Virasoro belongs to (the removed main.tex WIP text said 13; the M-row value is 250/3 - resolve which class Vir is in and whether 13 is the C-row or the Vir sum). (c) K3xE witness (0,0,3,5,24): the IV test compute/tests/test_kappa_tuple_iv.py passes (70 assertions, five disjoint routes) - audit the TEST itself for rubber-stamping (do the five routes genuinely compute, or assert literals?). (d) The conductor K^kappa = 8 = ord(H_1) on class B with hbar^2 K^kappa = -1: find the definition sites, check arithmetic consistency (what is hbar here; hbar^2 = -1/8 of what). (e) Cross-check 10 kappa values used in Vol II chapters against Vol I landscape_census.tex - disagreements are deliverables. (f) beta_4 = 13 numerology: 12(1/2+1/3+1/4) = 13 and rho_K contains 13 - does any chapter identify these? If so, attack the identification; if not, note the coincidence is unflagged.`,
  },
  {
    key: 'curved-dunn',
    brief: `CURVED-DUNN, MODULAR COMPOSITION, CLASS M AMBIENT. Files: locate via grep -rl 'curved-dunn\\|curved Dunn' chapters/ standalone/; standalone/curved_dunn_two_complex_bridge.tex; chapters/connections/relative_feynman_transform.tex (287-line WIP diff); the chapters holding thm:curved-dunn-H2-vanishing-all-genera, thm:irregular-singular-kzb-regularity, prop:bv-bar-class-m-weight-completed (a prior pass was reading 'class M files and the bv_brst proposition' when killed - bv_brst.tex has a 558-line WIP diff, attack it too).
Attack: (a) thm:curved-dunn-H2-vanishing-all-genera: what cochain complex is H^2 computed in, what does vanishing at g >= 2 buy, is the proof a computation or an assertion? The g=1 boundary: thm:irregular-singular-kzb-regularity replaces KZ Stokes via Jimbo-Miwa - does Jimbo-Miwa 1981 isomonodromy actually cover the claimed pole orders? (b) prop:bv-bar-class-m-weight-completed and 'genuinely false in Ch(Vect)': does the cited computation prove FALSITY (nonzero obstruction class) or non-convergence of one construction? (c) Fay trisecant usage in the F1 heal path - actual computation anywhere? (d) relative_feynman_transform.tex WIP diff - attack what changed. (e) Mittag-Leffler antighost-commutativity at higher genus - stated/used/proved where?`,
  },
  {
    key: 'gravity',
    brief: `THE 3D GRAVITY CLIMAX. File: chapters/connections/3d_gravity.tex (13256 lines; 1760-line WIP diff - read the diff FIRST; a prior pass was reading the period-2-parity recursion machinery for the m_4/m_5 derivation conventions when killed - include prop:gravity-m4 (line ~699) and prop:gravity-m5 (~772) and their recursion conventions). Also chapters/connections/universal_holography_functor.tex, chapters/connections/programme_climax_platonic.tex (thm:universal-holography-master at ~119 - KNOWN from main thread: it is properly ClaimStatusConditional with scope, licensing in prose not macros; and its internal Open Frontier numbering F1=W(p)-tempering, F3=class-M-original-complex DIFFERS from FRONTIER.md's F1/F3 - do not re-derive, but check the chapter's other cross-references).
Attack: (a) the four delta-hypotheses (hypBHdict, modular invariance, vacuum dominance, saddle dominance) at every algebra-to-physics crossing in 3d_gravity.tex - find crossings without them (the ~62 untagged theorems flagged by the licensing gate are candidates). (b) prop:gravity-page-curve: WIP recasts as raw-transseries model under sectorial Borel + Stokes + hypotheses - check the body matches and no downstream site uses the old unconditional reading. (c) scalar genus series radius 4pi^2 in hbar^2 - verify against A_scalar = (2pi)^2. (d) prop:pole-order-classification (~219), cor:gauge-gravity-dichotomy (~412), prop:formality-depth-discriminant (~512) - attack proof bodies. (e) alpha_T^mix pole-order material (recent commits acb7090) - internal consistency. (f) conj:gravity-line-hall-borcherds-comparison ~8429 - conjecture status respected everywhere? (g) Brown-Henneaux factor c = 3l/(2G_N) written consistently?`,
  },
  {
    key: 'thqg-a',
    brief: `TOPOLOGICAL-HOLOGRAPHY QG EXTENSIONS, ANALYTIC GROUP. Files: chapters/connections/thqg_perturbative_finiteness.tex (2732-line WIP diff), thqg_nonperturbative.tex, thqg_fredholm_partition_functions.tex, thqg_modular_bootstrap.tex.
LEAD FROM A PRIOR PASS (verify first): 'a likely critical issue in the HS-sewing verification' in thqg_perturbative_finiteness.tex (it had read regions 39-1405 when it found this; the computations it started: (i) Heisenberg normalized structure constants C(p,q,r=p+q) = r! k^r / (sqrt(p! k^p) sqrt(q! k^q) sqrt(r! k^r)) and an exponential-growth counterexample - does the normalized 3-point growth break the claimed Hilbert-Schmidt/sewing bound?; (ii) an F_g table through genus 10 of shape F_g = (2^{2g-1}-1)/2^{2g-1} * |B_{2g}|/(2g)!). Recompute both with sympy, pin the exact claim lines in the chapter, and determine whether the HS-sewing verification survives. Note C(p,q,p+q) = sqrt(binom(p+q,p)) after simplification - if the chapter claims boundedness or summability of these, the binomial growth sqrt(binom(2n,n)) ~ 2^n/n^{1/4} is the counterexample shape.
Then: (a) perturbative finiteness claims - what is precisely proved vs imported folklore; (b) Fredholm determinant identities - verify numerically; (c) thqg_modular_bootstrap: Cardy usage carries modular-invariance + vacuum-dominance? respects the 'not the Maloney-Witten orbit sum' boundary declared in the WIP main.tex? Rademacher/Farey convergence claims vs the actual literature; (d) consistency with 3d_gravity.tex licensing of the same material.`,
  },
  {
    key: 'thqg-b',
    brief: `TOPOLOGICAL-HOLOGRAPHY QG EXTENSIONS, HOLOGRAPHIC GROUP. Files: chapters/connections/thqg_holographic_reconstruction.tex (478-line WIP diff), thqg_gravitational_complexity.tex, thqg_gravitational_s_duality.tex, thqg_gravitational_yangian.tex, twisted_holography_quantum_gravity.tex (412-line WIP diff), thqg_3d_gravity_movements_vi_x.tex.
LEAD FROM A PRIOR PASS (CONFIRMED false entries - your first job is to pin and complete it): the Virasoro shadow table in thqg_holographic_reconstruction.tex region ~1403-1633, against the closed form S_r^cub = (-1)^(r-4) * 6^(r-5) * 240 / (c^(r-3) * (5c+22) * r) at c = 26 for r >= 5: recompute EVERY tabulated row as exact fractions, list each false entry with its line number, determine whether the closed form or the table is correct by checking the upstream derivation (where is S_r^cub derived?), and give the corrected values. Also check the N(N+1)/2 count appearing nearby, and the regions ~1824-2024, ~2356-2446, ~3002-3469 (the prior pass read those).
Then: (a) holographic reconstruction: bulk = Z^der_ch(A) (legitimate) vs dynamical-metric bulk (forbidden slogan 13) at every reconstruction claim; (b) thqg_gravitational_complexity: any ProvedHere about 'complexity' - is complexity given a mathematical definition?; (c) thqg_gravitational_s_duality: precise statement? verify any modular-transformation computation; (d) thqg_gravitational_yangian ~630-760: wedge-vs-full w_{1+infinity} conflation hunt (the celestial soft algebra is the WEDGE; loop corrections deform it); (e) twisted_holography_quantum_gravity: Costello-Gaiotto transplant licensed as analogy or construction?`,
  },
  {
    key: 'thqg-c',
    brief: `TOPOLOGICAL-HOLOGRAPHY QG EXTENSIONS, STRUCTURAL GROUP. Files: chapters/connections/thqg_anomaly_extensions.tex, thqg_bv_construction_extensions.tex, thqg_bv_ht_extensions.tex (216-line WIP diff), thqg_celestial_holography_extensions.tex, thqg_concordance_supplement.tex, thqg_critical_string_dichotomy.tex, thqg_fm_calculus_extensions.tex, thqg_ht_bbl_extensions.tex (247-line WIP diff), thqg_line_operators_extensions.tex, thqg_soft_graviton_theorems.tex, thqg_spectral_braiding_extensions.tex, thqg_symplectic_polarization.tex.
LEADS FROM A PRIOR PASS (verify first): (1) LIVE SIGN CONTRADICTION between thqg_anomaly_extensions.tex ~1139-1214 and anomaly_completed_core.tex ~780-840 on d(eta^2): quote both formulas, derive the correct sign from the stated parity of eta via graded Leibniz d(eta^2) = (d eta) eta + (-1)^{|eta|} eta (d eta), name which file is wrong and the propagation set. The gravitational-algebra definition needed to adjudicate sits at thqg_anomaly_extensions.tex ~737-847. (2) A 'suspicious genus-1 Heisenberg edit' in the thqg_bv_ht_extensions.tex WIP diff (first ~250 diff lines) - adjudicate. (3) The prior pass was checking thqg_symplectic_polarization.tex ~785-965 PTVV shifted-degree bookkeeping - finish that check ((-1)/0/1/2-shifted degrees consistent?).
Then: (a) soft graviton attributions and loop-correction claims in thqg_soft_graviton_theorems.tex (~90-150, ~520-580, ~774-934 were read; HMPS15/CachazoStrominger/BernDaviesNohle citation accuracy; subleading soft receives loop corrections - does the file claim exactness anywhere?); (b) critical string dichotomy c=26 vs Vir complementarity 13+13; (c) anomaly duplication across the four anomaly files (same polynomial, different coefficients = critical); (d) thqg_concordance_supplement label existence sample.`,
  },
  {
    key: 'celestial',
    brief: `CELESTIAL HOLOGRAPHY. Files: chapters/connections/universal_celestial_holography.tex (1167 lines), celestial_holography.tex, celestial_holography_core.tex, celestial_holography_frontier.tex, soft_graviton_mellin_shadow_bridge_platonic.tex (193-line WIP diff), celestial_moonshine_bridge.tex.
Attack: (a) thm:uch-main 'celestial OPE = chiral factorisation homology on P^1_cel': isomorphism of which precisely-defined objects, in which ambient; is the tree-level restriction IN the theorem statement (loop corrections break holomorphy)? (b) The four coverage rows: self-dual gauge -> affine KM; gauge+matter -> DS; gravity -> Virasoro + w_{1+infinity}; YM -> 'Beem-Rastelli chi-functor'. The chi-functor maps 4d N=2 SCFTs to VOAs; plain YM is not N=2 superconformal - category error or licensed construction? (c) wedge-vs-full w_{1+infinity} and classical-vs-quantum (W_{1+infinity}) conflation at every site (the Strominger soft algebra is the wedge; Costello-Paquette form-factor corrections deform the algebra at loop level). (d) Mellin/shadow conventions in the bridge file: verify one explicit transform (principal series Delta = 1 + i lambda; shadow weight 2 - Delta on the 2d celestial sphere; spin handling) by direct computation. (e) conj:uch-gravity-chain-level respected downstream? (f) celestial_moonshine_bridge licensing.`,
  },
  {
    key: 'qg',
    brief: `THE UNIFIED CHIRAL QUANTUM GROUP AND BRAIDING. Files: chapters/theory/unified_chiral_quantum_group.tex (1985 lines), chapters/connections/shifted_rtt_duality_orthogonal_coideals.tex (NEW uncommitted file), chapters/connections/typeA_baxter_rees_theta.tex (173-line WIP diff), chapters/connections/spectral-braiding.tex, spectral-braiding-core.tex (277-line diff), spectral-braiding-frontier.tex.
LEADS FROM A PRIOR PASS (CONFIRMED by computation - pin and complete): (1) a sign error in the shift convention of shifted_rtt_duality_orthogonal_coideals.tex: from the file's own relations [z,x] = hbar x and [z,y] = -hbar y one derives p(z) x = x p(z + hbar) (since z x = x(z + hbar)) - the file states the opposite shift direction somewhere; quote the wrong line, give the corrected convention, and propagate (which downstream formulas inherit the sign). (2) The same file contains the construct Y_{\\hbar}^{\\mathrm{RTT},\\vee}_{\\ge\\mu} - an illegal LaTeX double subscript that cannot compile - yet the June-9 build succeeded: check whether the file is \\input in main.tex at all (grep). If NOT input: the file is a dead new artifact (finding: half-finished work not wired in; assess whether its mathematics deserves completion or quarantine). If input: find how it compiles.
Then: (a) Q_g^{k,f,mu} nine specialisation fibres: enumerate as defined; each actually constructed vs named; completeness argument for 'cover'? (b) typeA_baxter_rees_theta: verify one theta-identity numerically. (c) spectral-braiding: 'stable envelopes = R-matrices (Maulik-Okounkov), Hall envelopes = chiral coproducts' discipline consistent; MO R-matrix as gluing residue carries cocycle/unitarity epsilon-license? (d) CoHA(C^3) = Y^+(gl_1) discipline (forbidden slogan 8) at every CoHA site.`,
  },
  {
    key: 'barcobar',
    brief: `THE BAR-COBAR CORE. Files: chapters/connections/bar-cobar-review.tex (534-line WIP diff - note: part of this WIP was authored by a prior Codex pass propagating Vol I conditional-status corrections, e.g. renaming cor:lagrangian-unconditional to cor:lagrangian-conditional-standard-landscape; audit that propagation's mathematical correctness), chapters/connections/ordered_associative_chiral_kd.tex + _core.tex (thm:opposite at ~451) + _frontier.tex, chapters/connections/brace.tex, appendices/brace-signs.tex, plus the main.tex preamble WIP change (reduced bar B-bar = sum_{n>=1} (s^{-1} A-bar)^{otimes n} inside B = k + B-bar).
Attack: (a) reduced-vs-full bar consistency across all usage sites (coaugmentation conventions; Omega(Bar(A)) = A needs the coaugmented version; silent mixing breaks Theorem A). (b) Theorem A on the chiral lane: does the proof exist as more than 'the classical proof transfers'? The BD chiral pseudo-tensor structure is not naively symmetric monoidal - where is the transfer justified? (c) thm:opposite + R_A sign (-1)^{sum_{i<j}(|a_i|+1)(|a_j|+1)}: verify by hand on length-2,3 words AND run compute/lib/f9_verdier_ordered_engine.py + its tests. (d) brace-signs appendix vs a from-scratch length-2 brace computation. (e) thm:two-color-master proof body. (f) NOTE: compute/tests/test_ordered_chiral_kd_engine.py::test_foundations_source_contains_ordered_bar_skeleton FAILS - the test pins a sentence in the foundations source that the WIP reworded; read the test and the manuscript site, decide which side regressed, and give the reconciliation.`,
  },
  {
    key: 'pva',
    brief: `THE PVA / QUANTIZATION DISCIPLINE. Files: chapters/theory/pva-descent-repaired.tex (1949 lines; input in main.tex:2228) vs chapters/theory/pva-descent.tex (1023 lines; dead, labels stripped with '% label removed' = transplant evidence), chapters/connections/modular_pva_quantization.tex + _core.tex + _frontier.tex, chapters/theory/pva-expanded-repaired.tex (live) vs appendices/pva-expanded.tex (dead).
LEADS FROM A PRIOR PASS (verify first, anchor, extend): (1) CITATION FRAUD CONFIRMED: bibitem RNW19 cites arXiv:1910.12006, which is 'New Geiger-Nuttall law for proton radioactivity' (nuclear physics), not convolution L-infinity algebras; the intended reference is Robert-Nicoud--Wierstra on convolution homotopy Lie algebras - find the bibitem in main.tex (inline thebibliography), give its line, the sites citing RNW19, and the correct reference (Daniel Robert-Nicoud and Felix Wierstra, 'Homotopy morphisms between convolution homotopy Lie algebras', J. Noncommut. Geom. 13 (2019); arXiv:1712.00794 - verify this ID yourself via WebFetch if available). (2) MISQUOTED CONSTANT: Khan-Zeng arXiv:2502.13227 eq. 3.77 has 32/(5c); a site in modular_pva_quantization*.tex quotes a different value (grep '22+5c', '{32}{5c}', '64' sites) - pin the exact site, wrong value, and corrected formula. (3) A 'major sign-consistency problem' in pva-descent-repaired.tex region ~1022-1952, likely around eq:sl2-lambda-brackets (~line 14xx) or the sesquilinearity axioms - find it, verify with sympy (the Virasoro lambda-Jacobi {T_l T} = (d+2l)T + (c/12) l^3 is the control case), and give corrected signs. (4) Semantic-diff the two dead/live pairs for LOST content (the repo forbids cutting content): pva-descent vs -repaired, pva-expanded vs -repaired; also FRONTIER.md and the architecture notes still cite the dead filenames. (5) Check labels thm:raviolo-PVA, prop:AOS-implies-Jacobi, hca-boundary-cech and bibitems PasolZagier2013, KZ25/KhanZeng2025, Mok25, DNP25, MS24 - existence and accuracy.
Then the discipline sweep: classical-only (forbidden slogan 16) at every quantum-flavoured site in the modular trio; the four quantum hypotheses named AT statements.`,
  },
  {
    key: 'examples',
    brief: `THE EXAMPLES LANDSCAPE AND beta_N. Files: chapters/examples/w-algebras*.tex (6 files; w-algebras-w3.tex has a 535-line WIP diff), examples-complete*.tex, examples-computing.tex, examples-worked.tex (395-line diff), rosetta_stone.tex (357-line diff), chapters/theory/beta_N_closed_form_all_platonic.tex, chapters/theory/wn_tempered_closure_platonic.tex.
KNOWN FROM MAIN THREAD (extend, do not re-derive): (i) beta_N = 12(H_N - 1) is integral exactly at N in {2,3,4} (values 6, 10, 13) and non-integral for all 5 <= N <= 200 by direct computation; the all-N >= 5 proof is Bertrand (a prime p in (N/2, N], p coprime to 12, v_p(beta_N) = -1). (ii) thm:beta-N-closed-form-proved-all-N is tagged ClaimStatusConjectured in the manuscript while its LABEL says 'proved-all-N' and FRONTIER.md lists it under closures - a three-way status inconsistency. (iii) compute/tests/test_w4_beta_direct.py FAILS on two text-pins ('They do not prove it.' sentence removed; w4-missing-bridge status pin) - read the test + chapter, decide which side regressed. Your jobs: state precisely what IS proved vs conjectured in the beta_N chapter (the identification of the W_N invariant beta_N with 12(H_N-1) vs the trivial harmonic arithmetic), whether the Bertrand strengthening is absent (if absent it is a missing-artifact finding with the proof supplied in the heal), and the exact status reconciliation.
Then: (a) run the three uncommitted new tests (test_w3_pva.py, test_wn_tempered_closure.py, test_wN_tensor_arakelov_weight_distribution.py) and attack what they pin against the manuscript values. (b) w-algebras-w3.tex WIP diff: attack changed material; verify W_3 structure constants against Zamolodchikov normalisation (C_33^4 ~ 16/(22+5c) structure). (c) rosetta_stone.tex: sample 8 rows vs Vol I landscape_census.tex. (d) test_deletion_ledger.py::test_thqg_faber_pandharipande_scalar_lane_uses_kappa_ch_hodge FAILS - reconcile (which side regressed).`,
  },
  {
    key: 'crossvolume',
    brief: `CROSS-VOLUME ANCHOR VERIFICATION. Files: chapters/connections/concordance.tex (120-line WIP diff), chapters/connections/dnp_identification_master.tex, plus /Users/raeez/chiral-bar-cobar (Vol I) and /Users/raeez/calabi-yau-quantum-groups (Vol III).
KNOWN FROM MAIN THREAD: the kappa IV test asserts c_N(0)/2 = (5, 2, 1, 1, 1/2) at N = (1,2,3,4,6) and c_1(0) = 10 checks against the Gritsenko-Nikulin Delta_5 seed phi_{0,1}; FRONTIER.md's trace-lane spread is w(N) in {5,2,1,1,1/2,1,1/4,0}. Verify the full N-to-weight table against Gritsenko 1999 / Gritsenko-Nikulin (which Delta_t has weight 1? is the N=4 -> 1, N=6 -> 1/2 assignment right, or is it t=3 -> 1, t=4 -> 1/2?) - an off-by-one in the (N,t) indexing is a finding.
Then: (a) Vol III fusion conjecture A_X = A_b 'verified at C^3, local P^2, conifold, K3xE' - find those verifications in Vol III, check each is a real computation, check Vol II quotes scope honestly. (b) The Universal Trace Identity {5,5,5} at N=1: locate the THREE computations (Vol II ghost trace; Vol II Pentagon trace; Vol III Borcherds weight) - genuinely independent routes or one computation cited thrice? (c) Vol I anchors: thm:mc5-class-m-chain-level-pro-ambient exists? Vol I Theorem B statement matches Vol II's P4 gloss? Vol I landscape census kappa conventions match Vol II usage? (d) concordance.tex: sample 12 rows, verify labels exist in the named volumes (grep both repos). (e) Vol III thm:G-eq-D-Yplus (cited by cor:open-cy-lane-unification at sc_chtop_heptagon.tex:1642): exists at theorem grade in Vol III, or still prose (the v2 critique says it NEEDS elevation - check current state)?`,
  },
  {
    key: 'foundations',
    brief: `FOUNDATIONS AND THE TYPED-CARRIER DISCIPLINE. Files: chapters/theory/foundations.tex (5069 lines), chapters/theory/axioms.tex (WIP diff region), chapters/frame/preface.tex (276-line diff).
Attack: (a) def:log-SC-algebra: self-contained, non-circular, strong enough to derive thm:FM-calculus and thm:physics-bridge? Read both proofs; check they use ONLY the definition. (b) Orphaned (H1)-(H4): grep -rn '(H1)' chapters/ appendices/ - known hit files: axioms.tex, concordance.tex, pva-preview.tex, examples-complete.tex, programme_climax_platonic.tex (the last is a LOCAL analytic package defined near line 1918 - distinct from the eliminated standing hypotheses; name collision between packages is itself a finding); classify every hit. (c) def:six-typed-carriers (C1)-(C6) + thm:typed-boundary-holographic-realisation: pairwise distinctions with named comparison morphisms; attack the realisation proof; FRONTIER F15 usage matches the numbering? (d) Five-objects-never-conflated (A_b, Bar, A^i, A^!, Z^der): grep for 'A^! = H' / 'Bar.*Koszul dual' violations in foundations/axioms/preface. (e) preface.tex WIP diff: hedged-in-chapter but unhedged-in-preface claims (the preface is the most collapse-prone surface). (f) test_ordered_chiral_kd_engine.py pins a foundations sentence about the ordered bar skeleton that the WIP broke - locate the site, reconcile.`,
  },
  {
    key: 'structure',
    brief: `MECHANICAL AND STRUCTURAL INTEGRITY. Scope: main.tex (110 inputs, 7 parts), build/test/gate layer, dead files, stale pointers, bibliography integrity.
Attack: (a) Orphan .tex inventory: files in chapters/ + appendices/ not referenced by main.tex; for each, theorems absent from any input file = lost mathematics (critical); known dead/live pairs: chapters/theory/pva-descent.tex (labels stripped) vs pva-descent-repaired.tex; appendices/pva-expanded.tex vs chapters/theory/pva-expanded-repaired.tex; check also preface_trimmed.tex, bar-cobar-review.tex (input or not?), the _core/_frontier pairs vs parents (which of celestial_holography{,_core,_frontier}.tex are input? duplication?), and chapters/connections/shifted_rtt_duality_orthogonal_coideals.tex (NEW file - input or dead?). (b) part_viii_synthesis.tex input at main.tex:2702 inside Part VII: does it read as a closing section or announce itself as 'Part VIII' (defect)? (c) BIBLIOGRAPHY INTEGRITY (inline thebibliography in main.tex): one fabricated bibitem already confirmed (RNW19 -> arXiv:1910.12006 is nuclear physics). Sample 25+ bibitems weighted toward recent arXiv IDs (2024-2026) and unusual venues; verify each against your knowledge (use WebFetch on arxiv.org/abs/<id> if available); list every suspicious entry with its line. (d) make verify-independence: run it, report. (e) The 8 IV files (test_bar_neq_bulk_iv.py etc.): disjoint routes or rubber stamps - read each. (f) metadata/claims.jsonl + theorem_registry.md: sample 10, verify labels. (g) The licensing gate scripts/verify_licensing.sh now counts hyp-tag misses (127 warnings): break down which of the 127 are genuinely-missing-delta-hypotheses vs prose-tagged theorems (sample 15 and classify). (h) test_deletion_ledger.py failure - reconcile pin vs manuscript.`,
  },
]

const VERIFY_PREAMBLE = `You are a skeptical referee on the Vol II monograph repo /Users/raeez/chiral-bar-cobar-vol2 (Vol I at /Users/raeez/chiral-bar-cobar, Vol III at /Users/raeez/calabi-yau-quantum-groups). Below is one finding from an adversarial audit, as JSON. Your job: try to REFUTE it. Re-read every cited anchor (file:line) in the actual files with surrounding context; re-run any cited computation (compute/.venv/bin/python, sympy available); check whether the attacked claim is actually what the manuscript says - misquoting and ignoring an in-context qualifier are the two most common false-finding modes; check whether the proposed heal is itself correct mathematics. Note the tree has uncommitted WIP - check the CURRENT file state. Set isReal=true ONLY if the defect is real as stated, the anchors check out, and the severity is justified. If the finding is directionally right but mis-anchored or wrong in detail, set isReal=true and put the corrected statement in assessment and corrected_heal. If you cannot confirm it after genuinely checking, set isReal=false. Do not edit any file. FINDING: `

phase('Attack')
log('Relaunching 17 attack agents with salvaged leads injected and death-proof journaling')

const results = await pipeline(
  AXES,
  a => agent(PREAMBLE.replace('AXISKEY', a.key) + a.brief, { label: 'attack:' + a.key, phase: 'Attack', schema: FINDINGS_SCHEMA }),
  (res, a) => {
    if (!res) return null
    const nonMinor = res.findings.filter(f => f.severity !== 'minor')
    const toVerify = nonMinor.slice(0, 8)
    if (nonMinor.length > 8) log(a.key + ': verifying top 8 of ' + nonMinor.length + ' non-minor findings; ' + (nonMinor.length - 8) + ' deferred unverified')
    return parallel(toVerify.map(f => () =>
      agent(VERIFY_PREAMBLE + JSON.stringify(f), { label: 'verify:' + a.key + ':' + f.title.slice(0, 40), phase: 'Verify', schema: VERDICT_SCHEMA })
        .then(v => ({ finding: f, referee: v }))
    )).then(checked => ({
      axis: a.key,
      summary: res.summary,
      sound: res.sound,
      minors: res.findings.filter(f => f.severity === 'minor'),
      unverified: nonMinor.slice(8),
      confirmed: checked.filter(Boolean).filter(x => x.referee && x.referee.isReal),
      refuted: checked.filter(Boolean).filter(x => x.referee && !x.referee.isReal),
    }))
  }
)

const out = results.filter(Boolean)
const totals = {
  axes: out.length,
  confirmed: out.reduce((n, r) => n + r.confirmed.length, 0),
  refuted: out.reduce((n, r) => n + r.refuted.length, 0),
  minors: out.reduce((n, r) => n + r.minors.length, 0),
}
log('Attack complete: ' + totals.confirmed + ' confirmed findings, ' + totals.refuted + ' refuted by referees, ' + totals.minors + ' minors across ' + totals.axes + ' axes')
return { totals, results: out }
