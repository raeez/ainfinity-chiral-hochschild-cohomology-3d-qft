# Symmetric Taylor maps on supported quantum currents

Collision002 is a new bounded construction. Candidate001 and its manifest remain unchanged. Its exact quantum contraction PASS is the starting chain comparison. The full ordered interval comparison remains the original target and is not claimed here.

## Question and deciding calculation

The first actual quantum Taylor equation must include the second-order Leibniz defect. For the native smooth observable algebra B with positive D=d+hbar Delta_mu, set

`l2(a,b)=D(ab)-(Da)b-(-1)^|a| a(Db)`.

The Cauchy chain isomorphism N=exp(-hbar K) has the binary coefficient

`F2(a,b)=N(ab)-N(a)N(b)`.

Direct differentiation gives `d F2-F2(Da,b)-(-1)^|a|F2(a,Db)=N l2(a,b)`. This is a higher map for the symmetric BV descendant brackets. It does not define an ordered bar map on an overlapping-support DGA, since that DGA does not exist.

The three-input calculation was derived before generalization and before the literature search. For distinct current representatives x_i, with p_ij=kappa/(a_i-a_j)^2, one gets `F3(x1,x2,x3)=0` but `F3(x1^2,x2,x3)=2 hbar^2 p12 p13`. Expanding the Gaussian exponential gives this coefficient directly. The three two-block terms remove all single contractions. The remaining matching connects all three input groups.

## Complete construction

The source is the actual compact smooth-kernel observable algebra, polynomial in field number and hbar. Functions have degree minus one, forms degree zero. Coinvariants use the graded symmetry signs. Compact support and the residue normalization 1/(2 pi i) remain fixed.

On the finite-word coaugmented symmetric coalgebra S^c(B), let E have multiplication as each Taylor coefficient. Its inverse coefficients are (-1)^(n-1)(n-1)! times multiplication, with Koszul signs. Partition inversion proves that E^{-1} Dhat E has only l1=D and l2 as above. Its square vanishes because D^2=0. The corresponding conjugation of the classical derivation d has only the unary term.

The coalgebra map `F=E^{-1} S^c(N) E` therefore intertwines those coderivations. Its displayed partition formula supplies every Taylor coefficient and its inverse. This proves an actual L-infinity[1] isomorphism, with graded symmetric degree-one source operations and degree-zero Taylor maps. Every identity is on the specified smooth complex. No bulk/centre or interacting equivalence follows.

Every contraction acts on disjoint individual field variables, even when its input groups contain many fields. Thus the construction uses tensor products of well-defined Cauchy contractions. It never multiplies singular distributions in a shared variable. Continuity on fixed compact supports and extension-by-zero naturality follow from the inherited Cauchy estimate. Each polynomial input has finitely many contractions.

## Binary collision domain and exact limit

The distribution `T(z,w)=-partial_z(1/(z-w))` exists because the Cauchy kernel is locally integrable. On compact smooth two-variable kernels it extends the off-diagonal double pole. Relative to area measure,

`partial_bar_z T=-pi partial_z delta^(2)(z-w)`.

The pairing is `P(u,v)=kappa/pi^2 <T,a(w)b(z)>` for u=a(w)dbar w and v=b(z)dbar z. The first formula differentiates the smooth input before integration and therefore works at overlapping supports.

For a compact function f and form v, `F2(dbar f,v)=-hbar mu(f,v)` while `F2(f,v)=0`. This proves the first contact homotopy, with the actual positive quantum differential. The normalization agrees with the inherited negative extracted level.

A limit of point-labelled currents is a different question. For nonzero kappa, `F2(u_a,u_b)=-hbar kappa/(a-b)^2` has no finite limit as b approaches a. Keeping the supports disjoint forces shrinking cutoffs and unbounded smooth seminorms. The proved continuity therefore gives no coincident point-distribution product. This is an exact failed route, not a renormalization prescription.

## Native consumer and source freeze

All mathematical sources are under `research-candidates/vol2_boundary026/collision002/`.

- New full proof: `chapters/current-taylor-maps.tex`, SHA-256 `2891afdc32ff69b413ffd00fd424cbadb1365d9c10bb68c22c4c44db515bd5fe`.
- Native higher-operation consumer: `chapters/Volume_II_Free_Field_Higher_Operations.tex`, SHA-256 `892dda9e1528424391e0770ca54d8aad42e9b1c1052accd22b110a0e9928c80b`.
- Source aggregate: `8659e9bac8bc4ffa80baad4511950320c3bfbf0d6fbb4a5f13da52ece657e2d3`.
- Manifest: `reports/research/vol2_boundary026/collision002/manifest.json`, SHA-256 `6745ca116625d29e1a3e5da90eb2ab555d191cb0b7280c060e021048efc7a9d5`.

The full native consumer includes the new section at the end of the existing higher-operation chapter. The final paragraph states the remaining ordered-coalgebra and collision-chain obligations. The BF and compact-current chapters are byte-identical copies of candidate001. All seven original candidate001 source hashes were rechecked after this construction.

The derivative kernel proof is internal. The general cumulant construction belongs to existing mathematics: the primary literature search found Ruggero Bandiera, *Cumulants, Koszul brackets and homological perturbation theory for commutative BV-infinity and IBL-infinity algebras*, https://arxiv.org/abs/2012.14812. The abstract identifies that subject. No unverified theorem locator from it is used as a proof dependency, and no novelty claim is made.

## Checks

`verify_taylor.py` checks 81 binary and 729 ternary morphism identities over exact rational arithmetic, including odd inputs. It independently encodes the new cumulant formulas and symmetric insertion signs. It reuses the polynomial arithmetic helper definitions from the previous calculation, with the helper hash recorded. It does not rerun the accepted contraction tests. The three-current examples and contact cancellation pass. Python is version 3.14.6.

The complete internal source builds to 37 pages. Consecutive deterministic pdfLaTeX passes have identical PDF and AUX hashes. There are zero undefined references, overfull boxes, underfull boxes, or duplicate labels. The remaining epstopdf warning states that shell escape is disabled. All 325 recorder inputs are archived and hashed. Physical pages 35–37, covering all new prose and its transition, were inspected as PNGs. No layout or firewall violation was found. No standalone PDF was opened.

The first build failed at the legacy `cal` font command in new notation. It was replaced by the standard `mathcal` command before the final freeze. No template or typography workaround was introduced.

## Scope, custody, and remaining work

This construction proves the symmetric Taylor comparison and the binary distributional contact homotopy on smooth test kernels. It does not complete the requested ordered A-infinity interval comparison, construct three-point diagonal extension domains, or define coincident point-supported products. It does not solve an interacting quantum master equation.

The coordinator approved two independent next constructions. Their active handles are `/root/dispatch_resume026/vol2_boundary026/ordered_taylor026` and `/root/dispatch_resume026/vol2_boundary026/collision_domain026`. They own separate prepared worktrees. The first owns ordered F2/F3 and the second owns the genuine three-point distributional domain. Their results are not incorporated into this frozen source or certified by this report.

Fresh independent review of collision002 remains required. The explicit runtime requirement is gpt-6-astra with ultra. Independently observed controls remain unavailable and unverified. No staging, commit, push, publication, central reader write, or modification of candidate001 occurred.
